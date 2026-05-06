"""
bacTRAP-to-HypoMap Mapping Tool

Streamlit application for mapping bacTRAP bulk RNA-seq data onto the
murine HypoMap single-cell atlas (Steuernagel et al., Nature Metabolism 2022).

Run with: streamlit run app.py
"""

import io
import logging
import re
import zipfile

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ---------------------------------------------------------------------------
# Logging setup — writes to bactrap_hypomap.log alongside app.py
# ---------------------------------------------------------------------------
_LOG_FILE = Path(__file__).parent / "bactrap_hypomap.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    handlers=[
        logging.FileHandler(_LOG_FILE, mode="a"),
        logging.StreamHandler(),          # also print to terminal
    ],
    force=True,
)
logger = logging.getLogger(__name__)
# Suppress noisy third-party loggers
logging.getLogger("fontTools").setLevel(logging.WARNING)
logging.getLogger("matplotlib").setLevel(logging.WARNING)
logger.info("="*60)
logger.info("App startup / Streamlit rerun")

st.set_page_config(
    page_title="bacTRAP → HypoMap Mapping",
    page_icon="🧬",
    layout="wide",
)

# ---------------------------------------------------------------------------
# Sidebar: file inputs and parameters
# ---------------------------------------------------------------------------

st.sidebar.title("bacTRAP → HypoMap Mapping")
st.sidebar.markdown("---")

st.sidebar.subheader("Input Files")

bactrap_file = st.sidebar.text_input(
    "bacTRAP FPKM table (.xlsx)",
    value="",
    placeholder="/path/to/IPvsInput_deg.xlsx",
    help="Path to the bacTRAP DESeq2 results Excel file.",
)

hypomap_file = st.sidebar.text_input(
    "HypoMap atlas (.h5ad)",
    value="",
    placeholder="/path/to/hypomap.h5ad",
    help="Path to the HypoMap AnnData h5ad file.",
)

st.sidebar.markdown("---")
st.sidebar.subheader("Parameters")

padj_cutoff = st.sidebar.slider(
    "Adjusted p-value cutoff", 0.001, 0.1, 0.05, 0.001, format="%.3f",
    help=(
        "Benjamini–Hochberg-adjusted p-value threshold for a bacTRAP gene to "
        "count as 'enriched'. Works in concert with the log₂FC cutoff — "
        "both must be satisfied. The standard DESeq2 threshold of 0.05 is "
        "a reasonable default; lower it (e.g. 0.01) for a more stringent "
        "signature, raise it toward 0.1 if your bacTRAP is underpowered."
    ),
)
log2fc_cutoff = st.sidebar.slider(
    "log₂FC cutoff", 0.0, 5.0, 1.0, 0.25,
    help=(
        "Minimum log₂ fold change (IP vs Input) for a bacTRAP gene to "
        "count as 'enriched'. 1.0 (2-fold) is the convention in the "
        "translational-profiling literature. Pair with the Min IP "
        "expression filter below — raising log₂FC alone admits "
        "pseudocount artefacts from genes with near-zero Input."
    ),
)
min_ip_expression = st.sidebar.slider(
    "Min IP expression (baseMean)", 0.0, 200.0, 10.0, 5.0,
    help=(
        "Minimum mean IP expression (DESeq2 baseMean-style, from the 'IP' "
        "column) required for a gene to pass the enrichment filter. "
        "DESeq2 inflates log₂FC for low-count genes with near-zero Input "
        "(e.g. 28 IP reads vs 0 Input → log₂FC≈7), so the default ranking "
        "otherwise pulls these artefacts above real Cre-driver signal. "
        "Set to 0 to disable. A typical value of 10 removes zero-Input "
        "low-count noise without dropping genuinely cell-type-specific genes."
    ),
)
ranking_metric_label = st.sidebar.selectbox(
    "Gene ranking metric",
    ["π-score (|log₂FC| × -log₁₀padj)", "log₂FoldChange", "-log₁₀(padj)"],
    index=0,
    help=(
        "How the app picks the top-N genes for AUCell, dot-plot, and "
        "heat-map. π-score (Xiao et al. 2014) balances effect size with "
        "significance so a strongly significant moderate-FC gene (e.g. "
        "Pnoc, FC=2, padj=1e-17) outranks a low-count pseudocount artefact "
        "(FC=7, padj=1e-4). Choose 'log₂FoldChange' for the historical "
        "behaviour."
    ),
)
_ranking_metric_map = {
    "π-score (|log₂FC| × -log₁₀padj)": "pi_score",
    "log₂FoldChange": "log2fc",
    "-log₁₀(padj)": "padj",
}
ranking_metric = _ranking_metric_map[ranking_metric_label]
top_n_genes = st.sidebar.slider(
    "Top N genes for scoring", 10, 500, 50, 10,
    help=(
        "Size of the bacTRAP signature passed to AUCell, the dot plot, and "
        "the heatmap (after the padj / log₂FC / min-IP-expression filters, "
        "ranked by the metric above). 50 is a good balance for rank-based "
        "scoring — small enough to stay inside the AUCell 5 %% top-ranked "
        "window (≈ 1,500 genes on HypoMap), large enough to be robust to "
        "dropout in any individual cell. If this exceeds "
        "~τ × total-genes, the AUCell window is automatically widened and "
        "the tab will warn you (see Methods, *AUCell scoring*)."
    ),
)
aucell_top_fraction = st.sidebar.slider(
    "AUCell top-ranked fraction", 0.01, 0.20, 0.05, 0.01, format="%.2f",
    help=(
        "Fraction of genes (ranked by expression within each cell) in which "
        "the bacTRAP signature must appear to contribute to the AUCell score. "
        "Smaller = more stringent: at 0.01 only cells where the signature "
        "concentrates in the top 1% of expressed genes score highly. "
        "AUCell default is 0.05 (Aibar et al. 2017). Use 0.01–0.03 for a "
        "more conservative call on bacTRAP-target identity."
    ),
)
n_markers_per_cluster = st.sidebar.slider(
    "Marker genes per cluster", 20, 500, 100, 10,
    help=(
        "Number of top marker genes retrieved per HypoMap cluster (Wilcoxon "
        "or t-test ranking via `scanpy.tl.rank_genes_groups`). These marker "
        "sets feed Fisher's exact overlap test (Suppl. S5) and preranked "
        "GSEA (Suppl. S9/S10). 100 is a common default; smaller sets "
        "emphasise the very top markers per cluster, larger sets give "
        "Fisher's test more statistical power at the cost of specificity. "
        "Does not affect AUCell."
    ),
)
min_cells_per_cluster = st.sidebar.slider(
    "Min cells per cluster (markers)", 1, 100, 10, 1,
    help=(
        "Minimum cell count required for a cluster to enter the marker / "
        "correlation / Fisher / GSEA pipelines. Small clusters yield "
        "unreliable Wilcoxon ranks. Distinct from 'Min cells for AUCell "
        "top-N ranking' below, which only gates AUCell figure rankings."
    ),
)
min_cells_for_rank = st.sidebar.slider(
    "Min cells for AUCell top-N ranking", 1, 200, 20, 1,
    help=(
        "Clusters with fewer than this many cells are excluded from the "
        "top-ranked set shown in main figures 1b/1c and supplementary S2. "
        "Mean AUCell for very small clusters is dominated by shrinkage "
        "variance, so a "
        "2-cell cluster with a slightly above-average mean can otherwise "
        "claim a top slot purely by chance (fix #4 from the critical "
        "evaluation). Raw per-cluster CSV is unaffected and still contains "
        "every cluster."
    ),
)
_marker_method_label = st.sidebar.selectbox(
    "Marker test",
    ["Wilcoxon (robust, slower)", "t-test overestim_var (faster)"],
    index=0,
    help=(
        "Wilcoxon is non-parametric and recommended for publication-grade "
        "p-values, but takes ~15 min on HypoMap-scale data. t-test "
        "overestim_var produces a comparable ranking in a fraction of the "
        "time — fine when downstream steps use the ranking, not the "
        "nominal p-values."
    ),
)
marker_method = "wilcoxon" if _marker_method_label.startswith("Wilcoxon") else "t-test_overestim_var"
umap_subsample = st.sidebar.slider(
    "UMAP subsample (cells)", 10000, 200000, 50000, 5000,
    help=(
        "How many cells to render in the UMAP panels (Figures 1a/1b, "
        "Suppl. S4). Subsampling only affects rendering speed and PDF "
        "file size — all cells are used for AUCell scoring and every "
        "statistical test. 50 k is a readable balance for HypoMap's "
        "~385 k cells; drop below 20 k for faster previews, raise above "
        "100 k if you need rare clusters to survive the random sample."
    ),
)
hide_unassigned = st.sidebar.checkbox(
    "Hide Unassigned / Mixed clusters in rankings",
    value=False,
    help=(
        "Filter clusters whose label contains 'Unassigned' or 'Mixed' from "
        "ranking tables and selection lists.  Useful because HypoMap's "
        "uncurated mixed clusters can artificially top GSEA and other "
        "rankings due to heterogeneous membership.  Underlying computations "
        "still include all clusters; only the displayed rankings are filtered."
    ),
)

st.sidebar.markdown("---")
st.sidebar.subheader("Cre-driver Sanity Check")
sanity_gene = st.sidebar.text_input(
    "Cre-driver gene",
    value="Pnoc",
    help=(
        "Marker gene used as a confidence check on the mapping — typically "
        "the Cre-driver of the bacTRAP line (e.g. Pnoc for Pnoc-Cre;NuTRAP). "
        "Top-ranked clusters that also express this gene are high-confidence; "
        "those that don't may reflect developmental Cre lineage tracing, "
        "snRNA-seq dropout, or background. Note: snRNA-seq dropout for "
        "neuropeptides means 'not detected' ≠ 'not expressed'."
    ),
)
sanity_fraction_threshold = st.sidebar.slider(
    "Expression threshold (fraction)",
    0.0, 0.5, 0.05, 0.01, format="%.2f",
    help=(
        "**Display-only.** Minimum fraction of cells expressing the "
        "Cre-driver gene for a cluster to be flagged as 'expressing' in "
        "the Cre-driver Check sanity table. Does NOT change any ranking "
        "— use the baseline slider below for that. 5% is a permissive "
        "default that tolerates dropout."
    ),
)
sanity_baseline_mean_expr = st.sidebar.slider(
    "Baseline Cre-driver mean expression (log-norm)",
    0.0, 1.5, 0.0, 0.01, format="%.2f",
    help=(
        "**Ranking filter.** When > 0, drop clusters whose mean "
        "(log-normalized) Cre-driver expression falls below this floor "
        "from every cluster-level ranking: correlation / Fisher / NNLS / "
        "GSEA / composite consensus (survivors are re-ranked against one "
        "another), AUCell cluster figures (1b, S2, 1c), and the heatmap "
        "(S5, with z-scores recomputed against the filtered reference). "
        "Per-cell panels (AUCell UMAP fig 1a, per-cell CSVs) are "
        "unaffected — they carry no cluster identity. Filtered CSV "
        "downloads are suffixed with the filter signature "
        "(e.g. `composite_ranking_pnoc_ge0p05.csv`). Set to 0 to disable. "
        "Caveat: snRNA-seq dropout for neuropeptides means 'not detected' "
        "≠ 'not expressed' — a strict floor can discard genuine positives. "
        "Start at ~0.05 and inspect the sanity-check table to tune."
    ),
)

st.sidebar.markdown("---")
st.sidebar.subheader("Figure Settings")
fig_width_mode = st.sidebar.radio(
    "Figure width", ["Single column (89mm)", "Double column (183mm)"],
    index=0,
    help=(
        "Target print width for the exported PDF / SVG figures. Nature's "
        "column widths are 89 mm (single) and 183 mm (double). Choose "
        "single for individual panels that will be placed in a one-column "
        "slot; double for figures that will span the page. Some functions "
        "force a layout regardless (e.g. two-panel UMAPs always render at "
        "double-column width)."
    ),
)
double_column = "Double" in fig_width_mode

# Annotation column selector — populated after data load
annotation_col_key = "annotation_col"

st.sidebar.markdown("---")
run_button = st.sidebar.button("Run Analysis", type="primary", use_container_width=True)

# ---------------------------------------------------------------------------
# Session state initialization
# ---------------------------------------------------------------------------

if "analysis_done" not in st.session_state:
    st.session_state.analysis_done = False

# ---------------------------------------------------------------------------
# Main area
# ---------------------------------------------------------------------------

st.title("bacTRAP-to-HypoMap Mapping Tool")
st.caption(
    "Map preoptic area (PoA) bacTRAP translational profiling data onto the "
    "murine HypoMap single-cell atlas to identify matching cell populations."
)

# Check file inputs — must be existing files, not directories
bt_path = Path(bactrap_file.strip()) if bactrap_file.strip() else None
hm_path = Path(hypomap_file.strip()) if hypomap_file.strip() else None

files_ready = (
    bt_path is not None
    and hm_path is not None
    and bt_path.is_file()
    and hm_path.is_file()
)

if not files_ready and (bt_path or hm_path):
    if bt_path and not bt_path.exists():
        st.warning(f"bacTRAP file not found: `{bt_path}`")
    elif bt_path and bt_path.is_dir():
        st.warning(f"bacTRAP path is a directory, not a file: `{bt_path}`")
    if hm_path and not hm_path.exists():
        st.warning(f"HypoMap file not found: `{hm_path}`")
    elif hm_path and hm_path.is_dir():
        st.warning(f"HypoMap path is a directory, not a file: `{hm_path}`")

if not files_ready:
    st.info("Enter file paths in the sidebar and click **Run Analysis** to begin.")
    st.stop()

# ---------------------------------------------------------------------------
# Data Loading (always loads when files are ready, uses caching)
# ---------------------------------------------------------------------------

from data_loading import (
    load_hypomap,
    load_bactrap,
    get_annotation_columns,
    match_genes,
    compute_cluster_mean_expression,
    compute_fraction_expressing,
    compute_single_gene_cluster_stats,
    get_gene_names_from_adata,
    _detect_gene_column,
    _build_adata_gene_lookup,
    _looks_like_ensembl,
)
from analysis import (
    compute_enrichment_correlation,
    get_enriched_genes,
    rank_enriched_genes,
    compute_marker_genes,
    load_precomputed_markers,
    fisher_overlap_test,
    compute_zscore_heatmap_data,
    compute_nnls_deconvolution,
    compute_gsea_enrichment,
    compute_aucell_scores,
    validate_aucell_input,
    compute_cluster_enrichment_stats,
    compute_composite_ranking,
)
from figures import (
    setup_nature_style,
    figure_bactrap_volcano,
    figure_correlation_barplot,
    figure_umap_enrichment,
    figure_dotplot,
    figure_volcano_enrichment,
    figure_heatmap,
    figure_nnls_barplot,
    figure_gsea_curves,
    figure_gsea_barplot,
    figure_aucell_umap,
    figure_celltype_umap,
    figure_aucell_cluster_barplot,
    figure_aucell_violins,
    figure_aucell_histogram,
    figure_composite_ranking,
    figure_marker_gene_diagnostic,
    fig_to_bytes,
)

# Load data with caching
with st.spinner("Loading data..."):
    bactrap_df = load_bactrap(bactrap_file.strip())
    adata = load_hypomap(hypomap_file.strip())

# Early UMAP sanity check — fail now rather than after minutes of analysis
if not any(key in adata.obsm for key in ("X_umap", "X_UMAP")):
    st.error(
        "No UMAP coordinates found in HypoMap `.obsm` (expected `X_umap`). "
        "UMAP projection is required for Figures 1a and S2. "
        "Please provide an atlas that includes precomputed UMAP coordinates."
    )
    st.stop()

# Early sanity check on bacTRAP required columns
for _req_col in ("padj", "log2FoldChange"):
    if _req_col not in bactrap_df.columns:
        st.error(
            f"Required column **`{_req_col}`** not found in the bacTRAP file. "
            f"Available columns: `{list(bactrap_df.columns)}`. "
            f"Expected DESeq2-style output with `padj` and `log2FoldChange` columns."
        )
        st.stop()

# Annotation column selection
ann_cols = get_annotation_columns(adata)
if len(ann_cols) == 0:
    st.error("No suitable annotation columns found in HypoMap .obs.")
    st.stop()

default_idx = 0
# Prefer C185_named (best resolution for bacTRAP mapping), then C66_named
for preferred in ["C185_named", "C66_named"]:
    if preferred in ann_cols:
        default_idx = ann_cols.index(preferred)
        break
else:
    # Fallback: look for any cell_type or cluster column
    for i, col in enumerate(ann_cols):
        col_lower = col.lower()
        if "cell_type" in col_lower or "celltype" in col_lower or "cluster" in col_lower:
            default_idx = i
            break

annotation_col = st.sidebar.selectbox(
    "Annotation column",
    ann_cols,
    index=default_idx,
    key=annotation_col_key,
    help=(
        "HypoMap cell-type annotation level used for every cluster-level "
        "analysis and figure (`C7_named` … `C465_named`). Finer levels "
        "(higher numbers) give more granular clusters but smaller cell "
        "counts per cluster — which may push populations below the "
        "'Min cells for AUCell top-N ranking' threshold. `C185_named` is "
        "the HypoMap default and a reasonable starting point for "
        "hypothalamic neuropeptide signatures."
    ),
)

# Gene column selection for bacTRAP data
# Build HypoMap lookup once for auto-detection
_adata_lookup, _adata_gnames, _adata_has_raw = _build_adata_gene_lookup(adata)
auto_gene_col = _detect_gene_column(bactrap_df)

# Build candidate list: auto-detected first, then other string/object columns + index
_gene_col_candidates = []
if auto_gene_col != "_index" and auto_gene_col in bactrap_df.columns:
    _gene_col_candidates.append(auto_gene_col)
for col in bactrap_df.columns:
    if col not in _gene_col_candidates:
        # Include string-like columns and columns with Ensembl-looking values
        dtype = bactrap_df[col].dtype
        if dtype == object or dtype.name in ("string", "category"):
            _gene_col_candidates.append(col)
_gene_col_candidates.append("(use row index)")

_default_gene_idx = 0
# Try auto-selecting the column with the best match count
_best_matches = 0
for i, col in enumerate(_gene_col_candidates):
    if col == "(use row index)":
        _vals = bactrap_df.index.astype(str)
    else:
        _vals = bactrap_df[col].astype(str)
    _n = sum(1 for v in _vals if str(v).strip().lower() in _adata_lookup)
    if _n > _best_matches:
        _best_matches = _n
        _default_gene_idx = i

gene_col_selection = st.sidebar.selectbox(
    "bacTRAP gene column",
    _gene_col_candidates,
    index=_default_gene_idx,
    help=(
        "Column in the bacTRAP file containing gene identifiers. "
        "Auto-detected by testing each candidate column against the "
        "HypoMap symbol + Ensembl lookup and picking the highest match "
        "rate. Both gene symbols and Ensembl IDs work. `(use row index)` "
        "lets you use the DataFrame index when gene IDs live there. "
        "Override only if the auto-pick looks wrong — the Gene Matching "
        "Diagnostics section of the Data Overview tab shows the per-"
        "column match rate that drove the auto-selection."
    ),
)
# Map UI selection to the internal value expected by match_genes
_gene_col_for_matching = "_index" if gene_col_selection == "(use row index)" else gene_col_selection

# Early validation: catch a stale / invalid gene-column selection here, not
# 100 lines later deep inside match_genes. Columns can disappear between
# the initial widget render and a rerun (e.g. after the user edits the
# bacTRAP file). The `_index` sentinel is always valid because every
# DataFrame has an index.
if (
    _gene_col_for_matching != "_index"
    and _gene_col_for_matching not in bactrap_df.columns
):
    st.error(
        f"Selected bacTRAP gene column **`{_gene_col_for_matching}`** is "
        f"not present in the bacTRAP file. Available columns: "
        f"`{list(bactrap_df.columns)}`. Pick a different column in the "
        f"sidebar."
    )
    st.stop()

# ---------------------------------------------------------------------------
# Run analysis
# ---------------------------------------------------------------------------

# Progress bar placeholder — rendered above tabs so it's always visible
progress_placeholder = st.empty()

tab1, tab_aucell, tab_sanity, tab3, tab4, tab5, tab6, tab7, tab8, tab_export = st.tabs([
    "📊 Data Overview",
    "⭐ AUCell (Main Figure)",
    "🔍 Cre-driver Check",
    "📈 Correlation (Suppl.)",
    "🗺️ UMAP Projection (Suppl.)",
    "🔬 Marker Overlap (Suppl.)",
    "🔥 Heatmap (Suppl.)",
    "⚖️ NNLS (Suppl.)",
    "📶 GSEA (Suppl.)",
    "📥 Export",
])

# Build a fingerprint of all analysis parameters so we can skip recomputation
# on Streamlit reruns when nothing changed.
_analysis_params = (
    bactrap_file.strip(), hypomap_file.strip(),
    _gene_col_for_matching, annotation_col,
    padj_cutoff, log2fc_cutoff, min_ip_expression, ranking_metric,
    top_n_genes, aucell_top_fraction,
    n_markers_per_cluster, min_cells_per_cluster, marker_method,
    umap_subsample,
)

if run_button or st.session_state.analysis_done:

    # Reuse cached results on rerun if parameters haven't changed
    _cached = st.session_state.get("_analysis_cache")
    _need_recompute = run_button or _cached is None or _cached.get("params") != _analysis_params

    if _need_recompute:
        # ---- Log user inputs (record-keeping) ----
        logger.info("-" * 60)
        logger.info("User inputs for this run:")
        logger.info("  bacTRAP file: %s", bactrap_file.strip())
        logger.info("  HypoMap file: %s", hypomap_file.strip())
        logger.info("  gene column: %s", _gene_col_for_matching)
        logger.info("  annotation column: %s", annotation_col)
        logger.info("  padj cutoff: %.3f", padj_cutoff)
        logger.info("  log2FC cutoff: %.2f", log2fc_cutoff)
        logger.info("  min IP expression: %.1f", min_ip_expression)
        logger.info("  ranking metric: %s", ranking_metric)
        logger.info("  top N genes: %d", top_n_genes)
        logger.info("  AUCell top fraction: %.2f", aucell_top_fraction)
        logger.info("  marker genes per cluster: %d", n_markers_per_cluster)
        logger.info("  min cells per cluster (markers): %d", min_cells_per_cluster)
        logger.info("  min cells for AUCell top-N ranking: %d", min_cells_for_rank)
        logger.info("  marker method: %s", marker_method)
        logger.info("  UMAP subsample: %d", umap_subsample)
        logger.info("  hide Unassigned/Mixed: %s", hide_unassigned)
        logger.info("  Cre-driver gene: %s", sanity_gene)
        logger.info("  Cre-driver expression fraction threshold: %.2f", sanity_fraction_threshold)
        logger.info("  Cre-driver baseline mean expression (log-norm): %.2f", sanity_baseline_mean_expr)
        logger.info("  figure width mode: %s", fig_width_mode)
        logger.info("-" * 60)

        # ---- Gene matching ----
        progress = progress_placeholder.progress(0, text="Matching genes...")

        bactrap_matched, matched_genes, gene_to_idx, matched_in_raw = match_genes(
            bactrap_df, adata, gene_col=_gene_col_for_matching,
            _prebuilt_lookup=(_adata_lookup, _adata_gnames, _adata_has_raw),
        )

        if len(matched_genes) == 0:
            progress.empty()
            st.error(
                f"**No genes matched.** The selected gene column "
                f"**`{_gene_col_for_matching}`** produced zero matches against "
                f"the HypoMap atlas lookup (size ≈ {len(_adata_lookup):,}).\n\n"
                f"**How to fix:**\n"
                f"1. Open the **Data Overview** tab and check **Gene Matching "
                f"Diagnostics** — it shows sample gene IDs from both sides so "
                f"you can see whether you're passing symbols where Ensembl IDs "
                f"are expected (or vice versa).\n"
                f"2. Override **bacTRAP gene column** in the sidebar and try "
                f"the column with the highest match rate.\n"
                f"3. If no column works, the atlas may be in a different "
                f"species or use an unusual symbol convention — check the "
                f"sample atlas keys logged to `bactrap_hypomap.log`."
            )
            st.stop()

        progress.progress(10, text="Genes matched. Computing cluster means...")

        # ---- Cluster mean expression ----
        # Explicit normalize=True: the raw layer holds counts, and downstream
        # correlation / NNLS require log-normalized means in comparable space.
        # Leaving it on auto is unsafe when the gene subset is sparse (the
        # heuristic samples only the selected genes and can falsely decide
        # the data is already normalized — see the second call below).
        gene_indices = [gene_to_idx[g] for g in matched_genes]
        cluster_mean_expr = compute_cluster_mean_expression(
            adata, gene_indices, annotation_col, min_cells=min_cells_per_cluster,
            indices_in_raw=matched_in_raw,
            normalize=True,
        )
        progress.progress(25, text="Cluster means computed. Identifying enriched genes...")

        # ---- Enriched genes ----
        # Pass the min-IP-expression filter through; it guards against the
        # DESeq2 low-count log₂FC-inflation artefact (see rank_enriched_genes
        # docstring for details).
        enriched_df = get_enriched_genes(
            bactrap_matched, padj_cutoff, log2fc_cutoff,
            min_ip_expression=min_ip_expression, ip_col="IP",
        )
        enriched_genes_list = enriched_df["_hypomap_gene_name"].tolist()

        # Rank enriched genes by the user-chosen metric (default π-score,
        # which balances effect size with statistical significance so
        # low-count / zero-Input pseudocount artefacts don't dominate).
        enriched_sorted = rank_enriched_genes(enriched_df, metric=ranking_metric)
        top_enriched_genes = enriched_sorted["_hypomap_gene_name"].tolist()[:top_n_genes]

        progress.progress(30, text="Computing enrichment correlation...")

        # ---- Correlation analysis ----
        # Use only enriched genes (positive FC + significant padj) — negative FC
        # genes are from non-target cell populations and would dilute the signal.
        corr_df = compute_enrichment_correlation(enriched_df, cluster_mean_expr)
        progress.progress(30, text="Computing marker gene overlap...")

        # ---- Marker gene overlap ----
        # Markers only depend on the atlas file, annotation column, n_genes,
        # min_cells, and the chosen test — NOT on padj/log2FC/top_n/umap
        # subsample.  They're cached separately from the full analysis
        # fingerprint so twiddling an unrelated slider doesn't trigger the
        # ~15-min Wilcoxon recompute.  Pre-computed markers from adata.uns
        # are reused when the cluster labels match the selected annotation.
        _markers_params = (
            hypomap_file.strip(), annotation_col,
            n_markers_per_cluster, min_cells_per_cluster, marker_method,
        )
        _markers_cached = st.session_state.get("_markers_cache")
        if _markers_cached is not None and _markers_cached.get("params") == _markers_params:
            markers = _markers_cached["markers"]
            logger.info("reusing cached markers (params unchanged)")
        else:
            markers = load_precomputed_markers(adata)
            if markers is not None:
                current_clusters = set(adata.obs[annotation_col].unique().astype(str))
                marker_clusters = set(markers.keys())
                overlap_ratio = len(current_clusters & marker_clusters) / max(len(current_clusters), 1)
                if overlap_ratio < 0.5:
                    markers = None  # mismatch — recompute for the selected annotation
            if markers is None:
                with st.spinner("Computing marker genes (this may take several minutes)..."):
                    try:
                        markers = compute_marker_genes(
                            adata, annotation_col,
                            n_genes=n_markers_per_cluster,
                            min_cells=min_cells_per_cluster,
                            method=marker_method,
                        )
                    except Exception as e:
                        # Common on pathologically small atlases or when
                        # min_cells_per_cluster is too aggressive for the
                        # selected annotation level. Fisher / GSEA / dotplot
                        # will all degrade gracefully on an empty marker dict.
                        logger.exception("Marker gene computation failed")
                        st.warning(
                            f"Marker gene computation failed: {e}. "
                            f"Fisher overlap, GSEA, dotplot and heatmap panels "
                            f"will be unavailable. Try lowering "
                            f"'Min cells per cluster (markers)' or selecting "
                            f"a coarser annotation level."
                        )
                        markers = {}
            st.session_state["_markers_cache"] = {
                "params": _markers_params, "markers": markers,
            }
        progress.progress(45, text="Running Fisher's exact test...")

        # Universe = bacTRAP-matched atlas genes. Pass the explicit set so
        # cluster markers (drawn from the full atlas) are intersected with
        # it before the 2×2 contingency table is built — otherwise
        # n_markers can exceed universe_size and the "neither" cell goes
        # negative, breaking Fisher's exact (audit fix #3).
        universe_size = len(matched_genes)
        fisher_df = fisher_overlap_test(
            enriched_genes_list, markers, universe_size,
            gene_universe=matched_genes,
        )
        progress.progress(50, text="Computing UMAP enrichment scores...")

        # ---- UMAP enrichment score ----
        if len(top_enriched_genes) == 0:
            st.warning(
                f"No genes pass enrichment thresholds (padj < {padj_cutoff}, "
                f"log₂FC > {log2fc_cutoff}). AUCell and all downstream "
                "scores will be zero — try relaxing the cutoffs."
            )
        progress.progress(55, text="Preparing figures...")

        # ---- Fraction expressing for dotplot ----
        enriched_gene_indices = [gene_to_idx[g] for g in top_enriched_genes if g in gene_to_idx]
        n_dropped = len(top_enriched_genes) - len(enriched_gene_indices)
        if n_dropped > 0:
            st.warning(f"{n_dropped} enriched gene(s) could not be mapped back to HypoMap indices and were excluded from the dot plot.")
        # Prefer clusters significant by both Pearson and Spearman for downstream
        # displays; fall back to all clusters if too few pass the dual filter.
        if len(corr_df) > 0 and "both_significant" in corr_df.columns:
            sig_clusters = corr_df.loc[corr_df["both_significant"], "cluster"].tolist()
            top_clusters_corr = sig_clusters[:15] if len(sig_clusters) >= 5 else corr_df["cluster"].tolist()[:15]
        else:
            top_clusters_corr = corr_df["cluster"].tolist()[:15] if len(corr_df) > 0 else []

        frac_expr = compute_fraction_expressing(
            adata, enriched_gene_indices, annotation_col,
            min_cells=min_cells_per_cluster,
            indices_in_raw=matched_in_raw,
        )
        # Explicit normalize=True — top_enriched_genes are typically sparse,
        # cell-type-specific markers (e.g. neuropeptides), so the auto-detect
        # heuristic sees a low sample max and skips normalisation, leaving
        # the dotplot in raw-count space while correlation/NNLS use log-norm.
        enriched_mean_expr = compute_cluster_mean_expression(
            adata, enriched_gene_indices, annotation_col,
            min_cells=min_cells_per_cluster,
            indices_in_raw=matched_in_raw,
            normalize=True,
        )

        # ---- Z-score heatmap data ----
        if len(corr_df) > 0 and "both_significant" in corr_df.columns:
            sig_clusters_hm = corr_df.loc[corr_df["both_significant"], "cluster"].tolist()
            top_clusters_heatmap = sig_clusters_hm[:20] if len(sig_clusters_hm) >= 5 else corr_df["cluster"].tolist()[:20]
        else:
            top_clusters_heatmap = corr_df["cluster"].tolist()[:20] if len(corr_df) > 0 else []
        top_genes_heatmap = top_enriched_genes[:30]
        zscore_df = compute_zscore_heatmap_data(
            cluster_mean_expr, top_genes_heatmap, top_clusters_heatmap,
        )

        progress.progress(60, text="Running NNLS deconvolution...")

        # ---- NNLS deconvolution ----
        # Use only enriched genes — negative FC genes from non-target populations
        # would distort the deconvolution by fitting against unwanted signal.
        nnls_df = compute_nnls_deconvolution(enriched_df, cluster_mean_expr)

        progress.progress(65, text="Running GSEA enrichment...")

        # ---- GSEA enrichment ----
        try:
            gsea_df, gsea_running_scores, gsea_ranked_genes = compute_gsea_enrichment(
                bactrap_matched, markers, n_perm=1000,
            )
        except Exception as e:
            # Full traceback to the log file so the operator can diagnose
            # post-hoc; UI gets a short message.
            logger.exception("GSEA computation failed")
            st.warning(
                f"GSEA computation failed: {e}. See `bactrap_hypomap.log` "
                f"for the full traceback. Analysis continues without GSEA."
            )
            gsea_df = pd.DataFrame()
            gsea_running_scores = {}
            gsea_ranked_genes = np.array([])

        progress.progress(78, text="Validating AUCell input layer...")

        # ---- AUCell input-layer validation (fix #2: guard against
        # non-raw-count layers silently being fed into AUCell) ----
        aucell_qc = validate_aucell_input(adata, use_raw=True)

        progress.progress(80, text="Computing AUCell scores...")

        # ---- AUCell scoring ----
        aucell_run_info: dict = {}
        aucell_scores = compute_aucell_scores(
            adata, top_enriched_genes, top_fraction=aucell_top_fraction,
            info_out=aucell_run_info,
            prebuilt_lookup=(_adata_lookup, _adata_gnames, _adata_has_raw),
        )

        # Merge scoring diagnostics into the QC report so the UI surfaces
        # match rate + n_top bumps (fixes #5 / #6) alongside the input-layer
        # check (fix #2).
        if aucell_run_info.get("unmatched"):
            _u = aucell_run_info["unmatched"]
            aucell_qc.setdefault("warnings", []).append(
                f"{len(_u)} / {aucell_run_info['n_query_requested']} signature "
                f"genes did not match the atlas gene-name lookup — first 20: "
                f"{_u[:20]}. These genes contribute nothing to the score."
            )
        if aucell_run_info.get("n_top_bumped"):
            aucell_qc.setdefault("warnings", []).append(
                f"Signature ({aucell_run_info['n_query_matched']} genes) is larger "
                f"than the requested top_fraction = {aucell_run_info['requested_top_fraction']:.2%}, "
                f"so n_top was bumped to {aucell_run_info['n_top']} "
                f"({aucell_run_info['effective_top_fraction']:.2%} of genes). "
                f"The AUCell window is wider than the slider suggests — "
                f"scores are less stringent than intended. Reduce 'Top N genes "
                f"for scoring' or raise 'AUCell top-ranked fraction' to align "
                f"the two."
            )
        aucell_qc.setdefault("info", []).append(
            f"AUCell window: n_top = {aucell_run_info.get('n_top', '?')} genes "
            f"({aucell_run_info.get('effective_top_fraction', 0) * 100:.2f}% of "
            f"{adata.raw.n_vars if adata.raw is not None else adata.n_vars}); "
            f"signature matched {aucell_run_info.get('n_query_matched', '?')} / "
            f"{aucell_run_info.get('n_query_requested', '?')} genes."
        )

        progress.progress(83, text="Computing per-cluster enrichment significance...")

        # ---- AUCell result tables (raw data underlying figures 1a–1c + S2/S3) ----
        _cell_labels_arr = adata.obs[annotation_col].values.astype(str)
        aucell_per_cell_df = pd.DataFrame({
            "cell_id": adata.obs_names.astype(str),
            "cluster": _cell_labels_arr,
            "aucell_score": aucell_scores,
        })

        # Per-cluster significance (fix #3: Welch's one-sided t-test
        # cluster-vs-rest with BH-FDR so users can separate "truly enriched"
        # from "small cluster with a slightly above-average mean")
        aucell_cluster_stats_df = compute_cluster_enrichment_stats(
            aucell_scores, _cell_labels_arr, min_cells=10, alpha=0.05,
        )

        # Preserve the original ordering (sorted by mean descending) for the
        # figures, but merge in the significance columns so the downloadable
        # table is the authoritative reference.
        _grp = aucell_per_cell_df.groupby("cluster")["aucell_score"]
        aucell_per_cluster_df = pd.DataFrame({
            "n_cells": _grp.count(),
            "mean": _grp.mean(),
            "median": _grp.median(),
            "std": _grp.std(ddof=1),
        })
        aucell_per_cluster_df["sem"] = (
            aucell_per_cluster_df["std"] / np.sqrt(aucell_per_cluster_df["n_cells"])
        )
        aucell_per_cluster_df = (
            aucell_per_cluster_df.sort_values("mean", ascending=False)
            .reset_index()
        )
        if not aucell_cluster_stats_df.empty:
            aucell_per_cluster_df = aucell_per_cluster_df.merge(
                aucell_cluster_stats_df[
                    ["cluster", "t_stat", "pvalue", "qvalue", "significant"]
                ],
                on="cluster", how="left",
            )

        progress.progress(85, text="Computing composite ranking...")

        # ---- Composite ranking ----
        composite_df = compute_composite_ranking(
            corr_df, fisher_df, nnls_df,
            gsea_df if len(gsea_df) > 0 else None,
        )

        progress.progress(95, text="Generating figures...")

        # ---- Subsample for UMAP ----
        sub_indices = None
        if adata.n_obs > umap_subsample:
            rng = np.random.default_rng(42)
            sub_indices = np.sort(rng.choice(adata.n_obs, size=umap_subsample, replace=False))

        # Get UMAP coordinates (existence already verified at load time)
        umap_key = "X_umap" if "X_umap" in adata.obsm else "X_UMAP"
        umap_coords = adata.obsm[umap_key]
        cell_labels = adata.obs[annotation_col].values.astype(str)

        progress.progress(100, text="Analysis complete!")
        progress_placeholder.empty()

        # ---- Cache all analysis results in session state ----
        st.session_state._analysis_cache = {
            "params": _analysis_params,
            "bactrap_matched": bactrap_matched,
            "matched_genes": matched_genes,
            "gene_to_idx": gene_to_idx,
            "matched_in_raw": matched_in_raw,
            "gene_indices": gene_indices,
            "cluster_mean_expr": cluster_mean_expr,
            "enriched_df": enriched_df,
            "enriched_genes_list": enriched_genes_list,
            "top_enriched_genes": top_enriched_genes,
            "corr_df": corr_df,
            "markers": markers,
            "fisher_df": fisher_df,
            "enriched_gene_indices": enriched_gene_indices,
            "top_clusters_corr": top_clusters_corr,
            "frac_expr": frac_expr,
            "enriched_mean_expr": enriched_mean_expr,
            "zscore_df": zscore_df,
            "nnls_df": nnls_df,
            "gsea_df": gsea_df,
            "gsea_running_scores": gsea_running_scores,
            "gsea_ranked_genes": gsea_ranked_genes,
            "aucell_scores": aucell_scores,
            "aucell_qc": aucell_qc,
            "aucell_per_cell_df": aucell_per_cell_df,
            "aucell_per_cluster_df": aucell_per_cluster_df,
            "aucell_cluster_stats_df": aucell_cluster_stats_df,
            "composite_df": composite_df,
            "sub_indices": sub_indices,
            "umap_coords": umap_coords,
            "cell_labels": cell_labels,
            "enriched_sorted": enriched_sorted,
            "top_genes_heatmap": top_genes_heatmap,
            "top_clusters_heatmap": top_clusters_heatmap,
        }
    else:
        # ---- Restore cached results (no recomputation needed) ----
        _c = _cached
        bactrap_matched = _c["bactrap_matched"]
        matched_genes = _c["matched_genes"]
        gene_to_idx = _c["gene_to_idx"]
        matched_in_raw = _c["matched_in_raw"]
        gene_indices = _c["gene_indices"]
        cluster_mean_expr = _c["cluster_mean_expr"]
        enriched_df = _c["enriched_df"]
        enriched_genes_list = _c["enriched_genes_list"]
        top_enriched_genes = _c["top_enriched_genes"]
        corr_df = _c["corr_df"]
        markers = _c["markers"]
        fisher_df = _c["fisher_df"]
        enriched_gene_indices = _c["enriched_gene_indices"]
        top_clusters_corr = _c["top_clusters_corr"]
        frac_expr = _c["frac_expr"]
        enriched_mean_expr = _c["enriched_mean_expr"]
        zscore_df = _c["zscore_df"]
        nnls_df = _c["nnls_df"]
        gsea_df = _c["gsea_df"]
        gsea_running_scores = _c["gsea_running_scores"]
        gsea_ranked_genes = _c["gsea_ranked_genes"]
        aucell_scores = _c["aucell_scores"]
        aucell_qc = _c.get("aucell_qc", {"warnings": [], "info": []})
        aucell_per_cell_df = _c["aucell_per_cell_df"]
        aucell_per_cluster_df = _c["aucell_per_cluster_df"]
        aucell_cluster_stats_df = _c.get("aucell_cluster_stats_df", pd.DataFrame())
        composite_df = _c["composite_df"]
        sub_indices = _c["sub_indices"]
        umap_coords = _c["umap_coords"]
        cell_labels = _c["cell_labels"]
        enriched_sorted = _c["enriched_sorted"]
        top_genes_heatmap = _c["top_genes_heatmap"]
        top_clusters_heatmap = _c["top_clusters_heatmap"]
        progress_placeholder.empty()

    st.session_state.analysis_done = True

    # ---- Cre-driver sanity stats (up-front so the baseline filter below can
    # feed the composite ranking, and the sanity tab can reuse the cache) ----
    # Floor at min(markers gate, AUCell gate) so the baseline filter never
    # drops an AUCell-eligible cluster just because it fell below sanity's
    # own size gate (user raising min_cells_per_cluster above
    # min_cells_for_rank would otherwise silently tighten AUCell too).
    _sanity_min_cells = min(min_cells_per_cluster, min_cells_for_rank)
    sanity_cache_key = (
        hypomap_file.strip(), annotation_col,
        _sanity_min_cells, sanity_gene.strip().lower(),
    )
    _sanity_cached = st.session_state.get("_sanity_cache")
    if _sanity_cached is not None and _sanity_cached.get("key") == sanity_cache_key:
        sanity_stats = _sanity_cached["stats"]
    else:
        with st.spinner(f"Computing per-cluster {sanity_gene} expression..."):
            sanity_stats = compute_single_gene_cluster_stats(
                adata, sanity_gene.strip(), annotation_col,
                adata_gene_lookup=_adata_lookup,
                has_raw=_adata_has_raw,
                min_cells=_sanity_min_cells,
                normalize=True,
            )
        st.session_state["_sanity_cache"] = {
            "key": sanity_cache_key, "stats": sanity_stats,
        }

    # ---- Optional baseline Cre-driver expression filter ----
    # When the slider is > 0, drop clusters whose mean Cre-driver expression
    # falls below the floor BEFORE re-running the composite vote so
    # correlation / Fisher / NNLS / GSEA all exclude them and the survivors
    # are re-ranked against each other.  Applied post-cache so toggling the
    # slider doesn't invalidate the expensive parts of the analysis.
    baseline_allowed = None
    baseline_filter_state = "disabled"  # "disabled" | "active" | "broken_missing_gene" | "broken_empty"
    if sanity_baseline_mean_expr > 0:
        if sanity_stats is None or len(sanity_stats) == 0:
            baseline_filter_state = "broken_missing_gene"
        else:
            baseline_allowed = set(
                sanity_stats.index[
                    sanity_stats["mean_expr"] >= sanity_baseline_mean_expr
                ].astype(str)
            )
            if len(baseline_allowed) == 0:
                baseline_allowed = None
                baseline_filter_state = "broken_empty"
            else:
                baseline_filter_state = "active"

    # ---- Optional display-time filter for Unassigned / Mixed clusters ----
    # HypoMap's "Unassigned" and "Mixed" clusters are uncurated aggregates
    # that can artificially top rankings (especially GSEA) because of
    # heterogeneous membership.  Filtering is applied after the analysis
    # cache so toggling the checkbox doesn't invalidate computations; it
    # only changes what the tabs display.
    _display_filter_active = hide_unassigned or baseline_allowed is not None
    if _display_filter_active:
        _excl_pat = re.compile(r"Unassigned|Mixed", re.IGNORECASE)

        def _filter_by_cluster(df, col="cluster"):
            if df is None or len(df) == 0 or col not in df.columns:
                return df
            mask = pd.Series(True, index=df.index)
            if hide_unassigned:
                mask &= ~df[col].astype(str).str.contains(_excl_pat, na=False)
            if baseline_allowed is not None:
                mask &= df[col].astype(str).isin(baseline_allowed)
            return df[mask].reset_index(drop=True)

        corr_df = _filter_by_cluster(corr_df)
        fisher_df = _filter_by_cluster(fisher_df)
        nnls_df = _filter_by_cluster(nnls_df)
        gsea_df = _filter_by_cluster(gsea_df)

        # When the baseline filter is active, re-run the composite vote on
        # the survivors so their percentiles reflect the restricted universe
        # (the user-requested behavior: "correlation/Fisher/NNLS/GSEA all
        # exclude those clusters" before the consensus is formed).  For
        # hide_unassigned alone, filtering the cached composite by cluster
        # is equivalent and cheaper.
        if baseline_allowed is not None:
            composite_df = compute_composite_ranking(
                corr_df, fisher_df, nnls_df,
                gsea_df if gsea_df is not None and len(gsea_df) > 0 else None,
            )
        else:
            composite_df = _filter_by_cluster(composite_df)

        # Recompute selection lists that drive dotplot / heatmap cluster sets
        # so they stay consistent with the filtered rankings.
        if len(corr_df) > 0 and "both_significant" in corr_df.columns:
            _sig = corr_df.loc[corr_df["both_significant"], "cluster"].tolist()
            top_clusters_corr = _sig[:15] if len(_sig) >= 5 else corr_df["cluster"].tolist()[:15]
            top_clusters_heatmap = _sig[:20] if len(_sig) >= 5 else corr_df["cluster"].tolist()[:20]
        else:
            top_clusters_corr = corr_df["cluster"].tolist()[:15] if len(corr_df) > 0 else []
            top_clusters_heatmap = corr_df["cluster"].tolist()[:20] if len(corr_df) > 0 else []

        # zscore_df is (genes × clusters); when the filter is active the
        # cached matrix still holds the original unfiltered top-20 columns,
        # and the heatmap (Suppl. S5) + heatmap_zscores.csv would otherwise
        # show clusters that fail the filter.  Recompute against the
        # filtered top_clusters_heatmap and the full cluster_mean_expr.
        zscore_df = compute_zscore_heatmap_data(
            cluster_mean_expr, top_genes_heatmap, top_clusters_heatmap,
        )

    # Apply the baseline filter to every AUCell output — the user's intent is
    # to remove Cre-driver-negative clusters from the AUCell analysis, so the
    # cluster-level panels (barplot S2, violin 1c, cell-type UMAP highlight
    # 1b, per-cluster CSV) AND the per-cell panels (AUCell UMAP fig 1a,
    # per-cell CSV) all drop cells belonging to filtered-out clusters.
    sub_indices_filtered = sub_indices
    if baseline_allowed is not None:
        aucell_per_cluster_df = aucell_per_cluster_df[
            aucell_per_cluster_df["cluster"].astype(str).isin(baseline_allowed)
        ].reset_index(drop=True)
        _baseline_cell_mask = np.isin(cell_labels, list(baseline_allowed))
        aucell_per_cell_df = aucell_per_cell_df[_baseline_cell_mask].reset_index(drop=True)
        if sub_indices is None:
            sub_indices_filtered = np.where(_baseline_cell_mask)[0]
        else:
            sub_indices_filtered = sub_indices[_baseline_cell_mask[sub_indices]]

    # Sidebar-side mirror of the filter state so the user sees the
    # effective setting even if they've scrolled past the global banner.
    if baseline_filter_state == "active":
        st.sidebar.caption(
            f"Baseline {sanity_gene} filter: **{len(baseline_allowed)}"
            f"/{len(sanity_stats)}** clusters pass mean ≥ "
            f"{sanity_baseline_mean_expr:.2f}."
        )
    elif baseline_filter_state == "broken_missing_gene":
        st.sidebar.caption(
            f":warning: `{sanity_gene}` not in atlas — filter ignored."
        )
    elif baseline_filter_state == "broken_empty":
        st.sidebar.caption(
            f":warning: threshold {sanity_baseline_mean_expr:.2f} rejects "
            f"every cluster — filter ignored."
        )

    # Cache figure bytes so the Export tab doesn't regenerate them.
    # Only reset when a new analysis run is triggered (run_button pressed
    # or parameters changed), not on every Streamlit rerun.
    if _need_recompute:
        st.session_state.fig_bytes = {}
        st.session_state.table_bytes = {}

    if "fig_bytes" not in st.session_state:
        st.session_state.fig_bytes = {}
    if "table_bytes" not in st.session_state:
        st.session_state.table_bytes = {}

    # Pre-encode AUCell result tables.  Both the per-cell (~400k rows) and
    # per-cluster (~185 rows) CSVs shrink when the baseline filter drops
    # clusters, so re-serialise whenever the filter signature changes.
    _baseline_sig = (
        tuple(sorted(baseline_allowed)) if baseline_allowed is not None else None
    )
    if st.session_state.table_bytes.get("_aucell_per_cell_sig") != _baseline_sig:
        st.session_state.table_bytes["aucell_per_cell"] = (
            aucell_per_cell_df.to_csv(index=False).encode()
        )
        st.session_state.table_bytes["_aucell_per_cell_sig"] = _baseline_sig
    if st.session_state.table_bytes.get("_aucell_per_cluster_sig") != _baseline_sig:
        st.session_state.table_bytes["aucell_per_cluster"] = (
            aucell_per_cluster_df.to_csv(index=False).encode()
        )
        st.session_state.table_bytes["_aucell_per_cluster_sig"] = _baseline_sig

    def _cache_fig(name, fig):
        st.session_state.fig_bytes[name] = {
            "pdf": fig_to_bytes(fig, "pdf"),
            "svg": fig_to_bytes(fig, "svg"),
        }

    # Filename suffix applied to filter-affected CSV downloads so a
    # collaborator who opens a 40-row composite_ranking.csv can tell from
    # the filename alone that it's a Pnoc-ge-0.05 subset, not the full 185.
    # Dots are replaced with 'p' (safe on every filesystem).
    def _filtered_name(basename: str) -> str:
        if baseline_allowed is None:
            return basename
        stem, _, ext = basename.rpartition(".")
        tag = f"{sanity_gene.lower()}_ge{sanity_baseline_mean_expr:.2f}".replace(".", "p")
        return f"{stem}_{tag}.{ext}"

    # ---- Global baseline-filter status banner (above the tab group) ----
    # Rendered once so every tab — not just AUCell — makes the filter state
    # obvious.  broken_* states distinguish "slider > 0 but doing nothing"
    # from "slider at 0" so a user whose threshold discards everything
    # doesn't silently see an unfiltered dashboard.
    if baseline_filter_state == "active":
        st.info(
            f"**Baseline {sanity_gene} filter active** — "
            f"{len(baseline_allowed)} / {len(sanity_stats)} clusters pass "
            f"mean {sanity_gene} ≥ {sanity_baseline_mean_expr:.2f}. "
            f"Applies to every cluster-level ranking (correlation, "
            f"Fisher, NNLS, GSEA, composite, AUCell cluster figs 1b / S2 "
            f"/ 1c, heatmap S5) **and** to the per-cell AUCell panel "
            f"(fig 1a) and per-cell CSV, which drop cells belonging to "
            f"filtered-out clusters."
        )
    elif baseline_filter_state == "broken_missing_gene":
        st.warning(
            f"**Baseline filter ignored** — `{sanity_gene}` not found in "
            f"the HypoMap atlas (slider at {sanity_baseline_mean_expr:.2f}). "
            f"Tabs show unfiltered data. Check gene-symbol casing "
            f"(e.g. `Pnoc`, not `PNOC`)."
        )
    elif baseline_filter_state == "broken_empty":
        st.error(
            f"**Baseline filter ignored** — threshold "
            f"mean {sanity_gene} ≥ {sanity_baseline_mean_expr:.2f} "
            f"discards every cluster. Tabs show unfiltered data. "
            f"Lower the slider."
        )

    # ======================================================================
    # TAB 1: Data Overview
    # ======================================================================
    with tab1:
        st.header("Data Overview")

        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("bacTRAP genes", len(bactrap_df))
        with col2:
            st.metric("HypoMap cells", f"{adata.n_obs:,}")
        with col3:
            n_clusters = adata.obs[annotation_col].nunique()
            st.metric(f"Clusters ({annotation_col})", n_clusters)

        st.markdown("---")

        col4, col5, col6 = st.columns(3)
        with col4:
            st.metric("Genes matched", len(matched_genes))
        with col5:
            st.metric("Enriched genes (padj + FC)", len(enriched_df))
        with col6:
            match_pct = (len(matched_genes) / len(bactrap_df) * 100) if len(bactrap_df) > 0 else 0.0
            st.metric("Match rate", f"{match_pct:.1f}%")

        st.subheader("Gene Matching Diagnostics")

        # Show which column was used and sample gene names from each dataset
        diag_col1, diag_col2 = st.columns(2)
        with diag_col1:
            st.markdown(f"**bacTRAP gene column:** `{gene_col_selection}`")
            if gene_col_selection == "(use row index)":
                _sample_bt = [str(x) for x in bactrap_df.index[:10]]
            else:
                _sample_bt = bactrap_df[gene_col_selection].dropna().head(10).astype(str).tolist()
            st.markdown("Sample bacTRAP gene IDs:")
            st.code("\n".join(_sample_bt))
        with diag_col2:
            _sample_hm = [str(x) for x in _adata_gnames[:10]]
            st.markdown(f"**HypoMap gene names** (resolved, n={len(_adata_gnames)}):")
            st.code("\n".join(_sample_hm))
            # Show raw var_names too if different
            if _adata_has_raw and adata.raw is not None:
                _raw_vn = [str(x) for x in adata.raw.var_names[:5]]
                if _raw_vn != [str(x) for x in _adata_gnames[:5]]:
                    st.markdown("Raw var_names (first 5):")
                    st.code("\n".join(_raw_vn))

        # Bar chart of matched vs unmatched
        match_data = pd.DataFrame({
            "Category": ["Matched", "Unmatched"],
            "Count": [len(matched_genes), len(bactrap_df) - len(matched_genes)],
        })
        st.bar_chart(match_data.set_index("Category"))

        st.subheader("bacTRAP Data Preview")
        st.dataframe(bactrap_df.head(20), use_container_width=True)

        st.subheader("Top Enriched Genes")
        st.caption(
            f"Filters: padj < {padj_cutoff:.3g}, log₂FC > {log2fc_cutoff:.2f}, "
            f"IP ≥ {min_ip_expression:.0f}. "
            f"Ranked by **{ranking_metric_label}**."
        )
        if len(enriched_df) > 0:
            display_cols = ["_hypomap_gene_name", "log2FoldChange", "padj", "IP", "Input"]
            available_cols = [c for c in display_cols if c in enriched_df.columns]
            st.dataframe(
                enriched_sorted[available_cols].head(30).reset_index(drop=True),
                use_container_width=True,
            )
        else:
            st.warning("No genes pass the current enrichment thresholds.")

        st.subheader("Figure: bacTRAP Volcano Plot")
        fig_volcano = figure_bactrap_volcano(
            bactrap_matched,
            highlight_genes=["Pnoc"],
            padj_cutoff=padj_cutoff,
            log2fc_cutoff=log2fc_cutoff,
            double_column=double_column,
        )
        st.pyplot(fig_volcano)
        _cache_fig("fig_volcano_bactrap", fig_volcano)

        col_pdf, col_svg = st.columns(2)
        with col_pdf:
            st.download_button(
                "Download PDF",
                st.session_state.fig_bytes["fig_volcano_bactrap"]["pdf"],
                "fig_volcano_bactrap.pdf", "application/pdf",
                key="dl_fig_volcano_bt_pdf",
            )
        with col_svg:
            st.download_button(
                "Download SVG",
                st.session_state.fig_bytes["fig_volcano_bactrap"]["svg"],
                "fig_volcano_bactrap.svg", "image/svg+xml",
                key="dl_fig_volcano_bt_svg",
            )
        plt.close(fig_volcano)

    # ======================================================================
    # TAB: AUCell (MAIN FIGURE)
    # ======================================================================
    with tab_aucell:
        st.header("Main Figure: AUCell Enrichment Analysis")
        st.markdown(
            "AUCell (rank-based Area Under the Curve) is the primary cell-type "
            "mapping method — it is **normalization-insensitive**, **threshold-free**, "
            "and quantifies per-cell enrichment of the bacTRAP gene set. "
            f"Scores computed from the top **{len(top_enriched_genes)}** enriched genes."
        )

        # Input-layer QC (fix #2) — warn loudly when the layer fed into
        # AUCell does not look like raw counts.
        _qc_warnings = aucell_qc.get("warnings", []) if aucell_qc else []
        if _qc_warnings:
            for _w in _qc_warnings:
                st.warning(_w)
        else:
            _info = aucell_qc.get("info", []) if aucell_qc else []
            if _info:
                with st.expander("AUCell input QC", expanded=False):
                    for _line in _info:
                        st.caption(_line)

        # Per-cluster significance summary (fix #3)
        if aucell_cluster_stats_df is not None and not aucell_cluster_stats_df.empty:
            _n_tested = len(aucell_cluster_stats_df)
            _n_sig = int(aucell_cluster_stats_df["significant"].sum())
            st.markdown(
                f"**Enrichment significance (Welch's t, BH-FDR):** "
                f"{_n_sig} / {_n_tested} clusters pass **q < 0.05** "
                f"(cluster AUCell distribution vs. rest of atlas). "
                f"Full per-cluster statistics — `t_stat`, `pvalue`, `qvalue`, "
                f"`significant` — are appended to the downloadable "
                f"`aucell_per_cluster.csv`."
            )

        # Figure 1a: AUCell UMAP
        st.subheader("Figure 1a: AUCell Enrichment UMAP")
        _fig_1a_n_cells = int(len(aucell_per_cell_df))
        _fig_1a_filter_clause = (
            f" restricted to clusters passing the baseline {sanity_gene} "
            f"≥ {sanity_baseline_mean_expr:.2f} filter"
            if baseline_allowed is not None else ""
        )
        st.markdown(
            f"**Figure 1a.** AUCell enrichment score for the bacTRAP signature "
            f"projected onto the HypoMap UMAP embedding"
            f"{_fig_1a_filter_clause} "
            f"(n = {_fig_1a_n_cells:,} cells). For each cell, the area under the "
            f"recovery curve for the top "
            f"**{len(top_enriched_genes)}** π-score-ranked bacTRAP-enriched "
            f"genes was computed within the top "
            f"**{aucell_top_fraction:.0%}** of genes by expression rank in "
            f"that cell and normalised to its theoretical maximum, so values "
            f"lie on [0, 1]. Colour encodes AUCell score (magma colormap); "
            f"the scale is clipped at the 2nd and 98th percentiles to "
            f"suppress outlier saturation. Bottom-left arrows mark the "
            f"UMAP1 / UMAP2 axes. See Methods (*AUCell scoring*) for the "
            f"full derivation."
        )
        fig_1a = figure_aucell_umap(
            umap_coords, aucell_scores,
            double_column=double_column,
            subsample_idx=sub_indices_filtered,
        )
        st.pyplot(fig_1a)
        _cache_fig("fig_1a_aucell_umap", fig_1a)

        col_pdf, col_svg, col_csv = st.columns(3)
        with col_pdf:
            st.download_button(
                "Download PDF",
                st.session_state.fig_bytes["fig_1a_aucell_umap"]["pdf"],
                "fig_1a_aucell_umap.pdf", "application/pdf",
                key="dl_fig_1a_pdf",
            )
        with col_svg:
            st.download_button(
                "Download SVG",
                st.session_state.fig_bytes["fig_1a_aucell_umap"]["svg"],
                "fig_1a_aucell_umap.svg", "image/svg+xml",
                key="dl_fig_1a_svg",
            )
        with col_csv:
            st.download_button(
                "Download CSV (per-cell)",
                st.session_state.table_bytes["aucell_per_cell"],
                _filtered_name("aucell_per_cell.csv"), "text/csv",
                key="dl_fig_1a_csv",
                help=(
                    "cell_id, cluster, aucell_score — filtered to "
                    "clusters passing the baseline filter when active."
                ),
            )
        plt.close(fig_1a)

        # Figure 1b: Cell-type annotation UMAP (top-15 AUCell clusters)
        st.subheader("Figure 1b: Cell-type Annotation UMAP (top-15 AUCell clusters)")
        _filter_clause_1b = (
            f" and passing the baseline {sanity_gene} filter "
            f"(mean ≥ {sanity_baseline_mean_expr:.2f})"
            if baseline_allowed is not None else ""
        )
        st.markdown(
            f"**Figure 1b.** Same UMAP layout as (a), coloured by HypoMap "
            f"cell-type annotation (`{annotation_col}`). Only the 15 "
            f"clusters with the highest mean AUCell score (among clusters "
            f"with ≥ **{min_cells_for_rank}** cells{_filter_clause_1b}) "
            f"are drawn in colour and overplotted on a light-grey "
            f"background of all remaining cells; this keeps small but "
            f"highly enriched populations visible while conveying overall "
            f"atlas topology. The legend lists the highlighted clusters "
            f"in rank order (top entry = highest mean AUCell). The 20-cell "
            f"floor excludes small populations whose cluster mean is "
            f"dominated by shrinkage variance and would otherwise claim "
            f"top slots by chance — they remain in `aucell_per_cluster.csv`."
        )
        _size_filtered_ranked = (
            aucell_per_cluster_df[
                aucell_per_cluster_df["n_cells"] >= min_cells_for_rank
            ]
            .sort_values("mean", ascending=False)
        )
        _top15_aucell_clusters_ct = (
            _size_filtered_ranked.head(15)["cluster"].astype(str).tolist()
        )
        fig_1a_ct = figure_celltype_umap(
            umap_coords=umap_coords,
            cell_labels=cell_labels,
            highlight_clusters=_top15_aucell_clusters_ct,
            double_column=double_column,
            subsample_idx=sub_indices,
        )
        st.pyplot(fig_1a_ct)
        _cache_fig("fig_1a_celltype_umap", fig_1a_ct)

        col_pdf, col_svg, col_csv = st.columns(3)
        with col_pdf:
            st.download_button(
                "Download PDF",
                st.session_state.fig_bytes["fig_1a_celltype_umap"]["pdf"],
                "fig_1a_celltype_umap.pdf", "application/pdf",
                key="dl_fig_1a_ct_pdf",
            )
        with col_svg:
            st.download_button(
                "Download SVG",
                st.session_state.fig_bytes["fig_1a_celltype_umap"]["svg"],
                "fig_1a_celltype_umap.svg", "image/svg+xml",
                key="dl_fig_1a_ct_svg",
            )
        with col_csv:
            st.download_button(
                "Download CSV (per-cluster mean)",
                st.session_state.table_bytes["aucell_per_cluster"],
                _filtered_name("aucell_per_cluster.csv"), "text/csv",
                key="dl_fig_1a_ct_csv",
                help="Top-15 rows are the highlighted clusters in this panel.",
            )
        plt.close(fig_1a_ct)

        st.markdown("---")
        st.caption(
            "Panels below are supplementary to the main figure (1a–c); they "
            "are published as supplementary figures in the Methods document."
        )

        # Supplementary S2: AUCell Cluster Barplot (was Figure 1b)
        st.subheader("Supplementary Figure S2: Mean AUCell Score per Cluster")
        _filter_clause_s2 = (
            f" and also filtered to clusters with mean {sanity_gene} ≥ "
            f"{sanity_baseline_mean_expr:.2f}"
            if baseline_allowed is not None else ""
        )
        st.markdown(
            f"**Supplementary Figure S2.** Horizontal barplot of mean AUCell "
            f"score per HypoMap cluster (top 25 by mean, clusters with "
            f"< **{min_cells_for_rank}** cells excluded from ranking"
            f"{_filter_clause_s2}). Error bars are standard error of the "
            f"mean. Bar colour encodes the cluster mean (magma colormap). "
            f"Source table: `aucell_per_cluster.csv`."
        )
        fig_1b = figure_aucell_cluster_barplot(
            aucell_scores, cell_labels,
            top_n=25, double_column=double_column,
            min_cluster_cells=min_cells_for_rank,
            allowed_clusters=baseline_allowed,
        )
        st.pyplot(fig_1b)
        _cache_fig("fig_1b_aucell_barplot", fig_1b)

        col_pdf, col_svg, col_csv = st.columns(3)
        with col_pdf:
            st.download_button(
                "Download PDF",
                st.session_state.fig_bytes["fig_1b_aucell_barplot"]["pdf"],
                "fig_1b_aucell_barplot.pdf", "application/pdf",
                key="dl_fig_1b_pdf",
            )
        with col_svg:
            st.download_button(
                "Download SVG",
                st.session_state.fig_bytes["fig_1b_aucell_barplot"]["svg"],
                "fig_1b_aucell_barplot.svg", "image/svg+xml",
                key="dl_fig_1b_svg",
            )
        with col_csv:
            st.download_button(
                "Download CSV (per-cluster)",
                st.session_state.table_bytes["aucell_per_cluster"],
                _filtered_name("aucell_per_cluster.csv"), "text/csv",
                key="dl_fig_1b_csv",
                help=(
                    "cluster, n_cells, mean, median, std, sem, t_stat, pvalue, "
                    "qvalue, significant — sorted by mean descending. Significance "
                    "is Welch's one-sided t (cluster > rest) with BH-FDR."
                ),
            )
        plt.close(fig_1b)

        # Figure 1c: AUCell Violin Plots
        st.subheader("Figure 1c: AUCell Score Distributions (Top-15 Clusters)")
        _filter_clause_1c = (
            f" Top-15 is taken over clusters with mean {sanity_gene} ≥ "
            f"{sanity_baseline_mean_expr:.2f}."
            if baseline_allowed is not None else ""
        )
        st.markdown(
            f"**Figure 1c.** Violin plots of the full AUCell score "
            f"distribution within each of the 15 top-ranked clusters from "
            f"(b), ordered from highest (top) to lowest (bottom) cluster "
            f"mean. The short **solid black** vertical bar inside each "
            f"violin marks the cluster mean; the **dashed grey** bar marks "
            f"the cluster median (the two nearly coincide when the "
            f"distribution is symmetric, in which case they read as a "
            f"single I-shape — see the on-figure legend). Violin fill "
            f"colour encodes the cluster mean (magma colormap). "
            f"Per-cluster means, medians, SEMs and Welch's one-sided "
            f"*t*-test *p*/*q*-values against the rest of the atlas are "
            f"available in `aucell_per_cluster.csv`.{_filter_clause_1c}"
        )
        fig_1c = figure_aucell_violins(
            aucell_scores, cell_labels,
            top_n=15, double_column=True,
            min_cluster_cells=min_cells_for_rank,
            allowed_clusters=baseline_allowed,
        )
        st.pyplot(fig_1c)
        _cache_fig("fig_1c_aucell_violins", fig_1c)

        col_pdf, col_svg, col_mean, col_cell = st.columns(4)
        with col_pdf:
            st.download_button(
                "Download PDF",
                st.session_state.fig_bytes["fig_1c_aucell_violins"]["pdf"],
                "fig_1c_aucell_violins.pdf", "application/pdf",
                key="dl_fig_1c_pdf",
            )
        with col_svg:
            st.download_button(
                "Download SVG",
                st.session_state.fig_bytes["fig_1c_aucell_violins"]["svg"],
                "fig_1c_aucell_violins.svg", "image/svg+xml",
                key="dl_fig_1c_svg",
            )
        with col_mean:
            st.download_button(
                "Download CSV (per-cluster mean)",
                st.session_state.table_bytes["aucell_per_cluster"],
                _filtered_name("aucell_per_cluster.csv"), "text/csv",
                key="dl_fig_1c_mean_csv",
                help=(
                    "Plotted quantities: cluster, n_cells, mean (black bar), "
                    "median (grey dashed bar), std, sem, plus Welch's t-test "
                    "significance columns (t_stat, pvalue, qvalue, significant). "
                    "Sorted by mean descending — the top 15 rows are the "
                    "clusters shown in the violin."
                ),
            )
        with col_cell:
            st.download_button(
                "Download CSV (per-cell)",
                st.session_state.table_bytes["aucell_per_cell"],
                "aucell_per_cell.csv", "text/csv",
                key="dl_fig_1c_csv",
                help="Per-cell AUCell scores — use to reconstruct the full violin shape.",
            )
        plt.close(fig_1c)

        # Supplementary S3: AUCell Score Histogram (was Figure 1d)
        st.subheader("Supplementary Figure S3: Global AUCell Score Distribution")
        st.markdown(
            "**Supplementary Figure S3.** Global histogram of per-cell "
            "AUCell scores across the atlas. Dashed vertical lines mark "
            "the 90th, 95th and 99th percentiles as well as the mean. "
            "Cells above the 95th percentile form the candidate pool "
            "for belonging to the bacTRAP target population."
        )
        fig_1d = figure_aucell_histogram(
            aucell_scores, double_column=double_column,
        )
        st.pyplot(fig_1d)
        _cache_fig("fig_1d_aucell_histogram", fig_1d)

        col_pdf, col_svg, col_csv = st.columns(3)
        with col_pdf:
            st.download_button(
                "Download PDF",
                st.session_state.fig_bytes["fig_1d_aucell_histogram"]["pdf"],
                "fig_1d_aucell_histogram.pdf", "application/pdf",
                key="dl_fig_1d_pdf",
            )
        with col_svg:
            st.download_button(
                "Download SVG",
                st.session_state.fig_bytes["fig_1d_aucell_histogram"]["svg"],
                "fig_1d_aucell_histogram.svg", "image/svg+xml",
                key="dl_fig_1d_svg",
            )
        with col_csv:
            st.download_button(
                "Download CSV (per-cell)",
                st.session_state.table_bytes["aucell_per_cell"],
                "aucell_per_cell.csv", "text/csv",
                key="dl_fig_1d_csv",
                help="Per-cell AUCell scores — source data for the histogram.",
            )
        plt.close(fig_1d)

        # Supplementary S11: Composite Consensus Ranking (was Figure 1e)
        st.markdown("---")
        st.subheader("Supplementary Figure S11: Composite Consensus Ranking")
        st.markdown(
            "**Supplementary Figure S11.** Composite consensus ranking "
            "across Spearman correlation, Fisher's exact marker overlap, "
            "NNLS deconvolution, and preranked GSEA (top 20 clusters). "
            "Each method's cluster scores are converted to percentile "
            "ranks (0–1) and averaged (`nanmean`, so methods with missing "
            "data for a given cluster are excluded rather than penalised "
            "with zero). Validation panel — confirms the AUCell-ranked "
            "clusters are also prioritised by orthogonal methods."
        )
        if len(composite_df) > 0:
            fig_1e = figure_composite_ranking(composite_df, top_n=20, double_column=True)
            st.pyplot(fig_1e)
            _cache_fig("fig_1e_composite_ranking", fig_1e)

            col_pdf, col_svg = st.columns(2)
            with col_pdf:
                st.download_button(
                    "Download PDF",
                    st.session_state.fig_bytes["fig_1e_composite_ranking"]["pdf"],
                    "fig_1e_composite_ranking.pdf", "application/pdf",
                    key="dl_fig_1e_pdf",
                )
            with col_svg:
                st.download_button(
                    "Download SVG",
                    st.session_state.fig_bytes["fig_1e_composite_ranking"]["svg"],
                    "fig_1e_composite_ranking.svg", "image/svg+xml",
                    key="dl_fig_1e_svg",
                )
            plt.close(fig_1e)
        else:
            st.warning("No composite ranking data available.")

    # ======================================================================
    # TAB: Cre-driver Sanity Check
    # ======================================================================
    with tab_sanity:
        st.header(f"Cre-driver Sanity Check: {sanity_gene}")
        st.markdown(
            f"In a **{sanity_gene}-Cre;NuTRAP** experiment we'd expect the "
            f"top-ranked HypoMap clusters to express *{sanity_gene}*. This panel "
            f"shows mean expression and fraction of cells expressing "
            f"*{sanity_gene}* across the top-ranked clusters from the "
            f"composite consensus."
        )
        st.info(
            "**Caveats — read before interpreting:** "
            "(1) Cre-lox is permanent lineage tracing — any cell that ever "
            f"expressed *{sanity_gene}* during development is labelled, even "
            f"if current mRNA is undetectable. (2) HypoMap is single-nucleus "
            f"data; neuropeptides like *{sanity_gene}* are notoriously prone "
            "to dropout. (3) Atlas expression is a snapshot; *Pnoc* is "
            "state-dependent (feeding, stress, estrous). Use this as a "
            "**confidence weight**, not a hard filter."
        )

        # sanity_stats is computed up-front (see post-cache block earlier)
        # so the baseline-expression slider can feed the composite ranking.
        # The tab just consumes the already-cached result.
        if sanity_stats is None or len(sanity_stats) == 0:
            st.error(
                f"**`{sanity_gene}`** was not found in the HypoMap atlas. "
                "Check the spelling and capitalization (mouse symbols are "
                "title-cased: `Pnoc`, not `PNOC` or `pnoc`)."
            )
        else:
            # ---- Cluster ordering: composite ranking, fall back to Pnoc mean ----
            _filter_suffix = (
                f" (after {sanity_gene} ≥ {sanity_baseline_mean_expr:.2f} filter)"
                if baseline_allowed is not None else ""
            )
            if len(composite_df) > 0:
                ranked_clusters = composite_df["cluster"].astype(str).tolist()
                rank_source = f"composite consensus{_filter_suffix}"
            elif len(corr_df) > 0:
                ranked_clusters = corr_df["cluster"].astype(str).tolist()
                rank_source = f"Spearman correlation{_filter_suffix}"
            else:
                ranked_clusters = sanity_stats.sort_values(
                    "mean_expr", ascending=False,
                ).index.tolist()
                rank_source = f"{sanity_gene} expression (no mapping ranking available)"

            top_n_sanity = st.slider(
                "Top N clusters to display", 5, 50, 20, 1,
                key="sanity_top_n",
                help=(
                    "How many of the composite-consensus top clusters to "
                    "show on the Cre-driver sanity panel. Widens or "
                    "narrows the barplot; does not change any computation."
                ),
            )
            top_clusters_sanity = ranked_clusters[:top_n_sanity]

            # Build the merged display table
            sanity_table = sanity_stats.loc[
                [c for c in top_clusters_sanity if c in sanity_stats.index]
            ].copy()
            # Add rank column from composite (or whichever source we used)
            sanity_table["rank"] = [
                ranked_clusters.index(c) + 1 if c in ranked_clusters else np.nan
                for c in sanity_table.index
            ]
            sanity_table = sanity_table.sort_values("rank")
            sanity_table["passes_threshold"] = (
                sanity_table["fraction_expressing"] >= sanity_fraction_threshold
            )

            # ---- Summary metric ----
            n_pass = int(sanity_table["passes_threshold"].sum())
            n_total = len(sanity_table)
            atlas_median_frac = float(sanity_stats["fraction_expressing"].median())

            mcol1, mcol2, mcol3 = st.columns(3)
            with mcol1:
                st.metric(
                    f"Top-{n_total} clusters expressing {sanity_gene}",
                    f"{n_pass} / {n_total}",
                    help=(
                        f"Clusters where ≥{sanity_fraction_threshold*100:.0f}% "
                        "of cells have non-zero counts for the Cre-driver gene."
                    ),
                )
            with mcol2:
                st.metric(
                    f"Atlas-wide median fraction expressing",
                    f"{atlas_median_frac*100:.1f}%",
                    help=(
                        "Median across ALL clusters in the atlas — useful "
                        "baseline for judging dropout."
                    ),
                )
            with mcol3:
                top_frac = float(sanity_table["fraction_expressing"].max()) if n_total else 0.0
                st.metric(
                    f"Highest fraction in top-{n_total}",
                    f"{top_frac*100:.1f}%",
                )

            st.caption(f"Cluster ordering: **{rank_source}**.")

            # ---- Diagnostic figure ----
            st.subheader(f"Figure: {sanity_gene} expression across top-ranked clusters")
            fig_sanity = figure_marker_gene_diagnostic(
                sanity_stats, top_clusters_sanity,
                gene_name=sanity_gene,
                fraction_threshold=sanity_fraction_threshold,
                double_column=True,
            )
            st.pyplot(fig_sanity)
            _cache_fig("fig_sanity_check", fig_sanity)

            col_pdf, col_svg = st.columns(2)
            with col_pdf:
                st.download_button(
                    "Download PDF",
                    st.session_state.fig_bytes["fig_sanity_check"]["pdf"],
                    f"fig_sanity_{sanity_gene.lower()}.pdf", "application/pdf",
                    key="dl_fig_sanity_pdf",
                )
            with col_svg:
                st.download_button(
                    "Download SVG",
                    st.session_state.fig_bytes["fig_sanity_check"]["svg"],
                    f"fig_sanity_{sanity_gene.lower()}.svg", "image/svg+xml",
                    key="dl_fig_sanity_svg",
                )
            plt.close(fig_sanity)

            # ---- Detailed table ----
            st.subheader(f"Per-cluster {sanity_gene} expression (top {n_total})")

            display_table = sanity_table[[
                "rank", "mean_expr", "fraction_expressing", "passes_threshold",
            ]].copy()
            display_table.columns = [
                "Composite rank", f"Mean {sanity_gene} (log-norm)",
                f"Fraction expressing {sanity_gene}",
                f"≥ {sanity_fraction_threshold*100:.0f}% threshold",
            ]

            def _highlight_fail(row):
                col = f"≥ {sanity_fraction_threshold*100:.0f}% threshold"
                if col in row.index and not row[col]:
                    return ["background-color: #fff3cd"] * len(row)
                return [""] * len(row)

            st.dataframe(
                display_table.style.apply(_highlight_fail, axis=1).format({
                    "Composite rank": "{:.0f}",
                    f"Mean {sanity_gene} (log-norm)": "{:.3f}",
                    f"Fraction expressing {sanity_gene}": "{:.1%}",
                }),
                use_container_width=True,
            )

            st.caption(
                "Highlighted rows = top-ranked clusters that **do not** pass "
                "the expression threshold. Investigate these: developmental "
                "lineage tracing, dropout, or potential artefact."
            )

            # ---- Download full table ----
            full_table = sanity_stats.copy()
            full_table["rank_in_mapping"] = [
                ranked_clusters.index(c) + 1 if c in ranked_clusters else np.nan
                for c in full_table.index
            ]
            full_table = full_table.reset_index()
            st.download_button(
                f"Download full {sanity_gene} per-cluster table (CSV)",
                full_table.to_csv(index=False).encode(),
                f"sanity_check_{sanity_gene.lower()}.csv", "text/csv",
                key="dl_sanity_csv",
            )

    # ======================================================================
    # TAB 3: Correlation Analysis (Supplementary)
    # ======================================================================
    with tab3:
        st.header("Supplementary: Enrichment Correlation Analysis")
        st.markdown(
            "Pearson and Spearman correlation between the bacTRAP log₂FC enrichment "
            "profile and mean expression per HypoMap cluster."
        )

        if len(corr_df) > 0:
            st.subheader("Top Correlated Clusters")

            # Highlight rows where Pearson is not significant (less credible)
            if "both_significant" in corr_df.columns:
                n_sig = int(corr_df["both_significant"].sum())
                n_total = len(corr_df)
                st.caption(
                    f"**{n_sig}/{n_total}** clusters significant by both Pearson and "
                    f"Spearman (BH-FDR < 0.05 across clusters). Rows failing the dual "
                    f"test are highlighted — concordance between parametric and rank-"
                    f"based correlations under multiple-testing correction is stronger "
                    f"evidence than either nominal p-value alone."
                )

                def _highlight_nonsig(row):
                    if "both_significant" in row.index and not row["both_significant"]:
                        return ["background-color: #fff3cd"] * len(row)
                    return [""] * len(row)

                styled = corr_df.style.apply(_highlight_nonsig, axis=1).format({
                    "pearson_r": "{:.4f}",
                    "pearson_pval": "{:.2e}",
                    "pearson_padj": "{:.2e}",
                    "spearman_r": "{:.4f}",
                    "spearman_pval": "{:.2e}",
                    "spearman_padj": "{:.2e}",
                })
            else:
                styled = corr_df.style.format({
                    "pearson_r": "{:.4f}",
                    "pearson_pval": "{:.2e}",
                    "spearman_r": "{:.4f}",
                    "spearman_pval": "{:.2e}",
                })

            st.dataframe(styled, use_container_width=True)

            st.subheader("Supplementary Figure S1: Correlation Barplot")
            fig_a = figure_correlation_barplot(corr_df, top_n=20, double_column=double_column)
            st.pyplot(fig_a)
            _cache_fig("fig_s1_correlation_barplot", fig_a)

            col_pdf, col_svg = st.columns(2)
            with col_pdf:
                st.download_button(
                    "Download PDF",
                    st.session_state.fig_bytes["fig_s1_correlation_barplot"]["pdf"],
                    "fig_s1_correlation.pdf", "application/pdf",
                    key="dl_fig_s1_pdf",
                )
            with col_svg:
                st.download_button(
                    "Download SVG",
                    st.session_state.fig_bytes["fig_s1_correlation_barplot"]["svg"],
                    "fig_s1_correlation.svg", "image/svg+xml",
                    key="dl_fig_s1_svg",
                )
            plt.close(fig_a)

            st.download_button(
                "Download correlation table (CSV)",
                corr_df.to_csv(index=False).encode(),
                _filtered_name("correlation_results.csv"), "text/csv",
                key="dl_corr_csv_tab2",
            )
        else:
            st.warning("No correlation results to display.")

    # ======================================================================
    # TAB 4: UMAP Projection (Supplementary)
    # ======================================================================
    with tab4:
        st.header("Supplementary: UMAP AUCell Projection")
        st.markdown(
            f"AUCell enrichment score projected onto the HypoMap UMAP "
            f"alongside the cell-type annotation for side-by-side comparison. "
            f"The left panel highlights the top-15 AUCell-ranked clusters "
            f"(matching main figures 1b/1c); all other cells are greyed out. "
            f"Score computed from the top **{len(top_enriched_genes)}** enriched "
            f"genes (padj < {padj_cutoff}, log₂FC > {log2fc_cutoff})."
        )

        st.subheader("Supplementary Figure S2: UMAP AUCell Map")
        _top15_aucell_clusters = (
            aucell_per_cluster_df[
                aucell_per_cluster_df["n_cells"] >= min_cells_for_rank
            ]
            .sort_values("mean", ascending=False)
            .head(15)["cluster"].astype(str).tolist()
        )
        fig_b = figure_umap_enrichment(
            umap_coords=umap_coords,
            cell_labels=cell_labels,
            enrichment_scores=aucell_scores,
            double_column=True,
            point_size=0.3,
            subsample_idx=sub_indices,
            score_title="AUCell enrichment score",
            score_label="AUCell score",
            highlight_clusters=_top15_aucell_clusters,
        )
        st.pyplot(fig_b)
        _cache_fig("fig_s2_umap_enrichment", fig_b)

        col_pdf, col_svg = st.columns(2)
        with col_pdf:
            st.download_button(
                "Download PDF",
                st.session_state.fig_bytes["fig_s2_umap_enrichment"]["pdf"],
                "fig_s2_umap.pdf", "application/pdf",
                key="dl_fig_s2_pdf",
            )
        with col_svg:
            st.download_button(
                "Download SVG",
                st.session_state.fig_bytes["fig_s2_umap_enrichment"]["svg"],
                "fig_s2_umap.svg", "image/svg+xml",
                key="dl_fig_s2_svg",
            )
        plt.close(fig_b)

    # ======================================================================
    # TAB 5: Marker Overlap (Supplementary)
    # ======================================================================
    with tab5:
        st.header("Supplementary: Marker Gene Overlap Analysis")
        st.markdown(
            "One-sided Fisher's exact test for overlap between bacTRAP-enriched genes "
            "and HypoMap cluster markers."
        )

        if len(fisher_df) > 0:
            st.subheader("Fisher's Exact Test Results")
            display_fisher = fisher_df[[
                "cluster", "overlap_count", "n_enriched", "n_markers",
                "odds_ratio", "pvalue", "padj", "overlap_genes",
            ]].copy()
            st.dataframe(
                display_fisher.style.format({
                    "odds_ratio": "{:.2f}",
                    "pvalue": "{:.2e}",
                    "padj": "{:.2e}",
                }),
                use_container_width=True,
            )

            st.subheader("Supplementary Figure S3: Enrichment Volcano Plot")
            fig_d = figure_volcano_enrichment(
                fisher_df, pval_threshold=padj_cutoff,
                double_column=double_column,
            )
            st.pyplot(fig_d)
            _cache_fig("fig_s3_volcano_enrichment", fig_d)

            col_pdf, col_svg = st.columns(2)
            with col_pdf:
                st.download_button(
                    "Download PDF",
                    st.session_state.fig_bytes["fig_s3_volcano_enrichment"]["pdf"],
                    "fig_s3_volcano.pdf", "application/pdf",
                    key="dl_fig_s3_pdf",
                )
            with col_svg:
                st.download_button(
                    "Download SVG",
                    st.session_state.fig_bytes["fig_s3_volcano_enrichment"]["svg"],
                    "fig_s3_volcano.svg", "image/svg+xml",
                    key="dl_fig_s3_svg",
                )
            plt.close(fig_d)

            # Dotplot — use correlation-ranked clusters (preferring dual-significant)
            st.subheader("Supplementary Figure S4: Enriched Gene Dot Plot")
            if len(corr_df) > 0 and "both_significant" in corr_df.columns:
                _sig_dp = corr_df.loc[corr_df["both_significant"], "cluster"].tolist()
                dotplot_clusters = _sig_dp[:15] if len(_sig_dp) >= 5 else corr_df["cluster"].tolist()[:15]
            else:
                dotplot_clusters = corr_df["cluster"].tolist()[:15] if len(corr_df) > 0 else []
            dotplot_genes = top_enriched_genes[:20]
            fig_c = figure_dotplot(
                enriched_mean_expr, frac_expr,
                dotplot_genes, dotplot_clusters,
                double_column=True,
            )
            st.pyplot(fig_c)
            _cache_fig("fig_s4_dotplot", fig_c)

            col_pdf, col_svg = st.columns(2)
            with col_pdf:
                st.download_button(
                    "Download PDF",
                    st.session_state.fig_bytes["fig_s4_dotplot"]["pdf"],
                    "fig_s4_dotplot.pdf", "application/pdf",
                    key="dl_fig_s4_pdf",
                )
            with col_svg:
                st.download_button(
                    "Download SVG",
                    st.session_state.fig_bytes["fig_s4_dotplot"]["svg"],
                    "fig_s4_dotplot.svg", "image/svg+xml",
                    key="dl_fig_s4_svg",
                )
            plt.close(fig_c)

            st.download_button(
                "Download Fisher's test results (CSV)",
                fisher_df.to_csv(index=False).encode(),
                _filtered_name("fisher_test_results.csv"), "text/csv",
                key="dl_fisher_csv_tab4",
            )
        else:
            st.warning("No marker overlap results to display.")

    # ======================================================================
    # TAB 6: Gene Heatmap (Supplementary)
    # ======================================================================
    with tab6:
        st.header("Supplementary: Gene Expression Heatmap")
        st.markdown(
            "Z-scored mean expression of top bacTRAP-enriched genes across "
            "the highest-correlating HypoMap clusters."
            + (
                f" When the baseline {sanity_gene} filter is active, the "
                f"cluster columns and the *z*-scoring reference are "
                f"recomputed against the filtered cluster set so the "
                f"displayed *z*-scores stay internally consistent with "
                f"the shown cluster panel."
                if baseline_allowed is not None else ""
            )
        )

        if not zscore_df.empty:
            st.subheader("Supplementary Figure S5: Heatmap")
            fig_e = figure_heatmap(zscore_df, double_column=double_column)
            st.pyplot(fig_e)
            _cache_fig("fig_s5_heatmap", fig_e)

            col_pdf, col_svg = st.columns(2)
            with col_pdf:
                st.download_button(
                    "Download PDF",
                    st.session_state.fig_bytes["fig_s5_heatmap"]["pdf"],
                    "fig_s5_heatmap.pdf", "application/pdf",
                    key="dl_fig_s5_pdf",
                )
            with col_svg:
                st.download_button(
                    "Download SVG",
                    st.session_state.fig_bytes["fig_s5_heatmap"]["svg"],
                    "fig_s5_heatmap.svg", "image/svg+xml",
                    key="dl_fig_s5_svg",
                )
            plt.close(fig_e)
        else:
            st.warning("No heatmap data available with current parameters.")

    # ======================================================================
    # TAB 7: NNLS Deconvolution (Supplementary)
    # ======================================================================
    with tab7:
        st.header("Supplementary: NNLS Deconvolution")
        st.markdown(
            "Non-negative least squares: find cluster weights that best "
            "reconstruct the bacTRAP enrichment profile from cluster-level "
            "expression signatures."
        )

        if len(nnls_df) > 0:
            nonzero = nnls_df[nnls_df["weight"] > 1e-6]
            st.metric("Clusters with non-zero weight", len(nonzero))

            st.subheader("NNLS Weights")
            st.dataframe(
                nnls_df[nnls_df["weight"] > 1e-6][["cluster", "weight", "weight_norm"]].style.format({
                    "weight": "{:.4f}",
                    "weight_norm": "{:.4f}",
                }),
                use_container_width=True,
            )

            st.subheader("Supplementary Figure S6: NNLS Deconvolution")
            fig_f = figure_nnls_barplot(nnls_df, top_n=20, double_column=double_column)
            st.pyplot(fig_f)
            _cache_fig("fig_s6_nnls", fig_f)

            col_pdf, col_svg = st.columns(2)
            with col_pdf:
                st.download_button(
                    "Download PDF",
                    st.session_state.fig_bytes["fig_s6_nnls"]["pdf"],
                    "fig_s6_nnls.pdf", "application/pdf",
                    key="dl_fig_s6_pdf",
                )
            with col_svg:
                st.download_button(
                    "Download SVG",
                    st.session_state.fig_bytes["fig_s6_nnls"]["svg"],
                    "fig_s6_nnls.svg", "image/svg+xml",
                    key="dl_fig_s6_svg",
                )
            plt.close(fig_f)

            st.download_button(
                "Download NNLS results (CSV)",
                nnls_df.to_csv(index=False).encode(),
                _filtered_name("nnls_results.csv"), "text/csv",
                key="dl_nnls_csv",
            )
        else:
            st.warning("NNLS deconvolution produced no results.")

    # ======================================================================
    # TAB 8: GSEA (Supplementary)
    # ======================================================================
    with tab8:
        st.header("Supplementary: GSEA Preranked Enrichment")
        st.markdown(
            "All matched genes are ranked by bacTRAP log₂FC. For each cluster's "
            "marker gene set, a running enrichment score is computed — more "
            "powerful than Fisher's binary overlap test because it uses the "
            "full ranking."
        )

        if len(gsea_df) > 0:
            st.subheader("Enrichment Results")
            st.dataframe(
                gsea_df[["cluster", "ES", "NES", "pvalue", "padj", "n_hits"]].style.format({
                    "ES": "{:.4f}",
                    "NES": "{:.4f}",
                    "pvalue": "{:.2e}",
                    "padj": "{:.2e}",
                }),
                use_container_width=True,
            )

            st.subheader("Supplementary Figure S7: Enrichment Curves (top 5)")
            fig_g = figure_gsea_curves(
                gsea_df, gsea_running_scores, gsea_ranked_genes,
                top_n=5, double_column=double_column,
            )
            st.pyplot(fig_g)
            _cache_fig("fig_s7_gsea_curves", fig_g)

            col_pdf, col_svg = st.columns(2)
            with col_pdf:
                st.download_button(
                    "Download PDF",
                    st.session_state.fig_bytes["fig_s7_gsea_curves"]["pdf"],
                    "fig_s7_gsea_curves.pdf", "application/pdf",
                    key="dl_fig_s7_pdf",
                )
            with col_svg:
                st.download_button(
                    "Download SVG",
                    st.session_state.fig_bytes["fig_s7_gsea_curves"]["svg"],
                    "fig_s7_gsea_curves.svg", "image/svg+xml",
                    key="dl_fig_s7_svg",
                )
            plt.close(fig_g)

            st.subheader("Supplementary Figure S8: NES Barplot")
            fig_h = figure_gsea_barplot(gsea_df, top_n=20, double_column=double_column)
            st.pyplot(fig_h)
            _cache_fig("fig_s8_gsea_barplot", fig_h)

            col_pdf, col_svg = st.columns(2)
            with col_pdf:
                st.download_button(
                    "Download PDF",
                    st.session_state.fig_bytes["fig_s8_gsea_barplot"]["pdf"],
                    "fig_s8_gsea_barplot.pdf", "application/pdf",
                    key="dl_fig_s8_pdf",
                )
            with col_svg:
                st.download_button(
                    "Download SVG",
                    st.session_state.fig_bytes["fig_s8_gsea_barplot"]["svg"],
                    "fig_s8_gsea_barplot.svg", "image/svg+xml",
                    key="dl_fig_s8_svg",
                )
            plt.close(fig_h)

            st.download_button(
                "Download GSEA results (CSV)",
                gsea_df.to_csv(index=False).encode(),
                _filtered_name("gsea_results.csv"), "text/csv",
                key="dl_gsea_csv",
            )
        else:
            st.warning("No GSEA results to display.")

    # ======================================================================
    # TAB: Export
    # ======================================================================
    with tab_export:
        st.header("Export All Results")

        st.subheader("Figures")
        st.markdown("Download all figures as a ZIP archive (PDF + SVG).")

        # Build ZIP from cached figure bytes (no regeneration needed)
        cached_bytes = st.session_state.get("fig_bytes", {})
        cached_tables = st.session_state.get("table_bytes", {})
        if cached_bytes:
            buf = io.BytesIO()
            with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
                for name, fmt_dict in cached_bytes.items():
                    zf.writestr(f"{name}.pdf", fmt_dict["pdf"])
                    zf.writestr(f"{name}.svg", fmt_dict["svg"])
                # Include AUCell raw-data tables alongside the figures.
                # Skip private metadata entries (e.g. cache signatures) that
                # share this dict but aren't serialised CSV payloads.
                for tbl_name, tbl_bytes in cached_tables.items():
                    if tbl_name.startswith("_") or not isinstance(
                        tbl_bytes, (bytes, bytearray, memoryview)
                    ):
                        continue
                    zf.writestr(f"{tbl_name}.csv", tbl_bytes)
                # Include log file (best-effort; skip if unreadable).
                if _LOG_FILE.is_file():
                    try:
                        zf.writestr(
                            "bactrap_hypomap.log",
                            _LOG_FILE.read_text(errors="replace"),
                        )
                    except OSError:
                        logger.exception(
                            "Failed to include log file in export ZIP",
                        )
            buf.seek(0)

            st.download_button(
                "Download All Figures (ZIP)",
                buf.getvalue(),
                "bactrap_hypomap_figures.zip",
                "application/zip",
                use_container_width=True,
                key="dl_all_figs_zip",
            )
        else:
            st.info("Run the analysis first to generate figures.")

        st.markdown("---")
        st.subheader("Tables")

        col_t1, col_t2 = st.columns(2)
        with col_t1:
            st.download_button(
                "Correlation results (CSV)",
                corr_df.to_csv(index=False).encode(),
                _filtered_name("correlation_results.csv"), "text/csv",
                key="dl_corr_csv_export",
            )
        with col_t2:
            if len(fisher_df) > 0:
                st.download_button(
                    "Fisher's test results (CSV)",
                    fisher_df.to_csv(index=False).encode(),
                    _filtered_name("fisher_test_results.csv"), "text/csv",
                    key="dl_fisher_csv_export",
                )

        col_t3, col_t4 = st.columns(2)
        with col_t3:
            st.download_button(
                "Matched genes (CSV)",
                bactrap_matched.to_csv(index=False).encode(),
                "matched_genes.csv", "text/csv",
                key="dl_matched_csv",
            )
        with col_t4:
            if len(enriched_df) > 0:
                st.download_button(
                    "Enriched genes (CSV)",
                    enriched_df.to_csv(index=False).encode(),
                    "enriched_genes.csv", "text/csv",
                    key="dl_enriched_csv",
                )

        if not zscore_df.empty:
            st.download_button(
                "Z-score heatmap data (CSV)",
                zscore_df.to_csv().encode(),
                _filtered_name("zscore_heatmap.csv"), "text/csv",
                key="dl_zscore_csv",
            )

        col_t5, col_t6 = st.columns(2)
        with col_t5:
            if len(nnls_df) > 0:
                st.download_button(
                    "NNLS results (CSV)",
                    nnls_df.to_csv(index=False).encode(),
                    _filtered_name("nnls_results.csv"), "text/csv",
                    key="dl_nnls_csv_export",
                )
        with col_t6:
            if len(gsea_df) > 0:
                st.download_button(
                    "GSEA results (CSV)",
                    gsea_df.to_csv(index=False).encode(),
                    _filtered_name("gsea_results.csv"), "text/csv",
                    key="dl_gsea_csv_export",
                )

        if len(composite_df) > 0:
            st.download_button(
                "Composite ranking (CSV)",
                composite_df.to_csv(index=False).encode(),
                _filtered_name("composite_ranking.csv"), "text/csv",
                key="dl_composite_csv_export",
            )

        # AUCell raw data — pre-encoded in table_bytes (populated during run)
        aucell_table_bytes = st.session_state.get("table_bytes", {})
        if "aucell_per_cell" in aucell_table_bytes:
            col_a1, col_a2 = st.columns(2)
            with col_a1:
                st.download_button(
                    "AUCell per-cell scores (CSV)",
                    aucell_table_bytes["aucell_per_cell"],
                    _filtered_name("aucell_per_cell.csv"), "text/csv",
                    key="dl_aucell_per_cell_export",
                    help=(
                        "cell_id, cluster, aucell_score — raw data for figures 1a, 1c and "
                        "supplementary S3. Restricted to clusters passing the baseline "
                        "filter when active."
                    ),
                )
            with col_a2:
                st.download_button(
                    "AUCell per-cluster summary (CSV)",
                    aucell_table_bytes["aucell_per_cluster"],
                    _filtered_name("aucell_per_cluster.csv"), "text/csv",
                    key="dl_aucell_per_cluster_export",
                    help=(
                        "cluster, n_cells, mean, median, std, sem, t_stat, pvalue, "
                        "qvalue, significant — raw data for figures 1b/1c and supplementary S2."
                    ),
                )

        st.markdown("---")
        st.subheader("Diagnostics")
        if _LOG_FILE.is_file():
            try:
                _log_bytes = _LOG_FILE.read_text(errors="replace").encode()
            except OSError as e:
                # File exists but can't be read (permissions, lock, disk
                # fault) — surface the reason rather than handing the user
                # an empty download.
                logger.exception("Failed to read log file %s", _LOG_FILE)
                st.info(
                    f"Log file exists at `{_LOG_FILE}` but could not be "
                    f"read: {e}. Check file permissions."
                )
            else:
                st.download_button(
                    "Download log file",
                    _log_bytes,
                    "bactrap_hypomap.log", "text/plain",
                    use_container_width=True,
                    key="dl_log_file",
                )
        else:
            st.info("No log file generated yet.")

else:
    st.info("Click **Run Analysis** in the sidebar to start.")
