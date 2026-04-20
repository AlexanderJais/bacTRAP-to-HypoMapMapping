"""
bacTRAP-to-HypoMap Mapping Tool

Streamlit application for mapping bacTRAP bulk RNA-seq data onto the
murine HypoMap single-cell atlas (Steuernagel et al., Nature Metabolism 2022).

Run with: streamlit run app.py
"""

import io
import logging
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
)
log2fc_cutoff = st.sidebar.slider(
    "log₂FC cutoff", 0.0, 5.0, 1.0, 0.25,
)
top_n_genes = st.sidebar.slider(
    "Top N genes for scoring", 10, 500, 50, 10,
)
n_markers_per_cluster = st.sidebar.slider(
    "Marker genes per cluster", 20, 500, 100, 10,
)
min_cells_per_cluster = st.sidebar.slider(
    "Min cells per cluster", 1, 100, 10, 1,
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
    help="Subsample cells for UMAP visualization to reduce rendering time.",
)

st.sidebar.markdown("---")
st.sidebar.subheader("Figure Settings")
fig_width_mode = st.sidebar.radio(
    "Figure width", ["Single column (89mm)", "Double column (183mm)"],
    index=0,
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
    get_gene_names_from_adata,
    _detect_gene_column,
    _build_adata_gene_lookup,
    _looks_like_ensembl,
)
from analysis import (
    compute_enrichment_correlation,
    get_enriched_genes,
    compute_marker_genes,
    load_precomputed_markers,
    fisher_overlap_test,
    compute_enrichment_score,
    compute_zscore_heatmap_data,
    compute_nnls_deconvolution,
    compute_gsea_enrichment,
    compute_aucell_scores,
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
    figure_aucell_cluster_barplot,
    figure_aucell_violins,
    figure_aucell_histogram,
    figure_composite_ranking,
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
    help="Select the cell-type annotation level from HypoMap .obs.",
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
        "Auto-detected based on which column yields the most matches against HypoMap."
    ),
)
# Map UI selection to the internal value expected by match_genes
_gene_col_for_matching = "_index" if gene_col_selection == "(use row index)" else gene_col_selection

# ---------------------------------------------------------------------------
# Run analysis
# ---------------------------------------------------------------------------

# Progress bar placeholder — rendered above tabs so it's always visible
progress_placeholder = st.empty()

tab1, tab_aucell, tab3, tab4, tab5, tab6, tab7, tab8, tab_export = st.tabs([
    "📊 Data Overview",
    "⭐ AUCell (Main Figure)",
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
    padj_cutoff, log2fc_cutoff, top_n_genes,
    n_markers_per_cluster, min_cells_per_cluster, marker_method,
    umap_subsample,
)

if run_button or st.session_state.analysis_done:

    # Reuse cached results on rerun if parameters haven't changed
    _cached = st.session_state.get("_analysis_cache")
    _need_recompute = run_button or _cached is None or _cached.get("params") != _analysis_params

    if _need_recompute:
        # ---- Gene matching ----
        progress = progress_placeholder.progress(0, text="Matching genes...")

        bactrap_matched, matched_genes, gene_to_idx, matched_in_raw = match_genes(
            bactrap_df, adata, gene_col=_gene_col_for_matching,
            _prebuilt_lookup=(_adata_lookup, _adata_gnames, _adata_has_raw),
        )

        if len(matched_genes) == 0:
            progress.empty()
            st.error("No genes could be matched between the bacTRAP data and HypoMap. "
                     "Check that gene_name symbols in your bacTRAP file correspond to "
                     "gene names in the HypoMap atlas.")
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
        enriched_df = get_enriched_genes(bactrap_matched, padj_cutoff, log2fc_cutoff)
        enriched_genes_list = enriched_df["_hypomap_gene_name"].tolist()

        # Top N enriched genes ranked by both significance and effect size
        enriched_sorted = enriched_df.sort_values("log2FoldChange", ascending=False).reset_index(drop=True)
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
                    markers = compute_marker_genes(
                        adata, annotation_col,
                        n_genes=n_markers_per_cluster,
                        min_cells=min_cells_per_cluster,
                        method=marker_method,
                    )
            st.session_state["_markers_cache"] = {
                "params": _markers_params, "markers": markers,
            }
        progress.progress(45, text="Running Fisher's exact test...")

        universe_size = len(matched_genes)
        fisher_df = fisher_overlap_test(enriched_genes_list, markers, universe_size)
        progress.progress(50, text="Computing UMAP enrichment scores...")

        # ---- UMAP enrichment score ----
        if len(top_enriched_genes) == 0:
            st.warning(
                f"No genes pass enrichment thresholds (padj < {padj_cutoff}, "
                f"log₂FC > {log2fc_cutoff}). UMAP enrichment score will be zero. "
                "Try relaxing the cutoffs."
            )
        enrichment_scores = compute_enrichment_score(adata, top_enriched_genes)
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
            st.warning(f"GSEA computation failed: {e}")
            gsea_df = pd.DataFrame()
            gsea_running_scores = {}
            gsea_ranked_genes = np.array([])

        progress.progress(80, text="Computing AUCell scores...")

        # ---- AUCell scoring ----
        aucell_scores = compute_aucell_scores(adata, top_enriched_genes)

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
            "enrichment_scores": enrichment_scores,
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
        enrichment_scores = _c["enrichment_scores"]
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
        composite_df = _c["composite_df"]
        sub_indices = _c["sub_indices"]
        umap_coords = _c["umap_coords"]
        cell_labels = _c["cell_labels"]
        enriched_sorted = _c["enriched_sorted"]
        top_genes_heatmap = _c["top_genes_heatmap"]
        top_clusters_heatmap = _c["top_clusters_heatmap"]
        progress_placeholder.empty()

    st.session_state.analysis_done = True

    # Cache figure bytes so the Export tab doesn't regenerate them.
    # Only reset when a new analysis run is triggered (run_button pressed
    # or parameters changed), not on every Streamlit rerun.
    if _need_recompute:
        st.session_state.fig_bytes = {}

    if "fig_bytes" not in st.session_state:
        st.session_state.fig_bytes = {}

    def _cache_fig(name, fig):
        st.session_state.fig_bytes[name] = {
            "pdf": fig_to_bytes(fig, "pdf"),
            "svg": fig_to_bytes(fig, "svg"),
        }

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

        # Figure 1a: AUCell UMAP
        st.subheader("Figure 1a: AUCell Enrichment UMAP")
        fig_1a = figure_aucell_umap(
            umap_coords, aucell_scores,
            double_column=double_column,
            subsample_idx=sub_indices,
        )
        st.pyplot(fig_1a)
        _cache_fig("fig_1a_aucell_umap", fig_1a)

        col_pdf, col_svg = st.columns(2)
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
        plt.close(fig_1a)

        # Figure 1b: AUCell Cluster Barplot
        st.subheader("Figure 1b: AUCell Score per Cluster")
        st.markdown(
            "Mean AUCell score per cluster (error bars = SEM). "
            "Clusters where cells consistently express the enriched gene set "
            "rank highest."
        )
        fig_1b = figure_aucell_cluster_barplot(
            aucell_scores, cell_labels,
            top_n=25, double_column=double_column,
        )
        st.pyplot(fig_1b)
        _cache_fig("fig_1b_aucell_barplot", fig_1b)

        col_pdf, col_svg = st.columns(2)
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
        plt.close(fig_1b)

        # Figure 1c: AUCell Violin Plots
        st.subheader("Figure 1c: AUCell Score Distributions (Top Clusters)")
        st.markdown(
            "Violin plots showing the full distribution of AUCell scores within "
            "each top-ranked cluster. Solid line = mean, dashed = median."
        )
        fig_1c = figure_aucell_violins(
            aucell_scores, cell_labels,
            top_n=15, double_column=True,
        )
        st.pyplot(fig_1c)
        _cache_fig("fig_1c_aucell_violins", fig_1c)

        col_pdf, col_svg = st.columns(2)
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
        plt.close(fig_1c)

        # Figure 1d: AUCell Score Histogram
        st.subheader("Figure 1d: AUCell Score Distribution (All Cells)")
        st.markdown(
            "Global distribution of AUCell scores with percentile markers. "
            "Cells above the 95th percentile are most likely part of the "
            "bacTRAP target population."
        )
        fig_1d = figure_aucell_histogram(
            aucell_scores, double_column=double_column,
        )
        st.pyplot(fig_1d)
        _cache_fig("fig_1d_aucell_histogram", fig_1d)

        col_pdf, col_svg = st.columns(2)
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
        plt.close(fig_1d)

        # Figure 1e: Composite Consensus Ranking
        st.markdown("---")
        st.subheader("Figure 1e: Composite Consensus Ranking")
        st.markdown(
            "Validation: consensus ranking across all methods confirms AUCell "
            "results. Each method's scores are converted to percentile ranks "
            "and averaged."
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
                "correlation_results.csv", "text/csv",
                key="dl_corr_csv_tab2",
            )
        else:
            st.warning("No correlation results to display.")

    # ======================================================================
    # TAB 4: UMAP Projection (Supplementary)
    # ======================================================================
    with tab4:
        st.header("Supplementary: UMAP Enrichment Projection")
        st.markdown(
            f"bacTRAP enrichment score projected onto the HypoMap UMAP. "
            f"Score computed from the top **{len(top_enriched_genes)}** enriched "
            f"genes (padj < {padj_cutoff}, log₂FC > {log2fc_cutoff})."
        )

        st.subheader("Supplementary Figure S2: UMAP Enrichment Map")
        fig_b = figure_umap_enrichment(
            umap_coords=umap_coords,
            cell_labels=cell_labels,
            enrichment_scores=enrichment_scores,
            double_column=True,
            point_size=0.3,
            subsample_idx=sub_indices,
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
                "fisher_test_results.csv", "text/csv",
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
                "nnls_results.csv", "text/csv",
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
                "gsea_results.csv", "text/csv",
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
        if cached_bytes:
            buf = io.BytesIO()
            with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
                for name, fmt_dict in cached_bytes.items():
                    zf.writestr(f"{name}.pdf", fmt_dict["pdf"])
                    zf.writestr(f"{name}.svg", fmt_dict["svg"])
                # Include log file
                if _LOG_FILE.is_file():
                    zf.writestr("bactrap_hypomap.log", _LOG_FILE.read_text(errors="replace"))
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
                "correlation_results.csv", "text/csv",
                key="dl_corr_csv_export",
            )
        with col_t2:
            if len(fisher_df) > 0:
                st.download_button(
                    "Fisher's test results (CSV)",
                    fisher_df.to_csv(index=False).encode(),
                    "fisher_test_results.csv", "text/csv",
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
                "zscore_heatmap.csv", "text/csv",
                key="dl_zscore_csv",
            )

        col_t5, col_t6 = st.columns(2)
        with col_t5:
            if len(nnls_df) > 0:
                st.download_button(
                    "NNLS results (CSV)",
                    nnls_df.to_csv(index=False).encode(),
                    "nnls_results.csv", "text/csv",
                    key="dl_nnls_csv_export",
                )
        with col_t6:
            if len(gsea_df) > 0:
                st.download_button(
                    "GSEA results (CSV)",
                    gsea_df.to_csv(index=False).encode(),
                    "gsea_results.csv", "text/csv",
                    key="dl_gsea_csv_export",
                )

        if len(composite_df) > 0:
            st.download_button(
                "Composite ranking (CSV)",
                composite_df.to_csv(index=False).encode(),
                "composite_ranking.csv", "text/csv",
                key="dl_composite_csv_export",
            )

        st.markdown("---")
        st.subheader("Diagnostics")
        if _LOG_FILE.is_file():
            st.download_button(
                "Download log file",
                _LOG_FILE.read_text(errors="replace").encode(),
                "bactrap_hypomap.log", "text/plain",
                use_container_width=True,
                key="dl_log_file",
            )
        else:
            st.info("No log file generated yet.")

else:
    st.info("Click **Run Analysis** in the sidebar to start.")
