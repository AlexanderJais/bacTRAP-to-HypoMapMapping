"""
bacTRAP-to-HypoMap Mapping Tool

Streamlit application for mapping bacTRAP bulk RNA-seq data onto the
murine HypoMap single-cell atlas (Steuernagel et al., Nature Metabolism 2022).

Run with: streamlit run app.py
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

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

# Check file inputs
files_ready = (
    bactrap_file.strip() != ""
    and hypomap_file.strip() != ""
    and Path(bactrap_file.strip()).exists()
    and Path(hypomap_file.strip()).exists()
)

if not files_ready and (bactrap_file.strip() or hypomap_file.strip()):
    missing = []
    if bactrap_file.strip() and not Path(bactrap_file.strip()).exists():
        missing.append(f"bacTRAP file not found: `{bactrap_file.strip()}`")
    if hypomap_file.strip() and not Path(hypomap_file.strip()).exists():
        missing.append(f"HypoMap file not found: `{hypomap_file.strip()}`")
    if missing:
        for m in missing:
            st.warning(m)

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
)
from analysis import (
    compute_enrichment_correlation,
    get_enriched_genes,
    compute_marker_genes,
    load_precomputed_markers,
    fisher_overlap_test,
    compute_enrichment_score,
    compute_zscore_heatmap_data,
)
from figures import (
    setup_nature_style,
    figure_correlation_barplot,
    figure_umap_enrichment,
    figure_dotplot,
    figure_volcano_enrichment,
    figure_heatmap,
    fig_to_bytes,
)

# Load data with caching
with st.spinner("Loading data..."):
    bactrap_df = load_bactrap(bactrap_file.strip())
    adata = load_hypomap(hypomap_file.strip())

# Annotation column selection
ann_cols = get_annotation_columns(adata)
if len(ann_cols) == 0:
    st.error("No suitable annotation columns found in HypoMap .obs.")
    st.stop()

default_idx = 0
# Try to find a reasonable default
for i, col in enumerate(ann_cols):
    col_lower = col.lower()
    if "cell_type" in col_lower or "celltype" in col_lower or "cluster" in col_lower:
        default_idx = i
        break

annotation_col = st.sidebar.selectbox(
    "Annotation column",
    ann_cols,
    index=default_idx,
    help="Select the cell-type annotation level from HypoMap .obs.",
)

# ---------------------------------------------------------------------------
# Tabs
# ---------------------------------------------------------------------------

tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
    "📊 Data Overview",
    "📈 Correlation Analysis",
    "🗺️ UMAP Projection",
    "🔬 Marker Overlap",
    "🔥 Gene Heatmap",
    "📥 Export",
])


# ---------------------------------------------------------------------------
# Run analysis
# ---------------------------------------------------------------------------

if run_button or st.session_state.analysis_done:

    # ---- Gene matching ----
    progress = st.progress(0, text="Matching genes...")

    bactrap_matched, matched_genes, gene_to_idx = match_genes(bactrap_df, adata)

    if len(matched_genes) == 0:
        progress.empty()
        st.error("No genes could be matched between the bacTRAP data and HypoMap. "
                 "Check that gene_name symbols in your bacTRAP file correspond to "
                 "gene names in the HypoMap atlas.")
        st.stop()

    progress.progress(10, text="Genes matched. Computing cluster means...")

    # ---- Cluster mean expression ----
    gene_indices = [gene_to_idx[g] for g in matched_genes]
    cluster_mean_expr = compute_cluster_mean_expression(
        adata, gene_indices, annotation_col, min_cells=min_cells_per_cluster,
    )
    progress.progress(25, text="Cluster means computed. Identifying enriched genes...")

    # ---- Enriched genes ----
    enriched_df = get_enriched_genes(bactrap_matched, padj_cutoff, log2fc_cutoff)
    enriched_genes_list = enriched_df["_hypomap_gene_name"].tolist()

    # Top N enriched genes ranked by both significance and effect size
    enriched_sorted = enriched_df.sort_values("log2FoldChange", ascending=False)
    top_enriched_genes = enriched_sorted["_hypomap_gene_name"].tolist()[:top_n_genes]

    progress.progress(30, text="Computing enrichment correlation...")

    # ---- Correlation analysis ----
    corr_df = compute_enrichment_correlation(bactrap_matched, cluster_mean_expr)
    progress.progress(45, text="Computing marker gene overlap...")

    # ---- Marker gene overlap ----
    # Try pre-computed markers first
    markers = load_precomputed_markers(adata)
    if markers is None:
        with st.spinner("Computing marker genes (this may take several minutes)..."):
            markers = compute_marker_genes(
                adata, annotation_col,
                n_genes=n_markers_per_cluster,
                min_cells=min_cells_per_cluster,
            )
    progress.progress(65, text="Running Fisher's exact test...")

    universe_size = len(matched_genes)
    fisher_df = fisher_overlap_test(enriched_genes_list, markers, universe_size)
    progress.progress(75, text="Computing UMAP enrichment scores...")

    # ---- UMAP enrichment score ----
    if len(top_enriched_genes) == 0:
        st.warning(
            f"No genes pass enrichment thresholds (padj < {padj_cutoff}, "
            f"log₂FC > {log2fc_cutoff}). UMAP enrichment score will be zero. "
            "Try relaxing the cutoffs."
        )
    enrichment_scores = compute_enrichment_score(adata, top_enriched_genes)
    progress.progress(85, text="Preparing figures...")

    # ---- Fraction expressing for dotplot ----
    enriched_gene_indices = [gene_to_idx[g] for g in top_enriched_genes if g in gene_to_idx]
    n_dropped = len(top_enriched_genes) - len(enriched_gene_indices)
    if n_dropped > 0:
        st.warning(f"{n_dropped} enriched gene(s) could not be mapped back to HypoMap indices and were excluded from the dot plot.")
    top_clusters_corr = corr_df["cluster"].tolist()[:15] if len(corr_df) > 0 else []

    frac_expr = compute_fraction_expressing(
        adata, enriched_gene_indices, annotation_col,
        min_cells=min_cells_per_cluster,
    )
    enriched_mean_expr = compute_cluster_mean_expression(
        adata, enriched_gene_indices, annotation_col,
        min_cells=min_cells_per_cluster,
    )

    # ---- Z-score heatmap data ----
    top_clusters_heatmap = corr_df["cluster"].tolist()[:20] if len(corr_df) > 0 else []
    top_genes_heatmap = top_enriched_genes[:30]
    zscore_df = compute_zscore_heatmap_data(
        cluster_mean_expr, top_genes_heatmap, top_clusters_heatmap,
    )

    progress.progress(95, text="Generating figures...")

    # ---- Subsample for UMAP ----
    sub_indices = None
    if adata.n_obs > umap_subsample:
        rng = np.random.default_rng(42)
        sub_indices = np.sort(rng.choice(adata.n_obs, size=umap_subsample, replace=False))

    # Get UMAP coordinates
    umap_key = None
    for key in ["X_umap", "X_UMAP"]:
        if key in adata.obsm:
            umap_key = key
            break
    if umap_key is None:
        st.error("No UMAP coordinates found in HypoMap .obsm. Expected 'X_umap'.")
        st.stop()

    umap_coords = adata.obsm[umap_key]
    cell_labels = adata.obs[annotation_col].values.astype(str)

    progress.progress(100, text="Analysis complete!")
    st.session_state.analysis_done = True

    # Cache figure bytes so the Export tab doesn't regenerate them.
    # Only reset when a new analysis run is triggered (run_button pressed),
    # not on every Streamlit rerun.
    if run_button:
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

        st.subheader("Gene Matching Summary")

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
            display_cols = ["gene_name", "log2FoldChange", "padj", "IP", "Input"]
            available_cols = [c for c in display_cols if c in enriched_df.columns]
            st.dataframe(
                enriched_sorted[available_cols].head(30).reset_index(drop=True),
                use_container_width=True,
            )
        else:
            st.warning("No genes pass the current enrichment thresholds.")

    # ======================================================================
    # TAB 2: Correlation Analysis
    # ======================================================================
    with tab2:
        st.header("Enrichment Correlation Analysis")
        st.markdown(
            "Pearson and Spearman correlation between the bacTRAP log₂FC enrichment "
            "profile and mean expression per HypoMap cluster."
        )

        if len(corr_df) > 0:
            st.subheader("Top Correlated Clusters")
            st.dataframe(
                corr_df.style.format({
                    "pearson_r": "{:.4f}",
                    "pearson_pval": "{:.2e}",
                    "spearman_r": "{:.4f}",
                    "spearman_pval": "{:.2e}",
                }),
                use_container_width=True,
            )

            st.subheader("Figure A: Correlation Barplot")
            fig_a = figure_correlation_barplot(corr_df, top_n=20, double_column=double_column)
            st.pyplot(fig_a)
            _cache_fig("fig_a_correlation_barplot", fig_a)

            col_pdf, col_svg = st.columns(2)
            with col_pdf:
                st.download_button(
                    "Download PDF",
                    st.session_state.fig_bytes["fig_a_correlation_barplot"]["pdf"],
                    "fig_a_correlation.pdf", "application/pdf",
                    key="dl_fig_a_pdf",
                )
            with col_svg:
                st.download_button(
                    "Download SVG",
                    st.session_state.fig_bytes["fig_a_correlation_barplot"]["svg"],
                    "fig_a_correlation.svg", "image/svg+xml",
                    key="dl_fig_a_svg",
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
    # TAB 3: UMAP Projection
    # ======================================================================
    with tab3:
        st.header("UMAP Enrichment Projection")
        st.markdown(
            f"bacTRAP enrichment score projected onto the HypoMap UMAP. "
            f"Score computed from the top **{len(top_enriched_genes)}** enriched "
            f"genes (padj < {padj_cutoff}, log₂FC > {log2fc_cutoff})."
        )

        st.subheader("Figure B: UMAP Enrichment Map")
        fig_b = figure_umap_enrichment(
            umap_coords=umap_coords,
            cell_labels=cell_labels,
            enrichment_scores=enrichment_scores,
            double_column=True,
            point_size=0.3,
            subsample_idx=sub_indices,
        )
        st.pyplot(fig_b)
        _cache_fig("fig_b_umap_enrichment", fig_b)

        col_pdf, col_svg = st.columns(2)
        with col_pdf:
            st.download_button(
                "Download PDF",
                st.session_state.fig_bytes["fig_b_umap_enrichment"]["pdf"],
                "fig_b_umap.pdf", "application/pdf",
                key="dl_fig_b_pdf",
            )
        with col_svg:
            st.download_button(
                "Download SVG",
                st.session_state.fig_bytes["fig_b_umap_enrichment"]["svg"],
                "fig_b_umap.svg", "image/svg+xml",
                key="dl_fig_b_svg",
            )
        plt.close(fig_b)

    # ======================================================================
    # TAB 4: Marker Overlap
    # ======================================================================
    with tab4:
        st.header("Marker Gene Overlap Analysis")
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

            st.subheader("Figure D: Enrichment Volcano Plot")
            fig_d = figure_volcano_enrichment(
                fisher_df, pval_threshold=padj_cutoff,
                double_column=double_column,
            )
            st.pyplot(fig_d)
            _cache_fig("fig_d_volcano_enrichment", fig_d)

            col_pdf, col_svg = st.columns(2)
            with col_pdf:
                st.download_button(
                    "Download PDF",
                    st.session_state.fig_bytes["fig_d_volcano_enrichment"]["pdf"],
                    "fig_d_volcano.pdf", "application/pdf",
                    key="dl_fig_d_pdf",
                )
            with col_svg:
                st.download_button(
                    "Download SVG",
                    st.session_state.fig_bytes["fig_d_volcano_enrichment"]["svg"],
                    "fig_d_volcano.svg", "image/svg+xml",
                    key="dl_fig_d_svg",
                )
            plt.close(fig_d)

            # Dotplot
            st.subheader("Figure C: Marker Overlap Dot Plot")
            top_fisher_clusters = fisher_df["cluster"].tolist()[:10]
            dotplot_genes = top_enriched_genes[:20]
            fig_c = figure_dotplot(
                enriched_mean_expr, frac_expr,
                dotplot_genes, top_fisher_clusters,
                double_column=True,
            )
            st.pyplot(fig_c)
            _cache_fig("fig_c_dotplot", fig_c)

            col_pdf, col_svg = st.columns(2)
            with col_pdf:
                st.download_button(
                    "Download PDF",
                    st.session_state.fig_bytes["fig_c_dotplot"]["pdf"],
                    "fig_c_dotplot.pdf", "application/pdf",
                    key="dl_fig_c_pdf",
                )
            with col_svg:
                st.download_button(
                    "Download SVG",
                    st.session_state.fig_bytes["fig_c_dotplot"]["svg"],
                    "fig_c_dotplot.svg", "image/svg+xml",
                    key="dl_fig_c_svg",
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
    # TAB 5: Gene Heatmap
    # ======================================================================
    with tab5:
        st.header("Gene Expression Heatmap")
        st.markdown(
            "Z-scored mean expression of top bacTRAP-enriched genes across "
            "the highest-correlating HypoMap clusters."
        )

        if not zscore_df.empty:
            st.subheader("Figure E: Heatmap")
            fig_e = figure_heatmap(zscore_df, double_column=double_column)
            st.pyplot(fig_e)
            _cache_fig("fig_e_heatmap", fig_e)

            col_pdf, col_svg = st.columns(2)
            with col_pdf:
                st.download_button(
                    "Download PDF",
                    st.session_state.fig_bytes["fig_e_heatmap"]["pdf"],
                    "fig_e_heatmap.pdf", "application/pdf",
                    key="dl_fig_e_pdf",
                )
            with col_svg:
                st.download_button(
                    "Download SVG",
                    st.session_state.fig_bytes["fig_e_heatmap"]["svg"],
                    "fig_e_heatmap.svg", "image/svg+xml",
                    key="dl_fig_e_svg",
                )
            plt.close(fig_e)
        else:
            st.warning("No heatmap data available with current parameters.")

    # ======================================================================
    # TAB 6: Export
    # ======================================================================
    with tab6:
        st.header("Export All Results")

        st.subheader("Figures")
        st.markdown("Download all figures as a ZIP archive (PDF + SVG).")

        # Build ZIP from cached figure bytes (no regeneration needed)
        cached_bytes = st.session_state.get("fig_bytes", {})
        if cached_bytes:
            import io, zipfile
            buf = io.BytesIO()
            with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
                for name, fmt_dict in cached_bytes.items():
                    zf.writestr(f"{name}.pdf", fmt_dict["pdf"])
                    zf.writestr(f"{name}.svg", fmt_dict["svg"])
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

else:
    st.info("Click **Run Analysis** in the sidebar to start.")
