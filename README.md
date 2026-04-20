# bacTRAP-to-HypoMap Mapping Tool

A Streamlit application for mapping bacTRAP (Translating Ribosome Affinity Purification) bulk RNA-seq data onto the murine [HypoMap](https://doi.org/10.1038/s42255-022-00657-y) single-cell atlas (Steuernagel et al., *Nature Metabolism* 2022). Identifies which hypothalamic cell types best match the translational profile captured by a bacTRAP pulldown using five complementary methods, and produces publication-ready, Nature-grade figures exportable as PDF/SVG.

---

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Launch the app
streamlit run app.py
```

Then open the URL shown in your terminal (default: `http://localhost:8501`).

---

## Input Files

You provide two files via the app's sidebar:

### 1. bacTRAP FPKM Table (`.xlsx`)

An Excel file containing DESeq2 differential expression results from a bacTRAP experiment. Required columns:

| Column | Description |
|---|---|
| `log2FoldChange` | DESeq2 log2 fold change (IP vs Input) |
| `padj` | Benjamini-Hochberg adjusted p-value |

The gene identifier column is **auto-detected**. The app tests all candidate columns (including `gene_id`, `gene_name`, and the DataFrame index) against the HypoMap gene lookup and selects the column yielding the highest match rate. Both Ensembl IDs and gene symbols are supported. You can override the auto-selection via the "bacTRAP gene column" dropdown in the sidebar.

Additional columns that will be used if present:

| Column | Description |
|---|---|
| `gene_id` | Ensembl gene IDs (e.g. `ENSMUSG00000074604`) |
| `gene_name` | Gene symbols (e.g. `Mgst2`) |
| `IP` | Mean IP expression across replicates |
| `Input` | Mean Input expression across replicates |
| `gene_biotype` | Gene biotype annotation |
| `gene_description` | Gene description |

### 2. HypoMap Atlas (`.h5ad`)

The HypoMap AnnData object (~384,925 cells). Download from [CellxGene](https://cellxgene.cziscience.com/). It must contain:

- Expression matrix in `.X` (or `.raw.X`)
- Cell-type annotations at one or more hierarchical levels in `.obs`
- UMAP coordinates in `.obsm['X_umap']`

The app inspects `.obs` columns on load and lets you select which annotation level to use (e.g. `C7_named`, `C25_named`, `C66_named`, `C185_named`, `C286_named`, `C465_named`).

**Gene name handling:** When atlas `var_names` are Ensembl IDs, the app resolves gene symbols from the `feature_name` column in `var` metadata (or `raw.var`). If no symbol column is found in `raw.var`, it falls back to building an Ensembl-to-symbol map from `adata.var`. Both gene symbols and Ensembl IDs are indexed in the lookup, so matching works regardless of the bacTRAP identifier format.

---

## Application Layout

### Sidebar Controls

| Control | Description | Default |
|---|---|---|
| bacTRAP file path | Path to the `.xlsx` file | -- |
| HypoMap file path | Path to the `.h5ad` file | -- |
| Annotation column | Cell-type label column from `.obs` | Auto-detected (prefers `C185_named`) |
| bacTRAP gene column | Column containing gene identifiers | Auto-detected (highest match rate) |
| padj cutoff | Adjusted p-value threshold for enriched genes | 0.05 |
| log2FC cutoff | Minimum log2 fold change for enriched genes | 1.0 |
| Top N genes | Number of top enriched genes for scoring | 50 |
| Markers per cluster | Number of marker genes per cluster for overlap test | 100 |
| Min cells per cluster | Minimum cells required to include a cluster | 10 |
| UMAP subsample | Number of cells to subsample for UMAP rendering | 50,000 |
| Cre-driver gene | Gene used for the Cre-driver Check sanity panel | `Pnoc` |
| Expression threshold (fraction) | Min fraction of cells expressing the Cre-driver gene to count a cluster as "expressing" | 0.05 |
| Figure width | Single column (89 mm) or double column (183 mm) | Single |

### Tabs

| Tab | Contents |
|---|---|
| **Data Overview** | Gene/cell/cluster counts, match rate, enriched gene list, bacTRAP volcano plot, gene matching diagnostics |
| **AUCell (Main Figure)** | Main Figure 1a-1e: AUCell UMAP, per-cluster barplot, violin distributions, score histogram, and composite consensus ranking |
| **Cre-driver Check** | Sanity check: per-cluster expression of the Cre-driver gene (default `Pnoc`) across the top-ranked clusters from the composite consensus, with a configurable "fraction expressing" threshold to flag mapping hits that may reflect Cre lineage tracing rather than current expression |
| **Correlation (Suppl.)** | Ranked cluster table + correlation barplot (Supplementary Figure S1) |
| **UMAP Projection (Suppl.)** | Two-panel UMAP: cell types + enrichment score (Supplementary Figure S2) |
| **Marker Overlap (Suppl.)** | Fisher's test table + Fisher volcano plot + dot plot (Supplementary Figures S3, S4) |
| **Heatmap (Suppl.)** | Z-scored heatmap of top genes across top clusters (Supplementary Figure S5) |
| **NNLS (Suppl.)** | NNLS weights table + weight barplot (Supplementary Figure S6) |
| **GSEA (Suppl.)** | GSEA results table + enrichment curves + NES barplot (Supplementary Figures S7, S8) |
| **Export** | Download all figures (ZIP of PDF+SVG), all result tables (CSV), and diagnostic log file |

---

## Analysis Pipeline

### 1. Data Loading & Gene Matching

- Loads both files with Streamlit caching (`@st.cache_resource` for the large h5ad).
- Builds a comprehensive bidirectional gene lookup from the HypoMap raw layer (51,216 genes), indexing both resolved gene symbols and Ensembl IDs.
- Auto-detects the best bacTRAP gene column by testing all candidates against the lookup. Achieves ~95% match rate with typical DESeq2 output.
- Deduplicates matched genes when multiple Ensembl IDs resolve to the same symbol.
- Reports overlap statistics, match diagnostics (sample gene names from both datasets), and enriched gene counts.

### 2. Enrichment Correlation

- Computes mean expression per cluster in HypoMap for all overlapping genes.
- Filters NaN values from the enrichment vector before computing correlations.
- Correlates each cluster's expression profile against the bacTRAP `log2FoldChange` vector using both Pearson and Spearman correlation.
- Output: clusters ranked by Spearman rho, with r-values and p-values.

### 3. Marker Gene Overlap (Fisher's Exact Test)

- Defines "enriched genes" as those passing both the padj and log2FC cutoffs.
- Computes marker genes per cluster via `scanpy.tl.rank_genes_groups` (Wilcoxon test), or loads pre-computed markers from `.uns['rank_genes_groups']` if available and matching the selected annotation column.
- Marker gene names are converted from Ensembl IDs to gene symbols when raw data uses Ensembl notation.
- Tests overlap using a one-sided Fisher's exact test per cluster, with Benjamini-Hochberg FDR correction.
- Output: clusters ranked by enrichment p-value, with odds ratios, overlap counts, and overlapping gene names.

### 4. UMAP Enrichment Projection

- Computes a per-cell "bacTRAP enrichment score" as the z-scored mean expression of the top N enriched genes.
- Projects this score onto the HypoMap UMAP alongside the cell-type annotation for side-by-side comparison.
- Cell-type legend placed below the annotation panel to avoid overlap.
- Subsamples cells for rendering performance (configurable, default 50k).

### 5. Gene Expression Heatmap

- Z-scores mean expression of the top enriched genes across the highest-correlating clusters.
- Rows (genes) are hierarchically clustered (Ward's linkage); columns (clusters) are kept in correlation-ranked order.

### 6. NNLS Deconvolution

- Non-Negative Least Squares: finds non-negative cluster weights **w** that minimize `||A @ w - b||` where **A** is the cluster mean expression matrix and **b** is the bacTRAP log2FC vector.
- Produces quantitative contribution weights per cell type, answering not just *which* clusters match but *how much* each contributes to the enrichment profile.
- Output: clusters ranked by normalized weight.

### 7. GSEA Preranked Enrichment

- Ranks all matched genes by bacTRAP log2FC (no hard cutoff).
- For each cluster's marker gene set, computes a running enrichment score (Subramanian et al., PNAS 2005).
- Significance assessed via 1,000 permutations (gene labels and enrichment values permuted in sync), with FDR correction.
- Normalized Enrichment Score (NES) enables cross-cluster comparison.
- More powerful than Fisher's binary test because it uses the full ranking.
- Output: clusters ranked by NES, with ES, p-values, and FDR-adjusted p-values.

### 8. AUCell Scoring

- Computes per-cell enrichment scores using the Area Under the recovery Curve method (Aibar et al., *Nature Methods* 2017).
- For each cell: ranks all genes by expression, then measures how quickly the bacTRAP-enriched gene set is recovered in the top-ranked genes.
- Advantages over simple z-scored mean: rank-based (normalization-insensitive), focuses on highly expressed genes, threshold-free.
- Vectorized implementation processes cells in chunks with `np.cumsum`-based AUC for performance on large atlases.
- Output: per-cell AUCell scores projected onto the HypoMap UMAP.

---

## Figures

All figures follow Nature journal specifications:

| Property | Value |
|---|---|
| Font | Arial / Helvetica (sans-serif) |
| Axis labels | 7 pt |
| Tick labels | 6 pt |
| Panel titles | 8 pt bold |
| Axis line width | 0.5 pt |
| Plot element line width | 0.75 pt |
| Spines | Bottom and left only |
| Background | White, no gridlines |
| Rasterized elements | 300 DPI |
| Color palettes | Colorblind-friendly (`viridis`, `magma`, `tab10`/`tab20`) |
| Text in exports | Editable (TrueType in PDF, native text in SVG) |
| Dimensions | Single column: 89 mm (3.5 in), Double column: 183 mm (7.2 in) |

### Figure Descriptions

**Main Figure 1 — AUCell cell-type mapping** (rank-based, normalization-insensitive, threshold-free; the primary mapping method):

| Panel | Type | Description |
|---|---|---|
| **1a** | UMAP | AUCell enrichment scores projected onto HypoMap UMAP (magma colormap) |
| **1b** | Horizontal barplot | Mean AUCell score per cluster with SEM error bars (top 25) |
| **1c** | Violin plots | Per-cluster AUCell score distributions for top 15 clusters |
| **1d** | Histogram | Global AUCell score distribution with 90/95/99th percentile markers |
| **1e** | Heatmap | Composite consensus ranking across all methods (percentile scores, YlOrRd) |

**Supplementary Figures** — complementary analysis approaches:

| Figure | Type | Description |
|---|---|---|
| **Volcano** | Scatter plot | bacTRAP gene-level volcano (log2FC vs -log10 padj), with Pnoc and top enriched genes labeled (Data Overview tab) |
| **S1** | Horizontal barplot | Top 20 clusters by Spearman correlation, colored by rho |
| **S2** | Two-panel UMAP | Left: cell-type annotation (legend below), Right: bacTRAP enrichment score (magma) |
| **S3** | Volcano plot | log2(odds ratio) vs -log10(p-value) from Fisher's test, top hits labeled |
| **S4** | Dot plot | Top enriched genes vs correlation-ranked clusters (size = % expressing, color = mean expression) |
| **S5** | Heatmap | Z-scored expression, genes clustered by Ward's linkage, diverging RdBu_r colormap |
| **S6** | Horizontal barplot | NNLS deconvolution weights per cluster (magma colormap) |
| **S7** | Line plot | Running GSEA enrichment score curves for top 5 clusters (legend right of plot) |
| **S8** | Horizontal barplot | Normalized Enrichment Scores with FDR significance coloring |

---

## Project Structure

```
bacTRAP-to-HypoMapMapping/
├── app.py              # Main Streamlit application (UI, 9 tabs, orchestration)
├── data_loading.py     # Data I/O, gene matching, cluster expression computation
├── analysis.py         # All analysis methods (correlation, Fisher, NNLS, GSEA, AUCell)
├── figures.py          # Nature-grade figure generation and export
├── METHODS.md          # Publication-ready Materials & Methods and figure legends
├── requirements.txt    # Python dependencies
└── README.md           # This file
```

### Module Responsibilities

**`data_loading.py`**
- `load_hypomap()` -- cached h5ad loading with sparse matrix enforcement; logs structure details
- `load_bactrap()` -- cached Excel loading; logs columns, dtypes, sample values
- `get_annotation_columns()` -- discovers categorical/string `.obs` columns
- `_build_adata_gene_lookup()` -- comprehensive bidirectional lookup (symbols + Ensembl IDs)
- `_auto_select_gene_col()` -- tests all bacTRAP columns, picks best match rate
- `match_genes()` -- gene matching with auto-detection and deduplication
- `_map_var_indices_to_raw()` -- translates adata.var indices to adata.raw.var space
- `_extract_gene_submatrix()` -- memory-efficient column extraction with chunk processing
- `compute_cluster_mean_expression()` -- mean expression per cluster (deduplicates index)
- `compute_fraction_expressing()` -- fraction of cells expressing per cluster (deduplicates index)

**`analysis.py`**
- `compute_enrichment_correlation()` -- Pearson + Spearman correlation per cluster (NaN-safe, dedup-safe)
- `get_enriched_genes()` -- filter by padj and log2FC thresholds (with column validation)
- `compute_marker_genes()` -- Wilcoxon-based marker detection via scanpy (Ensembl-to-symbol conversion)
- `load_precomputed_markers()` -- attempts to read `.uns['rank_genes_groups']`
- `fisher_overlap_test()` -- one-sided Fisher's exact test with FDR correction
- `compute_enrichment_score()` -- per-cell z-scored mean expression
- `compute_zscore_heatmap_data()` -- z-scored expression matrix for heatmaps
- `compute_nnls_deconvolution()` -- non-negative least squares decomposition
- `compute_gsea_enrichment()` -- preranked GSEA with permutation-based p-values
- `compute_aucell_scores()` -- rank-based AUCell scoring (vectorized)
- `compute_composite_ranking()` -- percentile-averaged consensus across all methods

**`figures.py`**
- `setup_nature_style()` -- global matplotlib configuration for Nature specs
- `figure_bactrap_volcano()` -- bacTRAP gene-level volcano with highlight support
- `figure_correlation_barplot()` -- correlation barplot
- `figure_umap_enrichment()` -- two-panel UMAP (cell types + enrichment)
- `figure_dotplot()` -- enriched gene dot plot
- `figure_volcano_enrichment()` -- Fisher's test volcano
- `figure_heatmap()` -- z-scored expression heatmap
- `figure_nnls_barplot()` -- NNLS weight barplot
- `figure_gsea_curves()` -- running enrichment score curves
- `figure_gsea_barplot()` -- NES barplot
- `figure_aucell_umap()` -- AUCell UMAP projection
- `fig_to_bytes()` -- convert figure to PDF/SVG bytes

---

## Logging

The application writes a detailed log to `bactrap_hypomap.log` in the app directory. The log captures:

- HypoMap structure: var_names, raw.var_names, all column names, sample values
- bacTRAP structure: all columns, dtypes, sample values from every string column
- Gene name resolution: Ensembl detection, symbol column search, fallback attempts
- Lookup construction: final size, sample entries
- Auto-detection: per-column match counts for gene column selection
- Match results: matched/unmatched gene samples, deduplication counts, final match rate

The log file can be downloaded from the **Export** tab (standalone or included in the figures ZIP).

---

## Memory & Performance Notes

- The HypoMap atlas is ~3.9 GB. The app enforces sparse matrix representation and caches the loaded object so it is not reloaded on interaction.
- Expression submatrices are extracted in column chunks (200 genes per chunk) to avoid materializing the full sparse matrix as dense.
- When `adata.raw` has more genes than `adata.var` (common after gene filtering), indices are remapped via gene name lookup rather than assuming positional correspondence.
- UMAP visualization subsamples cells (default 50k) while statistical analyses (correlation, Fisher's, NNLS, GSEA) use all cells.
- AUCell scoring is vectorized with `np.cumsum`-based AUC computation and processes cells in 5,000-cell chunks.
- GSEA uses 1,000 permutations per cluster; this is the most time-consuming step for atlases with many clusters.
- Scatter points in UMAP figures use `rasterized=True` to keep exported PDF/SVG file sizes manageable.
- Figure bytes are cached in Streamlit session state after first render to avoid regeneration.
- A progress bar above the tabs tracks all long-running steps.

---

## Dependencies

- Python 3.10+
- streamlit >= 1.28
- scanpy >= 1.9.6
- anndata >= 0.10
- pandas >= 2.0
- numpy >= 1.24
- scipy >= 1.11
- matplotlib >= 3.7
- seaborn >= 0.12
- openpyxl >= 3.1
- adjustText >= 0.8
- scikit-learn >= 1.3
- statsmodels >= 0.14
- h5py >= 3.9

---

## Context

The bacTRAP IP sample is already enriched for a specific Cre-expressing neuronal population via ribosomal tagging in the preoptic area (PoA). This is **not** a standard bulk deconvolution problem -- the goal is to identify which HypoMap cell types/subtypes best match the translational profile captured by the bacTRAP pulldown. All labels and metrics in the app reflect this framing (e.g., "bacTRAP enrichment score", not "cell-type proportion").

The tool applies five complementary approaches (correlation, Fisher's overlap, NNLS deconvolution, GSEA, and AUCell) to provide a robust, multi-method consensus on which cell populations are captured by the bacTRAP pulldown.

---

## References

- Steuernagel, L. et al. HypoMap -- a unified single-cell gene expression atlas of the murine hypothalamus. *Nature Metabolism* **4**, 1402--1419 (2022). [DOI: 10.1038/s42255-022-00657-y](https://doi.org/10.1038/s42255-022-00657-y)
- Heiman, M. et al. A translational profiling approach for the molecular characterization of CNS cell types. *Cell* **135**, 738--748 (2008).
- Subramanian, A. et al. Gene set enrichment analysis: a knowledge-based approach. *PNAS* **102**, 15545--15550 (2005).
- Aibar, S. et al. SCENIC: single-cell regulatory network inference and clustering. *Nature Methods* **14**, 1083--1086 (2017).
- Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* **15**, 550 (2014).
- Wolf, F. A., Angerer, P. & Theis, F. J. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology* **19**, 15 (2018).
