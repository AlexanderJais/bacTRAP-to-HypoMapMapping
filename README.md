# bacTRAP-to-HypoMap Mapping Tool

A Streamlit application for mapping bacTRAP (Translating Ribosome Affinity Purification) bulk RNA-seq data onto the murine [HypoMap](https://doi.org/10.1038/s42255-022-00657-y) single-cell atlas (Steuernagel et al., *Nature Metabolism* 2022). Identifies which hypothalamic cell types best match the translational profile captured by a bacTRAP pulldown, and produces publication-ready, Nature-grade figures exportable as PDF/SVG.

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

An Excel file containing DESeq2 differential expression results from a bacTRAP experiment. Expected columns:

| Column | Description |
|---|---|
| `gene_id` | Ensembl gene IDs (e.g. `ENSMUSG00000074604`) |
| `gene_name` | Gene symbols (e.g. `Mgst2`) |
| `PoA_IP1_fpkm` ... `PoA_IP3_fpkm` | Per-replicate FPKM for IP samples |
| `PoA_Input1_fpkm` ... `PoA_Input3_fpkm` | Per-replicate FPKM for Input samples |
| `PoA_IP1_count` ... `PoA_Input3_count` | Per-replicate raw counts |
| `IP` | Mean IP expression across replicates |
| `Input` | Mean Input expression across replicates |
| `log2FoldChange` | DESeq2 log2 fold change (IP vs Input) |
| `pvalue` | DESeq2 raw p-value |
| `padj` | Benjamini-Hochberg adjusted p-value |
| `gene_biotype` | Gene biotype annotation |
| `gene_description` | Gene description |
| `gene_chr` | Chromosome |
| `gene_length` | Gene length in bp |

### 2. HypoMap Atlas (`.h5ad`)

The HypoMap AnnData object (~384,925 cells). Download from [CellxGene](https://cellxgene.cziscience.com/). It must contain:

- Expression matrix in `.X` (or `.raw.X`)
- Cell-type annotations at one or more hierarchical levels in `.obs`
- UMAP coordinates in `.obsm['X_umap']`

The app inspects `.obs` columns on load and lets you select which annotation level to use.

---

## Application Layout

### Sidebar Controls

| Control | Description | Default |
|---|---|---|
| bacTRAP file path | Path to the `.xlsx` file | — |
| HypoMap file path | Path to the `.h5ad` file | — |
| Annotation column | Cell-type label column from `.obs` | Auto-detected |
| padj cutoff | Adjusted p-value threshold for enriched genes | 0.05 |
| log2FC cutoff | Minimum log2 fold change for enriched genes | 1.0 |
| Top N genes | Number of top enriched genes for scoring | 50 |
| Markers per cluster | Number of marker genes per cluster for overlap test | 100 |
| Min cells per cluster | Minimum cells required to include a cluster | 10 |
| UMAP subsample | Number of cells to subsample for UMAP rendering | 50,000 |
| Figure width | Single column (89 mm) or double column (183 mm) | Single |

### Tabs

| Tab | Contents |
|---|---|
| **Data Overview** | Gene/cell/cluster counts, match rate, enriched gene list, data preview |
| **Correlation Analysis** | Ranked cluster table + Figure A (correlation barplot) |
| **UMAP Projection** | Figure B (two-panel UMAP: cell types + enrichment score) |
| **Marker Overlap** | Fisher's test table + Figure C (dot plot) + Figure D (volcano plot) |
| **Gene Heatmap** | Figure E (z-scored heatmap of top genes across top clusters) |
| **Export** | Download all figures (ZIP of PDF+SVG) and all result tables (CSV) |

---

## Analysis Pipeline

### 1. Data Loading & QC

- Loads both files with Streamlit caching (`@st.cache_resource` for the large h5ad).
- Matches `gene_name` from the bacTRAP table to HypoMap gene names using case-insensitive matching. If HypoMap `var_names` are Ensembl IDs, it falls back to `var['gene_name']` or similar columns.
- Reports overlap statistics and enriched gene counts.

### 2. Enrichment Correlation

- Computes mean expression per cluster in HypoMap for all overlapping genes.
- Correlates each cluster's expression profile against the bacTRAP `log2FoldChange` vector using both Pearson and Spearman correlation.
- Output: clusters ranked by Spearman rho, with r-values and p-values.

### 3. Marker Gene Overlap (Fisher's Exact Test)

- Defines "enriched genes" as those passing both the padj and log2FC cutoffs.
- Computes marker genes per cluster via `scanpy.tl.rank_genes_groups` (Wilcoxon test), or loads pre-computed markers from `.uns['rank_genes_groups']` if available.
- Tests overlap using a one-sided Fisher's exact test per cluster, with Benjamini-Hochberg FDR correction.
- Output: clusters ranked by enrichment p-value, with odds ratios, overlap counts, and overlapping gene names.

### 4. UMAP Enrichment Projection

- Computes a per-cell "bacTRAP enrichment score" as the z-scored mean expression of the top N enriched genes (using `scanpy.tl.score_genes`).
- Projects this score onto the HypoMap UMAP alongside the cell-type annotation for side-by-side comparison.
- Subsamples cells for rendering performance (configurable, default 50k).

### 5. Gene Expression Heatmap

- Z-scores mean expression of the top enriched genes across the highest-correlating clusters.
- Rows (genes) are hierarchically clustered; columns (clusters) are kept in correlation-ranked order.

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

| Figure | Type | Description |
|---|---|---|
| **A** | Horizontal barplot | Top 20 clusters by Spearman correlation, colored by rho |
| **B** | Two-panel UMAP | Left: cell-type annotation, Right: bacTRAP enrichment score (magma) |
| **C** | Dot plot | Top enriched genes vs top clusters (size = % expressing, color = mean expression) |
| **D** | Volcano plot | log2(odds ratio) vs -log10(p-value) from Fisher's test, top hits labeled |
| **E** | Heatmap | Z-scored expression, genes clustered by Ward's linkage, diverging RdBu_r colormap |

---

## Project Structure

```
bacTRAP-to-HypoMapMapping/
├── app.py              # Main Streamlit application (UI, tabs, orchestration)
├── data_loading.py     # Data I/O, gene matching, cluster expression computation
├── analysis.py         # Correlation, Fisher's test, enrichment scoring, marker genes
├── figures.py          # Nature-grade figure generation and export utilities
├── requirements.txt    # Python dependencies
└── README.md           # This file
```

### Module Responsibilities

**`data_loading.py`**
- `load_hypomap()` — cached h5ad loading with sparse matrix enforcement
- `load_bactrap()` — cached Excel loading
- `get_annotation_columns()` — discovers categorical/string `.obs` columns
- `match_genes()` — case-insensitive gene symbol matching between datasets
- `compute_cluster_mean_expression()` — mean expression per cluster for a gene set
- `compute_fraction_expressing()` — fraction of cells expressing per cluster
- `subsample_adata()` — random subsampling for visualization

**`analysis.py`**
- `compute_enrichment_correlation()` — Pearson + Spearman correlation per cluster
- `get_enriched_genes()` — filter by padj and log2FC thresholds
- `compute_marker_genes()` — Wilcoxon-based marker detection via scanpy
- `load_precomputed_markers()` — attempts to read `.uns['rank_genes_groups']`
- `fisher_overlap_test()` — one-sided Fisher's exact test with FDR correction
- `compute_enrichment_score()` — per-cell enrichment score via `sc.tl.score_genes`
- `compute_zscore_heatmap_data()` — z-scored expression matrix for heatmaps

**`figures.py`**
- `setup_nature_style()` — global matplotlib configuration for Nature specs
- `figure_correlation_barplot()` — Figure A
- `figure_umap_enrichment()` — Figure B
- `figure_dotplot()` — Figure C
- `figure_volcano_enrichment()` — Figure D
- `figure_heatmap()` — Figure E
- `fig_to_bytes()` — convert figure to PDF/SVG bytes
- `create_all_figures_zip()` — bundle all figures into a ZIP archive

---

## Memory & Performance Notes

- The HypoMap atlas is ~3.9 GB. The app enforces sparse matrix representation and caches the loaded object so it is not reloaded on interaction.
- Cluster mean expressions are computed directly from sparse matrices without full dense conversion.
- UMAP visualization subsamples cells (default 50k) while statistical analyses (correlation, Fisher's test) use all cells.
- Scatter points in UMAP figures use `rasterized=True` to keep exported PDF/SVG file sizes manageable.
- A progress bar tracks long-running steps.

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

The bacTRAP IP sample is already enriched for a specific Cre-expressing neuronal population via ribosomal tagging in the preoptic area (PoA). This is **not** a standard bulk deconvolution problem — the goal is to identify which HypoMap cell types/subtypes best match the translational profile captured by the bacTRAP pulldown. All labels and metrics in the app reflect this framing (e.g., "bacTRAP enrichment score", not "cell-type proportion").

---

## References

- Steuernagel, L. et al. HypoMap — a unified single-cell gene expression atlas of the murine hypothalamus. *Nature Metabolism* **4**, 1402–1419 (2022). [DOI: 10.1038/s42255-022-00657-y](https://doi.org/10.1038/s42255-022-00657-y)
- Heiman, M. et al. A translational profiling approach for the molecular characterization of CNS cell types. *Cell* **135**, 738–748 (2008).
