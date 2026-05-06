# bacTRAP-to-HypoMap Mapping Tool

A Streamlit application for mapping bacTRAP (Translating Ribosome Affinity Purification) bulk RNA-seq data onto the murine [HypoMap](https://doi.org/10.1038/s42255-022-00657-y) single-cell atlas (Steuernagel et al., *Nature Metabolism* 2022). The primary mapping method is AUCell (Aibar et al., *Nat. Methods* 2017) — rank-based, threshold-free, normalisation-insensitive — cross-checked against preranked GSEA (Subramanian et al., *PNAS* 2005) on cluster marker genes. The two methods are combined into a single per-cluster SCORE (mean of AUCell and GSEA percentile ranks). Earlier versions of this tool also reported Spearman correlation, a Fisher's-exact marker-overlap test, and NNLS deconvolution; on the data we tested those produced unreliable rankings and were dropped. The output is a three-panel main figure (AUCell UMAP, cell-type UMAP with top-15 clusters highlighted, per-cluster violin distributions) plus the composite-score panel and supplementary AUCell / GSEA diagnostics, exported as publication-ready PDF / SVG.

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
| Min IP expression | Minimum mean IP count (DESeq2 baseMean) — suppresses log₂FC inflation for genes with near-zero Input | 10 |
| Gene ranking metric | How top-N genes are chosen: π-score (default), log₂FoldChange, or -log₁₀(padj) | π-score |
| Top N genes for scoring | Size of the bacTRAP signature used for AUCell (and for all other scoring methods) | 50 |
| AUCell top-ranked fraction | Fraction of genes per cell considered for the recovery-curve AUC (AUCell τ; Aibar 2017 default 5 %) | 0.05 |
| Min cells for AUCell top-N ranking | Cluster-size floor for the top-N figures (main 1b/1c + supplementary S2). Excludes small clusters whose mean is dominated by shrinkage variance | 20 |
| Markers per cluster | Number of marker genes per cluster used as gene sets in preranked GSEA | 100 |
| Min cells per cluster | Minimum cells required to include a cluster in marker analysis | 10 |
| UMAP subsample | Number of cells to subsample for UMAP rendering (scoring uses all cells) | 50,000 |
| Cre-driver gene | Gene used for the Cre-driver Check sanity panel and the baseline filter below | `Pnoc` |
| Expression threshold (fraction) | Display-only: min fraction of cells expressing the Cre-driver gene to flag a cluster as "expressing" in the sanity table; does **not** change any ranking | 0.05 |
| Baseline Cre-driver mean expression (log-norm) | **Ranking filter.** When > 0, drop clusters whose mean (log-normalised) Cre-driver expression falls below this floor before GSEA, the AUCell cluster rankings (figs 1b / S2 / 1c) and the composite (AUCell + GSEA) score; survivors are re-ranked against one another. Per-cell AUCell scores (Fig 1a) are unaffected. Set to 0 to disable | 0.00 |
| Figure width | Single column (89 mm) or double column (183 mm) | Single |

### Tabs

| Tab | Contents |
|---|---|
| **Data Overview** | Gene/cell/cluster counts, match rate, enriched gene list, bacTRAP volcano plot (Pnoc highlighted), gene matching diagnostics |
| **AUCell (Main Figure)** | Main **Figure 1a–c** (AUCell UMAP, cell-type UMAP with top-15 clusters highlighted, violin distributions), the supplementary AUCell barplot (S1) and global histogram (S2), and the composite **Figure 1e** SCORE (AUCell + GSEA percentile mean) |
| **Cre-driver Check** | Sanity check + optional baseline filter. Per-cluster expression of the Cre-driver gene (default `Pnoc`) across the top-ranked clusters from the composite SCORE, with a configurable "fraction expressing" threshold that highlights (but does not drop) suspicious hits. The sidebar's "Baseline Cre-driver mean expression" slider additionally enforces a minimum expression floor that feeds back into the cluster-level rankings (GSEA, AUCell cluster mean, composite SCORE); filtered CSV exports are suffixed with the filter signature (e.g. `composite_ranking_pnoc_ge0p05.csv`) |
| **UMAP Projection (Suppl.)** | Two-panel UMAP: cell-type annotation + AUCell score, side-by-side (Supplementary Figure S3) |
| **GSEA (Suppl.)** | GSEA results table + enrichment curves (Suppl. S4) + NES barplot (Suppl. S5) |
| **Export** | Download all figures (ZIP of PDF+SVG), all result tables (CSV), and diagnostic log file |

---

## Analysis Pipeline

### 1. Data Loading & Gene Matching

- Loads both files with Streamlit caching (`@st.cache_resource` for the large h5ad).
- Builds a comprehensive bidirectional gene lookup from the HypoMap raw layer (51,216 genes), indexing both resolved gene symbols and Ensembl IDs.
- Auto-detects the best bacTRAP gene column by testing all candidates against the lookup. Achieves ~95% match rate with typical DESeq2 output.
- Deduplicates matched genes when multiple Ensembl IDs resolve to the same symbol.
- Reports overlap statistics, match diagnostics (sample gene names from both datasets), and enriched gene counts.

### 2. Cluster Marker Genes (input to GSEA)

- Computes marker genes per cluster via `scanpy.tl.rank_genes_groups` (Wilcoxon test), or loads pre-computed markers from `.uns['rank_genes_groups']` if available and matching the selected annotation column.
- Marker gene names are converted from Ensembl IDs to gene symbols when raw data uses Ensembl notation.
- Output: a cluster → top-N gene-symbol dictionary fed to preranked GSEA below.

### 3. UMAP Enrichment Projection

- Projects the per-cell AUCell score (see §5) onto the HypoMap UMAP alongside the cell-type annotation for side-by-side comparison.
- Cell-type legend placed below the annotation panel to avoid overlap.
- Subsamples cells for rendering performance (configurable, default 50k).

### 4. GSEA Preranked Enrichment

- Ranks all matched genes by bacTRAP log2FC (no hard cutoff).
- For each cluster's marker gene set, computes a running enrichment score (Subramanian et al., PNAS 2005).
- Significance assessed via 1,000 permutations (gene labels and enrichment values permuted in sync), with FDR correction.
- Normalized Enrichment Score (NES) enables cross-cluster comparison.
- Output: clusters ranked by NES, with ES, p-values, and FDR-adjusted p-values.

### 5. AUCell Scoring

Per-cell enrichment of the bacTRAP signature, computed with a faithful Python implementation of the AUCell method (Aibar et al., *Nature Methods* 2017; R reference: `AUCell::AUCell_calcAUC`). See METHODS.md for the formal description; in summary:

- **Input-layer QC** — a 500-cell sample of the chosen layer (`adata.raw.X` by default) is checked to verify integer counts with max > 50; mismatches raise a user-visible warning (AUCell's tie structure is formulated on raw counts).
- **Signature matching** — uses the shared `_build_adata_gene_lookup` so symbols *and* Ensembl IDs map correctly regardless of namespace. Unmatched signature genes are logged and surfaced in the tab.
- **Per-cell ranking with random tie-breaking** — a per-cell uniform jitter strictly smaller than the smallest distinct expression gap (0.49 for integer counts) is added before partition/sort. This preserves the order of distinct values while randomising ties, equivalent to R's `ties.method = "random"`. Without this, `numpy.argpartition`/`argsort` would break ties by matrix-column position and bias scores toward signature genes that sit at low gene indices.
- **AUC of the recovery curve** — for each cell, the discrete area under the recovery curve is computed within the top *k* ranked genes (*k* = max(⌈τ·*G*⌉, *n*<sub>query</sub>), τ user-configurable, default 5 %) and normalised by the theoretical maximum, giving scores on [0, 1]. When the signature is wider than the τ window, *k* is raised so every signature gene can contribute; the bump is logged and surfaced as a warning.
- **Reproducible and vectorised** — cells are processed in 5 000-cell chunks with a seeded RNG (default `seed=0`). `argpartition` finds the top-*k* gene indices in linear time per cell; `np.cumsum` on the hit indicator evaluates the AUC.
- **Cluster-level summaries and significance** — the per-cluster table (`aucell_per_cluster.csv`) contains mean, median, SD, SEM, and a one-sided Welch's *t*-test of "cluster > rest-of-atlas" with Benjamini–Hochberg *q*-values. Clusters with fewer than 20 cells are excluded from all top-N figure rankings (configurable via the "Min cells for AUCell top-N ranking" slider) because very small clusters' means are dominated by shrinkage variance; the raw CSV is unaffected and still contains every cluster.
- **Output** — per-cell AUCell scores projected onto the HypoMap UMAP (Figure 1a), the same UMAP with the 15 highest-mean eligible clusters highlighted by cell type (Figure 1b), and violin distributions for those same 15 clusters (Figure 1c). Supplementary panels S1 (barplot), S2 (histogram) and S3 (paired UMAP) expose additional views.

### 6. Composite SCORE (AUCell + GSEA)

- For each method (AUCell cluster mean, GSEA NES) cluster scores are converted to percentile ranks (0 → 1).
- The composite SCORE is the arithmetic mean of the two percentiles (`nanmean`, so a cluster scored by only one method is not penalised with a zero).
- Clusters are ranked by descending SCORE and rendered as a heatmap (**Figure 1e**) with composite values annotated to the right.
- When the optional baseline Cre-driver filter is active, percentile conversion runs on the restricted cluster universe so the displayed scores compare survivors against one another rather than against the full atlas.

### 7. Cre-driver Expression Check and Baseline Filter

Two related knobs that share the Cre-driver gene set via the sidebar:

- **Diagnostic (always active).** Reports per-cluster log-normalised mean expression and fraction expressing for the Cre-driver gene itself (default `Pnoc`; configurable). Clusters in the composite-SCORE top-N that fall below the fraction threshold are highlighted in the sanity table so the user can spot hits that may reflect Cre lineage tracing or snRNA-seq dropout rather than ongoing transcription. This highlighting does *not* change any ranking.
- **Baseline filter (optional).** The **Baseline Cre-driver mean expression** slider applies a minimum mean-expression floor. When set above 0, clusters below the floor are dropped *before* percentile conversion in GSEA, the AUCell cluster mean, and the composite SCORE — so the surviving clusters re-rank against each other rather than inheriting their global positions. The AUCell cluster figures (1b / S1 / 1c) and `aucell_per_cluster.csv` are filtered the same way. Per-cell AUCell scores, the per-cell UMAP (Fig 1a), and `aucell_per_cell.csv` are intentionally left unfiltered because they carry no cluster identity.
- **Edge cases:** if the Cre-driver gene is absent from the atlas, or the threshold discards every cluster, the filter is disabled and a persistent warning/error banner is rendered above the tab group; all tabs revert to unfiltered behaviour and CSV filenames are the default.
- **Filenames of filtered exports** are suffixed with the filter signature (e.g. `gsea_results_pnoc_ge0p05.csv`) so a collaborator opening a trimmed table from the Export tab can tell at a glance that it's a subset, not the full atlas.
- **Output:** diagnostic barplot of fraction-expressing per cluster with the user-defined threshold marked, plus a per-cluster sanity table and CSV; when the baseline filter is active, every ranking-related tab acts on the restricted cluster universe.

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

**Main Figure 1 — AUCell maps the bacTRAP signature onto the HypoMap atlas.** Five panels:

| Panel | Type | Description | Responds to baseline filter? |
|---|---|---|---|
| **1a** | UMAP | AUCell enrichment score projected onto the HypoMap UMAP (magma colormap, 2nd/98th-percentile clip, bottom-left axis arrows) | No (per-cell) |
| **1b** | UMAP | Same UMAP layout coloured by cell-type annotation, with the **top-15 AUCell-ranked clusters** (≥ 20 cells, passing the baseline filter when active) highlighted over a grey "Other" background | Yes |
| **1c** | Violin plots | Per-cluster AUCell distributions for the same top-15 clusters, ordered by cluster mean (top = highest) | Yes |
| **1e** | Heatmap | Composite SCORE (AUCell mean + GSEA NES, percentile mean) for the top 20 clusters, YlOrRd colormap, composite values annotated | Yes |

**Supplementary Figures** — AUCell diagnostics + GSEA. Cluster-level panels (S1, S4, S5) honour the optional baseline Cre-driver expression filter and recompute their rankings over the surviving clusters when the slider is > 0; per-cell panels (S2, S3, Fig 1a) do not carry cluster identity and are unaffected.

| Figure | Type | Description | Baseline filter |
|---|---|---|---|
| **Volcano** (Data Overview) | Scatter plot | bacTRAP gene-level volcano (log₂FC vs −log₁₀ padj); only the user-specified highlight gene (default Pnoc) is labelled | No (gene-level) |
| **S1** | Horizontal barplot | Mean AUCell score per cluster (top 25, ≥ 20 cells), SEM error bars, magma colormap | Yes |
| **S2** | Histogram | Global AUCell score distribution with 90/95/99th-percentile markers | No (per-cell) |
| **S3** | Two-panel UMAP | Cell-type annotation (left) and per-cell AUCell score (right), paired for side-by-side comparison | No (per-cell) |
| **S4** | Line plot | Running GSEA enrichment curves for top 5 clusters (legend right of plot) | Yes (curves drawn for filtered top-5) |
| **S5** | Horizontal barplot | GSEA normalised enrichment scores (NES), significance coloured | Yes |
| **Sanity** (Cre-driver tab) | Two-panel barplot | Per-cluster mean expression and fraction of cells expressing the Cre-driver gene (default Pnoc) for the top-ranked clusters; diagnostic only — does not modify any ranking | Cluster ordering follows the active filter |

---

## Project Structure

```
bacTRAP-to-HypoMapMapping/
├── app.py              # Main Streamlit application (UI, 6 tabs, orchestration)
├── data_loading.py     # Data I/O, gene matching, single-gene cluster stats
├── analysis.py         # Analysis methods (markers, GSEA, AUCell, composite ranking)
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
- `compute_cluster_mean_expression()` -- mean expression per cluster (used by `compute_single_gene_cluster_stats`)
- `compute_fraction_expressing()` -- fraction of cells expressing per cluster (used by `compute_single_gene_cluster_stats`)
- `compute_single_gene_cluster_stats()` -- per-cluster mean expression + fraction expressing for one gene (Cre-driver sanity check)

**`analysis.py`**
- `get_enriched_genes()` -- filter by padj, log2FC, and optional minimum IP expression
- `rank_enriched_genes()` -- rank enriched genes by π-score (default), log₂FC, or -log₁₀(padj)
- `compute_marker_genes()` -- Wilcoxon-based marker detection via scanpy (Ensembl-to-symbol conversion); fed into GSEA
- `load_precomputed_markers()` -- attempts to read `.uns['rank_genes_groups']`
- `compute_gsea_enrichment()` -- preranked GSEA with permutation-based p-values
- `validate_aucell_input()` -- sanity-checks the AUCell input layer (raw counts, integer-ness, max value) and returns a machine-readable QC report
- `compute_aucell_scores()` -- rank-based AUCell scoring with per-cell random-jitter tie-breaking (Aibar 2017); optional `info_out` dict returns match rate, `n_top` and top-fraction-bump diagnostics
- `compute_cluster_enrichment_stats()` -- per-cluster Welch's one-sided *t*-test vs. rest of atlas, with Benjamini–Hochberg *q*-values
- `compute_composite_ranking()` -- percentile mean of AUCell cluster mean and GSEA NES (the two-method SCORE)

**`figures.py`**
- `setup_nature_style()` -- global matplotlib configuration for Nature specs
- `_add_umap_axis_arrows()` -- shared helper drawing the compact "UMAP1 / UMAP2" arrows in the bottom-left corner of any dimensionality-reduction panel
- `figure_bactrap_volcano()` -- bacTRAP gene-level volcano; only the user-supplied highlight gene is labelled
- `figure_umap_enrichment()` -- two-panel UMAP, cell types + enrichment (Suppl. S3)
- `figure_aucell_umap()` -- AUCell score projected onto HypoMap UMAP (**Figure 1a**)
- `figure_celltype_umap()` -- cell-type annotation UMAP with caller-supplied top-N clusters highlighted over a grey "Other" background (**Figure 1b**)
- `figure_aucell_cluster_barplot()` -- per-cluster AUCell means with SEM, size-filtered ranking (Suppl. S1)
- `figure_aucell_violins()` -- per-cluster AUCell distributions, size-filtered ranking (**Figure 1c**)
- `figure_aucell_histogram()` -- global AUCell score distribution (Suppl. S2)
- `figure_gsea_curves()` -- running enrichment score curves (Suppl. S4)
- `figure_gsea_barplot()` -- NES barplot (Suppl. S5)
- `figure_composite_ranking()` -- AUCell + GSEA percentile heatmap (**Figure 1e**)
- `figure_marker_gene_diagnostic()` -- Cre-driver gene per-cluster expression / fraction (sanity tab)
- `fig_to_bytes()` -- convert figure to PDF/SVG bytes

---

## Logging

The application writes a detailed log to `bactrap_hypomap.log` in the app directory. The log captures:

- User inputs: every sidebar parameter (cutoffs, thresholds, ranking metric, Cre-driver settings, figure width, etc.) is logged once per analysis run for reproducibility
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
- UMAP visualization subsamples cells (default 50k) while AUCell scoring and GSEA use all matched genes / all cells.
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

The tool combines two orthogonal cluster-mapping methods — AUCell (per-cell rank-sum on the bacTRAP signature) and preranked GSEA (preranked test on cluster markers) — into a single composite SCORE, with a Cre-driver expression sanity check that flags hits attributable to lineage tracing or snRNA-seq dropout rather than ongoing transcription.

---

## References

- Steuernagel, L. et al. HypoMap -- a unified single-cell gene expression atlas of the murine hypothalamus. *Nature Metabolism* **4**, 1402--1419 (2022). [DOI: 10.1038/s42255-022-00657-y](https://doi.org/10.1038/s42255-022-00657-y)
- Heiman, M. et al. A translational profiling approach for the molecular characterization of CNS cell types. *Cell* **135**, 738--748 (2008).
- Subramanian, A. et al. Gene set enrichment analysis: a knowledge-based approach. *PNAS* **102**, 15545--15550 (2005).
- Aibar, S. et al. SCENIC: single-cell regulatory network inference and clustering. *Nature Methods* **14**, 1083--1086 (2017).
- Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* **15**, 550 (2014).
- Wolf, F. A., Angerer, P. & Theis, F. J. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology* **19**, 15 (2018).
- Xiao, Y. et al. A novel significance score for gene selection and ranking. *Bioinformatics* **30**, 801--807 (2014).
