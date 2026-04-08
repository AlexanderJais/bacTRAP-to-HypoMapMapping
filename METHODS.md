# Materials and Methods: bacTRAP-to-HypoMap Transcriptomic Projection

*The text below is written in the style of a Nature journal Methods section. Adjust tense, length, and detail to match your target journal's requirements.*

---

## bacTRAP translational profiling

Translating ribosome affinity purification (bacTRAP) was performed on preoptic area (PoA) tissue from mice expressing a Cre-dependent EGFP-tagged ribosomal subunit (EGFP-L10a) in a cell-type-specific manner. Immunoprecipitated (IP) and total input mRNA were profiled by RNA-seq. Reads were aligned and quantified to obtain per-gene FPKM values across three biological replicates per condition. Differential expression between IP and Input was performed using DESeq2 (Love et al., 2014), yielding log2 fold changes and Benjamini-Hochberg adjusted p-values (padj) for each gene. Genes with padj < 0.05 and log2FC > 1 were classified as significantly enriched in the bacTRAP translational profile.

## Single-cell reference atlas

The HypoMap single-cell RNA-seq atlas of the murine hypothalamus (Steuernagel et al., 2022) was used as the reference dataset (~384,925 cells). The atlas was obtained as an AnnData h5ad file from the CellxGene repository. Cell-type annotations at the desired hierarchical level were selected from the `.obs` metadata. Expression data were accessed from the `.raw` slot when available (preserving all genes prior to feature selection) and maintained in sparse matrix format throughout to manage memory.

## Gene matching

Gene symbols from the bacTRAP dataset were matched to the HypoMap atlas using case-insensitive string matching. When the atlas used Ensembl gene identifiers as primary variable names, gene symbols were retrieved from the `gene_name` column of the variable metadata. When the atlas raw layer contained a different number of genes than the filtered variable set (common after highly variable gene selection), gene indices were remapped from filtered to raw space via name-based lookup to ensure correct expression values were retrieved.

## Enrichment correlation analysis

Mean expression was computed per cell-type cluster in HypoMap for all genes present in both datasets. Pearson and Spearman correlation coefficients were calculated between each cluster's mean expression profile and the bacTRAP log2 fold change vector. Genes with missing (NaN) enrichment values were excluded prior to correlation computation. A minimum of three overlapping genes was required. Clusters were ranked by Spearman rho.

## Marker gene overlap analysis

Cluster-specific marker genes were identified using the Wilcoxon rank-sum test via scanpy's `rank_genes_groups` function (Wolf et al., 2018), retaining the top 100 markers per cluster (user-configurable). When pre-computed marker genes were available in the atlas `.uns` slot and matched the selected annotation level (>50% cluster name overlap), these were used directly. Overlap between bacTRAP-enriched genes and each cluster's marker set was assessed using a one-sided Fisher's exact test. P-values were corrected for multiple testing using the Benjamini-Hochberg procedure (Benjamini and Hochberg, 1995). Clusters with fewer than 10 cells (user-configurable) were excluded.

## UMAP enrichment projection

A per-cell bacTRAP enrichment score was computed as the z-scored mean expression of the top N significantly enriched genes (default N = 50, ranked by log2FC among genes with padj < 0.05 and log2FC > 1). This score was projected onto the HypoMap UMAP embedding. For visualization, cells were randomly subsampled to 50,000 (user-configurable) while all cells were used for statistical analyses.

## Non-negative least squares deconvolution

To quantify the contribution of each cell type to the bacTRAP enrichment profile, non-negative least squares (NNLS) deconvolution was performed. The mean expression matrix **A** (genes x clusters) and the bacTRAP log2FC vector **b** were constructed from overlapping genes. The optimization problem min ||**A** **w** - **b**||₂ subject to **w** ≥ 0 was solved using `scipy.optimize.nnls` (Lawson and Hanson, 1995). The resulting weight vector **w** was normalized to sum to 1, providing fractional contribution estimates per cluster. Genes with missing enrichment values were excluded.

## Gene set enrichment analysis (GSEA)

A preranked GSEA approach (Subramanian et al., 2005) was applied to assess whether cluster-specific marker genes were enriched among the most bacTRAP-enriched genes. All matched genes were ranked by bacTRAP log2FC in descending order. For each cluster's marker gene set, a weighted running enrichment score was computed, where hit weights were proportional to the absolute enrichment value at each gene's rank. Significance was assessed by permutation testing (1,000 permutations), in which gene labels and their associated enrichment values were permuted in concert to preserve the rank-value structure while randomizing gene set membership. The enrichment score was normalized (NES) by dividing by the mean absolute enrichment score across permutations to enable comparison across gene sets of different sizes. P-values were corrected using the Benjamini-Hochberg procedure.

## AUCell scoring

Per-cell enrichment was additionally quantified using the AUCell method (Aibar et al., 2017). For each cell, all genes were ranked by expression level in descending order. The area under the recovery curve (AUC) was computed for the bacTRAP-enriched gene set within the top 5% of ranked genes (user-configurable). The AUC was normalized to the maximum possible value (number of query genes x number of top-ranked genes considered). This rank-based approach is robust to differences in normalization and sequencing depth, and does not require an expression threshold. AUCell scores were computed in chunks of 5,000 cells using vectorized cumulative sum operations.

## Composite consensus ranking

To integrate results across all methods, a composite ranking was generated. For each method (Spearman correlation, Fisher's exact test, NNLS, and GSEA NES), cluster scores were converted to percentile ranks scaled from 0 to 1. The composite score for each cluster was computed as the arithmetic mean of its percentile ranks across all methods that returned results for that cluster (methods with missing data for a cluster were excluded via `nanmean` rather than penalized with a zero score). Clusters were ranked by descending composite score.

## Figure generation

All figures were generated using matplotlib (v3.7+) following Nature journal specifications: Arial font (sans-serif), 7 pt axis labels, 6 pt tick labels, 8 pt bold panel titles, 0.5 pt axis line width, white background without gridlines, and removal of top and right spines. Colorblind-friendly palettes were used throughout (viridis, magma for sequential data; tab10/tab20 for categorical). UMAP scatter plots were rasterized at 300 DPI to manage file size while text elements remained as editable vector objects in PDF (TrueType embedding, fonttype 42) and SVG (native text, fonttype none) exports. Single-column figures were 89 mm (3.5 in) wide; double-column figures were 183 mm (7.2 in) wide.

## Software

All analyses were implemented in Python 3.10+ using a custom Streamlit application. Key dependencies: scanpy 1.9.6+ (Wolf et al., 2018), anndata 0.10+, scipy 1.11+ (Virtanen et al., 2020), pandas 2.0+, numpy 1.24+, matplotlib 3.7+ (Hunter, 2007), statsmodels 0.14+, and adjustText 0.8+. The HypoMap atlas was handled in sparse matrix format throughout to accommodate its size (~3.9 GB). Expression submatrices were extracted in chunks of 200 genes to limit peak memory usage. The application source code is available at [repository URL].

---

## References

Aibar, S. et al. SCENIC: single-cell regulatory network inference and clustering. *Nat. Methods* **14**, 1083–1086 (2017).

Benjamini, Y. & Hochberg, Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J. R. Stat. Soc. B* **57**, 289–300 (1995).

Hunter, J. D. Matplotlib: a 2D graphics environment. *Comput. Sci. Eng.* **9**, 90–95 (2007).

Lawson, C. L. & Hanson, R. J. *Solving Least Squares Problems*. (SIAM, 1995).

Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol.* **15**, 550 (2014).

Steuernagel, L. et al. HypoMap—a unified single-cell gene expression atlas of the murine hypothalamus. *Nat. Metab.* **4**, 1402–1419 (2022).

Subramanian, A. et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proc. Natl Acad. Sci. USA* **102**, 15545–15550 (2005).

Virtanen, P. et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nat. Methods* **17**, 261–272 (2020).

Wolf, F. A., Angerer, P. & Theis, F. J. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol.* **19**, 15 (2018).
