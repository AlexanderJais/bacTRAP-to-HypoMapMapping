# Materials and Methods: bacTRAP-to-HypoMap Transcriptomic Projection

*The text below is written in the style of a Nature journal Methods section. Adjust tense, length, and detail to match your target journal's requirements.*

---

## bacTRAP translational profiling

Translating ribosome affinity purification (bacTRAP) was performed on preoptic area (PoA) tissue from mice expressing a Cre-dependent EGFP-tagged ribosomal subunit (EGFP-L10a) in a cell-type-specific manner. Immunoprecipitated (IP) and total input mRNA were profiled by RNA-seq. Reads were aligned and quantified to obtain per-gene FPKM values across three biological replicates per condition. Differential expression between IP and Input was performed using DESeq2 (Love et al., 2014), yielding log2 fold changes and Benjamini-Hochberg adjusted p-values (padj) for each gene. Genes with padj < 0.05 and log2FC > 1 were classified as significantly enriched in the bacTRAP translational profile.

To suppress the DESeq2 low-count inflation artefact — in which genes with near-zero Input reads receive spuriously large fold changes from pseudocount division — an additional minimum mean IP expression filter (default baseMean ≥ 10) was applied. Signature gene sets used by the downstream scoring methods were ranked by the π-score (|log₂FC| × −log₁₀(padj); Xiao et al., 2014), which balances effect size against statistical significance; raw log₂FC and −log₁₀(padj) were available as alternative ranking metrics.

## Single-cell reference atlas

The HypoMap single-cell RNA-seq atlas of the murine hypothalamus (Steuernagel et al., 2022) was used as the reference dataset (~384,925 cells). The atlas was obtained as an AnnData h5ad file from the CellxGene repository. Cell-type annotations at the desired hierarchical level were selected from the `.obs` metadata. Expression data were accessed from the `.raw` slot when available (preserving all genes prior to feature selection) and maintained in sparse matrix format throughout to manage memory.

## Gene matching

Gene identifiers from the bacTRAP dataset were matched to the HypoMap atlas using a comprehensive bidirectional lookup. The atlas raw layer (containing all 51,216 genes prior to feature selection) was used as the primary matching target. When atlas variable names were Ensembl identifiers (ENSMUSG), gene symbols were resolved from the `feature_name` column of the variable metadata. Both gene symbols and Ensembl IDs were indexed in the lookup dictionary, enabling matching regardless of whether the bacTRAP data uses symbols or Ensembl identifiers. The bacTRAP gene column was auto-detected by testing all candidate columns against the atlas lookup and selecting the column yielding the highest match rate. Duplicate mappings (multiple Ensembl IDs resolving to the same gene symbol) were resolved by retaining the first occurrence. This procedure achieved a 94.8% match rate (26,426 of 27,884 bacTRAP genes).

## Cluster marker genes

Cluster-specific marker genes were identified using the Wilcoxon rank-sum test via scanpy's `rank_genes_groups` function (Wolf et al., 2018), retaining the top 100 markers per cluster (user-configurable). When pre-computed marker genes were available in the atlas `.uns` slot and matched the selected annotation level (>50% cluster name overlap), these were used directly. Clusters with fewer than 10 cells (user-configurable) were excluded. Marker gene names were converted from Ensembl IDs to gene symbols when the raw layer was indexed by Ensembl notation. The marker sets feed the preranked GSEA described below.

## AUCell scoring

Per-cell enrichment of the bacTRAP signature was quantified with AUCell (Aibar et al., 2017). AUCell is a rank-based, threshold-free, normalisation-insensitive measure of whether a gene set is concentrated among the most highly expressed genes of each cell; it has become a standard tool for projecting bulk or cell-type-specific signatures onto large single-cell atlases. The implementation used here follows the AUCell definition of the original paper and the R reference implementation (`AUCell::AUCell_calcAUC`), with the additional steps described below.

**Expression layer and QC.** Scores were computed from the HypoMap raw-count layer (`adata.raw.X`) so that the rank structure of the low-count tail — where most tied values live — is preserved. Prior to scoring, the chosen layer was inspected on a random sample of 500 cells to confirm that values were integers with a maximum above 50; layers failing this check raised a user-visible warning, since the AUCell definition (and in particular its tie-breaking behaviour) is formulated on raw counts.

**Signature matching.** The bacTRAP signature (top-*N* genes ranked by π-score after the DESeq2 enrichment filters, *N* user-configurable, default 50) was resolved against the atlas variable names through the shared bidirectional lookup described under *Gene matching*. The lookup indexes both gene symbols and Ensembl identifiers against the raw-layer column indices, so signature genes match regardless of the namespace used by either side of the mapping. Unmatched signature genes (if any) are logged and surfaced in the application; they contribute nothing to the score.

**Per-cell rank computation.** For each cell *c*, genes were ranked in descending order of raw expression. Ties were broken at random: a small uniform jitter was added to each cell independently, with the jitter magnitude chosen to be strictly smaller than the smallest distinct gap between expression values (0.49 for integer counts, 0.49 × smallest-nonzero-value otherwise). The order of *distinct* values is therefore preserved while the order of *tied* values is randomised. This reproduces the behaviour of the R implementation's `ties.method = "random"` option and eliminates the systematic bias that would otherwise arise from the deterministic tie-breaking of `numpy.argpartition`/`argsort`, which otherwise resolves ties by matrix-column position and therefore biases scores toward genes appearing earlier in the gene-name vector. The jitter RNG was seeded (default `seed = 0`) to make scores reproducible across reruns.

**AUC of the recovery curve.** Let *S* denote the signature gene set of size |*S*| = *n*<sub>query</sub>, *G* the total number of genes in the raw layer, and *k* = max(⌈τ · *G*⌉, *n*<sub>query</sub>) the number of top-ranked genes considered per cell, with τ the user-selected top-ranked fraction (default τ = 0.05, matching the AUCell default). For each cell the recovery curve *C*(*r*) = |{*g* ∈ *S* : rank<sub>*c*</sub>(*g*) ≤ *r*}| (0 ≤ *r* ≤ *k*) was implicitly constructed via a cumulative-sum over an indicator vector on the sorted top-*k* genes, and the discrete area under the curve AUC<sub>*c*</sub> = Σ<sub>*r*=1..*k*</sub> *C*(*r*) was computed. The AUC was normalised by its theoretical maximum — attained when all *n*<sub>query</sub> signature genes occupy the top *n*<sub>query</sub> positions — giving

&nbsp;&nbsp;&nbsp;&nbsp;*A*<sub>c</sub> = AUC<sub>*c*</sub> / [*n*<sub>query</sub> · (*k* − (*n*<sub>query</sub> − 1)/2)]

so that *A*<sub>c</sub> ∈ [0, 1] and comparisons across datasets or parameter choices are on a common scale. When *n*<sub>query</sub> exceeded ⌈τ · *G*⌉ — i.e. the signature was larger than the requested top-ranked fraction — *k* was raised to *n*<sub>query</sub> so that every signature gene could in principle contribute; this silent window-widening was logged and surfaced to the user as a warning, since the effective top-fraction then exceeds the τ slider.

**Computational details.** Cells were processed in chunks of 5 000 to keep the dense float32 expression buffer below ≈0.6 GB on 30 000-gene atlases. Within each chunk `argpartition` was used to identify top-*k* gene indices in linear time per cell; these were then sorted on the jittered expression and passed through the cumulative-sum / normalisation described above.

## Cluster-level enrichment inference

Per-cluster summaries of the AUCell scores were computed and used both for ranking and for formal inference.

**Cluster-level mean and shrinkage-aware ranking.** For each HypoMap cluster *c* with *n*<sub>*c*</sub> cells, the mean AUCell score *μ̂*<sub>*c*</sub>, median, standard deviation, and standard error of the mean (SEM) were tabulated. Clusters with fewer than 20 cells (user-configurable) were excluded from all *top-N* figures (UMAP cluster highlighting, per-cluster barplot, violin plot) because the cluster mean of a 2-cell population is dominated by shrinkage variance and can otherwise claim a top slot purely by chance; the excluded clusters remain available in the per-cluster CSV. This size filter is distinct from, and stricter than, the minimum-cell threshold used during marker analysis. When the optional baseline Cre-driver expression filter was active (see below), clusters failing that filter were additionally removed from the cluster-level barplot (Fig. S2), the violin panel (Fig. 1c), the cell-type UMAP top-15 highlight (Fig. 1b), and the exported per-cluster CSV; the per-cell AUCell scores and the per-cell UMAP (Fig. 1a) were left unmodified as they are not cluster-indexed.

**Significance testing.** For each eligible cluster, the one-sided null hypothesis that its AUCell scores are drawn from the same distribution as the rest of the atlas was tested using Welch's unequal-variance *t*-test (`scipy.stats.ttest_ind(..., equal_var=False)`), converting the two-sided *p*-value to a one-sided *p*-value in the upper tail. Multiple testing across clusters was controlled with the Benjamini–Hochberg procedure (Benjamini and Hochberg, 1995) to yield a per-cluster *q*-value; clusters with *q* < 0.05 were flagged as significant. The resulting per-cluster table (`aucell_per_cluster.csv`) contains `cluster`, `n_cells`, `mean`, `median`, `std`, `sem`, `t_stat`, `pvalue`, `qvalue`, `significant`.

## UMAP enrichment projection

Per-cell AUCell scores were projected onto the HypoMap UMAP embedding. To make small but highly enriched populations visible, a second UMAP panel coloured by cell-type annotation was rendered with only the top-15 AUCell-ranked clusters highlighted (in descending order of cluster mean, after the 20-cell size filter above); all other cells were drawn in light grey beneath the highlighted layer. For all UMAP panels, the conventional x/y axes were replaced with a pair of small axis arrows in the bottom-left corner — the compact convention common in single-cell publications — drawn in axes-fraction coordinates so they remain in the corner regardless of data range. Scatter points were rasterised at 300 DPI while all text elements remained as editable vector objects. For visualisation only, cells were randomly subsampled to 50 000 (user-configurable); all cells were used for scoring and statistical analyses.

## Gene set enrichment analysis (GSEA)

A preranked GSEA approach (Subramanian et al., 2005) was applied to assess whether cluster-specific marker genes were enriched among the most bacTRAP-enriched genes. All matched genes were ranked by bacTRAP log2FC in descending order. For each cluster's marker gene set, a weighted running enrichment score was computed, where hit weights were proportional to the absolute enrichment value at each gene's rank. Significance was assessed by permutation testing (1,000 permutations), in which gene labels and their associated enrichment values were permuted in concert to preserve the rank-value structure while randomizing gene set membership. The enrichment score was normalized (NES) by dividing by the mean absolute enrichment score across permutations to enable comparison across gene sets of different sizes. P-values were corrected using the Benjamini-Hochberg procedure.

## Composite cluster score (AUCell + GSEA)

Two methods are retained for cluster-level mapping: the per-cluster mean AUCell score (capturing per-cell signature recovery rolled up to the cluster) and the preranked GSEA NES (capturing the position of cluster markers in the bacTRAP enrichment ranking). They are orthogonal — AUCell is a per-cell rank-sum on the bacTRAP signature, GSEA is a preranked test on cluster markers — so a cluster prioritised by both is more credible than one prioritised by either alone. Earlier versions of this tool also reported Spearman correlation, a one-sided Fisher's exact marker-overlap test, and non-negative least squares deconvolution; on the data we tested, those rankings were noise-dominated and disagreed with both biology and the AUCell / GSEA result, so they were dropped.

For each method (AUCell mean, GSEA NES), cluster scores were converted to percentile ranks scaled from 0 to 1 (highest score → percentile 1.0). The composite score for each cluster was computed as the arithmetic mean of its two percentile ranks (`nanmean`, so a cluster scored by only one method is not penalised with a zero). Clusters were ranked by descending composite score. When the optional baseline Cre-driver expression filter was active, each method's cluster universe was restricted to the surviving clusters **before** percentile conversion, so the consensus percentiles reflect the restricted population rather than compressed ranks inherited from the full atlas.

## Cre-driver expression sanity check and baseline-expression filter

Per-cluster log-normalised mean expression and fraction expressing were computed for the Cre-driver gene itself (default `Pnoc` for Pnoc-Cre lines; user-configurable). These statistics served two distinct roles.

**Diagnostic (always active).** In the Cre-driver Check tab, clusters ranked by the composite consensus (see *Composite consensus ranking*, below) were annotated with their Cre-driver expression and flagged when the fraction of cells with non-zero raw counts fell below a user-defined threshold (default 5%). This view is intended to discriminate mapping hits that reflect current Cre-driver transcription from hits that may instead reflect Cre lineage tracing in cells that no longer express the driver (a known confound for bacTRAP-NuTRAP lines) or snRNA-seq dropout of low-copy neuropeptide transcripts. The diagnostic does not modify any ranking; flagged clusters remain visible and are merely highlighted.

**Baseline-expression filter (optional).** A separate slider applies a minimum mean (log-normalised) Cre-driver expression threshold. When set above zero, every cluster whose mean Cre-driver expression falls below the floor is dropped from the input of preranked GSEA, the AUCell cluster-level summaries (figs 1b, S2, 1c), the composite (AUCell + GSEA) score, and all corresponding CSV exports before percentile conversion — so the surviving clusters are re-ranked against one another rather than inheriting their global positions. Per-cell AUCell scores and all per-cell outputs (Fig. 1a, `aucell_per_cell.csv`) are unaffected because they carry no cluster identity. To avoid a silent coupling with the marker-test size gate, the cell-count floor used to compute the filter's reference statistics was set to `min(min_cells_per_cluster, min_cells_for_rank)`; this ensures that every cluster that enters the AUCell top-N figures has a Cre-driver statistic available, so the filter discriminates only on expression and never on sanity-side size mismatches. When the threshold is set high enough to discard every cluster, or when the Cre-driver symbol cannot be resolved against the atlas, the filter is disabled and the application renders a persistent warning banner; in those states all tabs revert to the unfiltered behaviour and CSV filenames use their default names. Filtered exports are suffixed with the filter signature (e.g. `composite_ranking_pnoc_ge0p05.csv`) so downstream users can tell a subset from a full table at a glance.

## Figure generation

All figures were generated using matplotlib (v3.7+) following Nature journal specifications: Arial font (sans-serif), 7 pt axis labels, 6 pt tick labels, 8 pt bold panel titles, 0.5 pt axis line width, white background without gridlines, and removal of top and right spines. Colorblind-friendly palettes were used throughout (viridis, magma for sequential data; tab10/tab20 for categorical). UMAP and other embedding panels omit the conventional x/y axes; instead, two short axis arrows labelled "UMAP1" and "UMAP2" are drawn in the bottom-left corner in axes-fraction coordinates — the compact convention now standard in single-cell publications. Scatter plots were rasterised at 300 DPI to manage file size while text elements remained as editable vector objects in PDF (TrueType embedding, fonttype 42) and SVG (native text, fonttype none) exports. Single-column figures were 89 mm (3.5 in) wide; double-column figures were 183 mm (7.2 in) wide.

## Software

All analyses were implemented in Python 3.10+ using a custom Streamlit application. Key dependencies: scanpy 1.9.6+ (Wolf et al., 2018), anndata 0.10+, scipy 1.11+ (Virtanen et al., 2020), pandas 2.0+, numpy 1.24+, matplotlib 3.7+ (Hunter, 2007), statsmodels 0.14+, and adjustText 0.8+. The HypoMap atlas was handled in sparse matrix format throughout to accommodate its size (~3.9 GB). Expression submatrices were extracted in chunks of 200 genes to limit peak memory usage. The application source code is available at https://github.com/AlexanderJais/bacTRAP-to-HypoMapMapping.

---

## Figure Legends

**Figure 1. AUCell enrichment of the bacTRAP signature across the HypoMap hypothalamic cell atlas.**
**(a)** AUCell enrichment score projected onto the HypoMap UMAP embedding (*n* = 384,925 cells scored; *n* = 50,000 rendered). For each cell, the area under the recovery curve for the bacTRAP signature (top *N* π-score-ranked genes after the padj, log₂FC and minimum-IP-expression filters; default *N* = 50; see *bacTRAP translational profiling*) was computed within the top 5 % of genes by expression rank in that cell, then normalised to its theoretical maximum so that scores lie on [0, 1] (Methods, *AUCell scoring*). Colour encodes the normalised AUCell score (magma colormap); the scale is clipped at the 2nd and 98th percentiles to suppress outlier saturation. The bottom-left arrows mark the UMAP1 / UMAP2 axes. Per-cell scores are computed on the full atlas irrespective of the baseline Cre-driver filter.
**(b)** Same UMAP layout as (a), coloured by HypoMap cell-type annotation at the selected hierarchical level. Only the 15 clusters with the highest mean AUCell score (among clusters with ≥ 20 cells, further restricted to clusters passing the optional baseline Cre-driver expression filter when active) are drawn in colour and overplotted on a light-grey background of all remaining cells; this keeps small but highly enriched populations visible while conveying the overall atlas topology. The in-figure legend lists the highlighted clusters in rank order.
**(c)** Violin plots of the full AUCell score distribution within each of the 15 top-ranked clusters from (b) (subject to the same size and baseline-expression filters), ordered from highest (top) to lowest (bottom) cluster mean. Violin colour encodes the cluster mean (magma colormap). Score calls and significance against the rest of the atlas (Welch's one-sided *t*-test, Benjamini–Hochberg FDR) are available in `aucell_per_cluster.csv` (filename suffixed with the filter signature when active).
**(e)** Composite cluster score across AUCell mean and GSEA NES, top 20 clusters. Each method's cluster scores were converted to percentile ranks (0–1) and averaged (`nanmean`). Colour intensity (YlOrRd) represents percentile score; composite scores are annotated to the right of each row. When the baseline filter is active, percentiles are computed on the restricted cluster universe so the displayed scores compare surviving clusters against one another rather than against the full atlas.

**Supplementary Figures:** Cluster-level supplementary panels honour the optional baseline Cre-driver expression filter: when active, the cluster universe on which the ranking is computed is restricted to clusters passing the filter before the top-N selection, and the corresponding CSV download is saved with a suffixed filename that records the filter signature (e.g. `composite_ranking_pnoc_ge0p05.csv`). Per-cell panels are unaffected because they carry no cluster identity.
**(S1)** Horizontal barplot of mean AUCell score per HypoMap cluster (top 25 clusters by mean, clusters with < 20 cells excluded from the ranking). Error bars are standard error of the mean. Bar colour encodes the mean AUCell score (magma colormap). Source table: `aucell_per_cluster.csv`.
**(S2)** Global histogram of per-cell AUCell scores across the atlas. Dashed vertical lines mark the 90th, 95th, and 99th percentiles of the distribution as well as the mean. Cells above the 95th percentile are the candidate pool for belonging to the bacTRAP target population.
**(S3)** Cell-type annotation (left) and AUCell enrichment score (right) rendered as a paired UMAP panel for side-by-side comparison. This is the two-panel version of Figure 1a/b with both panels sharing a colourbar layout suitable for supplementary presentation.
**(S4)** Preranked GSEA running enrichment score curves for the top 5 HypoMap clusters. All matched genes were ranked by bacTRAP log₂ fold change (descending); normalised enrichment scores (NES) are indicated in the legend. Significance was assessed by 1 000-permutation testing with Benjamini–Hochberg FDR correction.
**(S5)** Horizontal barplot of GSEA normalised enrichment scores (NES) for the top 20 clusters. Red bars indicate clusters with significant enrichment (*q* < 0.05); grey bars indicate non-significant clusters.
**(S6)** Cre-driver sanity panel: per-cluster mean log-normalised expression (left) and fraction of cells with non-zero counts (right) for the user-specified Cre-driver gene (default *Pnoc*), shown for the top-ranked clusters from the composite score. The fraction panel marks the user-defined "expressing" threshold (default 5%); clusters meeting both "ranked top" AND "fraction ≥ threshold" are filled in saturated colour, others in grey. This figure is diagnostic only and does not modify any ranking.

---

## References

Aibar, S. et al. SCENIC: single-cell regulatory network inference and clustering. *Nat. Methods* **14**, 1083–1086 (2017).

Benjamini, Y. & Hochberg, Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J. R. Stat. Soc. B* **57**, 289–300 (1995).

Hunter, J. D. Matplotlib: a 2D graphics environment. *Comput. Sci. Eng.* **9**, 90–95 (2007).

Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol.* **15**, 550 (2014).

Steuernagel, L. et al. HypoMap—a unified single-cell gene expression atlas of the murine hypothalamus. *Nat. Metab.* **4**, 1402–1419 (2022).

Subramanian, A. et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proc. Natl Acad. Sci. USA* **102**, 15545–15550 (2005).

Virtanen, P. et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nat. Methods* **17**, 261–272 (2020).

Wolf, F. A., Angerer, P. & Theis, F. J. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol.* **19**, 15 (2018).

Xiao, Y. et al. A novel significance score for gene selection and ranking. *Bioinformatics* **30**, 801–807 (2014).
