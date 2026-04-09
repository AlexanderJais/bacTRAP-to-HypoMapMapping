"""
Core analysis functions for bacTRAP-to-HypoMap mapping.

Includes: enrichment correlation, marker gene overlap (Fisher's exact),
UMAP enrichment scoring, marker gene computation, NNLS deconvolution,
GSEA-style enrichment, and AUCell scoring.
"""

import logging
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats, sparse
from scipy.optimize import nnls
from typing import Tuple, List, Dict, Optional
import warnings

logger = logging.getLogger(__name__)

from data_loading import (
    compute_cluster_mean_expression,
    compute_fraction_expressing,
    get_gene_names_from_adata,
    _looks_like_ensembl,
    _find_symbol_column,
)


def compute_enrichment_correlation(
    bactrap_matched: pd.DataFrame,
    cluster_mean_expr: pd.DataFrame,
    enrichment_col: str = "log2FoldChange",
) -> pd.DataFrame:
    """
    Compute Pearson and Spearman correlation between bacTRAP enrichment
    profile and each cluster's mean expression profile.

    Args:
        bactrap_matched: matched bacTRAP data with _hypomap_gene_name column
        cluster_mean_expr: DataFrame (genes x clusters) of mean expression
        enrichment_col: column in bactrap_matched to use as enrichment metric

    Returns:
        DataFrame with columns: cluster, pearson_r, pearson_pval,
        spearman_r, spearman_pval, sorted by spearman_r descending.
    """
    logger.info("compute_enrichment_correlation: %d bacTRAP genes, %d expr genes, %d clusters",
                len(bactrap_matched), len(cluster_mean_expr), len(cluster_mean_expr.columns))
    # Align genes between bacTRAP and cluster expression
    bt_genes = bactrap_matched["_hypomap_gene_name"].values
    expr_genes = cluster_mean_expr.index.values
    common = np.intersect1d(bt_genes, expr_genes)
    logger.info("  common genes: %d", len(common))

    if len(common) < 3:
        logger.warning("  <3 common genes — skipping correlation")
        return pd.DataFrame(columns=[
            "cluster", "pearson_r", "pearson_pval", "spearman_r", "spearman_pval",
            "n_genes",
        ])

    # Build aligned enrichment vector (deduplicate bacTRAP genes — keep first)
    bt_dedup = bactrap_matched.drop_duplicates(subset="_hypomap_gene_name", keep="first")
    bt_lookup = dict(zip(bt_dedup["_hypomap_gene_name"], bt_dedup[enrichment_col]))
    enrichment = np.array([bt_lookup[g] for g in common], dtype=float)

    # Guard against duplicate index in cluster_mean_expr
    expr_sub = cluster_mean_expr.loc[common]
    if expr_sub.index.duplicated().any():
        expr_sub = expr_sub.groupby(expr_sub.index).mean()
        # Re-align after dedup
        common = np.intersect1d(list(bt_lookup.keys()), expr_sub.index.values)
        enrichment = np.array([bt_lookup[g] for g in common], dtype=float)
        expr_sub = expr_sub.loc[common]

    # Remove genes with NaN enrichment values (e.g. NaN log2FoldChange)
    valid_mask = np.isfinite(enrichment)
    n_nan = int((~valid_mask).sum())
    if n_nan > 0:
        logger.info("  removed %d genes with NaN enrichment values", n_nan)
        enrichment = enrichment[valid_mask]
        expr_sub = expr_sub[valid_mask]
    logger.info("  genes for correlation: %d, clusters: %d", len(enrichment), len(expr_sub.columns))
    if len(enrichment) < 3:
        logger.warning("  <3 valid genes after NaN removal — skipping correlation")
        return pd.DataFrame(columns=[
            "cluster", "pearson_r", "pearson_pval", "spearman_r", "spearman_pval",
            "n_genes",
        ])

    results = []
    for cluster in expr_sub.columns:
        cluster_expr = expr_sub[cluster].values.astype(float)

        # Skip if no variance
        if np.std(enrichment) == 0 or np.std(cluster_expr) == 0:
            continue

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pr, pp = stats.pearsonr(enrichment, cluster_expr)
            sr, sp = stats.spearmanr(enrichment, cluster_expr)

        results.append({
            "cluster": cluster,
            "pearson_r": pr,
            "pearson_pval": pp,
            "spearman_r": sr,
            "spearman_pval": sp,
            "n_genes": len(enrichment),
        })

    df = pd.DataFrame(results)
    if len(df) > 0:
        df = df.sort_values("spearman_r", ascending=False).reset_index(drop=True)
        # Flag clusters where both Pearson and Spearman are significant —
        # concordance between parametric and rank-based tests is stronger evidence.
        df["both_significant"] = (df["pearson_pval"] < 0.05) & (df["spearman_pval"] < 0.05)
        n_both = int(df["both_significant"].sum())
        logger.info("  correlation result: %d clusters, top spearman_r=%.4f (%s), "
                    "%d/%d significant by both Pearson and Spearman (p<0.05)",
                    len(df), df["spearman_r"].iloc[0], df["cluster"].iloc[0],
                    n_both, len(df))
    else:
        logger.warning("  correlation result: 0 clusters (all had zero variance)")
    return df


def get_enriched_genes(
    bactrap_df: pd.DataFrame,
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 1.0,
) -> pd.DataFrame:
    """
    Filter bacTRAP data for significantly enriched genes.

    Returns subset of bactrap_df passing both padj and log2FC thresholds.
    """
    for required_col in ("padj", "log2FoldChange"):
        if required_col not in bactrap_df.columns:
            raise ValueError(
                f"Required column '{required_col}' not found in bacTRAP data. "
                f"Available columns: {list(bactrap_df.columns)}"
            )
    n_total = len(bactrap_df)
    has_padj = bactrap_df["padj"].notna().sum()
    passes_padj = (bactrap_df["padj"] < padj_cutoff).sum()
    passes_fc = (bactrap_df["log2FoldChange"] > log2fc_cutoff).sum()
    mask = (
        (bactrap_df["padj"].notna())
        & (bactrap_df["padj"] < padj_cutoff)
        & (bactrap_df["log2FoldChange"] > log2fc_cutoff)
    )
    n_enriched = mask.sum()
    logger.info("get_enriched_genes: %d total, %d with padj, %d pass padj<%.3f, "
                "%d pass log2FC>%.2f, %d pass both",
                n_total, has_padj, passes_padj, padj_cutoff, passes_fc, log2fc_cutoff, n_enriched)
    return bactrap_df[mask].copy()


def compute_marker_genes(
    adata, annotation_col: str, n_genes: int = 100, min_cells: int = 10
) -> Dict[str, List[str]]:
    """
    Compute marker genes per cluster using scanpy's rank_genes_groups.

    Returns dict mapping cluster name to list of marker gene names.
    """
    logger.info("compute_marker_genes: annotation_col=%s, n_genes=%d, min_cells=%d",
                annotation_col, n_genes, min_cells)
    # Filter to clusters with enough cells
    cluster_counts = adata.obs[annotation_col].value_counts()
    valid_clusters = cluster_counts[cluster_counts >= min_cells].index.tolist()
    logger.info("  %d/%d clusters pass min_cells=%d", len(valid_clusters), len(cluster_counts), min_cells)

    if len(valid_clusters) == 0:
        logger.warning("  no clusters pass min_cells filter — returning empty markers")
        return {}

    adata_sub = adata[adata.obs[annotation_col].isin(valid_clusters)].copy()

    # Use raw if available for marker gene computation
    if adata_sub.raw is not None:
        # Work with raw counts for proper statistical testing
        adata_work = adata_sub.raw.to_adata()
        adata_work.obs = adata_sub.obs.copy()
    else:
        adata_work = adata_sub.copy()

    # Ensure the data is log-normalized for rank_genes_groups
    sc.pp.filter_genes(adata_work, min_cells=1)

    # Check if data needs normalization by sampling a small subset to
    # avoid the cost of computing max() on the full sparse matrix.
    # Seeded for reproducibility so the normalize/log1p decision is stable.
    sample_size = min(1000, adata_work.n_obs)
    sample_idx = np.random.default_rng(42).choice(
        adata_work.n_obs, sample_size, replace=False
    )
    X_sample = adata_work.X[sample_idx, :]
    if sparse.issparse(X_sample):
        sample_max = X_sample.max()
    else:
        sample_max = np.max(X_sample)
    if sample_max > 50:
        logger.info("  sample_max=%.1f > 50 — applying normalize_total + log1p", float(sample_max))
        sc.pp.normalize_total(adata_work, target_sum=1e4)
        sc.pp.log1p(adata_work)
    else:
        logger.info("  sample_max=%.1f <= 50 — data appears already normalized", float(sample_max))

    sc.tl.rank_genes_groups(
        adata_work,
        groupby=annotation_col,
        method="wilcoxon",
        n_genes=n_genes,
        use_raw=False,
    )

    # Build a var_name -> gene_symbol mapping so that returned marker names
    # are gene symbols (matching how match_genes identifies genes), not
    # Ensembl IDs which adata_work.var_names may contain when built from raw.
    work_gene_names = get_gene_names_from_adata(adata_work)
    varname_to_symbol = {}
    for vn, gn in zip(adata_work.var_names, work_gene_names):
        varname_to_symbol[str(vn)] = str(gn)

    markers = {}
    n_failed = 0
    for cluster in valid_clusters:
        try:
            raw_names = sc.get.rank_genes_groups_df(
                adata_work, group=str(cluster)
            )["names"].tolist()[:n_genes]
            # Convert var_names to gene symbols where possible
            markers[str(cluster)] = [
                varname_to_symbol.get(g, g) for g in raw_names
            ]
        except Exception:
            n_failed += 1
            continue

    logger.info("  marker genes computed for %d clusters (%d failed)", len(markers), n_failed)
    return markers


def load_precomputed_markers(adata) -> Optional[Dict[str, List[str]]]:
    """Try to load pre-computed marker genes from adata.uns.

    If the stored gene names look like Ensembl IDs, they are converted to
    gene symbols using ``adata.var`` (or ``adata.raw.var``) metadata so that
    downstream comparisons (Fisher's, GSEA) work correctly against
    symbol-based enriched-gene lists.
    """
    logger.info("load_precomputed_markers: checking adata.uns for rank_genes_groups")
    if "rank_genes_groups" not in adata.uns:
        logger.info("  not found in adata.uns")
        return None
    try:
        groups = adata.uns["rank_genes_groups"]["names"].dtype.names
        markers = {}
        for group in groups:
            genes = adata.uns["rank_genes_groups"]["names"][group].tolist()
            markers[str(group)] = [str(g) for g in genes]
        logger.info("  loaded %d groups, %d genes/group (first group)", len(markers),
                    len(next(iter(markers.values()))) if markers else 0)
    except Exception as e:
        logger.warning("  failed to parse rank_genes_groups: %s", e)
        return None

    # ---- Ensembl → symbol conversion if needed ----
    sample_genes = []
    for glist in markers.values():
        sample_genes.extend(glist[:5])
        if len(sample_genes) >= 10:
            break
    if any(str(g).startswith(("ENSMUSG", "ENSG")) for g in sample_genes):
        logger.info("  marker genes appear to be Ensembl IDs — converting to symbols")
        # Build lookup from var (and raw.var if available)
        ensembl_to_symbol: Dict[str, str] = {}
        for source_var in ([adata.raw.var] if adata.raw is not None else []) + [adata.var]:
            sym_col = _find_symbol_column(source_var)
            if sym_col is not None:
                for ens_id, sym in zip(source_var.index, source_var[sym_col]):
                    sym_str = str(sym).strip()
                    if sym_str and sym_str.lower() != "nan":
                        ensembl_to_symbol[str(ens_id)] = sym_str

        if ensembl_to_symbol:
            logger.info("  Ensembl→symbol map: %d entries", len(ensembl_to_symbol))
            markers = {
                cluster: [ensembl_to_symbol.get(g, g) for g in genes]
                for cluster, genes in markers.items()
            }
        else:
            logger.warning("  no symbol mapping found — markers remain as Ensembl IDs")
    else:
        logger.info("  marker genes are already symbols")

    return markers


def fisher_overlap_test(
    enriched_genes: List[str],
    cluster_markers: Dict[str, List[str]],
    universe_size: int,
) -> pd.DataFrame:
    """
    For each cluster, compute overlap between bacTRAP enriched genes
    and cluster markers using a one-sided Fisher's exact test.

    Args:
        enriched_genes: list of bacTRAP-enriched gene names
        cluster_markers: dict of cluster -> marker gene list
        universe_size: total number of genes in the overlap universe

    Returns:
        DataFrame with columns: cluster, overlap_count, overlap_genes,
        odds_ratio, pvalue, neg_log10_pval, sorted by pvalue.
    """
    logger.info("fisher_overlap_test: %d enriched genes, %d clusters, universe=%d",
                len(enriched_genes), len(cluster_markers), universe_size)
    if len(enriched_genes) == 0 or len(cluster_markers) == 0:
        logger.warning("  empty input — returning empty results")
        return pd.DataFrame(columns=[
            "cluster", "overlap_count", "n_enriched", "n_markers",
            "overlap_genes", "odds_ratio", "pvalue", "neg_log10_pval",
            "log2_odds_ratio", "padj",
        ])

    enriched_set = set(g.lower() for g in enriched_genes)
    n_enriched = len(enriched_set)
    # Build case-mapping once outside the loop
    enriched_original = {g.lower(): g for g in enriched_genes}

    results = []
    for cluster, markers in cluster_markers.items():
        marker_set = set(g.lower() for g in markers)
        n_markers = len(marker_set)

        overlap = enriched_set & marker_set
        n_overlap = len(overlap)

        # Contingency table for Fisher's exact test
        a = n_overlap
        b = n_enriched - n_overlap
        c = n_markers - n_overlap
        d = max(universe_size - n_enriched - n_markers + n_overlap, 0)

        table = np.array([[a, b], [c, d]])
        odds_ratio, pvalue = stats.fisher_exact(table, alternative="greater")

        overlap_names = sorted([enriched_original.get(g, g) for g in overlap])

        results.append({
            "cluster": cluster,
            "overlap_count": n_overlap,
            "n_enriched": n_enriched,
            "n_markers": n_markers,
            "overlap_genes": ", ".join(overlap_names),
            "odds_ratio": odds_ratio if np.isfinite(odds_ratio) else 999.0,
            "pvalue": pvalue,
            "neg_log10_pval": -np.log10(max(pvalue, 1e-300)),
            "log2_odds_ratio": np.log2(odds_ratio) if odds_ratio > 0 and np.isfinite(odds_ratio) else 10.0,
        })

    df = pd.DataFrame(results)
    if len(df) > 0:
        # FDR correction
        from statsmodels.stats.multitest import multipletests
        _, padj, _, _ = multipletests(df["pvalue"].values, method="fdr_bh")
        df["padj"] = padj
        df = df.sort_values("pvalue").reset_index(drop=True)
        n_sig = (df["padj"] < 0.05).sum()
        logger.info("  Fisher result: %d clusters tested, %d significant (padj<0.05)", len(df), n_sig)
    return df


def compute_enrichment_score(
    adata,
    gene_names: List[str],
    score_name: str = "bacTRAP_enrichment",
    use_raw: bool = True,
) -> np.ndarray:
    """
    Compute a bacTRAP enrichment score for each cell.

    This is a z-scored mean expression of the given gene set across all cells.
    Gene lookup uses the **raw** layer when available so that all genes (not
    just the HVG-filtered subset in ``adata.var``) are considered.
    """
    # Determine expression source
    if use_raw and adata.raw is not None:
        X = adata.raw.X
        source_gene_names = get_gene_names_from_adata(adata, use_raw=True)
    else:
        X = adata.X
        source_gene_names = get_gene_names_from_adata(adata)

    # Build case-insensitive lookup: gene symbol -> column index
    gene_lower_to_idx = {}
    for i, g in enumerate(source_gene_names):
        key = str(g).strip().lower()
        if key and key != "nan":
            gene_lower_to_idx[key] = i

    # Map input genes to expression matrix column indices
    gene_idx = []
    for g in gene_names:
        g_lower = str(g).strip().lower()
        if g_lower in gene_lower_to_idx:
            gene_idx.append(gene_lower_to_idx[g_lower])

    logger.info("compute_enrichment_score: %d/%d genes found in atlas, %d cells",
                len(gene_idx), len(gene_names), adata.n_obs)
    if len(gene_idx) == 0:
        logger.warning("  no genes found — returning zero scores")
        return np.zeros(adata.n_obs)

    X_sub = X[:, gene_idx]
    if sparse.issparse(X_sub):
        X_sub = np.asarray(X_sub.toarray())
    else:
        X_sub = np.asarray(X_sub)

    scores = X_sub.mean(axis=1).flatten()

    # Z-score
    std = np.std(scores)
    if std > 0:
        scores = (scores - np.mean(scores)) / std

    return scores


def compute_zscore_heatmap_data(
    cluster_mean_expr: pd.DataFrame,
    top_genes: List[str],
    top_clusters: List[str],
) -> pd.DataFrame:
    """
    Extract and z-score expression data for a heatmap.

    Args:
        cluster_mean_expr: (genes x clusters) mean expression
        top_genes: gene names for rows
        top_clusters: cluster names for columns

    Returns:
        Z-scored DataFrame (genes x clusters).
    """
    # Filter to available genes and clusters
    available_genes = [g for g in top_genes if g in cluster_mean_expr.index]
    available_clusters = [c for c in top_clusters if c in cluster_mean_expr.columns]

    if len(available_genes) == 0 or len(available_clusters) == 0:
        return pd.DataFrame()

    sub = cluster_mean_expr.loc[available_genes, available_clusters].copy()

    # Z-score across clusters (rows)
    row_means = sub.mean(axis=1)
    row_stds = sub.std(axis=1)
    row_stds = row_stds.replace(0, 1)  # avoid division by zero
    zscored = sub.subtract(row_means, axis=0).divide(row_stds, axis=0)

    return zscored


# =========================================================================
# NNLS Deconvolution
# =========================================================================

def compute_nnls_deconvolution(
    bactrap_matched: pd.DataFrame,
    cluster_mean_expr: pd.DataFrame,
    enrichment_col: str = "log2FoldChange",
) -> pd.DataFrame:
    """
    Non-negative least squares deconvolution of the bacTRAP enrichment
    profile onto a cluster *specificity* signature matrix.

    The bacTRAP enrichment vector ``b`` (log2FC, IP vs Input) is a
    ratio/deviation quantity.  Fitting it against a raw cluster mean
    expression matrix is ill-posed because ``b`` and a raw-counts ``A``
    live in different spaces — clusters with high per-cell library depth
    (ependymal, endothelial, stromal) acquire inflated column norms and
    dominate the solution regardless of biology.

    To put ``A`` and ``b`` in comparable spaces we:

    1. Expect ``cluster_mean_expr`` to already be log-normalized (the
       upstream helper applies normalize_total + log1p).
    2. Row-standardize ``A`` across clusters per gene (z-score) so each
       row encodes *how cluster-specific* the gene is — the same kind of
       deviation quantity as ``log2FC``.
    3. Center ``b`` around its mean so NNLS fits relative enrichment
       rather than the absolute positive offset (all enriched genes have
       positive log2FC, which would otherwise bias every cluster up).

    Solves: ``min ||A_z @ w - b_c||_2``  subject to ``w >= 0``.

    Returns:
        DataFrame with columns: cluster, weight, weight_norm (0–1 scaled),
        sorted by weight descending.
    """
    logger.info("compute_nnls_deconvolution: %d bacTRAP genes, %d expr genes, %d clusters",
                len(bactrap_matched), len(cluster_mean_expr), len(cluster_mean_expr.columns))
    bt_dedup = bactrap_matched.drop_duplicates(subset="_hypomap_gene_name", keep="first")
    bt_genes = bt_dedup["_hypomap_gene_name"].values
    expr_genes = cluster_mean_expr.index.values
    common = np.intersect1d(bt_genes, expr_genes)
    logger.info("  common genes: %d", len(common))

    if len(common) < 5:
        logger.warning("  <5 common genes — skipping NNLS")
        return pd.DataFrame(columns=["cluster", "weight", "weight_norm"])

    bt_lookup = dict(zip(bt_dedup["_hypomap_gene_name"], bt_dedup[enrichment_col]))
    enrichment = np.array([bt_lookup[g] for g in common], dtype=float)
    valid = np.isfinite(enrichment)
    enrichment = enrichment[valid]
    common = common[valid]

    if len(enrichment) < 5:
        return pd.DataFrame(columns=["cluster", "weight", "weight_norm"])

    expr_sub = cluster_mean_expr.loc[common]
    if expr_sub.index.duplicated().any():
        expr_sub = expr_sub.groupby(expr_sub.index).mean()
        common = np.intersect1d(common, expr_sub.index.values)
        enrichment = np.array([bt_lookup[g] for g in common], dtype=float)
        expr_sub = expr_sub.loc[common]

    A_raw = expr_sub.values.astype(float)  # (genes, clusters)

    # ---- Row-wise z-score (cluster specificity per gene) -------------------
    # For each gene, compute how many std-devs each cluster is from the
    # gene's across-cluster mean.  Genes with no variance across clusters
    # (uniformly expressed or uniformly absent) are dropped — they carry no
    # information for deconvolution.
    row_mean = A_raw.mean(axis=1, keepdims=True)
    row_std = A_raw.std(axis=1, keepdims=True)
    informative = (row_std.ravel() > 1e-8)
    n_dropped = int((~informative).sum())
    if n_dropped > 0:
        logger.info("  dropping %d/%d genes with zero cross-cluster variance",
                    n_dropped, len(informative))
    A_raw = A_raw[informative]
    row_mean = row_mean[informative]
    row_std = row_std[informative]
    enrichment = enrichment[informative]
    common = common[informative]

    if len(common) < 5:
        logger.warning("  <5 informative genes — skipping NNLS")
        return pd.DataFrame(columns=["cluster", "weight", "weight_norm"])

    A = (A_raw - row_mean) / row_std  # z-score across clusters, per gene

    # ---- Centre b so NNLS fits relative (not offset) enrichment ------------
    # All values in b (log2FC for enriched genes) are positive; subtracting
    # the mean removes the global offset that NNLS would otherwise try to
    # absorb by spreading weight across many clusters.
    b = enrichment.astype(float) - float(enrichment.mean())

    logger.info("  NNLS input: A=%s (z-scored), b=%s (centered), "
                "b range=[%.3f, %.3f]",
                A.shape, b.shape, float(b.min()), float(b.max()))
    w, residual = nnls(A, b)
    n_nonzero = int((w > 0).sum())
    logger.info("  NNLS result: residual=%.4f, %d/%d clusters with nonzero weight, max_weight=%.4f",
                residual, n_nonzero, len(w), float(w.max()) if len(w) > 0 else 0)

    clusters = expr_sub.columns.tolist()
    df = pd.DataFrame({
        "cluster": clusters,
        "weight": w,
    })
    total = df["weight"].sum()
    df["weight_norm"] = df["weight"] / total if total > 0 else 0.0
    df["global_residual_norm"] = residual
    df = df.sort_values("weight", ascending=False).reset_index(drop=True)
    return df


# =========================================================================
# GSEA-style Preranked Enrichment
# =========================================================================

def _running_enrichment_score_vec(
    hit_mask: np.ndarray,
    abs_enrichment: np.ndarray,
) -> Tuple[float, np.ndarray]:
    """Vectorised running enrichment score (Subramanian et al., PNAS 2005).

    Args:
        hit_mask: boolean array — True where the ranked gene belongs to the set
        abs_enrichment: absolute enrichment values aligned with the ranked list

    Returns:
        (ES, running_scores) — peak enrichment score and the full running curve
    """
    N = len(hit_mask)
    n_hit = hit_mask.sum()
    n_miss = N - n_hit

    if n_hit == 0 or n_miss == 0:
        return 0.0, np.zeros(N)

    hit_weights = np.where(hit_mask, abs_enrichment, 0.0)
    weight_sum = hit_weights.sum()
    if weight_sum == 0:
        weight_sum = 1.0

    miss_penalty = 1.0 / n_miss
    increments = np.where(hit_mask, hit_weights / weight_sum, -miss_penalty)
    running = np.cumsum(increments)

    max_pos = running.max()
    max_neg = running.min()
    es = max_pos if abs(max_pos) >= abs(max_neg) else max_neg
    return es, running


def _es_from_mask(hit_mask: np.ndarray, abs_enrichment: np.ndarray) -> float:
    """Fast ES computation (no running curve returned)."""
    n_hit = hit_mask.sum()
    n_miss = len(hit_mask) - n_hit
    if n_hit == 0 or n_miss == 0:
        return 0.0
    hit_weights = np.where(hit_mask, abs_enrichment, 0.0)
    weight_sum = hit_weights.sum()
    if weight_sum == 0:
        weight_sum = 1.0
    increments = np.where(hit_mask, hit_weights / weight_sum, -1.0 / n_miss)
    running = np.cumsum(increments)
    max_pos = running.max()
    max_neg = running.min()
    return max_pos if abs(max_pos) >= abs(max_neg) else max_neg


def compute_gsea_enrichment(
    bactrap_matched: pd.DataFrame,
    cluster_markers: Dict[str, List[str]],
    enrichment_col: str = "log2FoldChange",
    n_perm: int = 1000,
) -> pd.DataFrame:
    """
    Preranked GSEA: rank all matched genes by bacTRAP enrichment, then
    compute an enrichment score for each cluster's marker gene set.

    Uses vectorised numpy operations for the running-score computation
    and batch-generates all permutation indices up front, so the
    185-cluster × 1 000-permutation workload stays tractable.

    Returns:
        (DataFrame, running_scores_dict, ranked_genes)
    """
    empty_result = (
        pd.DataFrame(columns=["cluster", "ES", "NES", "pvalue", "padj", "n_hits"]),
        {},
        np.array([]),
    )
    logger.info("compute_gsea_enrichment: %d genes, %d clusters, n_perm=%d",
                len(bactrap_matched), len(cluster_markers), n_perm)
    if len(bactrap_matched) == 0 or len(cluster_markers) == 0:
        logger.warning("  empty input — returning empty GSEA results")
        return empty_result

    # Rank genes by enrichment (descending)
    df_sorted = bactrap_matched.dropna(subset=[enrichment_col]).sort_values(
        enrichment_col, ascending=False,
    )
    ranked_genes = df_sorted["_hypomap_gene_name"].values
    enrichment_vals = df_sorted[enrichment_col].values.astype(float)
    abs_enrichment = np.abs(enrichment_vals)
    logger.info("  ranked %d genes (after dropna), FC range: [%.2f, %.2f]",
                len(ranked_genes), float(enrichment_vals[-1]), float(enrichment_vals[0]))

    N = len(ranked_genes)
    # Pre-compute lowercase names once (avoids repeated .lower() in loops)
    ranked_lower = np.array([g.lower() for g in ranked_genes])

    # Pre-generate all permutation indices at once
    rng = np.random.default_rng(42)
    perm_indices = np.empty((n_perm, N), dtype=np.intp)
    for p in range(n_perm):
        perm_indices[p] = rng.permutation(N)

    results = []
    running_scores_dict = {}

    for cluster, markers in cluster_markers.items():
        marker_set = set(g.lower() for g in markers)
        # Boolean hit mask — vectorised membership test
        hit_mask = np.array([g in marker_set for g in ranked_lower])
        n_hits = int(hit_mask.sum())

        es, running = _running_enrichment_score_vec(hit_mask, abs_enrichment)
        running_scores_dict[cluster] = running

        # Permutation null — reuse pre-generated indices
        null_es = np.empty(n_perm)
        for p in range(n_perm):
            perm_idx = perm_indices[p]
            null_es[p] = _es_from_mask(hit_mask[perm_idx], abs_enrichment[perm_idx])

        # p-value (one-sided)
        if es >= 0:
            pval = (np.sum(null_es >= es) + 1) / (n_perm + 1)
        else:
            pval = (np.sum(null_es <= es) + 1) / (n_perm + 1)

        null_mean = np.mean(np.abs(null_es))
        nes = es / null_mean if null_mean > 0 else 0.0

        results.append({
            "cluster": cluster,
            "ES": es,
            "NES": nes,
            "pvalue": pval,
            "n_hits": n_hits,
        })

    df = pd.DataFrame(results)
    if len(df) > 0:
        from statsmodels.stats.multitest import multipletests
        _, padj, _, _ = multipletests(df["pvalue"].values, method="fdr_bh")
        df["padj"] = padj
        df = df.sort_values("NES", ascending=False).reset_index(drop=True)
        n_sig = (df["padj"] < 0.05).sum()
        logger.info("  GSEA result: %d clusters, %d significant (padj<0.05), top NES=%.2f (%s)",
                    len(df), n_sig, df["NES"].iloc[0], df["cluster"].iloc[0])
    else:
        logger.warning("  GSEA result: 0 clusters")

    return df, running_scores_dict, ranked_genes


# =========================================================================
# AUCell Scoring
# =========================================================================

def compute_aucell_scores(
    adata,
    gene_names: List[str],
    use_raw: bool = True,
    top_fraction: float = 0.05,
) -> np.ndarray:
    """
    Compute AUCell scores for each cell.

    AUCell (Aibar et al., Nature Methods 2017) ranks genes by expression
    within each cell, then computes the Area Under the recovery Curve (AUC)
    for the gene set of interest within the top-ranked genes.

    This is more robust than simple mean expression because:
    - It's rank-based (insensitive to normalization differences)
    - It focuses on highly expressed genes per cell
    - It's threshold-free

    Args:
        adata: AnnData object
        gene_names: list of bacTRAP-enriched gene names
        use_raw: whether to use adata.raw for expression
        top_fraction: fraction of ranked genes to consider (default 5%)

    Returns:
        Array of AUCell scores, one per cell.
    """
    # Determine expression source and build gene name lookup
    if use_raw and adata.raw is not None:
        X = adata.raw.X
        source_gene_names = get_gene_names_from_adata(adata, use_raw=True)
    else:
        X = adata.X
        source_gene_names = get_gene_names_from_adata(adata)

    gene_lower_to_idx = {str(g).lower(): i for i, g in enumerate(source_gene_names)}

    # Map input genes to expression matrix column indices
    query_idx = []
    for g in gene_names:
        g_lower = str(g).lower()
        if g_lower in gene_lower_to_idx:
            query_idx.append(gene_lower_to_idx[g_lower])

    n_cells = X.shape[0]
    n_total_genes = X.shape[1]
    n_query = len(query_idx)

    logger.info("compute_aucell_scores: %d/%d genes found, %d cells, top_fraction=%.2f",
                n_query, len(gene_names), n_cells, top_fraction)
    if n_query == 0:
        logger.warning("  no genes found — returning zero AUCell scores")
        return np.zeros(n_cells)

    # Boolean mask for query genes (vectorized membership test)
    query_mask = np.zeros(n_total_genes, dtype=bool)
    query_mask[query_idx] = True

    # Number of top genes to consider per cell
    n_top = max(int(n_total_genes * top_fraction), n_query)
    n_top = min(n_top, n_total_genes)
    # Theoretical maximum of sum(cumsum(is_hit)) when all n_query query genes
    # occupy the top n_query positions:
    #   cumsum = [1, 2, ..., n_query, n_query, ..., n_query]  (n_top entries)
    #   sum    = n_query*(n_query+1)/2 + n_query*(n_top - n_query)
    #          = n_query * (n_top - (n_query - 1)/2)
    # Note: n_top >= n_query is guaranteed above, so this is always positive.
    max_auc = n_query * (n_top - (n_query - 1) / 2)

    # Process in cell chunks — vectorized within each chunk
    chunk_size = 5000
    scores = np.zeros(n_cells, dtype=np.float32)

    for start in range(0, n_cells, chunk_size):
        end = min(start + chunk_size, n_cells)
        X_chunk = X[start:end, :]
        if sparse.issparse(X_chunk):
            X_chunk = np.asarray(X_chunk.toarray())
        else:
            X_chunk = np.asarray(X_chunk)

        chunk_n = X_chunk.shape[0]

        # For each cell, get top-n gene indices via argpartition (O(n) per cell)
        # Then check which are query genes and compute cumulative AUC
        top_idx = np.argpartition(X_chunk, -n_top, axis=1)[:, -n_top:]

        for i in range(chunk_n):
            cell_top = top_idx[i]
            # Sort by expression descending
            order = np.argsort(X_chunk[i, cell_top])[::-1]
            sorted_top = cell_top[order]

            # Vectorized: check query membership and cumsum for AUC
            is_hit = query_mask[sorted_top]
            cumhits = np.cumsum(is_hit)
            auc = cumhits.sum()
            scores[start + i] = auc / max_auc if max_auc > 0 else 0.0

    logger.info("  AUCell scores: mean=%.4f, std=%.4f, min=%.4f, max=%.4f",
                float(scores.mean()), float(scores.std()), float(scores.min()), float(scores.max()))
    return scores


# =========================================================================
# Composite Ranking
# =========================================================================

def compute_composite_ranking(
    corr_df: pd.DataFrame,
    fisher_df: pd.DataFrame,
    nnls_df: pd.DataFrame,
    gsea_df: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """
    Combine multiple ranking methods into a single consensus ranking.

    For each method, ranks are converted to percentile scores (0–1),
    then averaged. This produces a robust ranking that doesn't depend
    on any single method's assumptions.

    Returns:
        DataFrame with columns: cluster, corr_rank, fisher_rank, nnls_rank,
        gsea_rank (if available), composite_score, sorted by composite_score.
    """
    # Collect rankings from each method
    rankings = {}

    if len(corr_df) > 0:
        rank_corr = corr_df[["cluster", "spearman_r"]].copy()
        rank_corr["corr_rank"] = rank_corr["spearman_r"].rank(ascending=False, method="min")
        n = len(rank_corr)
        rank_corr["corr_pctl"] = 1 - (rank_corr["corr_rank"] - 1) / max(n - 1, 1)
        rankings["corr"] = rank_corr.set_index("cluster")

    if len(fisher_df) > 0:
        rank_fisher = fisher_df[["cluster", "pvalue"]].copy()
        rank_fisher["fisher_rank"] = rank_fisher["pvalue"].rank(ascending=True, method="min")
        n = len(rank_fisher)
        rank_fisher["fisher_pctl"] = 1 - (rank_fisher["fisher_rank"] - 1) / max(n - 1, 1)
        rankings["fisher"] = rank_fisher.set_index("cluster")

    if len(nnls_df) > 0:
        rank_nnls = nnls_df[["cluster", "weight"]].copy()
        rank_nnls["nnls_rank"] = rank_nnls["weight"].rank(ascending=False, method="min")
        n = len(rank_nnls)
        rank_nnls["nnls_pctl"] = 1 - (rank_nnls["nnls_rank"] - 1) / max(n - 1, 1)
        rankings["nnls"] = rank_nnls.set_index("cluster")

    if gsea_df is not None and len(gsea_df) > 0:
        rank_gsea = gsea_df[["cluster", "NES"]].copy()
        rank_gsea["gsea_rank"] = rank_gsea["NES"].rank(ascending=False, method="min")
        n = len(rank_gsea)
        rank_gsea["gsea_pctl"] = 1 - (rank_gsea["gsea_rank"] - 1) / max(n - 1, 1)
        rankings["gsea"] = rank_gsea.set_index("cluster")

    if len(rankings) == 0:
        return pd.DataFrame(columns=["cluster", "composite_score"])

    # Merge on cluster
    all_clusters = set()
    for r in rankings.values():
        all_clusters.update(r.index)

    rows = []
    for cluster in all_clusters:
        row = {"cluster": cluster}
        pctls = []
        for method, rdf in rankings.items():
            if cluster in rdf.index:
                pctl = rdf.loc[cluster, f"{method}_pctl"]
                row[f"{method}_rank"] = rdf.loc[cluster, f"{method}_rank"]
                row[f"{method}_pctl"] = pctl
                pctls.append(pctl)
            else:
                row[f"{method}_rank"] = np.nan
                row[f"{method}_pctl"] = np.nan
                pctls.append(np.nan)
        # nanmean: only average over methods that have data for this cluster
        row["composite_score"] = float(np.nanmean(pctls))
        rows.append(row)

    df = pd.DataFrame(rows)
    df = df.sort_values("composite_score", ascending=False).reset_index(drop=True)
    return df
