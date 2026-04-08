"""
Core analysis functions for bacTRAP-to-HypoMap mapping.

Includes: enrichment correlation, marker gene overlap (Fisher's exact),
UMAP enrichment scoring, and marker gene computation.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats, sparse
from typing import Tuple, List, Dict, Optional
import warnings

from data_loading import (
    compute_cluster_mean_expression,
    compute_fraction_expressing,
    get_gene_names_from_adata,
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
    # Align genes between bacTRAP and cluster expression
    bt_genes = bactrap_matched["_hypomap_gene_name"].values
    expr_genes = cluster_mean_expr.index.values
    common = np.intersect1d(bt_genes, expr_genes)

    if len(common) < 3:
        return pd.DataFrame(columns=[
            "cluster", "pearson_r", "pearson_pval", "spearman_r", "spearman_pval",
            "n_genes",
        ])

    # Build aligned enrichment vector
    bt_lookup = dict(zip(
        bactrap_matched["_hypomap_gene_name"],
        bactrap_matched[enrichment_col],
    ))
    enrichment = np.array([bt_lookup[g] for g in common], dtype=float)
    expr_sub = cluster_mean_expr.loc[common]

    # Remove genes with NaN enrichment values (e.g. NaN log2FoldChange)
    valid_mask = np.isfinite(enrichment)
    if not np.all(valid_mask):
        enrichment = enrichment[valid_mask]
        expr_sub = expr_sub.iloc[valid_mask]
    if len(enrichment) < 3:
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
    mask = (
        (bactrap_df["padj"].notna())
        & (bactrap_df["padj"] < padj_cutoff)
        & (bactrap_df["log2FoldChange"] > log2fc_cutoff)
    )
    return bactrap_df[mask].copy()


def compute_marker_genes(
    adata, annotation_col: str, n_genes: int = 100, min_cells: int = 10
) -> Dict[str, List[str]]:
    """
    Compute marker genes per cluster using scanpy's rank_genes_groups.

    Returns dict mapping cluster name to list of marker gene names.
    """
    # Filter to clusters with enough cells
    cluster_counts = adata.obs[annotation_col].value_counts()
    valid_clusters = cluster_counts[cluster_counts >= min_cells].index.tolist()

    if len(valid_clusters) == 0:
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
    sample_size = min(1000, adata_work.n_obs)
    sample_idx = np.random.choice(adata_work.n_obs, sample_size, replace=False)
    X_sample = adata_work.X[sample_idx, :]
    if sparse.issparse(X_sample):
        sample_max = X_sample.max()
    else:
        sample_max = np.max(X_sample)
    if sample_max > 50:
        sc.pp.normalize_total(adata_work, target_sum=1e4)
        sc.pp.log1p(adata_work)

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
            continue

    return markers


def load_precomputed_markers(adata) -> Optional[Dict[str, List[str]]]:
    """Try to load pre-computed marker genes from adata.uns."""
    if "rank_genes_groups" in adata.uns:
        try:
            groups = adata.uns["rank_genes_groups"]["names"].dtype.names
            markers = {}
            for group in groups:
                genes = adata.uns["rank_genes_groups"]["names"][group].tolist()
                markers[str(group)] = [str(g) for g in genes]
            return markers
        except Exception:
            return None
    return None


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
    if len(enriched_genes) == 0 or len(cluster_markers) == 0:
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
    return df


def compute_enrichment_score(
    adata,
    gene_names: List[str],
    score_name: str = "bacTRAP_enrichment",
    use_raw: bool = True,
) -> np.ndarray:
    """
    Compute a bacTRAP enrichment score for each cell using scanpy's score_genes.

    This is a z-scored mean expression of the given gene set across all cells.
    """
    # Find which genes are present in adata
    adata_genes = get_gene_names_from_adata(adata)
    adata_genes_lower = {str(g).lower(): str(g) for g in adata_genes}

    present_genes = []
    for g in gene_names:
        g_lower = str(g).lower()
        if g_lower in adata_genes_lower:
            present_genes.append(adata_genes_lower[g_lower])

    if len(present_genes) == 0:
        return np.zeros(adata.n_obs)

    # Use var_names as they appear in the adata
    var_names_list = list(adata.var_names)
    adata_gene_names_arr = get_gene_names_from_adata(adata)

    # Map gene names to var_names indices
    gene_name_to_var = {}
    for i, (vn, gn) in enumerate(zip(var_names_list, adata_gene_names_arr)):
        gene_name_to_var[str(gn).lower()] = vn

    score_gene_list = []
    for g in present_genes:
        g_lower = g.lower()
        if g_lower in gene_name_to_var:
            score_gene_list.append(gene_name_to_var[g_lower])

    if len(score_gene_list) == 0:
        return np.zeros(adata.n_obs)

    # Manual z-scored mean — avoids adata.copy() which doubles memory for
    # the full atlas. sc.tl.score_genes requires a copy and uses more RAM
    # than we can afford with a ~3.9GB object.
    # Build index lookup for the correct expression source (var vs raw.var)
    if use_raw and adata.raw is not None:
        X = adata.raw.X
        source_var_names = list(adata.raw.var_names)
    else:
        X = adata.X
        source_var_names = var_names_list

    source_var_to_idx = {vn: i for i, vn in enumerate(source_var_names)}
    gene_idx = [source_var_to_idx[g] for g in score_gene_list if g in source_var_to_idx]
    if len(gene_idx) == 0:
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
