"""
Core analysis functions for bacTRAP-to-HypoMap mapping.

Includes: enrichment correlation, marker gene overlap (Fisher's exact),
UMAP enrichment scoring, marker gene computation, NNLS deconvolution,
GSEA-style enrichment, and AUCell scoring.
"""

import logging
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats, sparse
from scipy.optimize import nnls
from statsmodels.stats.multitest import multipletests
from typing import Tuple, List, Dict, Iterable, Optional

logger = logging.getLogger(__name__)

from data_loading import (
    compute_cluster_mean_expression,
    compute_fraction_expressing,
    get_gene_names_from_adata,
    _build_adata_gene_lookup,
    _looks_like_ensembl,
    _find_symbol_column,
)


def compute_enrichment_correlation(
    bactrap_matched: pd.DataFrame,
    cluster_mean_expr: pd.DataFrame,
    enrichment_col: str = "log2FoldChange",
) -> pd.DataFrame:
    """
    Compute Pearson and Spearman correlation between the bacTRAP enrichment
    profile and each cluster's *specificity* profile.

    ``cluster_mean_expr`` is row-standardised (z-scored across clusters per
    gene) before correlating.  Without this step, correlation against raw
    log-mean expression systematically goes negative: highly cell-type-
    specific markers have large log2FC but low across-cluster *mean*
    expression (because even the target cluster contains many cells that
    don't express them), while low-log2FC housekeeping genes have high mean
    everywhere.  Z-scoring per gene converts "how highly expressed" into
    "how cluster-specific", which is what the log2FC axis is measuring on
    the bacTRAP side.

    Per-cluster p-values are corrected across clusters using Benjamini–
    Hochberg FDR, and significance flags use the adjusted values so the
    "both_significant" column reflects multiple-testing-controlled calls.
    """
    # Stable output schema — every return (success or early exit) must have
    # these exact columns so downstream `"both_significant" in df.columns`
    # checks disappear and consumers (app.py, CSV export) don't need to
    # branch on schema variations.
    _CORR_COLUMNS = [
        "cluster", "pearson_r", "pearson_pval", "spearman_r", "spearman_pval",
        "n_genes", "pearson_padj", "spearman_padj", "both_significant",
    ]

    logger.info("compute_enrichment_correlation: %d bacTRAP genes, %d expr genes, %d clusters",
                len(bactrap_matched), len(cluster_mean_expr), len(cluster_mean_expr.columns))
    # Align genes between bacTRAP and cluster expression
    bt_genes = bactrap_matched["_hypomap_gene_name"].values
    expr_genes = cluster_mean_expr.index.values
    common = np.intersect1d(bt_genes, expr_genes)
    logger.info("  common genes: %d", len(common))

    if len(common) < 3:
        logger.warning("  <3 common genes — skipping correlation")
        return pd.DataFrame(columns=_CORR_COLUMNS)

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
        return pd.DataFrame(columns=_CORR_COLUMNS)

    # Row-standardise expression across clusters per gene so correlation
    # measures cluster *specificity* rather than absolute expression.
    expr_values = expr_sub.values.astype(float)
    row_mean = expr_values.mean(axis=1, keepdims=True)
    row_std = expr_values.std(axis=1, keepdims=True)
    informative = (row_std.ravel() > 1e-8)
    n_uninformative = int((~informative).sum())
    if n_uninformative > 0:
        logger.info("  dropping %d/%d genes with zero cross-cluster variance",
                    n_uninformative, len(informative))
        expr_values = expr_values[informative]
        row_mean = row_mean[informative]
        row_std = row_std[informative]
        enrichment = enrichment[informative]
    if len(enrichment) < 3:
        logger.warning("  <3 informative genes — skipping correlation")
        return pd.DataFrame(columns=_CORR_COLUMNS)
    expr_z = (expr_values - row_mean) / row_std
    expr_sub = pd.DataFrame(
        expr_z, index=expr_sub.index[informative], columns=expr_sub.columns,
    )

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
    if len(df) == 0:
        # Preserve the schema even when every cluster had zero variance —
        # downstream code can rely on the column set without guards.
        return pd.DataFrame(columns=_CORR_COLUMNS)

    # BH-FDR across clusters — without this, "significance" is inflated
    # because each cluster's p-value is independently computed over
    # ~hundreds of genes and many will cross p<0.05 under weak signal.
    _, pearson_padj, _, _ = multipletests(df["pearson_pval"].values, method="fdr_bh")
    _, spearman_padj, _, _ = multipletests(df["spearman_pval"].values, method="fdr_bh")
    df["pearson_padj"] = pearson_padj
    df["spearman_padj"] = spearman_padj
    df = df.sort_values("spearman_r", ascending=False).reset_index(drop=True)
    # Flag clusters where both Pearson and Spearman pass FDR — concordance
    # between parametric and rank-based tests under multiple-testing
    # correction is stronger evidence than either nominal p-value alone.
    df["both_significant"] = (df["pearson_padj"] < 0.05) & (df["spearman_padj"] < 0.05)
    n_both = int(df["both_significant"].sum())
    logger.info("  correlation result: %d clusters, top spearman_r=%.4f (%s), "
                "%d/%d significant by both Pearson and Spearman (FDR<0.05)",
                len(df), df["spearman_r"].iloc[0], df["cluster"].iloc[0],
                n_both, len(df))
    return df


def get_enriched_genes(
    bactrap_df: pd.DataFrame,
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 1.0,
    min_ip_expression: float = 0.0,
    ip_col: str = "IP",
) -> pd.DataFrame:
    """
    Filter bacTRAP data for significantly enriched genes.

    Returns subset of bactrap_df passing padj, log2FC, and (optionally)
    minimum IP expression thresholds.

    The *min_ip_expression* filter suppresses the DESeq2 low-count
    artefact in which a gene with near-zero Input counts gets an
    inflated log₂FC from pseudocount division (e.g. 28 IP reads vs
    0 Input → log₂FC ~7).  These genes would dominate any
    log₂FC-ranked top-N list despite being statistically weak.  When
    ``ip_col`` is not present in the DataFrame the filter is skipped
    with a warning.
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
    ip_filter_applied = False
    if min_ip_expression > 0.0:
        if ip_col in bactrap_df.columns:
            ip_values = pd.to_numeric(bactrap_df[ip_col], errors="coerce")
            ip_mask = ip_values.fillna(0) >= min_ip_expression
            n_before_ip = mask.sum()
            mask = mask & ip_mask
            ip_filter_applied = True
            logger.info(
                "get_enriched_genes: IP filter '%s' >= %.3g removed %d/%d genes "
                "(%d → %d after filter)",
                ip_col, min_ip_expression,
                n_before_ip - mask.sum(), n_before_ip,
                n_before_ip, mask.sum(),
            )
        else:
            logger.warning(
                "get_enriched_genes: min_ip_expression=%.3g requested but "
                "column '%s' not in bacTRAP data — filter SKIPPED. "
                "Available columns: %s",
                min_ip_expression, ip_col, list(bactrap_df.columns),
            )
    n_enriched = mask.sum()
    logger.info("get_enriched_genes: %d total, %d with padj, %d pass padj<%.3f, "
                "%d pass log2FC>%.2f, %d pass both%s",
                n_total, has_padj, passes_padj, padj_cutoff, passes_fc, log2fc_cutoff, n_enriched,
                f" (after min_IP>={min_ip_expression:.3g} filter)" if ip_filter_applied else "")
    return bactrap_df[mask].copy()


def filter_signature_genes_by_atlas(
    candidate_genes: List[str],
    cluster_mean_expr: pd.DataFrame,
    cell_detection_rate: pd.Series,
    *,
    min_detection_rate: float = 0.02,
    min_max_cluster_mean: float = 0.05,
    specificity_cluster_mean_thresh: float = 0.5,
    specificity_max_cluster_fraction: float = 0.5,
    apply_detectability: bool = True,
    apply_specificity: bool = True,
    logger=None,
) -> Tuple[List[str], pd.DataFrame]:
    """Pre-filter candidate signature genes against HypoMap atlas statistics.

    Two operationally-defined filters, intended to be applied to the bacTRAP
    candidate signature *after* the DESeq2 padj / log₂FC / min-IP filters and
    *before* π-score ranking / top-N selection:

    1. **HypoMap detectability** — drop genes that are essentially absent from
       the atlas: cell detection rate < ``min_detection_rate`` OR maximum
       per-cluster mean log-norm expression < ``min_max_cluster_mean``.  A
       gene that never enters any cell's top-τ window contributes nothing to
       AUCell and only adds noise to the signature.
    2. **Specificity** — drop genes that are highly expressed across most
       clusters: per-cluster mean log-norm expression exceeds
       ``specificity_cluster_mean_thresh`` in more than
       ``specificity_max_cluster_fraction`` of clusters.  Pan-neuronal genes
       (Snap25, Syt1, Stmn2, Tubb3, Map2, …) are enriched in any neuronal IP
       but carry no cell-type specificity; "expressed broadly across clusters"
       is a clean operational proxy that needs no curated list.

    Parameters
    ----------
    candidate_genes
        Candidate signature gene symbols, in atlas namespace.
    cluster_mean_expr
        Log-normalised mean expression as **genes (index) × clusters
        (columns)**.  If the candidate genes appear to live on the columns
        instead, the frame is transposed automatically.  Genes absent from
        the frame are treated as undetectable (max-cluster-mean = 0).
    cell_detection_rate
        Series mapping gene symbol → fraction of atlas cells with raw count
        > 0.  Genes absent from the series are treated as detection rate 0.
    apply_detectability, apply_specificity
        Toggle each filter independently (a disabled filter never drops a
        gene but its statistics still appear in the drop log).

    Returns
    -------
    (kept_genes, drop_log_df)
        ``kept_genes`` preserves the input order.  ``drop_log_df`` has one
        row per candidate gene with columns ``gene``, ``status`` ("kept" or a
        ``;``-joined drop-reason string), ``cell_detection_rate``,
        ``max_cluster_mean`` and ``frac_clusters_above_thresh``.
    """
    log = logger or logging.getLogger(__name__)
    cme = cluster_mean_expr
    if cme is not None and len(cme.index) and len(cme.columns):
        cand_set = {str(g) for g in candidate_genes}
        n_on_index = len(cand_set & {str(i) for i in cme.index})
        n_on_cols = len(cand_set & {str(c) for c in cme.columns})
        if n_on_cols > n_on_index:
            cme = cme.T

    if cell_detection_rate is None:
        cell_detection_rate = pd.Series(dtype=float)

    rows = []
    kept: List[str] = []
    for g in candidate_genes:
        det = float(cell_detection_rate.get(g, 0.0))
        if cme is not None and g in cme.index:
            row_vals = cme.loc[g]
            if isinstance(row_vals, pd.DataFrame):  # duplicate index guard
                row_vals = row_vals.mean(axis=0)
            row_arr = pd.to_numeric(row_vals, errors="coerce").to_numpy(dtype=float)
            finite = row_arr[np.isfinite(row_arr)]
            max_cluster_mean = float(finite.max()) if finite.size else 0.0
            n_clusters = int(finite.size)
            n_above = int((finite > specificity_cluster_mean_thresh).sum())
            frac_above = (n_above / n_clusters) if n_clusters else 0.0
        else:
            max_cluster_mean = 0.0
            frac_above = 0.0
        reasons: List[str] = []
        if apply_detectability:
            if det < min_detection_rate:
                reasons.append(f"low_detection_rate(<{min_detection_rate:g})")
            if max_cluster_mean < min_max_cluster_mean:
                reasons.append(f"low_max_cluster_mean(<{min_max_cluster_mean:g})")
        if apply_specificity and frac_above > specificity_max_cluster_fraction:
            reasons.append(
                f"broadly_expressed(>{specificity_cluster_mean_thresh:g} in "
                f">{specificity_max_cluster_fraction:.0%} of clusters)"
            )
        if reasons:
            status = "; ".join(reasons)
        else:
            status = "kept"
            kept.append(g)
        rows.append({
            "gene": g,
            "status": status,
            "cell_detection_rate": det,
            "max_cluster_mean": max_cluster_mean,
            "frac_clusters_above_thresh": frac_above,
        })

    drop_log_df = pd.DataFrame(rows, columns=[
        "gene", "status", "cell_detection_rate", "max_cluster_mean",
        "frac_clusters_above_thresh",
    ])
    n_dropped = len(candidate_genes) - len(kept)
    log.info(
        "filter_signature_genes_by_atlas: %d candidates -> %d kept, %d dropped "
        "(detectability=%s, specificity=%s; min_det=%.3g, min_max_mean=%.3g, "
        "spec_thresh=%.3g, spec_frac=%.3g)",
        len(candidate_genes), len(kept), n_dropped,
        apply_detectability, apply_specificity,
        min_detection_rate, min_max_cluster_mean,
        specificity_cluster_mean_thresh, specificity_max_cluster_fraction,
    )
    if n_dropped:
        _by_reason = drop_log_df.loc[drop_log_df["status"] != "kept", "status"].value_counts()
        log.info("  drop reasons: %s", _by_reason.to_dict())
    return kept, drop_log_df


def rank_enriched_genes(
    df: pd.DataFrame,
    metric: str = "pi_score",
) -> pd.DataFrame:
    """
    Sort an enriched-gene DataFrame by the chosen ranking metric,
    most-enriched first.

    Metrics
    -------
    ``"pi_score"`` (recommended default)
        ``|log₂FC| × -log₁₀(padj)`` — the π-score of Xiao et al.
        (Bioinformatics 2014). Combines effect size and significance so
        a gene with a modest but highly significant fold-change (e.g.
        Pnoc in a Pnoc-Cre line: FC=2, padj=1e-17) outranks a low-count
        zero-Input artefact (e.g. FC=7, padj=1e-4).
    ``"log2fc"``
        Raw ``log₂FoldChange``. Historical default; vulnerable to
        low-count inflation.
    ``"padj"``
        ``-log₁₀(padj)``. Ranks by statistical significance only; ignores
        effect size.

    Ties are broken by padj (ascending) then by log₂FC (descending).
    """
    if len(df) == 0:
        return df.copy()
    df = df.copy()
    eps = 1e-300  # avoid log10(0) for genes with padj==0

    if metric == "pi_score":
        score = np.abs(df["log2FoldChange"]) * -np.log10(df["padj"].clip(lower=eps))
    elif metric == "log2fc":
        score = df["log2FoldChange"].astype(float)
    elif metric == "padj":
        score = -np.log10(df["padj"].clip(lower=eps))
    else:
        raise ValueError(
            f"Unknown ranking metric '{metric}'. "
            f"Valid options: 'pi_score', 'log2fc', 'padj'."
        )

    df = df.assign(_rank_score=score)
    df = df.sort_values(
        ["_rank_score", "padj", "log2FoldChange"],
        ascending=[False, True, False],
        kind="mergesort",
    ).drop(columns="_rank_score").reset_index(drop=True)
    logger.info(
        "rank_enriched_genes: metric='%s', %d genes ranked. "
        "Top 3: %s",
        metric, len(df),
        df["_hypomap_gene_name"].head(3).tolist() if "_hypomap_gene_name" in df.columns else df.head(3).index.tolist(),
    )
    return df


def compute_marker_genes(
    adata,
    annotation_col: str,
    n_genes: int = 100,
    min_cells: int = 10,
    method: str = "wilcoxon",
) -> Dict[str, List[str]]:
    """
    Compute marker genes per cluster using scanpy's rank_genes_groups.

    The ``method`` argument is forwarded to ``sc.tl.rank_genes_groups``.
    ``"wilcoxon"`` is the default (robust, non-parametric) but scales
    poorly on very large atlases — on HypoMap (~385K cells × 185
    clusters) it takes ~15 minutes.  ``"t-test_overestim_var"`` produces
    a comparable ranking in a fraction of the time and is suitable when
    rankings — not exact p-values — are what downstream steps consume.

    Returns dict mapping cluster name to list of marker gene names.
    """
    logger.info("compute_marker_genes: annotation_col=%s, n_genes=%d, min_cells=%d, method=%s",
                annotation_col, n_genes, min_cells, method)
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
        method=method,
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
        d = universe_size - n_enriched - n_markers + n_overlap

        # When `d` (the "neither" cell) is < 5, Fisher's exact test is
        # essentially uninformative — odds_ratio explodes to infinity for
        # d=0 and the "greater" alternative is trivially satisfied. A
        # mis-set `universe_size` is the usual cause. Emit NaN so the
        # cluster drops out of downstream ranking / FDR rather than
        # silently topping the list with a sentinel odds-ratio.
        if d < 5:
            logger.warning(
                "fisher_overlap_test: cluster %s has d=%d < 5 "
                "(universe=%d, enriched=%d, markers=%d, overlap=%d); "
                "marking pvalue/odds_ratio as NaN. Verify universe_size.",
                cluster, d, universe_size, n_enriched, n_markers, n_overlap,
            )
            odds_ratio, pvalue = np.nan, np.nan
        else:
            table = np.array([[a, b], [c, d]])
            odds_ratio, pvalue = stats.fisher_exact(table, alternative="greater")

        overlap_names = sorted([enriched_original.get(g, g) for g in overlap])

        # Non-finite odds-ratio / log-OR become NaN (was 999.0 / 10.0
        # magic sentinels that looked like real values in tables).
        or_val = float(odds_ratio) if np.isfinite(odds_ratio) else np.nan
        log2_or = float(np.log2(or_val)) if np.isfinite(or_val) and or_val > 0 else np.nan
        neg_log10_p = (
            float(-np.log10(max(pvalue, 1e-300)))
            if np.isfinite(pvalue) else np.nan
        )

        results.append({
            "cluster": cluster,
            "overlap_count": n_overlap,
            "n_enriched": n_enriched,
            "n_markers": n_markers,
            "overlap_genes": ", ".join(overlap_names),
            "odds_ratio": or_val,
            "pvalue": pvalue,
            "neg_log10_pval": neg_log10_p,
            "log2_odds_ratio": log2_or,
        })

    df = pd.DataFrame(results)
    if len(df) > 0:
        # FDR correction — mask NaN pvalues (clusters that failed the
        # d>=5 guard) so they propagate as NaN padj instead of breaking
        # multipletests.
        mask = df["pvalue"].notna().values
        padj = np.full(len(df), np.nan)
        if mask.any():
            _, padj_valid, _, _ = multipletests(
                df.loc[mask, "pvalue"].values, method="fdr_bh",
            )
            padj[mask] = padj_valid
        df["padj"] = padj
        df = df.sort_values("pvalue", na_position="last").reset_index(drop=True)
        n_sig = int((df["padj"] < 0.05).sum())
        n_skipped = int((~mask).sum())
        logger.info(
            "  Fisher result: %d clusters tested (%d valid, %d skipped "
            "for d<5), %d significant (padj<0.05)",
            len(df), int(mask.sum()), n_skipped, n_sig,
        )
    return df


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
        DataFrame with columns ``cluster``, ``weight``, ``weight_norm``
        (sum-to-1 rescaling of ``weight``), and ``global_residual_norm``
        (the scalar ‖Aw − b‖₂ returned by scipy's NNLS, repeated on every
        row so the CSV / downstream pipelines see a stable schema).
        Sorted by ``weight`` descending.
    """
    # Stable output schema — referenced by app.py caching and CSV export.
    # Every early-exit path must return a DataFrame with these exact columns
    # (empty but typed) so consumers never branch on schema variations.
    _NNLS_COLUMNS = ["cluster", "weight", "weight_norm", "global_residual_norm"]

    logger.info("compute_nnls_deconvolution: %d bacTRAP genes, %d expr genes, %d clusters",
                len(bactrap_matched), len(cluster_mean_expr), len(cluster_mean_expr.columns))
    bt_dedup = bactrap_matched.drop_duplicates(subset="_hypomap_gene_name", keep="first")
    bt_genes = bt_dedup["_hypomap_gene_name"].values
    expr_genes = cluster_mean_expr.index.values
    common = np.intersect1d(bt_genes, expr_genes)
    logger.info("  common genes: %d", len(common))

    if len(common) < 5:
        logger.warning("  <5 common genes — skipping NNLS")
        return pd.DataFrame(columns=_NNLS_COLUMNS)

    bt_lookup = dict(zip(bt_dedup["_hypomap_gene_name"], bt_dedup[enrichment_col]))
    enrichment = np.array([bt_lookup[g] for g in common], dtype=float)
    valid = np.isfinite(enrichment)
    enrichment = enrichment[valid]
    common = common[valid]

    if len(enrichment) < 5:
        return pd.DataFrame(columns=_NNLS_COLUMNS)

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
        return pd.DataFrame(columns=_NNLS_COLUMNS)

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
    # Wrap in a pandas Series once so the per-cluster membership test can
    # use the vectorised C-level `isin` instead of an N-long Python list
    # comprehension. On HypoMap (~30k genes × ~185 clusters) this removes
    # a ~5.5M-iteration Python loop from the GSEA hot path.
    ranked_series = pd.Series(ranked_lower)

    # Pre-generate all permutation indices at once
    rng = np.random.default_rng(42)
    perm_indices = np.empty((n_perm, N), dtype=np.intp)
    for p in range(n_perm):
        perm_indices[p] = rng.permutation(N)

    results = []
    running_scores_dict = {}

    for cluster, markers in cluster_markers.items():
        marker_set = set(g.lower() for g in markers)
        # Boolean hit mask — vectorised membership test (pandas isin uses
        # a hash-based lookup under the hood; ~20–100× faster than the
        # previous Python-level `[g in marker_set for g in ranked_lower]`).
        hit_mask = ranked_series.isin(marker_set).to_numpy()
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

        # Sign-split NES normalization (Subramanian et al., PNAS 2005):
        # positive and negative enrichment scores live on asymmetric null
        # distributions, so each ES is normalized by the mean of null ESs
        # with the same sign.  Normalizing by mean(|null|) mixes the two
        # and can inflate |NES| when the null is dominated by one sign.
        if es >= 0:
            same_sign_null = null_es[null_es > 0]
            norm = float(same_sign_null.mean()) if same_sign_null.size else 0.0
        else:
            same_sign_null = null_es[null_es < 0]
            norm = float(-same_sign_null.mean()) if same_sign_null.size else 0.0
        # When the null has no values on the same side as the observed ES
        # (e.g. every permutation landed positive but observed ES was
        # slightly negative), NES is undefined. Previous behaviour was to
        # return NES=0, which collides with "ES exactly zero" and makes
        # unreachable-null clusters look like no-enrichment hits in the
        # NES-sorted output. Return NaN so downstream ranking / FDR
        # handles them explicitly.
        if norm > 0:
            nes = es / norm
        else:
            logger.warning(
                "compute_gsea_enrichment: cluster %s has empty same-sign "
                "null (ES=%.3f, n_perm=%d); NES set to NaN.",
                cluster, es, n_perm,
            )
            nes = np.nan

        results.append({
            "cluster": cluster,
            "ES": es,
            "NES": nes,
            "pvalue": pval,
            "n_hits": n_hits,
        })

    df = pd.DataFrame(results)
    if len(df) > 0:
        _, padj, _, _ = multipletests(df["pvalue"].values, method="fdr_bh")
        df["padj"] = padj
        # NaN NES (empty same-sign null) sorted to the end so real hits
        # appear first.
        df = df.sort_values(
            "NES", ascending=False, na_position="last",
        ).reset_index(drop=True)
        n_sig = int((df["padj"] < 0.05).sum())
        n_undef_nes = int(df["NES"].isna().sum())
        top_nes = df["NES"].iloc[0]
        top_cluster = df["cluster"].iloc[0]
        top_str = "NaN" if pd.isna(top_nes) else f"{top_nes:.2f}"
        logger.info(
            "  GSEA result: %d clusters, %d significant (padj<0.05), "
            "%d with undefined NES (empty same-sign null), top NES=%s (%s)",
            len(df), n_sig, n_undef_nes, top_str, top_cluster,
        )
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
    seed: int = 0,
    info_out: Optional[Dict] = None,
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

    Tie-breaking: in sparse single-cell data thousands of genes per cell are
    tied at low integer counts (especially 0). ``np.argpartition`` and
    ``np.argsort`` break ties by memory layout, so signature genes at low
    matrix indices would systematically win their ties (or lose, depending on
    layout) — biasing AUCell scores. We add per-cell uniform jitter smaller
    than the smallest gap between distinct expression values, which preserves
    the order of *distinct* values but randomises the order of ties. This
    matches R AUCell's ``ties.method = "random"`` (Aibar et al. 2017,
    Methods §"Building the rankings").

    Args:
        adata: AnnData object
        gene_names: list of bacTRAP-enriched gene names
        use_raw: whether to use adata.raw for expression
        top_fraction: fraction of ranked genes to consider (default 5%)
        seed: RNG seed for the per-cell tie-breaking jitter (reproducible)
        info_out: optional dict; if provided, augmented with scoring diagnostics
            (``n_query_matched``, ``unmatched``, ``n_top``, ``n_top_bumped``,
            ``effective_top_fraction``) so the caller can surface them to users.

    Returns:
        Array of AUCell scores, one per cell.
    """
    # Use the shared gene-name lookup built by data_loading — it resolves
    # both symbols and Ensembl IDs, matching the lookup used everywhere else
    # in the pipeline (fix #5: previously this function used a simpler
    # lowercase-only lookup that silently dropped genes when the matrix
    # layer used a different namespace than the signature).
    lookup, source_gene_names, is_raw_lookup = _build_adata_gene_lookup(
        adata, use_raw=use_raw,
    )
    if is_raw_lookup and adata.raw is not None:
        X = adata.raw.X
    else:
        X = adata.X

    query_idx: List[int] = []
    unmatched: List[str] = []
    for g in gene_names:
        g_lower = str(g).strip().lower()
        if g_lower in lookup:
            query_idx.append(lookup[g_lower][1])
        else:
            unmatched.append(str(g))

    n_cells = X.shape[0]
    n_total_genes = X.shape[1]
    n_query = len(query_idx)

    logger.info("compute_aucell_scores: %d/%d genes found, %d cells, top_fraction=%.2f, seed=%d",
                n_query, len(gene_names), n_cells, top_fraction, seed)
    if unmatched:
        _sample = unmatched[:20]
        logger.warning(
            "  %d/%d signature genes did not match the %s layer lookup "
            "(first 20: %s). Check gene-namespace / alias resolution between "
            "the bacTRAP DE table and the atlas.",
            len(unmatched), len(gene_names),
            "raw" if is_raw_lookup else "X",
            _sample,
        )
    if n_query == 0:
        logger.warning("  no genes found — returning zero AUCell scores")
        if info_out is not None:
            info_out.update({
                "n_query_matched": 0,
                "unmatched": unmatched,
                "n_top": 0,
                "n_top_bumped": False,
                "effective_top_fraction": 0.0,
            })
        return np.zeros(n_cells)

    # Boolean mask for query genes (vectorized membership test)
    query_mask = np.zeros(n_total_genes, dtype=bool)
    query_mask[query_idx] = True

    # Number of top genes to consider per cell. When the signature is larger
    # than top_fraction * n_total_genes, n_top is bumped to the signature
    # size — the canonical Aibar formula requires n_top >= n_query so every
    # hit can be recovered. Warn loudly in that case (fix #6): the user
    # thinks they're at "top 5%" but the window is actually wider, which
    # makes scores less stringent than the slider suggests.
    requested_n_top = int(n_total_genes * top_fraction)
    n_top = max(requested_n_top, n_query)
    n_top = min(n_top, n_total_genes)
    n_top_bumped = n_top > requested_n_top
    effective_top_fraction = float(n_top) / max(n_total_genes, 1)
    if n_top_bumped:
        logger.warning(
            "  n_top bumped from %d (%.2f%% of genes) to %d (%.2f%%) to fit "
            "the signature of %d genes. The AUCell window is wider than the "
            "requested top_fraction — effective threshold is %.2f%%.",
            requested_n_top, top_fraction * 100,
            n_top, effective_top_fraction * 100,
            n_query, effective_top_fraction * 100,
        )
    else:
        logger.info("  n_top=%d (%.2f%% of %d genes)",
                    n_top, effective_top_fraction * 100, n_total_genes)
    # Theoretical maximum of sum(cumsum(is_hit)) when all n_query query genes
    # occupy the top n_query positions:
    #   cumsum = [1, 2, ..., n_query, n_query, ..., n_query]  (n_top entries)
    #   sum    = n_query*(n_query+1)/2 + n_query*(n_top - n_query)
    #          = n_query * (n_top - (n_query - 1)/2)
    # Note: n_top >= n_query is guaranteed above, so this is always positive.
    max_auc = n_query * (n_top - (n_query - 1) / 2)

    # Two independent RNG streams so that the per-cell jitter is purely a
    # function of `seed` and not of how many cells were sampled for the
    # jitter-scale heuristic. Without the split, rng.choice(n_cells, ...)
    # consumed a variable amount of state before the jitter draws, so
    # otherwise-identical signature runs on atlases with different
    # n_cells would produce different AUCell numbers at the same seed.
    sampling_rng, jitter_rng = np.random.default_rng(seed).spawn(2)

    # Determine a safe jitter scale: must be smaller than the smallest gap
    # between distinct expression values, otherwise jitter could reorder
    # genuinely distinct values. For raw integer counts the smallest gap is 1
    # (so jitter < 0.5 is safe). For non-integer (e.g. log-normalised) data
    # we fall back on a sample-based estimate of half the smallest nonzero
    # value, which is a conservative proxy for the smallest distinct gap.
    sample_n = min(500, n_cells)
    sample_pick = (
        sampling_rng.choice(n_cells, sample_n, replace=False)
        if n_cells > sample_n else np.arange(n_cells)
    )
    sample_X = X[sample_pick, :]
    if sparse.issparse(sample_X):
        sample_X = np.asarray(sample_X.toarray())
    else:
        sample_X = np.asarray(sample_X)
    sample_max = float(sample_X.max()) if sample_X.size else 0.0
    sample_is_integer = (
        sample_X.size > 0 and bool(np.all(sample_X == np.round(sample_X)))
    )
    if sample_is_integer and sample_max > 5:
        jitter_scale = np.float32(0.49)
        logger.info("  input looks like raw integer counts; jitter_scale=0.49")
    else:
        nonzero_vals = sample_X[sample_X > 0]
        if nonzero_vals.size > 0:
            jitter_scale = np.float32(0.49 * float(nonzero_vals.min()))
        else:
            jitter_scale = np.float32(1e-6)
        logger.warning(
            "  input does not look like raw integer counts (max=%.2f, integer=%s); "
            "using scaled jitter (%.2e). AUCell is designed for raw counts "
            "(Aibar et al. 2017) — set use_raw=True against an integer-count layer "
            "for the cleanest behaviour.",
            sample_max, sample_is_integer, float(jitter_scale),
        )

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
        # Make a writable float32 copy so the in-place jitter add doesn't
        # touch the underlying adata buffer.
        X_chunk = X_chunk.astype(np.float32, copy=True)

        chunk_n = X_chunk.shape[0]

        # Per-cell uniform jitter ∈ [0, jitter_scale). Smaller than the
        # smallest distinct gap, so distinct values keep their order while
        # ties are randomised — equivalent to ties.method="random" in R.
        X_chunk += jitter_rng.random(
            (chunk_n, n_total_genes), dtype=np.float32,
        ) * jitter_scale

        # For each cell, get top-n gene indices via argpartition (O(n) per cell)
        # Then check which are query genes and compute cumulative AUC
        top_idx = np.argpartition(X_chunk, -n_top, axis=1)[:, -n_top:]

        for i in range(chunk_n):
            cell_top = top_idx[i]
            # Sort by jittered expression descending
            order = np.argsort(X_chunk[i, cell_top])[::-1]
            sorted_top = cell_top[order]

            # Vectorized: check query membership and cumsum for AUC
            is_hit = query_mask[sorted_top]
            cumhits = np.cumsum(is_hit)
            auc = cumhits.sum()
            scores[start + i] = auc / max_auc if max_auc > 0 else 0.0

    logger.info("  AUCell scores: mean=%.4f, std=%.4f, min=%.4f, max=%.4f",
                float(scores.mean()), float(scores.std()), float(scores.min()), float(scores.max()))

    if info_out is not None:
        info_out.update({
            "n_query_matched": int(n_query),
            "n_query_requested": int(len(gene_names)),
            "unmatched": list(unmatched),
            "n_top": int(n_top),
            "n_top_bumped": bool(n_top_bumped),
            "requested_top_fraction": float(top_fraction),
            "effective_top_fraction": float(effective_top_fraction),
            "source_layer": "raw" if is_raw_lookup else "X",
        })

    return scores


def _atlas_gene_mean_logexpr(adata, use_raw: bool = True) -> np.ndarray:
    """Atlas-wide mean expression per gene, log1p-scaled (raw-layer order).

    Used purely to bin genes by expression level for control-set matching, so
    a (monotone) raw-mean → log1p transform is sufficient — the bin
    assignments are rank-based and unaffected by the transform.  Cached on
    ``adata.uns['gene_mean_expr']`` so it is computed once per atlas load;
    the sparse column sums are cheap (no dense materialisation).
    """
    cache_key = "gene_mean_expr"
    if use_raw and adata.raw is not None:
        X = adata.raw.X
    else:
        X = adata.X
    n_genes_expected = X.shape[1]
    cached = adata.uns.get(cache_key)
    if cached is not None:
        cached = np.asarray(cached)
        if cached.shape[0] == n_genes_expected:
            return cached
    if sparse.issparse(X):
        gene_sum = np.asarray(X.sum(axis=0)).ravel().astype(np.float64)
    else:
        gene_sum = np.asarray(X, dtype=np.float64).sum(axis=0)
    gene_mean = np.log1p(gene_sum / max(X.shape[0], 1))
    adata.uns[cache_key] = gene_mean
    logger.info("Built bin lookup over %d atlas genes", len(gene_mean))
    return gene_mean


def _expression_bins(gene_mean_expr: np.ndarray, n_bins: int) -> np.ndarray:
    """Assign each gene to an expression quantile bin (0..n_bins-1).

    Ties (e.g. the large mass of never-expressed genes) collapse bins; the
    returned array still gives every gene a finite bin index.  Cached-friendly
    but cheap, so recomputed per call to follow the user's ``n_bins``.
    """
    s = pd.Series(np.asarray(gene_mean_expr, dtype=np.float64))
    try:
        bins = pd.qcut(s.rank(method="first"), q=int(max(n_bins, 1)),
                       labels=False, duplicates="drop")
    except ValueError:
        bins = pd.Series(np.zeros(len(s), dtype=int))
    return bins.fillna(0).astype(int).to_numpy()


def compute_empirical_null_aucell(
    adata,
    signature_genes: List[str],
    cluster_labels,
    compute_aucell_fn,
    *,
    n_control_sets: int = 100,
    n_bins: int = 5,
    seed: int = 0,
    top_fraction: float = 0.05,
    min_cluster_size: int = 20,
    use_raw: bool = True,
    logger=None,
    progress_callback=None,
) -> pd.DataFrame:
    """Empirical-null AUCell with expression-matched control gene sets.

    For each cluster, compares the bacTRAP-signature AUCell mean against the
    distribution of AUCell means obtained from ``n_control_sets`` random
    control gene sets matched to the signature in size and atlas-wide
    expression structure (Aibar et al., *Nat. Methods* 2017).  This removes
    the "baseline gene-rank-width" bias that inflates AUCell scores in
    broadly-active neuronal clusters.

    Parameters
    ----------
    adata
        Atlas the AUCell scoring runs against (neuronal subset when the
        neuronal-only mask is active).
    signature_genes
        The bacTRAP signature gene symbols (atlas namespace).
    cluster_labels
        Per-cell cluster labels aligned positionally with ``adata.obs_names``.
    compute_aucell_fn
        Callable with the ``compute_aucell_scores`` signature — injected so
        the function is testable and so the same scorer (and tie-breaking
        regime) is used for the signature and the controls.
    n_control_sets, n_bins, seed, top_fraction, min_cluster_size
        See module / sidebar docs.
    progress_callback
        Optional ``fn(i, n)`` called after scoring control set ``i`` of ``n``.

    Returns
    -------
    DataFrame indexed by cluster with columns
        ``null_mean``, ``null_sd``, ``z_empirical``, ``pvalue_empirical``,
        ``qvalue_empirical``, ``n_control_sets_used``.
    """
    log = logger or globals()["logger"]
    rng = np.random.default_rng(int(seed))

    # ---- Per-cell signature AUCell (scored once, here, so signature and
    # controls share the exact same scorer / jitter regime) ----
    sig_scores = np.asarray(
        compute_aucell_fn(adata, list(signature_genes),
                          top_fraction=top_fraction, seed=int(seed)),
        dtype=np.float64,
    )
    labels = pd.Series(np.asarray(cluster_labels).astype(str))
    if len(labels) != len(sig_scores):
        raise ValueError(
            f"cluster_labels ({len(labels)}) and AUCell scores ({len(sig_scores)}) "
            f"length mismatch"
        )
    cluster_sizes = labels.value_counts()
    eligible_clusters = cluster_sizes[cluster_sizes >= min_cluster_size].index.tolist()
    if not eligible_clusters:
        log.warning("compute_empirical_null_aucell: no cluster >= %d cells", min_cluster_size)
        return pd.DataFrame(columns=[
            "null_mean", "null_sd", "z_empirical", "pvalue_empirical",
            "qvalue_empirical", "n_control_sets_used",
        ])

    def _per_cluster_means(scores: np.ndarray) -> pd.Series:
        return pd.Series(scores).groupby(labels.values).mean().reindex(eligible_clusters)

    sig_means = _per_cluster_means(sig_scores)

    # ---- Expression-matched control gene sets ----
    lookup, gene_names, is_raw = _build_adata_gene_lookup(adata, use_raw=use_raw)
    gene_mean_expr = _atlas_gene_mean_logexpr(adata, use_raw=is_raw)
    gene_bins = _expression_bins(gene_mean_expr, n_bins)
    adata.uns["gene_expr_bins"] = gene_bins

    sig_idx: List[int] = []
    matched_genes: List[str] = []
    for g in signature_genes:
        key = str(g).strip().lower()
        if key in lookup:
            sig_idx.append(lookup[key][1])
            matched_genes.append(str(g))
    sig_idx = list(dict.fromkeys(sig_idx))  # de-dup, preserve order
    sig_idx_set = set(sig_idx)

    # Bin → candidate indices (excluding signature genes themselves)
    bin_to_candidates: Dict[int, np.ndarray] = {}
    for b in np.unique(gene_bins):
        cand = np.where(gene_bins == b)[0]
        cand = cand[~np.isin(cand, list(sig_idx_set))]
        bin_to_candidates[int(b)] = cand

    matchable_idx: List[int] = []
    dropped_idx: List[int] = []
    for idx in sig_idx:
        b = int(gene_bins[idx])
        if bin_to_candidates.get(b) is not None and bin_to_candidates[b].size > 0:
            matchable_idx.append(idx)
        else:
            dropped_idx.append(idx)
    k_total = len(sig_idx)
    k_matched = len(matchable_idx)
    k_dropped = len(dropped_idx)
    log.info(
        "Empirical null: enabled (N=%d, bins=%d, seed=%d)",
        int(n_control_sets), int(n_bins), int(seed),
    )
    log.info(
        "Signature size: %d, in-bin matchable: %d, dropped: %d",
        k_total, k_matched, k_dropped,
    )
    if k_dropped:
        log.warning(
            "  %d signature gene(s) have no in-bin control candidates and were "
            "excluded from control matching (signature scoring still uses all "
            "matched genes): %s",
            k_dropped, [gene_names[i] for i in dropped_idx][:20],
        )
    if k_matched == 0:
        log.warning("compute_empirical_null_aucell: no matchable signature genes — skipping")
        return pd.DataFrame(columns=[
            "null_mean", "null_sd", "z_empirical", "pvalue_empirical",
            "qvalue_empirical", "n_control_sets_used",
        ])

    import time as _time
    control_cluster_means = np.empty((int(n_control_sets), len(eligible_clusters)), dtype=np.float64)
    n_used = 0
    _t0 = _time.time()
    _per_set_times: List[float] = []
    for c in range(int(n_control_sets)):
        _ts = _time.time()
        ctrl_idx = []
        for idx in matchable_idx:
            cand = bin_to_candidates[int(gene_bins[idx])]
            ctrl_idx.append(int(rng.choice(cand)))
        ctrl_names = [str(gene_names[i]) for i in ctrl_idx]
        ctrl_scores = np.asarray(
            compute_aucell_fn(adata, ctrl_names, top_fraction=top_fraction, seed=int(seed)),
            dtype=np.float64,
        )
        control_cluster_means[c, :] = _per_cluster_means(ctrl_scores).to_numpy()
        n_used += 1
        _per_set_times.append(_time.time() - _ts)
        if progress_callback is not None:
            try:
                progress_callback(c + 1, int(n_control_sets))
            except Exception:
                pass
    _t_total = _time.time() - _t0
    log.info(
        "Mean control set scoring time: %.1fs; total null time: %.1fs",
        float(np.mean(_per_set_times)) if _per_set_times else 0.0, _t_total,
    )

    null_mean = control_cluster_means.mean(axis=0)
    null_sd = control_cluster_means.std(axis=0, ddof=1) if n_used > 1 else np.zeros(len(eligible_clusters))
    sig_vec = sig_means.to_numpy(dtype=np.float64)
    degenerate = ~np.isfinite(null_sd) | (null_sd == 0)
    z = np.where(degenerate, np.nan, (sig_vec - null_mean) / np.where(degenerate, 1.0, null_sd))
    # one-sided empirical p (add-1 smoothing): controls >= signature
    ge_counts = (control_cluster_means >= sig_vec[None, :]).sum(axis=0)
    pvals = (1.0 + ge_counts) / (n_used + 1.0)

    if degenerate.any():
        log.warning(
            "Clusters with degenerate null_sd (set to NaN z): %s",
            [eligible_clusters[i] for i in np.where(degenerate)[0]],
        )

    out = pd.DataFrame({
        "null_mean": null_mean,
        "null_sd": null_sd,
        "z_empirical": z,
        "pvalue_empirical": pvals,
        "n_control_sets_used": n_used,
    }, index=pd.Index(eligible_clusters, name="cluster"))
    valid_p = out["pvalue_empirical"].notna().to_numpy()
    qvals = np.full(len(out), np.nan)
    if valid_p.any():
        _, q, _, _ = multipletests(out.loc[valid_p, "pvalue_empirical"].values, method="fdr_bh")
        qvals[valid_p] = q
    out["qvalue_empirical"] = qvals
    out = out[[
        "null_mean", "null_sd", "z_empirical", "pvalue_empirical",
        "qvalue_empirical", "n_control_sets_used",
    ]]

    # Sanity check: rank correlation between signature mean and z (should be
    # positive but < ~0.95 — the null must re-order at least some clusters).
    try:
        valid = np.isfinite(z) & np.isfinite(sig_vec)
        if valid.sum() >= 3:
            rho = stats.spearmanr(sig_vec[valid], z[valid]).statistic
            log.info("Empirical null sanity: Spearman(mean, z_empirical) = %.3f over %d clusters",
                     float(rho), int(valid.sum()))
            ranked_by_mean = sig_means[valid].sort_values(ascending=False).index.tolist()
            z_series = pd.Series(z, index=eligible_clusters)
            ranked_by_z = z_series[valid].sort_values(ascending=False).index.tolist()
            for cl in ranked_by_mean[:10]:
                drop = ranked_by_z.index(cl) - ranked_by_mean.index(cl)
                if drop > 0:
                    log.info("  %s: rank by mean=%d -> by z=%d (drops %d)",
                             cl, ranked_by_mean.index(cl) + 1, ranked_by_z.index(cl) + 1, drop)
    except Exception:
        log.debug("Empirical null sanity check failed", exc_info=True)

    return out


def validate_aucell_input(adata, use_raw: bool = True, sample_size: int = 500) -> Dict:
    """Sanity-check the expression layer that will be fed into AUCell.

    AUCell is defined on gene-expression *ranks* per cell; the paper (Aibar
    et al. 2017) and the reference implementations build rankings from raw
    UMI counts. Feeding in pre-log-normalised data is not strictly wrong —
    the rank order of *distinct* values is preserved by any monotonic
    transform — but ties collapse very differently (most log-normalised
    matrices lose the fine-grained tie structure of low counts), so scores
    and their interpretation shift. This function inspects a random sample
    of the chosen layer and returns a machine-readable QC report that the
    caller can surface to the user.

    Returns:
        dict with keys
            layer: "raw" or "X"
            n_cells, n_genes, n_cells_sampled
            max_value, min_nonzero
            is_integer: bool — all sampled values are integers
            looks_like_counts: bool — is_integer AND max_value > 50
            warnings: list[str] — user-facing strings (empty when clean)
            info: list[str] — informational lines
    """
    report: Dict = {
        "layer": "raw" if (use_raw and adata.raw is not None) else "X",
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_cells_sampled": 0,
        "max_value": 0.0,
        "min_nonzero": float("nan"),
        "is_integer": False,
        "looks_like_counts": False,
        "warnings": [],
        "info": [],
    }
    if use_raw and adata.raw is not None:
        X = adata.raw.X
        report["n_genes"] = int(adata.raw.n_vars)
    else:
        X = adata.X
        if use_raw and adata.raw is None:
            report["warnings"].append(
                "use_raw=True was requested but adata.raw is None; falling "
                "back to adata.X. If adata.X is log-normalised, AUCell tie "
                "structure will differ from the canonical raw-count version."
            )

    n_cells = X.shape[0]
    sample_n = min(sample_size, n_cells)
    if n_cells > sample_n:
        rng = np.random.default_rng(42)
        pick = rng.choice(n_cells, sample_n, replace=False)
        sample = X[pick, :]
    else:
        sample = X
    if sparse.issparse(sample):
        sample = np.asarray(sample.toarray())
    else:
        sample = np.asarray(sample)
    report["n_cells_sampled"] = int(sample.shape[0])

    if sample.size == 0:
        report["warnings"].append("Selected layer is empty.")
        return report

    sample_max = float(sample.max())
    nonzero = sample[sample > 0]
    sample_min_nz = float(nonzero.min()) if nonzero.size > 0 else float("nan")
    is_int = bool(np.all(sample == np.round(sample)))
    looks_like_counts = is_int and sample_max > 50

    report["max_value"] = sample_max
    report["min_nonzero"] = sample_min_nz
    report["is_integer"] = is_int
    report["looks_like_counts"] = looks_like_counts
    report["info"].append(
        f"AUCell input layer: {report['layer']} ({sample.shape[0]} cells x "
        f"{sample.shape[1]} genes sampled, max={sample_max:.2f}, "
        f"min_nonzero={sample_min_nz:.3g}, integer={is_int})"
    )

    if not looks_like_counts:
        if is_int and sample_max <= 50:
            report["warnings"].append(
                f"AUCell input has integer values but max = {sample_max:.0f} "
                f"(<= 50) — unusually low for raw UMI counts. Verify that "
                f"adata.raw contains counts and not e.g. a binarised layer."
            )
        else:
            report["warnings"].append(
                f"AUCell input does NOT look like raw integer counts "
                f"(max={sample_max:.2f}, integer={is_int}). AUCell "
                f"(Aibar et al. 2017) is defined on raw counts; log-"
                f"normalised or scaled input will preserve distinct-value "
                f"rank order but collapses tie structure differently. "
                f"Consider pointing adata.raw at a raw-count layer before "
                f"scoring — compute_aucell_scores will still run and will "
                f"scale its tie-breaking jitter accordingly."
            )

    logger.info("validate_aucell_input: layer=%s, looks_like_counts=%s, "
                "max=%.2f, integer=%s, n_warnings=%d",
                report["layer"], looks_like_counts, sample_max, is_int,
                len(report["warnings"]))
    return report


def compute_cluster_enrichment_stats(
    aucell_scores: np.ndarray,
    cell_labels: np.ndarray,
    min_cells: int = 10,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Per-cluster enrichment significance from AUCell scores.

    For each cluster we test the one-sided null that its AUCell scores are
    drawn from the same distribution as the rest of the atlas, using
    Welch's t-test (unequal variance). Multiple testing across clusters is
    corrected with Benjamini-Hochberg to give a per-cluster q-value.

    This fills the gap the previous pipeline had: ranking clusters by mean
    AUCell alone cannot distinguish "strongly enriched" from "a small
    cluster that happens to have a slightly above-average mean" — the
    q-value addresses exactly that. Clusters with fewer than *min_cells*
    cells are dropped (Welch's t breaks down at very small n).

    Returns a DataFrame (sorted by qvalue, ascending) with columns:
        cluster, n_cells, mean, sem, std, t_stat, pvalue, qvalue, significant
    """
    scores = np.asarray(aucell_scores, dtype=np.float64)
    labels = np.asarray(cell_labels)
    if scores.shape[0] != labels.shape[0]:
        raise ValueError(
            f"aucell_scores ({scores.shape[0]}) and cell_labels "
            f"({labels.shape[0]}) length mismatch"
        )

    unique_labels = pd.unique(labels)
    rows = []
    for cl in unique_labels:
        mask = labels == cl
        n = int(mask.sum())
        if n < min_cells:
            continue
        in_scores = scores[mask]
        out_scores = scores[~mask]
        if in_scores.size < 2 or out_scores.size < 2:
            continue
        # Welch's one-sided t-test: cluster > rest
        t_stat, p_two = stats.ttest_ind(
            in_scores, out_scores, equal_var=False, nan_policy="omit",
        )
        # scipy's ttest_ind returns two-sided; convert to one-sided upper tail
        if np.isnan(t_stat):
            continue
        p_one = (p_two / 2.0) if t_stat > 0 else (1.0 - p_two / 2.0)
        rows.append({
            "cluster": str(cl),
            "n_cells": n,
            "mean": float(in_scores.mean()),
            "std": float(in_scores.std(ddof=1)),
            "sem": float(in_scores.std(ddof=1) / np.sqrt(n)),
            "t_stat": float(t_stat),
            "pvalue": float(p_one),
        })

    df = pd.DataFrame(rows)
    if df.empty:
        logger.warning("compute_cluster_enrichment_stats: no clusters passed "
                       "the min_cells=%d filter", min_cells)
        df["qvalue"] = []
        df["significant"] = []
        return df

    _, qvals, _, _ = multipletests(df["pvalue"].values, method="fdr_bh")
    df["qvalue"] = qvals
    df["significant"] = df["qvalue"] < alpha

    df = df.sort_values(["qvalue", "pvalue"], ascending=True).reset_index(drop=True)

    n_sig = int(df["significant"].sum())
    logger.info("compute_cluster_enrichment_stats: %d/%d clusters tested, "
                "%d significant at BH-FDR q < %.3f",
                len(df), len(unique_labels), n_sig, alpha)
    return df


# =========================================================================
# Composite Ranking
# =========================================================================

def compute_composite_ranking(
    corr_df: pd.DataFrame,
    fisher_df: pd.DataFrame,
    nnls_df: pd.DataFrame,
    gsea_df: Optional[pd.DataFrame] = None,
    allowed_clusters: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    """
    Combine multiple ranking methods into a single consensus ranking.

    For each method, ranks are converted to percentile scores (0–1),
    then averaged. This produces a robust ranking that doesn't depend
    on any single method's assumptions.

    If ``allowed_clusters`` is provided, each input table is restricted
    to that set BEFORE per-method ranks / percentiles are computed, so
    the survivors are re-ranked against each other rather than inheriting
    their global positions. Used to drop clusters that fail a Cre-driver
    baseline-expression check before the composite vote.

    Returns:
        DataFrame with columns: cluster, corr_rank, fisher_rank, nnls_rank,
        gsea_rank (if available), composite_score, sorted by composite_score.
    """
    if allowed_clusters is not None:
        allowed_set = {str(c) for c in allowed_clusters}

        def _restrict(df):
            if df is None or len(df) == 0 or "cluster" not in df.columns:
                return df
            mask = df["cluster"].astype(str).isin(allowed_set)
            return df.loc[mask].reset_index(drop=True)

        corr_df = _restrict(corr_df)
        fisher_df = _restrict(fisher_df)
        nnls_df = _restrict(nnls_df)
        gsea_df = _restrict(gsea_df)

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
