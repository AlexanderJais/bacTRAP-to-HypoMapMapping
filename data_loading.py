"""
Data loading and gene matching utilities for bacTRAP-to-HypoMap mapping.
"""

import logging
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import streamlit as st
from scipy import sparse
from typing import Tuple, Optional, Dict, List

logger = logging.getLogger(__name__)


@st.cache_resource(show_spinner="Loading HypoMap atlas (this may take a few minutes)...")
def load_hypomap(file_path: str) -> ad.AnnData:
    """Load the HypoMap h5ad atlas with memory-efficient settings."""
    logger.info("Loading HypoMap from %s", file_path)
    adata = sc.read_h5ad(file_path, backed=None)
    logger.info("HypoMap loaded: %d cells x %d vars", adata.n_obs, adata.n_vars)
    logger.info("  var_names (first 10): %s", list(adata.var_names[:10]))
    logger.info("  var.columns: %s", list(adata.var.columns))
    if adata.raw is not None:
        logger.info("  raw present: %d vars", adata.raw.n_vars)
        logger.info("  raw.var_names (first 10): %s", list(adata.raw.var_names[:10]))
        logger.info("  raw.var.columns: %s", list(adata.raw.var.columns))
        # Log sample values from each raw.var column
        for col in adata.raw.var.columns:
            sample = adata.raw.var[col].dropna().head(5).tolist()
            logger.info("  raw.var['%s'] sample: %s", col, sample)
    else:
        logger.info("  raw: None")
    logger.info("  obs.columns: %s", list(adata.obs.columns))
    # Ensure expression matrix is sparse
    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X)
    if adata.raw is not None and not sparse.issparse(adata.raw.X):
        # Rebuild raw layer with sparse matrix via the public API
        import anndata
        raw_adata = anndata.AnnData(
            X=sparse.csr_matrix(adata.raw.X),
            var=adata.raw.var,
        )
        adata.raw = raw_adata
    return adata


@st.cache_data(show_spinner="Loading bacTRAP data...")
def load_bactrap(file_path: str) -> pd.DataFrame:
    """Load the bacTRAP FPKM/DESeq2 results from an Excel file."""
    logger.info("Loading bacTRAP from %s", file_path)
    df = pd.read_excel(file_path, engine="openpyxl")
    logger.info("bacTRAP loaded: %d rows x %d columns", len(df), len(df.columns))
    logger.info("  columns: %s", list(df.columns))
    logger.info("  dtypes:\n%s", df.dtypes.to_string())
    logger.info("  index dtype: %s, sample: %s", df.index.dtype, list(df.index[:5]))
    for col in df.columns:
        if df[col].dtype == object or df[col].dtype.name in ("string", "category"):
            sample = df[col].dropna().head(5).tolist()
            logger.info("  column '%s' sample: %s", col, sample)
    return df


def get_annotation_columns(adata: ad.AnnData) -> List[str]:
    """Return candidate annotation columns from .obs (categorical or string dtypes)."""
    candidates = []
    for col in adata.obs.columns:
        dtype = adata.obs[col].dtype
        if dtype.name == "category" or dtype == object:
            nunique = adata.obs[col].nunique()
            # Likely an annotation if it has a reasonable number of unique values
            if 2 <= nunique <= 5000:
                candidates.append(col)
    return sorted(candidates)


def _looks_like_ensembl(gene_names: np.ndarray, sample_size: int = 20) -> bool:
    """Check whether gene names look like Ensembl IDs by sampling."""
    if len(gene_names) == 0:
        return False
    sample = gene_names[:sample_size]
    n_ens = sum(1 for g in sample if str(g).startswith(("ENSMUSG", "ENSG")))
    return n_ens > len(sample) * 0.5


def _find_symbol_column(var_df: pd.DataFrame) -> Optional[str]:
    """Find a gene-symbol column in a var DataFrame."""
    checked = []
    for col in ["gene_name", "gene_symbol", "symbol", "Gene", "gene_short_name",
                "external_gene_name", "mgi_symbol", "feature_name"]:
        if col in var_df.columns:
            # Verify this column actually has non-Ensembl values
            sample = var_df[col].dropna().head(20)
            if len(sample) > 0 and not _looks_like_ensembl(sample.values):
                logger.info("_find_symbol_column: found '%s', sample: %s", col, sample.tolist()[:5])
                return col
            checked.append((col, sample.tolist()[:5]))
    logger.info("_find_symbol_column: no symbol column found. Checked: %s. Available: %s",
                checked, list(var_df.columns))
    return None


def get_gene_names_from_adata(adata: ad.AnnData, use_raw: bool = False) -> np.ndarray:
    """Extract gene names from the AnnData object, trying multiple locations.

    When *use_raw* is True and ``adata.raw`` exists, gene names are taken from
    the raw layer (which typically contains the full pre-HVG-filtering gene set).

    If var_names look like Ensembl IDs, a gene-symbol column is searched in
    both the target var DataFrame AND ``adata.var`` as fallback (since
    ``adata.raw.var`` may lack annotation columns).
    """
    if use_raw and adata.raw is not None:
        gene_names = np.array(adata.raw.var_names)
        var_df = adata.raw.var
        logger.info("get_gene_names_from_adata: using raw layer (%d genes)", len(gene_names))
    else:
        gene_names = adata.var_names.values.copy()
        var_df = adata.var
        logger.info("get_gene_names_from_adata: using var layer (%d genes)", len(gene_names))

    logger.info("  var_names sample (first 10): %s", [str(g) for g in gene_names[:10]])
    is_ensembl = _looks_like_ensembl(gene_names)
    logger.info("  looks_like_ensembl: %s", is_ensembl)

    # Check if var_names look like Ensembl IDs; if so, look for a symbol column
    if is_ensembl:
        logger.info("  searching for symbol column in target var_df (columns: %s)", list(var_df.columns))
        sym_col = _find_symbol_column(var_df)
        if sym_col is not None:
            gene_names = var_df[sym_col].values.copy()
            logger.info("  resolved via column '%s', sample: %s", sym_col, [str(g) for g in gene_names[:10]])
        elif use_raw and adata.raw is not None:
            # Fallback: adata.var may have the symbol column even if raw.var doesn't.
            logger.info("  no symbol column in raw.var, trying adata.var fallback (columns: %s)",
                        list(adata.var.columns))
            sym_col_main = _find_symbol_column(adata.var)
            if sym_col_main is not None:
                ens_to_sym = {}
                for ens_id, sym in zip(adata.var.index, adata.var[sym_col_main]):
                    sym_str = str(sym).strip()
                    if sym_str and sym_str.lower() != "nan":
                        ens_to_sym[str(ens_id)] = sym_str
                logger.info("  built Ensembl→symbol map from adata.var['%s']: %d entries", sym_col_main, len(ens_to_sym))
                if ens_to_sym:
                    gene_names = np.array([
                        ens_to_sym.get(str(g), str(g)) for g in gene_names
                    ])
                    n_resolved = sum(1 for g in gene_names if not str(g).startswith(("ENSMUSG", "ENSG")))
                    logger.info("  after map: %d/%d resolved to symbols", n_resolved, len(gene_names))
                    logger.info("  sample after resolution: %s", [str(g) for g in gene_names[:10]])
            else:
                logger.warning("  NO symbol column found anywhere — gene names remain as Ensembl IDs!")
    else:
        logger.info("  var_names are not Ensembl IDs, using as-is")

    return gene_names


def _detect_gene_column(bactrap_df: pd.DataFrame) -> str:
    """Auto-detect the column in bacTRAP data that contains gene identifiers.

    Checks common column names and the DataFrame index.  Returns the best
    candidate column name (or ``"_index"`` if the index should be used).
    """
    # Priority-ordered list of likely gene-name column names
    candidates = [
        "gene_name", "gene_symbol", "symbol", "Gene", "GeneSymbol",
        "gene_id", "GeneID", "external_gene_name", "mgi_symbol",
        "SYMBOL", "gene_short_name", "feature_name", "Name",
    ]
    for col in candidates:
        if col in bactrap_df.columns:
            return col

    # Check if first column looks like gene names (common in DESeq2 output)
    first_col = bactrap_df.columns[0]
    sample_vals = bactrap_df[first_col].dropna().head(20).astype(str)
    if len(sample_vals) > 0:
        # If most values are non-numeric strings, likely gene names
        n_alpha = sum(1 for v in sample_vals if v and not v.replace(".", "").replace("-", "").replace("_", "").isdigit())
        if n_alpha > len(sample_vals) * 0.5:
            return first_col

    # Check the index
    if bactrap_df.index.dtype == object or bactrap_df.index.dtype.name == "string":
        sample_idx = [str(x) for x in bactrap_df.index[:20]]
        n_alpha = sum(1 for v in sample_idx if v and not v.replace(".", "").replace("-", "").replace("_", "").isdigit())
        if n_alpha > len(sample_idx) * 0.5:
            return "_index"

    # Last resort: return first column
    return first_col


@st.cache_resource(show_spinner=False)
def _build_adata_gene_lookup(
    adata: ad.AnnData, use_raw: bool = True,
) -> Tuple[Dict[str, Tuple[str, int]], np.ndarray, bool]:
    """Build a comprehensive gene lookup from an AnnData object.

    Returns:
        lookup: lowercase gene name -> (display_name, column_index)
        gene_names: array of resolved gene names
        is_raw: whether indices point into adata.raw.var
    """
    has_raw = adata.raw is not None and use_raw
    logger.info("_build_adata_gene_lookup: has_raw=%s", has_raw)
    gene_names = get_gene_names_from_adata(adata, use_raw=has_raw)

    lookup: Dict[str, Tuple[str, int]] = {}
    for idx, name in enumerate(gene_names):
        name_str = str(name).strip()
        key = name_str.lower()
        if key and key != "nan":
            lookup[key] = (name_str, idx)

    logger.info("  primary lookup size: %d (from resolved gene names)", len(lookup))

    # Also add Ensembl IDs as lookup keys (pointing to the same indices)
    # so that bacTRAP data with Ensembl IDs can still match.
    if has_raw and adata.raw is not None:
        raw_var_names = np.array(adata.raw.var_names)
    else:
        raw_var_names = adata.var_names.values

    raw_is_ensembl = _looks_like_ensembl(raw_var_names)
    names_are_ensembl = _looks_like_ensembl(gene_names)
    logger.info("  raw_var_names looks_like_ensembl: %s", raw_is_ensembl)
    logger.info("  resolved gene_names looks_like_ensembl: %s", names_are_ensembl)

    if raw_is_ensembl:
        n_added = 0
        for idx, ens_id in enumerate(raw_var_names):
            ens_str = str(ens_id).strip().lower()
            if ens_str and ens_str != "nan" and ens_str not in lookup:
                display = str(gene_names[idx]).strip()
                if display and display.lower() != "nan":
                    lookup[ens_str] = (display, idx)
                    n_added += 1
        logger.info("  added %d Ensembl ID keys to lookup", n_added)
    elif not names_are_ensembl:
        var_df = adata.raw.var if has_raw and adata.raw is not None else adata.var
        for col in ["gene_ids", "gene_id", "ensembl_id", "Ensembl", "ensembl"]:
            if col in var_df.columns:
                n_added = 0
                for idx, ens_id in enumerate(var_df[col]):
                    ens_str = str(ens_id).strip().lower()
                    if ens_str and ens_str != "nan" and ens_str not in lookup:
                        display = str(gene_names[idx]).strip()
                        if display and display.lower() != "nan":
                            lookup[ens_str] = (display, idx)
                            n_added += 1
                logger.info("  added %d Ensembl keys from var_df['%s']", n_added, col)
                break

    logger.info("  final lookup size: %d", len(lookup))
    # Log a few sample entries from the lookup
    sample_keys = list(lookup.keys())[:10]
    logger.info("  lookup sample: %s", {k: lookup[k][0] for k in sample_keys})

    return lookup, gene_names, has_raw


def match_genes(
    bactrap_df: pd.DataFrame,
    adata: ad.AnnData,
    gene_col: Optional[str] = None,
    _prebuilt_lookup: Optional[Tuple[Dict[str, Tuple[str, int]], np.ndarray, bool]] = None,
) -> Tuple[pd.DataFrame, List[str], Dict[str, int], bool]:
    """
    Match bacTRAP gene identifiers to HypoMap genes.

    When *gene_col* is None the column is auto-detected by trying common
    column names and selecting the one that yields the most matches.

    When ``adata.raw`` exists the lookup is built from the **raw** layer so
    that *all* genes are available for matching (not just the highly-variable
    subset stored in ``adata.var``).

    Parameters:
        _prebuilt_lookup: Optional pre-built lookup tuple from
            ``_build_adata_gene_lookup`` to avoid redundant recomputation.

    Returns:
        bactrap_matched: subset of bactrap_df with matched genes
        matched_gene_names: list of matched gene symbols (as they appear in HypoMap)
        gene_to_adata_idx: mapping from gene name to column index.
            Indices are in **raw** space when ``adata.raw`` exists,
            otherwise in ``adata.var`` space.
        matched_in_raw: True when indices refer to adata.raw.var space.
    """
    logger.info("match_genes called with gene_col=%s", gene_col)

    # Use pre-built lookup if provided, otherwise build it (cached)
    if _prebuilt_lookup is not None:
        adata_gene_lookup, adata_gene_names, has_raw = _prebuilt_lookup
    else:
        adata_gene_lookup, adata_gene_names, has_raw = _build_adata_gene_lookup(adata)

    # Determine which bacTRAP column to use for gene matching
    if gene_col is None:
        gene_col = _auto_select_gene_col(bactrap_df, adata_gene_lookup)

    logger.info("match_genes: using gene_col='%s'", gene_col)

    if gene_col == "_index":
        bt_gene_values = bactrap_df.index.astype(str)
    elif gene_col in bactrap_df.columns:
        bt_gene_values = bactrap_df[gene_col].astype(str)
    else:
        raise ValueError(
            f"Column '{gene_col}' not found in bacTRAP data. "
            f"Available columns: {list(bactrap_df.columns)}"
        )

    bt_sample = [str(v) for v in bt_gene_values[:10]]
    logger.info("  bacTRAP gene values sample (first 10): %s", bt_sample)

    # Match bacTRAP genes
    matched_rows = []
    matched_gene_names = []
    gene_to_adata_idx = {}
    unmatched_sample = []

    for i, bt_gene_raw in enumerate(bt_gene_values):
        bt_gene = bt_gene_raw.strip()
        bt_gene_lower = bt_gene.lower()

        if bt_gene_lower in adata_gene_lookup:
            original_name, adata_idx = adata_gene_lookup[bt_gene_lower]
            matched_rows.append(bactrap_df.iloc[i])
            matched_gene_names.append(original_name)
            gene_to_adata_idx[original_name] = adata_idx
        elif len(unmatched_sample) < 20:
            unmatched_sample.append(bt_gene)

    logger.info("match_genes result: %d/%d matched (%.1f%%)",
                len(matched_rows), len(bt_gene_values),
                100 * len(matched_rows) / max(len(bt_gene_values), 1))
    if matched_gene_names:
        logger.info("  matched sample: %s", matched_gene_names[:10])
    logger.info("  unmatched sample: %s", unmatched_sample)

    if len(matched_rows) == 0:
        empty_df = bactrap_df.iloc[:0].copy()
        empty_df["_hypomap_gene_name"] = pd.Series(dtype=str)
        return empty_df, [], {}, has_raw

    bactrap_matched = pd.DataFrame(matched_rows)
    bactrap_matched = bactrap_matched.reset_index(drop=True)
    bactrap_matched["_hypomap_gene_name"] = matched_gene_names

    # Deduplicate: keep first occurrence when multiple bacTRAP rows map to
    # the same HypoMap gene (e.g. multiple Ensembl IDs → same symbol).
    n_before = len(bactrap_matched)
    bactrap_matched = bactrap_matched.drop_duplicates(subset="_hypomap_gene_name", keep="first")
    bactrap_matched = bactrap_matched.reset_index(drop=True)
    matched_gene_names = bactrap_matched["_hypomap_gene_name"].tolist()
    gene_to_adata_idx = {g: gene_to_adata_idx[g] for g in matched_gene_names}
    if n_before != len(bactrap_matched):
        logger.info("  deduplicated %d → %d genes (removed %d duplicate HypoMap mappings)",
                     n_before, len(bactrap_matched), n_before - len(bactrap_matched))

    return bactrap_matched, matched_gene_names, gene_to_adata_idx, has_raw


def _auto_select_gene_col(
    bactrap_df: pd.DataFrame,
    adata_gene_lookup: Dict[str, Tuple[str, int]],
) -> str:
    """Try multiple candidate columns and pick the one with the most matches."""
    candidates = []

    # Try known column names
    for col in ["gene_name", "gene_symbol", "symbol", "Gene", "GeneSymbol",
                "gene_id", "GeneID", "external_gene_name", "mgi_symbol",
                "SYMBOL", "gene_short_name", "feature_name", "Name"]:
        if col in bactrap_df.columns:
            candidates.append(col)

    # Also try the first column and the index
    first_col = bactrap_df.columns[0]
    if first_col not in candidates:
        candidates.append(first_col)
    candidates.append("_index")

    best_col = candidates[0] if candidates else "_index"
    best_count = 0

    logger.info("_auto_select_gene_col: testing %d candidates: %s", len(candidates), candidates)

    for col in candidates:
        if col == "_index":
            values = bactrap_df.index.astype(str)
        else:
            values = bactrap_df[col].astype(str)

        n_matched = sum(
            1 for v in values
            if str(v).strip().lower() in adata_gene_lookup
        )
        sample = [str(v) for v in values[:5]]
        logger.info("  col='%s': %d matches, sample=%s", col, n_matched, sample)
        if n_matched > best_count:
            best_count = n_matched
            best_col = col

    logger.info("_auto_select_gene_col: selected '%s' with %d matches", best_col, best_count)
    return best_col


def _resolve_gene_names(
    adata: ad.AnnData,
    gene_indices: np.ndarray,
    from_raw: bool = False,
) -> List[str]:
    """Resolve gene indices to display names.

    When *from_raw* is True, indices are looked up in ``adata.raw.var``.
    """
    adata_gene_names = get_gene_names_from_adata(adata, use_raw=from_raw)
    return [str(adata_gene_names[i]) for i in gene_indices]


def _map_var_indices_to_raw(
    adata: ad.AnnData,
    gene_indices: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Map gene indices from adata.var space to adata.raw.var space.

    When adata.raw has more genes than adata (common after gene filtering),
    indices into adata.var do NOT correspond to the same columns in
    adata.raw.X. This function translates them via gene name lookup.

    Returns:
        raw_indices: array of indices into adata.raw.var
        survived_mask: boolean mask over gene_indices indicating which
            original indices were successfully mapped (for correct label alignment)
    """
    if adata.raw is None:
        return gene_indices, np.ones(len(gene_indices), dtype=bool)

    # Get gene names for the requested indices in adata.var
    var_gene_names = get_gene_names_from_adata(adata)
    query_names = [str(var_gene_names[i]).lower() for i in gene_indices]

    # Build lookup for raw var names — prefer gene symbols over Ensembl IDs
    # so that the lookup matches how match_genes() found these genes.
    raw_var_names = adata.raw.var_names
    raw_lookup = {}

    # First pass: raw var_names (may be Ensembl IDs or symbols)
    for i, name in enumerate(raw_var_names):
        raw_lookup[str(name).lower()] = i

    # Second pass: if raw var_names are Ensembl IDs, add gene symbol entries.
    # Gene symbols take precedence (overwrite) since queries use symbols.
    if _looks_like_ensembl(np.array(raw_var_names)):
        for col in ["gene_name", "gene_symbol", "symbol", "Gene", "gene_short_name"]:
            if col in adata.raw.var.columns:
                for i, name in enumerate(adata.raw.var[col]):
                    name_lower = str(name).strip().lower()
                    if name_lower and name_lower != "nan":
                        raw_lookup[name_lower] = i
                break

    raw_indices = []
    survived = []
    for j, name in enumerate(query_names):
        if name in raw_lookup:
            raw_indices.append(raw_lookup[name])
            survived.append(True)
        else:
            survived.append(False)

    return (
        np.array(raw_indices, dtype=int),
        np.array(survived, dtype=bool),
    )


def _extract_gene_submatrix(
    adata: ad.AnnData,
    gene_indices: np.ndarray,
    use_raw: bool = True,
    indices_in_raw: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract a dense (n_cells, n_genes) submatrix for the given gene indices.

    When *indices_in_raw* is False (legacy behaviour), gene_indices are in
    adata.var space and are remapped to adata.raw.var when *use_raw* is True.

    When *indices_in_raw* is True, gene_indices are already in adata.raw.var
    space and no remapping is performed.

    Returns:
        X_sub: dense array of shape (n_cells, n_survived_genes)
        survived_mask: boolean mask over gene_indices indicating which
            genes were successfully extracted (True = present in output)
    """
    survived_mask = np.ones(len(gene_indices), dtype=bool)

    if indices_in_raw and adata.raw is not None:
        # Indices already point into adata.raw.var — use directly
        X = adata.raw.X
    elif use_raw and adata.raw is not None:
        X = adata.raw.X
        raw_indices, survived_mask = _map_var_indices_to_raw(adata, gene_indices)
        gene_indices = raw_indices
    else:
        X = adata.X

    n_cells = X.shape[0]
    n_genes = len(gene_indices)

    if n_genes == 0:
        return np.empty((n_cells, 0), dtype=np.float32), survived_mask

    # Validate indices are within bounds
    max_idx = X.shape[1]
    valid = gene_indices < max_idx
    if not np.all(valid):
        # survived_positions maps each entry in gene_indices (post-mapping)
        # back to its position in the original gene_indices array
        survived_positions = np.where(survived_mask)[0]
        for i, is_valid in enumerate(valid):
            if not is_valid:
                survived_mask[survived_positions[i]] = False
        gene_indices = gene_indices[valid]
        n_genes = len(gene_indices)

    if n_genes == 0:
        return np.empty((n_cells, 0), dtype=np.float32), survived_mask

    # For small gene sets, direct column slicing is fine
    if n_genes <= 500:
        X_sub = X[:, gene_indices]
        if sparse.issparse(X_sub):
            result = np.asarray(X_sub.toarray())
        else:
            result = np.asarray(X_sub)
        return result, survived_mask

    # For larger sets, process in chunks to limit peak memory
    chunk_size = 200
    out = np.empty((n_cells, n_genes), dtype=np.float32)
    for start in range(0, n_genes, chunk_size):
        end = min(start + chunk_size, n_genes)
        chunk_idx = gene_indices[start:end]
        X_chunk = X[:, chunk_idx]
        if sparse.issparse(X_chunk):
            out[:, start:end] = np.asarray(X_chunk.toarray())
        else:
            out[:, start:end] = np.asarray(X_chunk)
    return out, survived_mask


def compute_cluster_mean_expression(
    adata: ad.AnnData,
    gene_indices: List[int],
    annotation_col: str,
    min_cells: int = 10,
    use_raw: bool = True,
    indices_in_raw: bool = False,
) -> pd.DataFrame:
    """
    Compute mean expression per cluster for a set of genes.

    When *indices_in_raw* is False (legacy), gene_indices refer to
    positions in ``adata.var`` and are remapped to ``adata.raw.var``
    internally.  When True, they already point into ``adata.raw.var``.

    Returns a DataFrame with shape (n_survived_genes, n_clusters).
    """
    gene_indices_arr = np.array(gene_indices)
    labels = adata.obs[annotation_col].values

    unique_labels, counts = np.unique(labels, return_counts=True)
    valid_labels = unique_labels[counts >= min_cells]

    X_genes, survived_mask = _extract_gene_submatrix(
        adata, gene_indices_arr, use_raw=use_raw, indices_in_raw=indices_in_raw,
    )

    result = {}
    for label in valid_labels:
        mask = labels == label
        result[str(label)] = X_genes[mask, :].mean(axis=0)

    # Use survived_mask to pick the correct gene names
    gene_names = _resolve_gene_names(
        adata, gene_indices_arr[survived_mask], from_raw=indices_in_raw,
    )
    df = pd.DataFrame(result, index=gene_names)
    # Deduplicate: multiple indices can resolve to the same gene symbol
    if df.index.duplicated().any():
        logger.info("compute_cluster_mean_expression: dropping %d duplicate gene names",
                     df.index.duplicated().sum())
        df = df.groupby(df.index).mean()
    return df


def compute_fraction_expressing(
    adata: ad.AnnData,
    gene_indices: List[int],
    annotation_col: str,
    min_cells: int = 10,
    use_raw: bool = True,
    threshold: float = 0.0,
    indices_in_raw: bool = False,
) -> pd.DataFrame:
    """
    Compute fraction of cells expressing each gene (>threshold) per cluster.

    Returns a DataFrame with shape (n_survived_genes, n_clusters).
    """
    gene_indices_arr = np.array(gene_indices)
    labels = adata.obs[annotation_col].values

    unique_labels, counts = np.unique(labels, return_counts=True)
    valid_labels = unique_labels[counts >= min_cells]

    X_genes, survived_mask = _extract_gene_submatrix(
        adata, gene_indices_arr, use_raw=use_raw, indices_in_raw=indices_in_raw,
    )

    result = {}
    for label in valid_labels:
        mask = labels == label
        result[str(label)] = (X_genes[mask, :] > threshold).mean(axis=0)

    gene_names = _resolve_gene_names(
        adata, gene_indices_arr[survived_mask], from_raw=indices_in_raw,
    )
    df = pd.DataFrame(result, index=gene_names)
    if df.index.duplicated().any():
        df = df.groupby(df.index).mean()
    return df


