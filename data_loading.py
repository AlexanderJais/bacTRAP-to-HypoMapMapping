"""
Data loading and gene matching utilities for bacTRAP-to-HypoMap mapping.
"""

import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import streamlit as st
from scipy import sparse
from typing import Tuple, Optional, Dict, List


@st.cache_resource(show_spinner="Loading HypoMap atlas (this may take a few minutes)...")
def load_hypomap(file_path: str) -> ad.AnnData:
    """Load the HypoMap h5ad atlas with memory-efficient settings."""
    adata = sc.read_h5ad(file_path, backed=None)
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
    df = pd.read_excel(file_path, engine="openpyxl")
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


def get_gene_names_from_adata(adata: ad.AnnData) -> np.ndarray:
    """Extract gene names from the AnnData object, trying multiple locations."""
    # First try var_names directly
    gene_names = adata.var_names.values.copy()

    # Check if var_names look like Ensembl IDs; if so, look for a symbol column
    if len(gene_names) > 0 and str(gene_names[0]).startswith("ENSMUSG"):
        for col in ["gene_name", "gene_symbol", "symbol", "Gene", "gene_short_name"]:
            if col in adata.var.columns:
                gene_names = adata.var[col].values.copy()
                break

    return gene_names


def match_genes(
    bactrap_df: pd.DataFrame,
    adata: ad.AnnData,
    gene_col: str = "gene_name",
) -> Tuple[pd.DataFrame, List[str], Dict[str, int]]:
    """
    Match bacTRAP gene symbols to HypoMap var_names.

    Returns:
        bactrap_matched: subset of bactrap_df with matched genes
        matched_gene_names: list of matched gene symbols (as they appear in HypoMap)
        gene_to_adata_idx: mapping from gene name to index in adata.var
    """
    if gene_col not in bactrap_df.columns:
        raise ValueError(f"Column '{gene_col}' not found in bacTRAP data.")

    # Get gene names from HypoMap
    adata_gene_names = get_gene_names_from_adata(adata)

    # Build case-insensitive lookup: lowercase -> (original_name, index)
    adata_gene_lookup = {}
    for idx, name in enumerate(adata_gene_names):
        name_str = str(name).strip()
        adata_gene_lookup[name_str.lower()] = (name_str, idx)

    # Match bacTRAP genes
    matched_rows = []
    matched_gene_names = []
    gene_to_adata_idx = {}

    for _, row in bactrap_df.iterrows():
        bt_gene = str(row[gene_col]).strip()
        bt_gene_lower = bt_gene.lower()

        if bt_gene_lower in adata_gene_lookup:
            original_name, adata_idx = adata_gene_lookup[bt_gene_lower]
            matched_rows.append(row)
            matched_gene_names.append(original_name)
            gene_to_adata_idx[original_name] = adata_idx

    if len(matched_rows) == 0:
        empty_df = bactrap_df.iloc[:0].copy()
        empty_df["_hypomap_gene_name"] = pd.Series(dtype=str)
        return empty_df, [], {}

    bactrap_matched = pd.DataFrame(matched_rows)
    bactrap_matched = bactrap_matched.reset_index(drop=True)
    bactrap_matched["_hypomap_gene_name"] = matched_gene_names

    return bactrap_matched, matched_gene_names, gene_to_adata_idx


def _resolve_gene_names(adata: ad.AnnData, gene_indices: np.ndarray) -> List[str]:
    """Resolve gene indices in adata.var to display names."""
    adata_gene_names = get_gene_names_from_adata(adata)
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
    if len(raw_var_names) > 0 and str(raw_var_names[0]).startswith("ENSMUSG"):
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
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract a dense (n_cells, n_genes) submatrix for the given gene indices.

    gene_indices are always in adata.var space. When use_raw=True and
    adata.raw exists, they are remapped to adata.raw.var space.

    Returns:
        X_sub: dense array of shape (n_cells, n_survived_genes)
        survived_mask: boolean mask over gene_indices indicating which
            genes were successfully extracted (True = present in output)
    """
    survived_mask = np.ones(len(gene_indices), dtype=bool)

    if use_raw and adata.raw is not None:
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
) -> pd.DataFrame:
    """
    Compute mean expression per cluster for a set of genes.

    gene_indices refer to positions in adata.var. When use_raw=True,
    they are remapped internally to adata.raw.var.

    Returns a DataFrame with shape (n_survived_genes, n_clusters).
    """
    gene_indices_arr = np.array(gene_indices)
    labels = adata.obs[annotation_col].values

    unique_labels, counts = np.unique(labels, return_counts=True)
    valid_labels = unique_labels[counts >= min_cells]

    X_genes, survived_mask = _extract_gene_submatrix(
        adata, gene_indices_arr, use_raw=use_raw,
    )

    result = {}
    for label in valid_labels:
        mask = labels == label
        result[str(label)] = X_genes[mask, :].mean(axis=0)

    # Use survived_mask to pick the correct gene names
    gene_names = _resolve_gene_names(adata, gene_indices_arr[survived_mask])
    return pd.DataFrame(result, index=gene_names)


def compute_fraction_expressing(
    adata: ad.AnnData,
    gene_indices: List[int],
    annotation_col: str,
    min_cells: int = 10,
    use_raw: bool = True,
    threshold: float = 0.0,
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
        adata, gene_indices_arr, use_raw=use_raw,
    )

    result = {}
    for label in valid_labels:
        mask = labels == label
        result[str(label)] = (X_genes[mask, :] > threshold).mean(axis=0)

    gene_names = _resolve_gene_names(adata, gene_indices_arr[survived_mask])
    return pd.DataFrame(result, index=gene_names)


