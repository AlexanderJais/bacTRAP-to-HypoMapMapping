"""Synthetic-data unit tests for the AUCell-false-positive-reduction changes.

Fast (< 30 s total) and self-contained — they never touch the real HypoMap
atlas.  Run with::

    python -m pytest tests/test_changes.py -q
"""

import numpy as np
import pandas as pd
import pytest

from analysis import (
    filter_signature_genes_by_atlas,
    compute_aucell_scores,
    compute_aucell_scores_multi,
    compute_empirical_null_aucell,
)
from data_loading import get_neuronal_cell_mask


# ---------------------------------------------------------------------------
# Change 2 — signature refinement filters
# ---------------------------------------------------------------------------

def test_filter_signature_genes_by_atlas_planted_scenarios():
    clusters = ["c1", "c2", "c3", "c4"]
    # Three planted genes:
    #  - Specific:    high in one cluster, near-zero elsewhere, decent detection
    #  - Broad:       high in every cluster (pan-neuronal-like)
    #  - Undetectable: ~0 everywhere, near-zero detection rate
    cluster_mean_expr = pd.DataFrame(
        {
            "c1": [2.0, 1.0, 0.0],
            "c2": [0.01, 1.1, 0.0],
            "c3": [0.0, 0.9, 0.005],
            "c4": [0.0, 1.3, 0.0],
        },
        index=["Specific", "Broad", "Undetectable"],
    )
    cell_detection_rate = pd.Series(
        {"Specific": 0.35, "Broad": 0.7, "Undetectable": 0.001}
    )
    kept, drop_log = filter_signature_genes_by_atlas(
        ["Specific", "Broad", "Undetectable"],
        cluster_mean_expr, cell_detection_rate,
    )
    assert kept == ["Specific"]
    status = dict(zip(drop_log["gene"], drop_log["status"]))
    assert status["Specific"] == "kept"
    assert "broadly_expressed" in status["Broad"]
    assert ("low_detection_rate" in status["Undetectable"]
            or "low_max_cluster_mean" in status["Undetectable"])


def test_filter_signature_genes_by_atlas_disabled_is_noop():
    cluster_mean_expr = pd.DataFrame(
        {"c1": [5.0], "c2": [5.0]}, index=["G"],
    )
    cell_detection_rate = pd.Series({"G": 0.0})
    kept, drop_log = filter_signature_genes_by_atlas(
        ["G"], cluster_mean_expr, cell_detection_rate,
        apply_detectability=False, apply_specificity=False,
    )
    assert kept == ["G"]
    assert (drop_log["status"] == "kept").all()


def test_filter_signature_genes_by_atlas_transpose_autodetect():
    # cluster x gene orientation should be auto-transposed
    cluster_by_gene = pd.DataFrame(
        {"Specific": [2.0, 0.0, 0.0], "Broad": [1.0, 1.1, 0.9]},
        index=["c1", "c2", "c3"],
    )
    det = pd.Series({"Specific": 0.3, "Broad": 0.7})
    kept, _ = filter_signature_genes_by_atlas(
        ["Specific", "Broad"], cluster_by_gene, det,
    )
    assert kept == ["Specific"]


# ---------------------------------------------------------------------------
# Change 3 — neuronal-only atlas mask
# ---------------------------------------------------------------------------

class _FakeAdata:
    """Minimal stand-in exposing only what get_neuronal_cell_mask reads."""

    def __init__(self, obs):
        self.obs = obs
        self.obs_names = obs.index


def test_get_neuronal_cell_mask_basic():
    obs = pd.DataFrame(
        {"C7_named": [
            "GABA-1", "GLU-2", "Oligodendrocytes", "Astrocytes",
            "Immune", "Endothelial", "GABA-3", "Tanycytes", "ParsTuber",
        ]},
        index=[f"cell{i}" for i in range(9)],
    )
    mask = get_neuronal_cell_mask(_FakeAdata(obs))
    assert mask.dtype == bool
    assert list(mask) == [True, True, False, False, False, False, True, False, False]
    assert mask.sum() == 3


def test_get_neuronal_cell_mask_fallback_to_c25():
    obs = pd.DataFrame(
        {"C25_named": ["GABA-1", "Microglia", "GLU-1"]},
        index=["a", "b", "c"],
    )
    mask = get_neuronal_cell_mask(_FakeAdata(obs))
    assert list(mask) == [True, False, True]


def test_get_neuronal_cell_mask_missing_columns_returns_all_true():
    obs = pd.DataFrame({"something_else": ["x", "y", "z"]}, index=["a", "b", "c"])
    mask = get_neuronal_cell_mask(_FakeAdata(obs))
    assert mask.all()
    assert len(mask) == 3


# ---------------------------------------------------------------------------
# Change 1 — empirical-null AUCell reproducibility
# ---------------------------------------------------------------------------

def _toy_adata(n_cells=200, n_genes=500, n_clusters=5, seed=1):
    import anndata as ad
    from scipy import sparse

    rng = np.random.default_rng(seed)
    # Negative-binomial-ish raw counts, mostly zeros
    X = rng.poisson(0.4, size=(n_cells, n_genes)).astype(np.float32)
    # Make a handful of "signature" genes co-vary with one cluster
    clusters = np.array([f"cl{i % n_clusters}" for i in range(n_cells)])
    sig_genes = [f"g{i}" for i in range(20)]
    for gi in range(20):
        boost = (clusters == "cl0").astype(np.float32) * rng.poisson(3, size=n_cells)
        X[:, gi] += boost
    var = pd.DataFrame(index=[f"g{i}" for i in range(n_genes)])
    obs = pd.DataFrame({"cluster": clusters}, index=[f"c{i}" for i in range(n_cells)])
    adata = ad.AnnData(X=sparse.csr_matrix(X), obs=obs, var=var)
    adata.raw = adata
    return adata, sig_genes


def test_compute_aucell_scores_multi_matches_single():
    adata, sig_genes = _toy_adata()
    # A few signatures of varying composition (all small relative to the τ
    # window, so the shared n_top equals each one's single-scorer n_top).
    sigs = [
        sig_genes,                         # the planted 20-gene signature
        [f"g{i}" for i in range(30, 55)],  # 25 unrelated genes
        ["g0", "g1", "g2", "MISSING_GENE"],  # tiny + an unmatched name
    ]
    multi = compute_aucell_scores_multi(adata, sigs, seed=0)
    assert multi.shape == (adata.n_obs, len(sigs))
    for j, s in enumerate(sigs):
        single = compute_aucell_scores(adata, s, seed=0)
        assert np.array_equal(multi[:, j], single), f"signature {j} mismatch"
    # empty list → (n_cells, 0)
    assert compute_aucell_scores_multi(adata, [], seed=0).shape == (adata.n_obs, 0)


def test_compute_empirical_null_aucell_reproducible():
    adata, sig_genes = _toy_adata()
    labels = adata.obs["cluster"]

    def run():
        # fresh uns each time so the gene-mean cache doesn't leak between runs
        adata.uns = {}
        return compute_empirical_null_aucell(
            adata, sig_genes, labels, compute_aucell_scores,
            n_control_sets=10, n_bins=5, seed=0,
            top_fraction=0.05, min_cluster_size=5,
        )

    out1 = run()
    out2 = run()
    assert not out1.empty
    assert list(out1.columns) == [
        "null_mean", "null_sd", "z_empirical", "pvalue_empirical",
        "qvalue_empirical", "n_control_sets_used",
    ]
    pd.testing.assert_frame_equal(out1, out2)
    # The cluster the signature was planted in should have the highest z.
    z = out1["z_empirical"].dropna()
    if len(z):
        assert z.idxmax() == "cl0"


def test_compute_empirical_null_aucell_different_seed_changes_z():
    adata, sig_genes = _toy_adata()
    labels = adata.obs["cluster"]
    adata.uns = {}
    a = compute_empirical_null_aucell(
        adata, sig_genes, labels, compute_aucell_scores,
        n_control_sets=10, n_bins=5, seed=0, min_cluster_size=5,
    )
    adata.uns = {}
    b = compute_empirical_null_aucell(
        adata, sig_genes, labels, compute_aucell_scores,
        n_control_sets=10, n_bins=5, seed=123, min_cluster_size=5,
    )
    # Different control draws → at least some z-scores differ.
    assert not np.allclose(
        a["z_empirical"].fillna(0).to_numpy(),
        b["z_empirical"].fillna(0).to_numpy(),
    )


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
