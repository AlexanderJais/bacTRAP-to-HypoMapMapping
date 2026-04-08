"""
Nature-grade publication figure generation for bacTRAP-to-HypoMap mapping.

All figures follow Nature journal specifications:
- Sans-serif font (Arial/Helvetica), editable text in PDF/SVG
- Axis labels 7pt, tick labels 6pt, panel titles 8pt bold
- Line width 0.5pt axes, 0.75pt plot elements
- Single column 89mm (3.5in) or double column 183mm (7.2in)
- 300 DPI for rasterized elements
- Colorblind-friendly palettes
- White background, no gridlines, no top/right spines
"""

import io
import zipfile
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import Normalize
from matplotlib import cm
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from adjustText import adjust_text
from typing import Optional, List, Tuple, Dict


# ---------------------------------------------------------------------------
# Global style configuration
# ---------------------------------------------------------------------------

def setup_nature_style():
    """Configure matplotlib for Nature-grade figures."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.titleweight": "bold",
        "axes.labelsize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "axes.linewidth": 0.5,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size": 2,
        "ytick.major.size": 2,
        "lines.linewidth": 0.75,
        "patch.linewidth": 0.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.transparent": False,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "pdf.fonttype": 42,  # TrueType fonts in PDF (editable)
        "ps.fonttype": 42,
        "svg.fonttype": "none",  # editable text in SVG
        # Do NOT enable constrained_layout globally — it conflicts with
        # bbox_to_anchor legends and manual colorbar pad placement.
        "figure.constrained_layout.use": False,
    })


def get_figure_width(double_column: bool = False) -> float:
    """Return figure width in inches for Nature format."""
    return 7.2 if double_column else 3.5


def get_qualitative_palette(n: int) -> List[str]:
    """Return a colorblind-friendly qualitative palette."""
    if n <= 10:
        return list(sns.color_palette("tab10", n).as_hex())
    elif n <= 20:
        return list(sns.color_palette("tab20", n).as_hex())
    else:
        base = list(sns.color_palette("tab20", 20).as_hex())
        extra = list(sns.color_palette("Set3", min(n - 20, 12)).as_hex())
        return (base + extra)[:n]


# ---------------------------------------------------------------------------
# Figure A: Correlation Barplot
# ---------------------------------------------------------------------------

def figure_correlation_barplot(
    corr_df: pd.DataFrame,
    top_n: int = 20,
    double_column: bool = False,
) -> plt.Figure:
    """
    Horizontal barplot of top clusters ranked by Spearman correlation.
    Bars colored by correlation strength (sequential colormap).
    """
    setup_nature_style()
    width = get_figure_width(double_column)

    df = corr_df.head(top_n).copy()

    if len(df) == 0:
        fig, ax = plt.subplots(figsize=(width, 2))
        ax.text(0.5, 0.5, "No data available", ha="center", va="center",
                transform=ax.transAxes)
        return fig

    df = df.iloc[::-1]  # reverse for bottom-to-top plotting

    n_bars = len(df)
    height = max(width * 0.6, n_bars * 0.18 + 0.8)

    fig, ax = plt.subplots(figsize=(width, height))

    # Color by Spearman r
    norm = Normalize(
        vmin=min(df["spearman_r"].min(), 0),
        vmax=df["spearman_r"].max(),
    )
    cmap = plt.colormaps["viridis"]
    colors = [cmap(norm(v)) for v in df["spearman_r"]]

    bars = ax.barh(
        range(n_bars),
        df["spearman_r"],
        color=colors,
        edgecolor="none",
        height=0.7,
    )

    # Add r-value text
    for i, (_, row) in enumerate(df.iterrows()):
        r_val = row["spearman_r"]
        offset = 0.005 if r_val >= 0 else -0.005
        ha = "left" if r_val >= 0 else "right"
        ax.text(
            r_val + offset, i, f"{r_val:.3f}",
            va="center", ha=ha, fontsize=5,
        )

    # Vertical dashed line at r=0
    ax.axvline(x=0, color="black", linestyle="--", linewidth=0.5)

    ax.set_yticks(range(n_bars))
    ax.set_yticklabels(df["cluster"], fontsize=6)
    ax.set_xlabel("Spearman correlation (ρ)")
    ax.set_title("bacTRAP enrichment correlation with HypoMap clusters")

    # Colorbar
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.6, aspect=20, pad=0.02)
    cbar.set_label("Spearman ρ", fontsize=6)
    cbar.ax.tick_params(labelsize=5)

    return fig


# ---------------------------------------------------------------------------
# Figure B: UMAP Enrichment Map
# ---------------------------------------------------------------------------

def figure_umap_enrichment(
    umap_coords: np.ndarray,
    cell_labels: np.ndarray,
    enrichment_scores: np.ndarray,
    double_column: bool = True,
    point_size: float = 0.3,
    subsample_idx: Optional[np.ndarray] = None,
    max_legend_items: int = 20,
) -> plt.Figure:
    """
    Two-panel UMAP: left colored by cell-type annotation, right by enrichment score.
    """
    setup_nature_style()
    width = get_figure_width(double_column=True)  # always double for two panels
    height = width * 0.45

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(width, height))

    if subsample_idx is not None:
        umap_coords = umap_coords[subsample_idx]
        cell_labels = cell_labels[subsample_idx]
        enrichment_scores = enrichment_scores[subsample_idx]

    # Shuffle points for fair overlapping
    rng = np.random.default_rng(42)
    order = rng.permutation(len(umap_coords))
    umap_coords = umap_coords[order]
    cell_labels = cell_labels[order]
    enrichment_scores = enrichment_scores[order]

    # --- Left panel: cell-type annotation ---
    unique_labels = np.unique(cell_labels)
    n_labels = len(unique_labels)

    other_label = None
    if n_labels > max_legend_items:
        # Keep top N by frequency, rest grouped as a catch-all
        from collections import Counter
        counts = Counter(cell_labels)
        top_labels = [label for label, _ in counts.most_common(max_legend_items)]
        top_set = set(top_labels)
        other_label = "Other (grouped)"
        plot_labels = np.array([l if l in top_set else other_label for l in cell_labels])
        unique_plot = sorted(set(plot_labels) - {other_label}) + [other_label]
    else:
        plot_labels = cell_labels
        unique_plot = sorted(unique_labels)

    palette = get_qualitative_palette(len(unique_plot))
    color_map = dict(zip(unique_plot, palette))
    if other_label is not None and other_label in color_map:
        color_map[other_label] = "#d3d3d3"
    point_colors = [color_map[l] for l in plot_labels]

    ax1.scatter(
        umap_coords[:, 0], umap_coords[:, 1],
        c=point_colors, s=point_size, alpha=0.6,
        edgecolors="none", rasterized=True,
    )
    ax1.set_xlabel("UMAP1")
    ax1.set_ylabel("UMAP2")
    ax1.set_title("Cell-type annotation")
    ax1.set_xticks([])
    ax1.set_yticks([])
    for spine in ax1.spines.values():
        spine.set_visible(False)

    # Legend outside
    handles = [
        plt.Line2D([0], [0], marker="o", color="w",
                    markerfacecolor=color_map[l], markersize=3, label=l)
        for l in unique_plot
    ]
    ncol = max(1, -(-len(unique_plot) // 15))  # ceiling division by 15
    ax1.legend(
        handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5),
        fontsize=4, frameon=False, ncol=ncol,
        handletextpad=0.2, columnspacing=0.5,
    )

    # --- Right panel: enrichment score ---
    vmin = np.percentile(enrichment_scores, 2)
    vmax = np.percentile(enrichment_scores, 98)

    sc = ax2.scatter(
        umap_coords[:, 0], umap_coords[:, 1],
        c=enrichment_scores, cmap="magma", s=point_size, alpha=0.7,
        edgecolors="none", rasterized=True,
        vmin=vmin, vmax=vmax,
    )
    ax2.set_xlabel("UMAP1")
    ax2.set_ylabel("UMAP2")
    ax2.set_title("bacTRAP enrichment score (PoA)")
    ax2.set_xticks([])
    ax2.set_yticks([])
    for spine in ax2.spines.values():
        spine.set_visible(False)

    cbar = fig.colorbar(sc, ax=ax2, shrink=0.7, aspect=20, pad=0.02)
    cbar.set_label("Enrichment score", fontsize=6)
    cbar.ax.tick_params(labelsize=5)

    return fig


# ---------------------------------------------------------------------------
# Figure C: Marker Overlap Dot Plot
# ---------------------------------------------------------------------------

def figure_dotplot(
    mean_expr: pd.DataFrame,
    frac_expr: pd.DataFrame,
    top_genes: List[str],
    top_clusters: List[str],
    double_column: bool = True,
) -> plt.Figure:
    """
    Dot plot: dot size = fraction expressing, dot color = mean expression.
    Rows = genes, columns = clusters.
    """
    setup_nature_style()

    genes = [g for g in top_genes if g in mean_expr.index and g in frac_expr.index]
    clusters = [c for c in top_clusters if c in mean_expr.columns and c in frac_expr.columns]

    if len(genes) == 0 or len(clusters) == 0:
        fig, ax = plt.subplots(figsize=(3.5, 2))
        ax.text(0.5, 0.5, "No data available", ha="center", va="center")
        return fig

    mean_sub = mean_expr.loc[genes, clusters]
    frac_sub = frac_expr.loc[genes, clusters]

    n_genes = len(genes)
    n_clusters = len(clusters)

    width = get_figure_width(double_column)
    height = max(2, n_genes * 0.2 + 1.5)

    fig, ax = plt.subplots(figsize=(width, height))

    # Scale dot sizes
    max_dot_size = 80
    min_dot_size = 5

    # Vectorized: build coordinate/size/color arrays for a single scatter call
    col_idx, row_idx = np.meshgrid(np.arange(n_clusters), np.arange(n_genes))
    x_flat = col_idx.ravel().astype(float)
    y_flat = row_idx.ravel().astype(float)
    frac_flat = frac_sub.values.ravel()
    mean_flat = mean_sub.values.ravel()
    sizes = min_dot_size + frac_flat * (max_dot_size - min_dot_size)

    vmin = np.nanmin(mean_flat) if np.any(np.isfinite(mean_flat)) else 0
    vmax = np.nanmax(mean_flat) if np.any(np.isfinite(mean_flat)) else 1

    ax.scatter(
        x_flat, y_flat, s=sizes, c=mean_flat, cmap="viridis",
        vmin=vmin, vmax=vmax,
        edgecolors="black", linewidths=0.3, zorder=3,
    )

    ax.set_xticks(range(n_clusters))
    ax.set_xticklabels(clusters, rotation=45, ha="right", fontsize=5)
    ax.set_yticks(range(n_genes))
    ax.set_yticklabels(genes, fontsize=5)

    ax.set_xlim(-0.5, n_clusters - 0.5)
    ax.set_ylim(-0.5, n_genes - 0.5)

    # Light grid
    for i in range(n_genes):
        ax.axhline(i, color="#eeeeee", linewidth=0.3, zorder=0)
    for j in range(n_clusters):
        ax.axvline(j, color="#eeeeee", linewidth=0.3, zorder=0)

    ax.set_title("bacTRAP-enriched gene expression in top HypoMap clusters")

    # Colorbar for mean expression
    sm = cm.ScalarMappable(
        cmap="viridis",
        norm=Normalize(vmin=vmin, vmax=vmax),
    )
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.5, aspect=15, pad=0.02)
    cbar.set_label("Mean expression", fontsize=6)
    cbar.ax.tick_params(labelsize=5)

    # Size legend
    for frac_val, label in [(0.25, "25%"), (0.50, "50%"), (0.75, "75%"), (1.0, "100%")]:
        ax.scatter(
            [], [], s=min_dot_size + frac_val * (max_dot_size - min_dot_size),
            c="gray", edgecolors="black", linewidths=0.3, label=label,
        )
    ax.legend(
        title="% expressing", loc="upper left", bbox_to_anchor=(1.15, 1.0),
        fontsize=5, title_fontsize=5, frameon=False, handletextpad=0.3,
    )

    return fig


# ---------------------------------------------------------------------------
# Figure D: Volcano-style Enrichment Plot (Fisher's test results)
# ---------------------------------------------------------------------------

def figure_volcano_enrichment(
    fisher_df: pd.DataFrame,
    pval_threshold: float = 0.05,
    top_n_labels: int = 10,
    double_column: bool = False,
) -> plt.Figure:
    """
    Volcano-style plot: x = log2(odds ratio), y = -log10(p-value).
    Label top hits.
    """
    setup_nature_style()
    width = get_figure_width(double_column)
    height = width * 0.8

    fig, ax = plt.subplots(figsize=(width, height))

    df = fisher_df.copy()
    # Cap extreme values for plotting
    df["log2_odds_ratio"] = df["log2_odds_ratio"].clip(-10, 10)
    df["neg_log10_pval"] = df["neg_log10_pval"].clip(0, 50)

    sig_mask = df["pvalue"] < pval_threshold
    nonsig = df[~sig_mask]
    sig = df[sig_mask]

    # Non-significant points
    ax.scatter(
        nonsig["log2_odds_ratio"], nonsig["neg_log10_pval"],
        c="#bbbbbb", s=15, alpha=0.6, edgecolors="none",
        label="Not significant",
    )

    # Significant points
    ax.scatter(
        sig["log2_odds_ratio"], sig["neg_log10_pval"],
        c="#d62728", s=20, alpha=0.8, edgecolors="none",
        label=f"p < {pval_threshold}",
    )

    # Significance threshold line
    thresh_y = -np.log10(pval_threshold)
    ax.axhline(y=thresh_y, color="black", linestyle="--", linewidth=0.5, alpha=0.5)
    ax.text(
        ax.get_xlim()[1] * 0.95, thresh_y + 0.3,
        f"p = {pval_threshold}", fontsize=5, ha="right", va="bottom",
    )

    # Label top hits
    top = df.nsmallest(top_n_labels, "pvalue")
    texts = []
    for _, row in top.iterrows():
        texts.append(
            ax.text(
                row["log2_odds_ratio"], row["neg_log10_pval"],
                row["cluster"], fontsize=4.5,
            )
        )
    if len(texts) > 0:
        adjust_text(
            texts, ax=ax,
            arrowprops=dict(arrowstyle="-", color="gray", lw=0.3),
        )

    ax.set_xlabel(r"$\log_2$(odds ratio)")
    ax.set_ylabel(r"$-\log_{10}$(p-value)")
    ax.set_title("Marker gene overlap enrichment (Fisher's exact test)")

    ax.legend(fontsize=5, frameon=False, loc="upper left")

    return fig


# ---------------------------------------------------------------------------
# Figure E: Heatmap
# ---------------------------------------------------------------------------

def figure_heatmap(
    zscore_df: pd.DataFrame,
    double_column: bool = True,
) -> plt.Figure:
    """
    Heatmap of z-scored mean expression for top enriched genes across clusters.
    Rows clustered by hierarchical clustering; columns in ranked order.
    """
    setup_nature_style()

    if zscore_df.empty:
        fig, ax = plt.subplots(figsize=(3.5, 2))
        ax.text(0.5, 0.5, "No data available", ha="center", va="center")
        return fig

    width = get_figure_width(double_column)
    n_genes = len(zscore_df)
    height = max(2.5, n_genes * 0.15 + 1.5)

    # Hierarchical clustering of rows (genes)
    data = zscore_df.values
    if data.shape[0] > 1:
        try:
            Z = linkage(data, method="ward", metric="euclidean")
            row_order = leaves_list(Z)
        except Exception:
            row_order = np.arange(data.shape[0])
    else:
        row_order = np.arange(data.shape[0])

    ordered_df = zscore_df.iloc[row_order]

    fig, ax = plt.subplots(figsize=(width, height))

    vals = ordered_df.values
    finite_vals = vals[np.isfinite(vals)]
    if len(finite_vals) > 0:
        vmax = min(max(abs(finite_vals.min()), abs(finite_vals.max())), 3.0)
    else:
        vmax = 3.0

    im = ax.imshow(
        ordered_df.values,
        aspect="auto",
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        interpolation="nearest",
    )

    ax.set_xticks(range(len(ordered_df.columns)))
    ax.set_xticklabels(ordered_df.columns, rotation=45, ha="right", fontsize=5)
    ax.set_yticks(range(len(ordered_df.index)))
    ax.set_yticklabels(ordered_df.index, fontsize=4.5)

    ax.set_title("Z-scored expression of bacTRAP-enriched genes across top clusters")

    cbar = fig.colorbar(im, ax=ax, shrink=0.6, aspect=20, pad=0.02)
    cbar.set_label("Z-score", fontsize=6)
    cbar.ax.tick_params(labelsize=5)

    return fig


# ---------------------------------------------------------------------------
# Figure F: NNLS Deconvolution Barplot
# ---------------------------------------------------------------------------

def figure_nnls_barplot(
    nnls_df: pd.DataFrame,
    top_n: int = 20,
    double_column: bool = False,
) -> plt.Figure:
    """
    Horizontal barplot of NNLS deconvolution weights per cluster.
    Only shows clusters with non-zero weights.
    """
    setup_nature_style()
    width = get_figure_width(double_column)

    df = nnls_df[nnls_df["weight"] > 1e-6].head(top_n).copy()

    if len(df) == 0:
        fig, ax = plt.subplots(figsize=(width, 2))
        ax.text(0.5, 0.5, "No clusters with non-zero NNLS weights",
                ha="center", va="center", transform=ax.transAxes)
        return fig

    df = df.iloc[::-1]
    n_bars = len(df)
    height = max(width * 0.5, n_bars * 0.18 + 0.8)

    fig, ax = plt.subplots(figsize=(width, height))

    norm = Normalize(vmin=0, vmax=df["weight_norm"].max())
    cmap = plt.colormaps["magma"]
    colors = [cmap(norm(v)) for v in df["weight_norm"]]

    ax.barh(
        range(n_bars), df["weight_norm"], color=colors,
        edgecolor="none", height=0.7,
    )

    for i, (_, row) in enumerate(df.iterrows()):
        ax.text(
            row["weight_norm"] + 0.005, i,
            f"{row['weight_norm']:.3f}",
            va="center", ha="left", fontsize=5,
        )

    ax.set_yticks(range(n_bars))
    ax.set_yticklabels(df["cluster"], fontsize=6)
    ax.set_xlabel("NNLS weight (normalized)")
    ax.set_title("NNLS deconvolution of bacTRAP enrichment profile")

    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.6, aspect=20, pad=0.02)
    cbar.set_label("Normalized weight", fontsize=6)
    cbar.ax.tick_params(labelsize=5)

    return fig


# ---------------------------------------------------------------------------
# Figure G: GSEA Enrichment Curves
# ---------------------------------------------------------------------------

def figure_gsea_curves(
    gsea_df: pd.DataFrame,
    running_scores: Dict[str, np.ndarray],
    ranked_genes: np.ndarray,
    top_n: int = 5,
    double_column: bool = False,
) -> plt.Figure:
    """
    Running enrichment score curves for the top GSEA-enriched clusters.
    """
    setup_nature_style()
    width = get_figure_width(double_column)

    top = gsea_df.head(top_n)
    if len(top) == 0:
        fig, ax = plt.subplots(figsize=(width, 2))
        ax.text(0.5, 0.5, "No GSEA results", ha="center", va="center",
                transform=ax.transAxes)
        return fig

    n_genes = len(ranked_genes)
    height = width * 0.6
    fig, ax = plt.subplots(figsize=(width, height))

    palette = get_qualitative_palette(len(top))
    for i, (_, row) in enumerate(top.iterrows()):
        cluster = row["cluster"]
        if cluster in running_scores:
            curve = running_scores[cluster]
            ax.plot(
                np.arange(len(curve)), curve,
                color=palette[i], linewidth=0.8,
                label=f"{cluster} (NES={row['NES']:.2f})",
            )

    ax.axhline(0, color="black", linewidth=0.3, linestyle="-")
    ax.set_xlabel("Gene rank (by bacTRAP enrichment)")
    ax.set_ylabel("Running enrichment score")
    ax.set_title("GSEA: cluster marker enrichment in bacTRAP-ranked genes")
    ax.legend(fontsize=5, frameon=False, loc="upper right")
    ax.set_xlim(0, n_genes)

    return fig


# ---------------------------------------------------------------------------
# Figure H: GSEA NES Barplot
# ---------------------------------------------------------------------------

def figure_gsea_barplot(
    gsea_df: pd.DataFrame,
    top_n: int = 20,
    double_column: bool = False,
) -> plt.Figure:
    """Horizontal barplot of top clusters by Normalized Enrichment Score."""
    setup_nature_style()
    width = get_figure_width(double_column)

    df = gsea_df.head(top_n).copy()
    if len(df) == 0:
        fig, ax = plt.subplots(figsize=(width, 2))
        ax.text(0.5, 0.5, "No GSEA results", ha="center", va="center",
                transform=ax.transAxes)
        return fig

    df = df.iloc[::-1]
    n_bars = len(df)
    height = max(width * 0.5, n_bars * 0.18 + 0.8)
    fig, ax = plt.subplots(figsize=(width, height))

    sig_mask = df["padj"] < 0.05
    colors = ["#d62728" if s else "#bbbbbb" for s in sig_mask]

    ax.barh(range(n_bars), df["NES"], color=colors, edgecolor="none", height=0.7)
    ax.axvline(0, color="black", linewidth=0.5, linestyle="--")

    for i, (_, row) in enumerate(df.iterrows()):
        label = f"{row['NES']:.2f}"
        if row["padj"] < 0.05:
            label += "*"
        offset = 0.02 if row["NES"] >= 0 else -0.02
        ha = "left" if row["NES"] >= 0 else "right"
        ax.text(row["NES"] + offset, i, label, va="center", ha=ha, fontsize=5)

    ax.set_yticks(range(n_bars))
    ax.set_yticklabels(df["cluster"], fontsize=6)
    ax.set_xlabel("Normalized Enrichment Score (NES)")
    ax.set_title("GSEA: preranked enrichment of cluster markers")

    # Custom legend
    from matplotlib.patches import Patch
    handles = [
        Patch(facecolor="#d62728", label="FDR < 0.05"),
        Patch(facecolor="#bbbbbb", label="Not significant"),
    ]
    ax.legend(handles=handles, fontsize=5, frameon=False, loc="lower right")

    return fig


# ---------------------------------------------------------------------------
# Figure I: AUCell UMAP
# ---------------------------------------------------------------------------

def figure_aucell_umap(
    umap_coords: np.ndarray,
    aucell_scores: np.ndarray,
    double_column: bool = False,
    point_size: float = 0.3,
    subsample_idx: Optional[np.ndarray] = None,
) -> plt.Figure:
    """UMAP colored by AUCell enrichment scores."""
    setup_nature_style()
    width = get_figure_width(double_column)
    height = width * 0.8

    fig, ax = plt.subplots(figsize=(width, height))

    if subsample_idx is not None:
        umap_coords = umap_coords[subsample_idx]
        aucell_scores = aucell_scores[subsample_idx]

    # Shuffle for fair overlap
    rng = np.random.default_rng(42)
    order = rng.permutation(len(umap_coords))
    umap_coords = umap_coords[order]
    aucell_scores = aucell_scores[order]

    vmin = np.percentile(aucell_scores, 2)
    vmax = np.percentile(aucell_scores, 98)

    sc = ax.scatter(
        umap_coords[:, 0], umap_coords[:, 1],
        c=aucell_scores, cmap="magma", s=point_size, alpha=0.7,
        edgecolors="none", rasterized=True,
        vmin=vmin, vmax=vmax,
    )
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("AUCell enrichment score (PoA bacTRAP)")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    cbar = fig.colorbar(sc, ax=ax, shrink=0.7, aspect=20, pad=0.02)
    cbar.set_label("AUCell score", fontsize=6)
    cbar.ax.tick_params(labelsize=5)

    return fig


# ---------------------------------------------------------------------------
# Figure J: Composite Ranking Heatmap
# ---------------------------------------------------------------------------

def figure_composite_ranking(
    composite_df: pd.DataFrame,
    top_n: int = 20,
    double_column: bool = True,
) -> plt.Figure:
    """
    Heatmap showing percentile scores across all methods for top clusters.
    Columns = methods, rows = clusters, color = percentile (0–1).
    """
    setup_nature_style()
    width = get_figure_width(double_column)

    pctl_cols = [c for c in composite_df.columns if c.endswith("_pctl")]
    if len(pctl_cols) == 0 or len(composite_df) == 0:
        fig, ax = plt.subplots(figsize=(width, 2))
        ax.text(0.5, 0.5, "No composite ranking data", ha="center",
                va="center", transform=ax.transAxes)
        return fig

    df = composite_df.head(top_n).copy()
    n_clusters = len(df)
    n_methods = len(pctl_cols)

    height = max(2.5, n_clusters * 0.2 + 1.2)
    fig, ax = plt.subplots(figsize=(width, height))

    data = df[pctl_cols].values
    method_labels = [c.replace("_pctl", "").upper() for c in pctl_cols]

    im = ax.imshow(
        data, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1,
        interpolation="nearest",
    )

    ax.set_xticks(range(n_methods))
    ax.set_xticklabels(method_labels, fontsize=6)
    ax.set_yticks(range(n_clusters))
    ax.set_yticklabels(df["cluster"].values, fontsize=5)

    # Add score text in cells
    for i in range(n_clusters):
        for j in range(n_methods):
            val = data[i, j]
            text_color = "white" if val > 0.6 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=4, color=text_color)

    # Add composite score as right-side annotation
    for i, (_, row) in enumerate(df.iterrows()):
        ax.text(
            n_methods - 0.3, i, f"{row['composite_score']:.2f}",
            ha="left", va="center", fontsize=5, fontweight="bold",
        )

    ax.set_title("Consensus ranking across all methods")

    cbar = fig.colorbar(im, ax=ax, shrink=0.6, aspect=20, pad=0.02)
    cbar.set_label("Percentile score", fontsize=6)
    cbar.ax.tick_params(labelsize=5)

    return fig


# ---------------------------------------------------------------------------
# Export utilities
# ---------------------------------------------------------------------------

def fig_to_bytes(fig: plt.Figure, fmt: str = "pdf") -> bytes:
    """Convert a matplotlib figure to bytes in the specified format."""
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, dpi=300, bbox_inches="tight")
    buf.seek(0)
    return buf.getvalue()


def create_all_figures_zip(figures: Dict[str, plt.Figure]) -> bytes:
    """Create a ZIP archive containing all figures as PDF and SVG."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for name, fig in figures.items():
            pdf_bytes = fig_to_bytes(fig, "pdf")
            svg_bytes = fig_to_bytes(fig, "svg")
            zf.writestr(f"{name}.pdf", pdf_bytes)
            zf.writestr(f"{name}.svg", svg_bytes)
    buf.seek(0)
    return buf.getvalue()
