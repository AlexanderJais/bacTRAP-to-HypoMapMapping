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
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from adjustText import adjust_text
from typing import Optional, List, Dict

logger = logging.getLogger(__name__)


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


def _add_umap_axis_arrows(
    ax,
    x_label: str = "UMAP1",
    y_label: str = "UMAP2",
    length: float = 0.14,
    origin: tuple = (0.02, 0.02),
    linewidth: float = 0.9,
    fontsize: float = 6,
) -> None:
    """Draw two small axis arrows in the bottom-left corner of a UMAP panel.

    Replaces the conventional x/y axes on dimensionality-reduction plots with
    the compact convention common in single-cell publications: two arrows
    anchored at the bottom-left, labelled "UMAP1" / "UMAP2". *length* and
    *origin* are expressed in axes fraction coordinates, so the arrows scale
    with the panel and stay in the corner regardless of data range.

    Call after all data have been plotted so annotations sit on top.
    """
    x0, y0 = origin
    arrow_style = dict(
        arrowstyle="-|>,head_length=3,head_width=2",
        linewidth=linewidth,
        color="black",
        shrinkA=0, shrinkB=0,
    )
    ax.annotate(
        "", xy=(x0 + length, y0), xytext=(x0, y0),
        xycoords="axes fraction", textcoords="axes fraction",
        arrowprops=arrow_style,
    )
    ax.annotate(
        "", xy=(x0, y0 + length), xytext=(x0, y0),
        xycoords="axes fraction", textcoords="axes fraction",
        arrowprops=arrow_style,
    )
    ax.text(
        x0 + length + 0.005, y0, x_label,
        transform=ax.transAxes, ha="left", va="center", fontsize=fontsize,
    )
    ax.text(
        x0, y0 + length + 0.005, y_label,
        transform=ax.transAxes, ha="center", va="bottom", fontsize=fontsize,
        rotation=90, rotation_mode="anchor",
    )


# ---------------------------------------------------------------------------
# Supplementary Figure S1: Correlation Barplot
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

    # Determine significance status per bar (both Pearson & Spearman p < 0.05)
    has_sig_col = "both_significant" in df.columns
    is_sig = df["both_significant"].values if has_sig_col else np.ones(n_bars, dtype=bool)

    bars = ax.barh(
        range(n_bars),
        df["spearman_r"],
        color=colors,
        edgecolor="none",
        height=0.7,
    )

    # Hatch non-significant bars (Pearson p >= 0.05)
    for i, bar in enumerate(bars):
        if not is_sig[i]:
            bar.set_hatch("//")
            bar.set_edgecolor("grey")
            bar.set_alpha(0.6)

    # Add r-value text (with n.s. annotation for non-significant clusters)
    for i, (_, row) in enumerate(df.iterrows()):
        r_val = row["spearman_r"]
        offset = 0.005 if r_val >= 0 else -0.005
        ha = "left" if r_val >= 0 else "right"
        label = f"{r_val:.3f}"
        if has_sig_col and not row["both_significant"]:
            label += " (n.s.)"
        ax.text(
            r_val + offset, i, label,
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
# Supplementary Figure S2: UMAP Enrichment Map
# ---------------------------------------------------------------------------

def figure_umap_enrichment(
    umap_coords: np.ndarray,
    cell_labels: np.ndarray,
    enrichment_scores: np.ndarray,
    double_column: bool = True,
    point_size: float = 0.3,
    subsample_idx: Optional[np.ndarray] = None,
    max_legend_items: int = 20,
    score_title: str = "bacTRAP enrichment score (PoA)",
    score_label: str = "Enrichment score",
    highlight_clusters: Optional[List[str]] = None,
) -> plt.Figure:
    """
    Two-panel UMAP: left colored by cell-type annotation, right by enrichment score.

    *score_title* / *score_label* let the caller override the right-panel
    title and colorbar text so the same function can render either the
    z-scored mean signature or the AUCell score without duplicating code.

    When *highlight_clusters* is provided, those clusters are coloured in the
    left panel (in the order supplied) and all other cells are greyed out.
    This keeps the annotation panel consistent with per-cluster ranking
    figures (e.g. AUCell 1b/1c) — otherwise the default "top N by cell
    frequency" heuristic hides small, highly-enriched populations.
    """
    logger.info("figure_umap_enrichment: %d cells, %d unique labels, subsample=%s",
                len(umap_coords), len(np.unique(cell_labels)),
                len(subsample_idx) if subsample_idx is not None else "none")
    setup_nature_style()
    width = get_figure_width(double_column=True)  # always double for two panels
    height = width * 0.5

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(width, height),
                                    gridspec_kw={"wspace": 0.8})

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
    if highlight_clusters is not None and len(highlight_clusters) > 0:
        # Caller-driven selection: colour the supplied clusters in order,
        # grey out everything else. Used by the AUCell UMAP to match the
        # top-N cluster set shown in the ranking figures.
        top_labels = [c for c in highlight_clusters if c in set(unique_labels)]
        top_set = set(top_labels)
        other_label = "Other"
        plot_labels = np.array([l if l in top_set else other_label for l in cell_labels])
        unique_plot = top_labels + [other_label]
    elif n_labels > max_legend_items:
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

    palette = get_qualitative_palette(max(len(unique_plot) - (1 if other_label else 0), 1))
    # Preserve caller-supplied ordering when highlight_clusters is used so
    # the legend matches the 1b/1c rank order; otherwise zip normally.
    color_map = {}
    for i, l in enumerate(unique_plot):
        if l == other_label:
            continue
        color_map[l] = palette[i % len(palette)]
    if other_label is not None:
        color_map[other_label] = "#d3d3d3"
    point_colors = np.array([color_map[l] for l in plot_labels])

    if highlight_clusters is not None and other_label is not None:
        # Draw the grey "Other" layer first, then each highlighted cluster
        # on top — otherwise small populations (e.g. Chat.GABA-7) get buried
        # under ~380k grey cells from the rng-shuffled concatenation.
        is_other = plot_labels == other_label
        ax1.scatter(
            umap_coords[is_other, 0], umap_coords[is_other, 1],
            c=point_colors[is_other], s=point_size, alpha=0.4,
            edgecolors="none", rasterized=True,
        )
        for cl in top_labels:
            mask = plot_labels == cl
            if not mask.any():
                continue
            ax1.scatter(
                umap_coords[mask, 0], umap_coords[mask, 1],
                c=[color_map[cl]], s=point_size * 2.5, alpha=0.95,
                edgecolors="none", rasterized=True,
            )
    else:
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

    # Legend below the left panel — avoids overlapping with right panel
    handles = [
        plt.Line2D([0], [0], marker="o", color="w",
                    markerfacecolor=color_map[l], markersize=3, label=l)
        for l in unique_plot
    ]
    ncol = max(2, -(-len(unique_plot) // 10))  # spread across columns
    ax1.legend(
        handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.05),
        fontsize=3.5, frameon=False, ncol=ncol,
        handletextpad=0.1, columnspacing=0.3, labelspacing=0.2,
    )

    # --- Right panel: enrichment score ---
    vmin = np.nanpercentile(enrichment_scores, 2)
    vmax = np.nanpercentile(enrichment_scores, 98)

    sc = ax2.scatter(
        umap_coords[:, 0], umap_coords[:, 1],
        c=enrichment_scores, cmap="magma", s=point_size, alpha=0.7,
        edgecolors="none", rasterized=True,
        vmin=vmin, vmax=vmax,
    )
    ax2.set_xlabel("UMAP1")
    ax2.set_ylabel("UMAP2")
    ax2.set_title(score_title)
    ax2.set_xticks([])
    ax2.set_yticks([])
    for spine in ax2.spines.values():
        spine.set_visible(False)

    cbar = fig.colorbar(sc, ax=ax2, shrink=0.7, aspect=20, pad=0.02)
    cbar.set_label(score_label, fontsize=6)
    cbar.ax.tick_params(labelsize=5)

    return fig


# ---------------------------------------------------------------------------
# Supplementary Figure S4: Marker Overlap Dot Plot
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
    logger.info("figure_dotplot: %d genes requested, %d clusters requested, "
                "mean_expr=%s, frac_expr=%s",
                len(top_genes), len(top_clusters), mean_expr.shape, frac_expr.shape)
    setup_nature_style()

    genes = [g for g in top_genes if g in mean_expr.index and g in frac_expr.index]
    clusters = [c for c in top_clusters if c in mean_expr.columns and c in frac_expr.columns]
    logger.info("  after filtering: %d/%d genes, %d/%d clusters available",
                len(genes), len(top_genes), len(clusters), len(top_clusters))

    if len(genes) == 0 or len(clusters) == 0:
        logger.warning("  no genes or clusters available — empty dotplot")
        fig, ax = plt.subplots(figsize=(3.5, 2))
        ax.text(0.5, 0.5, "No data available", ha="center", va="center")
        return fig

    mean_sub = mean_expr.loc[genes, clusters]
    frac_sub = frac_expr.loc[genes, clusters]
    logger.info("  frac_expr range: [%.3f, %.3f], mean_expr range: [%.3f, %.3f]",
                float(frac_sub.values.min()), float(frac_sub.values.max()),
                float(mean_sub.values.min()), float(mean_sub.values.max()))

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
# Supplementary Figure S3: Volcano-style Enrichment Plot (Fisher's test results)
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

    if len(fisher_df) == 0:
        ax.text(0.5, 0.5, "No data available", ha="center", va="center",
                transform=ax.transAxes)
        return fig

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
# bacTRAP Gene Volcano Plot (Data Overview)
# ---------------------------------------------------------------------------

def figure_bactrap_volcano(
    bactrap_matched: pd.DataFrame,
    highlight_genes: Optional[List[str]] = None,
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 1.0,
    top_n_labels: int = 15,
    double_column: bool = False,
) -> plt.Figure:
    """
    Classic volcano plot of bacTRAP DESeq2 results.
    x = log2FoldChange, y = -log10(padj).
    Highlights significantly enriched genes and optionally labels specific genes.
    """
    setup_nature_style()
    width = get_figure_width(double_column)
    height = width * 0.8

    fig, ax = plt.subplots(figsize=(width, height))

    if len(bactrap_matched) == 0:
        ax.text(0.5, 0.5, "No data available", ha="center", va="center",
                transform=ax.transAxes)
        return fig

    df = bactrap_matched.copy()
    n_before = len(df)
    df = df.dropna(subset=["log2FoldChange", "padj"])
    logger.info("figure_bactrap_volcano: %d genes (%d dropped for NaN), highlight=%s",
                len(df), n_before - len(df), highlight_genes)
    df["neg_log10_padj"] = -np.log10(df["padj"].clip(lower=1e-300))
    df["neg_log10_padj"] = df["neg_log10_padj"].clip(upper=50)

    # Classify points
    sig_up = (df["padj"] < padj_cutoff) & (df["log2FoldChange"] > log2fc_cutoff)
    sig_down = (df["padj"] < padj_cutoff) & (df["log2FoldChange"] < -log2fc_cutoff)
    nonsig = ~sig_up & ~sig_down
    logger.info("  sig_up=%d, sig_down=%d, nonsig=%d", sig_up.sum(), sig_down.sum(), nonsig.sum())

    # Plot non-significant
    ax.scatter(
        df.loc[nonsig, "log2FoldChange"], df.loc[nonsig, "neg_log10_padj"],
        c="#bbbbbb", s=8, alpha=0.4, edgecolors="none", zorder=1,
    )
    # Plot significant down
    ax.scatter(
        df.loc[sig_down, "log2FoldChange"], df.loc[sig_down, "neg_log10_padj"],
        c="#4575b4", s=12, alpha=0.7, edgecolors="none", zorder=2,
        label="Down-regulated",
    )
    # Plot significant up (enriched in IP)
    ax.scatter(
        df.loc[sig_up, "log2FoldChange"], df.loc[sig_up, "neg_log10_padj"],
        c="#d62728", s=12, alpha=0.7, edgecolors="none", zorder=2,
        label="Enriched in IP",
    )

    # Threshold lines
    thresh_y = -np.log10(padj_cutoff)
    ax.axhline(y=thresh_y, color="black", linestyle="--", linewidth=0.4, alpha=0.4)
    ax.axvline(x=log2fc_cutoff, color="black", linestyle="--", linewidth=0.4, alpha=0.4)
    ax.axvline(x=-log2fc_cutoff, color="black", linestyle="--", linewidth=0.4, alpha=0.4)

    # Label highlight genes (e.g. Pnoc) — always label these regardless of significance
    if highlight_genes is None:
        highlight_genes = []
    highlight_set = set(g.lower() for g in highlight_genes)

    # Auto-label top enriched genes + forced highlights
    top_up = df[sig_up].nlargest(top_n_labels, "neg_log10_padj")
    genes_to_label = set(top_up["_hypomap_gene_name"].tolist())

    # Add highlight genes
    for _, row in df.iterrows():
        gname = str(row.get("_hypomap_gene_name", ""))
        if gname.lower() in highlight_set:
            genes_to_label.add(gname)

    texts = []
    for _, row in df.iterrows():
        gname = str(row.get("_hypomap_gene_name", ""))
        if gname in genes_to_label:
            is_highlight = gname.lower() in highlight_set
            texts.append(
                ax.text(
                    row["log2FoldChange"], row["neg_log10_padj"],
                    gname, fontsize=5 if is_highlight else 4.5,
                    fontweight="bold" if is_highlight else "normal",
                    color="#d62728" if is_highlight else "black",
                )
            )
            # Mark highlight genes with a ring
            if is_highlight:
                ax.scatter(
                    [row["log2FoldChange"]], [row["neg_log10_padj"]],
                    s=50, facecolors="none", edgecolors="#d62728",
                    linewidths=1.0, zorder=5,
                )

    if len(texts) > 0:
        adjust_text(
            texts, ax=ax,
            arrowprops=dict(arrowstyle="-", color="gray", lw=0.3),
        )

    ax.set_xlabel(r"$\log_2$(Fold Change)")
    ax.set_ylabel(r"$-\log_{10}$(adjusted p-value)")
    ax.set_title("bacTRAP translational profiling (PoA IP vs Input)")
    ax.legend(fontsize=5, frameon=False, loc="upper left")

    return fig


# ---------------------------------------------------------------------------
# Supplementary Figure S5: Heatmap
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
# Supplementary Figure S6: NNLS Deconvolution Barplot
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
# Supplementary Figure S7: GSEA Enrichment Curves
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
    ax.legend(fontsize=5, frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5))
    ax.set_xlim(0, n_genes)
    fig.subplots_adjust(right=0.7)

    return fig


# ---------------------------------------------------------------------------
# Supplementary Figure S8: GSEA NES Barplot
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
# Main Figure 1a: AUCell UMAP
# ---------------------------------------------------------------------------

def figure_aucell_umap(
    umap_coords: np.ndarray,
    aucell_scores: np.ndarray,
    double_column: bool = False,
    point_size: float = 0.3,
    subsample_idx: Optional[np.ndarray] = None,
) -> plt.Figure:
    """Publication-ready UMAP coloured by AUCell enrichment scores.

    Title and axis labels are omitted (belong in the figure caption) and the
    conventional x/y axes are replaced with two bottom-left arrows via
    `_add_umap_axis_arrows`, matching the single-cell publication convention.
    """
    setup_nature_style()
    width = get_figure_width(double_column)
    height = width * 0.9

    fig, ax = plt.subplots(figsize=(width, height))

    if subsample_idx is not None:
        umap_coords = umap_coords[subsample_idx]
        aucell_scores = aucell_scores[subsample_idx]

    # Shuffle for fair overlap
    rng = np.random.default_rng(42)
    order = rng.permutation(len(umap_coords))
    umap_coords = umap_coords[order]
    aucell_scores = aucell_scores[order]

    vmin = np.nanpercentile(aucell_scores, 2)
    vmax = np.nanpercentile(aucell_scores, 98)

    sc = ax.scatter(
        umap_coords[:, 0], umap_coords[:, 1],
        c=aucell_scores, cmap="magma", s=point_size, alpha=0.7,
        edgecolors="none", rasterized=True,
        vmin=vmin, vmax=vmax,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    _add_umap_axis_arrows(ax)

    cbar = fig.colorbar(sc, ax=ax, shrink=0.6, aspect=20, pad=0.02)
    cbar.set_label("AUCell score", fontsize=6)
    cbar.ax.tick_params(labelsize=5)
    cbar.outline.set_linewidth(0.4)

    return fig


# ---------------------------------------------------------------------------
# Main Figure: Cell-type annotation UMAP (top-N highlighted)
# ---------------------------------------------------------------------------

def figure_celltype_umap(
    umap_coords: np.ndarray,
    cell_labels: np.ndarray,
    highlight_clusters: List[str],
    double_column: bool = False,
    point_size: float = 0.3,
    subsample_idx: Optional[np.ndarray] = None,
) -> plt.Figure:
    """Publication-ready UMAP coloured by cell-type annotation.

    The clusters in *highlight_clusters* (in the supplied order — typically
    top-N by AUCell mean) are drawn in colour on top of a grey "Other" layer,
    so small but highly enriched populations remain visible. Axis labels and
    title are omitted; a legend sits to the right of the panel and the two
    bottom-left arrows mark UMAP1 / UMAP2.
    """
    setup_nature_style()
    width = get_figure_width(double_column)
    height = width * 0.9

    fig, ax = plt.subplots(figsize=(width, height))

    if subsample_idx is not None:
        umap_coords = umap_coords[subsample_idx]
        cell_labels = cell_labels[subsample_idx]

    rng = np.random.default_rng(42)
    order = rng.permutation(len(umap_coords))
    umap_coords = umap_coords[order]
    cell_labels = cell_labels[order]

    unique = set(np.unique(cell_labels).tolist())
    top_labels = [c for c in highlight_clusters if c in unique]
    top_set = set(top_labels)

    palette = get_qualitative_palette(max(len(top_labels), 1))
    color_map = {cl: palette[i % len(palette)] for i, cl in enumerate(top_labels)}
    other_color = "#d9d9d9"

    is_other = np.array([lbl not in top_set for lbl in cell_labels])
    ax.scatter(
        umap_coords[is_other, 0], umap_coords[is_other, 1],
        c=other_color, s=point_size, alpha=0.4,
        edgecolors="none", rasterized=True,
    )
    for cl in top_labels:
        mask = cell_labels == cl
        if not mask.any():
            continue
        ax.scatter(
            umap_coords[mask, 0], umap_coords[mask, 1],
            c=[color_map[cl]], s=point_size * 2.5, alpha=0.95,
            edgecolors="none", rasterized=True,
        )

    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    _add_umap_axis_arrows(ax)

    handles = [
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=color_map[cl], markersize=3.5, label=cl)
        for cl in top_labels
    ]
    handles.append(
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=other_color, markersize=3.5, label="Other")
    )
    ax.legend(
        handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5),
        fontsize=5, frameon=False, handletextpad=0.3,
        labelspacing=0.35, borderaxespad=0,
    )

    logger.info("figure_celltype_umap: %d highlighted clusters",
                len(top_labels))

    return fig


# ---------------------------------------------------------------------------
# Main Figure 1b: AUCell Cluster Barplot
# ---------------------------------------------------------------------------

def figure_aucell_cluster_barplot(
    aucell_scores: np.ndarray,
    cell_labels: np.ndarray,
    top_n: int = 25,
    double_column: bool = False,
) -> plt.Figure:
    """Horizontal barplot of mean AUCell score per cluster, ranked."""
    setup_nature_style()
    width = get_figure_width(double_column)

    # Compute mean AUCell score per cluster
    df = pd.DataFrame({"score": aucell_scores, "cluster": cell_labels})
    cluster_stats = df.groupby("cluster")["score"].agg(["mean", "std", "count"])
    cluster_stats = cluster_stats.sort_values("mean", ascending=False)
    cluster_stats = cluster_stats.head(top_n).iloc[::-1]  # reverse for bottom-to-top

    if len(cluster_stats) == 0:
        fig, ax = plt.subplots(figsize=(width, 2))
        ax.text(0.5, 0.5, "No data available", ha="center", va="center",
                transform=ax.transAxes)
        return fig

    n_bars = len(cluster_stats)
    height = max(width * 0.6, n_bars * 0.18 + 0.8)
    fig, ax = plt.subplots(figsize=(width, height))

    # Color by score
    norm = Normalize(vmin=0, vmax=cluster_stats["mean"].max())
    cmap = plt.colormaps["magma"]
    colors = [cmap(norm(v)) for v in cluster_stats["mean"]]

    ax.barh(
        range(n_bars),
        cluster_stats["mean"],
        xerr=cluster_stats["std"] / np.sqrt(cluster_stats["count"]),  # SEM
        color=colors,
        edgecolor="none",
        height=0.7,
        capsize=1.5,
        error_kw={"linewidth": 0.5},
    )

    for i, (_, row) in enumerate(cluster_stats.iterrows()):
        ax.text(row["mean"] + cluster_stats["mean"].max() * 0.02, i,
                f"{row['mean']:.4f}", va="center", ha="left", fontsize=5)

    ax.set_yticks(range(n_bars))
    ax.set_yticklabels(cluster_stats.index, fontsize=6)
    ax.set_xlabel("Mean AUCell score")
    ax.set_title("AUCell enrichment per cluster (top %d)" % top_n)

    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.6, aspect=20, pad=0.02)
    cbar.set_label("AUCell score", fontsize=6)
    cbar.ax.tick_params(labelsize=5)

    logger.info("figure_aucell_cluster_barplot: %d clusters, top=%s (%.4f)",
                n_bars, cluster_stats.index[-1], cluster_stats["mean"].iloc[-1])

    return fig


# ---------------------------------------------------------------------------
# Main Figure 1c: AUCell Violin Plot (top clusters)
# ---------------------------------------------------------------------------

def figure_aucell_violins(
    aucell_scores: np.ndarray,
    cell_labels: np.ndarray,
    top_n: int = 15,
    double_column: bool = True,
) -> plt.Figure:
    """Violin plots of AUCell score distributions for top clusters."""
    setup_nature_style()
    width = get_figure_width(double_column)

    df = pd.DataFrame({"score": aucell_scores, "cluster": cell_labels})
    cluster_means = df.groupby("cluster")["score"].mean().sort_values(ascending=False)
    top_clusters = cluster_means.head(top_n).index.tolist()
    df_top = df[df["cluster"].isin(top_clusters)].copy()

    if len(df_top) == 0:
        fig, ax = plt.subplots(figsize=(width, 2))
        ax.text(0.5, 0.5, "No data available", ha="center", va="center",
                transform=ax.transAxes)
        return fig

    # Order by mean score (descending)
    df_top["cluster"] = pd.Categorical(df_top["cluster"], categories=top_clusters, ordered=True)

    height = max(width * 0.5, 3.5)
    fig, ax = plt.subplots(figsize=(width, height))

    parts = ax.violinplot(
        [df_top.loc[df_top["cluster"] == c, "score"].values for c in top_clusters],
        positions=range(len(top_clusters)),
        vert=False,
        showmeans=True,
        showmedians=True,
        showextrema=False,
    )

    # Style violins
    cmap = plt.colormaps["magma"]
    norm = Normalize(vmin=0, vmax=cluster_means.iloc[0])
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(cmap(norm(cluster_means.iloc[i])))
        body.set_alpha(0.7)
        body.set_edgecolor("grey")
        body.set_linewidth(0.5)
    if "cmeans" in parts:
        parts["cmeans"].set_linewidth(0.8)
        parts["cmeans"].set_color("black")
    if "cmedians" in parts:
        parts["cmedians"].set_linewidth(0.5)
        parts["cmedians"].set_color("grey")
        parts["cmedians"].set_linestyle("--")

    ax.set_yticks(range(len(top_clusters)))
    ax.set_yticklabels(top_clusters, fontsize=6)
    ax.set_xlabel("AUCell score")
    ax.tick_params(axis="x", labelsize=6)

    # Highest-mean cluster at the top of the plot (top_clusters[0]) rather
    # than at y=0 which matplotlib renders at the bottom.
    ax.invert_yaxis()

    # Explicit legend for mean/median — without it the two vertical ticks
    # inside each horizontal violin read as an ambiguous "I"-shape.
    legend_handles = [
        plt.Line2D([0], [0], color="black", linewidth=0.8, label="mean"),
        plt.Line2D([0], [0], color="grey", linewidth=0.5, linestyle="--", label="median"),
    ]
    ax.legend(
        handles=legend_handles, loc="lower right", fontsize=5,
        frameon=False, handlelength=1.5, handletextpad=0.4,
    )

    logger.info("figure_aucell_violins: %d clusters shown", len(top_clusters))

    return fig


# ---------------------------------------------------------------------------
# Main Figure 1d: AUCell Score Histogram
# ---------------------------------------------------------------------------

def figure_aucell_histogram(
    aucell_scores: np.ndarray,
    double_column: bool = False,
) -> plt.Figure:
    """Histogram of AUCell scores across all cells, with percentile markers."""
    setup_nature_style()
    width = get_figure_width(double_column)
    fig, ax = plt.subplots(figsize=(width, width * 0.6))

    # Remove zero scores for cleaner visualization
    nonzero = aucell_scores[aucell_scores > 0]
    all_scores = aucell_scores

    ax.hist(all_scores, bins=100, color="steelblue", edgecolor="none",
            alpha=0.8, density=True)

    # Mark percentiles
    for pct, ls, lbl in [(90, "--", "90th"), (95, "-.", "95th"), (99, ":", "99th")]:
        val = np.percentile(all_scores, pct)
        ax.axvline(val, color="firebrick", linestyle=ls, linewidth=0.8, alpha=0.8)
        ax.text(val, ax.get_ylim()[1] * 0.95, f" {lbl}\n {val:.4f}",
                fontsize=5, color="firebrick", va="top")

    mean_val = np.mean(all_scores)
    ax.axvline(mean_val, color="black", linestyle="-", linewidth=0.8)
    ax.text(mean_val, ax.get_ylim()[1] * 0.80, f" mean\n {mean_val:.4f}",
            fontsize=5, color="black", va="top")

    ax.set_xlabel("AUCell score")
    ax.set_ylabel("Density")
    ax.set_title("AUCell score distribution (all cells)")

    n_zero = int((all_scores == 0).sum())
    n_total = len(all_scores)
    ax.text(0.98, 0.98,
            f"n = {n_total:,}\nzero = {n_zero:,} ({100*n_zero/n_total:.1f}%)",
            transform=ax.transAxes, fontsize=5, va="top", ha="right",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="grey", alpha=0.8))

    logger.info("figure_aucell_histogram: %d cells, mean=%.4f, 95th=%.4f",
                n_total, mean_val, np.percentile(all_scores, 95))

    return fig


# ---------------------------------------------------------------------------
# Main Figure 1e: Composite Consensus Ranking Heatmap
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
            n_methods + 0.3, i, f"{row['composite_score']:.2f}",
            ha="left", va="center", fontsize=5, fontweight="bold",
        )

    ax.set_title("Consensus ranking across all methods")

    cbar = fig.colorbar(im, ax=ax, shrink=0.6, aspect=20, pad=0.02)
    cbar.set_label("Percentile score", fontsize=6)
    cbar.ax.tick_params(labelsize=5)

    return fig


# ---------------------------------------------------------------------------
# Sanity-check figure: Cre-driver gene expression across top-ranked clusters
# ---------------------------------------------------------------------------

def figure_marker_gene_diagnostic(
    gene_stats: pd.DataFrame,
    cluster_order: List[str],
    gene_name: str = "Pnoc",
    fraction_threshold: float = 0.05,
    double_column: bool = True,
) -> plt.Figure:
    """
    Two-panel sanity-check plot for a Cre-driver / marker gene.

    For each cluster in *cluster_order* (typically the top hits from the
    composite ranking) two horizontal bars are drawn:

      * Left  — mean log-normalised expression of *gene_name*
      * Right — fraction of cells with non-zero counts for *gene_name*

    A vertical reference line on the fraction panel marks
    *fraction_threshold*; bars meeting both "ranked top" AND
    "fraction ≥ threshold" are filled in saturated colour, the rest are
    greyed-out — making it easy to spot top-ranked clusters that fail the
    Pnoc check (likely lineage-tracing artefacts, dropout, or background).

    Parameters
    ----------
    gene_stats : DataFrame
        Output of ``compute_single_gene_cluster_stats`` — index = cluster,
        columns include ``mean_expr`` and ``fraction_expressing``.
    cluster_order : list[str]
        Clusters to display (top-down, e.g. the top 20 from composite
        ranking).  Missing clusters are skipped silently.
    """
    setup_nature_style()
    width = get_figure_width(double_column)

    available = [c for c in cluster_order if c in gene_stats.index]
    if len(available) == 0:
        fig, ax = plt.subplots(figsize=(width, 2))
        ax.text(0.5, 0.5, f"{gene_name} not detected in selected clusters",
                ha="center", va="center", transform=ax.transAxes)
        return fig

    # Reverse so the top-ranked cluster sits at the top of the plot
    df = gene_stats.loc[available, ["mean_expr", "fraction_expressing"]].iloc[::-1]
    n_bars = len(df)

    height = max(2.0, n_bars * 0.22 + 1.0)
    fig, (ax_mean, ax_frac) = plt.subplots(
        1, 2, figsize=(width, height), sharey=True,
        gridspec_kw={"wspace": 0.08},
    )

    pass_mask = (df["fraction_expressing"] >= fraction_threshold).values
    color_pass = "#762a83"   # saturated purple, colorblind-safe
    color_fail = "#bdbdbd"   # neutral grey
    bar_colors = np.where(pass_mask, color_pass, color_fail)

    # --- Panel 1: mean expression ----------------------------------------
    ax_mean.barh(
        range(n_bars), df["mean_expr"].values,
        color=bar_colors, edgecolor="none", height=0.7,
    )
    ax_mean.set_yticks(range(n_bars))
    ax_mean.set_yticklabels(df.index, fontsize=5)
    ax_mean.invert_xaxis()                       # bars grow leftward
    ax_mean.yaxis.tick_right()                   # labels live in the gutter
    ax_mean.tick_params(axis="y", which="both", length=0, pad=2)
    ax_mean.set_xlabel(f"Mean {gene_name} expression\n(log-norm)")
    # Tighten x-limit so the label row reads as 0 → max
    mean_max = float(df["mean_expr"].max()) if n_bars else 0.0
    if mean_max > 0:
        ax_mean.set_xlim(mean_max * 1.05, 0)

    # --- Panel 2: fraction expressing ------------------------------------
    ax_frac.barh(
        range(n_bars), df["fraction_expressing"].values,
        color=bar_colors, edgecolor="none", height=0.7,
    )
    ax_frac.axvline(
        fraction_threshold, color="black", linewidth=0.5,
        linestyle="--", zorder=4,
    )
    ax_frac.set_xlabel(f"Fraction of cells\nexpressing {gene_name}")
    ax_frac.set_xlim(0, max(1.0, float(df["fraction_expressing"].max()) * 1.1))
    # Hide redundant y tick labels on the right panel — they live with
    # the left axis (rotated to the right side via yaxis.tick_right).
    ax_frac.tick_params(axis="y", which="both", length=0, labelleft=False)

    # Annotate fraction values
    for i, frac in enumerate(df["fraction_expressing"].values):
        ax_frac.text(
            frac + 0.01, i, f"{frac*100:.0f}%",
            va="center", ha="left", fontsize=5,
        )

    fig.suptitle(
        f"{gene_name} expression across top-ranked clusters "
        f"(threshold = {fraction_threshold*100:.0f}%)",
        fontsize=8, fontweight="bold", y=0.995,
    )

    # Legend: colour-coded pass / fail
    pass_patch = plt.Rectangle((0, 0), 1, 1, color=color_pass)
    fail_patch = plt.Rectangle((0, 0), 1, 1, color=color_fail)
    ax_frac.legend(
        [pass_patch, fail_patch],
        [f"≥ {fraction_threshold*100:.0f}% expressing", "below threshold"],
        loc="lower right", fontsize=5, frameon=False, handlelength=1.2,
    )

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


