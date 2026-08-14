import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.transforms as transforms
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as mticker
import statistics
import math
import sys
import subprocess

thresholds=sys.argv[1]
if thresholds!="-none" and thresholds!=("-high") and thresholds.startswith("-thresh:")==False:
    print('enter valid threshold arguments')
    sys.exit(1)

filter=sys.argv[2]
if filter!="-none" and filter.startswith("-f:")==False:
    print('enter valid filtering argument')
    sys.exit(1)
else:
    if filter.startswith("-f:"):
        filter_term=filter.split(":")[1]
    else:
        filter_term=filter

DNA_COLORS = {"A": "#2ecc71", "C": "#3498db", "G": "#f39c12", "T": "#e74c3c",
              "N": "#95a5a6"}
 
AA_COLORS = {
    "A": "#e67e22", "V": "#e67e22", "I": "#e67e22", "L": "#e67e22",
    "M": "#e67e22", "F": "#e67e22", "W": "#e67e22", "P": "#e67e22",
    "S": "#2ecc71", "T": "#2ecc71", "C": "#2ecc71", "Y": "#2ecc71",
    "H": "#2ecc71", "N": "#2ecc71", "Q": "#2ecc71",
    "K": "#3498db", "R": "#3498db",
    "D": "#e74c3c", "E": "#e74c3c",
    "G": "#9b59b6",
    "*": "#95a5a6",
}

def _detect_alphabet(sequences):
    dna_chars = set("ACGTNacgtn")
    return "DNA" if set("".join(sequences)) <= dna_chars else "AA"
 
def _build_pwm(sequences, weights=None):
    seq_len = len(sequences[0])
    alpha_type = _detect_alphabet(sequences)
    all_chars = sorted(set("".join(s.upper() for s in sequences)))
    if alpha_type == "DNA":
        alphabet = [c for c in ["A", "C", "G", "T"] if c in all_chars]
    else:
        alphabet = [c for c in all_chars if c not in ("*",)]
 
    n = len(alphabet)
    idx = {c: i for i, c in enumerate(alphabet)}
 
    if weights is None:
        weights = [1.0] * len(sequences)
    total_weight = sum(weights)
 
    freq = np.zeros((seq_len, n))
    for seq, w in zip(sequences, weights):
        for pos, ch in enumerate(seq.upper()):
            if ch in idx:
                freq[pos, idx[ch]] += w
 
    freq /= total_weight
    freq = (freq + 0.001) / (1 + n * 0.001)   # pseudocount
    # Normalise each row to sum to exactly 1.0 so every column
    # fills the same total height (frequency logo, not IC logo)
    pwm = freq / freq.sum(axis=1, keepdims=True)
    return alphabet, pwm, alpha_type
 
def _letter_patch(letter, x, y, w, h, color, fp):
    """Return a PathPatch of *letter* scaled to the rectangle (x,y,w,h)."""
    tp = TextPath((0, 0), letter, size=1, prop=fp)
    bb = tp.get_extents()
    if bb.width == 0 or bb.height == 0:
        return None
    sx = w / bb.width
    sy = h / bb.height
    tr = (transforms.Affine2D()
          .translate(-bb.x0, -bb.y0)
          .scale(sx, sy)
          .translate(x, y))
    return PathPatch(tr.transform_path(tp), color=color, lw=0)
 
def _draw_logo_clean(ax, alphabet, pwm, alpha_type):
    """Render MEME-style letters onto ax with zero decorations."""
    color_map = DNA_COLORS if alpha_type == "DNA" else AA_COLORS
    fp = FontProperties(family="monospace", weight="bold")
    seq_len = pwm.shape[0]
    max_bits = np.log2(len(alphabet))
 
    ax.set_xlim(0, seq_len)
    ax.set_ylim(0, 1.0)
    ax.set_aspect("auto")
 
    for pos in range(seq_len):
        order = np.argsort(pwm[pos])
        y = 0.0
        for i in order:
            ch = alphabet[i]
            h = pwm[pos, i]
            if h < 1e-6:
                continue
            color = color_map.get(ch, "#888888")
            patch = _letter_patch(ch, pos, y, 1.0, h, color, fp)
            if patch is not None:
                ax.add_patch(patch)
            y += h
 
    # Strip everything
    ax.set_axis_off()
    for spine in ax.spines.values():
        spine.set_visible(False)
 
def _make_logo_fig(sequences, weights=None, figsize=None):
    alphabet, pwm, alpha_type = _build_pwm(sequences, weights)
    seq_len = pwm.shape[0]
    if figsize is None:
        figsize = (seq_len * 0.9 + 0.2, 2.2)
    fig, ax = plt.subplots(figsize=figsize)
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    _draw_logo_clean(ax, alphabet, pwm, alpha_type)
    return fig, ax
 
def d_gene_histogram(
    data,
    title="D Gene Size Distribution",
    xlabel="D Gene Size (bp)",
    ylabel="Number of Occurrences",
    color="#3498db",
    edgecolor="#2c3e50",
    figsize=(14, 7),
    save_path=None,
    show=True,
):
    """
    Bar chart of D gene sizes vs. occurrence counts.
 
    data : [[d_gene_size, occurrences], ...]
    """
    sizes  = [row[0] for row in data]
    counts = [row[1] for row in data]
    order  = np.argsort(sizes)
    sizes  = [sizes[i]  for i in order]
    counts = [counts[i] for i in order]
 
    fig, ax = plt.subplots(figsize=figsize)
    bars = ax.bar(sizes, counts, color=color, edgecolor=edgecolor,
                  linewidth=0.8, width=0.7)
    for bar, cnt in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(counts) * 0.01,
                f"{cnt:,}", ha="center", va="bottom",
                fontsize=9, color="#2c3e50")
    ax.set_xticks(sizes)
    ax.set_xticklabels([str(s) for s in sizes], fontsize=10)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.yaxis.set_major_formatter(
        plt.matplotlib.ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig, ax
 
def meme_from_sequences(
    sequences,
    figsize=None,
    save_path=None,
    show=True,
):
    """
    Clean MEME-style sequence logo — letters only, no axes or legend.
 
    sequences : list[str]  – all equal length; auto-detects DNA vs AA
    """
    if not sequences:
        raise ValueError("sequences list is empty")
    seq_len = len(sequences[0])
    if any(len(s) != seq_len for s in sequences):
        raise ValueError("All sequences must be the same length")
 
    fig, ax = _make_logo_fig(sequences, figsize=figsize)
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight",
                    facecolor="white")
    return fig, ax
 
def meme_from_counts(
    seq_count_pairs,
    figsize=None,
    save_path=None,
    show=True,
):
    """
    Clean MEME-style sequence logo weighted by occurrence counts — letters only.
 
    seq_count_pairs : [[seq, count], ...]  – all seqs must be the same length
    """
    if not seq_count_pairs:
        raise ValueError("seq_count_pairs is empty")
    sequences = [p[0] for p in seq_count_pairs]
    weights   = [float(p[1]) for p in seq_count_pairs]
    seq_len   = len(sequences[0])
    if any(len(s) != seq_len for s in sequences):
        raise ValueError("All sequences must be the same length")
 
    fig, ax = _make_logo_fig(sequences, weights=weights, figsize=figsize)
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight",
                    facecolor="white")
    return fig, ax
 
def _parse(data):
    """
    data row: [heptamer, occurrences, orders, species, haplotypes,
               30-36bp, 33bp, 50-56bp, 53bp]
    Skips the header row if the first element is non-numeric.
    Returns dict of lists keyed by field name.
    """
    rows = [r for r in data if _is_numeric(r[1])]
    return {
        "seqs":   [r[0]        for r in rows],
        "occ":    [float(r[1]) for r in rows],
        "ord":    [float(r[2]) for r in rows],
        "spe":    [float(r[3]) for r in rows],
        "hap":    [float(r[4]) for r in rows],
        "b3036":  [float(r[5]) for r in rows],
        "b33":    [float(r[6]) for r in rows],
        "b5056":  [float(r[7]) for r in rows],
        "b53":    [float(r[8]) for r in rows],
    }


def _is_numeric(v):
    try:
        float(v)
        return True
    except (TypeError, ValueError):
        return False


def _base_ax(fig, ax, x, seqs, xlabel, ylabel, title, max_label_len):
    labels = [s if len(s) <= max_label_len else s[:max_label_len] + "…"
              for s in seqs]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=40, ha="right",
                       fontsize=9, fontfamily="monospace")
    ax.set_xlabel(xlabel, fontsize=12, labelpad=8)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
    ax.set_xlim(-0.6, len(seqs) - 0.4)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_formatter(
        plt.matplotlib.ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _overlay_bars(ax, x, *series_color_label, bar_width=0.7):
    """Draw bars tallest→shortest so smaller bars stay visible on top."""
    items = sorted(
        series_color_label,
        key=lambda t: max(t[0]),
        reverse=True,
    )
    for zorder, (values, color, label) in enumerate(items, start=2):
        ax.bar(x, values, width=bar_width, color=color,
               zorder=zorder, label=label, linewidth=0)


def _count_labels(ax, x, values, max_val):
    for i, v in enumerate(values):
        ax.text(x[i], v + max_val * 0.012, f"{int(v):,}",
                ha="center", va="bottom", fontsize=8.5,
                color="#1a1a2e", fontweight="bold")

C_33    = "#abd9e9"
C_3036  = "#2c7bb6"
C_53    = "#fdae61"
C_5056  = "#d94701"
C_OTHER = "#888780"
C_HAP   = "#1d9e75"
C_SPE   = "#78c679"
C_ORD   = "#d7191c"

LBL_KW = dict(ha="center", va="bottom", fontsize=8, color="#1a1a2e", fontweight="bold")


def occurrence_dgene_barchart(
    data,
    title="Occurrence by D-Gene Group",
    xlabel="Heptamer",
    ylabel="Count",
    figsize=None,
    bar_width=0.7,
    save_path=None,
    show=True,
    max_label_len=14,
):
    """
    Stacked bar chart — one bar per sequence, five stacked segments:
      33bp  |  30-36bp remainder  |  53bp  |  50-56bp remainder  |  other
    Labels (same font/color) on: 30-36bp total, 50-56bp total, other, and
    the grand total above each bar. 33bp and 53bp segments are unlabelled.

    data : [[heptamer, occurrences, orders, species, haplotypes,
              30-36bp, 33bp, 50-56bp, 53bp], ...]
    """
    def parse(data):
        rows = []
        for r in data:
            try:
                float(r[1])
                rows.append(r)
            except (TypeError, ValueError):
                pass
        seqs   = [r[0] for r in rows]
        occ    = [float(r[1]) for r in rows]
        b3036  = [float(r[5]) for r in rows]
        b33    = [float(r[6]) for r in rows]
        b5056  = [float(r[7]) for r in rows]
        b53    = [float(r[8]) for r in rows]
        b3rest = [max(0.0, a - b) for a, b in zip(b3036, b33)]
        b5rest = [max(0.0, a - b) for a, b in zip(b5056, b53)]
        other  = [max(0.0, o - a - b) for o, a, b in zip(occ, b3036, b5056)]
        return seqs, occ, b3036, b33, b3rest, b5056, b53, b5rest, other

    seqs, occ, b3036, b33, b3rest, b5056, b53, b5rest, other = parse(data)
    n = len(seqs)
    if n == 0:
        raise ValueError("data contains no numeric rows")

    x = np.arange(n)
    labels = [s if len(s) <= max_label_len else s[:max_label_len] + "…" for s in seqs]

    if figsize is None:
        figsize = (max(8, n * 1.1 + 2), 6)

    fig, ax = plt.subplots(figsize=figsize)

    # Stack order (bottom to top): 33, 30-36rest, 53, 50-56rest, other
    layers = [
        (b33,    C_33,    None),
        (b3rest, C_3036,  b3036),   # label shows the full 30-36 total
        (b53,    C_53,    None),
        (b5rest, C_5056,  b5056),   # label shows the full 50-56 total
        (other,  C_OTHER, other),
    ]

    bottoms = np.zeros(n)
    seg_mids = []   # (mid_value, label_value_or_None)
    for vals, color, label_vals in layers:
        vals = np.array(vals)
        ax.bar(x, vals, bottom=bottoms, width=bar_width,
               color=color, linewidth=0, zorder=2)
        seg_mids.append((bottoms + vals / 2, label_vals))
        bottoms += vals

    max_val = max(bottoms)
    off_sm  = max_val * 0.008

    # Segment labels (inside bar, white text centred in segment)
    for mid_arr, label_vals in seg_mids:
        if label_vals is None:
            continue
        for i, (mid, lv) in enumerate(zip(mid_arr, label_vals)):
            seg_h = bottoms[i]  # approximate — use segment height check
            # only label if segment is tall enough to read
            seg_val = label_vals[i] if hasattr(label_vals, '__getitem__') else lv
            if seg_val < max_val * 0.04:
                continue
            ax.text(x[i], mid, f"{int(seg_val):,}",
                    ha="center", va="center", fontsize=8,
                    color="black", fontweight="bold")

    # Grand total above each bar (same LBL_KW style)
    for i, total in enumerate(bottoms):
        ax.text(x[i], total + off_sm, f"{int(total):,}", **LBL_KW)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=40, ha="right",
                       fontsize=9, fontfamily="monospace")
    ax.set_xlabel(xlabel, fontsize=12, labelpad=8)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
    ax.set_xlim(-0.6, n - 0.4)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_formatter(
        plt.matplotlib.ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    handles = [
        mpatches.Patch(color=C_33,    label="33bp"),
        mpatches.Patch(color=C_3036,  label="30–36bp"),
        mpatches.Patch(color=C_53,    label="53bp"),
        mpatches.Patch(color=C_5056,  label="50–56bp"),
        mpatches.Patch(color=C_OTHER, label="Other"),
    ]
    ax.legend(handles=handles, loc="upper right",
              framealpha=0.85, fontsize=9, edgecolor="#cccccc")

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    return fig, ax


def taxonomy_barchart(
    data,
    title="Sequence Taxonomy Overview",
    xlabel="Heptamer",
    ylabel="Count",
    figsize=None,
    bar_width=0.7,
    save_path=None,
    show=True,
    max_label_len=14,
):
    """
    Overlaid bars (tallest → shortest so all are visible):
      haplotypes (main) → species → orders
    Labels (same font/color) above each bar's top.

    data : [[heptamer, occurrences, orders, species, haplotypes,
              30-36bp, 33bp, 50-56bp, 53bp], ...]
    """
    def parse(data):
        rows = []
        for r in data:
            try:
                float(r[1])
                rows.append(r)
            except (TypeError, ValueError):
                pass
        return (
            [r[0]        for r in rows],
            [float(r[2]) for r in rows],
            [float(r[3]) for r in rows],
            [float(r[4]) for r in rows],
        )

    seqs, ord_, spe, hap = parse(data)
    n = len(seqs)
    if n == 0:
        raise ValueError("data contains no numeric rows")

    x = np.arange(n)
    labels = [s if len(s) <= max_label_len else s[:max_label_len] + "…" for s in seqs]

    if figsize is None:
        figsize = (max(8, n * 1.1 + 2), 6)

    fig, ax = plt.subplots(figsize=figsize)

    ax.bar(x, hap,  width=bar_width, color=C_HAP, zorder=2, linewidth=0, label="Haplotypes")
    ax.bar(x, spe,  width=bar_width, color=C_SPE, zorder=3, linewidth=0, label="Species")
    ax.bar(x, ord_, width=bar_width, color=C_ORD, zorder=4, linewidth=0, label="Orders")

    max_val = max(hap)
    off     = max_val * 0.012

    for series in (hap, spe, ord_):
        for i, v in enumerate(series):
            ax.text(x[i], v + off, f"{int(v):,}", **LBL_KW)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=40, ha="right",
                       fontsize=9, fontfamily="monospace")
    ax.set_xlabel(xlabel, fontsize=12, labelpad=8)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
    ax.set_xlim(-0.6, n - 0.4)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_formatter(
        plt.matplotlib.ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    handles = [
        mpatches.Patch(color=C_HAP, label="Haplotypes"),
        mpatches.Patch(color=C_SPE, label="Species"),
        mpatches.Patch(color=C_ORD, label="Orders"),
    ]
    ax.legend(handles=handles, loc="upper right",
              framealpha=0.85, fontsize=10, edgecolor="#cccccc")

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    return fig, ax

def seq_taxonomy_barcharts(data, save_prefix=None, show=True, **kwargs):
    """
    Calls both chart functions and returns (fig1, ax1), (fig2, ax2).
    Pass save_prefix="out" to save as out_dgene.png and out_taxonomy.png.
    """
    sp1 = f"{save_prefix}_dgene.png"    if save_prefix else None
    sp2 = f"{save_prefix}_taxonomy.png" if save_prefix else None
    r1 = occurrence_dgene_barchart(data, save_path=sp1, show=show, **kwargs)
    r2 = taxonomy_barchart(data,         save_path=sp2, show=show, **kwargs)
    return r1, r2

def order_gene_barchart(
    data,
    title="Order Gene Number & Distribution",
    ylabel="Number of Genes",
    figsize=None,
    bar_width=0.72,
    save_path=None,
    show=True,
):
    """
    Stacked bar chart (bottom → top):
        Genes 30–36  (green)  – at-33 portion has a black overlay darkening it
        Genes 50–56  (red)    – at-53 portion has a black overlay darkening it
        Other genes  (purple, top)
 
    One colour per segment. The peak-bin subset is darkened by drawing a
    semi-transparent black bar over it — no second colour, no hatching.
 
    Input row: [order_name, n_genes, n_haps, n_species, mean_genes, variance,
                n_30_36, n_at_33, n_50_56, n_at_53, n_other]
    """
    if not data:
        raise ValueError("data is empty")
 
    orders    = [str(row[0])    for row in data]
    n_genes   = [float(row[1])  for row in data]
    n_haps    = [float(row[2])  for row in data]
    n_species = [float(row[3])  for row in data]
    mean_g    = [float(row[4])  for row in data]
    variance  = [float(row[5])  for row in data]
    n_30_36   = [float(row[6])  for row in data]
    n_at_33   = [float(row[7])  for row in data]
    n_50_56   = [float(row[8])  for row in data]
    n_at_53   = [float(row[9])  for row in data]
    n_other   = [float(row[10]) for row in data]
 
    n   = len(orders)
    x   = np.arange(n, dtype=float)
    W   = bar_width
 
    C_30_36 = "#2ecc71"   # green
    C_50_56 = "#e74c3c"   # red
    C_OTHER = "#7b3294"   # purple
    DARK_ALPHA = 0.38     # strength of the black darkening filter
 
    BG   = "#f7f7f7"
    GRID = "#e0e0e0"
 
    if figsize is None:
        figsize = (max(10, n * 2.1 + 2), 9)
 
    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor("white")
    ax.set_facecolor(BG)
 
    n3  = np.array(n_30_36)
    n5  = np.array(n_50_56)
    b50 = n3                 # base of 50-56 block
    bot = n3 + n5            # base of other block
 
    # ── 30-36 block (green), then dark filter over at-33 portion ─────────────
    ax.bar(x, n3, bottom=0, width=W,
           color=C_30_36, zorder=3, linewidth=0)
    ax.bar(x, n_at_33, bottom=0, width=W,
           color="black", alpha=DARK_ALPHA, zorder=4, linewidth=0)
 
    # ── 50-56 block (red), then dark filter over at-53 portion ───────────────
    ax.bar(x, n5, bottom=b50, width=W,
           color=C_50_56, zorder=3, linewidth=0)
    ax.bar(x, n_at_53, bottom=b50, width=W,
           color="black", alpha=DARK_ALPHA, zorder=4, linewidth=0)
 
    # ── other block (purple, top) ─────────────────────────────────────────────
    ax.bar(x, n_other, bottom=bot, width=W,
           color=C_OTHER, zorder=3, linewidth=0)
 
    # ── white dividers between main segments ──────────────────────────────────
    for i in range(n):
        for y_div in [b50[i], bot[i]]:
            ax.plot([x[i] - W/2, x[i] + W/2], [y_div, y_div],
                    color="white", lw=2.0, zorder=6, solid_capstyle="butt")
 
    # ── total count above each bar ────────────────────────────────────────────
    max_h = max(n_genes)
    for i, total in enumerate(n_genes):
        ax.text(x[i], total + max_h * 0.013, f"{int(total):,}",
                ha="center", va="bottom", fontsize=9,
                fontweight="bold", color="#1a1a2e")
 
    # ── x-axis vertical list labels ───────────────────────────────────────────
    labels = []
    for i in range(n):
        lines = [
            orders[i],
            f"{int(n_haps[i])} Haplotypes",
            f"{int(n_species[i])} Species",
            f"Mean {mean_g[i]:.1f} genes per haplotype",
            f"Variance: {variance[i]:.2f}",
        ]
        labels.append("\n".join(lines))
 
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8.5, linespacing=1.8,
                        multialignment="center")
    ax.tick_params(axis="x", length=0, pad=10)
 
    # ── grid & spines ─────────────────────────────────────────────────────────
    ax.yaxis.grid(True, color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#cccccc")
    ax.spines["bottom"].set_color("#cccccc")
 
    ax.set_xlim(-0.6, n - 0.4)
    ax.set_ylabel(ylabel, fontsize=12, labelpad=8)
    ax.set_title(title, fontsize=15, fontweight="bold", pad=16)
    ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    ax.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
 
    # ── legend ────────────────────────────────────────────────────────────────
    # Simulate the darkened look in legend by blending manually
    import matplotlib.colors as mcolors
    def darken(hex_col, amount=DARK_ALPHA):
        rgb = np.array(mcolors.to_rgb(hex_col))
        return tuple(rgb * (1 - amount))
 
    leg_handles = [
        mpatches.Patch(facecolor=C_30_36,           label="Genes 30–36"),
        mpatches.Patch(facecolor=darken(C_30_36),   label="  ↳ Exactly 33"),
        mpatches.Patch(facecolor=C_50_56,           label="Genes 50–56"),
        mpatches.Patch(facecolor=darken(C_50_56),   label="  ↳ Exactly 53"),
        mpatches.Patch(facecolor=C_OTHER,           label="Other genes"),
    ]
    ax.legend(handles=leg_handles, loc="upper right",
              framealpha=0.92, fontsize=9.5,
              edgecolor="#cccccc", ncol=1,
              handlelength=1.6, handleheight=1.3)
 
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight", facecolor="white")
    return fig, ax


d_genes=[]
with open("bird_d_genes.csv","r") as read:
    reader=csv.reader(read)
    header=next(reader)
    for row in reader:
        d_genes.append(row)
    read.close()


lengths=[]
for r in range(51):
    lengths.append([])
for d in d_genes:
    length=len(d[5])
    lengths[length-20].append([d[5],d[8],d[9],d[10],d[11]])

n=20
length_data=[]
seqs_33=[]
seqs_53=[]

up_heptamers_33=[]
up_nonamers_33=[]
down_heptamers_33=[]
down_nonamers_33=[]

up_heptamers_53=[]
up_nonamers_53=[]
down_heptamers_53=[]
down_nonamers_53=[]
for l in lengths:
    print(n,"bp: ",len(l))
    length_data.append([n,len(l)])
    for z in l:
        if n==33:
            seqs_33.append(z[0])
            up_heptamers_33.append(z[1])
            up_nonamers_33.append(z[2])
            down_heptamers_33.append(z[3])
            down_nonamers_33.append(z[4])
        if n==53:
            seqs_53.append(z[0])
            up_heptamers_53.append(z[1])
            up_nonamers_53.append(z[2])
            down_heptamers_53.append(z[3])
            down_nonamers_53.append(z[4])
    n+=1

d_gene_histogram(length_data,save_path="d_gene_lengths_histogram.png", show=False)

'''
# makes MEME plots for 33 and 53 bp genes
if os.path.isdir("33_53_motifs")==False:
    os.mkdir("33_53_motifs")
meme_from_sequences(seqs_33, save_path="33_53_motifs/33bp_d_gene_motif.png", show=False)
meme_from_sequences(up_heptamers_33, save_path="33_53_motifs/33bp_upstream_heptamer_motif.png", show=False)
meme_from_sequences(up_nonamers_33, save_path="33_53_motifs/33bp_upstream_nonamer_motif.png", show=False)
meme_from_sequences(down_heptamers_33, save_path="33_53_motifs/33bp_downstream_heptamer_motif.png", show=False)
meme_from_sequences(down_nonamers_33, save_path="33_53_motifs/33bp_downstream_nonamer_motif.png", show=False)

meme_from_sequences(seqs_53, save_path="33_53_motifs/53bp_d_gene_motif.png", show=False)
meme_from_sequences(up_heptamers_53, save_path="33_53_motifs/53bp_upstream_heptamer_motif.png", show=False)
meme_from_sequences(up_nonamers_53, save_path="33_53_motifs/53bp_upstream_nonamer_motif.png", show=False)
meme_from_sequences(down_heptamers_53, save_path="33_53_motifs/53bp_downstream_heptamer_motif.png", show=False)
meme_from_sequences(down_nonamers_53, save_path="33_53_motifs/53bp_downstream_nonamer_motif.png", show=False)
'''


up_heps=[]
up_nons=[]
down_heps=[]
down_nons=[]
for d in d_genes:
    found=False
    if up_heps!=[]:
        for uh in up_heps:
            if uh[0]==d[8]:
                found=True
                uh[1]+=1
                if d[0].split("/")[0] not in uh[2]:
                    uh[2].append(d[0].split("/")[0])
                if d[0].split("/")[1] not in uh[3]:
                    uh[3].append(d[0].split("/")[1])
                if d[0].split("/")[2] not in uh[4]:
                    uh[4].append(d[0].split("/")[2])
                if 36>=len(d[5])>=30:
                    uh[5]+=1
                    if len(d[5])==33:
                        uh[6]+=1
                elif 56>=len(d[5])>=50:
                    uh[7]+=1
                    if len(d[5])==53:
                        uh[8]+=1                
                break
    if found==False:
        if 36>=len(d[5])>=30:
            x30_36=1
            x50_56=0
            x53=0
            if len(d[5])==33:
                x33=1
            else:
                x33=0
        elif 56>=len(d[5])>=50:
            x30_36=0
            x50_56=1
            x33=0
            if len(d[5])==53:
                x53=1
            else:
                x53=0
        else:
            x30_36=0
            x50_56=0
            x53=0
            x33=0
        up_heps.append([d[8],1,[d[0].split("/")[0]],[d[0].split("/")[1]],[d[0].split("/")[2]],x30_36,x33,x50_56,x53])

    found=False
    if up_nons!=[]:
        for un in up_nons:
            if un[0]==d[9]:
                found=True
                un[1]+=1
                if d[0].split("/")[0] not in un[2]:
                    un[2].append(d[0].split("/")[0])
                if d[0].split("/")[1] not in un[3]:
                    un[3].append(d[0].split("/")[1])
                if d[0].split("/")[2] not in un[4]:
                    un[4].append(d[0].split("/")[2])
                if 36>=len(d[5])>=30:
                    un[5]+=1
                    if len(d[5])==33:
                        un[6]+=1
                elif 56>=len(d[5])>=50:
                    un[7]+=1
                    if len(d[5])==53:
                        un[8]+=1 
                break
    if found==False:
        if 36>=len(d[5])>=30:
            x30_36=1
            x50_56=0
            x53=0
            if len(d[5])==33:
                x33=1
            else:
                x33=0
        elif 56>=len(d[5])>=50:
            x30_36=0
            x50_56=1
            x33=0
            if len(d[5])==53:
                x53=1
            else:
                x53=0
        else:
            x30_36=0
            x50_56=0
            x53=0
            x33=0
        up_nons.append([d[9],1,[d[0].split("/")[0]],[d[0].split("/")[1]],[d[0].split("/")[2]],x30_36,x33,x50_56,x53])
    
    found=False
    if down_heps!=[]:
        for dh in down_heps:
            if dh[0]==d[10]:
                found=True
                dh[1]+=1
                if d[0].split("/")[0] not in dh[2]:
                    dh[2].append(d[0].split("/")[0])
                if d[0].split("/")[1] not in dh[3]:
                    dh[3].append(d[0].split("/")[1])
                if d[0].split("/")[2] not in dh[4]:
                    dh[4].append(d[0].split("/")[2])
                
                if 36>=len(d[5])>=30:
                    dh[5]+=1
                    if len(d[5])==33:
                        dh[6]+=1
                elif 56>=len(d[5])>=50:
                    dh[7]+=1
                    if len(d[5])==53:
                        dh[8]+=1     
                break
    if found==False:
        if 36>=len(d[5])>=30:
            x30_36=1
            x50_56=0
            x53=0
            if len(d[5])==33:
                x33=1
            else:
                x33=0
        elif 56>=len(d[5])>=50:
            x30_36=0
            x50_56=1
            x33=0
            if len(d[5])==53:
                x53=1
            else:
                x53=0
        else:
            x30_36=0
            x50_56=0
            x53=0
            x33=0
        down_heps.append([d[10],1,[d[0].split("/")[0]],[d[0].split("/")[1]],[d[0].split("/")[2]],x30_36,x33,x50_56,x53])
    
    found=False
    if down_nons!=[]:
        for dn in down_nons:
            if dn[0]==d[11]:
                found=True
                dn[1]+=1
                if d[0].split("/")[0] not in dn[2]:
                    dn[2].append(d[0].split("/")[0])
                if d[0].split("/")[1] not in dn[3]:
                    dn[3].append(d[0].split("/")[1])
                if d[0].split("/")[2] not in dn[4]:
                    dn[4].append(d[0].split("/")[2])
                if 36>=len(d[5])>=30:
                    dn[5]+=1
                    if len(d[5])==33:
                        dn[6]+=1
                elif 56>=len(d[5])>=50:
                    dn[7]+=1
                    if len(d[5])==53:
                        dn[8]+=1 
                break
    if found==False:
        if 36>=len(d[5])>=30:
            x30_36=1
            x50_56=0
            x53=0
            if len(d[5])==33:
                x33=1
            else:
                x33=0
        elif 56>=len(d[5])>=50:
            x30_36=0
            x50_56=1
            x33=0
            if len(d[5])==53:
                x53=1
            else:
                x53=0
        else:
            x30_36=0
            x50_56=0
            x53=0
            x33=0
        down_nons.append([d[11],1,[d[0].split("/")[0]],[d[0].split("/")[1]],[d[0].split("/")[2]],x30_36,x33,x50_56,x53])

meme_weight=4      # 4 to weight by haplotypes, 1 to weight by occurances, 2 to weight by orders, 3 to weight by speices

if thresholds=="-none":
    filter_type=4      # 4 to filter by haplotypes, 1 to filter by occurances, 2 to filter by orders, 3 to filter by speices
    sort_type=4        # 4 to sort by haplotypes, 1 to sort by occurances, 2 to sort by orders, 3 to sort by speices
    hep_threshold=0
    non_threshold=0
elif thresholds=="-high":
    filter_type=4      # 4 to filter by haplotypes, 1 to filter by occurances, 2 to filter by orders, 3 to filter by speices
    sort_type=4        # 4 to sort by haplotypes, 1 to sort by occurances, 2 to sort by orders, 3 to sort by speices
    hep_threshold=25
    non_threshold=10
else:
    sort_term = thresholds.split(":")[2]
    if sort_term=="hap":
        filter_type=4
        sort_type=4
    elif sort_term=="species":
        filter_type=3
        sort_type=3
    elif sort_term=="order":
        filter_type=2
        sort_type=2
    elif sort_term=="occur":
        filter_type=1
        sort_type=1
    else:
        print("enter valid sorting term (-hap, -species, -order, -occur)")
        sys.exit(1)
    hep_threshold=int(thresholds.split(":")[1].split("-")[0])
    non_threshold=int(thresholds.split(":")[1].split("-")[1])

with open("upstream_d_heptamers.csv","w",newline="") as write:
    writer=csv.writer(write)
    writer.writerow(["Heptamer","Number of Occurances","Number of Orders","Number of Species","Number of Haplotypes","Number Corresponding to 30-36bp D Genes","Number Corresponding to 33bp D Genes","Number Corresponding to 50-56bp D Genes","Number Corresponding to 53bp D Genes"])
    up_hep_nums=[]
    meme_data=[]
    for x in up_heps:
        if filter_type==1:
            if int(x[filter_type])>hep_threshold:
                up_hep_nums.append([x[0],x[1],len(x[2]),len(x[3]),len(x[4]),x[5],x[6],x[7],x[8]])
                meme_data.append([x[0],len(x[meme_weight])])
        else:
            if int(len(x[filter_type]))>hep_threshold:
                up_hep_nums.append([x[0],x[1],len(x[2]),len(x[3]),len(x[4]),x[5],x[6],x[7],x[8]])
                meme_data.append([x[0],len(x[meme_weight])])
    meme_from_counts(meme_data, save_path="upstream_heptamer_meme.png", show=False)
    up_hep_nums.sort(key=lambda x: x[sort_type], reverse=True)
    for x in up_hep_nums:
        writer.writerow(x)
    seq_taxonomy_barcharts(up_hep_nums,save_prefix="upstream_heptamer")
    write.close()

with open("upstream_d_nonamers.csv","w",newline="") as write:
    writer=csv.writer(write)
    writer.writerow(["Nonamer","Number of Occurances","Number of Orders","Number of Species","Number of Haplotypes","Number Corresponding to 30-36bp D Genes","Number Corresponding to 33bp D Genes","Number Corresponding to 50-56bp D Genes","Number Corresponding to 53bp D Genes"])
    up_non_nums=[]
    meme_data=[]
    for x in up_nons:
        if filter_type==1:
            if int(x[filter_type])>hep_threshold:
                up_non_nums.append([x[0],x[1],len(x[2]),len(x[3]),len(x[4]),x[5],x[6],x[7],x[8]])
                meme_data.append([x[0],len(x[meme_weight])])
        else:
            if int(len(x[filter_type]))>hep_threshold:
                up_non_nums.append([x[0],x[1],len(x[2]),len(x[3]),len(x[4]),x[5],x[6],x[7],x[8]])
                meme_data.append([x[0],len(x[meme_weight])])
    meme_from_counts(meme_data, save_path="upstream_nonamer_meme.png", show=False)
    up_non_nums.sort(key=lambda x: x[sort_type], reverse=True)
    for x in up_non_nums:
        writer.writerow(x)
    seq_taxonomy_barcharts(up_non_nums,save_prefix="upstream_nonamer")
    write.close()

with open("downstream_d_heptamers.csv","w",newline="") as write:
    writer=csv.writer(write)
    writer.writerow(["Heptamer","Number of Occurances","Number of Orders","Number of Species","Number of Haplotypes","Number Corresponding to 30-36bp D Genes","Number Corresponding to 33bp D Genes","Number Corresponding to 50-56bp D Genes","Number Corresponding to 53bp D Genes"])
    down_hep_nums=[]
    meme_data=[]
    for x in down_heps:
        if filter_type==1:
            if int(x[filter_type])>hep_threshold:
                down_hep_nums.append([x[0],x[1],len(x[2]),len(x[3]),len(x[4]),x[5],x[6],x[7],x[8]])
                meme_data.append([x[0],len(x[meme_weight])])
        else:
            if int(len(x[filter_type]))>hep_threshold:
                down_hep_nums.append([x[0],x[1],len(x[2]),len(x[3]),len(x[4]),x[5],x[6],x[7],x[8]])
                meme_data.append([x[0],len(x[meme_weight])])
    meme_from_counts(meme_data, save_path="downstream_heptamer_meme.png", show=False)
    down_hep_nums.sort(key=lambda x: x[sort_type], reverse=True)
    for x in down_hep_nums:
        writer.writerow(x)
    seq_taxonomy_barcharts(down_hep_nums,save_prefix="downstream_heptamer")
    write.close()

with open("downstream_d_nonamers.csv","w",newline="") as write:
    writer=csv.writer(write)
    writer.writerow(["Nonamer","Number of Occurances","Number of Orders","Number of Species","Number of Haplotypes","Number Corresponding to 30-36bp D Genes","Number Corresponding to 33bp D Genes","Number Corresponding to 50-56bp D Genes","Number Corresponding to 53bp D Genes"])
    down_non_nums=[]
    meme_data=[]
    for x in down_nons:
        if filter_type==1:
            if int(x[filter_type])>hep_threshold:
                down_non_nums.append([x[0],x[1],len(x[2]),len(x[3]),len(x[4]),x[5],x[6],x[7],x[8]])
                meme_data.append([x[0],len(x[meme_weight])])
        else:
            if int(len(x[filter_type]))>hep_threshold:
                down_non_nums.append([x[0],x[1],len(x[2]),len(x[3]),len(x[4]),x[5],x[6],x[7],x[8]])
                meme_data.append([x[0],len(x[meme_weight])])
                
    meme_from_counts(meme_data, save_path="downstream_nonamer_meme.png", show=False)
    down_non_nums.sort(key=lambda x: x[sort_type], reverse=True)
    for x in down_non_nums:
        writer.writerow(x)
    seq_taxonomy_barcharts(down_non_nums,save_prefix="downstream_nonamer")
    write.close()


order_names=["Cormorants","Cranes","Doves","Eagles","Falcons","Hummingbirds","Ibises","Landfowl","MiscBirds","Owls","Parrots","Plovers","Songbirds","Suboscines","Waterfowl","Woodpeckers"] #names of all bird order folders
data_dir="/local/storage/kav67/clean_birds"
all_species=[]
all_orders=[]
all_haps=[]
for f in os.listdir(data_dir):
    if f in order_names:
        for f1 in os.listdir(data_dir+"/"+f):
            if os.path.isdir(data_dir+"/"+f+"/"+f1) and f1!="patchworkplot":
                for f2 in os.listdir(data_dir+"/"+f+"/"+f1):
                    if os.path.isdir(data_dir+"/"+f+"/"+f1+"/"+f2):
                        if filter_term=="-none" or filter_term in (f+"/"+f1+"/"+f2):
                            if f not in all_orders:
                                all_orders.append([f])
                            if f1 not in all_species:
                                all_species.append([f,f1])
                            if f2 not in all_haps:
                                all_haps.append([f,f1,f2])

haps=[]
for d in d_genes:
    found=False
    for s in haps:
        if d[0].split("/")[2]==s[2]:
            found=True
            s[3]+=1
            if d[4]=="+":
                s[4]+=1
            elif d[4]=="-":
                s[5]+=1

    if found==False:
        if d[4]=="+":
            pos=1
            neg=0
        elif d[4]=="-":
            pos=0
            neg=1
        haps.append([d[0].split("/")[0],d[0].split("/")[1],d[0].split("/")[2],1,pos,neg])

v_genes=[]
with open("/local/storage/kav67/clean_birds/gene_list.csv","r") as read:
    reader=csv.reader(read)
    header=next(reader)
    for row in reader:
        if str(row[11])=="True" and str(row[7])=="IGH" and str(row[1])=="V":
            v_genes.append(row)
    read.close()

for a in all_haps:
    found=False
    for s in haps:
        if a[2]==s[2]:
            found=True
    if found==False:
        haps.append([a[0],a[1],a[2],0,0,0])

for s in haps:
    v_count=0
    for v in v_genes:
        if (s[0]+"/"+s[1]+"/"+s[2]) == v[0]:
            v_count+=1
    s.insert(3, v_count)


with open("haplotypes_d_genes.csv","w",newline="") as write:
    writer=csv.writer(write)
    writer.writerow(["Order","Species","Haplotype","Number of V Genes","Number of D Genes","Number +","Number -"])
    for s in haps:
        writer.writerow(s)
    write.close()

speices=[]
for d in d_genes:
    found=False
    for s in speices:
        if d[0].split("/")[1]==s[1]:
            found=True
            s[2]+=1
            if d[4]=="+":
                s[3]+=1
            elif d[4]=="-":
                s[4]+=1

    if found==False:
        if d[4]=="+":
            pos=1
            neg=0
        elif d[4]=="-":
            pos=0
            neg=1
        speices.append([d[0].split("/")[0],d[0].split("/")[1],1,int(pos),int(neg)])
for a in all_species:
    found=False
    for s in speices:
        if a[1]==s[1]:
            found=True
    if found==False:
        speices.append([a[0],a[1],0,0,0])

for s in speices:
    v_count=0
    for v in v_genes:
        if (s[0]+"/"+s[1]) in v[0]:
            v_count+=1
    s.insert(2, v_count)

with open("species_d_genes.csv","w",newline="") as write:
    writer=csv.writer(write)
    writer.writerow(["Order","Species","Number of V Genes","Number of D Genes","Number +","Number -"])
    for s in speices:
        writer.writerow(s)
    write.close()

orders=[]
for d in d_genes:
    found=False
    for s in orders:
        if d[0].split("/")[0]==s[0]:
            found=True
            s[1]+=1
            if d[4]=="+":
                s[2]+=1
            elif d[4]=="-":
                s[3]+=1

    if found==False:
        if d[4]=="+":
            pos=1
            neg=0
        elif d[4]=="-":
            pos=0
            neg=1
        orders.append([d[0].split("/")[0],1,pos,neg])

for s in orders:
    v_count=0
    for v in v_genes:
        if s[0] in v[0]:
            v_count+=1
    s.insert(1, v_count)

with open("order_d_genes.csv","w",newline="") as write:
    writer=csv.writer(write)
    writer.writerow(["Order","Number of V Genes","Number of D Genes","Number +","Number -"])
    for s in orders:
        writer.writerow(s)
    write.close()


if os.path.isdir("order_data")==False:
    os.mkdir("order_data")

order_data_list=[]
for order in order_names:
    order_genes=[]
    order_30_36=0
    order_33=0
    order_50_56=0
    order_53=0
    motif_33=[]
    motif_53=[]
    order_speices=[]
    order_haps=[]
    for d in d_genes:
        if order in d[0]:
            if d[0].split("/")[1] not in order_speices:
                order_speices.append(d[0].split("/")[1])
            if d[0].split("/")[2] not in order_haps:
                order_haps.append(d[0].split("/")[2])
            order_genes.append(d)
            if 36>=len(d[5])>=30:
                order_30_36+=1
                if len(d[5])==33:
                    order_33+=1
                    motif_33.append(d[5])
            if 56>=len(d[5])>=50:
                order_50_56+=1
                if len(d[5])==53:
                    order_53+=1
                    motif_53.append(d[5])
    hap_nums=[]
    for h in order_haps:
        hap_num=0
        for d in d_genes:
            if h in d[0]:
                hap_num+=1
        hap_nums.append(hap_num)
    try:
        varience = math.sqrt(statistics.variance(hap_nums))
    except:
        varience=0
    
    try:
        mean_genes = len(order_genes)/len(order_haps)
    except:
        mean_genes=0

    other = len(order_genes)-order_30_36-order_50_56
    order_data=[order,len(order_genes),len(order_haps),len(order_speices),mean_genes,varience,order_30_36,order_33,order_50_56,order_53,other]
    order_data_list.append(order_data)
order_data_list.sort(key=lambda x: x[1],reverse=True)
order_gene_barchart(order_data_list, save_path="order_d_gene_barchart.png", show=False)

subprocess.run(["python","double_gene_checker.py"])