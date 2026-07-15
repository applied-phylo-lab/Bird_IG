import csv
import sys
from collections import Counter
import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
from matplotlib.patches import FancyBboxPatch, ConnectionPatch
import copy

order_names=["Cormorants","Cranes","Doves","Eagles","Falcons","finches","house_finches","Hummingbirds","Ibises","Landfowl","MiscBirds","Owls","Parrots","Plovers","Songbirds","Suboscines","Waterfowl","Woodpeckers","Capuchino_Seedeaters"] #names of all bird order folders
data_dir="/local/storage/kav67/clean_birds" #location of igdetective data

delete=sys.argv[1]
if delete!="-n" and delete!="-d":
    print("Enter valid delete argument (-n to keep existing, -d to delete existing)")
    sys.exit(1)

genome_paths=[]
with open("/local/storage/kav67/clean_birds/bird_genome_paths.csv","r") as read:
    reader=csv.reader(read)
    header=next(reader)
    for row in reader:
        genome_paths.append(row)
    read.close()

contig_list=[]
with open("/local/storage/kav67/clean_birds/ig_contig_list.csv","r") as read:
    reader=csv.reader(read)
    header=next(reader)
    for row in reader:
        contig_list.append(row)
    read.close()

rss_list=[]
rss_fold="/local/storage/dhardesty/assemblies/#rss_scripts/#rss_procedure"
for dfold in os.listdir(rss_fold):
    if os.path.isdir(rss_fold+"/"+dfold) and "all_birds" in dfold:
        for dfile in os.listdir(rss_fold+"/"+dfold):
            if os.path.isfile(rss_fold+"/"+dfold+"/"+dfile) and "combined_rss_zones" in dfile:
                with open(rss_fold+"/"+dfold+"/"+dfile,"r") as dread:
                    reader=csv.reader(dread)
                    header=next(reader)
                    for row in reader:
                        rss_list.append(row)
                    dread.close()

rss_filter=sys.argv[2]
if rss_filter!="-n":
    with open(rss_filter,"r") as rf_read:
        reader=csv.reader(rf_read)
        valid_rss=[]
        for row in reader:
            valid_rss.append(str(row[0]))
        rf_read.close()
    temp_list=copy.deepcopy(rss_list)
    for r in temp_list:
        if str(r[9]) not in valid_rss:
            num=0
            for r1 in rss_list:
                if r==r1:
                    rss_list.pop(num)
                    break
                num+=1
else:
    print("No RSS filter, labelling all matches to reference")

contig_filter=sys.argv[3]
if contig_filter!='-none' and contig_filter.startswith("-contig:")==False:
    print('enter valid contig filter argument')
    sys.exit(1)

filter=sys.argv[4]
if filter=="-h":
    haplotype_filter=sys.argv[5]
    species_filter="None"
    order_filter="None"
elif filter=="-s":
    haplotype_filter="None"
    species_filter=sys.argv[5]
    order_filter="None"
elif filter=="-o":
    haplotype_filter="None"
    species_filter="None"
    order_filter=sys.argv[5]
elif filter=="-n":
    haplotype_filter="None"
    species_filter="None"
    order_filter="None"
else:
    print("enter valid filtering argument")

def is_low_complexity(seq, threshold=0.7):
    """
    Returns True if the sequence is considered low complexity.
    A sequence is low complexity if the most frequent base or the sum
    of the two most frequent bases exceed the threshold proportion.
    """
    seq = str(seq).upper()
    if len(seq) == 0:
        return True

    counts = Counter(seq)
    freqs = sorted(counts.values(), reverse=True)
    
    if freqs[0] / len(seq) >= threshold:
        return True
    if len(freqs) > 1 and sum(freqs[:2]) / len(seq) >= threshold:
        return True
    return False


def _spread_labels(positions, min_gap):
    pos = np.array(positions, dtype=float)
    for _ in range(500):
        moved = False
        for i in range(len(pos) - 1):
            gap = pos[i + 1] - pos[i]
            if gap < min_gap:
                push = (min_gap - gap) / 2.0
                pos[i]     -= push
                pos[i + 1] += push
                moved = True
        if not moved:
            break
    return pos
 
 
# ── locus color palettes ─────────────────────────────────────────────────────
# Each palette holds the 4 colors used to render a single locus_type:
#   V_with    -> V gene that is labeled (has RSS)
#   V_without -> V gene that is unlabeled (no RSS)
#   D         -> D gene
#   J         -> J gene
#
# Palette[0] is the ORIGINAL fixed color scheme from the base script. It is
# always what gets used when only one locus_type is present in the data, so
# single-locus plots are pixel-for-pixel identical to the old behavior.
_LOCUS_PALETTES = [
    {"V_with": "#172869", "V_without": "#EA7580", "D": "#088BBE", "J": "#F8CD9C"},  # original scheme
    {"V_with": "#1B5E20", "V_without": "#EF6C00", "D": "#00838F", "J": "#FFCA28"},  # locus 2
    {"V_with": "#4A148C", "V_without": "#D81B60", "D": "#00695C", "J": "#FFB300"},  # locus 3
    {"V_with": "#6A1B9A", "V_without": "#C62828", "D": "#0277BD", "J": "#F4511E"},  # locus 4
    {"V_with": "#263238", "V_without": "#AD1457", "D": "#00ACC1", "J": "#FDD835"},  # locus 5
    {"V_with": "#3E2723", "V_without": "#E65100", "D": "#006064", "J": "#FFF176"},  # locus 6
]
 
 
def _generate_palette(seed_rgb):
    """Derive a 4-color V/D/J palette from a base RGB color (0-1 floats)."""
    h, l, s = colorsys.rgb_to_hls(*seed_rgb[:3])
 
    def _shade(l_delta, s_delta=0.0):
        rr, gg, bb = colorsys.hls_to_rgb(
            h,
            min(max(l + l_delta, 0.05), 0.95),
            min(max(s + s_delta, 0.15), 1.0),
        )
        return "#{:02x}{:02x}{:02x}".format(int(rr * 255), int(gg * 255), int(bb * 255))
 
    return {
        "V_with":    _shade(-0.15),
        "V_without": _shade(0.20, -0.10),
        "D":         _shade(-0.05, 0.10),
        "J":         _shade(0.28, -0.20),
    }
 
 
def _build_locus_color_map(locus_types_present):
    """
    Map each locus_type -> a 4-color palette dict (V_with/V_without/D/J).
 
    - Exactly one locus_type present -> it gets the ORIGINAL palette
      (index 0), so single-locus behavior/coloring is unchanged.
    - More than one locus_type present -> each locus_type (sorted
      alphabetically, for deterministic output) gets its own distinct
      palette. If there are more loci than predefined palettes, extra
      palettes are generated procedurally from a colormap so the function
      never runs out.
    """
    loci_sorted = sorted(locus_types_present)
 
    if len(loci_sorted) <= 1:
        return {locus: _LOCUS_PALETTES[0] for locus in loci_sorted}
 
    palettes = list(_LOCUS_PALETTES)
    if len(loci_sorted) > len(palettes):
        cmap = cm.get_cmap("tab10")
        for i in range(len(palettes), len(loci_sorted)):
            base_rgb = cmap((i - len(_LOCUS_PALETTES)) % 10)[:3]
            palettes.append(_generate_palette(base_rgb))
 
    return {locus: palettes[i] for i, locus in enumerate(loci_sorted)}
 
 
def _gene_color(label, gene_type, locus_type, locus_color_map):
    palette = locus_color_map[locus_type]
    if gene_type == "D":
        return palette["D"]
    if gene_type == "J":
        return palette["J"]
    # V gene
    return palette["V_with"] if label else palette["V_without"]
 
 
def _render(contig_size, genes_subset, title, figsize, out_path, padding_bp, include_d):
 
    box_color = "#F6A1A5"
 
    # ── determine locus color mapping ────────────────────────────────────────
    locus_types_present = {g[6] for g in genes_subset}
    locus_color_map = _build_locus_color_map(locus_types_present)
 
    tick_height = 0.35
 
    def _tick_h(productive):
        return tick_height if productive else tick_height * 0.5
 
    zoom_genes = [g for g in genes_subset if g[4]]
 
    zoom_locs  = [g[0] for g in zoom_genes]
    region_min = min(zoom_locs)
    region_max = max(zoom_locs)
    zoom_min   = max(0,           region_min - padding_bp)
    zoom_max   = min(contig_size, region_max + padding_bp)
    locus_size = zoom_max - zoom_min
 
    fig, (ax_ov, ax_zm) = plt.subplots(
        2, 1, figsize=figsize,
        gridspec_kw={"height_ratios": [1, 2.2], "hspace": 0.55}
    )
    fig.suptitle(title, fontsize=13, fontweight="bold", y=0.99)
 
    def _draw_backbone(ax, x0, x1):
        ax.hlines(0, x0, x1, colors="black", linewidths=3, zorder=1)
        ax.vlines([x0, x1], -0.15, 0.15, colors="black", linewidths=2, zorder=2)
 
    _draw_backbone(ax_ov, 0, contig_size)
 
    box_h = 0.38
    ax_ov.add_patch(FancyBboxPatch(
        (zoom_min, -box_h), zoom_max - zoom_min, 2 * box_h,
        boxstyle="square,pad=0",
        linewidth=2, edgecolor=box_color, facecolor=box_color, alpha=0.25, zorder=2
    ))
 
    for location, label, strand, productive, in_zoom, gene_type, locus_type in genes_subset:
        color = _gene_color(label, gene_type, locus_type, locus_color_map)
        alpha = 0.8 if in_zoom else 0.35
        th    = _tick_h(productive) * 0.7
        y_min, y_max = (0, th) if strand == "+" else (-th, 0)
        ax_ov.vlines(location, y_min, y_max, colors=color, linewidths=1.4, zorder=3, alpha=alpha)
 
    ax_ov.set_xlim(-contig_size * 0.02, contig_size * 1.02)
    ax_ov.set_ylim(-0.65, 0.75)
    ax_ov.set_yticks([])
    ax_ov.set_xlabel("Position (bp)", fontsize=9)
    ax_ov.set_title("")
    for spine in ax_ov.spines.values():
        spine.set_visible(False)
 
    zm_pad     = (zoom_max - zoom_min) * 0.01
    zm_xlim_lo = zoom_min - zm_pad
    zm_xlim_hi = zoom_max + zm_pad
    zm_ylim    = (-0.95, 0.95)
 
    _draw_backbone(ax_zm, zoom_min, zoom_max)
 
    for location, label, strand, productive, in_zoom, gene_type, locus_type in zoom_genes:
        color = _gene_color(label, gene_type, locus_type, locus_color_map)
        th    = _tick_h(productive)
        y_min, y_max = (0, th) if strand == "+" else (-th, 0)
        ax_zm.vlines(location, y_min, y_max, colors=color, linewidths=1.8, zorder=3, alpha=0.85)
 
    ax_zm.set_xlim(zm_xlim_lo, zm_xlim_hi)
    ax_zm.set_ylim(*zm_ylim)
    ax_zm.set_yticks([])
    ax_zm.set_title("")
    ax_zm.set_xlabel("")
    ax_zm.tick_params(axis="x", direction="in", pad=-12, labelsize=8)
 
    for spine in ax_zm.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor(box_color)
        spine.set_linewidth(2)
 
    ax_zm.text(0.5, 1.0, f"{locus_size:,} bp",
               transform=ax_zm.transAxes, ha="center", va="bottom",
               fontsize=9, color=box_color, fontweight="bold", clip_on=False)
    ax_zm.text(zoom_min, -0.20, f"{int(zoom_min):,}", ha="center", va="top", fontsize=8, color="black")
    ax_zm.text(zoom_max, -0.20, f"{int(zoom_max):,}", ha="center", va="top", fontsize=8, color="black")
 
    # ── inset legend ─────────────────────────────────────────────────────────
    inset_color  = "black"
    x_span       = zm_xlim_hi - zm_xlim_lo
    x_right      = zm_xlim_hi - x_span * 0.045
    x_left       = x_right    - x_span * 0.075
    y_prod_ins   = zm_ylim[1] - 0.07
    y_unprod_ins = y_prod_ins - tick_height * 0.5
    y_base_ins   = y_prod_ins - tick_height
    label_y_ins  = y_base_ins - 0.03
 
    char_w        = x_span * 0.006
    left_extent   = x_left  - char_w * len("Productive")   / 2
    right_extent  = x_right + char_w * len("Unproductive") / 2
    text_depth    = 0.09
    pad_x         = x_span * 0.01
    pad_y         = 0.02
    inset_box_x0  = left_extent  - pad_x
    inset_box_x1  = right_extent + pad_x
 
    # ── spread labels ─────────────────────────────────────────────────────────
    char_width_bp = locus_size / 80.0
 
    for strand_side in ["+", "-"]:
        labelled = [(loc, lbl, prod, gtype, ltype)
                    for loc, lbl, strand, prod, in_zoom, gtype, ltype in zoom_genes
                    if lbl and strand == strand_side and not include_d]
        if not labelled:
            continue
 
        labelled.sort(key=lambda x: x[0])
        raw_x   = [x[0] for x in labelled]
        min_gap = max(char_width_bp * 5, locus_size * 0.04)
        spread_x = list(_spread_labels(raw_x, min_gap))
 
        if strand_side == "+":
            for i, lx in enumerate(spread_x):
                if lx > inset_box_x0:
                    spread_x[i] = inset_box_x0 - min_gap * 0.5
 
        gene_label_y = (tick_height + 0.22) if strand_side == "+" else -(tick_height + 0.22)
        va           = "bottom"              if strand_side == "+" else "top"
 
        for (orig_x, lbl, prod, gtype, ltype), lx in zip(labelled, spread_x):
            color = _gene_color(lbl, gtype, ltype, locus_color_map)
            ax_zm.annotate(
                "", xy=(orig_x, _tick_h(prod) if strand_side == "+" else -_tick_h(prod)),
                xytext=(lx, gene_label_y),
                arrowprops=dict(arrowstyle="-", color=color, lw=0.8, connectionstyle="arc3,rad=0.0"),
                zorder=4
            )
            ax_zm.text(lx, gene_label_y, lbl, ha="center", va=va, fontsize=7,
                       color=color, rotation=0, zorder=5)
 
    # draw inset box and contents
    ax_zm.add_patch(FancyBboxPatch(
        (inset_box_x0, label_y_ins - text_depth - pad_y),
        (inset_box_x1 - inset_box_x0),
        (zm_ylim[1] - pad_y) - (label_y_ins - text_depth - pad_y),
        boxstyle="round,pad=0.005",
        linewidth=0.8, edgecolor="#aaaaaa", facecolor="white", alpha=1.0, zorder=5
    ))
    ax_zm.vlines(x_left,  y_base_ins, y_prod_ins,   colors=inset_color, linewidths=2.0, zorder=6)
    ax_zm.vlines(x_right, y_base_ins, y_unprod_ins, colors=inset_color, linewidths=2.0, zorder=6)
    ax_zm.text(x_left,  label_y_ins, "Productive",   ha="center", va="top", fontsize=6.5, color=inset_color, zorder=7)
    ax_zm.text(x_right, label_y_ins, "Unproductive", ha="center", va="top", fontsize=6.5, color=inset_color, zorder=7)
 
    # ── legend ────────────────────────────────────────────────────────────────
    # Entries are labeled locus_type + gene_type (e.g. "IGHV", "IGHD", "IGHJ").
    # V genes keep the "with RSS" / "without RSS" distinction as a suffix,
    # but only when both variants are actually present for that locus.
    legend_handles = []
    loci_sorted = sorted(locus_types_present)
 
    for locus in loci_sorted:
        palette = locus_color_map[locus]
        locus_genes = [g for g in genes_subset if g[6] == locus]
        gene_types_here = {g[5] for g in locus_genes}
 
        if "V" in gene_types_here:
            has_with    = any(g[5] == "V" and g[1] for g in locus_genes)
            has_without = any(g[5] == "V" and not g[1] for g in locus_genes)
            v_label = f"{locus}V"
            if has_with:
                suffix = " (with RSS)" if has_without else ""
                legend_handles.append(mpatches.Patch(color=palette["V_with"], label=f"{v_label}{suffix}"))
            if has_without:
                suffix = " (without RSS)" if has_with else ""
                legend_handles.append(mpatches.Patch(color=palette["V_without"], label=f"{v_label}{suffix}"))
        if include_d and "D" in gene_types_here:
            legend_handles.append(mpatches.Patch(color=palette["D"], label=f"{locus}D"))
        if "J" in gene_types_here:
            legend_handles.append(mpatches.Patch(color=palette["J"], label=f"{locus}J"))
 
    ncol = min(len(legend_handles), 4) if legend_handles else 1
    fig.legend(handles=legend_handles, loc="lower center", ncol=ncol,
               frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.01))
 
    # ── connector lines ───────────────────────────────────────────────────────
    fig.canvas.draw()
    for ov_x, zm_xdata in [(zoom_min, zm_xlim_lo), (zoom_max, zm_xlim_hi)]:
        fig.add_artist(ConnectionPatch(
            xyA=(ov_x, -box_h), coordsA=ax_ov.transData,
            xyB=(zm_xdata, zm_ylim[1]), coordsB=ax_zm.transData,
            color=box_color, linewidth=1.6, linestyle="--",
            zorder=0, clip_on=False, arrowstyle="-"
        ))
 
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Figure saved to {out_path}")
    plt.close(fig)
    return fig, (ax_ov, ax_zm)
 
 
def plot_gene_positions(contig_size, genes, title, figsize, output, padding_bp,
                        d_gene_mode="with"):
    """
    Parameters
    ----------
    contig_size : int
    genes : list of [location, label, strand, productive, in_zoom, gene_type, locus_type]
        gene_type  : "V", "D", or "J"
        locus_type : e.g. "IGH", "IGK", "IGL", "TRA", "TRB", ...
                     Genes are colored by (locus_type, gene_type). If the
                     input data contains only a single locus_type, coloring
                     is identical to the original fixed color scheme. If
                     more than one locus_type is present, each additional
                     locus gets its own distinct color palette so loci are
                     visually distinguishable. Legend entries are labeled
                     locus_type + gene_type (e.g. "IGHV", "IGHD", "IGHJ").
    title, figsize, output, padding_bp : as before
    d_gene_mode : str
        "with"    - produce one plot including D genes  (default)
        "without" - produce one plot excluding D genes
        "both"    - produce two plots; output path gets '_with_D' / '_without_D' suffixes
    J genes are always included regardless of d_gene_mode.
    """
    if d_gene_mode not in ("with", "without", "both"):
        raise ValueError("d_gene_mode must be 'with', 'without', or 'both'")
 
    def _suffixed(path, suffix):
        base, ext = os.path.splitext(path)
        return f"{base}{suffix}{ext}"
 
    if d_gene_mode == "with":
        return _render(contig_size, genes, title, figsize, output, padding_bp, include_d=True)
 
    elif d_gene_mode == "without":
        genes_no_d = [g for g in genes if g[5] != "D"]
        return _render(contig_size, genes_no_d, title, figsize, output, padding_bp, include_d=False)
 
    else:  # both
        out_with    = _suffixed(output, "_with_D")
        out_without = _suffixed(output, "_without_D")
        genes_no_d  = [g for g in genes if g[5] != "D"]
        r1 = _render(contig_size, genes,      title, figsize, out_with,    padding_bp, include_d=True)
        r2 = _render(contig_size, genes_no_d, title, figsize, out_without, padding_bp, include_d=False)
        return r1, r2


def plot_contig(genome,contig,rss,v_genes,d_genes,j_genes,out):
    subprocess.run(["samtools","faidx",genome])
    with open(genome+".fai","r") as read:
        reader=csv.reader(read,delimiter="\t")
        for row in reader:
            if row[0]==contig:
                contig_size=int(row[1])
                break
        read.close()
    gene_positions=[]
    loci=[]
    for g in v_genes:
        if str(g[5])=="True":
            bo=True
        elif str(g[5])=="False":
            bo=False
        gene_positions.append([int(g[2]),False,str(g[3]),bo,True,"V",g[-1]])
        if g[-1] not in loci:
            loci.append(g[-1])
    for r in rss:
        for g in gene_positions:
            if int(r[3])==g[0] and r[7]==g[-1]:
                g[1]=str(r[9])

    if "IGH" in loci or "TRD" in loci or "TRB" in loci:
        d_gene_mode="both"
    else:
        d_gene_mode="without"

    if d_gene_mode=="both" or d_gene_mode=="with":
        temp_d=[]
        for d in d_genes:
            if str(d[6])=="True":
                bo=True
            elif str(d[6])=="False":
                bo=False
            else:
                bo=True
            temp_d.append([int(d[3]),False,d[4],bo,True,"D",d[7]])
            if d[7] not in loci:
                loci.append(d[7])
        temp_j=[]
        for j in j_genes:
            if str(j[6])=="True":
                bo=True
            elif str(j[6])=="False":
                bo=False
            else:
                bo=True
            temp_j.append([int(j[3]),False,j[4],bo,True,"J",j[7]])
            if j[7] not in loci:
                loci.append(j[7])
        for d in temp_d:
            gene_positions.append(d)
        for j in temp_j:
            gene_positions.append(j) 
    
    sorted_gene_postions = sorted(gene_positions, key=lambda x: x[0])
    locus_size=int(sorted_gene_postions[-1][0])-int(sorted_gene_postions[0][0])
    if locus_size!=0:
        padding=int(locus_size*0.1)
    else:
        padding=1000
    
    locus_types=[]
    plot_gene_positions(contig_size,sorted_gene_postions,out.split("/")[-1].replace("_plot_v2.svg",""),(16, 4),out,padding,d_gene_mode=d_gene_mode)


all_d_genes=[]
for d_file in os.listdir("#d_gene_scripts/final_bird_D"):
    with open("#d_gene_scripts/final_bird_D/"+d_file+"/bird_d_genes.csv","r") as d_read:
        reader=csv.reader(d_read)
        header=next(reader)
        for row in reader:
            all_d_genes.append(row)
        d_read.close()

all_j_genes=[]
for j_file in os.listdir("#j_gene_scripts/final_bird_J"):
    with open("#j_gene_scripts/final_bird_J/"+j_file+"/bird_j_genes.csv","r") as j_read:
        reader=csv.reader(j_read)
        header=next(reader)
        for row in reader:
            all_j_genes.append(row)
        j_read.close()

for f in os.listdir(data_dir):
    if f in order_names:
        if order_filter=="None" or f==order_filter:
            for f1 in os.listdir(data_dir+"/"+f):
                if os.path.isdir(data_dir+"/"+f+"/"+f1):
                    if species_filter=="None" or f1==species_filter:
                        for f2 in os.listdir(data_dir+"/"+f+"/"+f1):
                            if os.path.isdir(data_dir+"/"+f+"/"+f1+"/"+f2):
                                if haplotype_filter=="None" or f2==haplotype_filter:
                                    temp_rss=[]
                                    for r in rss_list:
                                        if (f+"/"+f1+"/"+f2) in r[0]:
                                            temp_rss.append(r)
                                    for g in genome_paths:
                                        if f2 in g[0]:
                                            genome = g[1]
                                    contigs=[]
                                    for c in contig_list:
                                        if f1+"/"+f2 in c[1]:
                                            contigs.append(c)

                                    temp_genes=[]
                                    for f3 in os.listdir(data_dir+"/"+f+"/"+f1+"/"+f2):
                                        if f3.startswith("combined_genes") and "clean" not in f3:
                                            with open(data_dir+"/"+f+"/"+f1+"/"+f2+"/"+f3,"r") as read_gene:
                                                reader=csv.reader(read_gene,delimiter="\t")
                                                header=next(reader)
                                                temp_locus_genes=[]
                                                for row in reader:
                                                    temp_locus_genes.append(row)
                                                if len(temp_locus_genes)>1:
                                                    for te in temp_locus_genes:
                                                        if is_low_complexity(te[4])==False and len(te[4])>=250:
                                                            temp_genes.append(te)
                                                else:
                                                    for te in temp_locus_genes:
                                                        temp_genes.append(te)
                                                read_gene.close()
                                    unique_contigs=[]
                                    for temp_con in contigs:
                                        if temp_con[0] not in unique_contigs:
                                            unique_contigs.append(temp_con[0])
                                    
                                    for un in unique_contigs:
                                        if contig_filter!="-none" and un!=contig_filter.split(":")[1]:
                                            continue
                                        cont=[]
                                        for contig in contigs:
                                            if contig[0]==un:
                                                cont.append(contig)
                                        
                                        current_genes=[]
                                        current_rss=[]
                                        d_genes=[]
                                        j_genes=[]
                                        loci=""
                                        for con in cont:                            
                                            if loci=="":
                                                loci=str(con[2])
                                            else:
                                                loci=loci+"-"+str(con[2])

                                            for temp in temp_genes:
                                                if con[0]==temp[1] and con[2]==temp[-1]:
                                                    current_genes.append(temp)

                                            for trss in temp_rss:
                                                if con[0]==trss[2] and con[2]==trss[7]:
                                                    current_rss.append(trss)

                                            for d in all_d_genes:
                                                if (f+"/"+f1+"/"+f2) == d[0] and con[0] == d[2]:
                                                    d_genes.append(d)
                                            
                                            for j in all_j_genes:
                                                if f+"/"+f1+"/"+f2 == j[0] and con[0] == j[2]:
                                                    j_genes.append(j)

                                        output=data_dir+"/"+f+"/"+f1+"/"+f2+"/"+un+"_"+loci+"_plot_v2.svg"
                                        if os.path.isfile(output) or os.path.isfile(output.replace(".svg","_with_D.svg")) or os.path.isfile(output.replace(".svg","_without_D.svg")):
                                            if delete=="-n":
                                                continue
                                            elif delete=="-d":
                                                if os.path.isfile(output):
                                                    os.remove(output)
                                                if os.path.isfile(output.replace(".svg","_with_D.svg")):
                                                    os.remove(output.replace(".svg","_with_D.svg"))
                                                if os.path.isfile(output.replace(".svg","_without_D.svg")):
                                                    os.remove(output.replace(".svg","_without_D.svg"))
                                        print("\ngenome ",genome)
                                        print("path ",data_dir+"/"+f+"/"+f1+"/"+f2)
                                        print("contig ",un)
                                        #print(current_genes)
                                        #print(current_rss)
                                        try:
                                            plot_contig(genome,un,current_rss,current_genes,d_genes,j_genes,output)
                                        except:
                                            print("Failure on "+f+"/"+f1+"/"+f2+" : "+un)
                                    