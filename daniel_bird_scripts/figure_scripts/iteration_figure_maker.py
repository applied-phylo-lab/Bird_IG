import os
import csv
import sys
import copy
from collections import Counter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import numpy as np

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

# ---------------------------------------------------------------------------
# Typography & colour constants
# ---------------------------------------------------------------------------
TITLE_SIZE     = 10
ROW_LABEL_SZ   = 9.5
AXIS_LABEL_SZ  = 8.5
TICK_SZ        = 8
XTICKLABEL_SZ  = 7.5
GROUP_LABEL_SZ = 12
HAPLO_LABEL_SZ = 7.5
LEGEND_SZ      = 7
 
PALETTE = ["#172869FF","#EA7580FF"]

GROUP_BG  = ['#EEF3FB', '#FBF0EE']
GRID_COL  = '#E8E8E8'
SPINE_COL = "#00000000"
SPINE_LW  = 0.9
PRODUCTIVE_LIGHTEN = 0.52
 
 
def _lighten(hex_col, amount=0.5):
    hex_col = hex_col.lstrip('#')
    r, g, b = [int(hex_col[i:i+2], 16) for i in (0, 2, 4)]
    return (f'#{int(r+(255-r)*amount):02X}'
            f'{int(g+(255-g)*amount):02X}'
            f'{int(b+(255-b)*amount):02X}')
 
 
def _fmt_large(x, pos=None):
    if x >= 1_000_000: return f'{x/1_000_000:.1f}M'
    if x >= 10_000:    return f'{x/1_000:.0f}k'
    if x >= 1_000:     return f'{x/1_000:.1f}k'
    return f'{int(x)}'
 
 
def make_gene_figure(data, output_path):
    """
    Three-row gene-statistics figure.
 
    Parameters
    ----------
    data : list[list]
        data[0]  IG iterations
        data[1]  TCR iterations
        Each iteration: [name, [[gene, total, productive, contigs],...],
                                [[gene, haplo_with, haplo_total],...]]
    output_path : str – destination SVG path.
    """
 
    # -----------------------------------------------------------------------
    # 1. Parse
    # -----------------------------------------------------------------------
    def parse_group(iterations):
        gene_set, parsed = [], []
        for entry in iterations:
            name  = entry[0]
            genes = {r[0]: (r[1], r[2], r[3]) for r in entry[1]}
            haplo = {r[0]: (r[1], r[2])        for r in entry[2]}
            for g in genes:
                if g not in gene_set:
                    gene_set.append(g)
            parsed.append({'name': name, 'genes': genes, 'haplo': haplo})
        return gene_set, parsed
 
    groups       = [parse_group(grp) for grp in data]
    group_labels = ['IG', 'TCR']
    n_groups     = len(groups)
 
    gene_colours = {"IGH":"#ABCBE7","IGL":"#92B09A","TRA":"#F28B87","TRB":"#F6A1A5","TRG":"#F7B7A0","TRD":"#F8CD9C"}
 
    # -----------------------------------------------------------------------
    # 2. Compute minimum total haplotypes per gene type per group
    # -----------------------------------------------------------------------
    group_haplo_info = []
    for gi, (gene_types, iters) in enumerate(groups):
        info = {}
        for g in gene_types:
            totals = [it['haplo'].get(g, (0, 0))[1] for it in iters
                      if it['haplo'].get(g, (0, 0))[1] > 0]
            info[g] = min(totals) if totals else 0
        group_haplo_info.append(info)
 
    # -----------------------------------------------------------------------
    # 3. Compute shared y-axis scales per row per group
    #    Row 0: max of total genes across all gene types in the group
    #    Row 1: max of contigs across all gene types in the group
    #    Row 2: always 0–108 (percentage), same for all
    # -----------------------------------------------------------------------
    group_row_ylims = []   # list of (ylim_row0, ylim_row1) per group
    for gi, (gene_types, iters) in enumerate(groups):
        max_total = 0
        max_contig = 0
        for g in gene_types:
            for it in iters:
                t, p, c = it['genes'].get(g, (0, 0, 0))
                max_total  = max(max_total,  t)
                max_contig = max(max_contig, c)
        group_row_ylims.append((max_total * 1.22, max_contig * 1.18))
 
    # -----------------------------------------------------------------------
    # 4. Layout constants (all in inches)
    # -----------------------------------------------------------------------
    N_ROWS   = 3
    ROW_H    = 2.8
    COL_W    = 2.9
 
    LEFT_PAD  = 1.45
    RIGHT_PAD = 0.50
    TOP_PAD   = 1.55
    BOT_PAD   = 1.20
 
    COL_GAP   = 0.55
    GROUP_GAP = 1.10
    ROW_GAP   = 0.85
 
    # IG box: both left and right sides use SIDE_PAD (= measured ylabel coverage)
    # TCR box: left = SIDE_PAD, right = SIDE_PAD  (unchanged from user's version)
    SIDE_PAD    = 0.79
    WHITE_SPACE = 0.28   # preserved gap between the two boxes
    # IG right pad = SIDE_PAD (same as left, per request)
    # White space between boxes = GROUP_GAP - IG_right - TCR_left
    #                           = 1.10 - 0.79 - 0.79 = -0.48  → boxes would overlap
    # So we widen GROUP_GAP to guarantee WHITE_SPACE:
    GROUP_GAP = SIDE_PAD + WHITE_SPACE + SIDE_PAD   # 0.79+0.28+0.79 = 1.86"
 
    PAD_BOT  = 0.95
    PAD_TOP  = 0.15
    BANNER_H = 0.80
 
    total_cols = sum(len(gt) for gt, _ in groups)
    n_col_gaps = total_cols - 1
    inner_w = (total_cols * COL_W
               + (n_col_gaps - (n_groups - 1)) * COL_GAP
               + (n_groups - 1) * GROUP_GAP)
 
    fig_w = LEFT_PAD + inner_w + RIGHT_PAD
    fig_h = TOP_PAD  + N_ROWS * ROW_H + (N_ROWS - 1) * ROW_GAP + BOT_PAD
 
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=300)
    fig.patch.set_facecolor('white')
 
    # -----------------------------------------------------------------------
    # 5. Column x-positions (inches from left)
    # -----------------------------------------------------------------------
    col_x = []
    x_cursor = LEFT_PAD
    for gi, (gene_types, _) in enumerate(groups):
        for li in range(len(gene_types)):
            col_x.append(x_cursor)
            is_last = li == len(gene_types) - 1
            gap = GROUP_GAP if (is_last and gi < n_groups - 1) else COL_GAP
            x_cursor += COL_W + gap
 
    def row_bottom_abs(ri):
        return BOT_PAD + (N_ROWS - 1 - ri) * (ROW_H + ROW_GAP)
 
    def fx(x_in): return x_in / fig_w
    def fy(y_in): return y_in / fig_h
 
    # -----------------------------------------------------------------------
    # 6. Create axes
    # -----------------------------------------------------------------------
    axes = {}
    global_col = 0
    for gi, (gene_types, _) in enumerate(groups):
        for li in range(len(gene_types)):
            for ri in range(N_ROWS):
                rect = [fx(col_x[global_col]), fy(row_bottom_abs(ri)),
                        fx(COL_W), fy(ROW_H)]
                ax = fig.add_axes(rect)
                axes[(gi, li, ri)] = ax
            global_col += 1
 
    # -----------------------------------------------------------------------
    # 7. Axis styling
    # -----------------------------------------------------------------------
    def style_ax(ax):
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
        for spine in ['left', 'bottom']:
            ax.spines[spine].set_color(SPINE_COL)
            ax.spines[spine].set_linewidth(SPINE_LW)
        ax.tick_params(axis='both', which='major',
                       labelsize=TICK_SZ, length=4, width=0.8,
                       color=SPINE_COL, pad=3)
        ax.yaxis.grid(True, color=GRID_COL, linewidth=0.7, zorder=0)
        ax.set_axisbelow(True)
 
    def set_y_formatter(ax):
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(_fmt_large))
        ax.yaxis.set_major_locator(
            ticker.MaxNLocator(nbins=4, integer=True, min_n_ticks=3))
 
    # -----------------------------------------------------------------------
    # 8. Plot rows
    # -----------------------------------------------------------------------
    BAR_W = 0.55
 
    for gi, (gene_types, iters) in enumerate(groups):
        n_iters = len(iters)
        iter_labels = [f'Iteration {i+1}' for i in range(n_iters)]
        ylim_row0, ylim_row1 = group_row_ylims[gi]
 
        for li, g in enumerate(gene_types):
            is_leftmost = (li == 0)
            col   = gene_colours[g]
            col_l = _lighten(col, PRODUCTIVE_LIGHTEN)
            xs    = np.arange(n_iters)
 
            # ROW 0 – productive (bottom, darker) + unproductive (top, lighter)
            ax0 = axes[(gi, li, 0)]
            totals      = [it['genes'].get(g, (0,0,0))[0] for it in iters]
            productives = [it['genes'].get(g, (0,0,0))[1] for it in iters]
            unprod      = [t - p for t, p in zip(totals, productives)]
 
            ax0.bar(xs, productives, color=col,   width=BAR_W, zorder=2, linewidth=0)
            ax0.bar(xs, unprod, bottom=productives, color=col_l, width=BAR_W,
                    zorder=2, linewidth=0)
 
            ax0.set_ylim(0, ylim_row0)
            set_y_formatter(ax0)
            if is_leftmost:
                ax0.set_ylabel('# of Genes', fontsize=AXIS_LABEL_SZ, labelpad=5)
            ax0.set_xticks(xs)
            ax0.set_xticklabels(iter_labels, rotation=38, ha='right',
                                fontsize=XTICKLABEL_SZ)
            style_ax(ax0)
 
            prod_patch   = mpatches.Patch(color=col,   label='Functional')
            unprod_patch = mpatches.Patch(color=col_l, label='Pseudo')
            ax0.legend(handles=[prod_patch, unprod_patch],
                       fontsize=LEGEND_SZ, loc='upper left',
                       frameon=True, framealpha=0.85,
                       edgecolor='#CCCCCC', borderpad=0.4,
                       handlelength=1.0, handleheight=0.8,
                       labelspacing=0.3, handletextpad=0.4)
 
            # ROW 1 – contigs
            ax1 = axes[(gi, li, 1)]
            contigs = [it['genes'].get(g, (0,0,0))[2] for it in iters]
            ax1.bar(xs, contigs, color=col, width=BAR_W, zorder=2, linewidth=0)
            ax1.set_ylim(0, ylim_row1)
            set_y_formatter(ax1)
            if is_leftmost:
                ax1.set_ylabel('# of Contigs', fontsize=AXIS_LABEL_SZ, labelpad=5)
            ax1.set_xticks(xs)
            ax1.set_xticklabels(iter_labels, rotation=38, ha='right',
                                fontsize=XTICKLABEL_SZ)
            style_ax(ax1)
 
            # ROW 2 – % haplotypes with genes
            ax2 = axes[(gi, li, 2)]
            pcts = []
            for it in iters:
                hw, ht = it['haplo'].get(g, (0, 1))
                pcts.append(100 * hw / ht if ht else 0)
            ax2.bar(xs, pcts, color=col, width=BAR_W, zorder=2, linewidth=0)
            ax2.set_ylim(0, 108)
            ax2.yaxis.set_major_locator(ticker.MultipleLocator(25))
            ax2.yaxis.set_major_formatter(
                ticker.FuncFormatter(lambda v, _: f'{v:.0f}%'))
            if is_leftmost:
                ax2.set_ylabel(f"% of total Haplotypes with genes",
                               fontsize=AXIS_LABEL_SZ, labelpad=5)
            ax2.set_xticks(xs)
            ax2.set_xticklabels(iter_labels, rotation=38, ha='right',
                                fontsize=XTICKLABEL_SZ)
            style_ax(ax2)
 
    # -----------------------------------------------------------------------
    # 9. Group background boxes + banners
    # -----------------------------------------------------------------------
    global_col = 0
    for gi, (gene_types, iters) in enumerate(groups):
        nc = len(gene_types)
        if nc == 0:
            global_col += nc
            continue
 
        x0_abs    = col_x[global_col]
        x1_abs    = col_x[global_col + nc - 1] + COL_W
        y_bot_abs = row_bottom_abs(N_ROWS - 1)
        y_top_abs = row_bottom_abs(0) + ROW_H
 
        # IG:  left = SIDE_PAD, right = SIDE_PAD  (symmetric, per request)
        # TCR: left = SIDE_PAD, right = SIDE_PAD  (unchanged from user's version)
        pad_l = SIDE_PAD
        pad_r = SIDE_PAD
 
        rx = fx(x0_abs - pad_l)
        ry = fy(y_bot_abs - PAD_BOT)
        rw = fx(x1_abs - x0_abs + pad_l + pad_r)
        rh = fy(y_top_abs - y_bot_abs + PAD_BOT + PAD_TOP)
 
        fig.add_artist(mpatches.FancyBboxPatch(
            (rx, ry), rw, rh,
            boxstyle='round,pad=0.005',
            linewidth=1.1, edgecolor='#C4C4C4',
            facecolor=GROUP_BG[gi % len(GROUP_BG)],
            zorder=0, transform=fig.transFigure, clip_on=False,
        ))
 
        bx = fx(x0_abs - pad_l)
        by = fy(y_top_abs + PAD_TOP)
        bw = fx(x1_abs - x0_abs + pad_l + pad_r)
        bh = fy(BANNER_H)
 
        fig.add_artist(mpatches.FancyBboxPatch(
            (bx, by), bw, bh,
            boxstyle='round,pad=0.005',
            linewidth=0, facecolor=PALETTE[gi % len(PALETTE)],
            zorder=1, transform=fig.transFigure, clip_on=False,
        ))
 
        banner_bot = y_top_abs + PAD_TOP
 
        haplo_info     = group_haplo_info[gi]
        gene_types_gi  = groups[gi][0]
        all_haplo_vals = [haplo_info.get(g, 0) for g in gene_types_gi
                          if haplo_info.get(g, 0) > 0]
        group_n = min(all_haplo_vals) if all_haplo_vals else 0
 
        group_y = banner_bot + BANNER_H * 0.65
        if group_labels[gi] == 'IG':
            group_str = f'{group_labels[gi]}\n(From VGP  |  n = {group_n:,})'
        else:
            group_str = f'{group_labels[gi]}\n(VGP + Additional Data  |  n = {group_n:,})'
 
        fig.text(fx((x0_abs + x1_abs) / 2), fy(group_y),
                 group_str,
                 ha='center', va='center',
                 fontsize=GROUP_LABEL_SZ, fontweight='bold',
                 color='white', transform=fig.transFigure, zorder=2)
 
        for li, g in enumerate(gene_types_gi):
            col_cx = col_x[global_col + li] + COL_W / 2
            gene_y = banner_bot + BANNER_H * 0.22
            fig.text(fx(col_cx), fy(gene_y), g,
                     ha='center', va='center',
                     fontsize=TITLE_SIZE, fontweight='bold',
                     color='white', transform=fig.transFigure, zorder=2)
 
        global_col += nc
 
    # -----------------------------------------------------------------------
    # 10. Row section labels (vertical, far left)
    # -----------------------------------------------------------------------
    row_section_labels = [
        'Total & Productive Genes',
        'Total Contigs',
        '% Haplotypes with Genes',
    ]
    SECTION_X_ABS = 0.25
    for ri, label in enumerate(row_section_labels):
        fig.text(fx(SECTION_X_ABS), fy(row_bottom_abs(ri) + ROW_H / 2),
                 label, ha='center', va='center',
                 fontsize=ROW_LABEL_SZ, fontweight='semibold',
                 color='#3A3A3A', rotation=90,
                 transform=fig.transFigure)
 
    # -----------------------------------------------------------------------
    # 11. Save
    # -----------------------------------------------------------------------
    fig.savefig(output_path, format='svg', bbox_inches='tight',
                dpi=300, facecolor='white')
    plt.close(fig)
    print(f'Saved → {output_path}')

loci_folders=[]
argn=0
for arg in sys.argv:
    if arg.startswith("-l:"):
        locus_data=[arg.replace("-l:",""),[]]
        for a in sys.argv[argn+1].split(","):
            locus_data[1].append(a)
        loci_folders.append(locus_data)
    argn+=1

order_names=["Cormorants","Cranes","Doves","Eagles","Falcons","Hummingbirds","Ibises","Landfowl","MiscBirds","Owls","Parrots","Plovers","Songbirds","Suboscines","Waterfowl","Woodpeckers"] #names of all bird order folders
tcr_iteration_data=[]
ig_iteration_data=[]
for gt in loci_folders:
    if gt[0]=="IG":
        initialize=[["IGH",0,0,0],["IGL",0,0,0]]
    elif gt[0]=="TCR":
        initialize=[["TRA",0,0,0],["TRB",0,0,0],["TRG",0,0,0],["TRD",0,0,0]]
    for folder in gt[1]:
        haps=[]
        species_list=[]
        orders=[]
        total_genes=copy.deepcopy(initialize)
        for order in os.listdir(folder):
            if os.path.isdir(folder+"/"+order) and order in order_names:
                order_genes=copy.deepcopy(initialize)
                for species in os.listdir(folder+"/"+order):
                    if os.path.isdir(folder+"/"+order+"/"+species):
                        species_genes=copy.deepcopy(initialize)
                        for hap in os.listdir(folder+"/"+order+"/"+species):
                            if os.path.isdir(folder+"/"+order+"/"+species+"/"+hap):
                                genes=copy.deepcopy(initialize)
                                for data in os.listdir(folder+"/"+order+"/"+species+"/"+hap):
                                    file=folder+"/"+order+"/"+species+"/"+hap+"/"+data
                                    if data.startswith("combined") and data.endswith(".txt"):
                                        with open(file,"r") as read:
                                            reader=csv.reader(read,delimiter="\t")
                                            header=next(reader)
                                            total_gene_count=0
                                            prod_gene_count=0
                                            contigs=[]
                                            for row in reader:
                                                if len(row[4])>=250 and is_low_complexity(row[4])==False:
                                                    total_gene_count+=1
                                                    if row[1] not in contigs:
                                                        contigs.append(row[1])
                                                    if str(row[5])=="True":
                                                        prod_gene_count+=1
                                            for g in genes:
                                                if str(data.replace("combined_genes_","").replace(".txt",""))==g[0]:
                                                    g[1]=total_gene_count
                                                    g[2]=prod_gene_count
                                                    g[3]=len(contigs)
                                haps.append([order+"/"+species+"/"+hap])
                                #print("\n\nHap")
                                #print(genes)
                                for ge in genes:
                                    haps[-1].append(ge)
                                for g in genes:
                                    for sg in species_genes:
                                        if g[0]==sg[0]:
                                            sg[1]=int(sg[1])+int(g[1])
                                            sg[2]=int(sg[2])+int(g[2])
                                            sg[3]=int(sg[3])+int(g[3])
                        species_list.append([order+"/"+species])
                        #print("\n\nSpecies")
                        #print(species_genes)
                        for ge in species_genes:
                            species_list[-1].append(ge)
                        for sg in species_genes:
                            for og in order_genes:
                                if sg[0]==og[0]:
                                    og[1]=int(og[1])+int(sg[1])
                                    og[2]=int(og[2])+int(sg[2])
                                    og[3]=int(og[3])+int(sg[3])
                orders.append([order])
                #print("\n\nOrder")
                #print(order_genes)
                for ge in order_genes:
                    orders[-1].append(ge)
                for og in order_genes:
                    for tg in total_genes:
                        if tg[0]==og[0]:
                            tg[1]=int(tg[1])+int(og[1])
                            tg[2]=int(tg[2])+int(og[2])
                            tg[3]=int(tg[3])+int(og[3])
                #print(total_genes)
        if gt[0]=="IG":
            ig_iteration_data.append([folder.split("/")[-1],total_genes,orders,species_list,haps])
        if gt[0]=="TCR":
            tcr_iteration_data.append([folder.split("/")[-1],total_genes,orders,species_list,haps])

if ig_iteration_data!=[]:
    order_inputs=[]
    total_input=[]
    for i in ig_iteration_data:
        total_hap_percent=[["IGH",0,0],["IGL",0,0]]
        for order in i[2]:
            hap_percent=[["IGH",0,0],["IGL",0,0]]
            for hap in i[-1]:
                if order[0] in hap[0]:
                    for locus in hap[1:]:
                        if locus[0]=="IGH":
                            hap_percent[0][2]+=1
                            if locus[1]!=0:
                                hap_percent[0][1]+=1
                        elif locus[0]=="IGL":
                            hap_percent[1][2]+=1
                            if locus[1]!=0:
                                hap_percent[1][1]+=1
            order.append(hap_percent)
            function_input=[i[0],order[1:3],order[3:]]
            for ord in order_inputs:
                if ord[0]==order[0]:
                    ord[1].append(function_input)
            for t in total_hap_percent:
                for h in hap_percent:
                    if t[0]==h[0]:
                        t[1]=t[1]+h[1]
                        t[2]=t[2]+h[2]
        i.append(total_hap_percent)
        function_input=[i[0],i[1],i[5]]
        total_input.append(function_input)
    ig_total_input=total_input
    ig_order_input=order_inputs

if tcr_iteration_data!=[]:
    order_inputs=[]
    total_input=[]
    for i in tcr_iteration_data:
        total_hap_percent=[["TRA",0,0],["TRB",0,0],["TRG",0,0],["TRD",0,0]]
        for order in i[2]:
            found=False
            for ord in order_inputs:
                if order[0]==ord[0]:
                    found=True
            if found==False:
                order_inputs.append([order[0],[]])
            hap_percent=[["TRA",0,0],["TRB",0,0],["TRG",0,0],["TRD",0,0]]
            for hap in i[-1]:
                if order[0] in hap[0]:
                    for locus in hap[1:]:
                        if locus[0]=="TRA":
                            hap_percent[0][2]+=1
                            if locus[1]!=0:
                                hap_percent[0][1]+=1
                        elif locus[0]=="TRB":
                            hap_percent[1][2]+=1
                            if locus[1]!=0:
                                hap_percent[1][1]+=1
                        elif locus[0]=="TRG":
                            hap_percent[2][2]+=1
                            if locus[1]!=0:
                                hap_percent[2][1]+=1
                        elif locus[0]=="TRD":
                            hap_percent[3][2]+=1
                            if locus[1]!=0:
                                hap_percent[3][1]+=1
            order.append(hap_percent)
            function_input=[i[0],order[1:5],order[5:]]
            for ord in order_inputs:
                if ord[0]==order[0]:
                    ord[1].append(function_input)
            for t in total_hap_percent:
                for h in hap_percent:
                    if t[0]==h[0]:
                        t[1]=t[1]+h[1]
                        t[2]=t[2]+h[2]
        i.append(total_hap_percent)
        function_input=[i[0],i[1],i[5]]
        total_input.append(function_input)

    tcr_total_input=total_input
    tcr_order_input=order_inputs

if os.path.isdir("iteration_figures")==False:
    os.mkdir("iteration_figures")
if tcr_iteration_data!=[] and ig_iteration_data!=[]:
    function_args=[ig_total_input,tcr_total_input]
    print(function_args)
    make_gene_figure(function_args,"iteration_figures/total_iteration_figure.svg")
    for tcr_order in tcr_order_input:
        for ig_order in ig_order_input:
            if ig_order[0]==tcr_order[0]:
                function_args=[ig_order[1],tcr_order[1]]
                make_gene_figure(function_args,"iteration_figures/"+ig_order[0]+"_iteration_figure.svg")