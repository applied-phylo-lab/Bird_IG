import os
import sys
import csv
csv.field_size_limit(10000000000) 
import subprocess
import shutil
from collections import Counter
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

folder=sys.argv[1]
delete=sys.argv[2]
if delete!="-d" and delete!="-n":
    print("enter valid delete argument")
    sys.exit(1)

plot_contigs=sys.argv[3]
if plot_contigs!="-p:k" and plot_contigs!="-p:d" and plot_contigs!="-n" :
    print("enter valid contig plot argument")
    sys.exit(1)

contigs=[]
with open("ig_contig_list.csv","r") as read:
    reader=csv.reader(read)
    header=next(reader)
    for row in reader:
        contigs.append(row)
    read.close()

genomes=[]
with open("bird_genome_paths.csv","r") as read:
    reader=csv.reader(read)
    header=next(reader)
    for row in reader:
        genomes.append(row)
    read.close()

all_d_genes=[]
for d_file in os.listdir("/local/storage/dhardesty/assemblies/#rss_scripts/#d_gene_scripts/final_bird_D"):
    with open("/local/storage/dhardesty/assemblies/#rss_scripts/#d_gene_scripts/final_bird_D/"+d_file+"/bird_d_genes.csv","r") as d_read:
        reader=csv.reader(d_read)
        header=next(reader)
        for row in reader:
            all_d_genes.append(row)
        d_read.close()

all_j_genes=[]
for j_file in os.listdir("/local/storage/dhardesty/assemblies/#rss_scripts/#j_gene_scripts/final_bird_J"):
    with open("/local/storage/dhardesty/assemblies/#rss_scripts/#j_gene_scripts/final_bird_J/"+j_file+"/bird_j_genes.csv","r") as j_read:
        reader=csv.reader(j_read)
        header=next(reader)
        for row in reader:
            all_j_genes.append(row)
        j_read.close()


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

def ig_detective_parse(source,contig,folder):
    order_list=["Cormorants","Cranes","Doves","Eagles","Falcons","Hummingbirds","Ibises","Landfowl","MiscBirds","Owls","Parrots","Plovers","Songbirds","Suboscines","Waterfowl","Woodpeckers"] #names of all bird order folders
    for fold in os.listdir(folder):
        if os.path.isdir(folder+"/"+fold) and fold in order_list:
            for fold1 in os.listdir(folder+"/"+fold):
                if os.path.isdir(folder+"/"+fold+"/"+fold1):
                    for fold2 in os.listdir(folder+"/"+fold+"/"+fold1):
                        if (fold+"/"+fold1+"/"+fold2)==source:
                            positions=[]
                            for file in os.listdir(folder+"/"+fold+"/"+fold1+"/"+fold2):
                                if file.startswith("combined_genes") and "clean" not in file:
                                    with open(folder+"/"+fold+"/"+fold1+"/"+fold2+"/"+file,"r") as read:
                                        reader=csv.reader(read,delimiter="\t")
                                        header=next(reader)
                                        for row in reader:
                                            if row[1]==contig and len(row[4])>=250 and is_low_complexity(row[4])==False:
                                                positions.append([int(row[2]),row[3],row[0],row[5],file.replace("combined_genes_","").replace(".txt","")])
                                        read.close()
                            if positions!=[]:
                                return positions, (folder+"/"+fold+"/"+fold1+"/"+fold2)
                            else:
                                for file in os.listdir(folder+"/"+fold+"/"+fold1+"/"+fold2):
                                    if file.startswith("combined_genes") and 'clean' not in file:
                                        with open(folder+"/"+fold+"/"+fold1+"/"+fold2+"/"+file,"r") as read:
                                            reader=csv.reader(read,delimiter="\t")
                                            header=next(reader)
                                            for row in reader:
                                                if row[1]==contig:
                                                    positions.append([int(row[2]),row[3],row[0],row[5],file.replace("combined_genes_","").replace(".txt","")])
                                        read.close()
                                return positions, (folder+"/"+fold+"/"+fold1+"/"+fold2)

def annotate_gepard_dotplot(positions, locus_start, locus_end, dotplot_path,
                             output_path, dpi=150, thickness_px=80):
    """
    Adds gene-position number lines to the bottom and right of the plot
    box in a Gepard dotplot image, plus a legend in the top right keyed
    by LOCUS_TYPE + GENE_TYPE (e.g. "IGHV", "IGHD", "IGHJ") on row 1,
    Productive/Non-productive on row 2, rendered at exactly the same
    font size as Gepard's own info text, without altering the original
    image contents in any way. Returns the output path.
 
    Everything (box detection, line rendering, legend placement, font
    calibration) lives inside this one function as nested helpers.
    """
 
    # ------------------------------------------------------------------
    # Color scheme, keyed by (locus_type, gene_type). Edit/extend this
    # dict freely -- it is the single source of truth for tick + legend
    # swatch colors. Any (locus_type, gene_type) combo not listed here
    # falls back to a neutral gray (see style_for()).
    # ------------------------------------------------------------------
    LOCUS_GENE_SHADES = {
        ('IGH', 'V'): "#ABCBE7",
        ('IGH', 'D'): "#5198D6",
        ('IGH', 'J'): "#D7E4EF",

        ('IGL', 'V'): "#92B09A",
        ('IGL', 'J'): "#CBE2D1",

        ('TRA', 'V'): "#F28B87",
        ('TRA', 'J'): "#F1D1D0",

        ('TRB', 'V'): "#F6A1A5",
        ('TRB', 'D'): "#E36378",
        ('TRB', 'J'): "#F2D9DD",

        ('TRG', 'V'): "#F7B7A0",
        ('TRG', 'J'): "#F2DED4",

        ('TRD', 'V'): "#F8CD9C",
        ('TRD', 'D'): "#E7AF5A",
        ('TRD', 'J'): "#F2E6D4",
    }
    DEFAULT_SHADE = "#888888"  # fallback for any combo not in the dict above
 
    # Fixed draw order for gene types within a locus:
    # V (variable) -> D (diversity) -> J (joining).
    GENE_TYPE_ORDER = {"V": 0, "D": 1, "J": 2}
 
    # ------------------------------------------------------------------
    # Locate the square plot box inside the Gepard PNG.
    # ------------------------------------------------------------------
    def detect_gepard_box(dark_thresh=128, rel_thresh=0.9):
        img = Image.open(dotplot_path).convert("L")
        arr = np.array(img)
        img_h, img_w = arr.shape
 
        dark = arr < dark_thresh
        row_counts = dark.sum(axis=1)
        col_counts = dark.sum(axis=0)
 
        if row_counts.max() == 0 or col_counts.max() == 0:
            raise ValueError(
                f"Could not find any dark pixels in {dotplot_path}; "
                "is this really a Gepard dotplot image?"
            )
 
        row_idxs = np.where(row_counts >= rel_thresh * row_counts.max())[0]
        col_idxs = np.where(col_counts >= rel_thresh * col_counts.max())[0]
 
        top, bottom = int(row_idxs.min()), int(row_idxs.max())
        left, right = int(col_idxs.min()), int(col_idxs.max())
 
        return {
            "left": left, "right": right, "top": top, "bottom": bottom,
            "width": right - left, "height": bottom - top,
            "image_width": img_w, "image_height": img_h,
        }
 
    # ------------------------------------------------------------------
    # Parse the loose (string-y) input format.
    # positions entries: (pos, strand, gene_type, productive, locus_type)
    # ------------------------------------------------------------------
    def parse_positions():
        parsed = []
        for entry in positions:
            pos, strand, gene_type, productive, locus_type = entry
            if isinstance(productive, str):
                productive = productive.strip().lower() == "true"
            parsed.append({
                "position": int(pos),
                "strand": str(strand).strip(),
                "gene_type": str(gene_type).strip().upper(),
                "productive": bool(productive),
                "locus_type": str(locus_type).strip().upper(),
            })
        return parsed
 
    def combo_key(gene):
        return (gene["locus_type"], gene["gene_type"])
 
    def combo_label(locus_type, gene_type):
        # e.g. "IGH" + "V" -> "IGHV"
        return f"{locus_type}{gene_type}"
 
    def style_for(gene):
        color = LOCUS_GENE_SHADES.get(combo_key(gene), DEFAULT_SHADE)
        return color, "-", 1.0
 
    def dash_frac_for(gene, base_dash_frac):
        # Productive genes get a full-length dash; non-productive half.
        return base_dash_frac if gene["productive"] else base_dash_frac / 2
 
    # ------------------------------------------------------------------
    # Horizontal number line (goes under the box).
    # '+' dashes point UP (toward the box), '-' dashes point DOWN (away).
    # ------------------------------------------------------------------
    def render_horizontal_line(genes, width_px, dash_frac=0.6,
                                linewidth=1.2, dash_width=1.6):
        fig = plt.figure(figsize=(width_px / dpi, thickness_px / dpi), dpi=dpi)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(locus_start, locus_end)
        ax.set_ylim(-1, 1)
        ax.axis("off")
 
        ax.axhline(0, color="black", linewidth=linewidth, zorder=1)
 
        for gene in genes:
            color, linestyle, alpha = style_for(gene)
            d = dash_frac_for(gene, dash_frac)
            y0, y1 = (0, d) if gene["strand"] == "+" else (-d, 0)
            ax.plot([gene["position"], gene["position"]], [y0, y1],
                    color=color, linestyle=linestyle, linewidth=dash_width,
                    alpha=alpha, solid_capstyle="butt", zorder=2)
 
        tmp_path = "_tmp_horizontal_line.png"
        fig.savefig(tmp_path, dpi=dpi, transparent=True)
        plt.close(fig)
 
        img = Image.open(tmp_path)
        if img.width != width_px:
            img = img.crop((0, 0, width_px, img.height))
        return img
 
    # ------------------------------------------------------------------
    # Vertical number line (goes right of the box).
    # '+' dashes point LEFT (toward the box), '-' dashes point RIGHT (away).
    # Top of image = locus_start, bottom = locus_end.
    # ------------------------------------------------------------------
    def render_vertical_line(genes, height_px, dash_frac=0.6,
                              linewidth=1.2, dash_width=1.6):
        fig = plt.figure(figsize=(thickness_px / dpi, height_px / dpi), dpi=dpi)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(-1, 1)
        ax.set_ylim(locus_end, locus_start)  # inverted so top=start, bottom=end
        ax.axis("off")
 
        ax.axvline(0, color="black", linewidth=linewidth, zorder=1)
 
        for gene in genes:
            color, linestyle, alpha = style_for(gene)
            d = dash_frac_for(gene, dash_frac)
            x0, x1 = (-d, 0) if gene["strand"] == "+" else (0, d)
            ax.plot([x0, x1], [gene["position"], gene["position"]],
                    color=color, linestyle=linestyle, linewidth=dash_width,
                    alpha=alpha, solid_capstyle="butt", zorder=2)
 
        tmp_path = "_tmp_vertical_line.png"
        fig.savefig(tmp_path, dpi=dpi, transparent=True)
        plt.close(fig)
 
        img = Image.open(tmp_path)
        if img.height != height_px:
            img = img.crop((0, 0, img.width, height_px))
        return img
 
    def group_text_bands(row_has_dark, gap_tolerance=3):
        bands = []
        start = None
        blank_run = 0
        last_dark = None
        for r, has in enumerate(row_has_dark):
            if has:
                if start is None:
                    start = r
                blank_run = 0
                last_dark = r
            elif start is not None:
                blank_run += 1
                if blank_run > gap_tolerance:
                    bands.append((start, last_dark + 1))
                    start = None
        if start is not None:
            bands.append((start, last_dark + 1))
        return bands
 
    def split_oversized_band(r0, r1, row_counts, max_height=20, edge=2):
        if r1 - r0 <= max_height:
            return [(r0, r1)]
        lo, hi = r0 + edge, r1 - edge
        if hi <= lo:
            return [(r0, r1)]
        valley = lo + int(np.argmin(row_counts[lo:hi]))
        if valley <= r0 or valley >= r1:
            return [(r0, r1)]
        return (split_oversized_band(r0, valley, row_counts, max_height, edge)
                + split_oversized_band(valley, r1, row_counts, max_height, edge))
 
    def find_legend_zone(box, dark_thresh=200, pad=10, wide_line_frac=0.75):
        img = Image.open(dotplot_path).convert("L")
        arr = np.array(img)
        img_w = arr.shape[1]
 
        header = arr[:box["top"], :]
        dark = header < dark_thresh
        row_has_dark = dark.any(axis=1)
 
        bands = group_text_bands(row_has_dark)
 
        text_right_edge = 0
        for r0, r1 in bands:
            band_dark = dark[r0:r1, :]
            cols = np.where(band_dark.any(axis=0))[0]
            if len(cols) == 0:
                continue
            right = int(cols.max())
            if right / img_w <= wide_line_frac:  # skip full-width lines
                text_right_edge = max(text_right_edge, right)
 
        left = min(text_right_edge + pad, box["right"] - 1)
        right = box["right"]
        top = pad
        bottom = max(box["top"] - pad, top + 1)
 
        return {"left": left, "right": right, "top": top, "bottom": bottom,
                "width": max(right - left, 1), "height": max(bottom - top, 1)}
 
    def measure_gepard_line_ink_height(box, dark_thresh=128, wide_line_frac=0.75):
        img = Image.open(dotplot_path).convert("L")
        arr = np.array(img)
        img_w = arr.shape[1]
 
        header = arr[:box["top"], :]
        dark = header < dark_thresh
        row_has_dark = dark.any(axis=1)
        row_counts = dark.sum(axis=1)
 
        raw_bands = group_text_bands(row_has_dark)
        line_bands = []
        for r0, r1 in raw_bands:
            line_bands.extend(split_oversized_band(r0, r1, row_counts))
 
        candidates = []
        for r0, r1 in line_bands:
            band_dark = dark[r0:r1, :]
            cols = np.where(band_dark.any(axis=0))[0]
            if len(cols) == 0:
                continue
            if cols.max() / img_w <= wide_line_frac:  # skip title/axis-label
                rows = np.where(band_dark.any(axis=1))[0]
                candidates.append(int(rows.max() - rows.min() + 1))
 
        if not candidates:
            return None
        candidates.sort()
        return candidates[len(candidates) // 2]  # median line height
 
    def calibrate_font_size_pt(target_px_height, sample_text="Zoom: 14 : 1",
                                trial_pt=10.0):
        """
        Renders `sample_text` at `trial_pt` and measures its ink height, then
        scales linearly to find the font size that reproduces
        `target_px_height` -- i.e. the exact size of Gepard's own on-image
        text. This is the ONLY place font size is computed; it is fixed
        once here and reused everywhere below (never re-derived/variable).
        """
        fig = plt.figure(figsize=(4, 1), dpi=dpi)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(0, 400)
        ax.set_ylim(0, dpi)
        ax.axis("off")
        ax.text(10, dpi / 2, sample_text, fontsize=trial_pt,
                 family="monospace", va="center", ha="left", color="black")
        tmp_path = "_tmp_calibration.png"
        fig.savefig(tmp_path, dpi=dpi, transparent=True)
        plt.close(fig)
 
        arr = np.array(Image.open(tmp_path).convert("LA"))
        alpha_rows = np.where((arr[:, :, 1] > 10).any(axis=1))[0]
        if len(alpha_rows) == 0:
            return trial_pt
        measured_px = int(alpha_rows.max() - alpha_rows.min() + 1)
        return trial_pt * target_px_height / measured_px
 
    def render_legend_row(entries, row_height_px, font_size, max_width_px,
                           item_gap_px=14, swatch_gap_px=5, swatch_len_px=14):
        """
        Renders ONE horizontal row of legend entries: a short vertical
        tick swatch followed by its label, all inline. Cropped tightly
        to the content actually used (both width and height).
        """
        fig = plt.figure(figsize=(max_width_px / dpi, row_height_px / dpi), dpi=dpi)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(0, max_width_px)
        ax.set_ylim(0, row_height_px)
        ax.axis("off")
 
        y_mid = row_height_px / 2.0
        renderer = fig.canvas.get_renderer()
 
        x = 2.0
        for label, color, length_frac in entries:
            half = (swatch_len_px / 2.0) * length_frac
            ax.plot([x + half, x + half], [y_mid - half, y_mid + half],
                    color=color, linewidth=2.0, solid_capstyle="butt")
            x += swatch_len_px + swatch_gap_px
 
            txt = ax.text(x, y_mid, label, fontsize=font_size,
                           family="monospace", va="center", ha="left",
                           color="black")
            fig.canvas.draw()
            bbox = txt.get_window_extent(renderer=renderer)
            x += bbox.width + item_gap_px
 
        content_width_px = int(np.ceil(x - item_gap_px + 2))
 
        tmp_path = "_tmp_legend_row.png"
        fig.savefig(tmp_path, dpi=dpi, transparent=True)
        plt.close(fig)
 
        img = Image.open(tmp_path)
        content_width_px = min(content_width_px, img.width)
        img = img.crop((0, 0, content_width_px, img.height))
 
        # Crop vertically to the actual ink (tick marks + glyphs) too, so
        # stacking two rows later doesn't inherit unused figure padding.
        arr = np.array(img.convert("LA"))
        ink_rows = np.where((arr[:, :, 1] > 10).any(axis=1))[0]
        if len(ink_rows) > 0:
            top, bottom = int(ink_rows.min()), int(ink_rows.max()) + 1
            img = img.crop((0, top, img.width, bottom))
        return img
 
    def render_legend(combos_present, height_px, font_size,
                       max_width_px=3000, row_gap_px=4):
        """
        Two-row legend:
          row 1 -> the (locus_type, gene_type) combos present, labeled
                   as LOCUS_TYPE+GENE_TYPE (e.g. "IGHV", "IGHD", "IGHJ")
          row 2 -> Productive / Non-productive
        Both rows share the exact same (fixed) font size, and are cropped
        tight to their own ink so only a small fixed gap separates them.
        """
        row1_entries = [
            (combo_label(locus_type, gene_type),
             LOCUS_GENE_SHADES.get((locus_type, gene_type), DEFAULT_SHADE),
             1.0)
            for locus_type, gene_type in combos_present
        ]
        row2_entries = [("Productive", "#000000", 1.0),
                         ("Non-productive", "#000000", 0.5)]
 
        row_height_px = max(int(height_px / 2), 1)
 
        row1_img = render_legend_row(row1_entries, row_height_px, font_size, max_width_px)
        row2_img = render_legend_row(row2_entries, row_height_px, font_size, max_width_px)
 
        width = max(row1_img.width, row2_img.width)
        total_h = row1_img.height + row_gap_px + row2_img.height
        stacked = Image.new("RGBA", (width, total_h), (0, 0, 0, 0))
        stacked.paste(row1_img, (0, 0), row1_img)
        stacked.paste(row2_img, (0, row1_img.height + row_gap_px), row2_img)
        return stacked
 
    # ------------------------------------------------------------------
    # Main body
    # ------------------------------------------------------------------
    box = detect_gepard_box()
    genes = parse_positions()
 
    original = Image.open(dotplot_path).convert("RGBA")
    orig_w, orig_h = original.size
 
    bottom_line = render_horizontal_line(genes, box["width"])
    right_line = render_vertical_line(genes, box["height"])
 
    # Only include (locus_type, gene_type) combos that are actually
    # present in the data, ordered by gene type (V -> D -> J) then by
    # locus type alphabetically.
    present_combos_set = {combo_key(g) for g in genes}
    present_combos = sorted(
        present_combos_set,
        key=lambda c: (GENE_TYPE_ORDER.get(c[1], 99), c[0])
    )
 
    # Font size is computed exactly once, by matching Gepard's own
    # on-image text pixel height -- it is a fixed constant from here on,
    # never recomputed or overridden anywhere else.
    target_ink_px = measure_gepard_line_ink_height(box)
    font_size = (calibrate_font_size_pt(target_ink_px)
                 if target_ink_px else 8.0)
 
    legend_zone = find_legend_zone(box)
    legend_img = render_legend(present_combos, legend_zone["height"], font_size)
 
    # If the legend needs more horizontal room than the blank space to
    # the right of Gepard's info text provides, widen the whole canvas
    # so the font size never has to shrink to fit.
    available_width = (orig_w - legend_zone["left"]) + thickness_px
    extra_width = max(0, legend_img.width - available_width)
    canvas_w = orig_w + thickness_px + extra_width
    canvas_h = orig_h + thickness_px
 
    # Vertically center the (two-row) legend within the blank header band.
    legend_y = legend_zone["top"] + (legend_zone["height"] - legend_img.height) // 2
 
    canvas = Image.new("RGBA", (canvas_w, canvas_h), (255, 255, 255, 255))
    canvas.paste(original, (0, 0))  # original pixels copied through as-is
    canvas.paste(bottom_line, (box["left"], orig_h), bottom_line)
    canvas.paste(right_line, (orig_w, box["top"]), right_line)
    canvas.paste(legend_img, (legend_zone["left"], legend_y), legend_img)
 
    canvas.convert("RGB").save(output_path)
    return output_path

haplotypes=[]
for contig in contigs:
    found=False
    for h in haplotypes:
        if contig[1]==h[0]:
            found=True
            h[1].append(contig)
    if found==False:
        haplotypes.append([contig[1],[contig]])

for hap in haplotypes:
    if os.path.isdir(folder+"/"+hap[0]+'/locus_data')==False:
        os.mkdir(folder+"/"+hap[0]+'/locus_data')
    if delete=="-d":
        shutil.rmtree(folder+"/"+hap[0]+'/locus_data')
        os.mkdir(folder+"/"+hap[0]+'/locus_data')
    
    unique_contigs=[]
    for temp_con in hap[1]:
        if temp_con[0] not in unique_contigs:
            unique_contigs.append(temp_con[0])

    for ucon in unique_contigs:
        contig=[]
        loci=""
        for hcon in hap[1]:
            if ucon==hcon[0]:
                contig.append(hcon)
                if loci=="":
                    loci=str(hcon[2])
                else:
                    loci=loci+"-"+str(hcon[2])

        for g in genomes:
            if hap[0]==g[0]:
                genome=g[1]
                break

        subprocess.run(["samtools","faidx",genome])
        with open(genome+".fai","r") as read:
            reader=csv.reader(read,delimiter="\t")
            for row in reader:
                if row[0]==ucon:
                    contig_size=int(row[1])
                    break
            read.close()

        positions, data_dir = ig_detective_parse(hap[0],ucon,folder)

        if os.path.isdir(data_dir+'/locus_data/'+str(ucon)+"_"+str(loci))==False:
            os.mkdir(data_dir+'/locus_data/'+str(ucon)+"_"+str(loci))

        if plot_contigs.startswith("-p"):
            if plot_contigs.endswith("d"):
                if os.path.isfile(data_dir+'/locus_data/'+str(ucon)+"_"+str(loci)+"/"+(str(ucon)+"_contig_dotplot.png")):
                    os.remove(data_dir+'/locus_data/'+str(ucon)+"_"+str(loci)+"/"+(str(ucon)+"_contig_dotplot.png"))
                contig_run=True
            elif plot_contigs.endswith("k"):
                if os.path.isfile(data_dir+'/locus_data/'+str(ucon)+"_"+str(loci)+"/"+(str(ucon)+"_contig_dotplot.png")):
                    contig_run=False
                else:
                    contig_run=True
            if contig_run==True:
                contig_fasta=str(ucon)+".fasta"
                contig_start=int(max(1,int(sort_positions[0][0])-7500000))
                contig_end=int(min(contig_size,int(sort_positions[-1][0])+7500000))
                contig_region=""
                if contig_start==1 and contig_end==contig_size:
                    subprocess.run(["samtools", "faidx", genome, ucon, "-o", contig_fasta])
                else:
                    contig_region=str(ucon)+":"+str(contig_start)+"-"+str(contig_end)
                    subprocess.run(["samtools", "faidx", genome, contig_region, "-o", contig_fasta])
                with open(contig_fasta,"r") as read:
                    reader=csv.reader(read)
                    contig_fasta_data=[]
                    for row in reader:
                        contig_fasta_data.append(row)
                    read.close()
                with open(contig_fasta,"w",newline="") as write:
                    writer=csv.writer(write)
                    for l in contig_fasta_data:
                        if str(l[0]).startswith(">"):
                            if contig_region!="":
                                l[0]=">"+str(ucon)+" ("+contig_region.replace(str(ucon)+":","")+") ("+str(round((int(contig_end)-int(contig_start))/int(contig_size)*100,1))+"% of Contig)"
                            else:
                                l[0]=">"+str(ucon)+" (Full Contig)"
                        writer.writerow(l)
                
                subprocess.run(["java","-cp","/home/dhardesty/miniconda3/envs/gepard/share/gepard/dist/Gepard-2.1.jar","org.gepard.client.cmdline.CommandLine","-seq",contig_fasta,contig_fasta,"-matrix","/home/dhardesty/miniconda3/envs/gepard/share/gepard/src/matrices/edna.mat","-outfile",str(ucon)+"_contig_dotplot.png"])
                contig_dotplot=str(ucon)+"_contig_dotplot.png"
                os.rename(str(ucon)+"_contig_dotplot0.png",contig_dotplot)
                shutil.move(contig_dotplot, data_dir+'/locus_data/'+str(ucon)+"_"+str(loci)+"/"+contig_dotplot)
                shutil.move(contig_fasta, data_dir+'/locus_data/'+str(ucon)+"_"+str(loci)+"/"+contig_fasta)

        def plot_locus(sort_positions, ucon, contig_size, loci, data_dir, output_name):
            true_locus_size=int(sort_positions[-1][0])-int(sort_positions[0][0])
            if true_locus_size!=0:
                start_pad = int(sort_positions[0][0])-int(true_locus_size*0.1)
                end_pad = int(sort_positions[-1][0])+int(true_locus_size*0.1)
            else:
                start_pad = int(sort_positions[0][0])-1000
                end_pad = int(sort_positions[-1][0])+1000
            locus_start=int(max(1,start_pad))
            locus_end=int(min(contig_size,end_pad))

            locus_region=str(ucon)+":"+str(locus_start)+"-"+str(locus_end)
            locus_fasta=str(ucon)+"_"+str(loci)+".fasta"
            subprocess.run(["samtools", "faidx", genome, locus_region, "-o", locus_fasta])
            with open(locus_fasta,"r") as read:
                reader=csv.reader(read)
                locus_fasta_data=[]
                for row in reader:
                    locus_fasta_data.append(row)
                read.close()
            with open(locus_fasta,"w",newline="") as write:
                writer=csv.writer(write)
                for l in locus_fasta_data:
                    if str(l[0]).startswith(">"):
                        l[0]=str(l[0]).replace(locus_region, locus_fasta.replace(".fasta",""))+" ("+locus_region.replace(str(ucon)+":","")+")"
                    writer.writerow(l)
            locus_dotplot=str(ucon)+"_"+str(loci)+"_locus_dotplot"+output_name+".png"

            if os.path.isfile(data_dir+'/locus_data/'+str(ucon)+"_"+str(loci)+"/"+locus_dotplot)==True:
                os.remove(data_dir+'/locus_data/'+str(ucon)+"_"+str(loci)+"/"+locus_dotplot)
            if os.path.isfile(data_dir+'/locus_data/'+str(ucon)+"_"+str(loci)+"/"+locus_fasta)==True:
                os.remove(data_dir+'/locus_data/'+str(ucon)+"_"+str(loci)+"/"+locus_fasta)
            
            subprocess.run(["java","-cp","/home/dhardesty/miniconda3/envs/gepard/share/gepard/dist/Gepard-2.1.jar","org.gepard.client.cmdline.CommandLine","-seq",locus_fasta,locus_fasta,"-matrix","/home/dhardesty/miniconda3/envs/gepard/share/gepard/src/matrices/edna.mat","-outfile",str(ucon)+"_"+str(loci)+"_locus_dotplot.png"])
            os.rename(str(ucon)+"_"+str(loci)+"_locus_dotplot0.png",locus_dotplot)
            
            print(contig)

            print(sort_positions)
            print(locus_start, locus_end)
            print(locus_dotplot)
            
            annotate_gepard_dotplot(sort_positions, locus_start, locus_end, locus_dotplot, locus_dotplot.replace(".png","_annotated.png"))
            os.remove("_tmp_calibration.png")
            os.remove("_tmp_horizontal_line.png")
            os.remove("_tmp_legend_row.png")
            os.remove("_tmp_vertical_line.png")

            shutil.move(locus_dotplot.replace(".png","_annotated.png"), data_dir+'/locus_data/'+str(ucon)+"_"+str(loci)+"/"+locus_dotplot)
            os.remove(locus_dotplot)
            shutil.move(locus_fasta, data_dir+'/locus_data/'+str(ucon)+"_"+str(loci)+"/"+locus_fasta)
        
        sort_positions=sorted(positions, key=lambda x: x[0])
        if "TRD" in loci or "IGH" in loci or "TRB" in loci:
            outname="_without_d"
        else:
            outname=""
        plot_locus(sort_positions, ucon, contig_size, loci, data_dir, outname)

        if "TRD" in loci or "IGH" in loci or "TRB" in loci:
            d_genes=[]
            for d in all_d_genes:
                if d[0]==hap[0] and d[2]==ucon and str(d[7]) in loci:
                    d_genes.append([int(d[3]),d[4],d[1],"True",d[7]])
            
            j_genes=[]
            for j in all_j_genes:
                if j[0]==hap[0] and j[2]==ucon and str(j[7]) in loci:
                    j_genes.append([int(j[3]),j[4],j[1],"True",j[7]])

            for d in d_genes:
                positions.append(d)
            for j in j_genes:
                positions.append(j)
            sort_positions=sorted(positions, key=lambda x: x[0])
            plot_locus(sort_positions, ucon, contig_size, loci, data_dir, "_with_d")