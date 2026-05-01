#!/usr/bin/env python3
"""
Plot IGH and IGL phylogenetic trees as colored SVGs.

- IGH: skip if SVG already exists
- IGL: always run
- gene_list.csv is filtered by Locus
"""

import argparse
import os
import sys
import subprocess
import tempfile
import pandas as pd
from multiprocessing import Pool

PLOT_R_SCRIPT = os.path.join(os.path.dirname(__file__), "plot_tree.R")

COLOR_MAP = {
    (True,  True):  "#172869",
    (True,  False): "#0076BB",
    (False, True):  "#EA7580",
    (False, False): "#F8CD9C",
}

UNMATCHED_COLOR = "#999999"


def strand_to_tree(strand):
    return "_" if strand == "+" else strand


def build_color_lookup(gene_list_df, order, species, haplotype, locus,prefix):
    source = f"{order}/{species}/{haplotype}"

    sub = gene_list_df[
        (gene_list_df["Source"] == source) &
        (gene_list_df["Locus"] == locus)
    ]
    #print(sub)
    lookup = {}
    for _, row in sub.iterrows():
        gene_type = str(row["GeneType"])
        strand    = strand_to_tree(str(row["Strand"]))

        productive_raw = str(row["Productive"]).strip().lower()
        productive     = productive_raw in ("true", "1", "yes")

        heptamer = str(row.get("Heptamer", "")).strip()
        has_hept = heptamer not in ("", "nan", "None")

        tip_name = (
            f"{prefix}{haplotype}.{row['Pos']}.{row['Contig']}"
            f".{gene_type}.{row['Productive']}.{strand}"
        )

        lookup[tip_name] = COLOR_MAP[(productive, has_hept)]

    return lookup


def process_row(args):
    row, gene_list_df = args

    input_dir = row["InputDir"]
    order     = row["Order"]
    species   = row["Species"]
    haplotype = row["Haplotype"]
    locus     = row["Locus"]

    label = f"{species}/{haplotype}/{locus}"

    output_dir = os.path.join(input_dir, order, species, haplotype, "tree")

    # --- File naming depending on locus ---
    prefix = "IGL_" if locus == "IGL" else ""

    treefile   = os.path.join(output_dir, f"{prefix}{haplotype}_tree.treefile")
    output_svg = os.path.join(output_dir, f"{prefix}{haplotype}_tree.svg")

    # --- Skip if already exists ---
    if locus=="IGH" and os.path.isfile(output_svg):
        #print(f"[{label}] SVG exists, skipping.", flush=True)
        return

    if not os.path.isfile(treefile):
        #print(f"[{label}][WARNING] Tree file not found, skipping: {treefile}",
        #      file=sys.stderr, flush=True)
        return

    color_lookup = build_color_lookup(
        gene_list_df, order, species, haplotype, locus,prefix
    )

    # Write temp color mapping
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", delete=False, prefix=f"{haplotype}_colors_"
    ) as fh:
        color_tsv = fh.name
        fh.write("tip_name\tcolor\n")
        for tip, color in color_lookup.items():
            fh.write(f"{tip}\t{color}\n")

    try:
        result = subprocess.run(
            ["Rscript", PLOT_R_SCRIPT, treefile, color_tsv, output_svg],
            text=True, capture_output=True,
        )
        if result.returncode != 0:
            print(f"[{label}][ERROR] R script failed:\n{result.stderr}",
                  file=sys.stderr, flush=True)
        else:
            print(f"[{label}] Saved: {output_svg}", flush=True)
            if result.stderr.strip():
                print(f"[{label}][R] {result.stderr.strip()}", flush=True)
    finally:
        os.unlink(color_tsv)


def main():
    parser = argparse.ArgumentParser(
        description="Plot IGH + IGL phylogenetic trees"
    )
    parser.add_argument("-i", "--input_dir", required=True)
    parser.add_argument("-s", "--summary",   required=True)
    parser.add_argument("-c", "--cores",     type=int, default=4)

    args = parser.parse_args()

    gene_list_path = os.path.join(args.input_dir, "gene_list.csv")
    if not os.path.isfile(gene_list_path):
        print(f"[ERROR] gene_list.csv not found: {gene_list_path}", file=sys.stderr)
        sys.exit(1)

    gene_list_df = pd.read_csv(gene_list_path)

    df = pd.read_csv(args.summary, sep="\t")
    df["InputDir"] = args.input_dir

    tasks = [(row, gene_list_df) for _, row in df.iterrows()]

    print(f"Processing {len(tasks)} trees with {args.cores} workers...\n")

    with Pool(args.cores) as pool:
        pool.map(process_row, tasks)

    print("\nDone.")


if __name__ == "__main__":
    main()