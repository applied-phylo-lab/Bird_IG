#!/usr/bin/env python3
"""
For each row in the summary table, build a tip-color mapping from gene_list.csv
and call plot_tree.R to render a ladderized rectangular SVG tree.
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
    """FASTA headers replace '+' with '_'."""
    return "_" if strand == "+" else strand


def build_color_lookup(gene_list_df, order, species, haplotype):
    source = f"{order}/{species}/{haplotype}"
    sub = gene_list_df[gene_list_df["Source"] == source]

    lookup = {}
    for _, row in sub.iterrows():
        gene_type = str(row["GeneType"])
        strand    = strand_to_tree(str(row["Strand"]))

        productive_raw = str(row["Productive"]).strip().lower()
        productive     = productive_raw in ("true", "1", "yes")

        heptamer = str(row.get("Heptamer", "")).strip()
        has_hept = heptamer not in ("", "nan", "None")

        tip_name = (
            f"{haplotype}.{row['Pos']}.{row['Contig']}"
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
    label     = f"{species}/{haplotype}"

    output_dir = os.path.join(input_dir, order, species, haplotype, "tree")
    treefile   = os.path.join(output_dir, f"{haplotype}_tree.treefile")

    if not os.path.isfile(treefile):
        print(f"[{label}][WARNING] Tree file not found, skipping: {treefile}",
              file=sys.stderr, flush=True)
        return

    color_lookup = build_color_lookup(gene_list_df, order, species, haplotype)
    output_svg   = os.path.join(output_dir, f"{haplotype}_tree.svg")

    # Write color mapping to a temp TSV for R to consume
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
        description="Plot IGH phylogenetic trees as colored SVGs"
    )
    parser.add_argument("-i", "--input_dir", required=True, help="Top-level input directory")
    parser.add_argument("-s", "--summary",   required=True, help="Path to summary_table_IGH.tsv")
    parser.add_argument("-c", "--cores",     type=int, default=4, help="Number of parallel processes")

    args = parser.parse_args()

    gene_list_path = os.path.join(args.input_dir, "gene_list.csv")
    if not os.path.isfile(gene_list_path):
        print(f"[ERROR] gene_list.csv not found: {gene_list_path}", file=sys.stderr)
        sys.exit(1)

    gene_list_df = pd.read_csv(gene_list_path)

    df = pd.read_csv(args.summary, sep="\t")
    df["InputDir"] = args.input_dir

    tasks = [(row, gene_list_df) for _, row in df.iterrows()]

    print(f"Processing {len(df)} trees with {args.cores} parallel workers...\n")

    with Pool(args.cores) as pool:
        pool.map(process_row, tasks)

    print("\nDone.")


if __name__ == "__main__":
    main()
