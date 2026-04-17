#!/usr/bin/env python3
"""
Pipeline script: for each row in a summary TSV, runs createFastaFromCSV.py,
aligns with clustalo, and builds a phylogenetic tree with iqtree2.
Rows are processed in parallel using multiprocessing.Pool.
"""

import argparse
import os
import subprocess
import sys
import pandas as pd
from multiprocessing import Pool


CREATE_FASTA_SCRIPT = "createFastaFromCSV.py"


def run_cmd(cmd, step_name, label=""):
    """Run a shell command, raising RuntimeError on failure."""
    tag = f"[{label}]" if label else ""
    print(f"{tag}[{step_name}] Running: {' '.join(cmd)}", flush=True)
    result = subprocess.run(cmd, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"{step_name} failed with return code {result.returncode} for {label}"
        )
    print(f"{tag}[{step_name}] Done.", flush=True)


def process_row(row):
    input_dir = row['InputDir']
    order     = row['Order']
    species   = row['Species']
    haplotype = row['Haplotype']

    label = f"{species}/{haplotype}"

    combined_genes_path = os.path.join(
        input_dir, order, species, haplotype,
        'combined_genes_IGH_clean.txt'
    )

    output_dir  = os.path.join(input_dir, order, species, haplotype, "tree")
    os.makedirs(output_dir, exist_ok=True)

    prefix      = haplotype
    fasta_out   = os.path.join(output_dir, f"{prefix}.fasta")
    aligned_out = os.path.join(output_dir, f"{prefix}_aligned.fasta")

    print(f"\n[{label}] Starting pipeline", flush=True)
    print(f"[{label}]   Input CSV : {combined_genes_path}", flush=True)
    print(f"[{label}]   Output dir: {output_dir}", flush=True)

    if not os.path.isfile(combined_genes_path):
        print(f"[{label}][WARNING] Input file not found, skipping: {combined_genes_path}",
              file=sys.stderr, flush=True)
        return

    try:
        # ------------------------------------------------------------------ #
        # Step 1 – createFastaFromCSV.py
        # ------------------------------------------------------------------ #
        run_cmd(
            ["python3", CREATE_FASTA_SCRIPT, combined_genes_path, output_dir, prefix],
            step_name="createFastaFromCSV",
            label=label,
        )

        if not os.path.isfile(fasta_out):
            raise FileNotFoundError(
                f"Expected FASTA not found after createFastaFromCSV: {fasta_out}"
            )

        # ------------------------------------------------------------------ #
        # Step 2 – Multiple sequence alignment with clustalo
        # ------------------------------------------------------------------ #
        run_cmd(
            ["clustalo", "-i", fasta_out, "-o", aligned_out, "-t", "DNA", "--force"],
            step_name="clustalo",
            label=label,
        )

        # ------------------------------------------------------------------ #
        # Step 3 – Phylogenetic tree with IQ-TREE 2
        # ------------------------------------------------------------------ #
        iqtree_prefix = os.path.join(output_dir, f"{prefix}_tree")
        run_cmd(
            ["iqtree2", "-s", aligned_out, "--prefix", iqtree_prefix, "--redo"],
            step_name="iqtree2",
            label=label,
        )

        print(f"[{label}] Pipeline complete.", flush=True)

    except Exception as e:
        print(f"[{label}][ERROR] {e}", file=sys.stderr, flush=True)


def main():
    parser = argparse.ArgumentParser(
        description="Create FASTA, align, and build phylogenetic trees for IGH loci."
    )
    parser.add_argument("-i", "--input_dir",  required=True,
                        help="Top-level input directory")
    parser.add_argument("-s", "--summary",    required=True,
                        help="Path to summary TSV file")
    parser.add_argument("-c", "--cores",      type=int, default=4,
                        help="Number of parallel processes (default: 4)")
    parser.add_argument("--create-fasta-script", default="tree_analyses/createFastaFromCSV.py",
                        dest="create_fasta_script",
                        help="Path to createFastaFromCSV.py (default: %(default)s)")

    args = parser.parse_args()

    # Allow overriding the script path at runtime
    global CREATE_FASTA_SCRIPT
    CREATE_FASTA_SCRIPT = args.create_fasta_script

    if not os.path.isfile(args.summary):
        print(f"[ERROR] Summary file not found: {args.summary}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(args.summary, sep='\t')
    df['InputDir'] = args.input_dir

    print(f"Processing {len(df)} rows with {args.cores} parallel workers...\n")

    with Pool(args.cores) as pool:
        pool.map(process_row, [row for _, row in df.iterrows()])

    print("\nAll rows processed.")


if __name__ == "__main__":
    main()