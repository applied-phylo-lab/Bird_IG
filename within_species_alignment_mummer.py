#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import itertools
import subprocess
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed

# -------------------------
# Argument parsing
# -------------------------
def parse_args():
    parser = argparse.ArgumentParser(
        description="Pairwise MUMmer alignments for IGH FASTAs with dot & synteny plots (parallelized)"
    )
    parser.add_argument("-i", "--input_dir", required=True)
    parser.add_argument("-s", "--summary_tsv", required=True)
    parser.add_argument("-o", "--output_dir", required=True)
    parser.add_argument("-c", "--cores", type=int, default=1,
                        help="Number of parallel processes (default: 1)")
    parser.add_argument("--nucmer_path", default="nucmer")
    parser.add_argument("-m","--mummerplot_path", default="mummerplot")
    return parser.parse_args()

# -------------------------
# Utility
# -------------------------
def run_cmd(cmd, stdout_file=None):
    print("Running:", " ".join(cmd))
    if stdout_file:
        with open(stdout_file, "w") as fh:
            subprocess.run(cmd, check=True, stdout=fh)
    else:
        subprocess.run(cmd, check=True)

# -------------------------
# Plotting
# -------------------------
def create_dotplot(delta_file, output_prefix, mummerplot_path):
    png_file = output_prefix + ".png"
    if os.path.exists(png_file):
        return

    cmd = [
        mummerplot_path,
        "--png",
        "--layout",
        "-p", output_prefix,
        delta_file
    ]
    run_cmd(cmd)


# -------------------------
# Core worker function
# -------------------------
def process_pair(args_tuple):
    (
        args,
        order,
        species,
        h1,
        h2
    ) = args_tuple

    fasta1 = os.path.join(
        args.input_dir, order, h1['Species'], h1['Haplotype'],
        'refined_ig_loci', 'igloci_fasta',
        f"IGH_{h1['Contig']}_{h1['NumV']}Vs.fasta"
    )
    fasta2 = os.path.join(
        args.input_dir, order, h2['Species'], h2['Haplotype'],
        'refined_ig_loci', 'igloci_fasta',
        f"IGH_{h2['Contig']}_{h2['NumV']}Vs.fasta"
    )

    if not os.path.exists(fasta1) or not os.path.exists(fasta2):
        return f"Skipping missing FASTA: {species} {h1['Haplotype']} vs {h2['Haplotype']}"

    out_prefix = os.path.join(
        args.output_dir,
        f"{species}_{h1['Haplotype']}_vs_{h2['Haplotype']}"
    )

    delta = out_prefix + ".delta"
    filtered_delta = out_prefix + ".filtered.delta"
    coords = out_prefix + ".coords"

    # nucmer
    if not os.path.exists(delta):
        run_cmd([
            args.nucmer_path,
            "--maxmatch",
            "-p", out_prefix,
            fasta1,
            fasta2
        ])

    # delta-filter
    if not os.path.exists(filtered_delta):
        run_cmd(
            ["delta-filter", "-r", "-q", delta],
            stdout_file=filtered_delta
        )

    # show-coords
    run_cmd(
            ["show-coords", "-H", "-T", "-r", "-c", "-l", filtered_delta],
            stdout_file=coords
        )

    # plots
    #create_dotplot(filtered_delta, out_prefix, args.mummerplot_path)
    #create_synteny_plot(coords, out_prefix)

    return f"Finished {species}: {h1['Haplotype']} vs {h2['Haplotype']}"

# -------------------------
# Main
# -------------------------
def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    df = pd.read_csv(args.summary_tsv, sep="\t")

    jobs = []
    for order, order_df in df.groupby("Order"):
        for species, sp_df in order_df.groupby("Species"):
            haplotypes = sp_df.to_dict(orient="records")
            for h1, h2 in itertools.combinations(haplotypes, 2):
                jobs.append((args, order, species, h1, h2))

    print(f"Total pairwise jobs: {len(jobs)}")
    print(f"Using {args.cores} cores")

    with ProcessPoolExecutor(max_workers=args.cores) as executor:
        futures = [executor.submit(process_pair, job) for job in jobs]
        for f in as_completed(futures):
            print(f.result())

    print("All jobs completed.")

if __name__ == "__main__":
    main()
