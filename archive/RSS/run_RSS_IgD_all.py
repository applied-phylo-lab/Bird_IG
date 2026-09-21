#!/usr/bin/env python3

import argparse
import os
import subprocess
import pandas as pd
from glob import glob
from multiprocessing import Pool

def pick_fasta_for_haplotype(fasta_dir, haplotype, species):
    """
    Determine which FASTA should be used for the given haplotype.
    Handles pri/alt naming rules.
    """
    fasta_files = glob(os.path.join(fasta_dir, "*.fna"))

    if len(fasta_files) == 0:
        fasta_files = glob(os.path.join(fasta_dir+"_", "*.fna"))
        if len(fasta_files) == 0:
            raise FileNotFoundError(f"No FASTA files found in {fasta_dir}")

    # If only one fasta → assume primary haplotype
    if len(fasta_files) == 1:
        return fasta_files[0]

    # Haplotype naming rules
    hap = haplotype.lower()

    # pri / primary / mat / hap1 / p / maternal
    pri_markers = ["pri", "primary", "mat", "hap1", "p", "maternal"]

    # alt / alternate / hap2 / pat / a / paternal
    alt_markers = ["alt", "alternate", "hap2", "pat", "a", "paternal"]

    if any(hap.endswith(m) for m in pri_markers):
        wanted = pri_markers
    elif any(hap.endswith(m) for m in alt_markers):
        wanted = alt_markers
    else:
        # unknown naming → assume primary
        wanted = pri_markers

    # Search for FASTA filenames containing any marker
    for f in fasta_files:
        name = os.path.basename(f).lower()
        if any(m in name for m in wanted):
            return f

    # If nothing matches markers → fallback
    return fasta_files[0]


def run_igd_task(task):
    """
    Execute IGDetective.py with correct arguments.
    """
    (igd_dir, fasta, out_path) = task

    os.makedirs(out_path, exist_ok=True)

    cmd = [
        "python",
        os.path.join(igd_dir, "py", "IGDetective.py"),
        "-i", fasta,
        "-o", out_path,
        "-m", "1",
        "-r",
        "-l", "IGH"
    ]

    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=igd_dir)


def main():
    parser = argparse.ArgumentParser(description="Batch-run IGDetective for many samples")
    parser.add_argument("-i", "--input_dir", required=True)
    parser.add_argument("-f", "--fasta_dir", required=True)
    parser.add_argument("-s", "--summary", required=True)
    parser.add_argument("-d", "--igd_dir", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-c", "--cores", type=int, default=4)

    args = parser.parse_args()

    df = pd.read_csv(args.summary, sep="\t")
    df['InputDir'] = args.input_dir
    tasks = []
    for idx, row in df.iterrows():

        order = row["Order"]
        species = row["Species"]
        haplotype = row["Haplotype"]
        input_subdir = row["InputDir"]
        contig = row["Contig"]
        numv = row["NumV"]

        # FASTA selection
        fasta_path = pick_fasta_for_haplotype(
            os.path.join(args.fasta_dir, "#"+order,species),
            haplotype,
            species
        )

        # Output directory
        out_path = os.path.join(
            args.output,
            order,
            species,
            haplotype,
            "RSS"
        )
        tasks.append((args.igd_dir, fasta_path, out_path))
    
    with Pool(processes=args.cores) as pool:
        pool.map(run_igd_task, tasks)
        
       

if __name__ == "__main__":
    main()
