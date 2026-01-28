#!/usr/bin/env python3

import os
import argparse
import subprocess
import pandas as pd
import gzip
import re
from concurrent.futures import ProcessPoolExecutor, as_completed

# -----------------------------
# RepeatMasker parameters
# -----------------------------
RM_ENGINE = "hmmer"
RM_SPECIES = "aves"

# -----------------------------
# Helpers
# -----------------------------
def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)

def infer_class_family(repeat_name):
    r = repeat_name.upper()
    if r.startswith("CR1"):
        return "LINE/CR1"
    if r.startswith("ERV"):
        return "LTR/ERV"
    if r.startswith("GYPSY"):
        return "LTR/Gypsy"
    if r.startswith("COPIA"):
        return "LTR/Copia"
    if r.startswith("HAT"):
        return "DNA/hAT"
    if r.startswith("TCMAR"):
        return "DNA/TcMar"
    if "SAT" in r:
        return "Satellite"
    if r.startswith("(") or "SIMPLE" in r:
        return "Simple_repeat"
    return "Other"

def gff_to_bed(gff_file, bed_file):
    with open_maybe_gzip(gff_file) as infile, open(bed_file, "w") as out:
        for line in infile:
            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")
            if len(fields) < 9:
                continue

            chrom, _, _, start, end, score, strand, _, attrs = fields
            bed_start = int(start) - 1
            bed_end = int(end)

            m = re.search(r'Target\s+"Motif:([^"]+)"', attrs)
            if not m:
                continue

            repeat_name = m.group(1)
            class_family = infer_class_family(repeat_name)

            out.write(
                f"{chrom}\t{bed_start}\t{bed_end}\t"
                f"{repeat_name}\t{score}\t{strand}\t{class_family}\n"
            )

# -----------------------------
# Worker function
# -----------------------------
def process_row(row, input_dir):
    order = row["Order"]
    species = row["Species"]
    hap = row["Haplotype"]
    contig = row["Contig"]
    numv = row["NumV"]

    fasta = os.path.join(
        input_dir,
        order,
        species,
        hap,
        "refined_ig_loci",
        "igloci_fasta",
        f"IGH_{contig}_{numv}Vs.fasta"
    )

    if not os.path.exists(fasta):
        return f"[missing fasta] {fasta}"

    outdir = os.path.join(
        input_dir,
        order,
        species,
        hap,
        "RepeatMasker"
    )
    os.makedirs(outdir, exist_ok=True)

    gff = os.path.join(outdir, os.path.basename(fasta) + ".out.gff")
    bed = os.path.join(outdir, os.path.basename(fasta) + ".repeats.bed")

    # -----------------------------
    # Run RepeatMasker if needed
    # -----------------------------
    if not os.path.exists(gff):
        cmd = [
            "RepeatMasker",
            "-pa", "1",                 # IMPORTANT: 1 thread per job
            "-engine", RM_ENGINE,
            "-species", RM_SPECIES,
            "-gff",
            "-nolow",
            "-dir", outdir,
            fasta
        ]
        subprocess.run(cmd, check=True)

    # -----------------------------
    # Convert GFF → BED
    # -----------------------------
    if os.path.exists(gff) and not os.path.exists(bed):
        gff_to_bed(gff, bed)

    return f"[done] {hap} {contig}"

# -----------------------------
# Main
# -----------------------------
def main(args):
    summary = pd.read_csv(args.summary, sep="\t")

    rows = [row for _, row in summary.iterrows()]

    with ProcessPoolExecutor(max_workers=args.cores) as executor:
        futures = [
            executor.submit(process_row, row, args.input_dir)
            for row in rows
        ]

        for f in as_completed(futures):
            try:
                print(f.result())
            except Exception as e:
                print(f"[error] {e}")

# -----------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parallel RepeatMasker on IGH FASTAs + GFF to BED"
    )
    parser.add_argument(
        "-s", "--summary", required=True,
        help="Summary TSV"
    )
    parser.add_argument(
        "-i", "--input_dir", required=True,
        help="Base input directory"
    )
    parser.add_argument(
        "-c", "--cores", type=int, default=4,
        help="Number of parallel RepeatMasker jobs"
    )

    args = parser.parse_args()
    main(args)
