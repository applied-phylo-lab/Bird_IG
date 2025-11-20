#!/usr/bin/env python3
"""
summarize_inversion_units_parallel_optimized.py

Calculate minimal repeated inversion units from IGH_self.tsv files.

For each species/haplotype:
  - Reads LASTZ self-alignment (IGH_self.tsv)
  - Identifies inversions (strand1 != strand2)
  - Computes intercept b for each inversion assuming slope m = -1
  - For each pair of inversions that overlap in both coordinates,
    computes perpendicular distance between their anti-diagonal lines:
        d = |b1 - b2| / sqrt(1 + m^2)
    where m = -1 → denominator = sqrt(2)
  - The minimal repeat unit is min(distance, smaller_inversion_length)
  - Reports the global minimal repeat unit across all pairs

Usage:
  python summarize_inversion_units_parallel_optimized.py \
      -i /top/input/dir \
      -s summary_table_IGH.tsv \
      -o out.tsv \
      -c 4
"""

import os
import argparse
import pandas as pd
import numpy as np
from multiprocessing import Pool
from itertools import combinations
import math


def read_tsv_with_header(path):
    """Read TSV even if header starts with #"""
    with open(path) as f:
        header = f.readline().strip()
        if header.startswith("#"):
            header = header[1:]
        columns = header.split("\t")
    df = pd.read_csv(path, sep="\t", comment="#", names=columns, skiprows=1)
    df.columns = [c.replace("#", "").replace("%", "").replace("+", "") for c in df.columns]
    return df


def inversion_line_b(inv):
    """
    Given a row (with numeric start1,end1,start2,end2),
    return slope m=-1 and intercept b for the line passing through
    points (start1, start2) and (end1, end2).

    To be robust to noise, compute b from both endpoints and average:
      b_start = y_start - m * x_start
      b_end   = y_end   - m * x_end
      b = (b_start + b_end) / 2

    Note: start2 is negated here so that inversions (anti-diagonals)
    map correctly in this coordinate system.
    """
    m = -1.0
    x_start = float(inv["start1"])
    y_start = float(inv["end2"])
    x_end = float(inv["end1"])
    y_end = float(inv["start2"])

    b_start = y_start - m * x_start
    b_end = y_end - m * x_end
    b = 0.5 * (b_start + b_end)
    return b


def blocks_overlap(a_start, a_end, b_start, b_end):
    """Return True if [a_start, a_end] and [b_start, b_end] overlap (inclusive)."""
    return not (a_end < b_start or b_end < a_start)


def process_row_dict(row_dict, top_input_dir):
    """
    Process one row of the summary table and compute minimal inversion unit.
    """
    order = str(row_dict.get("Order", "")).strip()
    species = str(row_dict.get("Species", "")).strip()
    haplotype = str(row_dict.get("Haplotype", "")).strip()
    input_dir_field = str(row_dict.get("InputDir", "")).strip()


    aln_path = os.path.join(top_input_dir, order, species, haplotype, "IGH_self.tsv")

    out = {
        "InputDir": input_dir_field,
        "Order": order,
        "Species": species,
        "Haplotype": haplotype,
        "AlnPath": aln_path,
        "NumInversions": 0,
        "MinimalUnit_bp": np.nan,
        "Note": ""
    }

    # --- check file existence ---
    if not os.path.exists(aln_path):
        out["Note"] = "missing_alignment_file"
        return out

    # --- read alignment ---
    try:
        df = read_tsv_with_header(aln_path)
    except Exception as e:
        out["Note"] = f"read_error:{e}"
        return out

    required_cols = {"start1", "end1", "start2", "end2", "strand1", "strand2"}
    if not required_cols.issubset(set(df.columns)):
        out["Note"] = "missing_required_columns"
        return out

    # --- clean numeric columns ---
    for c in ("start1", "end1", "start2", "end2"):
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # --- detect inversions ---
    df["inversion"] = df["strand1"] != df["strand2"]
    inv_df = df[df["inversion"]].copy()
    inv_df = inv_df.dropna(subset=["start1", "end1", "start2", "end2"])
    inv_df = inv_df.reset_index(drop=True)

    n_inv = len(inv_df)
    out["NumInversions"] = int(n_inv)
    if n_inv < 2:
        out["Note"] = "fewer_than_2_inversions"
        return out

    # --- compute inversion length + intercept b ---
    inv_df["inv_len"] = (inv_df["end1"] - inv_df["start1"]).abs() + \
                        (inv_df["end2"] - inv_df["start2"]).abs()
    inv_df["b"] = inv_df.apply(inversion_line_b, axis=1)

    # --- compute pairwise minimal units ---
    m = -1.0
    denom = math.sqrt(1 + m * m)  # = sqrt(2)

    pair_minimals = []

    for (i, A), (j, B) in combinations(inv_df.iterrows(), 2):
        # check for overlap
        overlap_first = blocks_overlap(A["start1"], A["end1"], B["start1"], B["end1"])
        overlap_second = blocks_overlap(A["start2"], A["end2"], B["start2"], B["end2"])
        if not (overlap_first and overlap_second):
            # skip non-overlapping inversions
            continue

        smaller_len = min(A["inv_len"], B["inv_len"])
        dist = abs(A["b"] - B["b"]) / denom
        pair_minimals.append(min(dist, smaller_len))

    if len(pair_minimals) == 0:
        out["Note"] = "no_overlapping_pairs"
        return out

    out["MinimalUnit_bp"] = float(min(pair_minimals))
    out["Note"] = f"computed_from_{len(pair_minimals)}_pairs"
    print(f"Processed: {order}/{species}/{haplotype}")
    return out


def main():
    parser = argparse.ArgumentParser(description="Summarize inversion stats from LASTZ alignments and BED files")
    parser.add_argument("-i", "--input_dir", required=True, help="Top-level input directory")
    parser.add_argument("-s", "--summary", required=True, help="Path to summary_table_IGH.tsv")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    parser.add_argument("-c", "--cores", type=int, default=4, help="Number of parallel workers")
    args = parser.parse_args()

    summary_df = pd.read_csv(args.summary, sep="\t", dtype=str)
    tasks = [(row[1].to_dict(), args.input_dir) for row in summary_df.iterrows()]

    with Pool(processes=args.cores) as pool:
        results = pool.starmap(process_row_dict, tasks)

    out_df = pd.DataFrame(results)
    out_df.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote: {args.output}")


if __name__ == "__main__":
    main()
