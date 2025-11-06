#!/usr/bin/env python3
"""
summarize_inversion_units_parallel.py

Calculate minimal repeated inversion units from IGH_self.tsv files.
Uses the geometric shortcut for anti-diagonal inversions (slope m = -1):
  treat each inversion as a line y = m*x + b, compute b, distance between lines:
    d = |b1 - b2| / sqrt(1 + m^2)  (with m = -1 -> denom = sqrt(2))

Usage:
  python summarize_inversion_units_parallel.py -i /top/input/dir -s summary_table_IGH.tsv -o out.tsv -c 4
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

def process_row_dict(row_dict, top_input_dir):
    order = str(row_dict.get("Order", "")).strip()
    species = str(row_dict.get("Species", "")).strip()
    haplotype = str(row_dict.get("Haplotype", "")).strip()
    input_dir_field = str(row_dict.get("InputDir", "")).strip()

    if input_dir_field and os.path.isabs(input_dir_field):
        base_dir = input_dir_field
    elif input_dir_field:
        base_dir = os.path.join(top_input_dir, input_dir_field)
    else:
        base_dir = os.path.join(top_input_dir, order, species, haplotype)

    aln_path = os.path.join(base_dir, order, species, haplotype, "IGH_self.tsv")

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

    if not os.path.exists(aln_path):
        out["Note"] = "missing_alignment_file"
        return out

    try:
        df = read_tsv_with_header(aln_path)
    except Exception as e:
        out["Note"] = f"read_error:{e}"
        return out

    required_cols = {"length1", "length2", "strand1", "strand2"}
    if not required_cols.issubset(set(df.columns)):
        out["Note"] = "missing_required_columns"
        return out

    df["inversion"] = df["strand1"] != df["strand2"]
    inv_df = df[df["inversion"]].copy()
    inv_df = inv_df.dropna(subset=["length1", "length2"])
    inv_df = inv_df.reset_index(drop=True)

    n_inv = len(inv_df)
    out["NumInversions"] = int(n_inv)
    if n_inv < 2:
        out["Note"] = "fewer_than_2_inversions"
        return out

    # compute b and inversion length
    inv_df["b"] = 0.5 * (inv_df["length1"] + inv_df["length2"])
    inv_df["inv_len"] = inv_df[["length1", "length2"]].min(axis=1)

    m = -1.0
    denom = math.sqrt(1 + m * m)  # = sqrt(2)

    pair_minimals = []
    for (i, rowA), (j, rowB) in combinations(inv_df.iterrows(), 2):
        smaller_len = min(rowA["inv_len"], rowB["inv_len"])
        dist = abs(rowA["b"] - rowB["b"]) / denom
        pair_minimals.append(min(dist, smaller_len))

    if len(pair_minimals) == 0:
        out["Note"] = "no_pairs_computed"
        return out

    out["MinimalUnit_bp"] = float(min(pair_minimals))
    out["Note"] = f"computed_from_{len(pair_minimals)}_pairs"
    return out

def main():
    parser = argparse.ArgumentParser(description="Summarize inversion stats from LASTZ alignments and BED files")
    parser.add_argument("-i", "--input_dir", required=True, help="Top-level input directory")
    parser.add_argument("-s", "--summary", required=True, help="Path to summary_table_IGH.tsv")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    parser.add_argument("-c", "--cores", type=int, default=4, help="Number of parallel workers")
    args = parser.parse_args()

    summary_df = pd.read_csv(args.summary, sep="\t", dtype=str)

    tasks = [ (row[1].to_dict(), args.input_dir) for row in summary_df.iterrows() ]

    with Pool(processes=args.cores) as pool:
        results = pool.starmap(process_row_dict, tasks)

    out_df = pd.DataFrame(results)
    out_df.to_csv(args.output, sep="\t", index=False)
    print("Wrote:", args.output)

if __name__ == "__main__":
    main()
