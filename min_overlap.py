#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import re
from collections import defaultdict
import numpy as np
from multiprocessing import Pool
from itertools import combinations

def read_tsv_with_header(path):
    """Read TSV even if header starts with #"""
    with open(path) as f:
        header = f.readline().strip()
        if header.startswith("#"):
            header = header[1:]  # remove leading #
        columns = header.split("\t")
    df = pd.read_csv(path, sep="\t", comment="#", names=columns, skiprows=1)
    df.columns = [c.replace("#", "").replace("%", "").replace("+", "") for c in df.columns]
    return df


def calc_overlap(a_start, a_end, b_start, b_end):
    """Return overlap length (0 if none)."""
    overlap = min(a_end, b_end) - max(a_start, b_start)
    return max(0, overlap)


def process_row(row):
    input_dir = row['InputDir']
    order = row['Order']
    species = row['Species']
    haplotype = row['Haplotype']

    aln_path = os.path.join(input_dir, order, species, haplotype, 'IGH_self.tsv')
    bed_path = os.path.join(input_dir, order, species, haplotype, 'IGH.bed')

    if not os.path.exists(aln_path) or not os.path.exists(bed_path):
        print(f"Skipping {order}/{species}/{haplotype}: missing files")
        return []

    try:
        df = read_tsv_with_header(aln_path)
        bed_df = pd.read_csv(bed_path, sep='\t', header=None)
        bed_df.columns = ['contig', 'start', 'end', 'na1', 'na2', 'strand', 'na3', 'na4', 'color']
    except Exception as e:
        print(f"Error reading {order}/{species}/{haplotype}: {e}")
        return []

    # Only consider inversions (opposite strands)
    inv = df[df['strand1'] != df['strand2']].copy()
    if inv.empty or len(inv) < 2:
        return []

    min_overlap = None
    min_info = None

    # Compare all pairs of inversions
    for (i, row1), (j, row2) in combinations(inv.iterrows(), 2):
        # Overlap in first coordinate space
        ov1 = calc_overlap(row1['start1'], row1['end1'], row2['start1'], row2['end1'])
        # Overlap in second coordinate space
        ov2 = calc_overlap(row1['start2'], row1['end2'], row2['start2'], row2['end2'])
        # If any overlap occurs
        # Check first coordinate space
        if ov1 > 20:
            if min_overlap is None or ov1 < min_overlap:
                min_overlap = ov1
                min_info = {
                    'Order': order,
                    'Species': species,
                    'Haplotype': haplotype,
                    'index1': i,
                    'index2': j,
                    'overlap_len': ov1,
                    'overlap_type': 'start1-end1',
                    'region1_start': row1['start1'],
                    'region1_end': row1['end1'],
                    'region2_start': row2['start1'],
                    'region2_end': row2['end1']
                }

        # Check second coordinate space separately
        if ov2 > 20:
            if min_overlap is None or ov2 < min_overlap:
                min_overlap = ov2
                min_info = {
                    'Order': order,
                    'Species': species,
                    'Haplotype': haplotype,
                    'index1': i,
                    'index2': j,
                    'overlap_len': ov2,
                    'overlap_type': 'start2-end2',
                    'region1_start': row1['start2'],
                    'region1_end': row1['end2'],
                    'region2_start': row2['start2'],
                    'region2_end': row2['end2']
                }

    print(f"Processed: {order}/{species}/{haplotype}")
    if min_info is None:
        return []
    return [min_info]


def main():
    parser = argparse.ArgumentParser(description="Calculate minimal overlap between inversions in LASTZ alignments")
    parser.add_argument('-i', '--input_dir', required=True, help="Top-level input directory")
    parser.add_argument('-s', '--summary', required=True, help="Path to summary_table_IGH.tsv")
    parser.add_argument('-o', '--output', required=True, help="Output TSV file")
    parser.add_argument('-c', '--cores', type=int, default=4, help="Number of parallel workers")

    args = parser.parse_args()
    df = pd.read_csv(args.summary, sep='\t')
    df['InputDir'] = args.input_dir

    with Pool(args.cores) as pool:
        results = pool.map(process_row, [row for _, row in df.iterrows()])

    # Flatten results and save
    flat = [item for sublist in results for item in sublist if sublist]
    if not flat:
        print("No overlaps found in any species.")
        return

    out_df = pd.DataFrame(flat)
    out_df.to_csv(args.output, sep='\t', index=False)
    print(f"Saved results to {args.output}")


if __name__ == "__main__":
    main()