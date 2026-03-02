#!/usr/bin/env python3
import pandas as pd
import argparse
import os

def read_inversions(bed_path):
    """
    Reads the new combined_bed file (header or not) and returns a DataFrame.
    """
    # Read with pandas; detect header
    with open(bed_path) as f:
        first = f.readline().strip().split("\t")
    has_header = not all(c.replace("-", "").isdigit() for c in first[2:4])  # crude check

    df = pd.read_csv(bed_path, sep="\t", header=0 if has_header else None)
    if not has_header:
        df.columns = [
            "name1","strand1","start1","end1","length1",
            "name2","strand2","start2","end2","length2","id"
        ]

    # Ensure numeric columns
    for c in ["start1","end1","start2","end2"]:
        df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")

    # Extract contig name (assuming name1 == name2 per inversion)
    df["chr"] = df["name1"]

    # Filter out incomplete rows
    df = df.dropna(subset=["start1","end1","start2","end2"])
    return df


def compute_union_length(intervals):
    """Compute the total covered length of potentially overlapping intervals."""
    if len(intervals) == 0:
        return 0

    # sort by start
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = [intervals[0]]

    for current in intervals[1:]:
        last = merged[-1]
        if current[0] <= last[1]:  # overlap
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            merged.append(current)

    total = sum(e - s for s, e in merged)
    return total


def compute_inversion_density(bed_path, window_size, output):
    df = read_inversions(bed_path)
    summary = []

    for contig, group in df.groupby("name1"):
        max_pos = max(group["end1"].max(), group["end2"].max())
        #windows = sliding_windows(contig_len, args.window, step)
        for window_start in range(0, max_pos, window_size):
            window_end = window_start + window_size

            # filter inversions where both ends are within this window
            subset = group[
                (group["start1"] >= window_start) & (group["end1"] <= window_end) &
                (group["start2"] >= window_start) & (group["end2"] <= window_end)
            ]

            num_inv = len(subset)

            # collect all coordinates (both halves of inversion)
            all_intervals = (
                subset[["start1", "end1"]].values.tolist() +
                subset[["start2", "end2"]].values.tolist()
            )
            covered_len = compute_union_length(all_intervals)
            covered_fraction = covered_len / window_size
            center = (window_start + window_end) // 2
            summary.append({
                "chr": contig,
                "window_start": window_start,
                "window_end": window_end,
                "num_inversions": num_inv,
                "covered_fraction": covered_fraction,
                "center": center
            })

    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(output, sep="\t", index=False)
    print(f"✅ Saved summary to {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute inversion density in sliding windows (accounting for overlaps)"
    )
    parser.add_argument("--bed", required=True, help="Path to combined inversions BED")
    parser.add_argument("--window","-w", type=int, required=True, help="Window size in bp")
    parser.add_argument("--out", default="inversion_window_summary.tsv",
                        help="Output summary TSV")
    args = parser.parse_args()

    compute_inversion_density(args.bed, args.window, args.out)
