#!/usr/bin/env python3

import os
import argparse
import pandas as pd
from itertools import combinations
from collections import defaultdict

# --------------------------------------------------
# PARAMETERS
# --------------------------------------------------

MIN_INV_LEN = 500        # minimum inversion size (bp)
MIN_IDENTITY = 75.0      # minimum percent identity
MIN_RECIP_OVERLAP = 0.7  # reciprocal overlap threshold for clustering

# --------------------------------------------------
# LASTZ PARSER
# --------------------------------------------------

def parse_lastz_txt(txt_file):
    """
    Parse LASTZ txt output.
    Returns DataFrame with:
        ref_start, ref_end, strand, length, identity
    """
    rows = []

    with open(txt_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split()

            ref_start = int(parts[2])
            ref_end   = int(parts[3])
            strand    = str(parts[6])
            identity  = float(parts[10].strip("%"))

            length = abs(ref_end - ref_start)

            rows.append({
                "ref_start": min(ref_start, ref_end),
                "ref_end": max(ref_start, ref_end),
                "strand": strand,
                "length": length,
                "identity": identity
            })

    return pd.DataFrame(rows)


# --------------------------------------------------
# EXTRACT HIGH-CONFIDENCE INVERSIONS
# --------------------------------------------------

def extract_inversions(df):
    """
    Filter LASTZ alignments for inversions.
    """
    if df.empty:
        return df

    inv = df[
        (df["strand"] == "-") &
        (df["length"] >= MIN_INV_LEN) &
        (df["identity"] >= MIN_IDENTITY)
    ].copy()

    return inv


# --------------------------------------------------
# CLUSTERING FUNCTION
# --------------------------------------------------

def reciprocal_overlap(a_start, a_end, b_start, b_end):
    """
    Compute reciprocal overlap between two intervals.
    """
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start))
    if overlap == 0:
        return 0

    a_len = a_end - a_start
    b_len = b_end - b_start

    return min(overlap / a_len, overlap / b_len)


def cluster_intervals(intervals):
    """
    Cluster intervals using reciprocal overlap threshold.
    Input:
        intervals = list of dicts with keys:
            start, end, hap1, hap2
    Returns:
        list of clusters (each cluster = list of intervals)
    """

    clusters = []

    for inv in intervals:
        placed = False

        for cluster in clusters:
            for c in cluster:
                ro = reciprocal_overlap(
                    inv["start"], inv["end"],
                    c["start"], c["end"]
                )

                if ro >= MIN_RECIP_OVERLAP:
                    cluster.append(inv)
                    placed = True
                    break

            if placed:
                break

        if not placed:
            clusters.append([inv])

    return clusters


# --------------------------------------------------
# MAIN
# --------------------------------------------------

def main(args):

    summary = pd.read_csv(args.summary, sep="\t")
    all_presence = []

    for species, species_df in summary.groupby("Species"):

        print(f"\nProcessing species: {species}")

        haps = sorted(species_df["Haplotype"].unique())
        pairwise_inversions = []

        # --------------------------------------------------
        # ALL-VS-ALL PAIRWISE ALIGNMENTS
        # --------------------------------------------------

        for h1, h2 in combinations(haps, 2):

            txt1 = os.path.join(args.lastz_dir, f"{h1}_{h2}.txt")
            txt2 = os.path.join(args.lastz_dir, f"{h2}_{h1}.txt")

            if os.path.exists(txt1):
                txt = txt1
            elif os.path.exists(txt2):
                txt = txt2
            else:
                print(f"  Missing alignment: {h1} vs {h2}")
                continue

            df = parse_lastz_txt(txt)
            inv = extract_inversions(df)

            if inv.empty:
                continue

            for _, r in inv.iterrows():
                pairwise_inversions.append({
                    "start": r.ref_start,
                    "end": r.ref_end,
                    "hap1": h1,
                    "hap2": h2
                })

        if not pairwise_inversions:
            print("  No inversions detected.")
            continue

        # --------------------------------------------------
        # CLUSTER INVERSION EVENTS
        # --------------------------------------------------

        clusters = cluster_intervals(pairwise_inversions)

        print(f"  Detected {len(clusters)} inversion clusters")

        # --------------------------------------------------
        # BUILD PRESENCE MATRIX
        # --------------------------------------------------

        for i, cluster in enumerate(clusters):

            inv_id = f"{species}_inv_{i+1}"

            supporting_haps = set()

            for inv in cluster:
                supporting_haps.add(inv["hap1"])
                supporting_haps.add(inv["hap2"])

            for hap in haps:
                all_presence.append({
                    "Species": species,
                    "Inversion": inv_id,
                    "Haplotype": hap,
                    "Present": 1 if hap in supporting_haps else 0
                })

    # --------------------------------------------------
    # SAVE OUTPUT
    # --------------------------------------------------

    out = pd.DataFrame(all_presence)

    if out.empty:
        print("\nNo inversion presence data generated.")
    else:
        out.to_csv(args.output, sep="\t", index=False)
        print(f"\nSaved: {args.output}")


# --------------------------------------------------
# ARGPARSE
# --------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Reference-free inversion presence detection using all-vs-all LASTZ alignments."
    )

    parser.add_argument("-s", "--summary", required=True,
                        help="TSV with Species and Haplotype columns")

    parser.add_argument("-l", "--lastz_dir", required=True,
                        help="Directory containing LASTZ txt files")

    parser.add_argument("-o", "--output",
                        default="inversion_presence_all_vs_all.tsv",
                        help="Output TSV file")

    args = parser.parse_args()

    main(args)