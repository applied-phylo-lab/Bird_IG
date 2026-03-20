#!/usr/bin/env python3

import os
import argparse
import pandas as pd
from itertools import combinations
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np

MIN_INV_LEN = 500
MIN_IDENTITY = 75.0
MIN_RECIP_OVERLAP = 0.7

# ----------------------------
# LASTZ PARSING
# ----------------------------
def parse_lastz_txt(txt_file):
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

def extract_inversions(df):
    if df.empty:
        return df
    return df[(df["strand"] == "-") &
              (df["length"] >= MIN_INV_LEN) &
              (df["identity"] >= MIN_IDENTITY)].copy()

def reciprocal_overlap(a_start, a_end, b_start, b_end):
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start))
    if overlap == 0:
        return 0
    a_len = a_end - a_start
    b_len = b_end - b_start
    return min(overlap / a_len, overlap / b_len)

def cluster_intervals(intervals):
    clusters = []
    for inv in intervals:
        placed = False
        for cluster in clusters:
            for c in cluster:
                ro = reciprocal_overlap(inv["start"], inv["end"], c["start"], c["end"])
                if ro >= MIN_RECIP_OVERLAP:
                    cluster.append(inv)
                    placed = True
                    break
            if placed:
                break
        if not placed:
            clusters.append([inv])
    return clusters

# ----------------------------
# PARALLEL LASTZ PROCESSING
# ----------------------------
def process_pair(args):
    h1, h2, lastz_dir = args
    txt1 = os.path.join(lastz_dir, f"{h1}_{h2}.txt")
    txt2 = os.path.join(lastz_dir, f"{h2}_{h1}.txt")

    if os.path.exists(txt1):
        txt = txt1
        anchor = h1
    elif os.path.exists(txt2):
        txt = txt2
        anchor = h2
    else:
        return []  # missing alignment

    df = parse_lastz_txt(txt)
    inv = extract_inversions(df)
    out = []
    for _, r in inv.iterrows():
        out.append({
            "start": r.ref_start,
            "end": r.ref_end,
            "length": r.length,
            "hap1": h1,
            "hap2": h2,
            "anchor": anchor
        })
    return out

# ----------------------------
# MAIN
# ----------------------------
def main(args):

    summary = pd.read_csv(args.summary, sep="\t")

    all_presence = []
    all_stats = []

    for species, species_df in summary.groupby("Species"):

        print(f"\nProcessing {species}")
        haps = sorted(species_df["Haplotype"].unique())
        reference_hap = haps[0]
        print(f"  Reporting coordinates relative to: {reference_hap}")

        # prepare pairwise arguments
        pair_args = [(h1, h2, args.lastz_dir) for h1, h2 in combinations(haps, 2)]

        pairwise_inversions = []

        # ----------------------------
        # PARALLEL EXECUTION
        # ----------------------------
        with ProcessPoolExecutor(max_workers=args.cores) as executor:
            futures = {executor.submit(process_pair, pa): pa for pa in pair_args}
            for future in as_completed(futures):
                result = future.result()
                pairwise_inversions.extend(result)

        if not pairwise_inversions:
            print("  No inversions detected.")
            continue

        clusters = cluster_intervals(pairwise_inversions)
        print(f"  Detected {len(clusters)} inversion clusters")

        for i, cluster in enumerate(clusters):
            inv_id = f"{species}_inv_{i+1}"
            supporting_haps = set()
            lengths, starts, ends = [], [], []

            for inv in cluster:
                supporting_haps.add(inv["hap1"])
                supporting_haps.add(inv["hap2"])
                lengths.append(inv["length"])
                starts.append(inv["start"])
                ends.append(inv["end"])

            mean_len = np.mean(lengths)
            mean_start = int(np.mean(starts))
            mean_end = int(np.mean(ends))
            freq = len(supporting_haps) / len(haps)

            # stats
            all_stats.append({
                "Species": species,
                "Inversion": inv_id,
                "RefHaplotype": reference_hap,
                "RefStart": mean_start,
                "RefEnd": mean_end,
                "MeanLength": round(mean_len, 1),
                "NumSupportingHaps": len(supporting_haps),
                "TotalHaps": len(haps),
                "Frequency": round(freq,3)
            })

            # presence
            for hap in haps:
                all_presence.append({
                    "Species": species,
                    "Inversion": inv_id,
                    "Haplotype": hap,
                    "Present": 1 if hap in supporting_haps else 0
                })

    pd.DataFrame(all_presence).to_csv(args.output_prefix + "_presence.tsv", sep="\t", index=False)
    pd.DataFrame(all_stats).to_csv(args.output_prefix + "_stats.tsv", sep="\t", index=False)
    print("\nDone.")

# ----------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--summary", required=True)
    parser.add_argument("-l", "--lastz_dir", required=True)
    parser.add_argument("-o", "--output_prefix", default="inversion")
    parser.add_argument("-c", "--cores", type=int, default=4, help="Number of cores for parallel processing")
    args = parser.parse_args()

    main(args)