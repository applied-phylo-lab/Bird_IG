#!/usr/bin/env python3

import os
import argparse
import pandas as pd
from itertools import combinations_with_replacement, combinations
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
    if not intervals:
        return []
    
    # Sort by start position
    sorted_inv = sorted(intervals, key=lambda x: x["start"])
    
    clusters = []
    current_cluster = [sorted_inv[0]]
    current_start = sorted_inv[0]["start"]
    current_end = sorted_inv[0]["end"]
    
    for inv in sorted_inv[1:]:
        ro = reciprocal_overlap(current_start, current_end, inv["start"], inv["end"])
        if ro >= MIN_RECIP_OVERLAP:
            current_cluster.append(inv)
            current_end = max(current_end, inv["end"])
        else:
            clusters.append(current_cluster)
            current_cluster = [inv]
            current_start = inv["start"]
            current_end = inv["end"]
    
    clusters.append(current_cluster)
    return clusters

def process_pairs_batch(batch_args):
    results = []
    for args in batch_args:
        results.extend(process_pair(args))
    return results


# ----------------------------
# PARALLEL LASTZ PROCESSING
# ----------------------------
def process_pair(args):
    k1, k2, lastz_dir = args
    is_self = (k1 == k2)

    txt1 = os.path.join(lastz_dir, f"{k1}_vs_{k2}.txt")
    txt2 = os.path.join(lastz_dir, f"{k2}_vs_{k1}.txt")

    if os.path.exists(txt1):
        txt = txt1
    elif os.path.exists(txt2):
        txt = txt2
    else:
        return []

    df = parse_lastz_txt(txt)
    inv = extract_inversions(df)
    out = []
    for _, r in inv.iterrows():
        out.append({
            "start": r.ref_start,
            "end": r.ref_end,
            "length": r.length,
            "key1": k1,
            "key2": k2,
            "is_self": is_self
        })
    return out

# ----------------------------
# RESOLVE WHICH KEYS SUPPORT AN INVERSION CLUSTER
# ----------------------------
def resolve_supporting_keys(cluster, self_inversion_intervals):
    """
    For each inversion cluster, determine which keys truly carry the inversion.

    - Self-alignment hits: the key definitively has the inversion.
    - Cross-alignment hits: check which of the two keys has a self-alignment
      inversion overlapping this cluster. If one does -> that key carries it.
      If both do -> both carry it. If neither does -> can't resolve, mark both.
    """
    supporting = set()
    lengths, starts, ends = [], [], []

    for inv in cluster:
        lengths.append(inv["length"])
        starts.append(inv["start"])
        ends.append(inv["end"])

        if inv["is_self"]:
            supporting.add(inv["key1"])
        else:
            k1_has_self = any(
                reciprocal_overlap(inv["start"], inv["end"], s["start"], s["end"]) >= MIN_RECIP_OVERLAP
                for s in self_inversion_intervals.get(inv["key1"], [])
            )
            k2_has_self = any(
                reciprocal_overlap(inv["start"], inv["end"], s["start"], s["end"]) >= MIN_RECIP_OVERLAP
                for s in self_inversion_intervals.get(inv["key2"], [])
            )

            if k1_has_self or k2_has_self:
                if k1_has_self:
                    supporting.add(inv["key1"])
                if k2_has_self:
                    supporting.add(inv["key2"])
            else:
                pass

    return supporting, lengths, starts, ends

# ----------------------------
# MAIN
# ----------------------------
def main(args):

    summary = pd.read_csv(args.summary, sep="\t")
    summary["Key"] = summary["Haplotype"] + "_" + summary["Contig"].astype(str)
    key_to_hap = dict(zip(summary["Key"], summary["Haplotype"]))

    all_presence = []
    all_stats = []

    for species, species_df in summary.groupby("Species"):

        print(f"\nProcessing {species}")

        keys = sorted(species_df["Key"].unique())
        haps = sorted(species_df["Haplotype"].unique())
        reference_key = keys[0]
        print(f"  Keys: {keys}")
        print(f"  Reporting coordinates relative to: {reference_key}")

        # All pairs including self-alignments
        pair_args = [(k1, k2, args.lastz_dir) for k1, k2 in combinations_with_replacement(keys, 2)]

        pairwise_inversions = []
        # Then batch into chunks
        chunk_size = max(1, len(pair_args) // (args.cores * 4))
        batches = [pair_args[i:i+chunk_size] for i in range(0, len(pair_args), chunk_size)]

        with ProcessPoolExecutor(max_workers=args.cores) as executor:
            futures = [executor.submit(process_pairs_batch, batch) for batch in batches]
            for future in as_completed(futures):
                pairwise_inversions.extend(future.result())

        if not pairwise_inversions:
            print("  No inversions detected.")
            continue

        # Separate self-alignment inversions for use in resolution
        self_inversion_intervals = {}
        for inv in pairwise_inversions:
            if inv["is_self"]:
                self_inversion_intervals.setdefault(inv["key1"], []).append(inv)

        # Cluster ALL inversions together (self + cross)
        clusters = cluster_intervals(pairwise_inversions)
        print(f"  Detected {len(clusters)} inversion clusters")

        for i, cluster in enumerate(clusters):
            inv_id = f"{species}_inv_{i+1}"

            supporting_keys, lengths, starts, ends = resolve_supporting_keys(
                cluster, self_inversion_intervals
            )
            supporting_haps = set(key_to_hap[k] for k in supporting_keys)

            mean_len = np.mean(lengths)
            mean_start = int(np.mean(starts))
            mean_end = int(np.mean(ends))
            freq = len(supporting_haps) / len(haps)

            all_stats.append({
                "Species": species,
                "Inversion": inv_id,
                "RefKey": reference_key,
                "RefStart": mean_start,
                "RefEnd": mean_end,
                "MeanLength": round(mean_len, 1),
                "SupportingKeys": ",".join(sorted(supporting_keys)),
                "NumSupportingHaps": len(supporting_haps),
                "TotalHaps": len(haps),
                "Frequency": round(freq, 3)
            })

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