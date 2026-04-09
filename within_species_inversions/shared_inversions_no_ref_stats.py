#!/usr/bin/env python3

import os
import argparse
import pandas as pd
from itertools import combinations_with_replacement
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
            ref_start   = int(parts[2])
            ref_end     = int(parts[3])
            strand      = str(parts[6])
            query_start = int(parts[7])
            query_end   = int(parts[8])
            identity    = float(parts[10].strip("%"))
            length      = abs(ref_end - ref_start)
            rows.append({
                "ref_start":   min(ref_start, ref_end),
                "ref_end":     max(ref_start, ref_end),
                "query_start": min(query_start, query_end),
                "query_end":   max(query_start, query_end),
                "strand":      strand,
                "length":      length,
                "identity":    identity
            })
    return pd.DataFrame(rows)


def extract_inversions(df):
    if df.empty:
        return df
    return df[
        (df["strand"]   == "-") &
        (df["length"]   >= MIN_INV_LEN) &
        (df["identity"] >= MIN_IDENTITY)
    ].copy()


def reciprocal_overlap(a_start, a_end, b_start, b_end):
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start))
    if overlap == 0:
        return 0.0
    a_len = a_end - a_start
    b_len = b_end - b_start
    return min(overlap / a_len, overlap / b_len)


# ----------------------------
# CLUSTERING  (O(n log n))
# ----------------------------
def cluster_intervals(intervals):
    if not intervals:
        return []
    sorted_inv      = sorted(intervals, key=lambda x: x["start"])
    clusters        = []
    current_cluster = [sorted_inv[0]]
    current_start   = sorted_inv[0]["start"]
    current_end     = sorted_inv[0]["end"]

    for inv in sorted_inv[1:]:
        ro = reciprocal_overlap(current_start, current_end, inv["start"], inv["end"])
        if ro >= MIN_RECIP_OVERLAP:
            current_cluster.append(inv)
            current_end = max(current_end, inv["end"])
        else:
            clusters.append(current_cluster)
            current_cluster = [inv]
            current_start   = inv["start"]
            current_end     = inv["end"]

    clusters.append(current_cluster)
    return clusters


# ----------------------------
# PARALLEL LASTZ PROCESSING
# ----------------------------
def process_pair(k1, k2, lastz_dir):
    is_self = (k1 == k2)
    txt1 = os.path.join(lastz_dir, f"{k1}_vs_{k2}.txt")
    txt2 = os.path.join(lastz_dir, f"{k2}_vs_{k1}.txt")

    if os.path.exists(txt1):
        txt = txt1
    elif os.path.exists(txt2):
        txt = txt2
    else:
        return []

    df  = parse_lastz_txt(txt)
    inv = extract_inversions(df)

    inversion_rows = []
    for _, r in inv.iterrows():
        inversion_rows.append({
            "start":       r.ref_start,
            "end":         r.ref_end,
            "query_start": r.query_start,
            "query_end":   r.query_end,
            "length":      r.length,
            "identity":    r.identity,
            "key1":        k1,
            "key2":        k2,
            "is_self":     is_self
        })
    return inversion_rows


def process_pairs_batch(batch_args):
    results = []
    for k1, k2, lastz_dir in batch_args:
        results.extend(process_pair(k1, k2, lastz_dir))
    return results


# ----------------------------
# RESOLVE WHICH KEYS SUPPORT A CLUSTER
# ----------------------------
def resolve_supporting_keys(cluster, self_inversion_intervals):
    supporting = set()
    for inv in cluster:
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
            if k1_has_self:
                supporting.add(inv["key1"])
            if k2_has_self:
                supporting.add(inv["key2"])
    return supporting


# ----------------------------
# PICK REPRESENTATIVE COORDS FOR DOTPLOT
# Priority: reference key self-alignment > any supporting key self-alignment
# Returns (ref_start, ref_end, query_start, query_end, rep_key)
# ----------------------------
def pick_representative(cluster, supporting_keys, reference_key):
    # Prefer reference key self-alignment block
    for inv in cluster:
        if inv["is_self"] and inv["key1"] == reference_key and reference_key in supporting_keys:
            return inv["start"], inv["end"], inv["query_start"], inv["query_end"], reference_key

    # Fall back to first supporting key alphabetically
    for key in sorted(supporting_keys):
        for inv in cluster:
            if inv["is_self"] and inv["key1"] == key:
                return inv["start"], inv["end"], inv["query_start"], inv["query_end"], key

    # Last resort: use cross-alignment coords involving reference
    for inv in cluster:
        if inv["key1"] == reference_key and inv["key2"] in supporting_keys:
            return inv["start"], inv["end"], inv["query_start"], inv["query_end"], inv["key2"]
        if inv["key2"] == reference_key and inv["key1"] in supporting_keys:
            return inv["query_start"], inv["query_end"], inv["start"], inv["end"], inv["key1"]

    return None


# ----------------------------
# MAIN
# ----------------------------
def main(args):

    summary = pd.read_csv(args.summary, sep="\t")
    summary["Key"] = summary["Haplotype"] + "_" + summary["Contig"].astype(str)
    key_to_hap = dict(zip(summary["Key"], summary["Haplotype"]))

    all_presence = []
    all_stats    = []
    all_coords   = []
    all_dotplot  = []

    for species, species_df in summary.groupby("Species"):

        print(f"\nProcessing {species}")

        keys = sorted(species_df["Key"].unique())
        haps = sorted(species_df["Haplotype"].unique())
        reference_key = keys[0]
        print(f"  Keys ({len(keys)}): {keys}")
        print(f"  Reference key: {reference_key}")

        pair_args = [
            (k1, k2, args.lastz_dir)
            for k1, k2 in combinations_with_replacement(keys, 2)
        ]

        chunk_size = max(1, len(pair_args) // (args.cores * 4))
        batches = [
            pair_args[i:i + chunk_size]
            for i in range(0, len(pair_args), chunk_size)
        ]

        pairwise_inversions = []
        with ProcessPoolExecutor(max_workers=args.cores) as executor:
            futures = [executor.submit(process_pairs_batch, b) for b in batches]
            for future in as_completed(futures):
                pairwise_inversions.extend(future.result())

        if not pairwise_inversions:
            print("  No inversions detected.")
            continue

        self_inversion_intervals = {}
        for inv in pairwise_inversions:
            if inv["is_self"]:
                self_inversion_intervals.setdefault(inv["key1"], []).append(inv)

        clusters = cluster_intervals(pairwise_inversions)
        print(f"  Detected {len(clusters)} inversion clusters")

        for i, cluster in enumerate(clusters):
            inv_id = f"{species}_inv_{i+1}"

            supporting_keys = resolve_supporting_keys(cluster, self_inversion_intervals)
            if not supporting_keys:
                continue

            supporting_haps = set(key_to_hap[k] for k in supporting_keys)
            freq = len(supporting_haps) / len(haps)

            key_coords = {}
            for inv in cluster:
                if inv["is_self"] and inv["key1"] in supporting_keys:
                    key_coords.setdefault(inv["key1"], []).append(inv)

            all_lengths = [inv["length"] for inv in cluster]

            all_stats.append({
                "Species":           species,
                "Inversion":         inv_id,
                "MeanLength":        round(np.mean(all_lengths), 1),
                "MinLength":         np.min(all_lengths),
                "MaxLength":         np.max(all_lengths),
                "SupportingKeys":    ",".join(sorted(supporting_keys)),
                "NumSupportingHaps": len(supporting_haps),
                "TotalHaps":         len(haps),
                "Frequency":         round(freq, 3)
            })

            for hap in haps:
                all_presence.append({
                    "Species":   species,
                    "Inversion": inv_id,
                    "Haplotype": hap,
                    "Present":   1 if hap in supporting_haps else 0
                })

            for key, inv_list in key_coords.items():
                hap = key_to_hap[key]
                for inv in inv_list:
                    all_coords.append({
                        "Species":   species,
                        "Inversion": inv_id,
                        "Haplotype": hap,
                        "Key":       key,
                        "Start":     inv["start"],
                        "End":       inv["end"],
                        "Length":    inv["length"]
                    })

            # --- Dotplot: one representative row per inversion ---
            rep = pick_representative(cluster, supporting_keys, reference_key)
            if rep is None:
                continue
            ref_start, ref_end, qry_start, qry_end, rep_key = rep

            all_dotplot.append({
                "Species":           species,
                "Inversion":         inv_id,
                "RepKey":            rep_key,
                "RefStart":          ref_start,
                "RefEnd":            ref_end,
                "QueryStart":        qry_start,
                "QueryEnd":          qry_end,
                "MeanLength":        round(np.mean(all_lengths), 1),
                "NumSupportingHaps": len(supporting_haps),
                "TotalHaps":         len(haps),
                "Frequency":         round(freq, 3)
            })

    pd.DataFrame(all_presence).to_csv(
        args.output_prefix + "_presence.tsv", sep="\t", index=False)
    pd.DataFrame(all_stats).to_csv(
        args.output_prefix + "_stats.tsv", sep="\t", index=False)
    pd.DataFrame(all_coords).to_csv(
        args.output_prefix + "_coords.tsv", sep="\t", index=False)
    pd.DataFrame(all_dotplot).to_csv(
        args.output_prefix + "_dotplot.tsv", sep="\t", index=False)

    print("\nDone.")
    print(f"  Outputs: {args.output_prefix}_{{presence,stats,coords,dotplot}}.tsv")


# ----------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--summary",       required=True)
    parser.add_argument("-l", "--lastz_dir",     required=True)
    parser.add_argument("-o", "--output_prefix", default="inversion")
    parser.add_argument("-c", "--cores",         type=int, default=4)
    args = parser.parse_args()

    main(args)