#!/usr/bin/env python3

import os
import argparse
import pandas as pd
from collections import defaultdict
from intervaltree import Interval, IntervalTree

# -----------------------------
# Parameters (tweak as needed)
# -----------------------------
MIN_INV_LEN = 500        # bp
MIN_IDENTITY = 75.0       # percent
MIN_RECIP_OVERLAP = 0.7   # for clustering
REFERENCE_OVERRIDES = {
    "house_finches": "NY_1_hap1",
    "A.coerulescensAC": "A.coerulescens_AC_1603_72872_pri",
    "A.insularisAI": "A.insularis_AI_1363_74563_pri",
    "A.woodhouseiiAW": "A.woodhouseii_AW_366499_pri",
}
# -----------------------------
# LASTZ parser
# -----------------------------
def parse_lastz_txt(txt_file):
    """
    Parse LASTZ txt file.
    Returns DataFrame with ref_start, ref_end, strand, length, identity
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

# -----------------------------
# Extract inversions
# -----------------------------
def extract_inversions(df, label=None):
    """
    Extract inversion blocks from LASTZ dataframe.
    Handles empty dataframes gracefully.
    """
    if df.empty:
        msg = "  empty LASTZ alignment"
        if label:
            msg += f" ({label})"
        print(msg)
        return df

    inv = df[
        (df["strand"] == "-") &
        (df["length"] >= MIN_INV_LEN) &
        (df["identity"] >= MIN_IDENTITY)
    ].copy()

    if inv.empty:
        msg = "  no inversions after filtering"
        if label:
            msg += f" ({label})"
        print(msg)

    return inv



# -----------------------------
# Main
# -----------------------------
def main(args):
    summary = pd.read_csv(args.summary, sep="\t")

    all_presence = []

    for species, species_df in summary.groupby("Species"):

        species_df = species_df.reset_index(drop=True)

        if species in REFERENCE_OVERRIDES:
            ref_hap = REFERENCE_OVERRIDES[species]

            if ref_hap not in set(species_df["Haplotype"]):
                raise ValueError(
                    f"Reference haplotype {ref_hap} not found for species {species}"
                )
        else:
            ref_hap = species_df.iloc[0]["Haplotype"]


        print(f"[{species}] reference haplotype: {ref_hap}")

        inv_by_hap = defaultdict(list)

        for hap in species_df["Haplotype"]:
           
            if hap == ref_hap:
                lastz_prefix = os.path.join(
                args.lastz_dir,
                f"{ref_hap}_self"
            )
            else:
                lastz_prefix = os.path.join(
                args.lastz_dir,
                f"{ref_hap}_{hap}"
                )

            txt = lastz_prefix + ".txt"

            if not os.path.exists(txt):
                lastz_prefix = os.path.join(
                args.lastz_dir,
                f"{hap}_{ref_hap}"
                )
                txt = lastz_prefix + ".txt"

            if not os.path.exists(txt):
                print(f"  missing: {txt}")
                continue

            df = parse_lastz_txt(txt)
            inv = extract_inversions(df, label=f"{ref_hap} vs {hap}")

            for _, r in inv.iterrows():
                inv_by_hap[hap].append((r.ref_start, r.ref_end))

        ref_comparison_hap = next(iter(inv_by_hap.keys()), None)

        if ref_comparison_hap is None:
             print(f"  {species}: no inversions found")
             continue

        reference_inversions = inv_by_hap[ref_comparison_hap]

        OVERLAP_FRAC = 0.8

        for i, (rs, re) in enumerate(reference_inversions):
            inv_id = f"{species}_inv_{i+1}"

            inv_len = re - rs

            for hap in species_df["Haplotype"]:
                present = 0

                for s, e in inv_by_hap.get(hap, []):
                    overlap = max(0, min(e, re) - max(s, rs))
                    if overlap / inv_len >= OVERLAP_FRAC:
                        present = 1
                        break

                all_presence.append({
                    "Species": species,
                    "Inversion": inv_id,
                    "RefStart": rs,
                    "RefEnd": re,
                    "Haplotype": hap,
                    "Present": present
                })


    out = pd.DataFrame(all_presence)
    out.to_csv(args.output, sep="\t", index=False)
    print(f"\nSaved: {args.output}")

# -----------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s","--summary", required=True)
    parser.add_argument("-l","--lastz_dir", required=True)
    parser.add_argument("-o","--output", default="inversion_presence.tsv")
    args = parser.parse_args()

    main(args)
