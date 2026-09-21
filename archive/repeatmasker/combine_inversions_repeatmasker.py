#!/usr/bin/env python3

import os
import argparse
import pandas as pd

# -----------------------------
# Filters
# -----------------------------
MIN_INV_LEN = 500        # bp
MIN_IDENTITY = 75.0      # percent

# -----------------------------
def read_lastz_self(tsv):
    """
    Read LASTZ self-alignment TSV with header and extract inversions.
    """
    df = pd.read_csv(tsv, sep="\t", comment="#", header=0)

    # Fix column names if needed
    #df.columns = ["name1","strand1","start1","end1","length1","name2","strand2","start2","end2","length2","id"]
    df.columns = ["name1","strand1","start1","end1","length1","text1","name2","strand2","start2","end2","length2","text2","id"]
    df["length1"]= df["length1"].astype(int)
    #remove % from identity
    df["id"] = df["id"].str.rstrip('%').astype(float)
    inv = df[
        (df["strand2"] == "-") &
        (df["length1"] >= MIN_INV_LEN) &
        (df["id"] >= MIN_IDENTITY)
    ].copy()
    
    if inv.empty:
        return inv

    inv["InvStart"] = inv[["start1","end1"]].min(axis=1)
    inv["InvEnd"]   = inv[["start1","end1"]].max(axis=1)
    inv["InvLength"] = inv["InvEnd"] - inv["InvStart"]

    return inv[["name1","InvStart","InvEnd","InvLength"]]

# -----------------------------
def intersect_inversions_repeats(inv_df, rep_df, hap):
    rows = []

    for _, inv in inv_df.iterrows():
        i_start = inv["InvStart"]
        i_end   = inv["InvEnd"]

        overlaps = rep_df[
            (rep_df["end"] > i_start) &
            (rep_df["start"] < i_end)
        ]

        for _, r in overlaps.iterrows():
            ov = min(i_end, r["end"]) - max(i_start, r["start"])
            if ov <= 0:
                continue

            rows.append({
                "Haplotype": hap,
                "Contig": inv["name1"],
                "InvStart": i_start,
                "InvEnd": i_end,
                "InvLength": inv["InvLength"],
                "RepeatStart": r["start"],
                "RepeatEnd": r["end"],
                "RepeatName": r["name"],
                "RepeatClass": r["class_family"],
                "OverlapBp": ov
            })

    return rows

# -----------------------------
def main(args):
    summary = pd.read_csv(args.summary_tsv, sep="\t")

    for _, row in summary.iterrows():
        order = row["Order"]
        species = row["Species"]
        hap = row["Haplotype"]

        hap_dir = os.path.join(args.input_dir, order, species, hap,"RepeatMasker")
        lastz = os.path.join(args.input_dir, order, species, hap, "IGH_self.tsv")

        if not os.path.exists(lastz):
            print(f"[skip] {hap}: no IGH_self.tsv")
            continue

        beds = [f for f in os.listdir(hap_dir) if f.endswith(".repeats.bed")]
        if not beds:
            print(f"[skip] {hap}: no repeats BED")
            continue

        bed = os.path.join(hap_dir, beds[0])

        print(f"[{hap}] processing")

        inv_df = read_lastz_self(lastz)
        if inv_df.empty:
            print("  no inversions found")
            continue

        rep_df = pd.read_csv(
            bed, sep="\t",
            names=["chrom","start","end","name","score","strand","class_family"]
        )

        rows = intersect_inversions_repeats(inv_df, rep_df, hap)

        if not rows:
            print("  no repeat overlaps")
            continue

        out_df = pd.DataFrame(rows)
        out_file = os.path.join(hap_dir, "inversion_repeats.tsv")
        out_df.to_csv(out_file, sep="\t", index=False)

        print(f"  wrote {out_file}")

# -----------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Intersect inversions with RepeatMasker repeats for all haplotypes"
    )
    parser.add_argument(
        "-s","--summary_tsv", required=True,
        help="Summary TSV listing all haplotypes"
    )
    parser.add_argument(
        "-i","--input_dir", required=True,
        help="Base input directory (Order/Species/Haplotype)"
    )

    args = parser.parse_args()
    main(args)
