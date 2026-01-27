#!/usr/bin/env python3

import os
import argparse
import pandas as pd
from collections import defaultdict
import subprocess

# -----------------------------
# Parameters
# -----------------------------
MIN_INV_LEN = 500
MIN_IDENTITY = 75.0
OVERLAP_FRAC = 0.8

SUBGROUP_REF_OVERRIDE = {
    "house_finches_W": "bHaeMex1_pri"
}

STATE_REF_OVERRIDE = {
    "CA": "bHaeMex1_pri"
}

LASTZ_FLAGS = [
    "--step=20",
    "--notransition",
    "--format=general:name1,strand1,start1,end1,length1,"
    "name2,strand2,start2+,end2+,length2,id%"
]

# -----------------------------
# FASTA path helper
# -----------------------------
def get_fasta_path(input_dir, row):
    return os.path.join(
        input_dir,
        row["Order"],
        row["Species"],
        row["Haplotype"],
        "refined_ig_loci",
        "igloci_fasta",
        f"IGH_{row['Contig']}_{row['NumV']}Vs.fasta"
    )

# -----------------------------
# LASTZ helpers
# -----------------------------
def run_lastz(fa1, fa2, out_txt):
    cmd = ["lastz", fa1, fa2] + LASTZ_FLAGS
    with open(out_txt, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)

def ensure_self_alignment(ref_row, input_dir, aln_dir):
    hap = ref_row["Haplotype"]
    txt = os.path.join(aln_dir, f"{hap}_self.txt")
    if os.path.exists(txt):
        return txt

    fa = get_fasta_path(input_dir, ref_row)
    print(f"  running self-alignment for {hap}")
    run_lastz(fa, fa, txt)
    return txt

def ensure_pairwise(ref_row, hap_row, input_dir, aln_dir):
    ref = ref_row["Haplotype"]
    hap = hap_row["Haplotype"]

    txt1 = os.path.join(aln_dir, f"{ref}_{hap}.txt")
    txt2 = os.path.join(aln_dir, f"{hap}_{ref}.txt")

    if os.path.exists(txt1):
        return txt1, False
    if os.path.exists(txt2):
        return txt2, True

    print(f"  missing pairwise alignment, running LASTZ: {ref} vs {hap}")
    fa1 = get_fasta_path(input_dir, ref_row)
    fa2 = get_fasta_path(input_dir, hap_row)
    run_lastz(fa1, fa2, txt1)
    return txt1, False

# -----------------------------
# LASTZ parser
# -----------------------------
def parse_lastz_txt(txt_file, reversed_pair=False):
    rows = []
    with open(txt_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()

            if not reversed_pair:
                start = int(parts[2])
                end   = int(parts[3])
                strand = parts[6]
            else:
                start = int(parts[7])
                end   = int(parts[8])
                strand = parts[6]

            rows.append({
                "ref_start": min(start, end),
                "ref_end": max(start, end),
                "strand": strand,
                "length": abs(end - start),
                "identity": float(parts[10].strip("%"))
            })

    return pd.DataFrame(rows)

# -----------------------------
# Extract inversions
# -----------------------------
def extract_inversions(df, label=None):
    if df.empty:
        if label:
            print(f"  empty LASTZ alignment ({label})")
        return df

    inv = df[
        (df["strand"] == "-") &
        (df["length"] >= MIN_INV_LEN) &
        (df["identity"] >= MIN_IDENTITY)
    ].copy()

    if inv.empty and label:
        print(f"  no inversions after filtering ({label})")

    return inv

# -----------------------------
# Core logic
# -----------------------------
def process_group(df, group_col, ref_override, input_dir, aln_dir):
    rows = []

    for group, gdf in df.groupby(group_col):
        gdf = gdf.reset_index(drop=True)

        if group in ref_override:
            ref_hap = ref_override[group]
            ref_row = gdf[gdf["Haplotype"] == ref_hap].iloc[0]
        else:
            ref_row = gdf.iloc[0]
            ref_hap = ref_row["Haplotype"]

        print(f"[{group_col}={group}] reference haplotype: {ref_hap}")

        inv_by_hap = defaultdict(list)

        # Self-alignment
        self_txt = ensure_self_alignment(ref_row, input_dir, aln_dir)
        df_self = parse_lastz_txt(self_txt)
        inv_self = extract_inversions(df_self, label=f"{ref_hap} self")

        for _, r in inv_self.iterrows():
            inv_by_hap[ref_hap].append((r.ref_start, r.ref_end))

        # Pairwise alignments
        for _, hap_row in gdf.iterrows():
            hap = hap_row["Haplotype"]
            if hap == ref_hap:
                continue

            txt, reversed_pair = ensure_pairwise(
                ref_row, hap_row, input_dir, aln_dir
            )

            df_pw = parse_lastz_txt(txt, reversed_pair=reversed_pair)
            inv_pw = extract_inversions(df_pw, label=f"{ref_hap} vs {hap}")

            for _, r in inv_pw.iterrows():
                inv_by_hap[hap].append((r.ref_start, r.ref_end))

        ref_inversions = inv_by_hap.get(ref_hap, [])
        if not ref_inversions:
            print(f"  no reference inversions for {group}")
            continue

        for i, (rs, re) in enumerate(ref_inversions):
            inv_id = f"{group}_inv_{i+1}"
            inv_len = re - rs

            for hap in gdf["Haplotype"]:
                present = 0
                for s, e in inv_by_hap.get(hap, []):
                    overlap = max(0, min(e, re) - max(s, rs))
                    if overlap / inv_len >= OVERLAP_FRAC:
                        present = 1
                        break

                rows.append({
                    group_col: group,
                    "Inversion": inv_id,
                    "RefStart": rs,
                    "RefEnd": re,
                    "InvLen": inv_len,
                    "Haplotype": hap,
                    "Present": present
                })

    return pd.DataFrame(rows)

# -----------------------------
# Main
# -----------------------------
def main(args):
    df = pd.read_csv(args.summary, sep="\t")
    os.makedirs(args.lastz_dir, exist_ok=True)
    os.makedirs(args.output_dir, exist_ok=True)

    print("\n=== Processing SubGroups ===")
    sub_df = process_group(
        df,
        group_col="SubGroup",
        ref_override=SUBGROUP_REF_OVERRIDE,
        input_dir=args.input_dir,
        aln_dir=args.lastz_dir
    )
    sub_out = os.path.join(args.output_dir, "housefinch_inversions_by_SubGroup.tsv")
    sub_df.to_csv(sub_out, sep="\t", index=False)
    print(f"Saved {sub_out}")

    print("\n=== Processing States ===")
    state_df = process_group(
        df,
        group_col="State",
        ref_override=STATE_REF_OVERRIDE,
        input_dir=args.input_dir,
        aln_dir=args.lastz_dir
    )
    state_out = os.path.join(args.output_dir, "housefinch_inversions_by_State.tsv")
    state_df.to_csv(state_out, sep="\t", index=False)
    print(f"Saved {state_out}")

# -----------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--summary", required=True)
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Base directory containing Order/Species/Haplotype/")
    parser.add_argument("-l", "--lastz_dir", required=True,
                        help="Directory with pairwise LASTZ alignments")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="Directory for output TSVs")
    args = parser.parse_args()

    main(args)
