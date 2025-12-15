#!/usr/bin/env python3

import os
import argparse
import pandas as pd
from multiprocessing import Pool
from Bio.Seq import Seq
import random

# ================================================================
# Read TSV with optional # header
# ================================================================
def read_tsv_with_header(path):
    """Read TSV even if header starts with #."""
    with open(path) as f:
        header = f.readline().strip()
        if header.startswith("#"):
            header = header[1:]
        columns = header.split("\t")

    try:   
        df = pd.read_csv(path, sep="\t", comment="#", names=columns, skiprows=1)
    except:
        print(f"Error reading {path}")
    df.columns = [c.replace("#", "").replace("%", "").replace("+", "")
                  for c in df.columns]
    return df


# ================================================================
# Reverse complement
# ================================================================
def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())


# ================================================================
# Identity between two aligned sequences
# Gaps are counted as aligned but not as matches
# ================================================================
def alignment_identity(seq1, seq2):
    """
    Compute the identity between two aligned sequences.
    Only columns where both seq1 and seq2 have valid A/C/G/T are considered.
    Returns percentage of matches among valid columns.
    """
    matches = 0
    valid = 0

    valid_bases = set("ACGT")

    for a, b in zip(seq1, seq2):
        if a not in valid_bases or b not in valid_bases:
            continue  # skip gaps or non-ACGT
        valid += 1
        if a.upper() == b.upper():
            matches += 1

    if valid == 0:
        return 0.0

    return 100 * matches / valid


def lastz_percent_identity(seq1, seq2):
    """
    Reproduce LASTZ's percent_identical() behavior exactly.

    Only positions where both seq1 and seq2 have valid A/C/G/T are counted.
    Identity = (200*numMatches + denom) // (2*denom)
    """

    # LASTZ's nuc_to_bits encoding:
    nuc_to_bits = {
        'A': 0, 'C': 1, 'G': 2, 'T': 3,
        'a': 0, 'c': 1, 'g': 2, 't': 3,
    }

    numMatches = 0.0
    denom = 0.0

    # iterate over aligned sequences
    for c1, c2 in zip(seq1, seq2):

        b1 = nuc_to_bits.get(c1, -1)
        b2 = nuc_to_bits.get(c2, -1)

        # must be >= 0 → valid DNA base
        if b1 >= 0 and b2 >= 0:
            denom += 1
            if b1 == b2:
                numMatches += 1

    if denom == 0.0:
        return 0.0

    # LASTZ formula: integer rounded identity
    return round(100.0 * numMatches / denom, 1)



# ================================================================
# Extract middle window from aligned sequences
# ================================================================
def extract_middle(seq1, seq2, bp=100):
    L = len(seq1)
    mid = L // 2
    half = bp // 2
    start = max(0, mid - half)
    end = min(L, mid + half)
    return seq1[start:end], seq2[start:end]

def extract_random_window(seq1, seq2, bp=50):
    L = len(seq1)
    if L <= bp:
        return seq1, seq2
    start = random.randint(0, L - bp)
    end = start + bp
    return seq1[start:end], seq2[start:end]


# ================================================================
# Process one row of summary table
# ================================================================
def process_row(row):
    input_dir = row['InputDir']
    order = row['Order']
    species = row['Species']
    haplotype = row['Haplotype']
    contig = row['Contig']

    # Path to TSV with aligned sequences
    tsv_path = os.path.join(input_dir, order, species, haplotype, "IGH_self.tsv")
    if not os.path.exists(tsv_path):
        return []

    df = read_tsv_with_header(tsv_path)

    # Filter diagonal inversions
    inv = df[(df["strand1"] != df["strand2"]) &
             (df["start1"] == df["start2"]) &
             (df["end1"] == df["end2"])]

    results = []
    for _, row2 in inv.iterrows():
        seq1 = row2["text1"].upper()
        seq2 = row2["text2"].upper()
        #seq2_rc = reverse_complement(seq2)
        id=row2["id"].replace("%", "")

        whole_id = lastz_percent_identity(seq1, seq2)

        mid1, mid2 = extract_middle(seq1, seq2, bp=50)
        mid_id_50 = lastz_percent_identity(mid1, mid2)

        mid1, mid2 = extract_middle(seq1, seq2, bp=20)
        mid_id_20 = lastz_percent_identity(mid1, mid2)

        mid1, mid2 = extract_middle(seq1, seq2, bp=10)
        mid_id_10 = lastz_percent_identity(mid1, mid2)

        mid1, mid2 = extract_middle(seq1, seq2, bp=15)
        mid_id_15 = lastz_percent_identity(mid1, mid2)

        rand1, rand2 = extract_random_window(seq1, seq2, bp=50)
        rand_id_50 = lastz_percent_identity(rand1, rand2)
        rand1, rand2 = extract_random_window(seq1, seq2, bp=10)
        rand_id_10 = lastz_percent_identity(rand1, rand2)
        rand1, rand2 = extract_random_window(seq1, seq2, bp=20)
        rand_id_20 = lastz_percent_identity(rand1, rand2)
        rand1, rand2 = extract_random_window(seq1, seq2, bp=15)
        rand_id_15 = lastz_percent_identity(rand1, rand2)

        results.append({
            "Order": order,
            "Species": species,
            "Haplotype": haplotype,
            "Contig": contig,
            "Start": row2["start1"],
            "End": row2["end1"],
            "ID_Lastz"  : id,
            "WholeIdentity": whole_id,
            "MiddleIdentity10bp": mid_id_10,
            "MiddleIdentity15bp": mid_id_15,
            "MiddleIdentity20bp": mid_id_20,
            "MiddleIdentity50bp": mid_id_50,
            "MiddleMinusWhole": whole_id-mid_id_50,
            "RandomIdentity10bp": rand_id_10,
            "RandomIdentity15bp": rand_id_15,
            "RandomIdentity20bp": rand_id_20,
            "RandomIdentity50bp": rand_id_50,
            "MiddleSeq1": mid1,
            "MiddleSeq2": mid2,
        })

    return results


# ================================================================
# Main script
# ================================================================
def main():
    parser = argparse.ArgumentParser(description="Summarize inversion stats from alignments in IGH_self.tsv")
    parser.add_argument('-i','--input_dir', required=True, help="Top-level input directory")
    parser.add_argument('-s','--summary', required=True, help="Path to summary_table_IGH.tsv")
    parser.add_argument('-o','--output', required=True, help="Output TSV file")
    parser.add_argument('-c','--cores', type=int, default=4, help="Number of parallel workers")

    args = parser.parse_args()

    df = pd.read_csv(args.summary, sep='\t')
    df['InputDir'] = args.input_dir

    # Run in parallel
    with Pool(args.cores) as pool:
        results = pool.map(process_row, [row for _, row in df.iterrows()])

    all_rows = [r for sub in results for r in sub]
    if not all_rows:
        print("No valid data found.")
        return

    out_df = pd.DataFrame(all_rows)
    out_df.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(out_df)} rows to {args.output}")


if __name__ == "__main__":
    main()
