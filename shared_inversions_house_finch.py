#!/usr/bin/env python3
import os
import argparse
import pandas as pd
from collections import defaultdict

def categorize_length(length):
    """Categorize inversion length into bins."""
    if 1000 <= length < 2500:
        return "1000-2500"
    elif 2500 <= length < 5000:
        return "2500-5000"
    elif 5000 <= length < 7500:
        return "5000-7500"
    elif 7500 <= length <= 10000:
        return "7500-10000"
    else:
        return None  # out of range

def load_inversions(file_path):
    """
    Load inversion data.
    Expects columns: start1, end1, start2+, end2+, id%
    and identifies inversions as regions where strand2 is '-'
    """
    try:
        df = pd.read_csv(file_path, sep='\t')
    except Exception as e:
        print(f"Could not read {file_path}: {e}")
        return pd.DataFrame()

    # try to find strand columns
    strand_cols = [c for c in df.columns if 'strand' in c.lower()]
    if len(strand_cols) < 2:
        print(f"No strand columns found in {file_path}")
        return pd.DataFrame()

    # Inversions typically have opposite strand orientation
    inv_df = df[df[strand_cols[1]] == '-'].copy()

    if 'length1' in df.columns:
        inv_df['length'] = df['length1']
    elif 'length' in df.columns:
        inv_df['length'] = df['length']
    else:
        inv_df['length'] = abs(df['end1'] - df['start1']) if {'end1','start1'}.issubset(df.columns) else 0

    inv_df['length_category'] = inv_df['length'].apply(categorize_length)
    inv_df = inv_df.dropna(subset=['length_category'])
    return inv_df

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", required=True, help="Top-level input directory")
    parser.add_argument("-s", "--summary", required=True, help="Path to IGH_filtered_table.tsv")

    args = parser.parse_args()

    # Load summary
    summary = pd.read_csv(args.summary, sep="\t")
    hf_haplotypes = summary.loc[summary['Species'] == "House_Finch", 'Haplotype'].unique()

    songbird_dir = os.path.join(args.input_dir, "Songbirds")
    pairwise_dir = os.path.join(songbird_dir, "pairwise_alignments")
    if not os.path.exists(pairwise_dir):
        print(f"No pairwise_alignment folder found in {songbird_dir}")
        return

    # Prepare results
    results = []

    for _, row in summary.iterrows():
        species = row['Species']
        haplotype = row['Haplotype']

        if species == "House_Finch":
            continue

        for hf_hap in hf_haplotypes:
            # Check both possible file naming
            file1 = os.path.join(pairwise_dir, f"{hf_hap}_{haplotype}.txt")
            file2 = os.path.join(pairwise_dir, f"{haplotype}_{hf_hap}.txt")

            file_path = file1 if os.path.exists(file1) else (file2 if os.path.exists(file2) else None)
            if file_path is None:
                print(f"No pairwise file for {hf_hap} and {haplotype}")
                continue
            inv_df = load_inversions(file_path)
            if inv_df.empty:
                continue

            # Count by length category
            counts = inv_df['length_category'].value_counts().to_dict()

            for cat, num in counts.items():
                results.append({
                    "Species1": "House_Finch",
                    "Species2": species,
                    "Haplotype1": hf_hap,
                    "Haplotype2": haplotype,
                    "Inversion_Length_Category": cat,
                    "num_inversions": num
                })

    # Save output
    output_file = os.path.join(songbird_dir, "House_Finch_shared_inversions.tsv")
    pd.DataFrame(results).to_csv(output_file, sep="\t", index=False)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()
