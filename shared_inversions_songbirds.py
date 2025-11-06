#!/usr/bin/env python3
import os
import argparse
import pandas as pd

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
    
def categorize_identity(identity):
    """Categorize identity percentage into bins."""
    if identity < 60:
        return "<60%"
    elif 60 <= identity < 70:
        return "60-70%"
    elif 70 <= identity < 80:
        return "70-80%"
    elif 80 <= identity < 90:
        return "80-90%"
    elif 90 <= identity <= 100:
        return "90-100%"
    else:
        return None

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

    id_cols = [c for c in df.columns if 'id' in c.lower() or 'identity' in c.lower()]
    if id_cols:
        id_col = id_cols[0]
        # Convert to numeric (handle strings like '95.4%')
        inv_df[id_col] = inv_df[id_col].astype(str).str.replace('%', '', regex=False)
        inv_df[id_col] = pd.to_numeric(inv_df[id_col], errors='coerce')
        inv_df['identity_category'] = inv_df[id_col].apply(categorize_identity)
    else:
        inv_df['identity_category'] = None


    inv_df['length_category'] = inv_df['length'].apply(categorize_length)
    inv_df = inv_df.dropna(subset=['length_category', 'identity_category'])
    return inv_df

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", required=True, help="Top-level input directory")
    parser.add_argument("-s", "--summary", required=True, help="Path to IGH_filtered_table.tsv")
    args = parser.parse_args()

    # Load summary
    summary = pd.read_csv(args.summary, sep="\t")

    # Filter for Songbirds and *_pri haplotypes
    songbirds_summary = summary[
        (summary['Order'].str.contains("Songbirds", case=False)) &
        (summary['Haplotype'].str.endswith("_pri"))
    ]

    songbird_dir = os.path.join(args.input_dir, "Songbirds")
    pairwise_dir = os.path.join(songbird_dir, "pairwise_alignments")
    if not os.path.exists(pairwise_dir):
        print(f"No pairwise_alignment folder found in {songbird_dir}")
        return

    # Prepare results
    results = []

    haplotypes = songbirds_summary[['Species', 'Haplotype']].values.tolist()
    print(f"Processing {len(haplotypes)} haplotypes from {len(songbirds_summary['Species'].unique())} species.")
    # Iterate over all pairwise combinations
    for i in range(len(haplotypes)):
        species1, hap1 = haplotypes[i]
        for j in range(i+1, len(haplotypes)):
            species2, hap2 = haplotypes[j]

            # Check both possible file naming
            file1 = os.path.join(pairwise_dir, f"{hap1}_{hap2}.txt")
            file2 = os.path.join(pairwise_dir, f"{hap2}_{hap1}.txt")
            file_path = file1 if os.path.exists(file1) else (file2 if os.path.exists(file2) else None)
            if file_path is None:
                print(f"No pairwise file for {hap1} and {hap2}")
                continue

            inv_df = load_inversions(file_path)
            if inv_df.empty:
                continue

            # Count by length category
            counts = inv_df.groupby(['length_category', 'identity_category']).size().reset_index(name='num_inversions')
            for _, row in counts.iterrows():
                results.append({
                    "Species1": species1,
                    "Species2": species2,
                    "Haplotype1": hap1,
                    "Haplotype2": hap2,
                    "Inversion_Length_Category": row['length_category'],
                    "Identity_Category": row['identity_category'],
                    "num_inversions": row['num_inversions']
                })

    # Save output
    output_file = os.path.join(songbird_dir, "Songbirds_shared_inversions.tsv")
    pd.DataFrame(results).to_csv(output_file, sep="\t", index=False)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()
