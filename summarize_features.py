import os
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description="Summarize main feature stats.")
    parser.add_argument("-i", "--indir",required=True, help="Input directory containing self*.tsv files")
    #parser.add_argument("-t", "--table", required=True, help="tsv file with VGP summary table")
    parser.add_argument("-o", "--out",default="summary.tsv", help="Output TSV file")
    args = parser.parse_args()

    cols_to_keep = ["Locus", "Contig", "Length", "NumV", "NumProdV", "FracProdV"]

    # Store results
    summary_rows = []

    # Walk through all subdirectories
    for root, dirs, files in os.walk(args.indir):
        if "summary.csv" in files and "refined_ig_loci" in root:
            summary_path = os.path.join(root, "summary.csv")

            # Split path to extract order, species, haplotype
            parts = root.split(os.sep)
            try:
                order = parts[-4]
                species = parts[-3]
                haplotype = parts[-2]
            except IndexError:
                print(f"Skipping {summary_path} (unexpected folder structure)")
                continue

            # Read CSV
            try:
                df = pd.read_csv(summary_path, sep=",", dtype=str)
            except Exception as e:
                print(f"Could not read {summary_path}: {e}")
                continue

            # Keep only Locus == IGH or IGL
            df = df[df["Locus"].isin(["IGH", "IGL"])]

            # Select relevant columns
            df = df[cols_to_keep]

            # Add metadata columns
            df.insert(0, "Order", order)
            df.insert(1, "Species", species)
            df.insert(2, "Haplotype", haplotype)

            summary_rows.append(df)

    # Combine all
    if summary_rows:
        combined = pd.concat(summary_rows, ignore_index=True)
        combined.to_csv(args.out, index=False)
        print(f"Combined summary written to: {args.out}")
    else:
        print("No summary.csv files found.")

if __name__ == "__main__":
    main()
