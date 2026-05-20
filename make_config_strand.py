#!/usr/bin/env python3
"""
Create config_strand_IGH.csv and config_strand_IGL.csv for a given Order,
compatible with patchworkplot.

Columns: SampleID, Label, Fasta, Annotation, Strand
SampleID and Label: {Haplotype}

When multiple contigs exist for the same haplotype×locus, only the one with
the highest NumV is used.
"""

import argparse
import os
import pandas as pd


def build_rows(df, locus, input_dir):
    sub = (df[df["Locus"] == locus]
           .sort_values("NumV", ascending=False)
           .groupby("Haplotype", as_index=False)
           .first())

    rows = []
    for _, row in sub.iterrows():
        order     = row["Order"]
        species   = row["Species"]
        haplotype = row["Haplotype"]
        contig    = row["Contig"]
        numv      = row["NumV"]

        fasta = os.path.join(
            input_dir, order, species, haplotype,
            "refined_ig_loci", "igloci_fasta",
            f"{locus}_{contig}_{numv}Vs.fasta"
        )
        annotation = os.path.join(
            input_dir, order, species, haplotype,
            f"{contig}_{locus}_strand.bed"
        )

        rows.append({
            "SampleID":   haplotype,
            "Label":      haplotype,
            "Fasta":      fasta,
            "Annotation": annotation,
            "Strand":     "",
        })

    return pd.DataFrame(rows, columns=["SampleID", "Label", "Fasta", "Annotation", "Strand"])


def main():
    parser = argparse.ArgumentParser(
        description="Generate config_strand_IGH.csv and config_strand_IGL.csv for patchworkplot"
    )
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Top-level input directory (e.g. /local/storage/kav67/clean_birds)")
    parser.add_argument("-s", "--summary", required=True,
                        help="Path to summary_features.csv")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="Directory to write the config CSV files")
    parser.add_argument("--order", required=True,
                        help="Taxonomic order to include (e.g. Doves)")

    args = parser.parse_args()

    df = pd.read_csv(args.summary)
    df = df[df["Order"] == args.order]

    if df.empty:
        print(f"No rows found for Order='{args.order}'")
        return

    os.makedirs(args.output_dir, exist_ok=True)

    for locus in ("IGH", "IGL"):
        out_df = build_rows(df, locus, args.input_dir)
        if out_df.empty:
            print(f"No {locus} rows found for Order='{args.order}', skipping.")
            continue
        out_path = os.path.join(args.output_dir, f"config_strand_{locus}.csv")
        out_df.to_csv(out_path, index=False)
        print(f"[{locus}] Written {len(out_df)} rows to {out_path}")


if __name__ == "__main__":
    main()
