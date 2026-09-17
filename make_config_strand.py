#!/usr/bin/env python3
"""
Create config_strand_IGH.csv and config_strand_IGL.csv for a given Order (or
an explicit list of Species), compatible with patchworkplot.

Columns: SampleID, Label, Fasta, Annotation, Strand
SampleID: {Haplotype}. Label: {Species} (English common name, underscores as
spaces) on its own line(s), plus a second line with the haplotype suffix
(pri/alt/HiC) when the haplotype name ends in one -- the Latin-epithet or
contig-accession part of the haplotype name (e.g. "Scandens", "Magnirostris",
"GCA_902806625_1") is dropped.

When multiple contigs exist for the same haplotype×locus, only the one with
the highest NumV is used.
"""

import argparse
import os
import re
import textwrap
import pandas as pd


def build_label(species, haplotype, width=13):
    species_lines = textwrap.wrap(species.replace("_", " "), width=width) or [species]
    m = re.search(r"_(pri|alt|hic)$", haplotype, re.IGNORECASE)
    if m:
        suffix = m.group(1).lower()
        species_lines.append("HiC" if suffix == "hic" else suffix)
    return "\n".join(species_lines)


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
            "Label":      build_label(species, haplotype),
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
    parser.add_argument("--order",
                        help="Taxonomic order to include (e.g. Doves)")
    parser.add_argument("--species", nargs="+",
                        help="One or more Species names to include instead of --order "
                             "(e.g. --species Large_ground_finch Small_tree_finch)")
    parser.add_argument("--locus", choices=["IGH", "IGL", "both"], default="both",
                        help="Which locus config(s) to write (default: both)")

    args = parser.parse_args()

    if not args.order and not args.species:
        parser.error("one of --order or --species is required")

    df = pd.read_csv(args.summary)
    if args.species:
        df = df[df["Species"].isin(args.species)]
        missing = set(args.species) - set(df["Species"].unique())
        if missing:
            print(f"Warning: no rows found for species {sorted(missing)}")
    else:
        df = df[df["Order"] == args.order]

    if df.empty:
        print(f"No rows found for the given filter")
        return

    os.makedirs(args.output_dir, exist_ok=True)

    loci = ("IGH", "IGL") if args.locus == "both" else (args.locus,)
    for locus in loci:
        out_df = build_rows(df, locus, args.input_dir)
        if out_df.empty:
            print(f"No {locus} rows found for Order='{args.order}', skipping.")
            continue
        out_path = os.path.join(args.output_dir, f"config_strand_{locus}.csv")
        out_df.to_csv(out_path, index=False)
        print(f"[{locus}] Written {len(out_df)} rows to {out_path}")


if __name__ == "__main__":
    main()
