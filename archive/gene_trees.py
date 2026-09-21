#!/usr/bin/env python3
import argparse
import os
import pandas as pd
from multiprocessing import Pool
from collections import defaultdict
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def load_summary(path):
    df = pd.read_csv(path, sep="\t")
    return df


def process_one_row(args):
    """
    For one summary row:
      - read combined_genes_IGH_clean.txt
      - find main contig
      - return filtered dataframe
    """
    row, input_dir = args

    order = row["Order"]
    species = row["Species"]
    hap = row["Haplotype"]

    genes_path = os.path.join(
        input_dir, order, species, hap, "combined_genes_IGH_clean.txt"
    )

    if not os.path.exists(genes_path):
        print(f"[WARN] Missing: {genes_path}")
        return None

    df = pd.read_csv(genes_path, sep="\t")

    # find main contig = contig with most entries
    main_contig = df["Contig"].value_counts().idxmax()

    filtered = df[df["Contig"] == main_contig].copy()
    filtered["Order"] = order
    filtered["Species"] = species
    filtered["Haplotype"] = hap

    return filtered


def write_order_fastas(order, rows, outdir):
    """
    rows: a list of filtered dataframes belonging to this order
    Writes:
      - genes.fasta
      - genes_annotation.tsv
    Then runs clustalo + iqtree
    """

    if order == "MiscBirds":  # skip
        return

    order_dir = os.path.join(outdir, order)
    os.makedirs(order_dir, exist_ok=True)

    fasta_path = os.path.join(order_dir, "genes.fasta")
    annot_path = os.path.join(order_dir, "genes_annotation.tsv")

    # combine all dfs
    df = pd.concat(rows, ignore_index=True)

    records = []
    annot_rows = []

    for _, r in df.iterrows():
        seq = r["Sequence"]
        strand = r["Strand"]

        # strand symbol
        strand_symbol = "_" if strand == "+" else "-"

        header = f"{r['Species']}_{r['Haplotype']}_{r['Contig']}_{r['Pos']}_{strand_symbol}_{r['Productive']}"

        records.append(
            SeqRecord(
                Seq(seq),
                id=header,
                description=""
            )
        )

        annot_rows.append({
            "taxa": header,
            "Species": r["Species"],
            "Haplotype": r["Haplotype"],
            "Strand": r["Strand"],
            "Productive": r["Productive"]
        })

    # write FASTA
    with open(fasta_path, "w") as f:
        SeqIO.write(records, f, "fasta")

    # write annotation
    pd.DataFrame(annot_rows).to_csv(annot_path, sep="\t", index=False)

    # run clustalo
    aligned = os.path.join(order_dir, "genes_aligned.fasta")
    cmd1 = [
        "clustalo",
        "-i", fasta_path,
        "-t", "DNA",
        "--outfmt=fasta",
        "-o", aligned
    ]
    subprocess.run(cmd1, check=True)

    # run iqtree
    cmd2 = [
        "iqtree",
        "-s", aligned,
        "--prefix", os.path.join(order_dir, "genes_aligned")
    ]
    subprocess.run(cmd2, check=True)


def main():
    parser = argparse.ArgumentParser(description="Extract main-contig genes and build order-level trees")
    parser.add_argument("-i", "--input_dir", required=True)
    parser.add_argument("-s", "--summary", required=True)
    parser.add_argument("-c", "--cores", type=int, default=4)
    parser.add_argument("-o", "--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    summary = load_summary(args.summary)

    # parallel extraction
    work = [(row, args.input_dir) for _, row in summary.iterrows()]

    with Pool(args.cores) as p:
        results = p.map(process_one_row, work)

    # collect valid dfs
    results = [r for r in results if r is not None]

    # group by Order
    order_dict = defaultdict(list)
    for df in results:
        order = df["Order"].iloc[0]
        order_dict[order].append(df)

    # parallel build per order
    order_jobs = [(order, dfs, args.output_dir) for order, dfs in order_dict.items()]

    with Pool(args.cores) as p:
        p.starmap(write_order_fastas, order_jobs)


if __name__ == "__main__":
    main()
