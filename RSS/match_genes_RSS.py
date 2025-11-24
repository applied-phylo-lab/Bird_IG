#!/usr/bin/env python3

import argparse
import os
import pandas as pd
from multiprocessing import Pool


def process_row(task):
    """
    Process one summary row:
    - Load combined_genes_IGH_clean.txt
    - Choose majority contig
    - Load RSS file and filter to that contig
    - Check 7-mer index within 100bp downstream from each gene
    - Return summary counts + gene-RSS pair table
    """
    (
        input_dir,
        order,
        species,
        haplotype
    ) = task

    # Paths
    genes_path = os.path.join(
        input_dir, order, species, haplotype,
        "combined_genes_IGH_clean.txt"
    )

    rss_path = os.path.join(
        input_dir, order, species, haplotype,
        "RSS", "rss_V.csv"
    )

    if not os.path.exists(genes_path):
        print(f"[WARN] Missing genes file: {genes_path}")
        return None

    if not os.path.exists(rss_path):
        print(f"[WARN] Missing RSS file: {rss_path}")
        return None

    # Load data
    genes = pd.read_csv(genes_path, sep="\t")
    rss = pd.read_csv(rss_path, sep="\t")

    # -----------------------------------------------------
    # 1. Select the contig with the most rows
    # -----------------------------------------------------
    contig_counts = genes["Contig"].value_counts()
    if contig_counts.empty:
        return None

    major_contig = contig_counts.idxmax()

    # Filter genes
    genes_major = genes[genes["Contig"] == major_contig]

    # Filter RSS
    rss_major = rss[rss["reference contig"] == major_contig]

    # -----------------------------------------------------
    # 2. Check RSS within 100 bp downstream of gene Pos
    # -----------------------------------------------------
    results = []

    # Build intervals quickly
    gene_positions = genes_major[["GeneType", "Contig", "Pos", "Strand"]]

    for _, rss_row in rss_major.iterrows():
        rss_pos = rss_row["7-mer index"]

        # Find genes where Pos < rss_pos < Pos + 100

        

        for _, g in genes_major.iterrows():
            gene_pos = g["Pos"] + len(g["Sequence"])
            strand = g["Strand"]

            if strand == "+" and gene_pos < rss_pos <= gene_pos + 100:
                downstream = True
            elif strand == "-" and gene_pos - 100 <= rss_pos < gene_pos:
                downstream = True
            else:
                downstream = False

            if downstream:
                results.append({
                    "Order": order,
                    "Species": species,
                    "Haplotype": haplotype,
                    "Contig": major_contig,
                    "GeneType": g["GeneType"],
                    "GenePos": gene_pos,
                    "GeneStrand": strand,
                    "RSS_7mer_index": rss_pos,
                    "RSS_9mer_index": rss_row["9-mer index"],
                    "RSS_7mer": rss_row["7-mer"],
                    "RSS_9mer": rss_row["9-mer"],
                    "RSS_strand": rss_row["strand"],
                })
    # Count genes retained
    retained_rss_count = len(results)

    return {
        "order": order,
        "species": species,
        "haplotype": haplotype,
        "contig_rss": len(rss_major),
        "retained_rss": retained_rss_count,
        "pairs": results,
    }


def main():
    parser = argparse.ArgumentParser(description="Filter genes and RSS per contig")
    parser.add_argument("-i", "--input_dir", required=True)
    parser.add_argument("-s", "--summary", required=True)
    parser.add_argument("-c", "--cores", type=int, default=4)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.summary, sep="\t")

    tasks = [
        (args.input_dir, row["Order"], row["Species"], row["Haplotype"])
        for _, row in df.iterrows()
    ]

    with Pool(processes=args.cores) as pool:
        outputs = pool.map(process_row, tasks)

    # Remove None tasks
    outputs = [o for o in outputs if o is not None]

    # -----------------------------------------------------
    # 1. Summary: gene counts per species/haplotype
    # -----------------------------------------------------
    summary_rows = [
        {
            "Order": o["order"],
            "Species": o["species"],
            "Haplotype": o["haplotype"],
            "ContigRSS": o["contig_rss"],
            "RetainedRSS": o["retained_rss"],
        }
        for o in outputs
    ]

    summary_df = pd.DataFrame(summary_rows)
    summary_file = os.path.join(args.output, "RSS_counts.tsv")
    summary_df.to_csv(summary_file, sep="\t", index=False)

    # -----------------------------------------------------
    # 2. Combined gene–RSS pairs
    # -----------------------------------------------------
    all_pairs = []
    for o in outputs:
        all_pairs.extend(o["pairs"])

    pairs_df = pd.DataFrame(all_pairs)
    pairs_file = os.path.join(args.output, "gene_rss_pairs.tsv")
    pairs_df.to_csv(pairs_file, sep="\t", index=False)

    print("Done.")
    print(f"RSS count summary saved to:   {summary_file}")
    print(f"Gene-RSS pairs saved to:       {pairs_file}")


if __name__ == "__main__":
    main()
