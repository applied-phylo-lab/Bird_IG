import os
import argparse
import pandas as pd
from multiprocessing import Pool
from Bio import SeqIO
from Bio.Seq import Seq

HEPT_LEN = 7
SPACER_LEN = 23
NONAMER_LEN = 9
RSS_LEN = HEPT_LEN + SPACER_LEN + NONAMER_LEN


def extract_rss_block(seq, start, end, strand):
    """Extract the exact 7+23+9 RSS block immediately behind the V gene."""

    if strand == "+":
        rss_start = end
        rss_end = end + RSS_LEN
        if rss_end > len(seq):
            return None
        block = seq[rss_start:rss_end]

    else:  # strand == "-"
        rss_start = start - RSS_LEN
        rss_end = start
        if rss_start < 0:
            return None
        block = seq[rss_start:rss_end].reverse_complement()

    hept = str(block[:HEPT_LEN]).upper()
    spacer = str(block[HEPT_LEN:HEPT_LEN + SPACER_LEN]).upper()
    nona = str(block[HEPT_LEN + SPACER_LEN:RSS_LEN]).upper()

    return {
        "7-mer index": rss_start,                      # <-- REAL genomic coordinate
        "9-mer index": rss_start + HEPT_LEN + SPACER_LEN,
        "7-mer": hept,
        "spacer": spacer,
        "9-mer": nona
    }



def process_row(row):
    """Process one haplotype and return a list of RSS hits (not writing files)."""
    input_dir = row['InputDir']
    order = row['Order']
    species = row['Species']
    haplotype = row['Haplotype']
    contig = row['Contig']
    numv = row['NumV']

    fasta_path = os.path.join(
        input_dir, order, species, haplotype,
        "refined_ig_loci", "igloci_fasta",
        f"IGH_{contig}_{numv}Vs.fasta"
    )
    bed_path = os.path.join(input_dir, order, species, haplotype, "IGH.bed")

    if not (os.path.exists(fasta_path) and os.path.exists(bed_path)):
        print(f"Missing files for {order}/{species}/{haplotype}")
        return []

    fasta_records = {rec.id: rec for rec in SeqIO.parse(fasta_path, "fasta")}
    print(fasta_records.keys())
    bed = pd.read_csv(
    bed_path, sep="\t", header=None,
        names=[
            "contig",   # 0
            "start",    # 1
            "end",      # 2
            "NA1",      # 3
            "NA2",      # 4
            "strand",   # 5
            "NA3",      # 6
            "NA4",      # 7
            "color"     # 8
        ]
    )


    collected = []

    # Only one contig in the FASTA
    rec = next(iter(fasta_records.values()))
    seq = rec.seq
    ctg = rec.id   # add this so contig column is correct

    for _, b in bed.iterrows():
        start = int(b["start"])
        end = int(b["end"])
        strand = b["strand"]

        rss = extract_rss_block(seq, start, end, strand)
        if rss is None:
            continue

        rss["contig"] = ctg
        rss["strand"] = strand
        rss["Order"] = order
        rss["Species"] = species
        rss["Haplotype"] = haplotype

        collected.append(rss)

    print(f"Processed RSS for {order}/{species}/{haplotype}")

    return collected


def main():
    parser = argparse.ArgumentParser(description="Extract RSS motifs from V genes across all haplotypes.")
    parser.add_argument("-i", "--input_dir", required=True)
    parser.add_argument("-s", "--summary", required=True)
    parser.add_argument("-c", "--cores", type=int, default=4)
    args = parser.parse_args()

    df = pd.read_csv(args.summary, sep="\t")
    df["InputDir"] = args.input_dir

    with Pool(args.cores) as pool:
        all_results = pool.map(process_row, [row for _, row in df.iterrows()])

    # Flatten list of lists
    flat = [item for sublist in all_results for item in sublist]

    # Save combined file
    out_path = os.path.join(args.input_dir, "all_RSS.tsv")

    if flat:
        pd.DataFrame(flat).to_csv(out_path, sep="\t", index=False)
    else:
        pd.DataFrame(columns=[
            "Order", "Species", "Haplotype",
            "contig", "strand",
            "7-mer index", "9-mer index",
            "7-mer", "spacer", "9-mer"
        ]).to_csv(out_path, sep="\t", index=False)

    print(f"Generated combined RSS file: {out_path}")


if __name__ == "__main__":
    main()
