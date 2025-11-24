import argparse
import os
import pandas as pd
from Bio import SeqIO

WINDOW = 100
BED_COLUMNS = ["contig", "start", "end", "NA1", "NA2", "strand", "NA3", "NA4", "color"]

def extract_window(seq, start, end, strand, window=WINDOW):
    seq_len = len(seq)

    if strand == "+":
        s = end
        e = min(end + window, seq_len)
        block = seq[s:e]
    else:
        s = max(0, start - window)
        e = start
        block = seq[s:e].reverse_complement()

    block = str(block).upper()
    if len(block) < window:
        block += "N" * (window - len(block))
    return block


def process_row(row):
    input_dir = row['InputDir']
    order = row['Order']
    species = row['Species']
    haplotype = row['Haplotype']
    contig = row['Contig']  # Only one contig in each fasta/bed
    numv = row['NumV']

    # Paths
    fasta_path = os.path.join(
        input_dir, order, species, haplotype,
        'refined_ig_loci', 'igloci_fasta',
        f'IGH_{contig}_{numv}Vs.fasta'
    )
    bed_path = os.path.join(input_dir, order, species, haplotype, "IGH.bed")
    if not (os.path.exists(fasta_path) and os.path.exists(bed_path)):
        print(f"Missing files for {order}/{species}/{haplotype}")
        return []
    
    # Load FASTA contig
    record = next(SeqIO.parse(fasta_path, "fasta"))
    seq = record.seq

    # Read BED
    windows = []
    with open(bed_path) as f:
        for line in f:
            if line.strip() == "" or line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            fields += [""] * (9 - len(fields))  # ensure 9 columns
            bed = dict(zip(BED_COLUMNS, fields))

            start = int(bed["start"])
            end = int(bed["end"])
            strand = bed["strand"]
            contig_name = bed["contig"]

            block = extract_window(seq, start, end, strand)

            seq_id = (
                f"{order}|{species}|{haplotype}|{contig_name}|"
                f"{start}-{end}|{strand}"
            )

            windows.append((seq_id, block))

    return windows


def main():
    parser = argparse.ArgumentParser(description="Extract 100bp downstream of V genes for all haplotypes")
    parser.add_argument("-i", "--input_dir", required=True)
    parser.add_argument("-s", "--summary", required=True)
    parser.add_argument("-o", "--output", required=True, help="Combined output FASTA file")
    parser.add_argument("-c", "--cores", type=int, default=4)
    args = parser.parse_args()

    input_dir = args.input_dir
    summary_file = args.summary
    n_cores = args.cores
    output_file = args.output

    df = pd.read_csv(summary_file, sep="\t")
    df["InputDir"] = input_dir

    from multiprocessing import Pool
    all_windows = []

    with Pool(n_cores) as pool:
        results = pool.map(process_row, [row for _, row in df.iterrows()])

    for res in results:
        all_windows.extend(res)

    with open(output_file, "w") as out:
        for seq_id, seq in all_windows:
            out.write(f">{seq_id}\n{seq}\n")


if __name__ == "__main__":
    main()
