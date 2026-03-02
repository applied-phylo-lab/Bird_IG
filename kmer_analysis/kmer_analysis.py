import sys
import os
from collections import defaultdict, Counter


# ---------------------------
# Reverse complement
# ---------------------------
def revcomp(seq):
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


# ---------------------------
# Load FASTA
# ---------------------------
def load_fasta(fasta_file):
    genome = {}
    with open(fasta_file) as f:
        name = None
        seq = []
        for line in f:
            if line.startswith(">"):
                if name:
                    genome[name] = "".join(seq).upper()
                name = line.strip().split()[0][1:]
                seq = []
            else:
                seq.append(line.strip())
        if name:
            genome[name] = "".join(seq).upper()
    return genome


# ---------------------------
# Count kmers (canonicalized)
# ---------------------------
def count_kmers(genome, k):
    counts = Counter()

    for chrom in genome:
        seq = genome[chrom]
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            if "N" in kmer:
                continue

            rc = revcomp(kmer)
            canonical = min(kmer, rc)
            counts[canonical] += 1

    return counts


# ---------------------------
# Find positions for selected kmers
# ---------------------------
def find_positions(genome, k, selected_kmers):

    positions = []

    for chrom in genome:
        seq = genome[chrom]
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            if "N" in kmer:
                continue

            rc = revcomp(kmer)
            canonical = min(kmer, rc)

            if canonical in selected_kmers:
                strand = "+"
                if kmer != canonical:
                    strand = "-"

                positions.append(
                    (canonical, chrom, i, i+k, strand)
                )

    return positions


# ---------------------------
# Main
# ---------------------------
def main(fasta_file, k, output_dir, top_n):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    genome = load_fasta(fasta_file)

    print("Counting kmers...")
    counts = count_kmers(genome, k)

    print("Selecting top kmers...")
    top_kmers = counts.most_common(top_n)

    abundance_path = os.path.join(output_dir, "kmer_abundance.tsv")
    positions_path = os.path.join(output_dir, "kmer_positions.tsv")

    with open(abundance_path, "w") as out:
        out.write("kmer\tcount\n")
        for kmer, count in top_kmers:
            out.write(f"{kmer}\t{count}\n")

    selected_set = set(kmer for kmer, _ in top_kmers)

    print("Locating kmer positions...")
    positions = find_positions(genome, k, selected_set)

    with open(positions_path, "w") as out:
        out.write("kmer\tchrom\tstart\tend\tstrand\n")
        for kmer, chrom, start, end, strand in positions:
            out.write(f"{kmer}\t{chrom}\t{start}\t{end}\t{strand}\n")

    print("Done.")
    print("Outputs:")
    print(abundance_path)
    print(positions_path)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py genome.fa k output_dir top_n")
        sys.exit(1)

    main(
        sys.argv[1],
        int(sys.argv[2]),
        sys.argv[3],
        int(sys.argv[4])
    )
