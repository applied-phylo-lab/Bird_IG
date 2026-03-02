import sys
import os
from collections import Counter


def revcomp(seq):
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


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


def find_positions(genome, k, selected):
    records = []

    for chrom in genome:
        seq = genome[chrom]
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            if "N" in kmer:
                continue
            rc = revcomp(kmer)
            canonical = min(kmer, rc)
            if canonical in selected:
                strand = "+" if kmer == canonical else "-"
                records.append((k, canonical, chrom, i, i+k, strand))

    return records


def main(genome_file, k_min, k_max, top_n, outdir):

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    genome = load_fasta(genome_file)

    abundance_file = os.path.join(outdir, "kmer_abundance_all.tsv")
    top_positions_file = os.path.join(outdir, "kmer_positions.tsv")

    abundance_out = open(abundance_file, "w")
    abundance_out.write("k\tkmer\tcount\n")

    pos_out = open(top_positions_file, "w")
    pos_out.write("k\tkmer\tchrom\tstart\tend\tstrand\n")

    for k in range(k_min, k_max + 1):
        print(f"Processing k={k}")
        counts = count_kmers(genome, k)

        for kmer, count in counts.items():
            abundance_out.write(f"{k}\t{kmer}\t{count}\n")

        top_kmers = [x[0] for x in counts.most_common(top_n)]
        positions = find_positions(genome, k, set(top_kmers))

        for rec in positions:
            pos_out.write("\t".join(map(str, rec)) + "\n")

    abundance_out.close()
    pos_out.close()

    print("Done.")
    print("Outputs:")
    print(abundance_file)
    print(top_positions_file)


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python kmer_sweep.py genome.fa k_min k_max top_n output_dir")
        sys.exit(1)

    main(
        sys.argv[1],
        int(sys.argv[2]),
        int(sys.argv[3]),
        int(sys.argv[4]),
        sys.argv[5]
    )
