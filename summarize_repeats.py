import sys
import os
from collections import defaultdict


# ---------------------------
# Load genome
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
# Parse ALL alignments
# ---------------------------
def parse_lastz_all(file, min_length):
    intervals = defaultdict(list)

    with open(file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split()
            if len(parts) < 9:
                continue

            chrom1 = parts[0]
            start1 = int(parts[2])
            end1   = int(parts[3])
            chrom2 = parts[5]
            start2 = int(parts[7])
            end2   = int(parts[8])

            length = end1 - start1
            if length < min_length:
                continue

            # add BOTH sides
            intervals[chrom1].append((start1, end1))
            intervals[chrom2].append((start2, end2))

    return intervals


# ---------------------------
# Iterative overlap refinement
# ---------------------------
def refine_intervals(intervals):

    intervals = sorted(intervals)

    changed = True
    while changed:
        changed = False
        new_intervals = []

        for i in range(len(intervals)):
            s1, e1 = intervals[i]

            for j in range(i + 1, len(intervals)):
                s2, e2 = intervals[j]

                if s1 < e2 and s2 < e1:
                    # intersect
                    ns = max(s1, s2)
                    ne = min(e1, e2)

                    if ne - ns > 0 and (ns, ne) != (s1, e1):
                        s1, e1 = ns, ne
                        changed = True

            new_intervals.append((s1, e1))

        intervals = sorted(set(new_intervals))

    return intervals


# ---------------------------
# Main
# ---------------------------
def main(lastz_file, fasta_file, output_dir, min_length):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    genome = load_fasta(fasta_file)
    raw_intervals = parse_lastz_all(lastz_file, min_length)

    seq_path = os.path.join(output_dir, "repeat_sequences.tsv")
    pos_path = os.path.join(output_dir, "repeat_positions.tsv")

    seq_out = open(seq_path, "w")
    pos_out = open(pos_path, "w")

    seq_out.write("repeat_id\tsequence\n")
    pos_out.write("repeat_id\tchrom\tstart\tend\tlength\n")

    repeat_sequences = {}
    repeat_id_counter = 1

    for chrom in raw_intervals:

        refined = refine_intervals(raw_intervals[chrom])

        for start, end in refined:

            length = end - start
            if length < min_length:
                continue

            seq = genome[chrom][start:end].upper()

            # Deduplicate identical sequences
            if seq in repeat_sequences:
                rep_id = repeat_sequences[seq]
            else:
                rep_id = f"rep_{repeat_id_counter}"
                repeat_sequences[seq] = rep_id
                seq_out.write(f"{rep_id}\t{seq}\n")
                repeat_id_counter += 1

            pos_out.write(f"{rep_id}\t{chrom}\t{start}\t{end}\t{length}\n")

    seq_out.close()
    pos_out.close()


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py lastz.tsv genome.fa output_dir min_length")
        sys.exit(1)

    main(
        sys.argv[1],
        sys.argv[2],
        sys.argv[3],
        int(sys.argv[4])
    )
