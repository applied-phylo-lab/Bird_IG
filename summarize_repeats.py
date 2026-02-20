import sys
import os
from collections import defaultdict


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


def parse_lastz(file, min_length):
    inversions = defaultdict(list)

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
            strand = parts[6]
            start2 = int(parts[7])
            end2   = int(parts[8])

            length = end1 - start1
            if length < min_length:
                continue

            if chrom1 == chrom2 and strand == "-" and start1 == start2 and end1 == end2:
                inversions[chrom1].append((start1, end1))

    return inversions


def group_overlapping(intervals):
    intervals = sorted(intervals)
    groups = []

    for start, end in intervals:
        placed = False

        for group in groups:
            for s, e in group:
                if start < e and s < end:
                    group.append((start, end))
                    placed = True
                    break
            if placed:
                break

        if not placed:
            groups.append([(start, end)])

    return groups


def tile_supported(tile_start, tile_end, intervals):
    for s, e in intervals:
        if tile_start >= s and tile_end <= e:
            return True
    return False


def main(lastz_file, fasta_file, output_dir, min_length):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    genome = load_fasta(fasta_file)
    inversions = parse_lastz(lastz_file, min_length)

    pos_path = os.path.join(output_dir, "palindrome_positions.tsv")
    seq_path = os.path.join(output_dir, "palindrome_sequences.tsv")

    pos_out = open(pos_path, "w")
    seq_out = open(seq_path, "w")

    pos_out.write("pal_id\tchrom\tstart\tend\tlength\tstrand\n")
    seq_out.write("pal_id\tsequence\n")

    repeat_id = 1

    for chrom in inversions:

        groups = group_overlapping(inversions[chrom])

        for group in groups:

            shortest = min(group, key=lambda x: x[1] - x[0])
            s_short, e_short = shortest
            unit_len = (e_short - s_short) // 2
            if unit_len == 0:
                continue

            sequence = genome[chrom][s_short:s_short+unit_len].upper()
            pal_id = f"pal_{repeat_id}"
            seq_out.write(f"{pal_id}\t{sequence}\n")

            # tile across each inversion separately
            for inv_start, inv_end in group:

                midpoint = (inv_start + inv_end) // 2

                pos = inv_start
                tile_index = 0

                while pos + unit_len <= inv_end:

                    tile_start = pos
                    tile_end   = pos + unit_len

                    if tile_supported(tile_start, tile_end, group):

                        # strand assignment
                        if tile_start < midpoint:
                            strand = "+"
                        else:
                            strand = "-"

                        pos_out.write(
                            f"{pal_id}\t{chrom}\t{tile_start}\t{tile_end}\t{unit_len}\t{strand}\n"
                        )

                    pos += unit_len
                    tile_index += 1

            repeat_id += 1

    pos_out.close()
    seq_out.close()


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
