#!/usr/bin/env python3
"""
Print the average length of inversions across all self TSVs in a patchworkplot folder.
Inversions are alignments where strand2 == "-". Each pair appears twice in self alignments,
so we deduplicate by keeping only rows where start1 <= start2.
"""

import sys
import os
import glob


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <patchworkplot_folder>", file=sys.stderr)
        sys.exit(1)

    folder = sys.argv[1]
    tsv_pattern = os.path.join(folder, "pairwise_alignments", "self_*.tsv")
    tsv_files = sorted(glob.glob(tsv_pattern))

    if not tsv_files:
        print(f"No self_*.tsv files found in {tsv_pattern}", file=sys.stderr)
        sys.exit(1)

    lengths = []

    for tsv in tsv_files:
        with open(tsv) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10:
                    continue
                strand2 = fields[6]
                if strand2 != "-":
                    continue
                start1 = int(fields[2])
                start2 = int(fields[7])
                if start1 > start2:
                    continue  # skip duplicate (keep only one of each pair)
                lengths.append(int(fields[4]))  # length1

    if not lengths:
        print("No inversions found.")
        return

    avg = sum(lengths) / len(lengths)
    print(f"Inversions found : {len(lengths)}")
    print(f"Average length   : {avg:.1f} bp")


if __name__ == "__main__":
    main()
