#!/usr/bin/env python3

import argparse
import gzip
import re

def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)

def infer_class_family(repeat_name):
    r = repeat_name.upper()

    if r.startswith("CR1"):
        return "LINE/CR1"
    if r.startswith("ERV"):
        return "LTR/ERV"
    if r.startswith("GYPSY"):
        return "LTR/Gypsy"
    if r.startswith("COPIA"):
        return "LTR/Copia"
    if r.startswith("HAT"):
        return "DNA/hAT"
    if r.startswith("TCMAR"):
        return "DNA/TcMar"
    if "SAT" in r:
        return "Satellite"
    if r.startswith("(") or "SIMPLE" in r:
        return "Simple_repeat"

    return "Other"

def parse_repeatmasker_gff(gff_file, bed_file):
    with open_maybe_gzip(gff_file) as infile, open(bed_file, "w") as out:
        for line in infile:
            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attrs = fields

            bed_start = int(start) - 1
            bed_end = int(end)

            # Extract Target repeat name
            m = re.search(r'Target\s+"Motif:([^"]+)"', attrs)
            if not m:
                continue

            repeat_name = m.group(1)
            class_family = infer_class_family(repeat_name)

            out.write(
                f"{chrom}\t{bed_start}\t{bed_end}\t"
                f"{repeat_name}\t{score}\t{strand}\t{class_family}\n"
            )

def main():
    parser = argparse.ArgumentParser(
        description="Convert RepeatMasker (Dfam/HMMER GFF) to BED"
    )
    parser.add_argument("-i", "--input", required=True,
                        help="RepeatMasker .out.gff or .gff.gz")
    parser.add_argument("-o", "--output", required=True,
                        help="Output BED file")

    args = parser.parse_args()
    parse_repeatmasker_gff(args.input, args.output)

if __name__ == "__main__":
    main()
