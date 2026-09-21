import sys
import os
import subprocess


def write_temp_fasta(seq_tsv, out_fasta):
    with open(seq_tsv) as f, open(out_fasta, "w") as out:
        next(f)  # skip header
        for line in f:
            pal_id, seq = line.strip().split("\t")
            out.write(f">{pal_id}\n{seq}\n")


def run_lastz(query_fa, genome_fa, output_file):
    cmd = [
        "lastz",
        genome_fa,
        query_fa,
        "--format=general:name1,start1,end1,strand1,name2,start2,end2,strand2,identity",
        "--strand=both",
        "--step=1"
    ]

    with open(output_file, "w") as out:
        subprocess.run(cmd, stdout=out)


def parse_lastz_alignment(lastz_file, output_tsv):

    with open(lastz_file) as f, open(output_tsv, "w") as out:

        out.write("pal_id\tchrom\tstart\tend\tstrand\tidentity\n")

        for line in f:
            parts = line.strip().split()
            if len(parts) < 9:
                continue

            chrom = parts[0]
            start = parts[1]
            end   = parts[2]
            strand = parts[7]
            pal_id = parts[4]
            identity = parts[8]

            out.write(f"{pal_id}\t{chrom}\t{start}\t{end}\t{strand}\t{identity}\n")


def main(seq_tsv, genome_fa, output_dir):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    temp_fasta = os.path.join(output_dir, "palindromes_temp.fa")
    lastz_out  = os.path.join(output_dir, "palindromes_lastz.txt")
    final_out  = os.path.join(output_dir, "palindrome_realigned_positions.tsv")

    print("Writing temporary FASTA...")
    write_temp_fasta(seq_tsv, temp_fasta)

    print("Running LASTZ...")
    run_lastz(temp_fasta, genome_fa, lastz_out)

    print("Parsing alignments...")
    parse_lastz_alignment(lastz_out, final_out)

    print("Done.")
    print("Output:", final_out)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python align_palindromes.py palindrome_sequences.tsv genome.fa output_dir")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3])
