#!/usr/bin/env python3
import os
import subprocess
import argparse
import tempfile
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO

def run_lastz(contig_name, fasta_path, out_dir):
    """Run LASTZ on a single contig vs itself."""
    tmp_fa = os.path.join(out_dir, f"{contig_name}.fa")
    SeqIO.write([r for r in SeqIO.parse(fasta_path, "fasta") if r.id == contig_name], tmp_fa, "fasta")

    out_tsv = os.path.join(out_dir, f"{contig_name}_self.tsv")
    cmd = [
        "lastz", f"{tmp_fa}", f"{tmp_fa}",
        "--step=20","--notransition",
        "--format=general:name1,strand1,start1,end1,length1,name2,strand2,start2+,end2+,length2,id%"
    ]
    with open(out_tsv, "w") as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.DEVNULL)
    return out_tsv

def read_tsv_with_header(path):
    """Read TSV even if header starts with #"""
    with open(path) as f:
        header = f.readline().strip()
        if header.startswith("#"):
            header = header[1:]  # remove leading #
        columns = header.split("\t")
    df = pd.read_csv(path, sep="\t", comment="#", names=columns, skiprows=1)
    return df

def extract_inversions(tsv_file, out_bed):
    """Filter inversions (reverse strand) and save as BED."""
    df = read_tsv_with_header(tsv_file)
    df.columns = [c.replace("#", "").replace("%", "").replace("+", "") for c in df.columns]
    df = df[df['strand2'] == '-']
    #df_bed = df[['name1', 'start1', 'end1']].copy()
    return df
    #df_bed.to_csv(out_bed, sep="\t", header=False, index=False)


def main():
    parser = argparse.ArgumentParser(description="Self-align each contig and summarize inversion density per window")
    parser.add_argument("-i", "--input", required=True, help="Input genome FASTA")
    parser.add_argument("-o", "--output", required=True, help="Output combined BED file (inversion density)")
    parser.add_argument("-w", "--window", type=int, default=100000, help="Window size (default: 100000 bp)")
    parser.add_argument("-c", "--cores", type=int, default=4, help="Number of cores")
    args = parser.parse_args()

    os.makedirs(args.output+"tmp_self_align", exist_ok=True)
    contigs = [rec.id for rec in SeqIO.parse(args.input, "fasta")]

    print(f"Aligning {len(contigs)} contigs using {args.cores} cores...")

    # Run LASTZ for each contig in parallel
    with ProcessPoolExecutor(max_workers=args.cores) as ex:
        futures = {ex.submit(run_lastz, c, args.input, args.output+"tmp_self_align"): c for c in contigs}
        for fut in as_completed(futures):
            contig = futures[fut]
            try:
                fut.result()
            except Exception as e:
                print(f"Error with contig {contig}: {e}")
    print("LASTZ alignments complete.")

    # Extract inversions and combine
    combined_bed = args.output+"all_inversions.bed"
    with open(combined_bed, "w") as outfile:
        # write header to combined file
        outfile.write("name1\tstrand1\tstart1\tend1\tlength1\tname2\tstrand2\tstart2\tend2\tlength2\tid\n") 
        for contig in contigs:
            tsv = os.path.join(args.output,"tmp_self_align", f"{contig}_self.tsv")
            bed = os.path.join(args.output,"tmp_self_align", f"{contig}_self.bed")
            if os.path.exists(tsv):
                cur_bed = extract_inversions(tsv, bed)
                # add to combined file without header
                cur_bed.to_csv(outfile, sep="\t", header=False, index=False, mode="a")


    # Clean up
    # for f in bed_files + [combined_bed, windows_bed]:
    #     os.remove(f)
    # for f in os.listdir("tmp_self_align"):
    #     if f.endswith(".fa") or f.endswith(".tsv"):
    #         os.remove(os.path.join("tmp_self_align", f))

    print(f"Done. Output written to {args.output}")

if __name__ == "__main__":
    main()
