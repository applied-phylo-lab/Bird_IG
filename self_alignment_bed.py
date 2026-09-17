#!/usr/bin/env python3
"""Self-align IG loci with LASTZ and write gene-position BED files.

Replaces the former IGH_self_alignment_bed.py / IGL_self_alignment_bed.py pair;
the locus is now a flag.

For every row of the summary table matching --locus, this:
  1. runs a LASTZ self-alignment of the extracted locus FASTA, writing
     {Contig}_{Locus}.tsv next to the other per-haplotype files
  2. writes gene positions (locus-relative coordinates) as
     {Contig}_{Locus}.bed          -- all genes black
     {Contig}_{Locus}_strand.bed   -- minus-strand genes grey

Existing outputs are left alone, so re-running only fills in what is missing.

Usage:
    python self_alignment_bed.py -i INPUT_DIR -s summary_features.csv --locus IGH -c 20
"""

import os
import argparse
import subprocess
from multiprocessing import Pool

import pandas as pd

LASTZ_PARAMS = ['--step=20', '--notransition']
LASTZ_FORMAT = ('--format=general:name1,strand1,start1,end1,length1,'
                'name2,strand2,start2+,end2+,length2,id%')


def read_summary(path, locus):
    """Read the summary table (.csv comma separated, .tsv tab separated) and
    keep only rows for `locus`."""
    sep = ',' if path.endswith('.csv') else '\t'
    df = pd.read_csv(path, sep=sep)
    if 'Locus' in df.columns:
        df = df[df['Locus'] == locus].copy()
    return df


def run_lastz(lastz_bin, fasta, output):
    subprocess.run([lastz_bin, fasta, fasta, *LASTZ_PARAMS, LASTZ_FORMAT,
                    f'--output={output}'], check=True)


def gene_table_path(hap_dir, locus):
    """Prefer the filtered gene table, fall back to the raw one."""
    clean = os.path.join(hap_dir, f'combined_genes_{locus}_clean.txt')
    if os.path.exists(clean):
        return clean
    return os.path.join(hap_dir, f'combined_genes_{locus}.txt')


def write_bed(path, contig, annot_df, color_col):
    with open(path, 'w') as fh:
        for _, row in annot_df.iterrows():
            fh.write(f"{contig}\t{row['Start']}\t{row['End']}\tNA\tNA\t"
                     f"{row['Strand']}\tNA\tNA\t{row[color_col]}\n")


def process_row(row):
    input_dir = row['InputDir']
    locus     = row['Locus']
    lastz_bin = row['LastzBin']
    order     = row['Order']
    species   = row['Species']
    haplotype = row['Haplotype']
    contig    = row['Contig']
    numv      = row['NumV']

    hap_dir = os.path.join(input_dir, order, species, haplotype)
    label = f"{order}/{species}/{haplotype}/{contig} [{locus}]"

    fasta_path = os.path.join(hap_dir, 'refined_ig_loci', 'igloci_fasta',
                              f'{locus}_{contig}_{numv}Vs.fasta')
    summary_csv_path = os.path.join(hap_dir, 'refined_ig_loci', 'summary.csv')
    combined_genes_path = gene_table_path(hap_dir, locus)

    # ---- 1. LASTZ self-alignment -----------------------------------------
    lastz_output = os.path.join(hap_dir, f'{contig}_{locus}.tsv')
    if not os.path.exists(fasta_path):
        print(f"FASTA not found, skipping: {fasta_path}")
        return
    if row['Force'] or not os.path.exists(lastz_output):
        try:
            run_lastz(lastz_bin, fasta_path, lastz_output)
        except subprocess.CalledProcessError as e:
            print(f"lastz failed for {label}: {e}")
            return

    # ---- 2. BED files ----------------------------------------------------
    bed_path = os.path.join(hap_dir, f'{contig}_{locus}.bed')
    strand_bed_path = os.path.join(hap_dir, f'{contig}_{locus}_strand.bed')
    if not row['Force'] and os.path.exists(bed_path) and os.path.exists(strand_bed_path):
        return

    if not os.path.exists(combined_genes_path):
        print(f"Gene table not found, skipping BED for {label}")
        return

    igdetect_df = pd.read_csv(combined_genes_path, sep='\t')
    igdetect_df = igdetect_df[igdetect_df['Contig'] == contig]

    locus_df = pd.read_csv(summary_csv_path)
    locus_df = locus_df[locus_df['Contig'] == contig]

    if locus_df.empty or igdetect_df.empty:
        print(f"Skipping {label} because no data found.")
        return

    locus_start = locus_df['StartPos'].iloc[0]

    annot_df = pd.DataFrame({
        'Start':  igdetect_df['Pos'] - locus_start,
        'End':    igdetect_df['Pos'] - locus_start + igdetect_df['Sequence'].str.len(),
        'Strand': igdetect_df['Strand'],
    })
    annot_df['Color'] = '0,0,0'
    # minus strand grey, plus strand black
    annot_df['StrandColor'] = annot_df['Strand'].apply(
        lambda x: '211,211,211' if x == '-' else '0,0,0')

    write_bed(bed_path, contig, annot_df, 'Color')
    write_bed(strand_bed_path, contig, annot_df, 'StrandColor')

    print(f"Processed {label}")


def main():
    parser = argparse.ArgumentParser(
        description="Self-align IG loci with LASTZ and create BED annotation files")
    parser.add_argument('-i', '--input_dir', required=True,
                        help="Top-level input directory")
    parser.add_argument('-s', '--summary', required=True,
                        help="Path to summary_features.csv (or any table with "
                             "Order/Species/Haplotype/Locus/Contig/NumV)")
    parser.add_argument('--locus', choices=['IGH', 'IGL'], default='IGH',
                        help="Which locus to process (default: IGH)")
    parser.add_argument('-c', '--cores', type=int, default=4,
                        help="Number of parallel processes")
    parser.add_argument('--lastz', default='lastz',
                        help="Path to the lastz executable")
    parser.add_argument('--force', action='store_true',
                        help="Regenerate outputs even if they already exist")
    # Single-unit selection, so a workflow manager can drive one row at a time.
    parser.add_argument('--order', help="Only process this Order")
    parser.add_argument('--species', help="Only process this Species")
    parser.add_argument('--haplotype', help="Only process this Haplotype")
    parser.add_argument('--contig', help="Only process this Contig")

    args = parser.parse_args()

    df = read_summary(args.summary, args.locus)
    for col, val in (('Order', args.order), ('Species', args.species),
                     ('Haplotype', args.haplotype), ('Contig', args.contig)):
        if val is not None:
            df = df[df[col] == val]
    if df.empty:
        print(f"No {args.locus} rows found in {args.summary} matching the given filters")
        return

    df['InputDir'] = args.input_dir
    df['Locus'] = args.locus
    df['LastzBin'] = args.lastz
    df['Force'] = args.force

    print(f"Processing {len(df)} {args.locus} rows with {args.cores} workers...")

    with Pool(args.cores) as pool:
        pool.map(process_row, [row for _, row in df.iterrows()])


if __name__ == "__main__":
    main()
