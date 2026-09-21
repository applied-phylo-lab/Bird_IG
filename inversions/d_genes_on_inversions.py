#!/usr/bin/env python3
"""
For each IGH haplotype in summary_features.csv, calculate the fraction of
D genes (from IGHD.csv) that fall within inverted regions (from the lastz
self-alignment {contig}_IGH.tsv, where strand2 == '-').

D gene positions are stored relative to the full contig; they are converted
to locus-relative (fasta) coordinates using the locus StartPos from the
per-haplotype refined_ig_loci/summary.csv, matching what IGH_self_alignment_bed.py does:
    gene_pos_fasta = gene_row['Pos'] - locus_start

Inversion coordinates in the self-alignment TSV are already in fasta space.

Usage:
    python d_genes_on_inversions.py \\
        -i /local/storage/kav67/clean_birds \\
        -s /local/storage/kav67/clean_birds/summary_features.csv \\
        -o /local/storage/kav67/clean_birds/D_inversions.tsv \\
        --min_inv_len 250
"""

import argparse
import os
import sys
import pandas as pd
from multiprocessing import Pool


def process_row(args_tuple):
    row, input_dir, min_inv_len = args_tuple

    order     = row['Order']
    species   = row['Species']
    haplotype = row['Haplotype']
    contig    = row['Contig']

    hap_dir = os.path.join(input_dir, order, species, haplotype)

    # ── 1. Load D genes for this haplotype ────────────────────────────────────
    ighd_path = os.path.join(hap_dir, 'IGHD.csv')
    if not os.path.exists(ighd_path):
        print(f"[SKIP] IGHD.csv not found: {ighd_path}", file=sys.stderr)
        return None

    ighd_df = pd.read_csv(ighd_path)
    ighd_df = ighd_df[ighd_df['Contig'] == contig].copy()

    if ighd_df.empty:
        print(f"[SKIP] No D genes on contig {contig} for {haplotype}", file=sys.stderr)
        return None

    # ── 2. Get locus start position to convert to fasta coordinates ───────────
    locus_summary_path = os.path.join(hap_dir, 'refined_ig_loci', 'summary.csv')
    if not os.path.exists(locus_summary_path):
        print(f"[SKIP] summary.csv not found: {locus_summary_path}", file=sys.stderr)
        return None

    locus_df = pd.read_csv(locus_summary_path)
    locus_row = locus_df[(locus_df['Contig'] == contig) & (locus_df['Locus'] == 'IGH')]

    if locus_row.empty:
        print(f"[SKIP] No IGH locus entry for contig {contig} in {locus_summary_path}",
              file=sys.stderr)
        return None

    locus_start = locus_row['StartPos'].iloc[0]

    # Convert D gene positions to fasta (locus-relative) coordinates
    ighd_df['pos_fasta'] = ighd_df['Pos'] - locus_start

    # ── 3. Load self-alignment and extract inversion intervals ────────────────
    lastz_path = os.path.join(hap_dir, f'{contig}_IGH.tsv')
    if not os.path.exists(lastz_path):
        print(f"[SKIP] lastz self-alignment not found: {lastz_path}", file=sys.stderr)
        return None

    aln_df = pd.read_csv(lastz_path, sep='\t', comment='#',
                         header=None,
                         names=['name1', 'strand1', 'start1', 'end1', 'length1',
                                'name2', 'strand2', 'start2', 'end2', 'length2', 'id_pct'])

    # Inversions: query aligns on the minus strand relative to the reference
    inv_df = aln_df[
        (aln_df['strand2'] == '-') &
        (aln_df['length2'] >= min_inv_len)
    ].copy()

    # ── 4. Check which D genes fall within any inversion interval ─────────────
    def on_inversion(pos):
        return any(
            (inv_df['start2'] <= pos) & (pos <= inv_df['end2'])
        )

    n_total = len(ighd_df)

    if inv_df.empty:
        n_on_inv = 0
    else:
        ighd_df['on_inversion'] = ighd_df['pos_fasta'].apply(on_inversion)
        n_on_inv = int(ighd_df['on_inversion'].sum())

    frac = n_on_inv / n_total if n_total > 0 else float('nan')

    print(f"[OK] {order}/{species}/{haplotype}/{contig}: "
          f"{n_on_inv}/{n_total} D genes on inversions ({frac:.2%})")

    return {
        'Order':            order,
        'Species':          species,
        'Haplotype':        haplotype,
        'Contig':           contig,
        'n_d_genes':        n_total,
        'n_on_inversion':   n_on_inv,
        'frac_on_inversion': frac,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Calculate the fraction of IGH D genes lying on inverted regions."
    )
    parser.add_argument('-i', '--input_dir', required=True,
                        help='Top-level input directory')
    parser.add_argument('-s', '--summary', required=True,
                        help='Path to summary_features.csv')
    parser.add_argument('-o', '--output', required=True,
                        help='Path for output TSV (e.g. D_inversions.tsv)')
    parser.add_argument('-c', '--cores', type=int, default=4,
                        help='Number of parallel processes (default: 4)')
    parser.add_argument('--min_inv_len', type=int, default=250,
                        help='Minimum inversion length to consider (default: 250)')

    args = parser.parse_args()

    df = pd.read_csv(args.summary)
    df = df[df['Locus'] == 'IGH']

    if df.empty:
        print('[ERROR] No IGH rows found in summary file.', file=sys.stderr)
        sys.exit(1)

    # Keep only the highest-NumV contig per haplotype (same logic as pipeline)
    df = (df.sort_values('NumV', ascending=False)
            .groupby(['Order', 'Species', 'Haplotype'], as_index=False)
            .first())

    print(f"Processing {len(df)} IGH haplotypes with {args.cores} workers "
          f"(min inversion length: {args.min_inv_len} bp)...")

    task_args = [(row, args.input_dir, args.min_inv_len)
                 for _, row in df.iterrows()]

    with Pool(args.cores) as pool:
        results = pool.map(process_row, task_args)

    results = [r for r in results if r is not None]

    if not results:
        print('[ERROR] No results produced.', file=sys.stderr)
        sys.exit(1)

    out_df = pd.DataFrame(results, columns=[
        'Order', 'Species', 'Haplotype', 'Contig',
        'n_d_genes', 'n_on_inversion', 'frac_on_inversion'
    ])
    out_df = out_df.sort_values(['Order', 'Species', 'Haplotype'])

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    out_df.to_csv(args.output, sep='\t', index=False)
    print(f"\nDone. Written {len(out_df)} rows to {args.output}")


if __name__ == '__main__':
    main()
