#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import itertools
import multiprocessing as mp
from functools import partial

class LastZPairwiseAligner:
    def __init__(self, config=None):
        self.config = config
        self.lastz_params = '--step=20 --notransition --format=general:name1,strand1,start1,end1,length1,name2,strand2,start2+,end2+,length2,id%'

    def AlignTwoFasta(self, fasta1, fasta2, output_fname):
        cmd = f'lastz {fasta1} {fasta2} {self.lastz_params} --output={output_fname}'
        os.system(cmd)

    def GetAlignedDF(self, output_fname):
        return pd.read_csv(output_fname, sep='\t')


def run_pairwise_alignment(aligner, fasta1, fasta2, output_file):
    """Run LASTZ alignment for a pair of FASTAs if not already done."""
    if os.path.exists(output_file):
        print(f"Skipping existing alignment: {output_file}")
        return
    print(f"Aligning {os.path.basename(fasta1)} vs {os.path.basename(fasta2)}")
    aligner.AlignTwoFasta(fasta1, fasta2, output_file)


def main():
    parser = argparse.ArgumentParser(description="Self-align IGH loci and create BED files")
    parser.add_argument("-i", '--input_dir', required=True, help="Top-level input directory")
    parser.add_argument("-s", '--summary', required=True, help="Path to summary_table_IGH.tsv")
    parser.add_argument("-c", '--cores', type=int, default=4, help="Number of parallel processes")
    args = parser.parse_args()

    summary = pd.read_csv(args.summary, sep='\t')
    aligner = LastZPairwiseAligner()

    # Group by Order AND Species
    for order, order_df in summary.groupby('Order'):
        print(f"\n=== Processing Order: {order} ===")

        output_dir = os.path.join(args.input_dir, order, 'pairwise_alignments')
        os.makedirs(output_dir, exist_ok=True)

        for species, species_df in order_df.groupby('Species'):
            print(f"  --- Species: {species} ---")

            # Build mapping from key -> fasta path (within this species only)
            fasta_paths = {}
            for _, row in species_df.iterrows():
                fasta_path = os.path.join(
                    args.input_dir, order, row['Species'], row['Haplotype'],
                    'refined_ig_loci', 'igloci_fasta',
                    f"IGH_{row['Contig']}_{row['NumV']}Vs.fasta"
                )
                key = f"{row['Haplotype']}_{row['Contig']}"
                if os.path.exists(fasta_path):
                    fasta_paths[key] = fasta_path
                else:
                    print(f"  Missing fasta for {key}: {fasta_path}")

            keys = sorted(fasta_paths.keys())
            pairs = list(itertools.combinations_with_replacement(keys, 2))
            if not pairs:
                print(f"  No pairs for {species}")
                continue

            jobs = []
            for k1, k2 in pairs:
                output_file = os.path.join(output_dir, f"{k1}_vs_{k2}.txt")
                jobs.append((fasta_paths[k1], fasta_paths[k2], output_file))

            with mp.Pool(processes=args.cores) as pool:
                pool.starmap(partial(run_pairwise_alignment, aligner), jobs)

            print(f"  Finished {species}")


if __name__ == "__main__":
    main()