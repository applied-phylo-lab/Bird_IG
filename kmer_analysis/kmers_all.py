import os
import argparse
import pandas as pd
from multiprocessing import Pool
import subprocess

def process_row(row):
    input_dir = row['InputDir']
    order = row['Order']
    species = row['Species']
    haplotype = row['Haplotype']
    contig = row['Contig']
    numv = row['NumV']
    script = row['script']

    # Paths
    fasta_path = os.path.join(input_dir, order, species, haplotype,
                              'refined_ig_loci', 'igloci_fasta',
                              f'IGH_{contig}_{numv}Vs.fasta')
    output_path = os.path.join(input_dir, order, species, haplotype,
                               'kmers.csv')

    # Ensure fasta exists
    if not os.path.isfile(fasta_path):
        print(f"FASTA not found: {fasta_path}")
        return

    # Run Perl script
    with open(output_path, 'w') as out_file:
        subprocess.run(['perl', script, '--kmerlength', '21', fasta_path],
                    check=True, stdout=out_file)


def main():
    parser = argparse.ArgumentParser(description="Self-align IGH loci and create BED files")
    parser.add_argument("-i", '--input_dir', required=True, help="Top-level input directory")
    parser.add_argument("-s", '--summary', required=True, help="Path to summary_table_IGH.tsv")
    parser.add_argument("-c", '--cores', type=int, default=4, help="Number of parallel processes")
    parser.add_argument("-p", '--pearl_script', required=True, help="path to pearl script")
    
    args = parser.parse_args()
    input_dir = args.input_dir
    summary_file = args.summary
    n_cores = args.cores
    script = args.pearl_script

    # read summary table
    df = pd.read_csv(summary_file, sep='\t')
    df['InputDir'] = input_dir  # pass input_dir to each row
    df['script'] = script  # pass input_dir to each row
    print(f"Processing {len(df)} rows with {n_cores} cores...")
    # multiprocessing
    with Pool(n_cores) as pool:
        pool.map(process_row, [row for _, row in df.iterrows()])

if __name__ == "__main__":
    main()
