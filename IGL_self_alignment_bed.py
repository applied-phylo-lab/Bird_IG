import os
import pandas as pd
from multiprocessing import Pool
import argparse

class LastZPairwiseAligner:
    def __init__(self, config=None):
        self.config = config
        self.lastz_params = '--step=20 --notransition --format=general:name1,strand1,start1,end1,length1,name2,strand2,start2+,end2+,length2,id%'

    def AlignTwoFasta(self, fasta1, fasta2, output_dir, lastz_output):
        os.system(f'lastz {fasta1} {fasta2} {self.lastz_params} --output={lastz_output}')

    def GetAlignedDF(self, output_fname):
        return pd.read_csv(output_fname, sep='\t')

def process_row(row):
    input_dir = row['InputDir']
    order = row['Order']
    species = row['Species']
    haplotype = row['Haplotype']
    contig = row['Contig']
    numv = row['NumV']

    # Paths
    fasta_path = os.path.join(input_dir, order, species, haplotype,
                              'refined_ig_loci', 'igloci_fasta',
                              f'IGL_{contig}_{numv}Vs.fasta')
    summary_csv_path = os.path.join(input_dir, order, species, haplotype,
                                    'refined_ig_loci', 'summary.csv')
    combined_genes_path = os.path.join(input_dir, order, species, haplotype,
                                       'combined_genes_IGL_clean.txt')
    output_dir = os.path.join(input_dir, order, species, haplotype)

    # lastz self-alignment
    aligner = LastZPairwiseAligner()
    lastz_output = os.path.join(output_dir, f'{contig}_IGL.tsv')
    if os.path.exists(fasta_path) and not os.path.exists(lastz_output):
        aligner.AlignTwoFasta(fasta_path, fasta_path, output_dir, lastz_output)
    elif not os.path.exists(fasta_path):
        print(f"FASTA not found: {fasta_path}")
        return

    # create bed file
    # read gene table
    if not os.path.exists(combined_genes_path):
        combined_genes_path = os.path.join(input_dir, order, species, haplotype,
                                       'combined_genes_IGL.txt')
    igdetect_df = pd.read_csv(combined_genes_path, sep='\t')
    igdetect_df = igdetect_df[igdetect_df['Contig'] == contig]

    # read locus info
    locus_df = pd.read_csv(summary_csv_path)
    locus_df = locus_df[locus_df['Contig'] == contig]

    if locus_df.empty or igdetect_df.empty:
        print(f"Skipping {order}/{species}/{haplotype}/{contig} because no data found.")
        return

    bed_path        = os.path.join(output_dir, f'{contig}_IGL.bed')
    strand_bed_path = os.path.join(output_dir, f'{contig}_IGL_strand.bed')

    if os.path.exists(bed_path):
        return

    locus_start = locus_df['StartPos'].iloc[0]

    annot_df = {'Start': [], 'End': [], 'Strand': [], 'Color': []}
    for i, gene_row in igdetect_df.iterrows():
        gene_start_pos = gene_row['Pos'] - locus_start
        annot_df['Start'].append(gene_start_pos)
        annot_df['End'].append(gene_start_pos + len(gene_row['Sequence']))
        annot_df['Strand'].append(gene_row['Strand'])
        annot_df['Color'].append('0,0,0')

    annot_df = pd.DataFrame(annot_df)
    contig_name = contig
    annot_df_strand = annot_df.copy()
    # if strand is negative, color grey; if strand is positive, color black
    annot_df_strand['Color'] = annot_df_strand['Strand'].apply(lambda x: '211,211,211' if x == '-' else '0,0,0')
    with open(bed_path, 'w') as fh:
        for i, row_bed in annot_df.iterrows():
            fh.write(f"{contig_name}\t{row_bed['Start']}\t{row_bed['End']}\tNA\tNA\t{row_bed['Strand']}\tNA\tNA\t{row_bed['Color']}\n")
    with open(strand_bed_path, 'w') as fh:
        for i, row_bed in annot_df_strand.iterrows():
            fh.write(f"{contig_name}\t{row_bed['Start']}\t{row_bed['End']}\tNA\tNA\t{row_bed['Strand']}\tNA\tNA\t{row_bed['Color']}\n")

    print(f"Processed {order}/{species}/{haplotype}/{contig}")


def main():
    parser = argparse.ArgumentParser(description="Self-align IGL loci and create BED files")
    parser.add_argument("-i", '--input_dir', required=True, help="Top-level input directory")
    parser.add_argument("-s", '--summary', required=True, help="Path to summary CSV file")
    parser.add_argument("-c", '--cores', type=int, default=4, help="Number of parallel processes")

    args = parser.parse_args()
    input_dir = args.input_dir
    summary_file = args.summary
    n_cores = args.cores

    # read summary table, keep IGL rows only
    df = pd.read_csv(summary_file)
    df = df[df['Locus'] == 'IGL']
    df['InputDir'] = input_dir  # pass input_dir to each row

    # multiprocessing
    with Pool(n_cores) as pool:
        pool.map(process_row, [row for _, row in df.iterrows()])

if __name__ == "__main__":
    main()
