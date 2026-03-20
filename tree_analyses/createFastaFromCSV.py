import os
import sys
from collections import defaultdict
import pandas as pd

def main(input_csv, output_path, prefix):
    # Detect delimiter and read the input file
    with open(input_csv, mode="r") as file:
        first_line = file.readline()
        delimiter = "\t" if "\t" in first_line else ","

    df = pd.read_csv(input_csv, delimiter=delimiter)
    df = df[df.columns].loc[df[df.columns[0]] != df.columns[0]]  # Remove duplicate header rows

    # Ensure 'Sequence' column exists
    if 'Sequence' not in df.columns:
        raise KeyError("The input CSV does not contain a 'Sequence' column.")

    # Remove duplicate sequences
    df = df.drop_duplicates(subset=['Sequence'])

    # Organize data by GeneType
    gene_type_sequences = defaultdict(list)

    for _, row in df.iterrows():
        gene_type = row['GeneType']
        header = f"{prefix}.{row['Pos']}.{row['Contig']}.{gene_type}.{row['Productive']}.{row['Strand']}"
        sequence = row['Sequence'].upper()

        gene_type_sequences[gene_type].append((header, sequence))

    # Write the sequences to separate FASTA files
    os.makedirs(output_path, exist_ok=True)

    for gene_type, sequences in gene_type_sequences.items():
        if gene_type == "V":
            fasta_file_path = os.path.join(output_path, f"{prefix}.fasta")
            with open(fasta_file_path, mode="w") as fasta_file:
                for header, sequence in sequences:
                    fasta_file.write(f">{header}\n{sequence}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_csv> <output_path> <prefix>")
        sys.exit(1)

    input_csv = sys.argv[1]
    output_path = sys.argv[2]
    prefix = sys.argv[3]

    main(input_csv, output_path, prefix)