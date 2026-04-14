import os
import argparse
import pandas as pd


def create_configs(tsv_file, data_dir, out_dir, output_name):
    # Read TSV
    df = pd.read_csv(tsv_file, sep="\t")

    # Iterate over each species
    for (order, species), sub_df in df.groupby(["Order", "Species"]):

        rows = []

        for _, r in sub_df.iterrows():
            haplotype = r["Haplotype"]
            contig = r["Contig"]
            numV = r["NumV"]

            sample_id = f"{haplotype}.{contig}"
            label = haplotype

            # Paths now built from data_dir
            fasta_path = os.path.join(
                data_dir,
                order,
                species,
                haplotype,
                "refined_ig_loci",
                "igloci_fasta",
                f"IGH_{contig}_{numV}Vs.fasta"
            )

            annotation_path = os.path.join(
                data_dir,
                order,
                species,
                haplotype,
                f"{contig}_IGH_strand.bed"
            )

            rows.append({
                "SampleID": sample_id,
                "Label": label,
                "Fasta": fasta_path,
                "Annotation": annotation_path,
                "Strand": ""
            })

        # Output directory now uses out_dir
        species_out_dir = os.path.join(out_dir, order, species)
        os.makedirs(species_out_dir, exist_ok=True)

        out_file = os.path.join(species_out_dir, output_name)
        pd.DataFrame(rows).to_csv(out_file, index=False)

        print(f"Written: {out_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Create per-species config.csv files from a TSV."
    )

    parser.add_argument(
        "-i", "--input-tsv",
        required=True,
        help="Path to input TSV file"
    )

    parser.add_argument(
        "-d", "--data-dir",
        required=True,
        help="Base directory where FASTA and BED files are stored"
    )

    parser.add_argument(
        "-o", "--out-dir",
        required=True,
        help="Directory where config files will be written"
    )

    parser.add_argument(
        "-n", "--output-name",
        default="config.csv",
        help="Name of the output config file (default: config.csv)"
    )

    args = parser.parse_args()

    create_configs(
        tsv_file=args.input_tsv,
        data_dir=args.data_dir,
        out_dir=args.out_dir,
        output_name=args.output_name
    )


if __name__ == "__main__":
    main()