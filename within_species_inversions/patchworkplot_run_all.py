import os
import argparse
import subprocess
import pandas as pd


def run_patchwork(tsv_file, input_dir, config_name, output_prefix, dry_run):
    # Read TSV
    df = pd.read_csv(tsv_file, sep="\t")

    # Loop over each species
    for (order, species), _ in df.groupby(["Order", "Species"]):

        species_dir = os.path.join(input_dir, order, species)
        config_path = os.path.join(species_dir, config_name)

        # Output prefix per species
        out_prefix = os.path.join(species_dir, output_prefix)

        cmd = [
            "python", "/home/kav67/PatchWorkPlot_new/PatchWorkPlot.py",
            "-i", config_path,
            "-o", out_prefix,
            "--show-annot",
            "--lower",
            "--lwidth", "2",
            "--min-len", "15000",
            "--cmap", "viridis",
            "--min-pi", "75",
            "--transparent"
        ]

        print(f"\nRunning for {order} / {species}")
        print(" ".join(cmd))

        if not dry_run:
            try:
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error running PatchWorkPlot for {species}: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Run PatchWorkPlot.py for each species based on TSV input."
    )

    parser.add_argument(
        "-i", "--input-tsv",
        required=True,
        help="Path to input TSV file"
    )

    parser.add_argument(
        "-d", "--input-dir",
        required=True,
        help="Base directory containing Order/Species folders"
    )

    parser.add_argument(
        "-c", "--config-name",
        default="config.csv",
        help="Name of config file (default: config.csv)"
    )

    parser.add_argument(
        "-o", "--output-prefix",
        default="patchworkplot",
        help="Output prefix for PatchWorkPlot (default: patchworkplot)"
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing"
    )

    args = parser.parse_args()

    run_patchwork(
        tsv_file=args.input_tsv,
        input_dir=args.input_dir,
        config_name=args.config_name,
        output_prefix=args.output_prefix,
        dry_run=args.dry_run
    )


if __name__ == "__main__":
    main()