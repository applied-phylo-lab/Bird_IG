#!/usr/bin/env python3
"""
add_inversion_density.py

Reads a gene annotation TSV, parses Source into Order/Species/Haplotype,
loads the corresponding IGH_self.tsv lastz self-alignment for each haplotype,
and computes per-gene inversion density (number of deduplicated inversions
overlapping the gene) as a new column.

Usage:
    python add_inversion_density.py \
        -i genes.tsv \
        -d /path/to/alignments \
        -o genes_with_inversion_density.tsv \
        -c 8
"""

import os
import argparse
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def read_tsv_with_header(path):
    """Read TSV even if header starts with #"""
    with open(path) as f:
        header = f.readline().strip()
        if header.startswith("#"):
            header = header[1:]  # remove leading #
        columns = header.split("\t")
    df = pd.read_csv(path, sep="\t", comment="#", names=columns, skiprows=1)
    return df


def deduplicate_inversions(inv_df):
    """Remove duplicate inversion pairs (A-B vs B-A)."""
    inv_df = inv_df.copy()
    inv_df["pair_key"] = inv_df.apply(
        lambda row: tuple(sorted([(row["start1"], row["end1"]),
                                  (row["start2"], row["end2"])])),
        axis=1
    )
    inv_df = inv_df.drop_duplicates(subset="pair_key").drop(columns="pair_key")
    return inv_df


def load_inversions(aln_path):
    """
    Load and deduplicate inversions from a lastz self-alignment TSV.

    lastz self-alignment TSVs typically have columns:
        score, name1, strand1, size1, zstart1, end1,
        name2, strand2, size2, zstart2, end2, identity, ...

    We keep only rows where strand1 != strand2 (i.e. inverted alignments)
    and rename zstart1->start1, zstart2->start2 for clarity.

    If your file uses different column names adjust the mapping below.
    """
    aln = read_tsv_with_header(aln_path)
    aln.columns = [c.replace("#", "").replace("%", "").replace("+", "") for c in aln.columns]
    # Normalise column names (strip whitespace)
    aln.columns = [c.strip() for c in aln.columns]

    # Identify coordinate columns robustly
    # lastz --format=general uses: zstart1, end1, zstart2, end2 (0-based half-open)
    col_map = {}
    for c in aln.columns:
        lc = c.lower()
        if lc in ("zstart1", "start1"):
            col_map["start1"] = c
        elif lc == "end1":
            col_map["end1"] = c
        elif lc in ("zstart2", "start2"):
            col_map["start2"] = c
        elif lc == "end2":
            col_map["end2"] = c
        elif lc == "strand1":
            col_map["strand1"] = c
        elif lc == "strand2":
            col_map["strand2"] = c

    required = {"start1", "end1", "start2", "end2", "strand1", "strand2"}
    missing = required - set(col_map.keys())
    if missing:
        raise ValueError(
            f"Alignment file {aln_path!r} is missing expected columns: {missing}.\n"
            f"Found columns: {list(aln.columns)}"
        )

    # Rename to canonical names
    rename = {v: k for k, v in col_map.items()}
    aln = aln.rename(columns=rename)

    # Convert coords to numeric
    for col in ("start1", "end1", "start2", "end2"):
        aln[col] = pd.to_numeric(aln[col], errors="coerce")

    aln = aln.dropna(subset=["start1", "end1", "start2", "end2"])

    # Keep only inverted alignments (opposite strands)
    inversions = aln[aln["strand1"] != aln["strand2"]].copy()

    inversions = deduplicate_inversions(inversions)
    return inversions


def count_overlapping_inversions(inversions, gene_start, gene_end):
    """
    Count how many deduplicated inversions overlap [gene_start, gene_end).

    An inversion (defined by either arm) overlaps the gene if either arm
    overlaps the gene window.  We check both arms.
    """
    if inversions.empty:
        return 0

    # Arm 1 overlaps gene
    arm1 = (inversions["start1"] < gene_end) & (inversions["end1"] > gene_start)
    # Arm 2 overlaps gene
    arm2 = (inversions["start2"] < gene_end) & (inversions["end2"] > gene_start)

    return int((arm1 | arm2).sum())


# ---------------------------------------------------------------------------
# Worker function (runs in subprocess)
# ---------------------------------------------------------------------------

def process_haplotype(args):
    """
    Load inversions for one haplotype and compute density for all its genes.

    Parameters
    ----------
    args : tuple
        (key, aln_path, gene_rows)
        key       : (order, species, haplotype)
        aln_path  : str  - path to IGH_self.tsv
        gene_rows : list of (original_index, gene_start, gene_end) tuples

    Returns
    -------
    tuple of (log_message, list of (original_index, inversion_density))
    """
    key, aln_path, gene_rows = args
    order, species, haplotype, contig = key

    empty = pd.DataFrame(
        columns=["start1", "end1", "start2", "end2", "strand1", "strand2"]
    )

    if not os.path.isfile(aln_path):
        msg = f"[WARNING] Alignment file not found, density will be 0: {aln_path}"
        inversions = empty
    else:
        try:
            inversions = load_inversions(aln_path)
            msg = f"[INFO] {aln_path}: {len(inversions)} deduplicated inversions loaded."
        except Exception as exc:
            msg = f"[WARNING] Could not load {aln_path}: {exc}. Density will be 0."
            inversions = empty

    results = []
    for idx, gene_start, gene_end in gene_rows:
        density = count_overlapping_inversions(inversions, gene_start, gene_end)
        results.append((idx, density))

    return msg, results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Annotate genes with inversion density from lastz self-alignments."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input TSV with gene annotations.",
    )
    parser.add_argument(
        "-d", "--input_dir", required=True,
        help="Root directory containing Order/Species/Haplotype subdirs.",
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output TSV path.",
    )
    parser.add_argument(
        "-c", "--cores", type=int, default=1,
        help="Number of parallel worker processes (default: 1).",
    )
    args = parser.parse_args()

    # -----------------------------------------------------------------------
    # Load gene table
    # -----------------------------------------------------------------------
    genes = pd.read_csv(args.input)
    original_columns = list(genes.columns)
    #filter genes to IGH only
    genes = genes[genes["Locus"] == "IGH"].copy()
    # Split Source -> Order, Species, Haplotype
    split_source = genes["Source"].str.split("/", expand=True)
    if split_source.shape[1] < 3:
        raise ValueError(
            "Source column does not contain at least 3 '/'-separated fields "
            "(expected Order/Species/Haplotype)."
        )
    genes["Order"] = split_source[0]
    genes["Species"] = split_source[1]
    genes["Haplotype"] = split_source[2]

    # Pre-compute gene coordinates once, avoiding repeated work in workers
    genes["_gene_start"] = genes["Pos"].astype(int)
    genes["_gene_end"] = genes["_gene_start"] + genes["Sequence"].astype(str).str.len()

    # -----------------------------------------------------------------------
    # Build one work unit per haplotype
    # -----------------------------------------------------------------------
    work_units = []
    for (order, species, haplotype, contig), group in genes.groupby(
    ["Order", "Species", "Haplotype", "Contig"], sort=False
    ):  
        key = (order, species, haplotype, contig)

        aln_path = os.path.join(
        args.input_dir, order, species, haplotype, f"{contig}_IGH.tsv"
        )
        # Pass lightweight tuples to subprocesses, not full DataFrames
        gene_rows = list(
            zip(group.index, group["_gene_start"], group["_gene_end"])
        )
        work_units.append((key, aln_path, gene_rows))

    print(
        f"[INFO] Processing {len(work_units)} haplotype(s) across "
        f"{len(genes)} gene(s) using {args.cores} core(s)."
    )

    # -----------------------------------------------------------------------
    # Parallel execution — one process per haplotype
    # -----------------------------------------------------------------------
    densities = {}  # original DataFrame index -> inversion_density

    with ProcessPoolExecutor(max_workers=args.cores) as pool:
        futures = {pool.submit(process_haplotype, wu): wu[0] for wu in work_units}

        for future in as_completed(futures):
            key = futures[future]
            order, species, haplotype, contig = key
            try:
                msg, results = future.result()
                print(msg)
                for idx, density in results:
                    densities[idx] = density
            except Exception as exc:
                print(
                    f"[ERROR] {order}/{species}/{haplotype}/{contig} failed: {exc}."
                    "Densities for its genes will be 0."
                )
                # Recover: fill zeros for all genes in the failed haplotype
                failed_unit = next(wu for wu in work_units if wu[0] == key)
                for idx, _, _ in failed_unit[2]:
                    densities[idx] = 0

    # -----------------------------------------------------------------------
    # Attach results and write output
    # -----------------------------------------------------------------------
    genes["inversion_density"] = genes.index.map(densities)

    out_cols = original_columns + ["inversion_density"]
    genes[out_cols].to_csv(args.output, sep="\t", index=False)

    print(f"[DONE] Output written to {args.output}  ({len(genes)} rows)")


if __name__ == "__main__":
    main()