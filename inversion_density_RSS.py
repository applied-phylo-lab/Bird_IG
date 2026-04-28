#!/usr/bin/env python3
"""
add_inversion_density.py

Reads a gene annotation TSV, parses Source into Order/Species/Haplotype,
loads the corresponding {contig}_IGH.tsv lastz self-alignment for each
(haplotype, contig) pair, and computes two inversion metrics as new columns:

  inversion_density     - raw count of deduplicated inversions overlapping
                          the gene (either arm), in locus coordinates
  inversion_density_pct - percentile rank of that count among all genes
                          on the same contig (0-100), enabling cross-species
                          comparison

Gene positions (whole-contig coordinates) are converted to locus coordinates
by subtracting the locus StartPos read from:
    <input_dir>/<order>/<species>/<haplotype>/refined_ig_loci/summary.csv

Alignment files are at:
    <input_dir>/<order>/<species>/<haplotype>/<contig>_IGH.tsv

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
from scipy.stats import percentileofscore
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

    rename = {v: k for k, v in col_map.items()}
    aln = aln.rename(columns=rename)

    for col in ("start1", "end1", "start2", "end2"):
        aln[col] = pd.to_numeric(aln[col], errors="coerce")

    aln = aln.dropna(subset=["start1", "end1", "start2", "end2"])

    inversions = aln[aln["strand1"] != aln["strand2"]].copy()
    inversions = deduplicate_inversions(inversions)
    return inversions


def load_locus_start(summary_path, contig):
    """
    Read refined_ig_loci/summary.csv and return the StartPos for the given
    contig.  Raises ValueError if the contig is not found.
    """
    summary = pd.read_csv(summary_path)
    row = summary[summary["Contig"] == contig]
    if row.empty:
        raise ValueError(
            f"Contig {contig!r} not found in summary file {summary_path!r}.\n"
            f"Available contigs: {summary['Contig'].tolist()}"
        )
    if len(row) > 1:
        raise ValueError(
            f"Contig {contig!r} matches {len(row)} rows in {summary_path!r}; "
            "expected exactly one."
        )
    return int(row.iloc[0]["StartPos"])


def count_overlapping_inversions(inversions, gene_start, gene_end):
    """
    Count deduplicated inversions overlapping [gene_start, gene_end).

    Both coordinates are in locus space.  Either arm of an inversion
    overlapping the gene window counts.
    """
    if inversions.empty:
        return 0

    arm1 = (inversions["start1"] < gene_end) & (inversions["end1"] > gene_start)
    arm2 = (inversions["start2"] < gene_end) & (inversions["end2"] > gene_start)

    return int((arm1 | arm2).sum())


# ---------------------------------------------------------------------------
# Worker function (runs in subprocess)
# ---------------------------------------------------------------------------

def process_contig(args):
    """
    Load inversions for one (haplotype, contig) pair, convert gene coordinates
    to locus space, and compute inversion density for each gene.

    Parameters
    ----------
    args : tuple
        (key, aln_path, summary_path, contig, gene_rows)
        key          : (order, species, haplotype, contig)
        aln_path     : str  - path to {contig}_IGH.tsv
        summary_path : str  - path to refined_ig_loci/summary.csv
        contig       : str  - contig name (used to look up StartPos)
        gene_rows    : list of (original_index, gene_start, gene_end) tuples
                       in whole-contig coordinates

    Returns
    -------
    tuple of (log_messages, list of (original_index, inversion_density))
    """
    key, aln_path, summary_path, contig, gene_rows = args
    order, species, haplotype, _ = key
    msgs = []

    empty = pd.DataFrame(
        columns=["start1", "end1", "start2", "end2", "strand1", "strand2"]
    )

    # -- Load locus StartPos for coordinate conversion ----------------------
    try:
        locus_start = load_locus_start(summary_path, contig)
        msgs.append(f"[INFO] {order}/{species}/{haplotype}/{contig}: locus StartPos={locus_start}")
    except Exception as exc:
        msgs.append(f"[WARNING] Could not read locus start for {contig} from "
                    f"{summary_path}: {exc}. Density will be 0.")
        return "\n".join(msgs), [(idx, 0) for idx, _, _ in gene_rows]

    # -- Load inversions (already in locus coordinates) ---------------------
    if not os.path.isfile(aln_path):
        msgs.append(f"[WARNING] Alignment file not found, density will be 0: {aln_path}")
        inversions = empty
    else:
        try:
            inversions = load_inversions(aln_path)
            msgs.append(f"[INFO] {aln_path}: {len(inversions)} deduplicated inversions loaded.")
        except Exception as exc:
            msgs.append(f"[WARNING] Could not load {aln_path}: {exc}. Density will be 0.")
            inversions = empty

    # -- Compute density for each gene -------------------------------------
    results = []
    for idx, contig_start, contig_end in gene_rows:
        # Convert whole-contig coordinates -> locus coordinates
        locus_gene_start = contig_start - locus_start
        locus_gene_end   = contig_end   - locus_start
        density = count_overlapping_inversions(inversions, locus_gene_start, locus_gene_end)
        results.append((idx, density))

    return "\n".join(msgs), results


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
    #filter genes to IGH only
    genes = genes[genes["Locus"] == "IGH"].copy()

    original_columns = list(genes.columns)

    # Split Source -> Order, Species, Haplotype
    split_source = genes["Source"].str.split("/", expand=True)
    if split_source.shape[1] < 3:
        raise ValueError(
            "Source column does not contain at least 3 '/'-separated fields "
            "(expected Order/Species/Haplotype)."
        )
    genes["Order"]     = split_source[0]
    genes["Species"]   = split_source[1]
    genes["Haplotype"] = split_source[2]

    # Pre-compute gene coordinates in whole-contig space
    genes["_gene_start"] = genes["Pos"].astype(int)
    genes["_gene_end"]   = genes["_gene_start"] + genes["Sequence"].astype(str).str.len()

    # -----------------------------------------------------------------------
    # Build one work unit per (haplotype, contig)
    # -----------------------------------------------------------------------
    work_units = []
    for (order, species, haplotype, contig), group in genes.groupby(
        ["Order", "Species", "Haplotype", "Contig"], sort=False
    ):
        key          = (order, species, haplotype, contig)
        aln_path     = os.path.join(args.input_dir, order, species, haplotype,
                                    f"{contig}_IGH.tsv")
        summary_path = os.path.join(args.input_dir, order, species, haplotype,
                                    "refined_ig_loci", "summary.csv")
        gene_rows = list(
            zip(group.index, group["_gene_start"], group["_gene_end"])
        )
        work_units.append((key, aln_path, summary_path, contig, gene_rows))

    print(
        f"[INFO] Processing {len(work_units)} (haplotype, contig) pair(s) across "
        f"{len(genes)} gene(s) using {args.cores} core(s)."
    )

    # -----------------------------------------------------------------------
    # Parallel execution — one process per (haplotype, contig)
    # -----------------------------------------------------------------------
    densities = {}  # original DataFrame index -> inversion_density

    with ProcessPoolExecutor(max_workers=args.cores) as pool:
        futures = {pool.submit(process_contig, wu): wu[0] for wu in work_units}

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
                    f"[ERROR] {order}/{species}/{haplotype}/{contig} failed: {exc}. "
                    "Densities for its genes will be 0."
                )
                failed_unit = next(wu for wu in work_units if wu[0] == key)
                for idx, _, _ in failed_unit[4]:
                    densities[idx] = 0

    # -----------------------------------------------------------------------
    # Attach raw density
    # -----------------------------------------------------------------------
    genes["inversion_density"] = genes.index.map(densities)

    # -----------------------------------------------------------------------
    # Compute within-contig percentile rank
    # Each contig has its own coordinate space and alignment, so percentile
    # rank is computed per (haplotype, contig) group.
    # percentileofscore(..., kind="rank"): % of scores <= query score.
    # Ties receive the same (highest) percentile.
    # -----------------------------------------------------------------------
    pct_ranks = {}
    for (order, species, haplotype, contig), group in genes.groupby(
        ["Order", "Species", "Haplotype", "Contig"], sort=False
    ):
        scores = group["inversion_density"].tolist()
        for idx, score in zip(group.index, scores):
            pct_ranks[idx] = percentileofscore(scores, score, kind="rank")

    genes["inversion_density_pct"] = genes.index.map(pct_ranks)

    # -----------------------------------------------------------------------
    # Write output
    # -----------------------------------------------------------------------
    out_cols = original_columns + ["inversion_density", "inversion_density_pct"]
    genes[out_cols].to_csv(args.output, sep="\t", index=False)

    print(f"[DONE] Output written to {args.output}  ({len(genes)} rows)")


if __name__ == "__main__":
    main()