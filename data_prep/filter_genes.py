#!/usr/bin/env python3
"""Filter combined_genes_IGH.txt and combined_genes_IGL.txt per index TSV.

Usage: python filter_genes.py --index index.tsv --input-dir /path/to/root

Writes cleaned files named combined_genes_IGH_clean.txt and combined_genes_IGL_clean.txt
next to the originals.

Filtering rules (implemented):
- keep sequences longer than `min_len` (default 250) (strictly greater)
- drop sequences where the most common nucleotide fraction >= `overrep` OR
  the sum of the two most common nucleotide fractions >= `overrep`.

Assumption: "1 or 2 nucleotides are overrepresented (70% or more)" means
either the top-1 or top-2 nucleotide frequency >= threshold (default 0.7).

"""
import argparse
import os
import sys
from typing import Iterable

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Filter combined_genes files by length and composition")
    p.add_argument("--index","-i", required=True, help="Path to index TSV with columns Order,Species,Haplotype,...")
    p.add_argument("--input-dir","-d", required=True, help="Root input directory that contains Order/Species/Haplotype/...")
    p.add_argument("--min-len", type=int, default=250, help="Minimum sequence length; keep sequences longer than this")
    p.add_argument("--overrep", type=float, default=0.7, help="Overrepresentation threshold (fraction) for 1 or 2 nucleotides")
    p.add_argument("--dry-run", action="store_true", help="Run quick internal test of filter functions and exit")
    return p.parse_args()


def seq_nt_fractions(seq: str) -> dict:
    """Return nucleotide fractions for A,C,G,T (uppercased). Non-ACGT characters are ignored for fraction denom.

    If the sequence contains no A/C/G/T, returns empty dict.
    """
    seq = (seq or "").upper()
    counts = {n: 0 for n in "ACGT"}
    total = 0
    for ch in seq:
        if ch in counts:
            counts[ch] += 1
            total += 1
    if total == 0:
        return {}
    return {n: counts[n] / total for n in counts}


def is_overrepresented(seq: str, threshold: float) -> bool:
    """Return True if seq is overrepresented by 1 or 2 nucleotides per threshold.

    Logic: compute fractions for A,C,G,T. If max fraction >= threshold OR
    sum of top two fractions >= threshold -> True
    """
    fr = seq_nt_fractions(seq)
    if not fr:
        return True  # no real bases -> treat as overrepresented / bad
    freqs = sorted(fr.values(), reverse=True)
    top1 = freqs[0]
    top2 = freqs[1] if len(freqs) > 1 else 0.0
    return (top1 >= threshold) or ((top1 + top2) >= threshold)


def filter_dataframe(df: pd.DataFrame, min_len: int, overrep: float) -> pd.DataFrame:
    """Apply length and overrepresentation filters to dataframe with a 'Sequence' column.

    Keeps rows where len(Sequence) > min_len and not is_overrepresented.
    """
    if "Sequence" not in df.columns:
        raise ValueError("DataFrame missing required 'Sequence' column")
    seqs = df["Sequence"].fillna("")
    lengths = seqs.str.len().astype(int)
    keep_len = lengths > min_len

    # vectorize overrepresentation check
    overrep_mask = seqs.apply(lambda s: is_overrepresented(s, overrep))
    keep = keep_len & (~overrep_mask)
    return df.loc[keep].reset_index(drop=True)


def process_row(index_row: pd.Series, input_dir: str, min_len: int, overrep: float) -> None:
    """Process a single row from the index TSV and write cleaned files if combined files exist."""
    order = str(index_row.get("Order") or index_row.get("order") or "")
    species = str(index_row.get("Species") or index_row.get("species") or "")
    hap = str(index_row.get("Haplotype") or index_row.get("haplotype") or index_row.get("Haplotype") or "")
    if not (order and species and hap):
        print(f"Skipping index row due to missing Order/Species/Haplotype: {index_row.to_dict()}")
        return

    base_dir = os.path.join(input_dir, order, species, hap)
    for locus in ("IGH", "IGL"):
        fname = f"combined_genes_{locus}.txt"
        in_path = os.path.join(base_dir, fname)
        if not os.path.isfile(in_path):
            # try alternative: sometimes files may be in a nested locus folder
            print(f"File not found, skipping: {in_path}")
            continue
        try:
            df = pd.read_csv(in_path, sep="\t", dtype=str, na_filter=False)
        except Exception as e:
            print(f"Failed to read {in_path}: {e}")
            continue

        try:
            cleaned = filter_dataframe(df, min_len=min_len, overrep=overrep)
        except Exception as e:
            print(f"Error filtering {in_path}: {e}")
            continue

        out_path = os.path.join(base_dir, f"combined_genes_{locus}_clean.txt")
        try:
            cleaned.to_csv(out_path, sep="\t", index=False)
            print(f"Wrote {len(cleaned)} rows -> {out_path}")
        except Exception as e:
            print(f"Failed to write {out_path}: {e}")


def main():
    args = parse_args()
    if args.dry_run:
        # small internal tests for the filter logic
        print("Dry run: testing filter logic")
        samples = [
            ("A" * 300, False),  # length >250 but all A -> overrepresented -> filtered
            ("ACGT" * 70, True),  # length 280, balanced -> kept
            ("N" * 300, False),  # no ACGT -> treated as bad -> filtered
            ("A" * 200 + "CG" * 30, False),  # length 260, but A dominates? A=200/260~=0.769 -> filtered
        ]
        for seq, expect_ok in samples:
            ok = (len(seq) > args.min_len) and (not is_overrepresented(seq, args.overrep))
            print(f"len={len(seq)} top-ok={ok} expected={expect_ok}")
        sys.exit(0)

    try:
        idx = pd.read_csv(args.index, sep="\t", dtype=str, na_filter=False)
    except Exception as e:
        print(f"Failed to read index file {args.index}: {e}")
        sys.exit(2)

    if not ("Order" in idx.columns or "order" in idx.columns):
        print("Index file missing 'Order' column")
    if not ("Species" in idx.columns or "species" in idx.columns):
        print("Index file missing 'Species' column")
    if not ("Haplotype" in idx.columns or "haplotype" in idx.columns or "Haplotype" in idx.columns):
        print("Index file missing 'Haplotype' column")

    for i, row in idx.iterrows():
        process_row(row, args.input_dir, min_len=args.min_len, overrep=args.overrep)


if __name__ == "__main__":
    main()
