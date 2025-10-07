#!/usr/bin/env python3
import argparse
import os
import re
import pandas as pd
import numpy as np

def read_tsv_with_header(path):
    """Read TSV even if header starts with #"""
    with open(path) as f:
        header = f.readline().strip()
        if header.startswith("#"):
            header = header[1:]  # remove leading #
        columns = header.split("\t")
    df = pd.read_csv(path, sep="\t", comment="#", names=columns, skiprows=1)
    return df

def merge_intervals(intervals):
    """Merge overlapping intervals and return total covered length."""
    if not intervals:
        return 0, []
    intervals.sort(key=lambda x: x[0])
    merged = []
    start, end = intervals[0]
    for s, e in intervals[1:]:
        if s <= end:  # overlap
            end = max(end, e)
        else:
            merged.append((start, end))
            start, end = s, e
    merged.append((start, end))
    total_len = sum(e - s for s, e in merged)
    return total_len, merged

def deduplicate_inversions(inv_df):
    """Remove duplicate inversion pairs (A-B vs B-A)."""
    # Make canonical keys by sorting coords
    inv_df = inv_df.copy()
    inv_df["pair_key"] = inv_df.apply(
        lambda row: tuple(sorted([(row["start1"], row["end1"]),
                                  (row["start2"], row["end2"])])),
        axis=1
    )
    # Drop duplicates based on canonical pair
    inv_df = inv_df.drop_duplicates(subset="pair_key").drop(columns="pair_key")
    return inv_df



def summarize_table(df, minlen,bed_df):
    """Compute summary statistics for one alignment table at given minlen."""
    # filter by minlen

    df = df[(df["length1"] >= minlen) & (df["length2"] >= minlen)].copy()

    # total sequence length (first row, length1)
   #df100 = df[df["id"] == "100%"]
    total_seq_length = df["length1"].max()

    # inversion flag
    df["inversion"] = df["strand1"] != df["strand2"]

    # all inversions
    inv = df[df["inversion"]]

    num_inversions = len(inv)

    # inversions on diagonal (same coords)
    inv_diag = inv[(inv["start1"] == inv["start2"]) & (inv["end1"] == inv["end2"])]
    num_inversions_diag = len(inv_diag)

    # lengths of inversions
    inv_lengths = inv["length1"].tolist()
    avg_inv_len = np.mean(inv_lengths) if inv_lengths else 0
    sd_len = pd.Series(inv_lengths).std(ddof=1) if inv_lengths else 0  # sample standard deviation

    # covered region by diagonal inversions (merge overlaps)
    intervals = [(row["start1"], row["end1"]) for _, row in inv_diag.iterrows()]
    inv_cov_len, merged_intervals = merge_intervals(intervals)

    # --- gene overlap ---
    genes_on_inv = 0
    genes_pos = 0
    genes_neg = 0
    total_genes = len(bed_df) if bed_df is not None else 0
    if bed_df is not None and not bed_df.empty:
        for _, gene in bed_df.iterrows():
            g_start, g_end, strand = gene["start"], gene["end"], gene["strand"]
            for s, e in merged_intervals:
                if not (g_end < s or g_start > e):  # overlap
                    genes_on_inv += 1
                    if strand == "+":
                        genes_pos += 1
                    elif strand == "-":
                        genes_neg += 1
                    break  # gene counted once, move on
    

    return {
        "minlen": minlen,
        "num_inversions": num_inversions,
        "num_inversions_diag": num_inversions_diag,
        "avg_inv_len": avg_inv_len,
        "sd_inv_len": sd_len,
        "total_seq_length": total_seq_length,
        "inv_cov_len": inv_cov_len,
        "total_genes": total_genes,
        "genes_on_inv": genes_on_inv,
        "genes_on_inv_pos": genes_pos,
        "genes_on_inv_neg": genes_neg,
        "frac_genes_on_inv": genes_on_inv / total_genes if total_genes > 0 else 0,
    }

def shorten_name(filename):
    """
    Extract the sample code between the species code and '_pri'.
    Handles prefixes like 'self_' or '14-GCA_'.
    Examples:
      14-GCA_048181955.1_A.insularis_AI_1833_00687_pri_1.0_genomic.tsv
        -> 'AI_1833_00687'
      7-GCA_048174325.1_A.woodhouseii_AW_366494_pri_1.0_genomic.tsv
        -> 'AW_366494'
      GCA_048182105.1_A.insularis_AI_1603_79203_pri_1.0_genomic.bed
        -> 'AI_1603_79203'
    """
    base = os.path.basename(filename)
    noext = re.sub(r"\.(tsv|bed)$", "", base, flags=re.IGNORECASE)

    # Capture the chunk between species code and _pri
    m = re.search(r'_([A-Z]{2}_[0-9_]+)_pri', noext)
    if m:
        return m.group(1)  # returns e.g. "AI_1603_79203"
    else:
        return noext

def find_bed_file(tsv_file, bed_dir):
    """Find the matching BED file for a given self*.tsv file."""
    base = os.path.basename(tsv_file)
    # Remove "self_*-" prefix
    no_prefix = re.sub(r"^self[^-]*-", "", base)
    bed_name = re.sub(r"\.tsv$", ".bed", no_prefix, flags=re.IGNORECASE)
    return os.path.join(bed_dir, bed_name)

def main():
    parser = argparse.ArgumentParser(description="Summarize self-alignment inversion stats.")
    parser.add_argument("-i", "--indir",required=True, help="Input directory containing self*.tsv files")
    parser.add_argument("-b", "--beddir", required=True, help="Directory containing BED files")
    parser.add_argument("-o", "--out",default="summary.tsv", help="Output TSV file")
    args = parser.parse_args()

    minlens = [1000, 2500, 5000, 7500, 10000, 12500, 15000]

    all_stats = []

    for fname in os.listdir(args.indir):
        if not fname.startswith("self") or not fname.endswith(".tsv"):
            continue
        fpath = os.path.join(args.indir, fname)

        # read table
        df = pd.read_csv(fpath, sep="\t")
        # clean column names
        df.columns = [c.replace("#", "").replace("%", "").replace("+", "") for c in df.columns]
        # read corresponding bed file
        bed_path = find_bed_file(fpath, args.beddir)
        if os.path.exists(bed_path):
            bed_df = pd.read_csv(bed_path, sep="\t", header=None, usecols=[1, 2, 5],
                                 names=["start", "end", "strand"])
        else:
            bed_df = pd.DataFrame(columns=["start", "end", "strand"])

        sample_name = shorten_name(fname)

        for ml in minlens:
            stats = summarize_table(df, ml, bed_df)
            stats["sample"] = sample_name
            all_stats.append(stats)

    outdf = pd.DataFrame(all_stats)
    outdf.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    main()
