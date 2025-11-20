import os
import argparse
import pandas as pd
import numpy as np
from multiprocessing import Pool


def read_tsv_with_header(path):
    """Read TSV even if header starts with #"""
    with open(path) as f:
        header = f.readline().strip()
        if header.startswith("#"):
            header = header[1:]  # remosve leading #
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
        if s <= end:
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



def summarize_table(df, minlen, bed_df):
    """Compute summary statistics for one alignment table at given minlen."""

    df = df[(df["length1"] >= minlen) & (df["length2"] >= minlen)].copy()

    if df.empty:
        return {
            "minlen": minlen, "num_inversions": 0, "num_inversions_diag": 0,
            "avg_inv_len": 0, "sd_inv_len": 0, "total_seq_length": 0,
            "inv_cov_len": 0, "total_genes": len(bed_df),
            "genes_on_inv": 0, "genes_on_inv_pos": 0,
            "genes_on_inv_neg": 0, "frac_genes_on_inv": 0
        }

    # total sequence length (use max)
    total_seq_length = df["length1"].max()

    # inversion flag
    df["inversion"] = df["strand1"] != df["strand2"]
    inv = df[df["inversion"]]
    inv = deduplicate_inversions(inv)

    num_inversions = len(inv)
    inv_diag = inv[(inv["start1"] == inv["start2"]) & (inv["end1"] == inv["end2"])]
    num_inversions_diag = len(inv_diag)

    inv_lengths = inv["length1"].tolist()
    avg_inv_len = np.mean(inv_lengths) if inv_lengths else 0
    sd_len = pd.Series(inv_lengths).std(ddof=1) if inv_lengths else 0

    intervals = [(row["start1"], row["end1"]) for _, row in inv_diag.iterrows()]
    inv_cov_len, merged_intervals = merge_intervals(intervals)

    # gene overlap
    genes_on_inv = genes_pos = genes_neg = 0
    total_genes = len(bed_df) if bed_df is not None else 0

    if bed_df is not None and not bed_df.empty:
        for _, gene in bed_df.iterrows():
            g_start, g_end, strand = gene["start"], gene["end"], gene["strand"]
            for s, e in merged_intervals:
                if not (g_end < s or g_start > e):
                    genes_on_inv += 1
                    if strand == "+":
                        genes_pos += 1
                    elif strand == "-":
                        genes_neg += 1
                    break
    
    # make per-inversion DF for minlen=250
    inv_df = pd.DataFrame()
    if minlen == 250 and not inv.empty:
        inv_df = pd.DataFrame({
            "length": inv["length1"],
            "diagonal": (inv["start1"] == inv["start2"]) & (inv["end1"] == inv["end2"])
        })

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
    }, inv_df 



def process_row(row):
    input_dir = row['InputDir']
    order = row['Order']
    species = row['Species']
    haplotype = row['Haplotype']

    aln_path = os.path.join(input_dir, order, species, haplotype, 'IGH_self.tsv')
    bed_path = os.path.join(input_dir, order, species, haplotype, 'IGH.bed')

    if not os.path.exists(aln_path) or not os.path.exists(bed_path):
        print(f"Skipping {order}/{species}/{haplotype}: missing files")
        return [], pd.DataFrame()

    try:
        df = read_tsv_with_header(aln_path)
        df.columns = [c.replace("#", "").replace("%", "").replace("+", "") for c in df.columns]
        bed_df = pd.read_csv(bed_path, sep='\t', header=None)
        bed_df.columns = ['contig', 'start', 'end', 'na1', 'na2', 'strand', 'na3', 'na4', 'color']
    except Exception as e:
        print(f"Error reading {order}/{species}/{haplotype}: {e}")
        return [], pd.DataFrame()

    sample_name = f"{order}_{species}_{haplotype}"
    minlens = [250, 1000, 2500, 5000, 7500, 10000, 12500, 15000]
    all_stats = []
    inv_details_all = pd.DataFrame() 
    
    for ml in minlens:
        # try summarize_table, catch errors and print sample_name   
        try: 
            stats, inv_df = summarize_table(df, ml, bed_df)
        except Exception as e:
            print(f"Error processing {sample_name} at minlen {ml}: {e}")
            continue
        stats["sample"] = sample_name
        stats["order"] = order
        stats["species"] = species
        stats["haplotype"] = haplotype
        all_stats.append(stats)
    
        if ml == 250 and not inv_df.empty:
            inv_df["sample"] = sample_name
            inv_df["order"] = order
            inv_df["species"] = species
            inv_df["haplotype"] = haplotype
            inv_details_all = pd.concat([inv_details_all, inv_df], ignore_index=True)


    print(f"Processed {sample_name}")
    return all_stats,inv_details_all



def main():
    parser = argparse.ArgumentParser(description="Summarize inversion stats from LASTZ alignments and BED files")
    parser.add_argument('-i','--input_dir', required=True, help="Top-level input directory")
    parser.add_argument('-s','--summary', required=True, help="Path to summary_table_IGH.tsv")
    parser.add_argument('-o','--output', required=True, help="Output TSV file")
    parser.add_argument('-d','--details', required=True, help="Output TSV file (per inversion)")
    parser.add_argument('-c','--cores', type=int, default=4, help="Number of parallel workers")

    args = parser.parse_args()

    df = pd.read_csv(args.summary, sep='\t')
    df['InputDir'] = args.input_dir

    # run in parallel
    with Pool(args.cores) as pool:
        results = pool.map(process_row, [row for _, row in df.iterrows()])
    all_stats = []
    all_inv_details = []
   
    for stats, inv_details in results:
        if stats:  # if non-empty list
            all_stats.extend(stats)
        if not inv_details.empty:
            all_inv_details.append(inv_details)
    # flatten list of lists
    #all_stats = [item for sublist in results for item in sublist if sublist]
    #all_inv_details = pd.concat([inv for _, inv in results if not inv.empty], ignore_index=True)
   
    if all_stats:
        out_df = pd.DataFrame(all_stats)
        out_df.to_csv(args.output, sep='\t', index=False)
        print(f"Saved summary results to {args.output}")
    else:
        print("No summary results to save.")

    if all_inv_details:
        all_inv_details = pd.concat(all_inv_details, ignore_index=True)
        all_inv_details = all_inv_details[
            ["sample", "order", "species", "haplotype", "length", "diagonal"]
        ]
        all_inv_details.to_csv(args.details, sep='\t', index=False)
        print(f"Saved inversion details to {args.details}")
    else:
        print("No inversion details to save.")


if __name__ == "__main__":
    main()
