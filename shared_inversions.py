#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import itertools
from concurrent.futures import ProcessPoolExecutor, as_completed

def categorize_length(length):
    """Categorize inversion length into bins."""
    if 1000 <= length < 2500:
        return "1000-2500"
    elif 2500 <= length < 5000:
        return "2500-5000"
    elif 5000 <= length < 7500:
        return "5000-7500"
    elif 7500 <= length <= 10000:
        return "7500-10000"
    else:
        return None  # out of range

def load_inversions(file_path):
    """
    Load inversion data.
    Expects columns: start1, end1, start2+, end2+, id%
    and identifies inversions as regions where strand2 is '-'
    """
    try:
        df = pd.read_csv(file_path, sep='\t')
    except Exception as e:
        print(f"Could not read {file_path}: {e}")
        return pd.DataFrame()

    # try to find strand columns
    strand_cols = [c for c in df.columns if 'strand' in c.lower()]
    if len(strand_cols) < 2:
        print(f"No strand columns found in {file_path}")
        return pd.DataFrame()

    # Inversions typically have opposite strand orientation
    inv_df = df[df[strand_cols[1]] == '-'].copy()

    if 'length1' in df.columns:
        inv_df['length'] = df['length1']
    elif 'length' in df.columns:
        inv_df['length'] = df['length']
    else:
        inv_df['length'] = abs(df['end1'] - df['start1']) if {'end1','start1'}.issubset(df.columns) else 0

    inv_df['length_category'] = inv_df['length'].apply(categorize_length)
    inv_df = inv_df.dropna(subset=['length_category'])
    return inv_df

def overlaps(a_start, a_end, b_start, b_end):
    """Check if two intervals overlap."""
    return max(a_start, b_start) <= min(a_end, b_end)

def interval_overlap(a_start, a_end, b_start, b_end):
    """Return overlap length (bp) between intervals; 0 if none."""
    o_start = max(a_start, b_start)
    o_end = min(a_end, b_end)
    return max(0, o_end - o_start + 1)  # inclusive coordinates

def reciprocal_overlap_len(a_start, a_end, b_start, b_end):
    """Compute reciprocal overlap fraction: overlap/lenA and overlap/lenB."""
    overlap = interval_overlap(a_start, a_end, b_start, b_end)
    lenA = a_end - a_start + 1
    lenB = b_end - b_start + 1
    if lenA <= 0 or lenB <= 0:
        return 0.0, 0.0
    return overlap / lenA, overlap / lenB

def compute_shared_metrics(hap_df, other_df, pair_df,
                           min_overlap_bp=100, ro_thresh=0.5):
    """
    Compute three measures comparing hap_df (self inversions of haplotype A)
    to other_df (self inversions of haplotype B) using pair_df (pairwise alignments):
      - ro_count_fraction: fraction of inversions in hap_df that have >= ro_thresh reciprocal overlap with some pairwise inversion
      - shared_bases_fraction: (sum bp overlapped between hap_df inversions and pairwise inversions) / (total bp of hap_df inversions)
      - jaccard: shared_bp / (total_bp_hap_df + total_bp_other_df - shared_bp)
    All inputs should already be restricted to the same length_category if you want per-category measures.
    """
    # quick guards
    if hap_df.empty:
        return {
            "ro_count_fraction": 0.0,
            "shared_bases_fraction": 0.0,
            "jaccard": 0.0
        }

    # required coords in hap_df and pair_df: start1,end1 for the hap1 coordinate space
    if not {'start1', 'end1'}.issubset(hap_df.columns) or not {'start1', 'end1'}.issubset(pair_df.columns):
        return {"ro_count_fraction": 0.0, "shared_bases_fraction": 0.0, "jaccard": 0.0}

    # Precompute total bp in hap_df and in other_df (for jaccard denominator)
    total_bp_hap = (hap_df['end1'] - hap_df['start1'] + 1).sum()
    total_bp_other = 0
    if not other_df.empty and {'start1', 'end1'}.issubset(other_df.columns):
        total_bp_other = (other_df['end1'] - other_df['start1'] + 1).sum()

    shared_bp_sum = 0  # total bp from hap_df that are overlapped by pairwise inversion(s)
    ro_matched_count = 0

    # iterate inversions in hap_df (species A)
    for _, row_h in hap_df.iterrows():
        a_s, a_e = int(row_h['start1']), int(row_h['end1'])
        len_a = a_e - a_s + 1
        if len_a <= 0:
            continue

        # search pair_df for overlapping entries in hap1 coordinates
        # only consider rows in pair_df where strand2 == '-' (inversion in pair context)
        # but pair_df may have multiple columns; be robust about numeric types
        pair_cands = pair_df  # already filtered by length_category by caller

        # accumulate overlapping bp for this hap inversion
        hap_overlapped_bp = 0
        ro_hit = False

        for _, row_p in pair_cands.iterrows():
            # skip if explicit strand2 is present and not '-'
            if 'strand2' in pair_cands.columns:
                if str(row_p.get('strand2')).strip() != '-':
                    continue
            # coordinates in pair_df refer to hap1 coordinate space as start1/end1
            p_s, p_e = int(row_p['start1']), int(row_p['end1'])
            ov = interval_overlap(a_s, a_e, p_s, p_e)
            if ov <= 0:
                continue
            # optional min overlap filter
            if ov < min_overlap_bp:
                # still count toward overlapped bases but may not trigger RO
                hap_overlapped_bp += ov
                continue
            # compute reciprocal fractions
            fracA, fracP = reciprocal_overlap_len(a_s, a_e, p_s, p_e)
            hap_overlapped_bp += ov
            if fracA >= ro_thresh and fracP >= ro_thresh:
                ro_hit = True
                # break if you want one RO match per hap inversion (reciprocal best)
                break

        if ro_hit:
            ro_matched_count += 1
        shared_bp_sum += hap_overlapped_bp

    # metrics
    ro_count_fraction = ro_matched_count / len(hap_df) if len(hap_df) > 0 else 0.0
    shared_bases_fraction = shared_bp_sum / total_bp_hap if total_bp_hap > 0 else 0.0

    # For jaccard: shared_bp w.r.t union of hap_df and other_df (use shared_bp_sum as shared)
    union_bp = total_bp_hap + total_bp_other - shared_bp_sum
    jaccard = shared_bp_sum / union_bp if union_bp > 0 else 0.0

    return {
        "ro_count_fraction": ro_count_fraction,
        "shared_bases_fraction": shared_bases_fraction,
        "jaccard": jaccard
    }



def process_pair(args):
    """Worker function for one haplotype pair."""
    hap1, hap2, hap_info, pairwise_dir = args
    pair_file = os.path.join(pairwise_dir, f"{hap1}_{hap2}.txt")
    if not os.path.exists(pair_file):
        alt_file = os.path.join(pairwise_dir, f"{hap2}_{hap1}.txt")
        if os.path.exists(alt_file):
            pair_file = alt_file
        else:
            print(f"No pairwise alignment found for {hap1}, {hap2}")
            return []

    pair_df = load_inversions(pair_file)
    hap1_df = load_inversions(hap_info[hap1]['self_path'])
    hap2_df = load_inversions(hap_info[hap2]['self_path'])

    results = []
    for category in ["1000-2500", "2500-5000", "5000-7500", "7500-10000"]:
        hap1_cat = hap1_df[hap1_df['length_category'] == category]
        hap2_cat = hap2_df[hap2_df['length_category'] == category]
        pair_cat = pair_df[pair_df['length_category'] == category]

        # shared fraction: how many inversions in hap1 are found in pairwise inversion regions
        # compute metrics for hap1 (how many hap1 inversions are found in the pairwise df)
        metrics_hap1 = compute_shared_metrics(hap1_cat, hap2_cat, pair_cat,
                                      min_overlap_bp=100, ro_thresh=0.5)
        # compute metrics for hap2 (how many hap2 inversions are found in the pairwise df)
        metrics_hap2 = compute_shared_metrics(hap2_cat, hap1_cat, pair_cat,
                                            min_overlap_bp=100, ro_thresh=0.5)

        # create result row(s) — choose averaging or report both
        results.append({
            "Species1": hap_info[hap1]["species"],
            "Species2": hap_info[hap2]["species"],
            "Haplotype1": hap1,
            "Haplotype2": hap2,
            "Inversion_Length_Category": category,
            "RO_fraction_hap1": metrics_hap1['ro_count_fraction'],
            "RO_fraction_hap2": metrics_hap2['ro_count_fraction'],
            "RO_fraction_avg": (metrics_hap1['ro_count_fraction'] + metrics_hap2['ro_count_fraction']) / 2,
            "Shared_bases_fraction_hap1": metrics_hap1['shared_bases_fraction'],
            "Shared_bases_fraction_hap2": metrics_hap2['shared_bases_fraction'],
            "Shared_bases_fraction_avg": (metrics_hap1['shared_bases_fraction'] + metrics_hap2['shared_bases_fraction']) / 2,
            "Jaccard": (metrics_hap1['jaccard'] + metrics_hap2['jaccard']) / 2  # symmetric-ish
        })

    return results

def main():
    parser = argparse.ArgumentParser(description="Calculate shared inversion fractions per order (parallelized)")
    parser.add_argument("-i", "--input_dir", required=True, help="Top-level input directory")
    parser.add_argument("-s", "--summary", required=True, help="Path to IGH_filtered_table.tsv")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of parallel threads")
    args = parser.parse_args()

    summary = pd.read_csv(args.summary, sep='\t')

    for order, subdf in summary.groupby("Order"):
        if order == "MiscBirds":
            continue
        print(f"\n=== Processing Order: {order} ===")
        order_dir = os.path.join(args.input_dir, order)
        pairwise_dir = os.path.join(order_dir, "pairwise_alignments")

        if not os.path.exists(pairwise_dir):
            print(f"Missing pairwise dir for {order}")
            continue

        # Gather haplotype -> self alignment info
        hap_info = {}
        for _, row in subdf.iterrows():
            self_path = os.path.join(args.input_dir, order, row["Species"], row["Haplotype"], "IGH_self.tsv")
            if os.path.exists(self_path):
                hap_info[row["Haplotype"]] = {"species": row["Species"], "self_path": self_path}
            else:
                print(f"Missing self file for {row['Haplotype']}: {self_path}")

        haplotypes = sorted(hap_info.keys())
        pairs = list(itertools.combinations(haplotypes, 2))
        if not pairs:
            print(f"No valid haplotype pairs found for {order}")
            continue

        results_all = []
        with ProcessPoolExecutor(max_workers=args.threads) as executor:
            futures = {
                executor.submit(process_pair, (hap1, hap2, hap_info, pairwise_dir)): (hap1, hap2)
                for hap1, hap2 in pairs
            }
            for future in as_completed(futures):
                hap1, hap2 = futures[future]
                try:
                    res = future.result()
                    if res:
                        results_all.extend(res)
                except Exception as e:
                    print(f"Error processing pair {hap1}, {hap2}: {e}")

        if results_all:
            df_out = pd.DataFrame(results_all)
            out_path = os.path.join(order_dir, f"{order}_shared_inversions.tsv")
            df_out.to_csv(out_path, sep='\t', index=False)
            print(f"Saved: {out_path}")
        else:
            print(f"No results generated for {order}")

if __name__ == "__main__":
    main()
