#!/usr/bin/env python3
"""Count inversions in pairwise locus alignments, binned by length (and
optionally identity).

Replaces the former shared_inversions_songbirds.py /
shared_inversions_house_finch.py / shared_inversions_within_species.py trio
(the latter two were byte-identical); the comparison scheme is now --mode.

Modes:
  all-pairs        every pair of haplotypes in the filtered set
                   (was shared_inversions_songbirds.py, with
                    --order Songbirds --primary-only --identity-bins)
  reference        one reference species against every other haplotype
                   (was shared_inversions_house_finch.py, with
                    --reference-species House_Finch)
  within-species   every pair of haplotypes belonging to the same species

This counts inversions per length/identity bin. It is deliberately separate
from shared_inversions.py, which computes reciprocal-overlap and Jaccard
metrics against the per-haplotype self-alignments.

Pairwise alignment files are looked up under --pairwise-dir, accepting both
the legacy "{hapA}_{hapB}.txt" naming and PatchWorkPlot's
"pair_{i}-{hapA}_{j}-{hapB}.tsv" naming, in either haplotype order.

Usage:
    python shared_inversions_counts.py -i INPUT_DIR -s IGH_filtered_table.tsv \
        --mode all-pairs --order Songbirds --primary-only --identity-bins \
        --pairwise-dir INPUT_DIR/Songbirds/patchworkplot/IGH/pairwise_alignments \
        -o INPUT_DIR/Songbirds/Songbirds_shared_inversions.tsv
"""

import os
import glob
import argparse
import itertools

import pandas as pd

LENGTH_BINS = [(1000, 2500, "1000-2500"),
               (2500, 5000, "2500-5000"),
               (5000, 7500, "5000-7500"),
               (7500, 10000, "7500-10000")]

IDENTITY_BINS = [(float('-inf'), 60, "<60%"),
                 (60, 70, "60-70%"),
                 (70, 80, "70-80%"),
                 (80, 90, "80-90%"),
                 (90, 100, "90-100%")]


def categorize(value, bins):
    """Bin a value; the last bin is closed on the right, the rest half-open."""
    if pd.isna(value):
        return None
    last_index = len(bins) - 1
    for i, (lo, hi, label) in enumerate(bins):
        upper_ok = value <= hi if i == last_index else value < hi
        if lo <= value and upper_ok:
            return label
    return None


def load_inversions(file_path, want_identity):
    """Load a pairwise LASTZ alignment and keep the inverted alignments
    (strand2 == '-'), binned by length and optionally identity."""
    try:
        df = pd.read_csv(file_path, sep='\t')
    except Exception as e:
        print(f"Could not read {file_path}: {e}")
        return pd.DataFrame()

    strand_cols = [c for c in df.columns if 'strand' in c.lower()]
    if len(strand_cols) < 2:
        print(f"No strand columns found in {file_path}")
        return pd.DataFrame()

    inv_df = df[df[strand_cols[1]] == '-'].copy()
    if inv_df.empty:
        return pd.DataFrame()

    if 'length1' in inv_df.columns:
        inv_df['length'] = inv_df['length1']
    elif 'length' in inv_df.columns:
        pass
    elif {'end1', 'start1'}.issubset(inv_df.columns):
        inv_df['length'] = (inv_df['end1'] - inv_df['start1']).abs()
    else:
        inv_df['length'] = 0

    inv_df['length_category'] = inv_df['length'].apply(
        lambda v: categorize(v, LENGTH_BINS))
    subset = ['length_category']

    if want_identity:
        id_cols = [c for c in inv_df.columns
                   if 'id' in c.lower() or 'identity' in c.lower()]
        if id_cols:
            id_col = id_cols[0]
            ident = pd.to_numeric(
                inv_df[id_col].astype(str).str.replace('%', '', regex=False),
                errors='coerce')
            inv_df['identity_category'] = ident.apply(
                lambda v: categorize(v, IDENTITY_BINS))
        else:
            inv_df['identity_category'] = None
        subset.append('identity_category')

    return inv_df.dropna(subset=subset)


def find_pairwise_file(pairwise_dir, hap1, hap2):
    """Resolve the alignment file for a haplotype pair, trying the legacy
    "{a}_{b}.txt" naming and PatchWorkPlot's "pair_*-{a}_*-{b}.tsv", in both
    haplotype orders."""
    for a, b in ((hap1, hap2), (hap2, hap1)):
        direct = os.path.join(pairwise_dir, f"{a}_{b}.txt")
        if os.path.exists(direct):
            return direct
        hits = glob.glob(os.path.join(pairwise_dir, f"pair_*-{a}_*-{b}.tsv"))
        if hits:
            return sorted(hits)[0]
    return None


def build_pairs(summary, mode, reference_species):
    """Return a list of ((species1, hap1), (species2, hap2)) to compare."""
    units = (summary[['Species', 'Haplotype']]
             .drop_duplicates()
             .itertuples(index=False, name=None))
    units = sorted(set(units))

    if mode == 'all-pairs':
        return list(itertools.combinations(units, 2))

    if mode == 'reference':
        refs = [u for u in units if u[0] == reference_species]
        if not refs:
            raise SystemExit(
                f"No haplotypes found for --reference-species {reference_species!r}")
        others = [u for u in units if u[0] != reference_species]
        return [(r, o) for r in refs for o in others]

    if mode == 'within-species':
        pairs = []
        for species in sorted({u[0] for u in units}):
            same = [u for u in units if u[0] == species]
            pairs.extend(itertools.combinations(same, 2))
        return pairs

    raise SystemExit(f"Unknown mode: {mode}")


def main():
    parser = argparse.ArgumentParser(
        description="Count inversions in pairwise locus alignments, binned by "
                    "length and optionally identity")
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Top-level input directory")
    parser.add_argument("-s", "--summary", required=True,
                        help="Path to IGH_filtered_table.tsv (or summary_features.csv)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output TSV file")
    parser.add_argument("--mode", required=True,
                        choices=["all-pairs", "reference", "within-species"],
                        help="Which haplotype pairs to compare")
    parser.add_argument("--reference-species",
                        help="Species used as the reference in --mode reference "
                             "(e.g. House_Finch)")
    parser.add_argument("--pairwise-dir",
                        help="Directory holding the pairwise alignment files. "
                             "Default: {input_dir}/{order}/pairwise_alignments, "
                             "which requires a single --order.")
    parser.add_argument("--order",
                        help="Restrict to a taxonomic order (substring match, "
                             "case-insensitive)")
    parser.add_argument("--species", nargs="+",
                        help="Restrict to these species")
    parser.add_argument("--primary-only", action="store_true",
                        help="Keep only haplotypes whose name ends in _pri")
    parser.add_argument("--identity-bins", action="store_true",
                        help="Also bin by percent identity (the songbirds "
                             "analysis used this; the house finch one did not)")

    args = parser.parse_args()

    if args.mode == "reference" and not args.reference_species:
        parser.error("--mode reference requires --reference-species")

    sep = ',' if args.summary.endswith('.csv') else '\t'
    summary = pd.read_csv(args.summary, sep=sep)

    if args.order:
        summary = summary[summary['Order'].str.contains(args.order, case=False,
                                                        na=False)]
    if args.species:
        summary = summary[summary['Species'].isin(args.species)]
    if args.primary_only:
        summary = summary[summary['Haplotype'].str.endswith("_pri")]

    if summary.empty:
        print("No haplotypes left after filtering.")
        return

    pairwise_dir = args.pairwise_dir
    if not pairwise_dir:
        orders = sorted(summary['Order'].unique())
        if len(orders) != 1:
            parser.error("--pairwise-dir is required when the filtered summary "
                         f"spans multiple orders ({', '.join(orders)})")
        pairwise_dir = os.path.join(args.input_dir, orders[0], "pairwise_alignments")

    if not os.path.isdir(pairwise_dir):
        print(f"No pairwise alignment folder found: {pairwise_dir}")
        return

    pairs = build_pairs(summary, args.mode, args.reference_species)
    print(f"Comparing {len(pairs)} haplotype pairs from "
          f"{summary['Species'].nunique()} species "
          f"(mode={args.mode}, pairwise_dir={pairwise_dir})")

    group_cols = ['length_category']
    if args.identity_bins:
        group_cols.append('identity_category')

    results = []
    missing = 0
    for (species1, hap1), (species2, hap2) in pairs:
        file_path = find_pairwise_file(pairwise_dir, hap1, hap2)
        if file_path is None:
            missing += 1
            print(f"No pairwise file for {hap1} and {hap2}")
            continue

        inv_df = load_inversions(file_path, args.identity_bins)
        if inv_df.empty:
            continue

        counts = (inv_df.groupby(group_cols).size()
                  .reset_index(name='num_inversions'))
        for _, row in counts.iterrows():
            rec = {
                "Species1": species1,
                "Species2": species2,
                "Haplotype1": hap1,
                "Haplotype2": hap2,
                "Inversion_Length_Category": row['length_category'],
            }
            if args.identity_bins:
                rec["Identity_Category"] = row['identity_category']
            rec["num_inversions"] = row['num_inversions']
            results.append(rec)

    if not results:
        print("No inversions counted -- nothing written.")
        return

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    pd.DataFrame(results).to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(results)} rows to {args.output} "
          f"({missing} pairs had no alignment file)")


if __name__ == "__main__":
    main()
