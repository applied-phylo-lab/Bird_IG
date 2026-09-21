#!/usr/bin/env python3
"""
Quantify hairpin (palindrome) signal at the centre of diagonal inversions.

For every IGH haplotype in the summary table, the diagonal inversions
(strand1 != strand2, start1 == start2, end1 == end2) are pulled out of the
LASTZ self-alignment and the sequence identity of the whole alignment is
compared with the identity of a small window at its centre (the putative
hairpin tip) and of a random window of the same size.

Note on input files:
  The pipeline alignment `{contig}_IGH.tsv` written by IGH_self_alignment_bed.py
  does NOT contain the aligned sequences (no text1/text2 columns), which this
  analysis needs. This script therefore runs LASTZ once more per locus with
  text1/text2 added to the output format and caches the result next to the
  other per-haplotype files as `{contig}_IGH_text.tsv`.

Usage:
    python hairpin.py \\
        -i /local/storage/kav67/clean_birds \\
        -s /local/storage/kav67/clean_birds/IGH_filtered_table.tsv \\
        -o /local/storage/kav67/clean_birds/palindromes.tsv \\
        -c 20
"""

import os
import glob
import shutil
import argparse
import subprocess
import pandas as pd
from multiprocessing import Pool
from Bio.Seq import Seq
import random

# LASTZ settings: same as IGH_self_alignment_bed.py, plus text1/text2 so the
# aligned (gapped) sequences are available.
LASTZ_PARAMS = ['--step=20', '--notransition']
LASTZ_FORMAT = ('--format=general:name1,strand1,start1,end1,length1,text1,'
                'name2,strand2,start2+,end2+,length2,text2,id%')


# ================================================================
# Read TSV with optional # header
# ================================================================
def read_tsv_with_header(path):
    """Read TSV even if header starts with #."""
    with open(path) as f:
        header = f.readline().strip()
        if header.startswith("#"):
            header = header[1:]
        columns = header.split("\t")

    df = pd.read_csv(path, sep="\t", comment="#", names=columns, skiprows=1)
    df.columns = [c.replace("#", "").replace("%", "").replace("+", "")
                  for c in df.columns]
    return df


def read_summary(path):
    """Read the summary table (.csv is comma separated, .tsv tab separated)."""
    sep = ',' if path.endswith('.csv') else '\t'
    df = pd.read_csv(path, sep=sep)
    if 'Locus' in df.columns:
        df = df[df['Locus'] == 'IGH'].copy()
    return df


# ================================================================
# Reverse complement
# ================================================================
def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())


# ================================================================
# Identity between two aligned sequences
# Gaps are counted as aligned but not as matches
# ================================================================
def alignment_identity(seq1, seq2):
    """
    Compute the identity between two aligned sequences.
    Only columns where both seq1 and seq2 have valid A/C/G/T are considered.
    Returns percentage of matches among valid columns.
    """
    matches = 0
    valid = 0

    valid_bases = set("ACGT")

    for a, b in zip(seq1, seq2):
        if a not in valid_bases or b not in valid_bases:
            continue  # skip gaps or non-ACGT
        valid += 1
        if a.upper() == b.upper():
            matches += 1

    if valid == 0:
        return 0.0

    return 100 * matches / valid


def lastz_percent_identity(seq1, seq2):
    """
    Reproduce LASTZ's percent_identical() behavior exactly.

    Only positions where both seq1 and seq2 have valid A/C/G/T are counted.
    Identity = (200*numMatches + denom) // (2*denom)
    """

    # LASTZ's nuc_to_bits encoding:
    nuc_to_bits = {
        'A': 0, 'C': 1, 'G': 2, 'T': 3,
        'a': 0, 'c': 1, 'g': 2, 't': 3,
    }

    numMatches = 0.0
    denom = 0.0

    # iterate over aligned sequences
    for c1, c2 in zip(seq1, seq2):

        b1 = nuc_to_bits.get(c1, -1)
        b2 = nuc_to_bits.get(c2, -1)

        # must be >= 0 → valid DNA base
        if b1 >= 0 and b2 >= 0:
            denom += 1
            if b1 == b2:
                numMatches += 1

    if denom == 0.0:
        return 0.0

    # LASTZ formula: integer rounded identity
    return round(100.0 * numMatches / denom, 1)



# ================================================================
# Extract middle window from aligned sequences
# ================================================================
def extract_middle(seq1, seq2, bp=100):
    L = len(seq1)
    mid = L // 2
    half = bp // 2
    start = max(0, mid - half)
    end = min(L, mid + half)
    return seq1[start:end], seq2[start:end]

def extract_random_window(seq1, seq2, bp=50):
    L = len(seq1)
    if L <= bp:
        return seq1, seq2
    start = random.randint(0, L - bp)
    end = start + bp
    return seq1[start:end], seq2[start:end]


# ================================================================
# LASTZ self-alignment including the aligned sequences (text1/text2)
# ================================================================
def find_locus_fasta(hap_dir, contig, numv):
    """Locate the IGH locus FASTA for this contig."""
    fasta_dir = os.path.join(hap_dir, 'refined_ig_loci', 'igloci_fasta')
    if numv is not None:
        exact = os.path.join(fasta_dir, f'IGH_{contig}_{numv}Vs.fasta')
        if os.path.exists(exact):
            return exact
    # NumV in the summary table may not match the file name, fall back to glob
    hits = sorted(glob.glob(os.path.join(fasta_dir, f'IGH_{contig}_*Vs.fasta')))
    return hits[0] if hits else None


def self_align_with_text(hap_dir, contig, numv, lastz_bin, force=False):
    """
    Return the path to `{contig}_IGH_text.tsv`, running LASTZ if needed.
    Returns (path, message); path is None if the alignment could not be made.
    """
    out_path = os.path.join(hap_dir, f'{contig}_IGH_text.tsv')
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0 and not force:
        return out_path, None

    fasta = find_locus_fasta(hap_dir, contig, numv)
    if fasta is None:
        return None, f"locus FASTA not found for {contig} in {hap_dir}"

    tmp_path = out_path + '.tmp'
    cmd = [lastz_bin, fasta, fasta] + LASTZ_PARAMS + [LASTZ_FORMAT,
                                                      f'--output={tmp_path}']
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except FileNotFoundError:
        return None, f"lastz executable not found: {lastz_bin}"
    except subprocess.CalledProcessError as e:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        return None, f"lastz failed for {fasta}: {e.stderr.strip()}"

    os.replace(tmp_path, out_path)
    return out_path, None


# ================================================================
# Process one row of summary table
# ================================================================
def process_row(row):
    input_dir = row['InputDir']
    order = row['Order']
    species = row['Species']
    haplotype = row['Haplotype']
    contig = row['Contig']
    numv = row.get('NumV', None)
    lastz_bin = row['LastzBin']
    force = row['Force']

    hap_dir = os.path.join(input_dir, order, species, haplotype)

    # Self-alignment including the aligned sequences (text1/text2)
    tsv_path, err = self_align_with_text(hap_dir, contig, numv, lastz_bin, force)
    if tsv_path is None:
        print(f"[SKIP] {order}/{species}/{haplotype}/{contig}: {err}", flush=True)
        return []

    try:
        df = read_tsv_with_header(tsv_path)
    except Exception as e:
        print(f"[SKIP] {order}/{species}/{haplotype}/{contig}: cannot read {tsv_path}: {e}", flush=True)
        return []

    missing = [c for c in ("text1", "text2", "strand1", "strand2") if c not in df.columns]
    if missing:
        print(f"[SKIP] {order}/{species}/{haplotype}/{contig}: "
              f"{tsv_path} is missing column(s) {missing}", flush=True)
        return []

    # Filter diagonal inversions
    inv = df[(df["strand1"] != df["strand2"]) &
             (df["start1"] == df["start2"]) &
             (df["end1"] == df["end2"])]

    results = []
    for _, row2 in inv.iterrows():
        seq1 = str(row2["text1"]).upper()
        seq2 = str(row2["text2"]).upper()
        #seq2_rc = reverse_complement(seq2)
        id = str(row2["id"]).replace("%", "")

        whole_id = lastz_percent_identity(seq1, seq2)

        mid1_50, mid2_50 = extract_middle(seq1, seq2, bp=50)
        mid_id_50 = lastz_percent_identity(mid1_50, mid2_50)

        mid1_20, mid2_20 = extract_middle(seq1, seq2, bp=20)
        mid_id_20 = lastz_percent_identity(mid1_20, mid2_20)

        mid1_10, mid2_10 = extract_middle(seq1, seq2, bp=10)
        mid_id_10 = lastz_percent_identity(mid1_10, mid2_10)

        mid1_15, mid2_15 = extract_middle(seq1, seq2, bp=15)
        mid_id_15 = lastz_percent_identity(mid1_15, mid2_15)

        rand1, rand2 = extract_random_window(seq1, seq2, bp=50)
        rand_id_50 = lastz_percent_identity(rand1, rand2)
        rand1, rand2 = extract_random_window(seq1, seq2, bp=10)
        rand_id_10 = lastz_percent_identity(rand1, rand2)
        rand1, rand2 = extract_random_window(seq1, seq2, bp=20)
        rand_id_20 = lastz_percent_identity(rand1, rand2)
        rand1, rand2 = extract_random_window(seq1, seq2, bp=15)
        rand_id_15 = lastz_percent_identity(rand1, rand2)

        results.append({
            "Order": order,
            "Species": species,
            "Haplotype": haplotype,
            "Contig": contig,
            "Start": row2["start1"],
            "End": row2["end1"],
            "ID_Lastz"  : id,
            "WholeIdentity": whole_id,
            "MiddleIdentity10bp": mid_id_10,
            "MiddleIdentity15bp": mid_id_15,
            "MiddleIdentity20bp": mid_id_20,
            "MiddleIdentity50bp": mid_id_50,
            "MiddleMinusWhole": whole_id-mid_id_50,
            "RandomIdentity10bp": rand_id_10,
            "RandomIdentity15bp": rand_id_15,
            "RandomIdentity20bp": rand_id_20,
            "RandomIdentity50bp": rand_id_50,
            "MiddleSeq1": mid1_15,
            "MiddleSeq2": mid2_15,
        })

    print(f"[OK] {order}/{species}/{haplotype}/{contig}: {len(results)} diagonal inversions", flush=True)
    return results


# ================================================================
# Main script
# ================================================================
def main():
    parser = argparse.ArgumentParser(
        description="Summarize hairpin identity of diagonal inversions from "
                    "the IGH self-alignments ({contig}_IGH.tsv naming)")
    parser.add_argument('-i','--input_dir', required=True, help="Top-level input directory")
    parser.add_argument('-s','--summary', required=True,
                        help="Path to summary table (e.g. IGH_filtered_table.tsv)")
    parser.add_argument('-o','--output', required=True, help="Output TSV file")
    parser.add_argument('-c','--cores', type=int, default=4, help="Number of parallel workers")
    parser.add_argument('--lastz', default=shutil.which('lastz') or 'lastz',
                        help="Path to the lastz executable")
    parser.add_argument('--force', action='store_true',
                        help="Re-run lastz even if {contig}_IGH_text.tsv already exists")
    parser.add_argument('--seed', type=int, default=None,
                        help="Random seed for the random control windows")

    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    df = read_summary(args.summary)
    if df.empty:
        print("No IGH rows found in the summary table.")
        return

    df['InputDir'] = args.input_dir
    df['LastzBin'] = args.lastz
    df['Force'] = args.force

    print(f"Processing {len(df)} IGH haplotype/contig entries with {args.cores} workers...")

    # Run in parallel
    with Pool(args.cores) as pool:
        results = pool.map(process_row, [row for _, row in df.iterrows()])

    n_ok = sum(1 for sub in results if sub)
    all_rows = [r for sub in results for r in sub]
    if not all_rows:
        print("No valid data found.")
        return

    out_df = pd.DataFrame(all_rows)
    out_df.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(out_df)} rows from {n_ok}/{len(df)} entries to {args.output}")


if __name__ == "__main__":
    main()
