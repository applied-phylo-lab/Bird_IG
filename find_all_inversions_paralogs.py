#!/usr/bin/env python3
import os
import argparse
import pandas as pd
from collections import defaultdict
from multiprocessing import Pool


def read_tsv_with_header(path):
    """Read TSV even if header starts with #"""
    with open(path) as f:
        header = f.readline().strip()
        if header.startswith("#"):
            header = header[1:]
        columns = header.split("\t")
    df = pd.read_csv(path, sep="\t", comment="#", names=columns, skiprows=1)
    df.columns = [c.replace("#", "").replace("%", "").replace("+", "") for c in df.columns]
    return df


class DSU:
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0]*n
    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x
    def union(self, x, y):
        rx, ry = self.find(x), self.find(y)
        if rx == ry:
            return
        if self.rank[rx] < self.rank[ry]:
            self.parent[rx] = ry
        elif self.rank[rx] > self.rank[ry]:
            self.parent[ry] = rx
        else:
            self.parent[ry] = rx
            self.rank[rx] += 1


def get_diagonal_inversions(df, minlen=1000):
    """Return only diagonal inversions (palindromes)."""
    df["inversion"] = df["strand1"] != df["strand2"]
    inv = df[df["inversion"]].copy()
    inv_diag = inv[(inv["start1"] == inv["start2"]) & (inv["end1"] == inv["end2"])]
    return inv_diag[inv_diag["length1"] >= minlen]


def get_non_diagonal_inversions(df, minlen=1000):
    """Return non-diagonal inversions (non-palindromic)."""
    df["inversion"] = df["strand1"] != df["strand2"]
    inv = df[df["inversion"]].copy()
    inv_non_diag = inv[(inv["start1"] != inv["start2"]) | (inv["end1"] != inv["end2"])]
    return inv_non_diag[inv_non_diag["length1"] >= minlen]


def connect_genes_within_inversion(dsu, bed, inv_start, inv_end):
    """
    For each gene on + strand inside the inversion, compute its mirrored position
    on the - strand and connect if a gene exists there (and vice versa).
    """
    starts = bed["start"].astype(int).values
    ends = bed["end"].astype(int).values
    strands = bed["strand"].astype(str).values

    n = len(bed)
    for i in range(n):
        g_start, g_end, g_strand = starts[i], ends[i], strands[i]
        # Only consider genes fully inside inversion
        if g_start < inv_start or g_end > inv_end:
            continue

        # Compute mirrored coordinates
        offset_start = g_start - inv_start
        offset_end   = g_end - inv_start
        sister_start = inv_end - offset_end
        sister_end   = inv_end - offset_start
        sister_strand = "-" if g_strand == "+" else "+"

        # Find candidate mirror genes
        for j in range(n):
            if i == j:
                continue
            s_start, s_end, s_strand = starts[j], ends[j], strands[j]
            if s_strand != sister_strand:
                continue
            # Check if gene overlaps the mirrored coordinates
            if s_start <= sister_end and s_end >= sister_start:
                dsu.union(i, j)


def find_sister_groups(bed_df, inv_diag, inv_non_diag, include_singletons=True):
    """
    Find sister gene groups with union-find.
    Handles both diagonal (palindromic) and non-diagonal inversions.
    """
    if bed_df is None or bed_df.empty:
        return [], {}

    bed = bed_df.copy().reset_index().rename(columns={"index": "orig_idx"})
    bed["gene_label"] = bed.apply(
        lambda r: f"{int(r['start'])}_{int(r['end'])}_{r['strand']}_{int(r['orig_idx'])}", axis=1
    )

    n = len(bed)
    label_to_idx = {lab: i for i, lab in enumerate(bed["gene_label"])}
    idx_to_label = {i: lab for lab, i in label_to_idx.items()}

    dsu = DSU(n)

    # Diagonal inversions (palindromic)
    for _, inv in inv_diag.iterrows():
        inv_start = int(inv["start1"])
        inv_end   = int(inv["end1"])
        if inv_end <= inv_start:
            continue
        connect_genes_within_inversion(dsu, bed, inv_start, inv_end)

    # Non-diagonal inversions
    for _, inv in inv_non_diag.iterrows():
        inv_start = int(inv["start1"])
        inv_end   = int(inv["end1"])
        if inv_end <= inv_start:
            continue
        connect_genes_within_inversion(dsu, bed, inv_start, inv_end)

    # Collect groups
    groups_map = defaultdict(list)
    for i in range(n):
        root = dsu.find(i)
        groups_map[root].append(i)

    groups = []
    group_members = {}
    for root, members in groups_map.items():
        if not include_singletons and len(members) == 1:
            continue
        labels = [idx_to_label[m] for m in sorted(members)]
        groups.append(sorted(labels))
        group_members[root] = [
            {
                "label": idx_to_label[m],
                "orig_idx": int(bed.loc[m, "orig_idx"]),
                "start": int(bed.loc[m, "start"]),
                "end": int(bed.loc[m, "end"]),
                "strand": bed.loc[m, "strand"]
            }
            for m in sorted(members)
        ]

    return groups, group_members


def process_row(row):
    input_dir = row['InputDir']
    order = row['Order']
    species = row['Species']
    haplotype = row['Haplotype']

    aln_path = os.path.join(input_dir, order, species, haplotype, 'IGH_self.tsv')
    bed_path = os.path.join(input_dir, order, species, haplotype, 'IGH.bed')

    if not os.path.exists(aln_path) or not os.path.exists(bed_path):
        print(f"Skipping {order}/{species}/{haplotype}: missing files")
        return []

    try:
        df = read_tsv_with_header(aln_path)
        bed_df = pd.read_csv(bed_path, sep='\t', header=None)
        bed_df.columns = ['contig', 'start', 'end', 'na1', 'na2', 'strand', 'na3', 'na4', 'color']
    except Exception as e:
        print(f"Error reading {order}/{species}/{haplotype}: {e}")
        return []

    sample_name = f"{order}_{species}_{haplotype}"

    inv_diag = get_diagonal_inversions(df, minlen=1000)
    inv_non_diag = get_non_diagonal_inversions(df, minlen=1000)

    groups, group_members = find_sister_groups(bed_df, inv_diag, inv_non_diag)

    out_rows = []
    for g in groups:
        out_rows.append({
            "Order": order,
            "Species": species,
            "Haplotype": haplotype,
            "sample": sample_name,
            "group_size": len(g),
            "group_members": ";".join(g)
        })
    print(f"Processed {sample_name}")
    return out_rows


def main():
    parser = argparse.ArgumentParser(description="Summarize inversion-based paralog groups from LASTZ alignments and BED files")
    parser.add_argument('-i','--input_dir', required=True, help="Top-level input directory")
    parser.add_argument('-s','--summary', required=True, help="Path to summary_table_IGH.tsv")
    parser.add_argument('-o','--output', required=True, help="Output TSV file")
    parser.add_argument('-c','--cores', type=int, default=4, help="Number of parallel workers")

    args = parser.parse_args()

    df = pd.read_csv(args.summary, sep='\t')
    df['InputDir'] = args.input_dir

    with Pool(args.cores) as pool:
        results = pool.map(process_row, [row for _, row in df.iterrows()])

    all_rows = [r for sub in results for r in sub]
    if not all_rows:
        print("No valid data found.")
        return

    out_df = pd.DataFrame(all_rows)
    out_df.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(out_df)} rows to {args.output}")


if __name__ == "__main__":
    main()
