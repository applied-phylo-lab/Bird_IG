#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import re
from collections import defaultdict

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


def get_diagonal_inversions(df, minlen=1000):
    """Return only diagonal inversions (palindromes)."""
    df["inversion"] = df["strand1"] != df["strand2"]
    inv = df[df["inversion"]].copy()
    inv_diag = inv[(inv["start1"] == inv["start2"]) & (inv["end1"] == inv["end2"])]
    return inv_diag[inv_diag["length1"] >= minlen]


def find_sister_groups(bed_df, inv_diag, include_singletons=True):
    """
    Find sister gene groups with union-find.
    Now also keeps singleton groups if include_singletons=True.
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
    starts = bed["start"].astype(int).values
    ends = bed["end"].astype(int).values
    strands = bed["strand"].astype(str).values

    for _, inv in inv_diag.iterrows():
        inv_start = int(inv["start1"])
        inv_end   = int(inv["end1"])
        if inv_end <= inv_start:
            continue

        inside_mask = (starts >= inv_start) & (ends <= inv_end)
        idxs_inside = [i for i, flag in enumerate(inside_mask) if flag]
        for i in idxs_inside:
            g_start, g_end, g_strand = starts[i], ends[i], strands[i]
            offset_start = g_start - inv_start
            offset_end   = g_end   - inv_start
            sister_start = inv_end - offset_end
            sister_end   = inv_end - offset_start

            cand_mask = (starts <= sister_end) & (ends >= sister_start) & (strands != g_strand)
            cand_idxs = [j for j, flag in enumerate(cand_mask) if flag]
            for j in cand_idxs:
                dsu.union(i, j)

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


def main():
    parser = argparse.ArgumentParser(description="Find sister gene groups from inversions")
    parser.add_argument("-i", "--indir", required=True,
                        help="Input directory containing self*.tsv files")
    parser.add_argument("-b", "--beddir", required=True,
                        help="Directory containing BED files")
    parser.add_argument("-o", "--out", default="sisters.tsv",
                        help="Output TSV file (sister gene groups)")
    args = parser.parse_args()

    out_rows = []
    for fname in os.listdir(args.indir):
        if not fname.startswith("self") or not fname.endswith(".tsv"):
            continue

        fpath = os.path.join(args.indir, fname)
        df = pd.read_csv(fpath, sep="\t")
        df.columns = [c.replace("#", "").replace("%", "").replace("+", "") for c in df.columns]

        bed_path = find_bed_file(fpath, args.beddir)
        if bed_path and os.path.exists(bed_path):
            bed_df = pd.read_csv(
                bed_path, sep="\t", header=None, usecols=[1, 2, 5],
                names=["start", "end", "strand"]
            )
        else:
            bed_df = pd.DataFrame(columns=["start", "end", "strand"])

        sample = shorten_name(fname)

        inv_diag = get_diagonal_inversions(df, minlen=1000)
        groups, group_members = find_sister_groups(bed_df, inv_diag)

        # write out groups for this sample
        
        for g in groups:
            out_rows.append({
                "sample": sample,
                "group_size": len(g),
                "group_members": ";".join(g)
            })


    pd.DataFrame(out_rows).to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
