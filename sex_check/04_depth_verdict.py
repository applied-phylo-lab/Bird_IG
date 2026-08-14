#!/usr/bin/env python3
"""Windowed HiFi depth on the primary assembly -> is the Z diploid?

  female (ZW): Z is single copy -> ~0.5x the autosomal depth, plus a W also at 0.5x
  male   (ZZ): Z is diploid     -> ~1.0x, and no W at all

Why windows and not a per-contig mean: hifiasm leaves retained haplotigs in the
primary, and where one exists the reads split between the two near-identical
copies, so both read ~0.5x. A whole-contig mean therefore drags a perfectly
diploid Z down toward 0.5x and fakes a female. In bAgePho2 exactly this happens
-- ptg000020l duplicates ~30 Mb of the Z in ptg000019l -- so the verdict is
taken from the UNDUPLICATED part of the Z (its upper-mode window depth), and the
duplicated windows are reported separately rather than being averaged in.

A blanket "which contigs sit at half depth" sweep does not work here for the
same reason: retained haplotigs and repeat arrays fill that band.
"""
import statistics
import subprocess
from collections import defaultdict
from pathlib import Path

OUT = Path("/local/storage/kav67/RedwingedBlackbird_new/sex_check")
BAM = OUT / "depth_subsample.bam"
SAMTOOLS = "/programs/samtools-1.20/bin/samtools"
SUPER_Z, SUPER_W = "CM119746.1", "CM119744.1"
WINDOW = 5_000_000
MIN_LEN = 20_000_000     # contigs big enough for a stable windowed profile


def assignments():
    """contig -> 'Z' / 'W' from the primary-vs-female PAF."""
    paf = OUT / "bAgePho2.pri_vs_female.paf"
    if not paf.exists():
        raise SystemExit(f"missing {paf} -- run 01_map_to_ZW.sh first")
    bp = defaultdict(lambda: defaultdict(int))
    with open(paf) as fh:
        for line in fh:
            f = line.split("\t")
            if int(f[10]) >= 5000:
                bp[f[0]][f[5]] += int(f[3]) - int(f[2])
    out = {}
    for contig, d in bp.items():
        top, n = max(d.items(), key=lambda kv: kv[1])
        if n >= 1_000_000 and top in (SUPER_Z, SUPER_W):
            out[contig] = "Z" if top == SUPER_Z else "W"
    return out


def contig_lengths():
    lens = {}
    with open(OUT / "primary_contig_coverage.tsv") as fh:
        fh.readline()
        for line in fh:
            f = line.split("\t")
            lens[f[0]] = int(f[2])
    return lens


def window_depths(contig, length):
    out = []
    for start in range(1, length, WINDOW):
        stop = min(start + WINDOW - 1, length)
        if stop - start < WINDOW // 2:
            continue        # skip a short trailing window
        r = subprocess.run(
            [SAMTOOLS, "coverage", "-r", f"{contig}:{start}-{stop}", str(BAM)],
            capture_output=True, text=True, check=True)
        out.append(float(r.stdout.splitlines()[1].split("\t")[6]))
    return out


def main():
    if not BAM.exists():
        raise SystemExit(f"missing {BAM} -- run 03_read_depth.sh first")
    lab = assignments()
    lens = contig_lengths()

    autosomes = [c for c, L in lens.items() if L >= MIN_LEN and c not in lab]
    autosomes.sort(key=lambda c: -lens[c])
    zs = [c for c, L in lens.items() if L >= 5_000_000 and lab.get(c) == "Z"]
    ws = [c for c, L in lens.items() if L >= 5_000_000 and lab.get(c) == "W"]

    print("autosomal baseline (largest contigs):")
    base = []
    for c in autosomes[:5]:
        w = window_depths(c, lens[c])
        base += w
        print(f"  {c:<14}{lens[c]:>13,}  median {statistics.median(w):5.2f}x")
    baseline = statistics.median(base)
    print(f"  -> autosomal 2n depth = {baseline:.2f}x\n")

    def profile(contig, tag):
        w = window_depths(contig, lens[contig])
        rel = sorted(d / baseline for d in w)
        upper = statistics.median(rel[len(rel) // 2:])   # unduplicated part
        print(f"  {contig:<14}{lens[contig]:>13,}  {tag}")
        print("     windows (rel): " +
              " ".join(f"{r:.2f}" for r in sorted(rel, reverse=True)))
        print(f"     upper-mode relative depth = {upper:.2f}")
        return upper

    print("Z-assigned contigs:")
    z_upper = max(profile(c, "Z") for c in zs) if zs else None
    print("\nW-assigned contigs:")
    if ws:
        for c in ws:
            profile(c, "W")
    else:
        print("  (none >= 5 Mb)")

    print("\n" + "=" * 62)
    if z_upper is None:
        print("VERDICT: no Z contig found -- check 01/02 output")
    else:
        print(f"Z relative depth (unduplicated portion): {z_upper:.2f}x autosomal")
        if z_upper >= 0.80:
            print("VERDICT: MALE (ZZ) -- Z carries full diploid coverage")
        elif z_upper <= 0.65:
            print("VERDICT: FEMALE (ZW) -- Z is hemizygous")
        else:
            print("VERDICT: AMBIGUOUS -- inspect the window profiles above")
    print("=" * 62)


if __name__ == "__main__":
    main()
