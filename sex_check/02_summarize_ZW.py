#!/usr/bin/env python3
"""Call the sex of bAgePho2 from its alignment to a conspecific female genome.

Readout, per haplotype, is the breadth of reference sequence covered:

  SUPER_Z  (CM119746.1, 96.0 Mb) - a male (ZZ) has Z in BOTH haplotypes,
                                   a female (ZW) in only one
  SUPER_W  (CM119744.1, 24.6 Mb) - only a female has this

Two traps this script deliberately avoids, both learned the hard way here:

1. The PAR. The pseudoautosomal region sits in the first ~2 Mb of SUPER_W and
   is present in both sexes, so raw "bp aligned to W" is never zero in a male.
   W breadth is therefore also reported excluding that region.

2. The unlocalized W scaffolds. The reference's W_unloc_1/W_unloc_2
   (JBPJSG010000026.1 / JBPJSG010000024.1, 22.7 Mb combined) are repeat arrays
   that a female assembly binned onto W but which exist in both sexes -- in
   this bird they show up at similar coverage in BOTH haplotypes, i.e. diploid,
   i.e. not W-specific. They are reported but excluded from the verdict.
"""
from collections import defaultdict
from pathlib import Path

OUT = Path("/local/storage/kav67/RedwingedBlackbird_new/sex_check")
MIN_ALN = 5000

SUPER_Z = ("CM119746.1", 95_998_702)
SUPER_W = ("CM119744.1", 24_620_248)
PAR_END = 2_000_000          # PAR occupies roughly the first 2 Mb of SUPER_W
W_UNLOC = {"JBPJSG010000024.1": 12_197_094, "JBPJSG010000026.1": 10_509_783}


def merged_len(intervals, clip_start=0):
    intervals = sorted((max(s, clip_start), e) for s, e in intervals
                       if e > clip_start)
    if not intervals:
        return 0
    m = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= m[-1][1]:
            m[-1][1] = max(m[-1][1], e)
        else:
            m.append([s, e])
    return sum(e - s for s, e in m)


def load(hap):
    paf = OUT / f"bAgePho2.{hap}_vs_female.paf"
    if not paf.exists():
        raise SystemExit(f"missing {paf} -- run 01_map_to_ZW.sh first")
    iv = defaultdict(list)
    contigs = defaultdict(lambda: defaultdict(int))
    with open(paf) as fh:
        for line in fh:
            f = line.split("\t")
            if int(f[10]) < MIN_ALN:
                continue
            iv[f[5]].append((int(f[7]), int(f[8])))
            contigs[f[5]][f[0]] += int(f[3]) - int(f[2])
    return iv, contigs


def report(hap):
    iv, contigs = load(hap)

    z = merged_len(iv[SUPER_Z[0]])
    w_all = merged_len(iv[SUPER_W[0]])
    w_nopar = merged_len(iv[SUPER_W[0]], clip_start=PAR_END)

    print(f"\n===== {hap} haplotype =====")
    print(f"  SUPER_Z              {z:>12,} / {SUPER_Z[1]:>12,}  "
          f"({100*z/SUPER_Z[1]:5.1f}%)")
    print(f"  SUPER_W (all)        {w_all:>12,} / {SUPER_W[1]:>12,}  "
          f"({100*w_all/SUPER_W[1]:5.1f}%)")
    print(f"  SUPER_W (excl. PAR)  {w_nopar:>12,} / "
          f"{SUPER_W[1]-PAR_END:>12,}  "
          f"({100*w_nopar/(SUPER_W[1]-PAR_END):5.1f}%)   <- the real test")
    for acc, L in W_UNLOC.items():
        c = merged_len(iv[acc])
        print(f"  {acc:<20} {c:>12,} / {L:>12,}  ({100*c/L:5.1f}%)"
              f"   [repeat array, both sexes - not diagnostic]")

    top = sorted(contigs[SUPER_Z[0]].items(), key=lambda kv: -kv[1])[:3]
    print("  top Z contigs: " + ", ".join(f"{c} ({bp:,} bp)" for c, bp in top))
    return z / SUPER_Z[1], w_nopar / (SUPER_W[1] - PAR_END)


def main():
    res = {hap: report(hap) for hap in ("pri", "alt")}

    z_both = all(z >= 0.50 for z, _ in res.values())
    w_max = max(w for _, w in res.values())

    print("\n" + "=" * 62)
    print(f"Z breadth   pri {res['pri'][0]*100:.1f}%   alt {res['alt'][0]*100:.1f}%"
          f"   -> Z in both haplotypes: {z_both}")
    print(f"W breadth (excl. PAR), best haplotype: {w_max*100:.1f}%")
    if w_max >= 0.40:
        verdict = "FEMALE (ZW) -- W chromosome present"
    elif w_max < 0.15 and z_both:
        verdict = ("MALE (ZZ) -- Z in both haplotypes, no W beyond the "
                   "pseudoautosomal region")
    else:
        verdict = "AMBIGUOUS -- check the read-depth test (03/04)"
    print(f"VERDICT: {verdict}")
    print("=" * 62)


if __name__ == "__main__":
    main()
