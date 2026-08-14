# Sex check for bAgePho2 (red-winged blackbird, Lovette-BB-20993)

Birds are ZW: **males ZZ, females ZW**. So "is this bird a male?" becomes
"does its genome contain a W chromosome, and is its Z diploid?".

Two independent tests, both run off data already on disk.

## Test 1 — W presence / Z ploidy

`GCA_051311805.1_bAgePho1.hap1` (already in `other_assemblies/`) is a
**conspecific female** red-winged blackbird, chromosome level, carrying
`SUPER_Z` (CM119746.1, 96.0 Mb) and `SUPER_W` (CM119744.1, 24.6 Mb) plus two
unlocalized W scaffolds. Same species, so no cross-species divergence to work
around.

```bash
bash 01_map_to_ZW.sh
python3 02_summarize_ZW.py
```

Expected:

| | Z breadth | W breadth (excl. PAR) |
|---|---|---|
| female (ZW) | high in one haplotype only | most of chrW |
| male (ZZ) | high in **both** haplotypes | ~0 |

### Three traps, all of which produce a false "female"

1. **Do not restrict the reference to Z+W.** The first pass here did, and every
   autosomal repeat took its best hit on Z or W because there was nowhere else
   to go — manufacturing 372 Mb of "W-linked" sequence against a 47 Mb
   chromosome. Map against the whole female genome so repeats go home.
   (That pass is archived under `ZWonly_firstpass/` as a cautionary tale.)

2. **The PAR.** The pseudoautosomal region occupies roughly the first 2 Mb of
   SUPER_W and is present in both sexes, so "bp aligned to W" is never zero in
   a male. `02_summarize_ZW.py` reports W breadth both with and without it; the
   PAR-excluded number is the one to read.

3. **The unlocalized W scaffolds.** `W_unloc_1`/`W_unloc_2` (22.7 Mb combined)
   are repeat arrays the female reference binned onto W, but they exist in both
   sexes. The tell is that they appear at similar coverage in *both* haplotypes
   of bAgePho2 — diploid, therefore not W. Reported but excluded from the call.

Identity filtering does **not** rescue any of this: without `-c`, minimap2 does
no base-level alignment and PAF column 10 is a crude estimate, so per-alignment
identity from these files is not trustworthy.

## Test 2 — HiFi read depth (reference-free confirmation)

This bird's own HiFi reads are in `closeread/hifi_bam/bAgePho2/`. A subsample is
mapped back to the **primary assembly only**, so autosomal contigs collect reads
from both homologues:

- **female**: Z and W are single copy → ~0.5x the autosomal depth
- **male**: Z is diploid → ~1.0x, and nothing sits at half depth

```bash
bash 03_read_depth.sh          # FRAC=0.05 by default (~5x)
python3 04_depth_verdict.py
```

This one needs no reference and is immune to all three traps above — but it has
a trap of its own:

4. **Retained haplotigs.** hifiasm leaves duplicate haplotigs in the primary,
   and where one exists the reads split between the two near-identical copies,
   so *both* sit at ~0.5x. A per-contig mean therefore drags a perfectly diploid
   Z toward 0.5x and fakes a female — which is exactly what happens here, since
   `ptg000020l` duplicates ~30 Mb of the Z in `ptg000019l`. For the same reason
   a blanket "which contigs are at half depth" sweep is useless: haplotigs and
   repeat arrays fill that band. `04_depth_verdict.py` therefore reads depth in
   5 Mb windows and takes the verdict from the *unduplicated* part of the Z.

   Also take the autosomal baseline from the largest contigs only. A median over
   all contigs >1 Mb is dragged down by small haplotigs (5.70x here vs the true
   7.39x), which inflates every relative depth by ~30%.

## Result

Both tests agree: **male (ZZ)**.

Test 1 — alignment to the conspecific female reference:

| | Z breadth | W breadth (excl. PAR) |
|---|---|---|
| pri | 95.9% (ptg000019l, 96.1 Mb ≈ full Z) | 0.4% |
| alt | 73.1% (atg000022l + atg000033l) | 3.5% |

Z is present at near-full length in both haplotypes; the real W is empty apart
from the PAR. The alt's residual W hits are concentrated in 0–2 Mb of SUPER_W,
which is exactly where the PAR is.

Test 2 — HiFi depth, autosomal 2n baseline 7.39x:

| contig | windows (relative depth) |
|---|---|
| ptg000019l (Z, 96 Mb) | 1.02 … 0.86 over 0–60 Mb, then 0.56–0.44 |
| ptg000020l (Z haplotig, 19.7 Mb) | 0.49 0.47 0.46 0.42 |

The Z's unduplicated portion sits at **0.99x autosomal** — fully diploid. The
half-depth stretch is precisely where `ptg000020l` duplicates it, and the two
copies sum back to ~1.0x. In a female the Z would read ~0.5x along its whole
length and a ~20 Mb W would sit alongside it at ~0.5x. No such contig exists.

Useful by-products for the IG work:
- `ptg000019l` is the Z chromosome; `ptg000020l` is a retained Z haplotig.
- `ptg000032l` (11 Mb) is a repeat array matching the reference's unlocalized W
  scaffolds — a W-based analysis would misassign it.

## Notes

- Outputs go to `/local/storage/kav67/RedwingedBlackbird_new/sex_check/`.
- Both scripts skip work already done, so they are safe to re-run.
- Field ID is the obvious external check: adult male red-winged blackbirds are
  black with red/yellow epaulets, females streaky brown. If the genomic call and
  the field call disagree, suspect a sample swap rather than an unusual bird.
