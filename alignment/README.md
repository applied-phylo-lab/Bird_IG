# alignment/

LASTZ self-alignment of each extracted IG locus — the first compute stage, and
what everything downstream reads.

| Script | What it does | In the workflow? |
|--------|-------------|------------------|
| `self_alignment_bed.py` | Self-aligns each locus FASTA with LASTZ and writes `{Contig}_{Locus}.tsv`, plus gene positions as `{Contig}_{Locus}.bed` and `{Contig}_{Locus}_strand.bed` (minus strand grey). One script for both loci via `--locus`. | yes — `self_align_bed` |

```bash
snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 20 -- align
```

Standalone it takes `-i INPUT_DIR`, `-s summary`, `--locus`, `-c N`, `--lastz`,
plus `--order/--species/--haplotype/--contig` to run a single unit and `--force`
to regenerate existing output.

`rules/align.smk` also holds `self_align_text`, a second LASTZ pass that keeps the
aligned sequences (`text1`/`text2`) because `hairpin.py` needs them and the
stage-2 alignment does not carry them.
