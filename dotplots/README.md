# dotplots/

PatchWorkPlot triangle dot plots: every haplotype in a group aligned against every
other, with gene positions overlaid.

| Script | What it does | In the workflow? |
|--------|-------------|------------------|
| `make_config_strand.py` | Builds the config PatchWorkPlot consumes — one row per haplotype giving its locus FASTA and strand BED, taking the highest-NumV contig where a haplotype has several. | yes — `dotplot_config` |

```bash
# all orders x both loci
snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 20 -- dotplots
# just the Figure 1C panel (5 hand-picked species, config/config.yaml: fig1c_samples)
snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 8  -- fig1c
```

PatchWorkPlot itself is not on `$PATH`; `config/config.yaml` points at the
checkout (`patchworkplot:`). Its `pairwise_alignments/` output is what a rewrite
of the shared-inversion analysis should read — see
`inversions/shared_inversions_counts.py`.
