# archive/within_species_patchworkplot/

Superseded by the `dotplots` and `fig1c` targets in `workflow/rules/dotplots.smk`.

| Script | What it did |
|--------|-------------|
| `patchworkplot_create_config.py` | Wrote a per-species `config.csv` for PatchWorkPlot from a TSV of haplotypes. `dotplots/make_config_strand.py` plus the `dotplot_config` rule do this now, resolving each haplotype's highest-NumV contig through the index. |
| `patchworkplot_run_all.py` | Looped over species calling PatchWorkPlot with a hardcoded path to the checkout. The `dotplot` rule does this, with the path in `config/config.yaml`. |

The workflow version also declares `pairwise_alignments/` as an output, so those
alignments can be depended on rather than rediscovered.
