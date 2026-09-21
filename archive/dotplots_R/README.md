# archive/dotplots_R/

Dot plots drawn directly in R from a LASTZ self-alignment TSV, rather than
through PatchWorkPlot. Superseded by the `dotplots` / `fig1c` workflow targets.

| Script | What it did |
|--------|-------------|
| `dotplot_1000.R` | Dot plot of one haplotype's self-alignment, filtered to inversions of at least 1000 bp. Reads `IGH_self.tsv` under `Bird_data/` or `within_species/`. |
| `plot_diag_inversions.R` | Plots diagonal (palindromic) inversions specifically, from a jay pairwise alignment. |

Both point at data roots that no longer exist (`Bird_data/`, `within_species/`,
`jays/`), so neither runs as written.
