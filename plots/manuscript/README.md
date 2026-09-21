# plots/manuscript/

The six scripts that produce manuscript figures. Everything else in `plots/` is
exploratory.

| Script | Figure | Key inputs |
|--------|--------|-----------|
| `Figure1AB.R` | Locus length and strand bias, birds vs other vertebrates | `all_species_stats_pruned_12052025.csv` (built outside this repo by `match_names_VGP.R`) |
| `figure_1c.R` | Locus length / V count / strand / inversion distributions + phylolm scatter | `annotation_iroki_strand.tsv` (from `data_prep/build_iroki_annotation.R`), `inversion_stats.tsv`, tree |
| `mindir_tree.R` | MinDir on the phylogeny | gene list (see below), tree, `IGH_VGP_table.tsv` |
| `manhattan_plot.R` | Genome-wide inversion density | window summaries from `manhattan_plot_inversion_coverage/` |
| `hairpin_identity.R` | Hairpin identity at inversion centres | `palindromes.tsv` |
| `rss_position_oriented.R` | RSS positional distributions, biologically oriented | `gene_list.csv`, `IGH_VGP_table.tsv`, per-haplotype `IGHD.csv` |

Plus the PatchWorkPlot dot plots, chiefly
`{INPUT_DIR}/patchworkplot/plots_fig1c/` — built from `config_fig1c.csv`
(5 species: house finch, tawny owl, European golden plover, chicken, ostrich).

## Two switches

Both live at the top of the file that uses them.

- **`plots/_dataset.R`** — `DATASET <- "v1"` or `"v2"`. Resolves `TREE_FILE`,
  `VGP_TABLE` and `FIG_DIR`, so figures land in `figures/v1/` or `figures/v2/`.
  Only `figure_1c.R`, `mindir_tree.R` and `rss_position_oriented.R` actually
  depend on the tree; the other three ignore it apart from the output directory.
- **`mindir_tree.R`'s `GENE_LIST`** — `"legacy"` (`gene_list-igl_h3_n7.csv`,
  56,463 rows, what the manuscript used) or `"current"` (`gene_list.csv`, 82,451
  rows). Output filenames carry the choice.

## Saving figures

Figures are exported **by hand** from the RStudio viewer for the manuscript —
that gives better control over font sizes and framing than `ggsave`. Every script
still ends with `save_fig()` calls so the whole set reproduces unattended; treat
that output as a correctness check, not as final artwork.

`save_fig()` routes SVGs through `svglite` when available, because the default
`svg` device silently drops the Greek letters in the phylolm annotations.

## Cross-script dependency still outstanding

`rss_position_oriented.R` draws two composite panels using `p_combined_simple`,
which is built by `RSS/rss_correlation.R`. Those two panels are skipped with a
message unless that script has been sourced first; the manuscript figure itself
does not need it. **This is why `RSS/` cannot simply be archived.**
