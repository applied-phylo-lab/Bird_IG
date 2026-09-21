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
| `rss_correlation.R` | **The RSS manuscript figure**: genes with RSS vs total genes (phylolm per locus) beside single- vs multiple-productive-RSS positional density, IGH above IGL | `gene_list.csv`, `IGH_VGP_table.tsv`, tree |
| `rss_position_oriented.R` | RSS positional distributions, biologically oriented. **Saves nothing** — kept as the analysis behind the figure above, and because it is useful interactively | `gene_list.csv`, `IGH_VGP_table.tsv`, per-haplotype `IGHD.csv` |

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
  rows). **Now set to `current`.**
- **`mindir_tree.R`'s `AGGREGATION`** — how per-contig MinDir values are combined
  when a locus spans several contigs. **Now set to `weighted`.**

**`current` + `weighted` is what the pipeline runs.** The other combinations stay
available for comparison — switch, re-run, and the output filename records the
choice — but they are not produced by a normal run.

Output filenames carry both switches, e.g. `mindir_tree_current_weighted.svg`.

Only the butterfly is saved. The two-panel IGH/IGL view (`p`) is still built and
prints to the viewer, but writing it every run cluttered `figures/` with a panel
that is not in the manuscript.

### Why `weighted`

MinDir is a per-contig strand fraction, and a contig carrying one gene scores 1.0
by construction — 40 of 904 contigs (4.4%) in `gene_list.csv` are exactly that.
An unweighted mean lets them count as much as a 200-gene locus. That is what moved
the Red-billed Tropicbird's IGH from 0.50 to 0.83 between the two gene lists: two
scrap contigs of 2 and 1 genes outvoting its real 18-gene locus 2:1.

Legacy vs current gene list, by rule:

| rule | mean(legacy) | mean(current) | n changed | mean abs diff | max abs diff |
|---|---|---|---|---|---|
| `mean` | 0.7472 | 0.7467 | 20 | 0.0655 | 0.3333 |
| **`weighted`** | 0.7467 | 0.7466 | 20 | **0.0116** | **0.0714** |
| `pooled` | 0.7451 | 0.7441 | 24 | 0.0195 | 0.0779 |
| `min_genes` (>=5) | 0.7427 | 0.7421 | 8 | 0.0355 | 0.1006 |

`weighted` keeps every gene while cutting the worst legacy/current discrepancy
from 0.333 to 0.071, so the choice of gene list is close to cosmetic. `pooled` is
arguably more principled but shifts more species; `min_genes` discards 470 of
2670 contigs on an arbitrary threshold.

### Butterfly axis

The butterfly mirrors IGH onto negative x. `axis.params$text` cannot relabel it —
`ggtreeExtra:::build_axis` only honours that field when the axis has a single
break — so `mirror_axis_labels()` rewrites the tick layer afterwards, using the
magnitude as the label and the signed value as the position. Tick size is
`AXIS_TEXT_SIZE` at the top of the script.

## Saving figures

Figures are exported **by hand** from the RStudio viewer for the manuscript —
that gives better control over font sizes and framing than `ggsave`. Every script
still ends with `save_fig()` calls so the whole set reproduces unattended; treat
that output as a correctness check, not as final artwork.

`save_fig()` routes SVGs through `svglite` when available, because the default
`svg` device silently drops the Greek letters in the phylolm annotations.

## Cross-script dependency still outstanding

`rss_position_oriented.R` draws two composite panels using `p_combined_simple`,
which is built by `plots/manuscript/rss_correlation.R`. Those two panels are skipped with a
message unless that script has been sourced first; the manuscript figure itself
does not need it. **This is why `RSS/` cannot simply be archived.**
