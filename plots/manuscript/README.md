# plots/manuscript/

The six scripts that produce manuscript figures. Everything else in `plots/` is
exploratory.

| Script | Figures it saves | Key inputs |
|--------|------------------|-----------|
| `Figure1AB.R` | `Figure1AB_resized.svg` (735x354pt), `Figure1_contig_len.svg`, `Figure1_contig_len_green.svg`, `IGL_histograms.svg`, `IGL_histograms_green.svg`, `IGH_IGL_contig_strand.svg` | `all_species_stats_pruned_12052025.csv` (built outside this repo by `data_prep/match_names_VGP.R`) |
| `figure_1c.R` | `Figure1C.svg` -- avg inversion length \| genes in inversion region \| phylolm scatter | `annotation_iroki_strand.tsv` (from `data_prep/build_iroki_annotation.R`), `inversion_stats.tsv`, tree |
| `mindir_tree.R` | `mindir_tree_current_weighted.svg` -- the butterfly | gene list, tree, `IGH_VGP_table.tsv` |
| `manhattan_plot.R` | `manhattan_inversions_bHaeMex1.svg`, `..._IGL.svg`, `manhattan_covered_fraction_bHaeMex1.svg` | window summaries from `manhattan_plot_inversion_coverage/` |
| `hairpin_identity.R` | `hairpin.svg` (7x5in), `hairpin_all_windows.svg` (9x5in) | `palindromes.tsv` |
| `rss_position_oriented.R` | `IGH_IGL_RSS_pos_oriented.svg` -- the RSS manuscript figure | `gene_list-igl_h3_n7.csv`, `IGH_VGP_table.tsv`, per-haplotype `IGHD.csv` |
| `rss_correlation.R` | **none** -- builds `p_combined_simple`, the left panel of the figure above, and is sourced by it | `gene_list-igl_h3_n7.csv`, `IGH_VGP_table.tsv`, tree |

14 figures per dataset version. A clean run of all six scripts exits 0 for both
v1 and v2.

Some panels are deliberately built but not saved: the density version of
Figure 1A/B (`figure_1ab`), the two-panel `mindir_tree` view, `p_locus` and
`p_strand` from Figure 1C, and all of `rss_position_oriented.R`'s other RSS
panels. They print to the viewer and can be exported by hand; the pipeline just
does not write figures nobody uses.

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

### Canvas and font

`save_fig()` defaults to **1472 x 472 px at 96 dpi** = 1103 x 354 pt, the size of
the hand-exported figures. Exceptions set their own: `Figure1AB_resized.svg`
(980 px wide), the hairpin violins (7x5 and 9x5 in) and the `mindir_tree`
butterfly (9x12 in, a tall tree figure).

SVGs go through the **cairo** device, the same one RStudio's Export -> SVG uses.
It draws text as glyph outlines, so the font is baked into the file. `svglite`
keeps text as text but writes a single `font-family` with no fallback
(`"Liberation Sans"`), and any viewer without that exact font drops to its default
serif -- which made the figures look like Times New Roman. The trade-off is that
text in the SVG is no longer selectable.

Colour bars use `guide_colourbar(raster = FALSE)`: cairo renders a rasterised bar
through `<filter>`/`<feImage>`/`<mask>`, which many viewers show as a
checkerboard.

## Cross-script dependency still outstanding

`rss_position_oriented.R` draws two composite panels using `p_combined_simple`,
which is built by `plots/manuscript/rss_correlation.R`. Those two panels are skipped with a
message unless that script has been sourced first; the manuscript figure itself
does not need it. **This is why `RSS/` cannot simply be archived.**

## Which gene list the RSS figure uses

`rss_correlation.R` originally read `gene_list-igl_h2_n6.csv`, which no longer
exists. The `igl_hN_nM` suffix records the **IGL RSS calling thresholds**
(heptamer <= N, nonamer <= M mismatches), and that choice dominates the figure:

| gene list | IGH with RSS | IGL with RSS |
|---|---|---|
| `gene_list-igl_h3_n7.csv` | 2303 / 41861 (5.50%) | **4713 / 14602 (32.3%)** |
| `gene_list.csv` | 2341 / 42431 (5.52%) | **2078 / 14870 (14.0%)** |

IGH is unaffected; IGL RSS detection differs 2.3-fold. `h3_n7` is the surviving
member of the same family as the original `h2_n6`, so it is the default;
`GENE_LIST <- "default"` switches to `gene_list.csv` for comparison.

`gene_list.csv` additionally carries 25,150 TRA/TRB/TRG/TRD rows. Nothing
downstream in this script filtered them -- `counts`, the per-locus fits and the
combined plot all key off `unique(counts$Locus)` -- so they silently became six
loci instead of two. The script now filters to IGH/IGL on load.

Note `rss_position_oriented.R` still reads `gene_list.csv`. It filters to IGH/IGL
explicitly at every use, so it never had the TCR problem, but its IGL RSS calls
come from the stricter thresholds and so will not match this figure.
