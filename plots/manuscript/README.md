# plots/manuscript/

The six scripts that produce manuscript figures. Everything else in `plots/` is
exploratory.

| Script | Saves |
|--------|-------|
| `Figure1AB.R` | `Figure1AB_resized.svg`, `Figure1_contig_len{,_green}.svg`, `IGL_histograms{,_green}.svg`, `IGH_IGL_contig_strand.svg` |
| `figure_1c.R` | `Figure1C.svg` — avg inversion length \| genes in inversion region \| phylolm scatter |
| `mindir_tree.R` | `mindir_tree_current_weighted.svg` — the butterfly |
| `manhattan_plot.R` | `manhattan_inversions_bHaeMex1{,_IGL}.svg`, `manhattan_covered_fraction_bHaeMex1.svg` |
| `hairpin_identity.R` | `hairpin.svg`, `hairpin_all_windows.svg` |
| `rss_position_oriented.R` | `IGH_IGL_RSS_pos_oriented.svg` |
| `rss_correlation.R` | **nothing** — builds `p_combined_simple`, the left panel of the figure above, and is sourced by it |

14 figures per dataset version; a clean run of all six exits 0 for both v1 and v2.

Plus the PatchWorkPlot dot plots, chiefly `{INPUT_DIR}/patchworkplot/plots_fig1c/`
(house finch, tawny owl, European golden plover, chicken, ostrich).

Several panels are built but deliberately not saved — the two-panel `mindir_tree`
view, `p_locus`/`p_strand` from Figure 1C, and the other RSS panels. They print to
the viewer and can be exported by hand; the pipeline just doesn't write figures
nobody uses.

## Switches

- **`plots/_dataset.R`** — `DATASET <- "v1"` or `"v2"`, which resolves the tree,
  the VGP table and the output directory (`figures/v1/` or `figures/v2/`).
- **`mindir_tree.R`** — `GENE_LIST` (`current` / `legacy`) and `AGGREGATION`
  (`weighted` / `mean` / `pooled` / `min_genes`). **`current` + `weighted` is what
  the pipeline runs**; the rest stay available for comparison and the choice is
  recorded in the output filename.

`weighted` matters because MinDir is a per-contig strand fraction and a contig
carrying one gene scores 1.0 by construction (40 of 904 contigs). Weighting by
gene count stops those outvoting real loci, and cuts the worst legacy-vs-current
discrepancy from 0.33 to 0.07 — which makes the choice of gene list close to
cosmetic.

## Saving figures

Figures are exported **by hand** from the RStudio viewer for the manuscript —
better control over font size and framing. Every script still ends with
`save_fig()` so the set reproduces unattended; treat that output as a correctness
check, not as final artwork.

`save_fig()` defaults to **1472 × 472 px at 96 dpi**, matching the hand-exported
figures. Exceptions set their own size: `Figure1AB_resized.svg` (980 px wide), the
hairpin violins, and the `mindir_tree` butterfly (9 × 12 in).

Two rendering notes, both the result of chasing visual bugs:

- SVGs go through the **cairo** device, the same one RStudio's Export → SVG uses.
  It bakes text into glyph outlines (so text isn't selectable), but `svglite`
  writes a single `font-family` with no fallback, and any viewer lacking that font
  falls back to serif — which made the figures look like Times New Roman.
- Colour bars need `guide_colourbar(raster = FALSE)`. Cairo renders a rasterised
  bar through `<filter>`/`<feImage>`/`<mask>`, which many viewers show as a
  checkerboard.

The butterfly's mirrored axis can't be relabelled through `axis.params$text` —
`ggtreeExtra:::build_axis` only honours that field when the axis has a single
break — so `mirror_axis_labels()` rewrites the tick layer afterwards.

## Gene lists

`rss_correlation.R` reads `gene_list-igl_h3_n7.csv`. The `igl_hN_nM` suffix
records the IGL RSS calling thresholds (heptamer ≤ N, nonamer ≤ M mismatches),
and that choice dominates the figure: IGL RSS detection differs 2.3-fold between
`h3_n7` (32.3%) and `gene_list.csv` (14.0%). IGH is unaffected. `h3_n7` is the
surviving member of the same family as the original `h2_n6`, which no longer
exists on disk.

`gene_list.csv` also carries 25,150 TRA/TRB/TRG/TRD rows, which nothing in this
script filtered — they silently became six loci instead of two. It now filters to
IGH/IGL on load.

`rss_position_oriented.R` reads `gene_list.csv` but filters to IGH/IGL at every
use, so it never had the TCR problem. Its IGL RSS calls come from the stricter
thresholds and so won't match `rss_correlation.R`.

## Cross-script dependency

`rss_position_oriented.R` draws two composite panels using `p_combined_simple`
from `rss_correlation.R`, and skips them with a message if it hasn't been sourced.
The manuscript figure itself doesn't need it.
