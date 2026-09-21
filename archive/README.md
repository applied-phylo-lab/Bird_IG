# archive/

Superseded scripts, kept for reference. Nothing in the pipeline reads them.

| Script | Superseded by | Why |
|--------|---------------|-----|
| `create_summary_features.R` | `data_prep/create_summary_tables_clean.R` | Built the index by walking the tree for `refined_ig_loci/summary.csv` files. Wrote to `summary_features_old.csv`, not the name anything reads, and had no filtering for locus or V gene count. |
| `summarize_features.py` | `data_prep/create_summary_tables_clean.R` | Same job again in Python, walking for `summary.csv` under `refined_ig_loci`. A third route to one table. |

The index is now built once, from Daniel's `ig_contig_list.csv`, by
`create_summary_tables_clean.R` — see README.md. Having three scripts that each
produced a slightly different `summary_features.csv` is exactly the ambiguity
that let a four-month-old index go unnoticed.

## Individually archived scripts

| Script | Why |
|--------|-----|
| `alignment_all_pairwise.py` | All-vs-all LASTZ per species. Wrote `{hap}_{contig}_vs_{hap}_{contig}.txt`, a naming nothing ever read; PatchWorkPlot now produces the pairwise alignments. |
| `gene_trees.py` | Cross-species V gene trees plus an annotation TSV. The pipeline builds per-haplotype trees instead (`workflow/rules/trees.smk`). |
| `shared_inversions.py` | Reciprocal-overlap / Jaccard metrics for shared inversions. Expects `{order}/pairwise_alignments/{hap1}_{hap2}.txt`, a layout that no longer exists, and reads the pre-`{contig}_` `IGH_self.tsv`. Superseded in spirit by `inversions/shared_inversions_counts.py`. |
| `tree_building_pipeline.py` | Ran createFasta → clustalo → iqtree2 inside a multiprocessing.Pool with its own skip logic. Now three rules in `workflow/rules/trees.smk`. |
| `within_species_alignment_mummer.py` | MUMmer pri-vs-alt dot and synteny plots. |
| `avg_inversion_length.py` | One-off: mean inversion length over a patchworkplot folder. |
| `overview_inversions.R` | One-off: read a single jay summary TSV. |
| `BM_analysis_old.R` | Superseded Brownian-motion analysis. |
| `inversion_tree_figure.R` | Had no data loading at all — it opened with `order_nodes <- order_nodes[-c(6,13,15),]` and relied entirely on globals left by `plots/inversion_stats_overview.R`. |
| `repeatmasker/summarize_repeats.py` | Repeat content per locus; filed with the other RepeatMasker scripts. |
| `RSS/inversion_density_RSS.py` | Built `RSS_stats.tsv`. Its only consumer is `plots/RSS_inversion_density.R`, which is exploratory — no manuscript figure reads it. |
| `human_bird_comparison/` | `human_vs_bird_comparison.py` moved in alongside the CSVs it produces, then the whole folder archived. |
