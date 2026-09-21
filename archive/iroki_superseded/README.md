# archive/iroki_superseded/

Replaced by `data_prep/build_iroki_annotation.R`, which produces both iroki
annotation tables in one headless script.

| Script | Why it was replaced |
|--------|--------------------|
| `create_iroki_annotation.R` | Wrote `iroki_annotation.tsv`, but needed `species_stats_pruned`, `tree_pruned` and `dir` to already exist in the session — globals left behind by `plots/inversion_stats_overview.R`. Also mixed in interactive phylolm exploration. |
| `iroki_annnotation_strand.R` | Wrote `annotation_iroki_strand.tsv`, needing `vgp_bird_table` and `annotation_iroki` from the script above, plus `bird_tree_saved`, which was **never defined anywhere in the repo**. That is why the file had gone missing and `plots/manuscript/figure_1c.R` died on its first line. |
