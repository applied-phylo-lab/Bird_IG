# data_prep/

Stage 1: the index and lookup tables everything else is built from. These run by
hand, before the Snakemake workflow — `summary_features.csv` is read at parse time
to build the wildcard lists, so it must exist first.

| Script | Writes | Notes |
|--------|--------|-------|
| `create_summary_tables_clean.R` | `summary_features.csv`, `IGH_filtered_table.tsv`, `filtered_table.tsv` | The workflow's index, from Daniel's `ig_contig_list.csv`. Filters to IGH/IGL (that file also carries TRA/TRB/TRD/TRG) and applies `config/excluded_haplotypes.csv`. |
| `build_vgp_tables.R` | `IGH_VGP_table.tsv`, `IGH_table.tsv`, `vgp_birds.nwk` | Resolves each haplotype to a scientific name, prunes the VGP tree to the species present. `--suffix _v2` builds the updated dataset alongside the frozen v1. |
| `build_iroki_annotation.R` | `iroki_annotation.tsv`, `annotation_iroki_strand.tsv` | Bar heights for the iroki tree annotation; the six-bar version is what `figure_1c.R` reads. |
| `filter_genes.py` | `combined_genes_{IGH,IGL}_clean.txt` per haplotype | Length and composition filtering. Takes `--order/--species/--haplotype` so the workflow can drive one haplotype at a time. |
| `match_names_VGP.R` | `all_species_stats_pruned_12052025.csv` | Builds `Figure1AB.R`'s input. Operates on `/local/storage/kav67/IG_annotation_VGP2026/`, a different project directory. |

Typical order:

```bash
conda activate bird_ig
Rscript data_prep/create_summary_tables_clean.R -i $INPUT_DIR \
    --exclude config/excluded_haplotypes.csv
Rscript data_prep/build_vgp_tables.R      -i $INPUT_DIR              # v1
Rscript data_prep/build_vgp_tables.R      -i $INPUT_DIR --suffix _v2 # v2
Rscript data_prep/build_iroki_annotation.R -i $INPUT_DIR
```

`build_vgp_tables.R` defaults to primaries only with `NumV > 2`: the tree needs
exactly one haplotype per species. Both filters are flags
(`--all-haplotypes`, `--min-numv`) — in the old `overview_features.R` the first
had been commented out and the second came from silently reading a different
input, which is how the two dataset versions drifted apart.
