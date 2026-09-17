# annotation_tables/

Scripts that produce the final, publication-facing annotation tables — the ones
that get shared or submitted, as opposed to the intermediate per-haplotype files
the analysis pipeline works from.

| Script | What it does |
|--------|-------------|
| `assign_haplotype_source.py` | Classifies each haplotype's provenance (VGP / CCGP / pangenome project / unpublished) and haplotype type (Maternal/Paternal, Hap1/Hap2, Primary/Alternate, ...) by querying the NCBI Datasets API for the assembly's BioProject lineage, cross-referenced against the VGP master sheet. Writes `{INPUT_DIR}/haplotype_sources.csv` and caches raw NCBI responses in `ncbi_bioproject_cache.json`. |
| `build_summary_tables.py` | Combines `summary_features.csv`, `gene_list.csv`, `inversion_stats.tsv`, `D_inversions.tsv`, `palindromes.tsv`, `bird_d_genes.csv` and `haplotype_sources.csv` into `{INPUT_DIR}/summary_tables/`: `main_table`, `species_summary`, `v_gene_table`, `d_gene_table`, each also in a `_published` variant, plus a generated README documenting every column. |

Run in order — `build_summary_tables.py` needs `haplotype_sources.csv`:

```bash
conda activate snakemake
python annotation_tables/assign_haplotype_source.py $INPUT_DIR   # hits the NCBI API
python annotation_tables/build_summary_tables.py   $INPUT_DIR
```

## Inputs from elsewhere

- `species_traits_avonet.csv`, from `trait_analyses/fetch_avonet_traits.py` — used
  **only** as a LatinName lookup, not for the traits.
- `IGH_VGP_table.tsv`, from `data_prep/overview_features.R` — LatinName lookup.
- `gene_list.csv` and `bird_d_genes.csv`, from Daniel's upstream pipeline
  (`daniel_bird_scripts/`, read-only).

## Not yet in the Snakemake workflow

These are stage 7 and still run by hand. Wiring them up mainly means turning the
implicit `path()` reads in `build_summary_tables.py` into declared `input:` so the
DAG can tell when a table is stale — which is how the three-month-old
`D_inversions.tsv` went unnoticed.
