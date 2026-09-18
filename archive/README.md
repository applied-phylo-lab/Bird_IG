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
