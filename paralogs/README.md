# paralogs/

Inversion-mediated gene duplication: genes sitting on opposite arms of the same
inversion are candidate paralogs.

| Script | What it does | In the workflow? |
|--------|-------------|------------------|
| `find_all_inversions_paralogs.py` | Groups genes into paralog clusters using **all** inversions at or above `--min_len` (default 250). Union-find over genes whose mirrored coordinates across an inversion land on another gene. Writes `IGH_paralogs_all.tsv`. | yes — `paralogs_all` |
| `find_inversion_paralogs.py` | Older variant using **only diagonal** (palindromic) inversions, default `--min_len 1000`. Writes `IGH_paralogs_diag.tsv`. | yes — `paralogs_diagonal` |
| `gene_tree_paralogs.R` | Per-haplotype V gene tree coloured by paralog group membership. | no |
| `paralog_fraction_plot.R` | Fraction of V genes that fall into a paralog group, across species. | no |

```bash
snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 8 -- paralogs
```

Both Python scripts take `-i INPUT_DIR`, `-s summary`, `-o OUTPUT.tsv`, `-c N`,
`--locus {IGH,IGL}` and `--min_len N`, and read the per-contig
`{Contig}_{Locus}.tsv` / `.bed` files written by `self_alignment_bed.py`.

The two R scripts still read pre-`clean_birds` paths and have not been repointed.
