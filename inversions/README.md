# inversions/

Everything derived from the self-alignments: inversion statistics, D genes on
inversions, and hairpin/palindrome structure.

| Script | What it does | In the workflow? |
|--------|-------------|------------------|
| `summarize_inversions.py` | Parses each `{Contig}_{Locus}.tsv`, calls inversions (`strand1 != strand2`), and writes `inversion_stats.tsv` (one row per haplotype × contig) and `inversion_details.tsv` (one row per inversion). | yes — `summarize_inversions` |
| `hairpin.py` | For each diagonal inversion, compares identity of the whole alignment against a window at its centre (the putative hairpin tip) and a random window of the same size. Writes `palindromes.tsv`. | yes — `hairpin` |
| `d_genes_on_inversions.py` | Fraction of IGH D genes falling inside an inverted region. Writes `D_inversions.tsv`. | rule exists, **not a default target** — see `{INPUT_DIR}/old_unused/README.md` for why the table was retired |
| `confirm_palindrome_positions.py` | **One-off check**, not part of the pipeline. Re-extracted candidate palindrome sequences and re-ran the alignment to confirm the positions `hairpin.py` reported were real. Kept next to `hairpin.py` because that is what it was checking. | no |
| `shared_inversions_counts.py` | Counts inversions in *pairwise* (cross-haplotype) alignments, binned by length and optionally identity. Replaces the former `shared_inversions_songbirds.py` / `_house_finch.py` / `_within_species.py` trio via `--mode`. Reads PatchWorkPlot's `pairwise_alignments/`. | no — nothing calls it yet |

```bash
snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 20 -- inversions
```

`--min_len` defaults to 250 everywhere; that is the only threshold this project
uses.
