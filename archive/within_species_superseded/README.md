# archive/within_species_superseded/

`shared_inversions_no_ref_stats1.py` — an earlier iteration of
`within_species_inversions/shared_inversions_no_ref_stats.py`.

The two take **identical** command lines (`-s/--summary`, `-l/--lastz_dir`,
`-o/--output_prefix`, `-c/--cores`), so there was no flag to merge them behind.
`_stats.py` is a strict superset:

- adds `pick_representative()`, choosing a representative inversion per cluster
- writes a third output, `{prefix}_coords.tsv`, alongside `_presence.tsv` and
  `_stats.tsv`

Nothing `_stats1.py` did is unavailable in `_stats.py`, so it was archived rather
than merged.
