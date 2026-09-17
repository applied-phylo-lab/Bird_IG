# trait_analyses/ — SHELVED

Ecological-trait and pathogen-exposure correlates of IGH structure. **Neither
line produced a usable signal**; both are kept for the record and are expected
to move to an `obsolete/` directory once that exists.

| Script | What it does |
|--------|-------------|
| `fetch_avonet_traits.py` | Joins the species list to AVONET (Tobias et al. 2022) for migration status, hand-wing index, mass, habitat and trophic traits. Writes `{INPUT_DIR}/species_traits_avonet.csv`. |
| `fetch_eid2.R` | Per-species pathogen richness from EID2 (host-pathogen associations, with `papersNo`/`sequencesNo` as sampling-effort proxies). Writes `{INPUT_DIR}/species_pathogens_eid2.csv`. |
| `fetch_malavi.R` | Per-species haemosporidian lineage diversity from MalAvi. Writes `{INPUT_DIR}/species_pathogens_malavi.csv`. |
| `migration_traits.R` | phylolm of NumV and inversion counts against migration and other AVONET traits. |
| `pathogen_exposure.R` | phylolm of NumV and inversion counts against pathogen richness. No signal once sampling effort is controlled for. |

## Why these are not in `data_prep/`

They build species-level *covariates* for specific hypotheses, not the index and
annotation tables the main pipeline depends on. Nothing in stages 2-7 reads their
output.

## One live dependency

`annotation_tables/build_summary_tables.py` reads `species_traits_avonet.csv`, but
**only as a LatinName lookup** (`load_latin_names()`), never for the traits
themselves. So `species_traits_avonet.csv` still needs to exist for the
publication tables; the analyses in this folder do not.

If these move to `obsolete/`, either keep `fetch_avonet_traits.py` reachable or
fold its LatinName mapping into the annotation-table step first.
