# Bird_IG — project context

Comparative genomics of bird immunoglobulin loci (**IGH** and **IGL**) across
VGP assemblies. Research questions: how inversions/palindromes shape V gene
copy number and diversity, RSS presence and orientation, D gene organisation,
and how these traits map onto the bird phylogeny. Output is a manuscript, so
most endpoints are figures + statistics, not software.

`README.md` is the authoritative pipeline documentation (step-by-step, per-script,
with arguments and the full data directory layout). **Read it before answering
questions about how a script fits into the workflow** — don't duplicate it here.

## Data lives outside the repo

Nothing but code is versioned. All data is under `INPUT_DIR`, normally
`/local/storage/kav67/clean_birds/`, organised as `{Order}/{Species}/{Haplotype}/`
(haplotypes are `*_pri` / `*_alt`). Key master tables at the top level:

- `summary_features.csv` — one row per haplotype × locus (Order, Species, Haplotype, Locus, Contig, NumV)
- `IGH_filtered_table.tsv`, `IGH_VGP_table.tsv` (adds `LatinName` for tree matching)
- `gene_list.csv` — all V genes with RSS annotations
- `inversion_stats.tsv`, `inversion_details.tsv`, `D_inversions.tsv`, `palindromes.tsv`
- `vgp_birds.nwk` — VGP species tree used for all phylogenetic corrections
- `species_traits_avonet.csv` — ecological traits (built by `data_prep/fetch_avonet_traits.py`)

Raw assemblies: `/local/storage/dhardesty/assemblies/`.

## Conventions

- **Python** = pipeline/compute (LASTZ self-alignments, inversion detection, tree
  building, RSS extraction). `argparse` with the shared flags `-i INPUT_DIR`,
  `-s summary_features.csv`, `-c cores`; parallelism via `multiprocessing`.
- **R** (`plots/`, `RSS/*.R`, `tree_analyses/*.R`) = statistics and figures, sourced
  interactively in RStudio. `INPUT_DIR` is hardcoded near the top — edit it there.
  New scripts should follow the existing style: `data.table`/`dplyr` + `ggplot2`,
  `ggsave` SVG into `figures/`.
- **Any cross-species comparison must be phylogenetically corrected** (`phylolm`
  with `vgp_birds.nwk`); traits are heavily clade-clustered, so raw correlations
  across species are not usable. See `plots/phylolm_tree.R`, `plots/migration_traits.R`.
- External tools assumed on `$PATH`: `lastz`, `clustalo`, `iqtree2`, `mummer`,
  `RepeatMasker`, `minimap2`/`samtools` (sex check).

## Conda environments

- `snakemake` — Python scripts
- `bird_ig` — R scripts (has `phylolm`, `ape`, `ggtree`)

Activate explicitly before running anything; the base env has neither.

## Directory notes

- `plots/` — one R script per figure/analysis; the manuscript figures are
  `Figure1AB.R`, `figure_1c.R`.
- `daniel_bird_scripts/` — collaborator (Daniel Hardesty) scripts for D gene
  search, RSS extraction, and contig evaluation. **Read-only** — never edit these.
- `sex_check/` — one-off numbered pipeline (ZW read-depth check for a given bird);
  has its own README.
- `human_bird_comparison/`, `jay_tests/`, `min_max_unit/`, `kmer_analysis/` —
  self-contained side analyses.
- `.ipynb_checkpoints/` — junk, ignore.

## Working preferences

- Long-running jobs (LASTZ, IQ-TREE over all haplotypes) are expensive — check for
  existing output before re-running; most pipeline scripts already skip completed work.
- **Never create a near-duplicate script.** If a variant of an existing analysis is
  needed, add a flag/argument to the existing script instead. The several
  `shared_inversions*.py` variants are the mistake to avoid repeating.
- When adding a script, add a row to the matching table in `README.md`.
