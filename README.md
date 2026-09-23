# Bird Immunoglobulin Locus Analysis

Comparative genomics of the **IGH** and **IGL** loci across VGP bird assemblies:
how inversions and palindromes shape V gene copy number, where RSS sit and which
way they point, and how all of it maps onto the bird phylogeny.

Preprint: https://www.biorxiv.org/content/10.64898/2026.09.05.749481v1

All data lives under `INPUT_DIR`, organised as `{Order}/{Species}/{Haplotype}/`
— see [docs/data_layout.md](docs/data_layout.md) for the full tree.

---

## Quick start

```bash
conda activate snakemake
snakemake -s workflow/Snakefile --configfile config/config.yaml -n      # dry run
snakemake -s workflow/Snakefile --configfile config/config.yaml --profile workflow/profiles/local
```

Individual stages, if you don't want the whole thing:

```bash
snakemake ... -- align      # self-alignments + BED
snakemake ... -- inversions # inversion stats, hairpins
snakemake ... -- paralogs
snakemake ... -- trees
snakemake ... -- dotplots   # expensive; not in the default run
snakemake ... -- tables     # publication tables; not in the default run
snakemake ... -- fig1c      # just the Figure 1C dot plot
```

**Note the `--`.** `--configfile` takes a value, so a bare target after it gets
swallowed as an argument.


---

## The pipeline

![rulegraph](docs/workflow_rulegraph.svg)

Everything is driven by `summary_features.csv`, one row per haplotype × locus ×
contig. It is a pipeline **input**, not a rule output — read at parse time to
build the wildcard lists, so it has to exist before Snakemake starts.

**Stage 1 — index and gene tables** (`data_prep/`, run by hand)
`create_summary_tables_clean.R` builds the index from the contig list,
applying `config/excluded_haplotypes.csv`. `filter_genes.py` cleans the raw gene
calls per haplotype. `build_vgp_tables.R` maps haplotypes to scientific names and
prunes the VGP tree to what's present.

**Stage 2 — self-alignment** (`self_align_bed`, `self_align_text`)
Each locus FASTA is aligned against itself with LASTZ. Inverted alignments
(`strand2 = -`) are the raw signal for everything downstream. Gene positions are
written alongside as BED, once plain and once coloured by strand.

**Stage 3 — inversions** (`summarize_inversions`, `hairpin`)
Parses the alignments into `inversion_stats.tsv` (per haplotype × contig) and
`inversion_details.tsv` (per inversion), at a minimum length of 250 bp. `hairpin`
re-runs LASTZ with sequence in the output to compare identity at the centre of
each diagonal inversion — the putative hairpin tip — against a random control
window, giving `palindromes.tsv`.

A haplotype's locus is sometimes split across contigs, so **one haplotype can
have several rows**. Aggregate accordingly.

`d_genes` is kept but off by default: the upstream `IGHD.csv` files became a
threshold-swept candidate list, so the counts aren't meaningful until a calling
threshold is agreed. The old output is in `{INPUT_DIR}/old_unused/`.

**Stage 4 — paralogs** (`paralogs/`)
Groups genes lying on opposite ends of the same inversion — evidence for
inversion-mediated duplication.

**Stage 5 — gene trees** (`tree_fasta` → `tree_align` → `tree_iqtree`)
A V gene phylogeny per haplotype: extract sequences, align with Clustal Omega,
build tree with IQ-TREE.

**Stages 6–7 — dot plots and tables** (not in the default run)
PatchWorkPlot dot plots per order, and the publication-facing tables under
`summary_tables/`.

Regenerate the graph after changing rules — the targets must be listed, since
`--rulegraph` only graphs what it's asked for:

```bash
snakemake -s workflow/Snakefile --configfile config/config.yaml --rulegraph \
  -- all dotplots tables fig1c d_genes | dot -Tsvg > docs/workflow_rulegraph.svg
```

---

## Figures

The six manuscript scripts are in `plots/manuscript/` — see the README there for
what each one saves and the switches it reads. Everything else in `plots/` is
exploratory and meant to be sourced in RStudio.

Figures exist in two versions, with the second version having added birds after using a better database:

| | tree | table | species in phylolm |
|---|---|---|---|
| **v1** (manuscript) | `vgp_birds.nwk` — 122 tips | `IGH_VGP_table.tsv` — 128 rows | 116 |
| **v2** (updated) | `vgp_birds_v2.nwk` — 137 tips | `IGH_VGP_table_v2.tsv` — 127 rows | 123 |

v1: Reference copies and the full explanation are in
`{INPUT_DIR}/frozen_v1_manuscript/`. 
v2 is a strict superset and is what `data_prep/build_vgp_tables.R` produces. The effect on the headline fit is small
(β 0.511 → 0.515, R² 0.641 → 0.645).

Switch with `DATASET <- "v1"` / `"v2"` at the top of `plots/_dataset.R`; output
lands in `figures/v1/` or `figures/v2/`.

---

## Repository layout

```
workflow/  config/     Snakemake pipeline (stages 2-7) + local/slurm profiles
data_prep/             Index tables, VGP tree/tables, gene filtering
alignment/  inversions/  paralogs/  dotplots/   Pipeline scripts, by stage
annotation_tables/     Publication-facing tables
plots/                 Exploratory R figures
plots/manuscript/      The six manuscript figure scripts
tree_analyses/         Tree-distance analyses
within_species_inversions/          Within-species comparison
manhattan_plot_inversion_coverage/  Window summaries for the manhattan figure
sex_check/             Self-contained ZW read-depth check
daniel_bird_scripts/   Collaborator scripts (read-only)
figures/v1/  figures/v2/            Rendered figures per dataset version
archive/               Superseded and shelved work, each folder with a README
```

---


## Dependencies

**Conda envs:** `snakemake` (Python scripts + `lastz`, `clustalo`, `iqtree2`),
`bird_ig` (R).

**Python:** pandas, numpy, biopython
**R:** tidyverse, data.table, ape, ggtree, ggtreeExtra, patchwork, phylolm, phytools, ggrepel, viridis
**External:** lastz, clustalo, iqtree2
