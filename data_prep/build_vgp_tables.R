#!/usr/bin/env Rscript
# Map each haplotype in summary_features.csv to a scientific name, prune the VGP
# species tree to the species we actually have, and write the three tables the
# phylogenetic analyses depend on.
#
# Writes:
#   IGH_VGP_table.tsv  IGH, highest-NumV contig per SPECIES, with LatinName
#                      -- one row per species; drives every phylolm analysis
#   IGH_table.tsv      IGH, highest-NumV contig per HAPLOTYPE (no LatinName)
#   vgp_birds.nwk      VGP tree pruned to species present in the data
#
# Usage:
#   Rscript data_prep/build_vgp_tables.R [-i INPUT_DIR] [--vgp-tree NWK]
#                                        [--vgp-table TSV] [--overrides CSV]
#
# Split out of the interactive overview_features.R: that file interleaved these
# three outputs with ggtree exploration, so the tables could not be rebuilt
# without also running the plots. The plotting half now lives in
# plots/overview_features_plots.R.

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(purrr)
  library(stringr); library(readr); library(ape)
})

INPUT_DIR <- "/local/storage/kav67/clean_birds"
VGP_TREE  <- "/local/storage/kav67/roadies_v1.1.4.nwk"
VGP_TABLE <- "/local/storage/kav67/VGP_phase1_copy062025.tsv"
OVERRIDES <- NULL
# IGH_VGP_table.tsv feeds the phylogenetic analyses, which need exactly one
# haplotype per species -- two haplotypes of the same species cannot both sit on
# one tree tip. The live table has always been primaries only with NumV > 2;
# in overview_features.R the _pri filter had been commented out and the NumV one
# came from silently reading IGH_filtered_table.tsv instead. Both are explicit now.
PRIMARY_ONLY <- TRUE
MIN_NUMV     <- 2
# Suffix on the three output filenames, so an updated dataset can be built
# alongside the frozen manuscript one instead of overwriting it.
# "" -> vgp_birds.nwk ; "_v2" -> vgp_birds_v2.nwk
SUFFIX       <- ""

.args <- commandArgs(trailingOnly = TRUE)
.i <- 1L
while (.i <= length(.args)) {
  .key <- .args[.i]
  if (.key %in% c("-h", "--help")) {
    cat("Usage: build_vgp_tables.R [-i INPUT_DIR] [--vgp-tree NWK]",
        "[--vgp-table TSV] [--overrides CSV] [--min-numv N]",
        "[--all-haplotypes] [--suffix _v2]\n"); quit(status = 0)
  }
  if (.i == length(.args)) stop("Missing value for argument: ", .key)
  .val <- .args[.i + 1L]
  if (.key %in% c("-i", "--input_dir")) INPUT_DIR <- .val
  else if (.key == "--vgp-tree")  VGP_TREE  <- .val
  else if (.key == "--vgp-table") VGP_TABLE <- .val
  else if (.key == "--overrides") OVERRIDES <- .val
  else if (.key == "--min-numv") MIN_NUMV <- as.integer(.val)
  else if (.key == "--suffix") SUFFIX <- .val
  else if (.key == "--all-haplotypes") { PRIMARY_ONLY <- FALSE; .i <- .i - 1L }
  else stop("Unknown argument: ", .key)
  .i <- .i + 2L
}
dir <- sub("/+$", "", INPUT_DIR)

for (f in c(file.path(dir, "summary_features.csv"), VGP_TREE, VGP_TABLE)) {
  if (!file.exists(f)) stop("Input not found: ", f)
}

summary_table <- fread(file.path(dir, "summary_features.csv"))
VGP_tree <- read.tree(VGP_TREE)
VGP_data <- fread(VGP_TABLE)

# ---- 1. relabel the VGP tree tips from assembly accession to scientific name --
# Tips come in as GCA/GCF accessions; both columns map to the same species.
lookup <- VGP_data %>%
  select(GCA = `Accession # for main haplotype`,
         GCF = `RefSeq annotation main haplotype`,
         ScientificName = `Scientific Name`) %>%
  pivot_longer(c(GCA, GCF), names_to = "Type", values_to = "Assembly") %>%
  filter(!is.na(Assembly), Assembly != "") %>%
  mutate(Assembly = str_trim(Assembly)) %>%
  select(Assembly, ScientificName) %>%
  distinct(Assembly, .keep_all = TRUE)

VGP_tree$tip.label <- str_trim(VGP_tree$tip.label)
relabelled <- lookup$ScientificName[match(VGP_tree$tip.label, lookup$Assembly)]
n_unmatched <- sum(is.na(relabelled))
VGP_tree$tip.label <- ifelse(is.na(relabelled), VGP_tree$tip.label, relabelled)
message(sprintf("Tree: relabelled %d of %d tips (%d kept their accession)",
                length(relabelled) - n_unmatched, length(relabelled), n_unmatched))

# ---- 2. resolve a LatinName for every haplotype ------------------------------
# Two independent routes, preferring the assembly accession over the name match.
normalize_species_name <- function(x) {
  if (is.na(x) || x == "") return(NA_character_)
  x <- str_trim(str_replace_all(tolower(x), "[^a-z0-9\\s\\-]", ""))
  words <- str_split(x, "\\s+")[[1]]
  paste(vapply(words, function(w) {
    paste0(str_to_title(str_split(w, "-", simplify = FALSE)[[1]]), collapse = "")
  }, character(1)), collapse = "_")
}

# English Name can be a comma-separated list of common names.
vgp_species_map <- VGP_data %>%
  transmute(AssemblyID = `Assembly ID`, ScientificName = `Scientific Name`,
            EnglishName = `English Name`) %>%
  mutate(EnglishParts = str_split(EnglishName, ",")) %>%
  unnest_longer(EnglishParts) %>%
  mutate(EnglishParts = str_trim(EnglishParts),
         SpeciesNorm  = map_chr(EnglishParts, normalize_species_name)) %>%
  filter(!is.na(SpeciesNorm)) %>%
  group_by(SpeciesNorm) %>% dplyr::slice(1) %>% ungroup() %>%
  select(SpeciesNorm, ScientificName_by_Species = ScientificName)

summary_filled <- summary_table %>%
  mutate(AssemblyID = str_extract(Haplotype, "^[^_]+")) %>%
  left_join(VGP_data %>%
              select(AssemblyID = `Assembly ID`,
                     LatinName_by_Assembly = `Scientific Name`) %>%
              distinct(AssemblyID, .keep_all = TRUE),
            by = "AssemblyID") %>%
  left_join(vgp_species_map, by = c("Species" = "SpeciesNorm")) %>%
  mutate(LatinName = coalesce(LatinName_by_Assembly, ScientificName_by_Species)) %>%
  select(-AssemblyID, -LatinName_by_Assembly, -ScientificName_by_Species)

# ---- 3. manual overrides -----------------------------------------------------
# Species with no VGP accession and no matching English name (mostly the jay
# pangenome). Kept in a file rather than inline so adding one needs no code edit.
DEFAULT_OVERRIDES <- data.table(
  Species = c("Florida_scrub_jay", "Yucatan_jay", "Island_scrub_jay",
              "Woodhouse_scrub_jay"),
  LatinName = c("Aphelocoma coerulescens", "Cyanocorax yucatanicus",
                "Aphelocoma insularis", "Aphelocoma woodhouseii"))
ov <- if (is.null(OVERRIDES)) DEFAULT_OVERRIDES else fread(OVERRIDES)
if (!all(c("Species", "LatinName") %in% names(ov))) {
  stop("Overrides file needs Species and LatinName columns: ", OVERRIDES)
}
hit <- match(summary_filled$Species, ov$Species)
summary_filled$LatinName[!is.na(hit)] <- ov$LatinName[hit[!is.na(hit)]]
message(sprintf("Applied %d name overrides covering %d rows",
                nrow(ov), sum(!is.na(hit))))

unmatched <- sort(unique(summary_filled$Species[is.na(summary_filled$LatinName) |
                                                summary_filled$LatinName == ""]))
if (length(unmatched)) {
  message("No LatinName resolved for ", length(unmatched), " species (dropped):")
  for (sp in unmatched) message("  - ", sp)
} else {
  message("All rows have a LatinName.")
}

summary_table_curated <- summary_filled %>% filter(!is.na(LatinName), LatinName != "")

n_before <- nrow(summary_table_curated)
if (PRIMARY_ONLY) {
  summary_table_curated <- summary_table_curated %>% filter(!str_detect(Haplotype, "_alt$"))
}
summary_table_curated <- summary_table_curated %>% filter(NumV > MIN_NUMV)
message(sprintf("Curated rows: %d -> %d (primary_only=%s, NumV > %d)",
                n_before, nrow(summary_table_curated), PRIMARY_ONLY, MIN_NUMV))

# ---- 4. prune the tree to species we have data for ---------------------------
tips_to_keep <- intersect(VGP_tree$tip.label, summary_table_curated$LatinName)
sub_tree <- drop.tip(VGP_tree, setdiff(VGP_tree$tip.label, tips_to_keep))

locus_counts_wide <- summary_table_curated %>%
  count(LatinName, Locus, name = "NumLoci") %>%
  pivot_wider(names_from = Locus, values_from = NumLoci, values_fill = 0) %>%
  dplyr::rename(label = LatinName)

tips_in_data    <- intersect(sub_tree$tip.label, locus_counts_wide$label)
sub_tree_pruned <- drop.tip(sub_tree, setdiff(sub_tree$tip.label, tips_in_data))

# ---- 5. write the tables -----------------------------------------------------
# One row per species: the highest-NumV IGH contig. This is what the phylolm
# analyses join against, so it must stay one-row-per-species.
summary_table_IGH <- summary_table_curated %>%
  filter(Locus == "IGH") %>%
  group_by(LatinName) %>% slice_max(order_by = NumV, n = 1, with_ties = FALSE) %>%
  ungroup()

# One row per haplotype, no LatinName -- includes species absent from the tree.
summary_table_IGH_all <- summary_table %>%
  filter(Locus == "IGH") %>%
  group_by(Haplotype) %>% slice_max(order_by = NumV, n = 1, with_ties = FALSE) %>%
  ungroup()

out_vgp  <- file.path(dir, sprintf("IGH_VGP_table%s.tsv", SUFFIX))
out_all  <- file.path(dir, sprintf("IGH_table%s.tsv", SUFFIX))
out_tree <- file.path(dir, sprintf("vgp_birds%s.nwk", SUFFIX))

write_tsv(summary_table_IGH,     out_vgp)
write_tsv(summary_table_IGH_all, out_all)
write.tree(sub_tree_pruned,      out_tree)

cat("\nDone. Wrote:\n")
cat(sprintf(" - %s (%d rows, %d species)\n", out_vgp,
            nrow(summary_table_IGH), n_distinct(summary_table_IGH$LatinName)))
cat(sprintf(" - %s (%d rows, %d haplotypes)\n", out_all,
            nrow(summary_table_IGH_all), n_distinct(summary_table_IGH_all$Haplotype)))
cat(sprintf(" - %s (%d tips)\n", out_tree, length(sub_tree_pruned$tip.label)))
