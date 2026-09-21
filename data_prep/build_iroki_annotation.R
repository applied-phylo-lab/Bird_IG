#!/usr/bin/env Rscript
# Build the two iroki annotation tables that plots/manuscript/figure_1c.R reads.
#
# Writes:
#   iroki_annotation.tsv        5 bars: locus length, NumV, avg inversion length,
#                               avg inversion count, fraction of genes on inversions
#   annotation_iroki_strand.tsv the same with a 6th bar inserted in position 3:
#                               percentage of IGH V genes on the + strand
#
# Replaces create_iroki_annotation.R + data_prep/iroki_annnotation_strand.R, which
# between them needed five objects that only existed if you had already sourced
# inversion_stats_overview.R in the same session (species_stats_pruned,
# tree_pruned, dir, vgp_bird_table, annotation_iroki) plus one -- bird_tree_saved --
# that was never defined anywhere in the repo. That is why
# annotation_iroki_strand.tsv had gone missing and figure_1c.R died on line 1.
#
# Usage:
#   Rscript data_prep/build_iroki_annotation.R [-i INPUT_DIR] [--suffix _v2]

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(readr); library(ape)
})

INPUT_DIR <- "/local/storage/kav67/clean_birds"
SUFFIX    <- ""

.args <- commandArgs(trailingOnly = TRUE)
.i <- 1L
while (.i <= length(.args)) {
  .key <- .args[.i]
  if (.key %in% c("-h", "--help")) {
    cat("Usage: build_iroki_annotation.R [-i INPUT_DIR] [--suffix _v2]\n"); quit(status = 0)
  }
  if (.i == length(.args)) stop("Missing value for argument: ", .key)
  .val <- .args[.i + 1L]
  if (.key %in% c("-i", "--input_dir")) INPUT_DIR <- .val
  else if (.key == "--suffix") SUFFIX <- .val
  else stop("Unknown argument: ", .key)
  .i <- .i + 2L
}
dir <- sub("/+$", "", INPUT_DIR)

tree_file  <- file.path(dir, sprintf("vgp_birds%s.nwk", SUFFIX))
vgp_file   <- file.path(dir, sprintf("IGH_VGP_table%s.tsv", SUFFIX))
inv_file   <- file.path(dir, "inversion_stats.tsv")
for (f in c(tree_file, vgp_file, inv_file)) {
  if (!file.exists(f)) stop("Input not found: ", f)
}

tree      <- read.tree(tree_file)
vgp_table <- fread(vgp_file)
inv       <- fread(inv_file)

# ---- per-species inversion statistics ---------------------------------------
# Averaged across contigs and haplotypes, reproducing what inversion_stats_overview.R
# computed and handed to create_iroki_annotation.R as `species_stats_pruned`. Note
# this averages fractions rather than recomputing them from summed numerators and
# denominators; kept as-is so the figure stays comparable to the published one.
species_stats <- merge(inv[minlen == 250],
                       unique(vgp_table[, .(Haplotype, LatinName)]),
                       by.x = "haplotype", by.y = "Haplotype", all.x = TRUE) %>%
  filter(!is.na(LatinName), !grepl("MiscBirds", LatinName)) %>%
  group_by(LatinName) %>%
  summarize(avg_inv_len        = mean(avg_inv_len, na.rm = TRUE),
            avg_num_inversions = mean(num_inversions, na.rm = TRUE),
            inv_cov_frac       = mean(inv_cov_len / total_seq_length, na.rm = TRUE),
            frac_genes_on_inv  = mean(frac_genes_on_inv, na.rm = TRUE),
            total_seq_length   = mean(total_seq_length, na.rm = TRUE),
            .groups = "drop")

tips <- intersect(tree$tip.label, gsub(" ", "_", species_stats$LatinName))
species_stats <- species_stats %>% filter(gsub(" ", "_", LatinName) %in% tips)
message(sprintf("Species with inversion stats and a tree tip: %d", nrow(species_stats)))

# ---- one row per species, with locus length ---------------------------------
one_per_species <- vgp_table %>%
  group_by(LatinName) %>% slice_max(NumV, n = 1, with_ties = FALSE) %>% ungroup() %>%
  inner_join(species_stats, by = "LatinName") %>%
  mutate(Length = total_seq_length)

# Two contigs whose recorded length is wrong in the source table; carried over
# from create_iroki_annotation.R, where they were patched by hand.
LENGTH_FIXES <- c("JAJHSZ010000035.1" = 545918, "CM106230.1" = 196575)
for (ctg in names(LENGTH_FIXES)) {
  if (any(one_per_species$Contig == ctg)) {
    one_per_species$Length[one_per_species$Contig == ctg] <- LENGTH_FIXES[[ctg]]
  }
}

# ---- five-bar table ---------------------------------------------------------
# iroki draws bars leftwards, so every height is negated.
annotation_iroki <- one_per_species %>%
  transmute(name        = gsub(" ", "_", LatinName),
            bar1_height = -Length,
            bar2_height = -NumV,
            bar3_height = -avg_inv_len,
            bar4_height = -avg_num_inversions,
            bar5_height = -frac_genes_on_inv,
            bar1_color  = "black",
            bar2_color  = "darkgrey",
            bar3_color  = "#003049",
            bar4_color  = "#669bbc",
            bar5_color  = "#90e0ef")

out1 <- file.path(dir, sprintf("iroki_annotation%s.tsv", SUFFIX))
write_tsv(annotation_iroki, out1)

# ---- six-bar table, with + strand percentage as bar 3 ------------------------
strand_pct <- vapply(seq_len(nrow(one_per_species)), function(i) {
  r <- one_per_species[i, ]
  f <- file.path(dir, r$Order, r$Species, r$Haplotype, "combined_genes_IGH_clean.txt")
  if (!file.exists(f)) return(NA_real_)
  g <- suppressWarnings(read_tsv(f, show_col_types = FALSE))
  g <- g[g$Contig == r$Contig, ]
  if (nrow(g) == 0) return(NA_real_)
  sum(g$Strand == "+") / nrow(g) * 100
}, numeric(1))

n_missing <- sum(is.na(strand_pct))
if (n_missing) message(sprintf("No strand data for %d species", n_missing))

annotation_strand <- annotation_iroki %>%
  dplyr::rename(bar4_height = bar3_height, bar5_height = bar4_height,
                bar6_height = bar5_height,
                bar4_color  = bar3_color,  bar5_color  = bar4_color,
                bar6_color  = bar5_color) %>%
  mutate(bar2_color  = "grey",
         bar3_height = -strand_pct,
         bar3_color  = "darkgrey") %>%
  select(name, bar1_height, bar2_height, bar3_height, bar4_height, bar5_height,
         bar6_height, bar1_color, bar2_color, bar3_color, bar4_color, bar5_color,
         bar6_color)

out2 <- file.path(dir, sprintf("annotation_iroki_strand%s.tsv", SUFFIX))
write_tsv(annotation_strand, out2)

cat("\nDone. Wrote:\n")
cat(sprintf(" - %s (%d species, 5 bars)\n", out1, nrow(annotation_iroki)))
cat(sprintf(" - %s (%d species, 6 bars incl. %% + strand)\n", out2, nrow(annotation_strand)))
