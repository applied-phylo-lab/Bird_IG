#!/usr/bin/env Rscript
# Build the index tables that drive the rest of the pipeline, from Daniel's
# ig_contig_list.csv (one row per haplotype x locus x contig).
#
# Writes:
#   summary_features.csv     every locus, every contig  -- the workflow's index
#   IGH_filtered_table.tsv   IGH only, NumV > min_numv  -- inversion analyses
#   filtered_table.tsv       every locus, tab separated
#
# Usage:
#   Rscript data_prep/create_summary_tables_clean.R \
#       [-i INPUT_DIR] [--input CSV] [--min-numv N]
#
# Defaults keep the previous hardcoded behaviour, so the file can still be
# sourced interactively in RStudio with no arguments.

suppressPackageStartupMessages(library(data.table))

# ---- argument parsing (base R, so no extra package is needed) --------------
INPUT_DIR <- "/local/storage/kav67/clean_birds"
INPUT_CSV <- NULL
MIN_NUMV  <- 2
# ig_contig_list.csv also carries TRA/TRB/TRD/TRG, which are leftover annotation
# from the shared upstream pipeline and out of scope here (see CLAUDE.md).
LOCI      <- c("IGH", "IGL")
# Haplotypes deliberately kept out of the cross-species analysis, with reasons.
EXCLUDE   <- NULL

# Parsed at top level on purpose: wrapping this in local() and using <<- for the
# counter assigns to the global i, never the loop's own, and spins forever.
.args <- commandArgs(trailingOnly = TRUE)
.i <- 1L
while (.i <= length(.args)) {
  .key <- .args[.i]
  if (.key %in% c("-h", "--help")) {
    cat("Usage: create_summary_tables_clean.R [-i INPUT_DIR]",
        "[--input CSV] [--min-numv N] [--loci IGH,IGL]",
        "[--exclude config/excluded_haplotypes.csv]\n")
    quit(status = 0)
  }
  if (.i == length(.args)) stop("Missing value for argument: ", .key)
  .val <- .args[.i + 1L]
  if (.key %in% c("-i", "--input_dir")) {
    INPUT_DIR <- .val
  } else if (.key == "--input") {
    INPUT_CSV <- .val
  } else if (.key == "--min-numv") {
    MIN_NUMV <- as.integer(.val)
  } else if (.key == "--loci") {
    LOCI <- strsplit(.val, ",", fixed = TRUE)[[1]]
  } else if (.key == "--exclude") {
    EXCLUDE <- .val
  } else {
    stop("Unknown argument: ", .key)
  }
  .i <- .i + 2L
}
rm(.args, .i, .key, .val)

dir <- sub("/+$", "", INPUT_DIR)
input_csv      <- if (is.null(INPUT_CSV)) file.path(dir, "ig_contig_list.csv") else INPUT_CSV
output_all     <- file.path(dir, "summary_features.csv")
output_igh     <- file.path(dir, "IGH_filtered_table.tsv")
output_all_tsv <- file.path(dir, "filtered_table.tsv")

if (!file.exists(input_csv)) stop("Contig list not found: ", input_csv)

dt <- fread(input_csv)

# Source looks like "#Order/Species/Haplotype"
dt[, Source := sub("^#", "", Source)]
dt[, c("Order", "Species", "Haplotype") := tstrsplit(Source, "/", fixed = TRUE)]
dt[, NumV := `Number of Genes (before filtering)`]

formatted <- dt[Locus %in% LOCI, .(Order, Species, Haplotype, Locus, Contig, NumV)]

# Drop deliberately excluded haplotypes (see config/excluded_haplotypes.csv).
# Reported individually so a silent drop can never go unnoticed.
if (!is.null(EXCLUDE)) {
  if (!file.exists(EXCLUDE)) stop("Exclusion list not found: ", EXCLUDE)
  ex <- fread(EXCLUDE)
  missing_cols <- setdiff(c("Order", "Species", "Haplotype"), names(ex))
  if (length(missing_cols)) {
    stop("Exclusion list is missing column(s): ", paste(missing_cols, collapse = ", "))
  }
  before <- nrow(formatted)
  keys   <- formatted[, paste(Order, Species, Haplotype, sep = "/")]
  ex_keys <- ex[, paste(Order, Species, Haplotype, sep = "/")]
  unmatched <- setdiff(ex_keys, unique(keys))
  formatted <- formatted[!keys %in% ex_keys]
  cat(sprintf("Excluded %d rows for %d haplotypes (%s):\n",
              before - nrow(formatted), length(intersect(ex_keys, unique(keys))), EXCLUDE))
  for (i in seq_len(nrow(ex))) {
    k <- ex_keys[i]
    if (k %in% unmatched) next
    cat(sprintf("  - %s\n      %s\n", k,
                if ("Reason" %in% names(ex)) ex$Reason[i] else "(no reason given)"))
  }
  # An entry that matches nothing usually means the upstream name changed again.
  for (k in unmatched) cat(sprintf("  ! not found in the contig list: %s\n", k))
}
formatted_igh <- formatted[Locus == "IGH"][NumV > MIN_NUMV]

# fwrite for all three -- the previous version called readr::write_tsv without
# attaching readr, so it only worked if the session happened to have it loaded.
fwrite(formatted, output_all)
fwrite(formatted_igh, output_igh, sep = "\t")
fwrite(formatted, output_all_tsv, sep = "\t")

cat(sprintf("Kept loci: %s\n", paste(LOCI, collapse = ", ")))
cat("Done. Wrote:\n")
cat(sprintf(" - %s (%d rows)\n", output_all, nrow(formatted)))
cat(sprintf(" - %s (%d rows, IGH with NumV > %d)\n",
            output_igh, nrow(formatted_igh), MIN_NUMV))
cat(sprintf(" - %s (%d rows)\n", output_all_tsv, nrow(formatted)))
