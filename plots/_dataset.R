# Shared dataset switch for the manuscript figures.
#
# The species tree and the VGP table exist in two versions (see README.md):
#
#   v1  the tree the manuscript was written against -- 122 tips, frozen,
#       not regenerable (pruned two months before its table was built)
#   v2  rebuilt by data_prep/build_vgp_tables.R -- 137 tips, strict superset
#
# Figures that touch the tree source this file and read TREE_FILE / VGP_TABLE
# instead of hardcoding a path, so switching versions is a one-line edit here
# and figures land in the matching figures/ subdirectory.
#
# Usage in a figure script:
#     source(file.path(dirname(dirname(sys.frame(1)$ofile)), "_dataset.R"))
# or, when sourcing interactively in RStudio, simply:
#     source("plots/_dataset.R")

# ---- switch this line ----
DATASET <- "v1"          # "v1" (manuscript) or "v2" (updated)
# --------------------------

if (!DATASET %in% c("v1", "v2")) stop("DATASET must be 'v1' or 'v2', got: ", DATASET)

INPUT_DIR <- "/local/storage/kav67/clean_birds"
REPO_DIR  <- "/home/kav67/Bird_IG"

.suffix    <- if (DATASET == "v2") "_v2" else ""
TREE_FILE  <- file.path(INPUT_DIR, sprintf("vgp_birds%s.nwk", .suffix))
VGP_TABLE  <- file.path(INPUT_DIR, sprintf("IGH_VGP_table%s.tsv", .suffix))
IGH_TABLE  <- file.path(INPUT_DIR, sprintf("IGH_table%s.tsv", .suffix))

# Figures are exported by hand for the manuscript (better control over font sizes
# and framing than ggsave gives). The ggsave calls in the figure scripts exist so
# the pipeline reproduces every panel unattended; treat their output as a check,
# not as the final artwork.
FIG_DIR <- file.path(REPO_DIR, "figures", DATASET)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

fig_path <- function(name) file.path(FIG_DIR, name)

# svglite handles UTF-8 (beta, lambda in the phylolm annotations); the default
# svg device silently drops them with a "conversion failure in mbcsToSbcs" warning.
# Standard canvas for the non-tree manuscript figures: 1472 x 472 px at 96 dpi,
# which is the size the hand-exported versions in figures/ already use
# (1103 x 354 pt = 1472 x 472 px). Tree figures are tall and set their own.
FIG_WIDTH_PX  <- 1472
FIG_HEIGHT_PX <- 472
FIG_DPI       <- 96

save_fig <- function(name, plot,
                     width = FIG_WIDTH_PX, height = FIG_HEIGHT_PX,
                     units = "px", dpi = FIG_DPI) {
  path <- fig_path(name)
  if (grepl("\\.svg$", name) && requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(path, plot, device = svglite::svglite,
                    width = width, height = height, units = units, dpi = dpi)
  } else {
    ggplot2::ggsave(path, plot,
                    width = width, height = height, units = units, dpi = dpi)
  }
  message(sprintf("[figure] wrote %s (%gx%g %s)", path, width, height, units))
  invisible(path)
}

message(sprintf("[dataset] %s  tree=%s  table=%s  figures -> %s",
                DATASET, basename(TREE_FILE), basename(VGP_TABLE), FIG_DIR))

for (f in c(TREE_FILE, VGP_TABLE)) {
  if (!file.exists(f)) stop("Missing dataset file: ", f,
                            "\nBuild it with data_prep/build_vgp_tables.R")
}
