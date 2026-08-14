suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(ape)
  library(phylolm)
})

# ── INPUTS ────────────────────────────────────────────────────────────────────
INPUT_DIR    <- "/local/storage/kav67/clean_birds/"
SUMMARY_FILE <- file.path(INPUT_DIR, "summary_features.csv")
INV_FILE     <- file.path(INPUT_DIR, "inversion_stats.tsv")
VGP_FILE     <- file.path(INPUT_DIR, "IGH_VGP_table.tsv")   # Species -> LatinName
TREE_FILE    <- file.path(INPUT_DIR, "vgp_birds.nwk")
FIG_DIR      <- "/home/kav67/Bird_IG/figures"
dir.create(FIG_DIR, showWarnings = FALSE)

# ── Load summary features — IGH only, one contig per haplotype ────────────────
meta <- read.csv(SUMMARY_FILE, stringsAsFactors = FALSE) %>%
  filter(Locus == "IGH") %>%
  group_by(Order, Species, Haplotype) %>%
  slice_max(NumV, n = 1, with_ties = FALSE) %>%
  ungroup()

# ── Read per-haplotype locus summary and derive full contig length ─────────────
# refined_ig_loci/summary.csv stores the locus end both as an absolute position
# (EndPos) and as a fraction of the full contig (RelEnd); ContigLength = EndPos/RelEnd.
read_contig_length <- function(order, species, haplotype, contig) {
  path <- file.path(INPUT_DIR, order, species, haplotype,
                    "refined_ig_loci", "summary.csv")
  if (!file.exists(path)) return(NA_real_)
  d <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(d)) return(NA_real_)
  row <- d %>% filter(Contig == contig, Locus == "IGH")
  if (nrow(row) == 0 || is.na(row$RelEnd[1]) || row$RelEnd[1] == 0) return(NA_real_)
  row$EndPos[1] / row$RelEnd[1]
}

contig_lengths <- meta %>%
  rowwise() %>%
  mutate(contig_length_bp = read_contig_length(Order, Species, Haplotype, Contig)) %>%
  ungroup() %>%
  filter(!is.na(contig_length_bp))

message(sprintf("Contig lengths computed for %d / %d haplotypes",
                nrow(contig_lengths), nrow(meta)))

# ── Load inversion stats and join ─────────────────────────────────────────────
# inversion_stats.tsv has one row per (haplotype x contig): when the IGH locus is
# split across contigs a haplotype has several rows (342 haplotypes -> 417 rows),
# distinguished only by total_seq_length (there is no contig-name column). The
# contig length above is measured on the MAIN contig (highest NumV), so keep the
# matching inversion row: the one with the largest locus (max total_seq_length).
# Without this, the join fans out and pairs a haplotype's main-contig length with
# inversion counts from its other contigs.
inv <- read.delim(INV_FILE, stringsAsFactors = FALSE) %>%
  filter(minlen == 250) %>%
  group_by(order, species, haplotype) %>%
  slice_max(total_seq_length, n = 1, with_ties = FALSE) %>%
  ungroup()

combined <- contig_lengths %>%
  inner_join(inv,
             by = c("Order" = "order", "Species" = "species", "Haplotype" = "haplotype")) %>%
  mutate(contig_length_mb = contig_length_bp / 1e6)

message(sprintf("%d haplotypes (%d species) after joining with inversion stats",
                nrow(combined), n_distinct(combined$Species)))

# ── Correlation helper ────────────────────────────────────────────────────────
# Spearman is the primary statistic: both variables are strongly right-skewed,
# non-normal even after log, contain outliers, and num_inversions has zeros.
# Rank correlation is robust to all of these; log-log Pearson is reported as a
# secondary check. Note both treat haplotypes as independent, which they are NOT
# (see phylogenetic section below) — so these p-values are anticonservative.
corr_label <- function(x, y) {
  sp <- suppressWarnings(cor.test(x, y, method = "spearman"))
  pc  <- min(y[y > 0]) / 2
  lp <- cor.test(log10(x), log10(y + pc))               # log-log Pearson
  sprintf("Spearman ρ = %.2f (p = %.1e, n = %d)   |   log–log Pearson r = %.2f",
          sp$estimate, sp$p.value, length(x), lp$estimate)
}

scatter_plot <- function(df, yvar, ylab, title) {
  lab <- corr_label(df$contig_length_mb, df[[yvar]])
  ggplot(df, aes(x = contig_length_mb, y = .data[[yvar]])) +
    geom_point(alpha = 0.7, size = 2, colour = "#87b4dc") +
    geom_smooth(method = "lm", se = TRUE, colour = "grey30", linewidth = 0.7) +
    scale_x_log10(labels = function(x) paste0(x, " Mb")) +
    scale_y_log10() +
    labs(x = "Contig length (Mb, log scale)", y = ylab,
         title = title, subtitle = lab) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 10, colour = "grey20"))
}

p1 <- scatter_plot(combined, "num_inversions",
                   "Number of inversions (log)",
                   "Contig length vs. number of inversions (IGH)")
p2 <- scatter_plot(combined, "genes_on_inv",
                   "Genes in inversion-associated regions (log)",
                   "Contig length vs. genes in inversion-associated regions (IGH)")

ggsave(file.path(FIG_DIR, "contig_length_vs_inversions.svg"), p1,
       width = 6.5, height = 5)
ggsave(file.path(FIG_DIR, "contig_length_vs_genes_on_inv.svg"), p2,
       width = 6.5, height = 5)

# ── Console summary of naive correlations ─────────────────────────────────────
sp <- suppressWarnings(cor.test(combined$contig_length_mb, combined$num_inversions,
                                method = "spearman"))
message(sprintf("Naive Spearman (haplotypes): rho = %.3f, p = %.3e, n = %d",
                sp$estimate, sp$p.value, nrow(combined)))

# ── Phylogenetic version ──────────────────────────────────────────────────────
# The haplotypes come from only ~163 species (many with several haplotypes) and
# species are phylogenetically correlated. Both violate the independence
# assumption behind the naive correlation. Mirror phylolm_tree.R, but instead of
# averaging over a species' haplotypes, keep ONE primary haplotype per species:
# prefer names ending in .pri/_pri, then the VGP reference assembly (b... IDs),
# with explicit overrides where wanted. Then log-transform and fit a
# Pagel's-lambda phylogenetic regression on the VGP tree.

# Explicit per-species haplotype picks (Species -> Haplotype). Population-sampled
# species carry many "_pri" haplotypes; pin the reference assembly here.
HAP_OVERRIDE <- c(House_Finch = "bHaeMex1.pri")

# Species that have IGH data (contig length + inversions) and a tip on the VGP
# tree, but are absent from IGH_VGP_table.tsv, so would otherwise be dropped at
# the LatinName join. LatinNames are taken from the tree tips / assembly IDs.
# Four of these have only an alt haplotype (no _pri); the haplotype-scoring below
# falls back to the alt automatically.
LATIN_FALLBACK <- tibble::tribble(
  ~Order,      ~Species,                ~LatinName,
  "Songbirds", "Nelsons_Sparrow",       "Ammospiza nelsoni",     # bAmmNel1_alt (alt only)
  "Songbirds", "Swamp_Sparrow",         "Melospiza georgiana",   # bMelGeo1_alt (alt only)
  "Songbirds", "WhiteThroated_Sparrow", "Zonotrichia albicollis",# bZonAlb1_pri
  "MiscBirds", "Speckled_Mousebird",    "Colius striatus",       # bColStr4_alt (alt only)
  "MiscBirds", "Common_Swift",          "Apus apus",             # bApuApu2_pri
  "Landfowl",  "King_Quail",            "Coturnix chinensis"     # bCotChi1_alt (alt only)
)

run_phylolm <- function(combined) {
  # LatinName is constant per species -> use the VGP table only for the mapping,
  # not for which haplotype to keep. Add the in-script fallbacks for species that
  # are on the tree but missing from IGH_VGP_table.tsv.
  sp_latin <- read.delim(VGP_FILE, stringsAsFactors = FALSE) %>%
    distinct(Order, Species, LatinName) %>%
    bind_rows(LATIN_FALLBACK) %>%
    distinct(Order, Species, .keep_all = TRUE)
  tree <- read.tree(TREE_FILE)
  tree_species <- gsub('"', '', tree$tip.label)

  joined <- combined %>% inner_join(sp_latin, by = c("Order", "Species"))

  # Score each haplotype; highest score wins, ties broken by name.
  ov  <- HAP_OVERRIDE[joined$Species]
  joined <- joined %>%
    mutate(.score = 1L * grepl("[._]pri$", Haplotype, ignore.case = TRUE) +
                    2L * grepl("^b[A-Za-z]+[0-9]", Haplotype) +
                    100L * (!is.na(ov) & ov == Haplotype))
  sel <- joined %>%
    group_by(Order, Species) %>%
    arrange(desc(.score), Haplotype, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()

  # Report the pick only for species that had >1 PRIMARY haplotype (the genuinely
  # ambiguous cases: population-sampled species). pri+alt pairs are unremarkable.
  multi_pri <- joined %>%
    filter(grepl("[._]pri$", Haplotype, ignore.case = TRUE)) %>%
    count(Order, Species) %>% filter(n > 1)
  picked <- sel %>% semi_join(multi_pri, by = c("Order", "Species"))
  if (nrow(picked) > 0) {
    message("Primary haplotype chosen for species with multiple '_pri' haplotypes:")
    for (i in seq_len(nrow(picked)))
      message(sprintf("  %-22s (%d pri available) -> %s",
                      picked$Species[i],
                      multi_pri$n[match(picked$Species[i], multi_pri$Species)],
                      picked$Haplotype[i]))
  }

  trait <- sel %>%
    transmute(LatinName,
              contig_length = log(contig_length_mb),
              num_inversions = log(num_inversions)) %>%
    filter(is.finite(contig_length), is.finite(num_inversions))

  trait <- trait[match(tree_species, gsub(" ", "_", trait$LatinName)), ]
  keep  <- !(is.na(trait$contig_length) | is.na(trait$num_inversions))
  trait <- trait[keep, ]
  tree  <- drop.tip(tree, tree$tip.label[!keep])
  rownames(trait) <- tree$tip.label

  message(sprintf("phylolm: %d species matched to tree", nrow(trait)))

  model <- phylolm(num_inversions ~ contig_length, data = trait,
                   phy = tree, model = "lambda")
  co <- summary(model)$coefficients
  slope <- co["contig_length", "Estimate"]
  pval  <- co["contig_length", "p.value"]
  lambda <- model$optpar
  r2 <- summary(model)$r.squared
  message(sprintf("phylolm num_inversions ~ contig_length (per species, log):"))
  message(sprintf("  slope = %.3f, lambda = %.3f, R2 = %.3f, p = %.3e",
                  slope, lambda, r2, pval))

  stars <- ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**",
                  ifelse(pval < 0.05, "*", "ns")))
  # Two aligned columns: right-aligned "name =" (so the '=' line up) and
  # left-aligned values. No surrounding box.
  names_lab <- "β =\nλ =\nR² =\np ="
  vals_lab  <- sprintf("%.3f%s\n%.3f\n%.3f\n%.2e", slope, stars, lambda, r2, pval)

  # Upper-right, nudged inward for buffer on both sides.
  xr <- range(trait$contig_length, na.rm = TRUE)
  yr <- range(trait$num_inversions, na.rm = TRUE)
  w  <- xr[2] - xr[1]; h <- yr[2] - yr[1]
  x_eq <- xr[1] + 0.82 * w          # x of the aligned '=' column
  x_val <- xr[1] + 0.845 * w        # values start just to its right
  ann_y <- yr[1] + 0.98 * h

  p <- ggplot(trait, aes(x = contig_length, y = num_inversions)) +
    geom_point(size = 2.4, alpha = 0.9, colour = "#87b4dc") +
    geom_abline(intercept = coef(model)[1], slope = slope,
                linetype = "dashed", linewidth = 0.6) +
    annotate("text", x = x_eq, y = ann_y, label = names_lab,
             hjust = 1, vjust = 1, size = 4.2, lineheight = 1.2) +
    annotate("text", x = x_val, y = ann_y, label = vals_lab,
             hjust = 0, vjust = 1, size = 4.2, lineheight = 1.2) +
    labs(x = "Contig length (log, Mb)", y = "# Inversions (log)") +
    theme_classic(base_size = 12) +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 10))
  ggsave(file.path(FIG_DIR, "contig_length_vs_inversions_phylolm.svg"), p,
         width = 6.5, height = 5)
  invisible(model)
  return(p)
}

if (file.exists(VGP_FILE) && file.exists(TREE_FILE)) {
  p<-run_phylolm(combined)
} else {
  message("VGP table or tree not found; skipping phylogenetic model.")
}
p
