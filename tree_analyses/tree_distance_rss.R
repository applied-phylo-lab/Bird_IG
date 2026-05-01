#!/usr/bin/env Rscript
# Tree distance analysis: IGH vs IGL, RSS-pair vs any-pair distances.
# Designed to be sourced or run section-by-section in RStudio.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ape)
  library(ggplot2)
  library(patchwork)
})

# ── INPUTS ────────────────────────────────────────────────────────────────────
INPUT_DIR  <- "/local/storage/kav67/clean_birds/"
SUMMARY    <- "/local/storage/kav67/clean_birds/filtered_table.tsv"
GENE_LIST  <- "/local/storage/kav67/clean_birds/gene_list.csv"
OUTPUT_DIR <- "/local/storage/kav67/clean_birds/trees/"

LOCUS_COLORS <- c(IGH = "#87b4dc", IGL = "#638E6E")

# ── READ INPUTS ───────────────────────────────────────────────────────────────
meta <- read_tsv(SUMMARY, show_col_types = FALSE) %>%
  group_by(Order, Species, Haplotype, Locus) %>%
  slice_max(NumV, n = 1, with_ties = FALSE) %>%
  ungroup()

gene_list <- read.csv(GENE_LIST, stringsAsFactors = FALSE, check.names = FALSE) %>%
  separate(Source, into = c("gl_Order", "gl_Species", "gl_Haplotype"),
           sep = "/", remove = FALSE) %>%
  mutate(has_rss = !is.na(Heptamer) & trimws(Heptamer) != "")

# Pre-build a lookup: haplotype × locus → integer vector of RSS positions
rss_lookup <- gene_list %>%
  filter(has_rss) %>%
  group_by(gl_Haplotype, Locus) %>%
  summarise(rss_pos = list(as.integer(Pos)), .groups = "drop")

# ── HELPERS ───────────────────────────────────────────────────────────────────
tree_prefix <- function(haplotype, locus) {
  if (locus == "IGL") paste0("IGL_", haplotype) else haplotype
}

# Tip labels: {prefix}.{Pos}.{Contig}.{GeneType}.{Productive}.{Strand}
# Strip the prefix (which may itself contain dots) by matching it as a literal
# prefix, then grab the first run of digits as the position.
tip_to_pos <- function(tips, prefix) {
  stripped <- sub(paste0("^", prefix, "\\."), "", tips)
  as.integer(str_extract(stripped, "^\\d+"))
}

get_rss_pos <- function(haplotype, locus) {
  idx <- rss_lookup$gl_Haplotype == haplotype & rss_lookup$Locus == locus
  if (!any(idx)) return(integer(0))
  rss_lookup$rss_pos[[which(idx)[1]]]
}

# ── MAIN LOOP ─────────────────────────────────────────────────────────────────
all_pairs    <- list()
summary_rows <- list()

for (i in seq_len(nrow(meta))) {
  row       <- meta[i, ]
  order     <- row$Order
  species   <- row$Species
  haplotype <- row$Haplotype
  locus     <- row$Locus
  contig    <- row$Contig

  pfx       <- tree_prefix(haplotype, locus)
  tree_file <- file.path(INPUT_DIR, order, species, haplotype,
                         "tree", paste0(pfx, "_tree.treefile"))

  if (!file.exists(tree_file)) {
    message("Skipping missing tree: ", tree_file)
    next
  }

  tree <- tryCatch(read.tree(tree_file), error = function(e) {
    message("Failed to read: ", tree_file); NULL
  })
  if (is.null(tree)) next

  # Restrict to the contig of interest
  keep_tips <- tree$tip.label[str_detect(tree$tip.label, fixed(contig))]
  if (length(keep_tips) < 2) {
    message("Too few tips on contig (", locus, "): ", haplotype)
    next
  }
  tree_sub <- keep.tip(tree, keep_tips)

  # Match tips → RSS status
  tip_pos    <- tip_to_pos(tree_sub$tip.label, pfx)
  rss_pos    <- get_rss_pos(haplotype, locus)
  has_rss    <- setNames(tip_pos %in% rss_pos, tree_sub$tip.label)

  # Pairwise cophenetic distances (upper triangle only)
  dist_mat   <- cophenetic.phylo(tree_sub)
  tips       <- rownames(dist_mat)

  pairs_df <- as.data.frame(as.table(dist_mat),
                             stringsAsFactors = FALSE) %>%
    dplyr::rename(Gene1 = Var1, Gene2 = Var2, dist = Freq) %>%
    filter(Gene1 < Gene2) %>%
    mutate(
      rss1      = has_rss[Gene1],
      rss2      = has_rss[Gene2],
      both_rss  = rss1 & rss2,
      Order     = order,
      Species   = species,
      Haplotype = haplotype,
      Locus     = locus
    )

  all_pairs[[length(all_pairs) + 1]] <- pairs_df

  any_dists <- pairs_df$dist
  rss_dists <- pairs_df$dist[pairs_df$both_rss]

  summary_rows[[length(summary_rows) + 1]] <- tibble(
    Order         = order,
    Species       = species,
    Haplotype     = haplotype,
    Locus         = locus,
    n_genes       = length(tips),
    n_rss_genes   = sum(has_rss),
    mean_dist_any = mean(any_dists),
    var_dist_any  = var(any_dists),
    n_any_pairs   = length(any_dists),
    mean_dist_rss = if (length(rss_dists) >= 1) mean(rss_dists) else NA_real_,
    var_dist_rss  = if (length(rss_dists) >= 2) var(rss_dists)  else NA_real_,
    n_rss_pairs   = length(rss_dists)
  )
}

pairs_all  <- bind_rows(all_pairs)
summary_df <- bind_rows(summary_rows)

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
write_tsv(summary_df, file.path(OUTPUT_DIR, "tree_distance_summary.tsv"))
write_tsv(pairs_all,  file.path(OUTPUT_DIR, "tree_distance_pairs.tsv"))
message("Saved: ", nrow(summary_df), " haplotype×locus rows")

# ── PLOT 1: Mean distance IGH vs IGL (per haplotype) ─────────────────────────
# select() before pivot_wider so only the id + value columns are kept;
# otherwise every numeric column becomes an id and pivot creates list columns.
scatter_df <- summary_df %>%
  select(Order, Species, Haplotype, Locus, mean_dist_any, mean_dist_rss) %>%
  pivot_wider(
    names_from  = Locus,
    values_from = c(mean_dist_any, mean_dist_rss)
  ) %>%
  filter(!is.na(mean_dist_any_IGH) & !is.na(mean_dist_any_IGL))

make_scatter <- function(df, x_col, y_col, x_lab, y_lab, title) {
  if (nrow(df) == 0) {
    message("No data for scatter: ", title)
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5, label = "No paired data yet",
                      size = 5, color = "grey50") +
             theme_void() + ggtitle(title))
  }
  vals <- c(df[[x_col]], df[[y_col]])
  lim  <- range(vals, na.rm = TRUE)
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color = Order)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "grey60", linewidth = 0.5) +
    geom_point(size = 2.5, alpha = 0.8) +
    coord_equal(xlim = lim, ylim = lim) +
    labs(x = x_lab, y = y_lab, title = title) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")
}

p_scatter_any <- make_scatter(
  scatter_df,
  "mean_dist_any_IGH", "mean_dist_any_IGL",
  "Mean pairwise distance – IGH", "Mean pairwise distance – IGL",
  "All-pair mean tree distance: IGH vs IGL"
)

scatter_rss <- scatter_df %>%
  filter(!is.na(mean_dist_rss_IGH) & !is.na(mean_dist_rss_IGL))

p_scatter_rss <- make_scatter(
  scatter_rss,
  "mean_dist_rss_IGH", "mean_dist_rss_IGL",
  "Mean RSS-pair distance – IGH", "Mean RSS-pair distance – IGL",
  "RSS-only pair mean tree distance: IGH vs IGL"
)

# ── PLOT 2: RSS-pair vs any-pair distances, per locus ─────────────────────────
# Use per-haplotype mean distances (not raw pairs) to avoid pseudoreplication.
dist_comparison <- summary_df %>%
  filter(!is.na(mean_dist_any)) %>%
  select(Order, Species, Haplotype, Locus, mean_dist_any, mean_dist_rss) %>%
  pivot_longer(
    cols      = c(mean_dist_any, mean_dist_rss),
    names_to  = "pair_type",
    values_to = "mean_dist"
  ) %>%
  mutate(
    pair_type = recode(pair_type,
                       mean_dist_any = "Any pair",
                       mean_dist_rss = "RSS pair"),
    Locus = factor(Locus, levels = c("IGH", "IGL"))
  ) %>%
  filter(!is.na(mean_dist))

p_dist_comparison <- ggplot(dist_comparison,
                            aes(x = pair_type, y = mean_dist,
                                fill = Locus, color = Locus)) +
  geom_violin(alpha = 0.5, trim = FALSE, linewidth = 0.4) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.8,
               color = "black", linewidth = 0.4) +
  geom_jitter(width = 0.06, size = 0.8, alpha = 0.4) +
  facet_wrap(~Locus, ncol = 2) +
  scale_fill_manual(values  = LOCUS_COLORS) +
  scale_color_manual(values = LOCUS_COLORS) +
  labs(
    x     = "Pair type",
    y     = "Mean pairwise tree distance (per haplotype)",
    title = "RSS-pair vs any-pair distances by locus"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "none",
    strip.text      = element_text(face = "bold"),
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 10)
  )

# ── PLOT 3: Difference (RSS-pair mean − any-pair mean) per haplotype ──────────
# Positive = RSS pairs are more distant; negative = RSS pairs are closer.
diff_df <- summary_df %>%
  filter(!is.na(mean_dist_any), !is.na(mean_dist_rss)) %>%
  mutate(
    dist_diff = mean_dist_rss - mean_dist_any,
    Locus     = factor(Locus, levels = c("IGH", "IGL"))
  )

p_dist_diff <- ggplot(diff_df,
                      aes(x = Locus, y = dist_diff,
                          fill = Locus, color = Locus)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  geom_violin(alpha = 0.5, trim = FALSE, linewidth = 0.4) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.8,
               color = "black", linewidth = 0.4) +
  geom_jitter(width = 0.06, size = 1.2, alpha = 0.5) +
  scale_fill_manual(values  = LOCUS_COLORS) +
  scale_color_manual(values = LOCUS_COLORS) +
  labs(
    x     = "Locus",
    y     = "Mean RSS-pair dist − Mean any-pair dist",
    title = "Difference: RSS pairs vs any pairs"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "none",
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 10)
  )

# ── PLOT 3b: dist_diff vs RSS fraction — checks for count-driven bias ─────────
# When many genes have RSS (high fraction), mean_dist_rss converges toward
# mean_dist_any by definition, pushing the difference toward 0 regardless of
# biology. This plot reveals whether the pattern in p_slopegraph is confounded
# by the number of RSS genes.
diff_df <- summary_df %>%
  filter(!is.na(mean_dist_any), !is.na(mean_dist_rss)) %>%
  mutate(
    rss_fraction = n_rss_genes / n_genes,
    dist_diff    = mean_dist_rss - mean_dist_any,
    Locus        = factor(Locus, levels = c("IGH", "IGL"))
  )

p_diff_vs_fraction <- ggplot(diff_df,
                             aes(x = rss_fraction, y = dist_diff,
                                 color = Locus)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.15) +
  facet_wrap(~Locus, ncol = 2) +
  scale_color_manual(values = LOCUS_COLORS) +
  scale_x_continuous(labels = scales::percent_format()) +
  labs(
    x     = "Fraction of genes with RSS",
    y     = "Mean RSS-pair dist − Mean any-pair dist",
    title = "Bias check: does RSS gene count drive the difference?"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "none",
    strip.text      = element_text(face = "bold"),
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 10)
  )

# ── PLOT 4: Mean distance IGH vs IGL (distribution, any-pair) ────────────────
# Shows whether one locus has systematically higher within-locus divergence.
p_locus_dist <- ggplot(
    summary_df %>% filter(!is.na(mean_dist_any)) %>%
      mutate(Locus = factor(Locus, levels = c("IGH", "IGL"))),
    aes(x = Locus, y = mean_dist_any, fill = Locus, color = Locus)
  ) +
  geom_violin(alpha = 0.5, trim = FALSE, linewidth = 0.4) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.8,
               color = "black", linewidth = 0.4) +
  geom_jitter(width = 0.06, size = 1.2, alpha = 0.5) +
  scale_fill_manual(values  = LOCUS_COLORS) +
  scale_color_manual(values = LOCUS_COLORS) +
  labs(
    x     = "Locus",
    y     = "Mean pairwise tree distance (per haplotype)",
    title = "Within-locus divergence: IGH vs IGL"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "none",
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 10)
  )

# ── PRINT ─────────────────────────────────────────────────────────────────────
print(p_scatter_any)
print(p_scatter_rss)
print(p_locus_dist)
print(p_dist_comparison)
print(p_dist_diff)
print(p_diff_vs_fraction)

# Combined layout
p_scatter_any | p_scatter_rss
(p_locus_dist | p_dist_comparison) /(p_dist_diff | p_diff_vs_fraction)

