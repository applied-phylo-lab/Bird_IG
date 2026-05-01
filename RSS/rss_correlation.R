#!/usr/bin/env Rscript
# Correlation between total genes and genes with RSS, one plot per locus.
# Uses phylolm (lambda model) to account for phylogenetic signal.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(ape)
  library(phylolm)
  library(ggrepel)
})

input_dir <- "/local/storage/kav67/clean_birds/"

# ── Load data ─────────────────────────────────────────────────────────────────

gene_list_path <- file.path(input_dir, "gene_list.csv")
if (!file.exists(gene_list_path)) stop("gene_list.csv not found: ", gene_list_path)

df <- read.csv(gene_list_path, stringsAsFactors = FALSE, check.names = FALSE)

df <- df %>%
  separate(Source, into = c("Order", "Species", "Haplotype"),
           sep = "/", remove = FALSE) %>%
  mutate(has_rss = !is.na(Heptamer) & trimws(Heptamer) != "")

# Add LatinName by joining with IGH_VGP_table.tsv (Species + Haplotype → LatinName)
vgp_table <- read.delim(file.path(input_dir, "IGH_VGP_table.tsv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
df <- df %>%
  left_join(vgp_table %>% select(Species, LatinName) %>% distinct(),
            by = c("Species"))

# ── Load and process tree ─────────────────────────────────────────────────────

bird_tree <- read.tree(file.path(input_dir, "vgp_birds.nwk"))

# Extract plain species names from tip labels (same logic as phylolm_tree.R)
tree_species <- gsub('.*"', '', bird_tree$tip.label)
tree_species <- gsub('"', '', tree_species)

# ── Aggregate per species × locus ─────────────────────────────────────────────
# Per haplotype first, then average per LatinName so each species is one point

counts <- df %>%
  mutate(productive = tolower(trimws(as.character(Productive))) %in% c("true", "1", "yes")) %>%
  group_by(Order, LatinName, Haplotype, Locus) %>%
  summarise(
    n_total        = n(),
    n_rss          = sum(has_rss),
    n_productive   = sum(productive),
    n_prod_rss     = sum(productive & has_rss),
    .groups = "drop"
  ) %>%
  group_by(Order, LatinName, Locus) %>%
  summarise(
    n_total      = mean(n_total),
    n_rss        = mean(n_rss),
    n_productive = mean(n_productive),
    n_prod_rss   = mean(n_prod_rss),
    .groups = "drop"
  )

# ── Phylolm helper ────────────────────────────────────────────────────────────

focal_orders  <- c("Waterfowl", "Landfowl")
order_palette <- c(Galloanserae = "#C4713B", Other = "#CCCCCC")

make_phylolm_label <- function(model) {
  s      <- summary(model)
  slope  <- coef(model)[2]
  p_val  <- s$coefficients[2, "p.value"]
  lambda <- model$optpar
  r2     <- s$r.squared
  stars  <- ifelse(p_val < 0.001, "***",
             ifelse(p_val < 0.01,  "**",
              ifelse(p_val < 0.05,  "*", "ns")))
  sprintf("β  = %-8s\nλ  = %-8.3f\nR² = %-8.3f",
          paste0(round(slope, 3), stars), lambda, r2)
}

# Shared helper: build trait_df and pruned tree for a given locus
prep_locus <- function(locus, counts, bird_tree, tree_species) {
  sub <- filter(counts, Locus == locus)
  sub$latin_tree <- gsub(" ", "_", sub$LatinName)
  idx         <- match(tree_species, sub$latin_tree)
  trait_df    <- as.data.frame(sub[idx, ])
  keep        <- !is.na(trait_df$n_total) & !is.na(trait_df$n_rss)
  trait_df    <- trait_df[keep, ]
  tree_pruned <- drop.tip(bird_tree, bird_tree$tip.label[!keep])
  rownames(trait_df) <- tree_pruned$tip.label
  if (nrow(trait_df) < 5) return(NULL)
  list(trait_df = trait_df, tree_pruned = tree_pruned)
}

# ── Version 1: original, colored by Order ─────────────────────────────────────
run_locus <- function(locus, counts, bird_tree, tree_species) {
  prep <- prep_locus(locus, counts, bird_tree, tree_species)
  if (is.null(prep)) {
    message("Skipping locus ", locus, ": too few species"); return(NULL)
  }
  trait_df    <- prep$trait_df
  tree_pruned <- prep$tree_pruned

  model <- phylolm(n_rss ~ n_total, data = trait_df, phy = tree_pruned, model = "lambda")

  ggplot(trait_df, aes(x = n_total, y = n_rss, color = Order)) +
    geom_point(size = 2.5, alpha = 0.9) +
    geom_abline(intercept = coef(model)[1], slope = coef(model)[2],
                linetype = "dashed", linewidth = 0.5, color = "black") +
    annotate("text", x = -Inf, y = Inf, label = make_phylolm_label(model),
             hjust = -0.05, vjust = 1.3, size = 3.2, family = "mono") +
    labs(title = paste("Locus:", locus),
         x = "Total genes (mean per species)",
         y = "Genes with RSS (mean per species)") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          axis.title = element_text(size = 12),
          axis.text  = element_text(size = 10))
}

# ── Version 2: Galloanserae highlighted, with group parameters ────────────────
run_locus_galloans <- function(locus, counts, bird_tree, tree_species) {
  prep <- prep_locus(locus, counts, bird_tree, tree_species)
  if (is.null(prep)) {
    message("Skipping locus ", locus, ": too few species"); return(NULL)
  }
  trait_df    <- prep$trait_df
  tree_pruned <- prep$tree_pruned

  trait_df$group <- ifelse(trait_df$Order %in% focal_orders, "Galloanserae", "Other")

  model_all <- phylolm(n_rss ~ n_total, data = trait_df, phy = tree_pruned, model = "lambda")

  p <- ggplot(trait_df, aes(x = n_total, y = n_rss, color = group)) +
    geom_point(size = 2.5, alpha = 0.9) +
    geom_abline(intercept = coef(model_all)[1], slope = coef(model_all)[2],
                color = "black", linetype = "dashed", linewidth = 0.6) +
    annotate("text", x = -Inf, y = Inf,
             label = paste0("All:\n", make_phylolm_label(model_all)),
             hjust = -0.05, vjust = 1.3, size = 3.2, color = "black", family = "mono") +
    scale_color_manual(values = order_palette,
                       breaks = c("Galloanserae", "Other"), name = "Group") +
    labs(title = paste("Locus:", locus, "— Galloanserae highlighted"),
         x = "Total genes (mean per species)",
         y = "Genes with RSS (mean per species)") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          axis.title = element_text(size = 12),
          axis.text  = element_text(size = 10))

  # Galloanserae sub-model
  grp_df <- trait_df[trait_df$group == "Galloanserae", , drop = FALSE]
  if (nrow(grp_df) >= 5) {
    tree_grp  <- drop.tip(tree_pruned, setdiff(tree_pruned$tip.label, rownames(grp_df)))
    model_grp <- tryCatch(
      phylolm(n_rss ~ n_total, data = grp_df, phy = tree_grp, model = "lambda"),
      error = function(e) { message("Galloanserae/", locus, " failed: ", e$message); NULL }
    )
    if (!is.null(model_grp)) {
      # Interaction test
      trait_df$is_focal <- as.integer(trait_df$group == "Galloanserae")
      model_int <- tryCatch(
        phylolm(n_rss ~ n_total * is_focal, data = trait_df, phy = tree_pruned, model = "lambda"),
        error = function(e) NULL
      )
      delta_label <- if (!is.null(model_int)) {
        s       <- summary(model_int)$coefficients
        int_row <- grep("n_total.*is_focal|is_focal.*n_total", rownames(s), value = TRUE)
        if (length(int_row)) {
          p_int  <- s[int_row, "p.value"]
          dslope <- s[int_row, "Estimate"]
          stars  <- ifelse(p_int < 0.001, "***", ifelse(p_int < 0.01, "**",
                     ifelse(p_int < 0.05, "*", "ns")))
          sprintf("Δβ vs rest = %s%s", round(dslope, 3), stars)
        } else ""
      } else ""

      grp_label <- paste0("Galloanserae:\n", make_phylolm_label(model_grp),
                          if (nchar(delta_label)) paste0("\n", delta_label) else "")

      p <- p +
        geom_abline(intercept = coef(model_grp)[1], slope = coef(model_grp)[2],
                    color = order_palette["Galloanserae"], linetype = "dashed", linewidth = 0.6) +
        annotate("text", x = -Inf, y = -Inf, label = grp_label,
                 hjust = -0.05, vjust = -0.3,
                 size = 3.2, color = order_palette["Galloanserae"], family = "mono")
    }
  } else {
    message("Skipping Galloanserae for locus ", locus, ": only ", nrow(grp_df), " species")
  }
  p
}

# ── Per-locus plots ───────────────────────────────────────────────────────────

loci        <- unique(counts$Locus)
plots       <- Filter(Negate(is.null), lapply(loci, run_locus,
                counts = counts, bird_tree = bird_tree, tree_species = tree_species))
plots_galloans <- Filter(Negate(is.null), lapply(loci, run_locus_galloans,
                counts = counts, bird_tree = bird_tree, tree_species = tree_species))
plots
plots_galloans
# ── Combined plot (all loci, color by Locus) ──────────────────────────────────

locus_colors <- c(IGH = "#87b4dc", IGL = "#638E6E")

# Keep only species present in tree
counts_in_tree <- counts %>%
  mutate(latin_tree = gsub(" ", "_", LatinName)) %>%
  filter(latin_tree %in% tree_species)

# Run phylolm per locus and collect regression lines + annotations
combined_fits <- lapply(unique(counts_in_tree$Locus), function(locus) {
  sub <- filter(counts_in_tree, Locus == locus)
  sub$latin_tree <- gsub(" ", "_", sub$LatinName)
  idx         <- match(tree_species, sub$latin_tree)
  trait_df    <- as.data.frame(sub[idx, ])
  keep        <- !is.na(trait_df$n_total) & !is.na(trait_df$n_rss)
  trait_df    <- trait_df[keep, ]
  tree_pruned <- drop.tip(bird_tree, bird_tree$tip.label[!keep])
  rownames(trait_df) <- tree_pruned$tip.label

  if (nrow(trait_df) < 5) return(NULL)

  model <- phylolm(n_rss ~ n_total, data = trait_df, phy = tree_pruned, model = "lambda")
  list(
    locus     = locus,
    intercept = coef(model)[1],
    slope     = coef(model)[2],
    label     = make_phylolm_label(model)
  )
})
combined_fits <- Filter(Negate(is.null), combined_fits)

# ── Galloanserae model on the combined data ────────────────────────────────────
galloans_fits <- lapply(unique(counts_in_tree$Locus), function(locus) {
  sub <- filter(counts_in_tree, Locus == locus)
  sub$latin_tree <- gsub(" ", "_", sub$LatinName)
  sub$group      <- ifelse(sub$Order %in% focal_orders, "Galloanserae", "Other")
  idx         <- match(tree_species, sub$latin_tree)
  trait_df    <- as.data.frame(sub[idx, ])
  keep        <- !is.na(trait_df$n_total) & !is.na(trait_df$n_rss)
  trait_df    <- trait_df[keep, ]
  tree_pruned <- drop.tip(bird_tree, bird_tree$tip.label[!keep])
  rownames(trait_df) <- tree_pruned$tip.label

  grp_df <- trait_df[trait_df$group == "Galloanserae", , drop = FALSE]
  if (nrow(grp_df) < 5) return(NULL)

  tree_grp  <- drop.tip(tree_pruned, setdiff(tree_pruned$tip.label, rownames(grp_df)))
  model_grp <- tryCatch(
    phylolm(n_rss ~ n_total, data = grp_df, phy = tree_grp, model = "lambda"),
    error = function(e) NULL
  )
  if (is.null(model_grp)) return(NULL)

  # Interaction: Galloanserae slope vs rest
  trait_df$is_focal <- as.integer(trait_df$group == "Galloanserae")
  model_int <- tryCatch(
    phylolm(n_rss ~ n_total * is_focal, data = trait_df, phy = tree_pruned, model = "lambda"),
    error = function(e) NULL
  )
  int_coef <- if (!is.null(model_int)) {
    s       <- summary(model_int)$coefficients
    int_row <- grep("n_total.*is_focal|is_focal.*n_total", rownames(s), value = TRUE)
    if (length(int_row)) s[int_row, , drop = FALSE] else NULL
  } else NULL

  delta_label <- if (!is.null(int_coef)) {
    p_int  <- int_coef[1, "p.value"]
    dslope <- int_coef[1, "Estimate"]
    stars  <- ifelse(p_int < 0.001, "***", ifelse(p_int < 0.01, "**",
               ifelse(p_int < 0.05, "*", "ns")))
    sprintf("Δβ vs rest = %s%s", round(dslope, 3), stars)
  } else ""

  list(
    locus     = locus,
    intercept = coef(model_grp)[1],
    slope     = coef(model_grp)[2],
    label     = paste0("Galloanserae:\n", make_phylolm_label(model_grp),
                       if (nchar(delta_label)) paste0("\n", delta_label) else "")
  )
})
galloans_fits <- Filter(Negate(is.null), galloans_fits)
galloans_fits
x_max <- max(counts_in_tree$n_total, na.rm = TRUE)

# Add Galloanserae group column for point shading
counts_in_tree <- counts_in_tree %>%
  mutate(group = ifelse(Order %in% focal_orders, "Galloanserae", "Other"))

# Species to highlight individually
species_colors <- c(
  "Gallus gallus"   = "#B22222",
  "Aythya fuligula" = "#6A0DAD"
)

highlighted <- counts_in_tree %>%
  filter(
    (LatinName == "Gallus gallus"   & Locus == "IGH") |
    (LatinName == "Aythya fuligula" & Locus == "IGH")
  ) %>%
  mutate(sp_color = species_colors[LatinName])

p_combined <- ggplot(counts_in_tree,
                     aes(x = n_total, y = n_rss, color = Locus, shape = group)) +
  geom_point(aes(alpha = group), size = 2.5) +
  geom_point(data = highlighted, aes(x = n_total, y = n_rss),
             color = highlighted$sp_color, shape = 18, size = 4, inherit.aes = FALSE) +
  ggrepel::geom_text_repel(data = highlighted,
             aes(x = n_total, y = n_rss, label = LatinName),
             color = highlighted$sp_color, size = 3, fontface = "italic",
             inherit.aes = FALSE) +
  scale_alpha_manual(values = c(Galloanserae = 1, Other = 0.5), guide = "none") +
  scale_shape_manual(values = c(Galloanserae = 17, Other = 16), name = "Group")

for (fit in combined_fits) {
  color <- locus_colors[fit$locus]
  label <- paste0(fit$locus, ":\n", fit$label)
  if (fit$locus == "IGL") {
    ann_x <- -Inf; ann_y <- Inf; ann_hjust <- -0.05; ann_vjust <- 1.3
  } else {
    ann_x <- x_max - 200; ann_y <- fit$intercept + fit$slope * x_max
    ann_hjust <- 0; ann_vjust <- -0.25
  }
  p_combined <- p_combined +
    geom_abline(intercept = fit$intercept, slope = fit$slope,
                color = color, linetype = "dashed", linewidth = 0.6) +
    annotate("text", x = ann_x, y = ann_y, label = label,
             hjust = ann_hjust, vjust = ann_vjust,
             size = 3.2, color = color, family = "mono")
}

# Galloanserae lines + annotations (bottom-left, stacked per locus)
vjust_galloans <- 1.3
for (fit in galloans_fits) {
  p_combined <- p_combined +
    geom_abline(intercept = fit$intercept, slope = fit$slope,
                color = order_palette["Galloanserae"], linetype = "dotted", linewidth = 0.6) +
    annotate("text", x = -Inf, y = -Inf, label = paste0("[", fit$locus, "] ", fit$label),
             hjust = -0.05, vjust = vjust_galloans,
             size = 3.0, color = order_palette["Galloanserae"], family = "mono")
  vjust_galloans <- vjust_galloans + 8
}

p_combined <- p_combined +
  scale_color_manual(values = locus_colors) +
  labs(
    x     = "Total genes",
    y     = "Genes with RSS",
    color = "Locus"
  ) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10))+
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "right"
  )
p_combined

# ── Simple combined plot: locus color only + phylolm ──────────────────────────

make_simple_combined <- function(counts_in_tree, bird_tree, tree_species,
                                  y_var, y_label) {
  run_fit <- function(locus) {
    sub <- filter(counts_in_tree, Locus == locus)
    sub$latin_tree <- gsub(" ", "_", sub$LatinName)
    idx        <- match(tree_species, sub$latin_tree)
    trait_df   <- as.data.frame(sub[idx, ])
    keep       <- !is.na(trait_df$n_total) & !is.na(trait_df[[y_var]])
    trait_df   <- trait_df[keep, ]
    tree_p     <- drop.tip(bird_tree, bird_tree$tip.label[!keep])
    rownames(trait_df) <- tree_p$tip.label
    if (nrow(trait_df) < 5) return(NULL)
    f     <- as.formula(paste(y_var, "~ n_total"))
    model <- phylolm(f, data = trait_df, phy = tree_p, model = "lambda")
    list(locus = locus, intercept = coef(model)[1], slope = coef(model)[2],
         label = make_phylolm_label(model))
  }

  fits  <- Filter(Negate(is.null), lapply(unique(counts_in_tree$Locus), run_fit))
  x_max <- max(counts_in_tree$n_total, na.rm = TRUE)

  p <- ggplot(counts_in_tree, aes(x = n_total, y = .data[[y_var]], color = Locus)) +
    geom_point(size = 2.5, alpha = 0.8)

  for (fit in fits) {
    color <- locus_colors[fit$locus]
    label <- paste0(fit$locus, ":\n", fit$label)
    if (fit$locus == "IGL") {
      ax <- -Inf; ay <- Inf; ah <- -0.05; av <- 1.3
    } else {
      ax <- x_max- 200; ay <- fit$intercept + fit$slope * x_max; ah <- -0.05; av <- 0.5
    }
    p <- p +
      geom_abline(intercept = fit$intercept, slope = fit$slope,
                  color = color, linetype = "dashed", linewidth = 0.6) +
      annotate("text", x = ax, y = ay, label = label,
               hjust = ah, vjust = av, size = 3.2, color = color, family = "mono")
  }

  p +
    scale_color_manual(values = locus_colors) +
    labs(x = "Total genes", y = y_label, color = "Locus") +
    theme_classic() +
    theme(axis.title = element_text(size = 14),
          axis.text  = element_text(size = 10),
          legend.position = "right")
}

p_combined_simple      <- make_simple_combined(counts_in_tree, bird_tree, tree_species,
                                               "n_rss",        "Genes with RSS")+theme(
                                                 legend.position = "none",
                                                 strip.text      = element_text(face = "bold"),
                                                 axis.title      = element_text(size = 14),
                                                 axis.text       = element_text(size = 10)
                                               )
p_combined_productive  <- make_simple_combined(counts_in_tree, bird_tree, tree_species,
                                               "n_productive", "Productive genes")
p_combined_prod_rss    <- make_simple_combined(counts_in_tree, bird_tree, tree_species,
                                               "n_prod_rss",   "Productive genes with RSS")

highlighted_prod_rss <- counts_in_tree %>%
  filter(
    (LatinName == "Gallus gallus"   & Locus == "IGH") |
    (LatinName == "Aythya fuligula" & Locus == "IGH")
  ) %>%
  mutate(sp_color = species_colors[LatinName])

p_combined_prod_rss_birds <- p_combined_prod_rss +
  geom_point(data = highlighted_prod_rss,
             aes(x = n_total, y = n_prod_rss),
             color = highlighted_prod_rss$sp_color, shape = 18, size = 4,
             inherit.aes = FALSE) +
  ggrepel::geom_text_repel(data = highlighted_prod_rss,
             aes(x = n_total, y = n_prod_rss, label = LatinName),
             color = highlighted_prod_rss$sp_color, size = 3, fontface = "italic",
             inherit.aes = FALSE)

highlighted_prod_rss_both <- counts_in_tree %>%
  filter(LatinName %in% c("Gallus gallus", "Aythya fuligula")) %>%
  mutate(sp_color = species_colors[LatinName])

p_combined_prod_rss_birds_both <- p_combined_prod_rss +
  geom_point(data = highlighted_prod_rss_both,
             aes(x = n_total, y = n_prod_rss),
             color = highlighted_prod_rss_both$sp_color, shape = 18, size = 4,
             inherit.aes = FALSE) +
  ggrepel::geom_text_repel(data = highlighted_prod_rss_both,
             aes(x = n_total, y = n_prod_rss, label = LatinName),
             color = highlighted_prod_rss_both$sp_color, size = 3, fontface = "italic",
             inherit.aes = FALSE)

p_combined_simple
p_combined_productive
p_combined_prod_rss
p_combined_prod_rss_birds
p_combined_prod_rss_birds_both

# ── Strand consistency of RSS genes ───────────────────────────────────────────
# Per haplotype × locus: which strand carries more RSS genes ("main strand")?
# Metric: % of RSS genes on the main strand (50% = perfectly split, 100% = all on one strand)

strand_stats <- df %>%
  filter(has_rss) %>%
  group_by(LatinName, Haplotype, Locus, Strand) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(LatinName, Haplotype, Locus) %>%
  summarise(
    n_rss      = sum(n),
    pct_main   = max(n) / sum(n) * 100,
    main_strand = Strand[which.max(n)],
    .groups = "drop"
  ) %>%
  filter(n_rss >= 2)   # need at least 2 RSS genes to have a meaningful ratio

p_strand <- ggplot(strand_stats, aes(x = Locus, y = pct_main, fill = Locus)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(aes(color = Locus), width = 0.08, size = 1.2, alpha = 0.5) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  scale_fill_manual(values  = locus_colors) +
  scale_color_manual(values = locus_colors) +
  scale_y_continuous(limits = c(45, 100), breaks = seq(50, 100, 10),
                     labels = function(x) paste0(x, "%")) +
  annotate("text", x = 0.55, y = 51, label = "50% (no bias)",
           hjust = 0, vjust = -0.3, size = 3, color = "grey50") +
  labs(
    x    = "Locus",
    y    = "RSS genes on main strand (%)",
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 10)
  )

# ── Position of all RSS genes: IGH vs IGL ────────────────────────────────────

all_rss_pos <- df %>%
  group_by(LatinName, Haplotype, Locus) %>%
  mutate(rel_pos = rank(Pos, ties.method = "first") / n()) %>%
  filter(has_rss) %>%
  mutate(folded_pos = pmin(rel_pos, 1 - rel_pos)) %>%
  ungroup()

p_all_rss_pos <- ggplot(all_rss_pos, aes(x = folded_pos, fill = Locus, color = Locus)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.025, alpha = 0.5, position = "identity") +
  geom_density(alpha = 0, linewidth = 0.8) +
  scale_fill_manual(values  = locus_colors) +
  scale_color_manual(values = locus_colors) +
  scale_x_continuous(limits = c(0, 0.5),
                     breaks = seq(0, 0.5, 0.1),
                     labels = function(x) paste0(round(x * 100), "%")) +
  labs(
    x     = "Folded relative position  (0% = at end, 50% = in middle)",
    y     = "Density",
    fill  = "Locus",
    color = "Locus"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 10),
    legend.position = "right"
  )

# ── Position of sole RSS gene in haplotypes with exactly one RSS gene ─────────
# Relative position = rank of RSS gene by Pos among all genes in haplotype × locus.
# Folded: pmin(rel, 1 - rel) so that 0.3 and 0.7 are equivalent (direction unknown).
# 0 = near an end, 0.5 = in the middle.

single_rss <- df %>%
  group_by(LatinName, Haplotype, Locus) %>%
  mutate(
    n_rss_hap = sum(has_rss),
    rel_pos   = rank(Pos, ties.method = "first") / n()
  ) %>%
  filter(n_rss_hap == 1, has_rss) %>%
  mutate(folded_pos = pmin(rel_pos, 1 - rel_pos)) %>%
  ungroup()

p_single_rss <- ggplot(single_rss, aes(x = folded_pos, fill = Locus, color = Locus)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.05, alpha = 0.5, position = "identity") +
  geom_density(alpha = 0, linewidth = 0.8) +
  geom_vline(xintercept = 0.25, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  scale_fill_manual(values  = locus_colors) +
  scale_color_manual(values = locus_colors) +
  scale_x_continuous(limits = c(0, 0.5),
                     breaks = seq(0, 0.5, 0.1),
                     labels = function(x) paste0(round(x * 100), "%")) +
  annotate("text", x = 0.255, y = Inf, label = "middle of\nend–middle range",
           hjust = -0.05, vjust = 1.4, size = 2.8, color = "grey50") +
  labs(
    x     = "Folded relative position  (0% = at end, 50% = in middle)",
    y     = "Density",
    fill  = "Locus",
    color = "Locus"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 10),
    legend.position = "right"
  )

# ── Position of the sole productive+RSS gene in haplotypes with multiple RSS ──
# Cases: >1 RSS gene total, but exactly 1 that is both productive and has RSS.
# Same folded relative position metric as above.

single_productive_rss <- df %>%
  mutate(productive = tolower(trimws(as.character(Productive))) %in% c("true", "1", "yes"),
         has_prod_rss = productive & has_rss) %>%
  group_by(LatinName, Haplotype, Locus) %>%
  mutate(
    n_rss_hap      = sum(has_rss),
    n_prod_rss_hap = sum(has_prod_rss),
    rel_pos        = rank(Pos, ties.method = "first") / n()
  ) %>%
  filter(n_rss_hap > 1, n_prod_rss_hap == 1, has_prod_rss) %>%
  mutate(folded_pos = pmin(rel_pos, 1 - rel_pos)) %>%
  ungroup()


p_single_productive_rss <- ggplot(single_productive_rss,
                                   aes(x = folded_pos, fill = Locus, color = Locus)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.05, alpha = 0.5, position = "identity") +
  geom_density(alpha = 0, linewidth = 0.8) +
  geom_vline(xintercept = 0.25, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  scale_fill_manual(values  = locus_colors) +
  scale_color_manual(values = locus_colors) +
  scale_x_continuous(limits = c(0, 0.5),
                     breaks = seq(0, 0.5, 0.1),
                     labels = function(x) paste0(round(x * 100), "%")) +
  labs(
    title = "Haplotypes with multiple RSS genes but only one productive+RSS",
    x     = "Folded relative position  (0% = at end, 50% = in middle)",
    y     = "Density",
    fill  = "Locus",
    color = "Locus"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 11),
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 10),
    legend.position = "right"
  )

# ── Save ──────────────────────────────────────────────────────────────────────

print(p_combined)
print(p_all_rss_pos)
print(p_single_rss)
# ── Variance in RSS gene counts: IGH vs IGL ───────────────────────────────────

rss_var <- counts %>%
  group_by(Locus) %>%
  summarise(
    mean_rss = mean(n_rss, na.rm = TRUE),
    var_rss  = var(n_rss,  na.rm = TRUE),
    sd_rss   = sd(n_rss,   na.rm = TRUE),
    n        = n(),
    .groups  = "drop"
  )

# F-test for equality of variances
igh_rss <- counts$n_rss[counts$Locus == "IGH"]
igl_rss <- counts$n_rss[counts$Locus == "IGL"]
ftest   <- var.test(igh_rss, igl_rss)

var_label <- sprintf(
  "F = %.2f\np = %.3g",
  ftest$statistic, ftest$p.value
)

var_annotations <- rss_var %>%
  mutate(label = sprintf("var = %.1f\nsd  = %.1f", var_rss, sd_rss))

p_rss_var <- ggplot(counts, aes(x = n_rss, fill = Locus)) +
  geom_histogram(binwidth = 1, alpha = 0.8, color = "white", linewidth = 0.3) +
  geom_text(data = var_annotations,
            aes(x = Inf, y = Inf,
                label = sprintf("mean = %.1f\nvar  = %.1f\nsd   = %.1f", mean_rss, var_rss, sd_rss)),
            hjust = 1.1, vjust = 1.4, size = 3.2, family = "mono", inherit.aes = FALSE) +
  annotate("text", x = -Inf, y = Inf,
           label = var_label,
           hjust = -0.1, vjust = 1.4, size = 3.2, color = "grey30", family = "mono") +
  facet_wrap(~Locus, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = locus_colors) +
  labs(
    x = "RSS genes per species (mean across haplotypes)",
    y = "Count"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text      = element_text(face = "bold"),
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 10)
  )

print(p_rss_var)
print(p_single_productive_rss)

# ── Positional bias of pseudogenes relative to most 3' RSS+ gene ──────────────
# H: in IGL, pseudogenes accumulate downstream of the most 3' RSS+ gene.
# Anchor = gene with highest Pos among has_rss=TRUE genes in each haplotype.
# Statistic: proportion of pseudogenes (across all haplotypes) that lie downstream
# of their haplotype's anchor. Null via permutation of productive labels within
# each haplotype (keeping RSS status fixed so the anchor is unchanged).

df_bias <- df %>%
  mutate(is_pseudo = !(tolower(trimws(as.character(Productive))) %in% c("true", "1", "yes")))

get_downstream_prop <- function(data) {
  data %>%
    group_by(LatinName, Haplotype, Locus) %>%
    filter(any(has_rss)) %>%
    mutate(anchor_pos = max(Pos[has_rss])) %>%
    filter(!(has_rss & Pos == anchor_pos)) %>%   # exclude the anchor gene itself
    mutate(downstream = Pos > anchor_pos) %>%
    summarise(
      n_pseudo_down  = sum(is_pseudo & downstream),
      n_pseudo_total = sum(is_pseudo),
      .groups = "drop"
    ) %>%
    filter(n_pseudo_total > 0) %>%
    group_by(Locus) %>%
    summarise(prop = sum(n_pseudo_down) / sum(n_pseudo_total), .groups = "drop")
}

obs_stats <- get_downstream_prop(df_bias)

set.seed(42)
n_perm <- 1000

perm_results <- bind_rows(
  lapply(seq_len(n_perm), function(i) {
    df_bias %>%
      group_by(LatinName, Haplotype, Locus) %>%
      mutate(is_pseudo = sample(is_pseudo)) %>%
      ungroup() %>%
      get_downstream_prop() %>%
      mutate(perm = i)
  })
)

# One-sided p-value: proportion of permutations >= observed (testing for downstream excess)
perm_summary <- perm_results %>%
  left_join(obs_stats %>% rename(obs_prop = prop), by = "Locus") %>%
  group_by(Locus) %>%
  summarise(
    obs_prop = first(obs_prop),
    p_val    = mean(prop >= obs_prop),
    .groups  = "drop"
  ) %>%
  mutate(label = sprintf("obs = %.3f\np   = %.3f", obs_prop, p_val))

p_bias <- ggplot(perm_results, aes(x = prop, fill = Locus)) +
  geom_histogram(binwidth = 0.01, alpha = 0.8, color = "white", linewidth = 0.3) +
  geom_vline(data = obs_stats, aes(xintercept = prop, color = Locus),
             linewidth = 1) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_text(data = perm_summary,
            aes(x = obs_prop, y = Inf, label = label, color = Locus),
            hjust = -0.1, vjust = 1.4, size = 3.2, family = "mono",
            inherit.aes = FALSE) +
  facet_wrap(~Locus, ncol = 1, scales = "free_y") +
  scale_fill_manual(values  = locus_colors) +
  scale_color_manual(values = locus_colors) +
  labs(
    x = "Proportion of pseudogenes downstream of most 3′ RSS+ gene",
    y = "Permutation count"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text      = element_text(face = "bold"),
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 10)
  )

print(p_bias)

# ── IGL RSS gene positions: single vs multiple productive RSS haplotypes ───────

make_prod_rss_pos_plot <- function(locus) {
  bins <- seq(0, 0.5, by = 0.025)

  pos_data <- df %>%
    filter(Locus == locus) %>%
    mutate(productive = tolower(trimws(as.character(Productive))) %in% c("true", "1", "yes")) %>%
    group_by(LatinName, Haplotype) %>%
    mutate(
      rel_pos    = rank(Pos, ties.method = "first") / n(),
      n_prod_rss = sum(has_rss & productive)
    ) %>%
    filter(has_rss) %>%
    mutate(
      folded_pos = pmin(rel_pos, 1 - rel_pos),
      group      = ifelse(n_prod_rss > 1, "Multiple productive RSS", "Single productive RSS")
    ) %>%
    ungroup() %>%
    filter(folded_pos >= 0, folded_pos <= 0.5) %>%
    mutate(bin = cut(folded_pos, breaks = bins, include.lowest = TRUE, right = FALSE)) %>%
    group_by(group, bin) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(group) %>%
    mutate(pct = count / sum(count) * 100) %>%
    ungroup() %>%
    mutate(bin_mid = bins[as.integer(bin)] + 0.0125,
           group   = factor(group, levels = c("Multiple productive RSS", "Single productive RSS")))

  group_colors <- c("Multiple productive RSS" = "#4E79A7", "Single productive RSS" = "#F28E2B")

  ggplot(pos_data, aes(x = bin_mid, y = pct, fill = group, color = group)) +
    geom_col(alpha = 0.5, position = "identity", width = 0.023) +
    scale_x_continuous(limits = c(0, 0.5),
                       breaks = seq(0, 0.5, 0.1),
                       labels = function(x) paste0(round(x * 100), "%")) +
    scale_fill_manual(values  = group_colors, name = NULL) +
    scale_color_manual(values = group_colors, name = NULL) +
    labs(
      title = locus,
      x     = "Folded relative position  (0% = at end, 50% = in middle)",
      y     = "% of genes in category"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold"),
      axis.title      = element_text(size = 12),
      axis.text       = element_text(size = 10),
      legend.position = "top"
    )
}

p_igl_prod_rss_pos <- make_prod_rss_pos_plot("IGL")
p_igh_prod_rss_pos <- make_prod_rss_pos_plot("IGH")

print(p_igl_prod_rss_pos)
print(p_igh_prod_rss_pos)

# Old-style density versions (geom_histogram + after_stat(density), binwidth 0.05)
make_prod_rss_pos_density <- function(locus) {
  group_colors <- c("Multiple productive RSS" = "#4E79A7", "Single productive RSS" = "#F28E2B")

  pos_data <- df %>%
    filter(Locus == locus) %>%
    mutate(productive = tolower(trimws(as.character(Productive))) %in% c("true", "1", "yes")) %>%
    group_by(LatinName, Haplotype) %>%
    mutate(
      rel_pos    = rank(Pos, ties.method = "first") / n(),
      n_prod_rss = sum(has_rss & productive)
    ) %>%
    filter(has_rss) %>%
    mutate(
      folded_pos = pmin(rel_pos, 1 - rel_pos),
      group      = ifelse(n_prod_rss > 1, "Multiple productive RSS", "Single productive RSS")
    ) %>%
    ungroup()

  ggplot(pos_data, aes(x = folded_pos, fill = group, color = group)) +
    geom_histogram(aes(y = after_stat(density)),
                   binwidth = 0.025, alpha = 0.5, position = "identity") +
    scale_x_continuous(limits = c(0, 0.5),
                       breaks = seq(0, 0.5, 0.1),
                       labels = function(x) paste0(round(x * 100), "%")) +
    scale_fill_manual(values  = group_colors, name = NULL) +
    scale_color_manual(values = group_colors, name = NULL) +
    labs(
      title = locus,
      x     = "Folded relative position  (0% = at end, 50% = in middle)",
      y     = "Density"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold"),
      axis.title      = element_text(size = 12),
      axis.text       = element_text(size = 10),
      legend.position = "top"
    )
}

p_igl_prod_rss_pos_density <- make_prod_rss_pos_density("IGL")
p_igh_prod_rss_pos_density <- make_prod_rss_pos_density("IGH")

print(p_igl_prod_rss_pos_density)
print(p_igh_prod_rss_pos_density)

# Unfolded version (0–100% relative position)
make_prod_rss_pos_unfolded <- function(locus,x_label=TRUE) {
  
  locus_color<-if(locus=="IGH"){"#87b4dc"}else{"#638E6E"}
  group_colors <- c("Multiple productive RSS" = "lightgrey", "Single productive RSS" = locus_color)

  pos_data <- df %>%
    filter(Locus == locus) %>%
    mutate(productive = tolower(trimws(as.character(Productive))) %in% c("true", "1", "yes")) %>%
    group_by(LatinName, Haplotype) %>%
    mutate(
      rel_pos    = rank(Pos, ties.method = "first") / n(),
      n_prod_rss = sum(has_rss & productive)
    ) %>%
    filter(has_rss) %>%
    mutate(group = ifelse(n_prod_rss > 1, "Multiple productive RSS", "Single productive RSS")) %>%
    ungroup()

  ggplot(pos_data, aes(x = rel_pos, fill = group, color = group)) +
    geom_histogram(aes(y = after_stat(density)),
                   binwidth = 0.05, alpha = 0.5, position = "identity") +
    scale_x_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, 0.2),
                       labels = function(x) paste0(round(x * 100), "%")) +
    scale_fill_manual(values  = group_colors, name = NULL) +
    scale_color_manual(values = group_colors, name = NULL) +
    labs(
      x     = if(x_label){"Relative position  (0% = start, 100% = end)"}else{""},
      y     = "Density"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold"),
      axis.title      = element_text(size = 12),
      axis.text       = element_text(size = 10),
      legend.position = "none"
    )
}

p_igl_prod_rss_pos_unfolded <- make_prod_rss_pos_unfolded("IGL")
p_igh_prod_rss_pos_unfolded <- make_prod_rss_pos_unfolded("IGH",FALSE)

print(p_igl_prod_rss_pos_unfolded)
print(p_igh_prod_rss_pos_unfolded)

# Unfolded + percentage y-axis
make_prod_rss_pos_unfolded_pct <- function(locus) {
  bins         <- seq(0, 1, by = 0.05)
  group_colors <- c("Multiple productive RSS" = "#4E79A7", "Single productive RSS" = "#F28E2B")

  pos_data <- df %>%
    filter(Locus == locus) %>%
    mutate(productive = tolower(trimws(as.character(Productive))) %in% c("true", "1", "yes")) %>%
    group_by(LatinName, Haplotype) %>%
    mutate(
      rel_pos    = rank(Pos, ties.method = "first") / n(),
      n_prod_rss = sum(has_rss & productive)
    ) %>%
    filter(has_rss) %>%
    mutate(group = ifelse(n_prod_rss > 1, "Multiple productive RSS", "Single productive RSS")) %>%
    ungroup() %>%
    mutate(bin = cut(rel_pos, breaks = bins, include.lowest = TRUE, right = FALSE)) %>%
    group_by(group, bin) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(group) %>%
    mutate(pct = count / sum(count) * 100) %>%
    ungroup() %>%
    mutate(bin_mid = bins[as.integer(bin)] + 0.025,
           group   = factor(group, levels = c("Multiple productive RSS", "Single productive RSS")))

  ggplot(pos_data, aes(x = bin_mid, y = pct, fill = group, color = group)) +
    geom_col(alpha = 0.5, position = "identity", width = 0.05) +
    scale_x_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, 0.2),
                       labels = function(x) paste0(round(x * 100), "%")) +
    scale_fill_manual(values  = group_colors, name = NULL) +
    scale_color_manual(values = group_colors, name = NULL) +
    labs(
      title = locus,
      x     = "Relative position  (0% = start, 100% = end)",
      y     = "% of genes in category"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold"),
      axis.title      = element_text(size = 12),
      axis.text       = element_text(size = 10),
      legend.position = "top"
    )
}

p_igl_prod_rss_pos_unfolded_pct <- make_prod_rss_pos_unfolded_pct("IGL")
p_igh_prod_rss_pos_unfolded_pct <- make_prod_rss_pos_unfolded_pct("IGH")

print(p_igl_prod_rss_pos_unfolded_pct)
print(p_igh_prod_rss_pos_unfolded_pct)



# ── MinDir vs gene counts (IGH + IGL pooled) ─────────────────────────────────
# Note: lm used instead of phylolm because pooling both loci gives two rows per
# species, which cannot be matched to a single-tip-per-species tree.
library(data.table)

mindir_raw <- fread(file.path(input_dir, "all_species_stats_pruned_12052025.csv"))
mindir_raw <- mindir_raw[mindir_raw$IGH_AnnotationLevel < 2, ]
mindir_raw <- mindir_raw[mindir_raw$VertClass == "birds", ]

mindir_raw$latin_tree <- tolower(gsub(" ", "_", mindir_raw$LatinName))

vgp_order_map <- vgp_table %>%
  select(LatinName, Order) %>%
  distinct() %>%
  mutate(latin_tree = tolower(gsub(" ", "_", LatinName)))

mindir_birds <- merge(
  as.data.frame(mindir_raw) %>% select(-LatinName),
  vgp_order_map %>% select(latin_tree, LatinName, Order),
  by = "latin_tree", all.x = TRUE
)

# Pivot IGH_MinDir and IGL_MinDir to long format
mindir_long <- mindir_birds %>%
  select(LatinName, latin_tree, IGH_MinDir, IGL_MinDir) %>%
  pivot_longer(
    cols      = c(IGH_MinDir, IGL_MinDir),
    names_to  = "Locus",
    values_to = "MinDir"
  ) %>%
  mutate(Locus = sub("_MinDir", "", Locus))

# Join with counts_in_tree on LatinName × Locus
mindir_pooled <- counts_in_tree %>%
  left_join(mindir_long, by = c("LatinName", "Locus")) %>%
  filter(!is.na(MinDir))

make_lm_label <- function(model) {
  s     <- summary(model)
  slope <- coef(model)[2]
  p_val <- s$coefficients[2, "Pr(>|t|)"]
  r2    <- s$r.squared
  stars <- ifelse(p_val < 0.001, "***",
           ifelse(p_val < 0.01,  "**",
           ifelse(p_val < 0.05,  "*", "ns")))
  sprintf("β  = %-8s\nR² = %-8.3f",
          paste0(round(slope, 3), stars), r2)
}

make_mindir_pooled_plot <- function(y_var, y_label) {
  dat <- mindir_pooled %>% filter(!is.na(.data[[y_var]]))
  mod <- lm(as.formula(paste(y_var, "~ MinDir")), data = dat)
  x_rng <- range(dat$MinDir, na.rm = TRUE)
  ggplot(dat, aes(x = MinDir, y = .data[[y_var]], color = Locus)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_abline(intercept = coef(mod)[1], slope = coef(mod)[2],
                linetype = "dashed", linewidth = 0.6, color = "black") +
    annotate("text",
             x = x_rng[1], y = Inf,
             label  = make_lm_label(mod),
             hjust  = -0.05, vjust = 1.3,
             size   = 3.2, family = "mono") +
    scale_color_manual(values = locus_colors) +
    labs(x = "MinDir (# genes in minority direction)", y = y_label) +
    theme_classic(base_size = 12) +
    theme(axis.title = element_text(size = 12),
          axis.text  = element_text(size = 10))
}

p_mindir_rss   <- make_mindir_pooled_plot("n_rss",   "# genes with RSS")
p_mindir_total <- make_mindir_pooled_plot("n_total", "# total genes")

print(p_mindir_rss)
print(p_mindir_total)

# ── n_total vs n_rss / n_prod_rss colored by MinDir ──────────────────────────

mindir_counts <- counts_in_tree %>%
  left_join(mindir_long, by = c("LatinName", "Locus")) %>%
  filter(!is.na(MinDir))

p_mindir_color_rss <- ggplot(mindir_counts,
                             aes(x = n_total, y = n_rss, color = MinDir)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_viridis_c(name = "MinDir", option = "C") +
  labs(x = "Total genes (mean per species)",
       y = "Genes with RSS (mean per species)") +
  theme_classic(base_size = 12) +
  theme(axis.title  = element_text(size = 12),
        axis.text   = element_text(size = 10),
        strip.text  = element_text(face = "bold"))

p_mindir_color_prod_rss <- ggplot(mindir_counts,
                                  aes(x = n_total, y = n_prod_rss, color = MinDir)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_viridis_c(name = "MinDir", option = "C") +
  labs(x = "Total genes (mean per species)",
       y = "Productive genes with RSS (mean per species)") +
  theme_classic(base_size = 12) +
  theme(axis.title  = element_text(size = 12),
        axis.text   = element_text(size = 10),
        strip.text  = element_text(face = "bold"))

print(p_mindir_color_rss)
print(p_mindir_color_prod_rss)

# ── Terminal gap analysis: isolated genes at locus ends ───────────────────────
# For each haplotype (≥4 genes), compute gaps between consecutive genes sorted
# by position. Normalize each gap by the haplotype's median gap. Compare the
# distribution of first, last, and internal normalized gaps to see if terminal
# genes tend to be physically separated from the main cluster.

gap_df <- df %>%
  group_by(LatinName, Haplotype, Locus) %>%
  arrange(Pos, .by_group = TRUE) %>%
  filter(n() >= 4) %>%
  mutate(
    gene_rank  = row_number(),
    n_genes    = n(),
    gap        = Pos - lag(Pos)
  ) %>%
  filter(!is.na(gap)) %>%
  mutate(
    gap_type   = case_when(
      gene_rank == 1         ~ "First",
      gene_rank == 2         ~ "First",
      gene_rank == n_genes   ~ "Last",
      gene_rank == n_genes-1   ~ "Last",
      TRUE                   ~ "Internal"
    ),
    norm_gap   = gap / median(gap)
  ) %>%
  ungroup() %>%
  mutate(gap_type = factor(gap_type, levels = c("First", "Internal", "Last")))

# Plot 1: distribution of normalised gaps by type, faceted by locus
p_terminal_gaps <- ggplot(gap_df, aes(x = gap_type, y = norm_gap, fill = Locus)) +
  geom_violin(alpha = 0.6, trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.8,
               position = position_dodge(0)) +
  facet_wrap(~ Locus) +
  scale_fill_manual(values = locus_colors) +
  scale_y_log10() +
  labs(
    x = "Gap position",
    y = "Normalised gap size (log scale, relative to haplotype median)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text      = element_text(face = "bold"),
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 10)
  )



print(p_terminal_gaps)


final_figure<-(contig_length_p+mind_dir_p)/(p_combined_simple+(p_igh_prod_rss_pos_unfolded/p_igl_prod_rss_pos_unfolded))
final_figure
(contig_length_p+mind_dir_p)
(p_combined_simple+(p_igh_prod_rss_pos_unfolded/p_igl_prod_rss_pos_unfolded))

