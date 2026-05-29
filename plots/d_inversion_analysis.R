suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(ape)
  library(ggtree)
  library(ggtreeExtra)
})

# ── INPUTS ────────────────────────────────────────────────────────────────────
INPUT_DIR   <- "/local/storage/kav67/clean_birds/"
D_INV_FILE  <- file.path(INPUT_DIR, "D_inversions.tsv")
VGP_TABLE   <- file.path(INPUT_DIR, "IGH_VGP_table.tsv")
TREE_FILE   <- file.path(INPUT_DIR, "vgp_birds.nwk")

# ── READ DATA ─────────────────────────────────────────────────────────────────
d_inv <- read_tsv(D_INV_FILE, show_col_types = FALSE)

vgp   <- read_tsv(VGP_TABLE, show_col_types = FALSE) %>%
  select(Haplotype, LatinName) %>%
  distinct()

d_inv <- d_inv %>%
  left_join(vgp, by = "Haplotype") %>%
  # Average across haplotypes within species so pri/alt don't double-count
  group_by(Order, Species, LatinName) %>%
  summarise(
    n_d_genes        = mean(n_d_genes,        na.rm = TRUE),
    n_on_inversion   = mean(n_on_inversion,   na.rm = TRUE),
    frac_on_inversion = mean(frac_on_inversion, na.rm = TRUE),
    n_haplotypes     = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(frac_on_inversion))

message("Loaded ", nrow(d_inv), " species")

# ── PLOT 1: Distribution of D-gene inversion fraction by Order ────────────────
p_by_order <- ggplot(d_inv,
                     aes(x = reorder(Order, frac_on_inversion, median),
                         y = frac_on_inversion,
                         fill = Order)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.8, linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5, color = "black") +
  scale_y_continuous(labels = scales::percent_format(),
                     limits = c(0, 1)) +
  labs(
    x     = "Order",
    y     = "Fraction of D genes on inversions",
    title = "D genes on inversions by taxonomic order"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold"),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


# ── PLOT 3: Stacked bar — absolute counts, sorted by frac_on_inversion ────────
bar_df <- d_inv %>%
  arrange(desc(frac_on_inversion)) %>%
  mutate(label = factor(LatinName, levels = unique(LatinName)),
         n_off_inversion = n_d_genes - n_on_inversion) %>%
  select(label, Order, n_on_inversion, n_off_inversion) %>%
  pivot_longer(c(n_on_inversion, n_off_inversion),
               names_to = "group", values_to = "count") %>%
  mutate(group = recode(group,
                        n_on_inversion  = "On inversion",
                        n_off_inversion = "Off inversion"))

p_bar <- ggplot(bar_df,
                aes(x = label, y = count, fill = group)) +
  geom_col(position = "stack", width = 0.8) +
  scale_fill_manual(values = c("On inversion" = "#c0392b",
                               "Off inversion" = "#bdc3c7"),
                    name = NULL) +
  labs(
    x     = NULL,
    y     = "Number of D genes",
    title = "D genes on vs off inversions per species"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold"),
    axis.text.x     = element_text(angle = 55, hjust = 1, size = 7,
                                   face = "italic"),
    legend.position = "top"
  )

# ── PLOT 4: Overall summary histogram ─────────────────────────────────────────
overall_mean <- mean(d_inv$frac_on_inversion, na.rm = TRUE)
overall_med  <- median(d_inv$frac_on_inversion, na.rm = TRUE)

p_hist <- ggplot(d_inv, aes(x = frac_on_inversion)) +
  geom_histogram(binwidth = 0.05, fill = "#87b4dc",
                 color = "white", linewidth = 0.3) +
  geom_vline(xintercept = overall_mean, linetype = "dashed",
             color = "#c0392b", linewidth = 0.7) +
  geom_vline(xintercept = overall_med, linetype = "dotted",
             color = "#2c3e50", linewidth = 0.7) +
  annotate("text", x = overall_mean, y = Inf,
           label = sprintf("mean = %.2f", overall_mean),
           hjust = -0.1, vjust = 1.5, color = "#c0392b", size = 3.5) +
  annotate("text", x = overall_med, y = Inf,
           label = sprintf("median = %.2f", overall_med),
           hjust = -0.01, vjust = 1.5, color = "#2c3e50", size = 3.5) +
  scale_x_continuous(labels = scales::percent_format()) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x     = "Fraction of D genes on inversions",
    y     = "Number of species",
    title = "Distribution across all species"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# ── PLOT 5: D-gene inversion fraction mapped onto the VGP phylogenetic tree ───
bird_tree <- read.tree(TREE_FILE)

# Tree tip labels encode Latin names as e.g. 'Prefix"Genus_species"'
tree_species <- gsub('.*"', '', bird_tree$tip.label)
tree_species <- gsub('"',   '', tree_species)

# Match d_inv to tree; LatinName uses spaces, tree uses underscores
d_inv_tree <- d_inv %>%
  mutate(latin_tree = gsub(" ", "_", LatinName)) %>%
  filter(latin_tree %in% tree_species)

# Prune tree to species that have D-inversion data
tree_pruned <- drop.tip(
  bird_tree,
  bird_tree$tip.label[!tree_species %in% d_inv_tree$latin_tree]
)

# Rebuild tip-species mapping on pruned tree (tip labels still contain encoded names)
pruned_species <- gsub('.*"', '', tree_pruned$tip.label)
pruned_species <- gsub('"',   '', pruned_species)

# Add tip label column for ggtree join
d_inv_tree <- d_inv_tree %>%
  mutate(tip_label = tree_pruned$tip.label[match(latin_tree, pruned_species)])

tree_height <- max(10, length(tree_pruned$tip.label) * 0.18)

p_tree_base <- ggtree(tree_pruned) +
  ylim(0, length(tree_pruned$tip.label) + 4) +
  theme_tree()

# Bar track: frac_on_inversion
# MRCA node per Order (≥2 tips needed for a clade node)
order_nodes <- d_inv_tree %>%
  filter(!is.na(Order), Order != "MiscBirds") %>%
  group_by(Order) %>%
  summarise(tips = list(tip_label), n = n(), .groups = "drop") %>%
  filter(n >= 2) %>%
  mutate(node = sapply(tips, function(t) {
    t_in <- t[t %in% tree_pruned$tip.label]
    if (length(t_in) < 2) NA_integer_
    else getMRCA(tree_pruned, t_in)
  })) %>%
  filter(!is.na(node))

p_tree_plot <- p_tree_base +
  geom_cladelab(
    data    = order_nodes,
    mapping = aes(node = node, label = Order),
    align   = TRUE,
    offset  = 0.002,
    fontsize = 2.5,
    hjust   = 0
  ) +
  geom_fruit(
    data    = d_inv_tree,
    geom    = geom_bar,
    mapping = aes(y = tip_label, x = frac_on_inversion, fill = frac_on_inversion),
    stat    = "identity",
    orientation = "y",
    offset  = 0.2,
    width   = 0.7
  ) +
  scale_fill_viridis_c(
    option = "magma", direction = -1,
    name   = "Fraction D genes\non inversions",
    labels = scales::percent_format()
  ) +
  labs(title = "D-gene inversion fraction across the VGP bird phylogeny") +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    legend.position = "right"
  )

# ── PRINT ─────────────────────────────────────────────────────────────────────
print(p_by_order)
print(p_bar)
print(p_hist)
print(p_tree_plot)


# ── SUMMARY STATS ─────────────────────────────────────────────────────────────
cat("\n── Overall ──────────────────────────────────────────\n")
cat(sprintf("  Species with any D genes on inversions : %d / %d\n",
            sum(d_inv$n_on_inversion > 0), nrow(d_inv)))
cat(sprintf("  Mean fraction on inversions            : %.3f\n", overall_mean))
cat(sprintf("  Median fraction on inversions          : %.3f\n", overall_med))

cat("\n── By Order ─────────────────────────────────────────\n")
d_inv %>%
  group_by(Order) %>%
  summarise(
    n_species     = n(),
    mean_frac     = mean(frac_on_inversion),
    median_frac   = median(frac_on_inversion),
    mean_n_d      = mean(n_d_genes),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_frac)) %>%
  print(n = Inf)
