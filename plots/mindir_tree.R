library(data.table)
library(dplyr)
library(ggplot2)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(viridis)
library(tidyr)

# ── Load data ──────────────────────────────────────────────────────────────────

# all_species_data <- fread("/local/storage/kav67/clean_birds/all_species_stats_pruned_12052025.csv")
# all_species_data$bird <- FALSE
# all_species_data[all_species_data$VertClass == "birds", ]$bird <- TRUE
# all_species_data <- all_species_data[all_species_data$IGH_AnnotationLevel < 2, ]
# all_species_data <- all_species_data[all_species_data$VertClass %in% c("reptiles", "mammals", "birds"), ]
# bird_only <- all_species_data[all_species_data$VertClass == "birds", ]

gene_list <- fread("/local/storage/kav67/clean_birds/gene_list-igl_h2_n6.csv")
gene_list[, c("GrpOrder", "Species", "Haplotype") := tstrsplit(Source, "/", fixed = TRUE)]

main_strand_frac <- function(strands) {
  tbl <- table(strands)
  max(tbl) / sum(tbl)
}

bird_only <- gene_list %>%
  group_by(GrpOrder, Species, Haplotype, Locus, Contig) %>%
  summarise(MinDir = main_strand_frac(Strand), .groups = "drop") %>%
  group_by(GrpOrder, Species, Haplotype, Locus) %>%
  summarise(MinDir = mean(MinDir), .groups = "drop") %>%
  group_by(Species, Locus) %>%
  summarise(MinDir = mean(MinDir), .groups = "drop") %>%
  pivot_wider(names_from = Locus, values_from = MinDir, names_glue = "{Locus}_MinDir") %>%
  mutate(LatinName = gsub("_", " ", Species),
         latin_tree = tolower(Species)) %>%
  as.data.frame()

bird_tree_pruned <- read.tree("/local/storage/kav67/clean_birds/vgp_birds.nwk")
tree_species <- gsub('.*"', '', bird_tree_pruned$tip.label)
tree_species <- gsub('"', '', tree_species)
bird_tree_pruned$species <- tolower(tree_species)

# Join Order and proper LatinName from VGP table
vgp_table <- fread("/local/storage/kav67/clean_birds/IGH_VGP_table.tsv")
order_map  <- vgp_table %>%
  select(LatinName, Order,Species) %>%
  distinct() %>%
  mutate(latin_tree = tolower(gsub(" ", "_", Species)))

# ── Match to tree ──────────────────────────────────────────────────────────────

# Join on latin_tree (lowercase underscore) — also brings in proper LatinName
bird_only <- merge(as.data.frame(bird_only %>% select(!LatinName)),
                   order_map %>% select(latin_tree, Order, LatinName),
                   by = "latin_tree", all.x = TRUE)

idx      <- match(tree_species, gsub(" ","_",bird_only$LatinName))
tip_meta <- bird_only[idx, ]

keep         <- !is.na(tip_meta$LatinName)
tip_meta     <- tip_meta[keep, ]
tree_pruned  <- drop.tip(bird_tree_pruned, bird_tree_pruned$tip.label[!keep])

# Rename tip labels to LatinName so geom_fruit can match on y = LatinName
latin_map              <- setNames(tip_meta$LatinName, tree_pruned$tip.label)
tree_pruned$tip.label  <- latin_map[tree_pruned$tip.label]

# ── Order MRCA nodes for highlighting ─────────────────────────────────────────

order_nodes <- tip_meta %>%
  filter(!is.na(Order), Order != "MiscBirds") %>%
  group_by(Order) %>%
  summarise(tips = list(LatinName), n = n(), .groups = "drop") %>%
  filter(n >= 2) %>%
  mutate(node = sapply(tips, function(t) {
    t_in <- t[t %in% tree_pruned$tip.label]
    if (length(t_in) < 2) NA_integer_
    else getMRCA(tree_pruned, t_in)
  })) %>%
  filter(!is.na(node))

order_colors <- setNames(
  scales::hue_pal()(nrow(order_nodes)),
  order_nodes$Order
)

# ── Base tree ─────────────────────────────────────────────────────────────────

p <- ggtree(tree_pruned, layout = "rectangular") +
  geom_tiplab(size = 1.8, fontface = "italic", offset = 0.001)

# Order highlighting
for (i in seq_len(nrow(order_nodes))) {
  p <- p + geom_hilight(
    node  = order_nodes$node[i],
    fill  = order_colors[order_nodes$Order[i]],
    alpha = 0.25
  )
}

p <- p + geom_cladelab(
  data    = order_nodes,
  mapping = aes(node = node, label = Order),
  align   = TRUE,
  offset  = 0.002,
  fontsize = 2.5,
  hjust   = 0
)

# ── IGH_MinDir panel ──────────────────────────────────────────────────────────

p <- p + new_scale_fill()

p <- p + geom_fruit(
  data    = tip_meta,
  geom    = geom_bar,
  mapping = aes(y = LatinName, x = IGH_MinDir, fill = IGH_MinDir),
  stat    = "identity",
  orientation = "y",
  offset  = 0.15,
  pwidth  = 0.25,
  axis.params = list(axis = "x", text.size = 2, title = "IGH MinDir",
                     title.size = 3, title.height = 0.02)
) +
  scale_fill_viridis_c(name = "IGH\nMinDir", option = "C", na.value = "grey90",
                       limits = c(0, 1))

# ── IGL_MinDir panel ──────────────────────────────────────────────────────────

p <- p + new_scale_fill()

p <- p + geom_fruit(
  data    = tip_meta,
  geom    = geom_bar,
  mapping = aes(y = LatinName, x = IGL_MinDir, fill = IGL_MinDir),
  stat    = "identity",
  orientation = "y",
  offset  = 0.05,
  pwidth  = 0.25,
  axis.params = list(axis = "x", text.size = 2, title = "IGL MinDir",
                     title.size = 3, title.height = 0.02)
) +
  scale_fill_viridis_c(name = "IGL\nMinDir", option = "D", na.value = "grey90",
                       limits = c(0, 1))

# ── Save (0–1 scale) ──────────────────────────────────────────────────────────
p

# ── Version 2: scale from 0.5–1.0 ────────────────────────────────────────────

p2 <- ggtree(tree_pruned, layout = "rectangular") #+
  #geom_tiplab(size = 1.8, fontface = "italic", offset = 0.001)

for (i in seq_len(nrow(order_nodes))) {
  p2 <- p2 + geom_hilight(
    node  = order_nodes$node[i],
    fill  = order_colors[order_nodes$Order[i]],
    alpha = 0.25
  )
}

p2 <- p2 + geom_cladelab(
  data    = order_nodes,
  mapping = aes(node = node, label = Order),
  align   = TRUE,
  offset  = 0.002,
  fontsize = 2.5,
  hjust   = 0
)

butterfly_data <- tip_meta %>%
  select(LatinName, IGH_MinDir, IGL_MinDir) %>%
  pivot_longer(c(IGH_MinDir, IGL_MinDir), names_to = "Locus", values_to = "MinDir") %>%
  mutate(x_pos = if_else(Locus == "IGH_MinDir", -MinDir, MinDir),
         bar_alpha = 1)

# Transparent ghost bar at x = -1 using a real tree species so ggtreeExtra
# includes it in scale computation, extending the IGH side to -1
ghost_sp <- tip_meta$LatinName[which.max(tip_meta$IGH_MinDir)]
butterfly_data <- butterfly_data %>%
  add_row(LatinName = ghost_sp, Locus = "IGH_MinDir", MinDir = 1,
          x_pos = -1, bar_alpha = 0)

p2 <- p2 + new_scale_fill()

p2 <- p2 + geom_fruit(
  data    = butterfly_data,
  geom    = geom_bar,
  mapping = aes(y = LatinName, x = x_pos, fill = MinDir, alpha = bar_alpha),
  stat    = "identity",
  orientation = "y",
  offset  = 0.7,
  pwidth  = 0.5,
  axis.params = list(axis = "x", text.size = 2, title = "                        ← IGH | IGL →",
                     title.size = 3, title.height = 0.02, limits = c(-1, 1.1))
) +
  scale_fill_viridis_c(name = "MinDir", option = "D", na.value = "grey90",direction=-1,
                       limits = c(0.5, 1)) +
  scale_alpha_identity(guide = "none")

p2
