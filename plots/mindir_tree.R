library(data.table)
library(dplyr)
library(ggplot2)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(viridis)

# ── Load data ──────────────────────────────────────────────────────────────────

all_species_data <- fread("/local/storage/kav67/clean_birds/all_species_stats_pruned_12052025.csv")
all_species_data$bird <- FALSE
all_species_data[all_species_data$VertClass == "birds", ]$bird <- TRUE
all_species_data <- all_species_data[all_species_data$IGH_AnnotationLevel < 2, ]
all_species_data <- all_species_data[all_species_data$VertClass %in% c("reptiles", "mammals", "birds"), ]
bird_only <- all_species_data[all_species_data$VertClass == "birds", ]

bird_tree_pruned <- read.tree("/local/storage/kav67/clean_birds/vgp_birds.nwk")
tree_species <- gsub('.*"', '', bird_tree_pruned$tip.label)
tree_species <- gsub('"', '', tree_species)
bird_tree_pruned$species <- tolower(tree_species)

# Join Order and proper LatinName from VGP table
vgp_table <- fread("/local/storage/kav67/clean_birds/IGH_VGP_table.tsv")
order_map  <- vgp_table %>%
  select(LatinName, Order) %>%
  distinct() %>%
  mutate(latin_tree = tolower(gsub(" ", "_", LatinName)))

# ── Match to tree ──────────────────────────────────────────────────────────────

bird_only$latin_tree <- tolower(gsub(" ", "_", bird_only$LatinName))

# Join on latin_tree (lowercase underscore) — also brings in proper LatinName
bird_only <- merge(as.data.frame(bird_only %>% select(!LatinName)),
                   order_map %>% select(latin_tree, Order, LatinName),
                   by = "latin_tree", all.x = TRUE)

idx      <- match(tolower(tree_species), bird_only$latin_tree)
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
  scale_fill_viridis_c(name = "IGH\nMinDir", option = "C", na.value = "grey90")

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
  scale_fill_viridis_c(name = "IGL\nMinDir", option = "D", na.value = "grey90")

# ── Save ──────────────────────────────────────────────────────────────────────
p
