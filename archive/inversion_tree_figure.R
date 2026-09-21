library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(viridis)
order_nodes<-order_nodes[-c(6,13,15),]
p <- ggtree(tree_pruned) +
  theme(
    strip.text = element_blank(),   # remove panel titles
    strip.background = element_blank()
  )

p_tree <- p +
  geom_hilight(
    data = order_nodes,
    aes(node = node),
    fill = "grey85",
    color = NA,
    alpha = 0.6,
    inherit.aes = FALSE
  )

p_tree <- p_tree +
  geom_cladelab(
    data = order_nodes,
    mapping=aes(node = node, label = type),
    align = TRUE,
    offset = 0,
    fontsize = 3
  )
p_tree
p_h1 <- p_tree +
  geom_fruit(
    data = species_stats_pruned,
    geom = geom_tile,
    mapping = aes(y = LatinName, x = 1, fill = avg_inv_len),
    width = 0.08,
    offset = 0.2,
  ) +
  scale_fill_viridis_c(option = "magma", name = NULL,direction=-1)


p_h2 <- p_h1 +
  new_scale_fill() +
  geom_fruit(
    data = species_stats_pruned,
    geom = geom_tile,
    mapping = aes(y = LatinName, x = 1, fill = avg_num_inversions),
    width = 0.08
  ) +
  scale_fill_viridis_c(option = "plasma", name = NULL,direction=-1)


p_h3 <- p_h2 +
  new_scale_fill() +
  geom_fruit(
    data = species_stats_pruned,
    geom = geom_tile,
    mapping = aes(y = LatinName, x = 1, fill = frac_genes_on_inv),
    width = 0.08
  ) +
  scale_fill_viridis_c(option = "cividis", name = NULL)
p_h3




p <- ggtree(tree_pruned) +
  theme_tree() +   # removes tree scale axis
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )+ylim(0, length(tree_pruned$tip.label) + 6)

p_tree <- p +
  geom_cladelab(
    data = order_nodes,
    mapping = aes(node = node, label = type),
    align = TRUE,
    offset = 0,
    fontsize = 4
  )
#geom_hilight(
#  data = order_nodes,
#  aes(node = node),
#  fill = "lightgrey",
#  color = NA,
#  alpha = 0.6,
#  inherit.aes = FALSE
#) +


p_b1 <- p_tree +
  geom_fruit(
    data = species_stats_pruned,
    geom = geom_col,
    mapping = aes(y = LatinName, x = avg_inv_len),
    orientation = "y",
    width = 0.5,
    offset = 0.2,
    fill = "#003049",
    axis.params = list(
      axis = "x",
      text.size = 4,
      title = "Avg. Inversion \nLength",
      title.size = 4,
      title.height = 0.045,
      nbreak =2,
      vjust = 1.5,
      line.size =0
    )
  ) 
p_b2 <- p_b1 +
  geom_fruit(
    data = species_stats_pruned,
    geom = geom_col,
    mapping = aes(y = LatinName, x = avg_num_inversions),
    orientation = "y",
    width = 0.5,
    offset = 0,
    fill = "#669bbc",
    axis.params = list(
      axis = "x",
      title = "Inversion Count",
      title.size = 4,
      text.size = 4,
      title.height = 0.045,
      nbreak =2,
      vjust = 1.5,
      line.size =0
    )
  ) 
p_b3 <- p_b2 +
  geom_fruit(
    data = species_stats_pruned,
    geom = geom_col,
    mapping = aes(y = LatinName, x = frac_genes_on_inv*100),
    orientation = "y",
    width = 0.5,
    offset = 0,
    fill = "#90e0ef",
    axis.params = list(
      axis = "x",
      title = "% Gene \non Inversions",
      title.size = 4,
      text.size = 4,
      title.height = 0.045,
      nbreak =2.5,
      vjust = 1.5,
      line.size =0
    )
  ) 
p_b3 +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.grid = element_blank()
  )

