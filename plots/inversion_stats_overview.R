library(ggplot2)
library(dplyr)
library(tidyr)
library(phytools)
library(tibble)
dir<-"/local/storage/kav67/Bird_data/"
inversion_stats<-fread(paste0(dir,"IGH_inversions_all.tsv"))
species_data<-inversion_stats
summary_table_IGH<-fread(paste0(dir,"IGH_table.tsv")) #IGH_filtered_table.tsv

# Summarize by Order (average across samples per order)
inversion_stats_long <- inversion_stats[inversion_stats$minlen==250,] %>%
  mutate(
    inv_cov_frac = inv_cov_len / total_seq_length
  ) %>%
  select(order, avg_inv_len, inv_cov_frac, frac_genes_on_inv) %>%
  pivot_longer(
    cols = c(avg_inv_len, inv_cov_frac, frac_genes_on_inv),
    names_to = "metric",
    values_to = "value"
  )

# Optional: nicer facet labels
metric_labels <- c(
  avg_inv_len = "Average inversion length",
  inv_cov_frac = "Fraction of inverted region",
  frac_genes_on_inv = "Fraction of genes on inversions"
)

# Plot
ggplot(inversion_stats_long, aes(x = order, y = value, fill = order)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.8) +
  facet_wrap(.~metric, scales = "free_y", labeller = as_labeller(metric_labels), ncol = 3) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.background = element_rect(fill = "grey90")
  ) +
  labs(x = "Order", y = "Value", title = "Inversion statistics across orders")

inversion_stats1000<-inversion_stats[inversion_stats$minlen==250,]
inversion_stats1000 <- merge(
  inversion_stats1000,
  summary_table_IGH[, c("Haplotype", "LatinName")],
  by.x = "haplotype",
  by.y = "Haplotype",
  all.x = TRUE
)

species_stats <- inversion_stats1000[!is.na(inversion_stats1000$LatinName),] %>%
  filter(!grepl("MiscBirds", LatinName))  %>%
  group_by(LatinName) %>%
  summarize(
    avg_inv_len = mean(avg_inv_len, na.rm = TRUE),
    inv_cov_frac = mean(inv_cov_len / total_seq_length, na.rm = TRUE),
    frac_genes_on_inv = mean(frac_genes_on_inv, na.rm = TRUE)
  )

#species_stats$LatinName<-gsub(" ","_",species_stats$LatinName)

tips_in_data <- intersect(VGP_tree$tip.label, species_stats$LatinName)

tree_pruned <- drop.tip(VGP_tree, setdiff(VGP_tree$tip.label, tips_in_data))

species_stats_pruned<-species_stats[species_stats$LatinName %in% tips_in_data,]
species_data <- inversion_stats1000 %>%
  filter(!grepl("MiscBirds", LatinName))

p <- ggtree(tree_pruned, layout = "rectangular")
tree_pruned_data <- as_tibble(tree_pruned)
tree_pruned_data<-tree_pruned_data %>%
  left_join(
      species_data%>%       
      select(LatinName, order),
      by = c("label"="LatinName")
  )

order_nodes <- tree_pruned_data %>%
  # Keep only tips with orders
  filter(!is.na(order) & !is.na(label)) %>%
  group_by(order) %>%
  summarize(
    node = MRCA(tree_pruned, label),   # MRCA() from ggtree
    .groups = "drop"
  ) %>%
  rename(type = order)  # rename column to type
order_nodes <- order_nodes[-9,]

p <- ggtree(tree_pruned)

# Add all barplots first
p1 <- facet_plot(
  p, panel = "Average Inversion Length",
  data = species_stats_pruned,
  geom = geom_barh,
  mapping = aes(x = avg_inv_len),
  stat = "identity"
)

p2 <- facet_plot(
  p1, panel = "Fractions of locus covered by Inversions",
  data = species_stats_pruned,
  geom = geom_barh,
  mapping = aes(x = inv_cov_frac),
  stat = "identity"
)

p3 <- facet_plot(
  p2, panel = "Fractions of Genes lying on inversions",
  data = species_stats_pruned,
  geom = geom_barh,
  mapping = aes(x = frac_genes_on_inv),
  stat = "identity"
)

# Now add highlights *after* all facets
p3 + geom_hilight(
  data = order_nodes, 
  aes(node = node, fill = type),
  type = "roundrect",
  inherit.aes = FALSE
)+theme_tree2()



species_data <- inversion_stats1000%>%
  filter(!grepl("MiscBirds", order))
#species_data$LatinName<-gsub(" ","_",species_data$LatinName)
species_data$LatinName<-gsub("_"," ",species_data$LatinName)
# Drop tips from tree
bird_tree_pruned<-sub_tree_pruned
tree_filtered <- drop.tip(bird_tree_pruned,
                          setdiff(bird_tree_pruned$tip.label, species_data$LatinName))
tips_in <- intersect(bird_tree_pruned$tip.label, species_data$LatinName)

species_data<-species_data[species_data$LatinName %in% tips_in,]
species_to_order <- species_data %>%
  distinct(LatinName, order) %>%
  deframe()

tree_df <- as_tibble(tree_filtered)
ntips <- Ntip(tree_filtered)

# Add a logical column indicating whether the row is a tip
tree_df <- tree_df %>%
  mutate(isTip = node <= ntips)
tree_df$label <- as.character(tree_df$label)

# Add order info only for tips
tree_df <- tree_df %>%
  mutate(order = ifelse(isTip, species_to_order[label], NA))

# Collapse tips by order
order_tree <- tree_filtered

# Loop over orders to collapse monophyletic species
orders <- unique(species_to_order)
for (ord in orders) {
  tips_ord <- names(species_to_order[species_to_order == ord])
  if (length(tips_ord) > 1) {
    # collapse tips into single tip
    order_tree <- collapse.singles(order_tree) # placeholder, see below
  }
}
for( ord in orders){
  tips_ord <- names(species_to_order[species_to_order == ord])
  mrca_node <- getMRCA(order_tree, tips_ord)
  
  # Collapse MRCA node to a single tip (replace with order name)
  order_tree <- bind.tip(order_tree, ord, where = mrca_node)
  order_tree <- drop.tip(order_tree, tips_ord)
}


plot_species_data <- species_data %>%
  select(order, LatinName, num_inversions, avg_inv_len, inv_cov_len, total_seq_length, frac_genes_on_inv) %>%
  mutate(inv_cov_frac = inv_cov_len / total_seq_length) %>% 
  rename(label = order)  


order_tree<-ladderize(order_tree, right=FALSE)
p <- ggtree(order_tree, layout = "rectangular") +
  geom_tiplab(size = 3, align = TRUE, linetype = NA, linesize = 0.5)


# Add boxplots per order

p2<-facet_plot(p+xlim_tree(9), panel = "Average Inversion Length",
           data = plot_species_data,
           geom_boxplot,
           mapping = aes(x = avg_inv_len, group = label, fill = label))

p3<-facet_plot(p2, panel = "Inversion Count",
               data = plot_species_data,
               geom_boxplot,
               mapping = aes(x = num_inversions, group = label, fill = label))
p4<-facet_plot(p3, panel = "Fraction of genes on inversions",
               data = plot_species_data,
               geom_boxplot,
               mapping = aes(x = frac_genes_on_inv, group = label, fill = label))

p5<-facet_plot(p4, panel = "Fraction of locus covered by inversions",
               data = plot_species_data,
               geom_boxplot,
               mapping = aes(x = inv_cov_frac, group = label, fill = label))
p5+theme_tree2()
