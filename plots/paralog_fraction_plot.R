library(ggplot2)
library(dplyr)
library(tidyr)
library(phytools)
paralogs<-fread("/local/storage/kav67/birds/IGH_paralogs_all.tsv")

summary_table_IGH<-fread("/local/storage/kav67/birds/IGH_table.tsv")

for( o in unique(paralogs$Order)){
  p<-ggplot(paralogs[paralogs$Order==o], aes(x=group_size))+
    geom_histogram()+
    theme_minimal()+
    labs(title =o)
  print(p)
}

ggplot(paralogs, aes(x = group_size, fill = Order)) +
  geom_bar(position = "stack") +
  labs(title ="Paralog Group Sizes",
       x = "Group Size",
       y = "Count") +
  theme_minimal()

ggplot(paralogs[paralogs$group_size>1,], aes(x = group_size, fill = Order)) +
  geom_bar(position = "stack") +
  labs(title ="Paralog Group Sizes",
       x = "Group Size",
       y = "Count") +
  theme_minimal()

fractions <- data.table(paralogs %>%
  group_by(Order, Species, Haplotype, sample) %>%
  summarise(
    total_genes = n(),
    grouped_genes = sum(group_size > 1),
    fraction_grouped = grouped_genes / total_genes,
    .groups = "drop"
  ))%>%
  rename(label = Order)  


p <- ggtree(order_tree, layout = "rectangular") +
  geom_tiplab(size = 3, align = TRUE, linetype = NA, linesize = 0.5)


# Add boxplots per order
p2<-facet_plot(p+xlim_tree(9), panel = "Fraction of genes in paralog groups (>1)",
               data = fractions,
               geom_boxplot,
               mapping = aes(x = fraction_grouped, fill = label))
p2+theme_tree2()

