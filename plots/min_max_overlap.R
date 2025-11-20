max_overlap<-fread("/local/storage/kav67/Bird_data/maximal_overlap.tsv")



overlap_all <- data.table(max_overlap %>%
                          group_by(Order, Species, Haplotype)) %>%
  rename(label = Order) 

p <- ggtree(order_tree, layout = "rectangular") +
  geom_tiplab(size = 3, align = TRUE, linetype = NA, linesize = 0.5)



# Add boxplots per order
p2<-facet_plot(p+xlim_tree(9), panel = "Maximum Overlap",
               data = overlap_all,
               geom_boxplot,
               mapping = aes(x = overlap_len, fill = label))
p2+theme_tree2()


min_overlap<-fread("/local/storage/kav67/Bird_data/minimal_overlap.tsv")



min_overlap_all <- data.table(min_overlap %>%
                            group_by(Order, Species, Haplotype)) %>%
  rename(label = Order) 


# Add boxplots per order
p2<-facet_plot(p+xlim_tree(9), panel = "Minimum Overlap",
               data = min_overlap_all,
               geom_boxplot,
               mapping = aes(x = overlap_len, fill = label))
p2+theme_tree2()
