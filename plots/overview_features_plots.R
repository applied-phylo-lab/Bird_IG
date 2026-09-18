# Data-overview plots that used to live in the bottom half of
# data_prep/overview_features.R: V gene counts per locus mapped onto the pruned
# VGP tree, and the same tree with taxonomic orders highlighted.
#
# Self-contained: reads the outputs of data_prep/build_vgp_tables.R rather than
# relying on objects left in the global environment by another script. Run
# build_vgp_tables.R first, then source this interactively in RStudio.

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ape)
library(tidytree)
library(ggplot2)
library(ggtree)
library(ggstance)

INPUT_DIR <- "/local/storage/kav67/clean_birds/"

summary_features <- fread(file.path(INPUT_DIR, "summary_features.csv"))
vgp_table        <- fread(file.path(INPUT_DIR, "IGH_VGP_table.tsv"))
tree             <- read.tree(file.path(INPUT_DIR, "vgp_birds.nwk"))

# Species -> LatinName, from the table the tree was pruned against.
sp2latin <- unique(vgp_table[, .(Species, LatinName)])

curated <- summary_features %>%
  inner_join(sp2latin, by = "Species") %>%
  filter(!str_detect(Haplotype, "_alt$"))

# ---- V gene counts per locus, as barplots beside the tree --------------------
locus_counts_wide <- curated %>%
  count(LatinName, Locus, name = "NumLoci") %>%
  pivot_wider(names_from = Locus, values_from = NumLoci, values_fill = 0) %>%
  dplyr::rename(label = LatinName) %>%
  filter(label %in% tree$tip.label)

p <- ggtree(tree, layout = "rectangular")

p_loci <- p %>%
  facet_plot(panel = "IGH", data = locus_counts_wide,
             geom = geom_barh, mapping = aes(x = IGH), stat = "identity") %>%
  facet_plot(panel = "IGL", data = locus_counts_wide,
             geom = geom_barh, mapping = aes(x = IGL), stat = "identity")
p_loci

# ---- the same tree with taxonomic orders highlighted ------------------------
# order_nodes is the MRCA of each taxonomic order's tips. Orders that are not
# monophyletic in this tree give a meaningless MRCA spanning unrelated clades;
# the original version dropped those by hardcoded row index
# (order_nodes[-9,], then [c(3,7,8,9,11,12),]), which silently broke whenever the
# tree changed. Filtering on monophyly instead makes the intent explicit.
tip_orders <- as_tibble(tree) %>%
  left_join(unique(curated[, .(LatinName, Order)]),
            by = c("label" = "LatinName")) %>%
  filter(!is.na(Order), !is.na(label))

order_nodes <- tip_orders %>%
  group_by(Order) %>%
  filter(n() > 1) %>%
  summarize(node = MRCA(tree, label),
            n_tips = n(),
            monophyletic = is.monophyletic(tree, label),
            .groups = "drop") %>%
  filter(monophyletic) %>%
  dplyr::rename(type = Order)

message(sprintf("Highlighting %d monophyletic orders; %d non-monophyletic dropped",
                nrow(order_nodes),
                n_distinct(tip_orders$Order) - nrow(order_nodes)))

p + geom_tiplab(size = 2, angle = 90, hjust = 1) +
  layout_dendrogram() +
  scale_y_reverse() +
  geom_hilight(data = order_nodes, aes(node = node, fill = type),
               type = "roundrect", inherit.aes = FALSE)
