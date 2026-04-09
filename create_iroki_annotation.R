library(phylolm)
library(ape)
library(dplyr)
vgp_bird_table<-fread("/local/storage/kav67/clean_birds/IGH_VGP_table.tsv")

vgp_bird_table<-merge(vgp_bird_table,species_stats_pruned, by = "LatinName",all=FALSE)


ggplot(vgp_bird_table, aes(x=avg_num_inversions,y=NumV))+
  geom_point()+
  theme_classic()
tree_pruned


ggplot(vgp_bird_table, aes(x=NumV,y=frac_genes_on_inv))+
  geom_point()+
  theme_classic()

ggplot(vgp_bird_table, aes(x=avg_inv_len,y=frac_genes_on_inv))+
  geom_point()+
  theme_classic()
# Keep only species present in tree
data_use <- vgp_bird_table %>%
  filter(LatinName %in% tree_pruned$tip.label)

# Drop species from tree not in data
tree_use <- drop.tip(
  tree_pruned,
  setdiff(tree_pruned$tip.label, data_use$LatinName)
)

data_use <- data_use[match(tree_use$tip.label, data_use$LatinName), ]
all(tree_use$tip.label == data_use$LatinName)
rownames(data_use)<-data_use$LatinName
model <- phylolm(
  avg_num_inversions ~ NumV,
  data = data_use,
  phy = tree_use,
  model = "lambda"
)

summary(model)

data_use$log_NumV <- log10(data_use$NumV + 1)
data_use$log_inv  <- log10(data_use$avg_num_inversions + 1)

model <- phylolm(
  log_inv ~ log_NumV,
  data = data_use,
  phy = tree_use,
  model = "lambda"
)
summary(model)
data_use$fitted_phylo <- fitted(model)

ggplot(data_use, aes(x = log_NumV, y = log_inv)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_line(aes(y = fitted_phylo), color = "#003049", size = 1) +
  theme_classic() +
  labs(
    x = "Number of V genes (log)",
    y = "Average number of inversions (log)",
    title = "Phylogenetic regression (lambda model)"
  )

species_stats_pruned
df<-fread("/local/storage/kav67/clean_birds/summary_features_old.csv")
df <- df %>%
  group_by(Species, Locus) %>%
  filter(if (n() > 1) str_detect(Haplotype, "pri") else TRUE) %>%
  ungroup()

#vgp_bird_table<-vgp_bird_table%>%group_by(LatinName) %>%
#  summarize(
#    NumV = mean(NumV, na.rm = TRUE)
#  )

vgp_bird_table <- vgp_bird_table %>%
  group_by(Species, Locus) %>%
  filter(if (n() > 1) str_detect(Haplotype, "pri") else TRUE) %>%
  ungroup()

vgp_bird_table <- vgp_bird_table %>%
  left_join(df %>% select(Species, Order, Haplotype, Locus, Contig, NumV, Length),
            by = c("Species", "Order", "Haplotype", "Locus", "Contig", "NumV"))

vgp_bird_table[vgp_bird_table$Contig=="JAJHSZ010000035.1",]$Length<-545918
vgp_bird_table[vgp_bird_table$Contig=="CM106230.1",]$Length<-196575
#vgp_bird_table<-vgp_bird_table[!is.na(vgp_bird_table$Length),]

annotation_iroki<-vgp_bird_table%>%select(LatinName,Length,NumV,avg_inv_len,avg_num_inversions,inv_cov_frac,frac_genes_on_inv)


colnames(annotation_iroki)[1]<-"name"
colnames(annotation_iroki)[2]<-"bar1_height" #Length
colnames(annotation_iroki)[3]<-"bar2_height" #NumV
colnames(annotation_iroki)[4]<-"bar3_height" #avg_inv_len
colnames(annotation_iroki)[5]<-"bar4_height" #avg_num_inversions
colnames(annotation_iroki)[7]<-"bar5_height" #frac_genes_on_inv
annotation_iroki$bar1_color<-"black"
annotation_iroki$bar2_color<-"darkgrey"
annotation_iroki$bar3_color<-"#003049"
annotation_iroki$bar4_color<-"#669bbc"
annotation_iroki$bar5_color<-"#90e0ef"

annotation_iroki<-annotation_iroki[,-6]
annotation_iroki$bar1_height <- annotation_iroki$bar1_height * -1
annotation_iroki$bar2_height <- annotation_iroki$bar2_height * -1
annotation_iroki$bar3_height <- annotation_iroki$bar3_height * -1
annotation_iroki$bar4_height <- annotation_iroki$bar4_height * -1
annotation_iroki$bar5_height <- annotation_iroki$bar5_height * -1
annotation_iroki$name<-gsub(" ","_",annotation_iroki$name)



write_tsv(annotation_iroki,paste0(dir,"iroki_annotation.tsv"))
