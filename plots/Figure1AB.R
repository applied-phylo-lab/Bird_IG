all_species_data<-fread("/local/storage/kav67/clean_birds/all_species_stats_pruned_12052025.csv")
all_species_data$bird<-FALSE
all_species_data[all_species_data$VertClass=="birds",]$bird<-TRUE
all_species_data<-all_species_data[all_species_data$IGH_AnnotationLevel<2,]

ggplot(all_species_data, aes(x = IGH_TotalLength / 1e6)) +
  geom_histogram(bins = 100) +
  labs(
    x = "Locus length (Mbp)",
    y = "Count"
  ) +
  theme_classic()

locus_l_comp_hist<-ggplot(all_species_data, aes(x = IGH_TotalLength / 1e6, fill = bird)) +
  geom_histogram(bins = 40, alpha = 0.8, position = "identity") +
  scale_fill_manual(
    values = c("TRUE" = "#87b4dc", "FALSE" = "grey"),
    labels = c("TRUE" = "Bird Species", "FALSE" = "Other Species"),
    name = NULL
  ) +
  scale_x_log10() +
  labs(
    x = "Locus length (Mbp, log scale)",
    y = "Count"
  ) +
  theme_classic()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        legend.position = "none")

locus_l_comp<-ggplot(all_species_data, aes(x = IGH_TotalLength / 1e6, fill = bird)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(
    values = c("TRUE" = "#87b4dc", "FALSE" = "grey"),
    labels = c("TRUE" = "Bird Species", "FALSE" = "Other Species"),
    name = NULL
  ) +
  scale_x_log10() +
  labs(x = "Locus length (Mbp, log scale)", y = "Density") +
  theme_classic()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        legend.position = "none")

locus_strand_comp<-ggplot(all_species_data, aes(x = IGH_MinDir, fill = bird)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(
    values = c("TRUE" = "#87b4dc", "FALSE" = "grey"),
    labels = c("TRUE" = "Bird Species", "FALSE" = "Other Species"),
    name = NULL
  ) +
  labs(x = "Fraction of genes located on the same strand", y = "Density") +
  theme_classic()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10))

locus_strand_comp_hist<-ggplot(all_species_data, aes(x = IGH_MinDir, fill = bird)) +
  geom_histogram(
    position = "identity",
    alpha = 0.8,
    bins = 40
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#87b4dc", "FALSE" = "grey"),
    labels = c("TRUE" = "Bird Species", "FALSE" = "Other Species"),
    name = NULL
  ) +
  labs(
    x = "Fraction of genes located on the same strand",
    y = "Count"
  ) +
  theme_classic()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10))


locus_l_comp|locus_strand_comp
locus_l_comp_hist|locus_strand_comp_hist

median(all_species_data[all_species_data$bird==TRUE,]$IGH_TotalLength)
median(all_species_data[all_species_data$bird==FALSE,]$IGH_TotalLength)

ggplot(all_species_data, aes(x = IGL_MinDir, fill = bird)) +
  geom_histogram(
    position = "identity",
    alpha = 0.8,
    bins = 40
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#87b4dc", "FALSE" = "grey"),
    labels = c("TRUE" = "Bird Species", "FALSE" = "Other Species"),
    name = NULL
  ) +
  labs(
    x = "Fraction of genes located on the same strand",
    y = "Count"
  ) +
  theme_classic()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10))
