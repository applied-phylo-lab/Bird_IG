all_species_data<-fread("/local/storage/kav67/clean_birds/all_species_stats_pruned_12052025.csv")
all_species_data$bird<-FALSE
all_species_data[all_species_data$VertClass=="birds",]$bird<-TRUE
all_species_data<-all_species_data[all_species_data$IGH_AnnotationLevel<2,]
all_species_data<-all_species_data[all_species_data$VertClass %in% c("reptiles","mammals","birds"),]

ggplot(all_species_data, aes(x = IGH_TotalLength / 1e6)) +
  geom_histogram(bins = 100) +
  labs(
    x = "Locus length (Mbp)",
    y = "Count"
  ) +
  theme_classic()

locus_l_comp_hist<-ggplot(all_species_data, aes(x = IGH_TotalLength / 1e6, fill = bird)) +
  geom_histogram(bins = 40, alpha = 0.8, position = "stack") +
  scale_fill_manual(
    values = c("TRUE" = "#87b4dc", "FALSE" = "grey"),
    labels = c("TRUE" = "Birds", "FALSE" = "Mammals and Reptiles"),
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
    position = "stack",
    alpha = 0.8,
    bins = 40
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#87b4dc", "FALSE" = "grey"),
    labels = c("TRUE" = "Birds", "FALSE" = "Mammals and Reptiles"),
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

IGL_strand<-ggplot(all_species_data, aes(x = IGL_MinDir, fill = bird)) +
  geom_histogram(
    position = "identity",
    alpha = 0.8,
    bins = 40
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#87b4dc", "FALSE" = "grey"),
    labels = c("TRUE" = "Birds", "FALSE" = "Mammals and Reptiles"),
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



IGL_length<-ggplot(all_species_data, aes(x = IGL_TotalLength / 1e6, fill = bird)) +
  geom_histogram(bins = 40, alpha = 0.8, position = "identity") +
  scale_fill_manual(
    values = c("TRUE" = "#87b4dc", "FALSE" = "grey"),
    labels = c("TRUE" = "Birds", "FALSE" = "Mammals and Reptiles"),
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


IGL_length+IGL_strand

bird_only<-all_species_data[all_species_data$VertClass=="birds",]

plot_data <- bird_only %>%
  select(IGH_MinDir, IGL_MinDir) %>%
  pivot_longer(cols = everything(),
               names_to = "Locus",
               values_to = "MinDir")

# Plot
ggplot(plot_data, aes(x = MinDir, fill = Locus)) +
  geom_histogram(alpha = 0.7, position = "identity", bins = 30,color = "black",linewidth = 0.25) +
  scale_fill_manual(values = c("IGH_MinDir" = "#87b4dc",
                               "IGL_MinDir" = "#638E6E"),
                    labels = c("IGH", "IGL")) +
  labs(x = "Fraction of genes located on the same strand", y = "Count", fill = "Locus") +
  theme_classic()+
  theme(axis.title = element_text(size = 14),
                       axis.text = element_text(size = 10),
                       legend.position = "none")


plot_data <- bird_only %>%
  mutate(IGH_ContigLength = (IGH_TotalLength / IGH_LocusFraction) ,
         IGL_ContigLength = (IGL_TotalLength / IGL_LocusFraction) ) %>%
  select(IGH_ContigLength, IGL_ContigLength) %>%
  pivot_longer(cols = everything(),
               names_to = "Locus",
               values_to = "ContigLength") %>%
  filter(!is.infinite(ContigLength), !is.nan(ContigLength), ContigLength > 0)

# Plot
ggplot(plot_data, aes(x = ContigLength, fill = Locus)) +
  geom_histogram(alpha = 0.7, position = "identity", bins = 40) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_fill_manual(values = c("IGH_ContigLength" = "#87b4dc",
                               "IGL_ContigLength" = "#638E6E"),
                    labels = c("IGH", "IGL")) +
  labs(x = "Contig Length (bp)", y = "Count", fill = "Locus") +
  theme_classic()

