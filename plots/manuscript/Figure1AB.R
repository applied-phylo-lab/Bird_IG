# Figure 1A/B -- IGH and IGL locus length and strand bias, birds vs other
# vertebrates.
#
# The input is built outside this repo, by match_names_VGP.R operating on
# /local/storage/kav67/IG_annotation_VGP2026/, and is not produced by the
# pipeline; it is read here as a fixed input.
#
# Figures are exported by hand for the manuscript -- the save_fig calls at the
# bottom are for unattended reproduction, not final artwork.

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

source(file.path("/home/kav67/Bird_IG", "plots", "_dataset.R"))

# Locus colours, matching rss_position_oriented.R and mind_dir_p below.
IGH_BLUE  <- "#87b4dc"
IGL_GREEN <- "#638E6E"
OTHER_GREY <- "grey"

# Legend key for the bird/non-bird fill. "Birds" is drawn as a square split on the
# diagonal -- green for IGL, blue for IGH -- so one key says that both coloured
# panels are birds; everything else keeps a plain grey square.
draw_key_bird <- function(data, params, size) {
  fill <- as.character(data$fill)[1]
  if (is.na(fill) || tolower(fill) %in% c("grey", "gray", "#bebebe")) {
    return(grid::rectGrob(gp = grid::gpar(fill = fill, col = NA)))
  }
  grid::grobTree(
    grid::polygonGrob(c(0, 0, 1), c(0, 1, 1), gp = grid::gpar(fill = IGL_GREEN, col = NA)),
    grid::polygonGrob(c(0, 1, 1), c(0, 0, 1), gp = grid::gpar(fill = IGH_BLUE,  col = NA))
  )
}

all_species_data <- fread(file.path(INPUT_DIR, "all_species_stats_pruned_12052025.csv"))
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
        legend.position = "none")


locus_l_comp|locus_strand_comp
locus_l_comp_hist|locus_strand_comp_hist

mean(all_species_data[all_species_data$bird==TRUE,]$IGH_TotalLength)
mean(all_species_data[all_species_data$bird==FALSE,]$IGH_TotalLength)

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
mind_dir_p<-ggplot(plot_data, aes(x = MinDir, fill = Locus)) +
  geom_histogram(alpha = 0.7, position = "identity", bins = 30) +
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
contig_length_p<-ggplot(plot_data, aes(x = ContigLength, fill = Locus)) +
  geom_histogram(alpha = 0.7, position = "identity", bins = 40) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_fill_manual(values = c("IGH_ContigLength" = "#87b4dc",
                               "IGL_ContigLength" = "#638E6E"),
                    labels = c("IGH", "IGL")) +
  labs(x = "Contig Length (bp)", y = "Count", fill = "Locus") +
  theme_classic()+
  theme(
    legend.position = "none",
    strip.text      = element_text(face = "bold"),
    axis.title      = element_text(size = 14),
    axis.text       = element_text(size = 10)
  )

igh_contig_all <- all_species_data %>%
  mutate(ContigLength = IGH_TotalLength / IGH_LocusFraction) %>%
  filter(!is.infinite(ContigLength), !is.nan(ContigLength), ContigLength > 0)

igh_contig_length_p <- ggplot(igh_contig_all, aes(x = ContigLength, fill = bird)) +
  geom_histogram(alpha = 0.7, position = "stack", bins = 40) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_fill_manual(
    values = c("TRUE" = "#87b4dc", "FALSE" = "grey"),
    labels = c("TRUE" = "Birds", "FALSE" = "Mammals and Reptiles"),
    name = NULL
  ) +
  labs(x = "IGH Contig Length (bp, log scale)", y = "Count") +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text  = element_text(size = 10),
        legend.position = "none")

igl_contig_all <- all_species_data %>%
  mutate(ContigLength = IGL_TotalLength / IGL_LocusFraction) %>%
  filter(!is.infinite(ContigLength), !is.nan(ContigLength), ContigLength > 0)

igl_contig_length_p <- ggplot(igl_contig_all, aes(x = ContigLength, fill = bird)) +
  geom_histogram(alpha = 0.7, position = "stack", bins = 40) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_fill_manual(
    values = c("TRUE" = "#87b4dc", "FALSE" = "grey"),
    labels = c("TRUE" = "Birds", "FALSE" = "Mammals and Reptiles"),
    name = NULL
  ) +
  labs(x = "IGL Contig Length (bp, log scale)", y = "Count") +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text  = element_text(size = 10),
        legend.text = element_text(size = 10))

figure_1ab      <- locus_l_comp | locus_strand_comp
figure_1ab_hist <- locus_l_comp_hist | locus_strand_comp_hist
figure_contig   <- igh_contig_length_p | igl_contig_length_p

figure_1ab
figure_1ab_hist
figure_contig

# IGL_length + IGL_strand, and contig_length_p + mind_dir_p, were both assembled
# as bare expressions and never saved -- hence figures/manuscript_Sept/
# IGL_histograms.svg and IGH_IGL_contig_strand.svg had no counterpart in v1/v2.
# Both objects are built above; these are the same composites, at the same
# 1103 x 354 pt canvas as the hand-exported versions.
figure_igl_hist          <- IGL_length + IGL_strand
figure_igh_igl_contig_sd <- contig_length_p + mind_dir_p

# ---- IGL in its own colour ---------------------------------------------------
# Rebuilt rather than patched. Adding a scale to a finished plot only triggers
# "Scale for fill is already present", and key_glyph cannot be set on an existing
# layer -- adding a second geom_histogram to change it silently double-plots.

bird_fill <- function(colour, key = NULL) {
  list(
    scale_fill_manual(values = c("TRUE" = colour, "FALSE" = OTHER_GREY),
                      labels = c("TRUE" = "Birds", "FALSE" = "Mammals and Reptiles"),
                      name = NULL),
    theme_classic(),
    theme(axis.title = element_text(size = 14),
          axis.text  = element_text(size = 10),
          legend.text = element_text(size = 10))
  )
}

IGL_length_green <- ggplot(all_species_data, aes(x = IGL_TotalLength / 1e6, fill = bird)) +
  geom_histogram(bins = 40, alpha = 0.8, position = "identity") +
  bird_fill(IGL_GREEN) + scale_x_log10() +
  labs(x = "Locus length (Mbp, log scale)", y = "Count") +
  theme(legend.position = "none")

IGL_strand_green <- ggplot(all_species_data, aes(x = IGL_MinDir, fill = bird)) +
  geom_histogram(position = "identity", alpha = 0.8, bins = 40) +
  bird_fill(IGL_GREEN) +
  labs(x = "Fraction of genes located on the same strand", y = "Count")

figure_igl_hist_green <- IGL_length_green + IGL_strand_green

# Contig length: IGH in blue, IGL in green. Only the IGL panel carries the legend,
# and its "Birds" key is the split green/blue square, saying that both coloured
# panels are birds while everything else is grey.
igh_contig_blue <- ggplot(igh_contig_all, aes(x = ContigLength, fill = bird)) +
  geom_histogram(bins = 40, alpha = 0.8, position = "identity") +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  bird_fill(IGH_BLUE) +
  labs(x = "IGH Contig Length (bp, log scale)", y = "Count") +
  theme(legend.position = "none")

igl_contig_green <- ggplot(igl_contig_all, aes(x = ContigLength, fill = bird)) +
  geom_histogram(bins = 40, alpha = 0.8, position = "identity",
                 key_glyph = draw_key_bird) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  bird_fill(IGL_GREEN) +
  labs(x = "IGL Contig Length (bp, log scale)", y = "Count")

figure_contig_green <- igh_contig_blue | igl_contig_green

figure_igl_hist
figure_igh_igl_contig_sd

save_fig("Figure1AB.svg",             figure_1ab)
# Narrower cut of the same figure. figures/manuscript_Sept/Figure1AB_resized.svg
# is 735 x 354 pt, i.e. 980 x 472 px at 96 dpi -- same height as the standard
# canvas, two thirds the width.
save_fig("Figure1AB_resized.svg",     figure_1ab, width = 980)
save_fig("Figure1AB_hist.svg",        figure_1ab_hist)
save_fig("Figure1_contig_len.svg",    figure_contig)
save_fig("IGL_histograms.svg",         figure_igl_hist)
save_fig("IGH_IGL_contig_strand.svg",  figure_igh_igl_contig_sd)
save_fig("IGL_histograms_green.svg",   figure_igl_hist_green)
save_fig("Figure1_contig_len_green.svg", figure_contig_green)

cat(sprintf("IGH contig length, birds:      %.0f bp\n",
            mean(igh_contig_all[igh_contig_all$bird == TRUE, ]$ContigLength)))
cat(sprintf("IGH contig length, non-birds:  %.0f bp\n",
            mean(igh_contig_all[igh_contig_all$bird == FALSE, ]$ContigLength)))
cat(sprintf("IGL contig length, birds:      %.0f bp\n",
            mean(igl_contig_all[igl_contig_all$bird == TRUE, ]$ContigLength)))
cat(sprintf("IGL contig length, non-birds:  %.0f bp\n",
            mean(igl_contig_all[igl_contig_all$bird == FALSE, ]$ContigLength)))

