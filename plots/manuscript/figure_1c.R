# Figure 1C -- distributions of locus length, V gene count, strand bias,
# inversion length and inversion coverage, beside the phylolm scatter.
#
# Run data_prep/build_iroki_annotation.R first; it writes the input table.
# Figures are exported by hand for the manuscript -- the ggsave call at the bottom
# is for unattended reproduction, not final artwork.

library(data.table)
library(dplyr)
library(ggplot2)
library(ape)
library(phylolm)
library(patchwork)

source(file.path("/home/kav67/Bird_IG", "plots", "_dataset.R"))

# Six-bar table from data_prep/build_iroki_annotation.R. The five-bar
# iroki_annotation.tsv has no strand column, so the strand panel needs this one.
annotation_iroki <- fread(file.path(INPUT_DIR,
                                    sprintf("annotation_iroki_strand%s.tsv", .suffix)))

colnames(annotation_iroki) <- c("name", "locus_length", "V_number", "strand",
                                "inversion_length", "inversion_number",
                                "fraction_inversions",
                                "bar1_color", "bar2_color", "bar3_color",
                                "bar4_color", "bar5_color", "bar6_color")

# iroki draws bars leftwards so the builder negates every height; flip them back.
annotation_iroki$locus_length        <- annotation_iroki$locus_length * -1
annotation_iroki$V_number            <- annotation_iroki$V_number * -1
annotation_iroki$strand              <- annotation_iroki$strand * -1
annotation_iroki$inversion_length    <- annotation_iroki$inversion_length * -1
annotation_iroki$inversion_number    <- annotation_iroki$inversion_number * -1
annotation_iroki$fraction_inversions <- annotation_iroki$fraction_inversions * -1 * 100


# Create plots

p_locus <- ggplot(annotation_iroki, aes(x = locus_length / 1e6)) +
  geom_histogram(bins = 60, fill = "grey") +
  #geom_vline(
  #  xintercept = mean(annotation_iroki$locus_length, na.rm = TRUE) / 1e6,
  #  linetype = "dashed",
  #  size = 0.5
  #) +
  scale_x_continuous(
    breaks = seq(0.5, 2.5, by = 0.5),
  ) +
  labs(
    x = "Locus length (Mbp)",
    y = "Count"
  ) +
  theme_classic()


p_strand <- ggplot(annotation_iroki, aes(x = strand)) +
  geom_histogram(bins = 50, fill = "grey") +
  #geom_vline(xintercept = mean(annotation_iroki$strand, na.rm = TRUE),
  #           linetype = "dashed", size = 0.5) +
  scale_x_continuous(breaks = c(25, 50, 75, 100)) +
  labs(
       x = "Genes on positive strand (%)",
       y = "") +
  theme_classic()


# Function to make histogram with mean line
p_inv_len <-ggplot(annotation_iroki, aes(x = inversion_length)) +
    geom_histogram(bins = 40, fill = "#87b4dc") +
    #geom_vline(xintercept = mean_val, linetype = "dashed", size = 0.5) +
    labs(
      x = "Avg. Inversion Length (bp)",
      y = "Count") +
    theme_classic()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10))


p_frac_inv <- ggplot(annotation_iroki, aes(x = fraction_inversions)) +
  geom_histogram(bins = 40, fill = "#87b4dc") +
  #geom_vline(xintercept = mean(annotation_iroki$fraction_inversions, na.rm = TRUE),
  #           linetype = "dashed", size = 0.5) +
  scale_x_continuous(
    breaks = c(25, 50, 75, 100)
  ) +coord_cartesian(xlim = c(75, 101))+
  labs(
       x = "Genes in Inversion region (%)",
       y = "Count"
       ) +
  theme_classic()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10))




p_scatter <- ggplot(annotation_iroki, 
                    aes(x = V_number, y = inversion_number)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "V_number vs Inversion_number",
       x = "V_number",
       y = "Inversion_number") +
  theme_classic()


# ---- phylolm panel ----------------------------------------------------------
# p_phylo used to arrive as a global left behind by plots/phylolm_tree.R, so this
# script only worked if that one had been sourced first, in the same session.
# Built here instead so the figure stands alone.
vgp  <- fread(VGP_TABLE)
inv  <- fread(file.path(INPUT_DIR, "inversion_stats.tsv"))
tree <- read.tree(TREE_FILE)
setnames(inv, c("haplotype", "contig"), c("Haplotype", "Contig"))

trait <- vgp %>%
  group_by(LatinName, Haplotype) %>% slice_max(NumV, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  inner_join(inv, by = c("Haplotype", "Contig")) %>%
  group_by(LatinName) %>%
  summarise(NumV = log(mean(NumV)),
            num_inversions = log(mean(num_inversions)), .groups = "drop")

tsp   <- gsub('"', "", tree$tip.label)
trait <- trait[match(tsp, gsub(" ", "_", trait$LatinName)), ]
keep  <- !(is.na(trait$NumV) | is.na(trait$num_inversions))
trait$tip <- tree$tip.label
trait <- as.data.frame(trait)[keep, ]
rownames(trait) <- trait$tip
phy   <- drop.tip(tree, tree$tip.label[!keep])
trait <- trait[phy$tip.label, ]

model <- phylolm(NumV ~ num_inversions, data = trait, phy = phy, model = "lambda")
co    <- summary(model)$coefficients
stars <- if (co[2, 4] < 0.001) "***" else if (co[2, 4] < 0.01) "**" else
         if (co[2, 4] < 0.05) "*" else "ns"
label_text <- sprintf("\u03b2  = %.3f%s\n\u03bb  = %.5f\nR\u00b2 = %.3f",
                      coef(model)[2], stars, model$optpar,
                      summary(model)$r.squared)

p_phylo <- ggplot(trait, aes(x = num_inversions, y = NumV)) +
  geom_point(size = 2, alpha = 0.9, color = "#87b4dc") +
  geom_abline(intercept = coef(model)[1], slope = coef(model)[2],
              linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = -Inf, y = Inf, label = label_text,
           hjust = -0.2, vjust = 1.4, size = 4) +
  labs(x = "# Inversions (log)", y = "# V genes (log)") +
  theme_classic() +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 10))

message(sprintf("[figure_1c] phylolm on %d species (%s)", nrow(trait), DATASET))

# ---- assembled figure -------------------------------------------------------
# Three panels side by side: average inversion length, genes in inversion region,
# and the phylolm scatter. This is the last of the layout lines the original
# script tried (`p_inv_len | p_frac_inv | p_phylo`). p_locus (locus length) and
# p_strand (genes on the positive strand) are still built above and print to the
# viewer, but are not part of this figure.
figure_1c <- p_inv_len | p_frac_inv | p_phylo
figure_1c

save_fig("Figure1C.svg", figure_1c)

cat(sprintf("mean fraction of genes in inversions: %.3f\n",
            mean(annotation_iroki$fraction_inversions, na.rm = TRUE)))
cat(sprintf("mean inversion length: %.1f bp\n",
            mean(annotation_iroki$inversion_length, na.rm = TRUE)))
