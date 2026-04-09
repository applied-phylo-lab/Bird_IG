annotation_iroki<-fread("/local/storage/kav67/clean_birds/annotation_iroki_strand.tsv")

colnames(annotation_iroki)<-c("name","locus_length","V_number","strand","inversion_length","inversion_number","fraction_inversions","bar1_color"  ,"bar2_color" , "bar3_color",  "bar4_color" , "bar5_color"  ,"bar6_color" )
annotation_iroki$locus_length<-annotation_iroki$locus_length*-1
annotation_iroki$V_number<-annotation_iroki$V_number*-1
annotation_iroki$strand<-annotation_iroki$strand*-1
annotation_iroki$inversion_length<-annotation_iroki$inversion_length*-1
annotation_iroki$inversion_number<-annotation_iroki$inversion_number*-1
annotation_iroki$fraction_inversions<-annotation_iroki$fraction_inversions*-1
annotation_iroki$fraction_inversions<-annotation_iroki$fraction_inversions*100

library(ggplot2)


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
    geom_histogram(bins = 50, fill = "#87b4dc") +
    #geom_vline(xintercept = mean_val, linetype = "dashed", size = 0.5) +
    labs(
      x = "Avg. Inversion Length (bp)",
      y = "Count") +
    theme_classic()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10))


p_frac_inv <- ggplot(annotation_iroki, aes(x = fraction_inversions)) +
  geom_histogram(bins = 50, fill = "#87b4dc") +
  #geom_vline(xintercept = mean(annotation_iroki$fraction_inversions, na.rm = TRUE),
  #           linetype = "dashed", size = 0.5) +
  scale_x_continuous(
    breaks = c(25, 50, 75, 100)
  ) +
  labs(
       x = "Genes on Inversions (%)",
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


((p_locus | p_strand) /
  (p_inv_len | p_frac_inv) )|
  p_phylo

p_locus |p_strand|p_inv_len | p_frac_inv|p_phylo
p_locus |p_strand
p_inv_len | p_frac_inv|p_phylo
