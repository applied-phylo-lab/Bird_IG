library(ggplot2)
library(dplyr)
library(readr)

# === INPUT ===
inversions <- read_tsv("/local/storage/kav67/clean_birds/Songbirds/House_Finch/bHaeMex1_pri/self_align_contigs/inversion_window_summary.tsv", show_col_types = FALSE)
#california_scrub_jay
# optional: order contigs by size or name
inversions <- inversions %>%
  group_by(chr) %>%
  mutate(chr_len = max(window_end)) %>%
  ungroup() %>%
  arrange(chr, window_start)

# === BUILD CONTIG OFFSETS (for Manhattan-style continuous x-axis) ===
offsets <- inversions %>%
  group_by(chr) %>%
  summarise(chr_len = max(window_end)) %>%
  mutate(chr_start = lag(cumsum(chr_len), default = 0)) %>%
  select(chr, chr_start)

inversions <- inversions %>%
  left_join(offsets, by = "chr") %>%
  mutate(global_center = center + chr_start)

# === PLOT SETTINGS ===
highlight_chr<-"JAPYKC010000061.1"
inversions <- inversions %>%
  mutate(color = ifelse(chr == highlight_chr, "IGH", "other"))

# === MANHATTAN PLOT ===
ggplot(inversions, aes(x = global_center, y = num_inversions, color = color)) +
  geom_point(aes(size = color)) +
  scale_size_manual(values = c("IGH" = 3, "other" = 1)) +
  scale_color_manual(values = c("IGH" = "red", "other" = "black")) +
  labs(
    x = "Genomic position",
    y = "Inversion count per window"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 16),
    axis.text.y = element_text(size = 14)
  )