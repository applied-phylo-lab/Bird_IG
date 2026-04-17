INPUT_DIR  <- "/local/storage/kav67/clean_birds/"        # top-level input directory
SUMMARY    <- "/local/storage/kav67/clean_birds/IGH_filtered_table.tsv"
OUTPUT_DIR <- "/local/storage/kav67/clean_birds/trees/"

# ── Theme ─────────────────────────────────────────────────────────────────────
theme_igh <- function() {
  theme_classic(base_size = 12) +
    theme(
      text             = element_text(family = "sans"),
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(size = 10, color = "grey40"),
      panel.grid.minor = element_blank(),
      strip.text       = element_text(face = "bold"),
      legend.position  = "right"
    )
}

meta <- fread(paste0(INPUT_DIR,"trees/all_pairwise_distances.tsv"))
cat(sprintf("Loaded %d rows from summary table\n", nrow(meta)))

ggplot(meta, aes(x = "All", y = genetic_dist)) +
  geom_violin(fill = "steelblue", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  labs(x = "", y = "Genetic distance") +
  theme_classic()

ggplot(meta, aes(x = Order, y = genetic_dist)) +
  geom_violin(fill = "steelblue", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(meta, aes(x = "All", y = genetic_dist)) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(x = "", y = "Genetic distance") +
  theme_classic()


ggplot(meta, aes(x = Order, y = genetic_dist)) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


order_dist <- meta %>%
  group_by(Order) %>%
  summarise(mean_dist = mean(genetic_dist, na.rm = TRUE))


ggplot(order_dist, aes(x = Order, y = 1, fill = mean_dist)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(x = "", y = "", fill = "Mean genetic dist") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank())




ggplot(meta, aes(x = physical_dist, y = genetic_dist)) +
  geom_point(alpha = 0.2, size = 0.5) +
  theme_classic() +
  labs(x = "Physical distance", y = "Genetic distance")

ggplot(meta, aes(x = physical_dist, y = genetic_dist)) +
  geom_point(alpha = 0.15, size = 0.4) +
  geom_smooth(method = "loess", color = "red") +
  theme_classic()



ggplot(meta, aes(x = physical_dist, y = genetic_dist)) +
  geom_point(alpha = 0.15, size = 0.4) +
  geom_smooth(method = "loess", color = "red") +
  facet_wrap(~Order) +
  theme_classic()


ggplot(meta, aes(x = physical_dist, y = genetic_dist)) +
  geom_point(alpha = 0.15, size = 0.4) +
  geom_smooth(method = "loess", color = "red") +
  scale_x_log10() +
  theme_classic()

