rss <- read_tsv("/local/storage/kav67/clean_birds/RSS_stats.tsv", show_col_types = FALSE)

rss <- rss %>%
  mutate(
    Heptamer = if_else(is.na(Heptamer), "NO_HEPTAMER", Heptamer)
  )

#rss<-rss%>% filter(!grepl("Songbirds", Source))

top_heptamers <- rss %>%
  filter(Heptamer != "NO_HEPTAMER",
         `Passes Filtering` ==TRUE) %>%
  count(Heptamer, sort = TRUE) %>%
  slice_head(n = 20) %>%
  pull(Heptamer)
top_heptamers<-c(top_heptamers,"GATAGTG")

df_top <- rss %>%
  filter(Heptamer %in% top_heptamers | Heptamer == "NO_HEPTAMER")

df_top <- df_top %>%
  mutate(
    Heptamer = fct_reorder(Heptamer, inversion_density, .fun = median),
    Heptamer = fct_relevel(Heptamer, "NO_HEPTAMER")
  )

ggplot(df_top, aes(x = Heptamer, y = inversion_density)) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Heptamer (top 20 + no heptamer)",
    y = "Inversion density",
    title = "Inversion density: genes with vs without RSS heptamers"
  )

ggplot(df_top, aes(x = Heptamer, y = inversion_density_pct)) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Heptamer (top 20 + no heptamer)",
    y = "Inversion density percentile",
    title = "Inversion density: genes with vs without RSS heptamers"
  )


df_binary <- rss %>%
  mutate(
    HeptamerGroup = if_else(Heptamer == "NO_HEPTAMER", "No RSS", "Has RSS")
  )

ggplot(df_binary, aes(x = HeptamerGroup, y = inversion_density)) +
  geom_violin(fill = "steelblue", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_classic() +
  labs(
    x = "",
    y = "Inversion density",
    title = "Inversion density: RSS vs no RSS"
  )

ggplot(df_binary, aes(x = HeptamerGroup, y = inversion_density_pct)) +
  geom_violin(fill = "steelblue", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_classic() +
  labs(
    x = "",
    y = "Inversion density percentile",
    title = "Inversion density: RSS vs no RSS"
  )


summary_df <- df_top %>%
  group_by(Heptamer) %>%
  summarise(
    total = n(),
    zero_inv = sum(inversion_density == 0, na.rm = TRUE),
    pct_zero = zero_inv / total,
    .groups = "drop"
  )

ggplot(summary_df,
       aes(x = fct_reorder(Heptamer, pct_zero), y = pct_zero)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_classic() +
  labs(
    x = "Heptamer",
    y = "Fraction of genes with inversion_density = 0",
    title = "Proportion of genes without inversions by heptamer"
  )

ggplot(summary_df,
       aes(x = fct_reorder(Heptamer, pct_zero), y = pct_zero)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = paste0("n=", total)),
            hjust = -0.1, size = 3) +
  coord_flip() +
  theme_classic() +
  labs(
    x = "Heptamer",
    y = "Fraction with inversion_density = 0"
  ) +
  ylim(0, 1.05)
