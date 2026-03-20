palindromes<-fread("/local/storage/kav67/Bird_data/palindromes.tsv")


df_long <- palindromes %>%
  pivot_longer(
    cols = c(WholeIdentity, MiddleIdentity50bp, MiddleIdentity20bp,RandomIdentity50bp,MiddleIdentity10bp,RandomIdentity10bp,RandomIdentity20bp),
    names_to = "Region",
    values_to = "Identity"
  )%>%
  mutate(
    Region = factor(
      Region,
      levels = c("RandomIdentity10bp","MiddleIdentity10bp","RandomIdentity20bp","MiddleIdentity20bp","RandomIdentity50bp","MiddleIdentity50bp", "WholeIdentity")
    )
  )



ggplot(df_long, aes(x = Region, y = Identity, fill = Region)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.2, color = "black", alpha = 0.7) +
  theme_bw() +
  scale_fill_manual(
    values = c(
      "RandomIdentity10bp" = "#a6cee3",
      "MiddleIdentity10bp" = "#1f78b4",
      "RandomIdentity20bp" = "#b2df8a",
      "MiddleIdentity20bp" = "#33a02c",
      "RandomIdentity50bp" = "#fb9a99",
      "MiddleIdentity50bp" = "#e31a1c",
      "WholeIdentity"      = "#999999"
    )
  ) +
  labs(title = "Distribution of Identity Across Regions")

df_long <- palindromes %>%
  pivot_longer(
    cols = c(WholeIdentity, MiddleIdentity15bp, RandomIdentity15bp),
    names_to = "Region",
    values_to = "Identity"
  )%>%
  mutate(
    Region = factor(
      Region,
      levels = c("RandomIdentity15bp","MiddleIdentity15bp", "WholeIdentity")
    )
  )
df_long <- df_long %>%
  mutate(Region = recode(Region, 
                         'RandomIdentity15bp' = 'Random 15bp', 
                         'MiddleIdentity15bp' = 'Hairpin 15bp', 
                         'WholeIdentity' = 'Whole Inversion'))


ggplot(df_long, aes(x = Region, y = Identity, fill = Region)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.2, color = "black", alpha = 0.7) +
  theme_classic(base_size = 14) +
  scale_fill_manual(
    values = c(
      "Random 15bp" = "#a6cee3",
      "Hairpin 15bp" = "#1f78b4",
      "Whole Inversion"      = "#999999"
    )
  ) +
  labs(title = "Distribution of Identity in across regions in Inversions")
