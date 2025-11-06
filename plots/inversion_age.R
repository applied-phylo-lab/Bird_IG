library(viridis)
sub_tree

shared_inversions<-fread("/local/storage/kav67/birds/Songbirds/Songbirds_shared_inversions.tsv")
shared_inversions<-fread("/local/storage/kav67/birds/Songbirds/House_Finch_shared_inversions.tsv")
summary_filled

dist_matrix <- cophenetic.phylo(sub_tree)

name_map <- summary_filled %>%
  select(Species, LatinName) %>%
  distinct()  # just in case multiple rows per species

# Merge mapping for Species1 and Species2
shared_inversions <- shared_inversions %>%
  left_join(name_map, by = c("Species1" = "Species")) %>%
  rename(Latin1 = LatinName) %>%
  left_join(name_map, by = c("Species2" = "Species")) %>%
  rename(Latin2 = LatinName)

# Add distance column to shared_inversions
shared_inversions <- shared_inversions %>%
  rowwise() %>%
  mutate(tree_distance = dist_matrix[Latin1, Latin2]) %>%
  ungroup()


ggplot(shared_inversions, aes(x = tree_distance, y = num_inversions, color = Inversion_Length_Category)) +
  geom_point(alpha = 0.7, size = 3) +
  theme_minimal() +
  labs(
    x = "Genetic distance (tree-based)",
    y = "Number of shared inversions",
    color = "Inversion length category",
    title = "Shared Inversions vs Genetic Distance"
  ) +
  theme(text = element_text(size = 14))+
  scale_color_viridis(discrete = TRUE, option = "D", direction=-1) 


ggplot(shared_inversions, aes(x = tree_distance, y = num_inversions, color = Identity_Category)) +
  geom_point(alpha = 0.7, size = 2) +
  facet_wrap(~ Inversion_Length_Category, scales = "free_y") +
  scale_color_viridis(discrete = TRUE, option = "D", direction=-1) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Genetic distance (tree-based)",
    y = "Number of shared inversions",
    color = "Identity (%) category",
    title = "Shared Inversions vs Genetic Distance"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )



