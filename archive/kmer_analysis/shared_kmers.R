library(data.table)
library(dplyr)
library(stringr)
library(tidyr)

# --- 1. List all kmers.csv files recursively ---
dir_path <- "/local/storage/kav67/within_species/"  # change this
files <- list.files(path = dir_path, pattern = "kmers\\.csv$", recursive = TRUE, full.names = TRUE)

# --- 2. Read each file and clean ---
all_kmers <- lapply(files, function(f) {
  # read file
  df <- fread(f, header = FALSE)
  
  # assign proper column names
  colnames(df) <- c("kmer","count","count_fw","count_rc","meanGap","medianGap","madGap")
  
  # store file path
  df$file <- f
  return(df)
})

# combine all
kmers_df <- rbindlist(all_kmers)

# --- 3. Extract metadata from file path ---
# assume structure: dir/Order/Species/Haplotype/.../kmers.csv
kmers_df <- kmers_df %>%
  mutate(
    path_parts = str_split(file, "/"),
    Order = sapply(path_parts, function(x) x[length(x)-3]),      # adjust based on actual depth
    Species = sapply(path_parts, function(x) x[length(x)-2]),
    Haplotype = sapply(path_parts, function(x) x[length(x)-1])
  ) %>%
  select(-path_parts,-file)

# --- 4. Check for shared kmers ---
# within species
shared_within_species <- kmers_df %>%
  group_by(Species, kmer) %>%
  summarise(n_haplotypes = n_distinct(Haplotype), .groups="drop") %>%
  filter(n_haplotypes > 1)

# across species
shared_across_species <- kmers_df %>%
  group_by(kmer) %>%
  summarise(n_species = n_distinct(Species), .groups="drop") %>%
  filter(n_species > 1)

# --- 5. Optional: inspect results ---
print(head(shared_within_species))
print(head(shared_across_species))


kmers_species <- kmers_df %>%
  group_by(kmer) %>%
  summarise(
    total_count = sum(count),
    n_species = n_distinct(Species),
    .groups = "drop"
  )

# Plot: number of species vs abundance
ggplot(kmers_species, aes(x = n_species, y = total_count)) +
  geom_jitter(width = 0.2, height = 0) +
  scale_y_log10() +
  labs(
    x = "Number of species sharing kmer",
    y = "Total abundance (sum of counts)",
    title = "Shared kmers across species"
  ) +
  theme_minimal()


# pick top 50 kmers by total abundance
top_kmers <- kmers_df %>%
  group_by(kmer) %>%
  summarise(total_count = sum(count), .groups="drop") %>%
  arrange(desc(total_count)) %>%
  slice(1:50) %>%
  pull(kmer)

heatmap_df <- kmers_df %>%
  filter(kmer %in% top_kmers) %>%
  group_by(Species, Haplotype, kmer) %>%
  summarise(count = sum(count), .groups="drop") %>%
  pivot_wider(names_from = kmer, values_from = count, values_fill = 0)

# Convert to long format for ggplot
heatmap_long <- heatmap_df %>%
  pivot_longer(cols = -c(Species, Haplotype), names_to = "kmer", values_to = "count")

ggplot(heatmap_long, aes(x = kmer, y = interaction(Species,Haplotype), fill = count)) +
  geom_tile() +
  scale_fill_viridis_c(option="C",direction=-1) +
  labs(x = "kmer", y = "Species-Haplotype", title = "Top 50 kmers across species and haplotypes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

