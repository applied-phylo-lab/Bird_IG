library(data.table)
library(ggplot2)

df <- fread("/local/storage/kav67/within_species/Songbirds/house_finches/bHaeMex1_pri/kmers_long/kmer_abundance_all.tsv")
df <- fread("/local/storage/kav67/within_species/Songbirds/house_finches/bHaeMex1_pri/kmers_400-450/kmer_abundance_all.tsv")

# For each k, get top kmer abundance
top_df <- df[, .SD[which.max(count)], by = k]

ggplot(top_df, aes(x = k, y = count)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_minimal(base_size = 14) +
  labs(
    x = "k-mer length",
    y = "Top k-mer abundance"
  )

