library(data.table)
library(ggplot2)

df <- fread("/local/storage/kav67/within_species/Songbirds/house_finches/bHaeMex1_pri/kmers/kmer_abundance_all.tsv")

# For each k, get top kmer abundance
top_df <- df[, .SD[which.max(count)], by = k]

ggplot(top_df, aes(x = k, y = count)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  theme_classic(base_size = 14) +
  labs(
    x = "k-mer length",
    y = "Top k-mer abundance"
  )
