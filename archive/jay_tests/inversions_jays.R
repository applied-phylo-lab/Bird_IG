library(data.table)
library(ggplot2)
library(viridis)

# ---- read the three summaries ----
AW <- fread("/local/storage/kav67/jays/inversion_summary/AW_summary.tsv")
AI <- fread("/local/storage/kav67/jays/inversion_summary/AI_summary.tsv")
AC <- fread("/local/storage/kav67/jays/inversion_summary/AC_summary.tsv")

# add source column
AW[, source := "AW"]
AI[, source := "AI"]
AC[, source := "AC"]

# combine
summary_all <- rbindlist(list(AW, AI, AC), use.names = TRUE, fill = TRUE)

# ---- quick summary stats ----
# ---- summary stats by min_len and source ----
summary_stats <- summary_all[, .(
  mean_avg_len = mean(avg_inv_len, na.rm = TRUE),
  sd_avg_len   = sd(avg_inv_len, na.rm = TRUE),
  mean_num_inv = mean(num_inversions, na.rm = TRUE),
  sd_num_inv   = sd(num_inversions, na.rm = TRUE),
  mean_frac_cov = mean(inv_cov_len / total_seq_length, na.rm = TRUE),
  sd_frac_cov   = sd(inv_cov_len / total_seq_length, na.rm = TRUE)
), by = .(source, minlen)]

print(summary_stats)

df_1000 <- summary_all %>% filter(minlen == 1000)
ggplot(df_1000, aes(x=source,y = avg_inv_len, fill = source)) +
  geom_boxplot() +
  labs(title = "Distribution of inversion length",
       x="",
       y = "Mean inversion length",
       fill = "Species") +
  theme_bw()

ggplot(summary_all, aes(x=source,y = avg_inv_len, fill = source)) +
  geom_boxplot() +
  facet_wrap(~minlen) +
  labs(title = "Distribution of inversion length",
       x="",
       y = "Mean inversion length",
       fill = "Species") +
  theme_bw()

ggplot(df_1000, aes(x=source,y = num_inversions, fill = source)) +
  geom_boxplot() +
  labs(title = "Number of inversions",
       x="",
       y = "Number of inversions",
       fill = "Species") +
  theme_bw()

ggplot(summary_all, aes(x=source,y = num_inversions, fill = source)) +
  geom_boxplot() +
  facet_wrap(~minlen) +
  labs(title = "Number of inversions",
       x="",
       y = "Number of inversions",
       fill = "Species") +
  theme_bw()

ggplot(df_1000, aes(x=source,y = inv_cov_len / total_seq_length, fill = source)) +
  geom_boxplot() +
  labs(title = "Fraction covered by inversions",
       x="",
       y = "fraction",
       fill = "Species") +
  theme_bw()

ggplot(summary_all, aes(x=source,y = inv_cov_len / total_seq_length, fill = source)) +
  geom_boxplot() +
  facet_wrap(~minlen) +
  labs(title = "Fraction covered by inversions",
       x="",
       y = "fraction",
       fill = "Species") +
  theme_bw()

ggplot(df_1000, aes(x=source,y = genes_on_inv/total_genes, fill = source)) +
  geom_boxplot() +
  labs(title = "Fraction of Genes covered by inversions",
       x="",
       y = "fraction",
       fill = "Species") +
  theme_bw()

df_1000 %>%
  select(source, genes_on_inv_pos, genes_on_inv_neg) %>%
  tidyr::pivot_longer(cols = c("genes_on_inv_pos", "genes_on_inv_neg"),
                      names_to = "strand", values_to = "count") %>%
  ggplot(aes(x = source, y = count, fill = strand)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Strand bias of genes in inversions",
       y = "# Genes",
       fill = "Strand") +
  theme_bw()


ggplot(summary_stats, aes(x = minlen, y = mean_avg_len, color = source)) +
  geom_line() + geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_avg_len - sd_avg_len,
                    ymax = mean_avg_len + sd_avg_len),
                width = 300) +
  theme_minimal() +
  labs(title = "Average inversion length vs min_len",
       x = "min_len threshold", y = "Mean ± SD inversion length")


ggplot(summary_stats, aes(x = minlen, y = mean_num_inv, color = source)) +
  geom_line() + geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_num_inv - sd_num_inv,
                    ymax = mean_num_inv + sd_num_inv),
                width = 300) +
  theme_minimal() +
  labs(title = "Average number of inversions vs min_len",
       x = "min_len threshold", y = "Mean ± SD number of inversions")

ggplot(summary_stats, aes(x = minlen, y = mean_frac_cov, color = source)) +
  geom_line() + geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_frac_cov - sd_frac_cov,
                    ymax = mean_frac_cov + sd_frac_cov),
                width = 300) +
  theme_minimal() +
  labs(title = "Fraction of sequence covered by inversions vs min_len",
       x = "min_len threshold", y = "Fraction covered by inversions")




# Size vs number of inversions
ggplot(df_1000, aes(x = total_seq_length, y = num_inversions, color = source)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Total sequence length vs Number of inversions (min_len=1000)",
       x = "Total sequence length", y = "Number of inversions")

# Size vs fraction covered
ggplot(df_1000, aes(x = total_seq_length, y = inv_cov_len / total_seq_length, color = source)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Total sequence length vs Fraction covered (min_len=1000)",
       x = "Total sequence length", y = "Fraction covered")

# Number of inversions vs Fraction covered
ggplot(df_1000, aes(x = num_inversions, y = inv_cov_len / total_seq_length, color = source)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Number of inversions vs Fraction covered (min_len=1000)",
       x = "Number of inversions", y = "Fraction covered")

# add sister gene features
AW_sister <- fread("/local/storage/kav67/jays/inversion_summary/AW_sisters.tsv")
AI_sister <- fread("/local/storage/kav67/jays/inversion_summary/AI_sisters.tsv")
AC_sister <- fread("/local/storage/kav67/jays/inversion_summary/AC_sisters.tsv")

sisters_all <- rbindlist(list(AW_sister,AI_sister,AC_sister), use.names = TRUE, fill = TRUE)

df_1000<-merge(df_1000,sisters_all,by="sample")

df_plot1 <- df_1000 %>%
  mutate(group_frac = group_size / total_genes,
         group_frac_on_inv = group_size / genes_on_inv)

df_biggest <- df_1000 %>%
  group_by(sample) %>%
  slice_max(order_by = group_size, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(group_frac = group_size / total_genes,
         group_frac_on_inv = group_size / genes_on_inv)


ggplot(df_biggest, aes(x = source, y = group_frac)) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkblue") +
  labs(x = "Sample",
       y = "Fraction of genes in group (group_size / total_genes)",
       title = "Paralog group size relative to total genes per sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(df_biggest, aes(x = source, y = group_frac,color = source)) +
  geom_boxplot(fill = "lightblue", alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "darkblue") +
  labs(x = "Source",
       y = "Fraction of genes in largest paralog group",
       title = "Largest paralog group relative to total genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


df_plot2 <- df_1000 %>%
  group_by(sample,source) %>%
  summarise(num_groups = n(), .groups = "drop")

ggplot(df_plot2, aes(x = source, y = num_groups)) +
  geom_boxplot(fill = "lightblue",alpha = 0.6) +
  labs(x = "Sample",
       y = "Number of paralog groups",
       title = "Number of paralog groups per sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
