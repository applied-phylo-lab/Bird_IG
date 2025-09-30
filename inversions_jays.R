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

# ---- plots ----
# Avg inversion length
ggplot(summary_stats, aes(x = source, y = mean_avg_len, fill = source)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mean_avg_len - sd_avg_len,
                    ymax = mean_avg_len + sd_avg_len),
                width = 0.25, position = position_dodge(0.9)) +
  theme_minimal() +
  facet_wrap(~minlen, scales = "free_y") +
  labs(title = "Average inversion length across species",
       x = "Source", y = "Mean ± SD inversion length")

# Number of inversions
ggplot(summary_stats, aes(x = source, y = mean_num_inv, fill = source)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mean_num_inv - sd_num_inv,
                    ymax = mean_num_inv + sd_num_inv),
                width = 0.25, position = position_dodge(0.9)) +
  theme_minimal() +
  facet_wrap(~minlen, scales = "free_y") +
  labs(title = "Number of inversions across species",
       x = "Source", y = "Mean ± SD number of inversions")

# Fraction covered
ggplot(summary_stats, aes(x = source, y = mean_frac_cov, fill = source)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mean_frac_cov - sd_frac_cov,
                    ymax = mean_frac_cov + sd_frac_cov),
                width = 0.25, position = position_dodge(0.9)) +
  theme_minimal() +
  facet_wrap(~minlen, scales = "free_y") +
  labs(title = "Fraction of sequence covered by inversions",
       x = "Source", y = "Mean ± SD fraction")


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
ggplot(df_1000, aes(x = total_seq_length, y = frac_cov, color = source)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Total sequence length vs Fraction covered (min_len=1000)",
       x = "Total sequence length", y = "Fraction covered")

# Number of inversions vs Fraction covered
ggplot(df_1000, aes(x = num_inversions, y = frac_cov, color = source)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Number of inversions vs Fraction covered (min_len=1000)",
       x = "Number of inversions", y = "Fraction covered")

