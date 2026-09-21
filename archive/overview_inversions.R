summary<-fread("/local/storage/kav67/jays/inversion_summary/AW_summary.tsv")

ggplot(summary, aes(x = minlen, y = num_inversions, color = sample)) +
  geom_line() +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "D") +
  theme_minimal() +
  labs(title = "Number of inversions vs. min_len",
       x = "Minimum inversion length",
       y = "Number of inversions")

ggplot(summary, aes(x = minlen, y = avg_inv_len, color = sample)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = avg_inv_len - sd_inv_len,
                    ymax = avg_inv_len + sd_inv_len),
                width = 200) +
  scale_color_viridis(discrete = TRUE, option = "C") +
  theme_minimal() +
  labs(title = "Average inversion length vs. min_len",
       x = "Minimum inversion length",
       y = "Average inversion length (bp)")


summary <- summary %>%
  mutate(frac_cov = inv_cov_len / total_seq_length)

ggplot(summary, aes(x = minlen, y = frac_cov, color = sample)) +
  geom_line() +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "B") +
  theme_minimal() +
  labs(title = "Fraction of sequence covered by diagonal inversions",
       x = "Minimum inversion length",
       y = "Fraction covered")

