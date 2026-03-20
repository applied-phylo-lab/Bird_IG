library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# Load data
df <- read_csv("/local/storage/kav67/within_species/Songbirds/house_finches/bHaeMex1_pri/repaver_kcount.csv",
               col_names = c("kmer",
                             "total",
                             "forward",
                             "reverse",
                             "mean_gap",
                             "median_gap",
                             "mad_gap"))

# Sort by abundance
df <- df %>% arrange(desc(total))

# ----------------------------
# 1️⃣ Abundance distribution
# ----------------------------

p1 <- ggplot(df, aes(x = total)) +
  geom_histogram(bins = 50) +
  scale_x_log10() +
  theme_bw() +
  labs(title = "K-mer Abundance Distribution",
       x = "Total Count (log scale)",
       y = "Number of K-mers")

print(p1)


# ----------------------------
# 2️⃣ Forward vs Reverse Bias
# ----------------------------

df$strand_bias <- df$forward / df$total

p2 <- ggplot(df, aes(x = strand_bias)) +
  geom_histogram(bins = 50) +
  theme_bw() +
  labs(title = "Strand Bias Distribution",
       x = "Forward / Total",
       y = "Number of K-mers")

print(p2)


# ----------------------------
# 3️⃣ Top 20 Most Abundant Kmers
# ----------------------------

top20 <- df %>% slice_max(total, n = 20)

top20_long <- top20 %>%
  select(kmer, forward, reverse) %>%
  pivot_longer(cols = c(forward, reverse),
               names_to = "orientation",
               values_to = "count")

p3 <- ggplot(top20_long,
             aes(x = reorder(kmer, count),
                 y = count,
                 fill = orientation)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  theme_bw() +
  labs(title = "Top 20 Most Abundant K-mers",
       x = "K-mer",
       y = "Count")

print(p3)


# ----------------------------
# 4️⃣ Abundance vs Regularity
# ----------------------------

p4 <- ggplot(df, aes(x = total, y = median_gap)) +
  geom_point(alpha = 0.6) +
  scale_x_log10() +
  theme_bw() +
  labs(title = "Repeat Regularity vs Abundance",
       x = "Total Count (log scale)",
       y = "Median of Gap Length")

print(p4)

library(dplyr)
library(ggplot2)

# Top 10% by abundance
threshold <- quantile(df$total, 0.90)

high_df <- df %>% filter(total >= threshold)

ggplot(high_df, aes(x = strand_bias)) +
  geom_histogram(bins = 30) +
  theme_bw() +
  labs(title = "Strand Bias (Top 10% Most Abundant 21-mers)",
       x = "Forward / Total",
       y = "Count")

ggplot(df, aes(x = total, y = strand_bias)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  theme_bw() +
  labs(title = "Abundance vs Strand Bias (IGH, k=21)",
       x = "Total Count (log scale)",
       y = "Forward / Total")

ggplot(df, aes(x = strand_bias, y = mad_gap)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  labs(title = "Strand Bias vs Gap Irregularity",
       x = "Forward / Total",
       y = "MAD of Gap Length")

high_df <- df %>% filter(total > 40)

ggplot(high_df, aes(x = mean_gap)) +
  geom_histogram(bins = 30) +
  theme_bw() +
  labs(title = "Mean Gap Distribution (High Abundance 21-mers)",
       x = "Mean Gap Length",
       y = "Count")

