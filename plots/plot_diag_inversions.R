# Load libraries
library(ggplot2)
library(viridis)

# ==== Input file ====
# Change this to your TSV file
input_file <- "/local/storage/kav67/jays/patchworkplot/IGH/output_AW/pairwise_alignments/self_9-GCA_048174275.1_A.woodhouseii_AW_366498_pri_1.0_genomic.tsv"
minlen <- 1000   
# ==== Read data ====
df <- read.table(input_file,
                 header = TRUE,
                 sep = "\t",
                 stringsAsFactors = FALSE,
                 comment.char = "",
                 check.names = FALSE)   # <-- preserve column names

# If you want to strip weird symbols for easier handling:
clean_names <- gsub("[#%+]", "", colnames(df))
colnames(df) <- clean_names
df$id <- as.numeric(gsub("%", "", df$id))

df <- subset(df, length1 >= minlen & length2 >= minlen)

# ==== Plot ====
df$plot_start1 <- df$start1
df$plot_end1   <- df$end1
df$plot_start2 <- df$start2
df$plot_end2   <- df$end2

# Flip coords if strand is negative
df$plot_start1[df$strand1 == "-"] <- df$end1[df$strand1 == "-"]
df$plot_end1[df$strand1 == "-"]   <- df$start1[df$strand1 == "-"]

df$plot_start2[df$strand2 == "-"] <- df$end2[df$strand2 == "-"]
df$plot_end2[df$strand2 == "-"]   <- df$start2[df$strand2 == "-"]
df$inversion <- df$strand1 != df$strand2
# ==== Plot ====
p <- ggplot(df) +
  geom_segment(aes(x = plot_start1, xend = plot_end1,
                   y = plot_start2, yend = plot_end2,
                   color = id)) +
  scale_color_viridis(option = "plasma", direction = -1) +
  coord_equal() +
  theme_minimal(base_size = 14) +
  labs(x = "Genome 1 position",
       y = "Genome 2 position",
       color = "Identity %")
p


crossing_inversions <- df %>%
  filter(inversion) %>%
  filter(start1 == start2 & end1 == end2)

crossing_inversions<-rbind(crossing_inversions,df[1,])

# ==== Plot only crossing inversions ====
p_cross <- ggplot(crossing_inversions) +
  geom_segment(aes(x = plot_start1, xend = plot_end1,
                   y = plot_start2, yend = plot_end2,
                   color = id)) +
  scale_color_viridis(option = "plasma", direction = -1) +
  coord_equal() +
  theme_minimal(base_size = 14) +
  labs(x = "Genome 1 position",
       y = "Genome 2 position",
       color = "Identity %") +
  ggtitle("Inversions Crossing the Diagonal")
p_cross


# ==== Compute stats ====
total_alignments <- nrow(df)
total_inversions <- sum(df$inversion)

# lengths of inverted vs total
locus_length <- df$length1[1]
total_len_inv <- sum(df$length1[df$inversion])

inversion_stats <- data.frame(
  minlen = minlen,
  locus_length = locus_length,
  total_alignments = total_alignments,
  total_inversions = total_inversions,
  frac_inversions = total_inversions / total_alignments,
  total_len_inversions = total_len_inv,
  frac_locus_inverted = total_len_inv / locus_length
)
print(inversion_stats)
