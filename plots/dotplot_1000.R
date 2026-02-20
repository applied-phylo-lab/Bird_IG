# Load libraries
library(ggplot2)
library(viridis)

# ==== Input file ====
input_file <- "/local/storage/kav67/Bird_data/Cormorants/European_Shag/bGulAri2_pri/IGH_self.tsv"
input_file <- "/local/storage/kav67/within_species/Songbirds/house_finches/bHaeMex1_pri/IGH_self.tsv"

# ==== Parameters ====
minlen <- 10000   # <-- change this threshold as needed

# ==== Read data ====
df <- read.table(input_file,
                 header = TRUE,
                 sep = "\t",
                 stringsAsFactors = FALSE,
                 comment.char = "",
                 check.names = FALSE)

# Clean up column names for easier handling
colnames(df) <- gsub("[#%+]", "", colnames(df))
df$id<-gsub("%","",df$id)
df$id<-as.numeric(df$id)
# Columns now:, name1, strand1, start1, end1, length1, name2, strand2, start2, end2, length2, id

# ==== Filter by minlen ====
df <- subset(df, length1 >= minlen & length2 >= minlen)

# ==== Adjust coordinates depending on strand ====
df$plot_start1 <- df$start1
df$plot_end1   <- df$end1
df$plot_start2 <- df$start2
df$plot_end2   <- df$end2

# Flip coords if strand is negative
df$plot_start1[df$strand1 == "-"] <- df$end1[df$strand1 == "-"]
df$plot_end1[df$strand1 == "-"]   <- df$start1[df$strand1 == "-"]

df$plot_start2[df$strand2 == "-"] <- df$end2[df$strand2 == "-"]
df$plot_end2[df$strand2 == "-"]   <- df$start2[df$strand2 == "-"]

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

print(p)
df_diagonal<-df[df$start1==df$start2&df$end1==df$end2,]
df_diagonal<-df_diagonal[df_diagonal$length1>10000,]
ggplot(df_diagonal) +
  geom_segment(aes(x = plot_start1, xend = plot_end1,
                   y = plot_start2, yend = plot_end2,
                   color = id)) +
  scale_color_viridis(option = "plasma", direction = -1) +
  coord_equal() +
  theme_minimal(base_size = 14) +
  labs(x = "Genome 1 position",
       y = "Genome 2 position",
       color = "Identity %")
