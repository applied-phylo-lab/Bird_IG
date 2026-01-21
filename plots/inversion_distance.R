library(tidyverse)

# Set directory
dir <- "/local/storage/kav67/Bird_data/"

# Find all IGH_self.tsv files recursively
files <- list.files(
  path = dir,
  pattern = "^IGH_self\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

# Read and combine all files
df <- files %>%
  set_names() %>%  # keep filenames if you want later
  map_dfr(~ read_tsv(.x, show_col_types = FALSE),
          .id = "source_file")

# Compute distance = start2+ - end1
df <- df %>%
  mutate(
    distance = abs(`start2+` - start1)
  )

df <- df %>%
  mutate(
    Order = sub(
      "^/local/storage/kav67/Bird_data/+([^/]+)/.*",
      "\\1",
      source_file
    )
  )

df$`id%`

# Plot histogram
ggplot(df, aes(x = distance / 1000)) +
  geom_histogram(bins = 50) +
  scale_y_log10()+
  theme_minimal() +
  labs(
    title = "Distribution of distances",
    x = "Distance (kb)",
    y = "Count"
  )



ggplot(df, aes(x = distance / 1000, fill = Order)) +
  geom_histogram(
    bins = 100,
    position = "fill"
  ) +
  theme_minimal() +
  labs(
    title = "Relative Order composition across distances",
    x = "Distance (kb)",
    y = "Proportion"
  )

df$`id%`<-as.numeric(gsub("%","",df$`id%`))

ggplot(df, aes(x = distance / 1000, fill = `id%`)) +
  geom_histogram(
    bins = 100,
    position = "fill"
  ) +
  theme_minimal() +
  labs(
    title = "Relative Order composition across distances",
    x = "Distance (kb)",
    y = "Proportion"
  )

ggplot(df, aes(x = distance / 1000, y = `id%`, color = Order)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +  # log scale often makes sense for distances
  theme_minimal() +
  labs(
    title = "Relationship between distance and identity",
    x = "Distance (kb, log scale)",
    y = "Identity (%)"
  )


library(ggplot2)
library(hexbin)

ggplot(df[df$Order!="Cranes",], aes(x = distance / 1000, y = `id%`)) +
  geom_hex(bins = 50) +
  scale_fill_viridis_c() +
  scale_x_log10() +
  theme_minimal() +
  labs(
    title = "Distance vs Identity density",
    x = "Distance (kb, log scale)",
    y = "Identity (%)",
    fill = "Count"
  )

ggplot(df[df$Order!="Cranes",], aes(x = distance , y = `id%`)) +
  geom_hex(bins = 50) +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(
    title = "Distance vs Identity density",
    x = "Distance (bb)",
    y = "Identity (%)",
    fill = "Count"
  )

library(scales) 

