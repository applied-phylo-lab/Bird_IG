library(tidyverse)

# -----------------------------
# Input
# -----------------------------
# directory that contains Order/Species/Haplotype/...
base_dir <- "/local/storage/kav67/within_species/Songbirds"
base_dir <- "/local/storage/kav67/Bird_data"

# -----------------------------
# Read all inversion_repeats.tsv
# -----------------------------
files <- list.files(
  base_dir,
  pattern = "inversion_repeats.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No inversion_repeats.tsv files found")
}


# -----------------------------
# Clean & derive fields
# -----------------------------

inv_rep <- files %>%
  set_names() %>%
  map_dfr(read_tsv, .id = "file") %>%
  mutate(
    Order = file %>%
      str_remove(paste0("^", normalizePath(base_dir), "/?")) %>%
      str_split("/", simplify = TRUE) %>%
      .[, 1],
    Species = file %>%
      str_remove(paste0("^", normalizePath(base_dir), "/?")) %>%
      str_split("/", simplify = TRUE) %>%
      .[, 2]
  )
inv_rep <- inv_rep %>%
  mutate(
    InversionID = paste(Contig, InvStart, InvEnd, sep = ":"),
    RepeatClassSimple = str_replace(RepeatClass, "/.*", "")
  )
#inv_rep$Species<-gsub(base_dir,"",inv_rep$Species)
# -----------------------------
# 1) Repeat classes across inversions
#    (count each inversion once per class)
# -----------------------------
inv_by_class <- inv_rep %>%
  distinct(InversionID, RepeatClass) %>%
  count(RepeatClass, sort = TRUE)

ggplot(inv_by_class, aes(x = reorder(RepeatClass, n), y = n)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    x = "Repeat class / family",
    y = "Number of inversions",
    title = "Repeat classes observed across inversions"
  ) +
  theme_minimal(base_size = 13)


if ("Species" %in% colnames(inv_rep)) {
  
  inv_species <- inv_rep %>%
    distinct(Species, InversionID, RepeatClass) %>%
    count(Species, RepeatClass)
  
  ggplot(inv_species,
         aes(x = RepeatClass, y = n, fill = RepeatClass)) +
    geom_col(show.legend = FALSE) +
    facet_wrap(~ Species, scales = "free_y") +
    coord_flip() +
    labs(
      x = "Repeat class / family",
      y = "Number of inversions",
      title = "Repeat classes across inversions by species"
    ) +
    theme_minimal(base_size = 12)
}


if ("Order" %in% colnames(inv_rep)) {
  
  inv_species <- inv_rep %>%
    distinct(Order, InversionID, RepeatClass) %>%
    count(Order, RepeatClass)
  
  ggplot(inv_species,
         aes(x = RepeatClass, y = n, fill = RepeatClass)) +
    geom_col(show.legend = FALSE) +
    facet_wrap(~ Order, scales = "free_y") +
    coord_flip() +
    labs(
      x = "Repeat class / family",
      y = "Number of inversions",
      title = "Repeat classes across inversions by species"
    ) +
    theme_minimal(base_size = 12)
}
