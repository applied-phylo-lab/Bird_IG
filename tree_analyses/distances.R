library(tidyverse)
library(ape)

# ---- INPUTS ----
INPUT_DIR  <- "/local/storage/kav67/clean_birds/"
SUMMARY    <- "/local/storage/kav67/clean_birds/IGH_filtered_table.tsv"
OUTPUT_DIR <- "/local/storage/kav67/clean_birds/trees/"
OUTPUT_FILE <- file.path(OUTPUT_DIR, "all_pairwise_distances.tsv")

# ---- READ SUMMARY ----
meta <- read_tsv(SUMMARY, show_col_types = FALSE)

# ---- FUNCTION: extract position from gene name ----
extract_pos <- function(gene_name, haplotype) {
  # Remove haplotype prefix (only at the beginning)
  stripped <- sub(paste0("^", haplotype, "[._]?"), "", gene_name)
  
  # Extract first number after prefix
  pos <- str_extract(stripped, "\\d+")
  
  return(as.numeric(pos))
}

# ---- MAIN LOOP ----
all_results <- list()

for (i in seq_len(nrow(meta))) {
  
  row <- meta[i, ]
  
  haplotype <- row$Haplotype
  species   <- row$Species
  order     <- row$Order
  
  tree_file <- file.path(INPUT_DIR, order, species, haplotype,
                         "tree", paste0(haplotype, "_tree.treefile"))
  
  if (!file.exists(tree_file)) {
    message("Skipping missing: ", tree_file)
    next
  }
  
  # ---- READ TREE ----
  tree <- tryCatch({
    read.tree(tree_file)
  }, error = function(e) {
    message("Failed to read: ", tree_file)
    return(NULL)
  })
  
  if (is.null(tree)) next
  
  # ---- GET GENETIC DISTANCES ----
  dist_mat <- cophenetic.phylo(tree)
  
  genes <- rownames(dist_mat)
  
  # ---- GET PHYSICAL POSITIONS ----
  positions <- extract_pos(genes,haplotype)
  
  # skip if positions failed
  if (all(is.na(positions))) {
    message("No valid positions: ", tree_file)
    next
  }
  
  # ---- BUILD PAIRWISE TABLE ----
  df <- as.data.frame(as.table(dist_mat)) %>%
    rename(Gene1 = Var1,
           Gene2 = Var2,
           genetic_dist = Freq) %>%
    filter(Gene1 != Gene2) %>% 
    mutate(
      pos1 = extract_pos(Gene1, haplotype),
      pos2 = extract_pos(Gene2, haplotype),
      physical_dist = abs(pos1 - pos2),
      physical_dist = abs(pos1 - pos2),
      Order = order,
      Species = species,
      Haplotype = haplotype
    ) %>%
    select(Order, Species, Haplotype,
           Gene1, Gene2,
           physical_dist, genetic_dist)
  
  all_results[[length(all_results) + 1]] <- df
}

# ---- COMBINE + WRITE ----
final_df <- bind_rows(all_results)

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

write_tsv(final_df, OUTPUT_FILE)

message("Done! Output written to: ", OUTPUT_FILE)