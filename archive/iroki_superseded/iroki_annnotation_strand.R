bird_species<-gsub("_"," ",bird_tree_saved$tip.label)

library(dplyr)
library(readr)

results <- lapply(bird_species, function(species) {
  
  # Get the bird's metadata from vgp_bird_table
  bird_info <- vgp_bird_table %>%
    filter(LatinName == species)
  
  if (nrow(bird_info) == 0) {
    message("No entry found in vgp_bird_table for: ", species)
    return(data.frame(name = species, percentage_positive = NA))
  }
  
  order    <- bird_info$Order[1]
  sp       <- bird_info$Species[1]
  haplotype <- bird_info$Haplotype[1]
  contig   <- bird_info$Contig[1]
  
  # Build file path
  file_path <- file.path("/local/storage/kav67/clean_birds/", order, sp, haplotype, "combined_genes_IGH_clean.txt")
  
  if (!file.exists(file_path)) {
    message("File not found for: ", species, " -> ", file_path)
    return(data.frame(name = species, percentage_positive = NA))
  }
  
  # Read and filter
  gene_data <- read_tsv(file_path, show_col_types = FALSE) %>%
    filter(Contig == contig)
  
  if (nrow(gene_data) == 0) {
    message("No rows after Contig filter for: ", species)
    return(data.frame(name = species, percentage_positive = NA))
  }
  
  # Calculate % on + strand
  pct_positive <- sum(gene_data$Strand == "+") / nrow(gene_data) * 100
  
  data.frame(name = species, percentage_positive = pct_positive)
})

strand_results <- bind_rows(results)
strand_results$percentage_positive<-strand_results$percentage_positive*-1
strand_results$name<-gsub(" ","_",strand_results$name)
#combine strand_results and annotation_iroki


annotation_iroki_new <- annotation_iroki %>%
  select(name, 
         bar1_height, bar2_height, bar4_height = bar3_height, bar5_height = bar4_height, bar6_height = bar5_height,
         bar1_color,  bar2_color,  bar4_color  = bar3_color,  bar5_color  = bar4_color, bar6_color = bar5_color) %>%
  mutate(bar2_color="grey")%>% 
  left_join(strand_results %>% 
              select(name, bar3_height= percentage_positive) %>%
              mutate(bar3_color = "darkgrey"), 
            by = "name") %>%
  select(name, bar1_height, bar2_height, bar3_height, bar4_height, bar5_height,bar6_height,
         bar1_color,  bar2_color,  bar3_color,  bar4_color,  bar5_color,bar6_color)


write_tsv(annotation_iroki_new,"/local/storage/kav67/clean_birds/annotation_iroki_strand.tsv")
