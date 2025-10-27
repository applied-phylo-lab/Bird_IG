library(ggplot2)
library(dplyr)
library(gridExtra)  # for arranging plots
library(stringr)

inversion_details<-fread("/local/storage/kav67/birds/IGH_inversion_details.tsv")
output_pdf <- "/local/storage/kav67/birds/inversion_length_histograms.pdf"
inversion_details <- inversion_details %>%
  mutate(
    species_clean = gsub("_", " ", species),
    # Remove "genomic_igdetective" if present
    sample_clean = gsub("_genomic_igdetective$", "", haplotype),
    sample_clean = gsub("_igdetective$", "", haplotype),
    # Remove the first two underscore-separated parts (e.g., "GCA_048174505.1_A.woodhouseii_AW_366498_pri_1.0")
    sample_clean = sub("^[^_]+_[^_]+_", "", sample_clean)
  )

multi_hap_species <- inversion_details %>%
  distinct(species_clean, haplotype) %>%
  count(species_clean) %>%
  filter(n > 1) %>%
  pull(species_clean)


# Function to plot histograms for one order
plot_order_histograms <- function(order_name, df) {
  df_order <- df %>% filter(order == order_name)
  
  # Consistent scale within order
  max_length <- max(df_order$length, na.rm = TRUE)
  max_count <- df_order %>%
    split(.$sample_clean) %>%
    map_dbl(function(df) {
      # Add small buffer to cover edge cases
      h <- hist(df$length, breaks = seq(0, max_length + 1000, by = 1000), plot = FALSE)
      max(h$counts)
    }) %>%
    max(na.rm = TRUE)
  
  # Create plots per (species, haplotype)
  plots <- df_order %>%
    group_split(species_clean, haplotype) %>%
    map(~{
      species <- unique(.x$species_clean)
      haplo <- unique(.x$haplotype)
      sample <- unique(.x$sample_clean)
      title <- if (species %in% multi_hap_species) {
        paste0(species, " — \n", sample)
      } else {
        species
      }
      
      ggplot(.x, aes(x = length)) +
        geom_histogram(binwidth = 1000, fill = "steelblue", color = "black") +
        ggtitle(title) +
        xlim(0, max_length) +
        scale_y_continuous(limits = c(0, max_count), expand = c(0, 0))+
        xlab("Inversion length") + ylab("Count") +
        theme_minimal(base_size = 10) +
        theme(plot.title = element_text(size = 9, face = "bold"))
    })
  
  # Return list of pages (8 per page)
  n_pages <- ceiling(length(plots) / 8)
  pages <- map(seq_len(n_pages), ~{
    gridExtra::grid.arrange(
      grobs = plots[((.x - 1) * 8 + 1):min(.x * 8, length(plots))],
      ncol = 4, nrow = 2,
      top = paste("Order:", order_name)
    )
  })
  return(pages)
}

# Collect all pages across orders
orders <- unique(inversion_details$order)
all_pages <- map(orders, ~plot_order_histograms(.x, inversion_details)) %>% flatten()

# Write them to one multi-page PDF
pdf("/local/storage/kav67/birds/inversion_length_histograms.pdf", width = 11, height = 8.5)
for (page in all_pages) {
  grid::grid.newpage()     # ensure each layout gets a new page
  grid::grid.draw(page)
}
dev.off()

inversion_details<-inversion_details[inversion_details$diagonal==TRUE,]
all_pages <- map(orders, ~plot_order_histograms(.x, inversion_details)) %>% flatten()

# Write them to one multi-page PDF
pdf("/local/storage/kav67/birds/inversion_length_histograms_diag.pdf", width = 11, height = 8.5)
for (page in all_pages) {
  grid::grid.newpage()     # ensure each layout gets a new page
  grid::grid.draw(page)
}
dev.off()
