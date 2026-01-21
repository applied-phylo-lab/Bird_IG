dir <- "/local/storage/kav67/within_species/"

files <- list.files(
  path = file.path(dir),
  pattern = paste0("summary.csv"),
  recursive = TRUE,
  full.names = TRUE
)



df <- files %>%
  set_names() %>%  # keep filenames as names
  map_dfr(~ read_csv(.x) %>% mutate(File = .x), .id = NULL)

# Extract Order, Species, and Haplotype from file path
df <- df %>%
  mutate(
    File = gsub(dir,"",File),
    parts = strsplit(File, "/"),
    Order = sapply(parts, `[`, 2),
    Species = sapply(parts, `[`, 3),
    Haplotype = sapply(parts, `[`, 4)
  ) %>%
  select(-parts)%>%
  select(Order,Species,Haplotype,Locus,Contig,Length,NumV,NumProdV,FracProdV)

write_csv(df,paste0(dir,"/summary_features.csv"))
