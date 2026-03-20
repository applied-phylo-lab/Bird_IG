library(purrr)
library(dplyr)
library(readr)
dir <- "/local/storage/kav67/clean_birds/"

files <- list.files(
  path = file.path(dir),
  pattern = paste0("summary.csv"),
  recursive = TRUE,
  full.names = TRUE
)



df <- files %>%
  set_names() %>%
  map_dfr(~ {
    df_temp <- read_csv(.x)
    if (nrow(df_temp) == 0) return(NULL)
    df_temp %>% mutate(File = .x)
  }, .id = NULL)

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

write_csv(df,paste0(dir,"/summary_features_old.csv"))
