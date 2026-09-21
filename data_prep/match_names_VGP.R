library(dplyr)
library(stringr)

vgp_table<-fread("/local/storage/kav67/VGP_phase1_copy062025.tsv")

folder<-"/local/storage/kav67/IG_annotation_VGP2026/"
all_species<-fread(paste0(folder,"all_species_stats_pruned_12052025.csv"))


# Prepare vgp_table with a normalized LatinName-style column
vgp_clean <- vgp_table %>%
  select(
    `Scientific Name`,
    `Assembly ID`,
    `Assembly ID main haplotype`,
    `Accession # for main haplotype`,
    `RefSeq annotation main haplotype`
  ) %>%
  mutate(
    latin_key = str_replace_all(str_to_lower(`Scientific Name`), " ", "_")
  )

# --- Match 1: TreeID vs Accession # for main haplotype ---
match1 <- all_species %>%
  left_join(
    vgp_clean %>% select(`Assembly ID`, `Assembly ID main haplotype`, `Accession # for main haplotype`) %>%
      dplyr::rename(AssemblyID_m1 = `Assembly ID`, AssemblyID_main_m1 = `Assembly ID main haplotype`),
    by = c("TreeID" = "Accession # for main haplotype")
  )

# --- Match 2: TreeID vs RefSeq annotation main haplotype ---
match2 <- all_species %>%
  left_join(
    vgp_clean %>% select(`Assembly ID`, `Assembly ID main haplotype`, `RefSeq annotation main haplotype`) %>%
      dplyr::rename(AssemblyID_m2 = `Assembly ID`, AssemblyID_main_m2 = `Assembly ID main haplotype`),
    by = c("TreeID" = "RefSeq annotation main haplotype")
  )

# --- Match 3: LatinName vs normalized Scientific Name ---
match3 <- all_species %>%
  left_join(
    vgp_clean %>% select(`Assembly ID`, `Assembly ID main haplotype`, latin_key) %>%
      dplyr::rename(AssemblyID_m3 = `Assembly ID`, AssemblyID_main_m3 = `Assembly ID main haplotype`),
    by = c("LatinName" = "latin_key")
  )

# --- Coalesce all matches in priority order ---
all_species_joined <- all_species %>%
  mutate(
    AssemblyID_m1      = match1$AssemblyID_m1,
    AssemblyID_m2      = match2$AssemblyID_m2,
    AssemblyID_m3      = match3$AssemblyID_m3,
    AssemblyID_main_m1 = match1$AssemblyID_main_m1,
    AssemblyID_main_m2 = match2$AssemblyID_main_m2,
    AssemblyID_main_m3 = match3$AssemblyID_main_m3,
    `Assembly ID`               = coalesce(AssemblyID_m1, AssemblyID_m2, AssemblyID_m3),
    `Assembly ID main haplotype` = coalesce(AssemblyID_main_m1, AssemblyID_main_m2, AssemblyID_main_m3)
  ) %>%
  select(-AssemblyID_m1, -AssemblyID_m2, -AssemblyID_m3,
         -AssemblyID_main_m1, -AssemblyID_main_m2, -AssemblyID_main_m3)

# --- Diagnostics ---
cat("Total rows:", nrow(all_species_joined), "\n")
cat("Assembly ID matched via Accession:", sum(!is.na(match1$AssemblyID_m1)), "\n")
cat("Assembly ID matched via RefSeq:   ", sum(!is.na(match2$AssemblyID_m2)), "\n")
cat("Assembly ID matched via LatinName:", sum(!is.na(match3$AssemblyID_m3)), "\n")
cat("Assembly ID still unmatched:      ", sum(is.na(all_species_joined$`Assembly ID`)), "\n")
cat("Main haplotype still unmatched:   ", sum(is.na(all_species_joined$`Assembly ID main haplotype`)), "\n")

write_csv(all_species_joined,paste0(folder,"all_species_stats_pruned_12052025.csv"))


files<-list.files(paste0(folder,"annotation_stats/"))
library(dplyr)
library(stringr)
library(purrr)


# Function to normalize a string to bare lowercase letters and underscores only
normalize_key <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("-", "_") %>%   # dashes to underscores
    str_replace_all("'", "") %>%    # remove apostrophes
    str_replace_all("\\s+", "_")    # spaces to underscores
}

vgp_lookup <- vgp_table %>%
  select(`English Name`, `Scientific Name`, `Assembly ID`, `Assembly ID main haplotype`) %>%
  mutate(
    key_english = normalize_key(`English Name`),
    key_latin   = normalize_key(`Scientific Name`)
  )
vgp_lookup[vgp_lookup$key_english=="bolin’s_lanternfish",]<-"bolins_lanternfish"
vgp_lookup[vgp_lookup$key_english=="japanese_puffer_(torafugu)",]<-"japanese_puffer"
vgp_lookup[vgp_lookup$key_english=="boeseman’s_rainbowfish",]<-"boesemans_rainbowfish"

df_list <- list()
for (f in files) {
  df <- read.csv(paste0(folder, "annotation_stats/", f)) %>%
    select(-Path) %>%
    mutate(
      species_key = SpeciesID %>%
        normalize_key() %>%
        str_replace_all("_+", "_")  # collapse any double underscores from dash removal
    )
  
  # Join on english
  df <- df %>%
    left_join(
      vgp_lookup %>% select(key_english, `English Name`, `Scientific Name`, `Assembly ID`, `Assembly ID main haplotype`),
      by = c("species_key" = "key_english")
    ) %>%
    left_join(
      vgp_lookup %>% select(key_latin, `English Name`, `Scientific Name`, `Assembly ID`, `Assembly ID main haplotype`),
      by = c("species_key" = "key_latin"),
      suffix = c("", "_latin")
    ) %>%
    mutate(
      EnglishName     = coalesce(`English Name`, `English Name_latin`),
      LatinName       = coalesce(`Scientific Name`, `Scientific Name_latin`),
      AssemblyID      = coalesce(`Assembly ID`, `Assembly ID_latin`),
      AssemblyID_main = coalesce(`Assembly ID main haplotype`, `Assembly ID main haplotype_latin`)
    ) %>%
    select(-`English Name`, -`Scientific Name`, -`Assembly ID`, -`Assembly ID main haplotype`,
           -`English Name_latin`, -`Scientific Name_latin`, -`Assembly ID_latin`,
           -`Assembly ID main haplotype_latin`)
  
  # For still-unmatched rows, try fuzzy match by removing all underscores and comparing bare strings
  still_unmatched <- df %>% filter(is.na(EnglishName)) %>% distinct(species_key) %>% pull()
  
  if (length(still_unmatched) > 0) {
    # Strip all underscores for a bare string match
    bare_lookup <- vgp_lookup %>%
      mutate(
        bare_english = str_remove_all(key_english, "_"),
        bare_latin   = str_remove_all(key_latin, "_")
      )
    
    for (sk in still_unmatched) {
      bare_sk <- str_remove_all(sk, "_")
      match_row <- bare_lookup %>%
        filter(bare_english == bare_sk | bare_latin == bare_sk) %>%
        dplyr::slice(1)
      
      if (nrow(match_row) > 0) {
        df <- df %>%
          mutate(
            EnglishName     = if_else(species_key == sk & is.na(EnglishName),     match_row$`English Name`,                 EnglishName),
            LatinName       = if_else(species_key == sk & is.na(LatinName),       match_row$`Scientific Name`,              LatinName),
            AssemblyID      = if_else(species_key == sk & is.na(AssemblyID),      match_row$`Assembly ID`,                  AssemblyID),
            AssemblyID_main = if_else(species_key == sk & is.na(AssemblyID_main), match_row$`Assembly ID main haplotype`,   AssemblyID_main)
          )
      }
    }
  }
  
  # Print remaining unmatched
  unmatched <- df %>% filter(is.na(EnglishName)) %>% distinct(SpeciesID)
  if (nrow(unmatched) > 0) {
    cat("File:", f, "\n")
    cat("  Unmatched SpeciesIDs:\n")
    print(unmatched$SpeciesID)
    cat("\n")
  }
  
  df <- df %>% select(EnglishName,LatinName,AssemblyID,AssemblyID_main,AnnotationLevel,NumGenes,
                      ContigInfo,NumContigs,LocusFraction,Locus)
  df<-df[df$EnglishName!="NA",]
  df_list <- append(df_list, list(df))
  write.csv(df,paste0(folder, "annotation_stats/", f),row.names = FALSE, quote = FALSE)
}

df_list


