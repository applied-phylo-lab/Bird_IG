library(tidytree)
library(ape)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(ggtree)
library(ggplot2)
library(ggstance)
library(data.table)
library(readr)
dir<-"/local/storage/kav67/birds/"
dir<-"/local/storage/kav67/mammals/"
summary_table<-fread(paste0(dir,"summary_features.csv"))
summary_table <- summary_table %>%
  mutate(AssemblyID = str_extract(Haplotype, "^[^_]+"))  # take text before first underscore
VGP_tree<-read.tree("/local/storage/kav67/roadies_v1.1.4.nwk")
vgp_tree<-read.tree("/local/storage/kav67/VGP_figure/vgp_species_tree_08082025.nwk")
VGP_data<-fread("/local/storage/kav67/VGP_phase1_copy062025.tsv")

lookup <- VGP_data %>%
  select(
    GCA = `Accession # for main haplotype`,
    GCF = `RefSeq annotation main haplotype`,
    ScientificName = `Scientific Name`
  ) %>%
  pivot_longer(cols = c(GCA, GCF), names_to = "Type", values_to = "Assembly") %>%
  filter(!is.na(Assembly)) %>%
  select(Assembly, ScientificName)

# Make sure IDs match exactly (trim spaces)
lookup$Assembly <- str_trim(lookup$Assembly)
VGP_tree$tip.label <- str_trim(VGP_tree$tip.label)

# Rename tip labels based on lookup
new_labels <- sapply(VGP_tree$tip.label, function(x) {
  match_row <- lookup %>% filter(Assembly == x)
  if (nrow(match_row) == 1) {
    match_row$ScientificName
  } else {
    x  # keep original if no match
  }
})

VGP_tree$tip.label <- new_labels



summary_table <- summary_table %>%
  left_join(
    VGP_data %>%
      select(`Assembly ID`, `Scientific Name`) %>%
      rename(AssemblyID = `Assembly ID`, LatinName = `Scientific Name`),
    by = "AssemblyID"
  )

normalize_species_name <- function(x) {
  if (is.na(x) || x == "") return(NA_character_)
  # keep letters/numbers/spaces/hyphens, remove other punctuation
  x <- tolower(x)
  x <- str_replace_all(x, "[^a-z0-9\\s\\-]", "")
  x <- str_trim(x)
  # split on spaces into words, treat hyphenated parts as subwords to camel-case and join
  words <- str_split(x, "\\s+")[[1]]
  words_proc <- vapply(words, function(w) {
    parts <- str_split(w, "-", simplify = FALSE)[[1]]
    # Capitalize each subpart, then join them with no separator (Black + Legged -> BlackLegged)
    parts_cap <- str_to_title(parts)
    paste0(parts_cap, collapse = "")
  }, FUN.VALUE = character(1))
  # join the processed words with underscores: King_Vulture or BlackLegged_Kittiwake
  paste(words_proc, collapse = "_")
}

# ---- create a long mapping table: one row per normalized English name -> scientific name ----
# This handles cases where `English Name` is a comma-separated list of names.
vgp_species_map <- VGP_data %>%
  # keep Assembly ID and Scientific Name; split English Name on commas
  transmute(AssemblyID = `Assembly ID`,
            ScientificName = `Scientific Name`,
            EnglishName = `English Name`) %>%
  mutate(EnglishParts = str_split(EnglishName, ",")) %>%
  unnest_longer(EnglishParts) %>%
  mutate(EnglishParts = str_trim(EnglishParts),
         SpeciesNorm = map_chr(EnglishParts, normalize_species_name)) %>%
  select(AssemblyID, SpeciesNorm, ScientificName) %>%
  # if the same SpeciesNorm appears multiple times with different ScientificName rows,
  # keep the first—adjust as needed for your data.
  group_by(SpeciesNorm) %>%
  slice(1) %>%
  ungroup()

# ---- First: try matching by AssemblyID (extracted from Haplotype before first underscore) ----
summary_with_asm <- summary_table %>%
  mutate(AssemblyID = sub("_.*$", "", Haplotype)) %>%
  left_join(
    VGP_data %>%
      select(AssemblyID = `Assembly ID`, LatinName_by_Assembly = `Scientific Name`),
    by = "AssemblyID"
  )

# ---- Second: for any rows still missing LatinName, try matching by Species -> SpeciesNorm ----
# We assume summary_table$Species is already in the normalized style (e.g. King_Vulture, Kagu etc.)
summary_filled <- summary_with_asm %>%
  left_join(
    vgp_species_map %>%
      select(SpeciesNorm, ScientificName_by_Species = ScientificName),
    by = c("Species" = "SpeciesNorm")
  ) %>%
  # coalesce: prefer the Assembly-based LatinName, otherwise the Species-based ScientificName
  mutate(LatinName = coalesce(LatinName_by_Assembly, ScientificName_by_Species)) %>%
  # drop helper columns if you don't want them
  select(-AssemblyID, -LatinName_by_Assembly, -ScientificName_by_Species)

# ---- Quick diagnostics: which summary rows remain unmatched? ----
unmatched_species <- summary_filled %>%
  filter(is.na(LatinName) | LatinName == "") %>%
  distinct(Species) %>%
  pull(Species)

if (length(unmatched_species) > 0) {
  message("The following Species in summary_table did not get a LatinName (inspect these):")
  print(unmatched_species)
} else {
  message("All rows have LatinName filled (either via AssemblyID or Species match).")
}

summary_table_curated <- summary_filled %>%
  filter(
    !is.na(LatinName),          # keep only rows with a valid LatinName
    !str_detect(Haplotype, "_alt")  # remove alt haplotypes
  )


#bird tree
birds_vgp<-VGP_data[VGP_data$Lineage=="Birds",]
drop_tips<-VGP_data[VGP_data$Lineage!="Birds",]$`English Name`
bird_tree<-drop.tip(vgp_tree,drop_tips)
plot(bird_tree)

#other trees
tips_to_keep <- intersect(VGP_tree$tip.label, summary_filled$LatinName)
sub_tree <- drop.tip(VGP_tree, setdiff(VGP_tree$tip.label, tips_to_keep))
plot(sub_tree)

locus_counts <- summary_table_curated %>%
  group_by(LatinName, Locus) %>%
  summarise(NumLoci = n(), .groups = "drop")
locus_counts_wide <- locus_counts %>%
  pivot_wider(names_from = Locus, values_from = NumLoci, values_fill = 0) %>%
  rename(label = LatinName)
#locus_counts_wide$label<-gsub(" ","_",locus_counts_wide$label)

tips_in_data <- intersect(bird_tree$tip.label, locus_counts_wide$label)
tips_in_data <- intersect(sub_tree$tip.label, locus_counts_wide$label)

# 2. Drop tips without data
bird_tree_pruned <- drop.tip(bird_tree, setdiff(bird_tree$tip.label, tips_in_data))
sub_tree_pruned <- drop.tip(sub_tree, setdiff(sub_tree$tip.label, tips_in_data))

# 3. Subset your data to only the ones still in the tree
locus_counts_wide_sub <- locus_counts_wide[locus_counts_wide$label %in% tips_in_data, ]

# base tree (horizontal)
p <- ggtree(bird_tree_pruned, layout = "rectangular")
p <- ggtree(sub_tree_pruned, layout = "rectangular")
# IGH barplot
p1 <- facet_plot(
  p, panel = "IGH",
  data = locus_counts_wide,
  geom = geom_barh,
  mapping = aes(x = IGH),
  stat = "identity"
)

# IGL barplot underneath
p2 <- facet_plot(
  p1, panel = "IGL",
  data = locus_counts_wide,
  geom = geom_barh,
  mapping = aes(x = IGL),
  stat = "identity"
)

p2

# get only max contig for IGH locus, remove if they have 0
summary_table_IGH <- summary_table_curated %>%
  filter(Locus == "IGH") %>%            # only IGH rows
  group_by(LatinName) %>%               # group by species
  slice_max(order_by = NumV, n = 1) %>% # keep row with highest NumV
  ungroup()

write_tsv(summary_table_IGH,paste0(dir,"IGH_filtered_table.tsv"))


summary_table_IGH_all <- summary_filled %>%
  filter(Locus == "IGH") %>%            # only IGH rows
  group_by(Haplotype) %>%               # group by species
  slice_max(order_by = NumV, n = 1) %>% # keep row with highest NumV
  ungroup()

write_tsv(summary_table_IGH_all,paste0(dir,"IGH_table.tsv"))
