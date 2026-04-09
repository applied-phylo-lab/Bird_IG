IGH_df<-df %>% filter(Locus=="IGH")
IGH_df$LatinName<-""
IGH_df[IGH_df$Species=="Florida_scrub_jay",]$LatinName<-"Aphelocoma coerulescens"
IGH_df[IGH_df$Species=="Island_scrub_jay",]$LatinName<-"Aphelocoma insularis"
IGH_df[IGH_df$Species=="Woodhouse_scrub_jay",]$LatinName<-"Aphelocoma woodhousii"
IGH_df[IGH_df$Species=="house_finches",]$LatinName<-"Haemorhous mexicanus"
IGH_df[IGH_df$Species=="Chestnut_seedeater",]$LatinName<-"Sporophila cinnamomea"
IGH_df[IGH_df$Species=="Dark_throated_seedeater",]$LatinName<-"Sporophila ruficollis"
IGH_df[IGH_df$Species=="Tawny_bellied_seedeater",]$LatinName<-"Sporophila hypoxantha"
IGH_df[IGH_df$Species=="Marsh_seedeater",]$LatinName<-"Sporophila palustris"
IGH_df[IGH_df$Species=="Ibera_seedeater",]$LatinName<-"Sporophila iberaensis"
#IGH_df[IGH_df$Species=="Pearly_bellied_seedeater",]$LatinName<-"Sporophila pileata"
#IGH_df[IGH_df$Species=="Grey_and_chestnut_seedeater",]$LatinName<-"Sporophila hypochroma"

IGH_df<-IGH_df[IGH_df$LatinName!="",]

write_tsv(IGH_df,paste0(dir,"IGH_table.tsv"))

hf_df<-IGH_df %>% filter(Species=="house_finches")

hf_df <- hf_df %>%
  mutate(
    State = case_when(
      grepl("^bHaeMex", Haplotype) ~ "CA",
      TRUE ~ str_extract(Haplotype, "^[A-Z]+")
    ),
    Coast = case_when(
      State %in% c("CA", "WA", "AZ","NM") ~ "West",
      State %in% c("AL", "MA", "NY", "OH") ~ "East",
      TRUE ~ NA_character_
    )
  )
write_tsv(hf_df,paste0(dir,"IGH_table_housefinches.tsv"))
