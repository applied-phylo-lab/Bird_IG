#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)
library(patchwork)
# --- Load data ---
dot <- read_tsv("/local/storage/kav67/within_species_updated/Songbirds/_dotplot.tsv")


# -------------------------------------------------------------------
# Function: make one dotplot for a single species
# -------------------------------------------------------------------
plot_inversion_dotplot <- function(df, species_name,title,
                                   show_legend = FALSE,
                                   show_x_label = FALSE,
                                   show_y_label = FALSE) {
  
  sp<- species_name
  df <- df %>%
    filter(Species == species_name)
  
  n_haps <- df$TotalHaps[1]
  print(n_haps)
  # Axis limits for the true diagonal
  ax_min <- min(c(df$RefStart, df$QueryStart))
  ax_max <- max(c(df$RefEnd,   df$QueryEnd))
  
  p<-ggplot(df) +
    
    # True diagonal for orientation
    geom_segment(
      aes(x = ax_min, xend = ax_max, y = ax_min, yend = ax_max),
      color     = "black",
      linewidth = 0.4,
      linetype  = "solid"
    ) +
    
    # Inversions: swap QueryEnd -> y, QueryStart -> yend so they go anti-diagonal
    geom_segment(
      aes(
        x     = RefStart,
        xend  = RefEnd,
        y     = QueryEnd,    # swapped
        yend  = QueryStart,  # swapped
        color = Frequency
      ),
      linewidth = 0.5,
      alpha     = 0.85
    ) +
    
    scale_color_viridis_c(
      option    = "viridis",
      direction = -1,
      name      = "Frequency",
      limits    = c(0, 1),
      breaks    = c(0, 0.25, 0.5, 0.75, 1),
      labels    = c("0\n(private)", "0.25", "0.5", "0.75", "1\n(all haplotypes)")
    ) +
    
    labs(
      title    = title,
      subtitle = paste0(n_haps, " haplotypes  |  ", nrow(df), " inversion events"),
      x = if (show_x_label) "Reference position (kb)" else NULL,
      y = if (show_y_label) "Reference position (kb)" else NULL
    ) +
    
    theme_classic(base_size = 12) +
    theme(
      legend.position   = "right",
      legend.key.height = unit(1.5, "cm"),
      plot.title        = element_text(face = "bold", size = 14),
      plot.subtitle     = element_text(color = "grey40", size = 9)
    )+
    scale_x_continuous(
      labels = function(x) paste0(x / 1000)
    ) +
    scale_y_continuous(
      labels = function(x) paste0(x / 1000)
    ) 
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  return(p)
}


plot_inversion_dotplot_mirrored <- function(df, species_name, title,
                                            show_legend = FALSE,
                                            show_x_label = FALSE,
                                            show_y_label = FALSE) {
  
  sp <- species_name
  df <- df %>%
    dplyr::filter(Species == species_name)
  
  n_haps <- length(unique(df$RepKey))
  print(n_haps)
  
  # Axis limits for the diagonal
  ax_min <- min(c(df$RefStart, df$QueryStart))
  ax_max <- max(c(df$RefEnd,   df$QueryEnd))
  
  #df[df$NumSupportingHaps==1,]$Frequency<-0
  df<-df[df$NumSupportingHaps!=1,]
  #df<-df[df$MeanLength>10000,]
  
  p <- ggplot(df) +
    
    # Mirrored diagonal (top-left to bottom-right after reversing y)
    geom_segment(
      aes(x = ax_min, xend = ax_max, y = ax_min, yend = ax_max),
      color     = "black",
      linewidth = 0.4,
      linetype  = "solid"
    ) +
    
    # Inversions (same mapping — flipping axis handles mirroring)
    geom_segment(
      aes(
        x     = RefStart,
        xend  = RefEnd,
        y     = QueryEnd,
        yend  = QueryStart,
        color = Frequency
      ),
      linewidth = 0.5,
      alpha     = 0.85
    ) +
    
    scale_color_viridis_c(
      option    = "viridis",
      direction = -1,
      name      = "Shared",
      limits    = c(0, 1),
      breaks    = c(0, 0.25, 0.5, 0.75, 1),
      labels    = c("0\n(private)", "0.25", "0.5", "0.75", "1\n(all haplotypes)")
    ) +
    
    labs(
      title    = title,
      subtitle = paste0(n_haps, " haplotypes  |  ", nrow(df), " inversion events"),
      x = if (show_x_label) "Reference position (kb)" else NULL,
      y = if (show_y_label) "Reference position (kb)" else NULL
    ) +
    
    theme_classic(base_size = 12) +
    theme(
      legend.position   = "right",
      legend.key.height = unit(1, "cm"),
      plot.title        = element_text(face = "bold", size = 14),
      plot.subtitle     = element_text(color = "grey40", size = 9)
    ) +
    
    scale_x_continuous(
      limits = c(ax_min, ax_max),
      labels = function(x) paste0(x / 1000)
    ) +
    
    # 🔑 THIS is the key change
    scale_y_reverse(
      limits = c(ax_max, ax_min),  # reversed limits
      labels = function(x) paste0(x / 1000)
    )
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}

# -------------------------------------------------------------------
# Build a named list of plots, one per species
# -------------------------------------------------------------------
dot_filtered<-dot#[dot$MeanLength>5000,]

plot_ac<-plot_inversion_dotplot_mirrored(dot_filtered,"Florida_scrub_jay","Florida Scrub Jay",show_y_label = TRUE,show_x_label = TRUE)
plot_ai<-plot_inversion_dotplot_mirrored(dot_filtered,"Island_scrub_jay","Island Scrub Jay",show_x_label = TRUE)
plot_aw<-plot_inversion_dotplot_mirrored(dot_filtered,"Woodhouse_scrub_jay","Woodhouse Scrub Jay",show_x_label = TRUE)
plot_hf<-plot_inversion_dotplot_mirrored(dot_filtered,"house_finches","House Finch",show_y_label = TRUE,show_x_label = TRUE)
plot_c<-plot_inversion_dotplot_mirrored(dot_filtered,"Chestnut_seedeater","Chestnut Seedeater",show_x_label = TRUE)
plot_d<-plot_inversion_dotplot_mirrored(dot_filtered,"Dark_throated_seedeater","Dark Throated Seedeater",show_x_label = TRUE,show_legend = TRUE)
plot_i<-plot_inversion_dotplot_mirrored(dot_filtered,"Ibera_seedeater","Ibera Seedeater",show_y_label = TRUE,show_x_label = TRUE)
plot_m<-plot_inversion_dotplot_mirrored(dot_filtered,"Marsh_seedeater","Marsh Seedeater",show_x_label = TRUE)
plot_tb<-plot_inversion_dotplot_mirrored(dot_filtered,"Tawny_bellied_seedeater","Tawny Bellied Seedeater",show_x_label = TRUE)

#(plot_ac+plot_ai)/(plot_aw+plot_hf)
#(plot_c+plot_d)/(plot_i+plot_m+plot_tb)

(plot_ac+plot_ai+plot_aw+plot_hf+plot_c+plot_d+plot_i+plot_m+plot_tb)+plot_layout(ncol = 3)




## State Wise and Coast Wise
plot_inversion_dotplot_geo <- function(df, species_name, title,
                                       group_var = NULL, group_value = NULL,
                                       show_legend = FALSE,
                                       show_x_label = FALSE,
                                       show_y_label = FALSE) {
  
  # ---------------------------------------------------------------------------
  # 1. Filter species
  # ---------------------------------------------------------------------------
  d <- df %>% dplyr::filter(Species == species_name)
  
  # ---------------------------------------------------------------------------
  # 2. Optional: filter by geography (State or Coast)
  # ---------------------------------------------------------------------------
  if (!is.null(group_var) && !is.null(group_value)) {
    d <- d %>% dplyr::filter(.data[[group_var]] == group_value)
  }
  
  # ---------------------------------------------------------------------------
  # 3. Recompute haplotype counts within this subset
  # ---------------------------------------------------------------------------
  n_haps <- n_distinct(d$hap)
  
  d_sum <- d %>%
    group_by(Inversion) %>%
    summarise(
      RefStart = first(RefStart),
      RefEnd   = first(RefEnd),
      QueryStart = first(QueryStart),
      QueryEnd   = first(QueryEnd),
      NumSupportingHaps = n_distinct(hap),
      .groups = "drop"
    ) %>%
    mutate(
      Frequency = NumSupportingHaps / n_haps
    )
  
  print(n_haps)
  
  # ---------------------------------------------------------------------------
  # 4. Axis limits
  # ---------------------------------------------------------------------------
  ax_min <- min(c(d_sum$RefStart, d_sum$QueryStart))
  ax_max <- max(c(d_sum$RefEnd,   d_sum$QueryEnd))
  
  # ---------------------------------------------------------------------------
  # 5. Plot
  # ---------------------------------------------------------------------------
  p <- ggplot(d_sum) +
    
    geom_segment(
      aes(x = ax_min, xend = ax_max, y = ax_min, yend = ax_max),
      color = "black",
      linewidth = 0.4
    ) +
    
    geom_segment(
      aes(
        x     = RefStart,
        xend  = RefEnd,
        y     = QueryEnd,
        yend  = QueryStart,
        color = Frequency
      ),
      linewidth = 0.5,
      alpha     = 0.85
    ) +
    
    scale_color_viridis_c(
      option    = "viridis",
      direction = -1,
      name      = "Shared",
      limits    = c(0, 1),
      breaks    = c(0, 0.25, 0.5, 0.75, 1),
      labels    = c("0\n(private)", "0.25", "0.5", "0.75", "1\n(all haplotypes)")
    ) +
    
    labs(
      title    = title,
      subtitle = paste0(n_haps, " haplotypes  |  ", nrow(d_sum), " inversions"),
      x = if (show_x_label) "Reference position (kb)" else NULL,
      y = if (show_y_label) "Reference position (kb)" else NULL
    ) +
    
    theme_classic(base_size = 12) +
    theme(
      legend.position   = ifelse(show_legend, "right", "none"),
      legend.key.height = unit(1, "cm"),
      plot.title        = element_text(face = "bold", size = 14),
      plot.subtitle     = element_text(color = "grey40", size = 9)
    ) +
    
    scale_x_continuous(
      limits = c(ax_min, ax_max),
      labels = function(x) paste0(x / 1000)
    ) +
    
    scale_y_reverse(
      limits = c(ax_max, ax_min),
      labels = function(x) paste0(x / 1000)
    )
  
  return(p)
}

# -----------------------------------------------------------------------------
hf_df <- stats[stats$Species=="house_finches",] %>%
  separate_rows(SupportingKeys, sep = ",") %>%
  
  # -----------------------------------------------------------------------------
# 2. Extract haplotype ID (first 3 underscore-separated fields)
# -----------------------------------------------------------------------------
mutate(
  hap = str_extract(SupportingKeys, "^[^_]+_[^_]+_[^_]+")
) %>%
  
  # -----------------------------------------------------------------------------
# 3. Remove duplicate hap–inversion combinations
# (same hap may appear multiple times due to multiple contigs)
# -----------------------------------------------------------------------------
distinct(Species, Inversion, hap, .keep_all = TRUE)



hf_df <- hf_df %>%
  mutate(
    State = case_when(
      grepl("^bHaeMex", SupportingKeys) ~ "CA",
      TRUE ~ str_extract(SupportingKeys, "^[A-Z]+")
    ),
    Coast = case_when(
      State %in% c("CA", "WA", "AZ","NM") ~ "West",
      State %in% c("AL", "MA", "NY", "OH") ~ "East",
      TRUE ~ NA_character_
    ),
    RepKey = case_when(
      grepl("^bHaeMex", SupportingKeys) ~ "CA_3",
      TRUE ~ str_extract(SupportingKeys, "^[A-Z]+_[1-3]_[pri|alt]")
    )
  )

dot_coords <- dot %>%
  select(Species, Inversion, RefStart, RefEnd, QueryStart, QueryEnd) %>%
  distinct()

hf_df <- hf_df %>%
  left_join(dot_coords, by = c("Species", "Inversion"))

coasts <- unique(hf_df$Coast)

coast_plots <- lapply(coasts, function(cs) {
  plot_inversion_dotplot_geo(
    hf_df,
    species_name = "house_finches",
    title = paste("Coast:", cs),
    group_var = "Coast",
    group_value = cs
  )
})


names(coast_plots) <- coasts
coast_plots$West+coast_plots$East

states <- unique(hf_df$State)

state_plots <- lapply(states, function(cs) {
  df_sub <- hf_df %>% filter(State == cs)
  
  plot_inversion_dotplot_mirrored_cs(
    df_sub,
    "house_finches",
    title = cs
  )
})

names(state_plots) <- states
(state_plots$WA+state_plots$NY)/
  (state_plots$CA+state_plots$MA)/
  (state_plots$AZ+state_plots$OH)/
  (state_plots$NM+state_plots$AL)


test_dot<-hf_df%>% filter(
  hap %in% c("CA_1_a","OH_2_p","NM_1_p","AL_2_p")
)

random<-plot_inversion_dotplot_mirrored_cs(
  test_dot,
  "house_finches",
  title = "random"
)

test_dot<-hf_df%>% filter(
  hap %in% c("CA_1_a","WA_1_p","NM_1_p","AZ_1_p")
)

west_4<-plot_inversion_dotplot_mirrored_cs(
  test_dot,
  "house_finches",
  title = "West"
)

random+west_4+state_plots$NM


hap_mat <- hf_df %>%
  mutate(present = 1) %>%
  distinct(hap, Inversion, .keep_all = TRUE) %>%
  select(hap, Inversion, present) %>%
  pivot_wider(names_from = Inversion, values_from = present, values_fill = 0)

hap_mat<-hap_mat[!is.na(hap_mat$hap),]

mat <- hap_mat %>%
  column_to_rownames("hap") %>%
  as.matrix()

# Jaccard distance → convert to similarity
dist_mat <- vegdist(mat, method = "jaccard")

sim_mat <- 1 - as.matrix(dist_mat)

meta <- hf_df %>%
  select(hap, State, Coast) %>%
  distinct()

sim_df <- melt(sim_mat, varnames = c("hap1", "hap2"), value.name = "similarity") %>%
  filter(hap1 != hap2)

sim_df <- sim_df %>%
  left_join(meta, by = c("hap1" = "hap")) %>%
  dplyr::rename(State1 = State, Coast1 = Coast) %>%
  left_join(meta, by = c("hap2" = "hap")) %>%
  dplyr::rename(State2 = State, Coast2 = Coast) %>%
  mutate(
    same_state = State1 == State2,
    same_coast = Coast1 == Coast2
  )

wilcox.test(similarity ~ same_state, data = sim_df)
wilcox.test(similarity ~ same_coast, data = sim_df)

ggplot(sim_df, aes(x = same_state, y = similarity)) +
  geom_boxplot(fill="steelblue") +
  labs(x = "Same state", y = "Haplotype similarity")+
  theme_classic()

ggplot(sim_df, aes(x = same_coast, y = similarity)) +
  geom_boxplot(fill="steelblue") +
  labs(x = "Same coast", y = "Haplotype similarity")+
  theme_classic()


sim_df <- sim_df %>%
  mutate(
    geo_group = case_when(
      State1 == State2 ~ "Same state",
      Coast1 == Coast2 ~ "Same coast, diff state",
      TRUE             ~ "Different coast"
    )
  )
sim_df$geo_group <- factor(
  sim_df$geo_group,
  levels = c("Same state", "Same coast, diff state", "Different coast")
)

ggplot(sim_df, aes(x = geo_group, y = similarity, fill = geo_group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
  labs(
    x = NULL,
    y = "Haplotype similarity",
    title = "Geographic structure of inversion sharing"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1)
  )
