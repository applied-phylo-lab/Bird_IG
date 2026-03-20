library(tidyverse)

# ----------------------------
# Load data
# ----------------------------

stats <- read_tsv("/local/storage/kav67/within_species/Songbirds/inversions_shared_stats.tsv")

# Make frequency numeric just in case
stats <- stats %>%
  mutate(
    Frequency = as.numeric(Frequency),
    MeanLength = as.numeric(MeanLength),
    Midpoint = (RefStart + RefEnd) / 2
  )

# ----------------------------
# 1) DOTPLOT PER SPECIES
# ----------------------------
plot_inversion_dotplot <- function(
    sp_df,
    species_name,
    title,
    show_legend = FALSE,
    show_x_label = FALSE,
    show_y_label = FALSE
) {
  sp_df <- sp_df %>%
    filter(Species == species_name)
  
  # determine plot limits
  lims <- range(c(sp_df$RefStart, sp_df$RefEnd), na.rm = TRUE)
  
  p <- ggplot(sp_df) +
    # reference diagonal
    geom_abline(
      slope = 1,
      intercept = 0,
      color = "black",
      linewidth = 0.6
    ) +
    
    # inversion diagonals
    geom_segment(
      aes(
        x = RefStart,
        y = RefEnd,
        xend = RefEnd,
        yend = RefStart,
        color = Frequency
      ),
      linewidth = 0.75,
      lineend = "round"
    ) +
    
    scale_color_viridis_c(
      name = "Shared",
      limits = c(0, 1),
      labels = scales::percent,
      direction = -1,
    ) +
    
    scale_x_continuous(
      labels = function(x) paste0(x / 1000)
    ) +
    scale_y_continuous(
      labels = function(x) paste0(x / 1000)
    ) +
    
    coord_equal(xlim = lims, ylim = lims, expand = FALSE) +
    
    labs(
      title = title,
      x = if (show_x_label) "Reference position (kb)" else NULL,
      y = if (show_y_label) "Reference position (kb)" else NULL
    ) +
    
    theme_classic(base_size = 13)
  
  # toggle legend
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  p
}
stats_filtered<-stats[stats$Frequency>0.5,]
plot_ac<-plot_inversion_dotplot(stats,"A.coerulescensAC","Florida Scrub Jay",show_y_label = TRUE,show_x_label = TRUE)
plot_ai<-plot_inversion_dotplot(stats,"A.insularisAI","Island Scrub Jay",show_x_label = TRUE)
plot_aw<-plot_inversion_dotplot(stats,"A.woodhouseiiAW","Woodhouse Scrub Jay",show_y_label = TRUE,show_x_label = TRUE)
plot_hf<-plot_inversion_dotplot(stats,"house_finches","House Finch",show_legend = TRUE,show_x_label = TRUE)

(plot_ac+plot_ai)/(plot_aw+plot_hf)
(plot_ac+plot_ai+plot_aw+plot_hf)+
  plot_layout(ncol = 4)
# ----------------------------
# 2) LENGTH vs SHARING
# ----------------------------
support_df <- stats %>%
  mutate(
    support_class = case_when(
      Frequency >= 0.9 ~ ">90%",
      Frequency >= 0.75 ~ "75-90%",
      Frequency >= 0.5 ~ "50-75%",
      Frequency >= 0.25 ~ "25-50%",
      Frequency >= 0.1 ~ "10-25%",
      TRUE ~ "<10%"
    ),
    support_class = factor(
      support_class,
      levels = c( ">90%","75-90%","50-75%","25-50%", "10-25%","<10%"),
      ordered = TRUE
    )
  )

ggplot(
  support_df,
  aes(
    x = support_class,
    y = MeanLength,
    fill = support_class
  )
) +
  geom_boxplot(outlier.alpha = 0.4) +
  facet_wrap(~ Species, scales = "free_y") +
  scale_fill_viridis_d(
    direction =1,
    option = "viridis",
    end = 1,
    begin=0.25
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = "Support class",
    y = "Inversion length (bp)",
    fill = "% Shared"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom"
  )

ggplot(
  support_df,
  aes(
    x = support_class,
    y = MeanLength,
    fill = support_class
  )
) +
  geom_boxplot(outlier.alpha = 0.4) +
  scale_fill_viridis_d(
    direction =1,
    option = "viridis",
    end = 1,
    begin=0.25
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = "Support class",
    y = "Inversion length (bp)",
    fill = "Shared"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom"
  )+scale_y_log10()

inv_counts <- support_df %>%
  mutate(
    len_bin = cut(
      MeanLength,
      breaks = c(0, 1e3, 2e3, 5e3, 1e4, Inf),
      labels = c("<1 kb", "1–2 kb", "2–5 kb", "5–10 kb", ">10 kb")
    )
  ) %>%
  count(Species, support_class, len_bin, name = "n_inversions")




ggplot(
  inv_counts,
  aes(
    x = len_bin,
    y = n_inversions,
    fill = support_class
  )
) +
  geom_col(position = "stack") +
  facet_wrap(~ Species, scales = "free_y") +
  scale_fill_viridis_d(direction =1,
                       option = "viridis",
                       end = 1,
                       begin=0.25) +
  labs(
    x = "Inversion length",
    y = "Number of inversions",
    fill = "Shared"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

inv_counts <- support_df %>%
  mutate(
    len_bin = cut(
      MeanLength,
      breaks = c(0, 1e3, 2e3, 5e3, 1e4, Inf),
      labels = c("<1", "1–2", "2–5", "5–10", ">10")
    )
  ) %>%
  count(Species, support_class, len_bin, name = "n_inversions")

ggplot(
  inv_counts,
  aes(
    x = len_bin,
    y = n_inversions,
    fill = support_class
  )
) +
  geom_col(position = "stack") +
  scale_fill_viridis_d(direction =1,
                       option = "viridis",
                       end = 1,
                       begin=0.25) +
  labs(
    x = "Inversion length",
    y = "Number of inversions",
    fill = "Shared"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

inv_counts_pct <- inv_counts %>%
  group_by(len_bin) %>%
  mutate(
    pct = n_inversions / sum(n_inversions) * 100
  ) %>%
  ungroup()

length_shared_p<-ggplot(
  inv_counts_pct,
  aes(
    x = len_bin,
    y = pct,
    fill = support_class
  )
) +
  geom_col(position = "stack") +
  scale_fill_viridis_d(
    direction = 1,
    option = "viridis",
    end = 1,
    begin = 0.25
  ) +
  labs(
    x = "Inversion length (kb)",
    y = "Percentage of inversions (%)",
    fill = "Shared"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right"
  )


plot_ac+plot_ai+plot_aw+plot_hf+
  length_shared_p+
  plot_layout(ncol = 5)

((plot_ac+plot_ai)/(plot_aw+plot_hf))+length_shared_p
design <- "
ABX
CDX
"

(
  plot_ac + plot_ai +
    plot_aw + plot_hf + length_shared_p+plot_spacer()
) +
  plot_layout(design = design)
