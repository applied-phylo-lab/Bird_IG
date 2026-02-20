library(tidyverse)
library(dplyr)
library(ggplot2)
library(viridis)

# -----------------------------
# Parameters
# -----------------------------
tsv <- "/local/storage/kav67/within_species/Songbirds/inversion_analysis/inversion_presence.tsv"
#outdir <- "/local/storage/kav67/within_species/Songbirds/inversion_analysis/inversion_plots"
dir.create(outdir, showWarnings = FALSE)

# -----------------------------
# Read data
# -----------------------------
df <- read_tsv(tsv, show_col_types = FALSE)



support_df <- df %>%
  group_by(Species, Inversion, RefStart, RefEnd) %>%
  summarise(
    n_haplotypes = n(),
    n_support = sum(Present),
    support_frac = n_support / n_haplotypes,
    inv_len = RefEnd - RefStart,
    .groups = "drop"
  )
support_df<-unique(support_df)

species_list <- unique(support_df$Species)


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
        color = support_frac
      ),
      linewidth = 0.75,
      lineend = "round"
    ) +
    
    scale_color_viridis_c(
      name = "Support",
      limits = c(0, 1),
      labels = scales::percent
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
    
    theme_minimal(base_size = 13)
  
  # toggle legend
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  p
}



plot_ac<-plot_inversion_dotplot(support_df,"A.coerulescensAC","AC",show_y_label = TRUE)
plot_ai<-plot_inversion_dotplot(support_df,"A.insularisAI","AI")
plot_aw<-plot_inversion_dotplot(support_df,"A.woodhouseiiAW","AW",show_x_label = TRUE,show_y_label = TRUE)
plot_hf<-plot_inversion_dotplot(support_df,"house_finches","House Finches",show_legend = TRUE,show_x_label = TRUE)

(plot_ac+plot_ai)/(plot_aw+plot_hf)





support_df <- support_df %>%
  mutate(
    support_class = case_when(
      support_frac == 1 ~ "Fixed (100%)",
      support_frac >= 0.75 ~ "High (≥75%)",
      support_frac >= 0.25 ~ "Intermediate",
      TRUE ~ "Rare (<25%)"
    )
  )

ggplot(
  support_df,
  aes(
    x = support_class,
    y = inv_len,
    fill = support_class
  )
) +
  geom_boxplot(outlier.alpha = 0.4) +
  facet_wrap(~ Species, scales = "free_y") +
  scale_fill_viridis_d(
    direction =-1,
    option = "viridis",
    end = 1,
    begin=0.25
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = "Support class",
    y = "Inversion length (bp)",
    fill = "Support class"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom"
  )



inv_counts <- support_df %>%
  mutate(
    len_bin = cut(
      inv_len,
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
  scale_fill_viridis_d(direction =-1,
                       option = "viridis",
                       end = 1,
                       begin=0.25) +
  labs(
    x = "Inversion length",
    y = "Number of inversions",
    fill = "Support class"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


