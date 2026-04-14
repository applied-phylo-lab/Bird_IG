#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)
library(patchwork)

# =============================================================================
# LOAD DATA
# =============================================================================
dot    <- read_tsv("/local/storage/kav67/within_species_updated/Songbirds/_dotplot.tsv")
coords <- read_tsv("/local/storage/kav67/within_species_updated/Songbirds/_coords.tsv")
stats  <- read_tsv("/local/storage/kav67/within_species_updated/Songbirds/_stats.tsv")

# =============================================================================
# SHARED REGION FILTER
# Use inversion RefStart/RefEnd in dotplot as proxy for haplotype coverage.
# Shared region = max(RefStart) to min(RefEnd) across all inversions per species.
# =============================================================================
ref_coverage <- dot %>%
  group_by(Species,RepKey) %>%
  summarise(
    shared_start = min(RefStart),
    shared_end   = max(RefEnd),
    .groups = "drop"
  )

ref_coverage_all <- ref_coverage %>%
  group_by(Species) %>%
  summarise(
    shared_start_all = max(shared_start),
    shared_end_all   = min(shared_end),
    .groups = "drop"
  )
ref_coverage_all<-ref_coverage_all%>%filter(Species!="house_finches")

dot_filtered <- dot %>%
  left_join(ref_coverage_all, by = "Species") %>%
  filter(RefStart >= shared_start_all, RefEnd <= shared_end_all) %>%
  select(-shared_start_all, -shared_end_all)

# =============================================================================
# SIZE + FREQUENCY SUMMARY TABLE
# Join stats (has Frequency + MeanLength) with ref_coverage to filter
# to inversions whose representative coords fall in the shared region.
# =============================================================================
size_freq <- stats %>%
  left_join(ref_coverage_all, by = "Species") %>%
  # keep only inversions represented in the filtered dotplot
  semi_join(dot_filtered, by = c("Species", "Inversion")) %>%
  mutate(
    SizeClass = cut(
      MeanLength,
      breaks = c(0, 1000, 2000, 5000, 10000, Inf),
      labels = c("<1 kb", "1–2 kb", "2–5 kb", "5–10 kb", ">10 kb"),
      right  = TRUE
    ),
    FreqBin = cut(
      Frequency,
      breaks = c(0, 0.1,0.25, 0.5, 0.75,0.9, 1.001),
      labels = c("<10%","10-25%", "25–50%", "50–75%","75-90%", ">90%")
    )%>%
      factor(levels = c(">90%", "75-90%", "50–75%", "25–50%", "10-25%", "<10%"),
             ordered = TRUE)
  )



size_levels <- c("<1 kb", "1–2 kb", "2–5 kb", "5–10 kb", ">10 kb")
size_freq   <- size_freq %>%
  mutate(SizeClass = factor(SizeClass, levels = size_levels))



# =============================================================================
# OPTION A — HEATMAP
# x = size class, y = frequency bin, fill = count
# =============================================================================
plot_heatmap <- function(df, sp = NULL) {
  d <- if (!is.null(sp)) filter(df, Species == sp) else df
  title <- if (!is.null(sp)) sp else "All species"
  
  d %>%
    count(SizeClass, FreqBin) %>%
    complete(SizeClass, FreqBin, fill = list(n = 0)) %>%
    ggplot(aes(x = SizeClass, y = FreqBin, fill = n)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(n > 0, n, "")), size = 3.5) +
    scale_fill_gradient(low = "white", high = "#3B528B", name = "Count") +
    labs(title = title, subtitle = "Heatmap",
         x = "Inversion size class", y = "Sharing frequency") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
}

# =============================================================================
# OPTION B — STACKED BAR
# x = size class, fill = frequency bin, y = count
# =============================================================================
plot_stacked_bar <- function(df, sp = NULL) {
  d <- if (!is.null(sp)) filter(df, Species == sp) else df
  title <- if (!is.null(sp)) sp else "All species"
  
  d %>%
    count(SizeClass, FreqBin) %>%
    ggplot(aes(x = SizeClass, y = n, fill = FreqBin)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_viridis_d(
      direction = 1,
      option = "viridis",
      name = "Shared"
    ) +
    labs(title = title,
         x = "Inversion size", y = "Number of inversions") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "right")
}

# =============================================================================
# OPTION C — VIOLIN / BOXPLOT
# x = size class, y = continuous frequency
# =============================================================================
plot_violin <- function(df, sp = NULL) {
  d <- if (!is.null(sp)) filter(df, Species == sp) else df
  title <- if (!is.null(sp)) sp else "All species"
  
  ggplot(d, aes(x = SizeClass, y = Frequency)) +
    geom_violin(fill = "#21908C", alpha = 0.4, color = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "#3B528B") +
    geom_jitter(width = 0.08, size = 1.5, alpha = 0.6, color = "#3B528B") +
    scale_y_continuous(limits = c(0, 1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("0", "25%", "50%", "75%", "100%")) +
    labs(title = title, subtitle = "Violin + boxplot",
         x = "Inversion size class", y = "Sharing frequency") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
}

# =============================================================================
# OPTION D — FACETED BAR (one facet per species)
# =============================================================================
plot_faceted_bar <- function(df) {
  df %>%
    count(Species, SizeClass, FreqBin) %>%
    ggplot(aes(x = SizeClass, y = n, fill = FreqBin)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_viridis_d(
      direction = -1,
      option = "viridis",
      name = "Frequency"
    ) +
    facet_wrap(~ Species, scales = "free_y") +
    labs(title = "All species",
         x = "Inversion size", y = "Number of inversions") +
    theme_classic(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey92", color = "grey60"),
      strip.text       = element_text(face = "bold"),
      legend.position  = "right",
      axis.text.x      = element_text(angle = 30, hjust = 1)
    )
}


plot_stacked_bar_pct <- function(df, sp = NULL) {
  d <- if (!is.null(sp)) filter(df, Species == sp) else df
  title <- if (!is.null(sp)) sp else "All species"
  
  d %>%
    count(SizeClass, FreqBin) %>%
    group_by(SizeClass) %>%
    mutate(pct = n / sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x = SizeClass, y = pct, fill = FreqBin)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_viridis_d(
      direction = -1,
      option = "viridis",
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = title, subtitle = "Stacked bar (percent within size class)",
         x = "Inversion size class", y = "Percentage of inversions") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "right")
}

# =============================================================================
# BUILD PLOT LISTS
# =============================================================================
species_list <- unique(size_freq$Species)

# Dotplots (filtered to shared region)
dotplots <- dot_filtered %>%
  group_by(Species) %>%
  group_map(~ plot_inversion_dotplot(.x), .keep = TRUE) %>%
  setNames(unique(dot_filtered$Species))

heatmaps <- lapply(species_list, function(sp) {
  plot_heatmap(size_freq, sp)
})

plot_heatmap(size_freq)

# Stacked bars
stacked_bars <- lapply(species_list, function(sp) {
  plot_stacked_bar(size_freq, sp)
})

plot_stacked_bar(size_freq)
plot_stacked_bar_pct(size_freq)

# Violins
violins <- lapply(species_list, function(sp) {
  plot_violin(size_freq, sp)
})
plot_violin(size_freq)

# Cross-species faceted bar (single plot)
faceted_bar  <- plot_faceted_bar(size_freq)



