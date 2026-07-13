suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(patchwork)
})

# ── INPUTS ────────────────────────────────────────────────────────────────────
input_dir      <- "/local/storage/kav67/clean_birds/"
gene_list_path <- file.path(input_dir, "gene_list.csv")

# Minimum fraction that must agree on one side to call orientation.
# Applied to both IGH (D-gene location) and IGL (majority strand).
# Haplotypes below this threshold are excluded as ambiguous.
ORIENTATION_THRESHOLD <- 0.70

LOCUS_COLORS <- c(IGH = "#87b4dc", IGL = "#638E6E")
GROUP_COLORS <- list(
  IGH = c("Multiple productive RSS" = "lightgrey", "Single productive RSS" = "#87b4dc"),
  IGL = c("Multiple productive RSS" = "lightgrey", "Single productive RSS" = "#638E6E")
)

# ── Load V gene data (same source as rss_correlation.R) ───────────────────────
df <- read.csv(gene_list_path, stringsAsFactors = FALSE, check.names = FALSE) %>%
  separate(Source, into = c("Order", "Species", "Haplotype"),
           sep = "/", remove = FALSE) %>%
  mutate(
    has_rss    = !is.na(Heptamer) & trimws(Heptamer) != "",
    productive = tolower(trimws(as.character(Productive))) %in% c("true", "1", "yes","TRUE")
  )

vgp_table <- read.delim(file.path(input_dir, "IGH_VGP_table.tsv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
df <- df %>%
  left_join(vgp_table %>% select(Species, LatinName) %>% distinct(),
            by = "Species")

# ── Determine D-gene orientation per IGH haplotype ────────────────────────────
# Two methods are computed and kept for comparison:
#
#  d_orientation_meta — uses the pre-computed "Location Relative to V-Cluster"
#    column from IGHD.csv (values: "upstream", "downstream", "v_cluster").
#
#  d_orientation_pos  — analogous to IGL: directly compares D gene Pos values
#    against the median V gene Pos for the same haplotype (from gene_list.csv).
#    D genes mostly at higher Pos → downstream; mostly at lower Pos → upstream.
#
# d_orientation (used for oriented_pos) is set to d_orientation_pos.
# ORIENTATION_THRESHOLD controls how strict "mostly on one side" is for both.

# Method 1: pre-computed metadata field
get_d_orientation_meta <- function(order, species, haplotype) {
  path <- file.path(input_dir, order, species, haplotype, "IGHD.csv")
  if (!file.exists(path)) return(NA_character_)

  d <- tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(d)) return(NA_character_)
  if (!"Location Relative to V-Cluster" %in% names(d)) return(NA_character_)

  locs <- na.omit(trimws(d[["Location Relative to V-Cluster"]]))
  locs <- locs[locs %in% c("upstream", "downstream")]  # exclude v_cluster
  if (length(locs) == 0) return(NA_character_)

  n_down <- sum(locs == "downstream")
  n_up   <- sum(locs == "upstream")
  frac   <- max(n_down, n_up) / length(locs)

  if (frac < ORIENTATION_THRESHOLD) return(NA_character_)
  if (n_down >= n_up) "downstream" else "upstream"
}

# Method 2: position-based (analogous to IGL majority-strand approach)
get_d_orientation_pos <- function(order, species, haplotype, v_median_pos) {
  path <- file.path(input_dir, order, species, haplotype, "IGHD.csv")
  if (!file.exists(path)) return(NA_character_)

  d <- tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(d) || !"Pos" %in% names(d)) return(NA_character_)

  d_pos <- na.omit(d[["Pos"]])
  if (length(d_pos) == 0) return(NA_character_)

  n_down <- sum(d_pos > v_median_pos)   # D genes at higher Pos than V cluster
  n_up   <- sum(d_pos < v_median_pos)   # D genes at lower  Pos than V cluster
  total  <- n_down + n_up
  if (total == 0) return(NA_character_)

  frac <- max(n_down, n_up) / total
  if (frac < ORIENTATION_THRESHOLD) return(NA_character_)
  if (n_down >= n_up) "downstream" else "upstream"
}

# Precompute median V gene position per IGH haplotype
igh_v_medians <- df %>%
  filter(Locus == "IGH") %>%
  group_by(Order, Species, Haplotype) %>%
  summarise(v_median_pos = median(Pos, na.rm = TRUE), .groups = "drop")

igh_haplotypes <- df %>%
  filter(Locus == "IGH") %>%
  distinct(Order, Species, Haplotype) %>%
  left_join(igh_v_medians, by = c("Order", "Species", "Haplotype"))

orientation_df <- igh_haplotypes %>%
  rowwise() %>%
  mutate(
    d_orientation_meta = get_d_orientation_meta(Order, Species, Haplotype),
    d_orientation_pos  = get_d_orientation_pos(Order, Species, Haplotype, v_median_pos),
    d_orientation      = d_orientation_pos   # primary method used for oriented_pos
  ) %>%
  ungroup()

# ── Agreement between the two methods ─────────────────────────────────────────
both_called <- orientation_df %>%
  filter(!is.na(d_orientation_meta), !is.na(d_orientation_pos))
n_agree    <- sum(both_called$d_orientation_meta == both_called$d_orientation_pos)
n_disagree <- sum(both_called$d_orientation_meta != both_called$d_orientation_pos)
message(sprintf(
  "Orientation methods agree: %d / %d haplotypes  (%d disagree)",
  n_agree, nrow(both_called), n_disagree
))
if (n_disagree > 0) {
  message("Disagreeing haplotypes:")
  both_called %>%
    filter(d_orientation_meta != d_orientation_pos) %>%
    select(Order, Species, Haplotype, d_orientation_meta, d_orientation_pos) %>%
    as.data.frame() %>%
    print()
}

n_down     <- sum(orientation_df$d_orientation == "downstream", na.rm = TRUE)
n_up       <- sum(orientation_df$d_orientation == "upstream",   na.rm = TRUE)
n_excluded <- sum(is.na(orientation_df$d_orientation))
message(sprintf(
  "D orientation (pos): %d downstream (keep), %d upstream (flip), %d excluded (mixed/missing)",
  n_down, n_up, n_excluded
))

orientation_df <- orientation_df %>% filter(!is.na(d_orientation))

# Also collect the majority D-gene strand per haplotype (needed for the
# same-strand vs opposite-strand plot below).
# D genes can sit on both strands, so we apply the same threshold.
get_d_majority_strand <- function(order, species, haplotype) {
  path <- file.path(input_dir, order, species, haplotype, "IGHD.csv")
  if (!file.exists(path)) return(NA_character_)
  d <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
                error = function(e) NULL)
  if (is.null(d) || !"Strand" %in% names(d)) return(NA_character_)
  strands <- na.omit(d[["Strand"]])
  if (length(strands) == 0) return(NA_character_)
  n_plus  <- sum(strands == "+")
  n_minus <- sum(strands == "-")
  frac    <- max(n_plus, n_minus) / length(strands)
  if (frac < ORIENTATION_THRESHOLD) return(NA_character_)
  if (n_plus >= n_minus) "+" else "-"
}

orientation_df <- orientation_df %>%
  rowwise() %>%
  mutate(d_maj_strand = get_d_majority_strand(Order, Species, Haplotype)) %>%
  ungroup()

# ── Oriented IGH positions ─────────────────────────────────────────────────────
df_igh <- df %>%
  filter(Locus == "IGH") %>%
  inner_join(orientation_df, by = c("Order", "Species", "Haplotype")) %>%
  group_by(LatinName, Haplotype) %>%
  mutate(
    rel_pos      = (rank(Pos, ties.method = "first") - 1) / (n() - 1),
    oriented_pos = if_else(d_orientation == "upstream", 1 - rel_pos, rel_pos),
    n_prod_rss   = sum(has_rss & productive),
    # Infer reference strand from D-gene orientation:
    # downstream D genes (higher Pos) → + strand locus; upstream → - strand
    ref_strand   = if_else(d_orientation == "downstream", "+", "-"),
    direction    = if_else(Strand == ref_strand,
                           "same direction", "opposite direction")
  ) %>%
  ungroup()

# ── Determine strand-based orientation per IGL haplotype ─────────────────────
# V genes point their 3' end toward the J-C region.
# Majority + strand → locus runs left-to-right → J genes at higher Pos → keep
# Majority - strand → locus runs right-to-left → J genes at lower Pos → flip
# Haplotypes below ORIENTATION_THRESHOLD are excluded as ambiguous.

igl_orientation_df <- df %>%
  filter(Locus == "IGL") %>%
  group_by(Order, Species, Haplotype) %>%
  summarise(
    n_total  = n(),
    n_plus   = sum(Strand == "+"),
    n_minus  = sum(Strand == "-"),
    frac_maj = max(n_plus, n_minus) / n_total,
    maj_strand = if_else(n_plus >= n_minus, "+", "-"),
    .groups = "drop"
  ) %>%
  mutate(igl_orientation = if_else(
    frac_maj >= ORIENTATION_THRESHOLD,
    if_else(maj_strand == "+", "plus", "minus"),
    NA_character_
  ))

n_igl_plus     <- sum(igl_orientation_df$igl_orientation == "plus",  na.rm = TRUE)
n_igl_minus    <- sum(igl_orientation_df$igl_orientation == "minus", na.rm = TRUE)
n_igl_excluded <- sum(is.na(igl_orientation_df$igl_orientation))
message(sprintf(
  "IGL strand orientation: %d plus-strand (keep), %d minus-strand (flip), %d excluded (ambiguous)",
  n_igl_plus, n_igl_minus, n_igl_excluded
))

igl_orientation_df <- igl_orientation_df %>% filter(!is.na(igl_orientation))

# ── Oriented IGL positions ────────────────────────────────────────────────────
df_igl <- df %>%
  filter(Locus == "IGL") %>%
  inner_join(igl_orientation_df %>% select(Order, Species, Haplotype, igl_orientation),
             by = c("Order", "Species", "Haplotype")) %>%
  group_by(LatinName, Haplotype) %>%
  mutate(
    rel_pos      = (rank(Pos, ties.method = "first") - 1) / (n() - 1),
    oriented_pos = if_else(igl_orientation == "minus", 1 - rel_pos, rel_pos),
    n_prod_rss   = sum(has_rss & productive),
    # Reference strand is the majority strand used to orient the locus
    ref_strand   = if_else(igl_orientation == "plus", "+", "-"),
    direction    = if_else(Strand == ref_strand,
                           "same direction", "opposite direction")
  ) %>%
  ungroup()

# ── Plot helper ───────────────────────────────────────────────────────────────
make_oriented_plot <- function(data, locus, x_label = TRUE, stat = "density", y_limits = NULL) {
  pos_data <- data %>%
    filter(has_rss, productive) %>%
    mutate(group = if_else(n_prod_rss > 1,
                           "Multiple productive RSS",
                           "Single productive RSS"))

  x_lab <- if (x_label) {
    "Oriented relative position  (0% = away from J genes, 100% = toward J genes)"
  } else ""

  y_lab <- if (stat == "count") "Count" else "Density"

  n_multi  <- pos_data %>%
    distinct(LatinName, Haplotype, n_prod_rss) %>%
    filter(n_prod_rss > 1) %>%
    nrow()
  n_single <- pos_data %>%
    distinct(LatinName, Haplotype, n_prod_rss) %>%
    filter(n_prod_rss == 1) %>%
    nrow()

  hist_layer <- if (stat == "count") {
    geom_histogram(aes(y = after_stat(count)),
                   binwidth = 0.05, alpha = 0.5, position = "identity")
  } else {
    geom_histogram(aes(y = after_stat(density)),
                   binwidth = 0.05, alpha = 0.5, position = "identity")
  }

  ggplot(pos_data, aes(x = oriented_pos, fill = group, color = group)) +
    hist_layer +
    scale_x_continuous(
      breaks = seq(0, 1, 0.2),
      labels = function(x) paste0(round(x * 100), "%")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = y_limits) +
    scale_fill_manual(values  = GROUP_COLORS[[locus]], name = NULL) +
    scale_color_manual(values = GROUP_COLORS[[locus]], name = NULL) +
    annotate(
      "text",
      x = 0, y = Inf,
      label = sprintf("Single RSS: %d haplotypes\nMultiple RSS: %d haplotypes",
                      n_single, n_multi),
      hjust = 0, vjust = 1.4,
      size  = 3.2, color = "grey30"
    ) +
    labs(
      title = locus,
      x     = x_lab,
      y     = y_lab
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold"),
      axis.title      = element_text(size = 14),
      axis.text       = element_text(size = 10),
      legend.position = "none"
    )
}

p_igh_oriented <- make_oriented_plot(df_igh, "IGH", x_label = FALSE)
p_igl_oriented <- make_oriented_plot(df_igl, "IGL", x_label = TRUE)

(p_combined_simple+(p_igh_oriented / p_igl_oriented))
# ── Strand-split versions: + strand up, - strand down ────────────────────────
# Density for genes on the + strand is plotted above the axis;
# density for genes on the - strand is plotted below.
# Y-axis labels are absolute values (no minus sign).

make_oriented_plot_strand <- function(data, locus, x_label = TRUE, stat = "density", y_limits = NULL) {
  pos_data <- data %>%
    filter(has_rss, productive) %>%
    mutate(group = if_else(n_prod_rss > 1,
                           "Multiple productive RSS",
                           "Single productive RSS"))

  x_lab <- if (x_label) {
    "Oriented relative position  (0% = away from J genes, 100% = toward J genes)"
  } else ""

  y_lab <- if (stat == "count") "Count" else "Density"

  n_multi  <- pos_data %>%
    distinct(LatinName, Haplotype, n_prod_rss) %>%
    filter(n_prod_rss > 1) %>%
    nrow()
  n_single <- pos_data %>%
    distinct(LatinName, Haplotype, n_prod_rss) %>%
    filter(n_prod_rss == 1) %>%
    nrow()

  # Precompute binned values manually so each group is normalised by its own
  # total — matching how geom_histogram(after_stat(density)) works per fill
  # group in the plain plot.  Both direction halves share the same denominator,
  # so the split bars add up to the same height as the unsplit histogram.
  bin_w   <- 0.05
  breaks  <- seq(0, 1, by = bin_w)
  mids    <- breaks[-length(breaks)] + bin_w / 2

  group_totals <- pos_data %>%
    count(group, name = "group_total")

  dir_bins <- pos_data %>%
    mutate(bin_mid = mids[findInterval(oriented_pos, breaks, rightmost.closed = TRUE)]) %>%
    count(group, direction, bin_mid) %>%
    left_join(group_totals, by = "group") %>%
    mutate(
      density_val = n / (group_total * bin_w),
      y_val = if_else(direction == "same direction", density_val, -density_val),
      y_count = if_else(direction == "same direction", n, -n)
    )

  plot_y <- if (stat == "count") "y_count" else "y_val"

  ggplot(dir_bins, aes(x = bin_mid, y = .data[[plot_y]],
                       fill = group, color = group)) +
    geom_col(position = "identity", alpha = 0.5, width = bin_w) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey60") +
    #annotate("text", x = 1, y =  Inf, label = "same direction",
    #         hjust = 1, vjust = 1.5, size = 3.2, color = "grey40") +
    #annotate("text", x = 1, y = -Inf, label = "opposite direction",
    #         hjust = 1, vjust = -0.5, size = 3.2, color = "grey40") +
    annotate(
      "text",
      x = 0, y = Inf,
      label = sprintf("Single RSS: %d haplotypes\nMultiple RSS: %d haplotypes",
                      n_single, n_multi),
      hjust = 0, vjust = 1.4, size = 3.2, color = "grey30"
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = function(x) paste0(round(x * 100), "%")
    ) +
    scale_y_continuous(
      labels = function(y) abs(y),
      expand = expansion(mult = c(0.05, 0.05)),
      limits = y_limits
    ) +
    scale_fill_manual(values  = GROUP_COLORS[[locus]], name = NULL) +
    scale_color_manual(values = GROUP_COLORS[[locus]], name = NULL) +
    labs(
      title = locus,
      x     = x_lab,
      y     = y_lab
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold"),
      axis.title      = element_text(size = 14),
      axis.text       = element_text(size = 10),
      legend.position = "none"
    )
}

p_igh_oriented_strand <- make_oriented_plot_strand(df_igh, "IGH", x_label = FALSE)
p_igl_oriented_strand <- make_oriented_plot_strand(df_igl, "IGL", x_label = TRUE)

print(p_igh_oriented_strand / p_igl_oriented_strand)

# ── Count versions ────────────────────────────────────────────────────────────
p_igh_oriented_count        <- make_oriented_plot(df_igh, "IGH", x_label = FALSE, stat = "count")
p_igl_oriented_count        <- make_oriented_plot(df_igl, "IGL", x_label = TRUE,  stat = "count")
p_igh_oriented_strand_count <- make_oriented_plot_strand(df_igh, "IGH", x_label = FALSE, stat = "count")
p_igl_oriented_strand_count <- make_oriented_plot_strand(df_igl, "IGL", x_label = TRUE,  stat = "count")

print(p_igh_oriented_count / p_igl_oriented_count)
print(p_igh_oriented_strand_count / p_igl_oriented_strand_count)


# Density, shared scale
ylim_dens   <- c(0,14)
p_igh_oriented_shared <- make_oriented_plot(df_igh, "IGH", x_label = FALSE, y_limits = ylim_dens)
p_igl_oriented_shared <- make_oriented_plot(df_igl, "IGL", x_label = TRUE,  y_limits = ylim_dens)
print(p_igh_oriented_shared / p_igl_oriented_shared)
(p_combined_simple+(p_igh_oriented_shared / p_igl_oriented_shared))

# Strand-split density, shared symmetric scale
ylim_strand_dens  <- c(-20,20)
p_igh_oriented_strand_shared <- make_oriented_plot_strand(df_igh, "IGH", x_label = FALSE, y_limits = ylim_strand_dens)
p_igl_oriented_strand_shared <- make_oriented_plot_strand(df_igl, "IGL", x_label = TRUE,  y_limits = ylim_strand_dens)
print(p_igh_oriented_strand_shared / p_igl_oriented_strand_shared)


# ── Sanity check: how many haplotypes were flipped vs kept ────────────────────
# ── IGH: single productive RSS — same vs opposite strand as majority D genes ──
# For haplotypes with exactly one productive RSS gene, check whether that gene's
# strand matches the majority D-gene strand. Same strand = functional orientation
# (pointing toward the D-J-C region); opposite strand = inverted relative to D.

strand_comparison_df <- df_igh %>%
  # Only haplotypes where majority D-gene strand is clear
  dplyr::filter(!is.na(d_maj_strand)) %>%
  # Only haplotypes with exactly one productive RSS gene
  filter(n_prod_rss == 1) %>%
  # Only the RSS-bearing genes themselves
  filter(has_rss) %>%
  mutate(
    strand_match = if_else(
      Strand == d_maj_strand,
      "Same strand as D genes",
      "Opposite strand to D genes"
    )
  )

n_same <- sum(strand_comparison_df$strand_match == "Same strand as D genes")
n_opp  <- sum(strand_comparison_df$strand_match == "Opposite strand to D genes")
message(sprintf(
  "Single productive RSS strand: %d same as D, %d opposite to D", n_same, n_opp
))

strand_colors <- c(
  "Same strand as D genes"     = "#87b4dc",
  "Opposite strand to D genes" = "#c0392b"
)

p_igh_strand <- ggplot(strand_comparison_df,
                       aes(x = oriented_pos,
                           fill = strand_match, color = strand_match)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.05, alpha = 0.5, position = "identity") +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = function(x) paste0(round(x * 100), "%")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values  = strand_colors, name = NULL) +
  scale_color_manual(values = strand_colors, name = NULL) +
  labs(
    title = sprintf(
      "IGH single productive RSS: strand relative to D genes  (same = %d, opposite = %d)",
      n_same, n_opp
    ),
    x = "Oriented relative position  (0% = away from J genes, 100% = toward J genes)",
    y = "Density"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 10),
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 10),
    legend.position = "top"
  )

print(p_igh_strand)
