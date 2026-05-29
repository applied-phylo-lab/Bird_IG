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
    productive = tolower(trimws(as.character(Productive))) %in% c("true", "1", "yes")
  )

vgp_table <- read.delim(file.path(input_dir, "IGH_VGP_table.tsv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
df <- df %>%
  left_join(vgp_table %>% select(Species, LatinName) %>% distinct(),
            by = "Species")

# ── Determine D-gene orientation per IGH haplotype ────────────────────────────
# upstream   → D genes at lower genomic pos → rel_pos near 0 → flip needed
# downstream → D genes at higher genomic pos → rel_pos near 1 → keep as-is
# ORIENTATION_THRESHOLD controls how strict "mostly on one side" is

get_d_orientation <- function(order, species, haplotype) {
  path <- file.path(input_dir, order, species, haplotype, "IGHD.csv")
  if (!file.exists(path)) return(NA_character_)

  d <- tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(d)) return(NA_character_)
  if (!"Location Relative to V-Cluster" %in% names(d)) return(NA_character_)

  locs <- na.omit(trimws(d[["Location Relative to V-Cluster"]]))
  # v_cluster genes are inside the V region → ambiguous for orientation
  locs <- locs[locs %in% c("upstream", "downstream")]
  if (length(locs) == 0) return(NA_character_)

  n_down <- sum(locs == "downstream")
  n_up   <- sum(locs == "upstream")
  frac   <- max(n_down, n_up) / length(locs)

  if (frac < ORIENTATION_THRESHOLD) return(NA_character_)  # too mixed → exclude
  if (n_down >= n_up) "downstream" else "upstream"
}

igh_haplotypes <- df %>%
  filter(Locus == "IGH") %>%
  distinct(Order, Species, Haplotype)

orientation_df <- igh_haplotypes %>%
  rowwise() %>%
  mutate(d_orientation = get_d_orientation(Order, Species, Haplotype)) %>%
  ungroup()

n_down    <- sum(orientation_df$d_orientation == "downstream", na.rm = TRUE)
n_up      <- sum(orientation_df$d_orientation == "upstream",   na.rm = TRUE)
n_excluded <- sum(is.na(orientation_df$d_orientation))
message(sprintf(
  "D orientation: %d downstream (keep), %d upstream (flip), %d excluded (mixed/missing)",
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
    rel_pos      = rank(Pos, ties.method = "first") / n(),
    oriented_pos = if_else(d_orientation == "upstream", 1 - rel_pos, rel_pos),
    n_prod_rss   = sum(has_rss & productive)
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
    rel_pos      = rank(Pos, ties.method = "first") / n(),
    oriented_pos = if_else(igl_orientation == "minus", 1 - rel_pos, rel_pos),
    n_prod_rss   = sum(has_rss & productive)
  ) %>%
  ungroup()

# ── Plot helper ───────────────────────────────────────────────────────────────
make_oriented_plot <- function(data, locus, x_label = TRUE) {
  pos_data <- data %>%
    filter(has_rss) %>%
    mutate(group = if_else(n_prod_rss > 1,
                           "Multiple productive RSS",
                           "Single productive RSS"))

  x_lab <- if (x_label) {
    "Oriented relative position  (0% = away from J genes, 100% = toward J genes)"
  } else ""

  ggplot(pos_data, aes(x = oriented_pos, fill = group, color = group)) +
    geom_histogram(aes(y = after_stat(density)),
                   binwidth = 0.05, alpha = 0.5, position = "identity") +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = function(x) paste0(round(x * 100), "%")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values  = GROUP_COLORS[[locus]], name = NULL) +
    scale_color_manual(values = GROUP_COLORS[[locus]], name = NULL) +
    labs(
      title = if (locus == "IGH") {
        sprintf("IGH  (n = %d haplotypes, %d flipped via D-gene location)",
                n_distinct(data$Haplotype), n_up)
      } else {
        sprintf("IGL  (n = %d haplotypes, %d flipped via majority strand)",
                n_distinct(data$Haplotype), n_igl_minus)
      },
      x = x_lab,
      y = "Density"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold"),
      axis.title      = element_text(size = 12),
      axis.text       = element_text(size = 10),
      legend.position = "top"
    )
}

p_igh_oriented <- make_oriented_plot(df_igh, "IGH", x_label = TRUE)
p_igl_oriented <- make_oriented_plot(df_igl, "IGL", x_label = TRUE)

print(p_igh_oriented)
print(p_igl_oriented)
print(p_igh_oriented / p_igl_oriented)

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
