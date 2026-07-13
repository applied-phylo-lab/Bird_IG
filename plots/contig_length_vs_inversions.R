suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
})

# ── INPUTS ────────────────────────────────────────────────────────────────────
INPUT_DIR    <- "/local/storage/kav67/clean_birds/"
SUMMARY_FILE <- file.path(INPUT_DIR, "summary_features.csv")
INV_FILE     <- file.path(INPUT_DIR, "inversion_stats.tsv")

# ── Load summary features — IGH only, one contig per haplotype ────────────────
meta <- read.csv(SUMMARY_FILE, stringsAsFactors = FALSE) %>%
  filter(Locus == "IGH") %>%
  group_by(Order, Species, Haplotype) %>%
  slice_max(NumV, n = 1, with_ties = FALSE) %>%
  ungroup()

# ── Read per-haplotype locus summary and derive full contig length ─────────────
# The per-haplotype refined_ig_loci/summary.csv stores the locus start/end as
# absolute positions (StartPos, EndPos) and as fractions of the full contig
# (RelStart, RelEnd).  ContigLength = EndPos / RelEnd.
# RelStart can be 0 when the locus begins at the contig origin, so only RelEnd
# is used here.

read_contig_length <- function(order, species, haplotype, contig) {
  path <- file.path(INPUT_DIR, order, species, haplotype,
                    "refined_ig_loci", "summary.csv")
  if (!file.exists(path)) return(NA_real_)

  d <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(d)) return(NA_real_)

  row <- d %>% filter(Contig == contig, Locus == "IGH")
  if (nrow(row) == 0 || is.na(row$RelEnd[1]) || row$RelEnd[1] == 0) return(NA_real_)

  row$EndPos[1] / row$RelEnd[1]
}

contig_lengths <- meta %>%
  rowwise() %>%
  mutate(contig_length_bp = read_contig_length(Order, Species, Haplotype, Contig)) %>%
  ungroup() %>%
  filter(!is.na(contig_length_bp))

message(sprintf("Contig lengths computed for %d / %d haplotypes",
                nrow(contig_lengths), nrow(meta)))

# ── Load inversion stats ──────────────────────────────────────────────────────
inv <- read.delim(INV_FILE, stringsAsFactors = FALSE) %>%
  filter(minlen == 250)

# ── Join ──────────────────────────────────────────────────────────────────────
combined <- contig_lengths %>%
  inner_join(inv,
             by = c("Order" = "order",
                    "Species" = "species",
                    "Haplotype" = "haplotype")) %>%
  mutate(contig_length_mb = contig_length_bp / 1e6)

message(sprintf("%d haplotypes after joining with inversion stats", nrow(combined)))


# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot(combined, aes(x = contig_length_mb, y = num_inversions)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.7) +
  scale_x_log10(labels = function(x) paste0(x, " Mb")) +
  scale_y_log10() +
  labs(
    x     = "Contig length (Mb, log scale)",
    y     = "Number of inversions (log scale, min length 250 bp)",
    title = "Contig length vs. number of inversions (IGH)"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggplot(combined, aes(x = contig_length_mb, y = genes_on_inv)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.7) +
  scale_x_log10(labels = function(x) paste0(x, " Mb")) +
  scale_y_log10() +
  labs(
    x     = "Contig length (Mb, log scale)",
    y     = "Genes in inversion-associated regions (log scale)",
    title = "Contig length vs. genes in inversion-associated regions (IGH)"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))


# ── Quick summary ─────────────────────────────────────────────────────────────
cor_test <- cor.test(combined$contig_length_mb, combined$num_inversions,
                     method = "spearman")
message(sprintf("Spearman r = %.3f,  p = %.3e",
                cor_test$estimate, cor_test$p.value))
