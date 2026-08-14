library(ggplot2)
library(dplyr)
library(readr)

# === INPUT ===
inversions <- read_tsv("/local/storage/kav67/clean_birds/Songbirds/House_Finch/bHaeMex1_pri/self_align_contigs/inversion_window_summary.tsv", show_col_types = FALSE)
inversions <- read_tsv("/local/storage/kav67/clean_birds/Songbirds/House_Finch/bHaeMex1.pri/self_align_contigs/window_summary_deduplicated.tsv", show_col_types = FALSE)
inversions <- read_tsv("/local/storage/kav67/zebrafinch/window_summary_deduplicated.tsv", show_col_types = FALSE)

#california_scrub_jay
# optional: order contigs by size or name
inversions <- inversions %>%
  group_by(chr) %>%
  mutate(chr_len = max(window_end)) %>%
  ungroup() %>%
  arrange(chr, window_start)

# === BUILD CONTIG OFFSETS (for Manhattan-style continuous x-axis) ===
offsets <- inversions %>%
  group_by(chr) %>%
  summarise(chr_len = max(window_end)) %>%
  mutate(chr_start = lag(cumsum(chr_len), default = 0)) %>%
  select(chr, chr_start)

inversions <- inversions %>%
  left_join(offsets, by = "chr") %>%
  mutate(global_center = center + chr_start)

# === PLOT SETTINGS ===
highlight_chr<-"JAPYKC010000061.1"
highlight_chr<-"CM109827.1"

inversions <- inversions %>%
  mutate(color = ifelse(chr == highlight_chr, "IGH", "other"))

# === MANHATTAN PLOT ===
ggplot(inversions, aes(x = global_center, y = num_inversions, color = color)) +
  geom_point(aes(size = color)) +
  scale_size_manual(values = c("IGH" =1, "other" = 1)) +
  scale_color_manual(values = c("IGH" = "red", "other" = "black")) +
  labs(
    x = "Genomic position",
    y = "Inversion count per window"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )


ggplot(inversions, aes(x = global_center, y = covered_fraction, color = color)) +
  geom_point(aes(size = color)) +
  scale_size_manual(values = c("IGH" =1, "other" = 1)) +
  scale_color_manual(values = c("IGH" = "red", "other" = "black")) +
  labs(
    x = "Genomic position",
    y = "Inversion covered fraction"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )


# =====================================================================
# VERSION WITH ALTERNATING CHROMOSOME COLORS AND CHROMOSOME LABELS
# =====================================================================
# The contig names in the window summaries are GenBank accessions, so the
# accession -> chromosome mapping is taken from the NCBI assembly report
# (files under assembly_reports/, downloaded from the NCBI FTP site):
#   zebra finch bTaeGut7.mat  GCA_048771995.1  (CM1097xx.1 accessions)
#   house finch bHaeMex1.pri  GCF_027477595.1  (CM0502xx.1 / JAPYKC0100000xx.1)

report_dir <- "/home/kav67/Bird_IG/manhattan_plot_inversion_coverage/assembly_reports"

# --- read an NCBI assembly report and return accession -> chromosome ---
read_assembly_report <- function(path) {
  read_tsv(
    path,
    comment = "#",
    col_names = c("seq_name", "seq_role", "assigned_molecule",
                  "molecule_type", "genbank_accn", "relationship",
                  "refseq_accn", "assembly_unit", "seq_length", "ucsc_name"),
    col_types = cols(.default = col_character(), seq_length = col_double())
  ) %>%
    mutate(across(where(is.character), ~ gsub("\r", "", .x))) %>%
    # unlocalized scaffolds are drawn with the chromosome they belong to;
    # unplaced scaffolds have no chromosome and go into a trailing block
    mutate(chrom = ifelse(seq_role == "unplaced-scaffold",
                          "Unplaced", assigned_molecule)) %>%
    select(genbank_accn, refseq_accn, chrom, seq_length)
}

# --- sort key: chromosomes run 1, 1A, 2, 3, 4, 4A, ..., Z, W, MT, Unplaced ---
chrom_sort_key <- function(chrom) {
  lead_num <- suppressWarnings(as.numeric(sub("^([0-9]+).*$", "\\1", chrom)))
  special  <- match(chrom, c("Z", "W", "MT", "Unplaced"))
  tibble(
    key    = ifelse(!is.na(lead_num), lead_num, 1000 + special),
    suffix = sub("^[0-9]+", "", chrom)
  )
}

# --- attach chromosome names and Manhattan x-coordinates ---
add_chrom_coords <- function(df, report, unplaced_gap_frac = 0.012) {
  # accessions in the data may be GenBank or RefSeq
  map <- bind_rows(
    report %>% transmute(accn = genbank_accn, chrom, seq_length),
    report %>% transmute(accn = refseq_accn,  chrom, seq_length)
  ) %>% filter(!is.na(accn), accn != "na")

  df <- df %>% left_join(map, by = c("chr" = "accn"))
  missing <- unique(df$chr[is.na(df$chrom)])
  if (length(missing) > 0) {
    warning("no chromosome assignment for: ", paste(missing, collapse = ", "))
    df <- df %>% mutate(chrom = ifelse(is.na(chrom), "Unplaced", chrom))
  }

  # lay the sequences end to end: chromosomes in order, contigs within a
  # chromosome (the assembled molecule, then its unlocalized scaffolds)
  # ordered largest first
  seq_layout <- df %>%
    group_by(chr, chrom) %>%
    summarise(seq_length = max(seq_length[1], max(window_end), na.rm = TRUE),
              .groups = "drop") %>%
    bind_cols(chrom_sort_key(.$chrom)) %>%
    arrange(key, suffix, desc(seq_length)) %>%
    mutate(seq_start = lag(cumsum(seq_length), default = 0))

  # the unplaced scaffolds only add up to a sliver, too narrow for its axis
  # label; push the block right so the label clears the last chromosome
  gap <- unplaced_gap_frac * max(seq_layout$seq_start + seq_layout$seq_length)
  seq_layout <- seq_layout %>%
    mutate(seq_start = seq_start + ifelse(chrom == "Unplaced", gap, 0))

  chrom_layout <- seq_layout %>%
    group_by(chrom) %>%
    summarise(chrom_start = min(seq_start),
              chrom_end   = max(seq_start + seq_length),
              key = key[1], suffix = suffix[1], .groups = "drop") %>%
    arrange(key, suffix) %>%
    mutate(chrom_center = (chrom_start + chrom_end) / 2,
           band = row_number() %% 2)

  list(
    data = df %>%
      left_join(seq_layout %>% select(chr, seq_start), by = "chr") %>%
      left_join(chrom_layout %>% select(chrom, band), by = "chrom") %>%
      mutate(global_center = center + seq_start),
    chrom_layout = chrom_layout
  )
}

# --- x-axis ticks: label wide chromosomes individually, and collapse runs of
#     consecutive microchromosomes into a single range label ("23-39") ---
build_x_ticks <- function(cl, min_label_frac) {
  genome_len <- max(cl$chrom_end)
  cl <- cl %>%
    # only numbered chromosomes get merged; Z, W, MT and Unplaced keep their own
    mutate(small = (chrom_end - chrom_start) / genome_len < min_label_frac &
                     key < 1000)

  runs <- rle(cl$small)
  cl$run <- rep(seq_along(runs$lengths), runs$lengths)

  bind_rows(
    cl %>%
      filter(!small) %>%
      transmute(pos = chrom_center, label = chrom),
    cl %>%
      filter(small) %>%
      group_by(run) %>%
      summarise(pos = (min(chrom_start) + max(chrom_end)) / 2,
                label = if (n() == 1) first(chrom)
                        else paste0(first(chrom), "-", last(chrom)),
                .groups = "drop") %>%
      select(pos, label)
  ) %>%
    arrange(pos)
}

# --- the plot itself ---
# min_label_frac: chromosomes narrower than this fraction of the genome are
# merged into a range label. Raise it to merge more, set it to 0 to label
# every chromosome separately.
manhattan_alt <- function(df, report, yvar, ylab, min_label_frac = 0.019) {
  lay <- add_chrom_coords(df, report)
  d   <- lay$data
  ticks <- build_x_ticks(lay$chrom_layout, min_label_frac)

  # "Unplaced" is a wide label sitting over a very narrow block. Tick labels
  # are centred on the tick, so pad it on the left: the visible text is pushed
  # into the right margin instead of colliding with the last chromosome.
  ticks <- ticks %>%
    mutate(label = ifelse(label == "Unplaced",
                          paste0(strrep(" ", 2 * nchar(label)), label), label))

  ggplot(d, aes(x = global_center, y = .data[[yvar]],
                color = factor(band))) +
    geom_point(size = 0.9) +
    scale_color_manual(values = c("0" = "grey25", "1" = "grey65")) +
    scale_x_continuous(breaks = ticks$pos, labels = ticks$label,
                       expand = expansion(mult = 0.01)) +
    labs(x = "Chromosome", y = ylab) +
    theme_classic() +
    theme(
      legend.position = "none",
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.text = element_text(size = 10),
      axis.ticks.x = element_blank(),
      axis.title = element_text(size = 14),
      # room for the "Unplaced" label, which overhangs its narrow block
      plot.margin = margin(5.5, 45, 5.5, 5.5)
    )
}

# === HOUSE FINCH (bHaeMex1.pri, GCF_027477595.1) ===
hf <- read_tsv("/local/storage/kav67/clean_birds/Songbirds/House_Finch/bHaeMex1.pri/self_align_contigs/window_summary_deduplicated.tsv",
               show_col_types = FALSE)
hf_report <- read_assembly_report(
  file.path(report_dir, "GCF_027477595.1_bHaeMex1.pri_assembly_report.txt"))

p_counts <- manhattan_alt(hf, hf_report, "num_inversions",
                          "Inversion count per window")
p_counts

p_frac <- manhattan_alt(hf, hf_report, "covered_fraction",
                        "Inversion covered fraction")
p_frac

# --- where does IGL sit on this axis? -------------------------------------
# IGL is at CM050242.1:5522938-5545644, i.e. chromosome 19. Its x position is
# the locus midpoint plus that sequence's offset in the end-to-end layout.
igl_x <- local({
  lay <- add_chrom_coords(hf, hf_report)
  off <- lay$data %>% filter(chr == "CM050242.1") %>% slice(1) %>% pull(seq_start)
  (5522938 + 5545644) / 2 + off
})   # 921,528,336 of an axis spanning 0-1,155,016,870

p_counts +
  annotate("segment", x = igl_x, xend = igl_x, y = -Inf, yend = Inf,
           colour = "red", linewidth = 0.4, linetype = "dashed") +
  annotate("text", x = igl_x, y = Inf, label = "IGL", colour = "red",
           vjust = 1.4, hjust = -0.2, size = 4)

# ggsave("/home/kav67/Bird_IG/figures/manhattan_inversions_bHaeMex1.svg",
#        p_counts, width = 11, height = 4)
# ggsave("/home/kav67/Bird_IG/figures/manhattan_covered_fraction_bHaeMex1.svg",
#        p_frac, width = 11, height = 4)

# === ZEBRA FINCH (bTaeGut7.mat, the data loaded at the top of the script) ===
# zf_report <- read_assembly_report(
#   file.path(report_dir, "GCA_048771995.1_bTaeGut7.mat_assembly_report.txt"))
# manhattan_alt(inversions, zf_report, "num_inversions",
#               "Inversion count per window")

