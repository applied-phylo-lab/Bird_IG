#!/usr/bin/env Rscript

library(ggplot2)
library(data.table)


# ---- overlap function ----
calc_overlap <- function(a_start, a_end, b_start, b_end) {
  ov <- min(a_end, b_end) - max(a_start, b_start)
  return(max(0, ov))
}

# ---- main script ----

infile <- "/local/storage/kav67/Bird_data/Hummingbirds/LongTailed_Hermit/bPhaSup1_alt/IGH_self.tsv"


# Read LASTZ TSV (header may start with '#')
lines <- readLines(infile)
header <- gsub("^#", "", lines[1])
colnames <- strsplit(header, "\\s+")[[1]]
colnames <- gsub("\\+","",colnames)
df <- fread(infile, skip = 1, header = FALSE)
setnames(df, colnames)

# Filter inversions
inv <- df[df$strand1 != df$strand2]

if (nrow(inv) < 2) {
  stop("Less than 2 inversions found — no overlaps can be computed.")
}

# ---- compute all overlaps ----
overlaps <- c()

for (i in 1:(nrow(inv)-1)) {
  for (j in (i+1):nrow(inv)) {
    
    # Overlap in first coordinate system
    ov1 <- calc_overlap(inv$start1[i], inv$end1[i],
                        inv$start1[j], inv$end1[j])
    
    # Overlap in second coordinate system
    ov2 <- calc_overlap(inv$start2[i], inv$end2[i],
                        inv$start2[j], inv$end2[j])
    
    # Save nonzero overlaps
    if (ov1 > 0) overlaps <- c(overlaps, ov1)
    if (ov2 > 0) overlaps <- c(overlaps, ov2)
  }
}

if (length(overlaps) == 0) {
  stop("No inversion overlaps found.")
}

overlap_df <- data.frame(overlap_len = overlaps)

# ---- histogram ----
p <- ggplot(overlap_df, aes(x = overlap_len)) +
  geom_histogram(bins = 100) +
  theme_bw() +
  labs(
    title = "Histogram of inversion overlap lengths",
    x = "Overlap length (bp)",
    y = "Count"
  )

p
