# MinDir figure -- fraction of V genes on the majority strand, per species and
# locus, mapped onto the bird phylogeny.
#
# GENE_LIST selects which gene table MinDir is computed from:
#   "current"  gene_list.csv            -- 82,451 rows (2026-08-06)
#   "legacy"   gene_list-igl_h3_n7.csv  -- 56,463 rows (2026-05-07), what the
#              manuscript figure was made from
# Identical headers, so the two are drop-in interchangeable. Figure filenames
# carry the choice so both can be compared side by side.
#
# Figures are exported by hand for the manuscript -- the save_fig calls are for
# unattended reproduction, not final artwork.

library(data.table)
library(dplyr)
library(ggplot2)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(viridis)
library(tidyr)

source(file.path("/home/kav67/Bird_IG", "plots", "_dataset.R"))

# ---- switch this line ----
# "current" is the pipeline setting; "legacy" reproduces the older figure.
GENE_LIST <- "current"       # "legacy" (manuscript) or "current"
# --------------------------

GENE_LIST_FILE <- switch(GENE_LIST,
  legacy  = file.path(INPUT_DIR, "gene_list-igl_h3_n7.csv"),
  current = file.path(INPUT_DIR, "gene_list.csv"),
  stop("GENE_LIST must be 'legacy' or 'current', got: ", GENE_LIST))
if (!file.exists(GENE_LIST_FILE)) stop("Gene list not found: ", GENE_LIST_FILE)
message(sprintf("[mindir] gene list: %s (%s)", basename(GENE_LIST_FILE), GENE_LIST))

# ── Load data ──────────────────────────────────────────────────────────────────

# all_species_data <- fread("/local/storage/kav67/clean_birds/all_species_stats_pruned_12052025.csv")
# all_species_data$bird <- FALSE
# all_species_data[all_species_data$VertClass == "birds", ]$bird <- TRUE
# all_species_data <- all_species_data[all_species_data$IGH_AnnotationLevel < 2, ]
# all_species_data <- all_species_data[all_species_data$VertClass %in% c("reptiles", "mammals", "birds"), ]
# bird_only <- all_species_data[all_species_data$VertClass == "birds", ]

gene_list <- fread(GENE_LIST_FILE)
gene_list[, c("GrpOrder", "Species", "Haplotype") := tstrsplit(Source, "/", fixed = TRUE)]

main_strand_frac <- function(strands) {
  tbl <- table(strands)
  max(tbl) / sum(tbl)
}

# ── How MinDir is aggregated across contigs ──────────────────────────────────
# A haplotype's locus can be split over several contigs. MinDir is a per-contig
# strand fraction, so how those are combined matters more than it looks: a contig
# carrying ONE gene scores 1.0 by construction, and 40 of 904 contigs (4.4%) in
# gene_list.csv are exactly that. Under "mean" they count as much as a 200-gene
# locus. That is what pushed the Red-billed Tropicbird IGH from 0.50 to 0.83 --
# two scrap contigs of 2 and 1 genes outvoting its real 18-gene locus 2:1.
#
#   mean       unweighted mean over contigs (original behaviour, and what the
#              manuscript figure used)
#   weighted   mean over contigs weighted by gene count
#   pooled     pool every gene of a haplotype x locus, then take one fraction
#   min_genes  drop contigs with fewer than MIN_GENES genes, then unweighted mean
# "weighted" is the pipeline setting. The other three stay available for
# comparison -- switch, re-run, and the output filename records the choice -- but
# a normal pipeline run should always be weighted.
AGGREGATION <- "weighted"  # "mean" | "weighted" | "pooled" | "min_genes"
MIN_GENES   <- 5           # only used when AGGREGATION == "min_genes"

# Size of the numeric tick labels under the bar panels.
AXIS_TEXT_SIZE <- 3.2

per_contig <- gene_list %>%
  group_by(GrpOrder, Species, Haplotype, Locus, Contig) %>%
  summarise(MinDir = main_strand_frac(Strand), nGenes = n(), .groups = "drop")

per_haplotype <- switch(AGGREGATION,
  mean = per_contig %>%
    group_by(GrpOrder, Species, Haplotype, Locus) %>%
    summarise(MinDir = mean(MinDir), .groups = "drop"),

  weighted = per_contig %>%
    group_by(GrpOrder, Species, Haplotype, Locus) %>%
    summarise(MinDir = weighted.mean(MinDir, nGenes), .groups = "drop"),

  pooled = gene_list %>%
    group_by(GrpOrder, Species, Haplotype, Locus) %>%
    summarise(MinDir = main_strand_frac(Strand), .groups = "drop"),

  min_genes = per_contig %>%
    filter(nGenes >= MIN_GENES) %>%
    group_by(GrpOrder, Species, Haplotype, Locus) %>%
    summarise(MinDir = mean(MinDir), .groups = "drop"),

  stop("AGGREGATION must be one of mean/weighted/pooled/min_genes, got: ",
       AGGREGATION))

message(sprintf("[mindir] aggregation: %s%s | contigs %d -> haplotype x locus %d",
                AGGREGATION,
                if (AGGREGATION == "min_genes") sprintf(" (>=%d genes)", MIN_GENES) else "",
                nrow(per_contig), nrow(per_haplotype)))

bird_only <- per_haplotype %>%
  group_by(Species, Locus) %>%
  summarise(MinDir = mean(MinDir), .groups = "drop") %>%
  pivot_wider(names_from = Locus, values_from = MinDir, names_glue = "{Locus}_MinDir") %>%
  mutate(LatinName = gsub("_", " ", Species),
         latin_tree = tolower(Species)) %>%
  as.data.frame()

bird_tree_pruned <- read.tree(TREE_FILE)
tree_species <- gsub('.*"', '', bird_tree_pruned$tip.label)
tree_species <- gsub('"', '', tree_species)
bird_tree_pruned$species <- tolower(tree_species)

# Join Order and proper LatinName from VGP table
vgp_table <- fread(VGP_TABLE)
order_map  <- vgp_table %>%
  select(LatinName, Order,Species) %>%
  distinct() %>%
  mutate(latin_tree = tolower(gsub(" ", "_", Species)))

# ── Match to tree ──────────────────────────────────────────────────────────────

# Join on latin_tree (lowercase underscore) — also brings in proper LatinName
bird_only <- merge(as.data.frame(bird_only %>% select(!LatinName)),
                   order_map %>% select(latin_tree, Order, LatinName),
                   by = "latin_tree", all.x = TRUE)

idx      <- match(tree_species, gsub(" ","_",bird_only$LatinName))
tip_meta <- bird_only[idx, ]

keep         <- !is.na(tip_meta$LatinName)
tip_meta     <- tip_meta[keep, ]
tree_pruned  <- drop.tip(bird_tree_pruned, bird_tree_pruned$tip.label[!keep])

# Rename tip labels to LatinName so geom_fruit can match on y = LatinName
latin_map              <- setNames(tip_meta$LatinName, tree_pruned$tip.label)
tree_pruned$tip.label  <- latin_map[tree_pruned$tip.label]

# ── Order MRCA nodes for highlighting ─────────────────────────────────────────

order_nodes <- tip_meta %>%
  filter(!is.na(Order), Order != "MiscBirds") %>%
  group_by(Order) %>%
  summarise(tips = list(LatinName), n = n(), .groups = "drop") %>%
  filter(n >= 2) %>%
  mutate(node = sapply(tips, function(t) {
    t_in <- t[t %in% tree_pruned$tip.label]
    if (length(t_in) < 2) NA_integer_
    else getMRCA(tree_pruned, t_in)
  })) %>%
  filter(!is.na(node))

# ── Base tree ─────────────────────────────────────────────────────────────────

p <- ggtree(tree_pruned, layout = "rectangular") +
  geom_tiplab(size = 1.8, fontface = "italic", offset = 0.001)

# No order highlighting: geom_cladelab below already writes the order beside the
# clade, so the colour band was redundant.

p <- p + geom_cladelab(
  data    = order_nodes,
  mapping = aes(node = node, label = Order),
  align   = TRUE,
  offset  = 0.002,
  fontsize = 2.5,
  hjust   = 0
)

# ── IGH_MinDir panel ──────────────────────────────────────────────────────────

p <- p + new_scale_fill()

p <- p + geom_fruit(
  data    = tip_meta,
  geom    = geom_bar,
  mapping = aes(y = LatinName, x = IGH_MinDir, fill = IGH_MinDir),
  stat    = "identity",
  orientation = "y",
  offset  = 0.15,
  pwidth  = 0.25,
  axis.params = list(axis = "x", text.size = AXIS_TEXT_SIZE, title = "IGH MinDir",
                     title.size = 3, title.height = 0.02)
) +
  scale_fill_viridis_c(name = "IGH\nMinDir", option = "D", na.value = "grey90",
                       direction = -1, limits = c(0.5, 1))

# ── IGL_MinDir panel ──────────────────────────────────────────────────────────

p <- p + new_scale_fill()

p <- p + geom_fruit(
  data    = tip_meta,
  geom    = geom_bar,
  mapping = aes(y = LatinName, x = IGL_MinDir, fill = IGL_MinDir),
  stat    = "identity",
  orientation = "y",
  offset  = 0.05,
  pwidth  = 0.25,
  axis.params = list(axis = "x", text.size = AXIS_TEXT_SIZE, title = "IGL MinDir",
                     title.size = 3, title.height = 0.02)
) +
  scale_fill_viridis_c(name = "IGL\nMinDir", option = "D", na.value = "grey90",
                       direction = -1, limits = c(0.5, 1))

# MinDir is max(tbl)/sum(tbl), so it cannot fall below 0.5 -- the observed range
# is exactly 0.5000 to 1.0000. Both panels use viridis "D" reversed over that
# range, so 1.0 is dark purple and 0.5 is yellow. IGH previously used option "C"
# (plasma) against c(0, 1), which both clashed with the IGL panel and spent half
# the colour scale on values MinDir can never take.
p

# ── Butterfly version: IGH and IGL mirrored about zero ───────────────────────
# This is the manuscript figure.

p2 <- ggtree(tree_pruned, layout = "rectangular") #+
  #geom_tiplab(size = 1.8, fontface = "italic", offset = 0.001)

# No geom_hilight here: the order is already written beside the clade by
# geom_cladelab below, so the colour band only added noise.

p2 <- p2 + geom_cladelab(
  data    = order_nodes,
  mapping = aes(node = node, label = Order),
  align   = TRUE,
  offset  = 0.002,
  fontsize = 2.5,
  hjust   = 0
)
p2

butterfly_data <- tip_meta %>%
  select(LatinName, IGH_MinDir, IGL_MinDir) %>%
  pivot_longer(c(IGH_MinDir, IGL_MinDir), names_to = "Locus", values_to = "MinDir") %>%
  mutate(x_pos = if_else(Locus == "IGH_MinDir", -MinDir, MinDir),
         bar_alpha = 1)

# Transparent ghost bar at x = -1 using a real tree species so ggtreeExtra
# includes it in scale computation, extending the IGH side to -1
ghost_sp <- tip_meta$LatinName[which.max(tip_meta$IGH_MinDir)]
butterfly_data <- butterfly_data %>%
  add_row(LatinName = ghost_sp, Locus = "IGH_MinDir", MinDir = 1,
          x_pos = -1, bar_alpha = 0)

p2 <- p2 + new_scale_fill()

p2 <- p2 + geom_fruit(
  data    = butterfly_data,
  geom    = geom_bar,
  mapping = aes(y = LatinName, x = x_pos, fill = MinDir, alpha = bar_alpha),
  stat    = "identity",
  orientation = "y",
  offset  = 0.7,
  pwidth  = 0.5,
  # The x values are signed to mirror the two loci, but the tick labels should
  # read as magnitudes -- MinDir has no negative values. `text` overrides the
  # labels; line.alpha = 0 hides the axis rule under the bars.
  # line.alpha = 0 hides the axis rule that ran under the bars. The tick labels
  # are fixed up afterwards by mirror_axis_labels(); axis.params$text only takes
  # effect when there is a single break (see ggtreeExtra:::build_axis), so it
  # cannot relabel a mirrored axis.
  axis.params = list(axis = "x", text.size = AXIS_TEXT_SIZE,
                     title = "                        \u2190 IGH | IGL \u2192",
                     title.size = 3, title.height = 0.02,
                     limits = c(-1, 1.1),
                     line.alpha = 0)
) +
  scale_fill_viridis_c(name = "MinDir", option = "D", na.value = "grey90",direction=-1,
                       limits = c(0.5, 1)) +
  scale_alpha_identity(guide = "none")

# The butterfly mirrors IGH onto negative x, so ggtreeExtra labels those ticks
# -1, -0.8 ... MinDir has no negative values, and the sign only encodes which
# locus a bar belongs to, which the axis title already says. This rewrites the
# tick layer to the breaks we want, labelled by magnitude.
#
# The axis layer carries the label in one column and the plotting position in
# another (aes(x = new_<xid>, label = <xid>)), so the two can be set
# independently: label with the magnitude, position with the signed value.
mirror_axis_labels <- function(plot, breaks = c(-1, -0.5, 0, 0.5, 1),
                               digits = 2) {
  for (i in seq_along(plot$layers)) {
    d <- plot$layers[[i]]$data
    if (!is.data.frame(d) || ncol(d) != 2) next
    nm <- names(d)
    pos_col <- grep("^new_", nm, value = TRUE)
    if (length(pos_col) != 1) next
    lab_col <- setdiff(nm, pos_col)
    nz <- d[[lab_col]] != 0
    if (!any(nz)) next
    slope <- mean(d[[pos_col]][nz] / d[[lab_col]][nz])
    new_d <- data.frame(abs(round(breaks, digits)), breaks * slope)
    names(new_d) <- c(lab_col, pos_col)
    plot$layers[[i]]$data <- new_d
    return(plot)
  }
  warning("mirror_axis_labels(): no axis layer found; labels left unchanged")
  plot
}

p2 <- mirror_axis_labels(p2)
p2

# ---- save -------------------------------------------------------------------
# Only the butterfly (p2) is saved -- that is the manuscript figure. `p`, the
# two-panel IGH/IGL view, is still built above and prints to the viewer, but the
# pipeline should not write a figure nobody uses. Export it by hand if wanted:
#   save_fig(sprintf("mindir_tree_panels%s.svg", .tag), p, width = 9, height = 12)
#
# The filename carries GENE_LIST and AGGREGATION so a figure made with a
# non-default combination can never be mistaken for the pipeline's own.
.tag <- sprintf("_%s_%s", GENE_LIST, AGGREGATION)
save_fig(sprintf("mindir_tree%s.svg", .tag), p2, width = 9, height = 12)

cat(sprintf("[mindir] species on tree: %d | gene list: %s | dataset: %s\n",
            length(tree_pruned$tip.label), GENE_LIST, DATASET))
