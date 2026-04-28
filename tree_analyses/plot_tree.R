#!/usr/bin/env Rscript
# Usage: Rscript plot_tree.R <treefile> <color_tsv> <output_svg>
# color_tsv: two-column TSV with headers tip_name, color

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(svglite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript plot_tree.R <treefile> <color_tsv> <output_svg>")
}

treefile   <- args[1]
color_tsv  <- args[2]
output_svg <- args[3]

tree   <- read.tree(treefile)
tree   <- ladderize(tree, right = FALSE)
colors <- read.delim(color_tsv, stringsAsFactors = FALSE)

color_map            <- setNames(colors$color, colors$tip_name)
tip_df               <- data.frame(label = tree$tip.label, stringsAsFactors = FALSE)
tip_df$color         <- color_map[tip_df$label]
tip_df$color[is.na(tip_df$color)] <- "#999999"

n_tips <- length(tree$tip.label)
height <- max(4, n_tips * 0.18)

p <- ggtree(tree, layout = "rectangular") %<+% tip_df +
  geom_tiplab(aes(color = color), size = 2, align = FALSE) +
  scale_color_identity() +
  theme_tree2() +
  theme(legend.position = "none")

ggsave(output_svg, plot = p, device = svglite::svglite,
       width = 10, height = height, limitsize = FALSE)
