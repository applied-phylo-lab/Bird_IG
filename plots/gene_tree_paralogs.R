library(dplyr)
library(ggtree)      # Bioconductor
library(ggplot2)
library(ggnewscale)  # allow multiple fill scales
library(RColorBrewer)
library(stringr)
library(tidyr)
library(seqinr)
library(ggtree)
library(ggnewscale)
bird_dir<-"/local/storage/kav67/birds/Hummingbirds/Annas_Hummingbird/bCalAnn1_pri/"
iggenes<-fread(paste0(bird_dir,"combined_genes_IGH.txt"))
igsummary<-fread(paste0(bird_dir,"refined_ig_loci/summary.csv"))


iggenes_main <- iggenes %>%
  group_by(Contig) %>%
  mutate(n_rows = n()) %>%
  ungroup() %>%
  filter(Contig == Contig[which.max(n_rows)]) %>%
  select(-n_rows) %>%
  
  # join to get StartPos, removing unmatched Contigs
  inner_join(igsummary %>% select(Contig, StartPos), by = "Contig") %>%
  
  # calculate start and end
  mutate(
    start = Pos - StartPos,
    end = start + nchar(Sequence)
  ) %>%
  select(-StartPos) %>%
  # create name column
  mutate(name = paste0(start, "_", end, "_", Strand))


iggenes_main$tree_label<-paste0(iggenes_main$Contig,"_",iggenes_main$Pos)

write.fasta(sequences = as.list(iggenes_main$Sequence), names = iggenes_main$tree_label, file.out = paste0(bird_dir,"IGH_genes.fasta"), open = "w", nbchar = 60)


inversion_d<-fread("/local/storage/kav67/birds/IGH_paralogs_all.tsv")


# 1. Filter inversion_d for the right haplotype
inv_filtered <- inversion_d %>%
  filter(Haplotype == "bCalAnn1_pri")

# 2. Expand group_members into individual rows
inv_long <- inv_filtered %>%
  mutate(group_id = row_number()) %>%
  separate_rows(group_members, sep = ";") %>%
  mutate(
    # remove trailing "_0" or "_1" etc to match the format in iggenes_main$name
    clean_name = str_replace(group_members, "_\\d+$", "")
  ) %>%
  select(clean_name, group_id)

# 3. Add group_id to iggenes_main by matching to name
iggenes_main <- iggenes_main %>%
  left_join(inv_long, by = c("name" = "clean_name"))

  
gene_tree<-read.tree(file = paste0(bird_dir,"IGH_genes.treefile"))

annot <- iggenes_main %>%
  select(tree_label, group_id, Strand, Productive) %>%
  distinct(tree_label, .keep_all = TRUE) %>%
  mutate(tree_label = as.character(tree_label))

# ensure tree labels exist in annot; order annot to match tree tip labels
# make rownames = tip labels as required by gheatmap
annot_df <- as.data.frame(annot)
rownames(annot_df) <- annot_df$tree_label

#reorder to match tree tip order (this will insert NA rows for missing tips)
annot_df <- annot_df[gene_tree$tip.label, , drop = FALSE]

# --- determine singletons directly from iggenes_main --------------------------
# count non-NA occurrences of each group_id in iggenes_main
group_counts <- table(iggenes_main$group_id, useNA = "no")
singleton_ids <- names(group_counts[group_counts == 1])

# create factor for plotting; keep levels consistent
annot_df$group_id <- factor(annot_df$group_id, levels = names(group_counts))



n_groups <- length(levels(annot_df$group_id))
non_singleton_ids <- setdiff(levels(annot_df$group_id), singleton_ids)

# generate enough colors for non-singleton groups
if (length(non_singleton_ids) > 0) {
  palette_colors <- colorRampPalette(brewer.pal(8, "Set2"))(length(non_singleton_ids))
  group_colors <- setNames(palette_colors, non_singleton_ids)
} else {
  group_colors <- character(0)
}

# set singletons to white
if (length(singleton_ids) > 0) {
  # ensure singletons are included in the named vector (white)
  group_colors[singleton_ids] <- "white"
}

# For safety: ensure all factor levels have an entry; set any missing to a default (light grey)
missing_levels <- setdiff(levels(annot_df$group_id), names(group_colors))
if (length(missing_levels) > 0) {
  group_colors[missing_levels] <- "grey90"
}

# --- prepare separate small data.frames for gheatmap --------------------------
df_group  <- data.frame(group_id = annot_df$group_id, row.names = rownames(annot_df))
df_strand <- data.frame(Strand   = factor(annot_df$Strand), row.names = rownames(annot_df))
df_prod   <- data.frame(Productive = factor(annot_df$Productive), row.names = rownames(annot_df))

# --- plotting: circular tree + 3 concentric heatmap rings ---------------------


# Start with circular tree
p <- ggtree(gene_tree, layout = "circular")

# ---- 1. group_id ring (hide legend, singletons white) ----
p <- gheatmap(
  p, df_group,
  aes(fill = group_id),
  offset = 0.02, width = 0.1,
  colnames = FALSE
) +
  scale_fill_manual(
    values = group_colors,
    na.value = "white",
    guide = "none"
  )
p
# ---- 2. Strand ring ----
p <- ggtree(gene_tree, layout = "circular")
p1 <- p + new_scale_fill()
p1 <- gheatmap(
  p, df_strand,
  offset = 0.13, width = 0.1,
  colnames = FALSE
) +
  scale_fill_manual(
    values = c("+" = "#1f78b4", "-" = "#e31a1c"),
    name = "Strand"
  )
p
p2 <- p1 + new_scale_fill()
# ---- 3. Productive ring ----
p <- ggtree(gene_tree, layout = "circular")
p <- gheatmap(
  p, df_prod,
  offset = 0.19, width = 0.1,
  colnames = FALSE
) +
  scale_fill_manual(
    values = c("TRUE" = "#33a02c", "FALSE" = "#b2df8a"),
    name = "Productive"
  )

p
