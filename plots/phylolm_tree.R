
library(dplyr)
library(phytools)
library(viridis)
library(ape)
library(phylolm)


summary<-fread("/local/storage/kav67/clean_birds/IGH_filtered_table.tsv")
inversions<-fread("/local/storage/kav67/clean_birds/IGH_inversions.tsv")
inversions<-inversions[inversions$minlen==1000,]

colnames(inversions)[16]<-"Haplotype"
summary<-left_join(summary,inversions, by="Haplotype")

summary<-unique(summary)

# Reorder to match tree
tree_species <- gsub('.*"', '', bird_tree_pruned$tip.label)
tree_species <- gsub('"', '', tree_species)
bird_tree_pruned$species <- tree_species
# Aggregate traits and keep Species name
trait_df <- summary %>%
  filter(LatinName != "Gallinula chloropus") %>%   # remove outlier
  group_by(LatinName) %>%
  summarise(
    NumV = mean(NumV, na.rm=TRUE),
    num_inversions_diag = mean(num_inversions_diag, na.rm=TRUE),
    Species = first(Species),
    .groups="drop"
  )

# Match order to tree
trait_df <- trait_df[match(tree_species, trait_df$LatinName), ]

# Remove species with missing values
keep <- !(is.na(trait_df$NumV) | is.na(trait_df$num_inversions_diag))

trait_df <- trait_df[keep, ]


trait_df$tip_label <- bird_tree_pruned$tip.label
rownames(trait_df) <- bird_tree_pruned$tip.label

tree_plot <- drop.tip(bird_tree_pruned, bird_tree_pruned$tip.label[!keep])

trait_df_model <- trait_df
rownames(trait_df_model) <- tree_plot$tip.label

# Run phylogenetic linear model (lambda model)
model <- phylolm(
  num_inversions_diag ~ NumV,
  data = trait_df_model,
  phy = tree_plot,
  model = "lambda"
)

# Show results
summary(model)

# Extract useful values
lambda_estimate <- model$optpar
p_value <- summary(model)$coefficients["NumV","p.value"]
slope <- summary(model)$coefficients["NumV","Estimate"]

cat("\nLambda:", lambda_estimate,
    "\nSlope:", slope,
    "\nP-value:", p_value, "\n")




# Named vectors for contMap
NumV_values <- setNames(trait_df$NumV, tree_plot$tip.label)
inv_values  <- setNames(trait_df$num_inversions_diag, tree_plot$tip.label)

# Clean species names for labels
species_labels <- gsub("_", " ", trait_df$Species)

# Layout
layout(matrix(1:3,1,3), widths=c(0.4,0.2,0.4))
par(cex=1)

# ----- NumV contMap -----
obj <- contMap(tree_plot, NumV_values, outline=FALSE, plot=FALSE)
contmap_obj_viridis <- setMap(obj, rev(viridis(100)))

plot(contmap_obj_viridis,
     ftype=c("off","reg"),
     leg.txt="",
     legend=20,
     mar=c(1.1,0.1,4.1,0.1))
title(main="# V genes")

# ----- Middle species labels -----
ylim <- c(1-0.12*(length(tree_plot$tip.label)-1), length(tree_plot$tip.label))

plot.new()
plot.window(xlim=c(-0.1,0.1), ylim=ylim)

text(rep(0,length(tree_plot$tip.label)),
     1:length(tree_plot$tip.label),
     species_labels,
     font=3)

# ----- inversion contMap -----
obj <- contMap(tree_plot, inv_values, outline=FALSE, plot=FALSE)
contmap_obj_viridis <- setMap(obj, rev(viridis(100)))

plot(contmap_obj_viridis,
     ftype=c("off","reg"),
     direction="leftwards",
     leg.txt="",
     legend=20,
     mar=c(1.1,0.1,4.1,0.1))
title(main="# inversions")
#dev.off()
