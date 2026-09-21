
library(data.table)
library(dplyr)
library(phytools)
library(viridis)
library(ape)
library(phylolm)


summary<-fread("/local/storage/kav67/clean_birds/IGH_VGP_table.tsv")
inversions<-fread("/local/storage/kav67/clean_birds/inversion_stats.tsv")
#inversions<-inversions[inversions$minlen==1000,]
bird_tree_pruned<-read.tree("/local/storage/kav67/clean_birds/vgp_birds.nwk")

# IGH_VGP_table.tsv keeps only the main (highest-NumV) contig per haplotype, while
# inversion_stats.tsv has one row per contig. Joining on Haplotype alone fans out
# (128 -> 140 rows) and then pairs a haplotype's main-contig V gene count with
# inversion counts averaged over contigs whose genes were never counted. Join on
# the contig as well so each haplotype contributes its main contig's V genes AND
# that same contig's inversions.
setnames(inversions, c("haplotype", "contig"), c("Haplotype", "Contig"))
summary <- summary %>%
  group_by(LatinName, Haplotype) %>%
  slice_max(NumV, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  inner_join(inversions, by = c("Haplotype", "Contig"))

# Reorder to match tree
tree_species <- gsub('.*"', '', bird_tree_pruned$tip.label)
tree_species <- gsub('"', '', tree_species)
bird_tree_pruned$species <- tree_species
# Aggregate traits and keep Species name
trait_df <- summary %>%
  group_by(LatinName) %>%
  summarise(
    NumV = log(mean(NumV, na.rm=TRUE)),
    num_inversions = log(mean(num_inversions, na.rm=TRUE)),
    .groups="drop"
  )
#filter(LatinName != "Gallinula chloropus") %>%   # remove outlier
#    Species = first(Species),
# Match order to tree
trait_df <- trait_df[match(tree_species, gsub(" ","_",trait_df$LatinName)), ]

# Remove species with missing values
keep <- !(is.na(trait_df$NumV) | is.na(trait_df$num_inversions))

# label tips BEFORE dropping, so the labels still line up with all 122 rows
trait_df$tip_label <- bird_tree_pruned$tip.label
trait_df <- as.data.frame(trait_df)[keep, ]
rownames(trait_df) <- trait_df$tip_label

tree_plot <- drop.tip(bird_tree_pruned, bird_tree_pruned$tip.label[!keep])

trait_df_model <- trait_df[tree_plot$tip.label, ]

# Run phylogenetic linear model (lambda model)
model <- phylolm(
  NumV ~ num_inversions,
  data = trait_df_model,
  phy = tree_plot,
  model = "lambda"
)

# Show results
summary(model)

# Extract useful values
lambda_estimate <- model$optpar
p_value <- summary(model)$coefficients["num_inversions","p.value"]
slope <- summary(model)$coefficients["num_inversions","Estimate"]

cat("\nLambda:", lambda_estimate,
    "\nSlope:", slope,
    "\nP-value:", p_value, "\n")




# Named vectors for contMap
NumV_values <- setNames(trait_df$NumV, tree_plot$tip.label)
inv_values  <- setNames(trait_df$num_inversions, tree_plot$tip.label)

# Clean species names for labels
# trait_df is built by summarise() and keeps only LatinName, NumV and
# num_inversions -- the `Species = first(Species)` line that would have carried a
# Species column is commented out just above. Using trait_df$Species here gave
# character(0) and text() then failed on the length mismatch; this was masked
# until phytools was installed, because the script died at library(phytools)
# first.
species_labels <- gsub("_", " ", trait_df$LatinName)

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




library(ggplot2)

# Extract coefficients
coef_vals <- coef(model)
intercept <- coef_vals[1]
slope <- coef_vals[2]

# Add fitted values (optional but useful)
trait_df_model$fitted <- fitted(model)

# Extract lambda
lambda_val <- model$optpar

# Get p-value
p_val <- summary(model)$coefficients[2,4]

# Compute R² (pseudo-R² using correlation of fitted vs observed)
r2 <- cor(trait_df_model$num_inversions,
          trait_df_model$fitted,
          use = "complete.obs")^2
r2_phylo <- summary(model)$r.squared
# Significance stars
stars <- ifelse(p_val < 0.001, "***",
                ifelse(p_val < 0.01, "**",
                       ifelse(p_val < 0.05, "*", "ns")))
slope <- coef(model)[2]
# Create label
label_text <- paste0("β = ", round(slope, 3), stars,
                     "\n λ  = ", round(lambda_val, 3),"    ",
                     "\nR² = ", round(r2_phylo, 3),"   ")

slope_label <- sprintf("%.3f%s", round(slope, 3), stars)

# Then align everything
label_text <- sprintf(
  "β  = %-8s\nλ  = %-8.5f\nR² = %-8.3f",
  slope_label,
  lambda_val,
  round(r2_phylo, 3)
)
# Plot
p_phylo <- ggplot(trait_df_model,
                  aes(x = num_inversions, y = NumV)) +
  geom_point(size = 2, alpha = 0.9, color="#87b4dc") +
  
  # Regression line from model
  geom_abline(intercept = intercept,
              slope = slope,
              linetype = "dashed",
              size = 0.5) +
  
  # Annotation
  annotate("text",
           x = 0, y = Inf,
           label = label_text,
           hjust = -0.5, vjust = 2,
           size = 4) +
  
  labs(
       x = "# Inversions (log)",
       y = "# V genes (log)") +
  
  theme_classic()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10))

p_phylo

