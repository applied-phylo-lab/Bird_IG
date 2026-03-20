summary_table_curated
sub_tree

library(geiger)
library(dplyr)
library(ape)
library(ggplot2)

# Prepare data frame
locus<-"IGH"

summary_table_IGH <- summary_table_curated %>%
  filter(Locus == locus) %>%            # only IGH rows
  group_by(LatinName) %>%               # group by species
  slice_max(order_by = NumV, n = 1) %>% # keep row with highest NumV
  ungroup()


trait_df <- summary_table_IGH %>%
  select(LatinName, Order, NumProdV)


# Named trait vector
trait_vec <- setNames(trait_df$NumProdV, trait_df$LatinName)

# Keep only species present in the tree
trait_vec <- trait_vec[names(trait_vec) %in% sub_tree$tip.label]

orders <- unique(trait_df$Order)


# Empty results table
bm_results <- data.frame(
  Order = character(),
  Nspecies = numeric(),
  sigma2 = numeric(),
  logLik = numeric(),
  stringsAsFactors = FALSE
)

# ---- Loop over Orders ----

for (ord in orders) {
  
  spp <- trait_df %>% filter(Order == ord) %>% pull(LatinName)
  spp <- intersect(spp, sub_tree$tip.label)
  
  # skip tiny clades
  if (length(spp) < 3) {
    message("Skipping ", ord, " (too few species)")
    next
  }
  
  # order-specific subtree
  st <- drop.tip(sub_tree, setdiff(sub_tree$tip.label, spp))
  
  # matching trait vector
  v <- trait_vec[spp]
  
  # geiger requires names match exactly
  v <- v[st$tip.label]
  
  # Fit Brownian Motion model
  fit <- fitContinuous(st, v, model = "BM")
  
  # Extract parameters
  sigma2 <- fit$opt$sigsq  # the BM rate
  logL   <- fit$opt$lnL
  
  # Save results
  bm_results <- rbind(
    bm_results,
    data.frame(
      Order = ord,
      Nspecies = length(spp),
      sigma2 = sigma2,
      logLik = logL
    )
  )
}

bm_results


ggplot(bm_results, aes(x = Order, y = sigma2)) +
  geom_point(size = 3) +
  geom_segment(aes(x = Order, xend = Order, y = 0, yend = sigma2),
               alpha = 0.4) +
  theme_bw() +
  labs(
    title = "Brownian Motion Rate (σ²) of NumV per Order",
    y = "BM Rate (σ²)",
    x = "Order"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# keep only orders that were actually fitted
sigma_vec <- setNames(bm_results$sigma2, bm_results$Order)
orders_in_tree <- order_tree$tip.label

# match rates to tree order
sigma_mapped <- sigma_vec[orders_in_tree]
trait_for_plot <- sigma_mapped
names(trait_for_plot) <- orders_in_tree

# contMap expects no NA's; set NA=0 or mean
trait_for_plot_clean <- trait_for_plot
trait_for_plot_clean[is.na(trait_for_plot_clean)] <- 0  # grey branches later

cm <- contMap(order_tree, trait_for_plot_clean, plot=FALSE)

# Now manually grey out branches belonging to missing orders
missing_orders <- names(trait_for_plot)[is.na(trait_for_plot)]
missing_edges <- which(order_tree$tip.label %in% missing_orders)

# Plot the tree
plot(cm, lwd=5, fsize=0.9)
title("BM Rate (σ²) per Order")

# Label missing orders on the tree
tiplabels(pch=19, col="grey40", cex=1,
          tip = which(order_tree$tip.label %in% missing_orders))



tree_df <- data.frame(
  label = order_tree$tip.label,
  stringsAsFactors = FALSE
)

# Join σ² onto the tree tip labels
tree_df <- tree_df %>%
  left_join(bm_results, by = c("label" = "Order"))


p <- ggtree(order_tree) %<+% tree_df +
  geom_tree() +
  geom_tippoint(aes(color = sigma2), size = 4) +
  scale_color_viridis_c(option = "plasma", na.value = "grey80") +
  theme_tree2() +
  labs(color = "BM rate (σ²)",
       title = "Brownian Motion Rates per Order")


tree_df[is.na(tree_df)] <- 0  # grey branches later

p <- ggtree(order_tree, layout = "rectangular") +
  geom_tiplab(size = 3, align = TRUE, linetype = NA, linesize = 0.5)

p2<-facet_plot(p+xlim_tree(9), panel = "Brownian Motion Rate",
               data = tree_df,
               geom_point,
               mapping = aes(x = sigma2, group = label))
p2
