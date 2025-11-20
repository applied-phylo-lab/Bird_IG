library(ggplot2)
library(dplyr)
library(tidyr)
library(phytools)
paralogs<-fread("/local/storage/kav67/Bird_data/IGH_paralogs_all.tsv")

summary_table_IGH<-fread("/local/storage/kav67/Bird_data/IGH_table.tsv")

for( o in unique(paralogs$Order)){
  p<-ggplot(paralogs[paralogs$Order==o], aes(x=group_size))+
    geom_histogram()+
    theme_minimal()+
    labs(title =o)
  print(p)
}

ggplot(paralogs, aes(x = group_size, fill = Order)) +
  geom_bar(position = "stack") +
  scale_y_log10()+
  labs(title ="Paralog Group Sizes",
       x = "Group Size",
       y = "Count") +
  theme_minimal()

ggplot(paralogs[paralogs$group_size>1,], aes(x = group_size, fill = Order)) +
  geom_bar(position = "stack") +
  scale_y_log10()+
  labs(title ="Paralog Group Sizes",
       x = "Group Size",
       y = "Count") +
  theme_minimal()

fractions <- data.table(paralogs %>%
  group_by(Order, Species, Haplotype, sample) %>%
  summarise(
    total_genes = n(),
    grouped_genes = sum(group_size > 1),
    fraction_grouped = grouped_genes / total_genes,
    .groups = "drop"
  ))%>%
  rename(label = Order)  


p <- ggtree(order_tree, layout = "rectangular") +
  geom_tiplab(size = 3, align = TRUE, linetype = NA, linesize = 0.5)


# Add boxplots per order
p2<-facet_plot(p+xlim_tree(9), panel = "Fraction of genes in paralog groups (>1)",
               data = fractions,
               geom_boxplot,
               mapping = aes(x = fraction_grouped, fill = label))
p2+theme_tree2()


paralog_groups <- paralogs %>% #[paralogs$group_size>1,]
  distinct(Order, Species, Haplotype, sample, group_members, group_size)


group_sizes_summary <- data.table(paralog_groups) %>%
  rename(label = Order)
group_sizes_summary$log_group_size<-log(group_sizes_summary$group_size)
p3 <- facet_plot(
  p + xlim_tree(9),
  panel = "Paralog group sizes (log)",
  data = group_sizes_summary,
  geom_boxplot,
  mapping = aes(x = log_group_size, fill = label)
) +
  theme_tree2()
p3



paralog_groups<-paralogs[paralogs$group_size!=1,]

median_distance<- function(row){
  group<-strsplit(row["group_members"], ";")
  split_list <- strsplit(group$group_members, split = "_")
  df <- as.data.frame(do.call(rbind, split_list))
  negative<-df[df$V3=="-",]
  positive<-df[df$V3=="+",]
  combinations_df_n <- expand.grid(start1 = as.numeric(negative$V1), start2 = as.numeric(negative$V1))
  combinations_df_p <- expand.grid(start1 = as.numeric(positive$V1), start2 = as.numeric(positive$V1))
  combinations_df<-rbind(combinations_df_n,combinations_df_p)
  combinations_df$distance<-abs(combinations_df$start2-combinations_df$start1)
  return(median(combinations_df$distance))
}

paralog_groups$median_distance <- apply(paralog_groups, 1, median_distance)

paralog_summary <- data.table(paralog_groups) %>%
  rename(label = Order)

p4 <- facet_plot(
  p + xlim_tree(9),
  panel = "Median Distance between paralogs on the same strand",
  data = paralog_summary,
  geom_boxplot,
  mapping = aes(x = median_distance, fill = label)
) +
  theme_tree2()
p4


# histogram for one:

# Filter for the species of interest
df_sub <- paralog_groups %>% 
  filter(Species == "House_Finch")

# Function to extract ALL same-strand pairwise distances for a group
get_distances <- function(row) {
  group <- strsplit(row["group_members"], ";")[[1]]
  
  # split each member into components
  split_list <- strsplit(group, "_")
  df <- as.data.frame(do.call(rbind, split_list), stringsAsFactors = FALSE)
  colnames(df) <- c("start", "end", "strand")
  
  df$start <- as.numeric(df$start)
  
  # separate ± strands
  negative <- df[df$strand == "-", ]
  positive <- df[df$strand == "+", ]
  
  # expand combinations within each strand (all pairwise combinations)
  comb_neg <- expand.grid(start1 = negative$start, start2 = negative$start)
  comb_pos <- expand.grid(start1 = positive$start, start2 = positive$start)
  
  # combine
  comb <- rbind(comb_neg, comb_pos)
  
  # compute distances
  comb$distance <- abs(comb$start2 - comb$start1)
  
  return(comb$distance)
}

# Apply to all groups for Common Myna
all_distances <- unlist(apply(df_sub, 1, get_distances))

# Remove zeros if you don't want self-distance
all_distances <- all_distances[all_distances > 0]

dist_df <- data.frame(distance = all_distances)

# Median
med_val <- median(dist_df$distance)

# Plot histogram with median line
ggplot(dist_df, aes(x = distance)) +
  geom_histogram(binwidth = 100, color = "black", fill = "gray70") +
  geom_vline(aes(xintercept = med_val), color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = med_val, y = Inf, label = paste("Median:", med_val),
           vjust = -0.5, hjust = 0.5, color = "red") +
  theme_minimal() +
  labs(title = "Pairwise same-strand paralog distances in House Finch",
       x = "Distance (bp)",
       y = "Count")

