library(dplyr)
library(readr)

tsv <- "/local/storage/kav67/within_species/Songbirds/inversion_analysis/housefinch_inversions_by_State.tsv"
#tsv <- "/local/storage/kav67/within_species/Songbirds/inversion_analysis/housefinch_inversions_by_SubGroup.tsv"

# plot_inversion_dotplot() is defined in within_species_inversions/dotplots.R;
# this script relied on that file having been sourced first in the session.
source("/home/kav67/Bird_IG/within_species_inversions/dotplots.R")

df <- read_tsv(tsv, show_col_types = FALSE)
colnames(df)[1]<-"Species"

support_df <- df %>%
  group_by(Species, Inversion, RefStart, RefEnd) %>%
  summarise(
    n_haplotypes = n(),
    n_support = sum(Present),
    support_frac = n_support / n_haplotypes,
    # RefStart/RefEnd are grouping keys, so inside summarise() they are the
    # group's whole vector, not a scalar; newer dplyr rejects a non-size-1 result.
    inv_len = first(RefEnd) - first(RefStart),
    .groups = "drop"
  )
support_df<-unique(support_df)



species_list <- unique(support_df$Species)

plot_east<-plot_inversion_dotplot(support_df,"house_finches_E","East")
plot_west<-plot_inversion_dotplot(support_df,"house_finches_W","West")
plot_west+plot_east

plot_AL<-plot_inversion_dotplot(support_df,"AL","AL")
plot_AZ<-plot_inversion_dotplot(support_df,"AZ","AZ")
plot_CA<-plot_inversion_dotplot(support_df,"CA","CA")
plot_MA<-plot_inversion_dotplot(support_df,"MA","MA")
plot_NM<-plot_inversion_dotplot(support_df,"NM","NM")
plot_NY<-plot_inversion_dotplot(support_df,"NY","NY")
plot_OH<-plot_inversion_dotplot(support_df,"OH","OH")
plot_WA<-plot_inversion_dotplot(support_df,"WA","WA")


(plot_WA+plot_NY)/
  (plot_CA+plot_MA)/
  (plot_AZ+plot_OH)/
  (plot_NM+plot_AL)







