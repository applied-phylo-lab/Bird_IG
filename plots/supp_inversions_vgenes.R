# Supplementary panel: how far the inversion -> V gene relationship holds up.
#
#   A  the relationship on log-log axes (one species, Common Moorhen, is an
#      extreme outlier on inversions, locus length and V gene count alike)
#   B  the inversion slope with and without a locus-length covariate, across
#      three subsets -- the raw slope is robust, the locus-controlled one is not
#
# Writes figures/supp_inversions_vgenes.png

suppressMessages({library(data.table); library(ape); library(phylolm)
                  library(ggplot2); library(patchwork)})

INPUT_DIR <- "/local/storage/kav67/clean_birds/"
OUT       <- "/home/kav67/Bird_IG/figures/supp_inversions_vgenes.png"

tree <- read.tree(file.path(INPUT_DIR, "vgp_birds.nwk"))
vgp  <- unique(fread(file.path(INPUT_DIR, "IGH_VGP_table.tsv"))[, .(Species, LatinName)])

# inversion_stats.tsv is one row per CONTIG -- sum to haplotype level so that
# counts and sequence lengths are aggregated the same way.
inv <- fread(file.path(INPUT_DIR, "inversion_stats.tsv"))[minlen == 250,
        .(num_inversions = sum(num_inversions), locus_len = sum(total_seq_length)),
        by = .(species, haplotype)]
sf <- fread(file.path(INPUT_DIR, "summary_features.csv"))[Locus == "IGH",
        .(NumV = sum(NumV), n_contigs = .N), by = .(Species, Haplotype)]

sp <- merge(sf, inv, by.x = c("Species","Haplotype"), by.y = c("species","haplotype"))[
        , .(NumV = mean(NumV), num_inversions = mean(num_inversions),
            locus_len = mean(locus_len), n_contigs = mean(n_contigs)), by = Species]
sp <- merge(sp, vgp, by = "Species")
sp[, `:=`(tip = gsub(" ", "_", LatinName), log_locus = log10(locus_len))]

d <- as.data.frame(sp[tip %in% tree$tip.label & NumV > 0 & num_inversions > 0])
phy <- keep.tip(tree, d$tip); rownames(d) <- d$tip; d <- d[phy$tip.label, ]
d$assembly <- ifelse(d$n_contigs == 1, "Single contig", "Multiple contigs")

# --- slope under each subset x model ----------------------------------------
subsets <- list("All species"        = d,
                "Excl. Common Moorhen" = d[d$Species != "Common_Moorhen", ],
                "Excl. top 2"        = d[order(-d$num_inversions), ][-(1:2), ])

co <- rbindlist(lapply(names(subsets), function(sn) {
  dat <- subsets[[sn]]; p <- keep.tip(phy, dat$tip); dat <- dat[p$tip.label, ]
  rbindlist(lapply(c("Inversions alone", "+ locus length"), function(mn) {
    f <- if (mn == "Inversions alone") NumV ~ num_inversions else NumV ~ num_inversions + log_locus
    m <- phylolm(f, dat, phy = p, model = "lambda"); cf <- summary(m)$coefficients
    tc <- qt(0.975, m$n - length(coef(m)))
    data.table(subset = sn, model = mn, est = cf["num_inversions","Estimate"],
               lo = cf["num_inversions","Estimate"] - tc*cf["num_inversions","StdErr"],
               hi = cf["num_inversions","Estimate"] + tc*cf["num_inversions","StdErr"],
               p  = cf["num_inversions","p.value"])
  }))
}))
co[, `:=`(subset = factor(subset, levels = rev(names(subsets))),
          model  = factor(model, levels = c("Inversions alone", "+ locus length")))]

BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; INK <- "#1a1a19"; MUTED <- "#6b6a63"
base <- theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#e8e8e4", linewidth = 0.3),
        axis.title = element_text(colour = INK), axis.text = element_text(colour = MUTED),
        plot.title = element_text(colour = INK, face = "bold", size = 10),
        plot.subtitle = element_text(colour = MUTED, size = 8),
        plot.tag = element_text(colour = INK, face = "bold", size = 11),
        legend.key.size = unit(9, "pt"), legend.text = element_text(size = 7.5))

# the fitted model is linear in raw units, so it renders as a curve on log axes
cf1 <- coef(phylolm(NumV ~ num_inversions, d, phy = phy, model = "lambda"))
fitline <- data.frame(x = exp(seq(log(min(d$num_inversions)), log(max(d$num_inversions)), length = 200)))
fitline$y <- cf1[1] + cf1[2]*fitline$x

pA <- ggplot(d, aes(num_inversions, NumV)) +
  geom_line(data = fitline, aes(x, y), colour = INK, linewidth = 0.7) +
  geom_point(aes(fill = assembly), shape = 21, size = 2.1, stroke = 0.35, colour = "white", alpha = 0.9) +
  annotate("text", x = max(d$num_inversions), y = max(d$NumV)*1.18, label = "Common Moorhen",
           hjust = 1, size = 2.5, colour = MUTED) +
  scale_fill_manual(values = c("Single contig" = BLUE, "Multiple contigs" = ORANGE), name = NULL) +
  scale_x_log10() + scale_y_log10() +
  labs(tag = "A", x = "Inversions per locus (log)", y = "IGH V genes (log)",
       title = "More inversions, more V genes",
       subtitle = "Line: phylogenetic regression, fitted in linear units") +
  base + theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
               legend.background = element_rect(fill = "white", colour = NA))

pB <- ggplot(co, aes(est, subset, colour = model)) +
  geom_vline(xintercept = 0, colour = MUTED, linewidth = 0.3, linetype = "22") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0,
                linewidth = 0.7, position = position_dodge(width = 0.55)) +
  geom_point(size = 2.4, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = c("Inversions alone" = INK, "+ locus length" = BLUE), name = NULL) +
  labs(tag = "B", x = "V genes per additional inversion (95% CI)", y = NULL,
       title = "Only the uncontrolled slope is robust",
       subtitle = "With locus length in, the effect rests on one species") +
  base + theme(panel.grid.major.y = element_blank(),
               legend.position = c(0.98, 0.03), legend.justification = c(1, 0),
               legend.background = element_rect(fill = "white", colour = NA))

ggsave(OUT, pA + pB + plot_layout(widths = c(1, 1.05)),
       width = 190, height = 80, units = "mm", dpi = 300, bg = "white")
cat("wrote", OUT, "\n")
print(co[, .(subset, model, est = round(est,4), lo = round(lo,4), hi = round(hi,4), p = signif(p,2))])
