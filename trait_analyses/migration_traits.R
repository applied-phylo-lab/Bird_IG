# Test whether IG locus features differ between migratory and sedentary birds.
# Migration status from AVONET (see fetch_avonet_traits.py in this folder).
#
# Migration is strongly clustered by clade (e.g. Anseriformes mostly migratory,
# Struthioniformes all sedentary), so every test here is phylogenetically corrected.

library(data.table); library(dplyr); library(ape); library(phylolm)

INPUT_DIR <- "/local/storage/kav67/clean_birds/"

traits <- fread(file.path(INPUT_DIR, "species_traits_avonet.csv"))
igh    <- fread(file.path(INPUT_DIR, "IGH_VGP_table.tsv"))
tree   <- read.tree(file.path(INPUT_DIR, "vgp_birds.nwk"))

# One value per species: average V gene count over haplotypes.
numv <- igh[Locus == "IGH", .(NumV = mean(NumV, na.rm = TRUE)), by = .(Species, LatinName)]

# as.data.frame is required: phylolm indexes by rownames, and data.table would
# reinterpret dat[phy$tip.label, ] as a keyed join and silently return all-NA rows.
dat <- as.data.frame(merge(numv, traits[, -"LatinName"], by = "Species", all.x = TRUE)) %>%
  mutate(tip = gsub(" ", "_", LatinName),
         migratory = as.integer(Migration == 3),           # 3 = migratory
         moves     = as.integer(Migration >= 2),           # partial or migratory
         logMass   = log10(as.numeric(Mass)),
         HWI       = as.numeric(`Hand-Wing.Index`)) %>%
  filter(!is.na(Migration))

# Keep only species present in the tree.
dat <- dat[dat$tip %in% tree$tip.label, ]
cat("species with trait + tree data:", nrow(dat), "\n")
phy <- keep.tip(tree, dat$tip)
rownames(dat) <- dat$tip
dat <- dat[phy$tip.label, ]

cat("\n=== NumV by migration status (raw means) ===\n")
print(dat %>% group_by(MigrationLabel) %>%
        summarise(n = n(), mean_NumV = round(mean(NumV), 1), sd = round(sd(NumV), 1)))

run <- function(formula, label) {
  m <- phylolm(formula, data = dat, phy = phy, model = "lambda")
  cat("\n---", label, "---\n"); print(summary(m)$coefficients)
  invisible(m)
}

# Body mass is a standard covariate; HWI is a continuous dispersal proxy.
run(NumV ~ migratory,           "NumV ~ migratory (binary)")
run(NumV ~ migratory + logMass, "NumV ~ migratory + log10(mass)")
run(NumV ~ Migration,           "NumV ~ migration (ordinal 1-3)")
run(NumV ~ HWI,                 "NumV ~ hand-wing index")

# Same question for inversion structure, if the table is available.
inv_path <- file.path(INPUT_DIR, "inversion_stats.tsv")
if (file.exists(inv_path)) {
  inv <- fread(inv_path)[minlen == 250,
          .(num_inversions = mean(num_inversions, na.rm = TRUE),
            frac_genes_on_inv = mean(frac_genes_on_inv, na.rm = TRUE)), by = species]
  d2 <- as.data.frame(merge(dat, as.data.frame(inv), by.x = "Species", by.y = "species"))
  d2 <- d2[d2$tip %in% phy$tip.label, ]
  rownames(d2) <- d2$tip
  phy2 <- keep.tip(phy, d2$tip); d2 <- d2[phy2$tip.label, ]
  cat("\n\n########## INVERSIONS (n =", nrow(d2), "species) ##########\n")
  for (resp in c("num_inversions", "frac_genes_on_inv")) {
    m <- phylolm(as.formula(paste(resp, "~ migratory")), data = d2, phy = phy2, model = "lambda")
    cat("\n---", resp, "~ migratory ---\n"); print(summary(m)$coefficients)
  }
}
