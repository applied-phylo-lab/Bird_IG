# Test the pathogen-exposure hypothesis against both IG locus responses.
#
#   Q1  does migration predict pathogen richness?
#   Q2  does pathogen richness predict IGH V gene count?
#   Q3  does pathogen richness predict inversion structure?
#
# Two independent pathogen datasets:
#   MalAvi  - haemosporidian (avian malaria) cyt-b lineages, vector-borne only
#   EID2    - viruses, bacteria, protozoa, helminths, fungi, ectoparasites
#
# Every model is phylogenetically corrected (phylolm, lambda). Two confounds are
# controlled throughout:
#   - sampling effort, which dominates raw pathogen richness in both databases
#   - locus assembly length, which mechanically drives how many inversions are callable

suppressMessages({library(data.table); library(ape); library(phylolm)})

INPUT_DIR <- "/local/storage/kav67/clean_birds/"
tree <- read.tree(file.path(INPUT_DIR, "vgp_birds.nwk"))

igh  <- fread(file.path(INPUT_DIR, "IGH_VGP_table.tsv"))
numv <- igh[Locus == "IGH", .(NumV = mean(NumV, na.rm = TRUE)), by = Species]

# inversions_stats.tsv is one row per CONTIG. Sum to haplotype level first so
# counts and lengths aggregate the same way, then average over haplotypes.
inv_hap <- fread(file.path(INPUT_DIR, "inversions_stats.tsv"))[minlen == 250,
             .(num_inversions = sum(num_inversions, na.rm = TRUE),
               genes_on_inv   = sum(genes_on_inv,   na.rm = TRUE),
               total_genes    = sum(total_genes,    na.rm = TRUE),
               locus_len      = sum(total_seq_length, na.rm = TRUE)), by = .(species, haplotype)]
inv <- inv_hap[, .(num_inversions    = mean(num_inversions),
                   frac_genes_on_inv = sum(genes_on_inv) / sum(total_genes),
                   locus_len         = mean(locus_len)), by = species]

ig <- merge(numv, inv, by.x = "Species", by.y = "species", all = TRUE)

run <- function(f, label, dat, phy) {
  m <- try(phylolm(f, data = dat, phy = phy, model = "lambda"), silent = TRUE)
  if (inherits(m, "try-error")) { cat("\n---", label, "--- FAILED\n"); return(invisible(NULL)) }
  co <- round(summary(m)$coefficients, 4)
  cat("\n---", label, "---\n"); print(co)
  cat("    lambda =", round(m$optpar, 3), " n =", nrow(dat), "\n")
}

analyse <- function(path_file, richness_col, effort_col, source_name) {
  cat("\n\n", strrep("#", 68), "\n### ", source_name,
      "   richness = ", richness_col, " | effort = ", effort_col, "\n",
      strrep("#", 68), "\n", sep = "")

  p <- fread(file.path(INPUT_DIR, path_file))
  d <- merge(p, ig, by = "Species")
  d[, `:=`(tip        = gsub(" ", "_", LatinName),
           migratory  = as.integer(Migration == 3),
           log_rich   = log10(get(richness_col) + 1),
           log_effort = log10(get(effort_col) + 1),
           log_locus  = log10(locus_len),
           logMass    = log10(as.numeric(Mass)))]
  d <- as.data.frame(d[tip %in% tree$tip.label & is.finite(log_rich)])
  phy <- keep.tip(tree, d$tip); rownames(d) <- d$tip; d <- d[phy$tip.label, ]

  cat("\nn =", nrow(d), " (sedentary", sum(d$Migration == 1),
      "/ partial", sum(d$Migration == 2), "/ migratory", sum(d$Migration == 3), ")\n")
  cat("effort confound: OLS R2 of log(richness) ~ log(effort) =",
      round(summary(lm(log_rich ~ log_effort, d))$r.squared, 3), "\n")

  cat("\n## Q1: migration -> pathogen richness")
  run(log_rich ~ migratory,              "richness ~ migratory (uncontrolled)", d, phy)
  run(log_rich ~ migratory + log_effort, "richness ~ migratory + effort",       d, phy)
  run(log_rich ~ Migration + log_effort, "richness ~ migration ordinal + effort", d, phy)

  cat("\n## Q2: pathogen richness -> V gene count")
  dv <- d[is.finite(d$NumV), ]; pv <- keep.tip(phy, dv$tip); dv <- dv[pv$tip.label, ]
  run(NumV ~ log_rich,                        "NumV ~ richness",                  dv, pv)
  run(NumV ~ log_rich + log_effort,           "NumV ~ richness + effort",         dv, pv)
  run(NumV ~ log_rich + log_effort + logMass, "NumV ~ richness + effort + mass",  dv, pv)

  cat("\n## Q3: pathogen richness -> inversion structure")
  di <- d[is.finite(d$num_inversions) & is.finite(d$log_locus), ]
  pi <- keep.tip(phy, di$tip); di <- di[pi$tip.label, ]
  cat("\n(n =", nrow(di), "species with inversion data)\n")
  run(num_inversions ~ log_rich,                          "inversions ~ richness",              di, pi)
  run(num_inversions ~ log_rich + log_locus,              "inversions ~ richness + locus length", di, pi)
  run(num_inversions ~ log_rich + log_locus + log_effort, "inversions ~ richness + locus + effort", di, pi)
  run(frac_genes_on_inv ~ log_rich + log_locus,           "frac genes on inv ~ richness + locus", di, pi)
  run(num_inversions ~ migratory + log_locus,             "inversions ~ migratory + locus length", di, pi)
}

analyse("species_pathogens_eid2.csv",   "n_pathogens", "n_papers",  "EID2 (all pathogens)")
analyse("species_pathogens_eid2.csv",   "n_virus",     "n_papers",  "EID2 (viruses only)")
analyse("species_pathogens_malavi.csv", "n_lineages",  "n_studies", "MalAvi (haemosporidians)")
