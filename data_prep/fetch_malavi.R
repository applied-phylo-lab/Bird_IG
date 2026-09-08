# Build a per-species avian haemosporidian (malaria) diversity table from MalAvi,
# the public database of Plasmodium / Haemoproteus / Leucocytozoon cytochrome b
# lineages in birds (Bensch et al. 2009).
#
# MalAvi currently has no permanent web address; the database is bundled inside
# the malaviR R package, so we pull the package tarball and read the bundled .rds.
#
# Writes {INPUT_DIR}/species_pathogens_malavi.csv, one row per matched host species.
#
# IMPORTANT: MalAvi records only POSITIVE detections plus the number of birds
# screened. A species absent from MalAvi has not been sampled -- it is missing
# data, not a true zero. Absent species are omitted here rather than coded 0.

suppressMessages({library(data.table); library(utils)})

INPUT_DIR <- "/local/storage/kav67/clean_birds/"
TARBALL   <- "https://github.com/vincenzoaellis/malaviR/archive/refs/heads/master.tar.gz"
CACHE     <- file.path(INPUT_DIR, "malavi_db.rds")

if (!file.exists(CACHE)) {
  message("downloading MalAvi (bundled in malaviR) -> ", CACHE)
  tmp <- tempfile(fileext = ".tar.gz"); dir <- tempfile(); dir.create(dir)
  download.file(TARBALL, tmp, quiet = TRUE)
  untar(tmp, exdir = dir)
  rds <- list.files(dir, pattern = "^malavi_db_.*\\.rds$", recursive = TRUE, full.names = TRUE)
  stopifnot(length(rds) == 1)
  file.copy(rds, CACHE)
}

db <- readRDS(CACHE)
message("MalAvi release: ", as.character(db$version))

# Wild birds only -- captive birds do not reflect natural pathogen exposure.
h <- as.data.table(db$hosts_and_sites)[HOST_ENVIRONMENT == "Wild" & PARASITE_GENUS != "N/A"]
h[, NUMBER_TESTED := suppressWarnings(as.numeric(NUMBER_TESTED))]

per_host <- h[, .(
  n_lineages      = uniqueN(LINEAGE_NAME),
  n_para_genera   = uniqueN(PARASITE_GENUS),
  n_plasmodium    = uniqueN(LINEAGE_NAME[PARASITE_GENUS == "Plasmodium"]),
  n_haemoproteus  = uniqueN(LINEAGE_NAME[PARASITE_GENUS == "Haemoproteus"]),
  n_leucocytozoon = uniqueN(LINEAGE_NAME[PARASITE_GENUS == "Leucocytozoon"]),
  n_records       = .N,
  n_studies       = uniqueN(REFERENCE_NAME),          # sampling-effort proxy
  n_countries     = uniqueN(COUNTRY_NAME),
  total_tested    = sum(NUMBER_TESTED, na.rm = TRUE)  # sampling-effort proxy
), by = .(malavi_species = SPECIES_NAME)]

traits <- fread(file.path(INPUT_DIR, "species_traits_avonet.csv"))

# MalAvi follows its own taxonomy, so try our name and the AVONET name.
traits[, malavi_species := fifelse(LatinName  %in% per_host$malavi_species, LatinName,
                          fifelse(AvonetName %in% per_host$malavi_species, AvonetName, NA_character_))]

out <- merge(traits[!is.na(malavi_species)], per_host, by = "malavi_species")
setcolorder(out, c("Species", "LatinName", "malavi_species"))
fwrite(out, file.path(INPUT_DIR, "species_pathogens_malavi.csv"))

cat(sprintf("\nmatched %d / %d species to MalAvi\n", nrow(out), nrow(traits)))
cat("\nlineage counts by migration status:\n")
print(out[, .(n = .N,
              median_lineages = median(n_lineages),
              median_studies  = median(n_studies),
              median_tested   = median(total_tested)), by = MigrationLabel])
cat("\nspecies with NO MalAvi data (unsampled, excluded):", nrow(traits) - nrow(out), "\n")
