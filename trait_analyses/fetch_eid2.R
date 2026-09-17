# Build a per-species pathogen richness table from EID2 (ENHanCEd Infectious
# Diseases Database, University of Liverpool), which records host-pathogen
# ("carrier"-"cargo") associations mined from NCBI sequence records and the
# published literature.
#
# EID2 spans viruses, bacteria, protozoa, helminths, fungi and ectoparasitic
# arthropods, so it covers far more of pathogen space than MalAvi's blood
# parasites -- and it carries its own effort proxies (papersNo, sequencesNo).
#
# Public 2022 dump: https://figshare.com/articles/dataset/22060922
# Writes {INPUT_DIR}/species_pathogens_eid2.csv, one row per matched host species.
#
# Same caveat as MalAvi: a species absent from EID2 is unstudied, not pathogen-free.

suppressMessages({library(data.table); library(utils)})

INPUT_DIR <- "/local/storage/kav67/clean_birds/"
EID2_URL  <- "https://ndownloader.figshare.com/files/39418972"
CACHE     <- file.path(INPUT_DIR, "eid2_organism_interactions_2022.csv")

if (!file.exists(CACHE)) {
  message("downloading EID2 (~50 MB) -> ", CACHE)
  download.file(EID2_URL, CACHE, quiet = TRUE)
}

e <- fread(CACHE)
b <- e[carrierclass == "aves"]

# EID2 lists some hosts at subspecies rank (e.g. "gallus gallus spadiceus");
# collapse everything to a binomial so it joins to our species table.
b[, host := tolower(sub("^([a-z]+ [a-z.]+).*$", "\\1", carrier))]
b <- b[!grepl(" sp\\.$", host)]          # drop genus-level "Genus sp." records

# uniqueN(cargo[...]) must be written inline: a helper defined outside the
# data.table call cannot see the `cargo` / `cargoclass` columns.
per_host <- b[, .(
  n_pathogens  = uniqueN(cargo),
  n_virus      = uniqueN(cargo[cargoclass == "virus"]),
  n_bacteria   = uniqueN(cargo[cargoclass == "bacteria"]),
  n_protozoa   = uniqueN(cargo[cargoclass == "protozoa"]),
  n_helminth   = uniqueN(cargo[cargoclass %in% c("helminth", "worms")]),
  n_fungi      = uniqueN(cargo[cargoclass == "fungi"]),
  n_arthropod  = uniqueN(cargo[cargoclass == "arthropods"]),  # ectoparasites: mites, lice, ticks
  n_records    = .N,
  n_papers     = sum(papersNo,    na.rm = TRUE),   # sampling-effort proxy
  n_sequences  = sum(sequencesNo, na.rm = TRUE)    # sampling-effort proxy
), by = host]

traits <- fread(file.path(INPUT_DIR, "species_traits_avonet.csv"))
traits[, host := fifelse(tolower(LatinName)  %in% per_host$host, tolower(LatinName),
               fifelse(tolower(AvonetName) %in% per_host$host, tolower(AvonetName), NA_character_))]

out <- merge(traits[!is.na(host)], per_host, by = "host")
setcolorder(out, c("Species", "LatinName", "host"))
fwrite(out, file.path(INPUT_DIR, "species_pathogens_eid2.csv"))

cat(sprintf("\nmatched %d / %d species to EID2\n", nrow(out), nrow(traits)))
cat("\nby migration status (medians):\n")
print(out[, .(n = .N, pathogens = median(n_pathogens), virus = median(n_virus),
              papers = median(n_papers), seqs = median(n_sequences)), by = MigrationLabel])
