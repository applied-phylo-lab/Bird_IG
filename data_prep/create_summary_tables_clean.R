library(data.table)
library(stringr)

# ==== USER INPUT ====
dir<-"/local/storage/kav67/clean_birds/"
#input_csv  <- paste0(dir,"db_creation_scripts/good_contigs.csv")
input_csv  <- paste0(dir,"ig_contig_list.csv")

output_all <- paste0(dir,"summary_features.csv")
output_igh <- paste0(dir,"IGH_filtered_table.tsv")
# ====================

# Read data
dt <- fread(input_csv)

# Remove leading "#" from Source
dt[, Source := sub("^#", "", Source)]

# Split Source into components
# Expected format: Order/Species/Haplotype
dt[, c("Order", "Species", "Haplotype") := tstrsplit(Source, "/", fixed=TRUE)]

# Clean the contig column
# Example format: "['OZ222767.1', 35]"
#dt[, Contig := str_extract(`[Contig Name, Number of Genes]`, "'([^']+)'")]
#dt[, Contig := gsub("'", "", Contig)]

dt[, NumV := as.integer(str_extract(`[Contig Name, Number of Genes]`, "\\d+(?=\\])"))]
dt[, NumV := `Number of Genes (no filtering)`]

# Select final columns
formatted <- dt[, .(Order, Species, Haplotype, Locus, Contig,NumV)]

# Write full formatted CSV
fwrite(formatted, output_all)

# Filter for IGH
formatted_igh <- formatted[Locus == "IGH"]
formatted_igh <- formatted_igh[NumV >2]

# Write IGH-only CSV
write_tsv(formatted_igh, output_igh)

cat("Done.\n")
cat("Wrote:\n")
cat(paste0(" - ", output_all, "\n"))
cat(paste0(" - ", output_igh, "\n"))