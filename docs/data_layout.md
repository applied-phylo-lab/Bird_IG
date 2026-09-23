# Data directory layout

Reference for what lives under `INPUT_DIR` (`/local/storage/kav67/clean_birds/`).
Nothing here is versioned. See the main README for how the pipeline uses it.

```
INPUT_DIR/
├── summary_features.csv          # Master table: one row per haplotype × locus
│                                 # Columns: Order, Species, Haplotype, Locus, Contig, NumV
├── IGH_filtered_table.tsv        # IGH-only filtered subset of summary_features
├── IGH_VGP_table.tsv             # IGH table with added LatinName column (for tree matching)
├── gene_list.csv                 # All V genes across all species with RSS annotations
│                                 # Columns: Source (Order/Species/Haplotype), GeneType,
│                                 #          Contig, Pos, Strand, Sequence, Productive,
│                                 #          Locus, Heptamer, Nonamer, ...
├── vgp_birds.nwk                 # VGP bird species phylogenetic tree (Newick)
├── inversion_stats.tsv          # Inversion summary stats, one row per haplotype x contig
├── inversion_details.tsv        # Per-inversion details (length, diagonal flag, etc.)
├── old_unused/D_inversions.tsv   # Retired -- see that folder's README
│
└── {Order}/                      # e.g. Doves/, Eagles/, Waterfowl/
    └── {Species}/                # e.g. Pink_Pigeon/
        └── {Haplotype}/          # e.g. bNesMay2_pri/  (pri = primary, alt = alternate)
            │
            ├── combined_genes_IGH_clean.txt   # Filtered IGH V genes (TSV)
            ├── combined_genes_IGL_clean.txt   # Filtered IGL V genes (TSV)
            │                                  # Columns: GeneType, Contig, Pos, Strand,
            │                                  #          Sequence, Productive, Locus
            ├── IGHD.csv                       # D genes for this haplotype
            │                                  # Columns: Source, GeneType, Contig, Pos,
            │                                  #          Strand, Sequence, Productive,
            │                                  #          Locus, ..., Location Relative to V-Cluster
            │
            ├── {Contig}_IGH.tsv              # LASTZ self-alignment output for IGH locus
            ├── {Contig}_IGH.bed              # Gene positions as BED (all genes, black)
            ├── {Contig}_IGH_strand.bed       # Gene positions as BED (- strand = grey)
            ├── {Contig}_IGL.tsv              # Same for IGL locus
            ├── {Contig}_IGL.bed
            ├── {Contig}_IGL_strand.bed
            │
            ├── refined_ig_loci/
            │   ├── summary.csv               # Locus boundaries: StartPos, EndPos per contig
            │   └── igloci_fasta/
            │       ├── IGH_{Contig}_{NumV}Vs.fasta   # Extracted IGH locus sequence
            │       └── IGL_{Contig}_{NumV}Vs.fasta   # Extracted IGL locus sequence
            │
            └── tree/
                ├── {Haplotype}.fasta                  # IGH V gene sequences (FASTA)
                ├── {Haplotype}_aligned.fasta          # IGH multiple sequence alignment
                ├── {Haplotype}_tree.treefile          # IGH IQ-TREE phylogeny
                ├── IGL_{Haplotype}.fasta              # IGL V gene sequences
                ├── IGL_{Haplotype}_aligned.fasta      # IGL multiple sequence alignment
                └── IGL_{Haplotype}_tree.treefile      # IGL IQ-TREE phylogeny
```
