# archive/RSS/

The obsolete half of the old `RSS/` folder. These produce RSS extraction outputs
that nothing in the current pipeline or in any manuscript figure reads.

| Script | What it did |
|--------|-------------|
| `extract_RSS.py` | Pulled the heptamer/nonamer flanking each V gene out of the locus FASTA and wrote a per-haplotype RSS table. Superseded by the RSS columns Daniel's pipeline now writes directly into `gene_list.csv`. |
| `extract_RSS_100bp.py` | Same idea, emitting 100 bp windows downstream of each gene as FASTA, for motif discovery with MEME. |
| `match_genes_RSS.py` | Joined extracted RSS calls back onto the gene table and wrote summary + pair tables. |
| `run_RSS_IgD_all.py` | Batch driver running the RSS search across every haplotype for IGH D genes. |
| `plot_meme_output.R` | Plotted a MEME motif result. |

`rss_correlation.R` is **not** here: it turned out to hold a manuscript figure
(genes with RSS vs total genes, beside the single- vs multiple-productive-RSS
positional densities), so it moved to `plots/manuscript/rss_correlation.R`. It
also supplies `p_combined_simple` to `rss_position_oriented.R`.

With it gone, the old `RSS/` folder is empty and has been removed.
