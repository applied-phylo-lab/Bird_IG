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

`RSS/rss_correlation.R` is **not** here: it is still needed. It builds
`p_combined_simple`, which `plots/manuscript/rss_position_oriented.R` uses for two
composite panels. That script writes no files of its own — it is analysis and
plotting only.
