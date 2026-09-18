"""Stage 7 -- publication-facing annotation tables.

Both scripts already took INPUT_DIR positionally; the point of wiring them up is
declaring what they read, so the DAG can tell when a table has gone stale. That
was the gap that let D_inversions.tsv sit three months out of date behind a
refreshed IGHD.csv.

assign_haplotype_source rule queries the NCBI Datasets API, so it is a localrule
with retries and is never submitted to the cluster.
"""

localrules: haplotype_sources


rule haplotype_sources:
    """Assembly provenance per haplotype, from the NCBI Datasets API.

    Cached in ncbi_bioproject_cache.json, so a re-run only fetches accessions it
    has not seen. Network-bound and rate-limited: keep it off the cluster.
    """
    input:
        index=_index_path,
    output:
        f"{INPUT_DIR}/haplotype_sources.csv",
    params:
        script=f"{REPO_DIR}/annotation_tables/assign_haplotype_source.py",
    retries: 3
    log:
        f"{INPUT_DIR}/logs/haplotype_sources.log",
    resources:
        mem_mb=2000,
        runtime=120,
    shell:
        "python {params.script} {INPUT_DIR} > {log} 2>&1"


rule summary_tables:
    """main_table / species_summary / v_gene_table / d_gene_table, each also in a
    _published variant that drops the unpublished Darwin's finches."""
    input:
        index=_index_path,
        sources=f"{INPUT_DIR}/haplotype_sources.csv",
        inv=f"{INPUT_DIR}/inversion_stats.tsv",
        palindromes=f"{INPUT_DIR}/palindromes.tsv",
        vgp=f"{INPUT_DIR}/IGH_VGP_table.tsv",
        traits=f"{INPUT_DIR}/species_traits_avonet.csv",
        genes=f"{INPUT_DIR}/gene_list.csv",
        d_genes=f"{INPUT_DIR}/bird_d_genes.csv",
    output:
        SUMMARY_TABLES,
        readme=f"{INPUT_DIR}/summary_tables/README.md",
    params:
        script=f"{REPO_DIR}/annotation_tables/build_summary_tables.py",
    log:
        f"{INPUT_DIR}/logs/summary_tables.log",
    resources:
        mem_mb=16000,
        runtime=120,
    shell:
        "python {params.script} {INPUT_DIR} > {log} 2>&1"
