"""Stage 3 -- inversion summaries, D genes on inversions, hairpin/palindromes.

These three are gather steps: each reads every IGH self-alignment and writes one
table. They are cheap next to LASTZ and IQ-TREE, so they stay whole-table rather
than being sharded per haplotype.
"""

IGH_ALN = unit_files(inv_units, ".tsv")
IGH_BED = unit_files(inv_units, ".bed")


rule summarize_inversions:
    """inversion_stats.tsv (one row per haplotype x contig) and
    inversion_details.tsv (one row per inversion)."""
    input:
        aln=IGH_ALN,
        bed=IGH_BED,
        index=_index_path,
    output:
        stats=f"{INPUT_DIR}/inversion_stats.tsv",
        details=f"{INPUT_DIR}/inversion_details.tsv",
    params:
        script=f"{REPO_DIR}/inversions/summarize_inversions.py",
        min_len=MIN_INV_LEN,
    threads: 8
    log:
        f"{INPUT_DIR}/logs/summarize_inversions.log",
    resources:
        mem_mb=16000,
        runtime=120,
    shell:
        "python {params.script} -i {INPUT_DIR} -s {input.index}"
        " --locus IGH --min_len {params.min_len}"
        " -o {output.stats} -d {output.details} -c {threads} > {log} 2>&1"


# D gene calls come from per-haplotype IGHD.csv files, which are regenerated
# upstream by Daniel's d_gene_procedure scripts. Declaring them means a refreshed
# D gene annotation invalidates D_inversions.tsv instead of silently leaving a
# stale table behind (which is exactly what happened between May and August).
IGHD_FILES = sorted({
    f"{hap_dir(r.Order, r.Species, r.Haplotype)}/IGHD.csv"
    for r in inv_units.itertuples()
    if os.path.exists(f"{hap_dir(r.Order, r.Species, r.Haplotype)}/IGHD.csv")
})


rule d_genes_on_inversions:
    """Fraction of D genes falling inside an inverted region, per haplotype."""
    input:
        aln=IGH_ALN,
        ighd=IGHD_FILES,
        index=_index_path,
    output:
        f"{INPUT_DIR}/D_inversions.tsv",
    params:
        script=f"{REPO_DIR}/inversions/d_genes_on_inversions.py",
        min_len=MIN_INV_LEN,
    threads: 8
    log:
        f"{INPUT_DIR}/logs/d_genes_on_inversions.log",
    resources:
        mem_mb=8000,
        runtime=120,
    shell:
        "python {params.script} -i {INPUT_DIR} -s {input.index}"
        " -o {output} -c {threads} --min_inv_len {params.min_len} > {log} 2>&1"


rule hairpin:
    """Hairpin signal at the centre of each diagonal inversion.

    Depends on the text alignments from rule self_align_text, so hairpin.py finds
    them cached and does no LASTZ work of its own. --seed keeps the random
    control windows reproducible.
    """
    input:
        text=TEXT_ALN,
        index=_index_path,
    output:
        f"{INPUT_DIR}/palindromes.tsv",
    params:
        script=f"{REPO_DIR}/inversions/hairpin.py",
        lastz=config["lastz"],
        seed=config["hairpin_seed"],
    threads: 8
    log:
        f"{INPUT_DIR}/logs/hairpin.log",
    resources:
        mem_mb=16000,
        runtime=240,
    shell:
        "python {params.script} -i {INPUT_DIR} -s {input.index}"
        " -o {output} -c {threads} --lastz {params.lastz}"
        " --seed {params.seed} > {log} 2>&1"
