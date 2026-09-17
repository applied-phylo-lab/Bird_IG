"""Stage 5 -- per-haplotype V gene trees.

tree_building_pipeline.py ran these three steps itself, inside a
multiprocessing.Pool, skipping haplotypes whose output already existed. Split
into three rules the scheduler can see: ~800 IQ-TREE jobs (IGH + IGL) stop being
capped at one node, and a failure names the haplotype that failed.

Legacy file naming is preserved: IGH trees are {hap}.*, IGL trees are IGL_{hap}.*
-- distances.R and tree_distance_rss.R read those names. The `prefix` wildcard
carries it; prefix_to_locus_hap() recovers the locus.
"""

def gene_table_for_prefix(wc):
    locus, hap = prefix_to_locus_hap(wc.prefix)
    d = hap_dir(wc.order, wc.species, hap)
    clean = f"{d}/combined_genes_{locus}_clean.txt"
    return clean if os.path.exists(clean) else f"{d}/combined_genes_{locus}.txt"


rule tree_fasta:
    """V gene sequences as FASTA. Tip labels: {prefix}.{Pos}.{Contig}.{GeneType}.{Productive}.{Strand}"""
    input:
        genes=gene_table_for_prefix,
    output:
        "{input_dir}/{order}/{species}/{hap}/tree/{prefix}.fasta",
    params:
        script=f"{REPO_DIR}/tree_analyses/createFastaFromCSV.py",
        outdir=lambda wc, output: os.path.dirname(output[0]),
    log:
        "{input_dir}/{order}/{species}/{hap}/logs/{prefix}_fasta.log",
    resources:
        mem_mb=2000,
        runtime=20,
    shell:
        # createFastaFromCSV.py only writes a file when V genes are present, so
        # make an empty one rather than failing the whole DAG on a V-less locus.
        "python {params.script} {input.genes} {params.outdir} {wildcards.prefix}"
        " > {log} 2>&1; touch {output}"


rule tree_align:
    input:
        fasta="{input_dir}/{order}/{species}/{hap}/tree/{prefix}.fasta",
    output:
        aln="{input_dir}/{order}/{species}/{hap}/tree/{prefix}_aligned.fasta",
    threads: 4
    log:
        "{input_dir}/{order}/{species}/{hap}/logs/{prefix}_clustalo.log",
    resources:
        mem_mb=8000,
        runtime=240,
    shell:
        "clustalo -i {input.fasta} -o {output.aln} -t DNA --force"
        " --threads {threads} > {log} 2>&1"


rule tree_iqtree:
    """IQ-TREE writes several files sharing the prefix; the treefile is the one
    downstream analyses read. --redo because Snakemake, not IQ-TREE, decides
    when a job re-runs."""
    input:
        aln="{input_dir}/{order}/{species}/{hap}/tree/{prefix}_aligned.fasta",
    output:
        treefile="{input_dir}/{order}/{species}/{hap}/tree/{prefix}_tree.treefile",
        iqtree="{input_dir}/{order}/{species}/{hap}/tree/{prefix}_tree.iqtree",
    params:
        prefix=lambda wc, output: output.treefile[: -len(".treefile")],
    threads: 4
    log:
        "{input_dir}/{order}/{species}/{hap}/logs/{prefix}_iqtree.log",
    resources:
        mem_mb=8000,
        runtime=720,
    shell:
        "iqtree2 -s {input.aln} --prefix {params.prefix} --redo"
        " -T {threads} > {log} 2>&1"
