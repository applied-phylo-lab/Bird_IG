"""Stage 2 -- LASTZ self-alignment and gene-position BED files."""

LOCUS_FASTA = ("{input_dir}/{order}/{species}/{hap}/refined_ig_loci/igloci_fasta/"
               "{locus}_{contig}_{numv}Vs.fasta")


def find_locus_fasta(order, species, hap, locus, contig):
    """The extracted locus FASTA. NumV is part of the filename, so look it up."""
    row = units[(units.Order == order) & (units.Species == species) &
                (units.Haplotype == hap) & (units.Locus == locus) &
                (units.Contig == contig)]
    if row.empty:
        raise WorkflowError(f"No index row for {order}/{species}/{hap} {locus} {contig}")
    return LOCUS_FASTA.format(input_dir=INPUT_DIR, order=order, species=species,
                              hap=hap, locus=locus, contig=contig,
                              numv=int(row.NumV.iloc[0]))


def locus_fasta(wc):
    return find_locus_fasta(wc.order, wc.species, wc.hap, wc.locus, wc.contig)


def gene_table(wc):
    """Filtered gene table, falling back to the raw one where it is missing."""
    d = hap_dir(wc.order, wc.species, wc.hap)
    clean = f"{d}/combined_genes_{wc.locus}_clean.txt"
    return clean if os.path.exists(clean) else f"{d}/combined_genes_{wc.locus}.txt"


rule self_align_bed:
    """One LASTZ self-alignment plus both BED files for a single contig.

    self_alignment_bed.py skips work whose output already exists, so --force is
    passed: Snakemake has already decided this job needs to run.
    """
    input:
        fasta=locus_fasta,
        genes=gene_table,
        locus_summary=lambda wc: f"{hap_dir(wc.order, wc.species, wc.hap)}/refined_ig_loci/summary.csv",
        index=_index_path,
    output:
        aln="{input_dir}/{order}/{species}/{hap}/{contig}_{locus}.tsv",
        bed="{input_dir}/{order}/{species}/{hap}/{contig}_{locus}.bed",
        strand_bed="{input_dir}/{order}/{species}/{hap}/{contig}_{locus}_strand.bed",
    params:
        script=f"{REPO_DIR}/self_alignment_bed.py",
        lastz=config["lastz"],
    log:
        "{input_dir}/{order}/{species}/{hap}/logs/{contig}_{locus}_self_align.log",
    resources:
        mem_mb=4000,
        runtime=120,
    shell:
        "python {params.script}"
        " -i {wildcards.input_dir}"
        " -s {input.index}"
        " --locus {wildcards.locus}"
        " --order {wildcards.order}"
        " --species {wildcards.species}"
        " --haplotype {wildcards.hap}"
        " --contig {wildcards.contig}"
        " --lastz {params.lastz}"
        " --force -c 1"
        " > {log} 2>&1"


rule self_align_text:
    """Second self-alignment carrying the aligned sequences (text1/text2), which
    hairpin.py needs and the stage 2 alignment does not contain. Kept as its own
    rule so these run under the scheduler instead of inside hairpin.py."""
    input:
        fasta=lambda wc: find_locus_fasta(wc.order, wc.species, wc.hap,
                                          "IGH", wc.contig),
    output:
        "{input_dir}/{order}/{species}/{hap}/{contig}_IGH_text.tsv",
    params:
        lastz=config["lastz"],
        fmt=("--format=general:name1,strand1,start1,end1,length1,text1,"
             "name2,strand2,start2+,end2+,length2,text2,id%"),
    log:
        "{input_dir}/{order}/{species}/{hap}/logs/{contig}_IGH_text.log",
    resources:
        mem_mb=8000,
        runtime=240,
    shell:
        "{params.lastz} {input.fasta} {input.fasta} --step=20 --notransition"
        " '{params.fmt}' --output={output} 2> {log}"
