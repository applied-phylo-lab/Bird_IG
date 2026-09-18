"""Stage 6 -- PatchWorkPlot dot plots, one per taxonomic order x locus.

make_config_strand.py writes the config PatchWorkPlot consumes (FASTA + strand
BED per haplotype, one row per haplotype's highest-NumV contig); PatchWorkPlot
then aligns every pair and draws the triangle plot.

PatchWorkPlot is not on $PATH -- config["patchworkplot"] points at the checkout.
Its pairwise_alignments/ output is what a future rewrite of shared_inversions.py
should read, so it is declared as a directory output here.
"""

rule dotplot_config:
    """config_strand_{locus}.csv for one order."""
    input:
        index=_index_path,
        beds=lambda wc: [f"{hap_dir(r.Order, r.Species, r.Haplotype)}/{r.Contig}_{r.Locus}_strand.bed"
                         for r in units[(units.Order == wc.order) &
                                        (units.Locus == wc.locus)].itertuples()],
    output:
        "{input_dir}/{order}/patchworkplot/config_strand_{locus}.csv",
    params:
        script=f"{REPO_DIR}/make_config_strand.py",
        outdir=lambda wc, output: os.path.dirname(output[0]),
    log:
        "{input_dir}/{order}/patchworkplot/logs/config_{locus}.log",
    resources:
        mem_mb=2000,
        runtime=20,
    shell:
        "python {params.script} -i {wildcards.input_dir} -s {input.index}"
        " -o {params.outdir} --order {wildcards.order}"
        " --locus {wildcards.locus} > {log} 2>&1"


rule dotplot:
    """PatchWorkPlot for one order x locus."""
    input:
        cfg="{input_dir}/{order}/patchworkplot/config_strand_{locus}.csv",
    output:
        pdf="{input_dir}/{order}/patchworkplot/{locus}/patchworkplot.pdf",
        stats="{input_dir}/{order}/patchworkplot/{locus}/alignment_stats.csv",
        alns=directory("{input_dir}/{order}/patchworkplot/{locus}/pairwise_alignments"),
    params:
        pwp=config["patchworkplot"],
        extra=config["patchworkplot_args"],
        outdir=lambda wc, output: os.path.dirname(output.pdf),
    threads: 4
    log:
        "{input_dir}/{order}/patchworkplot/logs/{locus}.log",
    resources:
        mem_mb=16000,
        runtime=720,
    shell:
        "python {params.pwp} -i {input.cfg} -o {params.outdir}"
        " {params.extra} > {log} 2>&1"
