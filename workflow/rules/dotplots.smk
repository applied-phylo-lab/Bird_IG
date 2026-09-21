"""Stage 6 -- PatchWorkPlot dot plots, one per taxonomic order x locus.

make_config_strand.py writes the config PatchWorkPlot consumes (FASTA + strand
BED per haplotype, one row per haplotype's highest-NumV contig); PatchWorkPlot
then aligns every pair and draws the triangle plot.

PatchWorkPlot is not on $PATH -- config["patchworkplot"] points at the checkout.
Its pairwise_alignments/ output is what a future rewrite of shared_inversions.py
should read, so it is declared as a directory output here.
"""

# find_locus_fasta() and hap_dir() come from rules/align.smk and the Snakefile,
# both included ahead of this file.

def _fig1c_row(sample):
    """Index row for a fig1c sample: its highest-NumV IGH contig."""
    rows = units[(units.Species == sample["species"]) &
                 (units.Haplotype == sample["haplotype"]) &
                 (units.Locus == "IGH")]
    if rows.empty:
        raise WorkflowError(
            f"fig1c sample not in the index: "
            f"{sample['species']}/{sample['haplotype']} IGH")
    return rows.sort_values("NumV", ascending=False).iloc[0]


def fig1c_fasta(sample):
    r = _fig1c_row(sample)
    return find_locus_fasta(r.Order, r.Species, r.Haplotype, "IGH", r.Contig)


def fig1c_bed(sample):
    r = _fig1c_row(sample)
    return f"{hap_dir(r.Order, r.Species, r.Haplotype)}/{r.Contig}_IGH_strand.bed"


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
        script=f"{REPO_DIR}/dotplots/make_config_strand.py",
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


# --------------------------------------------------------------------------- #
# Figure 1C dot plot
#
# A hand-picked cross-section (house finch, tawny owl, golden plover, chicken,
# ostrich) rather than a taxonomic order, so it does not fit the per-order rules
# above and would otherwise require rebuilding all 16 orders to refresh one panel.
# --------------------------------------------------------------------------- #

rule fig1c_config:
    """Write config_fig1c.csv from config["fig1c_samples"].

    Resolves each sample's highest-NumV IGH contig through the index, so the FASTA
    and BED paths stay correct if the annotation is rebuilt.
    """
    input:
        index=_index_path,
        beds=lambda wc: [fig1c_bed(s) for s in config["fig1c_samples"]],
        fastas=lambda wc: [fig1c_fasta(s) for s in config["fig1c_samples"]],
    output:
        f"{INPUT_DIR}/patchworkplot/config_fig1c.csv",
    run:
        import csv as _csv
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        with open(output[0], "w", newline="") as fh:
            w = _csv.writer(fh)
            w.writerow(["SampleID", "Label", "Fasta", "Annotation", "Strand"])
            for s in config["fig1c_samples"]:
                w.writerow([s["haplotype"], s["label"],
                            fig1c_fasta(s), fig1c_bed(s), ""])
        print(f"wrote {output[0]} ({len(config['fig1c_samples'])} samples)")


rule fig1c_dotplot:
    input:
        cfg=f"{INPUT_DIR}/patchworkplot/config_fig1c.csv",
    output:
        pdf=f"{INPUT_DIR}/patchworkplot/plots_fig1c/patchworkplot.pdf",
        stats=f"{INPUT_DIR}/patchworkplot/plots_fig1c/alignment_stats.csv",
    params:
        pwp=config["patchworkplot"],
        extra=config["patchworkplot_args"],
        outdir=lambda wc, output: os.path.dirname(output.pdf),
    threads: 4
    log:
        f"{INPUT_DIR}/patchworkplot/logs/fig1c.log",
    resources:
        mem_mb=16000,
        runtime=720,
    shell:
        "python {params.pwp} -i {input.cfg} -o {params.outdir}"
        " {params.extra} > {log} 2>&1"
