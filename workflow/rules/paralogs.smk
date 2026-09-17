"""Stage 4 -- inversion paralog groups.

Two variants of the same analysis, kept separate because they answer different
questions: find_all_inversions_paralogs.py uses every inversion at >=250 bp,
find_inversion_paralogs.py only diagonal (palindromic) ones at >=1000 bp.
"""

rule paralogs_all:
    input:
        aln=IGH_ALN,
        bed=IGH_BED,
        index=_index_path,
    output:
        f"{INPUT_DIR}/IGH_paralogs_all.tsv",
    params:
        script=f"{REPO_DIR}/find_all_inversions_paralogs.py",
        min_len=MIN_INV_LEN,
    threads: 8
    log:
        f"{INPUT_DIR}/logs/paralogs_all.log",
    resources:
        mem_mb=8000,
        runtime=120,
    shell:
        "python {params.script} -i {INPUT_DIR} -s {input.index}"
        " --locus IGH --min_len {params.min_len}"
        " -o {output} -c {threads} > {log} 2>&1"


rule paralogs_diagonal:
    input:
        aln=IGH_ALN,
        bed=IGH_BED,
        index=_index_path,
    output:
        f"{INPUT_DIR}/IGH_paralogs_diag.tsv",
    params:
        script=f"{REPO_DIR}/find_inversion_paralogs.py",
    threads: 8
    log:
        f"{INPUT_DIR}/logs/paralogs_diag.log",
    resources:
        mem_mb=8000,
        runtime=120,
    shell:
        "python {params.script} -i {INPUT_DIR} -s {input.index}"
        " --locus IGH --min_len 1000"
        " -o {output} -c {threads} > {log} 2>&1"
