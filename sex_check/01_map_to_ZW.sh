#!/usr/bin/env bash
# Test 1: does the bAgePho2 assembly contain W-chromosome sequence?
#
# Aligns both bAgePho2 haplotypes to a conspecific FEMALE red-winged blackbird
# reference (bAgePho1.hap1, chromosome level, carries both Z and W).
# Females carry W, males do not -- so W-aligned assembly sequence is the readout.
#
# IMPORTANT: the reference is the WHOLE female genome, not just Z+W. Restricting
# it to Z+W forces every autosomal repeat to take its best hit on Z or W, which
# manufactures tens of Mb of fake "W" sequence. With all chromosomes present,
# repeats go home to their real chromosome and only genuine W maps to W.
set -euo pipefail
source "$(dirname "$0")/config.sh"

mkdir -p "$OUT"
cd "$OUT"

[[ -s ${REF_FEMALE}.fai ]] || "$SAMTOOLS" faidx "$REF_FEMALE"

# accession -> chromosome name table, used by 02_summarize_ZW.py
if [[ ! -s acc2chrom.tsv ]]; then
    curl -s "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/051/311/805/GCA_051311805.1_bAgePho1.hap1/GCA_051311805.1_bAgePho1.hap1_assembly_report.txt" \
      | awk -F'\t' '!/^#/ {print $5"\t"$3"\t"$9}' > acc2chrom.tsv
fi

for hap in pri alt; do
    case $hap in
        pri) QUERY=$PRI ;;
        alt) QUERY=$ALT ;;
    esac
    if [[ -s bAgePho2.${hap}_vs_female.paf ]]; then
        echo "$(date '+%F %T') - $hap PAF exists, skipping"
        continue
    fi
    echo "$(date '+%F %T') - mapping $hap -> whole female genome"
    # asm5: same species, different individual (~0.1-0.5% divergence)
    # --secondary=no: one best chromosome per assembly region, so Z-W gametologs
    #                 and the PAR are not counted twice
    "$MINIMAP2" -x asm5 -t "$THREADS" --secondary=no \
        "$REF_FEMALE" "$QUERY" > bAgePho2.${hap}_vs_female.paf.partial \
        2> minimap2.${hap}.log
    mv bAgePho2.${hap}_vs_female.paf.partial bAgePho2.${hap}_vs_female.paf
    echo "$(date '+%F %T') - $hap done: $(wc -l < bAgePho2.${hap}_vs_female.paf) alignments"
done

echo "$(date '+%F %T') - all mappings complete; now run 02_summarize_ZW.py"
