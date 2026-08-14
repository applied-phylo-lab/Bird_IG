#!/usr/bin/env bash
# Test 2 (independent confirmation, needs no reference): HiFi read depth.
#
# Map a subsample of this bird's own HiFi reads back to the PRIMARY assembly
# only. Autosomal contigs collect reads from both homologues (2n depth).
#   female (ZW): Z and W are single copy -> ~0.5x the autosomal depth
#   male   (ZZ): Z is diploid           -> ~1.0x, nothing sits at half depth
set -euo pipefail
source "$(dirname "$0")/config.sh"

FRAC=${FRAC:-0.05}          # fraction of HiFi reads to use (~5-6x is plenty)
mkdir -p "$OUT"
cd "$OUT"

if [[ ! -s depth_subsample.bam ]]; then
    echo "$(date '+%F %T') - subsampling ${FRAC} of HiFi reads and mapping to primary"
    "$SAMTOOLS" cat "${HIFI_BAMS[@]}" \
      | "$SAMTOOLS" view -@ 4 -b --subsample "$FRAC" --subsample-seed 42 - \
      | "$SAMTOOLS" fastq -@ 4 - \
      | "$MINIMAP2" -ax map-hifi -t "$THREADS" --secondary=no "$PRI" - \
      2> minimap2.depth.log \
      | "$SAMTOOLS" sort -@ 8 -m 4G -o depth_subsample.bam.partial -
    mv depth_subsample.bam.partial depth_subsample.bam
    "$SAMTOOLS" index -@ 8 depth_subsample.bam
fi

echo "$(date '+%F %T') - computing per-contig depth"
"$SAMTOOLS" coverage depth_subsample.bam > primary_contig_coverage.tsv
echo "$(date '+%F %T') - wrote $OUT/primary_contig_coverage.tsv; now run 04_depth_verdict.py"
