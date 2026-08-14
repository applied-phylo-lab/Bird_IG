#!/usr/bin/env bash
# Shared paths / settings for the bAgePho2 sex check.

ASM_DIR=/local/storage/kav67/RedwingedBlackbird_new
OUT=$ASM_DIR/sex_check

PRI=$ASM_DIR/bAgePho2.pri.fasta
ALT=$ASM_DIR/bAgePho2.alt.fasta

# Conspecific FEMALE red-winged blackbird, chromosome level (GCA_051311805.1).
# Contains SUPER_W (CM119744.1) + 2 unlocalized W scaffolds, and SUPER_Z (CM119746.1).
REF_FEMALE=$ASM_DIR/other_assemblies/GCA_051311805.1_bAgePho1.hap1_genomic.fna
CHR_Z=CM119746.1
CHR_W="CM119744.1 JBPJSG010000024.1 JBPJSG010000026.1"

# HiFi reads from the same individual (Revio, barcode bc2067)
HIFI_BAMS=("$ASM_DIR"/closeread/hifi_bam/bAgePho2/*.hifi_reads.bc2067.bam)

MINIMAP2=/home/kav67/miniconda3/envs/alignment_env/bin/minimap2
SAMTOOLS=/programs/samtools-1.20/bin/samtools

THREADS=48
