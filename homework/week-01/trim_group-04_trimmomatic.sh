#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

SRR_FILE="group-04_srr-list.txt"
THREADS=8
OUT_DIR="qc-group-04"
FASTQ_DIR="${OUT_DIR}/fastq"
TRIM_DIR="${OUT_DIR}/trimmed"
LOG_DIR="${OUT_DIR}/logs"

# Resolve Illumina TruSeq3 paired-end adapters bundled with Trimmomatic in ngs
NGS_PREFIX=$(mamba run -n ngs printenv CONDA_PREFIX)
ADAPTERS="${NGS_PREFIX}/share/trimmomatic/adapters/TruSeq3-PE-2.fa"

mkdir -p "${TRIM_DIR}" "${LOG_DIR}"
mapfile -t SRR_IDS < "${SRR_FILE}"

total=${#SRR_IDS[@]}
i=0

for srr_id in "${SRR_IDS[@]}"; do
    printf "\n=== [%d/%d] Trimming %s ===\n" "$((++i))" "${total}" "${srr_id}"

    # Trimmomatic applies steps strictly in argument order:
    #   -Xmx4g              : raise JVM heap above Bioconda's 1 GB default for 8-thread NovaSeq runs
    #   1. ILLUMINACLIP     : remove Illumina adapters first (palindrome mode for PE read-through;
    #                         2 seed mismatches, 30 palindrome / 10 simple score, min 2 bp, keep both reads)
    #   2. LEADING/TRAILING : trim low-quality bases (< Q3) from 5' and 3' ends
    #   3. SLIDINGWINDOW    : cut when the 4-bp sliding-window mean drops below Q15
    #   4. MINLEN           : discard reads shortened below 36 bp (too short for unique mapping)
    mamba run -n ngs trimmomatic -Xmx4g PE \
        -threads "${THREADS}" \
        -phred33 \
        -summary "${LOG_DIR}/${srr_id}_trimmomatic.summary" \
        "${FASTQ_DIR}/${srr_id}_1.fastq" \
        "${FASTQ_DIR}/${srr_id}_2.fastq" \
        "${TRIM_DIR}/${srr_id}_1_paired.fastq.gz" \
        "${TRIM_DIR}/${srr_id}_1_unpaired.fastq.gz" \
        "${TRIM_DIR}/${srr_id}_2_paired.fastq.gz" \
        "${TRIM_DIR}/${srr_id}_2_unpaired.fastq.gz" \
        "ILLUMINACLIP:${ADAPTERS}:2:30:10:2:true" \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36
done

printf "\nTrimmed reads : %s/\n" "${TRIM_DIR}"
printf "Summary logs  : %s/\n" "${LOG_DIR}"
