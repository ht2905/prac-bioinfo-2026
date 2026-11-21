#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

THREADS=8
OUT_DIR="qc-group-04"
TRIM_DIR="${OUT_DIR}/trimmed"
FASTQC_TRIM_DIR="${OUT_DIR}/fastqc-trimmed"
MULTIQC_TRIM_DIR="${OUT_DIR}/multiqc-trimmed"
REPORT_NAME="group-04_trimmed_multiqc_report"
REPORT_TITLE="BT338IU Group 04 — Post-Trimmomatic Read QC (PRJNA964723)"

mkdir -p "${FASTQC_TRIM_DIR}" "${MULTIQC_TRIM_DIR}"

printf "\n=== Running second-pass FastQC on surviving paired reads ===\n"

# Evaluate only *_paired.fastq.gz (the synchronized mates used for alignment)
# to verify that 3' adapter read-through and low-quality tails were removed.
mamba run -n ngs fastqc \
    --threads "${THREADS}" \
    --outdir "${FASTQC_TRIM_DIR}" \
    "${TRIM_DIR}"/*_paired.fastq.gz

printf "\n=== Aggregating post-trimming FastQC reports with MultiQC ===\n"

# Write to a separate directory from pre-trimming QC so before/after reports
# and exported plots can be compared side by side.
mamba run -n ngs multiqc \
    --module fastqc \
    --title "${REPORT_TITLE}" \
    --filename "${REPORT_NAME}.html" \
    --outdir "${MULTIQC_TRIM_DIR}" \
    --export \
    --force \
    "${FASTQC_TRIM_DIR}"

printf "\nFastQC zips : %s/\n" "${FASTQC_TRIM_DIR}"
printf "HTML report : %s/%s.html\n" "${MULTIQC_TRIM_DIR}" "${REPORT_NAME}"
printf "Parsed data : %s/%s_data/\n" "${MULTIQC_TRIM_DIR}" "${REPORT_NAME}"
printf "Static plots: %s/%s_plots/\n" "${MULTIQC_TRIM_DIR}" "${REPORT_NAME}"
