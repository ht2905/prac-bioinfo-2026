#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

OUT_DIR="qc-group-04"
FASTQC_DIR="${OUT_DIR}/fastqc"
MULTIQC_DIR="${OUT_DIR}/multiqc"
REPORT_NAME="group-04_multiqc_report"
REPORT_TITLE="BT338IU Group 04 — Raw Read QC (PRJNA964723)"

mkdir -p "${MULTIQC_DIR}"

printf "\n=== Aggregating FastQC reports with MultiQC ===\n"

# MultiQC parses fastqc_data.txt inside each *_fastqc.zip archive (not raw FASTQs)
# so all 12 R1/R2 libraries can be compared on shared axes.
#   --module fastqc : limit auto-detection to FastQC outputs to avoid scanning other files
#   --filename      : pin clean output names instead of slugifying --title
#   --export        : render static PNG/SVG/PDF plots alongside the interactive HTML report
#   --force         : keep re-runs idempotent by overwriting instead of numbering outputs
mamba run -n ngs multiqc \
    --module fastqc \
    --title "${REPORT_TITLE}" \
    --filename "${REPORT_NAME}.html" \
    --outdir "${MULTIQC_DIR}" \
    --export \
    --force \
    "${FASTQC_DIR}"

printf "\nHTML report : %s/%s.html\n" "${MULTIQC_DIR}" "${REPORT_NAME}"
printf "Parsed data : %s/%s_data/\n" "${MULTIQC_DIR}" "${REPORT_NAME}"
printf "Static plots: %s/%s_plots/\n" "${MULTIQC_DIR}" "${REPORT_NAME}"
