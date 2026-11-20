#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

SRR_FILE="group-04_srr-list.txt"
THREADS=8
OUT_DIR="homework/week-01/qc-group-04"
SRA_DIR="${OUT_DIR}/sra"
FASTQ_DIR="${OUT_DIR}/fastq"
QC_DIR="${OUT_DIR}/fastqc"
TMP_DIR="${OUT_DIR}/.tmp"

mkdir -p "${SRA_DIR}" "${FASTQ_DIR}" "${QC_DIR}" "${TMP_DIR}"
trap 'rm -rf "${TMP_DIR}"' EXIT
mapfile -t SRR_IDS < "${SRR_FILE}"

total=${#SRR_IDS[@]} # Adding `#` counts elements in array
i=0

for srr_id in "${SRR_IDS[@]}"; do
    printf "\n=== [%d/%d] %s ===\n" "$((++i))" "${total}" "${srr_id}"

    sra_url=$(mamba run -n ngs srapath "${srr_id}")
    # --summary-interval=0 keeps only the live progress bar instead of periodic log lines
    aria2c -x 8 -s 8 -c --summary-interval=0 -d "${SRA_DIR}" -o "${srr_id}.sra" "${sra_url}"

    mamba run -n ngs fasterq-dump \
        --progress \
        --threads "${THREADS}" \
        --temp "${TMP_DIR}" \
        --outdir "${FASTQ_DIR}" \
        "${SRA_DIR}/${srr_id}.sra"

    rm -f "${SRA_DIR}/${srr_id}.sra"
done

printf "\n=== Running FastQC (%d samples) ===\n" "${total}"
mamba run -n ngs fastqc \
    --threads "${THREADS}" \
    --outdir "${QC_DIR}" \
    "${FASTQ_DIR}"/*.fastq
