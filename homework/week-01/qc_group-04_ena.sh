#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

SRR_FILE="group-04_srr-list.txt"
THREADS=8
OUT_DIR="homework/week-01/qc-group-04"
FASTQ_DIR="${OUT_DIR}/fastq"
QC_DIR="${OUT_DIR}/fastqc"

mkdir -p "${FASTQ_DIR}" "${QC_DIR}"
mapfile -t SRR_IDS < "${SRR_FILE}"

total=${#SRR_IDS[@]} # Adding `#` counts elements in array
i=0

for srr_id in "${SRR_IDS[@]}"; do
    printf "\n=== [%d/%d] %s ===\n" "$((++i))" "${total}" "${srr_id}"

    # --summary-interval=0 keeps only the live progress bar instead of periodic log lines
    curl -fsS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${srr_id}&result=read_run&fields=fastq_ftp" \
        | awk -F'\t' 'NR==2 { n = split($2, urls, ";"); for (k = 1; k <= n; k++) print "https://" urls[k] }' \
        | aria2c -x 8 -s 8 -c --summary-interval=0 -d "${FASTQ_DIR}" -i -
done

printf "\n=== Running FastQC (%d samples) ===\n" "${total}"
mamba run -n ngs fastqc \
    --threads "${THREADS}" \
    --outdir "${QC_DIR}" \
    "${FASTQ_DIR}"/*.fastq.gz
