#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

THREADS=8
OUT_DIR="qc-group-04"
REF_DIR="${OUT_DIR}/ref"
BAM_DIR="${OUT_DIR}/align"
COUNTS_DIR="${OUT_DIR}/counts"
GTF_FILE="${REF_DIR}/genomic.gtf"
OUT_COUNTS="${COUNTS_DIR}/group-04_featurecounts.tsv"

# Strand-specificity: 0 = unstranded, 1 = forward-stranded, 2 = reverse-stranded.
# PRJNA964723 is reverse-stranded (dUTP library: -s 2 assigns ~7.20M vs. 0.19M for -s 1 on SRR24376036).
STRAND=2

mkdir -p "${COUNTS_DIR}"

printf "\n=== Counting paired-end fragments per gene with featureCounts ===\n"

# featureCounts processes all BAMs in one pass to build a single gene-by-sample matrix.
#   -p --countReadPairs : count each paired-end cDNA fragment once instead of each mate
#   -B -C               : require both mates mapped (-B) and exclude chimeric pairs (-C)
#   -t gene -g gene_id  : prokaryotic RefSeq GTFs annotate protein-coding loci as gene/CDS
#                         (only rRNA/tRNA have exon rows), so -t exon would miss coding genes
#   -s                  : apply strand-specificity setting defined above
mamba run -n ngs featureCounts \
    -T "${THREADS}" \
    -p --countReadPairs \
    -B \
    -C \
    -s "${STRAND}" \
    -t gene \
    -g gene_id \
    -a "${GTF_FILE}" \
    -o "${OUT_COUNTS}" \
    "${BAM_DIR}"/*.bam

printf "\nCount matrix : %s\n" "${OUT_COUNTS}"
printf "Summary log  : %s.summary\n" "${OUT_COUNTS}"
