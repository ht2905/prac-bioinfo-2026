# Workspace Instructions — BT338IU Practice in Bioinformatics

## Context

This repository houses coursework, automated bioinformatics pipelines, and laboratory reports for **BT338IU – Practice in Bioinformatics** (International University – VNUHCM, School of Biotechnology; Instructor: Dr. Tran Minh Nhựt).

## Technical Stack & Environments

- **Python CLI orchestration (`uv`):** Run project scripts via `uv run python3 src/<script>.py`.
- **NGS & Phylogeny binaries (`mamba` `ngs`):** `fastqc`, `fastp`, `trimmomatic`, `hisat2`, `star`, `salmon`, `samtools`, `bcftools`, `subread` (`featureCounts`), `htseq`, `multiqc`, `sra-tools`, `mafft`, `fasttree`, `blast`.
- **R / Bioconductor DGE & Enrichment (`mamba` `ngs` + `renv`):** `DESeq2`, `edgeR`, `GEOquery`, `clusterProfiler`, `enrichplot`, `ggtree`, `ape`. Run via `mamba run -n ngs Rscript src/dge_analysis.R`.
- **Structural Bioinformatics (`mamba` `sbvs`):** `pymol-open-source`, `pdb2pqr`, `apbs`, `modeller` (`mod10.7`), `biopython`.
- **CLI Standards:** All tools in `src/` must follow `self-documenting-code` conventions (`--help` to `stdout`, errors to `stderr`, `-n/--dry-run` support).
