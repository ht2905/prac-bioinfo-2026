# Workspace Instructions — BT338IU Practice in Bioinformatics

## Context

This repository houses coursework, laboratory exercises, and reports for **BT338IU – Practice in Bioinformatics** (International University – VNUHCM, School of Biotechnology; Instructor: Dr. Tran Minh Nhựt).

## Learning & Scripting Posture

- **Start with 1 example first:** Write flat, minimal, readable scripts (`15–40` lines) or direct commands for a single sample, gene, or structure inside `lecture-0X/` (or `.scratch/`) before scaling to full datasets.
- **Explain the biology & tools:** Explain the biological, statistical, and algorithmic purpose of key flags, thresholds, and intermediate formats (`FASTQ`, `BAM`, `.pir`, `.pqr`, `.dx`, `DESeq2` dispersion/shrunken $\log_2\text{FC}$, `DOPE`/`GA341`, `pLDDT`/`PAE`).
- **No unprompted CLI ceremony:** Avoid wrapping coursework scripts in `argparse`, `--dry-run`, or multi-preset frameworks unless explicitly asked (simple `@dataclass`/classes for structured records are fine when helpful). Treat existing scripts in `src/` as an optional reference archive.

## Technical Stack & Environments

- **Python (`uv`):** Run scripts via `uv run python3 <script>.py`.
- **NGS & Phylogeny (`mamba` `ngs`):** `fastqc`, `fastp`, `trimmomatic`, `hisat2`, `star`, `salmon`, `samtools`, `bcftools`, `subread` (`featureCounts`), `htseq`, `multiqc`, `sra-tools`, `mafft`, `fasttree`, `blast`.
- **R / Bioconductor DGE & Enrichment (`mamba` `ngs` + `renv`):** `DESeq2`, `edgeR`, `GEOquery`, `clusterProfiler`, `enrichplot`, `ggtree`, `ape`. Run via `mamba run -n ngs Rscript <script>.R`.
- **Structural Bioinformatics (`mamba` `sbvs`):** `pymol-open-source`, `pdb2pqr`, `apbs`, `modeller` (`mod10.7`), `biopython`.
