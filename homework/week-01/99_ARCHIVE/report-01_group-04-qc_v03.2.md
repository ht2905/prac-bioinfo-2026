# Week 01 Lab Report: Raw Read Quality Control (Group 04)

**Course:** BT338IU – Practice in Bioinformatics  
**Instructor:** Dr. Tran Minh Nhựt  
**Group:** 04 | **BioProject:** `PRJNA964723`  
**Organism:** *Lactobacillus johnsonii* N6.2 (`GCF_000498675.1`)

## 1. Study Background and Dataset

Our assigned dataset comes from Beliakoff et al. (2024), *"Bile promotes Lactobacillus johnsonii N6.2 extracellular vesicle production with conserved immunomodulatory properties"* (*Scientific Reports*, 14:12272; DOI: `10.1038/s41598-024-62843-0`). In this study, the probiotic bacterium *L. johnsonii* N6.2 was grown in control broth (vdMRS) or with 0.2% bovine bile to mimic stress in the small intestine. Bile caused the bacteria to release ~100-fold more extracellular vesicles (EVs), and those EVs kept their immune-regulating activity in human pancreatic (βlox5) and macrophage (THP-1) cells.

To identify genes responding to bile, the authors extracted total RNA from three biological replicates per condition (`n = 3` control, `n = 3` bile) and sequenced paired-end libraries on an Illumina platform (deposited under BioProject `PRJNA964723` as 150-bp Illumina NovaSeq 6000 runs). Bile significantly changed the expression of 967 transcripts (`padj < 0.05`), upregulating cell-wall remodeling enzymes (`murV`, `sdp`, `pepD`) and fatty acid utilization/uptake genes (including the regulator `fakT` and DegV-domain genes `fakB4`, `fakB1`) while downregulating bile salt transporters (`cbsT1`, `cbsT2`) and bile salt hydrolase (`bsh`).

| Condition | Replicate | Run ID | Sample Alias | Read Pairs |
| :--- | :---: | :--- | :--- | ---: |
| Control (vdMRS) | 1 | `SRR24376036` | `Ljn_CNTRL_01` | 7,859,010 |
| Control (vdMRS) | 2 | `SRR24376035` | `Ljn_CNTRL_02` | 6,995,853 |
| Control (vdMRS) | 3 | `SRR24376033` | `Ljn_CNTRL_03` | 9,385,497 |
| 0.2% Bovine Bile | 1 | `SRR24376032` | `Ljn_BILE_01` | 7,157,456 |
| 0.2% Bovine Bile | 2 | `SRR24376034` | `Ljn_BILE_02` | 7,312,110 |
| 0.2% Bovine Bile | 3 | `SRR24376031` | `Ljn_BILE_03` | 9,266,307 |

## 2. FastQC Overview

We ran `FastQC` (v0.13.0; Andrews, 2010) on all 12 FASTQ files (forward `_1` and reverse `_2` reads for each run) and combined the diagnostics using `MultiQC` (v1.35; Ewels et al., 2016). Eleven of the 12 files passed 7 out of 10 FastQC checks, while `SRR24376031_2` passed 8 out of 10. Across the dataset, warnings or failures appeared only in `Per base sequence content` (PBSC), `Sequence Duplication Levels` (SDL), and `Overrepresented sequences` (OS), as summarized in the status heatmap and table below.

![FastQC module status heatmap across all 12 Group 04 FASTQ files](figures/fastqc_status-heatmap.png)

| File | Group | Total Sequences | Length | %GC | PBSC | SDL | OS |
| :--- | :--- | ---: | :---: | :---: | :---: | :---: | :---: |
| `SRR24376036_1` | Control (Fwd) | 7,859,010 | 150 bp | 37% | FAIL | FAIL | WARN |
| `SRR24376036_2` | Control (Rev) | 7,859,010 | 150 bp | 37% | WARN | FAIL | WARN |
| `SRR24376035_1` | Control (Fwd) | 6,995,853 | 150 bp | 37% | FAIL | FAIL | WARN |
| `SRR24376035_2` | Control (Rev) | 6,995,853 | 150 bp | 37% | WARN | FAIL | WARN |
| `SRR24376033_1` | Control (Fwd) | 9,385,497 | 150 bp | 37% | FAIL | FAIL | WARN |
| `SRR24376033_2` | Control (Rev) | 9,385,497 | 150 bp | 37% | WARN | FAIL | WARN |
| `SRR24376032_1` | Bile (Fwd) | 7,157,456 | 150 bp | 38% | FAIL | FAIL | WARN |
| `SRR24376032_2` | Bile (Rev) | 7,157,456 | 150 bp | 37% | WARN | FAIL | WARN |
| `SRR24376034_1` | Bile (Fwd) | 7,312,110 | 150 bp | 38% | FAIL | FAIL | WARN |
| `SRR24376034_2` | Bile (Rev) | 7,312,110 | 150 bp | 37% | WARN | FAIL | WARN |
| `SRR24376031_1` | Bile (Fwd) | 9,266,307 | 150 bp | 37% | WARN | FAIL | WARN |
| `SRR24376031_2` | Bile (Rev) | 9,266,307 | 150 bp | 36% | WARN | FAIL | PASS |

## 3. Understanding the QC Parameters

### Basic Statistics, Read Length, and N Content (`PASS`)

Each sample has between 7.0 and 9.4 million read pairs, and all reads are uniformly 150 bp long. Ambiguous bases (`N`, meaning the sequencer could not call A, T, G, or C) stay near 0% across all positions, so base calling ran cleanly.

### Per Base and Per Sequence Quality (`PASS`)

Quality scores use the Phred scale, where `Q20` means 99% base-call accuracy and `Q30` means 99.9% accuracy. Across all 12 files, `Per base sequence quality` mean scores stay inside the green high-quality zone (roughly Q33 to Q37) from base 1 to base 150, with only a small drop near the 3' end of reverse reads (`_2`). In `Per sequence quality scores`, reads peak sharply at Q36 across every library (with 95.9%–98.5% of reads achieving mean quality ≥ Q30), confirming the raw reads are accurate enough for alignment.

![Mean Phred quality scores across all 150 read positions](figures/fastqc_per-base-quality.png)

![Distribution of mean per-sequence Phred quality scores across all 12 files](figures/fastqc_per-sequence-quality.png)

### Per Base Sequence Content (`WARN` / `FAIL`)

This check expects A, T, G, and C proportions to stay balanced along the read. In our files, the first ~10 bases fluctuate noticeably before leveling off from base 15 to 150 (~30–33% A/T and ~16–21% G/C). This 5' wobble is normal in RNA-seq libraries because the random hexamer primers used during cDNA synthesis bind certain short sequences more easily than others. It does not mean the reads are low quality.

![Per base sequence content in SRR24376032_1 (representative forward-read profile) showing 5-prime random-hexamer priming bias](figures/fastqc_per-base-sequence-content-r1.png)

### Per Sequence GC Content (`PASS`)

Every library forms a single smooth peak centered around 37%–39% GC (overall `%GC` 36%–38%, distribution mode 38%–39%), slightly above the *L. johnsonii* N6.2 reference genome (34.5% GC; `GCF_000498675.1`) because transcribed coding sequences and structural RNAs in AT-rich Firmicutes are higher in GC than intergenic spacers. A single unimodal curve without extra bumps indicates our reads come from the target bacterium without contamination from another organism.

![Per sequence GC content percentages across all 12 files](figures/fastqc_per-sequence-gc.png)

### Sequence Duplication Levels (`FAIL`) and Overrepresented Sequences (`WARN`)

FastQC expects most reads to be unique, which is true for whole-genome DNA sequencing. In bacterial RNA-seq, however, highly active genes make many copies of the same transcript inside a small 1.89-Mb genome, so identical reads appear often (only ~13%–17% of reads would remain after deduplication). Likewise, the flagged overrepresented sequences each make up less than 0.6% of a file and mark `No Hit` in FastQC's built-in contaminant/adapter list because they map to abundant endogenous *L. johnsonii* N6.2 RNAs (such as RNase P RNA `rnpB` and tmRNA `ssrA`), not adapter dimers.

![Sequence duplication level curves across all 12 FASTQ files](figures/fastqc_duplication-levels.png)

### Adapter Content (`PASS`)

Sequencing adapters can appear at the 3' end if a cDNA fragment is shorter than 150 bp. Adapter levels stay near 0% for most of the read and reach at most ~0.6% in Control and ~2.4% in Bile at the 3' tail, well below FastQC's 5% warning cutoff.

![Adapter content across read positions for all 12 files](figures/fastqc_adapter-content.png)

## 4. Conclusion and Next Steps

The Group 04 reads are high quality and balanced across all six replicates. Before mapping them to the *L. johnsonii* N6.2 reference genome (`GCF_000498675.1`) with `HISAT2` or `Bowtie2`, we only need light 3' adapter and quality trimming (`Trimmomatic` or `fastp`), without removing duplicate reads so gene expression counts (`HTSeq` / `featureCounts`) stay accurate.

## References

1. Beliakoff, R. E., Gonzalez, C. F., & Lorca, G. L. (2024). Bile promotes *Lactobacillus johnsonii* N6.2 extracellular vesicle production with conserved immunomodulatory properties. *Scientific Reports*, 14(1), 12272. https://doi.org/10.1038/s41598-024-62843-0
2. National Center for Biotechnology Information (NCBI). BioProject `PRJNA964723` (SRA runs `SRR24376031`–`SRR24376036`): *Lactobacillus johnsonii* N6.2 bile-stress RNA-seq. https://www.ncbi.nlm.nih.gov/bioproject/PRJNA964723
3. National Center for Biotechnology Information (NCBI). RefSeq Assembly `GCF_000498675.1` (`ASM49867v1` / `NC_022909.1`): *Lactobacillus johnsonii* N6.2 complete genome. https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000498675.1/
4. Andrews, S. (2010). *FastQC: A Quality Control Tool for High Throughput Sequence Data* [Computer software; v0.13.0 used]. Babraham Bioinformatics. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
5. Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354
