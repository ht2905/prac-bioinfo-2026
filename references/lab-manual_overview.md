# Practice in Bioinformatics (`BT338IU`) — Lab Manual Overview

- **Source:** [lab-manual.pdf](file:///home/ht2905/Documents/undergrad-materials-2024-2028/practice-in-bioinformatics/lab-manual.pdf) (74 pages)
- **Institution:** International University – VNUHCM, School of Biotechnology
- **Instructor:** Tran Minh Nhựt, PhD *(Manual composed/edited by Thu Vo, PhD & Nhut Tran, PhD — Nov 2025)*

---

## 1. Course Snapshot

`BT338IU` is an **applied, workflow-driven bioinformatics lab course** centered on two pillars:
1. **Transcriptomics (RNA-Seq & Differential Gene Expression):** Going from raw sequencing reads (`FASTQ`) through QC, alignment, read counting, statistical differential expression (`DESeq2`), and biological pathway interpretation across real NCBI GEO case studies.
2. **Structural Bioinformatics:** Visualizing 3D protein structures (`PyMOL`), running template-based homology modeling (`MODELLER`), and running AI-based 3D structure prediction (`AlphaFold 3`).

---

## 2. Tools, Platforms & Databases Expected

### A. Core Hands-On Tools (You Will Actively Use These)

| Tool / Platform | Role in Course | Where Used | Setup / Access Notes |
| :--- | :--- | :--- | :--- |
| **Galaxy** (`usegalaxy.eu`, `.org`, `.org.au`, `.fr`) | Web GUI workflow platform for NGS/RNA-seq pipelines (histories, tool runs, workflow extraction) | Lectures 1–2 | Free web account; browser-based |
| **R & RStudio** + **`DESeq2`** *(optional `edgeR`)* | Statistical Differential Gene Expression (DGE) analysis from read counts | Lectures 2, 5 | Local install; Bioconductor package (`DESeq2`) |
| **PyMOL** (Educational v3.1.5.1 shown) | 3D protein visualization, structure levels, hydrophobicity, electrostatics, H-bonds, alignment | Lectures 3–4 | Free Educational license (requires registration & email license file); GUI + short CLI commands + `.py` scripts |
| **MODELLER** (`mod10.7`) | Template-based protein homology modeling via pre-written Python scripts | Lecture 4 | Local install; free academic license key via email; invoked from terminal/PowerShell |
| **AlphaFold 3 Server** (`alphafoldserver.com`) | AI-based 3D structure prediction (single chains & complexes) | Lecture 4 | Web interface; requires sequence < 2,000 aa |

### B. RNA-Seq Pipeline & Annotation Tools (Mostly Run Inside Galaxy or Web)

- **Pre-alignment QC & Trimming:** **FastQC** (Phred quality checks) and **Trimmomatic** (adapter & quality trimming).
- **Alignment & BAM Handling:** **HISAT2** (splice-aware alignment to reference genome) and **SAMtools** (`.sam` → sorted `.bam`).
- **Quantification:** **HTSeq** (`htseq-count`) and **featureCounts** (mentioned in Case Study 4, p. 70) to produce gene count tables.
- **Functional Enrichment / GO Annotation:** **Blast2GO** (cloud BLAST + InterProScan + GO mapping), **ShinyGO 0.85** (linked in Lecture 1), and **LINCS/CMap** (mentioned as a learning outcome in Case Study 6).
- **Transcriptome Variant Calling (Lecture 2):** **SAMtools** + **BCFtools** (→ `.vcf`) visualized in **IGV** (Integrative Genomics Viewer).

### C. Biological Databases & Web Servers

- **Sequence & Literature:** **NCBI Entrez**, **GenBank**, **RefSeq**, **NCBI-BLAST**, **NCBI GEO** (all Lecture 5 datasets come from GEO `GSE*` accessions), **UniProt** (SWISS-Prot / TrEMBL).
- **Structure & Electrostatics:** **RCSB PDB** (`www.rcsb.org`), **PDB2PQR Server** (`server.poissonboltzmann.org/pdb2pqr` — generates `.pqr` and `.dx` files for PyMOL electrostatic surfaces), **AlphaFold DB**.
- **Web Utilities (Lecture 1):** **EBI ClustalW**, **ExPASy**.

### D. Mentioned-Only / Background Tools (No Step-by-Step Lab Exercises)

- **Qiagen CLC Main/Genomics Workbench**, **Geneious Prime**, **MEGA** (listed in Lecture 1 software overview).
- **UCSF Chimera** and **VMD** (mentioned as alternatives alongside PyMOL; Chimera is also mentioned in Lecture 4 step 4).
- **Notepad++** (suggested on Windows for editing `.PIR` files and reading `.log` files).
- **Basic Unix/Linux CLI concepts** (introduced conceptually in Lecture 1, though Galaxy replaces most raw CLI work).

---

## 3. Curriculum & Deliverables Breakdown

| Unit | Pages | Core Topics | Required Practice / Deliverables |
| :--- | :---: | :--- | :--- |
| **Lecture 1:** Intro to Bioinformatics | 2–19 | Databases, pitfalls, Entrez/GenBank, FASTA/FASTQ/SAM/BAM/VCF/GFF/BED/PDB, Galaxy UI, FastQC, Trimmomatic | Galaxy setup + **Homework: Quality controls** |
| **Lecture 2:** Galaxy & RNA-Seq Pipeline | 20–26 | Transcriptomics pipeline (`FastQC → Trimmomatic → HISAT2 → SAMtools → HTSeq → DESeq2 → Blast2GO`), variant calling (`bcftools`, `IGV`) | **Group Project:** Species-specific differential gene expression analysis on Galaxy |
| **Lecture 3:** Protein Visualization | 27–50 | PyMOL GUI/CLI, 4 structural levels, hydrophobicity (`color_h.py`), electrostatics (`PDB2PQR`), interfaces (`InterfaceResidues.py`), H-bonds, 3D alignment | **9 PDB structures total:** 5 structures for 4-level analysis (`1A3N`, `1GZX`, `3NIR`, `2MNR`, `1QTI`) + 4 structures for properties/H-bonds/alignment (`1A2C`, `1HCL`, `1MOE`, `4HHB`) |
| **Lecture 4:** Structure Prediction | 51–63 | Homology modeling (`MODELLER`: `.PIR`, `build_profile.py`, `align2d.py`, `model-single.py`, DOPE/GA341) + AI prediction (`AlphaFold 3`: pLDDT, PAE) | **5 MODELLER cases** (logs, best `.pdb`, PyMOL screenshot, 100–150 word reflection each) + **5 AlphaFold 3 predictions** (summaries of length, pLDDT, low-confidence regions, structural insights) |
| **Lecture 5:** DGE Case Studies | 64–74 | 7 real-world RNA-Seq studies (iPSC reprogramming `GSE33110`, *E. coli* ciprofloxacin `GSE220559`, *Arabidopsis* immunity `GSE80448`, *P. aeruginosa* `GSE120602`, mouse wound healing `GSE154132`, human macrophage LPS `GSE193336`/`GSE19420`, curcumin in colon cancer `GSE229613`) | Mini-project(s): Run `DESeq2` in R, list **top 10 upregulated & top 10 downregulated genes**, write biological interpretation reports (~800 words max) |

---

## 4. How Deep / Heavy Is This Subject?

### Overall Difficulty Profile: **Light–Moderate Coding, Moderate Theory, Moderate–High Lab Volume**

1. **Programming & Math Depth — Introductory ("Tool-User" Level)**
   - You are **not** expected to derive dynamic programming alignment math, write NGS pipelines in Bash from scratch, or program bioinformatics algorithms in Python.
   - **Galaxy abstracts away the Unix command line** for heavy NGS steps (`FastQC`, `HISAT2`, `HTSeq`).
   - **MODELLER** runs pre-built tutorial scripts (`build_profile.py`, `align2d.py`, `model-single.py`)—you edit plain-text `.PIR` sequence files and inspect `.log` tables rather than coding Python.
   - **Where real technical friction exists:** **R / RStudio (`DESeq2`)**. The manual does *not* include copy-pasteable R scripts—it points to the official `DESeq2` vignette/tutorials. If you have never used R, installing Bioconductor packages, formatting count matrices/metadata tables, and filtering by `log2FoldChange` and `padj` will be the steepest learning curve.

2. **Biological & Analytical Depth — Moderate (Core Focus of Grading)**
   - The professor repeatedly emphasizes **critical interpretation over blind button-clicking**:
     - Lecture 1 warns about database redundancy, misannotations, and algorithmic errors.
     - Lecture 3 expects you to connect 3D features (disulfide bridges, hydrophobic patches, electrostatic surfaces, polar contacts) to function.
     - Lecture 4 specifically tests whether you understand *why* homology modeling degrades below ~40% sequence identity (`CYP2E1` at ~37% vs `DHFR` at ~99%) and how to read `DOPE`, `GA341`, `pLDDT`, and `PAE`.
     - Lecture 5 asks you to connect top DEGs to specific pathways (e.g., bacterial SOS response `recA`/`lexA`/`sulA`, plant `PR1`/`NPR1`/`WRKY`, macrophage `NF-κB`/`IRF`/`STAT`).

3. **Workload & Time Commitment — Moderate–High (Due to Repetitive Multi-Sample Runs)**
   - There are **many individual deliverables** (9 PyMOL structures + 5 MODELLER cases + 5 AlphaFold predictions + Galaxy group project + DGE case study report).
   - **Queue/compute bottlenecks:** Public **UseGalaxy** servers and **AlphaFold Server** can take anywhere from minutes to **several hours** for medium/multimeric proteins (e.g., LDH homotetramer) or full RNA-seq datasets. Leaving labs to the night before the deadline is risky.

---

## 5. Practical Caveats & Quirks in the Lab Manual

> [!WARNING]
> **Things to clarify or prepare ahead of time** based on inspecting the manual text:

- **External Scripts & Links Required:** The manual relies on external URLs and scripts not embedded in the PDF—specifically PyMOL helper scripts (`color_h.py`, `InterfaceResidues.py`), MODELLER tutorial scripts (`salilab.org/modeller/tutorial/basic.html`), Galaxy Training Network tutorials, and the `DESeq2` guide.
- **License Registration Lead Times:** Both **PyMOL Educational** and **MODELLER** require academic email registration to receive credentials/license keys before you can run the software.
- **Windows-Centric Instructions (Easily Adapted to Linux):** Lecture 4 shows Windows PowerShell (`Run as administrator`, `C:\Users\...`) and Notepad++, but **PyMOL**, **MODELLER**, **R/RStudio**, and **Galaxy/AlphaFold** all run natively on Linux.
- **Minor Manual Typos / Discrepancies to Confirm with the Professor:**
  1. **Duplicate Case Study Numbering:** There are two sections titled *"Case Study 4"*—one for *Pseudomonas aeruginosa* (`GSE120602`, p. 68) and one for *Mus musculus* wound healing (`GSE154132`, p. 70).
  2. **Dataset ID Mismatch in Case Study 5 (p. 72):** The background names **`GSE193336`**, whereas the mini-project prompt asks you to use **`GSE19420`**.
  3. **Minor Wording Typos:** Case Study 2 (pp. 65–66) writes *"antibodies"* where *"antibiotics"* (ciprofloxacin) is meant; PyMOL §5.5 (p. 48) has minor filename typos (`1hag-1.pdb` vs `1h4g-2.pdb`).
