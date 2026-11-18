# BT338IU Practice in Bioinformatics — Workflow Automation Master Plan

This playbook maps every exercise in [lab-manual.pdf](../references/lab-manual.pdf) and [slides_week-01.pdf](../lecture-01/slides_week-01.pdf) to fast, local, reproducible CLI commands built for your Python (`sbvs` / `uv`) and R (`ngs` / `renv`) stack.

---

## 1. Environment Architecture

1. **Python Orchestration (`uv`):** Manages lightweight API clients, spreadsheet parsing, plotting, and automated testing (`pyproject.toml` + `uv.lock`).
2. **Mamba `ngs` Environment:** Provides multi-threaded C/C++/Java NGS and phylogenetic binaries (`fasterq-dump`, `fastqc`, `fastp`, `trimmomatic`, `hisat2`, `star`, `salmon`, `samtools`, `bcftools`, `subread`/`featureCounts`, `htseq`, `multiqc`, `mafft`, `fasttree`, `blast`) plus pre-compiled R/Bioconductor binaries (`DESeq2`, `edgeR`, `GEOquery`, `clusterProfiler`, `enrichplot`, `ggtree`, `ape`).
3. **R `renv` Lockfile (`renv.lock`):** Snapshots exact R package versions inside the repository so `dge_analysis.R` stays reproducible across sessions.
4. **Mamba `sbvs` Environment:** Provides `pymol-open-source`, `pdb2pqr`, `apbs`, `modeller` (`mod10.7`), and `biopython` for headless 3D rendering, offline electrostatics, and homology modeling.

---

## 2. Lecture-by-Lecture Execution Playbook

### Lecture 1 — Biological Databases, BLAST & Phylogenetics

- **Database Exploration (Exercises 1 & 2):**
  ```bash
  # Run preset Week 1 Exercise 1 or Exercise 2 queries across NCBI, UniProt, PubChem & RCSB PDB
  uv run python3 src/db_miner.py --preset w1-ex1 -o lecture-01/ex01-db
  uv run python3 src/db_miner.py --preset w1-ex2 -o lecture-01/ex02-db

  # Single-target queries (pairs with ncbi-fetch on $PATH)
  uv run python3 src/db_miner.py -s pubchem -q "Curcumin" --stdout
  uv run python3 src/db_miner.py -s uniprot -q "EGFR human" --stdout
  uv run python3 src/db_miner.py -s pdb -q "GLP1R" --stdout
  ncbi-fetch -d nuccore -i "NR_024570.1" --format fasta --stdout
  ```
- **BLAST & Phylogenetic Tree (Exercise 3 & Tree Exercise):**
  ```bash
  # Fetch top 100 homologs via BLAST, align with MAFFT, build FastTree + ggtree PNG + MEGA-ready FASTA/Newick
  uv run python3 src/phylo_pipeline.py -a AAB72175.1 -p blastp -m 100 -o lecture-01/phylogeny-cytc
  uv run python3 src/phylo_pipeline.py -a NR_024570.1 -p blastn -m 100 -o lecture-01/phylogeny-16s
  ```

### Lecture 2 — RNA-Seq QC, Alignment, Quantification & Variant Calling

Skip 40–66 minute UseGalaxy queues by running locally across 8–12 threads (`references/group-datasets.xlsx`):

```bash
# Inspect all parsed group datasets from references/group-datasets.xlsx
uv run python3 src/rnaseq_pipeline.py --list-groups

# Fast modern run (fasterq-dump -> fastp/FastQC/MultiQC -> HISAT2/STAR/Salmon -> featureCounts)
# plus --emit-compat-logs to format Dr. Nhut's slide-ready Trimmomatic & HISAT2 summary blocks
uv run python3 src/rnaseq_pipeline.py --group DATASET_G1 --engine modern --emit-compat-logs -t 8 -o lecture-02/g1

# Exact manual compatibility run (Trimmomatic -> HISAT2 -> SAMtools sort -> HTSeq + bcftools VCF)
uv run python3 src/rnaseq_pipeline.py --group DATASET_G1 --engine compat --call-variants -t 8 -o lecture-02/g1-compat
```

### Lecture 3 — Protein Visualization & Structural Analysis (PyMOL + Offline APBS)

Automates all 9 required PDB exercises (`1A3N, 1GZX, 3NIR, 2MNR, 1QTI` and `1A2C, 1HCL, 1MOE, 4HHB`) headlessly, skipping manual PDB2PQR web uploads:

```bash
# Render 4 protein structure levels (primary, secondary, CYS disulfides, quaternary chains)
uv run python3 src/structure_viz.py --preset lecture3-levels -o lecture-03

# Render hydrophobicity (color_h.py), offline APBS electrostatics, H-bonds, and 3D monomer alignments
uv run python3 src/structure_viz.py --preset lecture3-props -o lecture-03
```

Open any generated `.pse` file in PyMOL (`mamba run -n sbvs pymol lecture-03/sessions/1A2C.pse`) anytime you want to rotate or capture custom angles manually.

### Lecture 4 — Protein Structure Prediction (MODELLER & AlphaFold 3)

- **Part A — Homology Modeling (`MODELLER`):**
  > Ensure your academic key is exported once in your shell or `~/.bashrc`: `export KEY_MODELLER="YOUR_KEY"`.
  ```bash
  # Run TvLDH tutorial or all 5 graded MODELLER cases (ex1..ex5)
  uv run python3 src/homology_model.py --case all --num-models 5 -o lecture-04/modeller
  ```
  Produces `.pir`, `.ali`, `align2d.log`, `model-single.log`, best `.pdb` ranked by DOPE/GA341, PyMOL superposition PNGs, and a reflection metrics table.

- **Part B — AI Prediction (`AlphaFold`):**
  ```bash
  # Fetch instant AFDB structures & PAE JSONs for P00698, P01308, P00338, P42212, P00004
  # Renders pLDDT-colored PyMOL PNGs, PAE heatmaps, Markdown report, and AF3 Server batch JSON
  uv run python3 src/alphafold_client.py --preset lecture4-af -o lecture-04/alphafold
  ```
  For multi-chain complexes (such as the `4× P00338` LDH homotetramer or cleaved `P01308` Insulin), upload the generated `lecture-04/alphafold/jobs/af3_server_batch.json` directly to [alphafoldserver.com](https://alphafoldserver.com/) with one click.

### Lecture 5 — Differential Gene Expression (DGE) & GO Enrichment Case Studies

Run parameterized `DESeq2` / `edgeR` + `clusterProfiler` analysis from local count matrices or directly from NCBI GEO:

```bash
# Run any of the 7 Lecture 5 case studies (--case cs1, cs2, cs3, cs4a, cs4b, cs5, cs6)
mamba run -n ngs Rscript src/dge_analysis.R --case cs2 -o lecture-05/cs2_ecoli
mamba run -n ngs Rscript src/dge_analysis.R --case cs6 -o lecture-05/cs6_curcumin

# Or pass your own featureCounts/HTSeq count matrix + sample metadata
mamba run -n ngs Rscript src/dge_analysis.R \
    -c lecture-02/g1/counts/counts_matrix.tsv \
    -m lecture-02/g1/metadata.tsv \
    --method deseq2 --lfc 1.0 --padj 0.05 -o lecture-02/g1/dge
```

### Converting Reports to Submission Word (`.docx`) Files

After drafting a report in Markdown (applying `/humanize`), render it to HCMIU School of Biotechnology formatting with your existing CLI tool:

```bash
format-docx -m render lecture-01/report-01_biological-databases.md -o lecture-01/report-01_biological-databases.docx -p hcmiu-sbt
```
