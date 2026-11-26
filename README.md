# BT338IU: Practice in Bioinformatics (`prac-bioinfo-2026`)

Local, reproducible Python and R automation suite for **BT338IU – Practice in Bioinformatics** at International University – VNUHCM (School of Biotechnology; Instructor: Tran Minh Nhựt, PhD).

Designed to replace slow web-queue bottlenecks (UseGalaxy, PDB2PQR server, AlphaFold Server) and repetitive GUI clicking with fast, customizable local pipelines whose outputs match Dr. Nhựt's lecture slides and submission checklists 1:1.

---

## Repository Layout

```text
prac-bioinfo-2026/
├── AGENTS.md                                  # Report tone & environment instructions
├── pyproject.toml                             # uv Python dependencies & tooling
├── renv.lock                                  # Locked R / Bioconductor environment
├── docs/
│   └── workflow-automation_master-plan.md     # Lecture-by-lecture execution playbook
├── homework/                                  # Weekly homework QC, trimming, and count deliverables
├── prelab/                                    # Pre-lab preparation notes, figures, and assignments
├── references/
│   ├── lab-manual.pdf                         # Official 74-page BT338IU lab manual
│   ├── lab-manual_overview.md                 # Structured summary of tools & workload
│   └── group-datasets.xlsx                    # Class RNA-seq SRA & reference genome sheets
├── lecture-01/                                # Database mining, BLAST & phylogeny work
├── lecture-02/                                # Group RNA-seq pipeline outputs
├── lecture-03/                                # PyMOL 3D protein visualization outputs
├── lecture-04/                                # MODELLER homology models & AlphaFold outputs
├── lecture-05/                                # NCBI GEO Differential Gene Expression case studies
├── src/                                       # Self-documenting Python & R CLI tools
└── tests/                                     # CLI contract & offline unit tests
```

---

## Tooling Matrix

| Lecture(s) | Manual Workflow | Local Accelerated CLI | Runtime Environment |
| :--- | :--- | :--- | :--- |
| **1** | Manual NCBI/UniProt/PubChem/PDB browsing | `uv run python3 src/db_miner.py` (+ `ncbi-fetch`) | `uv` + `$PATH` |
| **1** | Web BLAST + MEGA GUI | `uv run python3 src/phylo_pipeline.py` | `uv` + `ngs` (`mafft`, `FastTree`, `blast`) |
| **1–2** | UseGalaxy web queues (`FastQC → HISAT2 → HTSeq`) | `uv run python3 src/rnaseq_pipeline.py` | `uv` + `ngs` (`fastp`, `Salmon`/`STAR`/`HISAT2`, `featureCounts`) |
| **2, 5** | Manual RStudio `DESeq2` + commercial `Blast2GO` | `mamba run -n ngs Rscript src/dge_analysis.R` | `ngs` + `renv` (`DESeq2`, `GEOquery`, `clusterProfiler`) |
| **3** | PyMOL GUI + PDB2PQR web server | `uv run python3 src/structure_viz.py` | `uv` + `sbvs` (`pymol`, `pdb2pqr`, `apbs`) |
| **4A** | Windows PowerShell `mod10.7` (15 manual runs) | `uv run python3 src/homology_model.py` | `uv` + `sbvs` (`modeller`, `pymol`) |
| **4B** | `alphafoldserver.com` wait queues | `uv run python3 src/alphafold_client.py` | `uv` + `sbvs` (`pymol`) |

---

## Quickstart

```bash
# Inspect any CLI contract
uv run python3 src/db_miner.py --help
uv run python3 src/rnaseq_pipeline.py --list-groups
mamba run -n ngs Rscript src/dge_analysis.R --help

# Run verification tests
uv run ruff check src tests
uv run pytest tests/
```

See [workflow-automation_master-plan.md](docs/workflow-automation_master-plan.md) for the complete week-by-week command guide.
