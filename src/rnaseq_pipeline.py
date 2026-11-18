#!/usr/bin/env python3
"""Run local multi-threaded RNA-seq QC, alignment, quantification, and variant calling pipelines."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import subprocess
import sys
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

import openpyxl
import requests

DEFAULT_XLSX = Path("references/group-datasets.xlsx")
NCBI_DATASETS_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download"
SRR_PATTERN = re.compile(r"SRR\d+")
GCF_PATTERN = re.compile(r"GCF_\d+\.\d+|GCA_\d+\.\d+")
PRJNA_PATTERN = re.compile(r"PRJNA\d+")


@dataclass
class DatasetGroup:
    group_id: str
    sheet: str
    label: str
    bioproject: str = ""
    organism: str = ""
    genome_gcf: str = ""
    conditions: dict[str, list[str]] = field(default_factory=dict)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run local RNA-seq pipelines (fastp/Salmon/STAR/HISAT2/featureCounts) "
            "with Trimmomatic/HISAT2 compat logs."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  # List all parsed class dataset IDs from references/group-datasets.xlsx
  rnaseq_pipeline.py --list-groups

  # Run modern pipeline with Dr. Nhut compatibility logs for dataBI Group 1
  rnaseq_pipeline.py -g DATABI_G1 --engine modern --emit-compat-logs -t 8 -o lecture-02/g1

  # Run exact manual compatibility pipeline (Trimmomatic + HISAT2 + HTSeq + BCFtools VCF)
  rnaseq_pipeline.py -g DATABI_G1 --engine compat --call-variants -o lecture-02/g1-compat
""",
    )
    parser.add_argument(
        "-g",
        "--group",
        type=str,
        default=None,
        help="Group identifier from --list-groups (e.g. DATABI_G1, DATASET_G4).",
    )
    parser.add_argument(
        "--xlsx",
        type=Path,
        default=DEFAULT_XLSX,
        help="Path to group datasets Excel workbook (default: references/group-datasets.xlsx).",
    )
    parser.add_argument(
        "--list-groups",
        action="store_true",
        help="Parse --xlsx and print all available group IDs, organisms, and SRR counts.",
    )
    parser.add_argument(
        "--sra",
        type=str,
        default=None,
        help="Comma-separated SRR accessions (overrides --group SRRs).",
    )
    parser.add_argument(
        "--genome",
        type=str,
        default=None,
        help="NCBI RefSeq/GenBank assembly accession (e.g. GCF_000005845.2).",
    )
    parser.add_argument(
        "-e",
        "--engine",
        choices=("modern", "compat", "both"),
        default="modern",
        help=(
            "Backend: modern (fastp+featureCounts), compat (Trimmomatic+HISAT2+HTSeq), "
            "or both (default: modern)."
        ),
    )
    parser.add_argument(
        "--aligner",
        choices=("hisat2", "star", "salmon"),
        default="hisat2",
        help="Alignment/pseudo-alignment engine in modern mode (default: hisat2).",
    )
    parser.add_argument(
        "--quantifier",
        choices=("featurecounts", "htseq"),
        default="featurecounts",
        help="Read quantifier for BAM files (default: featurecounts).",
    )
    parser.add_argument(
        "--emit-compat-logs",
        action="store_true",
        help="Write Trimmomatic and HISAT2-formatted summary logs from modern runs.",
    )
    parser.add_argument(
        "--call-variants",
        action="store_true",
        help="Run samtools/bcftools variant calling to generate a VCF file.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=8,
        help="Worker threads for downloads, QC, and alignment (default: 8).",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("lecture-02"),
        help="Output directory for FASTQ, QC, alignments, counts, and logs (default: lecture-02).",
    )
    parser.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Preview pipeline steps and commands without downloading or running tools.",
    )
    args = parser.parse_args(argv)
    if not args.list_groups and not args.group and not (args.sra and args.genome):
        parser.error(
            "Provide --list-groups, -g/--group, or --sra + --genome. "
            "Run 'rnaseq_pipeline.py --help' for usage."
        )
    return args


def _clean_cell(val: object) -> str:
    return re.sub(r"\s+", " ", str(val or "")).strip()


def _extract_organism(text: str) -> str:
    text = text.replace("Organism:", "").strip()
    parts = [p.strip() for p in text.split(";") if p.strip()]
    for part in parts:
        if not PRJNA_PATTERN.search(part) and not GCF_PATTERN.search(part):
            return part
    return ""


def _slugify_condition(raw_label: str, existing: dict[str, list[str]], fallback_idx: int) -> str:
    cleaned = re.sub(r"^(?:GSM\d+:\s*|RNA-[Ss]eq\s+of\s+)", "", raw_label.strip())
    slug = re.sub(r"[^a-zA-Z0-9._-]+", "_", cleaned).strip("_")[:60] or f"Condition_{fallback_idx}"
    if slug in existing:
        slug = f"{slug}_{fallback_idx}"
    return slug


def parse_groups_xlsx(xlsx_path: Path) -> dict[str, DatasetGroup]:
    """Extract structured SRA & reference genome metadata across workbook sheets."""
    if not xlsx_path.exists():
        raise FileNotFoundError(f"Workbook not found: {xlsx_path}")

    wb = openpyxl.load_workbook(xlsx_path, read_only=True, data_only=True)
    groups: dict[str, DatasetGroup] = {}

    if "dataBI" in wb.sheetnames:
        ws = wb["dataBI"]
        rows = [[_clean_cell(c) for c in row] for row in ws.iter_rows(values_only=True)]
        current_meta: list[str] = []
        pending_label = ""
        for row in rows:
            if not any(row):
                continue
            col1 = row[0] if len(row) > 0 else ""
            col2 = row[1] if len(row) > 1 else ""
            row_has_srr = any(SRR_PATTERN.search(c) for c in row[2:])
            is_group_tag = bool(col1 and (col1.startswith("G") or col1.lower() == "demo"))

            if not row_has_srr and (
                is_group_tag
                or PRJNA_PATTERN.search(col2)
                or (len(row) > 2 and row[2] and not SRR_PATTERN.search(row[2]))
            ):
                if is_group_tag:
                    pending_label = col1
                current_meta = row
                continue

            if row_has_srr and (is_group_tag or pending_label):
                label = col1 if is_group_tag else pending_label
                gid = f"DATABI_{label.upper()}"
                prj = ""
                gcf = ""
                org = ""
                for text in [col2, *(current_meta[1:2] if current_meta else [])]:
                    if (m := PRJNA_PATTERN.search(text)) and not prj:
                        prj = m.group(0)
                    if (m := GCF_PATTERN.search(text)) and not gcf:
                        gcf = m.group(0)
                    if not org:
                        org = _extract_organism(text)
                conds: dict[str, list[str]] = {}
                for idx in range(2, max(len(row), len(current_meta))):
                    cell_srr = row[idx] if idx < len(row) else ""
                    srrs = SRR_PATTERN.findall(cell_srr)
                    if srrs:
                        cond_label = (
                            current_meta[idx]
                            if idx < len(current_meta) and current_meta[idx]
                            else f"Condition_{idx - 1}"
                        )
                        cond_slug = _slugify_condition(cond_label, conds, idx - 1)
                        conds[cond_slug] = srrs
                if conds:
                    groups[gid] = DatasetGroup(
                        group_id=gid,
                        sheet="dataBI",
                        label=label,
                        bioproject=prj,
                        organism=org,
                        genome_gcf=gcf,
                        conditions=conds,
                    )
                current_meta = []
                pending_label = ""

    if "DATASET" in wb.sheetnames:
        ws = wb["DATASET"]
        rows = [[_clean_cell(c) for c in row] for row in ws.iter_rows(values_only=True)]
        active: DatasetGroup | None = None
        last_labels: dict[int, str] = {}
        label_to_slug: dict[str, str] = {}
        for row in rows:
            if not any(row):
                continue
            col1 = row[0] if len(row) > 0 else ""
            col3 = row[2] if len(row) > 2 else ""
            if col1.startswith("G") or col1 == "CLASS DEMO":
                slug = re.sub(r"\s+", "_", col1.upper())
                gid = f"DATASET_{slug}"
                active = DatasetGroup(group_id=gid, sheet="DATASET", label=col1)
                groups[gid] = active
                last_labels = {}
                label_to_slug = {}
            elif active is None:
                continue

            if m := PRJNA_PATTERN.search(col3):
                active.bioproject = active.bioproject or m.group(0)
            if m := GCF_PATTERN.search(col3):
                active.genome_gcf = active.genome_gcf or m.group(0)
            if not active.organism:
                active.organism = _extract_organism(col3)

            for idx in range(4, len(row)):
                val = row[idx]
                if not val:
                    continue
                srrs = SRR_PATTERN.findall(val)
                if srrs:
                    lbl = last_labels.get(idx) or f"Condition_{idx - 3}"
                    if lbl not in label_to_slug:
                        label_to_slug[lbl] = _slugify_condition(lbl, active.conditions, idx - 3)
                    cond_slug = label_to_slug[lbl]
                    active.conditions.setdefault(cond_slug, []).extend(srrs)
                elif "GCF_" not in val and "minutes" not in val:
                    last_labels[idx] = val

    return groups


def format_groups_table(groups: dict[str, DatasetGroup]) -> str:
    lines = [
        "| Group ID | BioProject | Organism | Genome Assembly | Samples (SRRs) |",
        "| :--- | :--- | :--- | :--- | ---: |",
    ]
    for gid, grp in groups.items():
        n_srr = sum(len(v) for v in grp.conditions.values())
        org = grp.organism or "N/A"
        lines.append(
            f"| `{gid}` | `{grp.bioproject or 'N/A'}` | *{org}* | "
            f"`{grp.genome_gcf or 'N/A'}` | {n_srr} |"
        )
    return "\n".join(lines)


def mamba_run(args: list[str], capture: bool = False) -> subprocess.CompletedProcess[str]:
    cmd = ["mamba", "run", "-n", "ngs", *args]
    return subprocess.run(cmd, capture_output=capture, text=True, check=True)


def fetch_genome_bundle(accession: str, ref_dir: Path) -> tuple[Path, Path | None, Path | None]:
    """Download genome FASTA, GTF/GFF, and RNA FASTA via NCBI Datasets v2 REST API."""
    ref_dir.mkdir(parents=True, exist_ok=True)
    zip_path = ref_dir / f"{accession}_bundle.zip"
    if not zip_path.exists():
        url = NCBI_DATASETS_URL.format(accession=accession)
        print(f"[INFO] Downloading NCBI reference package {accession}...", file=sys.stderr)
        resp = requests.get(
            url,
            params={"include_annotation_type": "GENOME_FASTA,GENOME_GFF,GENOME_GTF,RNA_FASTA"},
            stream=True,
            timeout=180,
        )
        resp.raise_for_status()
        with zip_path.open("wb") as fh:
            for chunk in resp.iter_content(chunk_size=1024 * 256):
                fh.write(chunk)

    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(ref_dir)

    fna = next(ref_dir.rglob("*_genomic.fna"), None)
    gtf = next(ref_dir.rglob("*.gtf"), None)
    gff = next(ref_dir.rglob("*.gff"), None)
    rna = next(ref_dir.rglob("rna.fna"), None)
    if fna is None:
        raise RuntimeError(f"No genomic FASTA found inside downloaded bundle for {accession}.")
    return fna, gtf or gff, rna


def download_srr_fastqs(srr: str, fastq_dir: Path, threads: int) -> tuple[Path, Path | None]:
    """Download FASTQ reads for an SRR accession using fasterq-dump."""
    fastq_dir.mkdir(parents=True, exist_ok=True)
    pe1 = fastq_dir / f"{srr}_1.fastq"
    pe2 = fastq_dir / f"{srr}_2.fastq"
    se = fastq_dir / f"{srr}.fastq"
    if pe1.exists() and pe2.exists():
        return pe1, pe2
    if se.exists():
        return se, None

    print(f"[INFO] Downloading {srr} via fasterq-dump ({threads} threads)...", file=sys.stderr)
    mamba_run(["fasterq-dump", "--split-3", "--threads", str(threads), "-O", str(fastq_dir), srr])
    if pe1.exists() and pe2.exists():
        return pe1, pe2
    if se.exists():
        return se, None
    raise FileNotFoundError(f"fasterq-dump completed but no FASTQ files found for {srr}.")


def write_compat_summary(
    sample: str,
    fastp_json: Path | None,
    hisat2_log: Path | None,
    out_path: Path,
) -> None:
    """Format modern QC/alignment metrics into Dr. Nhut's slide/manual summary layout."""
    lines = [f"=== Compatibility Summary: {sample} ==="]
    if fastp_json and fastp_json.exists():
        data = json.loads(fastp_json.read_text(encoding="utf-8"))
        before = data.get("summary", {}).get("before_filtering", {})
        after = data.get("summary", {}).get("after_filtering", {})
        in_reads = int(before.get("total_reads", 0))
        surv_reads = int(after.get("total_reads", 0))
        pct = (surv_reads / in_reads * 100.0) if in_reads else 0.0
        q30_before = before.get("q30_rate", 0.0) * 100
        q30_after = after.get("q30_rate", 0.0) * 100
        lines.extend(
            [
                "[Trimmomatic-equivalent Read Survival Summary]",
                (
                    f"Input Reads: {in_reads} Surviving: {surv_reads} ({pct:.2f}%) "
                    f"Dropped: {in_reads - surv_reads} ({100.0 - pct:.2f}%)"
                ),
                f"Q30 Before: {q30_before:.2f}% -> Q30 After: {q30_after:.2f}%",
                "",
            ]
        )
    if hisat2_log and hisat2_log.exists():
        lines.append("[HISAT2 Alignment Summary]")
        lines.append(hisat2_log.read_text(encoding="utf-8").strip())
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_featurecounts(
    bams: list[Path],
    sample_names: list[str],
    ann_path: Path,
    out_matrix: Path,
    threads: int,
    paired: bool,
) -> None:
    raw_fc = out_matrix.parent / "featurecounts_raw.tsv"
    fc_cmd = ["featureCounts", "-T", str(threads), "-a", str(ann_path), "-o", str(raw_fc)]
    if ann_path.suffix == ".gff":
        fc_cmd.extend(["-t", "gene", "-g", "ID"])
    if paired:
        fc_cmd.append("-p")
    fc_cmd.extend(str(b) for b in bams)
    mamba_run(fc_cmd)

    # Convert raw featureCounts output to a clean gene_id + sample count matrix for DESeq2
    with (
        raw_fc.open("r", encoding="utf-8") as src,
        out_matrix.open("w", encoding="utf-8", newline="") as dst,
    ):
        writer = csv.writer(dst, delimiter="\t")
        header_seen = False
        for line in src:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if not header_seen:
                writer.writerow(["gene_id", *sample_names])
                header_seen = True
                continue
            writer.writerow([cols[0], *cols[6:]])


def run_htseq(bams: list[Path], sample_names: list[str], ann_path: Path, out_matrix: Path) -> None:
    """Run htseq-count in an isolated Python 3.13 uvx environment (compat mode)."""
    counts_by_gene: dict[str, list[str]] = {}
    genes_order: list[str] = []
    for idx, bam in enumerate(bams):
        cmd = [
            "uvx",
            "--python",
            "3.13",
            "--from",
            "HTSeq",
            "htseq-count",
            "-f",
            "bam",
            "-r",
            "pos",
            "-s",
            "no",
            str(bam),
            str(ann_path),
        ]
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
        for line in res.stdout.splitlines():
            if line.startswith("__") or "\t" not in line:
                continue
            gid, cnt = line.split("\t", 1)
            if gid not in counts_by_gene:
                counts_by_gene[gid] = ["0"] * len(bams)
                genes_order.append(gid)
            counts_by_gene[gid][idx] = cnt

    with out_matrix.open("w", encoding="utf-8", newline="") as dst:
        writer = csv.writer(dst, delimiter="\t")
        writer.writerow(["gene_id", *sample_names])
        for gid in genes_order:
            writer.writerow([gid, *counts_by_gene[gid]])


def write_metadata_tsv(samples: list[tuple[str, str]], dest: Path) -> None:
    with dest.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["sample", "condition"])
        for sample_id, cond in samples:
            writer.writerow([sample_id, cond])


def execute_pipeline(
    args: argparse.Namespace, samples: list[tuple[str, str]], genome_acc: str
) -> int:
    out_dir = args.output_dir
    for sub in ("fastq", "qc", "trimmed", "ref", "align", "counts", "logs", "variants"):
        (out_dir / sub).mkdir(parents=True, exist_ok=True)

    write_metadata_tsv(samples, out_dir / "metadata.tsv")
    fna_path, ann_path, rna_path = fetch_genome_bundle(genome_acc, out_dir / "ref")

    use_hisat = args.engine in ("compat", "both") or args.aligner == "hisat2"
    hisat_prefix = out_dir / "ref" / "hisat2_idx"
    if use_hisat and not (out_dir / "ref" / "hisat2_idx.1.ht2").exists():
        print("[INFO] Building HISAT2 genome index...", file=sys.stderr)
        mamba_run(["hisat2-build", "-p", str(args.threads), str(fna_path), str(hisat_prefix)])

    salmon_idx = out_dir / "ref" / "salmon_idx"
    build_salmon = (
        args.engine in ("modern", "both")
        and args.aligner == "salmon"
        and rna_path is not None
        and not salmon_idx.exists()
    )
    if build_salmon and rna_path:
        print("[INFO] Building Salmon transcript index...", file=sys.stderr)
        mamba_run(
            ["salmon", "index", "-t", str(rna_path), "-i", str(salmon_idx), "-p", str(args.threads)]
        )

    star_idx = out_dir / "ref" / "star_idx"
    build_star = (
        args.engine in ("modern", "both")
        and args.aligner == "star"
        and not (star_idx / "SA").exists()
    )
    if build_star:
        genome_len = fna_path.stat().st_size
        sa_nbases = min(14, max(4, int(math.log2(max(genome_len, 1024)) / 2) - 1))
        star_idx.mkdir(parents=True, exist_ok=True)
        print(
            f"[INFO] Building STAR genome index (genomeSAindexNbases={sa_nbases})...",
            file=sys.stderr,
        )
        star_cmd = [
            "STAR",
            "--runMode",
            "genomeGenerate",
            "--runThreadN",
            str(args.threads),
            "--genomeDir",
            str(star_idx),
            "--genomeFastaFiles",
            str(fna_path),
            "--genomeSAindexNbases",
            str(sa_nbases),
        ]
        if ann_path and ann_path.suffix == ".gtf":
            star_cmd.extend(["--sjdbGTFfile", str(ann_path)])
        mamba_run(star_cmd)

    sorted_bams: list[Path] = []
    sample_names: list[str] = []
    is_paired_run = False

    for srr, _cond in samples:
        sample_names.append(srr)
        r1, r2 = download_srr_fastqs(srr, out_dir / "fastq", args.threads)
        is_paired = r2 is not None
        is_paired_run = is_paired_run or is_paired

        # FastQC on raw reads
        qc_inputs = [str(r1), str(r2)] if is_paired else [str(r1)]
        mamba_run(["fastqc", "-t", str(args.threads), "-o", str(out_dir / "qc"), *qc_inputs])

        trim1 = out_dir / "trimmed" / f"{srr}_trim_1.fastq.gz"
        trim2 = out_dir / "trimmed" / f"{srr}_trim_2.fastq.gz" if is_paired else None
        fastp_json = out_dir / "qc" / f"{srr}_fastp.json"
        fastp_html = out_dir / "qc" / f"{srr}_fastp.html"

        if args.engine in ("modern", "both"):
            fp_cmd = [
                "fastp",
                "-i",
                str(r1),
                "-o",
                str(trim1),
                "-w",
                str(min(args.threads, 16)),
                "-j",
                str(fastp_json),
                "-h",
                str(fastp_html),
            ]
            if is_paired and r2 and trim2:
                fp_cmd.extend(["-I", str(r2), "-O", str(trim2)])
            mamba_run(fp_cmd)
        else:
            if is_paired and r2 and trim2:
                unp1 = out_dir / "trimmed" / f"{srr}_unpaired_1.fastq.gz"
                unp2 = out_dir / "trimmed" / f"{srr}_unpaired_2.fastq.gz"
                mamba_run(
                    [
                        "trimmomatic",
                        "PE",
                        "-threads",
                        str(args.threads),
                        str(r1),
                        str(r2),
                        str(trim1),
                        str(unp1),
                        str(trim2),
                        str(unp2),
                        "LEADING:3",
                        "TRAILING:3",
                        "SLIDINGWINDOW:4:15",
                        "MINLEN:36",
                    ]
                )
            else:
                mamba_run(
                    [
                        "trimmomatic",
                        "SE",
                        "-threads",
                        str(args.threads),
                        str(r1),
                        str(trim1),
                        "LEADING:3",
                        "TRAILING:3",
                        "SLIDINGWINDOW:4:15",
                        "MINLEN:36",
                    ]
                )

        if args.aligner == "salmon" and args.engine == "modern" and salmon_idx.exists():
            quant_dir = out_dir / "counts" / f"{srr}_salmon"
            sal_cmd = [
                "salmon",
                "quant",
                "-i",
                str(salmon_idx),
                "-l",
                "A",
                "-p",
                str(args.threads),
                "-o",
                str(quant_dir),
            ]
            if is_paired and trim2:
                sal_cmd.extend(["-1", str(trim1), "-2", str(trim2)])
            else:
                sal_cmd.extend(["-r", str(trim1)])
            mamba_run(sal_cmd)
            if args.emit_compat_logs:
                log_out = out_dir / "logs" / f"{srr}_compat-summary.log"
                write_compat_summary(srr, fastp_json, None, log_out)
            continue

        hisat_log = out_dir / "logs" / f"{srr}_hisat2.log"
        sam_path = out_dir / "align" / f"{srr}.sam"
        bam_path = out_dir / "align" / f"{srr}.sorted.bam"

        if args.aligner == "star" and args.engine == "modern" and (star_idx / "SA").exists():
            prefix = str(out_dir / "align" / f"{srr}_")
            star_run = [
                "STAR",
                "--runThreadN",
                str(args.threads),
                "--genomeDir",
                str(star_idx),
                "--readFilesCommand",
                "zcat",
                "--outSAMtype",
                "BAM",
                "SortedByCoordinate",
                "--outFileNamePrefix",
                prefix,
                "--readFilesIn",
                str(trim1),
            ]
            if is_paired and trim2:
                star_run.append(str(trim2))
            mamba_run(star_run)
            shutil.move(f"{prefix}Aligned.sortedByCoord.out.bam", bam_path)
            mamba_run(["samtools", "index", str(bam_path)])
        else:
            h_cmd = [
                "hisat2",
                "-p",
                str(args.threads),
                "--summary-file",
                str(hisat_log),
                "-x",
                str(hisat_prefix),
                "-S",
                str(sam_path),
            ]
            if is_paired and trim2:
                h_cmd.extend(["-1", str(trim1), "-2", str(trim2)])
            else:
                h_cmd.extend(["-U", str(trim1)])
            mamba_run(h_cmd)
            mamba_run(
                ["samtools", "sort", "-@", str(args.threads), "-o", str(bam_path), str(sam_path)]
            )
            mamba_run(["samtools", "index", str(bam_path)])
            sam_path.unlink(missing_ok=True)

        sorted_bams.append(bam_path)
        if args.emit_compat_logs:
            log_out = out_dir / "logs" / f"{srr}_compat-summary.log"
            write_compat_summary(srr, fastp_json, hisat_log, log_out)

    qc_dir = out_dir / "qc"
    mamba_run(["multiqc", str(qc_dir), str(out_dir / "logs"), "-o", str(qc_dir), "-f"])

    if sorted_bams and ann_path:
        matrix_path = out_dir / "counts" / "counts_matrix.tsv"
        if args.engine == "compat" or args.quantifier == "htseq":
            run_htseq(sorted_bams, sample_names, ann_path, matrix_path)
        else:
            run_featurecounts(
                sorted_bams, sample_names, ann_path, matrix_path, args.threads, is_paired_run
            )

    if args.call_variants and sorted_bams:
        vcf_path = out_dir / "variants" / "variants.vcf"
        mpileup_cmd = [
            "mamba",
            "run",
            "-n",
            "ngs",
            "bcftools",
            "mpileup",
            "-Ou",
            "-f",
            str(fna_path),
            *(str(b) for b in sorted_bams),
        ]
        mpileup = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE)
        subprocess.run(
            ["mamba", "run", "-n", "ngs", "bcftools", "call", "-mv", "-Ov", "-o", str(vcf_path)],
            stdin=mpileup.stdout,
            check=True,
        )
        if mpileup.stdout:
            mpileup.stdout.close()
        mpileup.wait()

    print(f"[OK] RNA-seq pipeline complete in {out_dir}", file=sys.stderr)
    return 0


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.list_groups:
        groups = parse_groups_xlsx(args.xlsx)
        print(format_groups_table(groups))
        return 0

    samples: list[tuple[str, str]] = []
    genome_acc = args.genome or ""
    if args.group:
        groups = parse_groups_xlsx(args.xlsx)
        key = args.group.upper()
        if key not in groups:
            print(
                f"Error: Unknown group '{args.group}'. Run --list-groups for valid IDs.",
                file=sys.stderr,
            )
            return 1
        grp = groups[key]
        genome_acc = genome_acc or grp.genome_gcf
        for cond, srr_list in grp.conditions.items():
            for srr in srr_list:
                samples.append((srr, cond))

    if args.sra:
        samples = [
            (s.strip(), f"sample_{i}") for i, s in enumerate(args.sra.split(","), 1) if s.strip()
        ]

    if not genome_acc:
        print("Error: No genome assembly accession found; pass --genome GCF_...", file=sys.stderr)
        return 1

    if args.dry_run:
        eng_info = f"{args.engine} (aligner={args.aligner}, quantifier={args.quantifier})"
        sample_info = ", ".join(f"{s}({c})" for s, c in samples)
        print(f"[DRY-RUN] Engine         : {eng_info}")
        print(f"[DRY-RUN] Reference      : {genome_acc}")
        print(f"[DRY-RUN] Samples ({len(samples)})  : {sample_info}")
        print(
            f"[DRY-RUN] Compat logs    : {args.emit_compat_logs} "
            f"| Variants: {args.call_variants}"
        )
        print(f"[DRY-RUN] Output dir     : {args.output_dir}")
        return 0

    return execute_pipeline(args, samples, genome_acc)


if __name__ == "__main__":
    raise SystemExit(main())
