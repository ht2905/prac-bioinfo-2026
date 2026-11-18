#!/usr/bin/env python3
"""Run NCBI BLAST, MAFFT alignment, FastTree inference, and ggtree plotting for Lecture 1."""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from xml.etree import ElementTree

import requests

BLAST_API_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
EUTILS_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
DEFAULT_POLL_SEC = 10
MAX_POLL_ATTEMPTS = 60


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fetch homologs via BLAST, align with MAFFT, build a FastTree phylogeny, "
            "and export MEGA/NCBI-ready files."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  # Build Arabidopsis Cytochrome C phylogeny (top 100 nr hits)
  phylo_pipeline.py -a AAB72175.1 -p blastp -m 100 -o lecture-01/phylogeny-cytc

  # Align an existing multi-FASTA file and infer a Newick + PNG tree
  phylo_pipeline.py -i sequences.fasta -p blastn -o lecture-01/phylogeny-16s
""",
    )
    parser.add_argument(
        "-i",
        "--input-fasta",
        type=Path,
        default=None,
        help="Input FASTA (single query to BLAST, or >=4 sequences with --skip-blast).",
    )
    parser.add_argument(
        "-a",
        "--accession",
        type=str,
        default=None,
        help="NCBI accession ID to fetch as query sequence (e.g. AAB72175.1 or NR_024570.1).",
    )
    parser.add_argument(
        "-p",
        "--program",
        choices=("blastn", "blastp"),
        default="blastp",
        help="BLAST search program (default: blastp).",
    )
    parser.add_argument(
        "--database",
        type=str,
        default=None,
        help="NCBI BLAST database (default: nr for blastp, nt for blastn).",
    )
    parser.add_argument(
        "-m",
        "--max-targets",
        type=int,
        default=100,
        help="Maximum top BLAST homologs to retrieve (default: 100).",
    )
    parser.add_argument(
        "--aligner",
        choices=("mafft", "clustalo"),
        default="mafft",
        help="Multiple sequence aligner in ngs env (default: mafft).",
    )
    parser.add_argument(
        "--skip-blast",
        action="store_true",
        help="Skip BLAST and align/build tree directly from --input-fasta.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=8,
        help="CPU threads for MAFFT/FastTree (default: 8).",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("lecture-01/phylogeny"),
        help=(
            "Directory for FASTA, aligned FASTA, Newick tree, and PNG plot "
            "(default: lecture-01/phylogeny)."
        ),
    )
    parser.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Preview pipeline commands without executing network searches or tools.",
    )
    args = parser.parse_args(argv)
    if not args.input_fasta and not args.accession:
        parser.error(
            "Provide either -i/--input-fasta or -a/--accession. "
            "Run 'phylo_pipeline.py --help' for usage."
        )
    return args


def fetch_query_fasta(accession: str, program: str, dest: Path) -> Path:
    db = "protein" if program == "blastp" else "nuccore"
    ncbi_fetch = shutil.which("ncbi-fetch")
    if ncbi_fetch:
        proc = subprocess.run(
            [
                ncbi_fetch,
                "-d",
                db,
                "-i",
                accession,
                "--format",
                "fasta",
                "-o",
                str(dest),
                "--quiet",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode == 0 and dest.exists() and dest.stat().st_size > 0:
            return dest

    resp = requests.get(
        EUTILS_EFETCH_URL,
        params={"db": db, "id": accession, "rettype": "fasta", "retmode": "text"},
        timeout=30,
    )
    resp.raise_for_status()
    dest.write_text(resp.text, encoding="utf-8")
    return dest


def run_remote_blast(query_seq: str, program: str, database: str, max_targets: int) -> list[str]:
    """Submit a query to NCBI QBLAST REST API and return top hit accession/GI IDs."""
    put_resp = requests.post(
        BLAST_API_URL,
        data={
            "CMD": "Put",
            "PROGRAM": program,
            "DATABASE": database,
            "QUERY": query_seq,
            "HITLIST_SIZE": str(max_targets),
            "FORMAT_TYPE": "XML",
        },
        timeout=60,
    )
    put_resp.raise_for_status()
    rid_match = re.search(r"RID = (\S+)", put_resp.text)
    if not rid_match:
        raise RuntimeError("Failed to parse BLAST Request ID (RID) from NCBI response.")
    rid = rid_match.group(1)
    print(
        f"[INFO] Submitted NCBI BLAST job RID={rid} ({program}/{database}, max={max_targets})...",
        file=sys.stderr,
    )

    for _ in range(MAX_POLL_ATTEMPTS):
        time.sleep(DEFAULT_POLL_SEC)
        status_resp = requests.get(
            BLAST_API_URL,
            params={"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid},
            timeout=30,
        )
        status_resp.raise_for_status()
        if "Status=READY" in status_resp.text:
            break
        if "Status=FAILED" in status_resp.text or "Status=UNKNOWN" in status_resp.text:
            raise RuntimeError(f"NCBI BLAST job {rid} failed or expired.")
    else:
        timeout_sec = MAX_POLL_ATTEMPTS * DEFAULT_POLL_SEC
        raise TimeoutError(f"NCBI BLAST job {rid} timed out after {timeout_sec}s.")

    xml_resp = requests.get(
        BLAST_API_URL,
        params={"CMD": "Get", "FORMAT_TYPE": "XML", "RID": rid},
        timeout=60,
    )
    xml_resp.raise_for_status()
    root = ElementTree.fromstring(xml_resp.text)
    accessions: list[str] = []
    seen: set[str] = set()
    for hit in root.findall(".//Hit"):
        acc = (hit.findtext("Hit_accession") or hit.findtext("Hit_id") or "").strip()
        if acc and acc not in seen:
            seen.add(acc)
            accessions.append(acc)
    return accessions[:max_targets]


def batch_fetch_sequences(accessions: list[str], program: str, dest: Path) -> Path:
    """Download FASTA records for hit accessions via E-utilities POST."""
    db = "protein" if program == "blastp" else "nuccore"
    chunks: list[str] = []
    for i in range(0, len(accessions), 100):
        batch = ",".join(accessions[i : i + 100])
        resp = requests.post(
            EUTILS_EFETCH_URL,
            data={"db": db, "id": batch, "rettype": "fasta", "retmode": "text"},
            timeout=60,
        )
        resp.raise_for_status()
        chunks.append(resp.text.strip())
        time.sleep(0.35)
    dest.write_text("\n".join(c for c in chunks if c) + "\n", encoding="utf-8")
    return dest


def render_tree_plot(newick_path: Path, png_path: Path, title: str) -> None:
    """Render phylogenetic tree PNG using R (ggtree or ape fallback)."""
    r_code = f"""
    suppressPackageStartupMessages({{
      has_ggtree <- requireNamespace("ggtree", quietly = TRUE)
      has_ape <- requireNamespace("ape", quietly = TRUE)
    }})
    if (!has_ape) quit(status = 0)
    tree <- ape::read.tree("{newick_path}")
    if (has_ggtree) {{
      p <- ggtree::ggtree(tree) +
        ggtree::geom_tiplab(size = 2.2) +
        ggplot2::ggtitle("{title}") +
        ggplot2::theme_tree2()
      h_in <- max(7, length(tree$tip.label) * 0.14)
      ggplot2::ggsave("{png_path}", p, width = 11, height = h_in, limitsize = FALSE, dpi = 300)
    }} else {{
      png("{png_path}", width = 1800, height = max(1200, length(tree$tip.label) * 26), res = 220)
      plot(tree, cex = 0.55, main = "{title}")
      dev.off()
    }}
    """
    subprocess.run(["mamba", "run", "-n", "ngs", "Rscript", "-e", r_code], check=False)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    database = args.database or ("nr" if args.program == "blastp" else "nt")
    out_dir = args.output_dir

    query_fasta = args.input_fasta or (out_dir / f"{args.accession}_query.fasta")
    hits_fasta = out_dir / "homologs_raw.fasta"
    aligned_fasta = out_dir / "homologs_aligned.fasta"
    tree_nwk = out_dir / "phylogeny.nwk"
    tree_png = out_dir / "phylogeny_tree.png"

    if args.dry_run:
        ft_flags = " -wag" if args.program == "blastp" else " -nt -gtr"
        print(f"[DRY-RUN] Program     : {args.program} (db={database}, max={args.max_targets})")
        print(f"[DRY-RUN] Query       : {args.input_fasta or args.accession}")
        print(f"[DRY-RUN] Aligner     : {args.aligner} (threads={args.threads}) -> {aligned_fasta}")
        print(f"[DRY-RUN] Tree build  : FastTree{ft_flags} -> {tree_nwk}")
        print(f"[DRY-RUN] Plot        : {tree_png}")
        return 0

    out_dir.mkdir(parents=True, exist_ok=True)
    if args.accession and not args.input_fasta:
        fetch_query_fasta(args.accession, args.program, query_fasta)

    if args.skip_blast:
        shutil.copy2(query_fasta, hits_fasta)
    else:
        query_text = query_fasta.read_text(encoding="utf-8")
        hit_ids = run_remote_blast(query_text, args.program, database, args.max_targets)
        if not hit_ids:
            print("Error: BLAST returned 0 homologs.", file=sys.stderr)
            return 1
        print(f"[INFO] Downloading {len(hit_ids)} hit sequences...", file=sys.stderr)
        batch_fetch_sequences(hit_ids, args.program, hits_fasta)

    align_cmd = [
        "mamba",
        "run",
        "-n",
        "ngs",
        "mafft",
        "--auto",
        "--thread",
        str(args.threads),
        str(hits_fasta),
    ]
    print(f"[INFO] Aligning sequences -> {aligned_fasta}", file=sys.stderr)
    with aligned_fasta.open("w", encoding="utf-8") as handle:
        subprocess.run(align_cmd, stdout=handle, check=True)

    ft_args = ["mamba", "run", "-n", "ngs", "FastTree"]
    if args.program == "blastn":
        ft_args.extend(["-nt", "-gtr"])
    ft_args.append(str(aligned_fasta))
    print(f"[INFO] Inferring phylogenetic tree -> {tree_nwk}", file=sys.stderr)
    with tree_nwk.open("w", encoding="utf-8") as handle:
        subprocess.run(ft_args, stdout=handle, check=True)

    title = f"Phylogenetic Tree ({args.accession or query_fasta.stem})"
    render_tree_plot(tree_nwk, tree_png, title)
    print(
        f"[OK] Ready for MEGA 12 / NCBI inspection:\n"
        f"     - Aligned FASTA: {aligned_fasta}\n"
        f"     - Newick Tree  : {tree_nwk}\n"
        f"     - Rendered PNG : {tree_png}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
