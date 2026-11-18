#!/usr/bin/env python3
"""Query NCBI, UniProt, PubChem, and RCSB PDB APIs to emit Markdown reports and sequences."""

from __future__ import annotations

import argparse
import csv
import io
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import quote

import requests

HTTP_TIMEOUT_SEC = 30
MAX_ACTIVE_TARGETS = 8
MAX_PDB_DETAILS = 15
PAINS_ACTIVE_PCT_THRESHOLD = 10.0

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"
PUBCHEM_PUG_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_GRAPHQL_URL = "https://data.rcsb.org/graphql"

CAS_PATTERN = re.compile(r"^\d{2,7}-\d{2}-\d$")


@dataclass(frozen=True)
class PresetTask:
    source: str
    query: str
    database: str = "nuccore"
    label: str = ""


PRESETS: dict[str, list[PresetTask]] = {
    "w1-ex1": [
        PresetTask(
            "ncbi",
            "Escherichia coli 16S rRNA AND refseq[filter]",
            "nuccore",
            "Task 01: E. coli 16S rRNA",
        ),
        PresetTask(
            "uniprot",
            "BRCA1 AND organism_id:9606 AND reviewed:true",
            label="Task 02: BRCA1 Human",
        ),
        PresetTask("pubchem", "curcumin", label="Task 03: Curcumin"),
        PresetTask(
            "ncbi",
            "SIRT1 AND Homo sapiens[Organism] AND refseq[filter]",
            "nuccore",
            "Task 04: Human SIRT1",
        ),
        PresetTask(
            "uniprot",
            "GLP1R AND organism_id:9606 AND reviewed:true",
            label="Task 05: GLP1R Human",
        ),
        PresetTask("pubchem", "acetaminophen", label="Task 06: Acetaminophen"),
    ],
    "w1-ex2": [
        PresetTask(
            "ncbi",
            "Staphylococcus aureus 16S rRNA AND refseq[filter]",
            "nuccore",
            "S. aureus 16S rRNA",
        ),
        PresetTask("uniprot", "EGFR AND organism_id:9606 AND reviewed:true", label="Human EGFR"),
        PresetTask("pubchem", "resveratrol", label="Resveratrol"),
        PresetTask("pubchem", "orforglipron", label="Orforglipron"),
        PresetTask("pdb", "GLP1R", label="RCSB PDB GLP1R Structures"),
    ],
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract Markdown reports and sequence files from NCBI, UniProt, PubChem, and RCSB."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  # Run all Lecture 1 Exercise 1 database tasks and save Markdown + FASTA files
  db_miner.py --preset w1-ex1 -o lecture-01/ex01-db

  # Inspect PubChem bioactivity summary for curcumin on stdout
  db_miner.py -s pubchem -q curcumin --stdout

  # Count and summarize GLP1R 3D structures on RCSB PDB
  db_miner.py -s pdb -q GLP1R --stdout
""",
    )
    parser.add_argument(
        "-s",
        "--source",
        choices=("ncbi", "uniprot", "pubchem", "pdb", "all"),
        default="all",
        help="Target biological database API (default: all).",
    )
    parser.add_argument(
        "-q",
        "--query",
        type=str,
        default=None,
        help="Search term, gene symbol, compound name, or accession.",
    )
    parser.add_argument(
        "-d",
        "--database",
        type=str,
        default="nuccore",
        help="NCBI Entrez database when --source is ncbi (default: nuccore).",
    )
    parser.add_argument(
        "--preset",
        choices=tuple(PRESETS.keys()),
        default=None,
        help="Run a predefined Lecture 1 exercise batch (default: none).",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("lecture-01"),
        help=(
            "Directory for Markdown summaries and downloaded FASTA/PDB files "
            "(default: lecture-01)."
        ),
    )
    parser.add_argument(
        "--stdout",
        action="store_true",
        help="Print Markdown report output to stdout in addition to saving.",
    )
    parser.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Preview API queries and target files without network calls.",
    )
    args = parser.parse_args(argv)
    if not args.query and not args.preset:
        parser.error("Provide either -q/--query or --preset. Run 'db_miner.py --help' for usage.")
    return args


def slugify(text: str) -> str:
    cleaned = re.sub(r"[^a-zA-Z0-9._-]+", "-", text.strip().lower()).strip("-")
    return cleaned[:60] or "query"


def _get_json(url: str, params: dict[str, Any] | None = None) -> dict[str, Any]:
    resp = requests.get(url, params=params, timeout=HTTP_TIMEOUT_SEC)
    resp.raise_for_status()
    return resp.json()


def query_ncbi(query: str, database: str, out_dir: Path) -> str:
    """Delegate Entrez lookups and FASTA downloads to the user's ncbi-fetch CLI."""
    ncbi_fetch = shutil.which("ncbi-fetch")
    if not ncbi_fetch:
        return (
            f"### NCBI ({database}): `{query}`\n\n"
            "> Error: `ncbi-fetch` was not found on `$PATH`.\n"
        )

    tsv_proc = subprocess.run(
        [
            ncbi_fetch,
            "-d",
            database,
            "-q",
            query,
            "--retmax",
            "5",
            "--sort",
            "relevance",
            "--format",
            "summary-tsv",
            "--stdout",
            "--quiet",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    if tsv_proc.returncode != 0 or not tsv_proc.stdout.strip():
        err_msg = tsv_proc.stderr.strip() or "No matching records returned."
        return f"### NCBI `{database}`: `{query}`\n\n> Lookup failed: {err_msg}\n"

    rows = list(csv.DictReader(io.StringIO(tsv_proc.stdout), delimiter="\t"))
    lines = [
        f"### NCBI `{database}`: `{query}`",
        "",
        "| Accession | Title | Organism | Length | Updated |",
        "| :--- | :--- | :--- | ---: | :--- |",
    ]
    top_acc = ""
    for row in rows:
        acc = row.get("accession", "")
        if not top_acc and acc:
            top_acc = acc
        title = row.get("title", "").replace("|", "/")[:90]
        org = row.get("organism", "")
        slen = row.get("slen", "")
        upd = row.get("updatedate", "")
        lines.append(f"| `{acc}` | {title} | *{org}* | {slen} | {upd} |")

    if top_acc and database in {"nuccore", "protein"}:
        fasta_dir = out_dir / "sequences"
        fasta_dir.mkdir(parents=True, exist_ok=True)
        fasta_path = fasta_dir / f"{top_acc}.fasta"
        subprocess.run(
            [
                ncbi_fetch,
                "-d",
                database,
                "-i",
                top_acc,
                "--format",
                "fasta",
                "-o",
                str(fasta_path),
                "--quiet",
            ],
            check=False,
        )
        if fasta_path.exists():
            lines.append(f"\n- **Top FASTA saved:** `{fasta_path}`")
    lines.append("")
    return "\n".join(lines)


def query_uniprot(query: str, out_dir: Path) -> str:
    """Query UniProtKB REST API and return structured Markdown metadata + FASTA."""
    search_query = query
    if "human" in query.lower() and "organism_id:" not in query:
        cleaned = re.sub(r"\bhuman\b", "", query, flags=re.IGNORECASE).strip()
        search_query = f"({cleaned}) AND (organism_id:9606)"

    data = _get_json(
        f"{UNIPROT_BASE}/search",
        params={"query": search_query, "format": "json", "size": 5},
    )
    results = data.get("results", [])
    if not results:
        return f"### UniProtKB: `{query}`\n\n> 0 results returned for `{search_query}`.\n"

    # Prefer Reviewed (Swiss-Prot) entry if present in top hits
    entry = next(
        (r for r in results if "Swiss-Prot" in r.get("entryType", "")),
        results[0],
    )
    acc = entry.get("primaryAccession", "")
    uni_id = entry.get("uniProtkbId", "")
    sec_accs = entry.get("secondaryAccessions", [])
    org = entry.get("organism", {})
    org_name = org.get("scientificName", "")
    taxon = org.get("taxonId", "")

    pdesc = entry.get("proteinDescription", {})
    rec_name = pdesc.get("recommendedName", {}).get("fullName", {}).get("value", "N/A")
    alt_names = [
        alt.get("fullName", {}).get("value", "")
        for alt in pdesc.get("alternativeNames", [])
        if alt.get("fullName", {}).get("value")
    ]
    genes = [
        g.get("geneName", {}).get("value", "")
        for g in entry.get("genes", [])
        if g.get("geneName", {}).get("value")
    ]

    seq_info = entry.get("sequence", {})
    seq_len = seq_info.get("length", 0)
    mol_wt = seq_info.get("molWeight", 0)
    seq_aa = seq_info.get("value", "")

    func_texts = []
    isoform_count = 0
    for cmt in entry.get("comments", []):
        ctype = cmt.get("commentType", "")
        if ctype == "FUNCTION":
            for txt in cmt.get("texts", []):
                func_texts.append(txt.get("value", ""))
        elif ctype == "ALTERNATIVE PRODUCTS":
            isoform_count = len(cmt.get("isoforms", []))

    variants = sum(1 for f in entry.get("features", []) if f.get("type") == "Natural variant")
    xrefs = entry.get("uniProtKBCrossReferences", [])
    pdb_count = sum(1 for x in xrefs if x.get("database") == "PDB")
    drugbank_count = sum(1 for x in xrefs if x.get("database") == "DrugBank")
    proteome_id = next((x.get("id", "") for x in xrefs if x.get("database") == "Proteomes"), "N/A")

    fasta_dir = out_dir / "sequences"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = fasta_dir / f"{acc}.fasta"
    if seq_aa:
        header = f">{acc}|{uni_id} {rec_name} [{org_name}]"
        fasta_path.write_text(f"{header}\n{seq_aa}\n", encoding="utf-8")

    sec_preview = ", ".join(sec_accs[:8]) or "None"
    lines = [
        f"### UniProtKB: `{query}` (`{acc}` / `{uni_id}`)",
        "",
        "| Field | Value |",
        "| :--- | :--- |",
        f"| **Entry Status** | {entry.get('entryType', 'Unknown')} |",
        f"| **Recommended Name** | {rec_name} |",
        f"| **Alternative Names** | {', '.join(alt_names[:4]) or 'None'} |",
        f"| **Gene Symbol(s)** | {', '.join(genes) or 'N/A'} |",
        f"| **Organism** | *{org_name}* (TaxID: `{taxon}`) |",
        f"| **Primary / Secondary Accessions** | `{acc}` / {sec_preview} ({len(sec_accs)} total) |",
        f"| **Proteome / Length / Mass** | `{proteome_id}` / {seq_len} aa / {mol_wt:,} Da |",
        f"| **Isoforms / Natural Variants** | "
        f"{isoform_count or 1} isoform(s) / {variants} natural variant(s) |",
        f"| **Cross-References** | {pdb_count} PDB structures, {drugbank_count} DrugBank entries |",
        f"| **Saved FASTA** | `{fasta_path}` |",
    ]
    if func_texts:
        lines.append(f"\n- **Function:** {func_texts[0]}")
    lines.append("")
    return "\n".join(lines)


def _fetch_pubchem_assay_summary(cid: int) -> tuple[dict[str, int], list[dict[str, str]]]:
    url = f"{PUBCHEM_PUG_BASE}/compound/cid/{cid}/assaysummary/CSV"
    resp = requests.get(url, timeout=HTTP_TIMEOUT_SEC)
    if resp.status_code == 404:
        return {}, []
    resp.raise_for_status()

    counts: dict[str, int] = {"Active": 0, "Inactive": 0, "Inconclusive": 0, "Unspecified": 0}
    unique_aids: set[str] = set()
    active_hits: list[dict[str, str]] = []

    reader = csv.DictReader(io.StringIO(resp.text))
    for row in reader:
        aid = row.get("AID", "").strip()
        if aid:
            unique_aids.add(aid)
        outcome = row.get("Activity Outcome", "Unspecified").strip() or "Unspecified"
        counts[outcome] = counts.get(outcome, 0) + 1
        if outcome == "Active" and len(active_hits) < MAX_ACTIVE_TARGETS:
            target = row.get("Target Name", "") or row.get("Assay Name", "")
            if target:
                active_hits.append(
                    {
                        "aid": aid,
                        "target": target[:70].replace("|", "/"),
                        "potency": row.get("Activity Value [uM]", "") or "N/A",
                        "type": row.get("Activity Name", "") or "Activity",
                    }
                )
    counts["_total_rows"] = sum(v for k, v in counts.items() if not k.startswith("_"))
    counts["_unique_aids"] = len(unique_aids)
    return counts, active_hits


def query_pubchem(query: str) -> str:
    """Query PubChem PUG-REST for chemical descriptors, CAS, and bioassay statistics."""
    cid_url = f"{PUBCHEM_PUG_BASE}/compound/name/{quote(query)}/cids/JSON"
    cid_resp = requests.get(cid_url, timeout=HTTP_TIMEOUT_SEC)
    if cid_resp.status_code == 404:
        return f"### PubChem: `{query}`\n\n> Compound not found in PubChem.\n"
    cid_resp.raise_for_status()
    cids = cid_resp.json().get("IdentifierList", {}).get("CID", [])
    if not cids:
        return f"### PubChem: `{query}`\n\n> No CIDs returned.\n"
    cid = int(cids[0])

    props_url = (
        f"{PUBCHEM_PUG_BASE}/compound/cid/{cid}/property/"
        "Title,MolecularFormula,MolecularWeight,IUPACName,ConnectivitySMILES,XLogP,TPSA,"
        "HBondDonorCount,HBondAcceptorCount/JSON"
    )
    props = _get_json(props_url)["PropertyTable"]["Properties"][0]

    cas_num = "N/A"
    syn_url = f"{PUBCHEM_PUG_BASE}/compound/cid/{cid}/synonyms/JSON"
    syn_resp = requests.get(syn_url, timeout=HTTP_TIMEOUT_SEC)
    if syn_resp.ok:
        info = syn_resp.json().get("InformationList", {}).get("Information", [{}])
        syns = info[0].get("Synonym", [])
        cas_num = next((s for s in syns if CAS_PATTERN.match(s)), "N/A")

    counts, active_hits = _fetch_pubchem_assay_summary(cid)
    total_rows = counts.get("_total_rows", 0)
    active_n = counts.get("Active", 0)
    active_pct = (active_n / total_rows * 100.0) if total_rows else 0.0

    formula = props.get("MolecularFormula", "N/A")
    mw = props.get("MolecularWeight", "N/A")
    xlogp = props.get("XLogP", "N/A")
    tpsa = props.get("TPSA", "N/A")
    hbd = props.get("HBondDonorCount", "N/A")
    hba = props.get("HBondAcceptorCount", "N/A")
    lines = [
        f"### PubChem: `{props.get('Title', query)}` (`CID {cid}`)",
        "",
        "| Descriptor | Value |",
        "| :--- | :--- |",
        f"| **Formula / MW** | `{formula}` / {mw} g/mol |",
        f"| **CAS / IUPAC** | `{cas_num}` / {props.get('IUPACName', 'N/A')} |",
        f"| **Lipophilicity / TPSA** | `XLogP = {xlogp}`, `TPSA = {tpsa} A^2` |",
        f"| **HBD / HBA** | {hbd} / {hba} |",
        f"| **SMILES** | `{props.get('ConnectivitySMILES', 'N/A')}` |",
    ]

    if total_rows:
        pains_flag = (
            f" *(High promiscuity > {PAINS_ACTIVE_PCT_THRESHOLD:.0f}% — potential PAINS)*"
            if active_pct > PAINS_ACTIVE_PCT_THRESHOLD
            else ""
        )
        unique_aids = counts.get("_unique_aids", 0)
        lines.extend(
            [
                "",
                f"**BioAssay Summary:** {total_rows:,} records across {unique_aids:,} "
                f"unique BioAssays — Active: **{active_n:,} ({active_pct:.1f}%)**{pains_flag}, "
                f"Inactive: {counts.get('Inactive', 0):,}, "
                f"Inconclusive: {counts.get('Inconclusive', 0):,}, "
                f"Unspecified: {counts.get('Unspecified', 0):,}.",
            ]
        )
        if active_hits:
            lines.extend(
                [
                    "",
                    "| AID | Target / Assay | Metric | Value (uM) |",
                    "| :--- | :--- | :--- | ---: |",
                ]
            )
            for hit in active_hits:
                lines.append(
                    f"| `AID {hit['aid']}` | {hit['target']} | "
                    f"{hit['type']} | {hit['potency']} |"
                )
    lines.append("")
    return "\n".join(lines)


def query_rcsb_pdb(query: str) -> str:
    """Query RCSB PDB Search v2 + GraphQL API for experimental structure breakdown."""
    payload = {
        "query": {"type": "terminal", "service": "full_text", "parameters": {"value": query}},
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": MAX_PDB_DETAILS},
            "return_counts": True,
        },
    }
    resp = requests.post(RCSB_SEARCH_URL, json=payload, timeout=HTTP_TIMEOUT_SEC)
    if resp.status_code == 204:
        return f"### RCSB PDB: `{query}`\n\n> 0 experimental structures found.\n"
    resp.raise_for_status()
    data = resp.json()
    total_count = data.get("total_count", 0)
    pdb_ids = [hit["identifier"] for hit in data.get("result_set", [])]
    if not pdb_ids:
        return (
            f"### RCSB PDB: `{query}`\n\n"
            f"> Total count: {total_count}, no identifiers returned.\n"
        )

    gql = """
    query($ids: [String!]!) {
      entries(entry_ids: $ids) {
        rcsb_id
        struct { title }
        exptl { method }
        rcsb_entry_info { resolution_combined experimental_method }
      }
    }
    """
    gql_resp = requests.post(
        RCSB_GRAPHQL_URL,
        json={"query": gql, "variables": {"ids": pdb_ids}},
        timeout=HTTP_TIMEOUT_SEC,
    )
    gql_resp.raise_for_status()
    entries = gql_resp.json().get("data", {}).get("entries", [])

    method_counts: dict[str, int] = {}
    lines = [
        f"### RCSB PDB: `{query}` ({total_count} total structures)",
        "",
        "| PDB ID | Method | Resolution (A) | Title |",
        "| :--- | :--- | ---: | :--- |",
    ]
    for entry in entries:
        info = entry.get("rcsb_entry_info", {})
        exptl = entry.get("exptl") or [{}]
        method = info.get("experimental_method") or exptl[0].get("method", "Unknown")
        method_counts[method] = method_counts.get(method, 0) + 1
        res_list = info.get("resolution_combined") or []
        res_str = f"{res_list[0]:.2f}" if res_list else "N/A"
        title = (entry.get("struct", {}).get("title") or "").replace("|", "/")[:80]
        lines.append(f"| `{entry.get('rcsb_id', '')}` | {method} | {res_str} | {title} |")

    breakdown = ", ".join(f"{k}: {v}" for k, v in sorted(method_counts.items()))
    lines.insert(2, f"**Top {len(entries)} Hits Method Breakdown:** {breakdown}\n")
    lines.append("")
    return "\n".join(lines)


def run_task(task: PresetTask, out_dir: Path) -> str:
    """Dispatch a single database query task to its REST/CLI handler."""
    header = f"## {task.label}\n\n" if task.label else ""
    if task.source == "ncbi":
        return header + query_ncbi(task.query, task.database, out_dir)
    if task.source == "uniprot":
        return header + query_uniprot(task.query, out_dir)
    if task.source == "pubchem":
        return header + query_pubchem(task.query)
    if task.source == "pdb":
        return header + query_rcsb_pdb(task.query)
    raise ValueError(f"Unsupported source: {task.source}")


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    tasks: list[PresetTask] = []
    if args.preset:
        tasks.extend(PRESETS[args.preset])
        report_stem = f"db-summary_{args.preset}"
    else:
        sources = ("ncbi", "uniprot", "pubchem", "pdb") if args.source == "all" else (args.source,)
        for src in sources:
            tasks.append(PresetTask(src, args.query, args.database))
        report_stem = f"db-summary_{slugify(args.query)}"

    if args.dry_run:
        print(f"[DRY-RUN] Output directory: {args.output_dir}")
        for idx, task in enumerate(tasks, start=1):
            print(f"[DRY-RUN] {idx}. source={task.source} db={task.database} query={task.query!r}")
        return 0

    args.output_dir.mkdir(parents=True, exist_ok=True)
    sections = [f"# Database Query Report: `{args.preset or args.query}`\n"]
    for task in tasks:
        sections.append(run_task(task, args.output_dir))

    report_md = "\n".join(sections)
    report_path = args.output_dir / f"{report_stem}.md"
    report_path.write_text(report_md, encoding="utf-8")
    print(f"[OK] Saved Markdown summary to {report_path}", file=sys.stderr)
    if args.stdout:
        print(report_md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
