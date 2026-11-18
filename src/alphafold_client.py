#!/usr/bin/env python3
"""Fetch EBI AlphaFold DB structures/PAE, compute pLDDT stats, and emit AF3 Server JSONs."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import matplotlib
import numpy as np
import requests
from Bio.PDB import PDBParser

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

AFDB_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot}"
UNIPROT_JSON_URL = "https://rest.uniprot.org/uniprotkb/{uniprot}.json"

LECTURE4_PROTEINS: tuple[tuple[str, str, int], ...] = (
    ("P00698", "Lysozyme", 1),
    ("P01308", "Insulin", 1),
    ("P00338", "Lactate_dehydrogenase_A", 4),
    ("P42212", "Green_fluorescent_protein", 1),
    ("P00004", "Cytochrome_c", 1),
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Automate Lecture 4 AlphaFold analysis via EBI AFDB API and build "
            "AlphaFold 3 Server batch JSON."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  # Run all 5 Lecture 4 proteins (Lysozyme, Insulin, LDH, GFP, Cytochrome c)
  alphafold_client.py --preset lecture4-af -o lecture-04/alphafold

  # Analyze a specific UniProt accession
  alphafold_client.py -u P00698 -o .scratch/af-lysozyme
""",
    )
    parser.add_argument(
        "-u",
        "--uniprot",
        type=str,
        default=None,
        help="Comma-separated UniProt accessions (e.g. P00698,P01308).",
    )
    parser.add_argument(
        "--preset",
        choices=("lecture4-af",),
        default=None,
        help="Run all 5 Lecture 4 Part B target proteins (default: none).",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("lecture-04/alphafold"),
        help=(
            "Directory for PDBs, PAE JSONs, plots, Markdown summary, and AF3 JSON "
            "(default: lecture-04/alphafold)."
        ),
    )
    parser.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Preview target UniProt IDs and output artifacts without network calls.",
    )
    args = parser.parse_args(argv)
    if not args.uniprot and not args.preset:
        parser.error(
            "Provide either -u/--uniprot or --preset. Run 'alphafold_client.py --help' for usage."
        )
    return args


def extract_plddt_from_pdb(pdb_path: Path) -> list[tuple[int, str, float]]:
    """Return ordered list of (resseq, resname, pLDDT) from AlphaFold PDB B-factors."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    model = next(structure.get_models())
    residues: list[tuple[int, str, float]] = []
    for chain in model:
        for res in chain:
            if res.id[0] == " " and "CA" in res:
                residues.append((int(res.id[1]), res.get_resname(), float(res["CA"].get_bfactor())))
        break
    return residues


def find_low_confidence_spans(
    scores: list[tuple[int, str, float]], threshold: float = 70.0
) -> list[tuple[int, int]]:
    """Group contiguous residue ranges with pLDDT < threshold."""
    spans: list[tuple[int, int]] = []
    start: int | None = None
    prev: int | None = None
    for resnum, _resn, score in scores:
        if score < threshold:
            if start is None:
                start = resnum
            prev = resnum
        elif start is not None and prev is not None:
            spans.append((start, prev))
            start = None
            prev = None
    if start is not None and prev is not None:
        spans.append((start, prev))
    return spans


def summarize_uniprot_features(uniprot_acc: str) -> dict[str, object]:
    """Fetch canonical sequence, protein name, and key structural/functional annotations."""
    resp = requests.get(UNIPROT_JSON_URL.format(uniprot=uniprot_acc), timeout=30)
    resp.raise_for_status()
    data = resp.json()
    pname = (
        data.get("proteinDescription", {})
        .get("recommendedName", {})
        .get("fullName", {})
        .get("value", uniprot_acc)
    )
    seq = data.get("sequence", {}).get("value", "")
    active_sites: list[str] = []
    binding_sites: list[str] = []
    chains: list[tuple[str, int, int]] = []
    disulfides = 0
    for feat in data.get("features", []):
        ftype = feat.get("type", "")
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        if ftype == "Active site" and start:
            desc = feat.get("description", "").strip()
            active_sites.append(f"{start} ({desc})" if desc else str(start))
        elif ftype == "Binding site" and start:
            binding_sites.append(f"{start}-{end or start}")
        elif ftype == "Disulfide bond":
            disulfides += 1
        elif ftype == "Chain" and start and end:
            chains.append((feat.get("description", "Chain"), int(start), int(end)))
    return {
        "name": pname,
        "sequence": seq,
        "active_sites": active_sites[:6],
        "binding_sites": binding_sites[:6],
        "chains": chains,
        "disulfides": disulfides,
    }


def render_plddt_and_pae_plot(
    acc: str,
    scores: list[tuple[int, str, float]],
    pae_json_path: Path | None,
    out_png: Path,
) -> None:
    """Plot per-residue pLDDT line profile and 2D PAE heatmap."""
    x_pos = [s[0] for s in scores]
    y_val = [s[2] for s in scores]

    has_pae = pae_json_path is not None and pae_json_path.exists()
    fig, axes = plt.subplots(1, 2 if has_pae else 1, figsize=(12 if has_pae else 7, 4.5), dpi=300)
    ax1 = axes[0] if has_pae else axes

    ax1.plot(x_pos, y_val, color="#0053D6", linewidth=1.6)
    ax1.axhspan(90, 100, color="#0053D6", alpha=0.08, label="Very high (>90)")
    ax1.axhspan(70, 90, color="#65CBF3", alpha=0.10, label="Confident (70-90)")
    ax1.axhspan(50, 70, color="#FFDB13", alpha=0.12, label="Low (50-70)")
    ax1.axhspan(0, 50, color="#FF7D45", alpha=0.12, label="Very low (<50)")
    ax1.set_ylim(0, 102)
    ax1.set_xlabel("Residue Position")
    ax1.set_ylabel("pLDDT Score")
    ax1.set_title(f"{acc} Per-Residue pLDDT")
    ax1.legend(loc="lower right", fontsize=7)

    if has_pae and pae_json_path:
        pae_raw = json.loads(pae_json_path.read_text(encoding="utf-8"))
        pae_mat = np.array(pae_raw[0]["predicted_aligned_error"], dtype=float)
        ax2 = axes[1]
        im = ax2.imshow(pae_mat, cmap="Greens_r", vmin=0, vmax=30, origin="upper")
        ax2.set_title(f"{acc} Predicted Aligned Error (PAE)")
        ax2.set_xlabel("Scored Residue")
        ax2.set_ylabel("Aligned Residue")
        fig.colorbar(im, ax=ax2, label="Expected position error (A)")

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def render_pymol_plddt(pdb_path: Path, out_png: Path) -> None:
    """Render 3D structure colored by AlphaFold pLDDT confidence bands in headless PyMOL."""
    pml = out_png.with_suffix(".pml")
    lines = [
        f"load {pdb_path.resolve()}, af_model",
        "bg_color white",
        "set ray_opaque_background, 1",
        "hide everything",
        "show cartoon, af_model",
        "set_color af_vhigh, [0.00, 0.33, 0.84]",
        "set_color af_conf,  [0.40, 0.80, 0.95]",
        "set_color af_low,   [1.00, 0.86, 0.07]",
        "set_color af_vlow,  [1.00, 0.49, 0.27]",
        "color af_vlow,  af_model and b < 50",
        "color af_low,   af_model and b > 49.99 and b < 70",
        "color af_conf,  af_model and b > 69.99 and b < 90",
        "color af_vhigh, af_model and b > 89.99",
        "orient af_model",
        "ray 1400, 1000",
        f"png {out_png.resolve()}, dpi=300",
    ]
    pml.write_text("\n".join(lines) + "\n", encoding="utf-8")
    subprocess.run(["mamba", "run", "-n", "sbvs", "pymol", "-cq", str(pml)], check=True)


def build_af3_server_job(
    acc: str, label: str, oligo_count: int, up_info: dict[str, object]
) -> dict[str, object]:
    """Construct an AlphaFold 3 Server JSON job entry (mature chains for Insulin P01308)."""
    seq = str(up_info["sequence"])
    if acc == "P01308":
        # Mature human insulin: B-chain (25-54) and A-chain (90-110) cleaved from precursor
        chains = [
            {"proteinChain": {"sequence": seq[24:54], "count": 1}},
            {"proteinChain": {"sequence": seq[89:110], "count": 1}},
        ]
        return {"name": f"{acc}_{label}_mature_dimer", "modelSeeds": [], "sequences": chains}

    return {
        "name": f"{acc}_{label}_{oligo_count}mer",
        "modelSeeds": [],
        "sequences": [{"proteinChain": {"sequence": seq, "count": oligo_count}}],
    }


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    targets: list[tuple[str, str, int]] = []
    if args.preset == "lecture4-af":
        targets.extend(LECTURE4_PROTEINS)
    if args.uniprot:
        for acc in [a.strip().upper() for a in args.uniprot.split(",") if a.strip()]:
            if not any(t[0] == acc for t in targets):
                targets.append((acc, acc, 1))

    if args.dry_run:
        target_str = ", ".join(f"{a} ({c}x)" for a, _, c in targets)
        print(f"[DRY-RUN] Targets ({len(targets)}): {target_str}")
        print(f"[DRY-RUN] Output dir: {args.output_dir}")
        return 0

    out_dir = args.output_dir
    for sub in ("structures", "pae", "figures", "tables", "jobs"):
        (out_dir / sub).mkdir(parents=True, exist_ok=True)

    summary_rows = [
        "# Lecture 4B — AlphaFold Structure Prediction & Confidence Summary",
        "",
        "| UniProt | Protein | Length (aa) | Min / Mean / Max pLDDT "
        "| Low-Confidence Regions (<70) | Functional / Structural Highlights |",
        "| :--- | :--- | ---: | :--- | :--- | :--- |",
    ]
    af3_jobs: list[dict[str, object]] = []

    for acc, label, oligo_count in targets:
        print(f"[INFO] Fetching AlphaFold DB + UniProt metadata for {acc}...", file=sys.stderr)
        af_resp = requests.get(AFDB_API_URL.format(uniprot=acc), timeout=30)
        af_resp.raise_for_status()
        entry = af_resp.json()[0]

        pdb_path = out_dir / "structures" / f"{acc}_afdb.pdb"
        pdb_path.write_text(requests.get(entry["pdbUrl"], timeout=45).text, encoding="utf-8")

        pae_path: Path | None = None
        if entry.get("paeDocUrl"):
            pae_path = out_dir / "pae" / f"{acc}_pae.json"
            pae_path.write_text(requests.get(entry["paeDocUrl"], timeout=45).text, encoding="utf-8")

        up_info = summarize_uniprot_features(acc)
        af3_jobs.append(build_af3_server_job(acc, label, oligo_count, up_info))

        scores = extract_plddt_from_pdb(pdb_path)
        vals = [s[2] for s in scores]
        low_spans = find_low_confidence_spans(scores, threshold=70.0)
        span_str = ", ".join(f"{a}-{b}" for a, b in low_spans[:6]) or "None"

        fig_pae = out_dir / "figures" / f"{acc}_plddt-pae.png"
        render_plddt_and_pae_plot(acc, scores, pae_path, fig_pae)
        render_pymol_plddt(pdb_path, out_dir / "figures" / f"{acc}_3d-plddt.png")

        highlights: list[str] = []
        if up_info["active_sites"]:
            highlights.append("Active: " + "; ".join(up_info["active_sites"]))  # type: ignore[arg-type]
        if up_info["binding_sites"]:
            highlights.append("Binding: " + ", ".join(up_info["binding_sites"]))  # type: ignore[arg-type]
        if up_info["disulfides"]:
            highlights.append(f"{up_info['disulfides']} S-S bonds")
        hl_text = "; ".join(highlights) or "See PDB structure"

        summary_rows.append(
            f"| `{acc}` | {up_info['name']} | {len(scores)} | "
            f"{min(vals):.1f} / {sum(vals) / len(vals):.1f} / {max(vals):.1f} | "
            f"{span_str} | {hl_text} |"
        )

    summary_path = out_dir / "tables" / "alphafold_summary.md"
    summary_path.write_text("\n".join(summary_rows) + "\n", encoding="utf-8")

    batch_path = out_dir / "jobs" / "af3_server_batch.json"
    batch_path.write_text(json.dumps(af3_jobs, indent=2) + "\n", encoding="utf-8")
    print(f"[OK] Saved AlphaFold summary -> {summary_path}", file=sys.stderr)
    print(f"[OK] Saved AF3 Server JSON -> {batch_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
