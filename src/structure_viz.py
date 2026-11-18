#!/usr/bin/env python3
"""Render Lecture 3 PyMOL 3D protein visualizations, offline APBS maps, H-bonds, and alignments."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import requests
from Bio.PDB import PDBParser

RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
PYMOL_HELPERS_DIR = Path(__file__).resolve().parent / "pymol"

PRESETS: dict[str, dict[str, tuple[str, ...]]] = {
    "lecture3-tutorials": {"levels": ("1LYZ", "1A3N"), "props": ("1HEW", "1H4G")},
    "lecture3-levels": {"levels": ("1A3N", "1GZX", "3NIR", "2MNR", "1QTI"), "props": ()},
    "lecture3-props": {"levels": (), "props": ("1A2C", "1HCL", "1MOE", "4HHB")},
    "all": {
        "levels": ("1A3N", "1GZX", "3NIR", "2MNR", "1QTI"),
        "props": ("1A2C", "1HCL", "1MOE", "4HHB"),
    },
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run headless PyMOL rendering and local PDB2PQR/APBS electrostatics for "
            "Lecture 3 structural exercises."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  # Render all 5 Lecture 3 Section 4.5 structure-level exercises (1A3N, 1GZX, 3NIR, 2MNR, 1QTI)
  structure_viz.py --preset lecture3-levels -o lecture-03

  # Render hydrophobicity, offline APBS electrostatics, H-bonds, and alignment for 1HEW
  structure_viz.py -p 1HEW -m properties -o .scratch/pymol-test
""",
    )
    parser.add_argument(
        "-p",
        "--pdb",
        type=str,
        default=None,
        help="Comma-separated PDB IDs (e.g. 1HEW,1A2C,4HHB).",
    )
    parser.add_argument(
        "--preset",
        choices=tuple(PRESETS.keys()),
        default=None,
        help="Predefined Lecture 3 PDB preset batch (default: none).",
    )
    parser.add_argument(
        "-m",
        "--mode",
        choices=("levels", "properties", "both"),
        default="both",
        help="Visualization tasks to render (default: both).",
    )
    parser.add_argument(
        "--align-chains",
        type=str,
        default=None,
        help="Comma-separated chain pair to align, e.g. A,B (default: first two distinct chains).",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("lecture-03"),
        help="Output dir for structures, figures, sessions, and tables (default: lecture-03).",
    )
    parser.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Preview PDB downloads, APBS runs, and PyMOL render targets without executing.",
    )
    args = parser.parse_args(argv)
    if not args.pdb and not args.preset:
        parser.error(
            "Provide either -p/--pdb or --preset. Run 'structure_viz.py --help' for usage."
        )
    return args


def download_pdb(pdb_id: str, dest_dir: Path) -> Path:
    dest_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = dest_dir / f"{pdb_id.upper()}.pdb"
    if not pdb_path.exists():
        resp = requests.get(RCSB_PDB_URL.format(pdb_id=pdb_id.upper()), timeout=45)
        resp.raise_for_status()
        pdb_path.write_text(resp.text, encoding="utf-8")
    return pdb_path


def inspect_structure_stats(pdb_path: Path) -> dict[str, object]:
    """Extract chains, residue counts, cysteines, and disulfide bonds via BioPython."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    model = next(structure.get_models())
    chains: list[str] = []
    aa_count = 0
    cys_count = 0
    for chain in model:
        poly_res = [r for r in chain if r.id[0] == " "]
        if poly_res:
            chains.append(chain.id)
            aa_count += len(poly_res)
            cys_count += sum(1 for r in poly_res if r.get_resname() == "CYS")

    ssbond_count = sum(
        1 for line in pdb_path.read_text(encoding="utf-8").splitlines() if line.startswith("SSBOND")
    )
    return {
        "pdb_id": pdb_path.stem,
        "chains": ",".join(chains) or "None",
        "num_chains": len(chains),
        "residues": aa_count,
        "cys_residues": cys_count,
        "disulfide_bonds": ssbond_count,
    }


def strip_for_pdb2pqr(raw_pdb: Path, clean_pdb: Path) -> None:
    """Keep standard ATOM protein records so local pdb2pqr + apbs run reliably offline."""
    lines = [
        line
        for line in raw_pdb.read_text(encoding="utf-8").splitlines()
        if line.startswith(("ATOM  ", "TER", "END"))
    ]
    clean_pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_offline_apbs(pdb_path: Path, work_dir: Path) -> tuple[Path | None, Path | None]:
    """Run local pdb2pqr + apbs inside sbvs environment to produce .pqr and .dx files."""
    work_dir.mkdir(parents=True, exist_ok=True)
    clean_pdb = work_dir / f"{pdb_path.stem}_prot.pdb"
    pqr_path = work_dir / f"{pdb_path.stem}.pqr"
    apbs_in = work_dir / f"{pdb_path.stem}.in"
    strip_for_pdb2pqr(pdb_path, clean_pdb)

    p2p = subprocess.run(
        [
            "mamba",
            "run",
            "-n",
            "sbvs",
            "pdb2pqr",
            "--ff=AMBER",
            "--with-ph=7.4",
            f"--apbs-input={apbs_in.name}",
            clean_pdb.name,
            pqr_path.name,
        ],
        cwd=work_dir,
        capture_output=True,
        text=True,
        check=False,
    )
    if p2p.returncode != 0 or not apbs_in.exists():
        err_tail = p2p.stderr.splitlines()[-1:]
        print(f"[WARN] pdb2pqr skipped for {pdb_path.stem}: {err_tail}", file=sys.stderr)
        return None, None

    apbs_proc = subprocess.run(
        ["mamba", "run", "-n", "sbvs", "apbs", apbs_in.name],
        cwd=work_dir,
        capture_output=True,
        text=True,
        check=False,
    )
    dx_files = sorted(work_dir.glob("*.dx"))
    if apbs_proc.returncode != 0 or not dx_files:
        print(
            f"[WARN] apbs skipped for {pdb_path.stem}; continuing with other renders.",
            file=sys.stderr,
        )
        return pqr_path if pqr_path.exists() else None, None
    return pqr_path, dx_files[0]


def run_pymol_levels(
    pdb_path: Path, stats: dict[str, object], fig_dir: Path, sess_dir: Path
) -> None:
    """Render secondary, tertiary CYS disulfides, and quaternary chains in headless PyMOL."""
    pdb_id = pdb_path.stem
    pml = sess_dir / f"{pdb_id}_levels.pml"
    chains = str(stats["chains"]).split(",")
    palette = ("red", "marine", "forest", "gold", "magenta", "orange", "cyan", "salmon")

    lines = [
        f"load {pdb_path.resolve()}, prot",
        "bg_color white",
        "set ray_opaque_background, 1",
        "set cartoon_fancy_helices, 1",
        "hide everything",
        "show cartoon, prot",
        # Secondary structure coloring (Lecture 3 §4.2)
        "color red, ss H",
        "color yellow, ss S",
        "color cyan, ss L+",
        "orient prot",
        "ray 1400, 1000",
        f"png {(fig_dir / f'{pdb_id}_01-secondary.png').resolve()}, dpi=300",
        # Tertiary disulfide bonds (Lecture 3 §4.3)
        "color gray80, prot",
        "select cys_res, resn CYS",
        "show sticks, cys_res",
        "show spheres, cys_res and name SG",
        "set sphere_scale, 0.45",
        "color orange, cys_res",
        "ray 1400, 1000",
        f"png {(fig_dir / f'{pdb_id}_02-tertiary-disulfide.png').resolve()}, dpi=300",
        # Quaternary organization (Lecture 3 §4.4)
        "hide sticks",
        "hide spheres",
    ]
    for idx, ch in enumerate(chains):
        if ch and ch != "None":
            lines.append(f"color {palette[idx % len(palette)]}, chain {ch}")
    lines.extend(
        [
            "orient prot",
            "ray 1400, 1000",
            f"png {(fig_dir / f'{pdb_id}_03-quaternary-chains.png').resolve()}, dpi=300",
            f"save {(sess_dir / f'{pdb_id}_levels.pse').resolve()}",
        ]
    )
    pml.write_text("\n".join(lines) + "\n", encoding="utf-8")
    subprocess.run(["mamba", "run", "-n", "sbvs", "pymol", "-cq", str(pml)], check=True)


def run_pymol_properties(
    pdb_path: Path,
    stats: dict[str, object],
    pqr_path: Path | None,
    dx_path: Path | None,
    align_pair: tuple[str, str] | None,
    fig_dir: Path,
    sess_dir: Path,
) -> float | None:
    """Render hydrophobicity, APBS electrostatics, polar H-bonds, and optional chain alignment."""
    pdb_id = pdb_path.stem
    pml = sess_dir / f"{pdb_id}_props.pml"
    color_h_script = (PYMOL_HELPERS_DIR / "color_h.py").resolve()
    iface_script = (PYMOL_HELPERS_DIR / "interface_residues.py").resolve()

    lines = [
        f"run {color_h_script}",
        f"run {iface_script}",
        f"load {pdb_path.resolve()}, prot",
        "bg_color white",
        "set ray_opaque_background, 1",
        "remove resn HOH",
        # 1. Hydrophobicity surface (Lecture 3 §5.1)
        "hide everything",
        "show surface, prot",
        "color_h prot",
        "orient prot",
        "ray 1400, 1000",
        f"png {(fig_dir / f'{pdb_id}_04-hydrophobicity.png').resolve()}, dpi=300",
    ]

    # 2. Electrostatic potential surface (Lecture 3 §5.2)
    if pqr_path and dx_path:
        lines.extend(
            [
                f"load {pqr_path.resolve()}, prot_pqr",
                f"load {dx_path.resolve()}, potential",
                "hide everything, prot",
                "show surface, prot_pqr",
                "ramp_new elec_potential, potential, [-10, 0, 10]",
                "set surface_color, elec_potential, prot_pqr",
                "orient prot_pqr",
                "ray 1400, 1000",
                f"png {(fig_dir / f'{pdb_id}_05-electrostatic.png').resolve()}, dpi=300",
                "delete prot_pqr",
                "delete elec_potential",
                "delete potential",
            ]
        )

    # 3. Polar contacts / Hydrogen bonds (Lecture 3 §5.4)
    lines.extend(
        [
            "hide everything, prot",
            "show cartoon, prot",
            "color gray85, prot",
            "distance hbonds, prot and polymer.protein, prot and not polymer.protein, mode=2",
            "distance internal_hb, prot and chain A, prot and not chain A, mode=2",
            "show sticks, byres (prot and not polymer.protein)",
            "color yellow, hbonds",
            "color cyan, internal_hb",
            "orient prot",
            "ray 1400, 1000",
            f"png {(fig_dir / f'{pdb_id}_06-hbonds.png').resolve()}, dpi=300",
        ]
    )

    # 4. 3D Monomer Alignment if >=2 chains (Lecture 3 §5.5)
    chains = [c for c in str(stats["chains"]).split(",") if c and c != "None"]
    pair = align_pair if align_pair else ((chains[0], chains[1]) if len(chains) >= 2 else None)
    if pair:
        c1, c2 = pair
        lines.extend(
            [
                f"create mon1, prot and chain {c1}",
                f"create mon2, prot and chain {c2}",
                "hide everything",
                "show cartoon, mon1 or mon2",
                "color marine, mon1",
                "color orange, mon2",
                "python",
                "aln = cmd.align('mon1', 'mon2')",
                "print(f'ALIGN_RMSD={aln[0]:.4f} ALN_ATOMS={aln[1]}')",
                "python end",
                "orient mon1",
                "ray 1400, 1000",
                f"png {(fig_dir / f'{pdb_id}_07-align_{c1}-{c2}.png').resolve()}, dpi=300",
            ]
        )

    lines.append(f"save {(sess_dir / f'{pdb_id}_props.pse').resolve()}")
    pml.write_text("\n".join(lines) + "\n", encoding="utf-8")

    proc = subprocess.run(
        ["mamba", "run", "-n", "sbvs", "pymol", "-cq", str(pml)],
        capture_output=True,
        text=True,
        check=True,
    )
    for out_line in proc.stdout.splitlines():
        if "ALIGN_RMSD=" in out_line:
            return float(out_line.split("ALIGN_RMSD=")[1].split()[0])
    return None


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    levels_pdbs: list[str] = []
    props_pdbs: list[str] = []
    if args.preset:
        levels_pdbs.extend(PRESETS[args.preset]["levels"])
        props_pdbs.extend(PRESETS[args.preset]["props"])
    if args.pdb:
        user_pdbs = [p.strip().upper() for p in args.pdb.split(",") if p.strip()]
        if args.mode in ("levels", "both"):
            levels_pdbs.extend(user_pdbs)
        if args.mode in ("properties", "both"):
            props_pdbs.extend(user_pdbs)

    align_pair: tuple[str, str] | None = None
    if args.align_chains:
        parts = [c.strip() for c in args.align_chains.split(",") if c.strip()]
        if len(parts) != 2:
            print("Error: --align-chains expects exactly two chains like A,B.", file=sys.stderr)
            return 1
        align_pair = (parts[0], parts[1])

    if args.dry_run:
        print(f"[DRY-RUN] Levels PDBs     : {', '.join(levels_pdbs) or 'None'}")
        print(f"[DRY-RUN] Properties PDBs : {', '.join(props_pdbs) or 'None'}")
        print(f"[DRY-RUN] Output dir      : {args.output_dir}")
        return 0

    out_dir = args.output_dir
    pdb_dir = out_dir / "structures"
    fig_dir = out_dir / "figures"
    sess_dir = out_dir / "sessions"
    apbs_dir = out_dir / "apbs"
    table_dir = out_dir / "tables"
    for d in (pdb_dir, fig_dir, sess_dir, apbs_dir, table_dir):
        d.mkdir(parents=True, exist_ok=True)

    summary_rows: list[str] = [
        "# Lecture 3 Structural Analysis Summary",
        "",
        "| PDB ID | Chains | Total Residues | Cys Residues | Disulfide Bonds "
        "| Monomer Align RMSD (A) | APBS Electrostatics |",
        "| :--- | :--- | ---: | ---: | ---: | ---: | :--- |",
    ]
    all_pdbs = list(dict.fromkeys([*levels_pdbs, *props_pdbs]))
    for pdb_id in all_pdbs:
        print(f"[INFO] Processing {pdb_id}...", file=sys.stderr)
        pdb_path = download_pdb(pdb_id, pdb_dir)
        stats = inspect_structure_stats(pdb_path)
        if pdb_id in levels_pdbs:
            run_pymol_levels(pdb_path, stats, fig_dir, sess_dir)

        rmsd: float | None = None
        apbs_status = "Skipped"
        if pdb_id in props_pdbs:
            pqr_path, dx_path = run_offline_apbs(pdb_path, apbs_dir / pdb_id)
            apbs_status = "Local APBS OK" if dx_path else "Fallback"
            rmsd = run_pymol_properties(
                pdb_path, stats, pqr_path, dx_path, align_pair, fig_dir, sess_dir
            )

        rmsd_str = f"{rmsd:.3f}" if rmsd is not None else "N/A"
        summary_rows.append(
            f"| `{pdb_id}` | `{stats['chains']}` | {stats['residues']} | "
            f"{stats['cys_residues']} | {stats['disulfide_bonds']} | {rmsd_str} | {apbs_status} |"
        )

    summary_path = table_dir / "lecture-03_structure-summary.md"
    summary_path.write_text("\n".join(summary_rows) + "\n", encoding="utf-8")
    print(f"[OK] Saved structural summary to {summary_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
