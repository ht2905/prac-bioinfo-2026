#!/usr/bin/env python3
"""Automate Lecture 4 MODELLER homology modeling (PIR, align2d, automodel, DOPE/GA341, PyMOL)."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import requests

UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"
RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
TVLDH_TUTORIAL_ALI_URL = "https://salilab.org/modeller/tutorial/basic-example/TvLDH.ali"


@dataclass(frozen=True)
class ModelingCase:
    case_id: str
    title: str
    target_code: str
    target_uniprot: str | None
    template_pdb: str
    template_chain: str
    expected_id: str
    default_mutation: str | None = None
    notes: str = ""


PREDEFINED_CASES: dict[str, ModelingCase] = {
    "tvldh": ModelingCase(
        case_id="tvldh",
        title="Trichomonas vaginalis Lactate Dehydrogenase (Tutorial)",
        target_code="TvLDH",
        target_uniprot="Q27892",
        template_pdb="1BDM",
        template_chain="A",
        expected_id="~45%",
        notes="Manual basic tutorial using TvMDH (1BDM:A) as template.",
    ),
    "ex1": ModelingCase(
        case_id="ex1",
        title="Human Aldose Reductase",
        target_code="HsAKR1B1",
        target_uniprot="P15121",
        template_pdb="1AH3",
        template_chain="A",
        expected_id="~68%",
        notes="High-identity conserved enzyme template (pig aldose reductase).",
    ),
    "ex2": ModelingCase(
        case_id="ex2",
        title="Human Carbonic Anhydrase IX (CA IX)",
        target_code="HsCA9",
        target_uniprot="Q16790",
        template_pdb="3HS4",
        template_chain="A",
        expected_id="~44%",
        notes="Moderate identity; inspect active-site Zn-histidine conservation.",
    ),
    "ex3": ModelingCase(
        case_id="ex3",
        title="Human Cytochrome P450 2E1 (CYP2E1)",
        target_code="HsCYP2E1",
        target_uniprot="P05181",
        template_pdb="3E6I",
        template_chain="A",
        expected_id="~37%",
        notes="Low identity (<40%); verify conserved heme-binding cysteine pocket.",
    ),
    "ex4": ModelingCase(
        case_id="ex4",
        title="Human Dihydrofolate Reductase (DHFR) Point Mutant",
        target_code="HsDHFR_mut",
        target_uniprot="P00374",
        template_pdb="1RX2",
        template_chain="A",
        expected_id="~99%",
        default_mutation="F32S",
        notes="Single-residue point mutation modeled against wild-type DHFR.",
    ),
    "ex5": ModelingCase(
        case_id="ex5",
        title="E. coli Glutamine Synthetase",
        target_code="EcGlnA",
        target_uniprot="P0A9S5",
        template_pdb="2GLS",
        template_chain="A",
        expected_id="~82%",
        notes="Subunit model of a 12-mer oligomeric enzyme (Salmonella template).",
    ),
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run MODELLER homology modeling pipelines and PyMOL superpositions for Lecture 4."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  # Run Exercise 1 (Human Aldose Reductase P15121 on template 1AH3:A)
  homology_model.py -c ex1 --num-models 5 -o lecture-04/modeller

  # Run custom target UniProt + template PDB with a point mutation
  homology_model.py --target-uniprot P00374 --template-pdb 1RX2 --mutation F32S -o .scratch/dhfr
""",
    )
    parser.add_argument(
        "-c",
        "--case",
        choices=(*PREDEFINED_CASES.keys(), "all"),
        default=None,
        help="Predefined Lecture 4 MODELLER case (tvldh, ex1..ex5, or all).",
    )
    parser.add_argument(
        "--target-uniprot",
        type=str,
        default=None,
        help="Custom target UniProt accession (overrides --case target).",
    )
    parser.add_argument(
        "--template-pdb",
        type=str,
        default=None,
        help="Custom template 4-letter PDB code (overrides --case template).",
    )
    parser.add_argument(
        "--template-chain",
        type=str,
        default="A",
        help="Template chain ID (default: A).",
    )
    parser.add_argument(
        "--mutation",
        type=str,
        default=None,
        help="Optional point mutation in format WT_POS_MUT, e.g. F32S.",
    )
    parser.add_argument(
        "--num-models",
        type=int,
        default=5,
        help="Number of candidate MODELLER structures to generate (default: 5).",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("lecture-04/modeller"),
        help=(
            "Output dir for PIR/ALI files, logs, PDB models, and renders "
            "(default: lecture-04/modeller)."
        ),
    )
    parser.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Preview modeling cases and commands without invoking MODELLER or network calls.",
    )
    args = parser.parse_args(argv)
    if not args.case and not (args.target_uniprot and args.template_pdb):
        parser.error(
            "Provide -c/--case or both --target-uniprot and --template-pdb. "
            "Run 'homology_model.py --help' for usage."
        )
    return args


def apply_point_mutation(seq: str, mutation: str) -> str:
    """Apply a 1-indexed <WT><POS><MUT> amino-acid substitution to seq."""
    m = re.fullmatch(r"([A-Za-z])(\d+)([A-Za-z])", mutation.strip())
    if not m:
        raise ValueError(f"Invalid mutation spec '{mutation}'; expected format like F32S.")
    wt, pos_str, mut = m.group(1).upper(), m.group(2), m.group(3).upper()
    idx = int(pos_str) - 1
    if idx < 0 or idx >= len(seq):
        raise ValueError(
            f"Mutation position {pos_str} out of bounds for sequence length {len(seq)}."
        )
    if seq[idx] != wt:
        print(
            f"[WARN] Residue {pos_str} is '{seq[idx]}' (expected '{wt}'); substituting '{mut}'.",
            file=sys.stderr,
        )
    return seq[:idx] + mut + seq[idx + 1 :]


def format_pir_entry(code: str, sequence: str) -> str:
    """Format a protein sequence in MODELLER .pir alignment input syntax."""
    clean = re.sub(r"[^A-Za-z]", "", sequence).upper()
    wrapped = "\n".join(clean[i : i + 75] for i in range(0, len(clean), 75))
    return f">P1;{code}\nsequence:{code}:::::::0.00: 0.00\n{wrapped}*\n"


def fetch_target_sequence(case: ModelingCase, mutation: str | None) -> str:
    if case.case_id == "tvldh":
        resp = requests.get(TVLDH_TUTORIAL_ALI_URL, timeout=30)
        if resp.ok and ">P1;" in resp.text:
            body_lines = [
                ln.strip()
                for ln in resp.text.splitlines()
                if ln.strip() and not ln.startswith((">P1;", "sequence:"))
            ]
            seq = "".join(body_lines).rstrip("*")
            if seq:
                return apply_point_mutation(seq, mutation) if mutation else seq

    if not case.target_uniprot:
        raise ValueError(f"No UniProt accession configured for case {case.case_id}.")
    resp = requests.get(UNIPROT_FASTA_URL.format(acc=case.target_uniprot), timeout=30)
    resp.raise_for_status()
    seq = "".join(line.strip() for line in resp.text.splitlines() if not line.startswith(">"))
    return apply_point_mutation(seq, mutation) if mutation else seq


def download_template_pdb(pdb_id: str, dest: Path) -> Path:
    if not dest.exists():
        resp = requests.get(RCSB_PDB_URL.format(pdb_id=pdb_id.upper()), timeout=45)
        resp.raise_for_status()
        dest.write_text(resp.text, encoding="utf-8")
    return dest


def check_modeller_license() -> str | None:
    """Return configured KEY_MODELLER from env or sbvs config.py, else None."""
    env_key = os.environ.get("KEY_MODELLER", "").strip()
    if env_key and env_key != "XXXX":
        return env_key
    sbvs_root = Path("/home/ht2905/miniforge3/envs/sbvs/lib")
    for cfg in sbvs_root.glob("modeller-*/modlib/modeller/config.py"):
        text = cfg.read_text(encoding="utf-8")
        m = re.search(r"license\s*=\s*r?['\"]([^'\"]+)['\"]", text)
        if m and m.group(1) != "XXXX":
            return m.group(1)
    return None


def parse_model_single_log(log_text: str) -> list[dict[str, str | float]]:
    """Extract filename, molpdf, DOPE, and GA341 rows from MODELLER summary log."""
    rows: list[dict[str, str | float]] = []
    in_summary = False
    for line in log_text.splitlines():
        if "Summary of successfully produced models:" in line:
            in_summary = True
            continue
        if in_summary and line.strip().startswith("-"):
            continue
        if in_summary and ".pdb" in line:
            parts = line.split()
            if len(parts) >= 4:
                rows.append(
                    {
                        "filename": parts[0],
                        "molpdf": float(parts[1]),
                        "dope": float(parts[2]),
                        "ga341": float(parts[3]),
                    }
                )
    return rows


def render_superposition(
    model_pdb: Path, template_pdb: Path, chain: str, out_png: Path
) -> float | None:
    """Superpose best MODELLER structure onto template in headless PyMOL."""
    pml = out_png.with_suffix(".pml")
    lines = [
        f"load {template_pdb.resolve()}, tmpl",
        f"load {model_pdb.resolve()}, mod",
        "bg_color white",
        "set ray_opaque_background, 1",
        "hide everything",
        f"show cartoon, tmpl and chain {chain}",
        "show cartoon, mod",
        f"color gray70, tmpl and chain {chain}",
        "color marine, mod",
        "python",
        f"aln = cmd.align('mod', 'tmpl and chain {chain}')",
        "print(f'SUP_RMSD={aln[0]:.4f}')",
        "python end",
        "orient mod",
        "ray 1400, 1000",
        f"png {out_png.resolve()}, dpi=300",
    ]
    pml.write_text("\n".join(lines) + "\n", encoding="utf-8")
    proc = subprocess.run(
        ["mamba", "run", "-n", "sbvs", "pymol", "-cq", str(pml)],
        capture_output=True,
        text=True,
        check=True,
    )
    for line in proc.stdout.splitlines():
        if "SUP_RMSD=" in line:
            return float(line.split("SUP_RMSD=")[1].strip())
    return None


def run_case(
    case: ModelingCase,
    num_models: int,
    mutation: str | None,
    out_dir: Path,
    license_key: str,
) -> dict[str, str]:
    case_dir = out_dir / case.case_id
    case_dir.mkdir(parents=True, exist_ok=True)
    effective_mut = mutation or case.default_mutation
    seq = fetch_target_sequence(case, effective_mut)

    pir_path = case_dir / f"{case.target_code}.ali"
    pir_path.write_text(format_pir_entry(case.target_code, seq), encoding="utf-8")

    tmpl_code = case.template_pdb.lower()
    tmpl_pdb = download_template_pdb(case.template_pdb, case_dir / f"{tmpl_code}.pdb")
    ali_out = f"{case.target_code}-{tmpl_code}{case.template_chain}.ali"
    pap_out = f"{case.target_code}-{tmpl_code}{case.template_chain}.pap"

    align_script = case_dir / "align2d.py"
    ch = case.template_chain
    align_script.write_text(
        f"""import os
import modeller.config
modeller.config.license = {license_key!r}
from modeller import Environ, Alignment, Model

env = Environ()
aln = Alignment(env)
mdl = Model(env, file='{tmpl_code}', model_segment=('FIRST:{ch}', 'LAST:{ch}'))
aln.append_model(mdl, align_codes='{tmpl_code}{ch}', atom_files='{tmpl_code}.pdb')
aln.append(file='{pir_path.name}', align_codes='{case.target_code}')
aln.align2d(max_gap_length=50)
aln.write(file='{ali_out}', alignment_format='PIR')
aln.write(file='{pap_out}', alignment_format='PAP')
""",
        encoding="utf-8",
    )

    model_script = case_dir / "model-single.py"
    model_script.write_text(
        f"""import modeller.config
modeller.config.license = {license_key!r}
from modeller import Environ
from modeller.automodel import AutoModel, assess

env = Environ()
a = AutoModel(
    env,
    alnfile='{ali_out}',
    knowns='{tmpl_code}{ch}',
    sequence='{case.target_code}',
    assess_methods=(assess.DOPE, assess.GA341),
)
a.starting_model = 1
a.ending_model = {num_models}
a.make()
""",
        encoding="utf-8",
    )

    print(
        f"[INFO] Running MODELLER align2d + automodel ({num_models} models) for {case.case_id}...",
        file=sys.stderr,
    )
    align_log = case_dir / "align2d.log"
    with align_log.open("w", encoding="utf-8") as fh:
        subprocess.run(
            ["mamba", "run", "-n", "sbvs", "python", align_script.name],
            cwd=case_dir,
            stdout=fh,
            stderr=subprocess.STDOUT,
            check=True,
        )

    model_log = case_dir / "model-single.log"
    with model_log.open("w", encoding="utf-8") as fh:
        subprocess.run(
            ["mamba", "run", "-n", "sbvs", "python", model_script.name],
            cwd=case_dir,
            stdout=fh,
            stderr=subprocess.STDOUT,
            check=True,
        )

    scores = parse_model_single_log(model_log.read_text(encoding="utf-8"))
    if not scores:
        raise RuntimeError(f"MODELLER produced no scored models in {model_log}")

    best = min(scores, key=lambda r: float(r["dope"]))
    best_pdb = case_dir / str(best["filename"])
    best_copy = case_dir / f"{case.case_id}_best.pdb"
    shutil.copy2(best_pdb, best_copy)

    fig_png = case_dir / f"{case.case_id}_superposition.png"
    rmsd = render_superposition(best_copy, tmpl_pdb, case.template_chain, fig_png)
    return {
        "case": case.case_id,
        "title": case.title,
        "template": f"{case.template_pdb}:{case.template_chain} ({case.expected_id})",
        "best_model": str(best["filename"]),
        "dope": f"{float(best['dope']):.2f}",
        "ga341": f"{float(best['ga341']):.3f}",
        "rmsd": f"{rmsd:.3f}" if rmsd is not None else "N/A",
    }


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    cases: list[ModelingCase] = []
    if args.case == "all":
        cases.extend(PREDEFINED_CASES.values())
    elif args.case:
        cases.append(PREDEFINED_CASES[args.case])
    else:
        cases.append(
            ModelingCase(
                case_id=args.target_uniprot.lower(),
                title=f"Custom {args.target_uniprot}",
                target_code=args.target_uniprot,
                target_uniprot=args.target_uniprot,
                template_pdb=args.template_pdb.upper(),
                template_chain=args.template_chain,
                expected_id="custom",
            )
        )

    if args.dry_run:
        print(f"[DRY-RUN] MODELLER cases ({len(cases)}), num_models={args.num_models}:")
        for c in cases:
            print(
                f"[DRY-RUN]   - {c.case_id}: {c.target_uniprot} -> "
                f"{c.template_pdb}:{c.template_chain}"
            )
        print(f"[DRY-RUN] Output dir: {args.output_dir}")
        return 0

    license_key = check_modeller_license()
    if not license_key:
        print(
            "Error: MODELLER academic license key is not configured.\n"
            "Register at https://salilab.org/modeller/registration.html and either:\n"
            "  export KEY_MODELLER='YOUR_LICENSE_KEY'\n"
            "or edit /home/ht2905/miniforge3/envs/sbvs/lib/modeller-10.8/modlib/modeller/config.py",
            file=sys.stderr,
        )
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    results = [
        run_case(c, args.num_models, args.mutation, args.output_dir, license_key) for c in cases
    ]

    md_lines = [
        "# Lecture 4A — MODELLER Homology Modeling Summary",
        "",
        "| Case | Protein | Template (%ID) | Best Model | DOPE | GA341 | Template RMSD (A) |",
        "| :--- | :--- | :--- | :--- | ---: | ---: | ---: |",
    ]
    for r in results:
        md_lines.append(
            f"| `{r['case']}` | {r['title']} | `{r['template']}` | "
            f"`{r['best_model']}` | {r['dope']} | {r['ga341']} | {r['rmsd']} |"
        )
    summary_path = args.output_dir / "modeller_summary.md"
    summary_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    print(f"[OK] Saved MODELLER summary to {summary_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
