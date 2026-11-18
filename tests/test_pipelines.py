"""Offline unit and dry-run tests for BT338IU pipeline scripts."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

import alphafold_client  # noqa: E402
import homology_model  # noqa: E402
import rnaseq_pipeline  # noqa: E402


def test_parse_groups_xlsx() -> None:
    xlsx_path = REPO_ROOT / "references" / "group-datasets.xlsx"
    groups = rnaseq_pipeline.parse_groups_xlsx(xlsx_path)
    assert "DATABI_G1" in groups
    assert "DATABI_DEMO" in groups
    assert len(groups["DATABI_G4"].conditions) == 2
    g1 = groups["DATABI_G1"]
    assert g1.genome_gcf == "GCF_016027375.1"
    assert sum(len(v) for v in g1.conditions.values()) == 6


def test_write_compat_summary(tmp_path: Path) -> None:
    fp_json = tmp_path / "fastp.json"
    fp_json.write_text(
        json.dumps(
            {
                "summary": {
                    "before_filtering": {"total_reads": 1000, "q30_rate": 0.91},
                    "after_filtering": {"total_reads": 950, "q30_rate": 0.96},
                }
            }
        ),
        encoding="utf-8",
    )
    h2_log = tmp_path / "hisat2.log"
    h2_log.write_text("950 reads; of these:\n  95.00% overall alignment rate\n", encoding="utf-8")
    out_log = tmp_path / "compat.log"
    rnaseq_pipeline.write_compat_summary("SRR_TEST", fp_json, h2_log, out_log)
    text = out_log.read_text(encoding="utf-8")
    assert "Surviving: 950 (95.00%)" in text
    assert "95.00% overall alignment rate" in text


def test_pir_format_and_mutation() -> None:
    mutated = homology_model.apply_point_mutation("MSEAAH", "A4V")
    assert mutated == "MSEVAH"
    pir = homology_model.format_pir_entry("TvLDH", mutated)
    assert pir.startswith(">P1;TvLDH\nsequence:TvLDH:::::::0.00: 0.00\nMSEVAH*\n")


def test_parse_modeller_log() -> None:
    sample_log = """\
>> Summary of successfully produced models:
Filename                          molpdf     DOPE score    GA341 score
----------------------------------------------------------------------
TvLDH.B99990001.pdb              1570.23     -42567.12      0.98
TvLDH.B99990002.pdb              1523.84     -42690.45      0.99
"""
    rows = homology_model.parse_model_single_log(sample_log)
    assert len(rows) == 2
    best = min(rows, key=lambda r: float(r["dope"]))
    assert best["filename"] == "TvLDH.B99990002.pdb"


def test_af3_server_job_generation_insulin() -> None:
    seq = "A" * 24 + "B" * 30 + "C" * 35 + "D" * 21
    job = alphafold_client.build_af3_server_job("P01308", "Insulin", 1, {"sequence": seq})
    assert len(job["sequences"]) == 2  # type: ignore[arg-type]
    assert len(job["sequences"][0]["proteinChain"]["sequence"]) == 30  # type: ignore[index]
    assert len(job["sequences"][1]["proteinChain"]["sequence"]) == 21  # type: ignore[index]


def test_low_confidence_spans() -> None:
    scores = [
        (1, "MET", 45.0),
        (2, "ALA", 62.0),
        (3, "GLY", 88.0),
        (4, "SER", 65.0),
        (5, "LEU", 92.0),
    ]
    spans = alphafold_client.find_low_confidence_spans(scores, threshold=70.0)
    assert spans == [(1, 2), (4, 4)]


def test_cli_dry_runs() -> None:
    cases = [
        [str(REPO_ROOT / "src" / "db_miner.py"), "--preset", "w1-ex1", "--dry-run"],
        [str(REPO_ROOT / "src" / "phylo_pipeline.py"), "-a", "AAB72175.1", "--dry-run"],
        [str(REPO_ROOT / "src" / "rnaseq_pipeline.py"), "-g", "DATABI_G1", "--dry-run"],
        [str(REPO_ROOT / "src" / "structure_viz.py"), "--preset", "lecture3-levels", "--dry-run"],
        [str(REPO_ROOT / "src" / "homology_model.py"), "-c", "ex1", "--dry-run"],
        [str(REPO_ROOT / "src" / "alphafold_client.py"), "--preset", "lecture4-af", "--dry-run"],
        [str(REPO_ROOT / "src" / "dge_analysis.R"), "--case", "cs2", "--dry-run"],
    ]
    for cmd in cases:
        res = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True, check=False)
        assert res.returncode == 0, f"Dry-run failed for {cmd}: {res.stderr}"
        assert "[DRY-RUN]" in res.stdout
