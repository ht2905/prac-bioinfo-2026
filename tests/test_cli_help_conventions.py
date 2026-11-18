"""Contract tests for standardized --help and CLI error-stream conventions."""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent

ALL_SCRIPTS = [
    REPO_ROOT / "src" / "db_miner.py",
    REPO_ROOT / "src" / "phylo_pipeline.py",
    REPO_ROOT / "src" / "rnaseq_pipeline.py",
    REPO_ROOT / "src" / "structure_viz.py",
    REPO_ROOT / "src" / "homology_model.py",
    REPO_ROOT / "src" / "alphafold_client.py",
    REPO_ROOT / "src" / "dge_analysis.R",
]


@pytest.mark.parametrize(
    "script_path",
    ALL_SCRIPTS,
    ids=[str(p.relative_to(REPO_ROOT)) for p in ALL_SCRIPTS],
)
def test_help_follows_standardized_convention(script_path: Path) -> None:
    res = subprocess.run(
        [str(script_path), "--help"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert res.returncode == 0, f"{script_path} --help failed: {res.stderr}"
    assert res.stderr == "", f"{script_path} --help wrote to stderr: {res.stderr}"

    out = res.stdout
    assert out.startswith(f"usage: {script_path.name} "), out[:120]
    assert "\noptions:\n" in out, f"Missing 'options:' header in {script_path}"
    assert "\nexamples:\n" in out, f"Missing 'examples:' header in {script_path}"

    options_section = out.split("\noptions:\n", 1)[1].split("\nexamples:\n", 1)[0]
    assert not re.search(r"<[a-zA-Z0-9_-]+>", options_section), (
        f"Found <angle-bracket> metavar in {script_path} options:\n{options_section}"
    )


@pytest.mark.parametrize(
    "script_path",
    ALL_SCRIPTS,
    ids=[str(p.relative_to(REPO_ROOT)) for p in ALL_SCRIPTS],
)
def test_unknown_flag_errors_go_to_stderr(script_path: Path) -> None:
    res = subprocess.run(
        [str(script_path), "--definitely-unknown-flag-xyz"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert res.returncode != 0
    assert res.stdout == "", f"Unexpected stdout on invalid flag for {script_path}: {res.stdout}"
    assert res.stderr.strip() != "", f"Expected stderr error message for {script_path}"
    assert "--help" in res.stderr or "usage:" in res.stderr


def test_homology_model_help_does_not_leak_license_key() -> None:
    env = os.environ.copy()
    env["KEY_MODELLER"] = "SENTINEL_MODELLER_KEY_999"
    res = subprocess.run(
        [str(REPO_ROOT / "src" / "homology_model.py"), "--help"],
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert res.returncode == 0
    assert "SENTINEL_MODELLER_KEY_999" not in res.stdout
