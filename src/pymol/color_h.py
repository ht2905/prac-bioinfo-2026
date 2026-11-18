"""Eisenberg hydrophobicity scale coloring commands for PyMOL."""

from __future__ import annotations

from pymol import cmd

EISENBERG_SCALE: dict[str, float] = {
    "ILE": 1.38,
    "PHE": 1.19,
    "VAL": 1.08,
    "LEU": 1.06,
    "TRP": 0.81,
    "MET": 0.64,
    "ALA": 0.62,
    "GLY": 0.48,
    "CYS": 0.29,
    "TYR": 0.26,
    "PRO": 0.12,
    "THR": -0.05,
    "SER": -0.18,
    "HIS": -0.40,
    "GLU": -0.74,
    "ASN": -0.78,
    "GLN": -0.85,
    "ASP": -0.90,
    "LYS": -1.50,
    "ARG": -2.53,
}


def _apply_hydrophobicity_palette(selection: str, green_mode: bool = False) -> None:
    """Color residues along the normalized Eisenberg hydrophobicity gradient."""
    cmd.set_color("hydro_s0", [1.0, 1.0, 1.0])
    for resn, value in EISENBERG_SCALE.items():
        # Normalized to [0, 1] relative to max hydrophobic residue ILE (1.38)
        scaled = max(0.0, min(1.0, value / 1.38)) if value > 0 else 0.0
        color_name = f"hydro_{resn.lower()}_{'g' if green_mode else 'r'}"
        if green_mode:
            rgb = [1.0 - scaled, 1.0, 1.0 - scaled]
        else:
            rgb = [1.0, 1.0 - scaled, 1.0 - scaled]
        cmd.set_color(color_name, rgb)
        cmd.color(color_name, f"({selection}) and resn {resn}")


def color_h(selection: str = "all") -> None:
    """Color hydrophobic residues red to polar residues white (Lecture 3 color_h)."""
    _apply_hydrophobicity_palette(selection, green_mode=False)


def color_h1(selection: str = "all") -> None:
    """Color hydrophobic residues green to polar residues white (Lecture 3 color_h1)."""
    _apply_hydrophobicity_palette(selection, green_mode=True)


cmd.extend("color_h", color_h)
cmd.extend("color_h1", color_h1)
