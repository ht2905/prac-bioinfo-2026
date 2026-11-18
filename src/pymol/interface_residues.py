"""Solvent-accessible surface area (dASA) interface residue finder for PyMOL."""

from __future__ import annotations

from pymol import cmd, stored


def interface_residues(
    cmpx: str,
    cA: str = "c. A",
    cB: str = "c. B",
    cutoff: float = 1.0,
    selName: str = "interface",
) -> list[tuple[str, str]]:
    """Identify and select interface residues between two chains/selections in cmpx."""
    cutoff_val = float(cutoff)
    temp_c = "__tmp_iface_c"
    temp_a = "__tmp_iface_a"
    temp_b = "__tmp_iface_b"
    temp_ab = "__tmp_iface_ab"

    prev_dot_solvent = cmd.get("dot_solvent")
    cmd.set("dot_solvent", 1)

    cmd.create(temp_c, f"({cmpx}) and (({cA}) or ({cB}))")
    cmd.create(temp_a, f"({cmpx}) and ({cA})")
    cmd.create(temp_b, f"({cmpx}) and ({cB})")

    cmd.get_area(temp_c, load_b=1)
    cmd.get_area(temp_a, load_b=1)
    cmd.get_area(temp_b, load_b=1)

    cmd.create(temp_ab, f"{temp_a} or {temp_b}")

    stored.complex_b: dict[tuple[str, str, str], float] = {}
    stored.split_b: dict[tuple[str, str, str], float] = {}
    cmd.iterate(temp_c, "stored.complex_b[(chain, resi, name)] = b")
    cmd.iterate(temp_ab, "stored.split_b[(chain, resi, name)] = b")

    interface_res: set[tuple[str, str]] = set()
    for key, area_split in stored.split_b.items():
        area_complex = stored.complex_b.get(key, area_split)
        if (area_split - area_complex) >= cutoff_val:
            interface_res.add((key[0], key[1]))

    cmd.delete(temp_c)
    cmd.delete(temp_a)
    cmd.delete(temp_b)
    cmd.delete(temp_ab)
    cmd.set("dot_solvent", prev_dot_solvent)

    if interface_res:
        sel_terms = [
            f"(chain {ch} and resi {resi})" if ch else f"(resi {resi})"
            for ch, resi in sorted(interface_res)
        ]
        cmd.select(selName, f"({cmpx}) and (" + " or ".join(sel_terms) + ")")
    else:
        cmd.select(selName, "none")

    return sorted(interface_res)


cmd.extend("interfaceResidues", interface_residues)
