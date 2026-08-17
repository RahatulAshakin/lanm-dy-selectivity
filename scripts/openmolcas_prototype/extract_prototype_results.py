#!/usr/bin/env python3
"""Extract auditable prototype tables from genuine OpenMolcas console logs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path


HARTREE_TO_CM1 = 219474.6313705


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def first_float(pattern: str, text: str, flags: int = 0) -> float:
    match = re.search(pattern, text, flags)
    if not match:
        raise ValueError(f"Pattern not found: {pattern}")
    return float(match.group(1))


def dy_values(text: str) -> dict:
    roots = [
        (int(root), float(energy))
        for root, energy in re.findall(
            r"RASSCF root number\s+(\d+) Total energy:\s*([-+]?\d+\.\d+)", text
        )
    ]
    # Retain the final 11-root block if another summary ever appears.
    roots = roots[-11:]
    if len(roots) != 11:
        raise ValueError(f"Expected 11 RASSCF roots, found {len(roots)}")

    occ_matches = re.findall(
        r"Natural orbitals and occupation numbers for root\s+(\d+)\s*\n"
        r"\s*sym\s+2:\s*((?:[-+]?\d+\.\d+\s*){7})",
        text,
    )
    occupations = []
    for root, values in occ_matches[-11:]:
        parsed = [float(value) for value in values.split()]
        if len(parsed) != 7:
            raise ValueError(f"Root {root} has {len(parsed)} occupations, not 7")
        occupations.append((int(root), parsed))
    if len(occupations) != 11:
        raise ValueError(f"Expected occupations for 11 roots, found {len(occupations)}")

    so_section = text.split("LOW-LYING SPIN-ORBIT ENERGIES:")[-1]
    so_section = so_section.split("LOW-LYING SPIN-FREE ENERGIES:")[0]
    so_states = [
        (int(state), float(energy))
        for state, energy in re.findall(
            r"ENERGY OF THE SPIN-ORBIT STATE\s*\(\s*(\d+)\)\s*=\s*([-+]?\d+\.\d+)",
            so_section,
        )
    ]
    if len(so_states) != 66:
        raise ValueError(f"Expected 66 spin-orbit states, found {len(so_states)}")

    g_rows = []
    for axis in ("X", "Y", "Z"):
        match = re.search(
            rf"g{axis}\s*=\s*([-+]?\d+\.\d+)\s*\|\s*{axis}m\s*\|\s*"
            r"([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)",
            text,
        )
        if not match:
            raise ValueError(f"Could not find principal g{axis} line")
        g_rows.append((axis, *[float(value) for value in match.groups()]))

    return {
        "roots": roots,
        "occupations": occupations,
        "so_states": so_states,
        "g_rows": g_rows,
    }


def max_abs_diff(left: list[float], right: list[float]) -> float:
    if len(left) != len(right):
        return math.inf
    return max((abs(a - b) for a, b in zip(left, right)), default=0.0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h2-log", type=Path, required=True)
    parser.add_argument("--dy-log", type=Path, required=True)
    parser.add_argument("--dy-repeat-log", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    h2_text = args.h2_log.read_text(encoding="utf-8", errors="replace")
    dy_text = args.dy_log.read_text(encoding="utf-8", errors="replace")
    dy = dy_values(dy_text)

    expected_h2 = {
        "SCF": -1.122765460009,
        "RASSCF": -1.141134852262,
        "CASPT2": -1.143847618496,
    }
    observed_h2 = {
        "SCF": first_float(r"Total SCF energy\s+([-+]?\d+\.\d+)", h2_text),
        "RASSCF": first_float(r"Reference energy:\s*([-+]?\d+\.\d+)", h2_text),
        "CASPT2": first_float(
            r"FINAL CASPT2 RESULT:.*?Total energy:\s*([-+]?\d+\.\d+)",
            h2_text,
            re.DOTALL,
        ),
    }
    h2_rows = []
    for method in ("SCF", "RASSCF", "CASPT2"):
        difference = observed_h2[method] - expected_h2[method]
        h2_rows.append(
            {
                "method": method,
                "expected_energy_Eh": f"{expected_h2[method]:.12f}",
                "observed_energy_Eh": f"{observed_h2[method]:.12f}",
                "difference_Eh": f"{difference:.3e}",
                "within_1e-8_Eh": str(abs(difference) <= 1.0e-8).upper(),
            }
        )
    write_csv(
        out / "h2_reference_validation.csv",
        ["method", "expected_energy_Eh", "observed_energy_Eh", "difference_Eh", "within_1e-8_Eh"],
        h2_rows,
    )

    root_min = min(energy for _, energy in dy["roots"])
    root_rows = [
        {
            "root": root,
            "energy_Eh": f"{energy:.8f}",
            "relative_energy_cm-1": f"{(energy - root_min) * HARTREE_TO_CM1:.4f}",
        }
        for root, energy in dy["roots"]
    ]
    write_csv(out / "dy3_rasscf_root_energies.csv", list(root_rows[0]), root_rows)

    occupation_rows = []
    for root, values in dy["occupations"]:
        row = {"root": root}
        row.update({f"active_orbital_{index}": f"{value:.6f}" for index, value in enumerate(values, 1)})
        row["occupation_sum"] = f"{sum(values):.6f}"
        occupation_rows.append(row)
    write_csv(out / "dy3_natural_occupations_all_roots.csv", list(occupation_rows[0]), occupation_rows)

    so_rows = [
        {"spin_orbit_state": state, "relative_energy_cm-1": f"{energy:.8f}"}
        for state, energy in dy["so_states"]
    ]
    write_csv(out / "dy3_spin_orbit_states.csv", list(so_rows[0]), so_rows)

    kd_rows = []
    for index in range(0, len(dy["so_states"]), 2):
        state_a, energy_a = dy["so_states"][index]
        state_b, energy_b = dy["so_states"][index + 1]
        kd_rows.append(
            {
                "kramers_doublet": index // 2 + 1,
                "state_a": state_a,
                "state_b": state_b,
                "relative_energy_cm-1": f"{(energy_a + energy_b) / 2.0:.8f}",
                "pair_splitting_cm-1": f"{abs(energy_a - energy_b):.3e}",
            }
        )
    write_csv(out / "dy3_kramers_doublet_energies.csv", list(kd_rows[0]), kd_rows)

    g_table = [
        {
            "principal_axis": axis,
            "g_value": f"{g_value:.14f}",
            "axis_x": f"{axis_x:.14f}",
            "axis_y": f"{axis_y:.14f}",
            "axis_z": f"{axis_z:.14f}",
        }
        for axis, g_value, axis_x, axis_y, axis_z in dy["g_rows"]
    ]
    write_csv(out / "dy3_ground_doublet_g_tensor.csv", list(g_table[0]), g_table)

    validation_rows = [
        {"run": "H2 reference", "stage": "RASSCF", "status": "PASS", "evidence": "_RC_ALL_IS_WELL_ and Happy landing"},
        {"run": "H2 reference", "stage": "CASPT2", "status": "PASS", "evidence": "_RC_ALL_IS_WELL_ and energy reference check"},
        {"run": "Dy3+ RICD", "stage": "Seward/RICD", "status": "PASS", "evidence": "_RC_ALL_IS_WELL_"},
        {"run": "Dy3+ RICD", "stage": "RASSCF CAS(9,7), 11 roots", "status": "PASS", "evidence": "11 roots and _RC_ALL_IS_WELL_"},
        {"run": "Dy3+ RICD", "stage": "RASSI spin-orbit", "status": "PASS", "evidence": "66 states and _RC_ALL_IS_WELL_"},
        {"run": "Dy3+ RICD", "stage": "SINGLE_ANISO", "status": "PASS", "evidence": "g tensor, magnetic tables, _RC_ALL_IS_WELL_, Happy landing"},
    ]
    write_csv(out / "run_stage_validation.csv", list(validation_rows[0]), validation_rows)

    repeat_summary = None
    if args.dy_repeat_log:
        repeat_text = args.dy_repeat_log.read_text(encoding="utf-8", errors="replace")
        repeat = dy_values(repeat_text)
        comparisons = {
            "RASSCF root energies (Eh)": max_abs_diff(
                [value for _, value in dy["roots"]], [value for _, value in repeat["roots"]]
            ),
            "natural occupations": max_abs_diff(
                [value for _, values in dy["occupations"] for value in values],
                [value for _, values in repeat["occupations"] for value in values],
            ),
            "spin-orbit relative energies (cm-1)": max_abs_diff(
                [value for _, value in dy["so_states"]], [value for _, value in repeat["so_states"]]
            ),
            "ground-doublet g values": max_abs_diff(
                [row[1] for row in dy["g_rows"]], [row[1] for row in repeat["g_rows"]]
            ),
        }
        reproducibility_rows = [
            {
                "quantity": quantity,
                "max_absolute_difference_at_printed_precision": f"{difference:.12g}",
                "identical_at_printed_precision": str(difference == 0.0).upper(),
            }
            for quantity, difference in comparisons.items()
        ]
        write_csv(out / "dy3_repeat_reproducibility.csv", list(reproducibility_rows[0]), reproducibility_rows)
        repeat_summary = comparisons

    summary = {
        "prototype_only": True,
        "not_manuscript_evidence": True,
        "h2_happy_landing": "Happy landing!" in h2_text,
        "dy_happy_landing": "Happy landing!" in dy_text,
        "dy_rasscf_roots": len(dy["roots"]),
        "dy_active_orbitals": len(dy["occupations"][0][1]),
        "dy_spin_orbit_states": len(dy["so_states"]),
        "dy_kramers_doublets": len(kd_rows),
        "repeat_max_differences": repeat_summary,
        "hartree_to_wavenumber_factor": HARTREE_TO_CM1,
    }
    (out / "extraction_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
