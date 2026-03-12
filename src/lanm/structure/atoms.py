"""Atom table parsing for local structure CSV exports."""

from __future__ import annotations

import csv
from pathlib import Path

from lanm.models import AtomRecord


def _parse_optional_float(value: str) -> float | None:
    stripped = value.strip()
    return float(stripped) if stripped else None


def _parse_structure_id(path: Path) -> str:
    return path.stem.split("_", 1)[0].upper()


def parse_atom_row(row: dict[str, str], structure_id: str) -> AtomRecord:
    """Parse a CSV row into a typed atom record."""
    return AtomRecord(
        structure_id=structure_id,
        record_type=str(row["record_type"]).strip().upper(),
        atom_serial=int(str(row["atom_serial"]).strip()),
        atom_name=str(row["atom_name"]).strip(),
        alt_loc=str(row.get("alt_loc", "")).strip(),
        residue_name=str(row["residue_name"]).strip().upper(),
        chain_id=str(row["chain_id"]).strip(),
        residue_seq=int(str(row["residue_seq"]).strip()),
        insertion_code=str(row.get("insertion_code", "")).strip(),
        x=float(str(row["x"]).strip()),
        y=float(str(row["y"]).strip()),
        z=float(str(row["z"]).strip()),
        occupancy=_parse_optional_float(str(row.get("occupancy", ""))),
        b_factor=_parse_optional_float(str(row.get("b_factor", ""))),
        element=str(row["element"]).strip().upper(),
        charge=str(row.get("charge", "")).strip(),
    )


def read_atom_records(path: Path, structure_id: str | None = None) -> list[AtomRecord]:
    """Read atom records from a structure CSV export."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        resolved_structure_id = structure_id or _parse_structure_id(path)
        return [parse_atom_row(row, resolved_structure_id) for row in reader]

