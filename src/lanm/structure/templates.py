"""Deterministic parsing helpers for template structure inputs."""

from __future__ import annotations

from pathlib import Path

import gemmi

from lanm.data.fetch import require_path
from lanm.models import AtomRecord

_MISSING_CIF_VALUE = {"", ".", "?"}
_ATOM_SITE_TAGS = [
    "_atom_site.group_PDB",
    "_atom_site.id",
    "_atom_site.label_atom_id",
    "_atom_site.label_alt_id",
    "_atom_site.auth_comp_id",
    "_atom_site.auth_asym_id",
    "_atom_site.auth_seq_id",
    "_atom_site.pdbx_PDB_ins_code",
    "_atom_site.Cartn_x",
    "_atom_site.Cartn_y",
    "_atom_site.Cartn_z",
    "_atom_site.occupancy",
    "_atom_site.B_iso_or_equiv",
    "_atom_site.type_symbol",
    "_atom_site.pdbx_formal_charge",
    "_atom_site.pdbx_PDB_model_num",
]


def _parse_optional_float(value: str) -> float | None:
    stripped = value.strip()
    return float(stripped) if stripped else None


def _normalize_cif_value(value: str) -> str:
    stripped = value.strip().strip("'")
    return "" if stripped in _MISSING_CIF_VALUE else stripped


def _parse_structure_id(path: Path) -> str:
    return path.stem.upper()


def _read_cif_block(path: Path) -> gemmi.cif.Block:
    require_path(path)
    return gemmi.cif.read_file(str(path)).sole_block()


def read_cif_experimental_method(path: Path) -> str:
    """Read the experimental method annotation from a mmCIF file."""
    methods = {
        _normalize_cif_value(value)
        for value in _read_cif_block(path).find_values("_exptl.method")
    }
    methods.discard("")
    return "; ".join(sorted(methods))


def parse_cif_atom_records(
    path: Path,
    structure_id: str | None = None,
    model_number: int = 1,
) -> list[AtomRecord]:
    """Parse mmCIF atom rows into AtomRecord values using author identifiers."""
    block = _read_cif_block(path)
    resolved_structure_id = structure_id or _parse_structure_id(path)
    selected_model = str(model_number)
    atoms: list[AtomRecord] = []
    for row in block.find(_ATOM_SITE_TAGS):
        (
            record_type,
            atom_serial,
            atom_name,
            alt_loc,
            residue_name,
            chain_id,
            residue_seq,
            insertion_code,
            x,
            y,
            z,
            occupancy,
            b_factor,
            element,
            charge,
            row_model_number,
        ) = row
        if _normalize_cif_value(row_model_number) not in {"", selected_model}:
            continue
        normalized_residue_seq = _normalize_cif_value(residue_seq)
        if not normalized_residue_seq:
            continue
        atoms.append(
            AtomRecord(
                structure_id=resolved_structure_id,
                record_type=_normalize_cif_value(record_type).upper(),
                atom_serial=int(_normalize_cif_value(atom_serial)),
                atom_name=_normalize_cif_value(atom_name),
                alt_loc=_normalize_cif_value(alt_loc),
                residue_name=_normalize_cif_value(residue_name).upper(),
                chain_id=_normalize_cif_value(chain_id),
                residue_seq=int(normalized_residue_seq),
                insertion_code=_normalize_cif_value(insertion_code),
                x=float(_normalize_cif_value(x)),
                y=float(_normalize_cif_value(y)),
                z=float(_normalize_cif_value(z)),
                occupancy=_parse_optional_float(_normalize_cif_value(occupancy)),
                b_factor=_parse_optional_float(_normalize_cif_value(b_factor)),
                element=_normalize_cif_value(element).upper(),
                charge=_normalize_cif_value(charge),
            )
        )
    return atoms
