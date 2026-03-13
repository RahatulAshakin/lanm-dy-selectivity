"""Residue-level helpers for polymer sequence extraction."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable

from lanm.models import AtomRecord

AMINO_ACID_CODES = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


@dataclass(frozen=True, slots=True)
class ResidueRecord:
    chain_id: str
    residue_seq: int
    insertion_code: str
    residue_name: str

    @property
    def one_letter_code(self) -> str:
        return AMINO_ACID_CODES.get(self.residue_name, "X")

    @property
    def residue_id(self) -> str:
        insertion = self.insertion_code or ""
        return f"{self.chain_id}:{self.residue_seq}{insertion}"


def collect_polymer_residues(atoms: Iterable[AtomRecord]) -> dict[str, tuple[ResidueRecord, ...]]:
    """Collect ordered polymer residues from ATOM records, grouped by chain."""
    residues: dict[str, list[ResidueRecord]] = defaultdict(list)
    seen: dict[str, set[tuple[int, str]]] = defaultdict(set)
    for atom in sorted(atoms, key=lambda item: (item.chain_id, item.residue_seq, item.insertion_code, item.atom_serial)):
        if atom.record_type != "ATOM":
            continue
        residue_key = (atom.residue_seq, atom.insertion_code)
        if residue_key in seen[atom.chain_id]:
            continue
        seen[atom.chain_id].add(residue_key)
        residues[atom.chain_id].append(
            ResidueRecord(
                chain_id=atom.chain_id,
                residue_seq=atom.residue_seq,
                insertion_code=atom.insertion_code,
                residue_name=atom.residue_name,
            )
        )
    return {chain_id: tuple(chain_residues) for chain_id, chain_residues in residues.items()}
