"""Deterministic PDB serialization helpers for ProteinMPNN inputs."""

from __future__ import annotations

from collections.abc import Iterable

from lanm.models import AtomRecord


def _preferred_alt_loc_sort_key(atom: AtomRecord) -> tuple[int, float, str, int]:
    occupancy = atom.occupancy if atom.occupancy is not None else -1.0
    return (
        0 if atom.alt_loc == "" else 1,
        -occupancy,
        atom.alt_loc,
        atom.atom_serial,
    )


def select_preferred_atom_conformers(atoms: Iterable[AtomRecord]) -> list[AtomRecord]:
    """Collapse alternate locations to a single deterministic conformer per atom site."""
    grouped: dict[tuple[str, str, int, str, str], list[AtomRecord]] = {}
    for atom in atoms:
        atom_key = (
            atom.chain_id,
            atom.residue_name,
            atom.residue_seq,
            atom.insertion_code,
            atom.atom_name,
        )
        grouped.setdefault(atom_key, []).append(atom)
    return sorted(
        (
            min(group, key=_preferred_alt_loc_sort_key)
            for group in grouped.values()
        ),
        key=lambda atom: (
            atom.chain_id,
            atom.residue_seq,
            atom.insertion_code,
            atom.atom_serial,
            atom.atom_name,
        ),
    )


def _format_atom_name(atom_name: str, element: str) -> str:
    if len(atom_name) == 4:
        return atom_name
    if len(element) == 1:
        return f" {atom_name:<3}"
    return f"{atom_name:<4}"


def _format_ter_line(serial: int, atom: AtomRecord) -> str:
    return (
        f"TER   {serial:>5}      {atom.residue_name:>3} {atom.chain_id:1}"
        f"{atom.residue_seq:>4}{(atom.insertion_code or ' '):1}"
    )


def _format_atom_line(serial: int, atom: AtomRecord) -> str:
    atom_name = _format_atom_name(atom.atom_name, atom.element_upper)
    occupancy = atom.occupancy if atom.occupancy is not None else 1.00
    b_factor = atom.b_factor if atom.b_factor is not None else 0.00
    return (
        f"{atom.record_type:<6}{serial:>5} {atom_name}{' ':1}{atom.residue_name:>3} "
        f"{atom.chain_id:1}{atom.residue_seq:>4}{(atom.insertion_code or ' '):1}   "
        f"{atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}"
        f"{occupancy:>6.2f}{b_factor:>6.2f}          "
        f"{atom.element_upper:>2}{'':>2}"
    )


def _ordered_atom_key(
    atom: AtomRecord,
    chain_order: dict[str, int],
) -> tuple[int, int, str, int, str]:
    return (
        chain_order.get(atom.chain_id, len(chain_order)),
        atom.residue_seq,
        atom.insertion_code,
        atom.atom_serial,
        atom.atom_name,
    )


def render_selected_structure_pdb(
    atoms: Iterable[AtomRecord],
    *,
    selected_chain_ids: tuple[str, ...],
    included_het_residue_keys: frozenset[tuple[str, str, int, str]] = frozenset(),
) -> str:
    """Render a deterministic PDB for selected chains plus chosen HETATM residues."""
    chain_order = {chain_id: index for index, chain_id in enumerate(selected_chain_ids)}
    filtered_atoms = [
        atom
        for atom in atoms
        if (
            atom.record_type == "ATOM"
            and atom.chain_id in chain_order
        )
        or (
            atom.record_type == "HETATM"
            and atom.residue_key in included_het_residue_keys
        )
    ]
    selected_atoms = select_preferred_atom_conformers(filtered_atoms)
    extra_chain_ids = tuple(
        sorted(
            {
                atom.chain_id
                for atom in selected_atoms
                if atom.chain_id not in chain_order
            }
        )
    )
    ordered_chain_ids = (*selected_chain_ids, *extra_chain_ids)
    ordered_chain_map = {chain_id: index for index, chain_id in enumerate(ordered_chain_ids)}

    polymer_atoms_by_chain: dict[str, list[AtomRecord]] = {chain_id: [] for chain_id in ordered_chain_ids}
    hetero_atoms_by_chain: dict[str, list[AtomRecord]] = {chain_id: [] for chain_id in ordered_chain_ids}
    for atom in sorted(selected_atoms, key=lambda item: _ordered_atom_key(item, ordered_chain_map)):
        if atom.record_type == "ATOM":
            polymer_atoms_by_chain[atom.chain_id].append(atom)
        else:
            hetero_atoms_by_chain[atom.chain_id].append(atom)

    lines: list[str] = []
    serial = 1
    for chain_id in ordered_chain_ids:
        chain_polymer_atoms = polymer_atoms_by_chain[chain_id]
        for atom in chain_polymer_atoms:
            lines.append(_format_atom_line(serial, atom))
            serial += 1
        if chain_polymer_atoms:
            lines.append(_format_ter_line(serial, chain_polymer_atoms[-1]))
            serial += 1
        for atom in hetero_atoms_by_chain[chain_id]:
            lines.append(_format_atom_line(serial, atom))
            serial += 1
    lines.append("END")
    return "\n".join(lines) + "\n"


def render_protein_pdb(
    atoms: Iterable[AtomRecord],
    *,
    selected_chain_ids: tuple[str, ...],
) -> str:
    """Render ATOM-only PDB text for the requested chains."""
    return render_selected_structure_pdb(
        atoms,
        selected_chain_ids=selected_chain_ids,
    )
