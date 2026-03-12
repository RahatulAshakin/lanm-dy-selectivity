"""Geometry utilities for metal-site annotation."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable

from lanm.models import AtomRecord, MetalSiteSummary, ShellAnnotation

WATER_RESIDUES = {"HOH", "WAT", "H2O"}


@dataclass(frozen=True, slots=True)
class AtomDistance:
    atom: AtomRecord
    distance_A: float


def euclidean_distance(atom_a: AtomRecord, atom_b: AtomRecord) -> float:
    """Compute the 3D atom-atom distance in Angstrom."""
    return math.dist((atom_a.x, atom_a.y, atom_a.z), (atom_b.x, atom_b.y, atom_b.z))


def collect_atoms_within_cutoff(
    center: AtomRecord,
    atoms: Iterable[AtomRecord],
    cutoff_A: float,
) -> list[AtomDistance]:
    """Collect atoms within a cutoff, sorted deterministically."""
    matches: list[AtomDistance] = []
    for atom in atoms:
        if atom.atom_serial == center.atom_serial and atom.structure_id == center.structure_id:
            continue
        distance_A = euclidean_distance(center, atom)
        if distance_A <= cutoff_A:
            matches.append(AtomDistance(atom=atom, distance_A=distance_A))
    return sorted(
        matches,
        key=lambda item: (
            round(item.distance_A, 6),
            item.atom.chain_id,
            item.atom.residue_seq,
            item.atom.atom_serial,
        ),
    )


def annotate_metal_site(
    metal_atom: AtomRecord,
    atoms: Iterable[AtomRecord],
    first_shell_cutoff_A: float,
    second_sphere_cutoff_A: float,
) -> tuple[MetalSiteSummary, list[ShellAnnotation]]:
    """Compute first-shell donor atoms and second-sphere residue annotations."""
    nearby_atoms = collect_atoms_within_cutoff(metal_atom, atoms, second_sphere_cutoff_A)
    donor_matches = [
        item for item in nearby_atoms
        if item.atom.element_upper in {"O", "N", "S", "SE"}
        and item.distance_A <= first_shell_cutoff_A
    ]
    first_shell_residue_keys = {item.atom.residue_key for item in donor_matches}
    residue_contacts: dict[tuple[str, str, int, str], AtomDistance] = {}
    for item in nearby_atoms:
        residue_key = item.atom.residue_key
        if residue_key == metal_atom.residue_key and item.atom.chain_id == metal_atom.chain_id:
            continue
        best = residue_contacts.get(residue_key)
        if best is None or item.distance_A < best.distance_A:
            residue_contacts[residue_key] = item
    annotations: list[ShellAnnotation] = []
    for item in donor_matches:
        annotations.append(
            ShellAnnotation(
                structure_id=metal_atom.structure_id,
                metal_site_id=metal_atom.site_id,
                metal_element=metal_atom.element_upper,
                shell_type="first_shell",
                chain_id=item.atom.chain_id,
                residue_name=item.atom.residue_name,
                residue_seq=item.atom.residue_seq,
                atom_name=item.atom.atom_name,
                atom_serial=item.atom.atom_serial,
                element=item.atom.element_upper,
                record_type=item.atom.record_type,
                distance_A=round(item.distance_A, 3),
                is_water=item.atom.residue_name in WATER_RESIDUES,
            )
        )
    for residue_key, item in sorted(
        residue_contacts.items(),
        key=lambda pair: (
            round(pair[1].distance_A, 6),
            pair[1].atom.chain_id,
            pair[1].atom.residue_seq,
            pair[1].atom.atom_serial,
        ),
    ):
        if residue_key in first_shell_residue_keys:
            continue
        annotations.append(
            ShellAnnotation(
                structure_id=metal_atom.structure_id,
                metal_site_id=metal_atom.site_id,
                metal_element=metal_atom.element_upper,
                shell_type="second_sphere",
                chain_id=item.atom.chain_id,
                residue_name=item.atom.residue_name,
                residue_seq=item.atom.residue_seq,
                atom_name=item.atom.atom_name,
                atom_serial=item.atom.atom_serial,
                element=item.atom.element_upper,
                record_type=item.atom.record_type,
                distance_A=round(item.distance_A, 3),
                is_water=item.atom.residue_name in WATER_RESIDUES,
            )
        )
    donor_residue_count = len(first_shell_residue_keys)
    donor_distances = [item.distance_A for item in donor_matches]
    summary = MetalSiteSummary(
        structure_id=metal_atom.structure_id,
        metal_site_id=metal_atom.site_id,
        metal_element=metal_atom.element_upper,
        chain_id=metal_atom.chain_id,
        residue_name=metal_atom.residue_name,
        residue_seq=metal_atom.residue_seq,
        atom_serial=metal_atom.atom_serial,
        charge=metal_atom.charge,
        x=metal_atom.x,
        y=metal_atom.y,
        z=metal_atom.z,
        donor_atom_count=len(donor_matches),
        donor_residue_count=donor_residue_count,
        residue_within_6a_count=len(residue_contacts),
        nearest_donor_distance_A=round(min(donor_distances), 3) if donor_distances else None,
        farthest_donor_distance_A=round(max(donor_distances), 3) if donor_distances else None,
    )
    return summary, annotations

