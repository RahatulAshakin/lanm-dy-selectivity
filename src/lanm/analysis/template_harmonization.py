"""Phase 2A/2B template harmonization summaries and report generation."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import gemmi

from lanm.configuration import load_project_config
from lanm.data.fetch import load_structure_manifest, require_path
from lanm.filesystem import atomic_write_text, write_csv_rows, write_yaml
from lanm.models import (
    AtomRecord,
    CrossTemplateAlignmentRow,
    DesignMaskCandidate,
    DONOR_ELEMENTS,
    ResidueRoleAssignment,
)
from lanm.paths import REPO_ROOT
from lanm.structure.atoms import read_atom_records
from lanm.structure.geometry import WATER_RESIDUES, collect_atoms_within_cutoff, euclidean_distance
from lanm.structure.residues import ResidueRecord, collect_polymer_residues
from lanm.structure.templates import parse_cif_atom_records, read_cif_experimental_method
from lanm.viz.plots import render_template_harmonization_overview

_METAL_ELEMENTS = {"Y", "ND", "LA", "DY"}
_INTERCHAIN_CONTACT_CUTOFF_A = 5.0
_INTERFACE_NEIGHBORHOOD_WINDOW = 3
_AM1_TEMPLATE_IDS = frozenset({"6MI5", "8FNS"})
_HANS_TEMPLATE_IDS = frozenset({"8DQ2", "8FNR"})
_PROTECTED_RESIDUE_NAMES = frozenset({"GLY", "PRO", "CYS"})


@dataclass(frozen=True, slots=True)
class TemplateInput:
    template_id: str
    source_type: str
    path: Path


@dataclass(frozen=True, slots=True)
class TemplateChainSummaryRow:
    template_id: str
    source_type: str
    chain_id: str
    residue_start: str
    residue_end: str
    residue_range: str
    residue_count: int
    metal_identity: str
    metal_site_count: int
    experimental_method: str


@dataclass(frozen=True, slots=True)
class TemplateSiteSummaryRow:
    template_id: str
    source_type: str
    chain_id: str
    site_index: int
    site_label: str
    residue_name: str
    residue_seq: int
    insertion_code: str
    metal_identity: str
    experimental_method: str


@dataclass(frozen=True, slots=True)
class TemplateSummary:
    template_id: str
    source_type: str
    input_path: str
    chain_ids: tuple[str, ...]
    residue_ranges: tuple[str, ...]
    metal_identity: str
    metal_site_count: int
    experimental_method: str
    chain_rows: tuple[TemplateChainSummaryRow, ...]
    site_rows: tuple[TemplateSiteSummaryRow, ...]


@dataclass(frozen=True, slots=True)
class TemplateContext:
    template_id: str
    source_type: str
    input_path: str
    experimental_method: str
    atoms: tuple[AtomRecord, ...]
    chain_residues: dict[str, tuple[ResidueRecord, ...]]
    summary: TemplateSummary


@dataclass(frozen=True, slots=True)
class Am1MatureReference:
    sequence: str
    note: str


@dataclass(frozen=True, slots=True)
class AlignmentSummary:
    am1_reference_length: int
    am1_representative_template_id: str
    am1_representative_chain_id: str
    hans_representative_template_id: str
    hans_representative_chain_id: str
    canonical_family_position_count: int
    hans_only_canonical_position_count: int


@dataclass(frozen=True, slots=True)
class ResidueAlignmentMapping:
    canonical_family_position: int | None
    am1_mature_position: int | None
    alignment_status: str


@dataclass(frozen=True, slots=True)
class TemplateHarmonizationArtifacts:
    template_summaries: tuple[TemplateSummary, ...]
    chain_rows: tuple[TemplateChainSummaryRow, ...]
    site_rows: tuple[TemplateSiteSummaryRow, ...]
    cross_template_rows: tuple[CrossTemplateAlignmentRow, ...]
    residue_role_rows: tuple[ResidueRoleAssignment, ...]
    design_mask_rows: tuple[DesignMaskCandidate, ...]
    design_masks_config: dict[str, Any]
    template_specific_protected_positions: dict[str, dict[str, tuple[int, ...]]]
    report_markdown: str
    sequence_record_count: int


TEMPLATE_INPUTS = (
    TemplateInput("6MI5", "cif", Path("data/raw/public/structures/6MI5.cif")),
    TemplateInput("8FNS", "csv_atom_table", Path("data/raw/local_bundle/8fns_atoms.csv")),
    TemplateInput("8DQ2", "csv_atom_table", Path("data/raw/local_bundle/8dq2_atoms.csv")),
    TemplateInput("8FNR", "cif", Path("data/raw/public/structures/8FNR.cif")),
)

_FAMILY_COMPARISONS = (
    ("AM1/Mex family", ("6MI5", "8FNS")),
    ("Hans family", ("8DQ2", "8FNR")),
)
_ROLE_SORT_ORDER = {
    "first_shell": 0,
    "second_sphere": 1,
    "solvent_contact": 2,
    "interchain_contact": 3,
}
_TEMPLATE_SORT_ORDER = {item.template_id: index for index, item in enumerate(TEMPLATE_INPUTS)}


def _display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def load_sequence_rows(path: Path) -> list[dict[str, str]]:
    """Load the sequence reference table and fail with the exact missing path."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return sorted(
            ({key: str(value) for key, value in row.items()} for row in reader),
            key=lambda row: (
                row.get("organism", ""),
                row.get("strain", ""),
                row.get("gene_symbol", ""),
            ),
        )


def identify_am1_mature_sequence_reference(sequence_rows: list[dict[str, str]]) -> Am1MatureReference:
    """Identify the mature AM1 lanmodulin sequence used as the canonical reference."""
    for row in sequence_rows:
        strain = row.get("strain", "").strip().upper()
        sequence = row.get("sequence", "").strip()
        status = row.get("sequence_status", "").strip().lower()
        if strain == "AM1" and sequence and "mature" in status:
            return Am1MatureReference(sequence=sequence, note=row.get("notes", "").strip())
    raise ValueError("Unable to identify AM1 mature sequence reference")


def _load_manifest_method_map(manifest_path: Path) -> dict[str, str]:
    return {entry.pdb_id: entry.method for entry in load_structure_manifest(manifest_path)}


def _load_manifest_bound_metal_map(manifest_path: Path) -> dict[str, str]:
    return {entry.pdb_id: entry.bound_metal_symbol for entry in load_structure_manifest(manifest_path)}


def _format_residue_position(residue_seq: int, insertion_code: str) -> str:
    return f"{residue_seq}{insertion_code}" if insertion_code else str(residue_seq)


def _summarize_template(
    template_id: str,
    source_type: str,
    input_path: str,
    experimental_method: str,
    atoms: tuple[AtomRecord, ...],
) -> TemplateSummary:
    polymer_residues: dict[str, set[tuple[int, str]]] = defaultdict(set)
    metal_atoms = sorted(
        [
            atom for atom in atoms
            if atom.record_type == "HETATM" and atom.element_upper in _METAL_ELEMENTS
        ],
        key=lambda atom: (atom.chain_id, atom.residue_seq, atom.insertion_code, atom.atom_serial),
    )
    for atom in atoms:
        if atom.record_type != "ATOM":
            continue
        polymer_residues[atom.chain_id].add((atom.residue_seq, atom.insertion_code))
    metal_identities = sorted({atom.element_upper for atom in metal_atoms})
    chain_ids = tuple(sorted(polymer_residues))
    chain_site_counts: dict[str, int] = defaultdict(int)
    site_rows: list[TemplateSiteSummaryRow] = []
    for site_index, atom in enumerate(metal_atoms, start=1):
        chain_site_counts[atom.chain_id] += 1
        site_rows.append(
            TemplateSiteSummaryRow(
                template_id=template_id,
                source_type=source_type,
                chain_id=atom.chain_id,
                site_index=site_index,
                site_label=f"{atom.chain_id}:{_format_residue_position(atom.residue_seq, atom.insertion_code)}",
                residue_name=atom.residue_name,
                residue_seq=atom.residue_seq,
                insertion_code=atom.insertion_code,
                metal_identity=atom.element_upper,
                experimental_method=experimental_method,
            )
        )
    chain_rows: list[TemplateChainSummaryRow] = []
    residue_ranges: list[str] = []
    for chain_id in chain_ids:
        ordered_residues = sorted(polymer_residues[chain_id], key=lambda item: (item[0], item[1]))
        start_seq, start_insertion = ordered_residues[0]
        end_seq, end_insertion = ordered_residues[-1]
        residue_start = _format_residue_position(start_seq, start_insertion)
        residue_end = _format_residue_position(end_seq, end_insertion)
        residue_range = f"{chain_id}:{residue_start}-{residue_end}"
        residue_ranges.append(residue_range)
        chain_rows.append(
            TemplateChainSummaryRow(
                template_id=template_id,
                source_type=source_type,
                chain_id=chain_id,
                residue_start=residue_start,
                residue_end=residue_end,
                residue_range=residue_range,
                residue_count=len(ordered_residues),
                metal_identity=", ".join(metal_identities),
                metal_site_count=chain_site_counts.get(chain_id, 0),
                experimental_method=experimental_method,
            )
        )
    return TemplateSummary(
        template_id=template_id,
        source_type=source_type,
        input_path=input_path,
        chain_ids=chain_ids,
        residue_ranges=tuple(residue_ranges),
        metal_identity=", ".join(metal_identities),
        metal_site_count=len(site_rows),
        experimental_method=experimental_method,
        chain_rows=tuple(chain_rows),
        site_rows=tuple(site_rows),
    )


def _load_template_context(
    template_input: TemplateInput,
    manifest_methods: dict[str, str],
) -> TemplateContext:
    path = REPO_ROOT / template_input.path
    require_path(path)
    if template_input.source_type == "cif":
        atoms = tuple(parse_cif_atom_records(path, structure_id=template_input.template_id))
        experimental_method = read_cif_experimental_method(path)
    else:
        atoms = tuple(read_atom_records(path, structure_id=template_input.template_id))
        experimental_method = manifest_methods[template_input.template_id]
    input_path = _display_path(path)
    return TemplateContext(
        template_id=template_input.template_id,
        source_type=template_input.source_type,
        input_path=input_path,
        experimental_method=experimental_method,
        atoms=atoms,
        chain_residues=collect_polymer_residues(atoms),
        summary=_summarize_template(
            template_id=template_input.template_id,
            source_type=template_input.source_type,
            input_path=input_path,
            experimental_method=experimental_method,
            atoms=atoms,
        ),
    )


def _render_markdown_table(headers: tuple[str, ...], rows: list[tuple[str, ...]]) -> str:
    header_line = "| " + " | ".join(headers) + " |"
    divider_line = "| " + " | ".join(["---"] * len(headers)) + " |"
    body_lines = ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join([header_line, divider_line, *body_lines])


def _site_distribution(summary: TemplateSummary) -> str:
    return ", ".join(
        f"{row.chain_id}:{row.metal_site_count}"
        for row in summary.chain_rows
    )


def _alignment_scoring() -> gemmi.AlignmentScoring:
    scoring = gemmi.AlignmentScoring()
    scoring.match = 2
    scoring.mismatch = -1
    scoring.gapo = -10
    scoring.gape = -1
    scoring.good_gapo = 0
    scoring.bad_gapo = 0
    return scoring


def _align_sequences(query_sequence: str, target_sequence: str) -> tuple[str, str]:
    result = gemmi.align_string_sequences(
        list(query_sequence),
        list(target_sequence),
        [0] * len(target_sequence),
        _alignment_scoring(),
    )
    return result.add_gaps(query_sequence, 1), result.add_gaps(target_sequence, 2)


def _sequence_positions_for_residues(
    query_residues: tuple[ResidueRecord, ...],
    target_sequence: str,
) -> dict[int, int | None]:
    query_sequence = "".join(residue.one_letter_code for residue in query_residues)
    if target_sequence and target_sequence in query_sequence:
        start = query_sequence.index(target_sequence)
        return {
            index: (index - start if start < index <= start + len(target_sequence) else None)
            for index in range(1, len(query_residues) + 1)
        }
    if query_sequence and query_sequence in target_sequence:
        start = target_sequence.index(query_sequence)
        return {
            index: start + index
            for index in range(1, len(query_residues) + 1)
        }
    aligned_query, aligned_target = _align_sequences(query_sequence, target_sequence)
    query_index = 0
    target_index = 0
    mapping: dict[int, int | None] = {}
    for query_letter, target_letter in zip(aligned_query, aligned_target, strict=True):
        if query_letter != "-":
            query_index += 1
        if target_letter != "-":
            target_index += 1
        if query_letter != "-":
            mapping[query_index] = target_index if target_letter != "-" else None
    return mapping


def _alignment_columns(query_sequence: str, target_sequence: str) -> list[tuple[int | None, int | None]]:
    aligned_query, aligned_target = _align_sequences(query_sequence, target_sequence)
    query_index = 0
    target_index = 0
    columns: list[tuple[int | None, int | None]] = []
    for query_letter, target_letter in zip(aligned_query, aligned_target, strict=True):
        query_position: int | None = None
        target_position: int | None = None
        if query_letter != "-":
            query_index += 1
            query_position = query_index
        if target_letter != "-":
            target_index += 1
            target_position = target_index
        columns.append((query_position, target_position))
    return columns


def _classify_alignment_status(
    residue: ResidueRecord,
    am1_sequence: str,
    canonical_family_position: int | None,
    am1_mature_position: int | None,
) -> str:
    if canonical_family_position is None:
        return "outside_am1_mature_reference"
    if am1_mature_position is None:
        return "aligned_to_am1_gap"
    if residue.one_letter_code == am1_sequence[am1_mature_position - 1]:
        return "aligned_match"
    return "aligned_substitution"


def _build_cross_template_alignment_rows(
    template_contexts: tuple[TemplateContext, ...],
    am1_reference: Am1MatureReference,
) -> tuple[
    tuple[CrossTemplateAlignmentRow, ...],
    dict[tuple[str, str, int, str], ResidueAlignmentMapping],
    AlignmentSummary,
]:
    context_map = {context.template_id: context for context in template_contexts}
    am1_context = context_map["6MI5"]
    am1_chain_id = am1_context.summary.chain_ids[0]
    am1_residues = am1_context.chain_residues[am1_chain_id]
    am1_to_mature_position = _sequence_positions_for_residues(am1_residues, am1_reference.sequence)
    am1_family_residues = tuple(
        residue
        for index, residue in enumerate(am1_residues, start=1)
        if am1_to_mature_position[index] is not None
    )

    hans_context = context_map["8DQ2"]
    hans_chain_id = sorted(
        hans_context.chain_residues,
        key=lambda chain_id: (-len(hans_context.chain_residues[chain_id]), chain_id),
    )[0]
    hans_residues = hans_context.chain_residues[hans_chain_id]

    canonical_columns = _alignment_columns(
        "".join(residue.one_letter_code for residue in hans_residues),
        "".join(residue.one_letter_code for residue in am1_family_residues),
    )
    am1_position_to_canonical: dict[int, int] = {}
    hans_position_to_canonical: dict[int, int] = {}
    canonical_to_am1_position: dict[int, int | None] = {}
    for canonical_position, (hans_position, am1_position) in enumerate(canonical_columns, start=1):
        canonical_to_am1_position[canonical_position] = am1_position
        if am1_position is not None:
            am1_position_to_canonical[am1_position] = canonical_position
        if hans_position is not None:
            hans_position_to_canonical[hans_position] = canonical_position

    alignment_lookup: dict[tuple[str, str, int, str], ResidueAlignmentMapping] = {}
    rows: list[CrossTemplateAlignmentRow] = []
    for context in template_contexts:
        if context.template_id == "6MI5":
            for index, residue in enumerate(context.chain_residues[am1_chain_id], start=1):
                am1_mature_position = am1_to_mature_position[index]
                canonical_family_position = (
                    am1_position_to_canonical[am1_mature_position]
                    if am1_mature_position is not None
                    else None
                )
                alignment_status = _classify_alignment_status(
                    residue=residue,
                    am1_sequence=am1_reference.sequence,
                    canonical_family_position=canonical_family_position,
                    am1_mature_position=am1_mature_position,
                )
                row = CrossTemplateAlignmentRow(
                    template_id=context.template_id,
                    chain_id=residue.chain_id,
                    template_residue_seq=residue.residue_seq,
                    template_residue_name=residue.residue_name,
                    canonical_family_position=canonical_family_position,
                    am1_mature_position=am1_mature_position,
                    alignment_status=alignment_status,
                )
                rows.append(row)
                alignment_lookup[(context.template_id, residue.chain_id, residue.residue_seq, residue.insertion_code)] = ResidueAlignmentMapping(
                    canonical_family_position=canonical_family_position,
                    am1_mature_position=am1_mature_position,
                    alignment_status=alignment_status,
                )
            continue

        if context.template_id == "8FNS":
            target_positions = _sequence_positions_for_residues(
                context.chain_residues["A"],
                "".join(residue.one_letter_code for residue in am1_family_residues),
            )
            for index, residue in enumerate(context.chain_residues["A"], start=1):
                am1_mature_position = target_positions[index]
                canonical_family_position = (
                    am1_position_to_canonical[am1_mature_position]
                    if am1_mature_position is not None
                    else None
                )
                alignment_status = _classify_alignment_status(
                    residue=residue,
                    am1_sequence=am1_reference.sequence,
                    canonical_family_position=canonical_family_position,
                    am1_mature_position=am1_mature_position,
                )
                row = CrossTemplateAlignmentRow(
                    template_id=context.template_id,
                    chain_id=residue.chain_id,
                    template_residue_seq=residue.residue_seq,
                    template_residue_name=residue.residue_name,
                    canonical_family_position=canonical_family_position,
                    am1_mature_position=am1_mature_position,
                    alignment_status=alignment_status,
                )
                rows.append(row)
                alignment_lookup[(context.template_id, residue.chain_id, residue.residue_seq, residue.insertion_code)] = ResidueAlignmentMapping(
                    canonical_family_position=canonical_family_position,
                    am1_mature_position=am1_mature_position,
                    alignment_status=alignment_status,
                )
            continue

        for chain_id in sorted(context.chain_residues):
            chain_residues = context.chain_residues[chain_id]
            if context.template_id == "8DQ2" and chain_id == hans_chain_id:
                query_to_hans = {index: index for index in range(1, len(chain_residues) + 1)}
            else:
                query_to_hans = _sequence_positions_for_residues(
                    chain_residues,
                    "".join(residue.one_letter_code for residue in hans_residues),
                )
            for index, residue in enumerate(chain_residues, start=1):
                hans_position = query_to_hans[index]
                canonical_family_position = (
                    hans_position_to_canonical[hans_position]
                    if hans_position is not None
                    else None
                )
                am1_mature_position = (
                    canonical_to_am1_position[canonical_family_position]
                    if canonical_family_position is not None
                    else None
                )
                alignment_status = _classify_alignment_status(
                    residue=residue,
                    am1_sequence=am1_reference.sequence,
                    canonical_family_position=canonical_family_position,
                    am1_mature_position=am1_mature_position,
                )
                row = CrossTemplateAlignmentRow(
                    template_id=context.template_id,
                    chain_id=residue.chain_id,
                    template_residue_seq=residue.residue_seq,
                    template_residue_name=residue.residue_name,
                    canonical_family_position=canonical_family_position,
                    am1_mature_position=am1_mature_position,
                    alignment_status=alignment_status,
                )
                rows.append(row)
                alignment_lookup[(context.template_id, residue.chain_id, residue.residue_seq, residue.insertion_code)] = ResidueAlignmentMapping(
                    canonical_family_position=canonical_family_position,
                    am1_mature_position=am1_mature_position,
                    alignment_status=alignment_status,
                )

    ordered_rows = tuple(
        sorted(
            rows,
            key=lambda row: (
                _TEMPLATE_SORT_ORDER[row.template_id],
                row.chain_id,
                row.template_residue_seq,
                row.template_residue_name,
            ),
        )
    )
    summary = AlignmentSummary(
        am1_reference_length=len(am1_reference.sequence),
        am1_representative_template_id="6MI5",
        am1_representative_chain_id=am1_chain_id,
        hans_representative_template_id=hans_context.template_id,
        hans_representative_chain_id=hans_chain_id,
        canonical_family_position_count=len(canonical_columns),
        hans_only_canonical_position_count=sum(
            1 for am1_position in canonical_to_am1_position.values() if am1_position is None
        ),
    )
    return ordered_rows, alignment_lookup, summary


def _select_bound_metal_atoms(atoms: tuple[AtomRecord, ...], expected_element: str) -> list[AtomRecord]:
    return sorted(
        [
            atom for atom in atoms
            if atom.record_type == "HETATM" and atom.element_upper == expected_element
        ],
        key=lambda atom: (atom.chain_id, atom.residue_seq, atom.insertion_code, atom.atom_serial),
    )


def _lookup_alignment_mapping(
    template_id: str,
    atom: AtomRecord,
    alignment_lookup: dict[tuple[str, str, int, str], ResidueAlignmentMapping],
) -> ResidueAlignmentMapping | None:
    return alignment_lookup.get((template_id, atom.chain_id, atom.residue_seq, atom.insertion_code))


def _build_site_local_role_rows(
    context: TemplateContext,
    expected_element: str,
    alignment_lookup: dict[tuple[str, str, int, str], ResidueAlignmentMapping],
    first_shell_cutoff_A: float,
    second_sphere_cutoff_A: float,
) -> list[ResidueRoleAssignment]:
    chain_site_counts: dict[str, int] = defaultdict(int)
    role_rows: list[ResidueRoleAssignment] = []
    for metal_atom in _select_bound_metal_atoms(context.atoms, expected_element):
        chain_site_counts[metal_atom.chain_id] += 1
        site_index = chain_site_counts[metal_atom.chain_id]
        site_label = f"{metal_atom.chain_id}:{_format_residue_position(metal_atom.residue_seq, metal_atom.insertion_code)}"
        nearby_atoms = collect_atoms_within_cutoff(metal_atom, context.atoms, second_sphere_cutoff_A)
        residue_matches: dict[tuple[str, str, int, str], list] = defaultdict(list)
        for item in nearby_atoms:
            if item.atom.residue_key == metal_atom.residue_key and item.atom.chain_id == metal_atom.chain_id:
                continue
            residue_matches[item.atom.residue_key].append(item)
        for residue_key, items in sorted(
            residue_matches.items(),
            key=lambda pair: (
                round(min(item.distance_A for item in pair[1]), 6),
                pair[1][0].atom.chain_id,
                pair[1][0].atom.residue_seq,
                pair[1][0].atom.atom_serial,
            ),
        ):
            representative_item = min(
                items,
                key=lambda item: (
                    round(item.distance_A, 6),
                    item.atom.chain_id,
                    item.atom.residue_seq,
                    item.atom.atom_serial,
                ),
            )
            donor_items = [
                item for item in items
                if item.atom.element_upper in DONOR_ELEMENTS and item.distance_A <= first_shell_cutoff_A
            ]
            residue_distance = min(item.distance_A for item in items)
            residue_atom = representative_item.atom
            is_nonpolymer = residue_atom.record_type != "ATOM" and not residue_atom.is_metal
            if is_nonpolymer:
                role = "solvent_contact"
                distance_A = residue_distance
                alignment = None
                if residue_atom.residue_name in WATER_RESIDUES:
                    note = "water_or_buffer"
                else:
                    note = "hetero_solvent_contact"
            elif donor_items:
                role = "first_shell"
                distance_A = min(item.distance_A for item in donor_items)
                alignment = _lookup_alignment_mapping(context.template_id, residue_atom, alignment_lookup)
                note = f"donor_atoms={','.join(sorted({item.atom.atom_name for item in donor_items}))}"
            elif residue_atom.record_type == "ATOM" and residue_distance > first_shell_cutoff_A:
                role = "second_sphere"
                distance_A = residue_distance
                alignment = _lookup_alignment_mapping(context.template_id, residue_atom, alignment_lookup)
                note = ""
            else:
                continue
            role_rows.append(
                ResidueRoleAssignment(
                    template_id=context.template_id,
                    chain_id=residue_atom.chain_id,
                    site_index=site_index,
                    site_label=site_label,
                    role=role,
                    residue_id=f"{residue_atom.chain_id}:{_format_residue_position(residue_atom.residue_seq, residue_atom.insertion_code)}",
                    residue_name=residue_atom.residue_name,
                    residue_seq=residue_atom.residue_seq,
                    distance_A=round(distance_A, 3),
                    canonical_family_position=alignment.canonical_family_position if alignment else None,
                    am1_mature_position=alignment.am1_mature_position if alignment else None,
                    alignment_status=alignment.alignment_status if alignment else "non_polymer",
                    note=note,
                )
            )
    return role_rows


def _build_interchain_role_rows(
    context: TemplateContext,
    alignment_lookup: dict[tuple[str, str, int, str], ResidueAlignmentMapping],
) -> list[ResidueRoleAssignment]:
    polymer_atoms = [atom for atom in context.atoms if atom.record_type == "ATOM"]
    best_distances: dict[tuple[str, str, int, str], float] = {}
    representative_atoms: dict[tuple[str, str, int, str], AtomRecord] = {}
    partner_chains: dict[tuple[str, str, int, str], set[str]] = defaultdict(set)
    for index, atom in enumerate(polymer_atoms):
        for other in polymer_atoms[index + 1:]:
            if atom.chain_id == other.chain_id:
                continue
            distance_A = euclidean_distance(atom, other)
            if distance_A > _INTERCHAIN_CONTACT_CUTOFF_A:
                continue
            for source, partner in ((atom, other), (other, atom)):
                residue_key = source.residue_key
                current = best_distances.get(residue_key)
                if current is None or distance_A < current:
                    best_distances[residue_key] = distance_A
                    representative_atoms[residue_key] = source
                partner_chains[residue_key].add(partner.chain_id)
    role_rows: list[ResidueRoleAssignment] = []
    for residue_key in sorted(
        representative_atoms,
        key=lambda key: (
            key[0],
            key[2],
            key[3],
            key[1],
        ),
    ):
        atom = representative_atoms[residue_key]
        alignment = _lookup_alignment_mapping(context.template_id, atom, alignment_lookup)
        role_rows.append(
            ResidueRoleAssignment(
                template_id=context.template_id,
                chain_id=atom.chain_id,
                site_index=None,
                site_label="",
                role="interchain_contact",
                residue_id=f"{atom.chain_id}:{_format_residue_position(atom.residue_seq, atom.insertion_code)}",
                residue_name=atom.residue_name,
                residue_seq=atom.residue_seq,
                distance_A=round(best_distances[residue_key], 3),
                canonical_family_position=alignment.canonical_family_position if alignment else None,
                am1_mature_position=alignment.am1_mature_position if alignment else None,
                alignment_status=alignment.alignment_status if alignment else "unmapped",
                note=f"partner_chains={','.join(sorted(partner_chains[residue_key]))}",
            )
        )
    return role_rows


def _build_residue_role_rows(
    template_contexts: tuple[TemplateContext, ...],
    manifest_metals: dict[str, str],
    alignment_lookup: dict[tuple[str, str, int, str], ResidueAlignmentMapping],
    first_shell_cutoff_A: float,
    second_sphere_cutoff_A: float,
) -> tuple[ResidueRoleAssignment, ...]:
    rows: list[ResidueRoleAssignment] = []
    for context in template_contexts:
        expected_element = manifest_metals[context.template_id]
        rows.extend(
            _build_site_local_role_rows(
                context=context,
                expected_element=expected_element,
                alignment_lookup=alignment_lookup,
                first_shell_cutoff_A=first_shell_cutoff_A,
                second_sphere_cutoff_A=second_sphere_cutoff_A,
            )
        )
        rows.extend(_build_interchain_role_rows(context=context, alignment_lookup=alignment_lookup))
    return tuple(
        sorted(
            rows,
            key=lambda row: (
                _TEMPLATE_SORT_ORDER[row.template_id],
                row.chain_id,
                row.site_index is None,
                row.site_index if row.site_index is not None else 999,
                _ROLE_SORT_ORDER[row.role],
                row.residue_seq,
                row.residue_name,
                row.note,
            ),
        )
    )


def _summarize_role_observations(
    rows: list[ResidueRoleAssignment],
    *,
    label: str,
) -> str:
    template_ids = ", ".join(sorted({row.template_id for row in rows}))
    distances = [row.distance_A for row in rows if row.distance_A is not None]
    if distances:
        return f"{label} observed {len(rows)}x across {template_ids} (closest {min(distances):.3f} A)"
    return f"{label} observed {len(rows)}x across {template_ids}"


def _build_design_mask_rationale(
    *,
    first_shell_rows: list[ResidueRoleAssignment],
    second_sphere_rows: list[ResidueRoleAssignment],
    hans_interface_rows: list[ResidueRoleAssignment],
    fixed_first_shell: bool,
    mutable_second_sphere: bool,
    mutable_interface: bool,
    interface_neighborhood: bool,
    protection_reasons: list[str],
) -> str:
    fragments: list[str] = []
    if fixed_first_shell:
        fragments.append(
            "fixed_first_shell because "
            + _summarize_role_observations(first_shell_rows, label="first-shell donor")
        )
    if mutable_second_sphere:
        fragments.append(
            "mutable_second_sphere because "
            + _summarize_role_observations(second_sphere_rows, label="second-sphere contact")
        )
    elif second_sphere_rows and not fixed_first_shell:
        fragments.append(
            "second-sphere evidence retained but excluded from mutable_second_sphere"
        )
    if mutable_interface:
        fragments.append(
            "mutable_interface because "
            + _summarize_role_observations(hans_interface_rows, label="Hans interchain contact")
            + f" and the position lies within +/-{_INTERFACE_NEIGHBORHOOD_WINDOW} canonical positions "
            "of a Hans metal-associated residue"
        )
    elif hans_interface_rows and not interface_neighborhood:
        fragments.append(
            f"Hans interchain contact observed but outside the +/-{_INTERFACE_NEIGHBORHOOD_WINDOW}-position "
            "metal-site neighborhood filter"
        )
    if protection_reasons:
        fragments.append("protected_positions because " + "; ".join(protection_reasons))
    return "; ".join(fragments)


def _build_design_mask_candidates(
    cross_template_rows: tuple[CrossTemplateAlignmentRow, ...],
    residue_role_rows: tuple[ResidueRoleAssignment, ...],
) -> tuple[DesignMaskCandidate, ...]:
    cross_rows_by_position: dict[int, list[CrossTemplateAlignmentRow]] = defaultdict(list)
    for row in cross_template_rows:
        if row.canonical_family_position is None:
            continue
        cross_rows_by_position[row.canonical_family_position].append(row)

    role_rows_by_key: dict[tuple[int, str], list[ResidueRoleAssignment]] = defaultdict(list)
    hans_metal_associated_positions: set[int] = set()
    for row in residue_role_rows:
        if row.canonical_family_position is None:
            continue
        role_rows_by_key[(row.canonical_family_position, row.role)].append(row)
        if row.template_id in _HANS_TEMPLATE_IDS and row.role in {"first_shell", "second_sphere"}:
            hans_metal_associated_positions.add(row.canonical_family_position)

    candidates: list[DesignMaskCandidate] = []
    for canonical_family_position in sorted(cross_rows_by_position):
        alignment_rows = cross_rows_by_position[canonical_family_position]
        first_shell_rows = role_rows_by_key.get((canonical_family_position, "first_shell"), [])
        second_sphere_rows = role_rows_by_key.get((canonical_family_position, "second_sphere"), [])
        hans_interface_rows = [
            row
            for row in role_rows_by_key.get((canonical_family_position, "interchain_contact"), [])
            if row.template_id in _HANS_TEMPLATE_IDS
        ]
        am1_positions = sorted(
            {
                row.am1_mature_position
                for row in alignment_rows
                if row.am1_mature_position is not None
            }
        )
        am1_reference_row = next(
            (
                row
                for row in alignment_rows
                if row.template_id == "6MI5" and row.am1_mature_position is not None
            ),
            None,
        )
        am1_reference_residue = am1_reference_row.template_residue_name if am1_reference_row else ""
        observed_residue_identities = sorted(
            {row.template_residue_name for row in alignment_rows}
        )
        alignment_statuses = {row.alignment_status for row in alignment_rows}
        am1_family_identities = {
            row.template_residue_name
            for row in alignment_rows
            if row.template_id in _AM1_TEMPLATE_IDS
        }
        hans_family_identities = {
            row.template_residue_name
            for row in alignment_rows
            if row.template_id in _HANS_TEMPLATE_IDS
        }

        protection_reasons: list[str] = []
        am1_mature_position = am1_positions[0] if len(am1_positions) == 1 else None
        if len(am1_positions) > 1:
            protection_reasons.append("ambiguous AM1 mature mapping")
        elif am1_mature_position is None:
            protection_reasons.append("outside AM1 mature numbering")
        if "aligned_to_am1_gap" in alignment_statuses:
            protection_reasons.append("unresolved or gap-only alignment")
        if am1_reference_residue in _PROTECTED_RESIDUE_NAMES:
            protection_reasons.append(f"AM1 reference residue is {am1_reference_residue}")
        if len(am1_family_identities) > 1 or len(hans_family_identities) > 1:
            protection_reasons.append("within-family residue identity conflict")

        interface_neighborhood = bool(hans_interface_rows) and any(
            abs(canonical_family_position - source_position) <= _INTERFACE_NEIGHBORHOOD_WINDOW
            for source_position in hans_metal_associated_positions
        )
        fixed_first_shell = bool(first_shell_rows)
        protected_positions = bool(protection_reasons)
        mutable_second_sphere = bool(second_sphere_rows) and not fixed_first_shell and not protected_positions
        mutable_interface = interface_neighborhood and not fixed_first_shell and not protected_positions

        if not any((fixed_first_shell, mutable_second_sphere, mutable_interface, protected_positions)):
            continue

        candidates.append(
            DesignMaskCandidate(
                canonical_family_position=canonical_family_position,
                am1_mature_position=am1_mature_position,
                am1_reference_residue=am1_reference_residue,
                observed_residue_identities=",".join(observed_residue_identities),
                template_coverage_count=len(alignment_rows),
                first_shell_observation_count=len(first_shell_rows),
                second_sphere_observation_count=len(second_sphere_rows),
                hans_interface_observation_count=len(hans_interface_rows),
                interface_neighborhood=interface_neighborhood,
                fixed_first_shell=fixed_first_shell,
                mutable_second_sphere=mutable_second_sphere,
                mutable_interface=mutable_interface,
                protected_positions=protected_positions,
                protection_reasons="; ".join(protection_reasons),
                rationale=_build_design_mask_rationale(
                    first_shell_rows=first_shell_rows,
                    second_sphere_rows=second_sphere_rows,
                    hans_interface_rows=hans_interface_rows,
                    fixed_first_shell=fixed_first_shell,
                    mutable_second_sphere=mutable_second_sphere,
                    mutable_interface=mutable_interface,
                    interface_neighborhood=interface_neighborhood,
                    protection_reasons=protection_reasons,
                ),
            )
        )
    return tuple(candidates)


def _collect_template_specific_protected_positions(
    cross_template_rows: tuple[CrossTemplateAlignmentRow, ...],
) -> dict[str, dict[str, tuple[int, ...]]]:
    protected_positions: dict[str, dict[str, list[int]]] = defaultdict(lambda: defaultdict(list))
    for row in cross_template_rows:
        if (
            row.template_id == "6MI5"
            and row.alignment_status == "outside_am1_mature_reference"
            and row.template_residue_name == "HIS"
        ):
            protected_positions[row.template_id][row.chain_id].append(row.template_residue_seq)
    return {
        template_id: {
            chain_id: tuple(sorted(residue_positions))
            for chain_id, residue_positions in sorted(chain_map.items())
        }
        for template_id, chain_map in sorted(protected_positions.items())
    }


def _mask_positions(
    design_mask_rows: tuple[DesignMaskCandidate, ...],
    attribute_name: str,
) -> list[int]:
    return [
        row.canonical_family_position
        for row in design_mask_rows
        if getattr(row, attribute_name)
    ]


def _mask_am1_positions(
    design_mask_rows: tuple[DesignMaskCandidate, ...],
    attribute_name: str,
) -> list[int]:
    return sorted(
        {
            row.am1_mature_position
            for row in design_mask_rows
            if getattr(row, attribute_name) and row.am1_mature_position is not None
        }
    )


def _build_design_masks_config(
    design_mask_rows: tuple[DesignMaskCandidate, ...],
    template_specific_protected_positions: dict[str, dict[str, tuple[int, ...]]],
) -> dict[str, Any]:
    noncanonical_positions = {
        template_id: {
            chain_id: list(positions)
            for chain_id, positions in chain_map.items()
        }
        for template_id, chain_map in template_specific_protected_positions.items()
    }
    return {
        "version": 1,
        "index": {
            "primary": "canonical_family_position",
            "am1_reference": "mature AM1 numbering",
            "interface_neighborhood_window": _INTERFACE_NEIGHBORHOOD_WINDOW,
        },
        "masks": {
            "fixed_first_shell": {
                "canonical_family_positions": _mask_positions(design_mask_rows, "fixed_first_shell"),
                "am1_mature_positions": _mask_am1_positions(design_mask_rows, "fixed_first_shell"),
            },
            "mutable_second_sphere": {
                "canonical_family_positions": _mask_positions(design_mask_rows, "mutable_second_sphere"),
                "am1_mature_positions": _mask_am1_positions(design_mask_rows, "mutable_second_sphere"),
            },
            "mutable_interface": {
                "canonical_family_positions": _mask_positions(design_mask_rows, "mutable_interface"),
                "am1_mature_positions": _mask_am1_positions(design_mask_rows, "mutable_interface"),
            },
            "protected_positions": {
                "canonical_family_positions": _mask_positions(design_mask_rows, "protected_positions"),
                "am1_mature_positions": _mask_am1_positions(design_mask_rows, "protected_positions"),
                "template_specific_noncanonical": noncanonical_positions,
            },
        },
    }


def _candidate_label(row: DesignMaskCandidate) -> str:
    if row.am1_mature_position is None:
        return f"{row.canonical_family_position} (gap)"
    return f"{row.canonical_family_position} / AM1 {row.am1_mature_position}"


def _format_residue_positions(positions: tuple[int, ...]) -> str:
    if not positions:
        return ""
    if len(positions) > 1 and positions[-1] - positions[0] + 1 == len(positions):
        return f"{positions[0]}-{positions[-1]}"
    return ", ".join(str(position) for position in positions)


def _top_mask_rows(
    design_mask_rows: tuple[DesignMaskCandidate, ...],
    *,
    attribute_name: str,
    observation_count_field: str,
    limit: int,
) -> list[DesignMaskCandidate]:
    return sorted(
        [row for row in design_mask_rows if getattr(row, attribute_name)],
        key=lambda row: (
            -getattr(row, observation_count_field),
            row.am1_mature_position is None,
            row.am1_mature_position if row.am1_mature_position is not None else 999,
            row.canonical_family_position,
        ),
    )[:limit]


def render_template_harmonization_markdown(
    template_summaries: tuple[TemplateSummary, ...],
    sequence_rows: list[dict[str, str]],
    cross_template_rows: tuple[CrossTemplateAlignmentRow, ...],
    residue_role_rows: tuple[ResidueRoleAssignment, ...],
    design_mask_rows: tuple[DesignMaskCandidate, ...],
    template_specific_protected_positions: dict[str, dict[str, tuple[int, ...]]],
    alignment_summary: AlignmentSummary,
) -> str:
    summary_map = {summary.template_id: summary for summary in template_summaries}
    lines = [
        "# Template Harmonization",
        "",
        "Phase 2A/2B/2C deterministic template parsing, cross-template alignment, residue-role summary, and design-mask generation for four lanmodulin templates.",
        "",
        f"Sequence reference records loaded: {len(sequence_rows)}",
        "",
    ]
    am1_note = identify_am1_mature_sequence_reference(sequence_rows).note
    if am1_note:
        lines.extend(
            [
                "AM1 sequence note:",
                "",
                f"- {am1_note}",
                "",
            ]
        )
    for family_label, template_ids in _FAMILY_COMPARISONS:
        rows: list[tuple[str, ...]] = []
        for template_id in template_ids:
            summary = summary_map[template_id]
            rows.append(
                (
                    summary.template_id,
                    summary.source_type,
                    ", ".join(summary.chain_ids),
                    "; ".join(summary.residue_ranges),
                    summary.metal_identity,
                    str(summary.metal_site_count),
                    summary.experimental_method,
                )
            )
        lines.extend(
            [
                f"## {family_label}",
                "",
                _render_markdown_table(
                    (
                        "template_id",
                        "source_type",
                        "chain_ids",
                        "residue_ranges",
                        "metal_identity",
                        "metal_site_count",
                        "experimental_method",
                    ),
                    rows,
                ),
                "",
            ]
        )
        if family_label == "AM1/Mex family":
            lines.extend(
                [
                    f"- `6MI5` uses model 1 of the NMR ensemble for deterministic counting and reports {_site_distribution(summary_map['6MI5'])} metal sites.",
                    f"- `8FNS` is a single-chain crystal template spanning {'; '.join(summary_map['8FNS'].residue_ranges)} with {_site_distribution(summary_map['8FNS'])} metal sites.",
                    "",
                ]
            )
        else:
            lines.extend(
                [
                    f"- `8DQ2` and `8FNR` both resolve chains {', '.join(summary_map['8DQ2'].chain_ids)} over the Hans family backbone.",
                    f"- Metal occupancy differs across the tetramer: `8DQ2` shows {_site_distribution(summary_map['8DQ2'])}, while `8FNR` shows {_site_distribution(summary_map['8FNR'])}.",
                    "- `8FNR` chain D retains residues 24-133 but the observed model omits residues 34-38, reducing the chain residue count relative to `8DQ2`.",
                    "",
                ]
            )

    fns_positions = sorted(
        row.am1_mature_position
        for row in cross_template_rows
        if row.template_id == "8FNS" and row.am1_mature_position is not None
    )
    his_tail_count = sum(
        1
        for row in cross_template_rows
        if row.template_id == "6MI5" and row.alignment_status == "outside_am1_mature_reference"
    )
    role_counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in residue_role_rows:
        role_counts[row.template_id][row.role] += 1
    alignment_counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in cross_template_rows:
        alignment_counts[row.template_id][row.alignment_status] += 1
    phase_2b_rows = [
        (
            template_id,
            str(sum(alignment_counts[template_id].values())),
            str(alignment_counts[template_id]["outside_am1_mature_reference"]),
            str(role_counts[template_id]["first_shell"]),
            str(role_counts[template_id]["second_sphere"]),
            str(role_counts[template_id]["solvent_contact"]),
            str(role_counts[template_id]["interchain_contact"]),
        )
        for template_id in [item.template_id for item in TEMPLATE_INPUTS]
    ]
    lines.extend(
        [
            "## Alignment and residue roles",
            "",
            f"- AM1 mature reference length: {alignment_summary.am1_reference_length} aa from `lanmodulin_sequences.csv`.",
            f"- `6MI5` chain `{alignment_summary.am1_representative_chain_id}` maps residues 23-133 onto AM1 mature positions 1-111; {his_tail_count} C-terminal His-tag residues remain outside the mature reference.",
            f"- `8FNS` chain `A` covers AM1 mature positions {fns_positions[0]}-{fns_positions[-1]} in the representative AM1 alignment.",
            f"- Hans representative chain: `{alignment_summary.hans_representative_template_id}` chain `{alignment_summary.hans_representative_chain_id}` aligned against `6MI5` chain `{alignment_summary.am1_representative_chain_id}` to define {alignment_summary.canonical_family_position_count} canonical family positions; {alignment_summary.hans_only_canonical_position_count} positions are Hans-only relative to AM1 mature numbering.",
            "",
            _render_markdown_table(
                (
                    "template_id",
                    "aligned_rows",
                    "outside_am1_reference",
                    "first_shell",
                    "second_sphere",
                    "solvent_contact",
                    "interchain_contact",
                ),
                phase_2b_rows,
            ),
            "",
        ]
    )
    mask_count_rows = [
        ("fixed_first_shell", str(sum(1 for row in design_mask_rows if row.fixed_first_shell))),
        ("mutable_second_sphere", str(sum(1 for row in design_mask_rows if row.mutable_second_sphere))),
        ("mutable_interface", str(sum(1 for row in design_mask_rows if row.mutable_interface))),
        (
            "protected_positions",
            str(sum(1 for row in design_mask_rows if row.protected_positions)),
        ),
    ]
    top_second_rows = [
        (
            _candidate_label(row),
            row.am1_reference_residue or "-",
            str(row.second_sphere_observation_count),
            row.observed_residue_identities,
            row.rationale,
        )
        for row in _top_mask_rows(
            design_mask_rows,
            attribute_name="mutable_second_sphere",
            observation_count_field="second_sphere_observation_count",
            limit=8,
        )
    ]
    top_interface_rows = [
        (
            _candidate_label(row),
            row.am1_reference_residue or "-",
            str(row.hans_interface_observation_count),
            row.observed_residue_identities,
            row.rationale,
        )
        for row in _top_mask_rows(
            design_mask_rows,
            attribute_name="mutable_interface",
            observation_count_field="hans_interface_observation_count",
            limit=8,
        )
    ]
    protected_rows = [
        (
            _candidate_label(row),
            row.am1_reference_residue or "-",
            row.observed_residue_identities,
            row.protection_reasons,
        )
        for row in sorted(
            [row for row in design_mask_rows if row.protected_positions],
            key=lambda row: (
                row.am1_mature_position is None,
                row.am1_mature_position if row.am1_mature_position is not None else 999,
                row.canonical_family_position,
            ),
        )
    ]
    lines.extend(
        [
            "## Design mask candidates",
            "",
            _render_markdown_table(
                ("category", "canonical_position_count"),
                mask_count_rows,
            ),
            "",
            f"- `fixed_first_shell` spans {_mask_positions(design_mask_rows, 'fixed_first_shell')}.",
            f"- `mutable_second_sphere` uses the top-ranked positions below after excluding fixed and protected sites.",
            f"- `mutable_interface` applies the Hans-family interchain-contact filter within +/-{_INTERFACE_NEIGHBORHOOD_WINDOW} canonical positions of metal-associated residues.",
            "",
            "### Top mutable_second_sphere positions",
            "",
            _render_markdown_table(
                (
                    "canonical_family_position",
                    "am1_reference_residue",
                    "second_sphere_obs",
                    "observed_residue_identities",
                    "rationale",
                ),
                top_second_rows,
            ) if top_second_rows else "_No mutable_second_sphere positions passed the Phase 2C filter._",
            "",
            "### Top mutable_interface positions",
            "",
            _render_markdown_table(
                (
                    "canonical_family_position",
                    "am1_reference_residue",
                    "interface_obs",
                    "observed_residue_identities",
                    "rationale",
                ),
                top_interface_rows,
            ) if top_interface_rows else "_No mutable_interface positions passed the Phase 2C filter._",
            "",
            "### Positions explicitly protected",
            "",
            _render_markdown_table(
                (
                    "canonical_family_position",
                    "am1_reference_residue",
                    "observed_residue_identities",
                    "protection_reasons",
                ),
                protected_rows,
            ),
            "",
        ]
    )
    if template_specific_protected_positions:
        noncanonical_notes = [
            f"`{template_id}` chain `{chain_id}` residues {_format_residue_positions(positions)}"
            for template_id, chain_map in template_specific_protected_positions.items()
            for chain_id, positions in chain_map.items()
        ]
        lines.extend(
            [
                "- Additional noncanonical protected positions: " + "; ".join(noncanonical_notes) + ".",
                "",
            ]
        )
    return "\n".join(lines).rstrip() + "\n"


def build_template_harmonization_artifacts(
    *,
    sequences_path: Path,
    manifest_path: Path,
    template_inputs: tuple[TemplateInput, ...] = TEMPLATE_INPUTS,
) -> TemplateHarmonizationArtifacts:
    sequence_rows = load_sequence_rows(sequences_path)
    am1_reference = identify_am1_mature_sequence_reference(sequence_rows)
    manifest_methods = _load_manifest_method_map(manifest_path)
    manifest_metals = _load_manifest_bound_metal_map(manifest_path)
    config = load_project_config()
    template_contexts = tuple(
        _load_template_context(template_input, manifest_methods)
        for template_input in template_inputs
    )
    template_summaries = tuple(context.summary for context in template_contexts)
    chain_rows = tuple(
        row
        for summary in template_summaries
        for row in summary.chain_rows
    )
    site_rows = tuple(
        row
        for summary in template_summaries
        for row in summary.site_rows
    )
    cross_template_rows, alignment_lookup, alignment_summary = _build_cross_template_alignment_rows(
        template_contexts=template_contexts,
        am1_reference=am1_reference,
    )
    residue_role_rows = _build_residue_role_rows(
        template_contexts=template_contexts,
        manifest_metals=manifest_metals,
        alignment_lookup=alignment_lookup,
        first_shell_cutoff_A=config.first_shell_cutoff_A,
        second_sphere_cutoff_A=config.second_sphere_cutoff_A,
    )
    design_mask_rows = _build_design_mask_candidates(
        cross_template_rows=cross_template_rows,
        residue_role_rows=residue_role_rows,
    )
    template_specific_protected_positions = _collect_template_specific_protected_positions(
        cross_template_rows
    )
    design_masks_config = _build_design_masks_config(
        design_mask_rows=design_mask_rows,
        template_specific_protected_positions=template_specific_protected_positions,
    )
    return TemplateHarmonizationArtifacts(
        template_summaries=template_summaries,
        chain_rows=chain_rows,
        site_rows=site_rows,
        cross_template_rows=cross_template_rows,
        residue_role_rows=residue_role_rows,
        design_mask_rows=design_mask_rows,
        design_masks_config=design_masks_config,
        template_specific_protected_positions=template_specific_protected_positions,
        report_markdown=render_template_harmonization_markdown(
            template_summaries=template_summaries,
            sequence_rows=sequence_rows,
            cross_template_rows=cross_template_rows,
            residue_role_rows=residue_role_rows,
            design_mask_rows=design_mask_rows,
            template_specific_protected_positions=template_specific_protected_positions,
            alignment_summary=alignment_summary,
        ),
        sequence_record_count=len(sequence_rows),
    )


def write_template_harmonization_outputs(
    *,
    artifacts: TemplateHarmonizationArtifacts,
    report_path: Path,
    chain_summary_path: Path,
    site_summary_path: Path,
    cross_template_alignment_path: Path,
    residue_role_map_path: Path,
    design_mask_candidates_path: Path,
    design_masks_path: Path,
    figure_path: Path,
) -> None:
    write_csv_rows(chain_summary_path, artifacts.chain_rows)
    write_csv_rows(site_summary_path, artifacts.site_rows)
    write_csv_rows(cross_template_alignment_path, artifacts.cross_template_rows)
    write_csv_rows(residue_role_map_path, artifacts.residue_role_rows)
    write_csv_rows(design_mask_candidates_path, artifacts.design_mask_rows)
    write_yaml(design_masks_path, artifacts.design_masks_config)
    render_template_harmonization_overview(
        cross_template_rows=artifacts.cross_template_rows,
        residue_role_rows=artifacts.residue_role_rows,
        design_mask_rows=artifacts.design_mask_rows,
        output_path=figure_path,
    )
    atomic_write_text(report_path, artifacts.report_markdown)
