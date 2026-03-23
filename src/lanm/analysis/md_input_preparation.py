"""Phase 6A1 OpenMM-ready MD input preparation for the selected validation panel."""

from __future__ import annotations

import csv
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import yaml

from lanm.analysis.md_panel_selection import MDPanelRow
from lanm.configuration import load_project_config
from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_text, write_csv_rows, write_yaml
from lanm.models import AtomRecord
from lanm.paths import REPO_ROOT
from lanm.structure.atoms import read_atom_records
from lanm.structure.geometry import WATER_RESIDUES, collect_atoms_within_cutoff
from lanm.structure.pdb import read_pdb_atom_records, render_selected_structure_pdb, select_preferred_atom_conformers
from lanm.structure.templates import parse_cif_atom_records

TARGET_METALS = ("Dy", "Nd", "Y", "Al", "Fe")
TARGET_METAL_PDB_LABELS = {
    "Dy": "DY",
    "Nd": "ND",
    "Y": "Y",
    "Al": "AL",
    "Fe": "FE",
}
REPLICATE_COUNT = 3
TARGET_TEMPERATURE_K = 298
SUGGESTED_COLLECTIVE_VARIABLES = (
    "coordination_number",
    "mean_metal_oxygen_distance",
)


@dataclass(frozen=True, slots=True)
class DesignBackboneSource:
    backbone_id: str
    structure_id: str
    source_kind: str
    source_path: Path


@dataclass(frozen=True, slots=True)
class MDStructureTemplate:
    panel_member: MDPanelRow
    selected_chain_ids: tuple[str, ...]
    protein_atoms: tuple[AtomRecord, ...]
    metal_site_atoms: tuple[AtomRecord, ...]
    nearby_solvent_atoms: tuple[AtomRecord, ...]
    source_structure: str
    metal_site_template_source: str


@dataclass(frozen=True, slots=True)
class PreparedMDStartingStructure:
    pdb_text: str
    source_structure: str
    metal_site_template_source: str
    metal_site_count: int
    preserved_solvent_residue_count: int


@dataclass(frozen=True, slots=True)
class MDInputManifestRow:
    panel_rank: int
    panel_member_id: str
    panel_member_type: str
    panel_role: str
    candidate_id: str
    backbone_id: str
    topology_class: str
    target_metal: str
    source_structure: str
    metal_site_template_source: str
    starting_structure_path: str
    system_metadata_path: str
    metal_site_count: int
    preserved_solvent_residue_count: int
    replicate_count: int
    target_temperature_K: int
    suggested_collective_variables: str


def _display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def _repo_path(path_text: str) -> Path:
    path = Path(path_text).expanduser()
    if path.is_absolute():
        return path
    return REPO_ROOT / path


def _parse_optional_int(value: str) -> int | None:
    stripped = value.strip()
    return int(stripped) if stripped else None


def _parse_optional_float(value: str) -> float | None:
    stripped = value.strip()
    return float(stripped) if stripped else None


def _load_md_panel_rows(path: Path) -> tuple[MDPanelRow, ...]:
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            MDPanelRow(
                panel_rank=int(str(row["panel_rank"]).strip()),
                panel_member_id=str(row["panel_member_id"]).strip(),
                panel_member_type=str(row["panel_member_type"]).strip(),
                panel_role=str(row["panel_role"]).strip(),
                candidate_id=str(row["candidate_id"]).strip(),
                campaign_id=str(row["campaign_id"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                topology_class=str(row["topology_class"]).strip(),
                preserved_metal_identity=str(row["preserved_metal_identity"]).strip().upper(),
                proteinmpnn_rank=_parse_optional_int(str(row["proteinmpnn_rank"])),
                ligandmpnn_representative_ligand_confidence=_parse_optional_float(
                    str(row["ligandmpnn_representative_ligand_confidence"])
                ),
                ligandmpnn_representative_overall_confidence=_parse_optional_float(
                    str(row["ligandmpnn_representative_overall_confidence"])
                ),
                ligandmpnn_mean_ligand_confidence=_parse_optional_float(str(row["ligandmpnn_mean_ligand_confidence"])),
                rosetta_score_rank_within_topology=_parse_optional_int(str(row["rosetta_score_rank_within_topology"])),
                rosetta_total_score=_parse_optional_float(str(row["rosetta_total_score"])),
                starting_structure_path=str(row["starting_structure_path"]).strip(),
                selection_reason=str(row["selection_reason"]).strip(),
            )
            for row in reader
        ]
    return tuple(sorted(rows, key=lambda row: (row.panel_rank, row.panel_member_id)))


def _load_yaml_document(path: Path) -> dict[str, Any]:
    require_path(path)
    with path.open("r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"Expected mapping payload in {path}")
    return payload


def discover_md_validation_panel(
    *,
    panel_path: Path,
    md_panel_config_path: Path,
) -> tuple[MDPanelRow, ...]:
    """Load and validate the selected MD panel against the Phase 5B YAML snapshot."""
    panel_rows = _load_md_panel_rows(panel_path)
    config_payload = _load_yaml_document(md_panel_config_path)
    config_rows = config_payload.get("panel")
    if not isinstance(config_rows, list):
        raise ValueError(f"Expected 'panel' list in {md_panel_config_path}")
    if len(config_rows) != len(panel_rows):
        raise ValueError(
            f"Panel row count mismatch between {panel_path} ({len(panel_rows)}) "
            f"and {md_panel_config_path} ({len(config_rows)})"
        )
    for csv_row, yaml_row in zip(panel_rows, config_rows):
        if not isinstance(yaml_row, dict):
            raise ValueError(f"Expected mapping entries in {md_panel_config_path} panel list")
        for field_name in (
            "panel_rank",
            "panel_member_id",
            "panel_member_type",
            "panel_role",
            "candidate_id",
            "backbone_id",
            "topology_class",
            "starting_structure_path",
        ):
            if getattr(csv_row, field_name) != yaml_row.get(field_name):
                raise ValueError(
                    f"Panel mismatch for {csv_row.panel_member_id}: "
                    f"CSV {field_name}={getattr(csv_row, field_name)!r} "
                    f"vs YAML {yaml_row.get(field_name)!r}"
                )
    return panel_rows


def load_design_backbone_sources(path: Path) -> dict[str, DesignBackboneSource]:
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return {
            str(row["backbone_id"]).strip(): DesignBackboneSource(
                backbone_id=str(row["backbone_id"]).strip(),
                structure_id=str(row["structure_id"]).strip(),
                source_kind=str(row["source_kind"]).strip(),
                source_path=_repo_path(str(row["source_path"]).strip()),
            )
            for row in reader
        }


def _selected_chain_ids(atoms: tuple[AtomRecord, ...]) -> tuple[str, ...]:
    chain_ids: list[str] = []
    seen: set[str] = set()
    for atom in atoms:
        if atom.record_type != "ATOM":
            continue
        if atom.chain_id in seen:
            continue
        seen.add(atom.chain_id)
        chain_ids.append(atom.chain_id)
    if not chain_ids:
        raise ValueError("Expected at least one ATOM chain in the starting structure")
    return tuple(chain_ids)


def _load_backbone_source_atoms(source: DesignBackboneSource) -> tuple[AtomRecord, ...]:
    if source.source_kind == "csv_atom_table":
        return tuple(read_atom_records(source.source_path, structure_id=source.structure_id))
    if source.source_kind == "cif":
        return tuple(parse_cif_atom_records(source.source_path, structure_id=source.structure_id))
    raise ValueError(f"Unsupported backbone source_kind for {source.backbone_id}: {source.source_kind}")


def _resolve_primary_structure_atoms(panel_member: MDPanelRow) -> tuple[AtomRecord, ...]:
    source_path = _repo_path(panel_member.starting_structure_path)
    require_path(source_path)
    return tuple(
        select_preferred_atom_conformers(
            read_pdb_atom_records(source_path, structure_id=panel_member.panel_member_id.upper())
        )
    )


def _metal_site_atoms(
    atoms: tuple[AtomRecord, ...],
    *,
    selected_chain_ids: tuple[str, ...],
    preserved_metal_identity: str,
) -> tuple[AtomRecord, ...]:
    return tuple(
        atom
        for atom in atoms
        if atom.record_type == "HETATM"
        and atom.chain_id in selected_chain_ids
        and atom.element_upper == preserved_metal_identity
    )


def _collect_nearby_solvent_atoms(
    *,
    metal_atoms: tuple[AtomRecord, ...],
    context_atoms: tuple[AtomRecord, ...],
    nearby_solvent_cutoff_A: float,
) -> tuple[AtomRecord, ...]:
    solvent_residue_keys = {
        item.atom.residue_key
        for metal_atom in metal_atoms
        for item in collect_atoms_within_cutoff(metal_atom, context_atoms, nearby_solvent_cutoff_A)
        if item.atom.record_type == "HETATM"
        and item.atom.residue_name in WATER_RESIDUES
        and not item.atom.is_metal
    }
    return tuple(
        atom
        for atom in context_atoms
        if atom.record_type == "HETATM"
        and atom.residue_key in solvent_residue_keys
    )


def resolve_panel_member_structure_template(
    *,
    panel_member: MDPanelRow,
    design_backbone_sources: dict[str, DesignBackboneSource],
    nearby_solvent_cutoff_A: float,
) -> MDStructureTemplate:
    """Resolve the structure template used to generate all target-metal systems for one panel member."""
    primary_atoms = _resolve_primary_structure_atoms(panel_member)
    selected_chain_ids = _selected_chain_ids(primary_atoms)
    protein_atoms = tuple(atom for atom in primary_atoms if atom.record_type == "ATOM")
    metal_site_atoms = _metal_site_atoms(
        primary_atoms,
        selected_chain_ids=selected_chain_ids,
        preserved_metal_identity=panel_member.preserved_metal_identity,
    )
    metal_site_template_source = _display_path(_repo_path(panel_member.starting_structure_path))
    solvent_context_atoms = primary_atoms
    if not metal_site_atoms:
        backbone_source = design_backbone_sources.get(panel_member.backbone_id)
        if backbone_source is None:
            raise ValueError(f"No design_backbone_manifest entry found for {panel_member.backbone_id}")
        source_atoms = tuple(select_preferred_atom_conformers(_load_backbone_source_atoms(backbone_source)))
        metal_site_atoms = _metal_site_atoms(
            source_atoms,
            selected_chain_ids=selected_chain_ids,
            preserved_metal_identity=panel_member.preserved_metal_identity,
        )
        metal_site_template_source = _display_path(backbone_source.source_path)
        solvent_context_atoms = source_atoms
    if not metal_site_atoms:
        raise ValueError(
            f"No {panel_member.preserved_metal_identity} metal sites found for panel member {panel_member.panel_member_id}"
        )
    nearby_solvent_atoms = _collect_nearby_solvent_atoms(
        metal_atoms=metal_site_atoms,
        context_atoms=solvent_context_atoms,
        nearby_solvent_cutoff_A=nearby_solvent_cutoff_A,
    )
    return MDStructureTemplate(
        panel_member=panel_member,
        selected_chain_ids=selected_chain_ids,
        protein_atoms=protein_atoms,
        metal_site_atoms=metal_site_atoms,
        nearby_solvent_atoms=nearby_solvent_atoms,
        source_structure=_display_path(_repo_path(panel_member.starting_structure_path)),
        metal_site_template_source=metal_site_template_source,
    )


def _relabel_metal_atom(atom: AtomRecord, *, target_metal: str) -> AtomRecord:
    pdb_label = TARGET_METAL_PDB_LABELS[target_metal]
    return replace(
        atom,
        atom_name=pdb_label,
        residue_name=pdb_label,
        element=pdb_label,
        charge="",
    )


def render_target_metal_starting_structure(
    template: MDStructureTemplate,
    *,
    target_metal: str,
) -> PreparedMDStartingStructure:
    """Render a deterministic starting PDB with all bound-metal sites relabeled to one target metal."""
    if target_metal not in TARGET_METAL_PDB_LABELS:
        raise ValueError(f"Unsupported target metal: {target_metal}")
    relabeled_metal_atoms = tuple(_relabel_metal_atom(atom, target_metal=target_metal) for atom in template.metal_site_atoms)
    combined_atoms = tuple(
        select_preferred_atom_conformers(
            (
                *template.protein_atoms,
                *relabeled_metal_atoms,
                *template.nearby_solvent_atoms,
            )
        )
    )
    included_het_residue_keys = frozenset(
        atom.residue_key
        for atom in (*relabeled_metal_atoms, *template.nearby_solvent_atoms)
    )
    return PreparedMDStartingStructure(
        pdb_text=render_selected_structure_pdb(
            combined_atoms,
            selected_chain_ids=template.selected_chain_ids,
            included_het_residue_keys=included_het_residue_keys,
        ),
        source_structure=template.source_structure,
        metal_site_template_source=template.metal_site_template_source,
        metal_site_count=len(relabeled_metal_atoms),
        preserved_solvent_residue_count=len({atom.residue_key for atom in template.nearby_solvent_atoms}),
    )


def _build_system_metadata(
    *,
    panel_member: MDPanelRow,
    target_metal: str,
    source_structure: str,
    metal_site_template_source: str,
    starting_structure_path: Path,
    prepared_structure: PreparedMDStartingStructure,
) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "panel_member_id": panel_member.panel_member_id,
        "candidate_id": panel_member.candidate_id,
        "topology_class": panel_member.topology_class,
        "target_metal": target_metal,
        "source_structure": source_structure,
        "starting_structure_path": _display_path(starting_structure_path),
        "replicate_count": REPLICATE_COUNT,
        "target_temperature_K": TARGET_TEMPERATURE_K,
        "suggested_collective_variables": list(SUGGESTED_COLLECTIVE_VARIABLES),
        "panel_member_type": panel_member.panel_member_type,
        "panel_role": panel_member.panel_role,
        "backbone_id": panel_member.backbone_id,
        "metal_site_count": prepared_structure.metal_site_count,
        "preserved_solvent_residue_count": prepared_structure.preserved_solvent_residue_count,
    }
    if metal_site_template_source != source_structure:
        payload["metal_site_template_source"] = metal_site_template_source
    return payload


def _build_md_protocol_config(
    *,
    panel_rows: tuple[MDPanelRow, ...],
    manifest_rows: tuple[MDInputManifestRow, ...],
    panel_path: Path,
    md_panel_config_path: Path,
) -> dict[str, Any]:
    systems_by_panel_member: dict[str, list[MDInputManifestRow]] = {}
    for row in manifest_rows:
        systems_by_panel_member.setdefault(row.panel_member_id, []).append(row)
    return {
        "version": 1,
        "phase": "6A1",
        "description": "Prepare OpenMM-ready multi-metal starting structures and metadata only; no systems or simulations are built in this phase.",
        "source_artifacts": {
            "md_validation_panel": _display_path(panel_path),
            "md_panel_config": _display_path(md_panel_config_path),
        },
        "target_metals": list(TARGET_METALS),
        "replicate_count": REPLICATE_COUNT,
        "target_temperature_K": TARGET_TEMPERATURE_K,
        "suggested_collective_variables": list(SUGGESTED_COLLECTIVE_VARIABLES),
        "panel_members": [
            {
                "panel_rank": panel_row.panel_rank,
                "panel_member_id": panel_row.panel_member_id,
                "candidate_id": panel_row.candidate_id,
                "panel_member_type": panel_row.panel_member_type,
                "panel_role": panel_row.panel_role,
                "topology_class": panel_row.topology_class,
                "systems": [
                    {
                        "target_metal": system_row.target_metal,
                        "starting_structure_path": system_row.starting_structure_path,
                        "system_metadata_path": system_row.system_metadata_path,
                    }
                    for system_row in systems_by_panel_member.get(panel_row.panel_member_id, [])
                ],
            }
            for panel_row in panel_rows
        ],
    }


def render_md_input_preparation_markdown(
    *,
    panel_rows: tuple[MDPanelRow, ...],
    manifest_rows: tuple[MDInputManifestRow, ...],
) -> str:
    """Render the Phase 6A1 preparation report."""
    lines = [
        "# MD Input Preparation",
        "",
        "Phase 6A1 prepares deterministic OpenMM-ready starting structures and metadata for the selected MD validation panel.",
        "",
        "## Scope",
        "",
        f"- panel members prepared: `{len(panel_rows)}`",
        f"- target metals prepared per panel member: `{', '.join(TARGET_METALS)}`",
        f"- total systems written: `{len(manifest_rows)}`",
        f"- replicate_count recorded for later production MD: `{REPLICATE_COUNT}`",
        f"- target temperature recorded for later production MD: `{TARGET_TEMPERATURE_K} K`",
        f"- suggested collective variables recorded for later metadynamics: `{', '.join(SUGGESTED_COLLECTIVE_VARIABLES)}`",
        "- excluded in this phase: `OpenMM System building`, `MD`, `metadynamics`, `QM`, and quantum steps",
        "",
        "## Prepared Systems",
        "",
        "| panel_member_id | target_metal | topology_class | metal_site_count | preserved_solvent_residue_count | source_structure |",
        "| --- | --- | --- | ---: | ---: | --- |",
    ]
    for row in manifest_rows:
        lines.append(
            "| "
            f"{row.panel_member_id} | "
            f"{row.target_metal} | "
            f"{row.topology_class} | "
            f"{row.metal_site_count} | "
            f"{row.preserved_solvent_residue_count} | "
            f"{row.source_structure} |"
        )
    wild_type_rows = [row for row in manifest_rows if row.panel_member_type == "wild_type_reference"]
    if wild_type_rows:
        lines.extend(
            [
                "",
                "## Wild-Type Reference Handling",
                "",
                "- Wild-type reference backbones remain the panel entry points, but missing metal-site context is restored deterministically from the original backbone source structure when needed.",
            ]
        )
        fallback_rows = [
            row
            for row in wild_type_rows
            if row.metal_site_template_source != row.source_structure
        ]
        reported_pairs: set[tuple[str, str, str]] = set()
        for row in fallback_rows:
            pair = (
                row.panel_member_id,
                row.source_structure,
                row.metal_site_template_source,
            )
            if pair in reported_pairs:
                continue
            reported_pairs.add(pair)
            lines.append(
                f"- `{row.panel_member_id}` uses `{row.metal_site_template_source}` as the metal/water template while preserving the protein backbone from `{row.source_structure}`."
            )
    lines.append("")
    return "\n".join(lines)


def prepare_md_inputs(
    *,
    panel_path: Path,
    md_panel_config_path: Path,
    design_backbone_manifest_path: Path,
    manifest_path: Path,
    report_path: Path,
    protocol_config_path: Path,
    md_inputs_root: Path,
) -> tuple[MDInputManifestRow, ...]:
    """Prepare deterministic multi-metal MD starting structures and metadata."""
    project_config = load_project_config()
    configured_metals = (project_config.target_metal, *project_config.competitors)
    if configured_metals != TARGET_METALS:
        raise ValueError(
            f"Expected project metals {TARGET_METALS}, found {configured_metals}"
        )
    if project_config.temperature_K != TARGET_TEMPERATURE_K:
        raise ValueError(
            f"Expected project temperature {TARGET_TEMPERATURE_K} K, found {project_config.temperature_K} K"
        )
    panel_rows = discover_md_validation_panel(
        panel_path=panel_path,
        md_panel_config_path=md_panel_config_path,
    )
    design_backbone_sources = load_design_backbone_sources(design_backbone_manifest_path)
    manifest_rows: list[MDInputManifestRow] = []
    for panel_member in panel_rows:
        template = resolve_panel_member_structure_template(
            panel_member=panel_member,
            design_backbone_sources=design_backbone_sources,
            nearby_solvent_cutoff_A=project_config.second_sphere_cutoff_A,
        )
        for target_metal in TARGET_METALS:
            prepared_structure = render_target_metal_starting_structure(
                template,
                target_metal=target_metal,
            )
            output_dir = md_inputs_root / panel_member.panel_member_id / target_metal
            starting_structure_path = output_dir / "starting_structure.pdb"
            metadata_path = output_dir / "system_metadata.yaml"
            atomic_write_text(starting_structure_path, prepared_structure.pdb_text)
            write_yaml(
                metadata_path,
                _build_system_metadata(
                    panel_member=panel_member,
                    target_metal=target_metal,
                    source_structure=prepared_structure.source_structure,
                    metal_site_template_source=prepared_structure.metal_site_template_source,
                    starting_structure_path=starting_structure_path,
                    prepared_structure=prepared_structure,
                ),
            )
            manifest_rows.append(
                MDInputManifestRow(
                    panel_rank=panel_member.panel_rank,
                    panel_member_id=panel_member.panel_member_id,
                    panel_member_type=panel_member.panel_member_type,
                    panel_role=panel_member.panel_role,
                    candidate_id=panel_member.candidate_id,
                    backbone_id=panel_member.backbone_id,
                    topology_class=panel_member.topology_class,
                    target_metal=target_metal,
                    source_structure=prepared_structure.source_structure,
                    metal_site_template_source=prepared_structure.metal_site_template_source,
                    starting_structure_path=_display_path(starting_structure_path),
                    system_metadata_path=_display_path(metadata_path),
                    metal_site_count=prepared_structure.metal_site_count,
                    preserved_solvent_residue_count=prepared_structure.preserved_solvent_residue_count,
                    replicate_count=REPLICATE_COUNT,
                    target_temperature_K=TARGET_TEMPERATURE_K,
                    suggested_collective_variables=",".join(SUGGESTED_COLLECTIVE_VARIABLES),
                )
            )
    manifest_tuple = tuple(manifest_rows)
    write_csv_rows(manifest_path, manifest_tuple)
    atomic_write_text(
        report_path,
        render_md_input_preparation_markdown(
            panel_rows=panel_rows,
            manifest_rows=manifest_tuple,
        ),
    )
    write_yaml(
        protocol_config_path,
        _build_md_protocol_config(
            panel_rows=panel_rows,
            manifest_rows=manifest_tuple,
            panel_path=panel_path,
            md_panel_config_path=md_panel_config_path,
        ),
    )
    return manifest_tuple
