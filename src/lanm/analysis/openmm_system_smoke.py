"""Phase 6A3 representative OpenMM system-build smoke testing."""

from __future__ import annotations

import csv
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import yaml

from lanm.analysis.md_system_build import DEFAULT_PROTEIN_FORCEFIELD, DEFAULT_WATER_MODEL, load_md_protocol_build_config
from lanm.analysis.metal_parameter_mapping import load_md_system_build_manifest
from lanm.analysis.metal_parameter_values import GENERIC_12_6_4_HIGHLY_CHARGED
from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_text, write_csv_rows
from lanm.md.openmm_system_smoke import OpenMMSerializedSystem, build_openmm_serialized_system
from lanm.models import AtomRecord
from lanm.paths import REPO_ROOT, RESULTS_MD_INPUTS_DIR
from lanm.structure.atoms import read_atom_records
from lanm.structure.geometry import WATER_RESIDUES, collect_atoms_within_cutoff
from lanm.structure.pdb import read_pdb_atom_records, render_selected_structure_pdb, select_preferred_atom_conformers
from lanm.structure.templates import parse_cif_atom_records

EXPECTED_REPRESENTATIVE_SYSTEMS = (
    ("am1_monomer", "am1_mex_ss_only_u02", "Dy"),
    ("hans_monomer", "hans_pocket_ss_only_u02", "Dy"),
    ("hans_interface_multichain", "hans_interface_ss_plus_if_u04", "Dy"),
)
EXPECTED_FORCEFIELD_FILES = (DEFAULT_PROTEIN_FORCEFIELD, DEFAULT_WATER_MODEL)
FALLBACK_TEMPLATE_WATER_CUTOFF_A = 6.0


@dataclass(frozen=True, slots=True)
class MDValidationPanelRow:
    panel_rank: int
    panel_member_id: str
    panel_member_type: str
    topology_class: str


@dataclass(frozen=True, slots=True)
class MetalBuildReadinessRow:
    panel_rank: int
    panel_member_id: str
    target_metal: str
    topology_class: str
    build_request_path: str
    baseline_parameter_family: str
    ready_for_openmm_system_build_baseline: bool
    ready_for_openmm_system_build_tuned: bool


@dataclass(frozen=True, slots=True)
class MetalParameterValuesAudit:
    supported_water_models: tuple[str, ...]
    covered_metals: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class RepresentativeOpenMMSystem:
    panel_rank: int
    panel_member_id: str
    target_metal: str
    topology_class: str
    build_request_path: str


@dataclass(frozen=True, slots=True)
class OpenMMSystemBuildRequest:
    panel_member_id: str
    target_metal: str
    topology_class: str
    source_starting_structure_pdb: str
    chosen_protein_forcefield: str
    chosen_water_model: str
    restore_template_metal_proximal_waters_before_solvation: bool
    template_metal_proximal_waters_available: bool
    comparison_panel_member_id: str
    comparison_source: str
    comparison_preserved_solvent_residue_count: int


@dataclass(frozen=True, slots=True)
class PreparedStructureInput:
    pdb_text: str
    restored_template_waters: int
    template_water_reference_path: str


@dataclass(frozen=True, slots=True)
class OpenMMSystemSmokeSummaryRow:
    panel_member_id: str
    target_metal: str
    topology_class: str
    success_status: str
    failure_reason: str
    forcefield_files: str
    custom_metal_family: str
    restored_template_waters: int
    atom_count: int
    residue_count: int
    system_xml_path: str


@dataclass(frozen=True, slots=True)
class ProtocolRuntimeConfig:
    target_temperature_K: int


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


def _parse_required_mapping(path: Path) -> dict[str, Any]:
    require_path(path)
    with path.open("r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"Expected mapping payload in {path}")
    return payload


def _parse_csv_bool(row: dict[str, str], field_name: str) -> bool:
    value = str(row[field_name]).strip().lower()
    if value == "true":
        return True
    if value == "false":
        return False
    raise ValueError(f"Expected boolean text for {field_name}, found {row[field_name]!r}")


def load_md_validation_panel_rows(path: Path) -> tuple[MDValidationPanelRow, ...]:
    """Load the selected MD validation panel rows needed for representative selection."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = tuple(
            sorted(
                (
                    MDValidationPanelRow(
                        panel_rank=int(str(row["panel_rank"]).strip()),
                        panel_member_id=str(row["panel_member_id"]).strip(),
                        panel_member_type=str(row["panel_member_type"]).strip(),
                        topology_class=str(row["topology_class"]).strip(),
                    )
                    for row in reader
                ),
                key=lambda row: (row.panel_rank, row.panel_member_id),
            )
        )
    if not rows:
        raise ValueError(f"No rows found in {path}")
    return rows


def load_metal_build_readiness_rows(path: Path) -> tuple[MetalBuildReadinessRow, ...]:
    """Load the baseline/tuned readiness audit needed for representative filtering."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = tuple(
            MetalBuildReadinessRow(
                panel_rank=int(str(row["panel_rank"]).strip()),
                panel_member_id=str(row["panel_member_id"]).strip(),
                target_metal=str(row["target_metal"]).strip(),
                topology_class=str(row["topology_class"]).strip(),
                build_request_path=str(row["build_request_path"]).strip(),
                baseline_parameter_family=str(row["baseline_parameter_family"]).strip(),
                ready_for_openmm_system_build_baseline=_parse_csv_bool(
                    row,
                    "ready_for_openmm_system_build_baseline",
                ),
                ready_for_openmm_system_build_tuned=_parse_csv_bool(
                    row,
                    "ready_for_openmm_system_build_tuned",
                ),
            )
            for row in reader
        )
    if not rows:
        raise ValueError(f"No rows found in {path}")
    return rows


def load_metal_parameter_values_audit(path: Path) -> MetalParameterValuesAudit:
    """Validate that the generic baseline family is present for the Dy smoke systems."""
    payload = _parse_required_mapping(path)
    families = payload.get("parameter_value_families")
    numeric_values = payload.get("numeric_parameter_values")
    if not isinstance(families, list):
        raise ValueError(f"Expected parameter_value_families list in {path}")
    if not isinstance(numeric_values, list):
        raise ValueError(f"Expected numeric_parameter_values list in {path}")
    generic_family = next(
        (
            family
            for family in families
            if isinstance(family, dict)
            and str(family.get("source_family_label", "")).strip() == GENERIC_12_6_4_HIGHLY_CHARGED
        ),
        None,
    )
    if generic_family is None:
        raise ValueError(
            f"Missing {GENERIC_12_6_4_HIGHLY_CHARGED!r} family in {path}"
        )
    if not bool(generic_family.get("numeric_values_present", False)):
        raise ValueError(f"{GENERIC_12_6_4_HIGHLY_CHARGED!r} must be numerically available in {path}")
    supported_water_models = tuple(str(item).strip() for item in generic_family.get("water_model_compatibility", []))
    if DEFAULT_WATER_MODEL not in supported_water_models:
        raise ValueError(
            f"{GENERIC_12_6_4_HIGHLY_CHARGED!r} must support {DEFAULT_WATER_MODEL!r} in {path}"
        )
    covered_metals = tuple(str(item).strip() for item in generic_family.get("covered_metals", []))
    if "Dy" not in covered_metals:
        raise ValueError(f"{GENERIC_12_6_4_HIGHLY_CHARGED!r} must cover Dy in {path}")
    if not any(
        isinstance(entry, dict)
        and str(entry.get("metal_identity", "")).strip() == "Dy"
        and str(entry.get("source_family_label", "")).strip() == GENERIC_12_6_4_HIGHLY_CHARGED
        and str(entry.get("water_model", "")).strip() == DEFAULT_WATER_MODEL
        for entry in numeric_values
    ):
        raise ValueError(
            f"Missing Dy/{GENERIC_12_6_4_HIGHLY_CHARGED}/{DEFAULT_WATER_MODEL} numeric entry in {path}"
        )
    return MetalParameterValuesAudit(
        supported_water_models=supported_water_models,
        covered_metals=covered_metals,
    )


def load_openmm_system_smoke_protocol_runtime(path: Path) -> ProtocolRuntimeConfig:
    """Load the protocol temperature used for the Phase 6A3 serialization-only integrator."""
    payload = _parse_required_mapping(path)
    target_temperature = int(payload.get("target_temperature_K", 0) or 0)
    if target_temperature <= 0:
        raise ValueError(f"Expected positive target_temperature_K in {path}")
    return ProtocolRuntimeConfig(target_temperature_K=target_temperature)


def select_representative_baseline_systems(
    *,
    panel_rows: tuple[MDValidationPanelRow, ...],
    build_manifest_path: Path,
    readiness_rows: tuple[MetalBuildReadinessRow, ...],
    md_protocol_path: Path,
    metal_parameter_values_path: Path,
) -> tuple[RepresentativeOpenMMSystem, ...]:
    """Select one Dy designed candidate per topology class for the baseline smoke test."""
    build_config = load_md_protocol_build_config(md_protocol_path)
    if (
        build_config.protein_forcefield,
        build_config.water_model,
    ) != EXPECTED_FORCEFIELD_FILES:
        raise ValueError(
            "Phase 6A3 requires amber19-all.xml + amber19/opc3.xml, "
            f"found {(build_config.protein_forcefield, build_config.water_model)!r}"
        )
    load_metal_parameter_values_audit(metal_parameter_values_path)
    manifest_rows = load_md_system_build_manifest(build_manifest_path)
    manifest_by_key = {(row.panel_member_id, row.target_metal): row for row in manifest_rows}
    readiness_by_key = {(row.panel_member_id, row.target_metal): row for row in readiness_rows}
    selected_rows: list[RepresentativeOpenMMSystem] = []
    seen_topologies: set[str] = set()
    for panel_row in panel_rows:
        if panel_row.panel_member_type != "designed_candidate":
            continue
        if panel_row.topology_class in seen_topologies:
            continue
        key = (panel_row.panel_member_id, "Dy")
        readiness_row = readiness_by_key.get(key)
        if readiness_row is None:
            raise ValueError(f"Missing metal_build_readiness row for {panel_row.panel_member_id} Dy")
        if not readiness_row.ready_for_openmm_system_build_baseline:
            continue
        if readiness_row.baseline_parameter_family != GENERIC_12_6_4_HIGHLY_CHARGED:
            raise ValueError(
                f"{panel_row.panel_member_id} Dy must use {GENERIC_12_6_4_HIGHLY_CHARGED!r}, "
                f"found {readiness_row.baseline_parameter_family!r}"
            )
        manifest_row = manifest_by_key.get(key)
        if manifest_row is None:
            raise ValueError(f"Missing md_system_build_manifest row for {panel_row.panel_member_id} Dy")
        if manifest_row.topology_class != panel_row.topology_class:
            raise ValueError(
                f"Topology mismatch for {panel_row.panel_member_id}: "
                f"{manifest_row.topology_class!r} vs {panel_row.topology_class!r}"
            )
        if manifest_row.chosen_protein_forcefield != build_config.protein_forcefield:
            raise ValueError(
                f"{panel_row.panel_member_id} Dy uses unexpected protein force field "
                f"{manifest_row.chosen_protein_forcefield!r}"
            )
        if manifest_row.chosen_water_model != build_config.water_model:
            raise ValueError(
                f"{panel_row.panel_member_id} Dy uses unexpected water model "
                f"{manifest_row.chosen_water_model!r}"
            )
        if not manifest_row.generic_protein_water_preparation_ready:
            raise ValueError(f"{panel_row.panel_member_id} Dy is not generic protein/water ready")
        selected_rows.append(
            RepresentativeOpenMMSystem(
                panel_rank=panel_row.panel_rank,
                panel_member_id=panel_row.panel_member_id,
                target_metal="Dy",
                topology_class=panel_row.topology_class,
                build_request_path=manifest_row.build_request_path,
            )
        )
        seen_topologies.add(panel_row.topology_class)
    selected_rows = sorted(selected_rows, key=lambda row: (row.panel_rank, row.panel_member_id))
    observed = tuple((row.topology_class, row.panel_member_id, row.target_metal) for row in selected_rows)
    if observed != EXPECTED_REPRESENTATIVE_SYSTEMS:
        raise ValueError(
            "Phase 6A3 representative selection drifted; expected "
            f"{EXPECTED_REPRESENTATIVE_SYSTEMS!r}, found {observed!r}"
        )
    return tuple(selected_rows)


def parse_openmm_system_build_request(path: Path) -> OpenMMSystemBuildRequest:
    """Load and validate the per-system build request YAML used for Phase 6A3."""
    payload = _parse_required_mapping(path)
    reference = payload.get("template_water_restoration_reference")
    if not isinstance(reference, dict):
        raise ValueError(f"Expected template_water_restoration_reference mapping in {path}")
    request = OpenMMSystemBuildRequest(
        panel_member_id=str(payload["panel_member_id"]).strip(),
        target_metal=str(payload["target_metal"]).strip(),
        topology_class=str(payload["topology_class"]).strip(),
        source_starting_structure_pdb=str(payload["source_starting_structure_pdb"]).strip(),
        chosen_protein_forcefield=str(payload["chosen_protein_forcefield"]).strip(),
        chosen_water_model=str(payload["chosen_water_model"]).strip(),
        restore_template_metal_proximal_waters_before_solvation=bool(
            payload["restore_template_metal_proximal_waters_before_solvation"]
        ),
        template_metal_proximal_waters_available=bool(payload["template_metal_proximal_waters_available"]),
        comparison_panel_member_id=str(reference.get("comparison_panel_member_id", "") or "").strip(),
        comparison_source=str(reference.get("comparison_source", "") or "").strip(),
        comparison_preserved_solvent_residue_count=int(
            reference.get("comparison_preserved_solvent_residue_count", 0) or 0
        ),
    )
    if (
        request.chosen_protein_forcefield,
        request.chosen_water_model,
    ) != EXPECTED_FORCEFIELD_FILES:
        raise ValueError(
            f"Phase 6A3 build request must use {EXPECTED_FORCEFIELD_FILES!r}, "
            f"found {(request.chosen_protein_forcefield, request.chosen_water_model)!r} in {path}"
        )
    if request.restore_template_metal_proximal_waters_before_solvation:
        if not request.template_metal_proximal_waters_available:
            raise ValueError(
                f"{path} requests template-water restoration without marking template waters available"
            )
        if not request.comparison_panel_member_id and not request.comparison_source:
            raise ValueError(f"{path} requests template-water restoration without a usable reference")
    require_path(_repo_path(request.source_starting_structure_pdb))
    return request


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
        raise ValueError("Expected at least one polymer chain in the prepared structure")
    return tuple(chain_ids)


def _load_supported_atom_records(path: Path) -> tuple[AtomRecord, ...]:
    suffix = path.suffix.lower()
    if suffix == ".pdb":
        atoms = read_pdb_atom_records(path, structure_id=path.stem.upper())
    elif suffix == ".csv":
        atoms = read_atom_records(path, structure_id=path.stem.upper())
    elif suffix == ".cif":
        atoms = parse_cif_atom_records(path, structure_id=path.stem.upper())
    else:
        raise ValueError(f"Unsupported template-water reference format: {path}")
    return tuple(select_preferred_atom_conformers(atoms))


def _collect_template_water_atoms_from_reference(
    *,
    reference_path: Path,
    selected_chain_ids: tuple[str, ...],
) -> tuple[AtomRecord, ...]:
    reference_atoms = _load_supported_atom_records(reference_path)
    if reference_path.suffix.lower() == ".pdb":
        return tuple(
            atom
            for atom in reference_atoms
            if atom.record_type == "HETATM"
            and atom.residue_name in WATER_RESIDUES
            and atom.chain_id in selected_chain_ids
        )
    metal_atoms = tuple(
        atom
        for atom in reference_atoms
        if atom.record_type == "HETATM"
        and atom.is_metal
        and atom.chain_id in selected_chain_ids
    )
    if not metal_atoms:
        raise ValueError(f"No template metal atoms found in {reference_path}")
    solvent_residue_keys = {
        item.atom.residue_key
        for metal_atom in metal_atoms
        for item in collect_atoms_within_cutoff(metal_atom, reference_atoms, FALLBACK_TEMPLATE_WATER_CUTOFF_A)
        if item.atom.record_type == "HETATM"
        and item.atom.residue_name in WATER_RESIDUES
        and not item.atom.is_metal
        and item.atom.chain_id in selected_chain_ids
    }
    return tuple(
        atom
        for atom in reference_atoms
        if atom.record_type == "HETATM"
        and atom.residue_key in solvent_residue_keys
    )


def _resolve_template_water_reference_path(request: OpenMMSystemBuildRequest) -> Path:
    if request.comparison_panel_member_id:
        candidate = RESULTS_MD_INPUTS_DIR / request.comparison_panel_member_id / request.target_metal / "starting_structure.pdb"
        if candidate.exists():
            return candidate
    if request.comparison_source:
        path = _repo_path(request.comparison_source)
        require_path(path)
        return path
    raise ValueError(
        f"No template-water reference path could be resolved for {request.panel_member_id} {request.target_metal}"
    )


def prepare_structure_input(request: OpenMMSystemBuildRequest) -> PreparedStructureInput:
    """Assemble the deterministic pre-hydrogenation PDB text for one smoke-test system."""
    source_path = _repo_path(request.source_starting_structure_pdb)
    source_atoms = tuple(
        select_preferred_atom_conformers(
            read_pdb_atom_records(source_path, structure_id=request.panel_member_id.upper())
        )
    )
    selected_chain_ids = _selected_chain_ids(source_atoms)
    restored_template_waters = 0
    template_water_reference_path = ""
    combined_atoms = source_atoms
    if request.restore_template_metal_proximal_waters_before_solvation:
        reference_path = _resolve_template_water_reference_path(request)
        template_water_reference_path = _display_path(reference_path)
        template_water_atoms = _collect_template_water_atoms_from_reference(
            reference_path=reference_path,
            selected_chain_ids=selected_chain_ids,
        )
        if not template_water_atoms:
            raise ValueError(
                f"Template-water restoration was requested for {request.panel_member_id} "
                f"but no template waters were resolved from {reference_path}"
            )
        combined_atoms = tuple(select_preferred_atom_conformers((*source_atoms, *template_water_atoms)))
        restored_template_waters = len({atom.residue_key for atom in template_water_atoms})
    included_het_residue_keys = frozenset(
        atom.residue_key
        for atom in combined_atoms
        if atom.record_type == "HETATM"
    )
    return PreparedStructureInput(
        pdb_text=render_selected_structure_pdb(
            combined_atoms,
            selected_chain_ids=selected_chain_ids,
            included_het_residue_keys=included_het_residue_keys,
        ),
        restored_template_waters=restored_template_waters,
        template_water_reference_path=template_water_reference_path,
    )


def render_openmm_system_smoke_report(results: tuple[OpenMMSystemSmokeSummaryRow, ...]) -> str:
    """Render the concise Phase 6A3 Markdown report."""
    success_count = sum(row.success_status == "success" for row in results)
    failure_count = len(results) - success_count
    lines = [
        "# OpenMM System Smoke",
        "",
        "Phase 6A3 validates one Dy baseline-ready representative per topology class using the generic OPC3 12-6-4 family path.",
        "",
        "## Build Summary",
        "",
        f"- representative systems attempted: `{len(results)}`",
        f"- successful OpenMM `System` serializations: `{success_count}`",
        f"- failures recorded: `{failure_count}`",
        f"- force-field files: `{EXPECTED_FORCEFIELD_FILES[0]}`, `{EXPECTED_FORCEFIELD_FILES[1]}`",
        f"- custom metal family audited for this phase: `{GENERIC_12_6_4_HIGHLY_CHARGED}`",
        "",
        "| panel_member_id | target_metal | topology_class | success_status | restored_template_waters | atom_count | residue_count | failure_reason |",
        "| --- | --- | --- | --- | ---: | ---: | ---: | --- |",
    ]
    for row in results:
        failure_reason = row.failure_reason or "-"
        lines.append(
            "| "
            f"{row.panel_member_id} | "
            f"{row.target_metal} | "
            f"{row.topology_class} | "
            f"{row.success_status} | "
            f"{row.restored_template_waters} | "
            f"{row.atom_count} | "
            f"{row.residue_count} | "
            f"{failure_reason} |"
        )
    lines.append("")
    return "\n".join(lines)


def _write_system_build_log(
    *,
    log_path: Path,
    request: OpenMMSystemBuildRequest,
    restored_template_waters: int,
    template_water_reference_path: str,
    atom_count: int,
    residue_count: int,
    status: str,
    failure_reason: str,
    system_xml_path: str,
    traceback_text: str = "",
) -> None:
    lines = [
        "Phase 6A3 OpenMM system-build smoke",
        f"panel_member_id: {request.panel_member_id}",
        f"target_metal: {request.target_metal}",
        f"topology_class: {request.topology_class}",
        f"source_starting_structure_pdb: {request.source_starting_structure_pdb}",
        f"forcefield_files: {EXPECTED_FORCEFIELD_FILES[0]}, {EXPECTED_FORCEFIELD_FILES[1]}",
        f"custom_metal_family: {GENERIC_12_6_4_HIGHLY_CHARGED}",
        (
            "restore_template_metal_proximal_waters_before_solvation: "
            f"{request.restore_template_metal_proximal_waters_before_solvation}"
        ),
        f"template_water_reference_path: {template_water_reference_path or '-'}",
        f"restored_template_waters: {restored_template_waters}",
        f"atom_count: {atom_count}",
        f"residue_count: {residue_count}",
        f"system_xml_path: {system_xml_path or '-'}",
        f"status: {status}",
    ]
    if failure_reason:
        lines.append(f"failure_reason: {failure_reason}")
    if traceback_text:
        lines.extend(("", "traceback:", traceback_text.rstrip()))
    atomic_write_text(log_path, "\n".join(lines) + "\n")


def run_openmm_system_smoke(
    *,
    md_validation_panel_path: Path,
    md_system_build_manifest_path: Path,
    metal_build_readiness_path: Path,
    md_protocol_path: Path,
    metal_parameter_values_path: Path,
    output_root: Path,
    summary_path: Path,
    report_path: Path,
    system_builder: Callable[..., OpenMMSerializedSystem] = build_openmm_serialized_system,
) -> tuple[OpenMMSystemSmokeSummaryRow, ...]:
    """Run the deterministic Phase 6A3 representative OpenMM system-build smoke test."""
    panel_rows = load_md_validation_panel_rows(md_validation_panel_path)
    readiness_rows = load_metal_build_readiness_rows(metal_build_readiness_path)
    representative_systems = select_representative_baseline_systems(
        panel_rows=panel_rows,
        build_manifest_path=md_system_build_manifest_path,
        readiness_rows=readiness_rows,
        md_protocol_path=md_protocol_path,
        metal_parameter_values_path=metal_parameter_values_path,
    )
    runtime_config = load_openmm_system_smoke_protocol_runtime(md_protocol_path)
    results: list[OpenMMSystemSmokeSummaryRow] = []
    for representative in representative_systems:
        output_dir = output_root / representative.panel_member_id / representative.target_metal
        prepared_structure_path = output_dir / "prepared_structure.pdb"
        system_xml_path = output_dir / "system.xml"
        integrator_xml_path = output_dir / "integrator.xml"
        log_path = output_dir / "system_build_log.txt"
        for artifact_path in (prepared_structure_path, system_xml_path, integrator_xml_path):
            artifact_path.unlink(missing_ok=True)
        build_request: OpenMMSystemBuildRequest | None = None
        restored_template_waters = 0
        template_water_reference_path = ""
        atom_count = 0
        residue_count = 0
        system_xml_display_path = ""
        try:
            build_request = parse_openmm_system_build_request(_repo_path(representative.build_request_path))
            prepared_input = prepare_structure_input(build_request)
            restored_template_waters = prepared_input.restored_template_waters
            template_water_reference_path = prepared_input.template_water_reference_path
            built_system = system_builder(
                prepared_structure_pdb_text=prepared_input.pdb_text,
                forcefield_files=EXPECTED_FORCEFIELD_FILES,
                target_temperature_K=runtime_config.target_temperature_K,
            )
            atom_count = built_system.atom_count
            residue_count = built_system.residue_count
            atomic_write_text(prepared_structure_path, built_system.prepared_structure_pdb_text)
            atomic_write_text(system_xml_path, built_system.system_xml)
            atomic_write_text(integrator_xml_path, built_system.integrator_xml)
            system_xml_display_path = _display_path(system_xml_path)
            _write_system_build_log(
                log_path=log_path,
                request=build_request,
                restored_template_waters=restored_template_waters,
                template_water_reference_path=template_water_reference_path,
                atom_count=atom_count,
                residue_count=residue_count,
                status="success",
                failure_reason="",
                system_xml_path=system_xml_display_path,
            )
            results.append(
                OpenMMSystemSmokeSummaryRow(
                    panel_member_id=representative.panel_member_id,
                    target_metal=representative.target_metal,
                    topology_class=representative.topology_class,
                    success_status="success",
                    failure_reason="",
                    forcefield_files=";".join(EXPECTED_FORCEFIELD_FILES),
                    custom_metal_family=GENERIC_12_6_4_HIGHLY_CHARGED,
                    restored_template_waters=restored_template_waters,
                    atom_count=atom_count,
                    residue_count=residue_count,
                    system_xml_path=system_xml_display_path,
                )
            )
        except Exception as exc:
            failure_reason = f"{type(exc).__name__}: {exc}"
            traceback_text = traceback.format_exc()
            if build_request is None:
                build_request = OpenMMSystemBuildRequest(
                    panel_member_id=representative.panel_member_id,
                    target_metal=representative.target_metal,
                    topology_class=representative.topology_class,
                    source_starting_structure_pdb="",
                    chosen_protein_forcefield=EXPECTED_FORCEFIELD_FILES[0],
                    chosen_water_model=EXPECTED_FORCEFIELD_FILES[1],
                    restore_template_metal_proximal_waters_before_solvation=False,
                    template_metal_proximal_waters_available=False,
                    comparison_panel_member_id="",
                    comparison_source="",
                    comparison_preserved_solvent_residue_count=0,
                )
            _write_system_build_log(
                log_path=log_path,
                request=build_request,
                restored_template_waters=restored_template_waters,
                template_water_reference_path=template_water_reference_path,
                atom_count=atom_count,
                residue_count=residue_count,
                status="failure",
                failure_reason=failure_reason,
                system_xml_path="",
                traceback_text=traceback_text,
            )
            results.append(
                OpenMMSystemSmokeSummaryRow(
                    panel_member_id=representative.panel_member_id,
                    target_metal=representative.target_metal,
                    topology_class=representative.topology_class,
                    success_status="failure",
                    failure_reason=failure_reason,
                    forcefield_files=";".join(EXPECTED_FORCEFIELD_FILES),
                    custom_metal_family=GENERIC_12_6_4_HIGHLY_CHARGED,
                    restored_template_waters=restored_template_waters,
                    atom_count=atom_count,
                    residue_count=residue_count,
                    system_xml_path="",
                )
            )
    result_tuple = tuple(results)
    write_csv_rows(summary_path, result_tuple)
    atomic_write_text(report_path, render_openmm_system_smoke_report(result_tuple))
    if any(row.success_status != "success" for row in result_tuple):
        failure_count = sum(row.success_status != "success" for row in result_tuple)
        raise RuntimeError(
            f"OpenMM system smoke failed for {failure_count}/{len(result_tuple)} representative systems; "
            f"see {_display_path(summary_path)}"
        )
    return result_tuple
