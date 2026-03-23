"""Phase 6A4 full-panel OpenMM system builds for baseline-ready systems."""

from __future__ import annotations

import csv
import traceback
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable

import yaml

from lanm.analysis.md_system_build import (
    DEFAULT_PROTEIN_FORCEFIELD,
    DEFAULT_WATER_MODEL,
    discover_md_input_manifest,
    load_md_protocol_build_config,
)
from lanm.analysis.metal_parameter_mapping import load_md_system_build_manifest
from lanm.analysis.metal_parameter_values import (
    GENERIC_12_6_4_HIGHLY_CHARGED,
    NumericParameterValueRecord,
    validate_numeric_parameter_registry,
)
from lanm.analysis.openmm_system_smoke import (
    OpenMMSystemBuildRequest,
    parse_openmm_system_build_request,
    prepare_structure_input,
)
from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_text, write_csv_rows, write_yaml
from lanm.md.openmm_system_build import (
    DEFAULT_CUTOFF_NM,
    DEFAULT_FRICTION_COEFF_PS,
    DEFAULT_HYDROGEN_MASS_REPARTITIONING,
    DEFAULT_NONBONDED_METHOD,
    DEFAULT_TIMESTEP_FS,
    OpenMMSerializedSystem,
    build_openmm_serialized_system,
)
from lanm.paths import REPO_ROOT

EXPECTED_BASELINE_SYSTEM_COUNT = 30


@dataclass(frozen=True, slots=True)
class MetalBuildReadinessRow:
    panel_rank: int
    panel_member_id: str
    target_metal: str
    topology_class: str
    build_request_path: str
    chosen_water_model: str
    generic_protein_water_preparation_ready: bool
    phase_6a2_custom_metal_parameters_still_required: bool
    baseline_parameter_family: str
    baseline_numeric_values_present: bool
    ready_for_openmm_system_build_baseline: bool
    tuned_parameter_family: str
    tuned_numeric_values_present: bool
    ready_for_openmm_system_build_tuned: bool
    readiness_summary: str


@dataclass(frozen=True, slots=True)
class OpenMMSystemBuildRuntimeConfig:
    replicate_count: int
    target_temperature_K: int
    friction_coeff_ps: float
    timestep_fs: float
    nonbonded_method: str
    cutoff_nm: float
    hydrogen_mass_repartitioning: bool


@dataclass(frozen=True, slots=True)
class BaselineReadyOpenMMSystem:
    panel_rank: int
    panel_member_id: str
    panel_member_type: str
    panel_role: str
    topology_class: str
    target_metal: str
    build_request_path: str
    source_starting_structure_pdb: str
    chosen_protein_forcefield: str
    chosen_water_model: str
    replicate_count: int
    target_temperature_K: int
    baseline_parameter_family: str
    metal_parameter_provenance: str
    rmin_half_A: float
    epsilon_kcal_per_mol: float
    c4_kcal_per_mol_A4: int


@dataclass(frozen=True, slots=True)
class SimulationConfig:
    panel_member_id: str
    target_metal: str
    topology_class: str
    replicate_count: int
    target_temperature_K: int
    friction_coeff_ps: float
    timestep_fs: float
    nonbonded_method: str
    cutoff_nm: float
    hydrogen_mass_repartitioning: bool
    protein_forcefield: str
    water_model: str
    metal_parameter_family: str
    metal_parameter_provenance: str


@dataclass(frozen=True, slots=True)
class OpenMMSystemBuildSummaryRow:
    panel_member_id: str
    target_metal: str
    topology_class: str
    success_status: str
    failure_reason: str
    atom_count: int
    residue_count: int
    system_xml_path: str
    simulation_config_path: str


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


def load_metal_build_readiness_rows(path: Path) -> tuple[MetalBuildReadinessRow, ...]:
    """Load the Phase 6A2c baseline/tuned readiness table for all panel systems."""
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
                chosen_water_model=str(row["chosen_water_model"]).strip(),
                generic_protein_water_preparation_ready=_parse_csv_bool(
                    row,
                    "generic_protein_water_preparation_ready",
                ),
                phase_6a2_custom_metal_parameters_still_required=_parse_csv_bool(
                    row,
                    "phase_6a2_custom_metal_parameters_still_required",
                ),
                baseline_parameter_family=str(row["baseline_parameter_family"]).strip(),
                baseline_numeric_values_present=_parse_csv_bool(row, "baseline_numeric_values_present"),
                ready_for_openmm_system_build_baseline=_parse_csv_bool(
                    row,
                    "ready_for_openmm_system_build_baseline",
                ),
                tuned_parameter_family=str(row["tuned_parameter_family"]).strip(),
                tuned_numeric_values_present=_parse_csv_bool(row, "tuned_numeric_values_present"),
                ready_for_openmm_system_build_tuned=_parse_csv_bool(
                    row,
                    "ready_for_openmm_system_build_tuned",
                ),
                readiness_summary=str(row["readiness_summary"]).strip(),
            )
            for row in reader
        )
    if not rows:
        raise ValueError(f"No rows found in {path}")
    return tuple(sorted(rows, key=lambda row: (row.panel_rank, row.panel_member_id, row.target_metal)))


def load_numeric_parameter_registry(path: Path) -> tuple[NumericParameterValueRecord, ...]:
    """Load and validate the audited numeric metal-parameter registry."""
    payload = _parse_required_mapping(path)
    numeric_values = payload.get("numeric_parameter_values")
    if not isinstance(numeric_values, list):
        raise ValueError(f"Expected numeric_parameter_values list in {path}")
    registry = tuple(
        NumericParameterValueRecord(
            metal_identity=str(entry["metal_identity"]).strip(),
            source_family_label=str(entry["source_family_label"]).strip(),
            water_model=str(entry["water_model"]).strip(),
            formal_charge=int(entry["formal_charge"]),
            rmin_half_A=float(entry["rmin_half_A"]),
            epsilon_kcal_per_mol=float(entry["epsilon_kcal_per_mol"]),
            c4_kcal_per_mol_A4=int(entry["c4_kcal_per_mol_A4"]),
            parameter_provenance=str(entry["parameter_provenance"]).strip(),
            notes=str(entry["notes"]).strip(),
        )
        for entry in numeric_values
    )
    return validate_numeric_parameter_registry(registry)


def load_openmm_system_build_runtime_config(path: Path) -> OpenMMSystemBuildRuntimeConfig:
    """Load the deterministic Phase 6A4 simulation defaults from the MD protocol."""
    payload = _parse_required_mapping(path)
    replicate_count = int(payload.get("replicate_count", 0) or 0)
    target_temperature_K = int(payload.get("target_temperature_K", 0) or 0)
    if replicate_count <= 0:
        raise ValueError(f"Expected positive replicate_count in {path}")
    if target_temperature_K <= 0:
        raise ValueError(f"Expected positive target_temperature_K in {path}")
    return OpenMMSystemBuildRuntimeConfig(
        replicate_count=replicate_count,
        target_temperature_K=target_temperature_K,
        friction_coeff_ps=DEFAULT_FRICTION_COEFF_PS,
        timestep_fs=DEFAULT_TIMESTEP_FS,
        nonbonded_method=DEFAULT_NONBONDED_METHOD,
        cutoff_nm=DEFAULT_CUTOFF_NM,
        hydrogen_mass_repartitioning=DEFAULT_HYDROGEN_MASS_REPARTITIONING,
    )


def discover_baseline_ready_openmm_systems(
    *,
    md_input_manifest_path: Path,
    md_system_build_manifest_path: Path,
    metal_build_readiness_path: Path,
    md_protocol_path: Path,
    metal_parameter_values_path: Path,
) -> tuple[BaselineReadyOpenMMSystem, ...]:
    """Resolve the full Phase 6A4 baseline-ready panel in deterministic order."""
    md_input_rows = discover_md_input_manifest(md_input_manifest_path)
    build_manifest_rows = load_md_system_build_manifest(md_system_build_manifest_path)
    readiness_rows = load_metal_build_readiness_rows(metal_build_readiness_path)
    build_config = load_md_protocol_build_config(md_protocol_path)
    runtime_config = load_openmm_system_build_runtime_config(md_protocol_path)
    numeric_parameter_registry = load_numeric_parameter_registry(metal_parameter_values_path)

    if build_config.protein_forcefield != DEFAULT_PROTEIN_FORCEFIELD:
        raise ValueError(
            f"Phase 6A4 requires {DEFAULT_PROTEIN_FORCEFIELD!r}, found {build_config.protein_forcefield!r}"
        )
    if build_config.water_model != DEFAULT_WATER_MODEL:
        raise ValueError(
            f"Phase 6A4 requires {DEFAULT_WATER_MODEL!r}, found {build_config.water_model!r}"
        )

    input_by_key = {(row.panel_member_id, row.target_metal): row for row in md_input_rows}
    build_by_key = {(row.panel_member_id, row.target_metal): row for row in build_manifest_rows}
    readiness_by_key = {(row.panel_member_id, row.target_metal): row for row in readiness_rows}
    parameter_by_key = {
        (record.source_family_label, record.metal_identity, record.water_model): record
        for record in numeric_parameter_registry
    }

    if set(input_by_key) != set(build_by_key):
        raise ValueError("md_input_manifest.csv and md_system_build_manifest.csv do not cover the same systems")
    if set(input_by_key) != set(readiness_by_key):
        raise ValueError("md_input_manifest.csv and metal_build_readiness.csv do not cover the same systems")

    selected_rows: list[BaselineReadyOpenMMSystem] = []
    for readiness_row in readiness_rows:
        if not readiness_row.ready_for_openmm_system_build_baseline:
            continue
        key = (readiness_row.panel_member_id, readiness_row.target_metal)
        input_row = input_by_key.get(key)
        build_row = build_by_key.get(key)
        if input_row is None:
            raise ValueError(f"Missing md_input_manifest row for {key}")
        if build_row is None:
            raise ValueError(f"Missing md_system_build_manifest row for {key}")
        if readiness_row.topology_class != input_row.topology_class:
            raise ValueError(
                f"Topology mismatch for {key}: readiness {readiness_row.topology_class!r} vs "
                f"input {input_row.topology_class!r}"
            )
        if build_row.topology_class != input_row.topology_class:
            raise ValueError(
                f"Topology mismatch for {key}: build {build_row.topology_class!r} vs "
                f"input {input_row.topology_class!r}"
            )
        if build_row.build_request_path != readiness_row.build_request_path:
            raise ValueError(
                f"Build-request mismatch for {key}: {build_row.build_request_path!r} vs "
                f"{readiness_row.build_request_path!r}"
            )
        if build_row.source_starting_structure_pdb != input_row.starting_structure_path:
            raise ValueError(
                f"Starting-structure mismatch for {key}: {build_row.source_starting_structure_pdb!r} vs "
                f"{input_row.starting_structure_path!r}"
            )
        if build_row.chosen_protein_forcefield != DEFAULT_PROTEIN_FORCEFIELD:
            raise ValueError(
                f"{key} must use {DEFAULT_PROTEIN_FORCEFIELD!r}, found {build_row.chosen_protein_forcefield!r}"
            )
        if build_row.chosen_water_model != DEFAULT_WATER_MODEL:
            raise ValueError(
                f"{key} must use {DEFAULT_WATER_MODEL!r}, found {build_row.chosen_water_model!r}"
            )
        if readiness_row.chosen_water_model != DEFAULT_WATER_MODEL:
            raise ValueError(
                f"{key} readiness row must use {DEFAULT_WATER_MODEL!r}, found {readiness_row.chosen_water_model!r}"
            )
        if not build_row.generic_protein_water_preparation_ready:
            raise ValueError(f"{key} is not generic protein/water preparation ready")
        if not readiness_row.generic_protein_water_preparation_ready:
            raise ValueError(f"{key} readiness row is not generic protein/water preparation ready")
        if readiness_row.baseline_parameter_family != GENERIC_12_6_4_HIGHLY_CHARGED:
            raise ValueError(
                f"{key} must use {GENERIC_12_6_4_HIGHLY_CHARGED!r}, "
                f"found {readiness_row.baseline_parameter_family!r}"
            )
        if not readiness_row.baseline_numeric_values_present:
            raise ValueError(f"{key} is baseline-ready without numeric values present")
        numeric_parameters = parameter_by_key.get(
            (
                readiness_row.baseline_parameter_family,
                readiness_row.target_metal,
                DEFAULT_WATER_MODEL,
            )
        )
        if numeric_parameters is None:
            raise ValueError(
                f"Missing numeric parameters for {key} under {readiness_row.baseline_parameter_family!r}"
            )
        if input_row.replicate_count != runtime_config.replicate_count:
            raise ValueError(
                f"{key} replicate_count {input_row.replicate_count} does not match protocol "
                f"{runtime_config.replicate_count}"
            )
        if input_row.target_temperature_K != runtime_config.target_temperature_K:
            raise ValueError(
                f"{key} target_temperature_K {input_row.target_temperature_K} does not match protocol "
                f"{runtime_config.target_temperature_K}"
            )
        selected_rows.append(
            BaselineReadyOpenMMSystem(
                panel_rank=input_row.panel_rank,
                panel_member_id=input_row.panel_member_id,
                panel_member_type=input_row.panel_member_type,
                panel_role=input_row.panel_role,
                topology_class=input_row.topology_class,
                target_metal=input_row.target_metal,
                build_request_path=build_row.build_request_path,
                source_starting_structure_pdb=input_row.starting_structure_path,
                chosen_protein_forcefield=build_row.chosen_protein_forcefield,
                chosen_water_model=build_row.chosen_water_model,
                replicate_count=input_row.replicate_count,
                target_temperature_K=input_row.target_temperature_K,
                baseline_parameter_family=readiness_row.baseline_parameter_family,
                metal_parameter_provenance=numeric_parameters.parameter_provenance,
                rmin_half_A=numeric_parameters.rmin_half_A,
                epsilon_kcal_per_mol=numeric_parameters.epsilon_kcal_per_mol,
                c4_kcal_per_mol_A4=numeric_parameters.c4_kcal_per_mol_A4,
            )
        )

    if len(selected_rows) != EXPECTED_BASELINE_SYSTEM_COUNT:
        raise ValueError(
            f"Expected {EXPECTED_BASELINE_SYSTEM_COUNT} baseline-ready systems, found {len(selected_rows)}"
        )
    return tuple(sorted(selected_rows, key=lambda row: (row.panel_rank, row.panel_member_id, row.target_metal)))


def build_simulation_config(
    *,
    system_row: BaselineReadyOpenMMSystem,
    runtime_config: OpenMMSystemBuildRuntimeConfig,
) -> SimulationConfig:
    """Create the deterministic per-system simulation config artifact."""
    return SimulationConfig(
        panel_member_id=system_row.panel_member_id,
        target_metal=system_row.target_metal,
        topology_class=system_row.topology_class,
        replicate_count=runtime_config.replicate_count,
        target_temperature_K=runtime_config.target_temperature_K,
        friction_coeff_ps=runtime_config.friction_coeff_ps,
        timestep_fs=runtime_config.timestep_fs,
        nonbonded_method=runtime_config.nonbonded_method,
        cutoff_nm=runtime_config.cutoff_nm,
        hydrogen_mass_repartitioning=runtime_config.hydrogen_mass_repartitioning,
        protein_forcefield=system_row.chosen_protein_forcefield,
        water_model=system_row.chosen_water_model,
        metal_parameter_family=system_row.baseline_parameter_family,
        metal_parameter_provenance=system_row.metal_parameter_provenance,
    )


def render_openmm_system_build_report(
    results: tuple[OpenMMSystemBuildSummaryRow, ...],
    runtime_config: OpenMMSystemBuildRuntimeConfig,
) -> str:
    """Render the concise Phase 6A4 Markdown build report."""
    success_count = sum(row.success_status == "success" for row in results)
    failure_count = len(results) - success_count
    lines = [
        "# OpenMM System Building",
        "",
        "Phase 6A4 builds deterministic baseline OpenMM `System` objects and simulation configs for all baseline-ready MD panel systems.",
        "",
        "## Build Summary",
        "",
        f"- baseline-ready systems attempted: `{len(results)}`",
        f"- successful OpenMM `System` serializations: `{success_count}`",
        f"- failures recorded: `{failure_count}`",
        f"- force-field files: `{DEFAULT_PROTEIN_FORCEFIELD}`, `{DEFAULT_WATER_MODEL}`",
        f"- baseline metal parameter family: `{GENERIC_12_6_4_HIGHLY_CHARGED}`",
        (
            "- integrator defaults: "
            f"`LangevinMiddleIntegrator`, `{runtime_config.target_temperature_K} K`, "
            f"`{runtime_config.friction_coeff_ps:.1f} ps^-1`, `{runtime_config.timestep_fs:.1f} fs`"
        ),
        (
            "- nonbonded defaults: "
            f"`{runtime_config.nonbonded_method}`, cutoff `{runtime_config.cutoff_nm:.1f} nm`, "
            f"`hydrogen_mass_repartitioning={runtime_config.hydrogen_mass_repartitioning}`"
        ),
        "",
        "| panel_member_id | target_metal | topology_class | success_status | atom_count | residue_count | failure_reason |",
        "| --- | --- | --- | --- | ---: | ---: | --- |",
    ]
    for row in results:
        lines.append(
            "| "
            f"{row.panel_member_id} | "
            f"{row.target_metal} | "
            f"{row.topology_class} | "
            f"{row.success_status} | "
            f"{row.atom_count} | "
            f"{row.residue_count} | "
            f"{row.failure_reason or '-'} |"
        )
    lines.append("")
    return "\n".join(lines)


def _write_system_build_log(
    *,
    log_path: Path,
    system_row: BaselineReadyOpenMMSystem,
    request: OpenMMSystemBuildRequest,
    runtime_config: OpenMMSystemBuildRuntimeConfig,
    restored_template_waters: int,
    template_water_reference_path: str,
    atom_count: int,
    residue_count: int,
    status: str,
    failure_reason: str,
    system_xml_path: str,
    simulation_config_path: str,
    traceback_text: str = "",
) -> None:
    lines = [
        "Phase 6A4 OpenMM system build",
        f"panel_member_id: {system_row.panel_member_id}",
        f"target_metal: {system_row.target_metal}",
        f"topology_class: {system_row.topology_class}",
        f"source_starting_structure_pdb: {system_row.source_starting_structure_pdb}",
        f"build_request_path: {system_row.build_request_path}",
        f"forcefield_files: {system_row.chosen_protein_forcefield}, {system_row.chosen_water_model}",
        f"baseline_parameter_family: {system_row.baseline_parameter_family}",
        f"metal_parameter_provenance: {system_row.metal_parameter_provenance}",
        (
            "metal_parameter_values: "
            f"rmin_half_A={system_row.rmin_half_A:.3f}, "
            f"epsilon_kcal_per_mol={system_row.epsilon_kcal_per_mol:.8f}, "
            f"c4_kcal_per_mol_A4={system_row.c4_kcal_per_mol_A4}"
        ),
        f"replicate_count: {runtime_config.replicate_count}",
        f"target_temperature_K: {runtime_config.target_temperature_K}",
        f"friction_coeff_ps: {runtime_config.friction_coeff_ps:.1f}",
        f"timestep_fs: {runtime_config.timestep_fs:.1f}",
        f"nonbonded_method: {runtime_config.nonbonded_method}",
        f"cutoff_nm: {runtime_config.cutoff_nm:.1f}",
        f"hydrogen_mass_repartitioning: {runtime_config.hydrogen_mass_repartitioning}",
        (
            "restore_template_metal_proximal_waters_before_solvation: "
            f"{request.restore_template_metal_proximal_waters_before_solvation}"
        ),
        f"template_water_reference_path: {template_water_reference_path or '-'}",
        f"restored_template_waters: {restored_template_waters}",
        f"atom_count: {atom_count}",
        f"residue_count: {residue_count}",
        f"system_xml_path: {system_xml_path or '-'}",
        f"simulation_config_path: {simulation_config_path or '-'}",
        f"status: {status}",
    ]
    if failure_reason:
        lines.append(f"failure_reason: {failure_reason}")
    if traceback_text:
        lines.extend(("", "traceback:", traceback_text.rstrip()))
    atomic_write_text(log_path, "\n".join(lines) + "\n")


def run_openmm_system_build(
    *,
    md_input_manifest_path: Path,
    md_system_build_manifest_path: Path,
    metal_build_readiness_path: Path,
    md_protocol_path: Path,
    metal_parameter_values_path: Path,
    output_root: Path,
    summary_path: Path,
    report_path: Path,
    system_builder: Callable[..., OpenMMSerializedSystem] = build_openmm_serialized_system,
) -> tuple[OpenMMSystemBuildSummaryRow, ...]:
    """Build deterministic baseline OpenMM systems for the full baseline-ready panel."""
    runtime_config = load_openmm_system_build_runtime_config(md_protocol_path)
    baseline_systems = discover_baseline_ready_openmm_systems(
        md_input_manifest_path=md_input_manifest_path,
        md_system_build_manifest_path=md_system_build_manifest_path,
        metal_build_readiness_path=metal_build_readiness_path,
        md_protocol_path=md_protocol_path,
        metal_parameter_values_path=metal_parameter_values_path,
    )

    results: list[OpenMMSystemBuildSummaryRow] = []
    for system_row in baseline_systems:
        output_dir = output_root / system_row.panel_member_id / system_row.target_metal
        prepared_structure_path = output_dir / "prepared_structure.pdb"
        system_xml_path = output_dir / "system.xml"
        integrator_xml_path = output_dir / "integrator.xml"
        simulation_config_path = output_dir / "simulation_config.yaml"
        log_path = output_dir / "system_build_log.txt"
        for artifact_path in (
            prepared_structure_path,
            system_xml_path,
            integrator_xml_path,
            simulation_config_path,
        ):
            artifact_path.unlink(missing_ok=True)

        request: OpenMMSystemBuildRequest | None = None
        restored_template_waters = 0
        template_water_reference_path = ""
        atom_count = 0
        residue_count = 0
        system_xml_display_path = ""
        simulation_config_display_path = ""
        try:
            request = parse_openmm_system_build_request(_repo_path(system_row.build_request_path))
            prepared_input = prepare_structure_input(request)
            restored_template_waters = prepared_input.restored_template_waters
            template_water_reference_path = prepared_input.template_water_reference_path
            built_system = system_builder(
                prepared_structure_pdb_text=prepared_input.pdb_text,
                forcefield_files=(system_row.chosen_protein_forcefield, system_row.chosen_water_model),
                target_temperature_K=runtime_config.target_temperature_K,
                friction_coeff_ps=runtime_config.friction_coeff_ps,
                timestep_fs=runtime_config.timestep_fs,
                nonbonded_method=runtime_config.nonbonded_method,
                cutoff_nm=runtime_config.cutoff_nm,
                hydrogen_mass_repartitioning=runtime_config.hydrogen_mass_repartitioning,
            )
            atom_count = built_system.atom_count
            residue_count = built_system.residue_count
            simulation_config = build_simulation_config(
                system_row=system_row,
                runtime_config=runtime_config,
            )
            atomic_write_text(prepared_structure_path, built_system.prepared_structure_pdb_text)
            atomic_write_text(system_xml_path, built_system.system_xml)
            atomic_write_text(integrator_xml_path, built_system.integrator_xml)
            write_yaml(simulation_config_path, asdict(simulation_config))
            system_xml_display_path = _display_path(system_xml_path)
            simulation_config_display_path = _display_path(simulation_config_path)
            _write_system_build_log(
                log_path=log_path,
                system_row=system_row,
                request=request,
                runtime_config=runtime_config,
                restored_template_waters=restored_template_waters,
                template_water_reference_path=template_water_reference_path,
                atom_count=atom_count,
                residue_count=residue_count,
                status="success",
                failure_reason="",
                system_xml_path=system_xml_display_path,
                simulation_config_path=simulation_config_display_path,
            )
            results.append(
                OpenMMSystemBuildSummaryRow(
                    panel_member_id=system_row.panel_member_id,
                    target_metal=system_row.target_metal,
                    topology_class=system_row.topology_class,
                    success_status="success",
                    failure_reason="",
                    atom_count=atom_count,
                    residue_count=residue_count,
                    system_xml_path=system_xml_display_path,
                    simulation_config_path=simulation_config_display_path,
                )
            )
        except Exception as exc:
            failure_reason = f"{type(exc).__name__}: {exc}"
            traceback_text = traceback.format_exc()
            if request is None:
                request = OpenMMSystemBuildRequest(
                    panel_member_id=system_row.panel_member_id,
                    target_metal=system_row.target_metal,
                    topology_class=system_row.topology_class,
                    source_starting_structure_pdb=system_row.source_starting_structure_pdb,
                    chosen_protein_forcefield=system_row.chosen_protein_forcefield,
                    chosen_water_model=system_row.chosen_water_model,
                    restore_template_metal_proximal_waters_before_solvation=False,
                    template_metal_proximal_waters_available=False,
                    comparison_panel_member_id="",
                    comparison_source="",
                    comparison_preserved_solvent_residue_count=0,
                )
            _write_system_build_log(
                log_path=log_path,
                system_row=system_row,
                request=request,
                runtime_config=runtime_config,
                restored_template_waters=restored_template_waters,
                template_water_reference_path=template_water_reference_path,
                atom_count=atom_count,
                residue_count=residue_count,
                status="failure",
                failure_reason=failure_reason,
                system_xml_path="",
                simulation_config_path="",
                traceback_text=traceback_text,
            )
            results.append(
                OpenMMSystemBuildSummaryRow(
                    panel_member_id=system_row.panel_member_id,
                    target_metal=system_row.target_metal,
                    topology_class=system_row.topology_class,
                    success_status="failure",
                    failure_reason=failure_reason,
                    atom_count=atom_count,
                    residue_count=residue_count,
                    system_xml_path="",
                    simulation_config_path="",
                )
            )

    result_tuple = tuple(results)
    write_csv_rows(summary_path, result_tuple)
    atomic_write_text(report_path, render_openmm_system_build_report(result_tuple, runtime_config))
    if any(row.success_status != "success" for row in result_tuple):
        failure_count = sum(row.success_status != "success" for row in result_tuple)
        raise RuntimeError(
            f"OpenMM system build failed for {failure_count}/{len(result_tuple)} baseline-ready systems; "
            f"see {_display_path(summary_path)}"
        )
    return result_tuple
