"""Phase 7G1 reduced OpenMM MD rescreen for the three round-2 finalists."""

from __future__ import annotations

import csv
import os
import traceback
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import median
from typing import Callable, Mapping

from lanm.analysis.md_input_preparation import (
    SUGGESTED_COLLECTIVE_VARIABLES,
    TARGET_METALS,
    TARGET_TEMPERATURE_K,
    render_target_metal_starting_structure,
    resolve_panel_member_structure_template,
)
from lanm.analysis.md_panel_selection import MDPanelRow
from lanm.analysis.md_system_build import DEFAULT_PROTEIN_FORCEFIELD, DEFAULT_WATER_MODEL, load_md_protocol_build_config
from lanm.analysis.metal_parameter_values import GENERIC_12_6_4_HIGHLY_CHARGED
from lanm.analysis.openmm_equilibration_smoke import build_equilibration_runtime_config, load_simulation_config
from lanm.analysis.openmm_screening import (
    DEFAULT_SCREENING_CPU_THREAD_COUNT,
    DEFAULT_SCREENING_WORKER_COUNT,
    OpenMMScreeningRuntimeConfig,
    ScreeningReplicateSpec,
    StateDataSummary,
    build_openmm_screening_replicate_schedule,
    build_openmm_screening_runtime_config,
    summarize_state_data_csv,
)
from lanm.analysis.openmm_system_build import (
    SimulationConfig,
    load_numeric_parameter_registry,
    load_openmm_system_build_runtime_config,
)
from lanm.configuration import load_project_config
from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_text, write_csv_rows, write_yaml
from lanm.md.openmm_equilibration_smoke import OpenMMEquilibrationSmokeArtifacts, run_openmm_equilibration_smoke_protocol
from lanm.md.openmm_screening import OpenMMScreeningArtifacts, run_openmm_screening_replicate
from lanm.md.openmm_system_build import OpenMMSerializedSystem, build_openmm_serialized_system
from lanm.paths import REPO_ROOT

EXPECTED_ROUND2_FINALIST_IDS = (
    "am1_mex_ss_only_u02_r2u02",
    "hans_pocket_ss_only_u02_r2u01",
    "hans_interface_ss_plus_if_u04_r2u01",
)


@dataclass(frozen=True, slots=True)
class Round2OpenMMScreeningCandidate:
    panel_rank: int
    panel_member_id: str
    panel_role: str
    candidate_id: str
    campaign_id: str
    backbone_id: str
    topology_class: str
    preserved_metal_identity: str
    rosetta_score_rank_within_topology: int
    representative_design_id: int
    representative_sequence_id: str
    representative_ligand_confidence: float
    representative_overall_confidence: float
    total_score: float
    representative_packed_pdb_path: Path
    selection_reason: str


@dataclass(frozen=True, slots=True)
class Round2OpenMMSystem:
    panel_rank: int
    candidate_id: str
    panel_role: str
    campaign_id: str
    backbone_id: str
    topology_class: str
    preserved_metal_identity: str
    representative_packed_pdb_path: Path
    target_metal: str
    output_dir: Path
    starting_structure_path: Path
    system_metadata_path: Path
    prepared_structure_path: Path
    system_xml_path: Path
    integrator_xml_path: Path
    simulation_config_path: Path
    system_build_log_path: Path
    equilibration_output_dir: Path
    minimized_structure_path: Path
    nvt_final_path: Path
    npt_final_path: Path
    equilibration_log_path: Path
    equilibration_final_state_path: Path


@dataclass(frozen=True, slots=True)
class Round2BuildContext:
    protein_forcefield: str
    water_model: str
    parameter_provenance_by_metal: Mapping[str, str]


@dataclass(frozen=True, slots=True)
class Round2OpenMMScreeningReplicateSummaryRow:
    candidate_id: str
    target_metal: str
    topology_class: str
    replicate_id: int
    seed: int
    success_status: str
    failure_reason: str
    mean_temperature_K: float | None
    final_temperature_K: float | None
    mean_box_volume_nm3: float | None
    mean_min_metal_oxygen_distance_A: float | None
    inner_sphere_occupancy_fraction: float | None
    state_data_path: str
    final_state_path: str
    final_structure_path: str
    screening_log_path: str
    trajectory_path: str


@dataclass(frozen=True, slots=True)
class Round2OpenMMSystemSummaryRow:
    candidate_id: str
    target_metal: str
    topology_class: str
    replicate_success_count: int
    median_mean_min_metal_oxygen_distance_A: float | None
    median_inner_sphere_occupancy_fraction: float | None
    screening_status: str
    failure_reason: str


@dataclass(frozen=True, slots=True)
class Round2OpenMMPanelStatusRow:
    candidate_id: str
    dy_status: str
    nd_status: str
    y_status: str
    al_status: str
    fe_status: str
    candidate_keep_for_metadynamics: bool


@dataclass(frozen=True, slots=True)
class Round2OpenMMScreeningReplicateTask:
    system: Round2OpenMMSystem
    runtime_config: OpenMMScreeningRuntimeConfig
    replicate_spec: ScreeningReplicateSpec
    cpu_thread_count: int = DEFAULT_SCREENING_CPU_THREAD_COUNT


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


def _parse_optional_float(text: str) -> float | None:
    stripped = str(text).strip()
    if not stripped:
        return None
    lowered = stripped.lower()
    if lowered in {"nan", "none", "null", "-"}:
        return None
    return float(stripped)


def _median_or_none(values: list[float | None]) -> float | None:
    present = [value for value in values if value is not None]
    if not present:
        return None
    return float(median(present))


def _format_metric(value: float | None) -> str:
    if value is None:
        return "-"
    return f"{value:.6f}"


def _path_if_exists(path: Path) -> str:
    if not path.exists():
        return ""
    return _display_path(path)


def _screening_artifact_paths(
    system: Round2OpenMMSystem,
    replicate_id: int,
) -> tuple[Path, Path, Path, Path, Path]:
    output_dir = system.output_dir / f"replicate_{replicate_id}"
    return (
        output_dir / "state_data.csv",
        output_dir / "final_state.xml",
        output_dir / "final_structure.pdb",
        output_dir / "screening_log.txt",
        output_dir / "trajectory.dcd",
    )


def discover_round2_openmm_screening_candidates(
    panel_path: Path,
) -> tuple[Round2OpenMMScreeningCandidate, ...]:
    """Load the reduced Phase 7F panel and validate the three expected finalists."""
    require_path(panel_path)
    with panel_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            Round2OpenMMScreeningCandidate(
                panel_rank=int(str(row["panel_rank"]).strip()),
                panel_member_id=str(row["panel_member_id"]).strip(),
                panel_role=str(row["panel_role"]).strip(),
                candidate_id=str(row["candidate_id"]).strip(),
                campaign_id=str(row["campaign_id"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                topology_class=str(row["topology_class"]).strip(),
                preserved_metal_identity=str(row["preserved_metal_identity"]).strip().upper(),
                rosetta_score_rank_within_topology=int(str(row["rosetta_score_rank_within_topology"]).strip()),
                representative_design_id=int(str(row["representative_design_id"]).strip()),
                representative_sequence_id=str(row["representative_sequence_id"]).strip(),
                representative_ligand_confidence=float(str(row["representative_ligand_confidence"]).strip()),
                representative_overall_confidence=float(str(row["representative_overall_confidence"]).strip()),
                total_score=float(str(row["total_score"]).strip()),
                representative_packed_pdb_path=_repo_path(str(row["representative_packed_pdb"]).strip()),
                selection_reason=str(row["selection_reason"]).strip(),
            )
            for row in reader
        ]
    rows = sorted(rows, key=lambda row: (row.panel_rank, row.candidate_id))
    if len(rows) != len(EXPECTED_ROUND2_FINALIST_IDS):
        raise ValueError(
            f"Expected {len(EXPECTED_ROUND2_FINALIST_IDS)} round-2 finalists in {panel_path}, found {len(rows)}"
        )
    observed_ids = tuple(row.candidate_id for row in rows)
    if observed_ids != EXPECTED_ROUND2_FINALIST_IDS:
        raise ValueError(
            f"Round-2 finalist set mismatch in {panel_path}: "
            f"expected {EXPECTED_ROUND2_FINALIST_IDS}, found {observed_ids}"
        )
    for row in rows:
        require_path(row.representative_packed_pdb_path)
    return tuple(rows)


def build_round2_openmm_systems(
    *,
    candidates: tuple[Round2OpenMMScreeningCandidate, ...],
    output_root: Path,
) -> tuple[Round2OpenMMSystem, ...]:
    """Expand the reduced finalist panel into the 15 target-metal systems."""
    systems: list[Round2OpenMMSystem] = []
    for candidate in candidates:
        for target_metal in TARGET_METALS:
            output_dir = output_root / candidate.candidate_id / target_metal
            equilibration_output_dir = output_dir / "equilibration_smoke"
            systems.append(
                Round2OpenMMSystem(
                    panel_rank=candidate.panel_rank,
                    candidate_id=candidate.candidate_id,
                    panel_role=candidate.panel_role,
                    campaign_id=candidate.campaign_id,
                    backbone_id=candidate.backbone_id,
                    topology_class=candidate.topology_class,
                    preserved_metal_identity=candidate.preserved_metal_identity,
                    representative_packed_pdb_path=candidate.representative_packed_pdb_path,
                    target_metal=target_metal,
                    output_dir=output_dir,
                    starting_structure_path=output_dir / "starting_structure.pdb",
                    system_metadata_path=output_dir / "system_metadata.yaml",
                    prepared_structure_path=output_dir / "prepared_structure.pdb",
                    system_xml_path=output_dir / "system.xml",
                    integrator_xml_path=output_dir / "integrator.xml",
                    simulation_config_path=output_dir / "simulation_config.yaml",
                    system_build_log_path=output_dir / "system_build_log.txt",
                    equilibration_output_dir=equilibration_output_dir,
                    minimized_structure_path=equilibration_output_dir / "minimized_structure.pdb",
                    nvt_final_path=equilibration_output_dir / "nvt_final.pdb",
                    npt_final_path=equilibration_output_dir / "npt_final.pdb",
                    equilibration_log_path=equilibration_output_dir / "equilibration_log.csv",
                    equilibration_final_state_path=equilibration_output_dir / "final_state.xml",
                )
            )
    return tuple(systems)


def resolve_round2_build_context(
    *,
    md_protocol_path: Path,
    metal_parameter_values_path: Path,
) -> Round2BuildContext:
    """Resolve the shared baseline OpenMM build configuration for Phase 7G1."""
    build_config = load_md_protocol_build_config(md_protocol_path)
    runtime_config = load_openmm_system_build_runtime_config(md_protocol_path)
    if runtime_config.target_temperature_K != TARGET_TEMPERATURE_K:
        raise ValueError(
            f"Expected target_temperature_K={TARGET_TEMPERATURE_K} in {md_protocol_path}, "
            f"found {runtime_config.target_temperature_K}"
        )
    if build_config.protein_forcefield != DEFAULT_PROTEIN_FORCEFIELD:
        raise ValueError(
            f"Round-2 screening requires {DEFAULT_PROTEIN_FORCEFIELD!r}, "
            f"found {build_config.protein_forcefield!r}"
        )
    if build_config.water_model != DEFAULT_WATER_MODEL:
        raise ValueError(
            f"Round-2 screening requires {DEFAULT_WATER_MODEL!r}, found {build_config.water_model!r}"
        )
    numeric_registry = load_numeric_parameter_registry(metal_parameter_values_path)
    provenance_by_metal: dict[str, str] = {}
    for record in numeric_registry:
        key = (record.source_family_label, record.metal_identity, record.water_model)
        if key != (GENERIC_12_6_4_HIGHLY_CHARGED, record.metal_identity, DEFAULT_WATER_MODEL):
            continue
        provenance_by_metal[record.metal_identity] = record.parameter_provenance
    if tuple(sorted(provenance_by_metal, key=TARGET_METALS.index)) != TARGET_METALS:
        raise ValueError(
            "Missing baseline numeric parameter provenance entries for one or more target metals"
        )
    return Round2BuildContext(
        protein_forcefield=build_config.protein_forcefield,
        water_model=build_config.water_model,
        parameter_provenance_by_metal=provenance_by_metal,
    )


def _candidate_to_md_panel_row(candidate: Round2OpenMMScreeningCandidate) -> MDPanelRow:
    return MDPanelRow(
        panel_rank=candidate.panel_rank,
        panel_member_id=candidate.candidate_id,
        panel_member_type="designed_candidate",
        panel_role=candidate.panel_role,
        candidate_id=candidate.candidate_id,
        campaign_id=candidate.campaign_id,
        backbone_id=candidate.backbone_id,
        topology_class=candidate.topology_class,
        preserved_metal_identity=candidate.preserved_metal_identity,
        proteinmpnn_rank=None,
        ligandmpnn_representative_ligand_confidence=candidate.representative_ligand_confidence,
        ligandmpnn_representative_overall_confidence=candidate.representative_overall_confidence,
        ligandmpnn_mean_ligand_confidence=None,
        rosetta_score_rank_within_topology=candidate.rosetta_score_rank_within_topology,
        rosetta_total_score=candidate.total_score,
        starting_structure_path=_display_path(candidate.representative_packed_pdb_path),
        selection_reason=candidate.selection_reason,
    )


def render_round2_system_metadata(
    *,
    system: Round2OpenMMSystem,
    source_structure: str,
    metal_site_template_source: str,
    metal_site_count: int,
    preserved_solvent_residue_count: int,
) -> dict[str, object]:
    payload: dict[str, object] = {
        "candidate_id": system.candidate_id,
        "panel_role": system.panel_role,
        "campaign_id": system.campaign_id,
        "backbone_id": system.backbone_id,
        "topology_class": system.topology_class,
        "target_metal": system.target_metal,
        "source_structure": source_structure,
        "starting_structure_path": _display_path(system.starting_structure_path),
        "replicate_count": 3,
        "target_temperature_K": TARGET_TEMPERATURE_K,
        "suggested_collective_variables": list(SUGGESTED_COLLECTIVE_VARIABLES),
        "metal_site_count": metal_site_count,
        "preserved_solvent_residue_count": preserved_solvent_residue_count,
    }
    if metal_site_template_source != source_structure:
        payload["metal_site_template_source"] = metal_site_template_source
    return payload


def prepare_round2_md_inputs(
    *,
    candidates: tuple[Round2OpenMMScreeningCandidate, ...],
    systems: tuple[Round2OpenMMSystem, ...],
) -> None:
    """Write deterministic round-2 multi-metal starting structures and metadata."""
    project_config = load_project_config()
    configured_metals = (project_config.target_metal, *project_config.competitors)
    if configured_metals != TARGET_METALS:
        raise ValueError(f"Expected project metals {TARGET_METALS}, found {configured_metals}")
    if project_config.temperature_K != TARGET_TEMPERATURE_K:
        raise ValueError(
            f"Expected project temperature {TARGET_TEMPERATURE_K} K, found {project_config.temperature_K} K"
        )

    templates_by_candidate: dict[str, object] = {}
    for candidate in candidates:
        templates_by_candidate[candidate.candidate_id] = resolve_panel_member_structure_template(
            panel_member=_candidate_to_md_panel_row(candidate),
            design_backbone_sources={},
            nearby_solvent_cutoff_A=project_config.second_sphere_cutoff_A,
        )

    for system in systems:
        template = templates_by_candidate[system.candidate_id]
        prepared_structure = render_target_metal_starting_structure(template, target_metal=system.target_metal)
        atomic_write_text(system.starting_structure_path, prepared_structure.pdb_text)
        write_yaml(
            system.system_metadata_path,
            render_round2_system_metadata(
                system=system,
                source_structure=prepared_structure.source_structure,
                metal_site_template_source=prepared_structure.metal_site_template_source,
                metal_site_count=prepared_structure.metal_site_count,
                preserved_solvent_residue_count=prepared_structure.preserved_solvent_residue_count,
            ),
        )


def build_round2_simulation_config(
    *,
    system: Round2OpenMMSystem,
    build_context: Round2BuildContext,
    md_protocol_path: Path,
) -> SimulationConfig:
    """Build the deterministic per-system simulation config used by build/equilibration/screening."""
    runtime_config = load_openmm_system_build_runtime_config(md_protocol_path)
    return SimulationConfig(
        panel_member_id=system.candidate_id,
        target_metal=system.target_metal,
        topology_class=system.topology_class,
        replicate_count=runtime_config.replicate_count,
        target_temperature_K=runtime_config.target_temperature_K,
        friction_coeff_ps=runtime_config.friction_coeff_ps,
        timestep_fs=runtime_config.timestep_fs,
        nonbonded_method=runtime_config.nonbonded_method,
        cutoff_nm=runtime_config.cutoff_nm,
        hydrogen_mass_repartitioning=runtime_config.hydrogen_mass_repartitioning,
        protein_forcefield=build_context.protein_forcefield,
        water_model=build_context.water_model,
        metal_parameter_family=GENERIC_12_6_4_HIGHLY_CHARGED,
        metal_parameter_provenance=build_context.parameter_provenance_by_metal[system.target_metal],
    )


def build_round2_openmm_replicate_tasks(
    systems: tuple[Round2OpenMMSystem, ...],
) -> tuple[Round2OpenMMScreeningReplicateTask, ...]:
    """Expand the prepared systems into deterministic 3-replicate screening tasks."""
    tasks: list[Round2OpenMMScreeningReplicateTask] = []
    for system in systems:
        runtime_config = build_openmm_screening_runtime_config(system.simulation_config_path)
        for replicate_spec in build_openmm_screening_replicate_schedule(runtime_config):
            tasks.append(
                Round2OpenMMScreeningReplicateTask(
                    system=system,
                    runtime_config=runtime_config,
                    replicate_spec=replicate_spec,
                )
            )
    return tuple(tasks)


def classify_round2_screening_status(
    *,
    target_metal: str,
    replicate_success_count: int,
    median_inner_sphere_occupancy_fraction: float | None,
) -> str:
    """Apply the round-2 Phase 7G1 system-level classification rules."""
    if replicate_success_count < 3:
        return "failed"
    occupancy = median_inner_sphere_occupancy_fraction or 0.0
    if target_metal in {"Dy", "Nd", "Y"}:
        if occupancy >= 0.5:
            return "stable_bound"
        return "weak_binding_flag"
    if target_metal in {"Al", "Fe"}:
        if occupancy >= 0.5:
            return "persistent_capture_flag"
        return "no_persistent_capture"
    raise ValueError(f"Unsupported target metal {target_metal!r}")


def build_round2_openmm_screening_system_summaries(
    *,
    systems: tuple[Round2OpenMMSystem, ...],
    replicate_rows: tuple[Round2OpenMMScreeningReplicateSummaryRow, ...],
    stage_failure_reasons: Mapping[tuple[str, str], str],
) -> tuple[Round2OpenMMSystemSummaryRow, ...]:
    """Aggregate the 45 replicate attempts into the required 15-system summary table."""
    replicate_rows_by_key: dict[tuple[str, str], list[Round2OpenMMScreeningReplicateSummaryRow]] = defaultdict(list)
    for row in replicate_rows:
        replicate_rows_by_key[(row.candidate_id, row.target_metal)].append(row)

    system_rows: list[Round2OpenMMSystemSummaryRow] = []
    for system in systems:
        key = (system.candidate_id, system.target_metal)
        stage_failure_reason = stage_failure_reasons.get(key, "")
        if stage_failure_reason:
            system_rows.append(
                Round2OpenMMSystemSummaryRow(
                    candidate_id=system.candidate_id,
                    target_metal=system.target_metal,
                    topology_class=system.topology_class,
                    replicate_success_count=0,
                    median_mean_min_metal_oxygen_distance_A=None,
                    median_inner_sphere_occupancy_fraction=None,
                    screening_status="failed",
                    failure_reason=stage_failure_reason,
                )
            )
            continue

        system_replicates = sorted(
            replicate_rows_by_key.get(key, []),
            key=lambda row: row.replicate_id,
        )
        successful_rows = [row for row in system_replicates if row.success_status == "success"]
        replicate_success_count = len(successful_rows)
        median_distance = _median_or_none(
            [row.mean_min_metal_oxygen_distance_A for row in successful_rows]
        )
        median_occupancy = _median_or_none(
            [row.inner_sphere_occupancy_fraction for row in successful_rows]
        )
        failure_reasons = [
            f"replicate_{row.replicate_id}: {row.failure_reason}"
            for row in system_replicates
            if row.success_status != "success" and row.failure_reason
        ]
        system_rows.append(
            Round2OpenMMSystemSummaryRow(
                candidate_id=system.candidate_id,
                target_metal=system.target_metal,
                topology_class=system.topology_class,
                replicate_success_count=replicate_success_count,
                median_mean_min_metal_oxygen_distance_A=median_distance,
                median_inner_sphere_occupancy_fraction=median_occupancy,
                screening_status=classify_round2_screening_status(
                    target_metal=system.target_metal,
                    replicate_success_count=replicate_success_count,
                    median_inner_sphere_occupancy_fraction=median_occupancy,
                ),
                failure_reason="; ".join(failure_reasons),
            )
        )
    return tuple(system_rows)


def build_round2_openmm_panel_status_rows(
    system_rows: tuple[Round2OpenMMSystemSummaryRow, ...],
) -> tuple[Round2OpenMMPanelStatusRow, ...]:
    """Aggregate the required candidate-level keep/drop signal for later metadynamics."""
    rows_by_candidate: dict[str, dict[str, Round2OpenMMSystemSummaryRow]] = defaultdict(dict)
    for row in system_rows:
        rows_by_candidate[row.candidate_id][row.target_metal] = row

    panel_rows: list[Round2OpenMMPanelStatusRow] = []
    for candidate_id in sorted(rows_by_candidate):
        statuses = {metal: rows_by_candidate[candidate_id][metal].screening_status for metal in rows_by_candidate[candidate_id]}
        panel_rows.append(
            Round2OpenMMPanelStatusRow(
                candidate_id=candidate_id,
                dy_status=statuses.get("Dy", "failed"),
                nd_status=statuses.get("Nd", "failed"),
                y_status=statuses.get("Y", "failed"),
                al_status=statuses.get("Al", "failed"),
                fe_status=statuses.get("Fe", "failed"),
                candidate_keep_for_metadynamics=(
                    statuses.get("Dy") == "stable_bound"
                    and statuses.get("Al") == "no_persistent_capture"
                    and statuses.get("Fe") == "no_persistent_capture"
                ),
            )
        )
    return tuple(panel_rows)


def render_round2_openmm_screening_report(
    *,
    candidates: tuple[Round2OpenMMScreeningCandidate, ...],
    replicate_rows: tuple[Round2OpenMMScreeningReplicateSummaryRow, ...],
    system_rows: tuple[Round2OpenMMSystemSummaryRow, ...],
    panel_rows: tuple[Round2OpenMMPanelStatusRow, ...],
    runtime_config: OpenMMScreeningRuntimeConfig,
) -> str:
    """Render the Phase 7G1 Markdown report."""
    successful_replicates = sum(row.success_status == "success" for row in replicate_rows)
    failed_replicates = len(replicate_rows) - successful_replicates
    lines = [
        "# Round-2 OpenMM Screening",
        "",
        "Phase 7G1 runs a reduced deterministic OpenMM MD rescreen across the three round-2 finalists to compare Dy retention against Nd/Y and persistent Al/Fe capture.",
        "",
        "## Panel",
        "",
    ]
    for candidate in candidates:
        lines.append(
            f"- `{candidate.candidate_id}` ({candidate.topology_class}) from `{_display_path(candidate.representative_packed_pdb_path)}`"
        )
    lines.extend(
        [
            "",
            "## Screening Protocol",
            "",
            f"- finalists screened: `{len(candidates)}`",
            f"- systems screened: `{len(system_rows)}`",
            f"- replicate attempts: `{len(replicate_rows)}`",
            f"- successful replicates: `{successful_replicates}`",
            f"- failed replicates: `{failed_replicates}`",
            "- target metals per finalist: `Dy`, `Nd`, `Y`, `Al`, `Fe`",
            "- replicate schedule: `1->101`, `2->102`, `3->103`",
            f"- target temperature: `{runtime_config.target_temperature_K} K`",
            f"- pressure target: `{runtime_config.pressure_atm:.1f} atm`",
            f"- timestep: `{runtime_config.timestep_fs:.1f} fs`",
            f"- production duration per replicate: `{runtime_config.production_duration_ps:.1f} ps`",
            f"- report stride: `{runtime_config.report_stride_steps}` steps",
            f"- inner-sphere cutoff: `{runtime_config.inner_sphere_cutoff_A:.1f} A`",
            "- excluded in this phase: `metadynamics`, `QM`, and quantum steps",
            "",
            "## System Summary",
            "",
            "| candidate_id | target_metal | topology_class | replicate_success_count | median_mean_min_metal_oxygen_distance_A | median_inner_sphere_occupancy_fraction | screening_status | failure_reason |",
            "| --- | --- | --- | ---: | ---: | ---: | --- | --- |",
        ]
    )
    for row in system_rows:
        lines.append(
            "| "
            f"{row.candidate_id} | "
            f"{row.target_metal} | "
            f"{row.topology_class} | "
            f"{row.replicate_success_count} | "
            f"{_format_metric(row.median_mean_min_metal_oxygen_distance_A)} | "
            f"{_format_metric(row.median_inner_sphere_occupancy_fraction)} | "
            f"{row.screening_status} | "
            f"{row.failure_reason or '-'} |"
        )
    lines.extend(
        [
            "",
            "## Candidate Panel Status",
            "",
            "| candidate_id | Dy | Nd | Y | Al | Fe | keep_for_metadynamics |",
            "| --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in panel_rows:
        lines.append(
            "| "
            f"{row.candidate_id} | "
            f"{row.dy_status} | "
            f"{row.nd_status} | "
            f"{row.y_status} | "
            f"{row.al_status} | "
            f"{row.fe_status} | "
            f"{row.candidate_keep_for_metadynamics} |"
        )
    failure_rows = [row for row in system_rows if row.failure_reason]
    if failure_rows:
        lines.extend(
            [
                "",
                "## Failures",
                "",
            ]
        )
        for row in failure_rows:
            lines.append(
                f"- `{row.candidate_id}` / `{row.target_metal}`: {row.failure_reason}"
            )
    lines.append("")
    return "\n".join(lines)


def _write_round2_build_log(
    *,
    system: Round2OpenMMSystem,
    simulation_config: SimulationConfig,
    build_context: Round2BuildContext,
    atom_count: int,
    residue_count: int,
    status: str,
    failure_reason: str,
    traceback_text: str = "",
) -> None:
    lines = [
        "Phase 7G1 round-2 OpenMM system build",
        f"candidate_id: {system.candidate_id}",
        f"target_metal: {system.target_metal}",
        f"topology_class: {system.topology_class}",
        f"source_starting_structure_pdb: {_display_path(system.starting_structure_path)}",
        f"forcefield_files: {build_context.protein_forcefield}, {build_context.water_model}",
        f"metal_parameter_family: {simulation_config.metal_parameter_family}",
        f"metal_parameter_provenance: {simulation_config.metal_parameter_provenance}",
        f"replicate_count: {simulation_config.replicate_count}",
        f"target_temperature_K: {simulation_config.target_temperature_K}",
        f"friction_coeff_ps: {simulation_config.friction_coeff_ps:.1f}",
        f"timestep_fs: {simulation_config.timestep_fs:.1f}",
        f"nonbonded_method: {simulation_config.nonbonded_method}",
        f"cutoff_nm: {simulation_config.cutoff_nm:.1f}",
        f"hydrogen_mass_repartitioning: {simulation_config.hydrogen_mass_repartitioning}",
        f"atom_count: {atom_count}",
        f"residue_count: {residue_count}",
        f"system_xml_path: {_path_if_exists(system.system_xml_path) or '-'}",
        f"simulation_config_path: {_display_path(system.simulation_config_path)}",
        f"status: {status}",
    ]
    if failure_reason:
        lines.append(f"failure_reason: {failure_reason}")
    if traceback_text:
        lines.extend(("", "traceback:", traceback_text.rstrip()))
    atomic_write_text(system.system_build_log_path, "\n".join(lines) + "\n")


def _write_round2_equilibration_failure_log(log_path: Path, failure_reason: str, traceback_text: str) -> None:
    write_csv_rows(
        log_path,
        (
            {
                "stage_name": "failure",
                "status": "failure",
                "attempt_index": 0,
                "minimization_max_iterations": 0,
                "step_count": 0,
                "duration_ps": 0.0,
                "potential_energy_kj_per_mol": None,
                "temperature_K": None,
                "box_volume_nm3": None,
                "notes": failure_reason if not traceback_text else f"{failure_reason}\n{traceback_text.rstrip()}",
            },
        ),
    )


def _write_round2_screening_log(
    *,
    system: Round2OpenMMSystem,
    runtime_config: OpenMMScreeningRuntimeConfig,
    replicate_spec: ScreeningReplicateSpec,
    state_data_summary: StateDataSummary | None,
    artifacts: OpenMMScreeningArtifacts | None,
    log_path: Path,
    status: str,
    failure_reason: str,
    traceback_text: str = "",
) -> None:
    lines = [
        "Phase 7G1 round-2 OpenMM screening replicate",
        f"candidate_id: {system.candidate_id}",
        f"target_metal: {system.target_metal}",
        f"topology_class: {system.topology_class}",
        f"replicate_id: {replicate_spec.replicate_id}",
        f"seed: {replicate_spec.seed}",
        f"target_temperature_K: {runtime_config.target_temperature_K}",
        f"pressure_atm: {runtime_config.pressure_atm:.1f}",
        f"timestep_fs: {runtime_config.timestep_fs:.1f}",
        f"production_duration_ps: {runtime_config.production_duration_ps:.1f}",
        f"report_stride_steps: {runtime_config.report_stride_steps}",
        f"inner_sphere_cutoff_A: {runtime_config.inner_sphere_cutoff_A:.1f}",
        f"status: {status}",
    ]
    if state_data_summary is not None:
        lines.extend(
            [
                f"state_data_frame_count: {state_data_summary.frame_count}",
                f"mean_temperature_K: {_format_metric(state_data_summary.mean_temperature_K)}",
                f"final_temperature_K: {_format_metric(state_data_summary.final_temperature_K)}",
                f"mean_box_volume_nm3: {_format_metric(state_data_summary.mean_box_volume_nm3)}",
            ]
        )
    if artifacts is not None:
        lines.extend(
            [
                f"platform_name: {artifacts.platform_name}",
                f"platform_properties: {artifacts.platform_properties}",
                f"metal_atom_count: {artifacts.metal_atom_count}",
                f"oxygen_atom_count: {artifacts.oxygen_atom_count}",
                f"screened_frame_count: {artifacts.frame_count}",
                f"mean_min_metal_oxygen_distance_A: {_format_metric(artifacts.mean_min_metal_oxygen_distance_A)}",
                f"inner_sphere_occupancy_fraction: {_format_metric(artifacts.inner_sphere_occupancy_fraction)}",
            ]
        )
    if failure_reason:
        lines.append(f"failure_reason: {failure_reason}")
    if traceback_text:
        lines.extend(("", "traceback:", traceback_text.rstrip()))
    atomic_write_text(log_path, "\n".join(lines) + "\n")


def _build_replicate_summary_row(
    *,
    system: Round2OpenMMSystem,
    replicate_spec: ScreeningReplicateSpec,
    state_data_summary: StateDataSummary | None,
    artifacts: OpenMMScreeningArtifacts | None,
    failure_reason: str,
) -> Round2OpenMMScreeningReplicateSummaryRow:
    state_data_path, final_state_path, final_structure_path, screening_log_path, trajectory_path = _screening_artifact_paths(
        system,
        replicate_spec.replicate_id,
    )
    if artifacts is None or state_data_summary is None:
        return Round2OpenMMScreeningReplicateSummaryRow(
            candidate_id=system.candidate_id,
            target_metal=system.target_metal,
            topology_class=system.topology_class,
            replicate_id=replicate_spec.replicate_id,
            seed=replicate_spec.seed,
            success_status="failure",
            failure_reason=failure_reason,
            mean_temperature_K=None,
            final_temperature_K=None,
            mean_box_volume_nm3=None,
            mean_min_metal_oxygen_distance_A=None,
            inner_sphere_occupancy_fraction=None,
            state_data_path=_path_if_exists(state_data_path),
            final_state_path=_path_if_exists(final_state_path),
            final_structure_path=_path_if_exists(final_structure_path),
            screening_log_path=_path_if_exists(screening_log_path),
            trajectory_path=_path_if_exists(trajectory_path),
        )
    return Round2OpenMMScreeningReplicateSummaryRow(
        candidate_id=system.candidate_id,
        target_metal=system.target_metal,
        topology_class=system.topology_class,
        replicate_id=replicate_spec.replicate_id,
        seed=replicate_spec.seed,
        success_status="success",
        failure_reason="",
        mean_temperature_K=state_data_summary.mean_temperature_K,
        final_temperature_K=state_data_summary.final_temperature_K,
        mean_box_volume_nm3=state_data_summary.mean_box_volume_nm3,
        mean_min_metal_oxygen_distance_A=artifacts.mean_min_metal_oxygen_distance_A,
        inner_sphere_occupancy_fraction=artifacts.inner_sphere_occupancy_fraction,
        state_data_path=_display_path(state_data_path),
        final_state_path=_display_path(final_state_path),
        final_structure_path=_display_path(final_structure_path),
        screening_log_path=_display_path(screening_log_path),
        trajectory_path=_display_path(trajectory_path),
    )


def _parse_key_value_log(path: Path) -> dict[str, str]:
    payload: dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        if raw_line == "traceback:":
            break
        if ":" not in raw_line:
            continue
        key, value = raw_line.split(":", 1)
        payload[key.strip()] = value.strip()
    return payload


def load_existing_round2_successful_replicate(
    *,
    system: Round2OpenMMSystem,
    runtime_config: OpenMMScreeningRuntimeConfig,
    replicate_spec: ScreeningReplicateSpec,
) -> Round2OpenMMScreeningReplicateSummaryRow | None:
    """Load a previously completed successful replicate so reruns can resume deterministically."""
    state_data_path, final_state_path, final_structure_path, screening_log_path, trajectory_path = _screening_artifact_paths(
        system,
        replicate_spec.replicate_id,
    )
    required_paths = (
        state_data_path,
        final_state_path,
        final_structure_path,
        screening_log_path,
        trajectory_path,
    )
    if not all(path.exists() for path in required_paths):
        return None
    payload = _parse_key_value_log(screening_log_path)
    if payload.get("status") != "success":
        return None
    if payload.get("seed") != str(replicate_spec.seed):
        return None
    state_data_summary = summarize_state_data_csv(state_data_path)
    artifacts = OpenMMScreeningArtifacts(
        platform_name=payload.get("platform_name", ""),
        platform_properties=payload.get("platform_properties", ""),
        metal_atom_count=int(payload.get("metal_atom_count", "0") or 0),
        oxygen_atom_count=int(payload.get("oxygen_atom_count", "0") or 0),
        frame_count=int(payload.get("screened_frame_count", "0") or 0),
        mean_min_metal_oxygen_distance_A=_parse_optional_float(
            payload.get("mean_min_metal_oxygen_distance_A", "")
        )
        or 0.0,
        inner_sphere_occupancy_fraction=_parse_optional_float(
            payload.get("inner_sphere_occupancy_fraction", "")
        )
        or 0.0,
    )
    _write_round2_screening_log(
        system=system,
        runtime_config=runtime_config,
        replicate_spec=replicate_spec,
        state_data_summary=state_data_summary,
        artifacts=artifacts,
        log_path=screening_log_path,
        status="success",
        failure_reason="",
    )
    return _build_replicate_summary_row(
        system=system,
        replicate_spec=replicate_spec,
        state_data_summary=state_data_summary,
        artifacts=artifacts,
        failure_reason="",
    )


def _run_round2_screening_task_with_runner(
    task: Round2OpenMMScreeningReplicateTask,
    screening_runner: Callable[..., OpenMMScreeningArtifacts],
) -> Round2OpenMMScreeningReplicateSummaryRow:
    state_data_path, final_state_path, final_structure_path, screening_log_path, trajectory_path = _screening_artifact_paths(
        task.system,
        task.replicate_spec.replicate_id,
    )
    for artifact_path in (
        state_data_path,
        final_state_path,
        final_structure_path,
        screening_log_path,
        trajectory_path,
    ):
        artifact_path.unlink(missing_ok=True)

    artifacts: OpenMMScreeningArtifacts | None = None
    state_data_summary: StateDataSummary | None = None
    try:
        artifacts = screening_runner(
            prepared_structure_pdb_text=task.system.prepared_structure_path.read_text(encoding="utf-8"),
            system_xml=task.system.system_xml_path.read_text(encoding="utf-8"),
            equilibrated_state_xml=task.system.equilibration_final_state_path.read_text(encoding="utf-8"),
            target_metal=task.system.target_metal,
            target_temperature_K=task.runtime_config.target_temperature_K,
            friction_coeff_ps=task.runtime_config.friction_coeff_ps,
            timestep_fs=task.runtime_config.timestep_fs,
            cutoff_nm=task.runtime_config.cutoff_nm,
            production_duration_ps=task.runtime_config.production_duration_ps,
            pressure_atm=task.runtime_config.pressure_atm,
            report_stride_steps=task.runtime_config.report_stride_steps,
            inner_sphere_cutoff_A=task.runtime_config.inner_sphere_cutoff_A,
            seed=task.replicate_spec.seed,
            state_data_path=state_data_path,
            final_state_path=final_state_path,
            final_structure_path=final_structure_path,
            screening_log_path=screening_log_path,
            trajectory_path=trajectory_path,
            cpu_thread_count=task.cpu_thread_count,
        )
        state_data_summary = summarize_state_data_csv(state_data_path)
        _write_round2_screening_log(
            system=task.system,
            runtime_config=task.runtime_config,
            replicate_spec=task.replicate_spec,
            state_data_summary=state_data_summary,
            artifacts=artifacts,
            log_path=screening_log_path,
            status="success",
            failure_reason="",
        )
        return _build_replicate_summary_row(
            system=task.system,
            replicate_spec=task.replicate_spec,
            state_data_summary=state_data_summary,
            artifacts=artifacts,
            failure_reason="",
        )
    except Exception as exc:
        failure_reason = f"{type(exc).__name__}: {exc}"
        traceback_text = traceback.format_exc()
        if state_data_path.exists():
            try:
                state_data_summary = summarize_state_data_csv(state_data_path)
            except Exception:
                state_data_summary = None
        _write_round2_screening_log(
            system=task.system,
            runtime_config=task.runtime_config,
            replicate_spec=task.replicate_spec,
            state_data_summary=state_data_summary,
            artifacts=artifacts,
            log_path=screening_log_path,
            status="failure",
            failure_reason=failure_reason,
            traceback_text=traceback_text,
        )
        return _build_replicate_summary_row(
            system=task.system,
            replicate_spec=task.replicate_spec,
            state_data_summary=None,
            artifacts=None,
            failure_reason=failure_reason,
        )


def _run_default_round2_screening_task(
    task: Round2OpenMMScreeningReplicateTask,
) -> Round2OpenMMScreeningReplicateSummaryRow:
    return _run_round2_screening_task_with_runner(task, run_openmm_screening_replicate)


def run_round2_openmm_screening(
    *,
    panel_path: Path,
    md_protocol_path: Path,
    metal_parameter_values_path: Path,
    output_root: Path,
    summary_path: Path,
    panel_status_path: Path,
    report_path: Path,
    system_builder: Callable[..., OpenMMSerializedSystem] = build_openmm_serialized_system,
    equilibrator: Callable[..., OpenMMEquilibrationSmokeArtifacts] = run_openmm_equilibration_smoke_protocol,
    screening_runner: Callable[..., OpenMMScreeningArtifacts] = run_openmm_screening_replicate,
) -> tuple[Round2OpenMMSystemSummaryRow, ...]:
    """Execute the full Phase 7G1 reduced build/equilibration/screening cycle."""
    candidates = discover_round2_openmm_screening_candidates(panel_path)
    systems = build_round2_openmm_systems(candidates=candidates, output_root=output_root)
    prepare_round2_md_inputs(candidates=candidates, systems=systems)
    build_context = resolve_round2_build_context(
        md_protocol_path=md_protocol_path,
        metal_parameter_values_path=metal_parameter_values_path,
    )

    stage_failure_reasons: dict[tuple[str, str], str] = {}
    screenable_systems: list[Round2OpenMMSystem] = []

    for system in systems:
        simulation_config = build_round2_simulation_config(
            system=system,
            build_context=build_context,
            md_protocol_path=md_protocol_path,
        )
        write_yaml(system.simulation_config_path, asdict(simulation_config))

        build_artifacts_exist = all(
            path.exists()
            for path in (
                system.prepared_structure_path,
                system.system_xml_path,
                system.integrator_xml_path,
                system.simulation_config_path,
            )
        )
        if not build_artifacts_exist:
            for artifact_path in (
                system.prepared_structure_path,
                system.system_xml_path,
                system.integrator_xml_path,
                system.system_build_log_path,
            ):
                artifact_path.unlink(missing_ok=True)
            try:
                built_system = system_builder(
                    prepared_structure_pdb_text=system.starting_structure_path.read_text(encoding="utf-8"),
                    forcefield_files=(build_context.protein_forcefield, build_context.water_model),
                    target_temperature_K=simulation_config.target_temperature_K,
                    friction_coeff_ps=simulation_config.friction_coeff_ps,
                    timestep_fs=simulation_config.timestep_fs,
                    nonbonded_method=simulation_config.nonbonded_method,
                    cutoff_nm=simulation_config.cutoff_nm,
                    hydrogen_mass_repartitioning=simulation_config.hydrogen_mass_repartitioning,
                )
                atomic_write_text(system.prepared_structure_path, built_system.prepared_structure_pdb_text)
                atomic_write_text(system.system_xml_path, built_system.system_xml)
                atomic_write_text(system.integrator_xml_path, built_system.integrator_xml)
                _write_round2_build_log(
                    system=system,
                    simulation_config=simulation_config,
                    build_context=build_context,
                    atom_count=built_system.atom_count,
                    residue_count=built_system.residue_count,
                    status="success",
                    failure_reason="",
                )
            except Exception as exc:
                failure_reason = f"{type(exc).__name__}: {exc}"
                _write_round2_build_log(
                    system=system,
                    simulation_config=simulation_config,
                    build_context=build_context,
                    atom_count=0,
                    residue_count=0,
                    status="failure",
                    failure_reason=failure_reason,
                    traceback_text=traceback.format_exc(),
                )
                stage_failure_reasons[(system.candidate_id, system.target_metal)] = failure_reason
                continue

        equilibration_artifacts_exist = all(
            path.exists()
            for path in (
                system.minimized_structure_path,
                system.nvt_final_path,
                system.npt_final_path,
                system.equilibration_log_path,
                system.equilibration_final_state_path,
            )
        )
        if not equilibration_artifacts_exist:
            for artifact_path in (
                system.minimized_structure_path,
                system.nvt_final_path,
                system.npt_final_path,
                system.equilibration_log_path,
                system.equilibration_final_state_path,
            ):
                artifact_path.unlink(missing_ok=True)
            try:
                loaded_simulation_config = load_simulation_config(system.simulation_config_path)
                equilibration_runtime_config = build_equilibration_runtime_config(loaded_simulation_config)
                equilibration_artifacts = equilibrator(
                    prepared_structure_pdb_text=system.prepared_structure_path.read_text(encoding="utf-8"),
                    system_xml=system.system_xml_path.read_text(encoding="utf-8"),
                    target_temperature_K=equilibration_runtime_config.target_temperature_K,
                    friction_coeff_ps=equilibration_runtime_config.friction_coeff_ps,
                    timestep_fs=equilibration_runtime_config.timestep_fs,
                    cutoff_nm=equilibration_runtime_config.cutoff_nm,
                    nvt_duration_ps=equilibration_runtime_config.nvt_duration_ps,
                    npt_duration_ps=equilibration_runtime_config.npt_duration_ps,
                    pressure_atm=equilibration_runtime_config.pressure_atm,
                )
                atomic_write_text(system.minimized_structure_path, equilibration_artifacts.minimized_structure_pdb_text)
                atomic_write_text(system.nvt_final_path, equilibration_artifacts.nvt_final_pdb_text)
                atomic_write_text(system.npt_final_path, equilibration_artifacts.npt_final_pdb_text)
                atomic_write_text(system.equilibration_final_state_path, equilibration_artifacts.final_state_xml)
                write_csv_rows(system.equilibration_log_path, equilibration_artifacts.log_rows)
            except Exception as exc:
                failure_reason = f"{type(exc).__name__}: {exc}"
                _write_round2_equilibration_failure_log(
                    system.equilibration_log_path,
                    failure_reason,
                    traceback.format_exc(),
                )
                stage_failure_reasons[(system.candidate_id, system.target_metal)] = failure_reason
                continue

        screenable_systems.append(system)

    runtime_reference: OpenMMScreeningRuntimeConfig | None = None
    existing_replicate_rows: list[Round2OpenMMScreeningReplicateSummaryRow] = []
    pending_tasks: list[Round2OpenMMScreeningReplicateTask] = []
    for task in build_round2_openmm_replicate_tasks(tuple(screenable_systems)):
        runtime_reference = runtime_reference or task.runtime_config
        existing_row = load_existing_round2_successful_replicate(
            system=task.system,
            runtime_config=task.runtime_config,
            replicate_spec=task.replicate_spec,
        )
        if existing_row is not None:
            existing_replicate_rows.append(existing_row)
        else:
            pending_tasks.append(task)

    if screening_runner is run_openmm_screening_replicate and pending_tasks:
        with ProcessPoolExecutor(max_workers=DEFAULT_SCREENING_WORKER_COUNT) as executor:
            new_replicate_rows = list(executor.map(_run_default_round2_screening_task, pending_tasks))
    else:
        new_replicate_rows = [
            _run_round2_screening_task_with_runner(task, screening_runner)
            for task in pending_tasks
        ]
    replicate_rows = tuple(
        sorted(
            (*existing_replicate_rows, *new_replicate_rows),
            key=lambda row: (
                next(index for index, system in enumerate(systems) if (system.candidate_id, system.target_metal) == (row.candidate_id, row.target_metal)),
                row.replicate_id,
            ),
        )
    )

    if runtime_reference is None:
        first_system = systems[0]
        runtime_reference = build_openmm_screening_runtime_config(first_system.simulation_config_path)

    system_rows = build_round2_openmm_screening_system_summaries(
        systems=systems,
        replicate_rows=replicate_rows,
        stage_failure_reasons=stage_failure_reasons,
    )
    panel_rows = build_round2_openmm_panel_status_rows(system_rows)
    write_csv_rows(summary_path, system_rows)
    write_csv_rows(panel_status_path, panel_rows)
    atomic_write_text(
        report_path,
        render_round2_openmm_screening_report(
            candidates=candidates,
            replicate_rows=replicate_rows,
            system_rows=system_rows,
            panel_rows=panel_rows,
            runtime_config=runtime_reference,
        ),
    )
    return system_rows
