"""Phase 6B2 short multi-replicate OpenMM screening across the MD panel."""

from __future__ import annotations

import csv
import os
import traceback
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Callable

from lanm.analysis.openmm_equilibration_smoke import (
    load_openmm_system_build_summary,
    load_simulation_config,
)
from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_text, write_csv_rows
from lanm.md.openmm_screening import (
    DEFAULT_INNER_SPHERE_CUTOFF_A,
    DEFAULT_PRESSURE_ATM,
    DEFAULT_REPORT_STRIDE_STEPS,
    DEFAULT_SCREENING_DURATION_PS,
    OpenMMScreeningArtifacts,
    run_openmm_screening_replicate,
)
from lanm.paths import REPO_ROOT

SCREENING_REPLICATE_IDS = (1, 2, 3)
SCREENING_REPLICATE_SEEDS = (101, 102, 103)
DEFAULT_SCREENING_CPU_THREAD_COUNT = 2
DEFAULT_SCREENING_WORKER_COUNT = max(
    1,
    min(8, max(1, (os.cpu_count() or 1) // DEFAULT_SCREENING_CPU_THREAD_COUNT)),
)


@dataclass(frozen=True, slots=True)
class OpenMMEquilibrationSummaryRow:
    panel_member_id: str
    target_metal: str
    topology_class: str
    success_status: str
    failure_reason: str
    minimized_potential_energy: float | None
    final_nvt_potential_energy: float | None
    final_npt_potential_energy: float | None
    final_temperature_K: float | None
    final_box_volume_nm3: float | None
    final_structure_path: str


@dataclass(frozen=True, slots=True)
class ScreenableOpenMMSystem:
    panel_member_id: str
    target_metal: str
    topology_class: str
    system_xml_path: Path
    simulation_config_path: Path
    prepared_structure_path: Path
    equilibration_final_state_path: Path


@dataclass(frozen=True, slots=True)
class OpenMMScreeningRuntimeConfig:
    replicate_count: int
    target_temperature_K: int
    friction_coeff_ps: float
    timestep_fs: float
    cutoff_nm: float
    production_duration_ps: float
    pressure_atm: float
    report_stride_steps: int
    inner_sphere_cutoff_A: float


@dataclass(frozen=True, slots=True)
class ScreeningReplicateSpec:
    replicate_id: int
    seed: int


@dataclass(frozen=True, slots=True)
class StateDataSummary:
    frame_count: int
    mean_temperature_K: float | None
    final_temperature_K: float | None
    mean_box_volume_nm3: float | None


@dataclass(frozen=True, slots=True)
class OpenMMScreeningReplicateSummaryRow:
    panel_member_id: str
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
class OpenMMScreeningSystemSummaryRow:
    panel_member_id: str
    target_metal: str
    topology_class: str
    replicate_success_count: int
    median_mean_min_metal_oxygen_distance_A: float | None
    median_inner_sphere_occupancy_fraction: float | None
    screening_status: str


@dataclass(frozen=True, slots=True)
class OpenMMScreeningPanelStatusRow:
    panel_member_id: str
    topology_class: str
    dy_status: str
    nd_status: str
    y_status: str
    al_status: str
    fe_status: str
    stable_bound_rare_earth_count: int
    weak_binding_flag_count: int
    persistent_capture_flag_count: int
    no_persistent_capture_count: int
    failed_system_count: int


@dataclass(frozen=True, slots=True)
class ScreeningReplicateTask:
    system_row: ScreenableOpenMMSystem
    runtime_config: OpenMMScreeningRuntimeConfig
    replicate_spec: ScreeningReplicateSpec
    output_root: Path
    cpu_thread_count: int


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
    if lowered in {"nan", "none", "null"}:
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


def load_openmm_equilibration_smoke_summary(path: Path) -> tuple[OpenMMEquilibrationSummaryRow, ...]:
    """Load the Phase 6B1 equilibration summary table."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = tuple(
            OpenMMEquilibrationSummaryRow(
                panel_member_id=str(row["panel_member_id"]).strip(),
                target_metal=str(row["target_metal"]).strip(),
                topology_class=str(row["topology_class"]).strip(),
                success_status=str(row["success_status"]).strip(),
                failure_reason=str(row["failure_reason"]).strip(),
                minimized_potential_energy=_parse_optional_float(str(row["minimized_potential_energy"])),
                final_nvt_potential_energy=_parse_optional_float(str(row["final_nvt_potential_energy"])),
                final_npt_potential_energy=_parse_optional_float(str(row["final_npt_potential_energy"])),
                final_temperature_K=_parse_optional_float(str(row["final_temperature_K"])),
                final_box_volume_nm3=_parse_optional_float(str(row["final_box_volume_nm3"])),
                final_structure_path=str(row["final_structure_path"]).strip(),
            )
            for row in reader
        )
    if not rows:
        raise ValueError(f"No rows found in {path}")
    return rows


def discover_successful_openmm_screening_systems(
    *,
    system_build_summary_path: Path,
    equilibration_summary_path: Path,
) -> tuple[ScreenableOpenMMSystem, ...]:
    """Resolve the successful Phase 6B1 systems that should enter Phase 6B2."""
    build_rows = load_openmm_system_build_summary(system_build_summary_path)
    equilibration_rows = load_openmm_equilibration_smoke_summary(equilibration_summary_path)
    equilibration_by_key = {
        (row.panel_member_id, row.target_metal): row
        for row in equilibration_rows
    }

    selected_rows: list[ScreenableOpenMMSystem] = []
    for build_row in build_rows:
        if build_row.success_status != "success":
            continue
        key = (build_row.panel_member_id, build_row.target_metal)
        equilibration_row = equilibration_by_key.get(key)
        if equilibration_row is None:
            raise ValueError(f"Missing Phase 6B1 equilibration row for {key}")
        if equilibration_row.topology_class != build_row.topology_class:
            raise ValueError(
                f"Topology mismatch for {key}: Phase 6A4 {build_row.topology_class!r} vs "
                f"Phase 6B1 {equilibration_row.topology_class!r}"
            )
        if equilibration_row.success_status != "success":
            continue
        system_xml_path = _repo_path(build_row.system_xml_path)
        simulation_config_path = _repo_path(build_row.simulation_config_path)
        prepared_structure_path = system_xml_path.parent / "prepared_structure.pdb"
        if equilibration_row.final_structure_path:
            equilibration_dir = _repo_path(equilibration_row.final_structure_path).parent
        else:
            equilibration_dir = REPO_ROOT / "results" / "openmm_equilibration_smoke" / build_row.panel_member_id / build_row.target_metal
        equilibration_final_state_path = equilibration_dir / "final_state.xml"
        require_path(system_xml_path)
        require_path(simulation_config_path)
        require_path(prepared_structure_path)
        require_path(equilibration_final_state_path)
        selected_rows.append(
            ScreenableOpenMMSystem(
                panel_member_id=build_row.panel_member_id,
                target_metal=build_row.target_metal,
                topology_class=build_row.topology_class,
                system_xml_path=system_xml_path,
                simulation_config_path=simulation_config_path,
                prepared_structure_path=prepared_structure_path,
                equilibration_final_state_path=equilibration_final_state_path,
            )
        )

    if not selected_rows:
        raise ValueError(
            "No successful Phase 6B1 systems were available for Phase 6B2 screening"
        )
    return tuple(selected_rows)


def build_openmm_screening_runtime_config(
    simulation_config_path: Path,
) -> OpenMMScreeningRuntimeConfig:
    """Translate the per-system Phase 6A4 config into the fixed Phase 6B2 runtime config."""
    simulation_config = load_simulation_config(simulation_config_path)
    if simulation_config.replicate_count != len(SCREENING_REPLICATE_IDS):
        raise ValueError(
            "Phase 6B2 requires replicate_count=3, "
            f"found {simulation_config.replicate_count} in {simulation_config_path}"
        )
    if simulation_config.target_temperature_K != 298:
        raise ValueError(
            f"Phase 6B2 requires target_temperature_K=298, found {simulation_config.target_temperature_K}"
        )
    if simulation_config.timestep_fs != 2.0:
        raise ValueError(
            f"Phase 6B2 requires timestep_fs=2.0, found {simulation_config.timestep_fs}"
        )
    return OpenMMScreeningRuntimeConfig(
        replicate_count=simulation_config.replicate_count,
        target_temperature_K=simulation_config.target_temperature_K,
        friction_coeff_ps=simulation_config.friction_coeff_ps,
        timestep_fs=simulation_config.timestep_fs,
        cutoff_nm=simulation_config.cutoff_nm,
        production_duration_ps=DEFAULT_SCREENING_DURATION_PS,
        pressure_atm=DEFAULT_PRESSURE_ATM,
        report_stride_steps=DEFAULT_REPORT_STRIDE_STEPS,
        inner_sphere_cutoff_A=DEFAULT_INNER_SPHERE_CUTOFF_A,
    )


def build_openmm_screening_replicate_schedule(
    runtime_config: OpenMMScreeningRuntimeConfig,
) -> tuple[ScreeningReplicateSpec, ...]:
    """Construct the fixed deterministic replicate schedule for Phase 6B2."""
    if runtime_config.replicate_count != len(SCREENING_REPLICATE_IDS):
        raise ValueError(
            f"Phase 6B2 requires exactly {len(SCREENING_REPLICATE_IDS)} replicates, "
            f"found {runtime_config.replicate_count}"
        )
    return tuple(
        ScreeningReplicateSpec(replicate_id=replicate_id, seed=seed)
        for replicate_id, seed in zip(SCREENING_REPLICATE_IDS, SCREENING_REPLICATE_SEEDS, strict=True)
    )


def summarize_state_data_csv(path: Path) -> StateDataSummary:
    """Summarize the coarse screening `StateDataReporter` output."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        temperature_field = next(
            (field for field in fieldnames if "temperature" in field.lower()),
            None,
        )
        volume_field = next((field for field in fieldnames if "volume" in field.lower()), None)
        temperatures: list[float] = []
        volumes: list[float] = []
        for row in reader:
            if temperature_field is not None:
                temperature_value = _parse_optional_float(str(row.get(temperature_field, "")))
                if temperature_value is not None:
                    temperatures.append(temperature_value)
            if volume_field is not None:
                volume_value = _parse_optional_float(str(row.get(volume_field, "")))
                if volume_value is not None:
                    volumes.append(volume_value)
    return StateDataSummary(
        frame_count=len(temperatures) if temperatures else len(volumes),
        mean_temperature_K=(sum(temperatures) / len(temperatures)) if temperatures else None,
        final_temperature_K=temperatures[-1] if temperatures else None,
        mean_box_volume_nm3=(sum(volumes) / len(volumes)) if volumes else None,
    )


def build_screening_replicate_summary_row(
    *,
    system_row: ScreenableOpenMMSystem,
    replicate_spec: ScreeningReplicateSpec,
    output_dir: Path,
    state_data_summary: StateDataSummary | None,
    artifacts: OpenMMScreeningArtifacts | None,
    failure_reason: str,
) -> OpenMMScreeningReplicateSummaryRow:
    """Create one summary row for a Phase 6B2 replicate attempt."""
    state_data_path = output_dir / "state_data.csv"
    final_state_path = output_dir / "final_state.xml"
    final_structure_path = output_dir / "final_structure.pdb"
    screening_log_path = output_dir / "screening_log.txt"
    trajectory_path = output_dir / "trajectory.dcd"
    if artifacts is None or state_data_summary is None:
        return OpenMMScreeningReplicateSummaryRow(
            panel_member_id=system_row.panel_member_id,
            target_metal=system_row.target_metal,
            topology_class=system_row.topology_class,
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
    return OpenMMScreeningReplicateSummaryRow(
        panel_member_id=system_row.panel_member_id,
        target_metal=system_row.target_metal,
        topology_class=system_row.topology_class,
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


def build_openmm_screening_system_summaries(
    replicate_rows: tuple[OpenMMScreeningReplicateSummaryRow, ...],
) -> tuple[OpenMMScreeningSystemSummaryRow, ...]:
    """Aggregate per-replicate rows into one Phase 6B2 row per panel system."""
    grouped_rows: dict[tuple[str, str, str], list[OpenMMScreeningReplicateSummaryRow]] = defaultdict(list)
    for row in replicate_rows:
        grouped_rows[(row.panel_member_id, row.target_metal, row.topology_class)].append(row)

    system_rows: list[OpenMMScreeningSystemSummaryRow] = []
    for key in sorted(grouped_rows):
        panel_member_id, target_metal, topology_class = key
        rows = sorted(grouped_rows[key], key=lambda row: row.replicate_id)
        successful_rows = [row for row in rows if row.success_status == "success"]
        replicate_success_count = len(successful_rows)
        median_distance = _median_or_none(
            [row.mean_min_metal_oxygen_distance_A for row in successful_rows]
        )
        median_occupancy = _median_or_none(
            [row.inner_sphere_occupancy_fraction for row in successful_rows]
        )
        screening_status = classify_openmm_screening_status(
            target_metal=target_metal,
            replicate_success_count=replicate_success_count,
            median_inner_sphere_occupancy_fraction=median_occupancy,
        )
        system_rows.append(
            OpenMMScreeningSystemSummaryRow(
                panel_member_id=panel_member_id,
                target_metal=target_metal,
                topology_class=topology_class,
                replicate_success_count=replicate_success_count,
                median_mean_min_metal_oxygen_distance_A=median_distance,
                median_inner_sphere_occupancy_fraction=median_occupancy,
                screening_status=screening_status,
            )
        )
    return tuple(system_rows)


def classify_openmm_screening_status(
    *,
    target_metal: str,
    replicate_success_count: int,
    median_inner_sphere_occupancy_fraction: float | None,
) -> str:
    """Apply the fixed Phase 6B2 per-system screening classification rules."""
    if replicate_success_count < len(SCREENING_REPLICATE_IDS):
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
    raise ValueError(f"Unsupported target metal for Phase 6B2: {target_metal!r}")


def build_openmm_screening_panel_status_rows(
    system_rows: tuple[OpenMMScreeningSystemSummaryRow, ...],
) -> tuple[OpenMMScreeningPanelStatusRow, ...]:
    """Aggregate per-system statuses into one row per panel member."""
    rows_by_panel: dict[str, dict[str, OpenMMScreeningSystemSummaryRow]] = defaultdict(dict)
    topology_by_panel: dict[str, str] = {}
    for row in system_rows:
        rows_by_panel[row.panel_member_id][row.target_metal] = row
        topology_by_panel[row.panel_member_id] = row.topology_class

    panel_rows: list[OpenMMScreeningPanelStatusRow] = []
    for panel_member_id in sorted(rows_by_panel):
        metal_rows = rows_by_panel[panel_member_id]
        statuses = {metal: metal_rows[metal].screening_status for metal in metal_rows}
        status_counter = Counter(statuses.values())
        panel_rows.append(
            OpenMMScreeningPanelStatusRow(
                panel_member_id=panel_member_id,
                topology_class=topology_by_panel[panel_member_id],
                dy_status=statuses.get("Dy", ""),
                nd_status=statuses.get("Nd", ""),
                y_status=statuses.get("Y", ""),
                al_status=statuses.get("Al", ""),
                fe_status=statuses.get("Fe", ""),
                stable_bound_rare_earth_count=status_counter.get("stable_bound", 0),
                weak_binding_flag_count=status_counter.get("weak_binding_flag", 0),
                persistent_capture_flag_count=status_counter.get("persistent_capture_flag", 0),
                no_persistent_capture_count=status_counter.get("no_persistent_capture", 0),
                failed_system_count=status_counter.get("failed", 0),
            )
        )
    return tuple(panel_rows)


def render_openmm_screening_report(
    *,
    replicate_rows: tuple[OpenMMScreeningReplicateSummaryRow, ...],
    system_rows: tuple[OpenMMScreeningSystemSummaryRow, ...],
    panel_rows: tuple[OpenMMScreeningPanelStatusRow, ...],
    runtime_config: OpenMMScreeningRuntimeConfig,
) -> str:
    """Render the concise Phase 6B2 Markdown report."""
    successful_replicates = sum(row.success_status == "success" for row in replicate_rows)
    failed_replicates = len(replicate_rows) - successful_replicates
    status_counts = Counter(row.screening_status for row in system_rows)
    stride_ps = (runtime_config.report_stride_steps * runtime_config.timestep_fs) / 1000.0
    lines = [
        "# OpenMM Screening",
        "",
        "Phase 6B2 runs short deterministic 3-replicate NPT screening trajectories for every successful Phase 6B1 system to score stable rare-earth binding and persistent Al/Fe capture.",
        "",
        "## Screening Protocol",
        "",
        f"- screened systems: `{len(system_rows)}`",
        f"- replicate attempts: `{len(replicate_rows)}`",
        f"- successful replicates: `{successful_replicates}`",
        f"- failed replicates: `{failed_replicates}`",
        f"- replicate schedule: `1->101`, `2->102`, `3->103`",
        f"- target temperature: `{runtime_config.target_temperature_K} K`",
        f"- pressure target: `{runtime_config.pressure_atm:.1f} atm`",
        f"- timestep: `{runtime_config.timestep_fs:.1f} fs`",
        f"- production duration per replicate: `{runtime_config.production_duration_ps:.1f} ps`",
        f"- reporter stride: `{runtime_config.report_stride_steps}` steps (`{stride_ps:.1f} ps`)",
        f"- inner-sphere cutoff: `{runtime_config.inner_sphere_cutoff_A:.1f} A`",
        "",
        "## Status Counts",
        "",
    ]
    for status_name in (
        "stable_bound",
        "weak_binding_flag",
        "persistent_capture_flag",
        "no_persistent_capture",
        "failed",
    ):
        lines.append(f"- `{status_name}`: `{status_counts.get(status_name, 0)}`")
    lines.extend(
        [
            "",
            "## Per-System Summary",
            "",
            "| panel_member_id | target_metal | topology_class | replicate_success_count | median_mean_min_metal_oxygen_distance_A | median_inner_sphere_occupancy_fraction | screening_status |",
            "| --- | --- | --- | ---: | ---: | ---: | --- |",
        ]
    )
    for row in system_rows:
        lines.append(
            "| "
            f"{row.panel_member_id} | "
            f"{row.target_metal} | "
            f"{row.topology_class} | "
            f"{row.replicate_success_count} | "
            f"{_format_metric(row.median_mean_min_metal_oxygen_distance_A)} | "
            f"{_format_metric(row.median_inner_sphere_occupancy_fraction)} | "
            f"{row.screening_status} |"
        )
    lines.extend(
        [
            "",
            "## Panel Status",
            "",
            "| panel_member_id | topology_class | Dy | Nd | Y | Al | Fe |",
            "| --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in panel_rows:
        lines.append(
            "| "
            f"{row.panel_member_id} | "
            f"{row.topology_class} | "
            f"{row.dy_status or '-'} | "
            f"{row.nd_status or '-'} | "
            f"{row.y_status or '-'} | "
            f"{row.al_status or '-'} | "
            f"{row.fe_status or '-'} |"
        )
    failure_rows = [row for row in replicate_rows if row.success_status != "success"]
    if failure_rows:
        lines.extend(
            [
                "",
                "## Replicate Failures",
                "",
                "| panel_member_id | target_metal | replicate_id | failure_reason |",
                "| --- | --- | ---: | --- |",
            ]
        )
        for row in failure_rows:
            lines.append(
                "| "
                f"{row.panel_member_id} | "
                f"{row.target_metal} | "
                f"{row.replicate_id} | "
                f"{row.failure_reason} |"
            )
    lines.append("")
    return "\n".join(lines)


def _write_screening_log(
    *,
    log_path: Path,
    system_row: ScreenableOpenMMSystem,
    replicate_spec: ScreeningReplicateSpec,
    runtime_config: OpenMMScreeningRuntimeConfig,
    state_data_summary: StateDataSummary | None,
    artifacts: OpenMMScreeningArtifacts | None,
    status: str,
    failure_reason: str,
    traceback_text: str = "",
) -> None:
    lines = [
        "Phase 6B2 OpenMM screening replicate",
        f"panel_member_id: {system_row.panel_member_id}",
        f"target_metal: {system_row.target_metal}",
        f"topology_class: {system_row.topology_class}",
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


def _run_default_screening_replicate_task(
    task: ScreeningReplicateTask,
) -> OpenMMScreeningReplicateSummaryRow:
    output_dir = (
        task.output_root
        / task.system_row.panel_member_id
        / task.system_row.target_metal
        / f"replicate_{task.replicate_spec.replicate_id}"
    )
    state_data_path = output_dir / "state_data.csv"
    final_state_path = output_dir / "final_state.xml"
    final_structure_path = output_dir / "final_structure.pdb"
    screening_log_path = output_dir / "screening_log.txt"
    trajectory_path = output_dir / "trajectory.dcd"
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
        artifacts = run_openmm_screening_replicate(
            prepared_structure_pdb_text=task.system_row.prepared_structure_path.read_text(encoding="utf-8"),
            system_xml=task.system_row.system_xml_path.read_text(encoding="utf-8"),
            equilibrated_state_xml=task.system_row.equilibration_final_state_path.read_text(encoding="utf-8"),
            target_metal=task.system_row.target_metal,
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
        _write_screening_log(
            log_path=screening_log_path,
            system_row=task.system_row,
            replicate_spec=task.replicate_spec,
            runtime_config=task.runtime_config,
            state_data_summary=state_data_summary,
            artifacts=artifacts,
            status="success",
            failure_reason="",
        )
        return build_screening_replicate_summary_row(
            system_row=task.system_row,
            replicate_spec=task.replicate_spec,
            output_dir=output_dir,
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
        _write_screening_log(
            log_path=screening_log_path,
            system_row=task.system_row,
            replicate_spec=task.replicate_spec,
            runtime_config=task.runtime_config,
            state_data_summary=state_data_summary,
            artifacts=artifacts,
            status="failure",
            failure_reason=failure_reason,
            traceback_text=traceback_text,
        )
        return build_screening_replicate_summary_row(
            system_row=task.system_row,
            replicate_spec=task.replicate_spec,
            output_dir=output_dir,
            state_data_summary=None,
            artifacts=None,
            failure_reason=failure_reason,
        )


def run_openmm_screening(
    *,
    system_build_summary_path: Path,
    equilibration_summary_path: Path,
    output_root: Path,
    replicate_summary_path: Path,
    system_summary_path: Path,
    panel_status_path: Path,
    report_path: Path,
    screening_runner: Callable[..., OpenMMScreeningArtifacts] = run_openmm_screening_replicate,
) -> tuple[OpenMMScreeningReplicateSummaryRow, ...]:
    """Run Phase 6B2 on every successful Phase 6B1 system."""
    systems = discover_successful_openmm_screening_systems(
        system_build_summary_path=system_build_summary_path,
        equilibration_summary_path=equilibration_summary_path,
    )
    runtime_reference: OpenMMScreeningRuntimeConfig | None = None
    tasks: list[ScreeningReplicateTask] = []

    for system_row in systems:
        runtime_config = build_openmm_screening_runtime_config(system_row.simulation_config_path)
        runtime_reference = runtime_reference or runtime_config
        for replicate_spec in build_openmm_screening_replicate_schedule(runtime_config):
            tasks.append(
                ScreeningReplicateTask(
                    system_row=system_row,
                    runtime_config=runtime_config,
                    replicate_spec=replicate_spec,
                    output_root=output_root,
                    cpu_thread_count=DEFAULT_SCREENING_CPU_THREAD_COUNT,
                )
            )

    if screening_runner is run_openmm_screening_replicate:
        with ProcessPoolExecutor(max_workers=DEFAULT_SCREENING_WORKER_COUNT) as executor:
            replicate_rows = tuple(executor.map(_run_default_screening_replicate_task, tasks))
    else:
        results: list[OpenMMScreeningReplicateSummaryRow] = []
        for task in tasks:
            output_dir = (
                task.output_root
                / task.system_row.panel_member_id
                / task.system_row.target_metal
                / f"replicate_{task.replicate_spec.replicate_id}"
            )
            state_data_path = output_dir / "state_data.csv"
            final_state_path = output_dir / "final_state.xml"
            final_structure_path = output_dir / "final_structure.pdb"
            screening_log_path = output_dir / "screening_log.txt"
            trajectory_path = output_dir / "trajectory.dcd"
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
                    prepared_structure_pdb_text=task.system_row.prepared_structure_path.read_text(encoding="utf-8"),
                    system_xml=task.system_row.system_xml_path.read_text(encoding="utf-8"),
                    equilibrated_state_xml=task.system_row.equilibration_final_state_path.read_text(encoding="utf-8"),
                    target_metal=task.system_row.target_metal,
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
                )
                state_data_summary = summarize_state_data_csv(state_data_path)
                _write_screening_log(
                    log_path=screening_log_path,
                    system_row=task.system_row,
                    replicate_spec=task.replicate_spec,
                    runtime_config=task.runtime_config,
                    state_data_summary=state_data_summary,
                    artifacts=artifacts,
                    status="success",
                    failure_reason="",
                )
                results.append(
                    build_screening_replicate_summary_row(
                        system_row=task.system_row,
                        replicate_spec=task.replicate_spec,
                        output_dir=output_dir,
                        state_data_summary=state_data_summary,
                        artifacts=artifacts,
                        failure_reason="",
                    )
                )
            except Exception as exc:
                failure_reason = f"{type(exc).__name__}: {exc}"
                traceback_text = traceback.format_exc()
                if state_data_path.exists():
                    try:
                        state_data_summary = summarize_state_data_csv(state_data_path)
                    except Exception:
                        state_data_summary = None
                _write_screening_log(
                    log_path=screening_log_path,
                    system_row=task.system_row,
                    replicate_spec=task.replicate_spec,
                    runtime_config=task.runtime_config,
                    state_data_summary=state_data_summary,
                    artifacts=artifacts,
                    status="failure",
                    failure_reason=failure_reason,
                    traceback_text=traceback_text,
                )
                results.append(
                    build_screening_replicate_summary_row(
                        system_row=task.system_row,
                        replicate_spec=task.replicate_spec,
                        output_dir=output_dir,
                        state_data_summary=None,
                        artifacts=None,
                        failure_reason=failure_reason,
                    )
                )
        replicate_rows = tuple(results)
    system_rows = build_openmm_screening_system_summaries(replicate_rows)
    panel_rows = build_openmm_screening_panel_status_rows(system_rows)
    write_csv_rows(replicate_summary_path, replicate_rows)
    write_csv_rows(system_summary_path, system_rows)
    write_csv_rows(panel_status_path, panel_rows)
    if runtime_reference is None:  # pragma: no cover - defensive
        runtime_reference = OpenMMScreeningRuntimeConfig(
            replicate_count=3,
            target_temperature_K=298,
            friction_coeff_ps=1.0,
            timestep_fs=2.0,
            cutoff_nm=1.0,
            production_duration_ps=DEFAULT_SCREENING_DURATION_PS,
            pressure_atm=DEFAULT_PRESSURE_ATM,
            report_stride_steps=DEFAULT_REPORT_STRIDE_STEPS,
            inner_sphere_cutoff_A=DEFAULT_INNER_SPHERE_CUTOFF_A,
        )
    atomic_write_text(
        report_path,
        render_openmm_screening_report(
            replicate_rows=replicate_rows,
            system_rows=system_rows,
            panel_rows=panel_rows,
            runtime_config=runtime_reference,
        ),
    )
    failure_count = sum(row.success_status != "success" for row in replicate_rows)
    if failure_count:
        raise RuntimeError(
            f"OpenMM screening failed for {failure_count}/{len(replicate_rows)} replicates; "
            f"see {_display_path(replicate_summary_path)}"
        )
    return replicate_rows
