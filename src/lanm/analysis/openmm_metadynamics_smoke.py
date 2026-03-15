"""Phase 6C1 reduced-panel well-tempered metadynamics triage."""

from __future__ import annotations

import csv
import traceback
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from lanm.analysis.openmm_equilibration_smoke import (
    load_openmm_system_build_summary,
    load_simulation_config,
)
from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_text, write_csv_rows
from lanm.md.openmm_metadynamics_smoke import (
    DEFAULT_COORDINATION_NUMBER_BIAS_WIDTH,
    DEFAULT_COORDINATION_NUMBER_MAX,
    DEFAULT_COORDINATION_NUMBER_MIN,
    DEFAULT_COORDINATION_SWITCH_DISTANCE_A,
    DEFAULT_COORDINATION_SWITCH_POWER,
    DEFAULT_DONOR_CUTOFF_A,
    DEFAULT_ESCAPE_COORDINATION_THRESHOLD,
    DEFAULT_ESCAPE_DISTANCE_THRESHOLD_A,
    DEFAULT_MEAN_DISTANCE_BIAS_WIDTH_A,
    DEFAULT_MEAN_DISTANCE_MAX_A,
    DEFAULT_MEAN_DISTANCE_MIN_A,
    DEFAULT_METADYNAMICS_BIAS_FACTOR,
    DEFAULT_METADYNAMICS_CPU_THREAD_COUNT,
    DEFAULT_METADYNAMICS_DEPOSITION_FREQUENCY_STEPS,
    DEFAULT_METADYNAMICS_DURATION_PS,
    DEFAULT_METADYNAMICS_GAUSSIAN_HEIGHT_KJ_PER_MOL,
    DEFAULT_METADYNAMICS_REPORT_STRIDE_STEPS,
    OpenMMMetadynamicsArtifacts,
    run_openmm_metadynamics_smoke_system,
)
from lanm.paths import REPO_ROOT

REDUCED_METADYNAMICS_PANEL_MEMBER_IDS = (
    "am1_mex_ss_only_u02",
    "am1_mex_wt_reference",
    "hans_pocket_ss_only_u02",
    "hans_interface_ss_plus_if_u04",
    "hans_interface_wt_reference",
)
TARGET_METALS = ("Dy", "Nd", "Y", "Al", "Fe")
EXPECTED_METADYNAMICS_SYSTEM_COUNT = len(REDUCED_METADYNAMICS_PANEL_MEMBER_IDS) * len(TARGET_METALS)
DEFAULT_METADYNAMICS_WORKER_COUNT = 1


@dataclass(frozen=True, slots=True)
class OpenMMScreeningPanelStatusRow:
    panel_member_id: str
    topology_class: str
    dy_status: str
    nd_status: str
    y_status: str
    al_status: str
    fe_status: str


@dataclass(frozen=True, slots=True)
class ReducedMetadynamicsOpenMMSystem:
    panel_member_id: str
    target_metal: str
    topology_class: str
    system_xml_path: Path
    simulation_config_path: Path
    equilibrated_final_state_path: Path
    equilibrated_structure_path: Path
    seed: int


@dataclass(frozen=True, slots=True)
class OpenMMMetadynamicsRuntimeConfig:
    target_temperature_K: int
    friction_coeff_ps: float
    timestep_fs: float
    cutoff_nm: float
    duration_ps: float
    report_stride_steps: int
    bias_factor: float
    gaussian_height_kj_per_mol: float
    deposition_frequency_steps: int
    donor_cutoff_A: float
    coordination_switch_distance_A: float
    coordination_switch_power: int
    coordination_number_min: float
    coordination_number_max: float
    coordination_number_bias_width: float
    mean_distance_min_A: float
    mean_distance_max_A: float
    mean_distance_bias_width_A: float
    escape_coordination_threshold: float
    escape_distance_threshold_A: float


@dataclass(frozen=True, slots=True)
class OpenMMMetadynamicsSummaryRow:
    panel_member_id: str
    target_metal: str
    success_status: str
    failure_reason: str
    minimum_coordination_number: float | None
    maximum_mean_metal_oxygen_distance_A: float | None
    escape_event_detected: str
    final_coordination_number: float | None
    final_mean_metal_oxygen_distance_A: float | None
    metadynamics_status: str


@dataclass(frozen=True, slots=True)
class OpenMMMetadynamicsPanelStatusRow:
    panel_member_id: str
    dy_status: str
    nd_status: str
    y_status: str
    al_status: str
    fe_status: str
    candidate_keep_for_qm: str


@dataclass(frozen=True, slots=True)
class MetadynamicsTask:
    system_row: ReducedMetadynamicsOpenMMSystem
    runtime_config: OpenMMMetadynamicsRuntimeConfig
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


def _format_metric(value: float | None) -> str:
    if value is None:
        return "-"
    return f"{value:.6f}"


def _format_system_label(system_row: ReducedMetadynamicsOpenMMSystem) -> str:
    return f"{system_row.panel_member_id}/{system_row.target_metal}"


def build_metadynamics_seed(panel_member_id: str, target_metal: str) -> int:
    panel_index = REDUCED_METADYNAMICS_PANEL_MEMBER_IDS.index(panel_member_id) + 1
    metal_index = TARGET_METALS.index(target_metal) + 1
    return 61000 + (panel_index * 100) + metal_index


def load_openmm_screening_panel_status(path: Path) -> tuple[OpenMMScreeningPanelStatusRow, ...]:
    """Load the Phase 6B2 panel-status table."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = tuple(
            OpenMMScreeningPanelStatusRow(
                panel_member_id=str(row["panel_member_id"]).strip(),
                topology_class=str(row["topology_class"]).strip(),
                dy_status=str(row["dy_status"]).strip(),
                nd_status=str(row["nd_status"]).strip(),
                y_status=str(row["y_status"]).strip(),
                al_status=str(row["al_status"]).strip(),
                fe_status=str(row["fe_status"]).strip(),
            )
            for row in reader
        )
    if not rows:
        raise ValueError(f"No rows found in {path}")
    return rows


def discover_reduced_openmm_metadynamics_systems(
    *,
    screening_panel_status_path: Path,
    system_build_summary_path: Path,
) -> tuple[ReducedMetadynamicsOpenMMSystem, ...]:
    """Resolve the exact 25-system reduced Phase 6C1 panel."""
    screening_rows = load_openmm_screening_panel_status(screening_panel_status_path)
    build_rows = load_openmm_system_build_summary(system_build_summary_path)
    screening_by_panel = {row.panel_member_id: row for row in screening_rows}
    build_by_key = {
        (row.panel_member_id, row.target_metal): row
        for row in build_rows
        if row.success_status == "success"
    }

    selected_rows: list[ReducedMetadynamicsOpenMMSystem] = []
    for panel_member_id in REDUCED_METADYNAMICS_PANEL_MEMBER_IDS:
        screening_row = screening_by_panel.get(panel_member_id)
        if screening_row is None:
            raise ValueError(f"Missing Phase 6B2 panel-status row for {panel_member_id}")
        for target_metal in TARGET_METALS:
            key = (panel_member_id, target_metal)
            build_row = build_by_key.get(key)
            if build_row is None:
                raise ValueError(f"Missing successful Phase 6A4 build row for {key}")
            if build_row.topology_class != screening_row.topology_class:
                raise ValueError(
                    f"Topology mismatch for {key}: Phase 6A4 {build_row.topology_class!r} vs "
                    f"Phase 6B2 {screening_row.topology_class!r}"
                )
            system_xml_path = _repo_path(build_row.system_xml_path)
            simulation_config_path = _repo_path(build_row.simulation_config_path)
            equilibration_dir = REPO_ROOT / "results" / "openmm_equilibration_smoke" / panel_member_id / target_metal
            equilibrated_final_state_path = equilibration_dir / "final_state.xml"
            equilibrated_structure_path = equilibration_dir / "npt_final.pdb"
            require_path(system_xml_path)
            require_path(simulation_config_path)
            require_path(equilibrated_final_state_path)
            require_path(equilibrated_structure_path)
            selected_rows.append(
                ReducedMetadynamicsOpenMMSystem(
                    panel_member_id=panel_member_id,
                    target_metal=target_metal,
                    topology_class=build_row.topology_class,
                    system_xml_path=system_xml_path,
                    simulation_config_path=simulation_config_path,
                    equilibrated_final_state_path=equilibrated_final_state_path,
                    equilibrated_structure_path=equilibrated_structure_path,
                    seed=build_metadynamics_seed(panel_member_id, target_metal),
                )
            )

    if len(selected_rows) != EXPECTED_METADYNAMICS_SYSTEM_COUNT:
        raise ValueError(
            "Phase 6C1 must resolve exactly "
            f"{EXPECTED_METADYNAMICS_SYSTEM_COUNT} systems, found {len(selected_rows)}"
        )
    return tuple(selected_rows)


def build_openmm_metadynamics_runtime_config(
    simulation_config_path: Path,
) -> OpenMMMetadynamicsRuntimeConfig:
    """Translate the Phase 6A4 per-system config into the fixed Phase 6C1 runtime config."""
    simulation_config = load_simulation_config(simulation_config_path)
    if simulation_config.target_temperature_K != 298:
        raise ValueError(
            "Phase 6C1 requires target_temperature_K=298, "
            f"found {simulation_config.target_temperature_K}"
        )
    if simulation_config.timestep_fs != 2.0:
        raise ValueError(
            f"Phase 6C1 requires timestep_fs=2.0, found {simulation_config.timestep_fs}"
        )
    return OpenMMMetadynamicsRuntimeConfig(
        target_temperature_K=simulation_config.target_temperature_K,
        friction_coeff_ps=simulation_config.friction_coeff_ps,
        timestep_fs=simulation_config.timestep_fs,
        cutoff_nm=simulation_config.cutoff_nm,
        duration_ps=DEFAULT_METADYNAMICS_DURATION_PS,
        report_stride_steps=DEFAULT_METADYNAMICS_REPORT_STRIDE_STEPS,
        bias_factor=DEFAULT_METADYNAMICS_BIAS_FACTOR,
        gaussian_height_kj_per_mol=DEFAULT_METADYNAMICS_GAUSSIAN_HEIGHT_KJ_PER_MOL,
        deposition_frequency_steps=DEFAULT_METADYNAMICS_DEPOSITION_FREQUENCY_STEPS,
        donor_cutoff_A=DEFAULT_DONOR_CUTOFF_A,
        coordination_switch_distance_A=DEFAULT_COORDINATION_SWITCH_DISTANCE_A,
        coordination_switch_power=DEFAULT_COORDINATION_SWITCH_POWER,
        coordination_number_min=DEFAULT_COORDINATION_NUMBER_MIN,
        coordination_number_max=DEFAULT_COORDINATION_NUMBER_MAX,
        coordination_number_bias_width=DEFAULT_COORDINATION_NUMBER_BIAS_WIDTH,
        mean_distance_min_A=DEFAULT_MEAN_DISTANCE_MIN_A,
        mean_distance_max_A=DEFAULT_MEAN_DISTANCE_MAX_A,
        mean_distance_bias_width_A=DEFAULT_MEAN_DISTANCE_BIAS_WIDTH_A,
        escape_coordination_threshold=DEFAULT_ESCAPE_COORDINATION_THRESHOLD,
        escape_distance_threshold_A=DEFAULT_ESCAPE_DISTANCE_THRESHOLD_A,
    )


def classify_metadynamics_status(
    *,
    target_metal: str,
    success_status: str,
    escape_event_detected: str,
) -> str:
    """Apply the fixed Phase 6C1 status rules."""
    if success_status != "success":
        return "failed"
    escaped = escape_event_detected == "TRUE"
    if target_metal in {"Dy", "Nd", "Y"}:
        return "escape_or_weak_binding" if escaped else "retained_bound"
    if target_metal in {"Al", "Fe"}:
        return "escaped_under_bias" if escaped else "persistent_capture"
    raise ValueError(f"Unsupported target metal for Phase 6C1: {target_metal!r}")


def summarize_existing_cv_timeseries(
    path: Path,
) -> tuple[float, float, str, float, float]:
    """Summarize a previously written Phase 6C1 CV time series for resume support."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if not rows:
        raise ValueError(f"No CV time-series rows found in {path}")
    coordination_values = [float(str(row["coordination_number"]).strip()) for row in rows]
    mean_distance_values = [
        float(str(row["mean_metal_oxygen_distance_A"]).strip())
        for row in rows
    ]
    escape_event_detected = "TRUE" if any(str(row["escape_condition_met"]).strip() == "TRUE" for row in rows) else "FALSE"
    return (
        min(coordination_values),
        max(mean_distance_values),
        escape_event_detected,
        coordination_values[-1],
        mean_distance_values[-1],
    )


def build_openmm_metadynamics_summary_row(
    *,
    system_row: ReducedMetadynamicsOpenMMSystem,
    artifacts: OpenMMMetadynamicsArtifacts | None,
    failure_reason: str,
) -> OpenMMMetadynamicsSummaryRow:
    """Create one Phase 6C1 per-system summary row."""
    if artifacts is None:
        return OpenMMMetadynamicsSummaryRow(
            panel_member_id=system_row.panel_member_id,
            target_metal=system_row.target_metal,
            success_status="failure",
            failure_reason=failure_reason,
            minimum_coordination_number=None,
            maximum_mean_metal_oxygen_distance_A=None,
            escape_event_detected="FALSE",
            final_coordination_number=None,
            final_mean_metal_oxygen_distance_A=None,
            metadynamics_status="failed",
        )
    escape_event_detected = "TRUE" if artifacts.escape_event_detected else "FALSE"
    return OpenMMMetadynamicsSummaryRow(
        panel_member_id=system_row.panel_member_id,
        target_metal=system_row.target_metal,
        success_status="success",
        failure_reason="",
        minimum_coordination_number=artifacts.minimum_coordination_number,
        maximum_mean_metal_oxygen_distance_A=artifacts.maximum_mean_metal_oxygen_distance_A,
        escape_event_detected=escape_event_detected,
        final_coordination_number=artifacts.final_coordination_number,
        final_mean_metal_oxygen_distance_A=artifacts.final_mean_metal_oxygen_distance_A,
        metadynamics_status=classify_metadynamics_status(
            target_metal=system_row.target_metal,
            success_status="success",
            escape_event_detected=escape_event_detected,
        ),
    )


def _build_summary_row_from_existing_outputs(
    *,
    system_row: ReducedMetadynamicsOpenMMSystem,
    cv_timeseries_path: Path,
    bias_state_path: Path,
    final_state_path: Path,
    final_structure_path: Path,
    log_path: Path,
) -> OpenMMMetadynamicsSummaryRow | None:
    required_paths = (
        cv_timeseries_path,
        bias_state_path,
        final_state_path,
        final_structure_path,
        log_path,
    )
    if not all(path.exists() for path in required_paths):
        return None
    minimum_coordination_number, maximum_mean_distance_A, escape_event_detected, final_coordination_number, final_mean_distance_A = summarize_existing_cv_timeseries(
        cv_timeseries_path
    )
    return OpenMMMetadynamicsSummaryRow(
        panel_member_id=system_row.panel_member_id,
        target_metal=system_row.target_metal,
        success_status="success",
        failure_reason="",
        minimum_coordination_number=minimum_coordination_number,
        maximum_mean_metal_oxygen_distance_A=maximum_mean_distance_A,
        escape_event_detected=escape_event_detected,
        final_coordination_number=final_coordination_number,
        final_mean_metal_oxygen_distance_A=final_mean_distance_A,
        metadynamics_status=classify_metadynamics_status(
            target_metal=system_row.target_metal,
            success_status="success",
            escape_event_detected=escape_event_detected,
        ),
    )


def _is_obviously_stronger_than_reference(
    *,
    candidate_row: OpenMMMetadynamicsSummaryRow,
    reference_row: OpenMMMetadynamicsSummaryRow,
) -> bool:
    if (
        candidate_row.minimum_coordination_number is None
        or candidate_row.maximum_mean_metal_oxygen_distance_A is None
        or reference_row.minimum_coordination_number is None
        or reference_row.maximum_mean_metal_oxygen_distance_A is None
    ):
        return False
    return (
        candidate_row.minimum_coordination_number > reference_row.minimum_coordination_number + 0.5
        and candidate_row.maximum_mean_metal_oxygen_distance_A
        < reference_row.maximum_mean_metal_oxygen_distance_A - 0.1
    )


def classify_candidate_keep_for_qm(
    metal_rows: dict[str, OpenMMMetadynamicsSummaryRow],
) -> str:
    """Apply the Phase 6C1 panel-level keep gate."""
    required_rows = [metal_rows.get(metal) for metal in TARGET_METALS]
    if any(row is None for row in required_rows):
        return "FALSE"
    dy_row = metal_rows["Dy"]
    nd_row = metal_rows["Nd"]
    y_row = metal_rows["Y"]
    al_row = metal_rows["Al"]
    fe_row = metal_rows["Fe"]
    if any(row.success_status != "success" for row in (dy_row, nd_row, y_row, al_row, fe_row)):
        return "FALSE"
    if dy_row.metadynamics_status != "retained_bound":
        return "FALSE"
    if al_row.metadynamics_status != "escaped_under_bias":
        return "FALSE"
    if fe_row.metadynamics_status != "escaped_under_bias":
        return "FALSE"
    for competitor_row in (nd_row, y_row):
        if competitor_row.metadynamics_status == "failed":
            return "FALSE"
        if competitor_row.metadynamics_status == "retained_bound" and _is_obviously_stronger_than_reference(
            candidate_row=competitor_row,
            reference_row=dy_row,
        ):
            return "FALSE"
    return "TRUE"


def build_openmm_metadynamics_panel_status_rows(
    system_rows: tuple[OpenMMMetadynamicsSummaryRow, ...],
) -> tuple[OpenMMMetadynamicsPanelStatusRow, ...]:
    """Aggregate per-system Phase 6C1 results into one row per panel member."""
    rows_by_panel: dict[str, dict[str, OpenMMMetadynamicsSummaryRow]] = defaultdict(dict)
    for row in system_rows:
        rows_by_panel[row.panel_member_id][row.target_metal] = row

    def fallback_summary_row(target_metal: str) -> OpenMMMetadynamicsSummaryRow:
        return OpenMMMetadynamicsSummaryRow(
            panel_member_id="",
            target_metal=target_metal,
            success_status="failure",
            failure_reason="",
            minimum_coordination_number=None,
            maximum_mean_metal_oxygen_distance_A=None,
            escape_event_detected="FALSE",
            final_coordination_number=None,
            final_mean_metal_oxygen_distance_A=None,
            metadynamics_status="failed",
        )

    ordered_panel_ids = [
        panel_member_id
        for panel_member_id in REDUCED_METADYNAMICS_PANEL_MEMBER_IDS
        if panel_member_id in rows_by_panel
    ]
    ordered_panel_ids.extend(
        sorted(panel_member_id for panel_member_id in rows_by_panel if panel_member_id not in ordered_panel_ids)
    )

    panel_rows: list[OpenMMMetadynamicsPanelStatusRow] = []
    for panel_member_id in ordered_panel_ids:
        metal_rows = rows_by_panel.get(panel_member_id, {})
        panel_rows.append(
            OpenMMMetadynamicsPanelStatusRow(
                panel_member_id=panel_member_id,
                dy_status=metal_rows.get("Dy", fallback_summary_row("Dy")).metadynamics_status,
                nd_status=metal_rows.get("Nd", fallback_summary_row("Nd")).metadynamics_status,
                y_status=metal_rows.get("Y", fallback_summary_row("Y")).metadynamics_status,
                al_status=metal_rows.get("Al", fallback_summary_row("Al")).metadynamics_status,
                fe_status=metal_rows.get("Fe", fallback_summary_row("Fe")).metadynamics_status,
                candidate_keep_for_qm=classify_candidate_keep_for_qm(metal_rows),
            )
        )
    return tuple(panel_rows)


def render_openmm_metadynamics_smoke_report(
    *,
    system_rows: tuple[OpenMMMetadynamicsSummaryRow, ...],
    panel_rows: tuple[OpenMMMetadynamicsPanelStatusRow, ...],
    runtime_config: OpenMMMetadynamicsRuntimeConfig,
) -> str:
    """Render the concise Phase 6C1 Markdown report."""
    status_counts = Counter(row.metadynamics_status for row in system_rows)
    success_count = sum(row.success_status == "success" for row in system_rows)
    failure_count = len(system_rows) - success_count
    stride_ps = (runtime_config.report_stride_steps * runtime_config.timestep_fs) / 1000.0
    deposition_ps = (runtime_config.deposition_frequency_steps * runtime_config.timestep_fs) / 1000.0
    lines = [
        "# OpenMM Metadynamics Smoke",
        "",
        "Phase 6C1 runs one short deterministic well-tempered metadynamics triage trajectory per reduced-panel system to test whether local pocket binding is retained under bias.",
        "",
        "## Protocol",
        "",
        f"- reduced-panel systems attempted: `{len(system_rows)}`",
        f"- successful systems: `{success_count}`",
        f"- failed systems: `{failure_count}`",
        f"- panel members: `{', '.join(REDUCED_METADYNAMICS_PANEL_MEMBER_IDS)}`",
        f"- target metals: `{', '.join(TARGET_METALS)}`",
        f"- duration per system: `{runtime_config.duration_ps:.1f} ps`",
        f"- temperature: `{runtime_config.target_temperature_K} K`",
        f"- timestep: `{runtime_config.timestep_fs:.1f} fs`",
        f"- reporter stride: `{runtime_config.report_stride_steps}` steps (`{stride_ps:.1f} ps`)",
        f"- Gaussian deposition stride: `{runtime_config.deposition_frequency_steps}` steps (`{deposition_ps:.1f} ps`)",
        f"- bias factor: `{runtime_config.bias_factor:.1f}`",
        f"- Gaussian height: `{runtime_config.gaussian_height_kj_per_mol:.2f} kJ/mol`",
        f"- local-pocket donor cutoff: `{runtime_config.donor_cutoff_A:.1f} A`",
        f"- escape heuristic: `coordination <= {runtime_config.escape_coordination_threshold:.1f}` and `mean distance >= {runtime_config.escape_distance_threshold_A:.1f} A`",
        "",
        "## Status Counts",
        "",
    ]
    for status_name in (
        "retained_bound",
        "escape_or_weak_binding",
        "persistent_capture",
        "escaped_under_bias",
        "failed",
    ):
        lines.append(f"- `{status_name}`: `{status_counts.get(status_name, 0)}`")
    lines.extend(
        [
            "",
            "## Per-System Summary",
            "",
            "| panel_member_id | target_metal | success_status | minimum_coordination_number | maximum_mean_metal_oxygen_distance_A | escape_event_detected | final_coordination_number | final_mean_metal_oxygen_distance_A | metadynamics_status |",
            "| --- | --- | --- | ---: | ---: | --- | ---: | ---: | --- |",
        ]
    )
    for row in system_rows:
        lines.append(
            "| "
            f"{row.panel_member_id} | "
            f"{row.target_metal} | "
            f"{row.success_status} | "
            f"{_format_metric(row.minimum_coordination_number)} | "
            f"{_format_metric(row.maximum_mean_metal_oxygen_distance_A)} | "
            f"{row.escape_event_detected} | "
            f"{_format_metric(row.final_coordination_number)} | "
            f"{_format_metric(row.final_mean_metal_oxygen_distance_A)} | "
            f"{row.metadynamics_status} |"
        )
    lines.extend(
        [
            "",
            "## Panel Status",
            "",
            "| panel_member_id | Dy | Nd | Y | Al | Fe | candidate_keep_for_qm |",
            "| --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in panel_rows:
        lines.append(
            "| "
            f"{row.panel_member_id} | "
            f"{row.dy_status} | "
            f"{row.nd_status} | "
            f"{row.y_status} | "
            f"{row.al_status} | "
            f"{row.fe_status} | "
            f"{row.candidate_keep_for_qm} |"
        )
    failure_rows = [row for row in system_rows if row.success_status != "success"]
    if failure_rows:
        lines.extend(
            [
                "",
                "## Failures",
                "",
                "| panel_member_id | target_metal | failure_reason |",
                "| --- | --- | --- |",
            ]
        )
        for row in failure_rows:
            lines.append(
                "| "
                f"{row.panel_member_id} | "
                f"{row.target_metal} | "
                f"{row.failure_reason} |"
            )
    lines.append("")
    return "\n".join(lines)


def _write_metadynamics_log(
    *,
    log_path: Path,
    system_row: ReducedMetadynamicsOpenMMSystem,
    runtime_config: OpenMMMetadynamicsRuntimeConfig,
    artifacts: OpenMMMetadynamicsArtifacts | None,
    status: str,
    failure_reason: str,
    traceback_text: str = "",
) -> None:
    lines = [
        "Phase 6C1 OpenMM metadynamics smoke",
        f"panel_member_id: {system_row.panel_member_id}",
        f"target_metal: {system_row.target_metal}",
        f"topology_class: {system_row.topology_class}",
        f"seed: {system_row.seed}",
        f"target_temperature_K: {runtime_config.target_temperature_K}",
        f"timestep_fs: {runtime_config.timestep_fs:.1f}",
        f"cutoff_nm: {runtime_config.cutoff_nm:.2f}",
        f"duration_ps: {runtime_config.duration_ps:.1f}",
        f"report_stride_steps: {runtime_config.report_stride_steps}",
        f"deposition_frequency_steps: {runtime_config.deposition_frequency_steps}",
        f"bias_factor: {runtime_config.bias_factor:.1f}",
        f"gaussian_height_kj_per_mol: {runtime_config.gaussian_height_kj_per_mol:.2f}",
        f"donor_cutoff_A: {runtime_config.donor_cutoff_A:.1f}",
        f"coordination_switch_distance_A: {runtime_config.coordination_switch_distance_A:.1f}",
        f"status: {status}",
    ]
    if artifacts is not None:
        lines.extend(
            [
                f"platform_name: {artifacts.platform_name}",
                f"platform_properties: {artifacts.platform_properties}",
                f"metal_atom_count: {artifacts.metal_atom_count}",
                f"pocket_oxygen_atom_count: {artifacts.pocket_oxygen_atom_count}",
                f"time_series_frame_count: {len(artifacts.cv_timeseries_rows)}",
                f"minimum_coordination_number: {_format_metric(artifacts.minimum_coordination_number)}",
                "maximum_mean_metal_oxygen_distance_A: "
                f"{_format_metric(artifacts.maximum_mean_metal_oxygen_distance_A)}",
                f"escape_event_detected: {'TRUE' if artifacts.escape_event_detected else 'FALSE'}",
                f"final_coordination_number: {_format_metric(artifacts.final_coordination_number)}",
                "final_mean_metal_oxygen_distance_A: "
                f"{_format_metric(artifacts.final_mean_metal_oxygen_distance_A)}",
            ]
        )
    if failure_reason:
        lines.append(f"failure_reason: {failure_reason}")
    if traceback_text:
        lines.extend(("", "traceback:", traceback_text.rstrip()))
    atomic_write_text(log_path, "\n".join(lines) + "\n")


def _run_metadynamics_task(
    task: MetadynamicsTask,
    metadynamics_runner: Callable[..., OpenMMMetadynamicsArtifacts],
) -> OpenMMMetadynamicsSummaryRow:
    system_label = _format_system_label(task.system_row)
    output_dir = task.output_root / task.system_row.panel_member_id / task.system_row.target_metal
    cv_timeseries_path = output_dir / "cv_timeseries.csv"
    bias_state_path = output_dir / "bias_state.xml"
    final_state_path = output_dir / "metadynamics_final_state.xml"
    final_structure_path = output_dir / "metadynamics_final_structure.pdb"
    log_path = output_dir / "metadynamics_log.txt"
    existing_summary_row = _build_summary_row_from_existing_outputs(
        system_row=task.system_row,
        cv_timeseries_path=cv_timeseries_path,
        bias_state_path=bias_state_path,
        final_state_path=final_state_path,
        final_structure_path=final_structure_path,
        log_path=log_path,
    )
    if existing_summary_row is not None:
        print(f"[openmm_metadynamics_smoke] skipping completed system {system_label}", flush=True)
        return existing_summary_row
    for artifact_path in (
        cv_timeseries_path,
        bias_state_path,
        final_state_path,
        final_structure_path,
        log_path,
    ):
        artifact_path.unlink(missing_ok=True)

    artifacts: OpenMMMetadynamicsArtifacts | None = None
    print(f"[openmm_metadynamics_smoke] starting system {system_label}", flush=True)
    try:
        artifacts = metadynamics_runner(
            equilibrated_structure_pdb_text=task.system_row.equilibrated_structure_path.read_text(encoding="utf-8"),
            system_xml=task.system_row.system_xml_path.read_text(encoding="utf-8"),
            equilibrated_state_xml=task.system_row.equilibrated_final_state_path.read_text(encoding="utf-8"),
            target_metal=task.system_row.target_metal,
            target_temperature_K=task.runtime_config.target_temperature_K,
            friction_coeff_ps=task.runtime_config.friction_coeff_ps,
            timestep_fs=task.runtime_config.timestep_fs,
            cutoff_nm=task.runtime_config.cutoff_nm,
            duration_ps=task.runtime_config.duration_ps,
            report_stride_steps=task.runtime_config.report_stride_steps,
            bias_factor=task.runtime_config.bias_factor,
            gaussian_height_kj_per_mol=task.runtime_config.gaussian_height_kj_per_mol,
            deposition_frequency_steps=task.runtime_config.deposition_frequency_steps,
            donor_cutoff_A=task.runtime_config.donor_cutoff_A,
            coordination_switch_distance_A=task.runtime_config.coordination_switch_distance_A,
            coordination_switch_power=task.runtime_config.coordination_switch_power,
            coordination_number_min=task.runtime_config.coordination_number_min,
            coordination_number_max=task.runtime_config.coordination_number_max,
            coordination_number_bias_width=task.runtime_config.coordination_number_bias_width,
            mean_distance_min_A=task.runtime_config.mean_distance_min_A,
            mean_distance_max_A=task.runtime_config.mean_distance_max_A,
            mean_distance_bias_width_A=task.runtime_config.mean_distance_bias_width_A,
            escape_coordination_threshold=task.runtime_config.escape_coordination_threshold,
            escape_distance_threshold_A=task.runtime_config.escape_distance_threshold_A,
            seed=task.system_row.seed,
            cv_timeseries_path=cv_timeseries_path,
            bias_state_path=bias_state_path,
            final_state_path=final_state_path,
            final_structure_path=final_structure_path,
            cpu_thread_count=task.cpu_thread_count,
        )
        _write_metadynamics_log(
            log_path=log_path,
            system_row=task.system_row,
            runtime_config=task.runtime_config,
            artifacts=artifacts,
            status="success",
            failure_reason="",
        )
        print(f"[openmm_metadynamics_smoke] finished system {system_label}", flush=True)
        return build_openmm_metadynamics_summary_row(
            system_row=task.system_row,
            artifacts=artifacts,
            failure_reason="",
        )
    except Exception as exc:
        failure_reason = f"{type(exc).__name__}: {exc}"
        _write_metadynamics_log(
            log_path=log_path,
            system_row=task.system_row,
            runtime_config=task.runtime_config,
            artifacts=artifacts,
            status="failure",
            failure_reason=failure_reason,
            traceback_text=traceback.format_exc(),
        )
        print(
            f"[openmm_metadynamics_smoke] failed system {system_label}: {failure_reason}",
            flush=True,
        )
        return build_openmm_metadynamics_summary_row(
            system_row=task.system_row,
            artifacts=None,
            failure_reason=failure_reason,
        )


def _run_default_metadynamics_task(task: MetadynamicsTask) -> OpenMMMetadynamicsSummaryRow:
    return _run_metadynamics_task(task, run_openmm_metadynamics_smoke_system)


def run_openmm_metadynamics_smoke(
    *,
    screening_panel_status_path: Path,
    system_build_summary_path: Path,
    output_root: Path,
    system_summary_path: Path,
    panel_status_path: Path,
    report_path: Path,
    metadynamics_runner: Callable[..., OpenMMMetadynamicsArtifacts] = run_openmm_metadynamics_smoke_system,
    max_workers: int = DEFAULT_METADYNAMICS_WORKER_COUNT,
) -> tuple[OpenMMMetadynamicsSummaryRow, ...]:
    """Run Phase 6C1 across the exact reduced metadynamics panel."""
    if max_workers < 1:
        raise ValueError(f"max_workers must be at least 1, found {max_workers}")
    systems = discover_reduced_openmm_metadynamics_systems(
        screening_panel_status_path=screening_panel_status_path,
        system_build_summary_path=system_build_summary_path,
    )
    runtime_reference: OpenMMMetadynamicsRuntimeConfig | None = None
    tasks: list[MetadynamicsTask] = []
    for system_row in systems:
        runtime_config = build_openmm_metadynamics_runtime_config(system_row.simulation_config_path)
        runtime_reference = runtime_reference or runtime_config
        tasks.append(
            MetadynamicsTask(
                system_row=system_row,
                runtime_config=runtime_config,
                output_root=output_root,
                cpu_thread_count=DEFAULT_METADYNAMICS_CPU_THREAD_COUNT,
            )
        )

    if max_workers == 1:
        if metadynamics_runner is run_openmm_metadynamics_smoke_system:
            system_rows = tuple(_run_default_metadynamics_task(task) for task in tasks)
        else:
            system_rows = tuple(_run_metadynamics_task(task, metadynamics_runner) for task in tasks)
    elif metadynamics_runner is run_openmm_metadynamics_smoke_system:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            system_rows = tuple(executor.map(_run_default_metadynamics_task, tasks))
    else:
        system_rows = tuple(_run_metadynamics_task(task, metadynamics_runner) for task in tasks)

    panel_rows = build_openmm_metadynamics_panel_status_rows(system_rows)
    write_csv_rows(system_summary_path, system_rows)
    write_csv_rows(panel_status_path, panel_rows)
    if runtime_reference is None:  # pragma: no cover - defensive
        runtime_reference = OpenMMMetadynamicsRuntimeConfig(
            target_temperature_K=298,
            friction_coeff_ps=1.0,
            timestep_fs=2.0,
            cutoff_nm=1.0,
            duration_ps=DEFAULT_METADYNAMICS_DURATION_PS,
            report_stride_steps=DEFAULT_METADYNAMICS_REPORT_STRIDE_STEPS,
            bias_factor=DEFAULT_METADYNAMICS_BIAS_FACTOR,
            gaussian_height_kj_per_mol=DEFAULT_METADYNAMICS_GAUSSIAN_HEIGHT_KJ_PER_MOL,
            deposition_frequency_steps=DEFAULT_METADYNAMICS_DEPOSITION_FREQUENCY_STEPS,
            donor_cutoff_A=DEFAULT_DONOR_CUTOFF_A,
            coordination_switch_distance_A=DEFAULT_COORDINATION_SWITCH_DISTANCE_A,
            coordination_switch_power=DEFAULT_COORDINATION_SWITCH_POWER,
            coordination_number_min=DEFAULT_COORDINATION_NUMBER_MIN,
            coordination_number_max=DEFAULT_COORDINATION_NUMBER_MAX,
            coordination_number_bias_width=DEFAULT_COORDINATION_NUMBER_BIAS_WIDTH,
            mean_distance_min_A=DEFAULT_MEAN_DISTANCE_MIN_A,
            mean_distance_max_A=DEFAULT_MEAN_DISTANCE_MAX_A,
            mean_distance_bias_width_A=DEFAULT_MEAN_DISTANCE_BIAS_WIDTH_A,
            escape_coordination_threshold=DEFAULT_ESCAPE_COORDINATION_THRESHOLD,
            escape_distance_threshold_A=DEFAULT_ESCAPE_DISTANCE_THRESHOLD_A,
        )
    atomic_write_text(
        report_path,
        render_openmm_metadynamics_smoke_report(
            system_rows=system_rows,
            panel_rows=panel_rows,
            runtime_config=runtime_reference,
        ),
    )
    failure_count = sum(row.success_status != "success" for row in system_rows)
    if failure_count:
        raise RuntimeError(
            "OpenMM metadynamics smoke failed for "
            f"{failure_count}/{len(system_rows)} systems; see {_display_path(system_summary_path)}"
        )
    return system_rows
