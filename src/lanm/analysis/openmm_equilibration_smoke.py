"""Phase 6B1 deterministic minimization/NVT/NPT smoke runs for baseline OpenMM systems."""

from __future__ import annotations

import csv
import traceback
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable

import yaml

from lanm.analysis.openmm_system_build import SimulationConfig
from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_text, write_csv_rows
from lanm.md.openmm_equilibration_smoke import (
    DEFAULT_NPT_DURATION_PS,
    DEFAULT_NVT_DURATION_PS,
    DEFAULT_PRESSURE_ATM,
    EquilibrationLogRow,
    OpenMMEquilibrationSmokeArtifacts,
    run_openmm_equilibration_smoke_protocol,
)
from lanm.paths import REPO_ROOT


@dataclass(frozen=True, slots=True)
class OpenMMBuildSummaryRow:
    panel_member_id: str
    target_metal: str
    topology_class: str
    success_status: str
    failure_reason: str
    atom_count: int
    residue_count: int
    system_xml_path: str
    simulation_config_path: str


@dataclass(frozen=True, slots=True)
class BaselineBuiltOpenMMSystem:
    panel_member_id: str
    target_metal: str
    topology_class: str
    system_xml_path: Path
    simulation_config_path: Path
    prepared_structure_path: Path


@dataclass(frozen=True, slots=True)
class OpenMMEquilibrationSmokeRuntimeConfig:
    smoke_replicate_count: int
    target_temperature_K: int
    friction_coeff_ps: float
    timestep_fs: float
    cutoff_nm: float
    nvt_duration_ps: float
    npt_duration_ps: float
    pressure_atm: float


@dataclass(frozen=True, slots=True)
class OpenMMEquilibrationSmokeSummaryRow:
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


def _format_metric(value: float | None) -> str:
    if value is None:
        return "-"
    return f"{value:.6f}"


def load_openmm_system_build_summary(path: Path) -> tuple[OpenMMBuildSummaryRow, ...]:
    """Load the Phase 6A4 OpenMM build summary table."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = tuple(
            OpenMMBuildSummaryRow(
                panel_member_id=str(row["panel_member_id"]).strip(),
                target_metal=str(row["target_metal"]).strip(),
                topology_class=str(row["topology_class"]).strip(),
                success_status=str(row["success_status"]).strip(),
                failure_reason=str(row["failure_reason"]).strip(),
                atom_count=int(str(row["atom_count"]).strip() or 0),
                residue_count=int(str(row["residue_count"]).strip() or 0),
                system_xml_path=str(row["system_xml_path"]).strip(),
                simulation_config_path=str(row["simulation_config_path"]).strip(),
            )
            for row in reader
        )
    if not rows:
        raise ValueError(f"No rows found in {path}")
    return rows


def discover_successful_baseline_openmm_systems(path: Path) -> tuple[BaselineBuiltOpenMMSystem, ...]:
    """Resolve the successful baseline-built systems that should enter Phase 6B1."""
    build_rows = load_openmm_system_build_summary(path)
    successful_systems: list[BaselineBuiltOpenMMSystem] = []
    for row in build_rows:
        if row.success_status != "success":
            continue
        system_xml_path = _repo_path(row.system_xml_path)
        simulation_config_path = _repo_path(row.simulation_config_path)
        prepared_structure_path = system_xml_path.parent / "prepared_structure.pdb"
        require_path(system_xml_path)
        require_path(simulation_config_path)
        require_path(prepared_structure_path)
        successful_systems.append(
            BaselineBuiltOpenMMSystem(
                panel_member_id=row.panel_member_id,
                target_metal=row.target_metal,
                topology_class=row.topology_class,
                system_xml_path=system_xml_path,
                simulation_config_path=simulation_config_path,
                prepared_structure_path=prepared_structure_path,
            )
        )
    if not successful_systems:
        raise ValueError(f"No successful Phase 6A4 OpenMM systems found in {path}")
    return tuple(successful_systems)


def load_simulation_config(path: Path) -> SimulationConfig:
    """Load a per-system Phase 6A4 simulation config for Phase 6B1 reuse."""
    payload = _parse_required_mapping(path)
    return SimulationConfig(
        panel_member_id=str(payload["panel_member_id"]).strip(),
        target_metal=str(payload["target_metal"]).strip(),
        topology_class=str(payload["topology_class"]).strip(),
        replicate_count=int(payload["replicate_count"]),
        target_temperature_K=int(payload["target_temperature_K"]),
        friction_coeff_ps=float(payload["friction_coeff_ps"]),
        timestep_fs=float(payload["timestep_fs"]),
        nonbonded_method=str(payload["nonbonded_method"]).strip(),
        cutoff_nm=float(payload["cutoff_nm"]),
        hydrogen_mass_repartitioning=bool(payload["hydrogen_mass_repartitioning"]),
        protein_forcefield=str(payload["protein_forcefield"]).strip(),
        water_model=str(payload["water_model"]).strip(),
        metal_parameter_family=str(payload["metal_parameter_family"]).strip(),
        metal_parameter_provenance=str(payload["metal_parameter_provenance"]).strip(),
    )


def build_equilibration_runtime_config(simulation_config: SimulationConfig) -> OpenMMEquilibrationSmokeRuntimeConfig:
    """Translate Phase 6A4 config values into the fixed Phase 6B1 smoke protocol."""
    if simulation_config.target_temperature_K != 298:
        raise ValueError(
            f"Phase 6B1 requires target_temperature_K=298, found {simulation_config.target_temperature_K}"
        )
    if simulation_config.timestep_fs != 2.0:
        raise ValueError(f"Phase 6B1 requires timestep_fs=2.0, found {simulation_config.timestep_fs}")
    return OpenMMEquilibrationSmokeRuntimeConfig(
        smoke_replicate_count=1,
        target_temperature_K=simulation_config.target_temperature_K,
        friction_coeff_ps=simulation_config.friction_coeff_ps,
        timestep_fs=2.0,
        cutoff_nm=simulation_config.cutoff_nm,
        nvt_duration_ps=DEFAULT_NVT_DURATION_PS,
        npt_duration_ps=DEFAULT_NPT_DURATION_PS,
        pressure_atm=DEFAULT_PRESSURE_ATM,
    )


def build_equilibration_summary_row(
    *,
    system_row: BaselineBuiltOpenMMSystem,
    output_dir: Path,
    artifacts: OpenMMEquilibrationSmokeArtifacts | None,
    failure_reason: str,
) -> OpenMMEquilibrationSmokeSummaryRow:
    """Create one summary-table row for a Phase 6B1 system attempt."""
    if artifacts is None:
        return OpenMMEquilibrationSmokeSummaryRow(
            panel_member_id=system_row.panel_member_id,
            target_metal=system_row.target_metal,
            topology_class=system_row.topology_class,
            success_status="failure",
            failure_reason=failure_reason,
            minimized_potential_energy=None,
            final_nvt_potential_energy=None,
            final_npt_potential_energy=None,
            final_temperature_K=None,
            final_box_volume_nm3=None,
            final_structure_path="",
        )
    return OpenMMEquilibrationSmokeSummaryRow(
        panel_member_id=system_row.panel_member_id,
        target_metal=system_row.target_metal,
        topology_class=system_row.topology_class,
        success_status="success",
        failure_reason="",
        minimized_potential_energy=artifacts.minimized_potential_energy,
        final_nvt_potential_energy=artifacts.final_nvt_potential_energy,
        final_npt_potential_energy=artifacts.final_npt_potential_energy,
        final_temperature_K=artifacts.final_temperature_K,
        final_box_volume_nm3=artifacts.final_box_volume_nm3,
        final_structure_path=_display_path(output_dir / "npt_final.pdb"),
    )


def render_openmm_equilibration_smoke_report(
    results: tuple[OpenMMEquilibrationSmokeSummaryRow, ...],
    runtime_config: OpenMMEquilibrationSmokeRuntimeConfig,
) -> str:
    """Render the concise Phase 6B1 Markdown smoke-equilibration report."""
    success_count = sum(row.success_status == "success" for row in results)
    failure_count = len(results) - success_count
    lines = [
        "# OpenMM Equilibration Smoke",
        "",
        "Phase 6B1 runs one deterministic minimization plus short NVT/NPT equilibration smoke replicate for every successful Phase 6A4 baseline OpenMM build.",
        "",
        "## Smoke Protocol",
        "",
        f"- baseline-built systems attempted: `{len(results)}`",
        f"- successful smoke runs: `{success_count}`",
        f"- failures recorded: `{failure_count}`",
        f"- deterministic smoke replicate count per system: `{runtime_config.smoke_replicate_count}`",
        f"- target temperature: `{runtime_config.target_temperature_K} K`",
        f"- pressure target: `{runtime_config.pressure_atm:.1f} atm`",
        f"- timestep: `{runtime_config.timestep_fs:.1f} fs`",
        f"- NVT duration: `{runtime_config.nvt_duration_ps:.1f} ps`",
        f"- NPT duration: `{runtime_config.npt_duration_ps:.1f} ps`",
        (
            "- NPT note: baseline systems are serialized as non-periodic `NoCutoff` systems, "
            "so Phase 6B1 promotes only the pressure-coupled leg to a deterministic orthorhombic periodic box."
        ),
        "",
        "| panel_member_id | target_metal | topology_class | success_status | minimized_potential_energy | final_nvt_potential_energy | final_npt_potential_energy | final_temperature_K | final_box_volume_nm3 | failure_reason |",
        "| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in results:
        lines.append(
            "| "
            f"{row.panel_member_id} | "
            f"{row.target_metal} | "
            f"{row.topology_class} | "
            f"{row.success_status} | "
            f"{_format_metric(row.minimized_potential_energy)} | "
            f"{_format_metric(row.final_nvt_potential_energy)} | "
            f"{_format_metric(row.final_npt_potential_energy)} | "
            f"{_format_metric(row.final_temperature_K)} | "
            f"{_format_metric(row.final_box_volume_nm3)} | "
            f"{row.failure_reason or '-'} |"
        )
    lines.append("")
    return "\n".join(lines)


def _write_failure_log(log_path: Path, failure_reason: str, traceback_text: str) -> None:
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


def run_openmm_equilibration_smoke(
    *,
    system_build_summary_path: Path,
    output_root: Path,
    summary_path: Path,
    report_path: Path,
    equilibrator: Callable[..., OpenMMEquilibrationSmokeArtifacts] = run_openmm_equilibration_smoke_protocol,
) -> tuple[OpenMMEquilibrationSmokeSummaryRow, ...]:
    """Run Phase 6B1 on every successful Phase 6A4 baseline OpenMM build."""
    systems = discover_successful_baseline_openmm_systems(system_build_summary_path)
    results: list[OpenMMEquilibrationSmokeSummaryRow] = []
    runtime_reference: OpenMMEquilibrationSmokeRuntimeConfig | None = None

    for system_row in systems:
        output_dir = output_root / system_row.panel_member_id / system_row.target_metal
        minimized_structure_path = output_dir / "minimized_structure.pdb"
        nvt_final_path = output_dir / "nvt_final.pdb"
        npt_final_path = output_dir / "npt_final.pdb"
        equilibration_log_path = output_dir / "equilibration_log.csv"
        final_state_path = output_dir / "final_state.xml"
        for artifact_path in (
            minimized_structure_path,
            nvt_final_path,
            npt_final_path,
            equilibration_log_path,
            final_state_path,
        ):
            artifact_path.unlink(missing_ok=True)

        artifacts: OpenMMEquilibrationSmokeArtifacts | None = None
        try:
            simulation_config = load_simulation_config(system_row.simulation_config_path)
            runtime_config = build_equilibration_runtime_config(simulation_config)
            runtime_reference = runtime_reference or runtime_config
            artifacts = equilibrator(
                prepared_structure_pdb_text=system_row.prepared_structure_path.read_text(encoding="utf-8"),
                system_xml=system_row.system_xml_path.read_text(encoding="utf-8"),
                target_temperature_K=runtime_config.target_temperature_K,
                friction_coeff_ps=runtime_config.friction_coeff_ps,
                timestep_fs=runtime_config.timestep_fs,
                cutoff_nm=runtime_config.cutoff_nm,
                nvt_duration_ps=runtime_config.nvt_duration_ps,
                npt_duration_ps=runtime_config.npt_duration_ps,
                pressure_atm=runtime_config.pressure_atm,
            )
            atomic_write_text(minimized_structure_path, artifacts.minimized_structure_pdb_text)
            atomic_write_text(nvt_final_path, artifacts.nvt_final_pdb_text)
            atomic_write_text(npt_final_path, artifacts.npt_final_pdb_text)
            atomic_write_text(final_state_path, artifacts.final_state_xml)
            write_csv_rows(equilibration_log_path, artifacts.log_rows)
            results.append(
                build_equilibration_summary_row(
                    system_row=system_row,
                    output_dir=output_dir,
                    artifacts=artifacts,
                    failure_reason="",
                )
            )
        except Exception as exc:
            failure_reason = f"{type(exc).__name__}: {exc}"
            traceback_text = traceback.format_exc()
            _write_failure_log(equilibration_log_path, failure_reason, traceback_text)
            results.append(
                build_equilibration_summary_row(
                    system_row=system_row,
                    output_dir=output_dir,
                    artifacts=None,
                    failure_reason=failure_reason,
                )
            )

    result_tuple = tuple(results)
    write_csv_rows(summary_path, result_tuple)
    if runtime_reference is None:
        runtime_reference = OpenMMEquilibrationSmokeRuntimeConfig(
            smoke_replicate_count=1,
            target_temperature_K=298,
            friction_coeff_ps=1.0,
            timestep_fs=2.0,
            cutoff_nm=1.0,
            nvt_duration_ps=DEFAULT_NVT_DURATION_PS,
            npt_duration_ps=DEFAULT_NPT_DURATION_PS,
            pressure_atm=DEFAULT_PRESSURE_ATM,
        )
    atomic_write_text(report_path, render_openmm_equilibration_smoke_report(result_tuple, runtime_reference))
    if any(row.success_status != "success" for row in result_tuple):
        failure_count = sum(row.success_status != "success" for row in result_tuple)
        raise RuntimeError(
            f"OpenMM equilibration smoke failed for {failure_count}/{len(result_tuple)} systems; "
            f"see {_display_path(summary_path)}"
        )
    return result_tuple
