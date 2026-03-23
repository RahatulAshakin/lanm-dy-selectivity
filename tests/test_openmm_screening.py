from __future__ import annotations

import csv
from pathlib import Path

import pytest
import yaml

from lanm.analysis.openmm_screening import (
    OpenMMScreeningReplicateSummaryRow,
    OpenMMScreeningRuntimeConfig,
    build_openmm_screening_panel_status_rows,
    build_openmm_screening_replicate_schedule,
    build_openmm_screening_runtime_config,
    build_openmm_screening_system_summaries,
    discover_successful_openmm_screening_systems,
    run_openmm_screening,
)
from lanm.md.openmm_screening import OpenMMScreeningArtifacts
from lanm.paths import (
    OPENMM_EQUILIBRATION_SMOKE_SUMMARY_PATH,
    OPENMM_SYSTEM_BUILD_SUMMARY_PATH,
)


def _write_fake_phase_6b1_system(
    tmp_path: Path,
    *,
    panel_member_id: str,
    target_metal: str,
    topology_class: str = "am1_monomer",
) -> tuple[dict[str, str], dict[str, str]]:
    build_dir = tmp_path / "results" / "openmm_system_build" / panel_member_id / target_metal
    build_dir.mkdir(parents=True, exist_ok=True)
    (build_dir / "system.xml").write_text("<System/>", encoding="utf-8")
    (build_dir / "prepared_structure.pdb").write_text(
        "ATOM      1  O   HOH A   1       0.0   0.0   0.0\n"
        "HETATM    2 DY   DY  A 201       2.5   0.0   0.0\n"
        "END\n",
        encoding="utf-8",
    )
    with (build_dir / "simulation_config.yaml").open("w", encoding="utf-8") as handle:
        yaml.safe_dump(
            {
                "panel_member_id": panel_member_id,
                "target_metal": target_metal,
                "topology_class": topology_class,
                "replicate_count": 3,
                "target_temperature_K": 298,
                "friction_coeff_ps": 1.0,
                "timestep_fs": 2.0,
                "nonbonded_method": "NoCutoff",
                "cutoff_nm": 1.0,
                "hydrogen_mass_repartitioning": False,
                "protein_forcefield": "amber19-all.xml",
                "water_model": "amber19/opc3.xml",
                "metal_parameter_family": "generic_12_6_4_highly_charged",
                "metal_parameter_provenance": "published_opc3_12_6_4_baseline",
            },
            handle,
            sort_keys=False,
        )

    equilibration_dir = tmp_path / "results" / "openmm_equilibration_smoke" / panel_member_id / target_metal
    equilibration_dir.mkdir(parents=True, exist_ok=True)
    (equilibration_dir / "final_state.xml").write_text("<State/>", encoding="utf-8")
    (equilibration_dir / "npt_final.pdb").write_text(
        "ATOM      1  O   HOH A   1       0.0   0.0   0.0\nEND\n",
        encoding="utf-8",
    )

    build_row = {
        "panel_member_id": panel_member_id,
        "target_metal": target_metal,
        "topology_class": topology_class,
        "success_status": "success",
        "failure_reason": "",
        "atom_count": "100",
        "residue_count": "10",
        "system_xml_path": str(build_dir / "system.xml"),
        "simulation_config_path": str(build_dir / "simulation_config.yaml"),
    }
    equilibration_row = {
        "panel_member_id": panel_member_id,
        "target_metal": target_metal,
        "topology_class": topology_class,
        "success_status": "success",
        "failure_reason": "",
        "minimized_potential_energy": "-10.0",
        "final_nvt_potential_energy": "-9.0",
        "final_npt_potential_energy": "-8.0",
        "final_temperature_K": "298.0",
        "final_box_volume_nm3": "110.0",
        "final_structure_path": str(equilibration_dir / "npt_final.pdb"),
    }
    return build_row, equilibration_row


def _write_build_summary(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "panel_member_id",
                "target_metal",
                "topology_class",
                "success_status",
                "failure_reason",
                "atom_count",
                "residue_count",
                "system_xml_path",
                "simulation_config_path",
            ),
        )
        writer.writeheader()
        writer.writerows(rows)


def _write_equilibration_summary(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "panel_member_id",
                "target_metal",
                "topology_class",
                "success_status",
                "failure_reason",
                "minimized_potential_energy",
                "final_nvt_potential_energy",
                "final_npt_potential_energy",
                "final_temperature_K",
                "final_box_volume_nm3",
                "final_structure_path",
            ),
        )
        writer.writeheader()
        writer.writerows(rows)


def _write_fake_state_data_csv(path: Path, *, temperature_offset: float, volume_offset: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("Step", "Time (ps)", "Temperature (K)", "Box Volume (nm^3)"),
        )
        writer.writeheader()
        writer.writerows(
            (
                {
                    "Step": 500,
                    "Time (ps)": 1.0,
                    "Temperature (K)": 298.0 + temperature_offset,
                    "Box Volume (nm^3)": 110.0 + volume_offset,
                },
                {
                    "Step": 1000,
                    "Time (ps)": 2.0,
                    "Temperature (K)": 299.0 + temperature_offset,
                    "Box Volume (nm^3)": 111.0 + volume_offset,
                },
            )
        )


def test_discover_successful_openmm_screening_systems_matches_full_panel() -> None:
    systems = discover_successful_openmm_screening_systems(
        system_build_summary_path=OPENMM_SYSTEM_BUILD_SUMMARY_PATH,
        equilibration_summary_path=OPENMM_EQUILIBRATION_SMOKE_SUMMARY_PATH,
    )

    assert len(systems) == 30
    assert systems[0].panel_member_id == "am1_mex_ss_only_u02"
    assert systems[0].target_metal == "Al"
    assert systems[-1].panel_member_id == "hans_interface_ss_plus_if_u02"
    assert systems[-1].target_metal == "Y"
    assert all(system.system_xml_path.exists() for system in systems)
    assert all(system.equilibration_final_state_path.exists() for system in systems)


def test_build_openmm_screening_replicate_schedule_uses_phase_6b2_seeds() -> None:
    systems = discover_successful_openmm_screening_systems(
        system_build_summary_path=OPENMM_SYSTEM_BUILD_SUMMARY_PATH,
        equilibration_summary_path=OPENMM_EQUILIBRATION_SMOKE_SUMMARY_PATH,
    )
    runtime_config = build_openmm_screening_runtime_config(systems[0].simulation_config_path)

    schedule = build_openmm_screening_replicate_schedule(runtime_config)

    assert runtime_config.production_duration_ps == 100.0
    assert runtime_config.pressure_atm == 1.0
    assert runtime_config.timestep_fs == 2.0
    assert [spec.replicate_id for spec in schedule] == [1, 2, 3]
    assert [spec.seed for spec in schedule] == [101, 102, 103]


def test_build_openmm_screening_system_summaries_applies_status_rules() -> None:
    replicate_rows = (
        OpenMMScreeningReplicateSummaryRow(
            panel_member_id="candidate_a",
            target_metal="Dy",
            topology_class="am1_monomer",
            replicate_id=1,
            seed=101,
            success_status="success",
            failure_reason="",
            mean_temperature_K=298.1,
            final_temperature_K=298.5,
            mean_box_volume_nm3=111.0,
            mean_min_metal_oxygen_distance_A=2.8,
            inner_sphere_occupancy_fraction=0.8,
            state_data_path="a",
            final_state_path="b",
            final_structure_path="c",
            screening_log_path="d",
            trajectory_path="e",
        ),
        OpenMMScreeningReplicateSummaryRow(
            panel_member_id="candidate_a",
            target_metal="Dy",
            topology_class="am1_monomer",
            replicate_id=2,
            seed=102,
            success_status="success",
            failure_reason="",
            mean_temperature_K=298.0,
            final_temperature_K=298.4,
            mean_box_volume_nm3=111.5,
            mean_min_metal_oxygen_distance_A=2.9,
            inner_sphere_occupancy_fraction=0.7,
            state_data_path="a",
            final_state_path="b",
            final_structure_path="c",
            screening_log_path="d",
            trajectory_path="e",
        ),
        OpenMMScreeningReplicateSummaryRow(
            panel_member_id="candidate_a",
            target_metal="Dy",
            topology_class="am1_monomer",
            replicate_id=3,
            seed=103,
            success_status="success",
            failure_reason="",
            mean_temperature_K=297.9,
            final_temperature_K=298.3,
            mean_box_volume_nm3=110.9,
            mean_min_metal_oxygen_distance_A=3.0,
            inner_sphere_occupancy_fraction=0.6,
            state_data_path="a",
            final_state_path="b",
            final_structure_path="c",
            screening_log_path="d",
            trajectory_path="e",
        ),
        OpenMMScreeningReplicateSummaryRow(
            panel_member_id="candidate_a",
            target_metal="Fe",
            topology_class="am1_monomer",
            replicate_id=1,
            seed=101,
            success_status="success",
            failure_reason="",
            mean_temperature_K=298.1,
            final_temperature_K=298.1,
            mean_box_volume_nm3=109.0,
            mean_min_metal_oxygen_distance_A=4.1,
            inner_sphere_occupancy_fraction=0.2,
            state_data_path="a",
            final_state_path="b",
            final_structure_path="c",
            screening_log_path="d",
            trajectory_path="e",
        ),
        OpenMMScreeningReplicateSummaryRow(
            panel_member_id="candidate_a",
            target_metal="Fe",
            topology_class="am1_monomer",
            replicate_id=2,
            seed=102,
            success_status="success",
            failure_reason="",
            mean_temperature_K=298.1,
            final_temperature_K=298.1,
            mean_box_volume_nm3=109.1,
            mean_min_metal_oxygen_distance_A=4.2,
            inner_sphere_occupancy_fraction=0.3,
            state_data_path="a",
            final_state_path="b",
            final_structure_path="c",
            screening_log_path="d",
            trajectory_path="e",
        ),
        OpenMMScreeningReplicateSummaryRow(
            panel_member_id="candidate_a",
            target_metal="Fe",
            topology_class="am1_monomer",
            replicate_id=3,
            seed=103,
            success_status="success",
            failure_reason="",
            mean_temperature_K=298.1,
            final_temperature_K=298.1,
            mean_box_volume_nm3=109.2,
            mean_min_metal_oxygen_distance_A=4.3,
            inner_sphere_occupancy_fraction=0.4,
            state_data_path="a",
            final_state_path="b",
            final_structure_path="c",
            screening_log_path="d",
            trajectory_path="e",
        ),
        OpenMMScreeningReplicateSummaryRow(
            panel_member_id="candidate_b",
            target_metal="Y",
            topology_class="hans_monomer",
            replicate_id=1,
            seed=101,
            success_status="success",
            failure_reason="",
            mean_temperature_K=298.1,
            final_temperature_K=298.1,
            mean_box_volume_nm3=120.0,
            mean_min_metal_oxygen_distance_A=3.5,
            inner_sphere_occupancy_fraction=0.4,
            state_data_path="a",
            final_state_path="b",
            final_structure_path="c",
            screening_log_path="d",
            trajectory_path="e",
        ),
        OpenMMScreeningReplicateSummaryRow(
            panel_member_id="candidate_b",
            target_metal="Y",
            topology_class="hans_monomer",
            replicate_id=2,
            seed=102,
            success_status="failure",
            failure_reason="RuntimeError: synthetic failure",
            mean_temperature_K=None,
            final_temperature_K=None,
            mean_box_volume_nm3=None,
            mean_min_metal_oxygen_distance_A=None,
            inner_sphere_occupancy_fraction=None,
            state_data_path="",
            final_state_path="",
            final_structure_path="",
            screening_log_path="d",
            trajectory_path="",
        ),
        OpenMMScreeningReplicateSummaryRow(
            panel_member_id="candidate_b",
            target_metal="Y",
            topology_class="hans_monomer",
            replicate_id=3,
            seed=103,
            success_status="success",
            failure_reason="",
            mean_temperature_K=298.1,
            final_temperature_K=298.1,
            mean_box_volume_nm3=121.0,
            mean_min_metal_oxygen_distance_A=3.7,
            inner_sphere_occupancy_fraction=0.3,
            state_data_path="a",
            final_state_path="b",
            final_structure_path="c",
            screening_log_path="d",
            trajectory_path="e",
        ),
    )

    system_rows = build_openmm_screening_system_summaries(replicate_rows)
    panel_rows = build_openmm_screening_panel_status_rows(system_rows)

    assert len(system_rows) == 3
    dy_row = next(row for row in system_rows if row.target_metal == "Dy")
    fe_row = next(row for row in system_rows if row.target_metal == "Fe")
    y_row = next(row for row in system_rows if row.target_metal == "Y")
    assert dy_row.replicate_success_count == 3
    assert dy_row.screening_status == "stable_bound"
    assert dy_row.median_inner_sphere_occupancy_fraction == 0.7
    assert fe_row.screening_status == "no_persistent_capture"
    assert y_row.replicate_success_count == 2
    assert y_row.screening_status == "failed"

    assert len(panel_rows) == 2
    candidate_a = next(row for row in panel_rows if row.panel_member_id == "candidate_a")
    assert candidate_a.dy_status == "stable_bound"
    assert candidate_a.fe_status == "no_persistent_capture"
    assert candidate_a.stable_bound_rare_earth_count == 1


def test_run_openmm_screening_continues_after_replicate_failure_and_writes_summaries(tmp_path: Path) -> None:
    build_row_a, equilibration_row_a = _write_fake_phase_6b1_system(
        tmp_path,
        panel_member_id="candidate_a",
        target_metal="Dy",
    )
    build_row_b, equilibration_row_b = _write_fake_phase_6b1_system(
        tmp_path,
        panel_member_id="candidate_b",
        target_metal="Nd",
    )
    build_summary_path = tmp_path / "results" / "tables" / "openmm_system_build_summary.csv"
    equilibration_summary_path = tmp_path / "results" / "tables" / "openmm_equilibration_smoke_summary.csv"
    replicate_summary_path = tmp_path / "results" / "tables" / "openmm_screening_replicates.csv"
    system_summary_path = tmp_path / "results" / "tables" / "openmm_screening_summary.csv"
    panel_status_path = tmp_path / "results" / "tables" / "openmm_screening_panel_status.csv"
    report_path = tmp_path / "results" / "reports" / "openmm_screening.md"
    _write_build_summary(build_summary_path, [build_row_a, build_row_b])
    _write_equilibration_summary(equilibration_summary_path, [equilibration_row_a, equilibration_row_b])

    call_order: list[tuple[str, int]] = []

    def fake_screening_runner(**kwargs: object) -> OpenMMScreeningArtifacts:
        target_metal = str(kwargs["target_metal"])
        seed = int(kwargs["seed"])
        state_data_path = Path(kwargs["state_data_path"])
        final_state_path = Path(kwargs["final_state_path"])
        final_structure_path = Path(kwargs["final_structure_path"])
        trajectory_path = Path(kwargs["trajectory_path"])
        call_order.append((state_data_path.parent.parent.name, seed))
        if "candidate_a/Dy/replicate_2" in state_data_path.as_posix():
            raise RuntimeError("synthetic screening failure")
        _write_fake_state_data_csv(
            state_data_path,
            temperature_offset=float(seed - 100),
            volume_offset=0.5 if target_metal == "Nd" else 0.0,
        )
        final_state_path.write_text("<State/>", encoding="utf-8")
        final_structure_path.write_text("ATOM      1  O   HOH A   1       0.0   0.0   0.0\nEND\n", encoding="utf-8")
        trajectory_path.write_bytes(b"DCD")
        return OpenMMScreeningArtifacts(
            platform_name="CPU",
            platform_properties="Threads=16",
            metal_atom_count=1,
            oxygen_atom_count=4,
            frame_count=2,
            mean_min_metal_oxygen_distance_A=2.5 if target_metal == "Nd" else 2.8,
            inner_sphere_occupancy_fraction=0.75 if target_metal == "Nd" else 0.65,
        )

    with pytest.raises(RuntimeError, match="failed for 1/6 replicates"):
        run_openmm_screening(
            system_build_summary_path=build_summary_path,
            equilibration_summary_path=equilibration_summary_path,
            output_root=tmp_path / "results" / "openmm_screening",
            replicate_summary_path=replicate_summary_path,
            system_summary_path=system_summary_path,
            panel_status_path=panel_status_path,
            report_path=report_path,
            screening_runner=fake_screening_runner,
        )

    assert call_order == [
        ("Dy", 101),
        ("Dy", 102),
        ("Dy", 103),
        ("Nd", 101),
        ("Nd", 102),
        ("Nd", 103),
    ]
    with replicate_summary_path.open("r", encoding="utf-8", newline="") as handle:
        replicate_rows = list(csv.DictReader(handle))
    assert len(replicate_rows) == 6
    failed_row = next(row for row in replicate_rows if row["success_status"] == "failure")
    assert failed_row["panel_member_id"] == "candidate_a"
    assert failed_row["target_metal"] == "Dy"
    assert failed_row["replicate_id"] == "2"
    assert failed_row["failure_reason"] == "RuntimeError: synthetic screening failure"

    with system_summary_path.open("r", encoding="utf-8", newline="") as handle:
        system_rows = list(csv.DictReader(handle))
    assert len(system_rows) == 2
    candidate_a_row = next(row for row in system_rows if row["panel_member_id"] == "candidate_a")
    candidate_b_row = next(row for row in system_rows if row["panel_member_id"] == "candidate_b")
    assert candidate_a_row["screening_status"] == "failed"
    assert candidate_b_row["screening_status"] == "stable_bound"

    with panel_status_path.open("r", encoding="utf-8", newline="") as handle:
        panel_rows = list(csv.DictReader(handle))
    assert len(panel_rows) == 2
    assert report_path.exists()
    assert "synthetic screening failure" in report_path.read_text(encoding="utf-8")
    failed_log_path = tmp_path / "results" / "openmm_screening" / "candidate_a" / "Dy" / "replicate_2" / "screening_log.txt"
    assert failed_log_path.exists()
    assert "RuntimeError: synthetic screening failure" in failed_log_path.read_text(encoding="utf-8")
