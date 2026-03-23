from __future__ import annotations

import csv
import sys
import types
from pathlib import Path

import pytest
import yaml

import lanm.analysis.openmm_metadynamics_smoke as metadynamics_analysis
import lanm.md.openmm_metadynamics_smoke as metadynamics_md
from lanm.analysis.openmm_metadynamics_smoke import (
    OpenMMMetadynamicsSummaryRow,
    build_openmm_metadynamics_panel_status_rows,
    discover_reduced_openmm_metadynamics_systems,
    run_openmm_metadynamics_smoke,
)
from lanm.md.openmm_metadynamics_smoke import (
    OpenMMMetadynamicsArtifacts,
    build_metadynamics_collective_variable_setup,
    discover_local_pocket_cv_atoms,
)
from lanm.paths import (
    OPENMM_SCREENING_PANEL_STATUS_PATH,
    OPENMM_SYSTEM_BUILD_SUMMARY_PATH,
)


def _write_completed_metadynamics_outputs(
    output_dir: Path,
    *,
    coordination_rows: tuple[tuple[int, float, float, float, str], ...] = (
        (0, 0.0, 8.0, 2.5, "FALSE"),
    ),
    log_text: str = "Phase 6C1 OpenMM metadynamics smoke\n",
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "cv_timeseries.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "step",
                "time_ps",
                "coordination_number",
                "mean_metal_oxygen_distance_A",
                "bias_potential_kj_per_mol",
                "escape_condition_met",
            ),
        )
        writer.writeheader()
        for step, time_ps, coordination_number, mean_distance_A, escape_condition_met in coordination_rows:
            writer.writerow(
                {
                    "step": step,
                    "time_ps": time_ps,
                    "coordination_number": coordination_number,
                    "mean_metal_oxygen_distance_A": mean_distance_A,
                    "bias_potential_kj_per_mol": 0.0,
                    "escape_condition_met": escape_condition_met,
                }
            )
    (output_dir / "bias_state.xml").write_text("<bias/>", encoding="utf-8")
    (output_dir / "metadynamics_final_state.xml").write_text("<State/>", encoding="utf-8")
    (output_dir / "metadynamics_final_structure.pdb").write_text("END\n", encoding="utf-8")
    (output_dir / "metadynamics_log.txt").write_text(log_text, encoding="utf-8")


def _build_fake_metadynamics_artifacts(target_metal: str) -> OpenMMMetadynamicsArtifacts:
    escape_event_detected = target_metal in {"Al", "Fe"}
    return OpenMMMetadynamicsArtifacts(
        cv_timeseries_rows=(),
        bias_state_xml="<bias/>",
        final_state_xml="<State/>",
        final_structure_pdb_text="END\n",
        platform_name="CUDA",
        platform_properties="Precision=mixed, DeviceIndex=0, DeterministicForces=true, UseCpuPme=false",
        metal_atom_count=1,
        pocket_oxygen_atom_count=4,
        minimum_coordination_number=3.5 if escape_event_detected else 8.0,
        maximum_mean_metal_oxygen_distance_A=4.2 if escape_event_detected else 2.7,
        escape_event_detected=escape_event_detected,
        final_coordination_number=3.1 if escape_event_detected else 7.9,
        final_mean_metal_oxygen_distance_A=4.0 if escape_event_detected else 2.6,
    )


def _write_fake_phase_6c1_system(
    tmp_path: Path,
    *,
    panel_member_id: str,
    target_metal: str,
    topology_class: str = "am1_monomer",
) -> dict[str, str]:
    build_dir = tmp_path / "results" / "openmm_system_build" / panel_member_id / target_metal
    build_dir.mkdir(parents=True, exist_ok=True)
    (build_dir / "system.xml").write_text("<System/>", encoding="utf-8")
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
        "\n".join(
            [
                "ATOM      1  O   ASP A  35       0.000   0.000   0.000  1.00  0.00           O  ",
                f"HETATM    2 {target_metal.upper():>2}   {target_metal.upper():>2} A 201       2.300   0.000   0.000  1.00  0.00          {target_metal:>2}  ",
                "END",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return {
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


def _write_screening_panel_status(path: Path, panel_member_ids: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "panel_member_id",
                "topology_class",
                "dy_status",
                "nd_status",
                "y_status",
                "al_status",
                "fe_status",
            ),
        )
        writer.writeheader()
        for panel_member_id in panel_member_ids:
            writer.writerow(
                {
                    "panel_member_id": panel_member_id,
                    "topology_class": "am1_monomer",
                    "dy_status": "stable_bound",
                    "nd_status": "stable_bound",
                    "y_status": "stable_bound",
                    "al_status": "persistent_capture_flag",
                    "fe_status": "persistent_capture_flag",
                }
            )


def _patch_small_metadynamics_panel(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    panel_member_ids: tuple[str, ...],
) -> None:
    monkeypatch.setattr(metadynamics_analysis, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(
        metadynamics_analysis,
        "REDUCED_METADYNAMICS_PANEL_MEMBER_IDS",
        panel_member_ids,
    )
    monkeypatch.setattr(
        metadynamics_analysis,
        "EXPECTED_METADYNAMICS_SYSTEM_COUNT",
        len(panel_member_ids) * 5,
    )


def test_discover_reduced_openmm_metadynamics_systems_matches_target_panel() -> None:
    systems = discover_reduced_openmm_metadynamics_systems(
        screening_panel_status_path=OPENMM_SCREENING_PANEL_STATUS_PATH,
        system_build_summary_path=OPENMM_SYSTEM_BUILD_SUMMARY_PATH,
    )

    assert len(systems) == 25
    assert systems[0].panel_member_id == "am1_mex_ss_only_u02"
    assert systems[0].target_metal == "Dy"
    assert systems[-1].panel_member_id == "hans_interface_wt_reference"
    assert systems[-1].target_metal == "Fe"
    assert all(system.system_xml_path.exists() for system in systems)
    assert all(system.equilibrated_structure_path.exists() for system in systems)
    assert not any(system.panel_member_id == "hans_interface_ss_plus_if_u02" for system in systems)


def test_build_metadynamics_collective_variable_setup_uses_local_pocket_oxygens() -> None:
    pdb_text = "\n".join(
        [
            "HETATM    1 DY    DY A 201       0.000   0.000   0.000  1.00  0.00          DY  ",
            "ATOM      2  O   ASP A  35       2.400   0.000   0.000  1.00  0.00           O  ",
            "ATOM      3  O   HOH A 301       3.100   0.000   0.000  1.00  0.00           O  ",
            "ATOM      4  O   HOH A 302       6.500   0.000   0.000  1.00  0.00           O  ",
            "HETATM    5 DY    DY A 202      10.000   0.000   0.000  1.00  0.00          DY  ",
            "ATOM      6  O   GLU A  46      12.300   0.000   0.000  1.00  0.00           O  ",
            "END",
            "",
        ]
    )

    atom_selection = discover_local_pocket_cv_atoms(
        equilibrated_structure_pdb_text=pdb_text,
        target_metal="Dy",
    )
    cv_setup = build_metadynamics_collective_variable_setup(atom_selection=atom_selection)

    assert atom_selection.metal_atom_indices == (0, 4)
    assert atom_selection.pocket_oxygen_atom_indices == (1, 2, 5)
    assert cv_setup.coordination_number_descriptor.name == "coordination_number"
    assert cv_setup.coordination_number_descriptor.minimum_value == 0.0
    assert cv_setup.coordination_number_descriptor.maximum_value == 12.0
    assert cv_setup.coordination_number_descriptor.grid_width > 0
    assert cv_setup.mean_metal_oxygen_distance_descriptor.name == "mean_metal_oxygen_distance_A"
    assert cv_setup.mean_metal_oxygen_distance_descriptor.minimum_value == 1.8
    assert cv_setup.mean_metal_oxygen_distance_descriptor.maximum_value == 6.5
    assert cv_setup.mean_metal_oxygen_distance_descriptor.grid_width > 0


def test_build_openmm_metadynamics_panel_status_rows_applies_keep_rule() -> None:
    system_rows = (
        OpenMMMetadynamicsSummaryRow(
            panel_member_id="candidate_keep",
            target_metal="Dy",
            success_status="success",
            failure_reason="",
            minimum_coordination_number=7.9,
            maximum_mean_metal_oxygen_distance_A=2.8,
            escape_event_detected="FALSE",
            final_coordination_number=8.1,
            final_mean_metal_oxygen_distance_A=2.5,
            metadynamics_status="retained_bound",
        ),
        OpenMMMetadynamicsSummaryRow(
            panel_member_id="candidate_keep",
            target_metal="Nd",
            success_status="success",
            failure_reason="",
            minimum_coordination_number=7.8,
            maximum_mean_metal_oxygen_distance_A=2.9,
            escape_event_detected="FALSE",
            final_coordination_number=8.0,
            final_mean_metal_oxygen_distance_A=2.6,
            metadynamics_status="retained_bound",
        ),
        OpenMMMetadynamicsSummaryRow(
            panel_member_id="candidate_keep",
            target_metal="Y",
            success_status="success",
            failure_reason="",
            minimum_coordination_number=6.9,
            maximum_mean_metal_oxygen_distance_A=3.5,
            escape_event_detected="TRUE",
            final_coordination_number=5.5,
            final_mean_metal_oxygen_distance_A=3.6,
            metadynamics_status="escape_or_weak_binding",
        ),
        OpenMMMetadynamicsSummaryRow(
            panel_member_id="candidate_keep",
            target_metal="Al",
            success_status="success",
            failure_reason="",
            minimum_coordination_number=3.2,
            maximum_mean_metal_oxygen_distance_A=4.4,
            escape_event_detected="TRUE",
            final_coordination_number=3.0,
            final_mean_metal_oxygen_distance_A=4.1,
            metadynamics_status="escaped_under_bias",
        ),
        OpenMMMetadynamicsSummaryRow(
            panel_member_id="candidate_keep",
            target_metal="Fe",
            success_status="success",
            failure_reason="",
            minimum_coordination_number=3.8,
            maximum_mean_metal_oxygen_distance_A=4.0,
            escape_event_detected="TRUE",
            final_coordination_number=3.5,
            final_mean_metal_oxygen_distance_A=3.9,
            metadynamics_status="escaped_under_bias",
        ),
        OpenMMMetadynamicsSummaryRow(
            panel_member_id="candidate_reject",
            target_metal="Dy",
            success_status="success",
            failure_reason="",
            minimum_coordination_number=7.4,
            maximum_mean_metal_oxygen_distance_A=3.0,
            escape_event_detected="FALSE",
            final_coordination_number=7.8,
            final_mean_metal_oxygen_distance_A=2.7,
            metadynamics_status="retained_bound",
        ),
        OpenMMMetadynamicsSummaryRow(
            panel_member_id="candidate_reject",
            target_metal="Nd",
            success_status="success",
            failure_reason="",
            minimum_coordination_number=8.2,
            maximum_mean_metal_oxygen_distance_A=2.7,
            escape_event_detected="FALSE",
            final_coordination_number=8.4,
            final_mean_metal_oxygen_distance_A=2.4,
            metadynamics_status="retained_bound",
        ),
        OpenMMMetadynamicsSummaryRow(
            panel_member_id="candidate_reject",
            target_metal="Y",
            success_status="success",
            failure_reason="",
            minimum_coordination_number=7.2,
            maximum_mean_metal_oxygen_distance_A=3.1,
            escape_event_detected="FALSE",
            final_coordination_number=7.5,
            final_mean_metal_oxygen_distance_A=2.8,
            metadynamics_status="retained_bound",
        ),
        OpenMMMetadynamicsSummaryRow(
            panel_member_id="candidate_reject",
            target_metal="Al",
            success_status="success",
            failure_reason="",
            minimum_coordination_number=3.0,
            maximum_mean_metal_oxygen_distance_A=4.2,
            escape_event_detected="TRUE",
            final_coordination_number=2.9,
            final_mean_metal_oxygen_distance_A=4.0,
            metadynamics_status="escaped_under_bias",
        ),
        OpenMMMetadynamicsSummaryRow(
            panel_member_id="candidate_reject",
            target_metal="Fe",
            success_status="success",
            failure_reason="",
            minimum_coordination_number=3.1,
            maximum_mean_metal_oxygen_distance_A=4.3,
            escape_event_detected="TRUE",
            final_coordination_number=3.0,
            final_mean_metal_oxygen_distance_A=4.1,
            metadynamics_status="escaped_under_bias",
        ),
    )

    panel_rows = build_openmm_metadynamics_panel_status_rows(system_rows)

    keep_row = next(row for row in panel_rows if row.panel_member_id == "candidate_keep")
    reject_row = next(row for row in panel_rows if row.panel_member_id == "candidate_reject")
    assert keep_row.candidate_keep_for_qm == "TRUE"
    assert reject_row.candidate_keep_for_qm == "FALSE"


def test_run_openmm_metadynamics_smoke_max_workers_1_runs_serially_without_pool(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    panel_member_ids = ("candidate_a",)
    build_rows: list[dict[str, str]] = []
    for target_metal in ("Dy", "Nd", "Y", "Al", "Fe"):
        build_rows.append(
            _write_fake_phase_6c1_system(
                tmp_path,
                panel_member_id="candidate_a",
                target_metal=target_metal,
            )
        )

    screening_panel_status_path = tmp_path / "results" / "tables" / "openmm_screening_panel_status.csv"
    system_build_summary_path = tmp_path / "results" / "tables" / "openmm_system_build_summary.csv"
    system_summary_path = tmp_path / "results" / "tables" / "openmm_metadynamics_summary.csv"
    panel_status_path = tmp_path / "results" / "tables" / "openmm_metadynamics_panel_status.csv"
    report_path = tmp_path / "results" / "reports" / "openmm_metadynamics_smoke.md"
    _write_screening_panel_status(screening_panel_status_path, panel_member_ids)
    _write_build_summary(system_build_summary_path, build_rows)
    _patch_small_metadynamics_panel(monkeypatch, tmp_path, panel_member_ids)

    class FailIfPoolUsed:
        def __init__(self, *args: object, **kwargs: object) -> None:
            raise AssertionError("ProcessPoolExecutor should not be created in serial mode")

    call_order: list[str] = []

    def fake_default_runner(**kwargs: object) -> OpenMMMetadynamicsArtifacts:
        target_metal = str(kwargs["target_metal"])
        call_order.append(target_metal)
        _write_completed_metadynamics_outputs(Path(str(kwargs["cv_timeseries_path"])).parent)
        return _build_fake_metadynamics_artifacts(target_metal)

    monkeypatch.setattr(metadynamics_analysis, "ProcessPoolExecutor", FailIfPoolUsed)
    monkeypatch.setattr(
        metadynamics_analysis,
        "run_openmm_metadynamics_smoke_system",
        fake_default_runner,
    )

    system_rows = run_openmm_metadynamics_smoke(
        screening_panel_status_path=screening_panel_status_path,
        system_build_summary_path=system_build_summary_path,
        output_root=tmp_path / "results" / "openmm_metadynamics_smoke",
        system_summary_path=system_summary_path,
        panel_status_path=panel_status_path,
        report_path=report_path,
        metadynamics_runner=metadynamics_analysis.run_openmm_metadynamics_smoke_system,
        max_workers=1,
    )

    assert [row.target_metal for row in system_rows] == ["Dy", "Nd", "Y", "Al", "Fe"]
    assert call_order == ["Dy", "Nd", "Y", "Al", "Fe"]
    assert capsys.readouterr().out.splitlines() == [
        "[openmm_metadynamics_smoke] starting system candidate_a/Dy",
        "[openmm_metadynamics_smoke] finished system candidate_a/Dy",
        "[openmm_metadynamics_smoke] starting system candidate_a/Nd",
        "[openmm_metadynamics_smoke] finished system candidate_a/Nd",
        "[openmm_metadynamics_smoke] starting system candidate_a/Y",
        "[openmm_metadynamics_smoke] finished system candidate_a/Y",
        "[openmm_metadynamics_smoke] starting system candidate_a/Al",
        "[openmm_metadynamics_smoke] finished system candidate_a/Al",
        "[openmm_metadynamics_smoke] starting system candidate_a/Fe",
        "[openmm_metadynamics_smoke] finished system candidate_a/Fe",
    ]


def test_run_openmm_metadynamics_smoke_skips_completed_outputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    panel_member_ids = ("candidate_a",)
    build_rows: list[dict[str, str]] = []
    for target_metal in ("Dy", "Nd", "Y", "Al", "Fe"):
        build_rows.append(
            _write_fake_phase_6c1_system(
                tmp_path,
                panel_member_id="candidate_a",
                target_metal=target_metal,
            )
        )

    screening_panel_status_path = tmp_path / "results" / "tables" / "openmm_screening_panel_status.csv"
    system_build_summary_path = tmp_path / "results" / "tables" / "openmm_system_build_summary.csv"
    system_summary_path = tmp_path / "results" / "tables" / "openmm_metadynamics_summary.csv"
    panel_status_path = tmp_path / "results" / "tables" / "openmm_metadynamics_panel_status.csv"
    report_path = tmp_path / "results" / "reports" / "openmm_metadynamics_smoke.md"
    output_root = tmp_path / "results" / "openmm_metadynamics_smoke"
    _write_screening_panel_status(screening_panel_status_path, panel_member_ids)
    _write_build_summary(system_build_summary_path, build_rows)
    _patch_small_metadynamics_panel(monkeypatch, tmp_path, panel_member_ids)

    existing_output_dir = output_root / "candidate_a" / "Dy"
    existing_log_text = "existing metadynamics log\n"
    _write_completed_metadynamics_outputs(
        existing_output_dir,
        coordination_rows=(
            (0, 0.0, 8.8, 2.4, "FALSE"),
            (500, 1.0, 8.4, 2.6, "FALSE"),
        ),
        log_text=existing_log_text,
    )

    call_order: list[str] = []

    def fake_metadynamics_runner(**kwargs: object) -> OpenMMMetadynamicsArtifacts:
        target_metal = str(kwargs["target_metal"])
        call_order.append(target_metal)
        _write_completed_metadynamics_outputs(Path(str(kwargs["cv_timeseries_path"])).parent)
        return _build_fake_metadynamics_artifacts(target_metal)

    system_rows = run_openmm_metadynamics_smoke(
        screening_panel_status_path=screening_panel_status_path,
        system_build_summary_path=system_build_summary_path,
        output_root=output_root,
        system_summary_path=system_summary_path,
        panel_status_path=panel_status_path,
        report_path=report_path,
        metadynamics_runner=fake_metadynamics_runner,
        max_workers=1,
    )

    dy_row = next(row for row in system_rows if row.target_metal == "Dy")
    assert call_order == ["Nd", "Y", "Al", "Fe"]
    assert dy_row.success_status == "success"
    assert dy_row.minimum_coordination_number == 8.4
    assert dy_row.maximum_mean_metal_oxygen_distance_A == 2.6
    assert dy_row.final_coordination_number == 8.4
    assert dy_row.final_mean_metal_oxygen_distance_A == 2.6
    assert (existing_output_dir / "metadynamics_log.txt").read_text(encoding="utf-8") == existing_log_text
    output_lines = capsys.readouterr().out.splitlines()
    assert "[openmm_metadynamics_smoke] skipping completed system candidate_a/Dy" in output_lines


def test_choose_metadynamics_cuda_platform_uses_requested_cuda_properties(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    requested_names: list[str] = []
    platform_token = object()

    class FakePlatform:
        @staticmethod
        def getPlatformByName(name: str) -> object:
            requested_names.append(name)
            return platform_token

    openmm_module = types.ModuleType("openmm")
    openmm_module.Platform = FakePlatform
    monkeypatch.setitem(sys.modules, "openmm", openmm_module)

    platform, properties, platform_name, properties_display = (
        metadynamics_md._choose_metadynamics_cuda_platform()
    )

    assert platform is platform_token
    assert requested_names == ["CUDA"]
    assert properties == {
        "Precision": "mixed",
        "DeviceIndex": "0",
        "DeterministicForces": "true",
        "UseCpuPme": "false",
    }
    assert platform_name == "CUDA"
    assert (
        properties_display
        == "Precision=mixed, DeviceIndex=0, DeterministicForces=true, UseCpuPme=false"
    )


def test_run_openmm_metadynamics_smoke_continues_after_failure_and_writes_summary(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    panel_member_ids = ("candidate_a", "candidate_b", "candidate_c")
    build_rows: list[dict[str, str]] = []
    for panel_member_id in panel_member_ids:
        for target_metal in ("Dy", "Nd", "Y", "Al", "Fe"):
            build_rows.append(
                _write_fake_phase_6c1_system(
                    tmp_path,
                    panel_member_id=panel_member_id,
                    target_metal=target_metal,
                )
            )

    screening_panel_status_path = tmp_path / "results" / "tables" / "openmm_screening_panel_status.csv"
    system_build_summary_path = tmp_path / "results" / "tables" / "openmm_system_build_summary.csv"
    system_summary_path = tmp_path / "results" / "tables" / "openmm_metadynamics_summary.csv"
    panel_status_path = tmp_path / "results" / "tables" / "openmm_metadynamics_panel_status.csv"
    report_path = tmp_path / "results" / "reports" / "openmm_metadynamics_smoke.md"
    _write_screening_panel_status(screening_panel_status_path, panel_member_ids)
    _write_build_summary(system_build_summary_path, build_rows)
    _patch_small_metadynamics_panel(monkeypatch, tmp_path, panel_member_ids)

    call_count = 0

    def fake_metadynamics_runner(**kwargs: object) -> OpenMMMetadynamicsArtifacts:
        nonlocal call_count
        call_count += 1
        if call_count == 2:
            raise RuntimeError("synthetic metadynamics failure")
        cv_timeseries_path = Path(str(kwargs["cv_timeseries_path"]))
        bias_state_path = Path(str(kwargs["bias_state_path"]))
        final_state_path = Path(str(kwargs["final_state_path"]))
        final_structure_path = Path(str(kwargs["final_structure_path"]))
        cv_timeseries_path.parent.mkdir(parents=True, exist_ok=True)
        with cv_timeseries_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=(
                    "step",
                    "time_ps",
                    "coordination_number",
                    "mean_metal_oxygen_distance_A",
                    "bias_potential_kj_per_mol",
                    "escape_condition_met",
                ),
            )
            writer.writeheader()
            writer.writerow(
                {
                    "step": 0,
                    "time_ps": 0.0,
                    "coordination_number": 8.0,
                    "mean_metal_oxygen_distance_A": 2.5,
                    "bias_potential_kj_per_mol": 0.0,
                    "escape_condition_met": "FALSE",
                }
            )
        bias_state_path.write_text("<bias/>", encoding="utf-8")
        final_state_path.write_text("<State/>", encoding="utf-8")
        final_structure_path.write_text("END\n", encoding="utf-8")
        escape_event_detected = str(kwargs["target_metal"]) in {"Al", "Fe"}
        return OpenMMMetadynamicsArtifacts(
            cv_timeseries_rows=(),
            bias_state_xml="<bias/>",
            final_state_xml="<State/>",
            final_structure_pdb_text="END\n",
            platform_name="CPU",
            platform_properties="Threads=1",
            metal_atom_count=1,
            pocket_oxygen_atom_count=4,
            minimum_coordination_number=3.5 if escape_event_detected else 8.0,
            maximum_mean_metal_oxygen_distance_A=4.2 if escape_event_detected else 2.7,
            escape_event_detected=escape_event_detected,
            final_coordination_number=3.1 if escape_event_detected else 7.9,
            final_mean_metal_oxygen_distance_A=4.0 if escape_event_detected else 2.6,
        )

    with pytest.raises(RuntimeError, match="failed for 1/15 systems"):
        run_openmm_metadynamics_smoke(
            screening_panel_status_path=screening_panel_status_path,
            system_build_summary_path=system_build_summary_path,
            output_root=tmp_path / "results" / "openmm_metadynamics_smoke",
            system_summary_path=system_summary_path,
            panel_status_path=panel_status_path,
            report_path=report_path,
            metadynamics_runner=fake_metadynamics_runner,
        )

    assert call_count == 15
    with system_summary_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 15
    assert rows[1]["success_status"] == "failure"
    assert rows[1]["failure_reason"] == "RuntimeError: synthetic metadynamics failure"
    assert rows[-1]["success_status"] == "success"
    assert rows[-1]["metadynamics_status"] in {"retained_bound", "escaped_under_bias"}
    assert panel_status_path.exists()
    assert report_path.exists()
    assert "synthetic metadynamics failure" in report_path.read_text(encoding="utf-8")
    assert (
        tmp_path / "results" / "openmm_metadynamics_smoke" / "candidate_b" / "Dy" / "metadynamics_log.txt"
    ).exists()
