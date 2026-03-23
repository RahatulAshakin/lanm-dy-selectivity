from __future__ import annotations

import csv
from pathlib import Path

import yaml

from lanm.analysis.round2_openmm_screening import (
    Round2BuildContext,
    Round2OpenMMScreeningReplicateSummaryRow,
    build_round2_openmm_panel_status_rows,
    build_round2_openmm_replicate_tasks,
    build_round2_openmm_screening_system_summaries,
    build_round2_openmm_systems,
    build_round2_simulation_config,
    discover_round2_openmm_screening_candidates,
    prepare_round2_md_inputs,
    resolve_round2_build_context,
    run_round2_openmm_screening,
)
from lanm.md.openmm_equilibration_smoke import EquilibrationLogRow, OpenMMEquilibrationSmokeArtifacts
from lanm.md.openmm_screening import OpenMMScreeningArtifacts
from lanm.md.openmm_system_build import OpenMMSerializedSystem
from lanm.paths import MD_PROTOCOL_CONFIG_PATH, METAL_PARAMETER_VALUES_PATH


def _write_test_packed_pdb(path: Path, *, metal_label: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    metal_line = (
        f"HETATM    4 {metal_label:>4} {metal_label:>3} A 201       "
        f"2.500   0.000   0.000  1.00  0.00          {metal_label:>2}  "
    )
    path.write_text(
        "\n".join(
            (
                "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N  ",
                "ATOM      2  CA  GLY A   1       1.200   0.000   0.000  1.00  0.00           C  ",
                "ATOM      3  O   GLY A   1       1.900   0.000   0.000  1.00  0.00           O  ",
                metal_line,
                "END",
                "",
            )
        ),
        encoding="utf-8",
    )


def _write_round2_panel(path: Path, packed_paths: dict[str, Path]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "panel_rank",
                "panel_member_id",
                "panel_role",
                "candidate_id",
                "campaign_id",
                "backbone_id",
                "topology_class",
                "preserved_metal_identity",
                "rosetta_score_rank_within_topology",
                "topology_candidate_count",
                "representative_design_id",
                "representative_sequence_id",
                "representative_ligand_confidence",
                "representative_overall_confidence",
                "total_score",
                "representative_packed_pdb",
                "selection_reason",
            ),
        )
        writer.writeheader()
        writer.writerows(
            (
                {
                    "panel_rank": 1,
                    "panel_member_id": "am1_mex_ss_only_u02_r2u02",
                    "panel_role": "top_am1_mex_round2_candidate",
                    "candidate_id": "am1_mex_ss_only_u02_r2u02",
                    "campaign_id": "am1_mex_ss_only",
                    "backbone_id": "am1_mex_8fns_chain_a",
                    "topology_class": "am1_monomer",
                    "preserved_metal_identity": "ND",
                    "rosetta_score_rank_within_topology": 1,
                    "topology_candidate_count": 1,
                    "representative_design_id": 1,
                    "representative_sequence_id": "am1_mex_ss_only_u02_r2u02_design_01",
                    "representative_ligand_confidence": 0.2213,
                    "representative_overall_confidence": 0.3023,
                    "total_score": 218.234,
                    "representative_packed_pdb": str(packed_paths["am1_mex_ss_only_u02_r2u02"]),
                    "selection_reason": "synthetic am1 finalist",
                },
                {
                    "panel_rank": 2,
                    "panel_member_id": "hans_pocket_ss_only_u02_r2u01",
                    "panel_role": "top_hans_pocket_round2_candidate",
                    "candidate_id": "hans_pocket_ss_only_u02_r2u01",
                    "campaign_id": "hans_pocket_ss_only",
                    "backbone_id": "hans_pocket_8fnr_chain_a",
                    "topology_class": "hans_monomer",
                    "preserved_metal_identity": "DY",
                    "rosetta_score_rank_within_topology": 1,
                    "topology_candidate_count": 1,
                    "representative_design_id": 1,
                    "representative_sequence_id": "hans_pocket_ss_only_u02_r2u01_design_01",
                    "representative_ligand_confidence": 1.0,
                    "representative_overall_confidence": 0.3399,
                    "total_score": 144.070,
                    "representative_packed_pdb": str(packed_paths["hans_pocket_ss_only_u02_r2u01"]),
                    "selection_reason": "synthetic hans pocket finalist",
                },
                {
                    "panel_rank": 3,
                    "panel_member_id": "hans_interface_ss_plus_if_u04_r2u01",
                    "panel_role": "best_hans_interface_aware_round2_candidate",
                    "candidate_id": "hans_interface_ss_plus_if_u04_r2u01",
                    "campaign_id": "hans_interface_ss_plus_if",
                    "backbone_id": "hans_interface_8fnr_a_b_c_d",
                    "topology_class": "hans_interface_multichain",
                    "preserved_metal_identity": "DY",
                    "rosetta_score_rank_within_topology": 1,
                    "topology_candidate_count": 2,
                    "representative_design_id": 2,
                    "representative_sequence_id": "hans_interface_ss_plus_if_u04_r2u01_design_02",
                    "representative_ligand_confidence": 0.3454,
                    "representative_overall_confidence": 0.3303,
                    "total_score": 530.712,
                    "representative_packed_pdb": str(packed_paths["hans_interface_ss_plus_if_u04_r2u01"]),
                    "selection_reason": "synthetic hans interface finalist",
                },
            )
        )


def _build_test_panel(tmp_path: Path) -> Path:
    packed_root = tmp_path / "packed"
    packed_paths = {
        "am1_mex_ss_only_u02_r2u02": packed_root / "am1_mex_ss_only_u02_r2u02.pdb",
        "hans_pocket_ss_only_u02_r2u01": packed_root / "hans_pocket_ss_only_u02_r2u01.pdb",
        "hans_interface_ss_plus_if_u04_r2u01": packed_root / "hans_interface_ss_plus_if_u04_r2u01.pdb",
    }
    _write_test_packed_pdb(packed_paths["am1_mex_ss_only_u02_r2u02"], metal_label="ND")
    _write_test_packed_pdb(packed_paths["hans_pocket_ss_only_u02_r2u01"], metal_label="DY")
    _write_test_packed_pdb(packed_paths["hans_interface_ss_plus_if_u04_r2u01"], metal_label="DY")
    panel_path = tmp_path / "results" / "tables" / "round2_md_rescreen_panel.csv"
    _write_round2_panel(panel_path, packed_paths)
    return panel_path


def _round2_build_context() -> Round2BuildContext:
    return resolve_round2_build_context(
        md_protocol_path=MD_PROTOCOL_CONFIG_PATH,
        metal_parameter_values_path=METAL_PARAMETER_VALUES_PATH,
    )


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


def test_discover_round2_openmm_screening_candidates_reads_reduced_panel(tmp_path: Path) -> None:
    panel_path = _build_test_panel(tmp_path)

    candidates = discover_round2_openmm_screening_candidates(panel_path)

    assert [candidate.candidate_id for candidate in candidates] == [
        "am1_mex_ss_only_u02_r2u02",
        "hans_pocket_ss_only_u02_r2u01",
        "hans_interface_ss_plus_if_u04_r2u01",
    ]
    assert candidates[0].preserved_metal_identity == "ND"
    assert candidates[1].topology_class == "hans_monomer"
    assert candidates[2].representative_packed_pdb_path.exists()


def test_build_round2_openmm_systems_and_replicate_tasks_cover_reduced_panel(tmp_path: Path) -> None:
    panel_path = _build_test_panel(tmp_path)
    candidates = discover_round2_openmm_screening_candidates(panel_path)
    systems = build_round2_openmm_systems(
        candidates=candidates,
        output_root=tmp_path / "results" / "round2_openmm_screening",
    )

    prepare_round2_md_inputs(candidates=candidates, systems=systems)
    build_context = _round2_build_context()
    for system in systems:
        simulation_config = build_round2_simulation_config(
            system=system,
            build_context=build_context,
            md_protocol_path=MD_PROTOCOL_CONFIG_PATH,
        )
        system.simulation_config_path.parent.mkdir(parents=True, exist_ok=True)
        with system.simulation_config_path.open("w", encoding="utf-8") as handle:
            yaml.safe_dump(
                {
                    "panel_member_id": simulation_config.panel_member_id,
                    "target_metal": simulation_config.target_metal,
                    "topology_class": simulation_config.topology_class,
                    "replicate_count": simulation_config.replicate_count,
                    "target_temperature_K": simulation_config.target_temperature_K,
                    "friction_coeff_ps": simulation_config.friction_coeff_ps,
                    "timestep_fs": simulation_config.timestep_fs,
                    "nonbonded_method": simulation_config.nonbonded_method,
                    "cutoff_nm": simulation_config.cutoff_nm,
                    "hydrogen_mass_repartitioning": simulation_config.hydrogen_mass_repartitioning,
                    "protein_forcefield": simulation_config.protein_forcefield,
                    "water_model": simulation_config.water_model,
                    "metal_parameter_family": simulation_config.metal_parameter_family,
                    "metal_parameter_provenance": simulation_config.metal_parameter_provenance,
                },
                handle,
                sort_keys=False,
            )

    tasks = build_round2_openmm_replicate_tasks(systems)

    assert len(systems) == 15
    assert len(tasks) == 45
    assert systems[0].candidate_id == "am1_mex_ss_only_u02_r2u02"
    assert systems[0].target_metal == "Dy"
    assert systems[-1].candidate_id == "hans_interface_ss_plus_if_u04_r2u01"
    assert systems[-1].target_metal == "Fe"
    assert " DY " in systems[0].starting_structure_path.read_text(encoding="utf-8")
    assert " ND " in systems[1].starting_structure_path.read_text(encoding="utf-8")
    assert [task.replicate_spec.seed for task in tasks[:3]] == [101, 102, 103]


def test_build_round2_openmm_summary_rows_apply_phase7g1_rules(tmp_path: Path) -> None:
    panel_path = _build_test_panel(tmp_path)
    candidates = discover_round2_openmm_screening_candidates(panel_path)
    systems = build_round2_openmm_systems(
        candidates=candidates,
        output_root=tmp_path / "results" / "round2_openmm_screening",
    )
    am1_dy = next(system for system in systems if system.candidate_id == "am1_mex_ss_only_u02_r2u02" and system.target_metal == "Dy")
    am1_al = next(system for system in systems if system.candidate_id == "am1_mex_ss_only_u02_r2u02" and system.target_metal == "Al")
    pocket_fe = next(system for system in systems if system.candidate_id == "hans_pocket_ss_only_u02_r2u01" and system.target_metal == "Fe")
    interface_y = next(system for system in systems if system.candidate_id == "hans_interface_ss_plus_if_u04_r2u01" and system.target_metal == "Y")

    replicate_rows = (
        Round2OpenMMScreeningReplicateSummaryRow(
            candidate_id=am1_dy.candidate_id,
            target_metal="Dy",
            topology_class=am1_dy.topology_class,
            replicate_id=1,
            seed=101,
            success_status="success",
            failure_reason="",
            mean_temperature_K=298.1,
            final_temperature_K=298.2,
            mean_box_volume_nm3=110.0,
            mean_min_metal_oxygen_distance_A=2.6,
            inner_sphere_occupancy_fraction=0.8,
            state_data_path="a",
            final_state_path="b",
            final_structure_path="c",
            screening_log_path="d",
            trajectory_path="e",
        ),
        Round2OpenMMScreeningReplicateSummaryRow(
            candidate_id=am1_dy.candidate_id,
            target_metal="Dy",
            topology_class=am1_dy.topology_class,
            replicate_id=2,
            seed=102,
            success_status="success",
            failure_reason="",
            mean_temperature_K=298.0,
            final_temperature_K=298.1,
            mean_box_volume_nm3=110.1,
            mean_min_metal_oxygen_distance_A=2.7,
            inner_sphere_occupancy_fraction=0.7,
            state_data_path="a",
            final_state_path="b",
            final_structure_path="c",
            screening_log_path="d",
            trajectory_path="e",
        ),
        Round2OpenMMScreeningReplicateSummaryRow(
            candidate_id=am1_dy.candidate_id,
            target_metal="Dy",
            topology_class=am1_dy.topology_class,
            replicate_id=3,
            seed=103,
            success_status="success",
            failure_reason="",
            mean_temperature_K=297.9,
            final_temperature_K=298.0,
            mean_box_volume_nm3=110.2,
            mean_min_metal_oxygen_distance_A=2.8,
            inner_sphere_occupancy_fraction=0.6,
            state_data_path="a",
            final_state_path="b",
            final_structure_path="c",
            screening_log_path="d",
            trajectory_path="e",
        ),
        Round2OpenMMScreeningReplicateSummaryRow(
            candidate_id=am1_al.candidate_id,
            target_metal="Al",
            topology_class=am1_al.topology_class,
            replicate_id=1,
            seed=101,
            success_status="success",
            failure_reason="",
            mean_temperature_K=298.1,
            final_temperature_K=298.2,
            mean_box_volume_nm3=109.0,
            mean_min_metal_oxygen_distance_A=4.1,
            inner_sphere_occupancy_fraction=0.3,
            state_data_path="a",
            final_state_path="b",
            final_structure_path="c",
            screening_log_path="d",
            trajectory_path="e",
        ),
        Round2OpenMMScreeningReplicateSummaryRow(
            candidate_id=am1_al.candidate_id,
            target_metal="Al",
            topology_class=am1_al.topology_class,
            replicate_id=2,
            seed=102,
            success_status="success",
            failure_reason="",
            mean_temperature_K=298.0,
            final_temperature_K=298.1,
            mean_box_volume_nm3=109.1,
            mean_min_metal_oxygen_distance_A=4.2,
            inner_sphere_occupancy_fraction=0.2,
            state_data_path="a",
            final_state_path="b",
            final_structure_path="c",
            screening_log_path="d",
            trajectory_path="e",
        ),
        Round2OpenMMScreeningReplicateSummaryRow(
            candidate_id=am1_al.candidate_id,
            target_metal="Al",
            topology_class=am1_al.topology_class,
            replicate_id=3,
            seed=103,
            success_status="success",
            failure_reason="",
            mean_temperature_K=297.9,
            final_temperature_K=298.0,
            mean_box_volume_nm3=109.2,
            mean_min_metal_oxygen_distance_A=4.3,
            inner_sphere_occupancy_fraction=0.1,
            state_data_path="a",
            final_state_path="b",
            final_structure_path="c",
            screening_log_path="d",
            trajectory_path="e",
        ),
        Round2OpenMMScreeningReplicateSummaryRow(
            candidate_id=pocket_fe.candidate_id,
            target_metal="Fe",
            topology_class=pocket_fe.topology_class,
            replicate_id=1,
            seed=101,
            success_status="success",
            failure_reason="",
            mean_temperature_K=298.1,
            final_temperature_K=298.2,
            mean_box_volume_nm3=108.0,
            mean_min_metal_oxygen_distance_A=2.0,
            inner_sphere_occupancy_fraction=0.7,
            state_data_path="a",
            final_state_path="b",
            final_structure_path="c",
            screening_log_path="d",
            trajectory_path="e",
        ),
        Round2OpenMMScreeningReplicateSummaryRow(
            candidate_id=pocket_fe.candidate_id,
            target_metal="Fe",
            topology_class=pocket_fe.topology_class,
            replicate_id=2,
            seed=102,
            success_status="success",
            failure_reason="",
            mean_temperature_K=298.0,
            final_temperature_K=298.1,
            mean_box_volume_nm3=108.1,
            mean_min_metal_oxygen_distance_A=2.1,
            inner_sphere_occupancy_fraction=0.8,
            state_data_path="a",
            final_state_path="b",
            final_structure_path="c",
            screening_log_path="d",
            trajectory_path="e",
        ),
        Round2OpenMMScreeningReplicateSummaryRow(
            candidate_id=pocket_fe.candidate_id,
            target_metal="Fe",
            topology_class=pocket_fe.topology_class,
            replicate_id=3,
            seed=103,
            success_status="success",
            failure_reason="",
            mean_temperature_K=297.9,
            final_temperature_K=298.0,
            mean_box_volume_nm3=108.2,
            mean_min_metal_oxygen_distance_A=2.2,
            inner_sphere_occupancy_fraction=0.9,
            state_data_path="a",
            final_state_path="b",
            final_structure_path="c",
            screening_log_path="d",
            trajectory_path="e",
        ),
        Round2OpenMMScreeningReplicateSummaryRow(
            candidate_id=interface_y.candidate_id,
            target_metal="Y",
            topology_class=interface_y.topology_class,
            replicate_id=1,
            seed=101,
            success_status="success",
            failure_reason="",
            mean_temperature_K=298.1,
            final_temperature_K=298.2,
            mean_box_volume_nm3=112.0,
            mean_min_metal_oxygen_distance_A=3.8,
            inner_sphere_occupancy_fraction=0.4,
            state_data_path="a",
            final_state_path="b",
            final_structure_path="c",
            screening_log_path="d",
            trajectory_path="e",
        ),
        Round2OpenMMScreeningReplicateSummaryRow(
            candidate_id=interface_y.candidate_id,
            target_metal="Y",
            topology_class=interface_y.topology_class,
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
        Round2OpenMMScreeningReplicateSummaryRow(
            candidate_id=interface_y.candidate_id,
            target_metal="Y",
            topology_class=interface_y.topology_class,
            replicate_id=3,
            seed=103,
            success_status="success",
            failure_reason="",
            mean_temperature_K=297.9,
            final_temperature_K=298.0,
            mean_box_volume_nm3=112.2,
            mean_min_metal_oxygen_distance_A=3.9,
            inner_sphere_occupancy_fraction=0.3,
            state_data_path="a",
            final_state_path="b",
            final_structure_path="c",
            screening_log_path="d",
            trajectory_path="e",
        ),
    )
    stage_failure_reasons = {
        ("hans_pocket_ss_only_u02_r2u01", "Dy"): "RuntimeError: synthetic build failure",
    }

    system_rows = build_round2_openmm_screening_system_summaries(
        systems=systems,
        replicate_rows=replicate_rows,
        stage_failure_reasons=stage_failure_reasons,
    )
    panel_rows = build_round2_openmm_panel_status_rows(system_rows)

    am1_dy_row = next(row for row in system_rows if row.candidate_id == am1_dy.candidate_id and row.target_metal == "Dy")
    am1_al_row = next(row for row in system_rows if row.candidate_id == am1_al.candidate_id and row.target_metal == "Al")
    pocket_dy_row = next(row for row in system_rows if row.candidate_id == "hans_pocket_ss_only_u02_r2u01" and row.target_metal == "Dy")
    pocket_fe_row = next(row for row in system_rows if row.candidate_id == pocket_fe.candidate_id and row.target_metal == "Fe")
    interface_y_row = next(row for row in system_rows if row.candidate_id == interface_y.candidate_id and row.target_metal == "Y")
    am1_panel_row = next(row for row in panel_rows if row.candidate_id == "am1_mex_ss_only_u02_r2u02")

    assert am1_dy_row.screening_status == "stable_bound"
    assert am1_dy_row.median_inner_sphere_occupancy_fraction == 0.7
    assert am1_al_row.screening_status == "no_persistent_capture"
    assert pocket_dy_row.screening_status == "failed"
    assert pocket_dy_row.failure_reason == "RuntimeError: synthetic build failure"
    assert pocket_fe_row.screening_status == "persistent_capture_flag"
    assert interface_y_row.screening_status == "failed"
    assert "replicate_2: RuntimeError: synthetic failure" in interface_y_row.failure_reason
    assert am1_panel_row.dy_status == "stable_bound"
    assert am1_panel_row.al_status == "no_persistent_capture"
    assert am1_panel_row.candidate_keep_for_metadynamics is False


def test_run_round2_openmm_screening_continues_after_system_failure_and_writes_outputs(tmp_path: Path) -> None:
    panel_path = _build_test_panel(tmp_path)
    summary_path = tmp_path / "results" / "tables" / "round2_openmm_screening_summary.csv"
    panel_status_path = tmp_path / "results" / "tables" / "round2_openmm_screening_panel_status.csv"
    report_path = tmp_path / "results" / "reports" / "round2_openmm_screening.md"

    build_calls: list[tuple[str, str]] = []
    screening_calls: list[tuple[str, str, int]] = []

    def fake_builder(
        *,
        prepared_structure_pdb_text: str,
        forcefield_files: tuple[str, ...],
        target_temperature_K: int,
        friction_coeff_ps: float,
        timestep_fs: float,
        nonbonded_method: str,
        cutoff_nm: float,
        hydrogen_mass_repartitioning: bool,
    ) -> OpenMMSerializedSystem:
        assert forcefield_files == ("amber19-all.xml", "amber19/opc3.xml")
        candidate_id = "synthetic"
        target_metal = "unknown"
        for line in prepared_structure_pdb_text.splitlines():
            if line.startswith("HETATM"):
                target_metal = line[76:78].strip() or line[12:16].strip()
                break
        build_calls.append((candidate_id, target_metal))
        if len(build_calls) == 1:
            raise RuntimeError("synthetic build failure")
        return OpenMMSerializedSystem(
            prepared_structure_pdb_text=prepared_structure_pdb_text,
            system_xml="<System/>",
            integrator_xml="<Integrator/>",
            atom_count=10,
            residue_count=2,
        )

    def fake_equilibrator(**kwargs: object) -> OpenMMEquilibrationSmokeArtifacts:
        return OpenMMEquilibrationSmokeArtifacts(
            minimized_structure_pdb_text="MIN\n",
            nvt_final_pdb_text="NVT\n",
            npt_final_pdb_text="NPT\n",
            final_state_xml="<State/>",
            log_rows=(
                EquilibrationLogRow(
                    stage_name="minimization",
                    status="success",
                    attempt_index=1,
                    minimization_max_iterations=200,
                    step_count=0,
                    duration_ps=0.0,
                    potential_energy_kj_per_mol=-10.0,
                    temperature_K=None,
                    box_volume_nm3=None,
                    notes="synthetic",
                ),
            ),
            minimized_potential_energy=-10.0,
            final_nvt_potential_energy=-9.0,
            final_npt_potential_energy=-8.0,
            final_temperature_K=298.0,
            final_box_volume_nm3=110.0,
            platform_name="CPU",
            platform_properties="Threads=1",
            minimization_max_iterations_used=200,
        )

    def fake_screening_runner(**kwargs: object) -> OpenMMScreeningArtifacts:
        target_metal = str(kwargs["target_metal"])
        seed = int(kwargs["seed"])
        state_data_path = Path(kwargs["state_data_path"])
        final_state_path = Path(kwargs["final_state_path"])
        final_structure_path = Path(kwargs["final_structure_path"])
        trajectory_path = Path(kwargs["trajectory_path"])
        screening_calls.append((state_data_path.parents[1].name, target_metal, seed))
        _write_fake_state_data_csv(
            state_data_path,
            temperature_offset=float(seed - 100),
            volume_offset=0.5 if target_metal in {"Dy", "Nd", "Y"} else 0.0,
        )
        final_state_path.write_text("<State/>", encoding="utf-8")
        final_structure_path.write_text("ATOM      1  O   HOH A   1       0.0   0.0   0.0\nEND\n", encoding="utf-8")
        trajectory_path.write_bytes(b"DCD")
        return OpenMMScreeningArtifacts(
            platform_name="CPU",
            platform_properties="Threads=1",
            metal_atom_count=1,
            oxygen_atom_count=4,
            frame_count=2,
            mean_min_metal_oxygen_distance_A=2.6 if target_metal in {"Dy", "Nd", "Y"} else 4.1,
            inner_sphere_occupancy_fraction=0.75 if target_metal in {"Dy", "Nd", "Y"} else 0.25,
        )

    rows = run_round2_openmm_screening(
        panel_path=panel_path,
        md_protocol_path=MD_PROTOCOL_CONFIG_PATH,
        metal_parameter_values_path=METAL_PARAMETER_VALUES_PATH,
        output_root=tmp_path / "results" / "round2_openmm_screening",
        summary_path=summary_path,
        panel_status_path=panel_status_path,
        report_path=report_path,
        system_builder=fake_builder,
        equilibrator=fake_equilibrator,
        screening_runner=fake_screening_runner,
    )

    assert len(rows) == 15
    assert summary_path.exists()
    assert panel_status_path.exists()
    assert report_path.exists()
    assert len(screening_calls) == 42

    with summary_path.open("r", encoding="utf-8", newline="") as handle:
        summary_rows = list(csv.DictReader(handle))
    failed_row = next(
        row
        for row in summary_rows
        if row["candidate_id"] == "am1_mex_ss_only_u02_r2u02" and row["target_metal"] == "Dy"
    )
    success_row = next(
        row
        for row in summary_rows
        if row["candidate_id"] == "hans_pocket_ss_only_u02_r2u01" and row["target_metal"] == "Dy"
    )
    assert failed_row["screening_status"] == "failed"
    assert failed_row["failure_reason"] == "RuntimeError: synthetic build failure"
    assert success_row["screening_status"] == "stable_bound"

    with panel_status_path.open("r", encoding="utf-8", newline="") as handle:
        panel_rows = list(csv.DictReader(handle))
    pocket_panel_row = next(row for row in panel_rows if row["candidate_id"] == "hans_pocket_ss_only_u02_r2u01")
    assert pocket_panel_row["candidate_keep_for_metadynamics"] == "True"
    assert "synthetic build failure" in report_path.read_text(encoding="utf-8")
