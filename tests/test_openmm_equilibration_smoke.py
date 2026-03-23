from __future__ import annotations

import csv
from pathlib import Path

import pytest
import yaml

from lanm.analysis.openmm_equilibration_smoke import (
    build_equilibration_summary_row,
    discover_successful_baseline_openmm_systems,
    load_simulation_config,
    run_openmm_equilibration_smoke,
)
from lanm.md.openmm_equilibration_smoke import EquilibrationLogRow, OpenMMEquilibrationSmokeArtifacts
from lanm.paths import OPENMM_SYSTEM_BUILD_SUMMARY_PATH


def _write_fake_phase_6a4_system(
    tmp_path: Path,
    *,
    panel_member_id: str,
    target_metal: str,
    topology_class: str = "am1_monomer",
) -> dict[str, str]:
    output_dir = tmp_path / "results" / "openmm_system_build" / panel_member_id / target_metal
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "system.xml").write_text("<System/>", encoding="utf-8")
    (output_dir / "prepared_structure.pdb").write_text("ATOM      1  N   GLY A   1       0.0   0.0   0.0\nEND\n", encoding="utf-8")
    with (output_dir / "simulation_config.yaml").open("w", encoding="utf-8") as handle:
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
    return {
        "panel_member_id": panel_member_id,
        "target_metal": target_metal,
        "topology_class": topology_class,
        "success_status": "success",
        "failure_reason": "",
        "atom_count": "100",
        "residue_count": "10",
        "system_xml_path": str(output_dir / "system.xml"),
        "simulation_config_path": str(output_dir / "simulation_config.yaml"),
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


def test_discover_successful_baseline_openmm_systems_matches_full_panel() -> None:
    systems = discover_successful_baseline_openmm_systems(OPENMM_SYSTEM_BUILD_SUMMARY_PATH)

    assert len(systems) == 30
    assert systems[0].panel_member_id == "am1_mex_ss_only_u02"
    assert systems[0].target_metal == "Al"
    assert systems[-1].panel_member_id == "hans_interface_ss_plus_if_u02"
    assert systems[-1].target_metal == "Y"
    assert all(system.system_xml_path.exists() for system in systems)
    assert all(system.prepared_structure_path.exists() for system in systems)


def test_load_simulation_config_reads_phase_6a4_defaults() -> None:
    first_system = discover_successful_baseline_openmm_systems(OPENMM_SYSTEM_BUILD_SUMMARY_PATH)[0]

    config = load_simulation_config(first_system.simulation_config_path)

    assert config.panel_member_id == "am1_mex_ss_only_u02"
    assert config.target_metal == "Al"
    assert config.topology_class == "am1_monomer"
    assert config.replicate_count == 3
    assert config.target_temperature_K == 298
    assert config.friction_coeff_ps == 1.0
    assert config.timestep_fs == 2.0
    assert config.nonbonded_method == "NoCutoff"
    assert config.cutoff_nm == 1.0


def test_build_equilibration_summary_row_uses_npt_final_structure_path(tmp_path: Path) -> None:
    summary_row = build_equilibration_summary_row(
        system_row=discover_successful_baseline_openmm_systems(OPENMM_SYSTEM_BUILD_SUMMARY_PATH)[0],
        output_dir=tmp_path / "results" / "openmm_equilibration_smoke" / "am1_mex_ss_only_u02" / "Al",
        artifacts=OpenMMEquilibrationSmokeArtifacts(
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
            final_temperature_K=298.2,
            final_box_volume_nm3=120.0,
            platform_name="CPU",
            platform_properties="Threads=16",
            minimization_max_iterations_used=200,
        ),
        failure_reason="",
    )

    assert summary_row.success_status == "success"
    assert summary_row.final_structure_path.endswith("npt_final.pdb")
    assert summary_row.final_npt_potential_energy == -8.0
    assert summary_row.final_temperature_K == 298.2


def test_run_openmm_equilibration_smoke_continues_after_failure_and_writes_summary(tmp_path: Path) -> None:
    system_one = _write_fake_phase_6a4_system(tmp_path, panel_member_id="candidate_a", target_metal="Dy")
    system_two = _write_fake_phase_6a4_system(tmp_path, panel_member_id="candidate_b", target_metal="Nd")
    system_three = _write_fake_phase_6a4_system(tmp_path, panel_member_id="candidate_c", target_metal="Y")
    summary_path = tmp_path / "results" / "tables" / "openmm_system_build_summary.csv"
    report_path = tmp_path / "results" / "reports" / "openmm_equilibration_smoke.md"
    smoke_summary_path = tmp_path / "results" / "tables" / "openmm_equilibration_smoke_summary.csv"
    _write_build_summary(summary_path, [system_one, system_two, system_three])

    call_count = 0

    def fake_equilibrator(**kwargs: object) -> OpenMMEquilibrationSmokeArtifacts:
        nonlocal call_count
        call_count += 1
        if call_count == 2:
            raise RuntimeError("synthetic equilibration failure")
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
                    potential_energy_kj_per_mol=-10.0 * call_count,
                    temperature_K=None,
                    box_volume_nm3=None,
                    notes="synthetic",
                ),
                EquilibrationLogRow(
                    stage_name="nvt",
                    status="success",
                    attempt_index=1,
                    minimization_max_iterations=200,
                    step_count=10000,
                    duration_ps=20.0,
                    potential_energy_kj_per_mol=-9.0 * call_count,
                    temperature_K=298.0,
                    box_volume_nm3=None,
                    notes="synthetic",
                ),
                EquilibrationLogRow(
                    stage_name="npt",
                    status="success",
                    attempt_index=1,
                    minimization_max_iterations=200,
                    step_count=10000,
                    duration_ps=20.0,
                    potential_energy_kj_per_mol=-8.0 * call_count,
                    temperature_K=297.5,
                    box_volume_nm3=110.0 + call_count,
                    notes="synthetic",
                ),
            ),
            minimized_potential_energy=-10.0 * call_count,
            final_nvt_potential_energy=-9.0 * call_count,
            final_npt_potential_energy=-8.0 * call_count,
            final_temperature_K=297.5,
            final_box_volume_nm3=110.0 + call_count,
            platform_name="CPU",
            platform_properties="Threads=16",
            minimization_max_iterations_used=200,
        )

    with pytest.raises(RuntimeError, match="failed for 1/3 systems"):
        run_openmm_equilibration_smoke(
            system_build_summary_path=summary_path,
            output_root=tmp_path / "results" / "openmm_equilibration_smoke",
            summary_path=smoke_summary_path,
            report_path=report_path,
            equilibrator=fake_equilibrator,
        )

    assert call_count == 3
    with smoke_summary_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 3
    assert rows[1]["success_status"] == "failure"
    assert rows[1]["failure_reason"] == "RuntimeError: synthetic equilibration failure"
    assert rows[2]["success_status"] == "success"
    assert rows[0]["final_structure_path"].endswith("candidate_a/Dy/npt_final.pdb")
    assert rows[2]["final_structure_path"].endswith("candidate_c/Y/npt_final.pdb")
    assert report_path.exists()
    assert "synthetic equilibration failure" in report_path.read_text(encoding="utf-8")
    assert (
        tmp_path / "results" / "openmm_equilibration_smoke" / "candidate_b" / "Nd" / "equilibration_log.csv"
    ).exists()
