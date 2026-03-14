from __future__ import annotations

import csv
from pathlib import Path

import pytest
import yaml

from lanm.analysis.openmm_system_build import (
    build_simulation_config,
    discover_baseline_ready_openmm_systems,
    load_openmm_system_build_runtime_config,
    run_openmm_system_build,
)
from lanm.analysis.openmm_system_smoke import parse_openmm_system_build_request, prepare_structure_input
from lanm.md.openmm_system_build import OpenMMSerializedSystem, build_openmm_serialized_system
from lanm.paths import (
    MD_INPUT_MANIFEST_PATH,
    MD_PROTOCOL_CONFIG_PATH,
    MD_SYSTEM_BUILD_MANIFEST_PATH,
    METAL_BUILD_READINESS_PATH,
    METAL_PARAMETER_VALUES_PATH,
)


def test_discover_baseline_ready_openmm_systems_matches_full_panel() -> None:
    rows = discover_baseline_ready_openmm_systems(
        md_input_manifest_path=MD_INPUT_MANIFEST_PATH,
        md_system_build_manifest_path=MD_SYSTEM_BUILD_MANIFEST_PATH,
        metal_build_readiness_path=METAL_BUILD_READINESS_PATH,
        md_protocol_path=MD_PROTOCOL_CONFIG_PATH,
        metal_parameter_values_path=METAL_PARAMETER_VALUES_PATH,
    )

    assert len(rows) == 30
    assert rows[0].panel_member_id == "am1_mex_ss_only_u02"
    assert rows[0].target_metal == "Al"
    assert rows[-1].panel_member_id == "hans_interface_ss_plus_if_u02"
    assert rows[-1].target_metal == "Y"
    assert all(row.baseline_parameter_family == "generic_12_6_4_highly_charged" for row in rows)
    assert all(row.chosen_protein_forcefield == "amber19-all.xml" for row in rows)
    assert all(row.chosen_water_model == "amber19/opc3.xml" for row in rows)
    assert {row.topology_class for row in rows} == {
        "am1_monomer",
        "hans_monomer",
        "hans_interface_multichain",
    }


def test_build_simulation_config_uses_phase_6a4_defaults() -> None:
    system_row = discover_baseline_ready_openmm_systems(
        md_input_manifest_path=MD_INPUT_MANIFEST_PATH,
        md_system_build_manifest_path=MD_SYSTEM_BUILD_MANIFEST_PATH,
        metal_build_readiness_path=METAL_BUILD_READINESS_PATH,
        md_protocol_path=MD_PROTOCOL_CONFIG_PATH,
        metal_parameter_values_path=METAL_PARAMETER_VALUES_PATH,
    )[0]
    runtime_config = load_openmm_system_build_runtime_config(MD_PROTOCOL_CONFIG_PATH)

    simulation_config = build_simulation_config(
        system_row=system_row,
        runtime_config=runtime_config,
    )

    assert simulation_config.panel_member_id == system_row.panel_member_id
    assert simulation_config.target_metal == system_row.target_metal
    assert simulation_config.topology_class == system_row.topology_class
    assert simulation_config.replicate_count == 3
    assert simulation_config.target_temperature_K == 298
    assert simulation_config.friction_coeff_ps == 1.0
    assert simulation_config.timestep_fs == 2.0
    assert simulation_config.nonbonded_method == "NoCutoff"
    assert simulation_config.cutoff_nm == 1.0
    assert simulation_config.hydrogen_mass_repartitioning is False


def test_run_openmm_system_build_iterates_full_panel_and_writes_configs(tmp_path: Path) -> None:
    call_count = 0

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
        nonlocal call_count
        call_count += 1
        assert forcefield_files == ("amber19-all.xml", "amber19/opc3.xml")
        assert target_temperature_K == 298
        assert friction_coeff_ps == 1.0
        assert timestep_fs == 2.0
        assert nonbonded_method == "NoCutoff"
        assert cutoff_nm == 1.0
        assert hydrogen_mass_repartitioning is False
        return OpenMMSerializedSystem(
            prepared_structure_pdb_text=prepared_structure_pdb_text,
            system_xml=f"<System id='{call_count}'/>",
            integrator_xml="<Integrator type='LangevinMiddleIntegrator'/>",
            atom_count=1000 + call_count,
            residue_count=100 + call_count,
        )

    summary_path = tmp_path / "results" / "tables" / "openmm_system_build_summary.csv"
    report_path = tmp_path / "results" / "reports" / "openmm_system_building.md"
    rows = run_openmm_system_build(
        md_input_manifest_path=MD_INPUT_MANIFEST_PATH,
        md_system_build_manifest_path=MD_SYSTEM_BUILD_MANIFEST_PATH,
        metal_build_readiness_path=METAL_BUILD_READINESS_PATH,
        md_protocol_path=MD_PROTOCOL_CONFIG_PATH,
        metal_parameter_values_path=METAL_PARAMETER_VALUES_PATH,
        output_root=tmp_path / "results" / "openmm_system_build",
        summary_path=summary_path,
        report_path=report_path,
        system_builder=fake_builder,
    )

    assert call_count == 30
    assert len(rows) == 30
    assert summary_path.exists()
    assert report_path.exists()
    first_config_path = tmp_path / "results" / "openmm_system_build" / "am1_mex_ss_only_u02" / "Al" / "simulation_config.yaml"
    assert first_config_path.exists()
    config_payload = yaml.safe_load(first_config_path.read_text(encoding="utf-8"))
    assert config_payload["replicate_count"] == 3
    assert config_payload["target_temperature_K"] == 298
    assert config_payload["friction_coeff_ps"] == 1.0
    assert config_payload["timestep_fs"] == 2.0
    assert config_payload["nonbonded_method"] == "NoCutoff"
    assert config_payload["cutoff_nm"] == 1.0
    assert config_payload["hydrogen_mass_repartitioning"] is False


def test_run_openmm_system_build_writes_failure_tolerant_summary(tmp_path: Path) -> None:
    call_count = 0

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
        nonlocal call_count
        call_count += 1
        if call_count == 7:
            raise ValueError("synthetic baseline build failure")
        return OpenMMSerializedSystem(
            prepared_structure_pdb_text=prepared_structure_pdb_text,
            system_xml=f"<System id='{call_count}'/>",
            integrator_xml="<Integrator type='LangevinMiddleIntegrator'/>",
            atom_count=2000 + call_count,
            residue_count=200 + call_count,
        )

    summary_path = tmp_path / "results" / "tables" / "openmm_system_build_summary.csv"
    report_path = tmp_path / "results" / "reports" / "openmm_system_building.md"
    with pytest.raises(RuntimeError, match="failed for 1/30 baseline-ready systems"):
        run_openmm_system_build(
            md_input_manifest_path=MD_INPUT_MANIFEST_PATH,
            md_system_build_manifest_path=MD_SYSTEM_BUILD_MANIFEST_PATH,
            metal_build_readiness_path=METAL_BUILD_READINESS_PATH,
            md_protocol_path=MD_PROTOCOL_CONFIG_PATH,
            metal_parameter_values_path=METAL_PARAMETER_VALUES_PATH,
            output_root=tmp_path / "results" / "openmm_system_build",
            summary_path=summary_path,
            report_path=report_path,
            system_builder=fake_builder,
        )

    with summary_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 30
    assert rows[6]["success_status"] == "failure"
    assert rows[6]["failure_reason"] == "ValueError: synthetic baseline build failure"
    assert report_path.exists()
    assert "synthetic baseline build failure" in report_path.read_text(encoding="utf-8")


def test_build_openmm_serialized_system_resolves_fe_template_ambiguity() -> None:
    request = parse_openmm_system_build_request(
        Path("results/md_system_build/am1_mex_ss_only_u02/Fe/build_request.yaml")
    )
    prepared = prepare_structure_input(request)

    serialized = build_openmm_serialized_system(
        prepared_structure_pdb_text=prepared.pdb_text,
        forcefield_files=("amber19-all.xml", "amber19/opc3.xml"),
        target_temperature_K=298,
        friction_coeff_ps=1.0,
        timestep_fs=2.0,
        nonbonded_method="NoCutoff",
        cutoff_nm=1.0,
        hydrogen_mass_repartitioning=False,
    )

    assert serialized.atom_count > 0
    assert serialized.residue_count > 0
    assert "<System" in serialized.system_xml
    assert "LangevinMiddleIntegrator" in serialized.integrator_xml
    assert " FE " in serialized.prepared_structure_pdb_text
