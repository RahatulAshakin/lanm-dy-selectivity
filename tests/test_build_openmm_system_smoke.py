from __future__ import annotations

import csv
from pathlib import Path

import pytest

from lanm.analysis.openmm_system_smoke import (
    OpenMMSystemSmokeSummaryRow,
    parse_openmm_system_build_request,
    prepare_structure_input,
    render_openmm_system_smoke_report,
    run_openmm_system_smoke,
    load_md_validation_panel_rows,
    load_metal_build_readiness_rows,
    select_representative_baseline_systems,
)
from lanm.md.openmm_system_smoke import OpenMMSerializedSystem, build_openmm_serialized_system
from lanm.paths import (
    MD_PROTOCOL_CONFIG_PATH,
    MD_SYSTEM_BUILD_MANIFEST_PATH,
    MD_VALIDATION_PANEL_PATH,
    METAL_BUILD_READINESS_PATH,
    METAL_PARAMETER_VALUES_PATH,
)


def test_select_representative_baseline_systems_matches_phase_6a3_targets() -> None:
    panel_rows = load_md_validation_panel_rows(MD_VALIDATION_PANEL_PATH)
    readiness_rows = load_metal_build_readiness_rows(METAL_BUILD_READINESS_PATH)

    representative_rows = select_representative_baseline_systems(
        panel_rows=panel_rows,
        build_manifest_path=MD_SYSTEM_BUILD_MANIFEST_PATH,
        readiness_rows=readiness_rows,
        md_protocol_path=MD_PROTOCOL_CONFIG_PATH,
        metal_parameter_values_path=METAL_PARAMETER_VALUES_PATH,
    )

    assert [(row.panel_member_id, row.target_metal) for row in representative_rows] == [
        ("am1_mex_ss_only_u02", "Dy"),
        ("hans_pocket_ss_only_u02", "Dy"),
        ("hans_interface_ss_plus_if_u04", "Dy"),
    ]
    assert [row.topology_class for row in representative_rows] == [
        "am1_monomer",
        "hans_monomer",
        "hans_interface_multichain",
    ]


def test_parse_openmm_system_build_request_reads_template_water_fields() -> None:
    request = parse_openmm_system_build_request(
        Path("results/md_system_build/am1_mex_ss_only_u02/Dy/build_request.yaml")
    )

    assert request.panel_member_id == "am1_mex_ss_only_u02"
    assert request.target_metal == "Dy"
    assert request.topology_class == "am1_monomer"
    assert request.chosen_protein_forcefield == "amber19-all.xml"
    assert request.chosen_water_model == "amber19/opc3.xml"
    assert request.restore_template_metal_proximal_waters_before_solvation is True
    assert request.template_metal_proximal_waters_available is True
    assert request.comparison_panel_member_id == "am1_mex_wt_reference"
    assert request.comparison_preserved_solvent_residue_count == 25


def test_build_openmm_serialized_system_writes_valid_xml_for_restored_am1_system() -> None:
    request = parse_openmm_system_build_request(
        Path("results/md_system_build/am1_mex_ss_only_u02/Dy/build_request.yaml")
    )
    prepared = prepare_structure_input(request)

    assert prepared.restored_template_waters == 25

    serialized = build_openmm_serialized_system(
        prepared_structure_pdb_text=prepared.pdb_text,
        forcefield_files=("amber19-all.xml", "amber19/opc3.xml"),
        target_temperature_K=298,
    )

    assert serialized.atom_count > 0
    assert serialized.residue_count > 0
    assert "<System" in serialized.system_xml
    assert "LangevinMiddleIntegrator" in serialized.integrator_xml
    assert " DY " in serialized.prepared_structure_pdb_text
    assert "HOH" in serialized.prepared_structure_pdb_text


def test_run_openmm_system_smoke_writes_mixed_summary_when_one_system_fails(tmp_path: Path) -> None:
    call_count = 0

    def fake_builder(
        *,
        prepared_structure_pdb_text: str,
        forcefield_files: tuple[str, ...],
        target_temperature_K: int,
    ) -> OpenMMSerializedSystem:
        nonlocal call_count
        call_count += 1
        if call_count == 2:
            raise ValueError("synthetic serialization failure")
        return OpenMMSerializedSystem(
            prepared_structure_pdb_text=prepared_structure_pdb_text,
            system_xml="<System/>",
            integrator_xml="<LangevinMiddleIntegrator/>",
            atom_count=111 + call_count,
            residue_count=22 + call_count,
        )

    summary_path = tmp_path / "results" / "tables" / "openmm_system_smoke_summary.csv"
    report_path = tmp_path / "results" / "reports" / "openmm_system_smoke.md"
    with pytest.raises(RuntimeError, match="failed for 1/3 representative systems"):
        run_openmm_system_smoke(
            md_validation_panel_path=MD_VALIDATION_PANEL_PATH,
            md_system_build_manifest_path=MD_SYSTEM_BUILD_MANIFEST_PATH,
            metal_build_readiness_path=METAL_BUILD_READINESS_PATH,
            md_protocol_path=MD_PROTOCOL_CONFIG_PATH,
            metal_parameter_values_path=METAL_PARAMETER_VALUES_PATH,
            output_root=tmp_path / "results" / "openmm_system_smoke",
            summary_path=summary_path,
            report_path=report_path,
            system_builder=fake_builder,
        )

    assert summary_path.exists()
    assert report_path.exists()
    with summary_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 3
    assert [row["success_status"] for row in rows] == ["success", "failure", "success"]
    assert rows[1]["failure_reason"] == "ValueError: synthetic serialization failure"
    report_text = report_path.read_text(encoding="utf-8")
    assert "synthetic serialization failure" in report_text
    assert "hans_pocket_ss_only_u02" in report_text


def test_render_openmm_system_smoke_report_handles_success_and_failure_rows() -> None:
    rows = (
        OpenMMSystemSmokeSummaryRow(
            panel_member_id="am1_mex_ss_only_u02",
            target_metal="Dy",
            topology_class="am1_monomer",
            success_status="success",
            failure_reason="",
            forcefield_files="amber19-all.xml;amber19/opc3.xml",
            custom_metal_family="generic_12_6_4_highly_charged",
            restored_template_waters=25,
            atom_count=1608,
            residue_count=134,
            system_xml_path="results/openmm_system_smoke/am1_mex_ss_only_u02/Dy/system.xml",
        ),
        OpenMMSystemSmokeSummaryRow(
            panel_member_id="hans_pocket_ss_only_u02",
            target_metal="Dy",
            topology_class="hans_monomer",
            success_status="failure",
            failure_reason="ValueError: synthetic serialization failure",
            forcefield_files="amber19-all.xml;amber19/opc3.xml",
            custom_metal_family="generic_12_6_4_highly_charged",
            restored_template_waters=12,
            atom_count=0,
            residue_count=0,
            system_xml_path="",
        ),
    )

    report = render_openmm_system_smoke_report(rows)

    assert "successful OpenMM `System` serializations: `1`" in report
    assert "failures recorded: `1`" in report
    assert "ValueError: synthetic serialization failure" in report
