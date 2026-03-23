from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from lanm.analysis.md_system_build import (
    DEFAULT_WATER_MODEL,
    MetalModelRecord,
    TARGET_METALS,
    build_default_metal_model_registry,
    discover_md_input_manifest,
    load_md_protocol_build_config,
    plan_template_water_restoration,
    prepare_md_system_build,
    validate_md_input_manifest_against_protocol,
    validate_metal_model_registry,
)
from lanm.paths import MD_INPUT_MANIFEST_PATH, MD_PROTOCOL_CONFIG_PATH


def test_discover_md_input_manifest_matches_phase_6a1_outputs() -> None:
    manifest_rows = discover_md_input_manifest(MD_INPUT_MANIFEST_PATH)

    assert len(manifest_rows) == 30
    assert tuple(sorted({row.target_metal for row in manifest_rows}, key=TARGET_METALS.index)) == TARGET_METALS
    assert manifest_rows[0].panel_member_id == "am1_mex_ss_only_u02"
    assert [row.panel_rank for row in manifest_rows[:5]] == [1, 1, 1, 1, 1]
    assert manifest_rows[-1].panel_member_id == "hans_interface_ss_plus_if_u02"
    validate_md_input_manifest_against_protocol(
        manifest_rows=manifest_rows,
        protocol_path=MD_PROTOCOL_CONFIG_PATH,
    )


def test_prepare_md_system_build_writes_requests_with_protocol_water_override(tmp_path: Path) -> None:
    protocol_payload = yaml.safe_load(MD_PROTOCOL_CONFIG_PATH.read_text(encoding="utf-8"))
    protocol_payload["system_build_defaults"] = {
        "protein_forcefield": "amber19-all.xml",
        "water_model": "amber19/opc.xml",
        "box_padding_nm": 1.0,
        "ionic_strength_M": 0.20,
        "neutralization_ion_policy": "monovalent_background_ions_only_excluding_panel_metals",
    }
    protocol_path = tmp_path / "config" / "md_protocol.yaml"
    protocol_path.parent.mkdir(parents=True, exist_ok=True)
    protocol_path.write_text(yaml.safe_dump(protocol_payload, sort_keys=False), encoding="utf-8")

    manifest_rows = prepare_md_system_build(
        md_input_manifest_path=MD_INPUT_MANIFEST_PATH,
        md_protocol_path=protocol_path,
        metal_model_registry_path=tmp_path / "config" / "metal_models.yaml",
        md_system_build_root=tmp_path / "results" / "md_system_build",
        md_system_build_manifest_path=tmp_path / "results" / "tables" / "md_system_build_manifest.csv",
        md_system_build_report_path=tmp_path / "results" / "reports" / "md_system_building_plan.md",
    )

    assert len(manifest_rows) == 30

    build_request_path = tmp_path / "results" / "md_system_build" / "am1_mex_ss_only_u02" / "Dy" / "build_request.yaml"
    build_request = yaml.safe_load(build_request_path.read_text(encoding="utf-8"))
    assert build_request["source_starting_structure_pdb"] == "results/md_inputs/am1_mex_ss_only_u02/Dy/starting_structure.pdb"
    assert build_request["topology_class"] == "am1_monomer"
    assert build_request["target_metal"] == "Dy"
    assert build_request["chosen_protein_forcefield"] == "amber19-all.xml"
    assert build_request["chosen_water_model"] == "amber19/opc.xml"
    assert build_request["box_padding_target_nm"] == 1.0
    assert build_request["ionic_strength_target_M"] == 0.20
    assert build_request["neutralization_ion_policy"] == "monovalent_background_ions_only_excluding_panel_metals"
    assert build_request["custom_metal_parameters_still_required"] is True
    assert build_request["restore_template_metal_proximal_waters_before_solvation"] is True
    assert build_request["template_metal_proximal_waters_available"] is True
    assert build_request["template_water_restoration_reference"]["comparison_panel_member_id"] == "am1_mex_wt_reference"

    wild_type_request_path = (
        tmp_path / "results" / "md_system_build" / "hans_interface_wt_reference" / "Dy" / "build_request.yaml"
    )
    wild_type_request = yaml.safe_load(wild_type_request_path.read_text(encoding="utf-8"))
    assert wild_type_request["restore_template_metal_proximal_waters_before_solvation"] is False
    assert wild_type_request["template_metal_proximal_waters_available"] is True

    registry_payload = yaml.safe_load((tmp_path / "config" / "metal_models.yaml").read_text(encoding="utf-8"))
    assert [entry["metal_identity"] for entry in registry_payload["metals"]] == list(TARGET_METALS)
    assert all(entry["custom_required"] for entry in registry_payload["metals"])
    assert "custom_parameters_not_present_in_repo" in (tmp_path / "results" / "reports" / "md_system_building_plan.md").read_text(
        encoding="utf-8"
    )


def test_validate_metal_model_registry_rejects_ambiguous_support_mode() -> None:
    valid_registry = validate_metal_model_registry(build_default_metal_model_registry())

    assert len(valid_registry) == 5
    assert [record.formal_charge for record in valid_registry] == [3, 3, 3, 3, 3]

    invalid_registry = (
        MetalModelRecord(
            metal_identity="Dy",
            formal_charge=3,
            intended_model_family="custom_bound_site",
            standard_forcefield_supported=False,
            custom_required=False,
            parameter_source_status="missing",
            notes="invalid support mode",
        ),
        *valid_registry[1:],
    )
    with pytest.raises(ValueError, match="exactly one support mode"):
        validate_metal_model_registry(invalid_registry)


def test_plan_template_water_restoration_uses_matching_template_family_for_zero_solvent_designs() -> None:
    manifest_rows = discover_md_input_manifest(MD_INPUT_MANIFEST_PATH)
    hans_pocket_row = next(
        row
        for row in manifest_rows
        if row.panel_member_id == "hans_pocket_ss_only_u02" and row.target_metal == "Dy"
    )

    plan = plan_template_water_restoration(
        system_row=hans_pocket_row,
        manifest_rows=manifest_rows,
    )

    assert plan.template_metal_proximal_waters_available is True
    assert plan.restore_template_metal_proximal_waters_before_solvation is True
    assert plan.comparison_panel_member_id == "hans_interface_wt_reference"
    assert plan.comparison_source == "data/raw/public/structures/8FNR.cif"
    assert plan.comparison_preserved_solvent_residue_count == 40


def test_default_protocol_water_model_falls_back_to_opc3() -> None:
    build_config = load_md_protocol_build_config(MD_PROTOCOL_CONFIG_PATH)

    assert build_config.protein_forcefield == "amber19-all.xml"
    assert build_config.water_model == DEFAULT_WATER_MODEL
    assert build_config.box_padding_nm > 0.0
    assert build_config.ionic_strength_M >= 0.0
    assert DEFAULT_WATER_MODEL == "amber19/opc3.xml"
