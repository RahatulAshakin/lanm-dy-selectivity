from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from lanm.analysis.metal_parameter_mapping import (
    CHELATOR_TUNED_12_6_4_LANMODULIN,
    DIRECT_SUPPORT_STATUS,
    GENERIC_12_6_4_HIGHLY_CHARGED,
    NEEDS_DERIVATION_STATUS,
    PHASE_6A2_REGISTRY_FAMILY,
    PROXY_SUPPORT_STATUS,
    TARGET_METALS,
    build_metal_parameter_records,
    build_parameter_source_family_registry,
    build_per_system_parameter_mapping,
    classify_direct_support_status,
    load_md_system_build_manifest,
    load_metal_model_registry,
    prepare_metal_parameter_mapping,
    validate_metal_model_registry,
)
from lanm.paths import MD_SYSTEM_BUILD_MANIFEST_PATH, METAL_MODEL_REGISTRY_PATH


def test_load_metal_model_registry_matches_phase_6a2_outputs() -> None:
    registry = load_metal_model_registry(METAL_MODEL_REGISTRY_PATH)

    assert len(registry) == 5
    assert [record.metal_identity for record in registry] == list(TARGET_METALS)
    assert all(record.custom_required for record in registry)
    assert all(record.intended_model_family == PHASE_6A2_REGISTRY_FAMILY for record in registry)


def test_validate_metal_model_registry_rejects_noncustom_entry() -> None:
    registry = load_metal_model_registry(METAL_MODEL_REGISTRY_PATH)
    invalid_registry = (
        registry[0],
        registry[1],
        registry[2],
        registry[3],
        registry[4].__class__(
            metal_identity="Fe",
            formal_charge=3,
            intended_model_family=PHASE_6A2_REGISTRY_FAMILY,
            standard_forcefield_supported=False,
            custom_required=False,
            parameter_source_status="custom_parameters_not_present_in_repo",
            notes="invalid",
        ),
    )

    with pytest.raises(ValueError, match="must remain custom_required"):
        validate_metal_model_registry(invalid_registry)


def test_classify_direct_support_status_handles_direct_proxy_and_needs_derivation() -> None:
    manifest_rows = load_md_system_build_manifest(MD_SYSTEM_BUILD_MANIFEST_PATH)
    family_registry = build_parameter_source_family_registry(manifest_rows)

    assert (
        classify_direct_support_status(
            metal_identity="Dy",
            intended_parameter_family=CHELATOR_TUNED_12_6_4_LANMODULIN,
            source_family_label=CHELATOR_TUNED_12_6_4_LANMODULIN,
            source_family_registry=family_registry,
        )
        == DIRECT_SUPPORT_STATUS
    )
    assert (
        classify_direct_support_status(
            metal_identity="Fe",
            intended_parameter_family=CHELATOR_TUNED_12_6_4_LANMODULIN,
            source_family_label=GENERIC_12_6_4_HIGHLY_CHARGED,
            source_family_registry=family_registry,
        )
        == PROXY_SUPPORT_STATUS
    )
    assert (
        classify_direct_support_status(
            metal_identity="Dy",
            intended_parameter_family=CHELATOR_TUNED_12_6_4_LANMODULIN,
            source_family_label=None,
            source_family_registry=family_registry,
        )
        == NEEDS_DERIVATION_STATUS
    )


def test_build_metal_parameter_records_assigns_expected_family_by_metal() -> None:
    registry = load_metal_model_registry(METAL_MODEL_REGISTRY_PATH)
    manifest_rows = load_md_system_build_manifest(MD_SYSTEM_BUILD_MANIFEST_PATH)

    _, metal_records = build_metal_parameter_records(
        registry=registry,
        manifest_rows=manifest_rows,
    )

    mapping = {record.metal_identity: record for record in metal_records}
    assert mapping["Dy"].source_family_label == CHELATOR_TUNED_12_6_4_LANMODULIN
    assert mapping["Nd"].source_family_label == CHELATOR_TUNED_12_6_4_LANMODULIN
    assert mapping["Y"].source_family_label == CHELATOR_TUNED_12_6_4_LANMODULIN
    assert mapping["Al"].source_family_label == GENERIC_12_6_4_HIGHLY_CHARGED
    assert mapping["Fe"].source_family_label == GENERIC_12_6_4_HIGHLY_CHARGED
    assert mapping["Dy"].direct_support_status == DIRECT_SUPPORT_STATUS
    assert mapping["Fe"].direct_support_status == PROXY_SUPPORT_STATUS
    assert mapping["Dy"].water_model_compatibility == ("amber19/opc3.xml",)


def test_build_per_system_parameter_mapping_marks_all_current_panel_rows_not_ready() -> None:
    registry = load_metal_model_registry(METAL_MODEL_REGISTRY_PATH)
    manifest_rows = load_md_system_build_manifest(MD_SYSTEM_BUILD_MANIFEST_PATH)
    family_registry, metal_records = build_metal_parameter_records(
        registry=registry,
        manifest_rows=manifest_rows,
    )

    system_rows = build_per_system_parameter_mapping(
        manifest_rows=manifest_rows,
        metal_parameter_records=metal_records,
        source_family_registry=family_registry,
    )

    assert len(system_rows) == 30
    assert not any(row.ready_for_openmm_system_build for row in system_rows)
    dy_row = next(row for row in system_rows if row.panel_member_id == "am1_mex_ss_only_u02" and row.target_metal == "Dy")
    fe_row = next(row for row in system_rows if row.panel_member_id == "am1_mex_ss_only_u02" and row.target_metal == "Fe")
    assert dy_row.chosen_parameter_family == CHELATOR_TUNED_12_6_4_LANMODULIN
    assert dy_row.direct_support_status == DIRECT_SUPPORT_STATUS
    assert "numeric parameters still need to be derived" in dy_row.rationale
    assert fe_row.chosen_parameter_family == GENERIC_12_6_4_HIGHLY_CHARGED
    assert fe_row.direct_support_status == PROXY_SUPPORT_STATUS
    assert "generic 12-6-4 proxy" in fe_row.rationale


def test_prepare_metal_parameter_mapping_writes_outputs(tmp_path: Path) -> None:
    output_rows = prepare_metal_parameter_mapping(
        metal_model_registry_path=METAL_MODEL_REGISTRY_PATH,
        md_system_build_manifest_path=MD_SYSTEM_BUILD_MANIFEST_PATH,
        metal_parameter_config_path=tmp_path / "config" / "metal_parameters.yaml",
        metal_parameter_mapping_path=tmp_path / "results" / "tables" / "metal_parameter_mapping.csv",
        metal_parameter_strategy_report_path=tmp_path / "results" / "reports" / "metal_parameter_strategy.md",
    )

    assert len(output_rows) == 30
    config_payload = yaml.safe_load((tmp_path / "config" / "metal_parameters.yaml").read_text(encoding="utf-8"))
    assert [entry["metal_identity"] for entry in config_payload["metals"]] == list(TARGET_METALS)
    assert config_payload["metals"][0]["source_family_label"] == CHELATOR_TUNED_12_6_4_LANMODULIN
    assert config_payload["metals"][-1]["source_family_label"] == GENERIC_12_6_4_HIGHLY_CHARGED
    assert (tmp_path / "results" / "tables" / "metal_parameter_mapping.csv").exists()
    report_text = (tmp_path / "results" / "reports" / "metal_parameter_strategy.md").read_text(encoding="utf-8")
    assert "Metals with only generic baseline support" in report_text
    assert "Metals with LanM-adjacent chelator-tuned support" in report_text
