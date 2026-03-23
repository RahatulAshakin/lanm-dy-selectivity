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
    build_metal_build_readiness_rows,
    build_metal_parameter_records,
    build_parameter_source_family_registry,
    build_per_system_parameter_mapping,
    classify_direct_support_status,
    load_md_system_build_manifest,
    load_metal_model_registry,
    prepare_metal_parameter_mapping,
    validate_metal_model_registry,
)
from lanm.analysis.metal_parameter_values import (
    DEFAULT_OPC3_WATER_MODEL,
    NumericParameterValueRecord,
    build_default_numeric_parameter_registry,
    validate_numeric_parameter_registry,
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


def test_build_default_numeric_parameter_registry_matches_published_opc3_values() -> None:
    registry = build_default_numeric_parameter_registry()

    assert [record.metal_identity for record in registry] == list(TARGET_METALS)
    assert all(record.source_family_label == GENERIC_12_6_4_HIGHLY_CHARGED for record in registry)
    assert all(record.water_model == DEFAULT_OPC3_WATER_MODEL for record in registry)

    mapping = {record.metal_identity: record for record in registry}
    assert mapping["Al"].rmin_half_A == 1.361
    assert mapping["Al"].epsilon_kcal_per_mol == 0.01031847
    assert mapping["Al"].c4_kcal_per_mol_A4 == 363
    assert mapping["Fe"].rmin_half_A == 1.455
    assert mapping["Fe"].epsilon_kcal_per_mol == 0.02662782
    assert mapping["Fe"].c4_kcal_per_mol_A4 == 429
    assert mapping["Y"].rmin_half_A == 1.626
    assert mapping["Y"].epsilon_kcal_per_mol == 0.09289608
    assert mapping["Y"].c4_kcal_per_mol_A4 == 192
    assert mapping["Nd"].rmin_half_A == 1.712
    assert mapping["Nd"].epsilon_kcal_per_mol == 0.14640930
    assert mapping["Nd"].c4_kcal_per_mol_A4 == 184
    assert mapping["Dy"].rmin_half_A == 1.632
    assert mapping["Dy"].epsilon_kcal_per_mol == 0.09620220
    assert mapping["Dy"].c4_kcal_per_mol_A4 == 183


def test_validate_numeric_parameter_registry_rejects_non_opc3_baseline_entry() -> None:
    registry = build_default_numeric_parameter_registry()
    invalid_registry = (
        NumericParameterValueRecord(
            metal_identity="Dy",
            source_family_label=GENERIC_12_6_4_HIGHLY_CHARGED,
            water_model="amber19/opc.xml",
            formal_charge=3,
            rmin_half_A=1.632,
            epsilon_kcal_per_mol=0.09620220,
            c4_kcal_per_mol_A4=183,
            parameter_provenance=registry[0].parameter_provenance,
            notes="invalid water model",
        ),
        *registry[1:],
    )

    with pytest.raises(ValueError, match="must use 'amber19/opc3.xml'"):
        validate_numeric_parameter_registry(invalid_registry)


def test_classify_direct_support_status_handles_direct_proxy_and_needs_derivation() -> None:
    manifest_rows = load_md_system_build_manifest(MD_SYSTEM_BUILD_MANIFEST_PATH)
    family_registry = build_parameter_source_family_registry(
        manifest_rows,
        build_default_numeric_parameter_registry(),
    )

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


def test_build_metal_parameter_records_distinguishes_baseline_and_tuned_families() -> None:
    registry = load_metal_model_registry(METAL_MODEL_REGISTRY_PATH)
    manifest_rows = load_md_system_build_manifest(MD_SYSTEM_BUILD_MANIFEST_PATH)
    family_registry, metal_records = build_metal_parameter_records(
        registry=registry,
        manifest_rows=manifest_rows,
        numeric_parameter_registry=build_default_numeric_parameter_registry(),
    )

    family_lookup = {record.source_family_label: record for record in family_registry}
    assert family_lookup[GENERIC_12_6_4_HIGHLY_CHARGED].repo_numeric_parameters_present is True
    assert family_lookup[CHELATOR_TUNED_12_6_4_LANMODULIN].repo_numeric_parameters_present is False

    mapping = {record.metal_identity: record for record in metal_records}
    assert mapping["Dy"].source_family_label == CHELATOR_TUNED_12_6_4_LANMODULIN
    assert mapping["Dy"].direct_support_status == DIRECT_SUPPORT_STATUS
    assert mapping["Dy"].baseline_parameter_family == GENERIC_12_6_4_HIGHLY_CHARGED
    assert mapping["Dy"].baseline_numeric_values_present is True
    assert mapping["Dy"].tuned_parameter_family == CHELATOR_TUNED_12_6_4_LANMODULIN
    assert mapping["Dy"].tuned_numeric_values_present is False
    assert mapping["Fe"].source_family_label == GENERIC_12_6_4_HIGHLY_CHARGED
    assert mapping["Fe"].direct_support_status == PROXY_SUPPORT_STATUS
    assert mapping["Fe"].baseline_numeric_values_present is True
    assert mapping["Fe"].tuned_numeric_values_present is False


def test_build_per_system_parameter_mapping_marks_baseline_ready_and_tuned_pending() -> None:
    registry = load_metal_model_registry(METAL_MODEL_REGISTRY_PATH)
    manifest_rows = load_md_system_build_manifest(MD_SYSTEM_BUILD_MANIFEST_PATH)
    numeric_registry = build_default_numeric_parameter_registry()
    _, metal_records = build_metal_parameter_records(
        registry=registry,
        manifest_rows=manifest_rows,
        numeric_parameter_registry=numeric_registry,
    )

    system_rows = build_per_system_parameter_mapping(
        manifest_rows=manifest_rows,
        metal_parameter_records=metal_records,
        numeric_parameter_registry=numeric_registry,
    )

    assert len(system_rows) == 30
    assert all(row.baseline_numeric_values_present for row in system_rows)
    assert all(row.ready_for_openmm_system_build_baseline for row in system_rows)
    assert not any(row.tuned_numeric_values_present for row in system_rows)
    assert not any(row.ready_for_openmm_system_build_tuned for row in system_rows)

    dy_row = next(
        row
        for row in system_rows
        if row.panel_member_id == "am1_mex_ss_only_u02" and row.target_metal == "Dy"
    )
    fe_row = next(
        row
        for row in system_rows
        if row.panel_member_id == "am1_mex_ss_only_u02" and row.target_metal == "Fe"
    )
    assert dy_row.chosen_parameter_family == CHELATOR_TUNED_12_6_4_LANMODULIN
    assert dy_row.direct_support_status == DIRECT_SUPPORT_STATUS
    assert dy_row.baseline_parameter_family == GENERIC_12_6_4_HIGHLY_CHARGED
    assert dy_row.ready_for_openmm_system_build_baseline is True
    assert dy_row.ready_for_openmm_system_build_tuned is False
    assert "published generic OPC3 12-6-4 coefficients" in dy_row.rationale
    assert fe_row.chosen_parameter_family == GENERIC_12_6_4_HIGHLY_CHARGED
    assert fe_row.direct_support_status == PROXY_SUPPORT_STATUS
    assert fe_row.ready_for_openmm_system_build_baseline is True
    assert fe_row.ready_for_openmm_system_build_tuned is False


def test_build_metal_build_readiness_rows_keeps_phase_6a2_flag_but_marks_baseline_ready() -> None:
    registry = load_metal_model_registry(METAL_MODEL_REGISTRY_PATH)
    manifest_rows = load_md_system_build_manifest(MD_SYSTEM_BUILD_MANIFEST_PATH)
    numeric_registry = build_default_numeric_parameter_registry()
    _, metal_records = build_metal_parameter_records(
        registry=registry,
        manifest_rows=manifest_rows,
        numeric_parameter_registry=numeric_registry,
    )
    system_rows = build_per_system_parameter_mapping(
        manifest_rows=manifest_rows,
        metal_parameter_records=metal_records,
        numeric_parameter_registry=numeric_registry,
    )

    readiness_rows = build_metal_build_readiness_rows(
        manifest_rows=manifest_rows,
        system_mapping_rows=system_rows,
    )

    assert len(readiness_rows) == 30
    dy_row = next(
        row
        for row in readiness_rows
        if row.panel_member_id == "am1_mex_ss_only_u02" and row.target_metal == "Dy"
    )
    assert dy_row.phase_6a2_custom_metal_parameters_still_required is True
    assert dy_row.baseline_numeric_values_present is True
    assert dy_row.ready_for_openmm_system_build_baseline is True
    assert dy_row.tuned_numeric_values_present is False
    assert dy_row.ready_for_openmm_system_build_tuned is False
    assert dy_row.readiness_summary == "baseline_ready_tuned_pending"


def test_prepare_metal_parameter_mapping_writes_outputs(tmp_path: Path) -> None:
    output_rows = prepare_metal_parameter_mapping(
        metal_model_registry_path=METAL_MODEL_REGISTRY_PATH,
        md_system_build_manifest_path=MD_SYSTEM_BUILD_MANIFEST_PATH,
        metal_parameter_config_path=tmp_path / "config" / "metal_parameters.yaml",
        metal_parameter_values_path=tmp_path / "config" / "metal_parameter_values.yaml",
        metal_parameter_mapping_path=tmp_path / "results" / "tables" / "metal_parameter_mapping.csv",
        metal_build_readiness_path=tmp_path / "results" / "tables" / "metal_build_readiness.csv",
        metal_parameter_strategy_report_path=tmp_path / "results" / "reports" / "metal_parameter_strategy.md",
    )

    assert len(output_rows) == 30

    config_payload = yaml.safe_load((tmp_path / "config" / "metal_parameters.yaml").read_text(encoding="utf-8"))
    assert [entry["metal_identity"] for entry in config_payload["metals"]] == list(TARGET_METALS)
    assert config_payload["parameter_source_families"][0]["repo_numeric_parameters_present"] is True
    assert config_payload["parameter_source_families"][1]["repo_numeric_parameters_present"] is False
    assert config_payload["metals"][0]["source_family_label"] == CHELATOR_TUNED_12_6_4_LANMODULIN
    assert config_payload["metals"][0]["baseline_numeric_values_present"] is True
    assert config_payload["metals"][0]["tuned_numeric_values_present"] is False

    values_payload = yaml.safe_load(
        (tmp_path / "config" / "metal_parameter_values.yaml").read_text(encoding="utf-8")
    )
    assert [entry["metal_identity"] for entry in values_payload["numeric_parameter_values"]] == list(
        TARGET_METALS
    )
    assert values_payload["numeric_parameter_values"][0]["water_model"] == DEFAULT_OPC3_WATER_MODEL

    assert (tmp_path / "results" / "tables" / "metal_parameter_mapping.csv").exists()
    readiness_text = (tmp_path / "results" / "tables" / "metal_build_readiness.csv").read_text(
        encoding="utf-8"
    )
    assert "baseline_ready_tuned_pending" in readiness_text

    report_text = (tmp_path / "results" / "reports" / "metal_parameter_strategy.md").read_text(
        encoding="utf-8"
    )
    assert "Why Phase 6A2b Was Conservative" in report_text
    assert "Published generic OPC3 12-6-4 baseline coefficients are now present" in report_text
    assert "future refinement path" in report_text
