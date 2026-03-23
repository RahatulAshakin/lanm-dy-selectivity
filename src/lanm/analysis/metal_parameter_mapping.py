"""Deterministic Phase 6A2c metal-parameter mapping and readiness outputs."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from lanm.analysis.metal_parameter_values import (
    CHELATOR_TUNED_12_6_4_LANMODULIN,
    DEFAULT_OPC3_WATER_MODEL,
    EXPLICIT_TUNED_NUMERIC_COEFFICIENTS,
    GENERIC_12_6_4_HIGHLY_CHARGED,
    PUBLISHED_OPC3_12_6_4_BASELINE,
    TARGET_METALS,
    NumericParameterValueRecord,
    build_default_numeric_parameter_registry,
    build_numeric_parameter_family_registry,
    render_metal_parameter_values_payload,
    validate_numeric_parameter_registry,
)
from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_text, write_csv_rows, write_yaml
from lanm.paths import REPO_ROOT

PHASE_6A2_REGISTRY_FAMILY = "custom_bound_site_12_6_4_lj_chelator_tuned"
PHASE_6A2_PARAMETER_SOURCE_STATUS = "custom_parameters_not_present_in_repo"
DIRECT_SUPPORT_STATUS = "direct"
PROXY_SUPPORT_STATUS = "proxy"
NEEDS_DERIVATION_STATUS = "needs_derivation"
VALID_SUPPORT_STATUSES = (
    DIRECT_SUPPORT_STATUS,
    PROXY_SUPPORT_STATUS,
    NEEDS_DERIVATION_STATUS,
)
LANM_ADJACENT_DIRECT_METALS = ("Dy", "Nd", "Y")
GENERIC_BASELINE_ONLY_METALS = ("Al", "Fe")


@dataclass(frozen=True, slots=True)
class MetalModelRegistryRecord:
    metal_identity: str
    formal_charge: int
    intended_model_family: str
    standard_forcefield_supported: bool
    custom_required: bool
    parameter_source_status: str
    notes: str


@dataclass(frozen=True, slots=True)
class MDSystemBuildManifestRow:
    panel_rank: int
    panel_member_id: str
    panel_member_type: str
    panel_role: str
    topology_class: str
    target_metal: str
    source_starting_structure_pdb: str
    build_request_path: str
    chosen_protein_forcefield: str
    chosen_water_model: str
    box_padding_target_nm: float
    ionic_strength_target_M: float
    neutralization_ion_policy: str
    custom_metal_parameters_still_required: bool
    generic_protein_water_preparation_ready: bool
    template_metal_proximal_waters_available: bool
    restore_template_metal_proximal_waters_before_solvation: bool
    template_water_reference_panel_member_id: str
    template_water_reference_source: str


@dataclass(frozen=True, slots=True)
class ParameterSourceFamilyRecord:
    source_family_label: str
    family_scope: str
    covered_metals: tuple[str, ...]
    water_model_compatibility: tuple[str, ...]
    repo_numeric_parameters_present: bool
    parameter_provenance: str
    notes: str


@dataclass(frozen=True, slots=True)
class MetalParameterRecord:
    metal_identity: str
    formal_charge: int
    intended_parameter_family: str
    source_family_label: str
    direct_support_status: str
    baseline_parameter_family: str
    baseline_numeric_values_present: bool
    tuned_parameter_family: str
    tuned_numeric_values_present: bool
    water_model_compatibility: tuple[str, ...]
    notes: str


@dataclass(frozen=True, slots=True)
class MetalParameterMappingRow:
    panel_member_id: str
    target_metal: str
    chosen_parameter_family: str
    direct_support_status: str
    baseline_parameter_family: str
    baseline_numeric_values_present: bool
    tuned_parameter_family: str
    tuned_numeric_values_present: bool
    ready_for_openmm_system_build_baseline: bool
    ready_for_openmm_system_build_tuned: bool
    rationale: str


@dataclass(frozen=True, slots=True)
class MetalBuildReadinessRow:
    panel_rank: int
    panel_member_id: str
    target_metal: str
    topology_class: str
    build_request_path: str
    chosen_water_model: str
    generic_protein_water_preparation_ready: bool
    phase_6a2_custom_metal_parameters_still_required: bool
    baseline_parameter_family: str
    baseline_numeric_values_present: bool
    ready_for_openmm_system_build_baseline: bool
    tuned_parameter_family: str
    tuned_numeric_values_present: bool
    ready_for_openmm_system_build_tuned: bool
    readiness_summary: str


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


def _parse_csv_int(row: dict[str, str], field_name: str) -> int:
    return int(str(row[field_name]).strip())


def _parse_csv_float(row: dict[str, str], field_name: str) -> float:
    return float(str(row[field_name]).strip())


def _parse_csv_bool(row: dict[str, str], field_name: str) -> bool:
    value = str(row[field_name]).strip().lower()
    if value == "true":
        return True
    if value == "false":
        return False
    raise ValueError(f"Expected boolean text for {field_name}, found {row[field_name]!r}")


def load_metal_model_registry(path: Path) -> tuple[MetalModelRegistryRecord, ...]:
    """Load and validate the Phase 6A2 metal-model registry."""
    payload = _parse_required_mapping(path)
    metals = payload.get("metals")
    if not isinstance(metals, list):
        raise ValueError(f"Expected metals list in {path}")
    registry = tuple(
        MetalModelRegistryRecord(
            metal_identity=str(entry["metal_identity"]).strip(),
            formal_charge=int(entry["formal_charge"]),
            intended_model_family=str(entry["intended_model_family"]).strip(),
            standard_forcefield_supported=bool(entry["standard_forcefield_supported"]),
            custom_required=bool(entry["custom_required"]),
            parameter_source_status=str(entry["parameter_source_status"]).strip(),
            notes=str(entry["notes"]).strip(),
        )
        for entry in metals
    )
    return validate_metal_model_registry(registry)


def validate_metal_model_registry(
    registry: tuple[MetalModelRegistryRecord, ...],
) -> tuple[MetalModelRegistryRecord, ...]:
    """Validate the Phase 6A2 registry before translating it into source families."""
    if not registry:
        raise ValueError("Metal-model registry is empty")
    observed_metals = tuple(record.metal_identity for record in registry)
    if observed_metals != TARGET_METALS:
        raise ValueError(f"Expected registry metals {TARGET_METALS}, found {observed_metals}")
    for record in registry:
        if record.formal_charge != 3:
            raise ValueError(f"Registry entry for {record.metal_identity} must keep formal_charge 3")
        if record.standard_forcefield_supported:
            raise ValueError(
                f"Registry entry for {record.metal_identity} unexpectedly claims standard force-field support"
            )
        if not record.custom_required:
            raise ValueError(f"Registry entry for {record.metal_identity} must remain custom_required")
        if record.intended_model_family != PHASE_6A2_REGISTRY_FAMILY:
            raise ValueError(
                f"Registry entry for {record.metal_identity} must use {PHASE_6A2_REGISTRY_FAMILY!r}"
            )
        if record.parameter_source_status != PHASE_6A2_PARAMETER_SOURCE_STATUS:
            raise ValueError(
                f"Registry entry for {record.metal_identity} must use "
                f"{PHASE_6A2_PARAMETER_SOURCE_STATUS!r}"
            )
        if not record.notes:
            raise ValueError(f"Registry entry for {record.metal_identity} is missing notes")
    return registry


def load_md_system_build_manifest(path: Path) -> tuple[MDSystemBuildManifestRow, ...]:
    """Load and validate the Phase 6A2 system-build manifest."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = tuple(
            MDSystemBuildManifestRow(
                panel_rank=_parse_csv_int(row, "panel_rank"),
                panel_member_id=str(row["panel_member_id"]).strip(),
                panel_member_type=str(row["panel_member_type"]).strip(),
                panel_role=str(row["panel_role"]).strip(),
                topology_class=str(row["topology_class"]).strip(),
                target_metal=str(row["target_metal"]).strip(),
                source_starting_structure_pdb=str(row["source_starting_structure_pdb"]).strip(),
                build_request_path=str(row["build_request_path"]).strip(),
                chosen_protein_forcefield=str(row["chosen_protein_forcefield"]).strip(),
                chosen_water_model=str(row["chosen_water_model"]).strip(),
                box_padding_target_nm=_parse_csv_float(row, "box_padding_target_nm"),
                ionic_strength_target_M=_parse_csv_float(row, "ionic_strength_target_M"),
                neutralization_ion_policy=str(row["neutralization_ion_policy"]).strip(),
                custom_metal_parameters_still_required=_parse_csv_bool(
                    row,
                    "custom_metal_parameters_still_required",
                ),
                generic_protein_water_preparation_ready=_parse_csv_bool(
                    row,
                    "generic_protein_water_preparation_ready",
                ),
                template_metal_proximal_waters_available=_parse_csv_bool(
                    row,
                    "template_metal_proximal_waters_available",
                ),
                restore_template_metal_proximal_waters_before_solvation=_parse_csv_bool(
                    row,
                    "restore_template_metal_proximal_waters_before_solvation",
                ),
                template_water_reference_panel_member_id=str(
                    row["template_water_reference_panel_member_id"]
                ).strip(),
                template_water_reference_source=str(row["template_water_reference_source"]).strip(),
            )
            for row in reader
        )
    if not rows:
        raise ValueError(f"No rows found in {path}")
    sorted_rows = tuple(
        sorted(
            rows,
            key=lambda row: (
                row.panel_rank,
                row.panel_member_id,
                TARGET_METALS.index(row.target_metal),
            ),
        )
    )
    observed_metals = tuple(sorted({row.target_metal for row in sorted_rows}, key=TARGET_METALS.index))
    if observed_metals != TARGET_METALS:
        raise ValueError(f"Expected manifest metals {TARGET_METALS}, found {observed_metals}")
    for row in sorted_rows:
        require_path(_repo_path(row.source_starting_structure_pdb))
        require_path(_repo_path(row.build_request_path))
    return sorted_rows


def build_parameter_source_family_registry(
    manifest_rows: tuple[MDSystemBuildManifestRow, ...],
    numeric_parameter_registry: tuple[NumericParameterValueRecord, ...],
) -> tuple[ParameterSourceFamilyRecord, ...]:
    """Build the Phase 6A2c source-family registry with numeric readiness annotations."""
    if not manifest_rows:
        raise ValueError("No manifest rows were provided")
    observed_water_models = {row.chosen_water_model for row in manifest_rows}
    family_registry = build_numeric_parameter_family_registry(
        validate_numeric_parameter_registry(numeric_parameter_registry)
    )
    output_rows: list[ParameterSourceFamilyRecord] = []
    for family in family_registry:
        water_model_compatibility = tuple(
            model for model in family.water_model_compatibility if model in observed_water_models
        )
        if not water_model_compatibility:
            water_model_compatibility = family.water_model_compatibility
        output_rows.append(
            ParameterSourceFamilyRecord(
                source_family_label=family.source_family_label,
                family_scope=family.family_scope,
                covered_metals=family.covered_metals,
                water_model_compatibility=water_model_compatibility,
                repo_numeric_parameters_present=family.numeric_values_present,
                parameter_provenance=family.parameter_provenance,
                notes=family.notes,
            )
        )
    return tuple(output_rows)


def _parameter_source_families_by_label(
    families: tuple[ParameterSourceFamilyRecord, ...],
) -> dict[str, ParameterSourceFamilyRecord]:
    return {family.source_family_label: family for family in families}


def classify_direct_support_status(
    *,
    metal_identity: str,
    intended_parameter_family: str,
    source_family_label: str | None,
    source_family_registry: tuple[ParameterSourceFamilyRecord, ...],
) -> str:
    """Classify a source family as direct, proxy, or derivation-needed."""
    if not source_family_label:
        return NEEDS_DERIVATION_STATUS
    family = _parameter_source_families_by_label(source_family_registry).get(source_family_label)
    if family is None:
        return NEEDS_DERIVATION_STATUS
    if metal_identity not in family.covered_metals:
        return NEEDS_DERIVATION_STATUS
    if source_family_label == intended_parameter_family:
        return DIRECT_SUPPORT_STATUS
    return PROXY_SUPPORT_STATUS


def _translate_intended_parameter_family(record: MetalModelRegistryRecord) -> str:
    if record.intended_model_family != PHASE_6A2_REGISTRY_FAMILY:
        raise ValueError(
            f"Cannot translate unexpected intended_model_family {record.intended_model_family!r} "
            f"for {record.metal_identity}"
        )
    return CHELATOR_TUNED_12_6_4_LANMODULIN


def _choose_source_family_label(metal_identity: str) -> str:
    if metal_identity in LANM_ADJACENT_DIRECT_METALS:
        return CHELATOR_TUNED_12_6_4_LANMODULIN
    if metal_identity in GENERIC_BASELINE_ONLY_METALS:
        return GENERIC_12_6_4_HIGHLY_CHARGED
    raise ValueError(f"Unsupported panel metal for source-family mapping: {metal_identity}")


def _numeric_parameter_records_by_key(
    records: tuple[NumericParameterValueRecord, ...],
) -> dict[tuple[str, str, str], NumericParameterValueRecord]:
    return {
        (record.source_family_label, record.metal_identity, record.water_model): record
        for record in records
    }


def _family_has_numeric_values(
    *,
    metal_identity: str,
    source_family_label: str,
    water_model: str | None,
    numeric_parameter_registry: tuple[NumericParameterValueRecord, ...],
) -> bool:
    numeric_lookup = _numeric_parameter_records_by_key(numeric_parameter_registry)
    if water_model is not None:
        return (source_family_label, metal_identity, water_model) in numeric_lookup
    return any(
        record.source_family_label == source_family_label and record.metal_identity == metal_identity
        for record in numeric_parameter_registry
    )


def _build_metal_parameter_notes(
    *,
    metal_identity: str,
    support_status: str,
    source_family_label: str,
    baseline_numeric_values_present: bool,
    tuned_numeric_values_present: bool,
) -> str:
    baseline_note = (
        "Published generic OPC3 12-6-4 baseline coefficients are present for deterministic baseline builds."
        if baseline_numeric_values_present
        else "Published generic OPC3 12-6-4 baseline coefficients are not yet present in the repo."
    )
    tuned_note = (
        "Explicit numeric LanM-tuned coefficients are present for this metal."
        if tuned_numeric_values_present
        else "Explicit numeric LanM-tuned coefficients are still absent."
    )
    if support_status == DIRECT_SUPPORT_STATUS:
        return (
            f"{metal_identity} retains the conceptual LanM-adjacent family {source_family_label!r} as the "
            f"direct Phase 6A2b mapping. {baseline_note} {tuned_note}"
        )
    if support_status == PROXY_SUPPORT_STATUS:
        return (
            f"{metal_identity} continues to map conceptually through the generic highly charged family "
            f"{source_family_label!r}. {baseline_note} {tuned_note}"
        )
    return (
        f"{metal_identity} does not yet have a usable conceptual source-family mapping. "
        f"{baseline_note} {tuned_note}"
    )


def build_metal_parameter_records(
    *,
    registry: tuple[MetalModelRegistryRecord, ...],
    manifest_rows: tuple[MDSystemBuildManifestRow, ...],
    numeric_parameter_registry: tuple[NumericParameterValueRecord, ...],
) -> tuple[tuple[ParameterSourceFamilyRecord, ...], tuple[MetalParameterRecord, ...]]:
    """Translate the Phase 6A2 registry into Phase 6A2c family mappings and readiness state."""
    source_family_registry = build_parameter_source_family_registry(
        manifest_rows=manifest_rows,
        numeric_parameter_registry=numeric_parameter_registry,
    )
    family_lookup = _parameter_source_families_by_label(source_family_registry)
    metal_records: list[MetalParameterRecord] = []
    for record in registry:
        intended_parameter_family = _translate_intended_parameter_family(record)
        source_family_label = _choose_source_family_label(record.metal_identity)
        support_status = classify_direct_support_status(
            metal_identity=record.metal_identity,
            intended_parameter_family=intended_parameter_family,
            source_family_label=source_family_label,
            source_family_registry=source_family_registry,
        )
        baseline_numeric_values_present = _family_has_numeric_values(
            metal_identity=record.metal_identity,
            source_family_label=GENERIC_12_6_4_HIGHLY_CHARGED,
            water_model=None,
            numeric_parameter_registry=numeric_parameter_registry,
        )
        tuned_numeric_values_present = _family_has_numeric_values(
            metal_identity=record.metal_identity,
            source_family_label=intended_parameter_family,
            water_model=None,
            numeric_parameter_registry=numeric_parameter_registry,
        )
        family = family_lookup[source_family_label]
        metal_records.append(
            MetalParameterRecord(
                metal_identity=record.metal_identity,
                formal_charge=record.formal_charge,
                intended_parameter_family=intended_parameter_family,
                source_family_label=source_family_label,
                direct_support_status=support_status,
                baseline_parameter_family=GENERIC_12_6_4_HIGHLY_CHARGED,
                baseline_numeric_values_present=baseline_numeric_values_present,
                tuned_parameter_family=intended_parameter_family,
                tuned_numeric_values_present=tuned_numeric_values_present,
                water_model_compatibility=family.water_model_compatibility,
                notes=_build_metal_parameter_notes(
                    metal_identity=record.metal_identity,
                    support_status=support_status,
                    source_family_label=source_family_label,
                    baseline_numeric_values_present=baseline_numeric_values_present,
                    tuned_numeric_values_present=tuned_numeric_values_present,
                ),
            )
        )
    return source_family_registry, tuple(metal_records)


def _metal_parameter_records_by_metal(
    records: tuple[MetalParameterRecord, ...],
) -> dict[str, MetalParameterRecord]:
    return {record.metal_identity: record for record in records}


def _ready_for_openmm_system_build_baseline(
    *,
    manifest_row: MDSystemBuildManifestRow,
    metal_record: MetalParameterRecord,
    numeric_parameter_registry: tuple[NumericParameterValueRecord, ...],
) -> bool:
    if not manifest_row.generic_protein_water_preparation_ready:
        return False
    return _family_has_numeric_values(
        metal_identity=manifest_row.target_metal,
        source_family_label=metal_record.baseline_parameter_family,
        water_model=manifest_row.chosen_water_model,
        numeric_parameter_registry=numeric_parameter_registry,
    )


def _ready_for_openmm_system_build_tuned(
    *,
    manifest_row: MDSystemBuildManifestRow,
    metal_record: MetalParameterRecord,
    numeric_parameter_registry: tuple[NumericParameterValueRecord, ...],
) -> bool:
    if not manifest_row.generic_protein_water_preparation_ready:
        return False
    return _family_has_numeric_values(
        metal_identity=manifest_row.target_metal,
        source_family_label=metal_record.tuned_parameter_family,
        water_model=manifest_row.chosen_water_model,
        numeric_parameter_registry=numeric_parameter_registry,
    )


def _build_system_rationale(
    *,
    manifest_row: MDSystemBuildManifestRow,
    metal_record: MetalParameterRecord,
    ready_for_openmm_system_build_baseline: bool,
    ready_for_openmm_system_build_tuned: bool,
) -> str:
    if metal_record.direct_support_status == DIRECT_SUPPORT_STATUS:
        baseline_phrase = (
            "baseline OpenMM builds are now numerically ready"
            if ready_for_openmm_system_build_baseline
            else "baseline numeric handoff is still incomplete"
        )
        tuned_phrase = (
            "explicit tuned coefficients also exist"
            if ready_for_openmm_system_build_tuned
            else "the LanM-tuned refinement path still lacks explicit numeric coefficients"
        )
        return (
            "Phase 6A2b conservatively kept this system custom because the direct "
            f"{metal_record.tuned_parameter_family!r} family was only conceptual at the time. "
            f"Phase 6A2c now records published generic OPC3 12-6-4 coefficients for "
            f"{manifest_row.target_metal}, so {baseline_phrase}, while {tuned_phrase}."
        )
    if metal_record.direct_support_status == PROXY_SUPPORT_STATUS:
        baseline_phrase = (
            "baseline OpenMM builds are now numerically ready"
            if ready_for_openmm_system_build_baseline
            else "baseline numeric handoff is still incomplete"
        )
        tuned_phrase = (
            "explicit tuned coefficients also exist"
            if ready_for_openmm_system_build_tuned
            else "no explicit LanM-tuned coefficients are present in the repo"
        )
        return (
            "Phase 6A2b already mapped this system through the generic 12-6-4 proxy family for "
            f"{manifest_row.target_metal}. Phase 6A2c now records the published OPC3 baseline "
            f"coefficients for that family, so {baseline_phrase}, while {tuned_phrase}."
        )
    return (
        f"{manifest_row.target_metal} still lacks a usable family mapping, so neither baseline nor tuned "
        "OpenMM build readiness can be claimed."
    )


def build_per_system_parameter_mapping(
    *,
    manifest_rows: tuple[MDSystemBuildManifestRow, ...],
    metal_parameter_records: tuple[MetalParameterRecord, ...],
    numeric_parameter_registry: tuple[NumericParameterValueRecord, ...],
) -> tuple[MetalParameterMappingRow, ...]:
    """Map each MD panel system to conceptual and numeric baseline-vs-tuned readiness state."""
    metal_lookup = _metal_parameter_records_by_metal(metal_parameter_records)
    rows: list[MetalParameterMappingRow] = []
    for manifest_row in manifest_rows:
        metal_record = metal_lookup[manifest_row.target_metal]
        ready_for_openmm_system_build_baseline = _ready_for_openmm_system_build_baseline(
            manifest_row=manifest_row,
            metal_record=metal_record,
            numeric_parameter_registry=numeric_parameter_registry,
        )
        ready_for_openmm_system_build_tuned = _ready_for_openmm_system_build_tuned(
            manifest_row=manifest_row,
            metal_record=metal_record,
            numeric_parameter_registry=numeric_parameter_registry,
        )
        rows.append(
            MetalParameterMappingRow(
                panel_member_id=manifest_row.panel_member_id,
                target_metal=manifest_row.target_metal,
                chosen_parameter_family=metal_record.source_family_label,
                direct_support_status=metal_record.direct_support_status,
                baseline_parameter_family=metal_record.baseline_parameter_family,
                baseline_numeric_values_present=_family_has_numeric_values(
                    metal_identity=manifest_row.target_metal,
                    source_family_label=metal_record.baseline_parameter_family,
                    water_model=manifest_row.chosen_water_model,
                    numeric_parameter_registry=numeric_parameter_registry,
                ),
                tuned_parameter_family=metal_record.tuned_parameter_family,
                tuned_numeric_values_present=_family_has_numeric_values(
                    metal_identity=manifest_row.target_metal,
                    source_family_label=metal_record.tuned_parameter_family,
                    water_model=manifest_row.chosen_water_model,
                    numeric_parameter_registry=numeric_parameter_registry,
                ),
                ready_for_openmm_system_build_baseline=ready_for_openmm_system_build_baseline,
                ready_for_openmm_system_build_tuned=ready_for_openmm_system_build_tuned,
                rationale=_build_system_rationale(
                    manifest_row=manifest_row,
                    metal_record=metal_record,
                    ready_for_openmm_system_build_baseline=ready_for_openmm_system_build_baseline,
                    ready_for_openmm_system_build_tuned=ready_for_openmm_system_build_tuned,
                ),
            )
        )
    return tuple(rows)


def _readiness_summary(
    *,
    ready_for_openmm_system_build_baseline: bool,
    ready_for_openmm_system_build_tuned: bool,
) -> str:
    if ready_for_openmm_system_build_baseline and ready_for_openmm_system_build_tuned:
        return "baseline_and_tuned_ready"
    if ready_for_openmm_system_build_baseline:
        return "baseline_ready_tuned_pending"
    if ready_for_openmm_system_build_tuned:
        return "baseline_pending_tuned_ready"
    return "baseline_and_tuned_pending"


def build_metal_build_readiness_rows(
    *,
    manifest_rows: tuple[MDSystemBuildManifestRow, ...],
    system_mapping_rows: tuple[MetalParameterMappingRow, ...],
) -> tuple[MetalBuildReadinessRow, ...]:
    """Join manifest scaffolding state with numeric family readiness for audit-friendly output."""
    mapping_lookup = {
        (row.panel_member_id, row.target_metal): row
        for row in system_mapping_rows
    }
    readiness_rows: list[MetalBuildReadinessRow] = []
    for manifest_row in manifest_rows:
        mapping_row = mapping_lookup[(manifest_row.panel_member_id, manifest_row.target_metal)]
        readiness_rows.append(
            MetalBuildReadinessRow(
                panel_rank=manifest_row.panel_rank,
                panel_member_id=manifest_row.panel_member_id,
                target_metal=manifest_row.target_metal,
                topology_class=manifest_row.topology_class,
                build_request_path=manifest_row.build_request_path,
                chosen_water_model=manifest_row.chosen_water_model,
                generic_protein_water_preparation_ready=manifest_row.generic_protein_water_preparation_ready,
                phase_6a2_custom_metal_parameters_still_required=(
                    manifest_row.custom_metal_parameters_still_required
                ),
                baseline_parameter_family=mapping_row.baseline_parameter_family,
                baseline_numeric_values_present=mapping_row.baseline_numeric_values_present,
                ready_for_openmm_system_build_baseline=(
                    mapping_row.ready_for_openmm_system_build_baseline
                ),
                tuned_parameter_family=mapping_row.tuned_parameter_family,
                tuned_numeric_values_present=mapping_row.tuned_numeric_values_present,
                ready_for_openmm_system_build_tuned=mapping_row.ready_for_openmm_system_build_tuned,
                readiness_summary=_readiness_summary(
                    ready_for_openmm_system_build_baseline=(
                        mapping_row.ready_for_openmm_system_build_baseline
                    ),
                    ready_for_openmm_system_build_tuned=(
                        mapping_row.ready_for_openmm_system_build_tuned
                    ),
                ),
            )
        )
    return tuple(readiness_rows)


def render_metal_parameter_config_payload(
    *,
    source_family_registry: tuple[ParameterSourceFamilyRecord, ...],
    metal_parameter_records: tuple[MetalParameterRecord, ...],
    metal_model_registry_path: Path,
    md_system_build_manifest_path: Path,
    metal_parameter_values_path: Path,
) -> dict[str, Any]:
    """Render the Phase 6A2c YAML payload."""
    return {
        "version": 1,
        "phase": "6A2c",
        "description": (
            "Deterministic metal-parameter family mapping plus audited baseline-versus-tuned "
            "numeric readiness for the MD validation panel."
        ),
        "source_artifacts": {
            "metal_model_registry": _display_path(metal_model_registry_path),
            "md_system_build_manifest": _display_path(md_system_build_manifest_path),
            "metal_parameter_values": _display_path(metal_parameter_values_path),
        },
        "parameter_source_families": [
            {
                "source_family_label": family.source_family_label,
                "family_scope": family.family_scope,
                "covered_metals": list(family.covered_metals),
                "water_model_compatibility": list(family.water_model_compatibility),
                "repo_numeric_parameters_present": family.repo_numeric_parameters_present,
                "parameter_provenance": family.parameter_provenance,
                "notes": family.notes,
            }
            for family in source_family_registry
        ],
        "metals": [
            {
                "metal_identity": record.metal_identity,
                "formal_charge": record.formal_charge,
                "intended_parameter_family": record.intended_parameter_family,
                "source_family_label": record.source_family_label,
                "direct_support_status": record.direct_support_status,
                "baseline_parameter_family": record.baseline_parameter_family,
                "baseline_numeric_values_present": record.baseline_numeric_values_present,
                "tuned_parameter_family": record.tuned_parameter_family,
                "tuned_numeric_values_present": record.tuned_numeric_values_present,
                "water_model_compatibility": list(record.water_model_compatibility),
                "notes": record.notes,
            }
            for record in metal_parameter_records
        ],
    }


def render_metal_parameter_strategy_markdown(
    *,
    source_family_registry: tuple[ParameterSourceFamilyRecord, ...],
    metal_parameter_records: tuple[MetalParameterRecord, ...],
    numeric_parameter_registry: tuple[NumericParameterValueRecord, ...],
    system_mapping_rows: tuple[MetalParameterMappingRow, ...],
    readiness_rows: tuple[MetalBuildReadinessRow, ...],
    manifest_rows: tuple[MDSystemBuildManifestRow, ...],
) -> str:
    """Render the Phase 6A2c report."""
    baseline_ready_rows = [
        row for row in readiness_rows if row.ready_for_openmm_system_build_baseline
    ]
    tuned_ready_rows = [
        row for row in readiness_rows if row.ready_for_openmm_system_build_tuned
    ]
    custom_required_rows = [
        row for row in manifest_rows if row.custom_metal_parameters_still_required
    ]
    baseline_ready_metals = [
        record.metal_identity for record in metal_parameter_records if record.baseline_numeric_values_present
    ]
    conceptual_lanm_metals = [
        record.metal_identity
        for record in metal_parameter_records
        if record.source_family_label == CHELATOR_TUNED_12_6_4_LANMODULIN
    ]
    proxy_metals = [
        record.metal_identity
        for record in metal_parameter_records
        if record.direct_support_status == PROXY_SUPPORT_STATUS
    ]
    lines = [
        "# Metal Parameter Strategy",
        "",
        "Phase 6A2c preserves the Phase 6A2b conceptual family mapping, ingests published OPC3 12-6-4 baseline coefficients, and updates per-system baseline-versus-tuned build readiness.",
        "It does not build OpenMM `System` objects, run minimization, or run MD.",
        "",
        "## Why Phase 6A2b Was Conservative",
        "",
        f"- All `{len(custom_required_rows)}` Phase 6A2 panel systems still carry `custom_metal_parameters_still_required=True` in `results/tables/md_system_build_manifest.csv` because that artifact predated any audited numeric metal registry.",
        "- Phase 6A2b intentionally stopped at conceptual family mapping in `config/metal_parameters.yaml`; it did not ingest explicit Dy/Nd/Y/Al/Fe coefficient sets into the repo.",
        "- With only family labels and no auditable coefficient registry, the conservative and correct status in Phase 6A2b was to keep every system out of build-ready state.",
        "",
        "## Published OPC3 Baseline Registry",
        "",
        f"- Published generic OPC3 12-6-4 baseline coefficients are now present for `{', '.join(baseline_ready_metals)}` in `config/metal_parameter_values.yaml`.",
        f"- The conceptual LanM-tuned family is still the direct refinement path for `{', '.join(conceptual_lanm_metals)}`, while `{', '.join(proxy_metals)}` remain conceptually generic/proxy metals.",
        "- Result: the generic baseline family is now numerically build-ready for Dy, Nd, Y, Al, and Fe under OPC3, while the LanM-tuned family remains a future refinement path until explicit tuned coefficients are added.",
        "",
        "| metal | source_family_label | water_model | rmin_half_A | epsilon_kcal_per_mol | c4_kcal_per_mol_A4 |",
        "| --- | --- | --- | ---: | ---: | ---: |",
    ]
    for record in numeric_parameter_registry:
        lines.append(
            "| "
            f"{record.metal_identity} | "
            f"{record.source_family_label} | "
            f"{record.water_model} | "
            f"{record.rmin_half_A:.3f} | "
            f"{record.epsilon_kcal_per_mol:.8f} | "
            f"{record.c4_kcal_per_mol_A4:d} |"
        )
    lines.extend(
        [
            "",
            "## Family Status",
            "",
            "| source_family_label | family_scope | covered_metals | repo_numeric_parameters_present | parameter_provenance | water_model_compatibility |",
            "| --- | --- | --- | --- | --- | --- |",
        ]
    )
    for family in source_family_registry:
        lines.append(
            "| "
            f"{family.source_family_label} | "
            f"{family.family_scope} | "
            f"{', '.join(family.covered_metals)} | "
            f"{family.repo_numeric_parameters_present} | "
            f"{family.parameter_provenance} | "
            f"{', '.join(family.water_model_compatibility)} |"
        )
    lines.extend(
        [
            "",
            "## Per-Metal Mapping",
            "",
            "| metal | intended_parameter_family | conceptual_source_family | direct_support_status | baseline_parameter_family | baseline_numeric_values_present | tuned_parameter_family | tuned_numeric_values_present |",
            "| --- | --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for record in metal_parameter_records:
        lines.append(
            "| "
            f"{record.metal_identity} | "
            f"{record.intended_parameter_family} | "
            f"{record.source_family_label} | "
            f"{record.direct_support_status} | "
            f"{record.baseline_parameter_family} | "
            f"{record.baseline_numeric_values_present} | "
            f"{record.tuned_parameter_family} | "
            f"{record.tuned_numeric_values_present} |"
        )
    lines.extend(
        [
            "",
            "## Per-System Readiness",
            "",
            f"- systems mapped: `{len(system_mapping_rows)}`",
            f"- systems baseline-ready for OpenMM `System` build: `{len(baseline_ready_rows)}`",
            f"- systems tuned-ready for OpenMM `System` build: `{len(tuned_ready_rows)}`",
            f"- per-system readiness audit written to `results/tables/metal_build_readiness.csv`.",
            "",
            "| panel_member_id | target_metal | chosen_parameter_family | baseline_numeric_values_present | ready_for_openmm_system_build_baseline | tuned_numeric_values_present | ready_for_openmm_system_build_tuned | readiness_summary |",
            "| --- | --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    readiness_lookup = {
        (row.panel_member_id, row.target_metal): row
        for row in readiness_rows
    }
    for row in system_mapping_rows:
        readiness_row = readiness_lookup[(row.panel_member_id, row.target_metal)]
        lines.append(
            "| "
            f"{row.panel_member_id} | "
            f"{row.target_metal} | "
            f"{row.chosen_parameter_family} | "
            f"{row.baseline_numeric_values_present} | "
            f"{row.ready_for_openmm_system_build_baseline} | "
            f"{row.tuned_numeric_values_present} | "
            f"{row.ready_for_openmm_system_build_tuned} | "
            f"{readiness_row.readiness_summary} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            f"- The published generic OPC3 12-6-4 baseline family `{GENERIC_12_6_4_HIGHLY_CHARGED}` is now numerically ready for all panel metals under `{DEFAULT_OPC3_WATER_MODEL}`.",
            f"- The LanM-tuned family `{CHELATOR_TUNED_12_6_4_LANMODULIN}` remains intentionally not build-ready because `{EXPLICIT_TUNED_NUMERIC_COEFFICIENTS}` are absent from the repo.",
            f"- The generic baseline values come from `{PUBLISHED_OPC3_12_6_4_BASELINE}` and are being handed off as the auditable baseline path; tuned coefficients remain a future refinement path rather than a claimed ready state.",
            "",
        ]
    )
    return "\n".join(lines)


def prepare_metal_parameter_mapping(
    *,
    metal_model_registry_path: Path,
    md_system_build_manifest_path: Path,
    metal_parameter_config_path: Path,
    metal_parameter_values_path: Path,
    metal_parameter_mapping_path: Path,
    metal_build_readiness_path: Path,
    metal_parameter_strategy_report_path: Path,
) -> tuple[MetalParameterMappingRow, ...]:
    """Generate the deterministic Phase 6A2c metal-parameter mapping outputs."""
    registry = load_metal_model_registry(metal_model_registry_path)
    manifest_rows = load_md_system_build_manifest(md_system_build_manifest_path)
    numeric_parameter_registry = build_default_numeric_parameter_registry()
    source_family_registry, metal_parameter_records = build_metal_parameter_records(
        registry=registry,
        manifest_rows=manifest_rows,
        numeric_parameter_registry=numeric_parameter_registry,
    )
    system_mapping_rows = build_per_system_parameter_mapping(
        manifest_rows=manifest_rows,
        metal_parameter_records=metal_parameter_records,
        numeric_parameter_registry=numeric_parameter_registry,
    )
    readiness_rows = build_metal_build_readiness_rows(
        manifest_rows=manifest_rows,
        system_mapping_rows=system_mapping_rows,
    )
    write_yaml(
        metal_parameter_config_path,
        render_metal_parameter_config_payload(
            source_family_registry=source_family_registry,
            metal_parameter_records=metal_parameter_records,
            metal_model_registry_path=metal_model_registry_path,
            md_system_build_manifest_path=md_system_build_manifest_path,
            metal_parameter_values_path=metal_parameter_values_path,
        ),
    )
    write_yaml(
        metal_parameter_values_path,
        render_metal_parameter_values_payload(
            numeric_parameter_registry=numeric_parameter_registry,
            metal_model_registry_path=metal_model_registry_path,
            md_system_build_manifest_path=md_system_build_manifest_path,
        ),
    )
    write_csv_rows(metal_parameter_mapping_path, system_mapping_rows)
    write_csv_rows(metal_build_readiness_path, readiness_rows)
    atomic_write_text(
        metal_parameter_strategy_report_path,
        render_metal_parameter_strategy_markdown(
            source_family_registry=source_family_registry,
            metal_parameter_records=metal_parameter_records,
            numeric_parameter_registry=numeric_parameter_registry,
            system_mapping_rows=system_mapping_rows,
            readiness_rows=readiness_rows,
            manifest_rows=manifest_rows,
        ),
    )
    return system_mapping_rows
