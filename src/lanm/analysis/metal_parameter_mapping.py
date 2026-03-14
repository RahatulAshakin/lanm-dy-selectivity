"""Deterministic Phase 6A2b metal-parameter source mapping for the MD panel."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_text, write_csv_rows, write_yaml
from lanm.paths import REPO_ROOT

TARGET_METALS = ("Dy", "Nd", "Y", "Al", "Fe")
PHASE_6A2_REGISTRY_FAMILY = "custom_bound_site_12_6_4_lj_chelator_tuned"
PHASE_6A2_PARAMETER_SOURCE_STATUS = "custom_parameters_not_present_in_repo"
GENERIC_12_6_4_HIGHLY_CHARGED = "generic_12_6_4_highly_charged"
CHELATOR_TUNED_12_6_4_LANMODULIN = "chelator_tuned_12_6_4_lanmodulin"
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
    notes: str


@dataclass(frozen=True, slots=True)
class MetalParameterRecord:
    metal_identity: str
    formal_charge: int
    intended_parameter_family: str
    source_family_label: str
    direct_support_status: str
    water_model_compatibility: tuple[str, ...]
    notes: str


@dataclass(frozen=True, slots=True)
class MetalParameterMappingRow:
    panel_member_id: str
    target_metal: str
    chosen_parameter_family: str
    direct_support_status: str
    ready_for_openmm_system_build: bool
    rationale: str


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
    sorted_rows = tuple(sorted(rows, key=lambda row: (row.panel_rank, row.panel_member_id, TARGET_METALS.index(row.target_metal))))
    observed_metals = tuple(sorted({row.target_metal for row in sorted_rows}, key=TARGET_METALS.index))
    if observed_metals != TARGET_METALS:
        raise ValueError(f"Expected manifest metals {TARGET_METALS}, found {observed_metals}")
    for row in sorted_rows:
        require_path(_repo_path(row.source_starting_structure_pdb))
        require_path(_repo_path(row.build_request_path))
    return sorted_rows


def build_parameter_source_family_registry(
    manifest_rows: tuple[MDSystemBuildManifestRow, ...],
) -> tuple[ParameterSourceFamilyRecord, ...]:
    """Build the deterministic conceptual source-family registry for Phase 6A2b."""
    water_models = tuple(sorted({row.chosen_water_model for row in manifest_rows}))
    if not water_models:
        raise ValueError("No water models observed in the build manifest")
    return (
        ParameterSourceFamilyRecord(
            source_family_label=GENERIC_12_6_4_HIGHLY_CHARGED,
            family_scope="generic_baseline",
            covered_metals=GENERIC_BASELINE_ONLY_METALS,
            water_model_compatibility=water_models,
            repo_numeric_parameters_present=False,
            notes=(
                "Conceptual generic 12-6-4 baseline family for highly charged trivalent competitors. "
                "The repo does not currently store numeric Al or Fe bound-site parameters in this family."
            ),
        ),
        ParameterSourceFamilyRecord(
            source_family_label=CHELATOR_TUNED_12_6_4_LANMODULIN,
            family_scope="lanm_adjacent_chelator_tuned",
            covered_metals=LANM_ADJACENT_DIRECT_METALS,
            water_model_compatibility=water_models,
            repo_numeric_parameters_present=False,
            notes=(
                "Conceptual LanM-adjacent chelator-tuned 12-6-4 family for bound-site lanthanide models. "
                "The repo does not currently store numeric Dy, Nd, or Y bound-site parameters in this family."
            ),
        ),
    )


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


def _build_metal_parameter_notes(
    *,
    metal_identity: str,
    support_status: str,
    source_family_label: str,
) -> str:
    if support_status == DIRECT_SUPPORT_STATUS:
        return (
            f"{metal_identity} maps directly onto the LanM-adjacent chelator-tuned source family "
            f"{source_family_label!r}, but the repo still lacks numeric parameters for this metal. "
            "Parameter derivation is still required before OpenMM System creation."
        )
    if support_status == PROXY_SUPPORT_STATUS:
        return (
            f"{metal_identity} currently maps only to the generic highly charged 12-6-4 baseline "
            f"{source_family_label!r}. Any OpenMM use would require an explicit proxy decision, and the "
            "repo still lacks numeric bound-site parameters."
        )
    return (
        f"{metal_identity} does not yet have a usable source family mapping. Parameter derivation remains "
        "required before OpenMM System creation."
    )


def build_metal_parameter_records(
    *,
    registry: tuple[MetalModelRegistryRecord, ...],
    manifest_rows: tuple[MDSystemBuildManifestRow, ...],
) -> tuple[tuple[ParameterSourceFamilyRecord, ...], tuple[MetalParameterRecord, ...]]:
    """Translate the Phase 6A2 registry into Phase 6A2b metal-parameter mappings."""
    source_family_registry = build_parameter_source_family_registry(manifest_rows)
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
        family = family_lookup[source_family_label]
        metal_records.append(
            MetalParameterRecord(
                metal_identity=record.metal_identity,
                formal_charge=record.formal_charge,
                intended_parameter_family=intended_parameter_family,
                source_family_label=source_family_label,
                direct_support_status=support_status,
                water_model_compatibility=family.water_model_compatibility,
                notes=_build_metal_parameter_notes(
                    metal_identity=record.metal_identity,
                    support_status=support_status,
                    source_family_label=source_family_label,
                ),
            )
        )
    return source_family_registry, tuple(metal_records)


def _metal_parameter_records_by_metal(
    records: tuple[MetalParameterRecord, ...],
) -> dict[str, MetalParameterRecord]:
    return {record.metal_identity: record for record in records}


def _ready_for_openmm_system_build(
    *,
    manifest_row: MDSystemBuildManifestRow,
    metal_record: MetalParameterRecord,
    source_family_registry: tuple[ParameterSourceFamilyRecord, ...],
) -> bool:
    family = _parameter_source_families_by_label(source_family_registry).get(metal_record.source_family_label)
    if family is None:
        return False
    if manifest_row.custom_metal_parameters_still_required:
        return False
    if metal_record.direct_support_status != DIRECT_SUPPORT_STATUS:
        return False
    return family.repo_numeric_parameters_present


def _build_system_rationale(
    *,
    manifest_row: MDSystemBuildManifestRow,
    metal_record: MetalParameterRecord,
) -> str:
    if metal_record.direct_support_status == DIRECT_SUPPORT_STATUS:
        return (
            "Phase 6A2 kept this system in the custom-parameter-required state because the chosen "
            f"{metal_record.source_family_label!r} family is only a conceptual LanM-adjacent mapping in "
            "the repo; numeric parameters still need to be derived."
        )
    if metal_record.direct_support_status == PROXY_SUPPORT_STATUS:
        return (
            "Phase 6A2 kept this system in the custom-parameter-required state because the current "
            f"mapping for {manifest_row.target_metal} is a generic 12-6-4 proxy rather than a direct "
            "LanM-tuned family. An explicit proxy decision plus numeric parameter handoff is still needed."
        )
    return (
        f"Phase 6A2 kept this system in the custom-parameter-required state because {manifest_row.target_metal} "
        "still lacks a mapped source family. Parameter derivation is still required."
    )


def build_per_system_parameter_mapping(
    *,
    manifest_rows: tuple[MDSystemBuildManifestRow, ...],
    metal_parameter_records: tuple[MetalParameterRecord, ...],
    source_family_registry: tuple[ParameterSourceFamilyRecord, ...],
) -> tuple[MetalParameterMappingRow, ...]:
    """Map each MD panel system to its chosen parameter family and readiness status."""
    metal_lookup = _metal_parameter_records_by_metal(metal_parameter_records)
    rows: list[MetalParameterMappingRow] = []
    for manifest_row in manifest_rows:
        metal_record = metal_lookup[manifest_row.target_metal]
        rows.append(
            MetalParameterMappingRow(
                panel_member_id=manifest_row.panel_member_id,
                target_metal=manifest_row.target_metal,
                chosen_parameter_family=metal_record.source_family_label,
                direct_support_status=metal_record.direct_support_status,
                ready_for_openmm_system_build=_ready_for_openmm_system_build(
                    manifest_row=manifest_row,
                    metal_record=metal_record,
                    source_family_registry=source_family_registry,
                ),
                rationale=_build_system_rationale(
                    manifest_row=manifest_row,
                    metal_record=metal_record,
                ),
            )
        )
    return tuple(rows)


def render_metal_parameter_config_payload(
    *,
    source_family_registry: tuple[ParameterSourceFamilyRecord, ...],
    metal_parameter_records: tuple[MetalParameterRecord, ...],
    metal_model_registry_path: Path,
    md_system_build_manifest_path: Path,
) -> dict[str, Any]:
    """Render the Phase 6A2b YAML payload."""
    return {
        "version": 1,
        "phase": "6A2b",
        "description": "Deterministic metal-parameter source mapping for the MD validation panel.",
        "source_artifacts": {
            "metal_model_registry": _display_path(metal_model_registry_path),
            "md_system_build_manifest": _display_path(md_system_build_manifest_path),
        },
        "parameter_source_families": [
            {
                "source_family_label": family.source_family_label,
                "family_scope": family.family_scope,
                "covered_metals": list(family.covered_metals),
                "water_model_compatibility": list(family.water_model_compatibility),
                "repo_numeric_parameters_present": family.repo_numeric_parameters_present,
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
    system_mapping_rows: tuple[MetalParameterMappingRow, ...],
    manifest_rows: tuple[MDSystemBuildManifestRow, ...],
) -> str:
    """Render the Phase 6A2b report."""
    direct_rows = [row for row in system_mapping_rows if row.direct_support_status == DIRECT_SUPPORT_STATUS]
    proxy_rows = [row for row in system_mapping_rows if row.direct_support_status == PROXY_SUPPORT_STATUS]
    ready_rows = [row for row in system_mapping_rows if row.ready_for_openmm_system_build]
    custom_required_rows = [row for row in manifest_rows if row.custom_metal_parameters_still_required]
    direct_metals = [record.metal_identity for record in metal_parameter_records if record.direct_support_status == DIRECT_SUPPORT_STATUS]
    proxy_metals = [record.metal_identity for record in metal_parameter_records if record.direct_support_status == PROXY_SUPPORT_STATUS]
    derivation_metals = [
        record.metal_identity
        for record in metal_parameter_records
        if record.direct_support_status in {DIRECT_SUPPORT_STATUS, NEEDS_DERIVATION_STATUS}
    ]
    lines = [
        "# Metal Parameter Strategy",
        "",
        "Phase 6A2b converts the Phase 6A2 metal-model registry and system-build manifest into an explicit parameter-source mapping.",
        "It does not create numeric metal parameters and it does not build OpenMM `System` objects.",
        "",
        "## Why Phase 6A2 Marked Every System Custom",
        "",
        f"- All `{len(custom_required_rows)}` Phase 6A2 panel systems still carry `custom_metal_parameters_still_required=True` in `results/tables/md_system_build_manifest.csv`.",
        "- Dy, Nd, Y, Al, and Fe are modeled as pre-bound LanM-site metals rather than ordinary bulk-solvent ions.",
        "- Phase 6A2 recorded only an intended custom bound-site family in `config/metal_models.yaml`; it did not add numeric Dy/Nd/Y/Al/Fe parameter files to the repo.",
        "- Result: the protein/water scaffolding is ready to audit, but full OpenMM `System` creation remains blocked until the metal-parameter source is turned into actual numeric inputs.",
        "",
        "## Source Family Overview",
        "",
        f"- Metals with only generic baseline support: `{', '.join(proxy_metals)}`",
        f"- Metals with LanM-adjacent chelator-tuned support: `{', '.join(direct_metals)}`",
        f"- Metals that still require derivation before OpenMM System creation: `{', '.join(derivation_metals)}`",
        f"- Metals that still require an explicit proxy decision before OpenMM System creation: `{', '.join(proxy_metals)}`",
        "",
        "| source_family_label | family_scope | covered_metals | repo_numeric_parameters_present | water_model_compatibility |",
        "| --- | --- | --- | --- | --- |",
    ]
    for family in source_family_registry:
        lines.append(
            "| "
            f"{family.source_family_label} | "
            f"{family.family_scope} | "
            f"{', '.join(family.covered_metals)} | "
            f"{family.repo_numeric_parameters_present} | "
            f"{', '.join(family.water_model_compatibility)} |"
        )
    lines.extend(
        [
            "",
            "## Per-Metal Mapping",
            "",
            "| metal | formal_charge | intended_parameter_family | source_family_label | direct_support_status | water_model_compatibility |",
            "| --- | ---: | --- | --- | --- | --- |",
        ]
    )
    for record in metal_parameter_records:
        lines.append(
            "| "
            f"{record.metal_identity} | "
            f"{record.formal_charge:+d} | "
            f"{record.intended_parameter_family} | "
            f"{record.source_family_label} | "
            f"{record.direct_support_status} | "
            f"{', '.join(record.water_model_compatibility)} |"
        )
    lines.extend(
        [
            "",
            "## Per-System Mapping",
            "",
            f"- systems mapped: `{len(system_mapping_rows)}`",
            f"- systems with direct family mappings: `{len(direct_rows)}`",
            f"- systems with proxy family mappings: `{len(proxy_rows)}`",
            f"- systems ready for OpenMM `System` build: `{len(ready_rows)}`",
            "",
            "| panel_member_id | target_metal | chosen_parameter_family | direct_support_status | ready_for_openmm_system_build | rationale |",
            "| --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in system_mapping_rows:
        lines.append(
            "| "
            f"{row.panel_member_id} | "
            f"{row.target_metal} | "
            f"{row.chosen_parameter_family} | "
            f"{row.direct_support_status} | "
            f"{row.ready_for_openmm_system_build} | "
            f"{row.rationale} |"
        )
    lines.append("")
    return "\n".join(lines)


def prepare_metal_parameter_mapping(
    *,
    metal_model_registry_path: Path,
    md_system_build_manifest_path: Path,
    metal_parameter_config_path: Path,
    metal_parameter_mapping_path: Path,
    metal_parameter_strategy_report_path: Path,
) -> tuple[MetalParameterMappingRow, ...]:
    """Generate the deterministic Phase 6A2b metal-parameter mapping outputs."""
    registry = load_metal_model_registry(metal_model_registry_path)
    manifest_rows = load_md_system_build_manifest(md_system_build_manifest_path)
    source_family_registry, metal_parameter_records = build_metal_parameter_records(
        registry=registry,
        manifest_rows=manifest_rows,
    )
    system_mapping_rows = build_per_system_parameter_mapping(
        manifest_rows=manifest_rows,
        metal_parameter_records=metal_parameter_records,
        source_family_registry=source_family_registry,
    )
    write_yaml(
        metal_parameter_config_path,
        render_metal_parameter_config_payload(
            source_family_registry=source_family_registry,
            metal_parameter_records=metal_parameter_records,
            metal_model_registry_path=metal_model_registry_path,
            md_system_build_manifest_path=md_system_build_manifest_path,
        ),
    )
    write_csv_rows(metal_parameter_mapping_path, system_mapping_rows)
    atomic_write_text(
        metal_parameter_strategy_report_path,
        render_metal_parameter_strategy_markdown(
            source_family_registry=source_family_registry,
            metal_parameter_records=metal_parameter_records,
            system_mapping_rows=system_mapping_rows,
            manifest_rows=manifest_rows,
        ),
    )
    return system_mapping_rows
