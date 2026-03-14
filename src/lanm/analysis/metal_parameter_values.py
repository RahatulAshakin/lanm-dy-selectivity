"""Numeric metal-parameter registries for Phase 6A2c."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from lanm.filesystem import write_yaml
from lanm.paths import REPO_ROOT

TARGET_METALS = ("Dy", "Nd", "Y", "Al", "Fe")
GENERIC_12_6_4_HIGHLY_CHARGED = "generic_12_6_4_highly_charged"
CHELATOR_TUNED_12_6_4_LANMODULIN = "chelator_tuned_12_6_4_lanmodulin"
LANM_ADJACENT_DIRECT_METALS = ("Dy", "Nd", "Y")
DEFAULT_OPC3_WATER_MODEL = "amber19/opc3.xml"
PUBLISHED_OPC3_12_6_4_BASELINE = "published_opc3_12_6_4_baseline"
EXPLICIT_TUNED_NUMERIC_COEFFICIENTS = "explicit_tuned_numeric_coefficients_in_repo"


@dataclass(frozen=True, slots=True)
class NumericParameterValueRecord:
    metal_identity: str
    source_family_label: str
    water_model: str
    formal_charge: int
    rmin_half_A: float
    epsilon_kcal_per_mol: float
    c4_kcal_per_mol_A4: int
    parameter_provenance: str
    notes: str


@dataclass(frozen=True, slots=True)
class NumericParameterFamilyRecord:
    source_family_label: str
    family_scope: str
    covered_metals: tuple[str, ...]
    water_model_compatibility: tuple[str, ...]
    numeric_values_present: bool
    parameter_provenance: str
    notes: str


def _display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def build_default_numeric_parameter_registry() -> tuple[NumericParameterValueRecord, ...]:
    """Return the published OPC3 12-6-4 baseline coefficients for the panel metals."""
    registry = (
        NumericParameterValueRecord(
            metal_identity="Dy",
            source_family_label=GENERIC_12_6_4_HIGHLY_CHARGED,
            water_model=DEFAULT_OPC3_WATER_MODEL,
            formal_charge=3,
            rmin_half_A=1.632,
            epsilon_kcal_per_mol=0.09620220,
            c4_kcal_per_mol_A4=183,
            parameter_provenance=PUBLISHED_OPC3_12_6_4_BASELINE,
            notes="Published OPC3 12-6-4 baseline coefficients ingested exactly for Dy3+.",
        ),
        NumericParameterValueRecord(
            metal_identity="Nd",
            source_family_label=GENERIC_12_6_4_HIGHLY_CHARGED,
            water_model=DEFAULT_OPC3_WATER_MODEL,
            formal_charge=3,
            rmin_half_A=1.712,
            epsilon_kcal_per_mol=0.14640930,
            c4_kcal_per_mol_A4=184,
            parameter_provenance=PUBLISHED_OPC3_12_6_4_BASELINE,
            notes="Published OPC3 12-6-4 baseline coefficients ingested exactly for Nd3+.",
        ),
        NumericParameterValueRecord(
            metal_identity="Y",
            source_family_label=GENERIC_12_6_4_HIGHLY_CHARGED,
            water_model=DEFAULT_OPC3_WATER_MODEL,
            formal_charge=3,
            rmin_half_A=1.626,
            epsilon_kcal_per_mol=0.09289608,
            c4_kcal_per_mol_A4=192,
            parameter_provenance=PUBLISHED_OPC3_12_6_4_BASELINE,
            notes="Published OPC3 12-6-4 baseline coefficients ingested exactly for Y3+.",
        ),
        NumericParameterValueRecord(
            metal_identity="Al",
            source_family_label=GENERIC_12_6_4_HIGHLY_CHARGED,
            water_model=DEFAULT_OPC3_WATER_MODEL,
            formal_charge=3,
            rmin_half_A=1.361,
            epsilon_kcal_per_mol=0.01031847,
            c4_kcal_per_mol_A4=363,
            parameter_provenance=PUBLISHED_OPC3_12_6_4_BASELINE,
            notes="Published OPC3 12-6-4 baseline coefficients ingested exactly for Al3+.",
        ),
        NumericParameterValueRecord(
            metal_identity="Fe",
            source_family_label=GENERIC_12_6_4_HIGHLY_CHARGED,
            water_model=DEFAULT_OPC3_WATER_MODEL,
            formal_charge=3,
            rmin_half_A=1.455,
            epsilon_kcal_per_mol=0.02662782,
            c4_kcal_per_mol_A4=429,
            parameter_provenance=PUBLISHED_OPC3_12_6_4_BASELINE,
            notes="Published OPC3 12-6-4 baseline coefficients ingested exactly for Fe3+.",
        ),
    )
    return validate_numeric_parameter_registry(registry)


def validate_numeric_parameter_registry(
    registry: tuple[NumericParameterValueRecord, ...],
) -> tuple[NumericParameterValueRecord, ...]:
    """Validate the numeric parameter registry before it is written or consumed."""
    if not registry:
        raise ValueError("Numeric parameter registry is empty")
    observed_keys: set[tuple[str, str, str]] = set()
    observed_generic_metals: list[str] = []
    for record in registry:
        key = (record.source_family_label, record.metal_identity, record.water_model)
        if key in observed_keys:
            raise ValueError(f"Duplicate numeric parameter entry for {key}")
        observed_keys.add(key)
        if record.metal_identity not in TARGET_METALS:
            raise ValueError(f"Unsupported panel metal in numeric parameter registry: {record.metal_identity}")
        if record.source_family_label not in {
            GENERIC_12_6_4_HIGHLY_CHARGED,
            CHELATOR_TUNED_12_6_4_LANMODULIN,
        }:
            raise ValueError(
                f"Unsupported parameter family in numeric parameter registry: {record.source_family_label}"
            )
        if record.formal_charge != 3:
            raise ValueError(
                f"Numeric parameter entry for {record.metal_identity} must keep formal_charge 3"
            )
        if record.rmin_half_A <= 0.0:
            raise ValueError(f"Numeric parameter entry for {record.metal_identity} must have rmin_half_A > 0")
        if record.epsilon_kcal_per_mol <= 0.0:
            raise ValueError(
                f"Numeric parameter entry for {record.metal_identity} must have epsilon_kcal_per_mol > 0"
            )
        if record.c4_kcal_per_mol_A4 <= 0:
            raise ValueError(
                f"Numeric parameter entry for {record.metal_identity} must have c4_kcal_per_mol_A4 > 0"
            )
        if not record.parameter_provenance:
            raise ValueError(
                f"Numeric parameter entry for {record.metal_identity} is missing parameter_provenance"
            )
        if not record.notes:
            raise ValueError(f"Numeric parameter entry for {record.metal_identity} is missing notes")
        if record.source_family_label == GENERIC_12_6_4_HIGHLY_CHARGED:
            if record.parameter_provenance != PUBLISHED_OPC3_12_6_4_BASELINE:
                raise ValueError(
                    f"Generic baseline entry for {record.metal_identity} must use "
                    f"{PUBLISHED_OPC3_12_6_4_BASELINE!r}"
                )
            if record.water_model != DEFAULT_OPC3_WATER_MODEL:
                raise ValueError(
                    f"Generic baseline entry for {record.metal_identity} must use "
                    f"{DEFAULT_OPC3_WATER_MODEL!r}"
                )
            observed_generic_metals.append(record.metal_identity)
    if tuple(observed_generic_metals) != TARGET_METALS:
        raise ValueError(
            "Published OPC3 baseline registry must include generic 12-6-4 entries for "
            f"{TARGET_METALS}, found {tuple(observed_generic_metals)}"
        )
    return registry


def build_numeric_parameter_family_registry(
    numeric_parameter_registry: tuple[NumericParameterValueRecord, ...],
) -> tuple[NumericParameterFamilyRecord, ...]:
    """Summarize family-level numeric availability for audited readiness logic."""
    validate_numeric_parameter_registry(numeric_parameter_registry)
    generic_metals = tuple(
        record.metal_identity
        for record in numeric_parameter_registry
        if record.source_family_label == GENERIC_12_6_4_HIGHLY_CHARGED
    )
    generic_water_models = tuple(
        sorted(
            {
                record.water_model
                for record in numeric_parameter_registry
                if record.source_family_label == GENERIC_12_6_4_HIGHLY_CHARGED
            }
        )
    )
    tuned_metals = tuple(
        record.metal_identity
        for record in numeric_parameter_registry
        if record.source_family_label == CHELATOR_TUNED_12_6_4_LANMODULIN
    )
    tuned_water_models = tuple(
        sorted(
            {
                record.water_model
                for record in numeric_parameter_registry
                if record.source_family_label == CHELATOR_TUNED_12_6_4_LANMODULIN
            }
        )
    )
    return (
        NumericParameterFamilyRecord(
            source_family_label=GENERIC_12_6_4_HIGHLY_CHARGED,
            family_scope="generic_baseline",
            covered_metals=generic_metals,
            water_model_compatibility=generic_water_models,
            numeric_values_present=bool(generic_metals),
            parameter_provenance=PUBLISHED_OPC3_12_6_4_BASELINE,
            notes=(
                "Published OPC3 12-6-4 baseline coefficients are present in the repo for "
                "Dy, Nd, Y, Al, and Fe, so this family is numerically auditable for baseline builds."
            ),
        ),
        NumericParameterFamilyRecord(
            source_family_label=CHELATOR_TUNED_12_6_4_LANMODULIN,
            family_scope="lanm_adjacent_chelator_tuned",
            covered_metals=tuned_metals if tuned_metals else LANM_ADJACENT_DIRECT_METALS,
            water_model_compatibility=tuned_water_models if tuned_water_models else (DEFAULT_OPC3_WATER_MODEL,),
            numeric_values_present=bool(tuned_metals),
            parameter_provenance=EXPLICIT_TUNED_NUMERIC_COEFFICIENTS,
            notes=(
                "Lanmodulin-adjacent chelator-tuned coefficients remain a future refinement path. "
                "No explicit tuned numeric coefficients are present in the repo yet."
            ),
        ),
    )


def render_metal_parameter_values_payload(
    *,
    numeric_parameter_registry: tuple[NumericParameterValueRecord, ...],
    metal_model_registry_path: Path,
    md_system_build_manifest_path: Path,
) -> dict[str, Any]:
    """Render the numeric parameter YAML payload."""
    family_registry = build_numeric_parameter_family_registry(numeric_parameter_registry)
    return {
        "version": 1,
        "phase": "6A2c",
        "description": (
            "Auditable numeric metal-parameter registry for published OPC3 12-6-4 baseline "
            "handoff and tuned-family readiness tracking."
        ),
        "source_artifacts": {
            "metal_model_registry": _display_path(metal_model_registry_path),
            "md_system_build_manifest": _display_path(md_system_build_manifest_path),
        },
        "parameter_value_families": [
            {
                "source_family_label": family.source_family_label,
                "family_scope": family.family_scope,
                "covered_metals": list(family.covered_metals),
                "water_model_compatibility": list(family.water_model_compatibility),
                "numeric_values_present": family.numeric_values_present,
                "parameter_provenance": family.parameter_provenance,
                "notes": family.notes,
            }
            for family in family_registry
        ],
        "numeric_parameter_values": [
            {
                "metal_identity": record.metal_identity,
                "source_family_label": record.source_family_label,
                "water_model": record.water_model,
                "formal_charge": record.formal_charge,
                "rmin_half_A": record.rmin_half_A,
                "epsilon_kcal_per_mol": record.epsilon_kcal_per_mol,
                "c4_kcal_per_mol_A4": record.c4_kcal_per_mol_A4,
                "parameter_provenance": record.parameter_provenance,
                "notes": record.notes,
            }
            for record in numeric_parameter_registry
        ],
    }


def write_metal_parameter_values_yaml(
    *,
    path: Path,
    numeric_parameter_registry: tuple[NumericParameterValueRecord, ...],
    metal_model_registry_path: Path,
    md_system_build_manifest_path: Path,
) -> None:
    """Write the numeric parameter registry as deterministic YAML."""
    write_yaml(
        path,
        render_metal_parameter_values_payload(
            numeric_parameter_registry=numeric_parameter_registry,
            metal_model_registry_path=metal_model_registry_path,
            md_system_build_manifest_path=md_system_build_manifest_path,
        ),
    )
