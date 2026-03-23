"""Phase 6A2 deterministic OpenMM system-build scaffolding for the MD panel."""

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
DEFAULT_PROTEIN_FORCEFIELD = "amber19-all.xml"
DEFAULT_WATER_MODEL = "amber19/opc3.xml"
ALLOWED_WATER_MODELS = ("amber19/opc3.xml", "amber19/opc.xml")
DEFAULT_BOX_PADDING_NM = 1.2
DEFAULT_IONIC_STRENGTH_M = 0.15
DEFAULT_NEUTRALIZATION_ION_POLICY = "monovalent_background_ions_only_excluding_panel_metals"
CUSTOM_MODEL_FAMILY = "custom_bound_site_12_6_4_lj_chelator_tuned"
CUSTOM_PARAMETER_SOURCE_STATUS = "custom_parameters_not_present_in_repo"


@dataclass(frozen=True, slots=True)
class MDInputSystemRow:
    panel_rank: int
    panel_member_id: str
    panel_member_type: str
    panel_role: str
    candidate_id: str
    backbone_id: str
    topology_class: str
    target_metal: str
    source_structure: str
    metal_site_template_source: str
    starting_structure_path: str
    system_metadata_path: str
    metal_site_count: int
    preserved_solvent_residue_count: int
    replicate_count: int
    target_temperature_K: int
    suggested_collective_variables: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class OpenMMBuildScaffoldingConfig:
    protein_forcefield: str
    water_model: str
    box_padding_nm: float
    ionic_strength_M: float
    neutralization_ion_policy: str


@dataclass(frozen=True, slots=True)
class TemplateWaterRestorationPlan:
    template_metal_proximal_waters_available: bool
    restore_template_metal_proximal_waters_before_solvation: bool
    comparison_panel_member_id: str | None
    comparison_source: str | None
    comparison_preserved_solvent_residue_count: int


@dataclass(frozen=True, slots=True)
class MetalModelRecord:
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


def _parse_suggested_collective_variables(value: str) -> tuple[str, ...]:
    return tuple(item.strip() for item in value.split(",") if item.strip())


def discover_md_input_manifest(path: Path) -> tuple[MDInputSystemRow, ...]:
    """Load and validate the Phase 6A1 MD input manifest."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            MDInputSystemRow(
                panel_rank=_parse_csv_int(row, "panel_rank"),
                panel_member_id=str(row["panel_member_id"]).strip(),
                panel_member_type=str(row["panel_member_type"]).strip(),
                panel_role=str(row["panel_role"]).strip(),
                candidate_id=str(row["candidate_id"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                topology_class=str(row["topology_class"]).strip(),
                target_metal=str(row["target_metal"]).strip(),
                source_structure=str(row["source_structure"]).strip(),
                metal_site_template_source=str(row["metal_site_template_source"]).strip(),
                starting_structure_path=str(row["starting_structure_path"]).strip(),
                system_metadata_path=str(row["system_metadata_path"]).strip(),
                metal_site_count=_parse_csv_int(row, "metal_site_count"),
                preserved_solvent_residue_count=_parse_csv_int(row, "preserved_solvent_residue_count"),
                replicate_count=_parse_csv_int(row, "replicate_count"),
                target_temperature_K=_parse_csv_int(row, "target_temperature_K"),
                suggested_collective_variables=_parse_suggested_collective_variables(
                    str(row["suggested_collective_variables"]).strip()
                ),
            )
            for row in reader
        ]
    rows = sorted(rows, key=lambda row: (row.panel_rank, row.panel_member_id, row.target_metal))
    if not rows:
        raise ValueError(f"No rows found in {path}")
    observed_metals = tuple(sorted({row.target_metal for row in rows}, key=TARGET_METALS.index))
    if observed_metals != TARGET_METALS:
        raise ValueError(f"Expected MD input metals {TARGET_METALS}, found {observed_metals}")
    for row in rows:
        require_path(_repo_path(row.starting_structure_path))
        require_path(_repo_path(row.system_metadata_path))
    return tuple(rows)


def load_md_protocol_build_config(path: Path) -> OpenMMBuildScaffoldingConfig:
    """Load the Phase 6A1 protocol and normalize Phase 6A2 build defaults."""
    payload = _parse_required_mapping(path)
    protocol_metals = tuple(str(item).strip() for item in payload.get("target_metals", []))
    if protocol_metals != TARGET_METALS:
        raise ValueError(f"Expected protocol target_metals {TARGET_METALS}, found {protocol_metals}")
    defaults = payload.get("system_build_defaults", {})
    if defaults is None:
        defaults = {}
    if not isinstance(defaults, dict):
        raise ValueError(f"Expected mapping in {path} under system_build_defaults")
    protein_forcefield = str(defaults.get("protein_forcefield", DEFAULT_PROTEIN_FORCEFIELD)).strip()
    water_model = str(defaults.get("water_model", DEFAULT_WATER_MODEL)).strip()
    box_padding_nm = float(defaults.get("box_padding_nm", DEFAULT_BOX_PADDING_NM))
    ionic_strength_M = float(defaults.get("ionic_strength_M", DEFAULT_IONIC_STRENGTH_M))
    neutralization_ion_policy = str(
        defaults.get("neutralization_ion_policy", DEFAULT_NEUTRALIZATION_ION_POLICY)
    ).strip()
    if protein_forcefield != DEFAULT_PROTEIN_FORCEFIELD:
        raise ValueError(f"Unsupported protein_forcefield in {path}: {protein_forcefield!r}")
    if water_model not in ALLOWED_WATER_MODELS:
        raise ValueError(
            f"Unsupported water_model in {path}: {water_model!r}; expected one of {ALLOWED_WATER_MODELS}"
        )
    if box_padding_nm <= 0.0:
        raise ValueError(f"box_padding_nm must be positive in {path}")
    if ionic_strength_M < 0.0:
        raise ValueError(f"ionic_strength_M must be non-negative in {path}")
    if not neutralization_ion_policy:
        raise ValueError(f"neutralization_ion_policy must be set in {path}")
    return OpenMMBuildScaffoldingConfig(
        protein_forcefield=protein_forcefield,
        water_model=water_model,
        box_padding_nm=box_padding_nm,
        ionic_strength_M=ionic_strength_M,
        neutralization_ion_policy=neutralization_ion_policy,
    )


def validate_md_input_manifest_against_protocol(
    *,
    manifest_rows: tuple[MDInputSystemRow, ...],
    protocol_path: Path,
) -> None:
    """Validate that the manifest still matches the protocol snapshot."""
    payload = _parse_required_mapping(protocol_path)
    panel_members = payload.get("panel_members")
    if not isinstance(panel_members, list):
        raise ValueError(f"Expected panel_members list in {protocol_path}")
    protocol_rows: dict[tuple[str, str], dict[str, str]] = {}
    for panel_member in panel_members:
        if not isinstance(panel_member, dict):
            raise ValueError(f"Expected mapping entries in {protocol_path} panel_members")
        panel_member_id = str(panel_member.get("panel_member_id", "")).strip()
        topology_class = str(panel_member.get("topology_class", "")).strip()
        panel_member_type = str(panel_member.get("panel_member_type", "")).strip()
        panel_role = str(panel_member.get("panel_role", "")).strip()
        systems = panel_member.get("systems")
        if not isinstance(systems, list):
            raise ValueError(f"Expected systems list for {panel_member_id} in {protocol_path}")
        for system in systems:
            if not isinstance(system, dict):
                raise ValueError(f"Expected system mappings for {panel_member_id} in {protocol_path}")
            target_metal = str(system.get("target_metal", "")).strip()
            protocol_rows[(panel_member_id, target_metal)] = {
                "starting_structure_path": str(system.get("starting_structure_path", "")).strip(),
                "system_metadata_path": str(system.get("system_metadata_path", "")).strip(),
                "topology_class": topology_class,
                "panel_member_type": panel_member_type,
                "panel_role": panel_role,
            }
    if len(protocol_rows) != len(manifest_rows):
        raise ValueError(
            f"Protocol/manifest system count mismatch: {len(protocol_rows)} vs {len(manifest_rows)}"
        )
    for manifest_row in manifest_rows:
        key = (manifest_row.panel_member_id, manifest_row.target_metal)
        protocol_row = protocol_rows.get(key)
        if protocol_row is None:
            raise ValueError(f"Missing protocol system for {manifest_row.panel_member_id} {manifest_row.target_metal}")
        for field_name, expected_value in (
            ("starting_structure_path", manifest_row.starting_structure_path),
            ("system_metadata_path", manifest_row.system_metadata_path),
            ("topology_class", manifest_row.topology_class),
            ("panel_member_type", manifest_row.panel_member_type),
            ("panel_role", manifest_row.panel_role),
        ):
            if protocol_row[field_name] != expected_value:
                raise ValueError(
                    f"Protocol mismatch for {manifest_row.panel_member_id} {manifest_row.target_metal}: "
                    f"{field_name} {protocol_row[field_name]!r} vs {expected_value!r}"
                )


def build_default_metal_model_registry() -> tuple[MetalModelRecord, ...]:
    """Create the deterministic metal-model registry for the MD validation panel."""
    return (
        MetalModelRecord(
            metal_identity="Dy",
            formal_charge=3,
            intended_model_family=CUSTOM_MODEL_FAMILY,
            standard_forcefield_supported=False,
            custom_required=True,
            parameter_source_status=CUSTOM_PARAMETER_SOURCE_STATUS,
            notes=(
                "Dy remains a pre-bound LanM-site species in the Phase 6A2 build scaffold. "
                "The repo does not yet contain numeric OpenMM/AMBER bound-site parameters."
            ),
        ),
        MetalModelRecord(
            metal_identity="Nd",
            formal_charge=3,
            intended_model_family=CUSTOM_MODEL_FAMILY,
            standard_forcefield_supported=False,
            custom_required=True,
            parameter_source_status=CUSTOM_PARAMETER_SOURCE_STATUS,
            notes=(
                "Nd is handled as the same bound-site metal-model family used for Dy screening. "
                "No numeric bound-site parameters are present in the repo yet."
            ),
        ),
        MetalModelRecord(
            metal_identity="Y",
            formal_charge=3,
            intended_model_family=CUSTOM_MODEL_FAMILY,
            standard_forcefield_supported=False,
            custom_required=True,
            parameter_source_status=CUSTOM_PARAMETER_SOURCE_STATUS,
            notes=(
                "Y is tracked as a bound-site competitor rather than an ordinary solvent ion. "
                "Phase 6A2 records only the registry entry because custom numeric parameters are absent."
            ),
        ),
        MetalModelRecord(
            metal_identity="Al",
            formal_charge=3,
            intended_model_family=CUSTOM_MODEL_FAMILY,
            standard_forcefield_supported=False,
            custom_required=True,
            parameter_source_status=CUSTOM_PARAMETER_SOURCE_STATUS,
            notes=(
                "Al is treated as an off-target bound-site competitor in the LanM pocket. "
                "No numeric bound-site parameters are stored in the repo at this phase."
            ),
        ),
        MetalModelRecord(
            metal_identity="Fe",
            formal_charge=3,
            intended_model_family=CUSTOM_MODEL_FAMILY,
            standard_forcefield_supported=False,
            custom_required=True,
            parameter_source_status=CUSTOM_PARAMETER_SOURCE_STATUS,
            notes=(
                "Fe is carried as an Fe(III)-style trivalent competitor assumption for the validation panel. "
                "Phase 6A2 does not invent any numeric bound-site parameters."
            ),
        ),
    )


def validate_metal_model_registry(registry: tuple[MetalModelRecord, ...]) -> tuple[MetalModelRecord, ...]:
    """Validate the registry structure before it is written to disk."""
    if not registry:
        raise ValueError("Metal-model registry is empty")
    observed_metals = tuple(item.metal_identity for item in registry)
    if observed_metals != TARGET_METALS:
        raise ValueError(f"Expected registry metals {TARGET_METALS}, found {observed_metals}")
    for record in registry:
        if record.standard_forcefield_supported == record.custom_required:
            raise ValueError(
                f"Registry entry for {record.metal_identity} must select exactly one support mode"
            )
        if record.formal_charge == 0:
            raise ValueError(f"Registry entry for {record.metal_identity} must include a non-zero formal charge")
        if not record.intended_model_family:
            raise ValueError(f"Registry entry for {record.metal_identity} is missing intended_model_family")
        if not record.parameter_source_status:
            raise ValueError(f"Registry entry for {record.metal_identity} is missing parameter_source_status")
        if not record.notes:
            raise ValueError(f"Registry entry for {record.metal_identity} is missing notes")
    return registry


def render_metal_model_registry_payload(
    *,
    registry: tuple[MetalModelRecord, ...],
    protocol_path: Path,
    manifest_path: Path,
) -> dict[str, Any]:
    """Render the metal-model registry YAML payload."""
    validate_metal_model_registry(registry)
    return {
        "version": 1,
        "phase": "6A2",
        "description": "Deterministic metal-model registry for OpenMM system-build scaffolding only.",
        "source_artifacts": {
            "md_protocol": _display_path(protocol_path),
            "md_input_manifest": _display_path(manifest_path),
        },
        "metals": [
            {
                "metal_identity": record.metal_identity,
                "formal_charge": record.formal_charge,
                "intended_model_family": record.intended_model_family,
                "standard_forcefield_supported": record.standard_forcefield_supported,
                "custom_required": record.custom_required,
                "parameter_source_status": record.parameter_source_status,
                "notes": record.notes,
            }
            for record in registry
        ],
    }


def _structure_family_key(row: MDInputSystemRow) -> str:
    for value in (row.backbone_id, row.metal_site_template_source, row.source_structure):
        lowered = value.lower()
        for token in ("8fns", "8fnr", "6mi5", "8dq2"):
            if token in lowered:
                return token
    return row.backbone_id.lower()


def plan_template_water_restoration(
    *,
    system_row: MDInputSystemRow,
    manifest_rows: tuple[MDInputSystemRow, ...],
) -> TemplateWaterRestorationPlan:
    """Determine whether template metal-proximal waters are available for restoration."""
    if system_row.preserved_solvent_residue_count > 0:
        return TemplateWaterRestorationPlan(
            template_metal_proximal_waters_available=True,
            restore_template_metal_proximal_waters_before_solvation=False,
            comparison_panel_member_id=system_row.panel_member_id,
            comparison_source=system_row.starting_structure_path,
            comparison_preserved_solvent_residue_count=system_row.preserved_solvent_residue_count,
        )
    if system_row.panel_member_type != "designed_candidate":
        return TemplateWaterRestorationPlan(
            template_metal_proximal_waters_available=False,
            restore_template_metal_proximal_waters_before_solvation=False,
            comparison_panel_member_id=None,
            comparison_source=None,
            comparison_preserved_solvent_residue_count=0,
        )
    family_key = _structure_family_key(system_row)
    candidates = [
        row
        for row in manifest_rows
        if row.target_metal == system_row.target_metal
        and row.preserved_solvent_residue_count > 0
        and _structure_family_key(row) == family_key
    ]
    candidates.sort(
        key=lambda row: (
            row.panel_member_type != "wild_type_reference",
            -row.preserved_solvent_residue_count,
            row.panel_rank,
            row.panel_member_id,
        )
    )
    if not candidates:
        return TemplateWaterRestorationPlan(
            template_metal_proximal_waters_available=False,
            restore_template_metal_proximal_waters_before_solvation=False,
            comparison_panel_member_id=None,
            comparison_source=None,
            comparison_preserved_solvent_residue_count=0,
        )
    reference_row = candidates[0]
    reference_source = (
        reference_row.metal_site_template_source
        if reference_row.metal_site_template_source != reference_row.source_structure
        else reference_row.starting_structure_path
    )
    return TemplateWaterRestorationPlan(
        template_metal_proximal_waters_available=True,
        restore_template_metal_proximal_waters_before_solvation=True,
        comparison_panel_member_id=reference_row.panel_member_id,
        comparison_source=reference_source,
        comparison_preserved_solvent_residue_count=reference_row.preserved_solvent_residue_count,
    )


def _registry_by_metal(registry: tuple[MetalModelRecord, ...]) -> dict[str, MetalModelRecord]:
    return {record.metal_identity: record for record in registry}


def build_system_build_request_payload(
    *,
    system_row: MDInputSystemRow,
    build_config: OpenMMBuildScaffoldingConfig,
    metal_model: MetalModelRecord,
    water_restoration_plan: TemplateWaterRestorationPlan,
    metal_model_registry_path: Path,
    manifest_path: Path,
    protocol_path: Path,
) -> dict[str, Any]:
    """Build the deterministic per-system YAML request payload."""
    return {
        "version": 1,
        "phase": "6A2",
        "panel_member_id": system_row.panel_member_id,
        "panel_member_type": system_row.panel_member_type,
        "panel_role": system_row.panel_role,
        "candidate_id": system_row.candidate_id,
        "source_starting_structure_pdb": system_row.starting_structure_path,
        "source_system_metadata": system_row.system_metadata_path,
        "source_artifacts": {
            "md_input_manifest": _display_path(manifest_path),
            "md_protocol": _display_path(protocol_path),
            "metal_model_registry": _display_path(metal_model_registry_path),
        },
        "topology_class": system_row.topology_class,
        "target_metal": system_row.target_metal,
        "chosen_protein_forcefield": build_config.protein_forcefield,
        "chosen_water_model": build_config.water_model,
        "box_padding_target_nm": build_config.box_padding_nm,
        "ionic_strength_target_M": build_config.ionic_strength_M,
        "neutralization_ion_policy": build_config.neutralization_ion_policy,
        "custom_metal_parameters_still_required": metal_model.custom_required,
        "generic_protein_water_preparation_ready": True,
        "restore_template_metal_proximal_waters_before_solvation": (
            water_restoration_plan.restore_template_metal_proximal_waters_before_solvation
        ),
        "template_metal_proximal_waters_available": water_restoration_plan.template_metal_proximal_waters_available,
        "template_water_restoration_reference": {
            "comparison_panel_member_id": water_restoration_plan.comparison_panel_member_id,
            "comparison_source": water_restoration_plan.comparison_source,
            "comparison_preserved_solvent_residue_count": (
                water_restoration_plan.comparison_preserved_solvent_residue_count
            ),
        },
        "metal_model": {
            "metal_identity": metal_model.metal_identity,
            "formal_charge": metal_model.formal_charge,
            "intended_model_family": metal_model.intended_model_family,
            "standard_forcefield_supported": metal_model.standard_forcefield_supported,
            "custom_required": metal_model.custom_required,
            "parameter_source_status": metal_model.parameter_source_status,
            "notes": metal_model.notes,
        },
        "notes": [
            "This Phase 6A2 artifact plans deterministic OpenMM system-building inputs only.",
            "No OpenMM System object, minimization, MD, metadynamics, QM, or quantum step is executed here.",
        ],
    }


def render_md_system_building_plan_markdown(
    *,
    manifest_rows: tuple[MDSystemBuildManifestRow, ...],
    registry: tuple[MetalModelRecord, ...],
    build_config: OpenMMBuildScaffoldingConfig,
) -> str:
    """Render the Phase 6A2 MD system-building report."""
    direct_as_is_rows = [row for row in manifest_rows if not row.restore_template_metal_proximal_waters_before_solvation]
    restore_plan_rows = [row for row in manifest_rows if row.restore_template_metal_proximal_waters_before_solvation]
    blocked_rows = [row for row in manifest_rows if row.custom_metal_parameters_still_required]
    lines = [
        "# MD System Building Plan",
        "",
        "Phase 6A2 creates deterministic OpenMM build requests and a metal-model registry for the selected MD validation panel.",
        "",
        "## Build Defaults",
        "",
        f"- protein force field scaffold: `{build_config.protein_forcefield}`",
        f"- water model scaffold: `{build_config.water_model}`",
        f"- box padding target: `{build_config.box_padding_nm:.2f} nm`",
        f"- ionic strength target: `{build_config.ionic_strength_M:.2f} M`",
        f"- neutralization ion policy: `{build_config.neutralization_ion_policy}`",
        "",
        "## Why These Metals Stay Custom",
        "",
        "- Dy, Nd, Y, Al, and Fe are treated as custom bound species because they are carried in the starting structures as pre-bound LanM-site cofactors, not as freely exchanged bulk-solvent ions.",
        "- The project brief explicitly calls for a custom 12-6-4 / chelator-tuned treatment of the metal center, and the repo does not yet contain numeric bound-site parameters for any panel metal.",
        "- Phase 6A2 therefore keeps generic protein and water setup separate from the missing metal-site parameterization step needed for full OpenMM `System` creation.",
        "",
        "## Registry Summary",
        "",
        "| metal | formal_charge | intended_model_family | standard_forcefield_supported | custom_required | parameter_source_status |",
        "| --- | ---: | --- | --- | --- | --- |",
    ]
    for record in registry:
        lines.append(
            "| "
            f"{record.metal_identity} | "
            f"{record.formal_charge:+d} | "
            f"{record.intended_model_family} | "
            f"{record.standard_forcefield_supported} | "
            f"{record.custom_required} | "
            f"{record.parameter_source_status} |"
        )
    lines.extend(
        [
            "",
            "## System Readiness",
            "",
            f"- systems planned: `{len(manifest_rows)}`",
            f"- systems that can proceed directly with generic protein/water preparation as written: `{len(direct_as_is_rows)}`",
            f"- systems with optional template-water restoration flagged before solvation: `{len(restore_plan_rows)}`",
            f"- systems still blocked from full OpenMM `System` creation pending custom metal parameters: `{len(blocked_rows)}`",
            "",
            "| panel_member_id | target_metal | topology_class | generic_protein_water_preparation_ready | restore_template_waters | custom_parameters_still_required | template_water_reference |",
            "| --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in manifest_rows:
        reference = row.template_water_reference_panel_member_id or "-"
        lines.append(
            "| "
            f"{row.panel_member_id} | "
            f"{row.target_metal} | "
            f"{row.topology_class} | "
            f"{row.generic_protein_water_preparation_ready} | "
            f"{row.restore_template_metal_proximal_waters_before_solvation} | "
            f"{row.custom_metal_parameters_still_required} | "
            f"{reference} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- Wild-type reference systems already preserve their metal-proximal waters and can move straight into generic protein/water setup planning.",
            "- Designed systems with zero preserved solvent are compared against the matching wild-type/template family so optional metal-proximal waters can be restored deterministically before bulk solvation.",
            "- All systems still require custom metal parameters before an actual OpenMM `System` object can be created.",
            "",
        ]
    )
    return "\n".join(lines)


def prepare_md_system_build(
    *,
    md_input_manifest_path: Path,
    md_protocol_path: Path,
    metal_model_registry_path: Path,
    md_system_build_root: Path,
    md_system_build_manifest_path: Path,
    md_system_build_report_path: Path,
) -> tuple[MDSystemBuildManifestRow, ...]:
    """Generate deterministic Phase 6A2 OpenMM system-build scaffolding."""
    manifest_rows = discover_md_input_manifest(md_input_manifest_path)
    validate_md_input_manifest_against_protocol(
        manifest_rows=manifest_rows,
        protocol_path=md_protocol_path,
    )
    build_config = load_md_protocol_build_config(md_protocol_path)
    registry = validate_metal_model_registry(build_default_metal_model_registry())
    write_yaml(
        metal_model_registry_path,
        render_metal_model_registry_payload(
            registry=registry,
            protocol_path=md_protocol_path,
            manifest_path=md_input_manifest_path,
        ),
    )
    registry_by_metal = _registry_by_metal(registry)
    output_rows: list[MDSystemBuildManifestRow] = []
    for system_row in manifest_rows:
        water_restoration_plan = plan_template_water_restoration(
            system_row=system_row,
            manifest_rows=manifest_rows,
        )
        metal_model = registry_by_metal[system_row.target_metal]
        build_request_path = md_system_build_root / system_row.panel_member_id / system_row.target_metal / "build_request.yaml"
        write_yaml(
            build_request_path,
            build_system_build_request_payload(
                system_row=system_row,
                build_config=build_config,
                metal_model=metal_model,
                water_restoration_plan=water_restoration_plan,
                metal_model_registry_path=metal_model_registry_path,
                manifest_path=md_input_manifest_path,
                protocol_path=md_protocol_path,
            ),
        )
        output_rows.append(
            MDSystemBuildManifestRow(
                panel_rank=system_row.panel_rank,
                panel_member_id=system_row.panel_member_id,
                panel_member_type=system_row.panel_member_type,
                panel_role=system_row.panel_role,
                topology_class=system_row.topology_class,
                target_metal=system_row.target_metal,
                source_starting_structure_pdb=system_row.starting_structure_path,
                build_request_path=_display_path(build_request_path),
                chosen_protein_forcefield=build_config.protein_forcefield,
                chosen_water_model=build_config.water_model,
                box_padding_target_nm=build_config.box_padding_nm,
                ionic_strength_target_M=build_config.ionic_strength_M,
                neutralization_ion_policy=build_config.neutralization_ion_policy,
                custom_metal_parameters_still_required=metal_model.custom_required,
                generic_protein_water_preparation_ready=True,
                template_metal_proximal_waters_available=(
                    water_restoration_plan.template_metal_proximal_waters_available
                ),
                restore_template_metal_proximal_waters_before_solvation=(
                    water_restoration_plan.restore_template_metal_proximal_waters_before_solvation
                ),
                template_water_reference_panel_member_id=water_restoration_plan.comparison_panel_member_id or "",
                template_water_reference_source=water_restoration_plan.comparison_source or "",
            )
        )
    manifest_tuple = tuple(output_rows)
    write_csv_rows(md_system_build_manifest_path, manifest_tuple)
    atomic_write_text(
        md_system_build_report_path,
        render_md_system_building_plan_markdown(
            manifest_rows=manifest_tuple,
            registry=registry,
            build_config=build_config,
        ),
    )
    return manifest_tuple
