"""Phase 3A deterministic design-campaign export and ProteinMPNN preparation."""

from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

try:
    import yaml  # type: ignore
except ModuleNotFoundError:  # pragma: no cover - optional dependency in local envs
    yaml = None

from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_text, write_csv_rows, write_jsonl, write_yaml
from lanm.models import (
    CrossTemplateAlignmentRow,
    DesignBackboneManifestRow,
    DesignCampaignManifestRow,
    DesignCampaignPositionRow,
    DesignMaskCandidate,
)
from lanm.paths import REPO_ROOT
from lanm.structure.atoms import read_atom_records
from lanm.structure.pdb import render_protein_pdb
from lanm.structure.residues import collect_polymer_residues
from lanm.structure.templates import parse_cif_atom_records

_STRUCTURE_SOURCES: dict[str, tuple[str, Path]] = {
    "8FNS": ("csv_atom_table", Path("data/raw/local_bundle/8fns_atoms.csv")),
    "8FNR": ("cif", Path("data/raw/public/structures/8FNR.cif")),
}
_DESIGN_SET_DESCRIPTIONS = {
    "hard_fixed": "fixed_first_shell union protected_positions",
    "campaign_ss_only": "mutable_second_sphere minus hard_fixed",
    "campaign_ss_plus_if": "union(mutable_second_sphere, mutable_interface) minus hard_fixed",
    "campaign_if_only": "mutable_interface minus hard_fixed",
}


@dataclass(frozen=True, slots=True)
class TemplateChainRecord:
    template_id: str
    source_type: str
    chain_id: str
    residue_count: int
    metal_site_count: int
    experimental_method: str


@dataclass(frozen=True, slots=True)
class DesignSetDefinition:
    name: str
    canonical_positions: tuple[int, ...]
    am1_positions: tuple[int, ...]
    description: str


@dataclass(frozen=True, slots=True)
class BackboneDefinition:
    backbone_id: str
    label: str
    structure_id: str
    source_kind: str
    source_path: Path
    experimental_method: str
    selected_chains: tuple[str, ...]
    design_chains: tuple[str, ...]
    fixed_context_chains: tuple[str, ...]
    selected_residue_count: int
    design_chain_residue_count: int
    design_chain_metal_site_count: int
    representative_note: str


@dataclass(frozen=True, slots=True)
class BackboneResiduePosition:
    chain_id: str
    sequence_index: int
    residue_seq: int
    insertion_code: str
    residue_name: str
    canonical_family_position: int | None
    am1_mature_position: int | None


@dataclass(frozen=True, slots=True)
class CampaignDefinition:
    campaign_id: str
    label: str
    backbone_id: str
    design_set_name: str
    designed_chains: tuple[str, ...]
    fixed_context_chains: tuple[str, ...]
    designable_canonical_positions: tuple[int, ...]
    designable_am1_positions: tuple[int, ...]
    designable_residue_count: int
    fixed_residue_count: int
    chain_assignment_payload: dict[str, object]
    fixed_positions_payload: dict[str, object]


@dataclass(frozen=True, slots=True)
class DesignCampaignArtifacts:
    design_sets: tuple[DesignSetDefinition, ...]
    backbones: tuple[BackboneDefinition, ...]
    campaigns: tuple[CampaignDefinition, ...]
    campaign_position_rows: tuple[DesignCampaignPositionRow, ...]
    backbone_pdb_text_by_id: dict[str, str]


def _display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def _positions_text(positions: tuple[int, ...]) -> str:
    return ",".join(str(position) for position in positions)


def _parse_optional_int(value: str) -> int | None:
    stripped = value.strip()
    return int(stripped) if stripped else None


def _parse_bool(value: str) -> bool:
    return value.strip().lower() == "true"


def _load_yaml_mapping(path: Path) -> dict[str, Any]:
    require_path(path)
    if yaml is None:
        raise RuntimeError("PyYAML is required to read nested design mask configuration")
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Unexpected YAML payload in {path}")
    return payload


def load_design_mask_candidates(path: Path) -> tuple[DesignMaskCandidate, ...]:
    """Load Phase 2C design-mask candidates from CSV."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            DesignMaskCandidate(
                canonical_family_position=int(str(row["canonical_family_position"]).strip()),
                am1_mature_position=_parse_optional_int(str(row["am1_mature_position"])),
                am1_reference_residue=str(row["am1_reference_residue"]).strip(),
                observed_residue_identities=str(row["observed_residue_identities"]).strip(),
                template_coverage_count=int(str(row["template_coverage_count"]).strip()),
                first_shell_observation_count=int(str(row["first_shell_observation_count"]).strip()),
                second_sphere_observation_count=int(str(row["second_sphere_observation_count"]).strip()),
                hans_interface_observation_count=int(str(row["hans_interface_observation_count"]).strip()),
                interface_neighborhood=_parse_bool(str(row["interface_neighborhood"])),
                fixed_first_shell=_parse_bool(str(row["fixed_first_shell"])),
                mutable_second_sphere=_parse_bool(str(row["mutable_second_sphere"])),
                mutable_interface=_parse_bool(str(row["mutable_interface"])),
                protected_positions=_parse_bool(str(row["protected_positions"])),
                protection_reasons=str(row["protection_reasons"]).strip(),
                rationale=str(row["rationale"]).strip(),
            )
            for row in reader
        ]
    return tuple(sorted(rows, key=lambda row: row.canonical_family_position))


def load_cross_template_alignment_rows(path: Path) -> tuple[CrossTemplateAlignmentRow, ...]:
    """Load the cross-template residue alignment table from CSV."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            CrossTemplateAlignmentRow(
                template_id=str(row["template_id"]).strip(),
                chain_id=str(row["chain_id"]).strip(),
                template_residue_seq=int(str(row["template_residue_seq"]).strip()),
                template_residue_name=str(row["template_residue_name"]).strip(),
                canonical_family_position=_parse_optional_int(str(row["canonical_family_position"])),
                am1_mature_position=_parse_optional_int(str(row["am1_mature_position"])),
                alignment_status=str(row["alignment_status"]).strip(),
            )
            for row in reader
        ]
    return tuple(rows)


def load_template_chain_records(path: Path) -> tuple[TemplateChainRecord, ...]:
    """Load the Phase 2 template-chain summary table."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            TemplateChainRecord(
                template_id=str(row["template_id"]).strip(),
                source_type=str(row["source_type"]).strip(),
                chain_id=str(row["chain_id"]).strip(),
                residue_count=int(str(row["residue_count"]).strip()),
                metal_site_count=int(str(row["metal_site_count"]).strip()),
                experimental_method=str(row["experimental_method"]).strip(),
            )
            for row in reader
        ]
    return tuple(
        sorted(
            rows,
            key=lambda row: (row.template_id, row.chain_id),
        )
    )


def _mask_positions_from_rows(
    candidate_rows: tuple[DesignMaskCandidate, ...],
    attribute_name: str,
) -> tuple[int, ...]:
    return tuple(
        row.canonical_family_position
        for row in candidate_rows
        if getattr(row, attribute_name)
    )


def _mask_am1_positions_from_rows(
    candidate_rows: tuple[DesignMaskCandidate, ...],
    attribute_name: str,
) -> tuple[int, ...]:
    return tuple(
        sorted(
            {
                row.am1_mature_position
                for row in candidate_rows
                if getattr(row, attribute_name) and row.am1_mature_position is not None
            }
        )
    )


def _validate_mask_config(
    mask_config: dict[str, Any],
    candidate_rows: tuple[DesignMaskCandidate, ...],
) -> None:
    expected_masks = {
        "fixed_first_shell": "fixed_first_shell",
        "mutable_second_sphere": "mutable_second_sphere",
        "mutable_interface": "mutable_interface",
        "protected_positions": "protected_positions",
    }
    payload_masks = mask_config.get("masks")
    if not isinstance(payload_masks, dict):
        raise ValueError("design_masks.yaml is missing a top-level masks mapping")
    for mask_name, attribute_name in expected_masks.items():
        mask_payload = payload_masks.get(mask_name)
        if not isinstance(mask_payload, dict):
            raise ValueError(f"design_masks.yaml is missing masks.{mask_name}")
        yaml_positions = tuple(
            int(position)
            for position in mask_payload.get("canonical_family_positions", [])
        )
        csv_positions = _mask_positions_from_rows(candidate_rows, attribute_name)
        if yaml_positions != csv_positions:
            raise ValueError(
                f"Mask mismatch for {mask_name}: YAML has {yaml_positions}, CSV has {csv_positions}"
            )
        yaml_am1_positions = tuple(
            int(position)
            for position in mask_payload.get("am1_mature_positions", [])
        )
        csv_am1_positions = _mask_am1_positions_from_rows(candidate_rows, attribute_name)
        if yaml_am1_positions != csv_am1_positions:
            raise ValueError(
                f"AM1 mask mismatch for {mask_name}: YAML has {yaml_am1_positions}, CSV has {csv_am1_positions}"
            )


def _positions_to_am1(
    positions: set[int],
    candidate_by_position: dict[int, DesignMaskCandidate],
) -> tuple[int, ...]:
    return tuple(
        sorted(
            {
                candidate_by_position[position].am1_mature_position
                for position in positions
                if position in candidate_by_position
                and candidate_by_position[position].am1_mature_position is not None
            }
        )
    )


def _build_design_sets(
    *,
    mask_config: dict[str, Any],
    candidate_by_position: dict[int, DesignMaskCandidate],
) -> tuple[DesignSetDefinition, ...]:
    masks = mask_config["masks"]
    fixed_first_shell = {
        int(position)
        for position in masks["fixed_first_shell"]["canonical_family_positions"]
    }
    protected_positions = {
        int(position)
        for position in masks["protected_positions"]["canonical_family_positions"]
    }
    mutable_second_sphere = {
        int(position)
        for position in masks["mutable_second_sphere"]["canonical_family_positions"]
    }
    mutable_interface = {
        int(position)
        for position in masks["mutable_interface"]["canonical_family_positions"]
    }
    hard_fixed = fixed_first_shell | protected_positions
    derived_sets = {
        "hard_fixed": tuple(sorted(hard_fixed)),
        "campaign_ss_only": tuple(sorted(mutable_second_sphere - hard_fixed)),
        "campaign_ss_plus_if": tuple(sorted((mutable_second_sphere | mutable_interface) - hard_fixed)),
        "campaign_if_only": tuple(sorted(mutable_interface - hard_fixed)),
    }
    return tuple(
        DesignSetDefinition(
            name=name,
            canonical_positions=positions,
            am1_positions=_positions_to_am1(set(positions), candidate_by_position),
            description=_DESIGN_SET_DESCRIPTIONS[name],
        )
        for name, positions in derived_sets.items()
    )


def _load_structure_atoms(structure_id: str) -> tuple:
    source_kind, relative_path = _STRUCTURE_SOURCES[structure_id]
    source_path = REPO_ROOT / relative_path
    require_path(source_path)
    if source_kind == "cif":
        return tuple(parse_cif_atom_records(source_path, structure_id=structure_id))
    return tuple(read_atom_records(source_path, structure_id=structure_id))


def _select_hans_representative_chain(
    chain_records: tuple[TemplateChainRecord, ...],
) -> TemplateChainRecord:
    hans_rows = [row for row in chain_records if row.template_id == "8FNR"]
    if not hans_rows:
        raise ValueError("template_chain_summary.csv does not contain 8FNR rows")
    return sorted(
        hans_rows,
        key=lambda row: (-row.metal_site_count, -row.residue_count, row.chain_id),
    )[0]


def _build_backbones(
    chain_records: tuple[TemplateChainRecord, ...],
) -> tuple[BackboneDefinition, ...]:
    chain_map: dict[tuple[str, str], TemplateChainRecord] = {
        (row.template_id, row.chain_id): row
        for row in chain_records
    }
    am1_row = chain_map.get(("8FNS", "A"))
    if am1_row is None:
        raise ValueError("template_chain_summary.csv does not contain 8FNS chain A")
    hans_row = _select_hans_representative_chain(chain_records)
    hans_chain_ids = tuple(
        sorted(row.chain_id for row in chain_records if row.template_id == "8FNR")
    )
    hans_selected_residue_count = sum(
        chain_map[("8FNR", chain_id)].residue_count
        for chain_id in hans_chain_ids
    )
    return (
        BackboneDefinition(
            backbone_id="am1_mex_8fns_chain_a",
            label="AM1/Mex default",
            structure_id="8FNS",
            source_kind=_STRUCTURE_SOURCES["8FNS"][0],
            source_path=_STRUCTURE_SOURCES["8FNS"][1],
            experimental_method=am1_row.experimental_method,
            selected_chains=("A",),
            design_chains=("A",),
            fixed_context_chains=(),
            selected_residue_count=am1_row.residue_count,
            design_chain_residue_count=am1_row.residue_count,
            design_chain_metal_site_count=am1_row.metal_site_count,
            representative_note="Phase 3A default backbone specified by the workflow: 8FNS chain A.",
        ),
        BackboneDefinition(
            backbone_id=f"hans_pocket_8fnr_chain_{hans_row.chain_id.lower()}",
            label="Hans pocket-focused",
            structure_id="8FNR",
            source_kind=_STRUCTURE_SOURCES["8FNR"][0],
            source_path=_STRUCTURE_SOURCES["8FNR"][1],
            experimental_method=hans_row.experimental_method,
            selected_chains=(hans_row.chain_id,),
            design_chains=(hans_row.chain_id,),
            fixed_context_chains=(),
            selected_residue_count=hans_row.residue_count,
            design_chain_residue_count=hans_row.residue_count,
            design_chain_metal_site_count=hans_row.metal_site_count,
            representative_note=(
                "Deterministic 8FNR representative chain selected by descending metal_site_count, "
                f"descending residue_count, then ascending chain_id; this chooses chain {hans_row.chain_id}."
            ),
        ),
        BackboneDefinition(
            backbone_id=f"hans_interface_8fnr_{'_'.join(chain_id.lower() for chain_id in hans_chain_ids)}",
            label="Hans interface-aware",
            structure_id="8FNR",
            source_kind=_STRUCTURE_SOURCES["8FNR"][0],
            source_path=_STRUCTURE_SOURCES["8FNR"][1],
            experimental_method=hans_row.experimental_method,
            selected_chains=hans_chain_ids,
            design_chains=(hans_row.chain_id,),
            fixed_context_chains=tuple(
                chain_id
                for chain_id in hans_chain_ids
                if chain_id != hans_row.chain_id
            ),
            selected_residue_count=hans_selected_residue_count,
            design_chain_residue_count=hans_row.residue_count,
            design_chain_metal_site_count=hans_row.metal_site_count,
            representative_note=(
                f"Full 8FNR assembly exported with representative design chain {hans_row.chain_id} "
                "and the remaining chains retained as fixed interface context."
            ),
        ),
    )


def _build_backbone_residue_positions(
    *,
    structure_id: str,
    atoms: tuple,
    design_chains: tuple[str, ...],
    alignment_rows: tuple[CrossTemplateAlignmentRow, ...],
) -> dict[str, tuple[BackboneResiduePosition, ...]]:
    alignment_lookup = {
        (row.template_id, row.chain_id, row.template_residue_seq): row
        for row in alignment_rows
    }
    selected_atoms = [
        atom
        for atom in atoms
        if atom.record_type == "ATOM" and atom.chain_id in design_chains
    ]
    chain_residues = collect_polymer_residues(selected_atoms)
    residue_positions: dict[str, tuple[BackboneResiduePosition, ...]] = {}
    for chain_id in design_chains:
        residues = chain_residues.get(chain_id)
        if residues is None:
            raise ValueError(f"{structure_id} is missing design chain {chain_id}")
        chain_rows: list[BackboneResiduePosition] = []
        for index, residue in enumerate(residues, start=1):
            alignment = alignment_lookup.get((structure_id, residue.chain_id, residue.residue_seq))
            chain_rows.append(
                BackboneResiduePosition(
                    chain_id=residue.chain_id,
                    sequence_index=index,
                    residue_seq=residue.residue_seq,
                    insertion_code=residue.insertion_code,
                    residue_name=residue.residue_name,
                    canonical_family_position=alignment.canonical_family_position if alignment else None,
                    am1_mature_position=alignment.am1_mature_position if alignment else None,
                )
            )
        residue_positions[chain_id] = tuple(chain_rows)
    return residue_positions


def _build_campaign_position_rows(
    *,
    campaign_id: str,
    design_set_name: str,
    backbone_id: str,
    structure_id: str,
    residue_positions: dict[str, tuple[BackboneResiduePosition, ...]],
    designable_canonical_positions: set[int],
    hard_fixed_positions: set[int],
    candidate_by_position: dict[int, DesignMaskCandidate],
) -> tuple[DesignCampaignPositionRow, ...]:
    rows: list[DesignCampaignPositionRow] = []
    for chain_id in sorted(residue_positions):
        for residue in residue_positions[chain_id]:
            candidate = (
                candidate_by_position[residue.canonical_family_position]
                if residue.canonical_family_position in candidate_by_position
                else None
            )
            designable = residue.canonical_family_position in designable_canonical_positions
            hard_fixed = residue.canonical_family_position in hard_fixed_positions
            if designable:
                position_state = "designable"
                position_reason = "campaign_designable"
            elif hard_fixed:
                position_state = "fixed"
                position_reason = "hard_fixed"
            elif residue.canonical_family_position is None:
                position_state = "fixed"
                position_reason = "unmapped_backbone_position"
            else:
                position_state = "fixed"
                position_reason = "not_in_campaign_design_set"
            rows.append(
                DesignCampaignPositionRow(
                    campaign_id=campaign_id,
                    design_set_name=design_set_name,
                    backbone_id=backbone_id,
                    structure_id=structure_id,
                    chain_id=residue.chain_id,
                    sequence_index=residue.sequence_index,
                    residue_seq=residue.residue_seq,
                    insertion_code=residue.insertion_code,
                    residue_name=residue.residue_name,
                    canonical_family_position=residue.canonical_family_position,
                    am1_mature_position=residue.am1_mature_position,
                    fixed_first_shell=candidate.fixed_first_shell if candidate else False,
                    protected_positions=candidate.protected_positions if candidate else False,
                    mutable_second_sphere=candidate.mutable_second_sphere if candidate else False,
                    mutable_interface=candidate.mutable_interface if candidate else False,
                    hard_fixed=hard_fixed,
                    designable=designable,
                    position_state=position_state,
                    position_reason=position_reason,
                    rationale=candidate.rationale if candidate else "",
                )
            )
    return tuple(rows)


def _build_campaigns(
    *,
    design_sets: tuple[DesignSetDefinition, ...],
    backbones: tuple[BackboneDefinition, ...],
    alignment_rows: tuple[CrossTemplateAlignmentRow, ...],
    candidate_by_position: dict[int, DesignMaskCandidate],
) -> tuple[tuple[CampaignDefinition, ...], tuple[DesignCampaignPositionRow, ...]]:
    design_set_map = {design_set.name: design_set for design_set in design_sets}
    hard_fixed_positions = set(design_set_map["hard_fixed"].canonical_positions)
    backbone_map = {backbone.backbone_id: backbone for backbone in backbones}
    structure_atoms = {
        backbone.structure_id: _load_structure_atoms(backbone.structure_id)
        for backbone in backbones
    }
    residue_positions_by_backbone = {
        backbone.backbone_id: _build_backbone_residue_positions(
            structure_id=backbone.structure_id,
            atoms=structure_atoms[backbone.structure_id],
            design_chains=backbone.design_chains,
            alignment_rows=alignment_rows,
        )
        for backbone in backbones
    }
    campaign_specs = (
        (
            "am1_mex_ss_only",
            "AM1/Mex second-sphere only",
            "am1_mex_8fns_chain_a",
            "campaign_ss_only",
        ),
        (
            "hans_pocket_ss_only",
            "Hans pocket-focused second-sphere only",
            next(
                backbone.backbone_id
                for backbone in backbones
                if backbone.label == "Hans pocket-focused"
            ),
            "campaign_ss_only",
        ),
        (
            "hans_interface_ss_plus_if",
            "Hans interface-aware second-sphere plus interface",
            next(
                backbone.backbone_id
                for backbone in backbones
                if backbone.label == "Hans interface-aware"
            ),
            "campaign_ss_plus_if",
        ),
        (
            "hans_interface_if_only",
            "Hans interface-aware interface only",
            next(
                backbone.backbone_id
                for backbone in backbones
                if backbone.label == "Hans interface-aware"
            ),
            "campaign_if_only",
        ),
    )
    campaigns: list[CampaignDefinition] = []
    all_position_rows: list[DesignCampaignPositionRow] = []
    for campaign_id, label, backbone_id, design_set_name in campaign_specs:
        backbone = backbone_map[backbone_id]
        design_set = design_set_map[design_set_name]
        residue_positions = residue_positions_by_backbone[backbone_id]
        backbone_mapped_positions = {
            residue.canonical_family_position
            for residues in residue_positions.values()
            for residue in residues
            if residue.canonical_family_position is not None
        }
        missing_positions = set(design_set.canonical_positions) - backbone_mapped_positions
        if missing_positions:
            raise ValueError(
                f"{campaign_id} is missing canonical positions on {backbone_id}: {sorted(missing_positions)}"
            )
        position_rows = _build_campaign_position_rows(
            campaign_id=campaign_id,
            design_set_name=design_set_name,
            backbone_id=backbone_id,
            structure_id=backbone.structure_id,
            residue_positions=residue_positions,
            designable_canonical_positions=set(design_set.canonical_positions),
            hard_fixed_positions=hard_fixed_positions,
            candidate_by_position=candidate_by_position,
        )
        fixed_positions_payload = {
            campaign_id: {
                chain_id: [
                    row.sequence_index
                    for row in position_rows
                    if row.chain_id == chain_id and not row.designable
                ]
                for chain_id in backbone.design_chains
            }
        }
        designable_residue_count = sum(1 for row in position_rows if row.designable)
        fixed_residue_count = sum(1 for row in position_rows if not row.designable)
        all_position_rows.extend(position_rows)
        campaigns.append(
            CampaignDefinition(
                campaign_id=campaign_id,
                label=label,
                backbone_id=backbone_id,
                design_set_name=design_set_name,
                designed_chains=backbone.design_chains,
                fixed_context_chains=backbone.fixed_context_chains,
                designable_canonical_positions=design_set.canonical_positions,
                designable_am1_positions=design_set.am1_positions,
                designable_residue_count=designable_residue_count,
                fixed_residue_count=fixed_residue_count,
                chain_assignment_payload={
                    campaign_id: [list(backbone.design_chains), list(backbone.fixed_context_chains)]
                },
                fixed_positions_payload=fixed_positions_payload,
            )
        )
    return tuple(campaigns), tuple(all_position_rows)


def build_design_campaign_artifacts(
    *,
    design_masks_path: Path,
    design_mask_candidates_path: Path,
    cross_template_alignment_path: Path,
    template_chain_summary_path: Path,
) -> DesignCampaignArtifacts:
    """Build deterministic Phase 3A design-campaign artifacts from Phase 2 outputs."""
    mask_config = _load_yaml_mapping(design_masks_path)
    candidate_rows = load_design_mask_candidates(design_mask_candidates_path)
    _validate_mask_config(mask_config, candidate_rows)
    candidate_by_position = {
        row.canonical_family_position: row
        for row in candidate_rows
    }
    design_sets = _build_design_sets(
        mask_config=mask_config,
        candidate_by_position=candidate_by_position,
    )
    alignment_rows = load_cross_template_alignment_rows(cross_template_alignment_path)
    chain_records = load_template_chain_records(template_chain_summary_path)
    backbones = _build_backbones(chain_records)
    structure_atoms = {
        backbone.structure_id: _load_structure_atoms(backbone.structure_id)
        for backbone in backbones
    }
    backbone_pdb_text_by_id = {
        backbone.backbone_id: render_protein_pdb(
            structure_atoms[backbone.structure_id],
            selected_chain_ids=backbone.selected_chains,
        )
        for backbone in backbones
    }
    campaigns, campaign_position_rows = _build_campaigns(
        design_sets=design_sets,
        backbones=backbones,
        alignment_rows=alignment_rows,
        candidate_by_position=candidate_by_position,
    )
    return DesignCampaignArtifacts(
        design_sets=design_sets,
        backbones=backbones,
        campaigns=campaigns,
        campaign_position_rows=campaign_position_rows,
        backbone_pdb_text_by_id=backbone_pdb_text_by_id,
    )


def _render_markdown_table(headers: tuple[str, ...], rows: list[tuple[str, ...]]) -> str:
    header_line = "| " + " | ".join(headers) + " |"
    divider_line = "| " + " | ".join(["---"] * len(headers)) + " |"
    body_lines = ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join([header_line, divider_line, *body_lines])


def _build_output_path_map(
    *,
    backbones: tuple[BackboneDefinition, ...],
    campaigns: tuple[CampaignDefinition, ...],
    proteinmpnn_root: Path,
) -> tuple[dict[str, Path], dict[str, dict[str, Path]]]:
    backbone_paths = {
        backbone.backbone_id: proteinmpnn_root / "backbones" / f"{backbone.backbone_id}.pdb"
        for backbone in backbones
    }
    campaign_paths = {
        campaign.campaign_id: {
            "directory": proteinmpnn_root / "campaigns" / campaign.campaign_id,
            "pdb": proteinmpnn_root / "campaigns" / campaign.campaign_id / f"{campaign.campaign_id}.pdb",
            "chain_assignment": proteinmpnn_root / "campaigns" / campaign.campaign_id / "chain_id.jsonl",
            "fixed_positions": proteinmpnn_root / "campaigns" / campaign.campaign_id / "fixed_positions.jsonl",
        }
        for campaign in campaigns
    }
    return backbone_paths, campaign_paths


def _build_backbone_manifest_rows(
    *,
    backbones: tuple[BackboneDefinition, ...],
    backbone_paths: dict[str, Path],
) -> tuple[DesignBackboneManifestRow, ...]:
    return tuple(
        DesignBackboneManifestRow(
            backbone_id=backbone.backbone_id,
            structure_id=backbone.structure_id,
            source_kind=backbone.source_kind,
            source_path=_display_path(REPO_ROOT / backbone.source_path),
            experimental_method=backbone.experimental_method,
            selected_chains=",".join(backbone.selected_chains),
            design_chains=",".join(backbone.design_chains),
            fixed_context_chains=",".join(backbone.fixed_context_chains),
            selected_residue_count=backbone.selected_residue_count,
            design_chain_residue_count=backbone.design_chain_residue_count,
            design_chain_metal_site_count=backbone.design_chain_metal_site_count,
            representative_note=backbone.representative_note,
            backbone_pdb_path=_display_path(backbone_paths[backbone.backbone_id]),
        )
        for backbone in backbones
    )


def _build_campaign_manifest_rows(
    *,
    campaigns: tuple[CampaignDefinition, ...],
    campaign_paths: dict[str, dict[str, Path]],
) -> tuple[DesignCampaignManifestRow, ...]:
    return tuple(
        DesignCampaignManifestRow(
            campaign_id=campaign.campaign_id,
            backbone_id=campaign.backbone_id,
            design_set_name=campaign.design_set_name,
            designed_chains=",".join(campaign.designed_chains),
            fixed_context_chains=",".join(campaign.fixed_context_chains),
            exported_pdb_path=_display_path(campaign_paths[campaign.campaign_id]["pdb"]),
            chain_assignment_path=_display_path(campaign_paths[campaign.campaign_id]["chain_assignment"]),
            fixed_positions_path=_display_path(campaign_paths[campaign.campaign_id]["fixed_positions"]),
            designable_residue_count=campaign.designable_residue_count,
            fixed_residue_count=campaign.fixed_residue_count,
            designable_canonical_positions=_positions_text(campaign.designable_canonical_positions),
            designable_am1_positions=_positions_text(campaign.designable_am1_positions),
        )
        for campaign in campaigns
    )


def _render_design_campaigns_markdown(
    *,
    artifacts: DesignCampaignArtifacts,
    backbone_manifest_rows: tuple[DesignBackboneManifestRow, ...],
    campaign_manifest_rows: tuple[DesignCampaignManifestRow, ...],
) -> str:
    lines = [
        "# Design Campaigns",
        "",
        "Phase 3A deterministic ProteinMPNN input preparation from the existing Phase 2 masks.",
        "",
        "## Derived design sets",
        "",
        _render_markdown_table(
            ("design_set", "count", "canonical_positions", "am1_positions", "definition"),
            [
                (
                    design_set.name,
                    str(len(design_set.canonical_positions)),
                    _positions_text(design_set.canonical_positions),
                    _positions_text(design_set.am1_positions),
                    design_set.description,
                )
                for design_set in artifacts.design_sets
            ],
        ),
        "",
        "## Backbone manifest",
        "",
        _render_markdown_table(
            (
                "backbone_id",
                "structure_id",
                "selected_chains",
                "design_chains",
                "fixed_context_chains",
                "design_chain_residue_count",
                "note",
            ),
            [
                (
                    row.backbone_id,
                    row.structure_id,
                    row.selected_chains,
                    row.design_chains,
                    row.fixed_context_chains or "-",
                    str(row.design_chain_residue_count),
                    row.representative_note,
                )
                for row in backbone_manifest_rows
            ],
        ),
        "",
        "## Campaign manifest",
        "",
        _render_markdown_table(
            (
                "campaign_id",
                "backbone_id",
                "design_set",
                "designed_chains",
                "fixed_context_chains",
                "designable",
                "fixed",
            ),
            [
                (
                    row.campaign_id,
                    row.backbone_id,
                    row.design_set_name,
                    row.designed_chains,
                    row.fixed_context_chains or "-",
                    str(row.designable_residue_count),
                    str(row.fixed_residue_count),
                )
                for row in campaign_manifest_rows
            ],
        ),
        "",
        "ProteinMPNN note: `fixed_positions.jsonl` uses 1-based sequence indices on the exported design chains. Use `design_campaign_positions.csv` to translate those indices back to residue numbers and canonical family positions.",
        "",
    ]
    position_rows_by_campaign: dict[str, list[DesignCampaignPositionRow]] = defaultdict(list)
    for row in artifacts.campaign_position_rows:
        if row.designable:
            position_rows_by_campaign[row.campaign_id].append(row)
    for campaign in artifacts.campaigns:
        lines.extend(
            [
                f"## {campaign.campaign_id}",
                "",
                f"- Backbone: `{campaign.backbone_id}`",
                f"- Design set: `{campaign.design_set_name}`",
                f"- Designed chains: `{','.join(campaign.designed_chains)}`",
                f"- Fixed context chains: `{','.join(campaign.fixed_context_chains) if campaign.fixed_context_chains else '-'}`",
                "",
                _render_markdown_table(
                    (
                        "chain",
                        "seq_index",
                        "residue_seq",
                        "residue_name",
                        "canonical",
                        "am1",
                        "why_designable",
                    ),
                    [
                        (
                            row.chain_id,
                            str(row.sequence_index),
                            str(row.residue_seq),
                            row.residue_name,
                            str(row.canonical_family_position),
                            str(row.am1_mature_position),
                            row.rationale,
                        )
                        for row in position_rows_by_campaign[campaign.campaign_id]
                    ],
                ),
                "",
            ]
        )
    return "\n".join(lines).rstrip() + "\n"


def _build_design_campaigns_config(
    *,
    artifacts: DesignCampaignArtifacts,
    backbone_manifest_rows: tuple[DesignBackboneManifestRow, ...],
    campaign_manifest_rows: tuple[DesignCampaignManifestRow, ...],
) -> dict[str, Any]:
    return {
        "version": 1,
        "phase": "3A",
        "source_artifacts": {
            "design_masks": "config/design_masks.yaml",
            "design_mask_candidates": "results/tables/design_mask_candidates.csv",
            "cross_template_alignment": "results/tables/cross_template_residue_alignment.csv",
            "template_chain_summary": "results/tables/template_chain_summary.csv",
        },
        "design_sets": {
            design_set.name: {
                "canonical_family_positions": list(design_set.canonical_positions),
                "am1_mature_positions": list(design_set.am1_positions),
                "definition": design_set.description,
            }
            for design_set in artifacts.design_sets
        },
        "backbones": {
            row.backbone_id: {
                "structure_id": row.structure_id,
                "source_kind": row.source_kind,
                "source_path": row.source_path,
                "experimental_method": row.experimental_method,
                "selected_chains": row.selected_chains.split(",") if row.selected_chains else [],
                "design_chains": row.design_chains.split(",") if row.design_chains else [],
                "fixed_context_chains": row.fixed_context_chains.split(",") if row.fixed_context_chains else [],
                "design_chain_residue_count": row.design_chain_residue_count,
                "design_chain_metal_site_count": row.design_chain_metal_site_count,
                "representative_note": row.representative_note,
                "backbone_pdb_path": row.backbone_pdb_path,
            }
            for row in backbone_manifest_rows
        },
        "campaigns": {
            row.campaign_id: {
                "backbone_id": row.backbone_id,
                "design_set_name": row.design_set_name,
                "designed_chains": row.designed_chains.split(",") if row.designed_chains else [],
                "fixed_context_chains": row.fixed_context_chains.split(",") if row.fixed_context_chains else [],
                "designable_residue_count": row.designable_residue_count,
                "fixed_residue_count": row.fixed_residue_count,
                "designable_canonical_positions": [
                    int(position)
                    for position in row.designable_canonical_positions.split(",")
                    if position
                ],
                "designable_am1_positions": [
                    int(position)
                    for position in row.designable_am1_positions.split(",")
                    if position
                ],
                "exported_pdb_path": row.exported_pdb_path,
                "chain_assignment_path": row.chain_assignment_path,
                "fixed_positions_path": row.fixed_positions_path,
            }
            for row in campaign_manifest_rows
        },
    }


def write_design_campaign_outputs(
    *,
    artifacts: DesignCampaignArtifacts,
    report_path: Path,
    campaign_positions_path: Path,
    backbone_manifest_path: Path,
    campaign_manifest_path: Path,
    config_path: Path,
    proteinmpnn_root: Path,
) -> None:
    """Write Phase 3A design campaign outputs to disk."""
    backbone_paths, campaign_paths = _build_output_path_map(
        backbones=artifacts.backbones,
        campaigns=artifacts.campaigns,
        proteinmpnn_root=proteinmpnn_root,
    )
    backbone_manifest_rows = _build_backbone_manifest_rows(
        backbones=artifacts.backbones,
        backbone_paths=backbone_paths,
    )
    campaign_manifest_rows = _build_campaign_manifest_rows(
        campaigns=artifacts.campaigns,
        campaign_paths=campaign_paths,
    )
    for backbone in artifacts.backbones:
        atomic_write_text(
            backbone_paths[backbone.backbone_id],
            artifacts.backbone_pdb_text_by_id[backbone.backbone_id],
        )
    backbone_map = {
        backbone.backbone_id: backbone
        for backbone in artifacts.backbones
    }
    for campaign in artifacts.campaigns:
        campaign_dir = campaign_paths[campaign.campaign_id]["directory"]
        campaign_dir.mkdir(parents=True, exist_ok=True)
        atomic_write_text(
            campaign_paths[campaign.campaign_id]["pdb"],
            artifacts.backbone_pdb_text_by_id[campaign.backbone_id],
        )
        write_jsonl(
            campaign_paths[campaign.campaign_id]["chain_assignment"],
            [campaign.chain_assignment_payload],
        )
        write_jsonl(
            campaign_paths[campaign.campaign_id]["fixed_positions"],
            [campaign.fixed_positions_payload],
        )
        backbone = backbone_map[campaign.backbone_id]
        if tuple(campaign.designed_chains) != tuple(backbone.design_chains):
            raise ValueError(
                f"Campaign {campaign.campaign_id} design chains do not match backbone {campaign.backbone_id}"
            )
    write_csv_rows(backbone_manifest_path, backbone_manifest_rows)
    write_csv_rows(campaign_manifest_path, campaign_manifest_rows)
    write_csv_rows(campaign_positions_path, artifacts.campaign_position_rows)
    write_yaml(
        config_path,
        _build_design_campaigns_config(
            artifacts=artifacts,
            backbone_manifest_rows=backbone_manifest_rows,
            campaign_manifest_rows=campaign_manifest_rows,
        ),
    )
    atomic_write_text(
        report_path,
        _render_design_campaigns_markdown(
            artifacts=artifacts,
            backbone_manifest_rows=backbone_manifest_rows,
            campaign_manifest_rows=campaign_manifest_rows,
        ),
    )
