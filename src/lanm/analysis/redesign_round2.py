"""Phase 7A focused round-2 redesign planning from metadynamics-informed seed scaffolds."""

from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

try:
    import yaml  # type: ignore
except ModuleNotFoundError:  # pragma: no cover - optional dependency in local envs
    yaml = None

from lanm.analysis.design_campaigns import load_design_mask_candidates
from lanm.analysis.md_panel_selection import load_ligandmpnn_shortlist_rows
from lanm.analysis.proteinmpnn_candidates import load_design_campaign_position_rows
from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_text, copy_if_changed, write_csv_rows, write_jsonl, write_yaml
from lanm.models import DesignCampaignPositionRow, DesignMaskCandidate, LigandMPNNShortlistRow
from lanm.paths import REPO_ROOT
from lanm.structure.residues import AMINO_ACID_CODES

ROUND2_SEED_IDS = (
    "am1_mex_ss_only_u02",
    "hans_pocket_ss_only_u02",
    "hans_interface_ss_plus_if_u04",
)
EXPECTED_DY_STATUS = "retained_bound"
EXPECTED_COMPETITOR_STATUS = "persistent_capture"


@dataclass(frozen=True, slots=True)
class OpenMMMetadynamicsSummaryRecord:
    panel_member_id: str
    target_metal: str
    success_status: str
    failure_reason: str
    minimum_coordination_number: float | None
    maximum_mean_metal_oxygen_distance_A: float | None
    escape_event_detected: bool
    final_coordination_number: float | None
    final_mean_metal_oxygen_distance_A: float | None
    metadynamics_status: str


@dataclass(frozen=True, slots=True)
class OpenMMMetadynamicsPanelStatusRecord:
    panel_member_id: str
    dy_status: str
    nd_status: str
    y_status: str
    al_status: str
    fe_status: str
    candidate_keep_for_qm: bool


@dataclass(frozen=True, slots=True)
class DesignCampaignManifestRecord:
    campaign_id: str
    backbone_id: str
    design_set_name: str
    designed_chains: tuple[str, ...]
    fixed_context_chains: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class MDValidationPanelRecord:
    panel_rank: int
    panel_member_id: str
    panel_member_type: str
    panel_role: str
    candidate_id: str
    campaign_id: str
    backbone_id: str
    topology_class: str
    starting_structure_path: str
    selection_reason: str


@dataclass(frozen=True, slots=True)
class SeedResolvedPosition:
    chain_id: str
    sequence_index: int
    residue_seq: int
    insertion_code: str
    residue_name: str
    seed_amino_acid: str
    canonical_family_position: int | None
    am1_mature_position: int | None
    fixed_first_shell: bool
    protected_positions: bool
    mutable_second_sphere: bool
    mutable_interface: bool

    @property
    def residue_id(self) -> str:
        return f"{self.chain_id}{self.residue_seq}{self.insertion_code}".strip()

    @property
    def native_amino_acid(self) -> str:
        return AMINO_ACID_CODES.get(self.residue_name, "X")


@dataclass(frozen=True, slots=True)
class Round2SeedContext:
    seed_rank: int
    candidate_id: str
    campaign_id: str
    backbone_id: str
    design_set_name: str
    topology_class: str
    designed_chains: tuple[str, ...]
    fixed_context_chains: tuple[str, ...]
    starting_structure_path: Path
    designed_chain_sequence: str
    panel_rank: int
    panel_role: str
    selection_reason: str
    dy_status: str
    nd_status: str
    y_status: str
    al_status: str
    fe_status: str
    candidate_keep_for_qm: bool
    seed_mutation_count: int
    seed_mutation_tokens: tuple[str, ...]
    seed_mutation_residue_ids: tuple[str, ...]
    seed_mutation_canonical_positions: tuple[int, ...]
    seed_mutation_am1_positions: tuple[int, ...]
    resolved_positions: tuple[SeedResolvedPosition, ...]


@dataclass(frozen=True, slots=True)
class Round2PositionRow:
    seed_rank: int
    candidate_id: str
    campaign_id: str
    backbone_id: str
    topology_class: str
    chain_id: str
    sequence_index: int
    residue_seq: int
    insertion_code: str
    residue_id: str
    residue_name: str
    native_amino_acid: str
    seed_amino_acid: str
    canonical_family_position: int
    am1_mature_position: int
    mutable_second_sphere: bool
    mutable_interface: bool
    near_seed_mutation: bool
    near_known_interface_site: bool
    is_seed_mutation_site: bool
    selection_reason: str


@dataclass(frozen=True, slots=True)
class Round2SeedRow:
    seed_rank: int
    candidate_id: str
    campaign_id: str
    backbone_id: str
    topology_class: str
    panel_rank: int
    panel_role: str
    designed_chains: str
    fixed_context_chains: str
    dy_status: str
    nd_status: str
    y_status: str
    al_status: str
    fe_status: str
    candidate_keep_for_qm: bool
    starting_structure_path: str
    exported_backbone_path: str
    exported_pdb_path: str
    chain_assignment_path: str
    fixed_positions_path: str
    seed_mutation_count: int
    seed_mutation_tokens: str
    seed_mutation_residue_ids: str
    seed_mutation_canonical_positions: str
    seed_mutation_am1_positions: str
    redesignable_position_count: int
    redesignable_residue_ids: str
    redesignable_canonical_positions: str
    redesignable_am1_positions: str
    selection_reason: str


@dataclass(frozen=True, slots=True)
class Round2SeedExportPaths:
    backbone_path: Path
    campaign_dir: Path
    exported_pdb_path: Path
    chain_assignment_path: Path
    fixed_positions_path: Path


@dataclass(frozen=True, slots=True)
class Round2Artifacts:
    interface_window_size: int
    seed_contexts: tuple[Round2SeedContext, ...]
    seed_rows: tuple[Round2SeedRow, ...]
    position_rows: tuple[Round2PositionRow, ...]
    config_payload: dict[str, Any]
    report_markdown: str


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


def _parse_bool(value: str) -> bool:
    return value.strip().lower() == "true"


def _parse_optional_float(value: str) -> float | None:
    stripped = value.strip()
    return float(stripped) if stripped else None


def _split_csv_list(value: str) -> tuple[str, ...]:
    if not value.strip():
        return ()
    return tuple(item.strip() for item in value.split(",") if item.strip())


def _positions_text(positions: Sequence[int]) -> str:
    return ",".join(str(position) for position in positions)


def _load_yaml_mapping(path: Path) -> dict[str, Any]:
    require_path(path)
    if yaml is None:
        raise RuntimeError("PyYAML is required to read redesign_round2 configuration inputs")
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Unexpected YAML payload in {path}")
    return payload


def load_openmm_metadynamics_summary_rows(path: Path) -> tuple[OpenMMMetadynamicsSummaryRecord, ...]:
    """Load Phase 6C1 metadynamics system-level summary rows."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            OpenMMMetadynamicsSummaryRecord(
                panel_member_id=str(row["panel_member_id"]).strip(),
                target_metal=str(row["target_metal"]).strip(),
                success_status=str(row["success_status"]).strip(),
                failure_reason=str(row["failure_reason"]).strip(),
                minimum_coordination_number=_parse_optional_float(str(row["minimum_coordination_number"])),
                maximum_mean_metal_oxygen_distance_A=_parse_optional_float(
                    str(row["maximum_mean_metal_oxygen_distance_A"])
                ),
                escape_event_detected=_parse_bool(str(row["escape_event_detected"])),
                final_coordination_number=_parse_optional_float(str(row["final_coordination_number"])),
                final_mean_metal_oxygen_distance_A=_parse_optional_float(
                    str(row["final_mean_metal_oxygen_distance_A"])
                ),
                metadynamics_status=str(row["metadynamics_status"]).strip(),
            )
            for row in reader
        ]
    return tuple(sorted(rows, key=lambda row: (row.panel_member_id, row.target_metal)))


def load_openmm_metadynamics_panel_status_rows(path: Path) -> tuple[OpenMMMetadynamicsPanelStatusRecord, ...]:
    """Load Phase 6C1 metadynamics panel-status rows."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            OpenMMMetadynamicsPanelStatusRecord(
                panel_member_id=str(row["panel_member_id"]).strip(),
                dy_status=str(row["dy_status"]).strip(),
                nd_status=str(row["nd_status"]).strip(),
                y_status=str(row["y_status"]).strip(),
                al_status=str(row["al_status"]).strip(),
                fe_status=str(row["fe_status"]).strip(),
                candidate_keep_for_qm=_parse_bool(str(row["candidate_keep_for_qm"])),
            )
            for row in reader
        ]
    return tuple(sorted(rows, key=lambda row: row.panel_member_id))


def load_design_campaign_manifest_rows(path: Path) -> tuple[DesignCampaignManifestRecord, ...]:
    """Load Phase 3A design campaign manifest rows needed for round-2 exports."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            DesignCampaignManifestRecord(
                campaign_id=str(row["campaign_id"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                design_set_name=str(row["design_set_name"]).strip(),
                designed_chains=_split_csv_list(str(row["designed_chains"])),
                fixed_context_chains=_split_csv_list(str(row["fixed_context_chains"])),
            )
            for row in reader
        ]
    return tuple(sorted(rows, key=lambda row: row.campaign_id))


def load_md_validation_panel_rows(path: Path) -> tuple[MDValidationPanelRecord, ...]:
    """Load Phase 5B MD validation panel rows."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            MDValidationPanelRecord(
                panel_rank=int(str(row["panel_rank"]).strip()),
                panel_member_id=str(row["panel_member_id"]).strip(),
                panel_member_type=str(row["panel_member_type"]).strip(),
                panel_role=str(row["panel_role"]).strip(),
                candidate_id=str(row["candidate_id"]).strip(),
                campaign_id=str(row["campaign_id"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                topology_class=str(row["topology_class"]).strip(),
                starting_structure_path=str(row["starting_structure_path"]).strip(),
                selection_reason=str(row["selection_reason"]).strip(),
            )
            for row in reader
        ]
    return tuple(sorted(rows, key=lambda row: row.panel_rank))


def validate_design_mask_config(
    *,
    design_mask_rows: Sequence[DesignMaskCandidate],
    design_masks_config: dict[str, Any],
) -> int:
    """Validate that design_masks.yaml matches design_mask_candidates.csv and return the interface window."""
    payload_masks = design_masks_config.get("masks")
    if not isinstance(payload_masks, dict):
        raise ValueError("design_masks.yaml is missing a top-level masks mapping")

    expected_masks = {
        "fixed_first_shell": "fixed_first_shell",
        "mutable_second_sphere": "mutable_second_sphere",
        "mutable_interface": "mutable_interface",
        "protected_positions": "protected_positions",
    }
    for mask_name, attribute_name in expected_masks.items():
        mask_payload = payload_masks.get(mask_name)
        if not isinstance(mask_payload, dict):
            raise ValueError(f"design_masks.yaml is missing masks.{mask_name}")
        yaml_positions = tuple(
            sorted(
                int(position)
                for position in mask_payload.get("canonical_family_positions", [])
            )
        )
        csv_positions = tuple(
            sorted(
                row.canonical_family_position
                for row in design_mask_rows
                if getattr(row, attribute_name)
            )
        )
        if yaml_positions != csv_positions:
            raise ValueError(
                f"design_masks.yaml masks.{mask_name}.canonical_family_positions did not match "
                "design_mask_candidates.csv"
            )

    index_payload = design_masks_config.get("index")
    if not isinstance(index_payload, dict):
        raise ValueError("design_masks.yaml is missing index settings")
    interface_window_size = int(index_payload.get("interface_neighborhood_window", 3))
    if interface_window_size != 3:
        raise ValueError(
            f"Phase 7A requires interface_neighborhood_window=3, found {interface_window_size}"
        )
    return interface_window_size


def _split_designed_sequence(sequence: str, designed_chains: Sequence[str]) -> tuple[str, ...]:
    if len(designed_chains) == 1:
        return (sequence.replace("/", ""),)
    parts = tuple(part.strip() for part in sequence.split("/"))
    if len(parts) != len(designed_chains):
        raise ValueError(
            f"Expected {len(designed_chains)} designed-chain sequence parts, found {len(parts)} in {sequence!r}"
        )
    return tuple(part.replace("/", "") for part in parts)


def _group_campaign_positions(
    position_rows: Sequence[DesignCampaignPositionRow],
    designed_chains: Sequence[str],
) -> dict[str, tuple[DesignCampaignPositionRow, ...]]:
    rows_by_chain: dict[str, list[DesignCampaignPositionRow]] = defaultdict(list)
    for row in position_rows:
        if row.chain_id not in designed_chains:
            continue
        rows_by_chain[row.chain_id].append(row)
    grouped: dict[str, tuple[DesignCampaignPositionRow, ...]] = {}
    for chain_id in designed_chains:
        chain_rows = tuple(sorted(rows_by_chain[chain_id], key=lambda row: row.sequence_index))
        if not chain_rows:
            raise ValueError(f"No campaign positions found for designed chain {chain_id}")
        grouped[chain_id] = chain_rows
    return grouped


def _resolve_seed_positions(
    *,
    designed_chains: Sequence[str],
    designed_chain_sequence: str,
    position_rows: Sequence[DesignCampaignPositionRow],
) -> tuple[tuple[SeedResolvedPosition, ...], tuple[str, ...], tuple[str, ...], tuple[int, ...], tuple[int, ...]]:
    grouped_rows = _group_campaign_positions(position_rows, designed_chains)
    sequence_parts = _split_designed_sequence(designed_chain_sequence, designed_chains)

    resolved_positions: list[SeedResolvedPosition] = []
    mutation_tokens: list[str] = []
    mutation_residue_ids: list[str] = []
    mutation_canonical_positions: list[int] = []
    mutation_am1_positions: list[int] = []

    for chain_id, sequence_part in zip(designed_chains, sequence_parts):
        chain_rows = grouped_rows[chain_id]
        if len(chain_rows) != len(sequence_part):
            raise ValueError(
                f"Designed sequence length {len(sequence_part)} did not match campaign positions "
                f"{len(chain_rows)} for chain {chain_id}"
            )
        for row, seed_amino_acid in zip(chain_rows, sequence_part):
            resolved_position = SeedResolvedPosition(
                chain_id=row.chain_id,
                sequence_index=row.sequence_index,
                residue_seq=row.residue_seq,
                insertion_code=row.insertion_code,
                residue_name=row.residue_name,
                seed_amino_acid=seed_amino_acid,
                canonical_family_position=row.canonical_family_position,
                am1_mature_position=row.am1_mature_position,
                fixed_first_shell=row.fixed_first_shell,
                protected_positions=row.protected_positions,
                mutable_second_sphere=row.mutable_second_sphere,
                mutable_interface=row.mutable_interface,
            )
            resolved_positions.append(resolved_position)
            native_amino_acid = resolved_position.native_amino_acid
            if native_amino_acid == seed_amino_acid:
                continue
            mutation_tokens.append(f"{native_amino_acid}{row.sequence_index}{seed_amino_acid}")
            mutation_residue_ids.append(resolved_position.residue_id)
            if (
                row.canonical_family_position is not None
                and row.canonical_family_position not in mutation_canonical_positions
            ):
                mutation_canonical_positions.append(row.canonical_family_position)
            if row.am1_mature_position is not None and row.am1_mature_position not in mutation_am1_positions:
                mutation_am1_positions.append(row.am1_mature_position)

    return (
        tuple(resolved_positions),
        tuple(mutation_tokens),
        tuple(mutation_residue_ids),
        tuple(mutation_canonical_positions),
        tuple(mutation_am1_positions),
    )


def discover_round2_seed_contexts(
    *,
    seed_candidate_ids: Sequence[str],
    panel_status_rows: Sequence[OpenMMMetadynamicsPanelStatusRecord],
    summary_rows: Sequence[OpenMMMetadynamicsSummaryRecord],
    ligand_shortlist_rows: Sequence[LigandMPNNShortlistRow],
    md_panel_rows: Sequence[MDValidationPanelRecord],
    campaign_manifest_rows: Sequence[DesignCampaignManifestRecord],
    design_campaign_position_rows: Sequence[DesignCampaignPositionRow],
) -> tuple[Round2SeedContext, ...]:
    """Resolve the requested round-2 seed scaffolds onto their validated structural contexts."""
    panel_status_by_id = {
        row.panel_member_id: row
        for row in panel_status_rows
    }
    summary_by_seed_metal = {
        (row.panel_member_id, row.target_metal): row
        for row in summary_rows
    }
    ligand_by_candidate_id = {
        row.candidate_id: row
        for row in ligand_shortlist_rows
    }
    md_panel_by_candidate_id = {
        row.candidate_id: row
        for row in md_panel_rows
    }
    campaign_by_id = {
        row.campaign_id: row
        for row in campaign_manifest_rows
    }
    positions_by_campaign: dict[str, list[DesignCampaignPositionRow]] = defaultdict(list)
    for row in design_campaign_position_rows:
        positions_by_campaign[row.campaign_id].append(row)

    seed_contexts: list[Round2SeedContext] = []
    for seed_rank, candidate_id in enumerate(seed_candidate_ids, start=1):
        panel_status = panel_status_by_id.get(candidate_id)
        if panel_status is None:
            raise ValueError(f"No openmm_metadynamics_panel_status.csv entry found for {candidate_id}")
        ligand_row = ligand_by_candidate_id.get(candidate_id)
        if ligand_row is None:
            raise ValueError(f"No ligandmpnn_shortlist.csv entry found for {candidate_id}")
        md_panel_row = md_panel_by_candidate_id.get(candidate_id)
        if md_panel_row is None:
            raise ValueError(f"No md_validation_panel.csv entry found for {candidate_id}")
        if md_panel_row.panel_member_type != "designed_candidate":
            raise ValueError(f"Round-2 seed {candidate_id} is not a designed_candidate in the MD panel")
        campaign_row = campaign_by_id.get(ligand_row.campaign_id)
        if campaign_row is None:
            raise ValueError(f"No design_campaign_manifest.csv entry found for {ligand_row.campaign_id}")
        if campaign_row.backbone_id != ligand_row.backbone_id:
            raise ValueError(
                f"Backbone mismatch for {candidate_id}: {campaign_row.backbone_id} vs {ligand_row.backbone_id}"
            )
        if md_panel_row.campaign_id != ligand_row.campaign_id or md_panel_row.backbone_id != ligand_row.backbone_id:
            raise ValueError(
                f"MD panel context mismatch for {candidate_id}: "
                f"{md_panel_row.campaign_id}/{md_panel_row.backbone_id} vs "
                f"{ligand_row.campaign_id}/{ligand_row.backbone_id}"
            )
        dy_summary_row = summary_by_seed_metal.get((candidate_id, "Dy"))
        al_summary_row = summary_by_seed_metal.get((candidate_id, "Al"))
        fe_summary_row = summary_by_seed_metal.get((candidate_id, "Fe"))
        if dy_summary_row is None or al_summary_row is None or fe_summary_row is None:
            raise ValueError(f"Missing Dy/Al/Fe metadynamics rows for {candidate_id}")
        if dy_summary_row.metadynamics_status != panel_status.dy_status:
            raise ValueError(f"Dy status mismatch between summary and panel status for {candidate_id}")
        if al_summary_row.metadynamics_status != panel_status.al_status:
            raise ValueError(f"Al status mismatch between summary and panel status for {candidate_id}")
        if fe_summary_row.metadynamics_status != panel_status.fe_status:
            raise ValueError(f"Fe status mismatch between summary and panel status for {candidate_id}")
        if panel_status.dy_status != EXPECTED_DY_STATUS:
            raise ValueError(f"{candidate_id} did not retain Dy under metadynamics")
        if panel_status.al_status != EXPECTED_COMPETITOR_STATUS:
            raise ValueError(f"{candidate_id} did not show persistent Al capture")
        if panel_status.fe_status != EXPECTED_COMPETITOR_STATUS:
            raise ValueError(f"{candidate_id} did not show persistent Fe capture")

        starting_structure_path = _repo_path(md_panel_row.starting_structure_path)
        require_path(starting_structure_path)
        if _display_path(starting_structure_path) != ligand_row.packed_pdb_path:
            raise ValueError(
                f"Packed PDB mismatch for {candidate_id}: "
                f"{_display_path(starting_structure_path)} vs {ligand_row.packed_pdb_path}"
            )

        position_rows = positions_by_campaign.get(ligand_row.campaign_id)
        if not position_rows:
            raise ValueError(f"No design_campaign_positions.csv rows found for {ligand_row.campaign_id}")
        resolved_positions, mutation_tokens, mutation_residue_ids, mutation_canonical_positions, mutation_am1_positions = (
            _resolve_seed_positions(
                designed_chains=campaign_row.designed_chains,
                designed_chain_sequence=ligand_row.designed_chain_sequence,
                position_rows=position_rows,
            )
        )

        seed_contexts.append(
            Round2SeedContext(
                seed_rank=seed_rank,
                candidate_id=candidate_id,
                campaign_id=ligand_row.campaign_id,
                backbone_id=ligand_row.backbone_id,
                design_set_name=campaign_row.design_set_name,
                topology_class=md_panel_row.topology_class,
                designed_chains=campaign_row.designed_chains,
                fixed_context_chains=campaign_row.fixed_context_chains,
                starting_structure_path=starting_structure_path,
                designed_chain_sequence=ligand_row.designed_chain_sequence,
                panel_rank=md_panel_row.panel_rank,
                panel_role=md_panel_row.panel_role,
                selection_reason=md_panel_row.selection_reason,
                dy_status=panel_status.dy_status,
                nd_status=panel_status.nd_status,
                y_status=panel_status.y_status,
                al_status=panel_status.al_status,
                fe_status=panel_status.fe_status,
                candidate_keep_for_qm=panel_status.candidate_keep_for_qm,
                seed_mutation_count=len(mutation_tokens),
                seed_mutation_tokens=mutation_tokens,
                seed_mutation_residue_ids=mutation_residue_ids,
                seed_mutation_canonical_positions=mutation_canonical_positions,
                seed_mutation_am1_positions=mutation_am1_positions,
                resolved_positions=resolved_positions,
            )
        )

    return tuple(seed_contexts)


def select_round2_redesign_positions(
    *,
    seed_context: Round2SeedContext,
    design_mask_rows: Sequence[DesignMaskCandidate],
    interface_window_size: int,
) -> tuple[Round2PositionRow, ...]:
    """Select the intentionally narrow round-2 redesign positions for one seed scaffold."""
    interface_anchor_positions = {
        row.canonical_family_position
        for row in design_mask_rows
        if row.mutable_interface
    }
    seed_anchor_positions = set(seed_context.seed_mutation_canonical_positions)

    selected_rows: list[Round2PositionRow] = []
    for position in seed_context.resolved_positions:
        canonical_position = position.canonical_family_position
        am1_position = position.am1_mature_position
        if canonical_position is None or am1_position is None:
            continue
        if position.fixed_first_shell or position.protected_positions:
            continue
        if not (position.mutable_second_sphere or position.mutable_interface):
            continue
        near_seed_mutation = any(
            abs(canonical_position - anchor_position) <= interface_window_size
            for anchor_position in seed_anchor_positions
        )
        near_interface_site = any(
            abs(canonical_position - anchor_position) <= interface_window_size
            for anchor_position in interface_anchor_positions
        )
        if not (near_seed_mutation or near_interface_site):
            continue
        if near_seed_mutation and near_interface_site:
            selection_reason = "seed_mutation_window,known_interface_window"
        elif near_seed_mutation:
            selection_reason = "seed_mutation_window"
        else:
            selection_reason = "known_interface_window"
        selected_rows.append(
            Round2PositionRow(
                seed_rank=seed_context.seed_rank,
                candidate_id=seed_context.candidate_id,
                campaign_id=seed_context.campaign_id,
                backbone_id=seed_context.backbone_id,
                topology_class=seed_context.topology_class,
                chain_id=position.chain_id,
                sequence_index=position.sequence_index,
                residue_seq=position.residue_seq,
                insertion_code=position.insertion_code,
                residue_id=position.residue_id,
                residue_name=position.residue_name,
                native_amino_acid=position.native_amino_acid,
                seed_amino_acid=position.seed_amino_acid,
                canonical_family_position=canonical_position,
                am1_mature_position=am1_position,
                mutable_second_sphere=position.mutable_second_sphere,
                mutable_interface=position.mutable_interface,
                near_seed_mutation=near_seed_mutation,
                near_known_interface_site=near_interface_site,
                is_seed_mutation_site=canonical_position in seed_anchor_positions,
                selection_reason=selection_reason,
            )
        )
    return tuple(sorted(selected_rows, key=lambda row: (row.seed_rank, row.chain_id, row.sequence_index)))


def _build_seed_export_paths(
    *,
    seed_context: Round2SeedContext,
    proteinmpnn_round2_root: Path,
) -> Round2SeedExportPaths:
    campaign_dir = proteinmpnn_round2_root / "campaigns" / seed_context.candidate_id
    return Round2SeedExportPaths(
        backbone_path=proteinmpnn_round2_root / "backbones" / f"{seed_context.candidate_id}.pdb",
        campaign_dir=campaign_dir,
        exported_pdb_path=campaign_dir / f"{seed_context.candidate_id}.pdb",
        chain_assignment_path=campaign_dir / "chain_id.jsonl",
        fixed_positions_path=campaign_dir / "fixed_positions.jsonl",
    )


def _build_seed_rows(
    *,
    seed_contexts: Sequence[Round2SeedContext],
    position_rows: Sequence[Round2PositionRow],
    proteinmpnn_round2_root: Path,
) -> tuple[Round2SeedRow, ...]:
    positions_by_candidate: dict[str, list[Round2PositionRow]] = defaultdict(list)
    for row in position_rows:
        positions_by_candidate[row.candidate_id].append(row)

    seed_rows: list[Round2SeedRow] = []
    for seed_context in seed_contexts:
        selected_positions = tuple(
            sorted(
                positions_by_candidate[seed_context.candidate_id],
                key=lambda row: (row.chain_id, row.sequence_index),
            )
        )
        export_paths = _build_seed_export_paths(
            seed_context=seed_context,
            proteinmpnn_round2_root=proteinmpnn_round2_root,
        )
        seed_rows.append(
            Round2SeedRow(
                seed_rank=seed_context.seed_rank,
                candidate_id=seed_context.candidate_id,
                campaign_id=seed_context.campaign_id,
                backbone_id=seed_context.backbone_id,
                topology_class=seed_context.topology_class,
                panel_rank=seed_context.panel_rank,
                panel_role=seed_context.panel_role,
                designed_chains=",".join(seed_context.designed_chains),
                fixed_context_chains=",".join(seed_context.fixed_context_chains),
                dy_status=seed_context.dy_status,
                nd_status=seed_context.nd_status,
                y_status=seed_context.y_status,
                al_status=seed_context.al_status,
                fe_status=seed_context.fe_status,
                candidate_keep_for_qm=seed_context.candidate_keep_for_qm,
                starting_structure_path=_display_path(seed_context.starting_structure_path),
                exported_backbone_path=_display_path(export_paths.backbone_path),
                exported_pdb_path=_display_path(export_paths.exported_pdb_path),
                chain_assignment_path=_display_path(export_paths.chain_assignment_path),
                fixed_positions_path=_display_path(export_paths.fixed_positions_path),
                seed_mutation_count=seed_context.seed_mutation_count,
                seed_mutation_tokens=",".join(seed_context.seed_mutation_tokens),
                seed_mutation_residue_ids=",".join(seed_context.seed_mutation_residue_ids),
                seed_mutation_canonical_positions=_positions_text(seed_context.seed_mutation_canonical_positions),
                seed_mutation_am1_positions=_positions_text(seed_context.seed_mutation_am1_positions),
                redesignable_position_count=len(selected_positions),
                redesignable_residue_ids=",".join(row.residue_id for row in selected_positions),
                redesignable_canonical_positions=_positions_text(
                    tuple(row.canonical_family_position for row in selected_positions)
                ),
                redesignable_am1_positions=_positions_text(tuple(row.am1_mature_position for row in selected_positions)),
                selection_reason=seed_context.selection_reason,
            )
        )
    return tuple(seed_rows)


def _render_markdown_table(headers: tuple[str, ...], rows: Sequence[tuple[str, ...]]) -> str:
    header_line = "| " + " | ".join(headers) + " |"
    divider_line = "| " + " | ".join(["---"] * len(headers)) + " |"
    body_lines = ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join([header_line, divider_line, *body_lines])


def _render_round2_report(
    *,
    seed_contexts: Sequence[Round2SeedContext],
    seed_rows: Sequence[Round2SeedRow],
    position_rows: Sequence[Round2PositionRow],
    interface_window_size: int,
) -> str:
    positions_by_candidate: dict[str, list[Round2PositionRow]] = defaultdict(list)
    for row in position_rows:
        positions_by_candidate[row.candidate_id].append(row)

    seed_rows_by_candidate = {
        row.candidate_id: row
        for row in seed_rows
    }
    lines = [
        "# Round-2 Focused Redesign Plan",
        "",
        "Phase 7A prepares deterministic ProteinMPNN-ready campaign inputs for a narrow redesign round centered on the three metadynamics-validated seed scaffolds.",
        "No ProteinMPNN, LigandMPNN, Rosetta, MD, QM, or quantum execution is started here.",
        "",
        "## Why QM Is Not Started Yet",
        "",
        "- Phase 6C1 still classifies all three selected seed scaffolds as `candidate_keep_for_qm=FALSE`, so there is no metadynamics-backed signal yet to escalate any design directly into QM.",
        "- The immediate failure mode is shared across the seeds: `Dy` remains `retained_bound`, while both `Al` and `Fe` remain `persistent_capture`, so the next efficient step is to narrow sequence space around the tolerated local neighborhoods before spending on QM.",
        "- This phase therefore only exports a deterministic round-2 redesign plan and leaves all expensive downstream methods for later phases.",
        "",
        "## Why These Three Seed Scaffolds",
        "",
        _render_markdown_table(
            (
                "seed",
                "topology",
                "panel_role",
                "Dy",
                "Al",
                "Fe",
                "seed_mutations",
                "why_kept",
            ),
            [
                (
                    seed.candidate_id,
                    seed.topology_class,
                    seed.panel_role,
                    seed.dy_status,
                    seed.al_status,
                    seed.fe_status,
                    seed_rows_by_candidate[seed.candidate_id].seed_mutation_canonical_positions or "-",
                    seed.selection_reason,
                )
                for seed in seed_contexts
            ],
        ),
        "",
        "These three seeds span the AM1/Mex monomer, the Hans pocket-focused monomer, and the Hans interface-aware multichain context while preserving the exact metadynamics phenotype that motivates round 2: Dy retention without relief of persistent Al/Fe capture.",
        "",
        "## Why The Scope Is Intentionally Narrow",
        "",
        f"- `fixed_first_shell` and `protected_positions` stay fixed for every seed, so the Dy-retaining first-shell architecture and numbering-protected positions are not reopened.",
        f"- Redesignable sites must already be marked `mutable_second_sphere` or `mutable_interface`, must remain inside mature AM1 numbering, and must lie within +/-{interface_window_size} canonical positions of a seed mutation or of a known interface-aware mutable site.",
        "- This keeps round 2 focused on local neighborhoods most likely to tune competitor capture without destabilizing the Dy-compatible scaffold cores.",
        "",
        "## Planned Round-2 Campaigns",
        "",
        _render_markdown_table(
            (
                "seed",
                "designable_count",
                "designable_canonical_positions",
                "designable_residue_ids",
            ),
            [
                (
                    row.candidate_id,
                    str(row.redesignable_position_count),
                    row.redesignable_canonical_positions or "-",
                    row.redesignable_residue_ids or "-",
                )
                for row in seed_rows
            ],
        ),
        "",
    ]
    for seed in seed_contexts:
        selected_positions = tuple(
            sorted(
                positions_by_candidate[seed.candidate_id],
                key=lambda row: (row.chain_id, row.sequence_index),
            )
        )
        lines.extend(
            [
                f"## {seed.candidate_id}",
                "",
                f"- Starting structure: `{_display_path(seed.starting_structure_path)}`",
                f"- Seed mutation canonical positions: `{seed_rows_by_candidate[seed.candidate_id].seed_mutation_canonical_positions or '-'}`",
                f"- Redesignable canonical positions: `{seed_rows_by_candidate[seed.candidate_id].redesignable_canonical_positions or '-'}`",
                "",
                _render_markdown_table(
                    (
                        "residue_id",
                        "canonical",
                        "am1",
                        "seed_aa",
                        "mutable_second_sphere",
                        "mutable_interface",
                        "reason",
                    ),
                    [
                        (
                            row.residue_id,
                            str(row.canonical_family_position),
                            str(row.am1_mature_position),
                            row.seed_amino_acid,
                            str(row.mutable_second_sphere),
                            str(row.mutable_interface),
                            row.selection_reason,
                        )
                        for row in selected_positions
                    ],
                ),
                "",
            ]
        )
    return "\n".join(lines).rstrip() + "\n"


def _build_round2_config_payload(
    *,
    seed_contexts: Sequence[Round2SeedContext],
    seed_rows: Sequence[Round2SeedRow],
    interface_window_size: int,
    source_artifacts: dict[str, str],
) -> dict[str, Any]:
    seed_rows_by_candidate = {
        row.candidate_id: row
        for row in seed_rows
    }
    return {
        "version": 1,
        "phase": "7A",
        "source_artifacts": source_artifacts,
        "selection_policy": {
            "seed_candidate_ids": [seed.candidate_id for seed in seed_contexts],
            "keep_fixed_first_shell_fixed": True,
            "keep_protected_positions_fixed": True,
            "interface_window_size": interface_window_size,
            "redesign_position_requirements": [
                "originally mutable_second_sphere or mutable_interface",
                "within +/-3 canonical positions of a seed mutation or a known interface-aware mutable site",
                "not outside AM1 mature numbering",
            ],
            "qm_started": False,
            "external_execution_started": False,
        },
        "seeds": {
            seed.candidate_id: {
                "campaign_id": seed.campaign_id,
                "backbone_id": seed.backbone_id,
                "design_set_name": seed.design_set_name,
                "topology_class": seed.topology_class,
                "panel_rank": seed.panel_rank,
                "panel_role": seed.panel_role,
                "designed_chains": list(seed.designed_chains),
                "fixed_context_chains": list(seed.fixed_context_chains),
                "starting_structure_path": seed_rows_by_candidate[seed.candidate_id].starting_structure_path,
                "dy_status": seed.dy_status,
                "nd_status": seed.nd_status,
                "y_status": seed.y_status,
                "al_status": seed.al_status,
                "fe_status": seed.fe_status,
                "candidate_keep_for_qm": seed.candidate_keep_for_qm,
                "seed_mutation_canonical_positions": list(seed.seed_mutation_canonical_positions),
                "seed_mutation_am1_positions": list(seed.seed_mutation_am1_positions),
                "redesignable_canonical_positions": [
                    int(position)
                    for position in seed_rows_by_candidate[seed.candidate_id].redesignable_canonical_positions.split(",")
                    if position
                ],
                "redesignable_am1_positions": [
                    int(position)
                    for position in seed_rows_by_candidate[seed.candidate_id].redesignable_am1_positions.split(",")
                    if position
                ],
                "exported_backbone_path": seed_rows_by_candidate[seed.candidate_id].exported_backbone_path,
                "exported_pdb_path": seed_rows_by_candidate[seed.candidate_id].exported_pdb_path,
                "chain_assignment_path": seed_rows_by_candidate[seed.candidate_id].chain_assignment_path,
                "fixed_positions_path": seed_rows_by_candidate[seed.candidate_id].fixed_positions_path,
            }
            for seed in seed_contexts
        },
    }


def build_redesign_round2_artifacts(
    *,
    openmm_metadynamics_summary_path: Path,
    openmm_metadynamics_panel_status_path: Path,
    design_mask_candidates_path: Path,
    design_masks_path: Path,
    ligandmpnn_shortlist_path: Path,
    md_validation_panel_path: Path,
    design_campaign_manifest_path: Path,
    design_campaign_positions_path: Path,
    proteinmpnn_round2_root: Path,
    seed_candidate_ids: Sequence[str] = ROUND2_SEED_IDS,
) -> Round2Artifacts:
    """Build deterministic Phase 7A round-2 redesign planning artifacts."""
    design_mask_rows = load_design_mask_candidates(design_mask_candidates_path)
    design_masks_config = _load_yaml_mapping(design_masks_path)
    interface_window_size = validate_design_mask_config(
        design_mask_rows=design_mask_rows,
        design_masks_config=design_masks_config,
    )
    panel_status_rows = load_openmm_metadynamics_panel_status_rows(openmm_metadynamics_panel_status_path)
    summary_rows = load_openmm_metadynamics_summary_rows(openmm_metadynamics_summary_path)
    ligand_shortlist_rows = load_ligandmpnn_shortlist_rows(ligandmpnn_shortlist_path)
    md_panel_rows = load_md_validation_panel_rows(md_validation_panel_path)
    campaign_manifest_rows = load_design_campaign_manifest_rows(design_campaign_manifest_path)
    design_campaign_position_rows = load_design_campaign_position_rows(design_campaign_positions_path)

    seed_contexts = discover_round2_seed_contexts(
        seed_candidate_ids=seed_candidate_ids,
        panel_status_rows=panel_status_rows,
        summary_rows=summary_rows,
        ligand_shortlist_rows=ligand_shortlist_rows,
        md_panel_rows=md_panel_rows,
        campaign_manifest_rows=campaign_manifest_rows,
        design_campaign_position_rows=design_campaign_position_rows,
    )
    position_rows = tuple(
        row
        for seed_context in seed_contexts
        for row in select_round2_redesign_positions(
            seed_context=seed_context,
            design_mask_rows=design_mask_rows,
            interface_window_size=interface_window_size,
        )
    )
    seed_rows = _build_seed_rows(
        seed_contexts=seed_contexts,
        position_rows=position_rows,
        proteinmpnn_round2_root=proteinmpnn_round2_root,
    )
    source_artifacts = {
        "openmm_metadynamics_summary": _display_path(openmm_metadynamics_summary_path),
        "openmm_metadynamics_panel_status": _display_path(openmm_metadynamics_panel_status_path),
        "design_mask_candidates": _display_path(design_mask_candidates_path),
        "design_masks": _display_path(design_masks_path),
        "ligandmpnn_shortlist": _display_path(ligandmpnn_shortlist_path),
        "md_validation_panel": _display_path(md_validation_panel_path),
        "design_campaign_manifest": _display_path(design_campaign_manifest_path),
        "design_campaign_positions": _display_path(design_campaign_positions_path),
    }
    config_payload = _build_round2_config_payload(
        seed_contexts=seed_contexts,
        seed_rows=seed_rows,
        interface_window_size=interface_window_size,
        source_artifacts=source_artifacts,
    )
    return Round2Artifacts(
        interface_window_size=interface_window_size,
        seed_contexts=seed_contexts,
        seed_rows=seed_rows,
        position_rows=position_rows,
        config_payload=config_payload,
        report_markdown=_render_round2_report(
            seed_contexts=seed_contexts,
            seed_rows=seed_rows,
            position_rows=position_rows,
            interface_window_size=interface_window_size,
        ),
    )


def write_redesign_round2_outputs(
    *,
    artifacts: Round2Artifacts,
    config_path: Path,
    seed_table_path: Path,
    position_table_path: Path,
    report_path: Path,
    proteinmpnn_round2_root: Path,
) -> None:
    """Write Phase 7A round-2 redesign planning outputs to disk."""
    write_csv_rows(seed_table_path, artifacts.seed_rows)
    write_csv_rows(position_table_path, artifacts.position_rows)
    write_yaml(config_path, artifacts.config_payload)
    atomic_write_text(report_path, artifacts.report_markdown)

    position_rows_by_candidate: dict[str, list[Round2PositionRow]] = defaultdict(list)
    for row in artifacts.position_rows:
        position_rows_by_candidate[row.candidate_id].append(row)

    for seed_context in artifacts.seed_contexts:
        export_paths = _build_seed_export_paths(
            seed_context=seed_context,
            proteinmpnn_round2_root=proteinmpnn_round2_root,
        )
        export_paths.campaign_dir.mkdir(parents=True, exist_ok=True)
        copy_if_changed(seed_context.starting_structure_path, export_paths.backbone_path)
        copy_if_changed(seed_context.starting_structure_path, export_paths.exported_pdb_path)
        write_jsonl(
            export_paths.chain_assignment_path,
            [
                {
                    seed_context.candidate_id: [
                        list(seed_context.designed_chains),
                        list(seed_context.fixed_context_chains),
                    ]
                }
            ],
        )
        selected_positions = {
            (row.chain_id, row.sequence_index)
            for row in position_rows_by_candidate[seed_context.candidate_id]
        }
        fixed_positions_payload = {
            seed_context.candidate_id: {
                chain_id: [
                    position.sequence_index
                    for position in seed_context.resolved_positions
                    if position.chain_id == chain_id and (position.chain_id, position.sequence_index) not in selected_positions
                ]
                for chain_id in seed_context.designed_chains
            }
        }
        write_jsonl(export_paths.fixed_positions_path, [fixed_positions_payload])

    required_paths = (
        config_path,
        seed_table_path,
        position_table_path,
        report_path,
    )
    for required_path in required_paths:
        if not required_path.exists():
            raise FileNotFoundError(f"Expected Phase 7A output was not written: {required_path}")
    for seed_context in artifacts.seed_contexts:
        export_paths = _build_seed_export_paths(
            seed_context=seed_context,
            proteinmpnn_round2_root=proteinmpnn_round2_root,
        )
        for required_path in (
            export_paths.backbone_path,
            export_paths.exported_pdb_path,
            export_paths.chain_assignment_path,
            export_paths.fixed_positions_path,
        ):
            if not required_path.exists():
                raise FileNotFoundError(f"Expected Phase 7A campaign export was not written: {required_path}")
