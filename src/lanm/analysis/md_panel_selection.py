"""Phase 5B integrated ranking and deterministic MD validation panel selection."""

from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

from lanm.configuration import ProjectConfig, load_project_config
from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_text, write_csv_rows, write_yaml
from lanm.models import (
    LigandMPNNSmokeSummaryRow,
    LigandMPNNShortlistRow,
    ProteinMPNNShortlistRow,
)
from lanm.paths import REPO_ROOT

PANEL_LIMIT = 6
LIGAND_CONFIDENCE_WEIGHT = 0.50
ROSETTA_WITHIN_BACKBONE_WEIGHT = 0.35
PROTEINMPNN_RANK_WEIGHT = 0.15
ROSETTA_COMPARISON_NOTE = (
    "Raw Rosetta total_score is not a valid cross-topology comparator; only within-backbone/topology "
    "rank and normalized score are used across AM1 monomer, Hans monomer, and Hans multichain assemblies."
)
REFERENCE_STARTING_STRUCTURES = {
    "am1_mex_8fns_chain_a": {
        "reference_id": "am1_mex_wt_reference",
        "campaign_id": "am1_mex_wild_type",
        "panel_role": "am1_mex_wild_type_reference",
        "reference_label": "AM1/Mex wild-type reference",
        "starting_structure_path": "results/design_inputs/proteinmpnn/backbones/am1_mex_8fns_chain_a.pdb",
    },
    "hans_pocket_8fnr_chain_a": {
        "reference_id": "hans_pocket_wt_reference",
        "campaign_id": "hans_wild_type",
        "panel_role": "hans_wild_type_reference",
        "reference_label": "Hans wild-type reference",
        "starting_structure_path": "results/design_inputs/proteinmpnn/backbones/hans_pocket_8fnr_chain_a.pdb",
    },
    "hans_interface_8fnr_a_b_c_d": {
        "reference_id": "hans_interface_wt_reference",
        "campaign_id": "hans_wild_type",
        "panel_role": "hans_wild_type_reference",
        "reference_label": "Hans wild-type reference",
        "starting_structure_path": "results/design_inputs/proteinmpnn/backbones/hans_interface_8fnr_a_b_c_d.pdb",
    },
}
PANEL_ROLE_ORDER = {
    "am1_mex_designed_candidate": 1,
    "am1_mex_wild_type_reference": 2,
    "hans_pocket_designed_candidate": 3,
    "hans_interface_aware_candidate": 4,
    "hans_wild_type_reference": 5,
    "best_remaining_designed_candidate": 6,
}


@dataclass(frozen=True, slots=True)
class RosettaScoreRankingRow:
    rosetta_rank: int
    shortlist_rank: int
    candidate_id: str
    campaign_id: str
    backbone_id: str
    preserved_metal_identity: str
    representative_design_id: int
    representative_ligand_confidence: float
    representative_overall_confidence: float
    representative_packed_pdb: str
    scorefile_path: str
    total_score: float


@dataclass(frozen=True, slots=True)
class LigandCandidateContext:
    shortlist_rank: int
    candidate_id: str
    campaign_id: str
    backbone_id: str
    preserved_metal_identity: str
    representative_design_id: int
    representative_ligand_confidence: float
    representative_overall_confidence: float
    representative_packed_pdb: str
    source: str


@dataclass(frozen=True, slots=True)
class IntegratedCandidateRow:
    integrated_rank: int
    candidate_id: str
    campaign_id: str
    backbone_id: str
    topology_class: str
    preserved_metal_identity: str
    proteinmpnn_rank: int
    ligandmpnn_shortlist_rank: int
    ligandmpnn_representative_design_id: int
    ligandmpnn_representative_ligand_confidence: float
    ligandmpnn_representative_overall_confidence: float
    ligandmpnn_mean_ligand_confidence: float | None
    ligandmpnn_max_ligand_confidence: float | None
    ligandmpnn_mean_overall_confidence: float | None
    ligandmpnn_max_overall_confidence: float | None
    rosetta_global_rank: int
    rosetta_score_rank_within_topology: int
    rosetta_total_score: float
    rosetta_total_score_delta_from_backbone_best: float
    rosetta_total_score_normalized_within_backbone: float
    integrated_priority_score: float
    direct_cross_topology_comparison_valid: bool
    cross_topology_comparison_note: str
    representative_packed_pdb: str
    ligand_context_source: str


@dataclass(frozen=True, slots=True)
class MDPanelRow:
    panel_rank: int
    panel_member_id: str
    panel_member_type: str
    panel_role: str
    candidate_id: str
    campaign_id: str
    backbone_id: str
    topology_class: str
    preserved_metal_identity: str
    proteinmpnn_rank: int | None
    ligandmpnn_representative_ligand_confidence: float | None
    ligandmpnn_representative_overall_confidence: float | None
    ligandmpnn_mean_ligand_confidence: float | None
    rosetta_score_rank_within_topology: int | None
    rosetta_total_score: float | None
    starting_structure_path: str
    selection_reason: str


@dataclass(frozen=True, slots=True)
class WildTypeReference:
    panel_member_id: str
    panel_role: str
    candidate_id: str
    campaign_id: str
    backbone_id: str
    topology_class: str
    preserved_metal_identity: str
    starting_structure_path: str
    selection_reason: str


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


def _normalize_rank(rank: int, ranks: Sequence[int]) -> float:
    if not ranks:
        return 0.0
    minimum_rank = min(ranks)
    maximum_rank = max(ranks)
    if minimum_rank == maximum_rank:
        return 0.0
    return (rank - minimum_rank) / (maximum_rank - minimum_rank)


def infer_topology_class(*, campaign_id: str, backbone_id: str) -> str:
    """Map campaign/backbone identifiers onto the Phase 5B topology buckets."""
    normalized_campaign_id = campaign_id.strip().lower()
    normalized_backbone_id = backbone_id.strip().lower()
    if normalized_campaign_id.startswith("am1_mex") or normalized_backbone_id.startswith("am1_mex"):
        return "am1_monomer"
    if normalized_campaign_id.startswith("hans_pocket") or normalized_backbone_id.startswith("hans_pocket"):
        return "hans_monomer"
    if normalized_campaign_id.startswith("hans_interface") or normalized_backbone_id.startswith("hans_interface"):
        return "hans_interface_multichain"
    raise ValueError(f"Could not infer topology class from campaign={campaign_id!r}, backbone={backbone_id!r}")


def classify_panel_role(row: IntegratedCandidateRow) -> str:
    """Classify a designed candidate into the mandatory Phase 5B panel buckets."""
    if row.topology_class == "am1_monomer":
        return "am1_mex_designed_candidate"
    if row.campaign_id == "hans_pocket_ss_only":
        return "hans_pocket_designed_candidate"
    if row.topology_class == "hans_interface_multichain":
        return "hans_interface_aware_candidate"
    if row.topology_class == "hans_monomer":
        return "hans_pocket_designed_candidate"
    raise ValueError(f"Unsupported panel classification for {row.candidate_id}")


def load_proteinmpnn_shortlist_rows(path: Path) -> tuple[ProteinMPNNShortlistRow, ...]:
    """Load the Phase 3C ProteinMPNN shortlist rows."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            ProteinMPNNShortlistRow(
                shortlist_rank=int(str(row["shortlist_rank"]).strip()),
                candidate_id=str(row["candidate_id"]).strip(),
                campaign_id=str(row["campaign_id"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                design_set_name=str(row["design_set_name"]).strip(),
                designed_sequence=str(row["designed_sequence"]).strip(),
                representative_sequence_id=str(row["representative_sequence_id"]).strip(),
                occurrence_count=int(str(row["occurrence_count"]).strip()),
                campaign_rank=int(str(row["campaign_rank"]).strip()),
                temperature=float(str(row["temperature"]).strip()),
                best_score=float(str(row["best_score"]).strip()),
                best_global_score=float(str(row["best_global_score"]).strip()),
                best_seq_recovery=float(str(row["best_seq_recovery"]).strip()),
                mutation_count=int(str(row["mutation_count"]).strip()),
                mutation_string=str(row["mutation_string"]).strip(),
                canonical_family_positions_mutated=str(row["canonical_family_positions_mutated"]).strip(),
                am1_mature_positions_mutated=str(row["am1_mature_positions_mutated"]).strip(),
                includes_second_sphere_position=_parse_bool(str(row["includes_second_sphere_position"])),
                includes_interface_position=_parse_bool(str(row["includes_interface_position"])),
                retention_reason=str(row["retention_reason"]).strip(),
            )
            for row in reader
        ]
    return tuple(sorted(rows, key=lambda row: row.shortlist_rank))


def load_ligandmpnn_shortlist_rows(path: Path) -> tuple[LigandMPNNShortlistRow, ...]:
    """Load the Phase 4C LigandMPNN shortlist rows when available."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            LigandMPNNShortlistRow(
                shortlist_rank=int(str(row["shortlist_rank"]).strip()),
                sequence_id=str(row["sequence_id"]).strip(),
                candidate_id=str(row["candidate_id"]).strip(),
                campaign_id=str(row["campaign_id"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                preserved_metal_identity=str(row["preserved_metal_identity"]).strip(),
                candidate_shortlist_rank=int(str(row["candidate_shortlist_rank"]).strip()),
                candidate_unique_rank=int(str(row["candidate_unique_rank"]).strip()),
                design_id=int(str(row["design_id"]).strip()),
                total_occurrence_count=int(str(row["total_occurrence_count"]).strip()),
                source_candidate_count=int(str(row["source_candidate_count"]).strip()),
                source_candidate_ids=str(row["source_candidate_ids"]).strip(),
                overall_confidence=float(str(row["overall_confidence"]).strip()),
                ligand_confidence=float(str(row["ligand_confidence"]).strip()),
                seq_recovery=float(str(row["seq_recovery"]).strip()),
                mutation_count=int(str(row["mutation_count"]).strip()),
                mutation_string=str(row["mutation_string"]).strip(),
                redesigned_residue_identifiers=str(row["redesigned_residue_identifiers"]).strip(),
                designed_chain_sequence=str(row["designed_chain_sequence"]).strip(),
                input_pdb_path=str(row["input_pdb_path"]).strip(),
                output_fasta_path=str(row["output_fasta_path"]).strip(),
                backbone_pdb_path=str(row["backbone_pdb_path"]).strip(),
                packed_pdb_path=str(row["packed_pdb_path"]).strip(),
                retention_reason=str(row["retention_reason"]).strip(),
            )
            for row in reader
        ]
    return tuple(sorted(rows, key=lambda row: row.shortlist_rank))


def load_ligandmpnn_smoke_summary_rows(path: Path) -> tuple[LigandMPNNSmokeSummaryRow, ...]:
    """Load the Phase 4B LigandMPNN smoke summary rows."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            LigandMPNNSmokeSummaryRow(
                shortlist_rank=int(str(row["shortlist_rank"]).strip()),
                candidate_id=str(row["candidate_id"]).strip(),
                campaign_id=str(row["campaign_id"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                source_structure_id=str(row["source_structure_id"]).strip(),
                designed_chains=str(row["designed_chains"]).strip(),
                preserved_metal_identity=str(row["preserved_metal_identity"]).strip(),
                redesigned_residues=str(row["redesigned_residues"]).strip(),
                generated_sequence_count=int(str(row["generated_sequence_count"]).strip()),
                unique_sequence_count=int(str(row["unique_sequence_count"]).strip()),
                unique_sequence_fraction=float(str(row["unique_sequence_fraction"]).strip()),
                sequence_length=int(str(row["sequence_length"]).strip()),
                mean_overall_confidence=_parse_optional_float(str(row["mean_overall_confidence"])),
                min_overall_confidence=_parse_optional_float(str(row["min_overall_confidence"])),
                max_overall_confidence=_parse_optional_float(str(row["max_overall_confidence"])),
                mean_ligand_confidence=_parse_optional_float(str(row["mean_ligand_confidence"])),
                min_ligand_confidence=_parse_optional_float(str(row["min_ligand_confidence"])),
                max_ligand_confidence=_parse_optional_float(str(row["max_ligand_confidence"])),
                mean_pairwise_identity=_parse_optional_float(str(row["mean_pairwise_identity"])),
                min_pairwise_identity=_parse_optional_float(str(row["min_pairwise_identity"])),
                max_pairwise_identity=_parse_optional_float(str(row["max_pairwise_identity"])),
                temperature=float(str(row["temperature"]).strip()),
                seed=int(str(row["seed"]).strip()),
                output_fasta_path=str(row["output_fasta_path"]).strip(),
            )
            for row in reader
        ]
    return tuple(sorted(rows, key=lambda row: row.shortlist_rank))


def load_rosetta_score_ranking_rows(path: Path) -> tuple[RosettaScoreRankingRow, ...]:
    """Load the Phase 5A1 Rosetta score ranking rows."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            RosettaScoreRankingRow(
                rosetta_rank=int(str(row["rosetta_rank"]).strip()),
                shortlist_rank=int(str(row["shortlist_rank"]).strip()),
                candidate_id=str(row["candidate_id"]).strip(),
                campaign_id=str(row["campaign_id"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                preserved_metal_identity=str(row["preserved_metal_identity"]).strip(),
                representative_design_id=int(str(row["representative_design_id"]).strip()),
                representative_ligand_confidence=float(str(row["representative_ligand_confidence"]).strip()),
                representative_overall_confidence=float(str(row["representative_overall_confidence"]).strip()),
                representative_packed_pdb=str(row["representative_packed_pdb"]).strip(),
                scorefile_path=str(row["scorefile_path"]).strip(),
                total_score=float(str(row["total_score"]).strip()),
            )
            for row in reader
        ]
    return tuple(sorted(rows, key=lambda row: row.rosetta_rank))


def resolve_ligand_candidate_contexts(
    *,
    ligandmpnn_shortlist_path: Path,
    rosetta_rows: Sequence[RosettaScoreRankingRow],
) -> tuple[LigandCandidateContext, ...]:
    """Resolve candidate-level LigandMPNN context, falling back to smoke/Rosetta metadata when needed."""
    shortlist_map: dict[str, LigandMPNNShortlistRow] = {}
    if ligandmpnn_shortlist_path.exists():
        shortlist_map = {row.candidate_id: row for row in load_ligandmpnn_shortlist_rows(ligandmpnn_shortlist_path)}

    contexts: list[LigandCandidateContext] = []
    for rosetta_row in rosetta_rows:
        shortlist_row = shortlist_map.get(rosetta_row.candidate_id)
        if shortlist_row is not None:
            _validate_candidate_identity(
                rosetta_row=rosetta_row,
                campaign_id=shortlist_row.campaign_id,
                backbone_id=shortlist_row.backbone_id,
                preserved_metal_identity=shortlist_row.preserved_metal_identity,
            )
            contexts.append(
                LigandCandidateContext(
                    shortlist_rank=shortlist_row.shortlist_rank,
                    candidate_id=shortlist_row.candidate_id,
                    campaign_id=shortlist_row.campaign_id,
                    backbone_id=shortlist_row.backbone_id,
                    preserved_metal_identity=shortlist_row.preserved_metal_identity,
                    representative_design_id=shortlist_row.design_id,
                    representative_ligand_confidence=shortlist_row.ligand_confidence,
                    representative_overall_confidence=shortlist_row.overall_confidence,
                    representative_packed_pdb=shortlist_row.packed_pdb_path,
                    source="ligandmpnn_shortlist.csv",
                )
            )
            continue

        contexts.append(
            LigandCandidateContext(
                shortlist_rank=rosetta_row.shortlist_rank,
                candidate_id=rosetta_row.candidate_id,
                campaign_id=rosetta_row.campaign_id,
                backbone_id=rosetta_row.backbone_id,
                preserved_metal_identity=rosetta_row.preserved_metal_identity,
                representative_design_id=rosetta_row.representative_design_id,
                representative_ligand_confidence=rosetta_row.representative_ligand_confidence,
                representative_overall_confidence=rosetta_row.representative_overall_confidence,
                representative_packed_pdb=rosetta_row.representative_packed_pdb,
                source="inferred_from_ligandmpnn_smoke_and_rosetta_score",
            )
        )

    return tuple(sorted(contexts, key=lambda row: row.shortlist_rank))


def normalize_rosetta_scores_by_backbone(
    rosetta_rows: Sequence[RosettaScoreRankingRow],
) -> dict[str, dict[str, float | int]]:
    """Normalize raw Rosetta total_score values within each backbone/topology cohort."""
    rows_by_backbone: dict[str, list[RosettaScoreRankingRow]] = defaultdict(list)
    for row in rosetta_rows:
        rows_by_backbone[row.backbone_id].append(row)

    normalized: dict[str, dict[str, float | int]] = {}
    for backbone_id, grouped_rows in rows_by_backbone.items():
        del backbone_id
        ordered_rows = sorted(grouped_rows, key=lambda row: (row.total_score, row.candidate_id))
        minimum_score = ordered_rows[0].total_score
        maximum_score = ordered_rows[-1].total_score
        score_span = maximum_score - minimum_score
        for rank, row in enumerate(ordered_rows, start=1):
            normalized[row.candidate_id] = {
                "rank_within_topology": rank,
                "delta_from_backbone_best": row.total_score - minimum_score,
                "normalized_within_backbone": 0.0 if score_span == 0 else (row.total_score - minimum_score) / score_span,
            }
    return normalized


def build_integrated_candidate_rows(
    *,
    protein_rows: Sequence[ProteinMPNNShortlistRow],
    ligand_context_rows: Sequence[LigandCandidateContext],
    smoke_rows: Sequence[LigandMPNNSmokeSummaryRow],
    rosetta_rows: Sequence[RosettaScoreRankingRow],
) -> tuple[IntegratedCandidateRow, ...]:
    """Merge ProteinMPNN, LigandMPNN, and Rosetta data into one topology-aware ranking table."""
    protein_by_candidate = {row.candidate_id: row for row in protein_rows}
    ligand_by_candidate = {row.candidate_id: row for row in ligand_context_rows}
    smoke_by_candidate = {row.candidate_id: row for row in smoke_rows}
    rosetta_normalization = normalize_rosetta_scores_by_backbone(rosetta_rows)
    protein_ranks = [row.shortlist_rank for row in protein_rows]

    pending_rows: list[IntegratedCandidateRow] = []
    for rosetta_row in rosetta_rows:
        protein_row = protein_by_candidate.get(rosetta_row.candidate_id)
        if protein_row is None:
            raise ValueError(f"Missing ProteinMPNN shortlist row for {rosetta_row.candidate_id}")
        ligand_row = ligand_by_candidate.get(rosetta_row.candidate_id)
        if ligand_row is None:
            raise ValueError(f"Missing LigandMPNN context for {rosetta_row.candidate_id}")
        smoke_row = smoke_by_candidate.get(rosetta_row.candidate_id)
        if smoke_row is None:
            raise ValueError(f"Missing LigandMPNN smoke summary row for {rosetta_row.candidate_id}")

        _validate_candidate_identity(
            rosetta_row=rosetta_row,
            campaign_id=protein_row.campaign_id,
            backbone_id=protein_row.backbone_id,
            preserved_metal_identity=ligand_row.preserved_metal_identity,
        )
        if smoke_row.backbone_id != rosetta_row.backbone_id or smoke_row.campaign_id != rosetta_row.campaign_id:
            raise ValueError(f"Smoke summary mismatch for {rosetta_row.candidate_id}")

        topology_class = infer_topology_class(
            campaign_id=rosetta_row.campaign_id,
            backbone_id=rosetta_row.backbone_id,
        )
        protein_rank_component = _normalize_rank(protein_row.shortlist_rank, protein_ranks)
        rosetta_metadata = rosetta_normalization[rosetta_row.candidate_id]
        rosetta_component = float(rosetta_metadata["normalized_within_backbone"])
        ligand_component = 1.0 - ligand_row.representative_ligand_confidence
        integrated_priority_score = (
            LIGAND_CONFIDENCE_WEIGHT * ligand_component
            + ROSETTA_WITHIN_BACKBONE_WEIGHT * rosetta_component
            + PROTEINMPNN_RANK_WEIGHT * protein_rank_component
        )

        pending_rows.append(
            IntegratedCandidateRow(
                integrated_rank=0,
                candidate_id=rosetta_row.candidate_id,
                campaign_id=rosetta_row.campaign_id,
                backbone_id=rosetta_row.backbone_id,
                topology_class=topology_class,
                preserved_metal_identity=ligand_row.preserved_metal_identity,
                proteinmpnn_rank=protein_row.shortlist_rank,
                ligandmpnn_shortlist_rank=ligand_row.shortlist_rank,
                ligandmpnn_representative_design_id=ligand_row.representative_design_id,
                ligandmpnn_representative_ligand_confidence=ligand_row.representative_ligand_confidence,
                ligandmpnn_representative_overall_confidence=ligand_row.representative_overall_confidence,
                ligandmpnn_mean_ligand_confidence=smoke_row.mean_ligand_confidence,
                ligandmpnn_max_ligand_confidence=smoke_row.max_ligand_confidence,
                ligandmpnn_mean_overall_confidence=smoke_row.mean_overall_confidence,
                ligandmpnn_max_overall_confidence=smoke_row.max_overall_confidence,
                rosetta_global_rank=rosetta_row.rosetta_rank,
                rosetta_score_rank_within_topology=int(rosetta_metadata["rank_within_topology"]),
                rosetta_total_score=rosetta_row.total_score,
                rosetta_total_score_delta_from_backbone_best=float(rosetta_metadata["delta_from_backbone_best"]),
                rosetta_total_score_normalized_within_backbone=rosetta_component,
                integrated_priority_score=integrated_priority_score,
                direct_cross_topology_comparison_valid=False,
                cross_topology_comparison_note=ROSETTA_COMPARISON_NOTE,
                representative_packed_pdb=ligand_row.representative_packed_pdb,
                ligand_context_source=ligand_row.source,
            )
        )

    ordered_rows = sorted(
        pending_rows,
        key=lambda row: (
            row.integrated_priority_score,
            row.rosetta_total_score_normalized_within_backbone,
            -row.ligandmpnn_representative_ligand_confidence,
            row.proteinmpnn_rank,
            row.candidate_id,
        ),
    )
    return tuple(
        IntegratedCandidateRow(
            integrated_rank=index,
            candidate_id=row.candidate_id,
            campaign_id=row.campaign_id,
            backbone_id=row.backbone_id,
            topology_class=row.topology_class,
            preserved_metal_identity=row.preserved_metal_identity,
            proteinmpnn_rank=row.proteinmpnn_rank,
            ligandmpnn_shortlist_rank=row.ligandmpnn_shortlist_rank,
            ligandmpnn_representative_design_id=row.ligandmpnn_representative_design_id,
            ligandmpnn_representative_ligand_confidence=row.ligandmpnn_representative_ligand_confidence,
            ligandmpnn_representative_overall_confidence=row.ligandmpnn_representative_overall_confidence,
            ligandmpnn_mean_ligand_confidence=row.ligandmpnn_mean_ligand_confidence,
            ligandmpnn_max_ligand_confidence=row.ligandmpnn_max_ligand_confidence,
            ligandmpnn_mean_overall_confidence=row.ligandmpnn_mean_overall_confidence,
            ligandmpnn_max_overall_confidence=row.ligandmpnn_max_overall_confidence,
            rosetta_global_rank=row.rosetta_global_rank,
            rosetta_score_rank_within_topology=row.rosetta_score_rank_within_topology,
            rosetta_total_score=row.rosetta_total_score,
            rosetta_total_score_delta_from_backbone_best=row.rosetta_total_score_delta_from_backbone_best,
            rosetta_total_score_normalized_within_backbone=row.rosetta_total_score_normalized_within_backbone,
            integrated_priority_score=row.integrated_priority_score,
            direct_cross_topology_comparison_valid=row.direct_cross_topology_comparison_valid,
            cross_topology_comparison_note=row.cross_topology_comparison_note,
            representative_packed_pdb=row.representative_packed_pdb,
            ligand_context_source=row.ligand_context_source,
        )
        for index, row in enumerate(ordered_rows, start=1)
    )


def select_md_validation_panel(
    integrated_rows: Sequence[IntegratedCandidateRow],
    *,
    limit: int = PANEL_LIMIT,
) -> tuple[MDPanelRow, ...]:
    """Select the deterministic Phase 5B MD validation panel with required scaffold coverage."""
    if limit <= 0:
        raise ValueError("MD panel limit must be positive")

    selected_rows: list[MDPanelRow] = []
    selected_candidate_ids: set[str] = set()

    role_predicates = (
        ("am1_mex_designed_candidate", lambda row: row.topology_class == "am1_monomer"),
        ("hans_pocket_designed_candidate", lambda row: row.campaign_id == "hans_pocket_ss_only"),
        ("hans_interface_aware_candidate", lambda row: row.topology_class == "hans_interface_multichain"),
    )
    for panel_role, predicate in role_predicates:
        candidate = _choose_best_integrated_row(integrated_rows, predicate=predicate, excluded_ids=selected_candidate_ids)
        if candidate is None:
            continue
        selected_candidate_ids.add(candidate.candidate_id)
        selected_rows.append(_panel_row_from_integrated(candidate, panel_role=panel_role))

    for reference in build_wild_type_references(integrated_rows):
        if len(selected_rows) >= limit:
            break
        selected_rows.append(
            MDPanelRow(
                panel_rank=0,
                panel_member_id=reference.panel_member_id,
                panel_member_type="wild_type_reference",
                panel_role=reference.panel_role,
                candidate_id=reference.candidate_id,
                campaign_id=reference.campaign_id,
                backbone_id=reference.backbone_id,
                topology_class=reference.topology_class,
                preserved_metal_identity=reference.preserved_metal_identity,
                proteinmpnn_rank=None,
                ligandmpnn_representative_ligand_confidence=None,
                ligandmpnn_representative_overall_confidence=None,
                ligandmpnn_mean_ligand_confidence=None,
                rosetta_score_rank_within_topology=None,
                rosetta_total_score=None,
                starting_structure_path=reference.starting_structure_path,
                selection_reason=reference.selection_reason,
            )
        )

    remaining_capacity = limit - len(selected_rows)
    if remaining_capacity > 0:
        for candidate in sorted_integrated_rows_for_panel(integrated_rows):
            if candidate.candidate_id in selected_candidate_ids:
                continue
            selected_candidate_ids.add(candidate.candidate_id)
            selected_rows.append(
                _panel_row_from_integrated(
                    candidate,
                    panel_role="best_remaining_designed_candidate",
                )
            )
            if len(selected_rows) >= limit:
                break

    ordered_panel = sorted(
        selected_rows,
        key=lambda row: (
            PANEL_ROLE_ORDER.get(row.panel_role, 999),
            row.candidate_id,
        ),
    )
    return tuple(
        MDPanelRow(
            panel_rank=index,
            panel_member_id=row.panel_member_id,
            panel_member_type=row.panel_member_type,
            panel_role=row.panel_role,
            candidate_id=row.candidate_id,
            campaign_id=row.campaign_id,
            backbone_id=row.backbone_id,
            topology_class=row.topology_class,
            preserved_metal_identity=row.preserved_metal_identity,
            proteinmpnn_rank=row.proteinmpnn_rank,
            ligandmpnn_representative_ligand_confidence=row.ligandmpnn_representative_ligand_confidence,
            ligandmpnn_representative_overall_confidence=row.ligandmpnn_representative_overall_confidence,
            ligandmpnn_mean_ligand_confidence=row.ligandmpnn_mean_ligand_confidence,
            rosetta_score_rank_within_topology=row.rosetta_score_rank_within_topology,
            rosetta_total_score=row.rosetta_total_score,
            starting_structure_path=row.starting_structure_path,
            selection_reason=row.selection_reason,
        )
        for index, row in enumerate(ordered_panel, start=1)
    )


def build_wild_type_references(integrated_rows: Sequence[IntegratedCandidateRow]) -> tuple[WildTypeReference, ...]:
    """Build the required AM1 and Hans wild-type references from the available scaffold set."""
    rows_by_backbone: dict[str, list[IntegratedCandidateRow]] = defaultdict(list)
    for row in integrated_rows:
        rows_by_backbone[row.backbone_id].append(row)

    references: list[WildTypeReference] = []
    am1_reference = _build_reference_for_backbone("am1_mex_8fns_chain_a", rows_by_backbone)
    if am1_reference is not None:
        references.append(am1_reference)

    hans_backbone = (
        "hans_interface_8fnr_a_b_c_d"
        if "hans_interface_8fnr_a_b_c_d" in rows_by_backbone
        else "hans_pocket_8fnr_chain_a"
    )
    hans_reference = _build_reference_for_backbone(hans_backbone, rows_by_backbone)
    if hans_reference is not None:
        references.append(hans_reference)
    return tuple(references)


def sorted_integrated_rows_for_panel(
    integrated_rows: Sequence[IntegratedCandidateRow],
) -> tuple[IntegratedCandidateRow, ...]:
    """Sort integrated rows by the same deterministic priority used in the integrated table."""
    return tuple(
        sorted(
            integrated_rows,
            key=lambda row: (
                row.integrated_rank,
                row.candidate_id,
            ),
        )
    )


def build_md_panel_config(
    *,
    panel_rows: Sequence[MDPanelRow],
    integrated_rows: Sequence[IntegratedCandidateRow],
    proteinmpnn_shortlist_path: Path,
    ligandmpnn_shortlist_path: Path,
    ligandmpnn_smoke_summary_path: Path,
    rosetta_score_ranking_path: Path,
    integrated_ranking_path: Path,
) -> dict[str, Any]:
    """Build the YAML payload for the selected MD validation panel."""
    return {
        "version": 1,
        "phase": "5B",
        "source_artifacts": {
            "proteinmpnn_shortlist": _display_path(proteinmpnn_shortlist_path),
            "ligandmpnn_shortlist": _display_path(ligandmpnn_shortlist_path),
            "ligandmpnn_smoke_summary": _display_path(ligandmpnn_smoke_summary_path),
            "rosetta_score_candidate_ranking": _display_path(rosetta_score_ranking_path),
            "integrated_candidate_ranking": _display_path(integrated_ranking_path),
        },
        "selection_rules": {
            "panel_limit": PANEL_LIMIT,
            "rosetta_cross_topology_policy": ROSETTA_COMPARISON_NOTE,
            "required_designed_coverage": {
                "am1_mex_designed_candidate": 1,
                "hans_pocket_designed_candidate": 1,
                "hans_interface_aware_candidate": 1,
            },
            "required_references": [
                "am1_mex_wild_type_reference",
                "hans_wild_type_reference",
            ],
        },
        "integrated_candidates": [
            {
                "integrated_rank": row.integrated_rank,
                "candidate_id": row.candidate_id,
                "campaign_id": row.campaign_id,
                "topology_class": row.topology_class,
                "integrated_priority_score": round(row.integrated_priority_score, 6),
                "ligandmpnn_representative_ligand_confidence": row.ligandmpnn_representative_ligand_confidence,
                "rosetta_score_rank_within_topology": row.rosetta_score_rank_within_topology,
            }
            for row in integrated_rows
        ],
        "panel": [
            {
                "panel_rank": row.panel_rank,
                "panel_member_id": row.panel_member_id,
                "panel_member_type": row.panel_member_type,
                "panel_role": row.panel_role,
                "candidate_id": row.candidate_id,
                "campaign_id": row.campaign_id,
                "backbone_id": row.backbone_id,
                "topology_class": row.topology_class,
                "preserved_metal_identity": row.preserved_metal_identity,
                "proteinmpnn_rank": row.proteinmpnn_rank,
                "ligandmpnn_representative_ligand_confidence": row.ligandmpnn_representative_ligand_confidence,
                "rosetta_score_rank_within_topology": row.rosetta_score_rank_within_topology,
                "rosetta_total_score": row.rosetta_total_score,
                "starting_structure_path": row.starting_structure_path,
                "selection_reason": row.selection_reason,
            }
            for row in panel_rows
        ],
    }


def render_md_panel_selection_markdown(
    *,
    project_config: ProjectConfig,
    integrated_rows: Sequence[IntegratedCandidateRow],
    smoke_rows: Sequence[LigandMPNNSmokeSummaryRow],
    panel_rows: Sequence[MDPanelRow],
) -> str:
    """Render the Phase 5B markdown report."""
    lines = [
        "# MD Panel Selection",
        "",
        "Phase 5B integrates ProteinMPNN prioritization, LigandMPNN confidence, and Rosetta score_jd2 screening to nominate a deterministic MD validation panel.",
        "",
        "## Scope",
        "",
        f"- target metal: `{project_config.target_metal}`",
        f"- competitors for later MD validation: `{', '.join(project_config.competitors)}`",
        f"- temperature: `{project_config.temperature_K} K`",
        "- excluded in this phase: `MD`, `QM`, and quantum steps",
        "",
        "## Why raw Rosetta total_score is not used across topologies",
        "",
        f"- {ROSETTA_COMPARISON_NOTE}",
    ]
    for topology_class, lengths in _topology_lengths(smoke_rows).items():
        lines.append(
            f"- `{topology_class}` LigandMPNN sequence lengths in the current set: `{', '.join(str(length) for length in lengths)}` residues."
        )
    lines.extend(
        [
            "- Phase 5B therefore compares Rosetta scores only within the same backbone/topology cohort and uses within-backbone rank plus normalized score in the integrated ranking.",
            "",
            "## Integrated ranking",
            "",
            "| integrated_rank | candidate_id | topology_class | proteinmpnn_rank | ligand_confidence | rosetta_rank_within_topology | rosetta_total_score |",
            "| ---: | --- | --- | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in integrated_rows:
        lines.append(
            "| "
            f"{row.integrated_rank} | "
            f"{row.candidate_id} | "
            f"{row.topology_class} | "
            f"{row.proteinmpnn_rank} | "
            f"{row.ligandmpnn_representative_ligand_confidence:.4f} | "
            f"{row.rosetta_score_rank_within_topology} | "
            f"{row.rosetta_total_score:.3f} |"
        )

    lines.extend(
        [
            "",
            "## MD panel",
            "",
            "The selected panel keeps scaffold coverage, includes explicit wild-type baselines, and preserves candidates that are best-supported by LigandMPNN confidence plus within-topology Rosetta standing.",
            "",
        ]
    )
    for row in panel_rows:
        lines.append(
            f"- `{row.panel_member_id}` ({row.panel_role}): {row.selection_reason}"
        )

    lines.extend(
        [
            "",
            "## Dy-vs-Nd/Y/Al/Fe objective support",
            "",
            f"- The panel spans both AM1/Mex and Hans scaffolds, so later MD can test whether `{project_config.target_metal}` preference emerges from scaffold identity or from specific redesigned residues.",
            "- The Hans pocket-focused and interface-aware designs decouple local pocket tuning from multichain/interface effects, which is important for later Dy-versus-Nd/Y discrimination and for checking whether interface stabilization helps or hurts selectivity.",
            "- The two wild-type references provide baseline behavior for off-target Al/Fe coordination and for Dy-versus-Nd/Y comparison before mutations are introduced.",
            "- Because this phase keeps one best remaining designed candidate after satisfying coverage constraints, the panel still allocates one slot to the strongest extra design signal rather than stopping at categorical minimums.",
            "",
        ]
    )
    return "\n".join(lines)


def run_md_panel_selection(
    *,
    proteinmpnn_shortlist_path: Path,
    ligandmpnn_shortlist_path: Path,
    ligandmpnn_smoke_summary_path: Path,
    rosetta_score_ranking_path: Path,
    integrated_ranking_path: Path,
    panel_path: Path,
    report_path: Path,
    md_panel_config_path: Path,
) -> tuple[tuple[IntegratedCandidateRow, ...], tuple[MDPanelRow, ...]]:
    """Run the complete deterministic Phase 5B ranking and panel selection workflow."""
    protein_rows = load_proteinmpnn_shortlist_rows(proteinmpnn_shortlist_path)
    smoke_rows = load_ligandmpnn_smoke_summary_rows(ligandmpnn_smoke_summary_path)
    rosetta_rows = load_rosetta_score_ranking_rows(rosetta_score_ranking_path)
    ligand_context_rows = resolve_ligand_candidate_contexts(
        ligandmpnn_shortlist_path=ligandmpnn_shortlist_path,
        rosetta_rows=rosetta_rows,
    )
    integrated_rows = build_integrated_candidate_rows(
        protein_rows=protein_rows,
        ligand_context_rows=ligand_context_rows,
        smoke_rows=smoke_rows,
        rosetta_rows=rosetta_rows,
    )
    panel_rows = select_md_validation_panel(integrated_rows)

    write_csv_rows(integrated_ranking_path, integrated_rows)
    write_csv_rows(panel_path, panel_rows)
    atomic_write_text(
        report_path,
        render_md_panel_selection_markdown(
            project_config=load_project_config(),
            integrated_rows=integrated_rows,
            smoke_rows=smoke_rows,
            panel_rows=panel_rows,
        ),
    )
    write_yaml(
        md_panel_config_path,
        build_md_panel_config(
            panel_rows=panel_rows,
            integrated_rows=integrated_rows,
            proteinmpnn_shortlist_path=proteinmpnn_shortlist_path,
            ligandmpnn_shortlist_path=ligandmpnn_shortlist_path,
            ligandmpnn_smoke_summary_path=ligandmpnn_smoke_summary_path,
            rosetta_score_ranking_path=rosetta_score_ranking_path,
            integrated_ranking_path=integrated_ranking_path,
        ),
    )
    return integrated_rows, panel_rows


def _validate_candidate_identity(
    *,
    rosetta_row: RosettaScoreRankingRow,
    campaign_id: str,
    backbone_id: str,
    preserved_metal_identity: str,
) -> None:
    if campaign_id != rosetta_row.campaign_id:
        raise ValueError(f"Campaign mismatch for {rosetta_row.candidate_id}: {campaign_id} vs {rosetta_row.campaign_id}")
    if backbone_id != rosetta_row.backbone_id:
        raise ValueError(f"Backbone mismatch for {rosetta_row.candidate_id}: {backbone_id} vs {rosetta_row.backbone_id}")
    if preserved_metal_identity != rosetta_row.preserved_metal_identity:
        raise ValueError(
            f"Preserved metal mismatch for {rosetta_row.candidate_id}: "
            f"{preserved_metal_identity} vs {rosetta_row.preserved_metal_identity}"
        )


def _choose_best_integrated_row(
    integrated_rows: Sequence[IntegratedCandidateRow],
    *,
    predicate: Any,
    excluded_ids: set[str],
) -> IntegratedCandidateRow | None:
    for row in sorted_integrated_rows_for_panel(integrated_rows):
        if row.candidate_id in excluded_ids:
            continue
        if predicate(row):
            return row
    return None


def _panel_row_from_integrated(row: IntegratedCandidateRow, *, panel_role: str) -> MDPanelRow:
    if panel_role == "am1_mex_designed_candidate":
        selection_reason = (
            f"Retained as the best AM1/Mex design by integrated Phase 5B ranking; "
            f"LigandMPNN ligand_confidence={row.ligandmpnn_representative_ligand_confidence:.4f} and "
            f"Rosetta within-topology rank={row.rosetta_score_rank_within_topology}."
        )
    elif panel_role == "hans_pocket_designed_candidate":
        selection_reason = (
            f"Retained as the best Hans pocket-focused design; "
            f"LigandMPNN ligand_confidence={row.ligandmpnn_representative_ligand_confidence:.4f} and "
            f"Rosetta within-topology rank={row.rosetta_score_rank_within_topology}."
        )
    elif panel_role == "hans_interface_aware_candidate":
        selection_reason = (
            f"Retained as the best Hans interface-aware design to preserve the multichain hypothesis; "
            f"LigandMPNN ligand_confidence={row.ligandmpnn_representative_ligand_confidence:.4f} and "
            f"Rosetta within-topology rank={row.rosetta_score_rank_within_topology}."
        )
    else:
        selection_reason = (
            f"Retained as the strongest remaining designed candidate after mandatory scaffold coverage; "
            f"LigandMPNN ligand_confidence={row.ligandmpnn_representative_ligand_confidence:.4f} and "
            f"Rosetta within-topology rank={row.rosetta_score_rank_within_topology}."
        )
    return MDPanelRow(
        panel_rank=0,
        panel_member_id=row.candidate_id,
        panel_member_type="designed_candidate",
        panel_role=panel_role,
        candidate_id=row.candidate_id,
        campaign_id=row.campaign_id,
        backbone_id=row.backbone_id,
        topology_class=row.topology_class,
        preserved_metal_identity=row.preserved_metal_identity,
        proteinmpnn_rank=row.proteinmpnn_rank,
        ligandmpnn_representative_ligand_confidence=row.ligandmpnn_representative_ligand_confidence,
        ligandmpnn_representative_overall_confidence=row.ligandmpnn_representative_overall_confidence,
        ligandmpnn_mean_ligand_confidence=row.ligandmpnn_mean_ligand_confidence,
        rosetta_score_rank_within_topology=row.rosetta_score_rank_within_topology,
        rosetta_total_score=row.rosetta_total_score,
        starting_structure_path=row.representative_packed_pdb,
        selection_reason=selection_reason,
    )


def _build_reference_for_backbone(
    backbone_id: str,
    rows_by_backbone: dict[str, list[IntegratedCandidateRow]],
) -> WildTypeReference | None:
    metadata = REFERENCE_STARTING_STRUCTURES.get(backbone_id)
    if metadata is None:
        return None
    source_rows = rows_by_backbone.get(backbone_id)
    if not source_rows:
        return None
    topology_class = source_rows[0].topology_class
    preserved_metal_identity = source_rows[0].preserved_metal_identity
    starting_structure_path = metadata["starting_structure_path"]
    if metadata["panel_role"] == "am1_mex_wild_type_reference":
        selection_reason = (
            "Kept as the AM1/Mex wild-type baseline so later MD can measure whether redesigned AM1 variants "
            "improve Dy selectivity relative to the unmodified monomer scaffold."
        )
    else:
        selection_reason = (
            "Kept as the Hans wild-type baseline so later MD can separate mutation-driven effects from the native "
            "Hans scaffold and interface context."
        )
    return WildTypeReference(
        panel_member_id=str(metadata["reference_id"]),
        panel_role=str(metadata["panel_role"]),
        candidate_id=str(metadata["reference_id"]),
        campaign_id=str(metadata["campaign_id"]),
        backbone_id=backbone_id,
        topology_class=topology_class,
        preserved_metal_identity=preserved_metal_identity,
        starting_structure_path=starting_structure_path,
        selection_reason=selection_reason,
    )


def _topology_lengths(
    smoke_rows: Sequence[LigandMPNNSmokeSummaryRow],
) -> dict[str, tuple[int, ...]]:
    lengths_by_topology: dict[str, set[int]] = defaultdict(set)
    for row in smoke_rows:
        topology_class = infer_topology_class(campaign_id=row.campaign_id, backbone_id=row.backbone_id)
        lengths_by_topology[topology_class].add(row.sequence_length)
    return {
        topology_class: tuple(sorted(lengths))
        for topology_class, lengths in sorted(lengths_by_topology.items())
    }
