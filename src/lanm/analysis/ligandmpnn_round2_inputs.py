"""Phase 7C shortlist generation and LigandMPNN round-2 input preparation."""

from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Iterable, Sequence

from lanm.analysis.proteinmpnn_round2_smoke import (
    ProteinMPNNRound2Campaign,
    discover_round2_proteinmpnn_campaigns,
    summarize_round2_proteinmpnn_output,
)
from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_text, copy_if_changed, write_csv_rows, write_yaml
from lanm.models import (
    LigandMPNNRound2InputManifestRow,
    LigandMPNNRound2RedesignPositionRow,
    ProteinMPNNRound2SequenceCatalogRow,
    ProteinMPNNRound2ShortlistRow,
    ProteinMPNNRound2UniqueSequenceRow,
)
from lanm.paths import REPO_ROOT

ROUND2_PER_SEED_SHORTLIST_LIMIT = 2
ROUND2_GLOBAL_SHORTLIST_LIMIT = 6


@dataclass(frozen=True, slots=True)
class Round2MutationSummary:
    mutation_count: int
    mutation_string: str
    redesigned_residue_ids: str
    redesigned_canonical_positions: str
    redesigned_am1_positions: str


@dataclass(frozen=True, slots=True)
class Round2CampaignUniqueSequenceResult:
    campaign: ProteinMPNNRound2Campaign
    raw_sequence_count: int
    unique_rows: tuple[ProteinMPNNRound2UniqueSequenceRow, ...]


@dataclass(frozen=True, slots=True)
class Round2ShortlistSelectionResult:
    shortlist_rows: tuple[ProteinMPNNRound2ShortlistRow, ...]
    skipped_cross_campaign_duplicates: int


@dataclass(frozen=True, slots=True)
class Round2PreparationOutputs:
    unique_rows: tuple[ProteinMPNNRound2UniqueSequenceRow, ...]
    shortlist_rows: tuple[ProteinMPNNRound2ShortlistRow, ...]
    manifest_rows: tuple[LigandMPNNRound2InputManifestRow, ...]
    redesign_rows: tuple[LigandMPNNRound2RedesignPositionRow, ...]


def _display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def load_round2_sequence_catalog_rows(path: Path) -> tuple[ProteinMPNNRound2SequenceCatalogRow, ...]:
    """Load the Phase 7B ProteinMPNN round-2 sequence catalog."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            ProteinMPNNRound2SequenceCatalogRow(
                seed_rank=int(str(row["seed_rank"]).strip()),
                candidate_id=str(row["candidate_id"]).strip(),
                campaign_id=str(row["campaign_id"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                topology_class=str(row["topology_class"]).strip(),
                design_set_name=str(row["design_set_name"]).strip(),
                designed_chains=str(row["designed_chains"]).strip(),
                temperature=float(str(row["temperature"]).strip()),
                sample_number=int(str(row["sample_number"]).strip()),
                score=float(str(row["score"]).strip()),
                global_score=float(str(row["global_score"]).strip()),
                seq_recovery=float(str(row["seq_recovery"]).strip()),
                sequence_id=str(row["sequence_id"]).strip(),
                designed_sequence=str(row["designed_sequence"]).strip(),
                sequence_length=int(str(row["sequence_length"]).strip()),
                mutation_count_vs_seed=int(str(row["mutation_count_vs_seed"]).strip()),
                mutation_string_vs_seed=str(row["mutation_string_vs_seed"]).strip(),
                mutated_residue_ids=str(row["mutated_residue_ids"]).strip(),
                output_fasta_path=str(row["output_fasta_path"]).strip(),
            )
            for row in reader
        ]
    return tuple(
        sorted(
            rows,
            key=lambda row: (
                row.seed_rank,
                row.candidate_id,
                row.temperature,
                row.sample_number,
                row.sequence_id,
            ),
        )
    )


def extract_round2_mutation_summary(
    *,
    campaign: ProteinMPNNRound2Campaign,
    designed_sequence: str,
) -> Round2MutationSummary:
    """Summarize how a round-2 designed sequence differs from its seed scaffold."""
    seed_parts = _split_designed_sequence(campaign.seed_scaffold_sequence, campaign.designed_chains)
    designed_parts = _split_designed_sequence(designed_sequence, campaign.designed_chains)
    redesignable_positions = {
        (position.chain_id, position.sequence_index): position
        for position in campaign.redesignable_positions
    }
    designed_residue_map = {
        (residue.chain_id, residue.sequence_index): residue
        for residue in campaign.designed_residues
    }

    mutation_tokens: list[str] = []
    redesigned_residue_ids: list[str] = []
    redesigned_canonical_positions: list[int] = []
    redesigned_am1_positions: list[int] = []
    for chain_id, seed_part, designed_part in zip(campaign.designed_chains, seed_parts, designed_parts):
        if len(seed_part) != len(designed_part):
            raise ValueError(
                f"Designed sequence length mismatch for {campaign.candidate_id} chain {chain_id}: "
                f"{len(designed_part)} vs seed {len(seed_part)}"
            )
        for sequence_index, (seed_aa, designed_aa) in enumerate(zip(seed_part, designed_part), start=1):
            if seed_aa == designed_aa:
                continue
            residue = designed_residue_map[(chain_id, sequence_index)]
            position = redesignable_positions.get((chain_id, sequence_index))
            if position is None:
                raise ValueError(
                    f"Unexpected mutation outside round-2 redesignable positions for {campaign.candidate_id}: "
                    f"{residue.residue_id}"
                )
            mutation_tokens.append(f"{residue.residue_id}:{seed_aa}>{designed_aa}")
            redesigned_residue_ids.append(residue.residue_id)
            if position.canonical_family_position not in redesigned_canonical_positions:
                redesigned_canonical_positions.append(position.canonical_family_position)
            if position.am1_mature_position not in redesigned_am1_positions:
                redesigned_am1_positions.append(position.am1_mature_position)

    return Round2MutationSummary(
        mutation_count=len(mutation_tokens),
        mutation_string=",".join(mutation_tokens),
        redesigned_residue_ids=",".join(redesigned_residue_ids),
        redesigned_canonical_positions=",".join(str(position) for position in redesigned_canonical_positions),
        redesigned_am1_positions=",".join(str(position) for position in redesigned_am1_positions),
    )


def deduplicate_round2_campaign_sequences(
    *,
    campaign: ProteinMPNNRound2Campaign,
    sequence_catalog_rows: Sequence[ProteinMPNNRound2SequenceCatalogRow],
) -> Round2CampaignUniqueSequenceResult:
    """Collapse duplicate round-2 ProteinMPNN outputs within one seed campaign."""
    if not sequence_catalog_rows:
        raise ValueError(f"No round-2 sequence catalog rows provided for {campaign.candidate_id}")

    grouped_rows: dict[str, list[ProteinMPNNRound2SequenceCatalogRow]] = defaultdict(list)
    for row in sequence_catalog_rows:
        if row.seed_rank != campaign.seed_rank:
            raise ValueError(f"Seed-rank mismatch for {campaign.candidate_id}: {row.seed_rank} vs {campaign.seed_rank}")
        if row.candidate_id != campaign.candidate_id:
            raise ValueError(
                f"Sequence catalog candidate mismatch for {campaign.candidate_id}: {row.candidate_id}"
            )
        if row.campaign_id != campaign.campaign_id:
            raise ValueError(
                f"Sequence catalog campaign mismatch for {campaign.candidate_id}: {row.campaign_id}"
            )
        if row.backbone_id != campaign.backbone_id:
            raise ValueError(
                f"Sequence catalog backbone mismatch for {campaign.candidate_id}: {row.backbone_id}"
            )
        if row.topology_class != campaign.topology_class:
            raise ValueError(
                f"Sequence catalog topology mismatch for {campaign.candidate_id}: {row.topology_class}"
            )
        if row.design_set_name != campaign.design_set_name:
            raise ValueError(
                f"Sequence catalog design_set mismatch for {campaign.candidate_id}: {row.design_set_name}"
            )
        if row.designed_chains != campaign.designed_chains_text:
            raise ValueError(
                f"Sequence catalog designed_chains mismatch for {campaign.candidate_id}: {row.designed_chains}"
            )
        grouped_rows[row.designed_sequence].append(row)

    pending_rows: list[ProteinMPNNRound2UniqueSequenceRow] = []
    for designed_sequence, rows in grouped_rows.items():
        mutation_summary = extract_round2_mutation_summary(
            campaign=campaign,
            designed_sequence=designed_sequence,
        )
        for row in rows:
            if row.mutation_count_vs_seed != mutation_summary.mutation_count:
                raise ValueError(
                    f"Mutation-count mismatch for {campaign.candidate_id} sequence {row.sequence_id}: "
                    f"{row.mutation_count_vs_seed} vs {mutation_summary.mutation_count}"
                )
            if row.mutation_string_vs_seed != mutation_summary.mutation_string:
                raise ValueError(
                    f"Mutation-string mismatch for {campaign.candidate_id} sequence {row.sequence_id}: "
                    f"{row.mutation_string_vs_seed!r} vs {mutation_summary.mutation_string!r}"
                )
            if row.mutated_residue_ids != mutation_summary.redesigned_residue_ids:
                raise ValueError(
                    f"Mutated residue mismatch for {campaign.candidate_id} sequence {row.sequence_id}: "
                    f"{row.mutated_residue_ids!r} vs {mutation_summary.redesigned_residue_ids!r}"
                )

        representative_row = min(rows, key=_catalog_row_rank_key)
        pending_rows.append(
            ProteinMPNNRound2UniqueSequenceRow(
                candidate_id="",
                seed_rank=campaign.seed_rank,
                seed_candidate_id=campaign.candidate_id,
                campaign_id=campaign.campaign_id,
                backbone_id=campaign.backbone_id,
                topology_class=campaign.topology_class,
                design_set_name=campaign.design_set_name,
                designed_chains=campaign.designed_chains_text,
                designed_sequence=designed_sequence,
                representative_sequence_id=representative_row.sequence_id,
                occurrence_count=len(rows),
                campaign_rank=0,
                temperature=representative_row.temperature,
                best_score=min(row.score for row in rows),
                best_global_score=min(row.global_score for row in rows),
                best_seq_recovery=max(row.seq_recovery for row in rows),
                mutation_count=mutation_summary.mutation_count,
                mutation_string=mutation_summary.mutation_string,
                redesigned_residue_ids=mutation_summary.redesigned_residue_ids,
                redesigned_canonical_positions=mutation_summary.redesigned_canonical_positions,
                redesigned_am1_positions=mutation_summary.redesigned_am1_positions,
            )
        )

    ranked_rows = sorted(pending_rows, key=_unique_sequence_rank_key)
    unique_rows = tuple(
        replace(
            row,
            candidate_id=f"{campaign.candidate_id}_r2u{index:02d}",
            campaign_rank=index,
        )
        for index, row in enumerate(ranked_rows, start=1)
    )
    return Round2CampaignUniqueSequenceResult(
        campaign=campaign,
        raw_sequence_count=len(sequence_catalog_rows),
        unique_rows=unique_rows,
    )


def build_round2_shortlist(
    unique_rows: Sequence[ProteinMPNNRound2UniqueSequenceRow],
    *,
    per_seed_limit: int = ROUND2_PER_SEED_SHORTLIST_LIMIT,
    total_limit: int = ROUND2_GLOBAL_SHORTLIST_LIMIT,
) -> Round2ShortlistSelectionResult:
    """Create a balanced shortlist from round-2 unique sequences."""
    pool_by_seed: dict[str, list[ProteinMPNNRound2UniqueSequenceRow]] = defaultdict(list)
    for row in unique_rows:
        pool_by_seed[row.seed_candidate_id].append(row)
    for seed_candidate_id in pool_by_seed:
        pool_by_seed[seed_candidate_id] = sorted(
            pool_by_seed[seed_candidate_id],
            key=lambda row: row.campaign_rank,
        )[:per_seed_limit]

    seed_order = tuple(
        seed_candidate_id
        for seed_candidate_id, _ in sorted(
            (
                (seed_candidate_id, min(row.seed_rank for row in rows))
                for seed_candidate_id, rows in pool_by_seed.items()
            ),
            key=lambda item: (item[1], item[0]),
        )
    )
    selected_rows: list[ProteinMPNNRound2UniqueSequenceRow] = []
    selected_ids: set[str] = set()
    selected_sequences: set[str] = set()
    duplicate_skip_keys: set[tuple[str, str]] = set()
    reason_by_candidate_id: dict[str, str] = {}

    for seed_candidate_id in seed_order:
        next_row = _next_available_round2_row(
            pool_by_seed[seed_candidate_id],
            selected_ids=selected_ids,
            selected_sequences=selected_sequences,
            duplicate_skip_keys=duplicate_skip_keys,
        )
        if next_row is None:
            continue
        selected_rows.append(next_row)
        selected_ids.add(next_row.candidate_id)
        selected_sequences.add(next_row.designed_sequence)
        reason_by_candidate_id[next_row.candidate_id] = "required seed scaffold coverage"

    round_number = 2
    while len(selected_rows) < total_limit:
        round_added = False
        for seed_candidate_id in seed_order:
            next_row = _next_available_round2_row(
                pool_by_seed[seed_candidate_id],
                selected_ids=selected_ids,
                selected_sequences=selected_sequences,
                duplicate_skip_keys=duplicate_skip_keys,
            )
            if next_row is None:
                continue
            selected_rows.append(next_row)
            selected_ids.add(next_row.candidate_id)
            selected_sequences.add(next_row.designed_sequence)
            reason_by_candidate_id[next_row.candidate_id] = f"balanced round-robin fill (round {round_number})"
            round_added = True
            if len(selected_rows) >= total_limit:
                break
        if not round_added:
            break
        round_number += 1

    shortlist_rows = tuple(
        ProteinMPNNRound2ShortlistRow(
            shortlist_rank=index,
            candidate_id=row.candidate_id,
            seed_rank=row.seed_rank,
            seed_candidate_id=row.seed_candidate_id,
            campaign_id=row.campaign_id,
            backbone_id=row.backbone_id,
            topology_class=row.topology_class,
            design_set_name=row.design_set_name,
            designed_chains=row.designed_chains,
            designed_sequence=row.designed_sequence,
            representative_sequence_id=row.representative_sequence_id,
            occurrence_count=row.occurrence_count,
            campaign_rank=row.campaign_rank,
            temperature=row.temperature,
            best_score=row.best_score,
            best_global_score=row.best_global_score,
            best_seq_recovery=row.best_seq_recovery,
            mutation_count=row.mutation_count,
            mutation_string=row.mutation_string,
            redesigned_residue_ids=row.redesigned_residue_ids,
            redesigned_canonical_positions=row.redesigned_canonical_positions,
            redesigned_am1_positions=row.redesigned_am1_positions,
            retention_reason=reason_by_candidate_id[row.candidate_id],
        )
        for index, row in enumerate(selected_rows, start=1)
    )
    return Round2ShortlistSelectionResult(
        shortlist_rows=shortlist_rows,
        skipped_cross_campaign_duplicates=len(duplicate_skip_keys),
    )


def build_round2_redesign_rows(
    *,
    campaign: ProteinMPNNRound2Campaign,
    shortlist_row: ProteinMPNNRound2ShortlistRow,
) -> tuple[LigandMPNNRound2RedesignPositionRow, ...]:
    """Resolve shortlisted mutations onto LigandMPNN residue identifiers."""
    seed_parts = _split_designed_sequence(campaign.seed_scaffold_sequence, campaign.designed_chains)
    designed_parts = _split_designed_sequence(shortlist_row.designed_sequence, campaign.designed_chains)
    redesignable_positions = {
        (position.chain_id, position.sequence_index): position
        for position in campaign.redesignable_positions
    }
    designed_residue_map = {
        (residue.chain_id, residue.sequence_index): residue
        for residue in campaign.designed_residues
    }

    redesign_rows: list[LigandMPNNRound2RedesignPositionRow] = []
    mutation_tokens: list[str] = []
    redesigned_residue_ids: list[str] = []
    for chain_id, seed_part, designed_part in zip(campaign.designed_chains, seed_parts, designed_parts):
        if len(seed_part) != len(designed_part):
            raise ValueError(
                f"Designed sequence length mismatch for {shortlist_row.candidate_id} chain {chain_id}: "
                f"{len(designed_part)} vs seed {len(seed_part)}"
            )
        for sequence_index, (seed_aa, designed_aa) in enumerate(zip(seed_part, designed_part), start=1):
            if seed_aa == designed_aa:
                continue
            residue = designed_residue_map[(chain_id, sequence_index)]
            position = redesignable_positions.get((chain_id, sequence_index))
            if position is None:
                raise ValueError(
                    f"Unexpected shortlist mutation outside redesignable positions for "
                    f"{shortlist_row.candidate_id}: {residue.residue_id}"
                )
            mutation_token = f"{residue.residue_id}:{seed_aa}>{designed_aa}"
            mutation_tokens.append(mutation_token)
            redesigned_residue_ids.append(residue.residue_id)
            redesign_rows.append(
                LigandMPNNRound2RedesignPositionRow(
                    shortlist_rank=shortlist_row.shortlist_rank,
                    candidate_id=shortlist_row.candidate_id,
                    seed_rank=shortlist_row.seed_rank,
                    seed_candidate_id=shortlist_row.seed_candidate_id,
                    campaign_id=shortlist_row.campaign_id,
                    backbone_id=shortlist_row.backbone_id,
                    chain_id=residue.chain_id,
                    sequence_index=residue.sequence_index,
                    residue_seq=residue.residue_seq,
                    insertion_code=residue.insertion_code,
                    ligandmpnn_residue_id=residue.residue_id,
                    residue_name=residue.residue_name,
                    seed_amino_acid=seed_aa,
                    designed_amino_acid=designed_aa,
                    mutation_token=mutation_token,
                    canonical_family_position=position.canonical_family_position,
                    am1_mature_position=position.am1_mature_position,
                )
            )

    if len(redesign_rows) != shortlist_row.mutation_count:
        raise ValueError(
            f"Mutation count mismatch for {shortlist_row.candidate_id}: "
            f"{len(redesign_rows)} vs {shortlist_row.mutation_count}"
        )
    if ",".join(mutation_tokens) != shortlist_row.mutation_string:
        raise ValueError(
            f"Mutation string mismatch for {shortlist_row.candidate_id}: "
            f"{','.join(mutation_tokens)!r} vs {shortlist_row.mutation_string!r}"
        )
    if ",".join(redesigned_residue_ids) != shortlist_row.redesigned_residue_ids:
        raise ValueError(
            f"Redesigned residue mismatch for {shortlist_row.candidate_id}: "
            f"{','.join(redesigned_residue_ids)!r} vs {shortlist_row.redesigned_residue_ids!r}"
        )
    return tuple(redesign_rows)


def render_proteinmpnn_round2_shortlist_markdown(
    *,
    campaign_results: Sequence[Round2CampaignUniqueSequenceResult],
    shortlist_result: Round2ShortlistSelectionResult,
) -> str:
    """Render the Phase 7C round-2 shortlist report."""
    lines = [
        "# ProteinMPNN Round-2 Shortlist",
        "",
        "Phase 7C deterministic deduplication and shortlist selection from Phase 7B round-2 ProteinMPNN smoke outputs.",
        "",
        f"- Per-seed prefilter: top `{ROUND2_PER_SEED_SHORTLIST_LIMIT}` unique sequences by best_score then best_global_score.",
        f"- Global shortlist limit: `{ROUND2_GLOBAL_SHORTLIST_LIMIT}` candidates.",
        "- At least one candidate per seed scaffold is retained when a non-duplicate sequence is available.",
        (
            "- Exact designed-sequence duplicates across seed scaffolds were excluded from the final shortlist: "
            f"`{shortlist_result.skipped_cross_campaign_duplicates}` skipped candidates."
        ),
        "",
        "## Seed Counts",
        "",
        "| seed_candidate_id | campaign_id | raw_sequences | unique_sequences | top2_considered | backbone_id |",
        "| --- | --- | ---: | ---: | ---: | --- |",
    ]
    for result in sorted(campaign_results, key=lambda result: (result.campaign.seed_rank, result.campaign.candidate_id)):
        lines.append(
            "| "
            f"{result.campaign.candidate_id} | "
            f"{result.campaign.campaign_id} | "
            f"{result.raw_sequence_count} | "
            f"{len(result.unique_rows)} | "
            f"{min(ROUND2_PER_SEED_SHORTLIST_LIMIT, len(result.unique_rows))} | "
            f"{result.campaign.backbone_id} |"
        )

    lines.extend(
        [
            "",
            "## Final Shortlist",
            "",
            "| rank | candidate_id | seed_scaffold | campaign_rank | best_score | best_global_score | best_seq_recovery | redesigned_residue_ids | retained_because |",
            "| ---: | --- | --- | ---: | ---: | ---: | ---: | --- | --- |",
        ]
    )
    for row in shortlist_result.shortlist_rows:
        redesigned_residue_ids = row.redesigned_residue_ids or "native"
        lines.append(
            "| "
            f"{row.shortlist_rank} | "
            f"{row.candidate_id} | "
            f"{row.seed_candidate_id} | "
            f"{row.campaign_rank} | "
            f"{row.best_score:.4f} | "
            f"{row.best_global_score:.4f} | "
            f"{row.best_seq_recovery:.4f} | "
            f"{redesigned_residue_ids} | "
            f"{row.retention_reason} |"
        )
    lines.append("")
    return "\n".join(lines)


def export_ligandmpnn_round2_inputs(
    *,
    shortlist_rows: Sequence[ProteinMPNNRound2ShortlistRow],
    campaigns_by_seed: dict[str, ProteinMPNNRound2Campaign],
    manifest_path: Path,
    redesign_positions_path: Path,
    report_path: Path,
    input_root: Path,
    config_path: Path,
) -> tuple[tuple[LigandMPNNRound2InputManifestRow, ...], tuple[LigandMPNNRound2RedesignPositionRow, ...]]:
    """Export LigandMPNN-ready round-2 inputs for shortlisted ProteinMPNN candidates."""
    manifest_rows: list[LigandMPNNRound2InputManifestRow] = []
    redesign_rows: list[LigandMPNNRound2RedesignPositionRow] = []
    shortlist_by_candidate = {
        row.candidate_id: row
        for row in shortlist_rows
    }

    for shortlist_row in sorted(shortlist_rows, key=lambda row: row.shortlist_rank):
        campaign = campaigns_by_seed.get(shortlist_row.seed_candidate_id)
        if campaign is None:
            raise ValueError(f"No round-2 campaign metadata found for {shortlist_row.seed_candidate_id}")
        candidate_redesign_rows = build_round2_redesign_rows(
            campaign=campaign,
            shortlist_row=shortlist_row,
        )
        redesign_rows.extend(candidate_redesign_rows)

        candidate_dir = input_root / shortlist_row.candidate_id
        pdb_path = candidate_dir / f"{shortlist_row.candidate_id}.pdb"
        redesigned_residues_path = candidate_dir / "redesigned_residues.txt"
        copy_if_changed(campaign.pdb_path, pdb_path)
        atomic_write_text(
            redesigned_residues_path,
            _render_redesigned_residues_text(candidate_redesign_rows),
        )
        manifest_rows.append(
            LigandMPNNRound2InputManifestRow(
                shortlist_rank=shortlist_row.shortlist_rank,
                candidate_id=shortlist_row.candidate_id,
                seed_rank=shortlist_row.seed_rank,
                seed_candidate_id=shortlist_row.seed_candidate_id,
                campaign_id=shortlist_row.campaign_id,
                backbone_id=shortlist_row.backbone_id,
                topology_class=shortlist_row.topology_class,
                source_pdb_path=_display_path(campaign.pdb_path),
                designed_chains=campaign.designed_chains_text,
                fixed_context_chains=campaign.fixed_context_chains_text,
                redesigned_residue_count=len(candidate_redesign_rows),
                redesigned_residue_ids=shortlist_row.redesigned_residue_ids,
                pdb_path=_display_path(pdb_path),
                redesigned_residues_path=_display_path(redesigned_residues_path),
            )
        )

    manifest_row_tuple = tuple(sorted(manifest_rows, key=lambda row: row.shortlist_rank))
    redesign_row_tuple = tuple(
        sorted(
            redesign_rows,
            key=lambda row: (
                row.shortlist_rank,
                row.chain_id,
                row.sequence_index,
                row.residue_seq,
                row.insertion_code,
            ),
        )
    )
    write_csv_rows(manifest_path, manifest_row_tuple)
    write_csv_rows(redesign_positions_path, redesign_row_tuple)
    atomic_write_text(
        report_path,
        _render_ligandmpnn_round2_inputs_markdown(
            manifest_rows=manifest_row_tuple,
            shortlist_by_candidate=shortlist_by_candidate,
        ),
    )
    write_yaml(
        config_path,
        _build_ligandmpnn_round2_config_payload(
            manifest_rows=manifest_row_tuple,
            shortlist_by_candidate=shortlist_by_candidate,
        ),
    )
    return manifest_row_tuple, redesign_row_tuple


def prepare_ligandmpnn_round2_inputs(
    *,
    redesign_round2_config_path: Path,
    round2_seed_table_path: Path,
    round2_position_table_path: Path,
    round2_campaigns_dir: Path,
    round2_smoke_output_root: Path,
    round2_summary_path: Path,
    round2_sequence_catalog_path: Path,
    unique_sequences_path: Path,
    shortlist_path: Path,
    shortlist_report_path: Path,
    ligandmpnn_round2_root: Path,
    manifest_path: Path,
    redesign_positions_path: Path,
    inputs_report_path: Path,
    config_path: Path,
) -> Round2PreparationOutputs:
    """Run the deterministic Phase 7C shortlist and input-prep workflow."""
    require_path(round2_summary_path)
    campaigns = discover_round2_proteinmpnn_campaigns(
        redesign_round2_config_path=redesign_round2_config_path,
        round2_seed_table_path=round2_seed_table_path,
        round2_position_table_path=round2_position_table_path,
        campaigns_dir=round2_campaigns_dir,
        output_root=round2_smoke_output_root,
    )
    campaigns_by_seed = {
        campaign.candidate_id: campaign
        for campaign in campaigns
    }
    catalog_rows = load_round2_sequence_catalog_rows(round2_sequence_catalog_path)
    catalog_rows_by_seed: dict[str, list[ProteinMPNNRound2SequenceCatalogRow]] = defaultdict(list)
    for row in catalog_rows:
        catalog_rows_by_seed[row.candidate_id].append(row)

    discovered_seed_ids = {campaign.candidate_id for campaign in campaigns}
    catalog_seed_ids = set(catalog_rows_by_seed)
    if discovered_seed_ids != catalog_seed_ids:
        raise ValueError(
            "Round-2 seed mismatch between redesign artifacts and sequence catalog: "
            f"{sorted(discovered_seed_ids)} vs {sorted(catalog_seed_ids)}"
        )

    campaign_results: list[Round2CampaignUniqueSequenceResult] = []
    for campaign in campaigns:
        fasta_path = campaign.output_dir / "seqs" / f"{campaign.candidate_id}.fa"
        _, generated_rows = summarize_round2_proteinmpnn_output(
            campaign=campaign,
            fasta_path=fasta_path,
        )
        _validate_round2_catalog_rows(
            campaign=campaign,
            loaded_rows=catalog_rows_by_seed[campaign.candidate_id],
            generated_rows=generated_rows,
        )
        campaign_results.append(
            deduplicate_round2_campaign_sequences(
                campaign=campaign,
                sequence_catalog_rows=generated_rows,
            )
        )

    unique_rows = tuple(row for result in campaign_results for row in result.unique_rows)
    shortlist_result = build_round2_shortlist(unique_rows)
    write_csv_rows(unique_sequences_path, unique_rows)
    write_csv_rows(shortlist_path, shortlist_result.shortlist_rows)
    atomic_write_text(
        shortlist_report_path,
        render_proteinmpnn_round2_shortlist_markdown(
            campaign_results=campaign_results,
            shortlist_result=shortlist_result,
        ),
    )

    manifest_rows, redesign_rows = export_ligandmpnn_round2_inputs(
        shortlist_rows=shortlist_result.shortlist_rows,
        campaigns_by_seed=campaigns_by_seed,
        manifest_path=manifest_path,
        redesign_positions_path=redesign_positions_path,
        report_path=inputs_report_path,
        input_root=ligandmpnn_round2_root,
        config_path=config_path,
    )
    return Round2PreparationOutputs(
        unique_rows=unique_rows,
        shortlist_rows=shortlist_result.shortlist_rows,
        manifest_rows=manifest_rows,
        redesign_rows=redesign_rows,
    )


def _render_redesigned_residues_text(
    redesign_rows: Sequence[LigandMPNNRound2RedesignPositionRow],
) -> str:
    if not redesign_rows:
        return ""
    return " ".join(row.ligandmpnn_residue_id for row in redesign_rows) + "\n"


def _render_ligandmpnn_round2_inputs_markdown(
    *,
    manifest_rows: Sequence[LigandMPNNRound2InputManifestRow],
    shortlist_by_candidate: dict[str, ProteinMPNNRound2ShortlistRow],
) -> str:
    lines = [
        "# LigandMPNN Round-2 Inputs",
        "",
        "Phase 7C deterministic LigandMPNN-ready input preparation from the shortlisted round-2 ProteinMPNN candidates.",
        "",
        "- Source PDBs were copied directly from the Phase 7A round-2 campaign exports.",
        "- No LigandMPNN, Rosetta, MD, QM, or quantum execution steps were run in this phase.",
        "",
        "| rank | candidate_id | seed_scaffold | source_pdb | redesigned_residue_ids | mutation_string | exported_pdb |",
        "| ---: | --- | --- | --- | --- | --- | --- |",
    ]
    for row in manifest_rows:
        shortlist_row = shortlist_by_candidate[row.candidate_id]
        redesigned_residue_ids = row.redesigned_residue_ids or "native"
        mutation_string = shortlist_row.mutation_string or "native"
        lines.append(
            "| "
            f"{row.shortlist_rank} | "
            f"{row.candidate_id} | "
            f"{row.seed_candidate_id} | "
            f"{row.source_pdb_path} | "
            f"{redesigned_residue_ids} | "
            f"{mutation_string} | "
            f"{row.pdb_path} |"
        )
    lines.append("")
    return "\n".join(lines)


def _build_ligandmpnn_round2_config_payload(
    *,
    manifest_rows: Sequence[LigandMPNNRound2InputManifestRow],
    shortlist_by_candidate: dict[str, ProteinMPNNRound2ShortlistRow],
) -> dict[str, Any]:
    return {
        "version": 1,
        "phase": "7C",
        "source_artifacts": {
            "redesign_round2_config": "config/redesign_round2.yaml",
            "redesign_round2_seeds": "results/tables/redesign_round2_seeds.csv",
            "redesign_round2_positions": "results/tables/redesign_round2_positions.csv",
            "proteinmpnn_round2_summary": "results/tables/proteinmpnn_round2_summary.csv",
            "proteinmpnn_round2_sequence_catalog": "results/tables/proteinmpnn_round2_sequence_catalog.csv",
            "proteinmpnn_round2_unique_sequences": "results/tables/proteinmpnn_round2_unique_sequences.csv",
            "proteinmpnn_round2_shortlist": "results/tables/proteinmpnn_round2_shortlist.csv",
        },
        "selection_policy": {
            "per_seed_campaign_limit": ROUND2_PER_SEED_SHORTLIST_LIMIT,
            "shortlist_total_limit": ROUND2_GLOBAL_SHORTLIST_LIMIT,
            "include_each_seed_scaffold_if_available": True,
            "avoid_exact_duplicate_sequences_across_campaigns": True,
            "external_execution_started": False,
        },
        "candidates": {
            row.candidate_id: {
                "shortlist_rank": row.shortlist_rank,
                "seed_rank": row.seed_rank,
                "seed_candidate_id": row.seed_candidate_id,
                "campaign_id": row.campaign_id,
                "backbone_id": row.backbone_id,
                "topology_class": row.topology_class,
                "source_pdb_path": row.source_pdb_path,
                "designed_chains": _split_csv_list(row.designed_chains),
                "fixed_context_chains": _split_csv_list(row.fixed_context_chains),
                "redesigned_residue_count": row.redesigned_residue_count,
                "redesigned_residue_ids": _split_csv_list(row.redesigned_residue_ids),
                "best_score": shortlist_by_candidate[row.candidate_id].best_score,
                "best_global_score": shortlist_by_candidate[row.candidate_id].best_global_score,
                "best_seq_recovery": shortlist_by_candidate[row.candidate_id].best_seq_recovery,
                "mutation_count": shortlist_by_candidate[row.candidate_id].mutation_count,
                "mutation_string": shortlist_by_candidate[row.candidate_id].mutation_string,
                "pdb_path": row.pdb_path,
                "redesigned_residues_path": row.redesigned_residues_path,
            }
            for row in manifest_rows
        },
    }


def _validate_round2_catalog_rows(
    *,
    campaign: ProteinMPNNRound2Campaign,
    loaded_rows: Sequence[ProteinMPNNRound2SequenceCatalogRow],
    generated_rows: Sequence[ProteinMPNNRound2SequenceCatalogRow],
) -> None:
    loaded_identity = tuple(_round2_catalog_identity(row) for row in loaded_rows)
    generated_identity = tuple(_round2_catalog_identity(row) for row in generated_rows)
    if loaded_identity != generated_identity:
        raise ValueError(
            f"Round-2 sequence catalog did not match FASTA-derived records for {campaign.candidate_id}"
        )


def _round2_catalog_identity(
    row: ProteinMPNNRound2SequenceCatalogRow,
) -> tuple[object, ...]:
    return (
        row.seed_rank,
        row.candidate_id,
        row.campaign_id,
        row.backbone_id,
        row.topology_class,
        row.design_set_name,
        row.designed_chains,
        row.temperature,
        row.sample_number,
        row.score,
        row.global_score,
        row.seq_recovery,
        row.sequence_id,
        row.designed_sequence,
        row.sequence_length,
        row.mutation_count_vs_seed,
        row.mutation_string_vs_seed,
        row.mutated_residue_ids,
        row.output_fasta_path,
    )


def _catalog_row_rank_key(
    row: ProteinMPNNRound2SequenceCatalogRow,
) -> tuple[float, float, float, float, int, str]:
    return (
        row.score,
        row.global_score,
        -row.seq_recovery,
        row.temperature,
        row.sample_number,
        row.sequence_id,
    )


def _unique_sequence_rank_key(
    row: ProteinMPNNRound2UniqueSequenceRow,
) -> tuple[float, float, float, int, float, str, str]:
    return (
        row.best_score,
        row.best_global_score,
        -row.best_seq_recovery,
        row.mutation_count,
        row.temperature,
        row.representative_sequence_id,
        row.designed_sequence,
    )


def _next_available_round2_row(
    pool: Iterable[ProteinMPNNRound2UniqueSequenceRow],
    *,
    selected_ids: set[str],
    selected_sequences: set[str],
    duplicate_skip_keys: set[tuple[str, str]],
) -> ProteinMPNNRound2UniqueSequenceRow | None:
    for row in pool:
        if row.candidate_id in selected_ids:
            continue
        if row.designed_sequence in selected_sequences:
            duplicate_skip_keys.add((row.seed_candidate_id, row.designed_sequence))
            continue
        return row
    return None


def _split_designed_sequence(sequence: str, designed_chains: Sequence[str]) -> tuple[str, ...]:
    if len(designed_chains) == 1:
        return (sequence.replace("/", ""),)
    parts = tuple(part.strip() for part in sequence.split("/"))
    if len(parts) != len(designed_chains):
        raise ValueError(
            f"Expected {len(designed_chains)} designed-chain sequence parts, found {len(parts)} in {sequence!r}"
        )
    return tuple(part.replace("/", "") for part in parts)


def _split_csv_list(value: str) -> tuple[str, ...]:
    if not value.strip():
        return ()
    return tuple(item.strip() for item in value.split(",") if item.strip())
