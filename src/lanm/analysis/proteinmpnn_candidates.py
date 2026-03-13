"""Phase 3C ProteinMPNN candidate deduplication and shortlist selection."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable, Sequence

from lanm.analysis.proteinmpnn_smoke import (
    ProteinMPNNGeneratedRecord,
    ProteinMPNNParsedOutput,
    ProteinMPNNSmokeCampaign,
    discover_proteinmpnn_smoke_campaigns,
    parse_proteinmpnn_fasta,
)
from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_text, write_csv_rows
from lanm.models import (
    DesignCampaignPositionRow,
    ProteinMPNNShortlistRow,
    ProteinMPNNUniqueSequenceRow,
)

SHORTLIST_CAMPAIGN_LIMIT = 5
SHORTLIST_TOTAL_LIMIT = 12
PRIORITY_CAMPAIGN_ID = "hans_interface_ss_plus_if"
PRIORITY_CAMPAIGN_MINIMUM = 2
DEFAULT_CAMPAIGN_MINIMUM = 1


@dataclass(frozen=True, slots=True)
class MutationSummary:
    mutation_count: int
    mutation_string: str
    canonical_family_positions_mutated: str
    am1_mature_positions_mutated: str
    includes_second_sphere_position: bool
    includes_interface_position: bool


@dataclass(frozen=True, slots=True)
class CampaignUniqueSequenceResult:
    campaign: ProteinMPNNSmokeCampaign
    raw_sequence_count: int
    unique_rows: tuple[ProteinMPNNUniqueSequenceRow, ...]


@dataclass(frozen=True, slots=True)
class ShortlistSelectionResult:
    shortlist_rows: tuple[ProteinMPNNShortlistRow, ...]
    skipped_cross_campaign_duplicates: int


def load_design_campaign_position_rows(path: Path) -> tuple[DesignCampaignPositionRow, ...]:
    """Load the Phase 3A design campaign position table from CSV."""
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            DesignCampaignPositionRow(
                campaign_id=str(row["campaign_id"]).strip(),
                design_set_name=str(row["design_set_name"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                structure_id=str(row["structure_id"]).strip(),
                chain_id=str(row["chain_id"]).strip(),
                sequence_index=int(str(row["sequence_index"]).strip()),
                residue_seq=int(str(row["residue_seq"]).strip()),
                insertion_code=str(row["insertion_code"]).strip(),
                residue_name=str(row["residue_name"]).strip(),
                canonical_family_position=_parse_optional_int(str(row["canonical_family_position"])),
                am1_mature_position=_parse_optional_int(str(row["am1_mature_position"])),
                fixed_first_shell=_parse_bool(str(row["fixed_first_shell"])),
                protected_positions=_parse_bool(str(row["protected_positions"])),
                mutable_second_sphere=_parse_bool(str(row["mutable_second_sphere"])),
                mutable_interface=_parse_bool(str(row["mutable_interface"])),
                hard_fixed=_parse_bool(str(row["hard_fixed"])),
                designable=_parse_bool(str(row["designable"])),
                position_state=str(row["position_state"]).strip(),
                position_reason=str(row["position_reason"]).strip(),
                rationale=str(row["rationale"]).strip(),
            )
            for row in reader
        ]
    return tuple(rows)


def extract_mutation_summary(
    *,
    native_sequence: str,
    designed_sequence: str,
    designed_chains: Sequence[str],
    position_rows: Sequence[DesignCampaignPositionRow],
) -> MutationSummary:
    """Extract a mutation summary versus the native designed-chain sequence."""
    native_parts = _split_designed_sequence(native_sequence, designed_chains)
    designed_parts = _split_designed_sequence(designed_sequence, designed_chains)
    position_rows_by_chain = _group_position_rows_by_chain(position_rows, designed_chains)

    mutation_tokens: list[str] = []
    canonical_positions: list[int] = []
    am1_positions: list[int] = []
    includes_second_sphere = False
    includes_interface = False
    mutation_count = 0

    for chain_id, native_part, designed_part in zip(designed_chains, native_parts, designed_parts):
        chain_rows = position_rows_by_chain[chain_id]
        if len(native_part) != len(chain_rows):
            raise ValueError(
                f"Native designed-chain length {len(native_part)} did not match campaign positions "
                f"{len(chain_rows)} for chain {chain_id}"
            )
        if len(designed_part) != len(chain_rows):
            raise ValueError(
                f"Designed-chain length {len(designed_part)} did not match campaign positions "
                f"{len(chain_rows)} for chain {chain_id}"
            )
        for native_aa, designed_aa, row in zip(native_part, designed_part, chain_rows):
            if native_aa == designed_aa:
                continue
            mutation_count += 1
            mutation_tokens.append(f"{native_aa}{row.sequence_index}{designed_aa}")
            if row.canonical_family_position is not None and row.canonical_family_position not in canonical_positions:
                canonical_positions.append(row.canonical_family_position)
            if row.am1_mature_position is not None and row.am1_mature_position not in am1_positions:
                am1_positions.append(row.am1_mature_position)
            includes_second_sphere = includes_second_sphere or row.mutable_second_sphere
            includes_interface = includes_interface or row.mutable_interface

    return MutationSummary(
        mutation_count=mutation_count,
        mutation_string=",".join(mutation_tokens),
        canonical_family_positions_mutated=",".join(str(position) for position in canonical_positions),
        am1_mature_positions_mutated=",".join(str(position) for position in am1_positions),
        includes_second_sphere_position=includes_second_sphere,
        includes_interface_position=includes_interface,
    )


def deduplicate_campaign_sequences(
    *,
    campaign: ProteinMPNNSmokeCampaign,
    parsed_output: ProteinMPNNParsedOutput,
    position_rows: Sequence[DesignCampaignPositionRow],
) -> CampaignUniqueSequenceResult:
    """Collapse duplicate ProteinMPNN outputs for one campaign into ranked unique sequences."""
    native_record = parsed_output.native_record
    if native_record.name != campaign.campaign_id:
        raise ValueError(
            f"ProteinMPNN FASTA name {native_record.name} did not match campaign {campaign.campaign_id}"
        )
    if native_record.designed_chains != campaign.designed_chains:
        raise ValueError(
            f"Designed chains {native_record.designed_chains} did not match manifest {campaign.designed_chains}"
        )

    grouped_records: dict[str, list[ProteinMPNNGeneratedRecord]] = defaultdict(list)
    for record in parsed_output.generated_records:
        grouped_records[record.sequence].append(record)

    pending_rows: list[ProteinMPNNUniqueSequenceRow] = []
    for designed_sequence, records in grouped_records.items():
        representative_record = min(records, key=_generated_record_rank_key)
        mutation_summary = extract_mutation_summary(
            native_sequence=native_record.sequence,
            designed_sequence=designed_sequence,
            designed_chains=campaign.designed_chains,
            position_rows=position_rows,
        )
        pending_rows.append(
            ProteinMPNNUniqueSequenceRow(
                candidate_id="",
                campaign_id=campaign.campaign_id,
                backbone_id=campaign.backbone_id,
                design_set_name=campaign.design_set_name,
                designed_sequence=designed_sequence,
                representative_sequence_id=_sequence_id(campaign.campaign_id, representative_record),
                occurrence_count=len(records),
                campaign_rank=0,
                temperature=representative_record.temperature,
                best_score=min(record.score for record in records),
                best_global_score=min(record.global_score for record in records),
                best_seq_recovery=max(record.seq_recovery for record in records),
                mutation_count=mutation_summary.mutation_count,
                mutation_string=mutation_summary.mutation_string,
                canonical_family_positions_mutated=mutation_summary.canonical_family_positions_mutated,
                am1_mature_positions_mutated=mutation_summary.am1_mature_positions_mutated,
                includes_second_sphere_position=mutation_summary.includes_second_sphere_position,
                includes_interface_position=mutation_summary.includes_interface_position,
            )
        )

    ranked_rows = sorted(pending_rows, key=_unique_sequence_rank_key)
    unique_rows = tuple(
        replace(
            row,
            candidate_id=f"{campaign.campaign_id}_u{index:02d}",
            campaign_rank=index,
        )
        for index, row in enumerate(ranked_rows, start=1)
    )
    return CampaignUniqueSequenceResult(
        campaign=campaign,
        raw_sequence_count=len(parsed_output.generated_records),
        unique_rows=unique_rows,
    )


def build_balanced_shortlist(
    unique_rows: Sequence[ProteinMPNNUniqueSequenceRow],
    *,
    per_campaign_limit: int = SHORTLIST_CAMPAIGN_LIMIT,
    total_limit: int = SHORTLIST_TOTAL_LIMIT,
    priority_campaign_id: str = PRIORITY_CAMPAIGN_ID,
) -> ShortlistSelectionResult:
    """Build a deterministic, balanced shortlist from ranked unique sequences."""
    pool_by_campaign: dict[str, list[ProteinMPNNUniqueSequenceRow]] = defaultdict(list)
    for row in unique_rows:
        pool_by_campaign[row.campaign_id].append(row)
    for campaign_id in pool_by_campaign:
        pool_by_campaign[campaign_id] = sorted(
            pool_by_campaign[campaign_id],
            key=lambda row: row.campaign_rank,
        )[:per_campaign_limit]

    campaign_order = _shortlist_campaign_order(pool_by_campaign, priority_campaign_id)
    selected_rows: list[ProteinMPNNUniqueSequenceRow] = []
    selected_ids: set[str] = set()
    selected_sequences: set[str] = set()
    reason_by_candidate_id: dict[str, str] = {}
    duplicate_skip_keys: set[tuple[str, str]] = set()

    required_counts = {
        campaign_id: (
            PRIORITY_CAMPAIGN_MINIMUM if campaign_id == priority_campaign_id else DEFAULT_CAMPAIGN_MINIMUM
        )
        for campaign_id in campaign_order
    }

    for campaign_id in campaign_order:
        target_count = required_counts[campaign_id]
        retained_count = 0
        for row in pool_by_campaign[campaign_id]:
            if len(selected_rows) >= total_limit:
                break
            if row.designed_sequence in selected_sequences:
                duplicate_skip_keys.add((campaign_id, row.designed_sequence))
                continue
            selected_rows.append(row)
            selected_ids.add(row.candidate_id)
            selected_sequences.add(row.designed_sequence)
            retained_count += 1
            reason_by_candidate_id[row.candidate_id] = (
                "required interface campaign coverage"
                if campaign_id == priority_campaign_id
                else "required campaign coverage"
            )
            if retained_count >= target_count:
                break

    round_number = 1
    while len(selected_rows) < total_limit:
        round_added = False
        for campaign_id in campaign_order:
            next_row = _next_available_row(
                pool_by_campaign[campaign_id],
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
        ProteinMPNNShortlistRow(
            shortlist_rank=index,
            candidate_id=row.candidate_id,
            campaign_id=row.campaign_id,
            backbone_id=row.backbone_id,
            design_set_name=row.design_set_name,
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
            canonical_family_positions_mutated=row.canonical_family_positions_mutated,
            am1_mature_positions_mutated=row.am1_mature_positions_mutated,
            includes_second_sphere_position=row.includes_second_sphere_position,
            includes_interface_position=row.includes_interface_position,
            retention_reason=reason_by_candidate_id[row.candidate_id],
        )
        for index, row in enumerate(selected_rows, start=1)
    )
    return ShortlistSelectionResult(
        shortlist_rows=shortlist_rows,
        skipped_cross_campaign_duplicates=len(duplicate_skip_keys),
    )


def render_proteinmpnn_shortlist_markdown(
    *,
    campaign_results: Sequence[CampaignUniqueSequenceResult],
    shortlist_result: ShortlistSelectionResult,
) -> str:
    """Render the Phase 3C shortlist report."""
    unique_rows = [row for result in campaign_results for row in result.unique_rows]
    mutation_distribution = Counter(row.mutation_count for row in unique_rows)

    lines = [
        "# ProteinMPNN Candidate Shortlist",
        "",
        "Deterministic Phase 3C deduplication of Phase 3B ProteinMPNN smoke outputs.",
        "",
        f"- Per-campaign prefilter: top `{SHORTLIST_CAMPAIGN_LIMIT}` unique sequences by best score.",
        f"- Global shortlist limit: `{SHORTLIST_TOTAL_LIMIT}` sequences.",
        (
            f"- Priority campaign `{PRIORITY_CAMPAIGN_ID}` contributes at least "
            f"`{PRIORITY_CAMPAIGN_MINIMUM}` sequences when available."
        ),
        (
            "- Exact designed-sequence duplicates across campaigns were excluded from the final shortlist: "
            f"`{shortlist_result.skipped_cross_campaign_duplicates}` skipped candidates."
        ),
        "",
        "## Campaign Counts",
        "",
        "| campaign_id | raw_sequences | unique_sequences | top5_considered | backbone_id | design_set_name |",
        "| --- | ---: | ---: | ---: | --- | --- |",
    ]
    for result in campaign_results:
        lines.append(
            "| "
            f"{result.campaign.campaign_id} | "
            f"{result.raw_sequence_count} | "
            f"{len(result.unique_rows)} | "
            f"{min(SHORTLIST_CAMPAIGN_LIMIT, len(result.unique_rows))} | "
            f"{result.campaign.backbone_id} | "
            f"{result.campaign.design_set_name} |"
        )

    lines.extend(
        [
            "",
            "## Mutation Count Distribution",
            "",
            "| mutation_count | unique_sequences |",
            "| --- | ---: |",
        ]
    )
    for mutation_count in sorted(mutation_distribution):
        lines.append(f"| {mutation_count} | {mutation_distribution[mutation_count]} |")

    lines.extend(
        [
            "",
            "## Final Shortlist",
            "",
            "| rank | candidate_id | campaign_id | best_score | best_global_score | best_seq_recovery | mutations | retained_because |",
            "| ---: | --- | --- | ---: | ---: | ---: | --- | --- |",
        ]
    )
    for row in shortlist_result.shortlist_rows:
        mutation_summary = row.mutation_string or "native"
        lines.append(
            "| "
            f"{row.shortlist_rank} | "
            f"{row.candidate_id} | "
            f"{row.campaign_id} | "
            f"{row.best_score:.4f} | "
            f"{row.best_global_score:.4f} | "
            f"{row.best_seq_recovery:.4f} | "
            f"{mutation_summary} | "
            f"{row.retention_reason} |"
        )
    lines.append("")
    return "\n".join(lines)


def render_proteinmpnn_shortlist_fasta(shortlist_rows: Sequence[ProteinMPNNShortlistRow]) -> str:
    """Render the shortlist FASTA in shortlist rank order."""
    lines: list[str] = []
    for row in shortlist_rows:
        mutation_summary = row.mutation_string or "native"
        lines.append(
            ">"
            f"{row.candidate_id}"
            f"|rank={row.shortlist_rank}"
            f"|campaign={row.campaign_id}"
            f"|score={row.best_score:.4f}"
            f"|global_score={row.best_global_score:.4f}"
            f"|seq_recovery={row.best_seq_recovery:.4f}"
            f"|mutations={mutation_summary}"
        )
        lines.append(row.designed_sequence)
    return "\n".join(lines) + ("\n" if lines else "")


def select_proteinmpnn_candidates(
    *,
    campaign_manifest_path: Path,
    backbone_manifest_path: Path,
    campaigns_dir: Path,
    smoke_output_root: Path,
    design_campaign_positions_path: Path,
    unique_sequences_path: Path,
    shortlist_path: Path,
    report_path: Path,
    shortlist_fasta_path: Path,
) -> tuple[tuple[ProteinMPNNUniqueSequenceRow, ...], tuple[ProteinMPNNShortlistRow, ...]]:
    """Run the deterministic Phase 3C selection workflow and write repo artifacts."""
    campaigns = discover_proteinmpnn_smoke_campaigns(
        campaign_manifest_path=campaign_manifest_path,
        backbone_manifest_path=backbone_manifest_path,
        campaigns_dir=campaigns_dir,
        output_root=smoke_output_root,
    )
    if not campaigns:
        raise ValueError(f"No campaigns discovered under {campaigns_dir}")

    position_rows = load_design_campaign_position_rows(design_campaign_positions_path)
    positions_by_campaign: dict[str, list[DesignCampaignPositionRow]] = defaultdict(list)
    for row in position_rows:
        positions_by_campaign[row.campaign_id].append(row)

    campaign_results: list[CampaignUniqueSequenceResult] = []
    for campaign in campaigns:
        fasta_path = campaign.output_dir / "seqs" / f"{campaign.campaign_id}.fa"
        parsed_output = parse_proteinmpnn_fasta(fasta_path)
        campaign_results.append(
            deduplicate_campaign_sequences(
                campaign=campaign,
                parsed_output=parsed_output,
                position_rows=positions_by_campaign[campaign.campaign_id],
            )
        )

    unique_rows = tuple(row for result in campaign_results for row in result.unique_rows)
    shortlist_result = build_balanced_shortlist(unique_rows)

    write_csv_rows(unique_sequences_path, unique_rows)
    write_csv_rows(shortlist_path, shortlist_result.shortlist_rows)
    atomic_write_text(
        report_path,
        render_proteinmpnn_shortlist_markdown(
            campaign_results=campaign_results,
            shortlist_result=shortlist_result,
        ),
    )
    atomic_write_text(
        shortlist_fasta_path,
        render_proteinmpnn_shortlist_fasta(shortlist_result.shortlist_rows),
    )
    return unique_rows, shortlist_result.shortlist_rows


def _parse_optional_int(value: str) -> int | None:
    stripped = value.strip()
    return int(stripped) if stripped else None


def _parse_bool(value: str) -> bool:
    return value.strip().lower() == "true"


def _split_designed_sequence(sequence: str, designed_chains: Sequence[str]) -> tuple[str, ...]:
    if len(designed_chains) == 1:
        return (sequence.replace("/", ""),)
    parts = tuple(part.strip() for part in sequence.split("/"))
    if len(parts) != len(designed_chains):
        raise ValueError(
            f"Expected {len(designed_chains)} designed-chain sequence parts, found {len(parts)} in {sequence!r}"
        )
    return parts


def _group_position_rows_by_chain(
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


def _generated_record_rank_key(record: ProteinMPNNGeneratedRecord) -> tuple[float, float, float, float, int]:
    return (
        record.score,
        record.global_score,
        -record.seq_recovery,
        record.temperature,
        record.sample_number,
    )


def _unique_sequence_rank_key(
    row: ProteinMPNNUniqueSequenceRow,
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


def _sequence_id(campaign_id: str, record: ProteinMPNNGeneratedRecord) -> str:
    return f"{campaign_id}_T{record.temperature:g}_sample_{record.sample_number:02d}"


def _shortlist_campaign_order(
    pool_by_campaign: dict[str, list[ProteinMPNNUniqueSequenceRow]],
    priority_campaign_id: str,
) -> tuple[str, ...]:
    ordered = [campaign_id for campaign_id in sorted(pool_by_campaign) if campaign_id != priority_campaign_id]
    if priority_campaign_id in pool_by_campaign:
        return (priority_campaign_id, *ordered)
    return tuple(ordered)


def _next_available_row(
    pool: Iterable[ProteinMPNNUniqueSequenceRow],
    *,
    selected_ids: set[str],
    selected_sequences: set[str],
    duplicate_skip_keys: set[tuple[str, str]],
) -> ProteinMPNNUniqueSequenceRow | None:
    for row in pool:
        if row.candidate_id in selected_ids:
            continue
        if row.designed_sequence in selected_sequences:
            duplicate_skip_keys.add((row.campaign_id, row.designed_sequence))
            continue
        return row
    return None
