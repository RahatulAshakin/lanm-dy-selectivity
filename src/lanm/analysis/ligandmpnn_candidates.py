"""Phase 4C LigandMPNN candidate deduplication and shortlist selection."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

from lanm.analysis.ligandmpnn_smoke import parse_ligandmpnn_fasta
from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_text, write_csv_rows, write_yaml
from lanm.models import (
    LigandMPNNInputManifestRow,
    LigandMPNNRedesignPositionRow,
    LigandMPNNShortlistRow,
    LigandMPNNUniqueSequenceRow,
    ProteinMPNNShortlistRow,
)
from lanm.paths import REPO_ROOT

PER_CANDIDATE_SHORTLIST_LIMIT = 2
FINAL_SHORTLIST_LIMIT = 8
INTERFACE_CAMPAIGN_ID = "hans_interface_ss_plus_if"
INTERFACE_CAMPAIGN_MINIMUM = 2
POCKET_CAMPAIGN_ID = "hans_pocket_ss_only"
POCKET_CAMPAIGN_MINIMUM = 1
AM1_CAMPAIGN_PREFIX = "am1_mex"
AM1_CAMPAIGN_MINIMUM = 1
PACKED_SUFFIX = "_packed"


@dataclass(frozen=True, slots=True)
class LigandMPNNSelectionCandidate:
    shortlist_rank: int
    candidate_id: str
    campaign_id: str
    backbone_id: str
    designed_chains: tuple[str, ...]
    preserved_chains: tuple[str, ...]
    preserved_metal_identity: str
    input_sequence: str
    redesign_rows: tuple[LigandMPNNRedesignPositionRow, ...]
    input_pdb_path: Path
    output_dir: Path

    @property
    def output_fasta_path(self) -> Path:
        return self.output_dir / "seqs" / f"{self.candidate_id}.fa"


@dataclass(frozen=True, slots=True)
class LigandMPNNMutationSummary:
    mutation_count: int
    mutation_string: str
    redesigned_residue_identifiers: str


@dataclass(frozen=True, slots=True)
class LigandMPNNGeneratedDesign:
    sequence_id: str
    candidate_id: str
    campaign_id: str
    backbone_id: str
    preserved_metal_identity: str
    candidate_shortlist_rank: int
    design_id: int
    overall_confidence: float
    ligand_confidence: float
    seq_recovery: float
    mutation_count: int
    mutation_string: str
    redesigned_residue_identifiers: str
    designed_chain_sequence: str
    input_pdb_path: str
    output_fasta_path: str
    backbone_pdb_path: str
    packed_pdb_path: str


@dataclass(frozen=True, slots=True)
class CandidateUniqueSequence:
    sequence_id: str
    candidate_id: str
    campaign_id: str
    backbone_id: str
    preserved_metal_identity: str
    candidate_shortlist_rank: int
    candidate_unique_rank: int
    design_id: int
    occurrence_count: int
    overall_confidence: float
    ligand_confidence: float
    seq_recovery: float
    mutation_count: int
    mutation_string: str
    redesigned_residue_identifiers: str
    designed_chain_sequence: str
    input_pdb_path: str
    output_fasta_path: str
    backbone_pdb_path: str
    packed_pdb_path: str


@dataclass(frozen=True, slots=True)
class CandidateUniqueSequenceResult:
    candidate: LigandMPNNSelectionCandidate
    raw_sequence_count: int
    unique_rows: tuple[CandidateUniqueSequence, ...]


@dataclass(frozen=True, slots=True)
class ShortlistSelectionResult:
    shortlist_rows: tuple[LigandMPNNShortlistRow, ...]


def _display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def _repo_path(path_text: str) -> Path:
    candidate = Path(path_text).expanduser()
    if candidate.is_absolute():
        return candidate
    return REPO_ROOT / candidate


def _parse_bool(value: str) -> bool:
    return value.strip().lower() == "true"


def _split_csv_list(value: str) -> tuple[str, ...]:
    return tuple(item.strip() for item in value.split(",") if item.strip())


def _split_designed_sequence(sequence: str, designed_chains: Sequence[str]) -> tuple[str, ...]:
    if len(designed_chains) == 1:
        return (sequence.replace("/", ""),)
    parts = tuple(part.strip() for part in sequence.split("/"))
    if len(parts) != len(designed_chains):
        raise ValueError(
            f"Expected {len(designed_chains)} designed-chain sequence parts, found {len(parts)} in {sequence!r}"
        )
    return parts


def _load_proteinmpnn_shortlist_rows(path: Path) -> tuple[ProteinMPNNShortlistRow, ...]:
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


def _load_ligandmpnn_input_manifest_rows(path: Path) -> tuple[LigandMPNNInputManifestRow, ...]:
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            LigandMPNNInputManifestRow(
                shortlist_rank=int(str(row["shortlist_rank"]).strip()),
                candidate_id=str(row["candidate_id"]).strip(),
                campaign_id=str(row["campaign_id"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                source_structure_id=str(row["source_structure_id"]).strip(),
                source_kind=str(row["source_kind"]).strip(),
                source_path=str(row["source_path"]).strip(),
                preserved_chains=str(row["preserved_chains"]).strip(),
                designed_chains=str(row["designed_chains"]).strip(),
                fixed_context_chains=str(row["fixed_context_chains"]).strip(),
                preserved_metal_identity=str(row["preserved_metal_identity"]).strip(),
                preserved_metal_site_count=int(str(row["preserved_metal_site_count"]).strip()),
                preserved_solvent_residue_count=int(str(row["preserved_solvent_residue_count"]).strip()),
                redesigned_residue_count=int(str(row["redesigned_residue_count"]).strip()),
                pdb_path=str(row["pdb_path"]).strip(),
                redesigned_residues_path=str(row["redesigned_residues_path"]).strip(),
            )
            for row in reader
        ]
    return tuple(sorted(rows, key=lambda row: row.shortlist_rank))


def _load_ligandmpnn_redesign_rows(path: Path) -> tuple[LigandMPNNRedesignPositionRow, ...]:
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            LigandMPNNRedesignPositionRow(
                shortlist_rank=int(str(row["shortlist_rank"]).strip()),
                candidate_id=str(row["candidate_id"]).strip(),
                campaign_id=str(row["campaign_id"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                chain_id=str(row["chain_id"]).strip(),
                sequence_index=int(str(row["sequence_index"]).strip()),
                residue_seq=int(str(row["residue_seq"]).strip()),
                insertion_code=str(row["insertion_code"]).strip(),
                ligandmpnn_residue_id=str(row["ligandmpnn_residue_id"]).strip(),
                residue_name=str(row["residue_name"]).strip(),
                native_amino_acid=str(row["native_amino_acid"]).strip(),
                designed_amino_acid=str(row["designed_amino_acid"]).strip(),
                mutation_token=str(row["mutation_token"]).strip(),
                canonical_family_position=(
                    int(str(row["canonical_family_position"]).strip())
                    if str(row["canonical_family_position"]).strip()
                    else None
                ),
                am1_mature_position=(
                    int(str(row["am1_mature_position"]).strip())
                    if str(row["am1_mature_position"]).strip()
                    else None
                ),
            )
            for row in reader
        ]
    return tuple(rows)


def discover_ligandmpnn_selection_candidates(
    *,
    proteinmpnn_shortlist_path: Path,
    manifest_path: Path,
    redesign_positions_path: Path,
    output_root: Path,
) -> tuple[LigandMPNNSelectionCandidate, ...]:
    """Resolve Phase 4A inputs and Phase 4B outputs into Phase 4C selection candidates."""
    shortlist_rows = _load_proteinmpnn_shortlist_rows(proteinmpnn_shortlist_path)
    shortlist_map = {row.candidate_id: row for row in shortlist_rows}
    manifest_rows = _load_ligandmpnn_input_manifest_rows(manifest_path)
    redesign_rows = _load_ligandmpnn_redesign_rows(redesign_positions_path)

    redesign_rows_by_candidate: dict[str, list[LigandMPNNRedesignPositionRow]] = defaultdict(list)
    for row in redesign_rows:
        redesign_rows_by_candidate[row.candidate_id].append(row)

    candidates: list[LigandMPNNSelectionCandidate] = []
    seen_candidate_ids: set[str] = set()
    for manifest_row in manifest_rows:
        if manifest_row.candidate_id in seen_candidate_ids:
            raise ValueError(f"Duplicate candidate_id in LigandMPNN input manifest: {manifest_row.candidate_id}")
        seen_candidate_ids.add(manifest_row.candidate_id)

        shortlist_row = shortlist_map.get(manifest_row.candidate_id)
        if shortlist_row is None:
            raise ValueError(f"Candidate {manifest_row.candidate_id} was missing from proteinmpnn_shortlist.csv")
        if shortlist_row.shortlist_rank != manifest_row.shortlist_rank:
            raise ValueError(
                f"Shortlist rank mismatch for {manifest_row.candidate_id}: "
                f"{shortlist_row.shortlist_rank} vs {manifest_row.shortlist_rank}"
            )
        if shortlist_row.campaign_id != manifest_row.campaign_id:
            raise ValueError(
                f"Campaign mismatch for {manifest_row.candidate_id}: "
                f"{shortlist_row.campaign_id} vs {manifest_row.campaign_id}"
            )
        if shortlist_row.backbone_id != manifest_row.backbone_id:
            raise ValueError(
                f"Backbone mismatch for {manifest_row.candidate_id}: "
                f"{shortlist_row.backbone_id} vs {manifest_row.backbone_id}"
            )

        candidate_redesign_rows = tuple(
            sorted(
                redesign_rows_by_candidate.get(manifest_row.candidate_id, []),
                key=lambda row: (row.chain_id, row.sequence_index, row.residue_seq, row.insertion_code),
            )
        )
        if len(candidate_redesign_rows) != manifest_row.redesigned_residue_count:
            raise ValueError(
                f"Redesign row count mismatch for {manifest_row.candidate_id}: "
                f"{len(candidate_redesign_rows)} vs {manifest_row.redesigned_residue_count}"
            )

        input_pdb_path = _repo_path(manifest_row.pdb_path)
        require_path(input_pdb_path)

        candidates.append(
            LigandMPNNSelectionCandidate(
                shortlist_rank=manifest_row.shortlist_rank,
                candidate_id=manifest_row.candidate_id,
                campaign_id=manifest_row.campaign_id,
                backbone_id=manifest_row.backbone_id,
                designed_chains=_split_csv_list(manifest_row.designed_chains),
                preserved_chains=_split_csv_list(manifest_row.preserved_chains),
                preserved_metal_identity=manifest_row.preserved_metal_identity,
                input_sequence=shortlist_row.designed_sequence,
                redesign_rows=candidate_redesign_rows,
                input_pdb_path=input_pdb_path,
                output_dir=output_root / manifest_row.candidate_id,
            )
        )
    return tuple(candidates)


def extract_designed_chain_sequence(
    sequence: str,
    *,
    preserved_chains: Sequence[str],
    designed_chains: Sequence[str],
) -> str:
    """Extract the designed-chain sequence from a LigandMPNN full output sequence."""
    parts = tuple(part.strip() for part in sequence.split(":"))
    if len(parts) != len(preserved_chains):
        raise ValueError(
            f"Expected {len(preserved_chains)} preserved-chain sequence parts, found {len(parts)} in {sequence!r}"
        )
    chain_to_part = {chain_id: part for chain_id, part in zip(preserved_chains, parts)}
    missing = [chain_id for chain_id in designed_chains if chain_id not in chain_to_part]
    if missing:
        raise ValueError(
            f"Designed chains {','.join(missing)} were not present in preserved chains {','.join(preserved_chains)}"
        )
    designed_parts = [chain_to_part[chain_id] for chain_id in designed_chains]
    if len(designed_parts) == 1:
        return designed_parts[0]
    return "/".join(designed_parts)


def summarize_mutations_vs_candidate_input(
    *,
    candidate_input_sequence: str,
    designed_chain_sequence: str,
    designed_chains: Sequence[str],
    redesign_rows: Sequence[LigandMPNNRedesignPositionRow],
) -> LigandMPNNMutationSummary:
    """Summarize LigandMPNN mutations relative to the Phase 4A candidate input sequence."""
    input_parts = _split_designed_sequence(candidate_input_sequence, designed_chains)
    designed_parts = _split_designed_sequence(designed_chain_sequence, designed_chains)
    redesign_rows_by_chain: dict[str, dict[int, LigandMPNNRedesignPositionRow]] = defaultdict(dict)
    for row in redesign_rows:
        redesign_rows_by_chain[row.chain_id][row.sequence_index] = row

    mutation_tokens: list[str] = []
    redesigned_ids: list[str] = []
    for chain_id, input_part, designed_part in zip(designed_chains, input_parts, designed_parts):
        if len(input_part) != len(designed_part):
            raise ValueError(
                f"Designed sequence length mismatch for chain {chain_id}: "
                f"{len(designed_part)} vs candidate input {len(input_part)}"
            )
        for index, (input_aa, output_aa) in enumerate(zip(input_part, designed_part), start=1):
            if input_aa == output_aa:
                continue
            redesign_row = redesign_rows_by_chain.get(chain_id, {}).get(index)
            if redesign_row is None:
                raise ValueError(
                    f"Unexpected mutation outside redesign positions on chain {chain_id} at sequence index {index}"
                )
            if redesign_row.designed_amino_acid != input_aa:
                raise ValueError(
                    f"Candidate input mismatch for {redesign_row.candidate_id} at {chain_id}{index}: "
                    f"{input_aa} vs redesign table {redesign_row.designed_amino_acid}"
                )
            mutation_tokens.append(f"{input_aa}{index}{output_aa}")
            redesigned_ids.append(redesign_row.ligandmpnn_residue_id)

    return LigandMPNNMutationSummary(
        mutation_count=len(mutation_tokens),
        mutation_string=",".join(mutation_tokens),
        redesigned_residue_identifiers=",".join(redesigned_ids),
    )


def parse_ligandmpnn_candidate_outputs(
    candidate: LigandMPNNSelectionCandidate,
) -> tuple[LigandMPNNGeneratedDesign, ...]:
    """Parse one candidate FASTA and extract designed-chain sequences plus mutation summaries."""
    parsed_output = parse_ligandmpnn_fasta(candidate.output_fasta_path)
    if parsed_output.native_record.name != candidate.candidate_id:
        raise ValueError(
            f"LigandMPNN FASTA name {parsed_output.native_record.name} did not match candidate {candidate.candidate_id}"
        )

    native_designed_chain_sequence = extract_designed_chain_sequence(
        parsed_output.native_record.sequence,
        preserved_chains=candidate.preserved_chains,
        designed_chains=candidate.designed_chains,
    )
    summarize_mutations_vs_candidate_input(
        candidate_input_sequence=candidate.input_sequence,
        designed_chain_sequence=native_designed_chain_sequence,
        designed_chains=candidate.designed_chains,
        redesign_rows=candidate.redesign_rows,
    )

    designs: list[LigandMPNNGeneratedDesign] = []
    for record in parsed_output.generated_records:
        if record.name != candidate.candidate_id:
            raise ValueError(
                f"LigandMPNN generated sequence name {record.name} did not match candidate {candidate.candidate_id}"
            )
        if record.seed != parsed_output.native_record.seed:
            raise ValueError(
                f"LigandMPNN seed mismatch for {candidate.candidate_id}: "
                f"{record.seed} vs native {parsed_output.native_record.seed}"
            )

        designed_chain_sequence = extract_designed_chain_sequence(
            record.sequence,
            preserved_chains=candidate.preserved_chains,
            designed_chains=candidate.designed_chains,
        )
        mutation_summary = summarize_mutations_vs_candidate_input(
            candidate_input_sequence=candidate.input_sequence,
            designed_chain_sequence=designed_chain_sequence,
            designed_chains=candidate.designed_chains,
            redesign_rows=candidate.redesign_rows,
        )
        backbone_pdb_path = candidate.output_dir / "backbones" / f"{candidate.candidate_id}_{record.design_id}.pdb"
        packed_pdb_path = (
            candidate.output_dir / "packed" / f"{candidate.candidate_id}{PACKED_SUFFIX}_{record.design_id}_1.pdb"
        )
        require_path(backbone_pdb_path)
        require_path(packed_pdb_path)

        designs.append(
            LigandMPNNGeneratedDesign(
                sequence_id=f"{candidate.candidate_id}_design_{record.design_id:02d}",
                candidate_id=candidate.candidate_id,
                campaign_id=candidate.campaign_id,
                backbone_id=candidate.backbone_id,
                preserved_metal_identity=candidate.preserved_metal_identity,
                candidate_shortlist_rank=candidate.shortlist_rank,
                design_id=record.design_id,
                overall_confidence=record.overall_confidence,
                ligand_confidence=record.ligand_confidence,
                seq_recovery=record.seq_recovery,
                mutation_count=mutation_summary.mutation_count,
                mutation_string=mutation_summary.mutation_string,
                redesigned_residue_identifiers=mutation_summary.redesigned_residue_identifiers,
                designed_chain_sequence=designed_chain_sequence,
                input_pdb_path=_display_path(candidate.input_pdb_path),
                output_fasta_path=_display_path(candidate.output_fasta_path),
                backbone_pdb_path=_display_path(backbone_pdb_path),
                packed_pdb_path=_display_path(packed_pdb_path),
            )
        )
    return tuple(designs)


def deduplicate_candidate_sequences(
    candidate: LigandMPNNSelectionCandidate,
    generated_designs: Sequence[LigandMPNNGeneratedDesign],
) -> CandidateUniqueSequenceResult:
    """Collapse duplicate designed-chain sequences within one candidate."""
    grouped_records: dict[str, list[LigandMPNNGeneratedDesign]] = defaultdict(list)
    for design in generated_designs:
        grouped_records[design.designed_chain_sequence].append(design)

    pending_rows: list[CandidateUniqueSequence] = []
    for designs in grouped_records.values():
        ranked_designs = sorted(designs, key=_generated_design_rank_key)
        representative = ranked_designs[0]
        pending_rows.append(
            CandidateUniqueSequence(
                sequence_id=representative.sequence_id,
                candidate_id=representative.candidate_id,
                campaign_id=representative.campaign_id,
                backbone_id=representative.backbone_id,
                preserved_metal_identity=representative.preserved_metal_identity,
                candidate_shortlist_rank=representative.candidate_shortlist_rank,
                candidate_unique_rank=0,
                design_id=representative.design_id,
                occurrence_count=len(designs),
                overall_confidence=representative.overall_confidence,
                ligand_confidence=representative.ligand_confidence,
                seq_recovery=representative.seq_recovery,
                mutation_count=representative.mutation_count,
                mutation_string=representative.mutation_string,
                redesigned_residue_identifiers=representative.redesigned_residue_identifiers,
                designed_chain_sequence=representative.designed_chain_sequence,
                input_pdb_path=representative.input_pdb_path,
                output_fasta_path=representative.output_fasta_path,
                backbone_pdb_path=representative.backbone_pdb_path,
                packed_pdb_path=representative.packed_pdb_path,
            )
        )

    ranked_rows = sorted(pending_rows, key=_candidate_unique_sequence_rank_key)
    unique_rows = tuple(
        CandidateUniqueSequence(
            sequence_id=row.sequence_id,
            candidate_id=row.candidate_id,
            campaign_id=row.campaign_id,
            backbone_id=row.backbone_id,
            preserved_metal_identity=row.preserved_metal_identity,
            candidate_shortlist_rank=row.candidate_shortlist_rank,
            candidate_unique_rank=index,
            design_id=row.design_id,
            occurrence_count=row.occurrence_count,
            overall_confidence=row.overall_confidence,
            ligand_confidence=row.ligand_confidence,
            seq_recovery=row.seq_recovery,
            mutation_count=row.mutation_count,
            mutation_string=row.mutation_string,
            redesigned_residue_identifiers=row.redesigned_residue_identifiers,
            designed_chain_sequence=row.designed_chain_sequence,
            input_pdb_path=row.input_pdb_path,
            output_fasta_path=row.output_fasta_path,
            backbone_pdb_path=row.backbone_pdb_path,
            packed_pdb_path=row.packed_pdb_path,
        )
        for index, row in enumerate(ranked_rows, start=1)
    )
    return CandidateUniqueSequenceResult(
        candidate=candidate,
        raw_sequence_count=len(generated_designs),
        unique_rows=unique_rows,
    )


def collapse_global_unique_sequences(
    candidate_results: Sequence[CandidateUniqueSequenceResult],
    *,
    per_candidate_limit: int = PER_CANDIDATE_SHORTLIST_LIMIT,
) -> tuple[LigandMPNNUniqueSequenceRow, ...]:
    """Collapse exact designed-chain duplicates across candidates into a global unique table."""
    grouped_rows: dict[str, list[CandidateUniqueSequence]] = defaultdict(list)
    for result in candidate_results:
        for row in result.unique_rows:
            grouped_rows[row.designed_chain_sequence].append(row)

    pending_rows: list[LigandMPNNUniqueSequenceRow] = []
    for duplicate_rows in grouped_rows.values():
        ranked_rows = sorted(
            duplicate_rows,
            key=lambda row: _cross_candidate_representative_key(row, per_candidate_limit=per_candidate_limit),
        )
        representative = ranked_rows[0]
        source_candidate_ids = tuple(
            row.candidate_id for row in sorted(duplicate_rows, key=_source_candidate_order_key)
        )
        pending_rows.append(
            LigandMPNNUniqueSequenceRow(
                unique_sequence_rank=0,
                sequence_id=representative.sequence_id,
                candidate_id=representative.candidate_id,
                campaign_id=representative.campaign_id,
                backbone_id=representative.backbone_id,
                preserved_metal_identity=representative.preserved_metal_identity,
                candidate_shortlist_rank=representative.candidate_shortlist_rank,
                candidate_unique_rank=representative.candidate_unique_rank,
                candidate_prefilter_selected=(representative.candidate_unique_rank <= per_candidate_limit),
                design_id=representative.design_id,
                total_occurrence_count=sum(row.occurrence_count for row in duplicate_rows),
                source_candidate_count=len({row.candidate_id for row in duplicate_rows}),
                source_candidate_ids=",".join(source_candidate_ids),
                overall_confidence=representative.overall_confidence,
                ligand_confidence=representative.ligand_confidence,
                seq_recovery=representative.seq_recovery,
                mutation_count=representative.mutation_count,
                mutation_string=representative.mutation_string,
                redesigned_residue_identifiers=representative.redesigned_residue_identifiers,
                designed_chain_sequence=representative.designed_chain_sequence,
                input_pdb_path=representative.input_pdb_path,
                output_fasta_path=representative.output_fasta_path,
                backbone_pdb_path=representative.backbone_pdb_path,
                packed_pdb_path=representative.packed_pdb_path,
            )
        )

    ranked_pending_rows = sorted(pending_rows, key=_global_unique_sequence_rank_key)
    return tuple(
        LigandMPNNUniqueSequenceRow(
            unique_sequence_rank=index,
            sequence_id=row.sequence_id,
            candidate_id=row.candidate_id,
            campaign_id=row.campaign_id,
            backbone_id=row.backbone_id,
            preserved_metal_identity=row.preserved_metal_identity,
            candidate_shortlist_rank=row.candidate_shortlist_rank,
            candidate_unique_rank=row.candidate_unique_rank,
            candidate_prefilter_selected=row.candidate_prefilter_selected,
            design_id=row.design_id,
            total_occurrence_count=row.total_occurrence_count,
            source_candidate_count=row.source_candidate_count,
            source_candidate_ids=row.source_candidate_ids,
            overall_confidence=row.overall_confidence,
            ligand_confidence=row.ligand_confidence,
            seq_recovery=row.seq_recovery,
            mutation_count=row.mutation_count,
            mutation_string=row.mutation_string,
            redesigned_residue_identifiers=row.redesigned_residue_identifiers,
            designed_chain_sequence=row.designed_chain_sequence,
            input_pdb_path=row.input_pdb_path,
            output_fasta_path=row.output_fasta_path,
            backbone_pdb_path=row.backbone_pdb_path,
            packed_pdb_path=row.packed_pdb_path,
        )
        for index, row in enumerate(ranked_pending_rows, start=1)
    )


def build_ligandmpnn_shortlist(
    unique_rows: Sequence[LigandMPNNUniqueSequenceRow],
    *,
    per_candidate_limit: int = PER_CANDIDATE_SHORTLIST_LIMIT,
    total_limit: int = FINAL_SHORTLIST_LIMIT,
) -> ShortlistSelectionResult:
    """Build a deterministic, balanced shortlist from globally unique LigandMPNN sequences."""
    eligible_rows = [
        row
        for row in unique_rows
        if row.candidate_prefilter_selected and row.candidate_unique_rank <= per_candidate_limit
    ]
    pool_by_campaign: dict[str, list[LigandMPNNUniqueSequenceRow]] = defaultdict(list)
    for row in eligible_rows:
        pool_by_campaign[row.campaign_id].append(row)
    for campaign_id in pool_by_campaign:
        pool_by_campaign[campaign_id] = sorted(pool_by_campaign[campaign_id], key=_shortlist_pool_row_key)

    selected_rows: list[LigandMPNNUniqueSequenceRow] = []
    selected_sequence_ids: set[str] = set()
    retention_reasons: dict[str, str] = {}

    _retain_required_rows(
        selected_rows=selected_rows,
        selected_sequence_ids=selected_sequence_ids,
        retention_reasons=retention_reasons,
        pool=pool_by_campaign.get(INTERFACE_CAMPAIGN_ID, []),
        limit=min(total_limit, INTERFACE_CAMPAIGN_MINIMUM),
        reason="required interface campaign coverage",
    )
    _retain_required_rows(
        selected_rows=selected_rows,
        selected_sequence_ids=selected_sequence_ids,
        retention_reasons=retention_reasons,
        pool=pool_by_campaign.get(POCKET_CAMPAIGN_ID, []),
        limit=min(max(total_limit - len(selected_rows), 0), POCKET_CAMPAIGN_MINIMUM),
        reason="required pocket campaign coverage",
    )
    _retain_required_rows(
        selected_rows=selected_rows,
        selected_sequence_ids=selected_sequence_ids,
        retention_reasons=retention_reasons,
        pool=sorted(
            [
                row
                for row in eligible_rows
                if row.campaign_id.startswith(AM1_CAMPAIGN_PREFIX)
            ],
            key=_shortlist_pool_row_key,
        ),
        limit=min(max(total_limit - len(selected_rows), 0), AM1_CAMPAIGN_MINIMUM),
        reason="required AM1/Mex coverage",
    )

    campaign_order = _shortlist_campaign_order(pool_by_campaign)
    round_number = 1
    while len(selected_rows) < total_limit:
        round_added = False
        for campaign_id in campaign_order:
            next_row = _next_available_row(
                pool_by_campaign.get(campaign_id, ()),
                selected_sequence_ids=selected_sequence_ids,
            )
            if next_row is None:
                continue
            selected_rows.append(next_row)
            selected_sequence_ids.add(next_row.sequence_id)
            retention_reasons[next_row.sequence_id] = f"balanced round-robin fill (round {round_number})"
            round_added = True
            if len(selected_rows) >= total_limit:
                break
        if not round_added:
            break
        round_number += 1

    shortlist_rows = tuple(
        LigandMPNNShortlistRow(
            shortlist_rank=index,
            sequence_id=row.sequence_id,
            candidate_id=row.candidate_id,
            campaign_id=row.campaign_id,
            backbone_id=row.backbone_id,
            preserved_metal_identity=row.preserved_metal_identity,
            candidate_shortlist_rank=row.candidate_shortlist_rank,
            candidate_unique_rank=row.candidate_unique_rank,
            design_id=row.design_id,
            total_occurrence_count=row.total_occurrence_count,
            source_candidate_count=row.source_candidate_count,
            source_candidate_ids=row.source_candidate_ids,
            overall_confidence=row.overall_confidence,
            ligand_confidence=row.ligand_confidence,
            seq_recovery=row.seq_recovery,
            mutation_count=row.mutation_count,
            mutation_string=row.mutation_string,
            redesigned_residue_identifiers=row.redesigned_residue_identifiers,
            designed_chain_sequence=row.designed_chain_sequence,
            input_pdb_path=row.input_pdb_path,
            output_fasta_path=row.output_fasta_path,
            backbone_pdb_path=row.backbone_pdb_path,
            packed_pdb_path=row.packed_pdb_path,
            retention_reason=retention_reasons[row.sequence_id],
        )
        for index, row in enumerate(selected_rows, start=1)
    )
    return ShortlistSelectionResult(shortlist_rows=shortlist_rows)


def render_ligandmpnn_shortlist_markdown(
    *,
    candidate_results: Sequence[CandidateUniqueSequenceResult],
    unique_rows: Sequence[LigandMPNNUniqueSequenceRow],
    shortlist_result: ShortlistSelectionResult,
) -> str:
    """Render the Phase 4C LigandMPNN shortlist report."""
    retained_counts = Counter(row.candidate_id for row in shortlist_result.shortlist_rows)
    raw_sequence_count = sum(result.raw_sequence_count for result in candidate_results)
    within_candidate_unique_count = sum(len(result.unique_rows) for result in candidate_results)
    within_candidate_duplicates = raw_sequence_count - within_candidate_unique_count
    cross_candidate_duplicates = within_candidate_unique_count - len(unique_rows)
    cross_candidate_duplicate_groups = sum(1 for row in unique_rows if row.source_candidate_count > 1)

    lines = [
        "# LigandMPNN Candidate Shortlist",
        "",
        "Deterministic Phase 4C deduplication of existing LigandMPNN smoke outputs into a Rosetta-ready shortlist.",
        "",
        f"- Per-candidate prefilter: top `{PER_CANDIDATE_SHORTLIST_LIMIT}` unique designed-chain sequences by ligand_confidence then overall_confidence.",
        f"- Final shortlist limit: `{FINAL_SHORTLIST_LIMIT}` designs.",
        (
            f"- `{INTERFACE_CAMPAIGN_ID}` contributes at least "
            f"`{INTERFACE_CAMPAIGN_MINIMUM}` designs when available."
        ),
        (
            f"- `{POCKET_CAMPAIGN_ID}` contributes at least "
            f"`{POCKET_CAMPAIGN_MINIMUM}` design when available."
        ),
        (
            f"- `{AM1_CAMPAIGN_PREFIX}*` campaigns contribute at least "
            f"`{AM1_CAMPAIGN_MINIMUM}` design when available."
        ),
        "",
        "## Candidate Counts",
        "",
        "| candidate_id | campaign_id | raw_designs | unique_designed_chain_sequences | top2_considered | retained_shortlist |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for result in candidate_results:
        lines.append(
            "| "
            f"{result.candidate.candidate_id} | "
            f"{result.candidate.campaign_id} | "
            f"{result.raw_sequence_count} | "
            f"{len(result.unique_rows)} | "
            f"{min(PER_CANDIDATE_SHORTLIST_LIMIT, len(result.unique_rows))} | "
            f"{retained_counts[result.candidate.candidate_id]} |"
        )

    lines.extend(
        [
            "",
            "## Duplicate Counts",
            "",
            "| scope | count |",
            "| --- | ---: |",
            f"| raw generated designs | {raw_sequence_count} |",
            f"| within-candidate unique designed-chain sequences | {within_candidate_unique_count} |",
            f"| globally unique designed-chain sequences | {len(unique_rows)} |",
            f"| within-candidate duplicates collapsed | {within_candidate_duplicates} |",
            f"| cross-candidate duplicates collapsed | {cross_candidate_duplicates} |",
            f"| cross-candidate duplicate groups | {cross_candidate_duplicate_groups} |",
            "",
            "## Retained Shortlist",
            "",
            "| rank | sequence_id | candidate_id | campaign_id | ligand_confidence | overall_confidence | seq_recovery | mutations_vs_candidate_input | redesigned_residue_ids | retained_because |",
            "| ---: | --- | --- | --- | ---: | ---: | ---: | --- | --- | --- |",
        ]
    )
    for row in shortlist_result.shortlist_rows:
        lines.append(
            "| "
            f"{row.shortlist_rank} | "
            f"{row.sequence_id} | "
            f"{row.candidate_id} | "
            f"{row.campaign_id} | "
            f"{row.ligand_confidence:.4f} | "
            f"{row.overall_confidence:.4f} | "
            f"{row.seq_recovery:.4f} | "
            f"{row.mutation_string or 'native'} | "
            f"{row.redesigned_residue_identifiers or 'native'} | "
            f"{row.retention_reason} |"
        )
    lines.append("")
    return "\n".join(lines)


def render_ligandmpnn_shortlist_fasta(shortlist_rows: Sequence[LigandMPNNShortlistRow]) -> str:
    """Render the Phase 4C shortlist FASTA in shortlist rank order."""
    lines: list[str] = []
    for row in shortlist_rows:
        lines.append(
            ">"
            f"{row.sequence_id}"
            f"|rank={row.shortlist_rank}"
            f"|candidate={row.candidate_id}"
            f"|campaign={row.campaign_id}"
            f"|ligand_confidence={row.ligand_confidence:.4f}"
            f"|overall_confidence={row.overall_confidence:.4f}"
            f"|seq_recovery={row.seq_recovery:.4f}"
            f"|mutations={row.mutation_string or 'native'}"
        )
        lines.append(row.designed_chain_sequence)
    return "\n".join(lines) + ("\n" if lines else "")


def select_ligandmpnn_candidates(
    *,
    proteinmpnn_shortlist_path: Path,
    manifest_path: Path,
    redesign_positions_path: Path,
    smoke_output_root: Path,
    unique_sequences_path: Path,
    shortlist_path: Path,
    report_path: Path,
    shortlist_fasta_path: Path,
    rosetta_shortlist_path: Path,
) -> tuple[tuple[LigandMPNNUniqueSequenceRow, ...], tuple[LigandMPNNShortlistRow, ...]]:
    """Run the deterministic Phase 4C selection workflow and write repo artifacts."""
    candidates = discover_ligandmpnn_selection_candidates(
        proteinmpnn_shortlist_path=proteinmpnn_shortlist_path,
        manifest_path=manifest_path,
        redesign_positions_path=redesign_positions_path,
        output_root=smoke_output_root,
    )
    if not candidates:
        raise ValueError(f"No LigandMPNN candidates discovered from {manifest_path}")

    candidate_results: list[CandidateUniqueSequenceResult] = []
    for candidate in candidates:
        require_path(candidate.output_fasta_path)
        generated_designs = parse_ligandmpnn_candidate_outputs(candidate)
        candidate_results.append(
            deduplicate_candidate_sequences(
                candidate=candidate,
                generated_designs=generated_designs,
            )
        )

    unique_rows = collapse_global_unique_sequences(candidate_results)
    shortlist_result = build_ligandmpnn_shortlist(unique_rows)

    write_csv_rows(unique_sequences_path, unique_rows)
    write_csv_rows(shortlist_path, shortlist_result.shortlist_rows)
    atomic_write_text(
        report_path,
        render_ligandmpnn_shortlist_markdown(
            candidate_results=candidate_results,
            unique_rows=unique_rows,
            shortlist_result=shortlist_result,
        ),
    )
    atomic_write_text(
        shortlist_fasta_path,
        render_ligandmpnn_shortlist_fasta(shortlist_result.shortlist_rows),
    )
    write_yaml(
        rosetta_shortlist_path,
        _build_rosetta_shortlist_payload(shortlist_result.shortlist_rows),
    )
    return unique_rows, shortlist_result.shortlist_rows


def _generated_design_rank_key(
    design: LigandMPNNGeneratedDesign,
) -> tuple[float, float, float, int, int, str]:
    return (
        -design.ligand_confidence,
        -design.overall_confidence,
        -design.seq_recovery,
        design.mutation_count,
        design.design_id,
        design.sequence_id,
    )


def _candidate_unique_sequence_rank_key(
    row: CandidateUniqueSequence,
) -> tuple[float, float, float, int, int, str, str]:
    return (
        -row.ligand_confidence,
        -row.overall_confidence,
        -row.seq_recovery,
        row.mutation_count,
        row.design_id,
        row.sequence_id,
        row.designed_chain_sequence,
    )


def _cross_candidate_representative_key(
    row: CandidateUniqueSequence,
    *,
    per_candidate_limit: int,
) -> tuple[bool, float, float, float, int, int, int, str, str]:
    return (
        row.candidate_unique_rank > per_candidate_limit,
        -row.ligand_confidence,
        -row.overall_confidence,
        -row.seq_recovery,
        row.candidate_unique_rank,
        row.candidate_shortlist_rank,
        row.design_id,
        row.sequence_id,
        row.candidate_id,
    )


def _global_unique_sequence_rank_key(
    row: LigandMPNNUniqueSequenceRow,
) -> tuple[bool, float, float, float, int, int, int, str, str]:
    return (
        not row.candidate_prefilter_selected,
        -row.ligand_confidence,
        -row.overall_confidence,
        -row.seq_recovery,
        row.candidate_unique_rank,
        row.candidate_shortlist_rank,
        row.design_id,
        row.sequence_id,
        row.candidate_id,
    )


def _shortlist_pool_row_key(
    row: LigandMPNNUniqueSequenceRow,
) -> tuple[int, float, float, float, int, int, str]:
    return (
        row.candidate_unique_rank,
        -row.ligand_confidence,
        -row.overall_confidence,
        -row.seq_recovery,
        row.candidate_shortlist_rank,
        row.design_id,
        row.sequence_id,
    )


def _shortlist_campaign_order(
    pool_by_campaign: dict[str, list[LigandMPNNUniqueSequenceRow]],
) -> tuple[str, ...]:
    am1_campaigns = sorted(
        campaign_id
        for campaign_id in pool_by_campaign
        if campaign_id.startswith(AM1_CAMPAIGN_PREFIX) and campaign_id != INTERFACE_CAMPAIGN_ID
    )
    other_campaigns = sorted(
        campaign_id
        for campaign_id in pool_by_campaign
        if campaign_id not in {INTERFACE_CAMPAIGN_ID, POCKET_CAMPAIGN_ID}
        and not campaign_id.startswith(AM1_CAMPAIGN_PREFIX)
    )

    ordered: list[str] = []
    if INTERFACE_CAMPAIGN_ID in pool_by_campaign:
        ordered.append(INTERFACE_CAMPAIGN_ID)
    if POCKET_CAMPAIGN_ID in pool_by_campaign:
        ordered.append(POCKET_CAMPAIGN_ID)
    ordered.extend(am1_campaigns)
    ordered.extend(other_campaigns)
    return tuple(ordered)


def _retain_required_rows(
    *,
    selected_rows: list[LigandMPNNUniqueSequenceRow],
    selected_sequence_ids: set[str],
    retention_reasons: dict[str, str],
    pool: Sequence[LigandMPNNUniqueSequenceRow],
    limit: int,
    reason: str,
) -> None:
    retained = 0
    for row in pool:
        if retained >= limit:
            break
        if row.sequence_id in selected_sequence_ids:
            continue
        selected_rows.append(row)
        selected_sequence_ids.add(row.sequence_id)
        retention_reasons[row.sequence_id] = reason
        retained += 1


def _next_available_row(
    pool: Iterable[LigandMPNNUniqueSequenceRow],
    *,
    selected_sequence_ids: set[str],
) -> LigandMPNNUniqueSequenceRow | None:
    for row in pool:
        if row.sequence_id in selected_sequence_ids:
            continue
        return row
    return None


def _source_candidate_order_key(row: CandidateUniqueSequence) -> tuple[int, str]:
    return row.candidate_shortlist_rank, row.candidate_id


def _build_rosetta_shortlist_payload(
    shortlist_rows: Sequence[LigandMPNNShortlistRow],
) -> dict[str, Any]:
    return {
        "version": 1,
        "phase": "4C",
        "source_artifacts": {
            "proteinmpnn_shortlist": "results/tables/proteinmpnn_shortlist.csv",
            "ligandmpnn_input_manifest": "results/tables/ligandmpnn_input_manifest.csv",
            "ligandmpnn_redesign_positions": "results/tables/ligandmpnn_redesign_positions.csv",
            "ligandmpnn_smoke_summary": "results/tables/ligandmpnn_smoke_summary.csv",
            "ligandmpnn_sequence_catalog": "results/tables/ligandmpnn_sequence_catalog.csv",
        },
        "shortlist_rules": {
            "per_candidate_limit": PER_CANDIDATE_SHORTLIST_LIMIT,
            "total_limit": FINAL_SHORTLIST_LIMIT,
            "required_campaign_minima": {
                INTERFACE_CAMPAIGN_ID: INTERFACE_CAMPAIGN_MINIMUM,
                POCKET_CAMPAIGN_ID: POCKET_CAMPAIGN_MINIMUM,
                f"{AM1_CAMPAIGN_PREFIX}*": AM1_CAMPAIGN_MINIMUM,
            },
        },
        "designs": [
            {
                "shortlist_rank": row.shortlist_rank,
                "sequence_id": row.sequence_id,
                "candidate_id": row.candidate_id,
                "campaign_id": row.campaign_id,
                "backbone_id": row.backbone_id,
                "preserved_metal_identity": row.preserved_metal_identity,
                "candidate_shortlist_rank": row.candidate_shortlist_rank,
                "candidate_unique_rank": row.candidate_unique_rank,
                "design_id": row.design_id,
                "source_candidate_count": row.source_candidate_count,
                "source_candidate_ids": row.source_candidate_ids.split(",") if row.source_candidate_ids else [],
                "ligand_confidence": row.ligand_confidence,
                "overall_confidence": row.overall_confidence,
                "seq_recovery": row.seq_recovery,
                "mutation_count": row.mutation_count,
                "mutation_string": row.mutation_string,
                "redesigned_residue_identifiers": (
                    row.redesigned_residue_identifiers.split(",") if row.redesigned_residue_identifiers else []
                ),
                "designed_chain_sequence": row.designed_chain_sequence,
                "ligandmpnn_input_pdb_path": row.input_pdb_path,
                "ligandmpnn_output_fasta_path": row.output_fasta_path,
                "ligandmpnn_backbone_pdb_path": row.backbone_pdb_path,
                "ligandmpnn_packed_pdb_path": row.packed_pdb_path,
                "rosetta_starting_pdb_path": row.packed_pdb_path,
                "retention_reason": row.retention_reason,
            }
            for row in shortlist_rows
        ],
    }
