"""Phase 7E LigandMPNN round-2 deduplication and reduced Rosetta triage shortlist."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

from lanm.analysis.ligandmpnn_round2_smoke import discover_round2_ligandmpnn_smoke_candidates
from lanm.analysis.ligandmpnn_smoke import SMOKE_PACKED_SUFFIX, parse_ligandmpnn_fasta
from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_text, write_csv_rows, write_yaml
from lanm.models import (
    LigandMPNNRound2RedesignPositionRow,
    LigandMPNNRound2ShortlistRow,
    LigandMPNNRound2UniqueSequenceRow,
    ProteinMPNNRound2ShortlistRow,
)
from lanm.paths import REPO_ROOT
from lanm.structure.pdb import read_pdb_atom_records

ROUND2_SHORTLIST_TOTAL_LIMIT = 4
AM1_CAMPAIGN_PREFIX = "am1_mex"
POCKET_CAMPAIGN_ID = "hans_pocket_ss_only"
INTERFACE_CAMPAIGN_ID = "hans_interface_ss_plus_if"


@dataclass(frozen=True, slots=True)
class LigandMPNNRound2SelectionCandidate:
    shortlist_rank: int
    candidate_id: str
    seed_rank: int
    seed_candidate_id: str
    campaign_id: str
    backbone_id: str
    topology_class: str
    designed_chains: tuple[str, ...]
    fixed_context_chains: tuple[str, ...]
    output_chain_order: tuple[str, ...]
    preserved_metal_identity: str
    redesign_rows: tuple[LigandMPNNRound2RedesignPositionRow, ...]
    input_pdb_path: Path
    output_dir: Path

    @property
    def output_fasta_path(self) -> Path:
        return self.output_dir / "seqs" / f"{self.candidate_id}.fa"


@dataclass(frozen=True, slots=True)
class LigandMPNNRound2MutationSummary:
    mutation_count: int
    mutation_string: str
    redesigned_residue_identifiers: str


@dataclass(frozen=True, slots=True)
class LigandMPNNRound2GeneratedDesign:
    sequence_id: str
    candidate_id: str
    seed_rank: int
    seed_candidate_id: str
    campaign_id: str
    backbone_id: str
    topology_class: str
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
class CandidateUniqueRound2Sequence:
    sequence_id: str
    candidate_id: str
    seed_rank: int
    seed_candidate_id: str
    campaign_id: str
    backbone_id: str
    topology_class: str
    preserved_metal_identity: str
    candidate_shortlist_rank: int
    candidate_unique_rank: int
    design_id: int
    occurrence_count: int
    mean_ligand_confidence: float
    mean_overall_confidence: float
    best_seq_recovery: float
    mutation_count: int
    mutation_string: str
    redesigned_residue_identifiers: str
    designed_chain_sequence: str
    input_pdb_path: str
    output_fasta_path: str
    backbone_pdb_path: str
    packed_pdb_path: str


@dataclass(frozen=True, slots=True)
class CandidateUniqueRound2SequenceResult:
    candidate: LigandMPNNRound2SelectionCandidate
    raw_sequence_count: int
    unique_rows: tuple[CandidateUniqueRound2Sequence, ...]


@dataclass(frozen=True, slots=True)
class Round2ShortlistSelectionResult:
    shortlist_rows: tuple[LigandMPNNRound2ShortlistRow, ...]


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


def _split_designed_sequence(sequence: str, designed_chains: Sequence[str]) -> tuple[str, ...]:
    if len(designed_chains) == 1:
        return (sequence.replace("/", ""),)
    parts = tuple(part.strip() for part in sequence.split("/"))
    if len(parts) != len(designed_chains):
        raise ValueError(
            f"Expected {len(designed_chains)} designed-chain sequence parts, found {len(parts)} in {sequence!r}"
        )
    return parts


def _polymer_chain_order(pdb_path: Path) -> tuple[str, ...]:
    atoms = read_pdb_atom_records(pdb_path)
    chain_ids: list[str] = []
    seen_chain_ids: set[str] = set()
    for atom in atoms:
        if atom.record_type != "ATOM":
            continue
        if atom.chain_id in seen_chain_ids:
            continue
        seen_chain_ids.add(atom.chain_id)
        chain_ids.append(atom.chain_id)
    if not chain_ids:
        raise ValueError(f"Expected at least one polymer chain in {pdb_path}")
    return tuple(chain_ids)


def _load_proteinmpnn_round2_shortlist_rows(path: Path) -> tuple[ProteinMPNNRound2ShortlistRow, ...]:
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            ProteinMPNNRound2ShortlistRow(
                shortlist_rank=int(str(row["shortlist_rank"]).strip()),
                candidate_id=str(row["candidate_id"]).strip(),
                seed_rank=int(str(row["seed_rank"]).strip()),
                seed_candidate_id=str(row["seed_candidate_id"]).strip(),
                campaign_id=str(row["campaign_id"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                topology_class=str(row["topology_class"]).strip(),
                design_set_name=str(row["design_set_name"]).strip(),
                designed_chains=str(row["designed_chains"]).strip(),
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
                redesigned_residue_ids=str(row["redesigned_residue_ids"]).strip(),
                redesigned_canonical_positions=str(row["redesigned_canonical_positions"]).strip(),
                redesigned_am1_positions=str(row["redesigned_am1_positions"]).strip(),
                retention_reason=str(row["retention_reason"]).strip(),
            )
            for row in reader
        ]
    return tuple(sorted(rows, key=lambda row: row.shortlist_rank))


def _load_ligandmpnn_round2_redesign_rows(path: Path) -> tuple[LigandMPNNRound2RedesignPositionRow, ...]:
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            LigandMPNNRound2RedesignPositionRow(
                shortlist_rank=int(str(row["shortlist_rank"]).strip()),
                candidate_id=str(row["candidate_id"]).strip(),
                seed_rank=int(str(row["seed_rank"]).strip()),
                seed_candidate_id=str(row["seed_candidate_id"]).strip(),
                campaign_id=str(row["campaign_id"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                chain_id=str(row["chain_id"]).strip(),
                sequence_index=int(str(row["sequence_index"]).strip()),
                residue_seq=int(str(row["residue_seq"]).strip()),
                insertion_code=str(row["insertion_code"]).strip(),
                ligandmpnn_residue_id=str(row["ligandmpnn_residue_id"]).strip(),
                residue_name=str(row["residue_name"]).strip(),
                seed_amino_acid=str(row["seed_amino_acid"]).strip(),
                designed_amino_acid=str(row["designed_amino_acid"]).strip(),
                mutation_token=str(row["mutation_token"]).strip(),
                canonical_family_position=int(str(row["canonical_family_position"]).strip()),
                am1_mature_position=int(str(row["am1_mature_position"]).strip()),
            )
            for row in reader
        ]
    return tuple(rows)


def discover_ligandmpnn_round2_selection_candidates(
    *,
    proteinmpnn_round2_shortlist_path: Path,
    manifest_path: Path,
    redesign_positions_path: Path,
    smoke_output_root: Path,
) -> tuple[LigandMPNNRound2SelectionCandidate, ...]:
    """Resolve round-2 LigandMPNN smoke outputs into Phase 7E selection candidates."""
    smoke_candidates = discover_round2_ligandmpnn_smoke_candidates(
        manifest_path=manifest_path,
        output_root=smoke_output_root,
    )
    shortlist_rows = _load_proteinmpnn_round2_shortlist_rows(proteinmpnn_round2_shortlist_path)
    shortlist_by_candidate = {row.candidate_id: row for row in shortlist_rows}
    redesign_rows = _load_ligandmpnn_round2_redesign_rows(redesign_positions_path)

    redesign_rows_by_candidate: dict[str, list[LigandMPNNRound2RedesignPositionRow]] = defaultdict(list)
    for row in redesign_rows:
        redesign_rows_by_candidate[row.candidate_id].append(row)

    candidates: list[LigandMPNNRound2SelectionCandidate] = []
    for smoke_candidate in smoke_candidates:
        shortlist_row = shortlist_by_candidate.get(smoke_candidate.candidate_id)
        if shortlist_row is None:
            raise ValueError(
                f"Round-2 LigandMPNN candidate {smoke_candidate.candidate_id} was missing from proteinmpnn_round2_shortlist.csv"
            )
        if shortlist_row.shortlist_rank != smoke_candidate.shortlist_rank:
            raise ValueError(
                f"Shortlist rank mismatch for {smoke_candidate.candidate_id}: "
                f"{shortlist_row.shortlist_rank} vs {smoke_candidate.shortlist_rank}"
            )
        if shortlist_row.seed_rank != smoke_candidate.seed_rank:
            raise ValueError(
                f"Seed rank mismatch for {smoke_candidate.candidate_id}: "
                f"{shortlist_row.seed_rank} vs {smoke_candidate.seed_rank}"
            )
        if shortlist_row.seed_candidate_id != smoke_candidate.seed_candidate_id:
            raise ValueError(
                f"Seed candidate mismatch for {smoke_candidate.candidate_id}: "
                f"{shortlist_row.seed_candidate_id} vs {smoke_candidate.seed_candidate_id}"
            )
        if shortlist_row.campaign_id != smoke_candidate.campaign_id:
            raise ValueError(
                f"Campaign mismatch for {smoke_candidate.candidate_id}: "
                f"{shortlist_row.campaign_id} vs {smoke_candidate.campaign_id}"
            )
        if shortlist_row.backbone_id != smoke_candidate.backbone_id:
            raise ValueError(
                f"Backbone mismatch for {smoke_candidate.candidate_id}: "
                f"{shortlist_row.backbone_id} vs {smoke_candidate.backbone_id}"
            )
        if shortlist_row.topology_class != smoke_candidate.topology_class:
            raise ValueError(
                f"Topology mismatch for {smoke_candidate.candidate_id}: "
                f"{shortlist_row.topology_class} vs {smoke_candidate.topology_class}"
            )

        candidate_redesign_rows = tuple(
            sorted(
                redesign_rows_by_candidate.get(smoke_candidate.candidate_id, []),
                key=lambda row: (row.chain_id, row.sequence_index, row.residue_seq, row.insertion_code),
            )
        )
        if len(candidate_redesign_rows) != len(smoke_candidate.redesigned_residue_ids):
            raise ValueError(
                f"Redesign row count mismatch for {smoke_candidate.candidate_id}: "
                f"{len(candidate_redesign_rows)} vs {len(smoke_candidate.redesigned_residue_ids)}"
            )

        output_chain_order = _polymer_chain_order(smoke_candidate.pdb_path)
        missing_designed_chains = [
            chain_id for chain_id in smoke_candidate.designed_chains if chain_id not in output_chain_order
        ]
        if missing_designed_chains:
            raise ValueError(
                f"Designed chains {','.join(missing_designed_chains)} were not present in {smoke_candidate.pdb_path}"
            )

        candidates.append(
            LigandMPNNRound2SelectionCandidate(
                shortlist_rank=smoke_candidate.shortlist_rank,
                candidate_id=smoke_candidate.candidate_id,
                seed_rank=smoke_candidate.seed_rank,
                seed_candidate_id=smoke_candidate.seed_candidate_id,
                campaign_id=smoke_candidate.campaign_id,
                backbone_id=smoke_candidate.backbone_id,
                topology_class=smoke_candidate.topology_class,
                designed_chains=smoke_candidate.designed_chains,
                fixed_context_chains=smoke_candidate.fixed_context_chains,
                output_chain_order=output_chain_order,
                preserved_metal_identity=smoke_candidate.preserved_metal_identity,
                redesign_rows=candidate_redesign_rows,
                input_pdb_path=smoke_candidate.pdb_path,
                output_dir=smoke_candidate.output_dir,
            )
        )
    return tuple(candidates)


def extract_designed_chain_sequence(
    sequence: str,
    *,
    output_chain_order: Sequence[str],
    designed_chains: Sequence[str],
) -> str:
    """Extract the designed-chain sequence from a LigandMPNN round-2 full output sequence."""
    parts = tuple(part.strip() for part in sequence.split(":"))
    if len(parts) != len(output_chain_order):
        raise ValueError(
            f"Expected {len(output_chain_order)} chain sequence parts, found {len(parts)} in {sequence!r}"
        )
    chain_to_part = {chain_id: part for chain_id, part in zip(output_chain_order, parts)}
    missing = [chain_id for chain_id in designed_chains if chain_id not in chain_to_part]
    if missing:
        raise ValueError(
            f"Designed chains {','.join(missing)} were not present in output chain order {','.join(output_chain_order)}"
        )
    designed_parts = [chain_to_part[chain_id] for chain_id in designed_chains]
    if len(designed_parts) == 1:
        return designed_parts[0]
    return "/".join(designed_parts)


def summarize_mutations_vs_round2_seed_input(
    *,
    round2_seed_input_sequence: str,
    designed_chain_sequence: str,
    designed_chains: Sequence[str],
    redesign_rows: Sequence[LigandMPNNRound2RedesignPositionRow],
) -> LigandMPNNRound2MutationSummary:
    """Summarize LigandMPNN mutations relative to the round-2 seed input sequence."""
    seed_parts = _split_designed_sequence(round2_seed_input_sequence, designed_chains)
    designed_parts = _split_designed_sequence(designed_chain_sequence, designed_chains)
    redesign_rows_by_chain: dict[str, dict[int, LigandMPNNRound2RedesignPositionRow]] = defaultdict(dict)
    for row in redesign_rows:
        redesign_rows_by_chain[row.chain_id][row.sequence_index] = row

    mutation_tokens: list[str] = []
    redesigned_ids: list[str] = []
    for chain_id, seed_part, designed_part in zip(designed_chains, seed_parts, designed_parts):
        if len(seed_part) != len(designed_part):
            raise ValueError(
                f"Designed sequence length mismatch for chain {chain_id}: "
                f"{len(designed_part)} vs round-2 seed input {len(seed_part)}"
            )
        for index, (seed_aa, designed_aa) in enumerate(zip(seed_part, designed_part), start=1):
            if seed_aa == designed_aa:
                continue
            redesign_row = redesign_rows_by_chain.get(chain_id, {}).get(index)
            if redesign_row is None:
                raise ValueError(
                    f"Unexpected mutation outside redesign positions on chain {chain_id} at sequence index {index}"
                )
            if redesign_row.seed_amino_acid != seed_aa:
                raise ValueError(
                    f"Round-2 seed input mismatch for {redesign_row.candidate_id} at "
                    f"{redesign_row.ligandmpnn_residue_id}: {seed_aa} vs redesign table {redesign_row.seed_amino_acid}"
                )
            mutation_tokens.append(f"{redesign_row.ligandmpnn_residue_id}:{seed_aa}>{designed_aa}")
            redesigned_ids.append(redesign_row.ligandmpnn_residue_id)

    return LigandMPNNRound2MutationSummary(
        mutation_count=len(mutation_tokens),
        mutation_string=",".join(mutation_tokens),
        redesigned_residue_identifiers=",".join(redesigned_ids),
    )


def parse_ligandmpnn_round2_candidate_outputs(
    candidate: LigandMPNNRound2SelectionCandidate,
) -> tuple[LigandMPNNRound2GeneratedDesign, ...]:
    """Parse one round-2 LigandMPNN FASTA and extract designed-chain sequences plus mutation summaries."""
    parsed_output = parse_ligandmpnn_fasta(candidate.output_fasta_path)
    if parsed_output.native_record.name != candidate.candidate_id:
        raise ValueError(
            f"LigandMPNN FASTA name {parsed_output.native_record.name} did not match candidate {candidate.candidate_id}"
        )

    native_designed_chain_sequence = extract_designed_chain_sequence(
        parsed_output.native_record.sequence,
        output_chain_order=candidate.output_chain_order,
        designed_chains=candidate.designed_chains,
    )

    designs: list[LigandMPNNRound2GeneratedDesign] = []
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
            output_chain_order=candidate.output_chain_order,
            designed_chains=candidate.designed_chains,
        )
        mutation_summary = summarize_mutations_vs_round2_seed_input(
            round2_seed_input_sequence=native_designed_chain_sequence,
            designed_chain_sequence=designed_chain_sequence,
            designed_chains=candidate.designed_chains,
            redesign_rows=candidate.redesign_rows,
        )
        backbone_pdb_path = candidate.output_dir / "backbones" / f"{candidate.candidate_id}_{record.design_id}.pdb"
        packed_pdb_path = (
            candidate.output_dir / "packed" / f"{candidate.candidate_id}{SMOKE_PACKED_SUFFIX}_{record.design_id}_1.pdb"
        )
        require_path(backbone_pdb_path)
        require_path(packed_pdb_path)

        designs.append(
            LigandMPNNRound2GeneratedDesign(
                sequence_id=f"{candidate.candidate_id}_design_{record.design_id:02d}",
                candidate_id=candidate.candidate_id,
                seed_rank=candidate.seed_rank,
                seed_candidate_id=candidate.seed_candidate_id,
                campaign_id=candidate.campaign_id,
                backbone_id=candidate.backbone_id,
                topology_class=candidate.topology_class,
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
    candidate: LigandMPNNRound2SelectionCandidate,
    generated_designs: Sequence[LigandMPNNRound2GeneratedDesign],
) -> CandidateUniqueRound2SequenceResult:
    """Collapse exact designed-chain duplicates within one round-2 candidate."""
    grouped_records: dict[str, list[LigandMPNNRound2GeneratedDesign]] = defaultdict(list)
    for design in generated_designs:
        grouped_records[design.designed_chain_sequence].append(design)

    pending_rows: list[CandidateUniqueRound2Sequence] = []
    for designs in grouped_records.values():
        ranked_designs = sorted(designs, key=_generated_design_rank_key)
        representative = ranked_designs[0]
        occurrence_count = len(designs)
        pending_rows.append(
            CandidateUniqueRound2Sequence(
                sequence_id=representative.sequence_id,
                candidate_id=representative.candidate_id,
                seed_rank=representative.seed_rank,
                seed_candidate_id=representative.seed_candidate_id,
                campaign_id=representative.campaign_id,
                backbone_id=representative.backbone_id,
                topology_class=representative.topology_class,
                preserved_metal_identity=representative.preserved_metal_identity,
                candidate_shortlist_rank=representative.candidate_shortlist_rank,
                candidate_unique_rank=0,
                design_id=representative.design_id,
                occurrence_count=occurrence_count,
                mean_ligand_confidence=_mean(design.ligand_confidence for design in designs),
                mean_overall_confidence=_mean(design.overall_confidence for design in designs),
                best_seq_recovery=max(design.seq_recovery for design in designs),
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
        CandidateUniqueRound2Sequence(
            sequence_id=row.sequence_id,
            candidate_id=row.candidate_id,
            seed_rank=row.seed_rank,
            seed_candidate_id=row.seed_candidate_id,
            campaign_id=row.campaign_id,
            backbone_id=row.backbone_id,
            topology_class=row.topology_class,
            preserved_metal_identity=row.preserved_metal_identity,
            candidate_shortlist_rank=row.candidate_shortlist_rank,
            candidate_unique_rank=index,
            design_id=row.design_id,
            occurrence_count=row.occurrence_count,
            mean_ligand_confidence=row.mean_ligand_confidence,
            mean_overall_confidence=row.mean_overall_confidence,
            best_seq_recovery=row.best_seq_recovery,
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
    return CandidateUniqueRound2SequenceResult(
        candidate=candidate,
        raw_sequence_count=len(generated_designs),
        unique_rows=unique_rows,
    )


def collapse_global_unique_sequences(
    candidate_results: Sequence[CandidateUniqueRound2SequenceResult],
) -> tuple[LigandMPNNRound2UniqueSequenceRow, ...]:
    """Collapse exact designed-chain duplicates across round-2 candidates into a global unique table."""
    grouped_rows: dict[str, list[CandidateUniqueRound2Sequence]] = defaultdict(list)
    for result in candidate_results:
        for row in result.unique_rows:
            grouped_rows[row.designed_chain_sequence].append(row)

    pending_rows: list[LigandMPNNRound2UniqueSequenceRow] = []
    for duplicate_rows in grouped_rows.values():
        ranked_rows = sorted(duplicate_rows, key=_cross_candidate_representative_key)
        representative = ranked_rows[0]
        total_occurrence_count = sum(row.occurrence_count for row in duplicate_rows)
        ordered_source_rows = sorted(duplicate_rows, key=_source_candidate_order_key)
        source_candidate_ids = tuple(dict.fromkeys(row.candidate_id for row in ordered_source_rows))
        pending_rows.append(
            LigandMPNNRound2UniqueSequenceRow(
                unique_sequence_rank=0,
                sequence_id=representative.sequence_id,
                candidate_id=representative.candidate_id,
                seed_rank=representative.seed_rank,
                seed_candidate_id=representative.seed_candidate_id,
                campaign_id=representative.campaign_id,
                backbone_id=representative.backbone_id,
                topology_class=representative.topology_class,
                preserved_metal_identity=representative.preserved_metal_identity,
                candidate_shortlist_rank=representative.candidate_shortlist_rank,
                candidate_unique_rank=representative.candidate_unique_rank,
                design_id=representative.design_id,
                total_occurrence_count=total_occurrence_count,
                source_candidate_count=len(source_candidate_ids),
                source_candidate_ids=",".join(source_candidate_ids),
                mean_ligand_confidence=(
                    sum(row.mean_ligand_confidence * row.occurrence_count for row in duplicate_rows)
                    / total_occurrence_count
                ),
                mean_overall_confidence=(
                    sum(row.mean_overall_confidence * row.occurrence_count for row in duplicate_rows)
                    / total_occurrence_count
                ),
                best_seq_recovery=max(row.best_seq_recovery for row in duplicate_rows),
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
        LigandMPNNRound2UniqueSequenceRow(
            unique_sequence_rank=index,
            sequence_id=row.sequence_id,
            candidate_id=row.candidate_id,
            seed_rank=row.seed_rank,
            seed_candidate_id=row.seed_candidate_id,
            campaign_id=row.campaign_id,
            backbone_id=row.backbone_id,
            topology_class=row.topology_class,
            preserved_metal_identity=row.preserved_metal_identity,
            candidate_shortlist_rank=row.candidate_shortlist_rank,
            candidate_unique_rank=row.candidate_unique_rank,
            design_id=row.design_id,
            total_occurrence_count=row.total_occurrence_count,
            source_candidate_count=row.source_candidate_count,
            source_candidate_ids=row.source_candidate_ids,
            mean_ligand_confidence=row.mean_ligand_confidence,
            mean_overall_confidence=row.mean_overall_confidence,
            best_seq_recovery=row.best_seq_recovery,
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


def build_ligandmpnn_round2_shortlist(
    unique_rows: Sequence[LigandMPNNRound2UniqueSequenceRow],
    *,
    total_limit: int = ROUND2_SHORTLIST_TOTAL_LIMIT,
) -> Round2ShortlistSelectionResult:
    """Build the reduced Phase 7E round-2 shortlist for Rosetta score-only triage."""
    selected_rows: list[LigandMPNNRound2UniqueSequenceRow] = []
    selected_sequence_ids: set[str] = set()
    retention_reasons: dict[str, str] = {}

    _retain_first_available(
        selected_rows=selected_rows,
        selected_sequence_ids=selected_sequence_ids,
        retention_reasons=retention_reasons,
        pool=sorted(
            [row for row in unique_rows if row.campaign_id.startswith(AM1_CAMPAIGN_PREFIX)],
            key=_global_unique_sequence_rank_key,
        ),
        reason="top AM1/Mex round-2 candidate",
        total_limit=total_limit,
    )
    _retain_first_available(
        selected_rows=selected_rows,
        selected_sequence_ids=selected_sequence_ids,
        retention_reasons=retention_reasons,
        pool=sorted(
            [row for row in unique_rows if row.campaign_id == POCKET_CAMPAIGN_ID],
            key=_global_unique_sequence_rank_key,
        ),
        reason="top Hans pocket round-2 candidate",
        total_limit=total_limit,
    )
    _retain_first_available(
        selected_rows=selected_rows,
        selected_sequence_ids=selected_sequence_ids,
        retention_reasons=retention_reasons,
        pool=sorted(
            [row for row in unique_rows if row.campaign_id == INTERFACE_CAMPAIGN_ID],
            key=_global_unique_sequence_rank_key,
        ),
        reason="top Hans interface-aware round-2 candidate",
        total_limit=total_limit,
    )

    remaining_pool = sorted(
        [row for row in unique_rows if row.sequence_id not in selected_sequence_ids],
        key=_global_unique_sequence_rank_key,
    )
    _retain_first_available(
        selected_rows=selected_rows,
        selected_sequence_ids=selected_sequence_ids,
        retention_reasons=retention_reasons,
        pool=remaining_pool,
        reason="best remaining unique candidate",
        total_limit=total_limit,
    )

    ranked_selected_rows = sorted(selected_rows, key=_global_unique_sequence_rank_key)
    shortlist_rows = tuple(
        LigandMPNNRound2ShortlistRow(
            shortlist_rank=index,
            sequence_id=row.sequence_id,
            candidate_id=row.candidate_id,
            seed_rank=row.seed_rank,
            seed_candidate_id=row.seed_candidate_id,
            campaign_id=row.campaign_id,
            backbone_id=row.backbone_id,
            topology_class=row.topology_class,
            preserved_metal_identity=row.preserved_metal_identity,
            candidate_shortlist_rank=row.candidate_shortlist_rank,
            candidate_unique_rank=row.candidate_unique_rank,
            design_id=row.design_id,
            total_occurrence_count=row.total_occurrence_count,
            source_candidate_count=row.source_candidate_count,
            source_candidate_ids=row.source_candidate_ids,
            mean_ligand_confidence=row.mean_ligand_confidence,
            mean_overall_confidence=row.mean_overall_confidence,
            best_seq_recovery=row.best_seq_recovery,
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
        for index, row in enumerate(ranked_selected_rows, start=1)
    )
    return Round2ShortlistSelectionResult(shortlist_rows=shortlist_rows)


def render_ligandmpnn_round2_shortlist_markdown(
    *,
    candidate_results: Sequence[CandidateUniqueRound2SequenceResult],
    unique_rows: Sequence[LigandMPNNRound2UniqueSequenceRow],
    shortlist_result: Round2ShortlistSelectionResult,
) -> str:
    """Render the Phase 7E reduced LigandMPNN round-2 shortlist report."""
    retained_counts = Counter(row.seed_candidate_id for row in shortlist_result.shortlist_rows)
    raw_sequence_count = sum(result.raw_sequence_count for result in candidate_results)
    within_candidate_unique_count = sum(len(result.unique_rows) for result in candidate_results)
    within_candidate_duplicates = raw_sequence_count - within_candidate_unique_count
    cross_candidate_duplicates = within_candidate_unique_count - len(unique_rows)

    raw_by_seed: Counter[str] = Counter()
    unique_sequences_by_seed: dict[str, set[str]] = defaultdict(set)
    campaign_by_seed: dict[str, str] = {}
    topology_by_seed: dict[str, str] = {}
    for result in candidate_results:
        seed_candidate_id = result.candidate.seed_candidate_id
        raw_by_seed[seed_candidate_id] += result.raw_sequence_count
        campaign_by_seed[seed_candidate_id] = result.candidate.campaign_id
        topology_by_seed[seed_candidate_id] = result.candidate.topology_class
        for row in result.unique_rows:
            unique_sequences_by_seed[seed_candidate_id].add(row.designed_chain_sequence)

    lines = [
        "# LigandMPNN Round-2 Shortlist",
        "",
        "Phase 7E deterministic deduplication of Phase 7D LigandMPNN smoke outputs into a reduced Rosetta score-only round-2 shortlist.",
        "",
        f"- Final shortlist limit: `{ROUND2_SHORTLIST_TOTAL_LIMIT}` unique designs.",
        f"- Top `{AM1_CAMPAIGN_PREFIX}*` round-2 candidate retained when available.",
        f"- Top `{POCKET_CAMPAIGN_ID}` round-2 candidate retained when available.",
        f"- Top `{INTERFACE_CAMPAIGN_ID}` round-2 candidate retained when available.",
        "- One additional best remaining unique candidate retained across all remaining candidates when available.",
        "- Ranking uses mean ligand confidence, then mean overall confidence, then smallest round-2 shortlist rank.",
        "- No Rosetta, MD, QM, or quantum steps were run in this phase.",
        "",
        "## Seed Counts",
        "",
        "| seed_candidate_id | campaign_id | topology_class | raw_designs | unique_designed_chain_sequences | duplicates_collapsed | retained_shortlist |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for seed_candidate_id, _ in sorted(
        ((seed_candidate_id, min(row.candidate.shortlist_rank for row in candidate_results if row.candidate.seed_candidate_id == seed_candidate_id))
         for seed_candidate_id in raw_by_seed),
        key=lambda item: (item[1], item[0]),
    ):
        unique_count = len(unique_sequences_by_seed[seed_candidate_id])
        raw_count = raw_by_seed[seed_candidate_id]
        lines.append(
            "| "
            f"{seed_candidate_id} | "
            f"{campaign_by_seed[seed_candidate_id]} | "
            f"{topology_by_seed[seed_candidate_id]} | "
            f"{raw_count} | "
            f"{unique_count} | "
            f"{raw_count - unique_count} | "
            f"{retained_counts[seed_candidate_id]} |"
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
            "",
            "## Retained Shortlist",
            "",
            "| rank | sequence_id | candidate_id | seed_scaffold | campaign_id | mean_ligand_confidence | mean_overall_confidence | best_seq_recovery | mutations_vs_round2_seed | redesigned_residue_ids | retained_because |",
            "| ---: | --- | --- | --- | --- | ---: | ---: | ---: | --- | --- | --- |",
        ]
    )
    for row in shortlist_result.shortlist_rows:
        lines.append(
            "| "
            f"{row.shortlist_rank} | "
            f"{row.sequence_id} | "
            f"{row.candidate_id} | "
            f"{row.seed_candidate_id} | "
            f"{row.campaign_id} | "
            f"{row.mean_ligand_confidence:.4f} | "
            f"{row.mean_overall_confidence:.4f} | "
            f"{row.best_seq_recovery:.4f} | "
            f"{row.mutation_string or 'native'} | "
            f"{row.redesigned_residue_identifiers or 'native'} | "
            f"{row.retention_reason} |"
        )
    lines.append("")
    return "\n".join(lines)


def render_ligandmpnn_round2_shortlist_fasta(shortlist_rows: Sequence[LigandMPNNRound2ShortlistRow]) -> str:
    """Render the Phase 7E reduced shortlist FASTA in shortlist rank order."""
    lines: list[str] = []
    for row in shortlist_rows:
        lines.append(
            ">"
            f"{row.sequence_id}"
            f"|rank={row.shortlist_rank}"
            f"|candidate={row.candidate_id}"
            f"|seed={row.seed_candidate_id}"
            f"|campaign={row.campaign_id}"
            f"|mean_ligand_confidence={row.mean_ligand_confidence:.4f}"
            f"|mean_overall_confidence={row.mean_overall_confidence:.4f}"
            f"|best_seq_recovery={row.best_seq_recovery:.4f}"
            f"|mutations={row.mutation_string or 'native'}"
        )
        lines.append(row.designed_chain_sequence)
    return "\n".join(lines) + ("\n" if lines else "")


def select_ligandmpnn_round2_candidates(
    *,
    proteinmpnn_round2_shortlist_path: Path,
    manifest_path: Path,
    redesign_positions_path: Path,
    smoke_output_root: Path,
    unique_sequences_path: Path,
    shortlist_path: Path,
    report_path: Path,
    shortlist_fasta_path: Path,
    rosetta_shortlist_path: Path,
) -> tuple[tuple[LigandMPNNRound2UniqueSequenceRow, ...], tuple[LigandMPNNRound2ShortlistRow, ...]]:
    """Run the deterministic Phase 7E selection workflow and write repo artifacts."""
    candidates = discover_ligandmpnn_round2_selection_candidates(
        proteinmpnn_round2_shortlist_path=proteinmpnn_round2_shortlist_path,
        manifest_path=manifest_path,
        redesign_positions_path=redesign_positions_path,
        smoke_output_root=smoke_output_root,
    )
    if not candidates:
        raise ValueError(f"No LigandMPNN round-2 candidates discovered from {manifest_path}")

    candidate_results: list[CandidateUniqueRound2SequenceResult] = []
    for candidate in candidates:
        require_path(candidate.output_fasta_path)
        generated_designs = parse_ligandmpnn_round2_candidate_outputs(candidate)
        candidate_results.append(
            deduplicate_candidate_sequences(
                candidate=candidate,
                generated_designs=generated_designs,
            )
        )

    unique_rows = collapse_global_unique_sequences(candidate_results)
    shortlist_result = build_ligandmpnn_round2_shortlist(unique_rows)

    write_csv_rows(unique_sequences_path, unique_rows)
    write_csv_rows(shortlist_path, shortlist_result.shortlist_rows)
    atomic_write_text(
        report_path,
        render_ligandmpnn_round2_shortlist_markdown(
            candidate_results=candidate_results,
            unique_rows=unique_rows,
            shortlist_result=shortlist_result,
        ),
    )
    atomic_write_text(
        shortlist_fasta_path,
        render_ligandmpnn_round2_shortlist_fasta(shortlist_result.shortlist_rows),
    )
    write_yaml(
        rosetta_shortlist_path,
        _build_rosetta_round2_shortlist_payload(shortlist_result.shortlist_rows),
    )
    return unique_rows, shortlist_result.shortlist_rows


def _mean(values: Iterable[float]) -> float:
    series = tuple(values)
    if not series:
        raise ValueError("Expected at least one value")
    return sum(series) / len(series)


def _generated_design_rank_key(
    design: LigandMPNNRound2GeneratedDesign,
) -> tuple[float, float, float, int, str]:
    return (
        -design.ligand_confidence,
        -design.overall_confidence,
        -design.seq_recovery,
        design.design_id,
        design.sequence_id,
    )


def _candidate_unique_sequence_rank_key(
    row: CandidateUniqueRound2Sequence,
) -> tuple[float, float, int, float, int, str, str]:
    return (
        -row.mean_ligand_confidence,
        -row.mean_overall_confidence,
        row.candidate_shortlist_rank,
        -row.best_seq_recovery,
        row.design_id,
        row.sequence_id,
        row.designed_chain_sequence,
    )


def _cross_candidate_representative_key(
    row: CandidateUniqueRound2Sequence,
) -> tuple[float, float, int, float, int, int, str, str]:
    return (
        -row.mean_ligand_confidence,
        -row.mean_overall_confidence,
        row.candidate_shortlist_rank,
        -row.best_seq_recovery,
        row.candidate_unique_rank,
        row.design_id,
        row.sequence_id,
        row.candidate_id,
    )


def _global_unique_sequence_rank_key(
    row: LigandMPNNRound2UniqueSequenceRow,
) -> tuple[float, float, int, float, int, int, str, str]:
    return (
        -row.mean_ligand_confidence,
        -row.mean_overall_confidence,
        row.candidate_shortlist_rank,
        -row.best_seq_recovery,
        row.candidate_unique_rank,
        row.design_id,
        row.sequence_id,
        row.candidate_id,
    )


def _retain_first_available(
    *,
    selected_rows: list[LigandMPNNRound2UniqueSequenceRow],
    selected_sequence_ids: set[str],
    retention_reasons: dict[str, str],
    pool: Sequence[LigandMPNNRound2UniqueSequenceRow],
    reason: str,
    total_limit: int,
) -> None:
    if len(selected_rows) >= total_limit:
        return
    for row in pool:
        if row.sequence_id in selected_sequence_ids:
            continue
        selected_rows.append(row)
        selected_sequence_ids.add(row.sequence_id)
        retention_reasons[row.sequence_id] = reason
        return


def _source_candidate_order_key(row: CandidateUniqueRound2Sequence) -> tuple[int, int, str]:
    return row.seed_rank, row.candidate_shortlist_rank, row.candidate_id


def _build_rosetta_round2_shortlist_payload(
    shortlist_rows: Sequence[LigandMPNNRound2ShortlistRow],
) -> dict[str, Any]:
    return {
        "version": 1,
        "phase": "7E",
        "source_artifacts": {
            "proteinmpnn_round2_shortlist": "results/tables/proteinmpnn_round2_shortlist.csv",
            "ligandmpnn_round2_input_manifest": "results/tables/ligandmpnn_round2_input_manifest.csv",
            "ligandmpnn_round2_redesign_positions": "results/tables/ligandmpnn_round2_redesign_positions.csv",
            "ligandmpnn_round2_smoke_summary": "results/tables/ligandmpnn_round2_smoke_summary.csv",
            "ligandmpnn_round2_sequence_catalog": "results/tables/ligandmpnn_round2_sequence_catalog.csv",
            "ligandmpnn_round2_unique_sequences": "results/tables/ligandmpnn_round2_unique_sequences.csv",
            "ligandmpnn_round2_shortlist": "results/tables/ligandmpnn_round2_shortlist.csv",
        },
        "shortlist_rules": {
            "total_limit": ROUND2_SHORTLIST_TOTAL_LIMIT,
            "required_categories": {
                f"{AM1_CAMPAIGN_PREFIX}*": 1,
                POCKET_CAMPAIGN_ID: 1,
                INTERFACE_CAMPAIGN_ID: 1,
            },
            "wildcard_remaining": 1,
            "ranking": [
                "mean_ligand_confidence_desc",
                "mean_overall_confidence_desc",
                "candidate_shortlist_rank_asc",
            ],
        },
        "designs": [
            {
                "shortlist_rank": row.shortlist_rank,
                "sequence_id": row.sequence_id,
                "candidate_id": row.candidate_id,
                "seed_rank": row.seed_rank,
                "seed_candidate_id": row.seed_candidate_id,
                "campaign_id": row.campaign_id,
                "backbone_id": row.backbone_id,
                "topology_class": row.topology_class,
                "preserved_metal_identity": row.preserved_metal_identity,
                "candidate_shortlist_rank": row.candidate_shortlist_rank,
                "candidate_unique_rank": row.candidate_unique_rank,
                "design_id": row.design_id,
                "total_occurrence_count": row.total_occurrence_count,
                "source_candidate_count": row.source_candidate_count,
                "source_candidate_ids": row.source_candidate_ids.split(",") if row.source_candidate_ids else [],
                "mean_ligand_confidence": row.mean_ligand_confidence,
                "mean_overall_confidence": row.mean_overall_confidence,
                "best_seq_recovery": row.best_seq_recovery,
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
