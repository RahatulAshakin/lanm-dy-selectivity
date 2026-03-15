"""Phase 7F deterministic Rosetta score-only triage for the round-2 shortlist."""

from __future__ import annotations

import csv
import logging
import shlex
import shutil
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

from lanm.analysis.rosetta_score_smoke import (
    ROSETTA_FA_ELEC_WEIGHT,
    ROSETTA_NSTRUCT,
    ROSETTA_SEED,
    SCORE_APPLICATION_NAME,
    build_rosetta_score_command as build_shared_rosetta_score_command,
    parse_rosetta_scorefile as parse_shared_rosetta_scorefile,
    resolve_rosetta_executable,
    validate_rosetta_paths,
)
from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_text, write_csv_rows
from lanm.models import LigandMPNNRound2SequenceCatalogRow, LigandMPNNRound2ShortlistRow
from lanm.paths import REPO_ROOT

LOGGER = logging.getLogger(__name__)

AM1_CAMPAIGN_PREFIX = "am1_mex"
HANS_POCKET_CAMPAIGN_ID = "hans_pocket_ss_only"
HANS_INTERFACE_TOPOLOGY_CLASS = "hans_interface_multichain"
TOPOLOGY_CLASS_ORDER = {
    "am1_monomer": 1,
    "hans_monomer": 2,
    "hans_interface_multichain": 3,
}
ROUND2_CROSS_TOPOLOGY_COMPARISON_NOTE = (
    "Raw Rosetta total_score is not a valid cross-topology comparator; Phase 7F ranks candidates only within "
    "am1_monomer, hans_monomer, and hans_interface_multichain topology classes."
)


@dataclass(frozen=True, slots=True)
class RosettaRound2ScoreSmokeCandidate:
    shortlist_rank: int
    shortlist_sequence_id: str
    candidate_id: str
    seed_rank: int
    seed_candidate_id: str
    campaign_id: str
    backbone_id: str
    topology_class: str
    preserved_metal_identity: str
    shortlist_design_id: int
    shortlist_retention_reason: str
    representative_design_id: int
    representative_sequence_id: str
    representative_ligand_confidence: float
    representative_overall_confidence: float
    representative_packed_pdb: Path
    output_dir: Path


@dataclass(frozen=True, slots=True)
class RosettaRound2ScoreSmokeResult:
    candidate: RosettaRound2ScoreSmokeCandidate
    scores: Mapping[str, float | str]
    scorefile_path: Path

    @property
    def total_score(self) -> float:
        value = self.scores.get("total_score")
        if not isinstance(value, float):
            raise ValueError(f"Missing numeric total_score in {self.scorefile_path}")
        return value


@dataclass(frozen=True, slots=True)
class Round2MDRescreenPanelRow:
    panel_rank: int
    panel_member_id: str
    panel_role: str
    candidate_id: str
    campaign_id: str
    backbone_id: str
    topology_class: str
    preserved_metal_identity: str
    rosetta_score_rank_within_topology: int
    topology_candidate_count: int
    representative_design_id: int
    representative_sequence_id: str
    representative_ligand_confidence: float
    representative_overall_confidence: float
    total_score: float
    representative_packed_pdb: str
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


def _load_round2_shortlist_rows(path: Path) -> tuple[LigandMPNNRound2ShortlistRow, ...]:
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            LigandMPNNRound2ShortlistRow(
                shortlist_rank=int(str(row["shortlist_rank"]).strip()),
                sequence_id=str(row["sequence_id"]).strip(),
                candidate_id=str(row["candidate_id"]).strip(),
                seed_rank=int(str(row["seed_rank"]).strip()),
                seed_candidate_id=str(row["seed_candidate_id"]).strip(),
                campaign_id=str(row["campaign_id"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                topology_class=str(row["topology_class"]).strip(),
                preserved_metal_identity=str(row["preserved_metal_identity"]).strip(),
                candidate_shortlist_rank=int(str(row["candidate_shortlist_rank"]).strip()),
                candidate_unique_rank=int(str(row["candidate_unique_rank"]).strip()),
                design_id=int(str(row["design_id"]).strip()),
                total_occurrence_count=int(str(row["total_occurrence_count"]).strip()),
                source_candidate_count=int(str(row["source_candidate_count"]).strip()),
                source_candidate_ids=str(row["source_candidate_ids"]).strip(),
                mean_ligand_confidence=float(str(row["mean_ligand_confidence"]).strip()),
                mean_overall_confidence=float(str(row["mean_overall_confidence"]).strip()),
                best_seq_recovery=float(str(row["best_seq_recovery"]).strip()),
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


def _load_round2_sequence_catalog_rows(path: Path) -> tuple[LigandMPNNRound2SequenceCatalogRow, ...]:
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            LigandMPNNRound2SequenceCatalogRow(
                shortlist_rank=int(str(row["shortlist_rank"]).strip()),
                candidate_id=str(row["candidate_id"]).strip(),
                seed_rank=int(str(row["seed_rank"]).strip()),
                seed_candidate_id=str(row["seed_candidate_id"]).strip(),
                campaign_id=str(row["campaign_id"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                topology_class=str(row["topology_class"]).strip(),
                designed_chains=str(row["designed_chains"]).strip(),
                fixed_context_chains=str(row["fixed_context_chains"]).strip(),
                preserved_metal_identity=str(row["preserved_metal_identity"]).strip(),
                redesigned_residues=str(row["redesigned_residues"]).strip(),
                design_id=int(str(row["design_id"]).strip()),
                temperature=float(str(row["temperature"]).strip()),
                seed=int(str(row["seed"]).strip()),
                overall_confidence=float(str(row["overall_confidence"]).strip()),
                ligand_confidence=float(str(row["ligand_confidence"]).strip()),
                seq_recovery=float(str(row["seq_recovery"]).strip()),
                sequence_id=str(row["sequence_id"]).strip(),
                designed_sequence=str(row["designed_sequence"]).strip(),
                sequence_length=int(str(row["sequence_length"]).strip()),
                output_fasta_path=str(row["output_fasta_path"]).strip(),
                backbone_pdb_path=str(row["backbone_pdb_path"]).strip(),
                packed_pdb_path=str(row["packed_pdb_path"]).strip(),
            )
            for row in reader
        ]
    return tuple(rows)


def choose_round2_representative_packed_pdb(
    sequence_rows: Sequence[LigandMPNNRound2SequenceCatalogRow],
) -> LigandMPNNRound2SequenceCatalogRow:
    """Choose the deterministic representative packed PDB for one round-2 candidate."""
    if not sequence_rows:
        raise ValueError("Cannot choose a representative packed PDB from an empty round-2 sequence catalog")
    return min(
        sequence_rows,
        key=lambda row: (
            -row.ligand_confidence,
            -row.overall_confidence,
            row.design_id,
            row.sequence_id,
            row.packed_pdb_path,
        ),
    )


def discover_rosetta_round2_score_smoke_candidates(
    *,
    shortlist_path: Path,
    sequence_catalog_path: Path,
    output_root: Path,
) -> tuple[RosettaRound2ScoreSmokeCandidate, ...]:
    """Resolve round-2 shortlisted candidates and select one representative packed PDB each."""
    shortlist_rows = _load_round2_shortlist_rows(shortlist_path)
    sequence_rows = _load_round2_sequence_catalog_rows(sequence_catalog_path)
    shortlist_by_candidate: dict[str, LigandMPNNRound2ShortlistRow] = {}
    for row in shortlist_rows:
        existing = shortlist_by_candidate.get(row.candidate_id)
        if existing is None or row.shortlist_rank < existing.shortlist_rank:
            shortlist_by_candidate[row.candidate_id] = row

    sequence_rows_by_candidate: dict[str, list[LigandMPNNRound2SequenceCatalogRow]] = defaultdict(list)
    for row in sequence_rows:
        sequence_rows_by_candidate[row.candidate_id].append(row)

    candidates: list[RosettaRound2ScoreSmokeCandidate] = []
    for shortlist_row in sorted(shortlist_by_candidate.values(), key=lambda row: row.shortlist_rank):
        candidate_rows = sequence_rows_by_candidate.get(shortlist_row.candidate_id)
        if not candidate_rows:
            raise ValueError(f"No round-2 LigandMPNN sequence catalog rows found for {shortlist_row.candidate_id}")
        representative_row = choose_round2_representative_packed_pdb(candidate_rows)
        if representative_row.campaign_id != shortlist_row.campaign_id:
            raise ValueError(
                f"Campaign mismatch for {shortlist_row.candidate_id}: "
                f"{representative_row.campaign_id} vs {shortlist_row.campaign_id}"
            )
        if representative_row.backbone_id != shortlist_row.backbone_id:
            raise ValueError(
                f"Backbone mismatch for {shortlist_row.candidate_id}: "
                f"{representative_row.backbone_id} vs {shortlist_row.backbone_id}"
            )
        if representative_row.topology_class != shortlist_row.topology_class:
            raise ValueError(
                f"Topology mismatch for {shortlist_row.candidate_id}: "
                f"{representative_row.topology_class} vs {shortlist_row.topology_class}"
            )
        if representative_row.preserved_metal_identity != shortlist_row.preserved_metal_identity:
            raise ValueError(
                f"Preserved metal mismatch for {shortlist_row.candidate_id}: "
                f"{representative_row.preserved_metal_identity} vs {shortlist_row.preserved_metal_identity}"
            )
        representative_packed_pdb = _repo_path(representative_row.packed_pdb_path)
        require_path(representative_packed_pdb)
        candidates.append(
            RosettaRound2ScoreSmokeCandidate(
                shortlist_rank=shortlist_row.shortlist_rank,
                shortlist_sequence_id=shortlist_row.sequence_id,
                candidate_id=shortlist_row.candidate_id,
                seed_rank=shortlist_row.seed_rank,
                seed_candidate_id=shortlist_row.seed_candidate_id,
                campaign_id=shortlist_row.campaign_id,
                backbone_id=shortlist_row.backbone_id,
                topology_class=shortlist_row.topology_class,
                preserved_metal_identity=shortlist_row.preserved_metal_identity,
                shortlist_design_id=shortlist_row.design_id,
                shortlist_retention_reason=shortlist_row.retention_reason,
                representative_design_id=representative_row.design_id,
                representative_sequence_id=representative_row.sequence_id,
                representative_ligand_confidence=representative_row.ligand_confidence,
                representative_overall_confidence=representative_row.overall_confidence,
                representative_packed_pdb=representative_packed_pdb,
                output_dir=output_root / shortlist_row.candidate_id,
            )
        )
    return tuple(candidates)


def build_round2_rosetta_score_command(
    *,
    rosetta_executable: Path,
    rosetta_database: Path,
    input_pdb: Path,
    scorefile_path: Path,
) -> tuple[str, ...]:
    """Build the deterministic Phase 7F Rosetta score_jd2 baseline command."""
    return build_shared_rosetta_score_command(
        rosetta_executable=rosetta_executable,
        rosetta_database=rosetta_database,
        input_pdb=input_pdb,
        scorefile_path=scorefile_path,
    )


def parse_round2_rosetta_scorefile(path: Path) -> dict[str, float | str]:
    """Parse a single-pose Phase 7F Rosetta scorefile."""
    return parse_shared_rosetta_scorefile(path)


def build_rosetta_round2_score_summary_rows(
    results: Sequence[RosettaRound2ScoreSmokeResult],
) -> tuple[dict[str, object], ...]:
    """Build deterministic Phase 7F summary rows."""
    ordered_score_terms = _ordered_score_terms(results)
    rows: list[dict[str, object]] = []
    for result in sorted(results, key=lambda item: item.candidate.shortlist_rank):
        candidate = result.candidate
        row: dict[str, object] = {
            "shortlist_rank": candidate.shortlist_rank,
            "shortlist_sequence_id": candidate.shortlist_sequence_id,
            "candidate_id": candidate.candidate_id,
            "seed_rank": candidate.seed_rank,
            "seed_candidate_id": candidate.seed_candidate_id,
            "campaign_id": candidate.campaign_id,
            "backbone_id": candidate.backbone_id,
            "topology_class": candidate.topology_class,
            "preserved_metal_identity": candidate.preserved_metal_identity,
            "shortlist_design_id": candidate.shortlist_design_id,
            "shortlist_retention_reason": candidate.shortlist_retention_reason,
            "representative_design_id": candidate.representative_design_id,
            "representative_sequence_id": candidate.representative_sequence_id,
            "representative_ligand_confidence": candidate.representative_ligand_confidence,
            "representative_overall_confidence": candidate.representative_overall_confidence,
            "representative_packed_pdb": _display_path(candidate.representative_packed_pdb),
            "scorefile_path": _display_path(result.scorefile_path),
            "total_score": result.total_score,
            "score_description": result.scores.get("description"),
        }
        for term in ordered_score_terms:
            row[term] = result.scores.get(term)
        rows.append(row)
    return tuple(rows)


def build_rosetta_round2_candidate_ranking_rows(
    results: Sequence[RosettaRound2ScoreSmokeResult],
) -> tuple[dict[str, object], ...]:
    """Build within-topology ranking rows for Phase 7F."""
    summary_by_candidate = {
        str(row["candidate_id"]): row
        for row in build_rosetta_round2_score_summary_rows(results)
    }
    grouped_results: dict[str, list[RosettaRound2ScoreSmokeResult]] = defaultdict(list)
    for result in results:
        grouped_results[result.candidate.topology_class].append(result)

    ranking_rows: list[dict[str, object]] = []
    for topology_class in _sorted_topology_classes(grouped_results.keys()):
        ordered_results = sorted(
            grouped_results[topology_class],
            key=lambda item: (
                item.total_score,
                -item.candidate.representative_ligand_confidence,
                -item.candidate.representative_overall_confidence,
                item.candidate.shortlist_rank,
                item.candidate.candidate_id,
            ),
        )
        best_total_score = ordered_results[0].total_score
        previous_total_score: float | None = None
        current_rank = 0
        for position, result in enumerate(ordered_results, start=1):
            if previous_total_score is None or result.total_score != previous_total_score:
                current_rank = position
                previous_total_score = result.total_score
            ranking_row = {
                "rosetta_score_rank_within_topology": current_rank,
                "topology_candidate_count": len(ordered_results),
                "topology_best_total_score": best_total_score,
                "total_score_delta_from_topology_best": result.total_score - best_total_score,
                "cross_topology_total_score_comparison_valid": False,
                "cross_topology_comparison_note": ROUND2_CROSS_TOPOLOGY_COMPARISON_NOTE,
            }
            ranking_row.update(summary_by_candidate[result.candidate.candidate_id])
            ranking_rows.append(ranking_row)
    return tuple(ranking_rows)


def select_round2_md_rescreen_panel(
    ranking_rows: Sequence[Mapping[str, object]],
) -> tuple[Round2MDRescreenPanelRow, ...]:
    """Select the reduced three-candidate Phase 7F MD/metadynamics rescreen panel."""
    ordered_rows = tuple(
        sorted(
            ranking_rows,
            key=lambda row: (
                TOPOLOGY_CLASS_ORDER.get(str(row["topology_class"]), 99),
                int(row["rosetta_score_rank_within_topology"]),
                -float(row["representative_ligand_confidence"]),
                -float(row["representative_overall_confidence"]),
                int(row["shortlist_rank"]),
                str(row["candidate_id"]),
            ),
        )
    )
    am1_rows = [row for row in ordered_rows if str(row["campaign_id"]).startswith(AM1_CAMPAIGN_PREFIX)]
    pocket_rows = [row for row in ordered_rows if str(row["campaign_id"]) == HANS_POCKET_CAMPAIGN_ID]
    interface_rows = [row for row in ordered_rows if str(row["topology_class"]) == HANS_INTERFACE_TOPOLOGY_CLASS]

    selected_rows: list[Round2MDRescreenPanelRow] = []
    panel_rank = 1
    if am1_rows:
        best_am1 = _best_panel_candidate(am1_rows)
        selected_rows.append(
            _build_panel_row(
                panel_rank=panel_rank,
                row=best_am1,
                panel_role="top_am1_mex_round2_candidate",
                selection_reason=(
                    "Retained as the top available AM1/Mex round-2 candidate by within-topology Rosetta score "
                    f"rank={int(best_am1['rosetta_score_rank_within_topology'])} with "
                    f"ligand_confidence={float(best_am1['representative_ligand_confidence']):.4f}."
                ),
            )
        )
        panel_rank += 1

    if pocket_rows:
        best_pocket = _best_panel_candidate(pocket_rows)
        selected_rows.append(
            _build_panel_row(
                panel_rank=panel_rank,
                row=best_pocket,
                panel_role="top_hans_pocket_round2_candidate",
                selection_reason=(
                    "Retained as the top available Hans pocket round-2 candidate by within-topology Rosetta score "
                    f"rank={int(best_pocket['rosetta_score_rank_within_topology'])} with "
                    f"ligand_confidence={float(best_pocket['representative_ligand_confidence']):.4f}."
                ),
            )
        )
        panel_rank += 1

    if interface_rows:
        ordered_interface_rows = sorted(interface_rows, key=_panel_selection_key)
        winner = ordered_interface_rows[0]
        if len(ordered_interface_rows) > 1:
            runner_up = ordered_interface_rows[1]
            selection_reason = (
                "Retained as the better Hans interface-aware round-2 candidate using within-topology Rosetta score "
                "rank first and ligand_confidence second; "
                f"selected {winner['candidate_id']} (rank={int(winner['rosetta_score_rank_within_topology'])}, "
                f"ligand_confidence={float(winner['representative_ligand_confidence']):.4f}) over "
                f"{runner_up['candidate_id']} (rank={int(runner_up['rosetta_score_rank_within_topology'])}, "
                f"ligand_confidence={float(runner_up['representative_ligand_confidence']):.4f})."
            )
        else:
            selection_reason = (
                "Retained as the only available Hans interface-aware round-2 candidate after within-topology "
                "Rosetta ranking."
            )
        selected_rows.append(
            _build_panel_row(
                panel_rank=panel_rank,
                row=winner,
                panel_role="best_hans_interface_aware_round2_candidate",
                selection_reason=selection_reason,
            )
        )

    return tuple(selected_rows)


def render_rosetta_round2_score_smoke_markdown(
    *,
    rosetta_bin_dir: Path,
    rosetta_database: Path,
    ranking_rows: Sequence[Mapping[str, object]],
    panel_rows: Sequence[Round2MDRescreenPanelRow],
) -> str:
    """Render the Phase 7F Rosetta score-only triage report."""
    lines = [
        "# Rosetta Round-2 Score Smoke",
        "",
        "Deterministic Phase 7F Rosetta `score_jd2` baseline scoring for the 4 round-2 shortlisted candidates.",
        "",
        f"- rosetta_bin_dir: `{rosetta_bin_dir}`",
        f"- rosetta_database: `{rosetta_database}`",
        f"- score application: `{SCORE_APPLICATION_NAME}`",
        f"- seed: `{ROSETTA_SEED}`",
        f"- nstruct: `{ROSETTA_NSTRUCT}`",
        f"- fa_elec weight override: `{ROSETTA_FA_ELEC_WEIGHT}`",
        "- representative selection: `highest ligand_confidence`, then `highest overall_confidence`, then `smallest design_id`",
        "- ranking scope: `within topology class only`",
        f"- cross-topology note: {ROUND2_CROSS_TOPOLOGY_COMPARISON_NOTE}",
        "- excluded in this phase: `relax`, `MD`, `QM`, and quantum steps",
        "",
        "## Within-Topology Ranking",
        "",
        "| topology_class | rosetta_rank_within_topology | candidate_id | representative_design_id | ligand_confidence | total_score |",
        "| --- | ---: | --- | ---: | ---: | ---: |",
    ]
    for row in ranking_rows:
        lines.append(
            "| "
            f"{row['topology_class']} | "
            f"{row['rosetta_score_rank_within_topology']} | "
            f"{row['candidate_id']} | "
            f"{row['representative_design_id']} | "
            f"{float(row['representative_ligand_confidence']):.4f} | "
            f"{float(row['total_score']):.3f} |"
        )
    lines.extend(
        [
            "",
            "## Reduced MD Rescreen Panel",
            "",
            "| panel_rank | panel_role | candidate_id | topology_class | rosetta_rank_within_topology | ligand_confidence | total_score |",
            "| ---: | --- | --- | --- | ---: | ---: | ---: |",
        ]
    )
    for row in panel_rows:
        lines.append(
            "| "
            f"{row.panel_rank} | "
            f"{row.panel_role} | "
            f"{row.candidate_id} | "
            f"{row.topology_class} | "
            f"{row.rosetta_score_rank_within_topology} | "
            f"{row.representative_ligand_confidence:.4f} | "
            f"{row.total_score:.3f} |"
        )
    lines.append("")
    return "\n".join(lines)


def run_rosetta_round2_score_smoke(
    *,
    rosetta_bin_dir: Path,
    rosetta_database: Path,
    shortlist_path: Path,
    sequence_catalog_path: Path,
    output_root: Path,
    summary_path: Path,
    ranking_path: Path,
    panel_path: Path,
    report_path: Path,
) -> tuple[tuple[dict[str, object], ...], tuple[dict[str, object], ...], tuple[Round2MDRescreenPanelRow, ...]]:
    """Run deterministic Phase 7F Rosetta score-only triage for the round-2 shortlist."""
    validated_bin_dir, validated_database = validate_rosetta_paths(rosetta_bin_dir, rosetta_database)
    score_executable = resolve_rosetta_executable(validated_bin_dir, SCORE_APPLICATION_NAME)
    candidates = discover_rosetta_round2_score_smoke_candidates(
        shortlist_path=shortlist_path,
        sequence_catalog_path=sequence_catalog_path,
        output_root=output_root,
    )
    if not candidates:
        raise ValueError(f"No Phase 7F candidates discovered from {shortlist_path}")

    _reset_output_directory(output_root)
    results: list[RosettaRound2ScoreSmokeResult] = []
    for candidate in candidates:
        LOGGER.info("Running Phase 7F Rosetta score_jd2 baseline for %s", candidate.candidate_id)
        candidate.output_dir.mkdir(parents=True, exist_ok=True)
        scorefile_path = candidate.output_dir / "score.sc"
        score_command = build_round2_rosetta_score_command(
            rosetta_executable=score_executable,
            rosetta_database=validated_database,
            input_pdb=candidate.representative_packed_pdb,
            scorefile_path=scorefile_path,
        )
        atomic_write_text(candidate.output_dir / "score_command.txt", shlex.join(score_command) + "\n")
        _run_logged_command(
            command=score_command,
            cwd=candidate.output_dir,
            stdout_path=candidate.output_dir / "score_stdout.txt",
            stderr_path=candidate.output_dir / "score_stderr.txt",
        )
        results.append(
            RosettaRound2ScoreSmokeResult(
                candidate=candidate,
                scores=parse_round2_rosetta_scorefile(scorefile_path),
                scorefile_path=scorefile_path,
            )
        )

    summary_rows = build_rosetta_round2_score_summary_rows(results)
    ranking_rows = build_rosetta_round2_candidate_ranking_rows(results)
    panel_rows = select_round2_md_rescreen_panel(ranking_rows)
    write_csv_rows(summary_path, summary_rows)
    write_csv_rows(ranking_path, ranking_rows)
    write_csv_rows(panel_path, panel_rows)
    atomic_write_text(
        report_path,
        render_rosetta_round2_score_smoke_markdown(
            rosetta_bin_dir=validated_bin_dir,
            rosetta_database=validated_database,
            ranking_rows=ranking_rows,
            panel_rows=panel_rows,
        ),
    )
    return summary_rows, ranking_rows, panel_rows


def _build_panel_row(
    *,
    panel_rank: int,
    row: Mapping[str, object],
    panel_role: str,
    selection_reason: str,
) -> Round2MDRescreenPanelRow:
    return Round2MDRescreenPanelRow(
        panel_rank=panel_rank,
        panel_member_id=str(row["candidate_id"]),
        panel_role=panel_role,
        candidate_id=str(row["candidate_id"]),
        campaign_id=str(row["campaign_id"]),
        backbone_id=str(row["backbone_id"]),
        topology_class=str(row["topology_class"]),
        preserved_metal_identity=str(row["preserved_metal_identity"]),
        rosetta_score_rank_within_topology=int(row["rosetta_score_rank_within_topology"]),
        topology_candidate_count=int(row["topology_candidate_count"]),
        representative_design_id=int(row["representative_design_id"]),
        representative_sequence_id=str(row["representative_sequence_id"]),
        representative_ligand_confidence=float(row["representative_ligand_confidence"]),
        representative_overall_confidence=float(row["representative_overall_confidence"]),
        total_score=float(row["total_score"]),
        representative_packed_pdb=str(row["representative_packed_pdb"]),
        selection_reason=selection_reason,
    )


def _best_panel_candidate(rows: Sequence[Mapping[str, object]]) -> Mapping[str, object]:
    return min(rows, key=_panel_selection_key)


def _ordered_score_terms(results: Sequence[RosettaRound2ScoreSmokeResult]) -> tuple[str, ...]:
    ordered_terms: list[str] = []
    seen_terms: set[str] = set()
    for result in results:
        for key, value in result.scores.items():
            if key in {"description", "total_score"} or key in seen_terms:
                continue
            if isinstance(value, float):
                seen_terms.add(key)
                ordered_terms.append(key)
    return tuple(ordered_terms)


def _panel_selection_key(row: Mapping[str, object]) -> tuple[object, ...]:
    return (
        int(row["rosetta_score_rank_within_topology"]),
        -float(row["representative_ligand_confidence"]),
        -float(row["representative_overall_confidence"]),
        int(row["shortlist_rank"]),
        str(row["candidate_id"]),
    )


def _reset_output_directory(output_dir: Path) -> None:
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)


def _run_logged_command(
    *,
    command: Sequence[str],
    cwd: Path,
    stdout_path: Path,
    stderr_path: Path,
) -> None:
    completed = subprocess.run(
        list(command),
        cwd=cwd,
        check=False,
        capture_output=True,
        text=True,
    )
    atomic_write_text(stdout_path, completed.stdout)
    atomic_write_text(stderr_path, completed.stderr)
    if completed.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {completed.returncode}: {shlex.join(command)}")


def _sorted_topology_classes(topology_classes: Sequence[str]) -> tuple[str, ...]:
    return tuple(sorted(topology_classes, key=lambda value: (TOPOLOGY_CLASS_ORDER.get(value, 99), value)))
