"""Phase 5A1 deterministic Rosetta score_jd2 baseline scoring."""

from __future__ import annotations

import csv
import logging
import shlex
import shutil
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_text, write_csv_rows
from lanm.models import LigandMPNNSequenceCatalogRow, LigandMPNNShortlistRow
from lanm.paths import REPO_ROOT

LOGGER = logging.getLogger(__name__)

ROSETTA_SEED = 37
ROSETTA_NSTRUCT = 1
ROSETTA_FA_ELEC_WEIGHT = 0.0
SCORE_APPLICATION_NAME = "score_jd2"
_ROSETTA_EXECUTABLE_SUFFIXES = (
    ".default.linuxgccrelease",
    ".linuxgccrelease",
    ".default.linuxclangrelease",
    ".linuxclangrelease",
    ".default.macosclangrelease",
    ".macosclangrelease",
)


@dataclass(frozen=True, slots=True)
class RosettaScoreSmokeCandidate:
    shortlist_rank: int
    candidate_id: str
    campaign_id: str
    backbone_id: str
    preserved_metal_identity: str
    representative_design_id: int
    representative_ligand_confidence: float
    representative_overall_confidence: float
    representative_packed_pdb: Path
    output_dir: Path


@dataclass(frozen=True, slots=True)
class RosettaScoreSmokeResult:
    candidate: RosettaScoreSmokeCandidate
    scores: dict[str, float | str]
    scorefile_path: Path

    @property
    def total_score(self) -> float:
        value = self.scores.get("total_score")
        if not isinstance(value, float):
            raise ValueError(f"Missing numeric total_score in {self.scorefile_path}")
        return value


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


def validate_rosetta_paths(rosetta_bin_dir: Path, rosetta_database: Path) -> tuple[Path, Path]:
    """Validate the Rosetta bin/database locations and return normalized paths."""
    normalized_bin_dir = rosetta_bin_dir.expanduser()
    normalized_database = rosetta_database.expanduser()
    if not normalized_bin_dir.is_absolute():
        raise ValueError("--rosetta-bin-dir must be an absolute path")
    if not normalized_database.is_absolute():
        raise ValueError("--rosetta-database must be an absolute path")
    require_path(normalized_bin_dir)
    require_path(normalized_database)
    return normalized_bin_dir, normalized_database


def resolve_rosetta_executable(rosetta_bin_dir: Path, application_name: str) -> Path:
    """Resolve a Rosetta executable deterministically from a bin directory."""
    require_path(rosetta_bin_dir)
    for suffix in _ROSETTA_EXECUTABLE_SUFFIXES:
        candidate = rosetta_bin_dir / f"{application_name}{suffix}"
        if candidate.exists():
            return candidate
    fallback_matches = sorted(path for path in rosetta_bin_dir.glob(f"{application_name}*") if path.is_file())
    if fallback_matches:
        return fallback_matches[0]
    raise FileNotFoundError(str(rosetta_bin_dir / f"{application_name}{_ROSETTA_EXECUTABLE_SUFFIXES[0]}"))


def _load_ligandmpnn_shortlist_rows(path: Path) -> tuple[LigandMPNNShortlistRow, ...]:
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


def _load_ligandmpnn_sequence_catalog_rows(path: Path) -> tuple[LigandMPNNSequenceCatalogRow, ...]:
    require_path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [
            LigandMPNNSequenceCatalogRow(
                shortlist_rank=int(str(row["shortlist_rank"]).strip()),
                candidate_id=str(row["candidate_id"]).strip(),
                campaign_id=str(row["campaign_id"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                source_structure_id=str(row["source_structure_id"]).strip(),
                designed_chains=str(row["designed_chains"]).strip(),
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


def choose_representative_packed_pdb(
    sequence_rows: Sequence[LigandMPNNSequenceCatalogRow],
) -> LigandMPNNSequenceCatalogRow:
    """Choose the deterministic representative packed PDB for one shortlisted candidate."""
    if not sequence_rows:
        raise ValueError("Cannot choose a representative packed PDB from an empty sequence catalog")
    return min(
        sequence_rows,
        key=lambda row: (
            -row.ligand_confidence,
            -row.overall_confidence,
            row.design_id,
            row.packed_pdb_path,
        ),
    )


def discover_rosetta_score_smoke_candidates(
    *,
    shortlist_path: Path,
    sequence_catalog_path: Path,
    output_root: Path,
) -> tuple[RosettaScoreSmokeCandidate, ...]:
    """Resolve shortlisted candidates and select one representative packed PDB for each."""
    shortlist_rows = _load_ligandmpnn_shortlist_rows(shortlist_path)
    sequence_rows = _load_ligandmpnn_sequence_catalog_rows(sequence_catalog_path)
    shortlist_by_candidate: dict[str, LigandMPNNShortlistRow] = {}
    for row in shortlist_rows:
        existing = shortlist_by_candidate.get(row.candidate_id)
        if existing is None or row.shortlist_rank < existing.shortlist_rank:
            shortlist_by_candidate[row.candidate_id] = row

    sequence_rows_by_candidate: dict[str, list[LigandMPNNSequenceCatalogRow]] = defaultdict(list)
    for row in sequence_rows:
        sequence_rows_by_candidate[row.candidate_id].append(row)

    candidates: list[RosettaScoreSmokeCandidate] = []
    for shortlist_row in sorted(shortlist_by_candidate.values(), key=lambda row: row.shortlist_rank):
        candidate_rows = sequence_rows_by_candidate.get(shortlist_row.candidate_id)
        if not candidate_rows:
            raise ValueError(f"No LigandMPNN sequence catalog rows found for {shortlist_row.candidate_id}")
        representative_row = choose_representative_packed_pdb(candidate_rows)
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
        if representative_row.preserved_metal_identity != shortlist_row.preserved_metal_identity:
            raise ValueError(
                f"Preserved metal mismatch for {shortlist_row.candidate_id}: "
                f"{representative_row.preserved_metal_identity} vs {shortlist_row.preserved_metal_identity}"
            )
        representative_packed_pdb = _repo_path(representative_row.packed_pdb_path)
        require_path(representative_packed_pdb)
        candidates.append(
            RosettaScoreSmokeCandidate(
                shortlist_rank=shortlist_row.shortlist_rank,
                candidate_id=shortlist_row.candidate_id,
                campaign_id=shortlist_row.campaign_id,
                backbone_id=shortlist_row.backbone_id,
                preserved_metal_identity=shortlist_row.preserved_metal_identity,
                representative_design_id=representative_row.design_id,
                representative_ligand_confidence=representative_row.ligand_confidence,
                representative_overall_confidence=representative_row.overall_confidence,
                representative_packed_pdb=representative_packed_pdb,
                output_dir=output_root / shortlist_row.candidate_id,
            )
        )
    return tuple(candidates)


def build_rosetta_score_command(
    *,
    rosetta_executable: Path,
    rosetta_database: Path,
    input_pdb: Path,
    scorefile_path: Path,
) -> tuple[str, ...]:
    """Build the deterministic Rosetta score_jd2 baseline command."""
    return (
        str(rosetta_executable),
        "-database",
        str(rosetta_database),
        "-in:file:fullatom",
        "-s",
        str(input_pdb),
        "-out:file:scorefile",
        str(scorefile_path),
        "-out:overwrite",
        "-nstruct",
        str(ROSETTA_NSTRUCT),
        "-ignore_waters",
        "false",
        "-in:auto_setup_metals",
        "true",
        "-run:preserve_header",
        "true",
        "-out:file:renumber_pdb",
        "false",
        "-run:constant_seed",
        "true",
        "-run:jran",
        str(ROSETTA_SEED),
        "-score:set_weights",
        "fa_elec",
        str(ROSETTA_FA_ELEC_WEIGHT),
    )


def parse_rosetta_scorefile(path: Path) -> dict[str, float | str]:
    """Parse a Rosetta text scorefile that contains exactly one SCORE row."""
    require_path(path)
    header: list[str] | None = None
    parsed_rows: list[dict[str, float | str]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped == "SEQUENCE:" or not stripped.startswith("SCORE:"):
            continue
        tokens = stripped.split()
        if header is None:
            header = tokens[1:]
            continue
        values = tokens[1:]
        if len(values) != len(header):
            raise ValueError(f"Scorefile row length mismatch in {path}: {len(values)} vs {len(header)}")
        parsed_rows.append({key: _coerce_score_value(value) for key, value in zip(header, values)})
    if header is None:
        raise ValueError(f"No SCORE header found in {path}")
    if len(parsed_rows) != 1:
        raise ValueError(f"Expected exactly one SCORE row in {path}, found {len(parsed_rows)}")
    return parsed_rows[0]


def rank_rosetta_score_candidates(
    results: Sequence[RosettaScoreSmokeResult],
) -> tuple[RosettaScoreSmokeResult, ...]:
    """Rank candidates by total_score with deterministic candidate-id tie-breaking."""
    return tuple(sorted(results, key=lambda result: (result.total_score, result.candidate.candidate_id)))


def build_rosetta_score_summary_rows(
    results: Sequence[RosettaScoreSmokeResult],
) -> tuple[dict[str, object], ...]:
    """Build deterministic baseline scoring summary rows for downstream ranking."""
    ordered_score_terms = _ordered_score_terms(results)
    rows: list[dict[str, object]] = []
    for result in sorted(results, key=lambda item: item.candidate.shortlist_rank):
        candidate = result.candidate
        row: dict[str, object] = {
            "shortlist_rank": candidate.shortlist_rank,
            "candidate_id": candidate.candidate_id,
            "campaign_id": candidate.campaign_id,
            "backbone_id": candidate.backbone_id,
            "preserved_metal_identity": candidate.preserved_metal_identity,
            "representative_design_id": candidate.representative_design_id,
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


def build_rosetta_score_ranking_rows(
    results: Sequence[RosettaScoreSmokeResult],
) -> tuple[dict[str, object], ...]:
    """Build total_score-ranked candidate rows for baseline Rosetta comparison."""
    summary_by_candidate = {
        str(row["candidate_id"]): row
        for row in build_rosetta_score_summary_rows(results)
    }
    ranking_rows: list[dict[str, object]] = []
    for rank, result in enumerate(rank_rosetta_score_candidates(results), start=1):
        ranking_row = {"rosetta_rank": rank}
        ranking_row.update(summary_by_candidate[result.candidate.candidate_id])
        ranking_rows.append(ranking_row)
    return tuple(ranking_rows)


def render_rosetta_score_smoke_markdown(
    *,
    rosetta_bin_dir: Path,
    rosetta_database: Path,
    ranked_results: Sequence[RosettaScoreSmokeResult],
) -> str:
    """Render the Phase 5A1 baseline Rosetta scoring report."""
    lines = [
        "# Rosetta Score Smoke",
        "",
        "Deterministic Phase 5A1 Rosetta `score_jd2` baseline scoring for the LigandMPNN shortlist.",
        "",
        f"- rosetta_bin_dir: `{rosetta_bin_dir}`",
        f"- rosetta_database: `{rosetta_database}`",
        f"- score application: `{SCORE_APPLICATION_NAME}`",
        f"- seed: `{ROSETTA_SEED}`",
        f"- nstruct: `{ROSETTA_NSTRUCT}`",
        f"- fa_elec weight override: `{ROSETTA_FA_ELEC_WEIGHT}`",
        "- representative selection: `highest ligand_confidence`, then `highest overall_confidence`, then `smallest design_id`",
        "- preserved input handling: original packed PDB scored directly with waters retained, auto metal setup enabled, headers preserved, and PDB renumbering disabled",
        "- excluded in this phase: `relax`, `fixbb`, `MD`, `QM`, and quantum steps",
        "",
        "| rosetta_rank | candidate_id | campaign_id | backbone_id | metal | representative_design_id | total_score |",
        "| ---: | --- | --- | --- | --- | ---: | ---: |",
    ]
    for rank, result in enumerate(ranked_results, start=1):
        candidate = result.candidate
        lines.append(
            "| "
            f"{rank} | "
            f"{candidate.candidate_id} | "
            f"{candidate.campaign_id} | "
            f"{candidate.backbone_id} | "
            f"{candidate.preserved_metal_identity} | "
            f"{candidate.representative_design_id} | "
            f"{result.total_score:.3f} |"
        )
    lines.append("")
    return "\n".join(lines)


def run_rosetta_score_smoke(
    *,
    rosetta_bin_dir: Path,
    rosetta_database: Path,
    shortlist_path: Path,
    sequence_catalog_path: Path,
    output_root: Path,
    summary_path: Path,
    ranking_path: Path,
    report_path: Path,
) -> tuple[tuple[dict[str, object], ...], tuple[dict[str, object], ...]]:
    """Run deterministic Rosetta score_jd2 baseline scoring for the LigandMPNN shortlist."""
    validated_bin_dir, validated_database = validate_rosetta_paths(rosetta_bin_dir, rosetta_database)
    score_executable = resolve_rosetta_executable(validated_bin_dir, SCORE_APPLICATION_NAME)
    candidates = discover_rosetta_score_smoke_candidates(
        shortlist_path=shortlist_path,
        sequence_catalog_path=sequence_catalog_path,
        output_root=output_root,
    )
    if not candidates:
        raise ValueError(f"No Rosetta score smoke candidates discovered from {shortlist_path}")

    _reset_output_directory(output_root)
    results: list[RosettaScoreSmokeResult] = []
    for candidate in candidates:
        LOGGER.info("Running Rosetta score_jd2 baseline for %s", candidate.candidate_id)
        candidate.output_dir.mkdir(parents=True, exist_ok=True)
        scorefile_path = candidate.output_dir / "score.sc"
        score_command = build_rosetta_score_command(
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
            RosettaScoreSmokeResult(
                candidate=candidate,
                scores=parse_rosetta_scorefile(scorefile_path),
                scorefile_path=scorefile_path,
            )
        )

    summary_rows = build_rosetta_score_summary_rows(results)
    ranking_rows = build_rosetta_score_ranking_rows(results)
    write_csv_rows(summary_path, summary_rows)
    write_csv_rows(ranking_path, ranking_rows)
    atomic_write_text(
        report_path,
        render_rosetta_score_smoke_markdown(
            rosetta_bin_dir=validated_bin_dir,
            rosetta_database=validated_database,
            ranked_results=rank_rosetta_score_candidates(results),
        ),
    )
    return summary_rows, ranking_rows


def _ordered_score_terms(results: Sequence[RosettaScoreSmokeResult]) -> tuple[str, ...]:
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


def _coerce_score_value(value: str) -> float | str:
    try:
        return float(value)
    except ValueError:
        return value


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
