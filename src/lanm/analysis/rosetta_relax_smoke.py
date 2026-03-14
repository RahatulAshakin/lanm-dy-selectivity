"""Phase 5A2 deterministic restrained Rosetta relax smoke refinement."""

from __future__ import annotations

import csv
import logging
import os
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Iterable, Mapping, Sequence

from lanm.analysis.rosetta_score_smoke import (
    ROSETTA_FA_ELEC_WEIGHT,
    ROSETTA_NSTRUCT,
    ROSETTA_SEED,
    parse_rosetta_scorefile,
    resolve_rosetta_executable,
    validate_rosetta_paths,
)
from lanm.data.fetch import require_path
from lanm.filesystem import atomic_write_text, write_csv_rows
from lanm.paths import REPO_ROOT

LOGGER = logging.getLogger(__name__)

ROSETTA_RELAX_APPLICATION_NAME = "relax"
TOP_RELAX_CANDIDATE_COUNT = 3
SINGLE_THREADED_ENVIRONMENT = {
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
}


@dataclass(frozen=True, slots=True)
class RosettaRelaxSmokeCandidate:
    rosetta_rank: int
    shortlist_rank: int
    candidate_id: str
    campaign_id: str
    backbone_id: str
    preserved_metal_identity: str
    representative_design_id: int
    representative_ligand_confidence: float
    representative_overall_confidence: float
    representative_packed_pdb: Path
    pre_relax_scorefile_source: Path
    pre_relax_scores: Mapping[str, float | str]
    output_dir: Path


@dataclass(frozen=True, slots=True)
class RosettaRelaxSmokeResult:
    candidate: RosettaRelaxSmokeCandidate
    pre_relax_scorefile_copy: Path
    relax_scorefile_path: Path
    relaxed_pdb_path: Path | None
    post_relax_scores: Mapping[str, float | str] | None
    success_status: str
    failure_status: str
    failure_reason: str

    @property
    def pre_relax_total_score(self) -> float:
        value = self.candidate.pre_relax_scores.get("total_score")
        if not isinstance(value, float):
            raise ValueError(f"Missing numeric pre-relax total_score for {self.candidate.candidate_id}")
        return value

    @property
    def post_relax_total_score(self) -> float | None:
        if self.post_relax_scores is None:
            return None
        value = self.post_relax_scores.get("total_score")
        if not isinstance(value, float):
            raise ValueError(f"Missing numeric post-relax total_score for {self.candidate.candidate_id}")
        return value

    @property
    def delta_total_score(self) -> float | None:
        post_relax_total_score = self.post_relax_total_score
        if post_relax_total_score is None:
            return None
        return post_relax_total_score - self.pre_relax_total_score

    @property
    def is_success(self) -> bool:
        return self.success_status == "success" and self.post_relax_scores is not None and self.relaxed_pdb_path is not None


def _display_path(path: Path | None) -> str:
    if path is None:
        return ""
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def _repo_path(path_text: str) -> Path:
    path = Path(path_text).expanduser()
    if path.is_absolute():
        return path
    return REPO_ROOT / path


def _coerce_score_value(value: str) -> float | str:
    try:
        return float(value)
    except ValueError:
        return value


def _parse_pre_relax_scores(row: Mapping[str, str]) -> dict[str, float | str]:
    metadata_fields = {
        "rosetta_rank",
        "shortlist_rank",
        "candidate_id",
        "campaign_id",
        "backbone_id",
        "preserved_metal_identity",
        "representative_design_id",
        "representative_ligand_confidence",
        "representative_overall_confidence",
        "representative_packed_pdb",
        "scorefile_path",
        "score_description",
    }
    scores: dict[str, float | str] = {
        "description": str(row["score_description"]).strip(),
        "total_score": float(str(row["total_score"]).strip()),
    }
    for key, raw_value in row.items():
        if key in metadata_fields or key == "total_score":
            continue
        value = str(raw_value).strip()
        if not value:
            continue
        scores[key] = _coerce_score_value(value)
    return scores


def select_top_rosetta_relax_candidates(
    *,
    ranking_path: Path,
    output_root: Path,
    top_n: int = TOP_RELAX_CANDIDATE_COUNT,
) -> tuple[RosettaRelaxSmokeCandidate, ...]:
    """Select the top score-ranked Phase 5A1 candidates for Phase 5A2 relax."""
    require_path(ranking_path)
    if top_n <= 0:
        raise ValueError("top_n must be positive")

    with ranking_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        ordered_rows = sorted(reader, key=lambda row: int(str(row["rosetta_rank"]).strip()))

    selected_rows = ordered_rows[:top_n]
    candidates: list[RosettaRelaxSmokeCandidate] = []
    for row in selected_rows:
        representative_packed_pdb = _repo_path(str(row["representative_packed_pdb"]).strip())
        pre_relax_scorefile_source = _repo_path(str(row["scorefile_path"]).strip())
        require_path(representative_packed_pdb)
        require_path(pre_relax_scorefile_source)
        candidates.append(
            RosettaRelaxSmokeCandidate(
                rosetta_rank=int(str(row["rosetta_rank"]).strip()),
                shortlist_rank=int(str(row["shortlist_rank"]).strip()),
                candidate_id=str(row["candidate_id"]).strip(),
                campaign_id=str(row["campaign_id"]).strip(),
                backbone_id=str(row["backbone_id"]).strip(),
                preserved_metal_identity=str(row["preserved_metal_identity"]).strip(),
                representative_design_id=int(str(row["representative_design_id"]).strip()),
                representative_ligand_confidence=float(str(row["representative_ligand_confidence"]).strip()),
                representative_overall_confidence=float(str(row["representative_overall_confidence"]).strip()),
                representative_packed_pdb=representative_packed_pdb,
                pre_relax_scorefile_source=pre_relax_scorefile_source,
                pre_relax_scores=_parse_pre_relax_scores(row),
                output_dir=output_root / str(row["candidate_id"]).strip(),
            )
        )
    return tuple(candidates)


def build_rosetta_relax_command(
    *,
    rosetta_executable: Path,
    rosetta_database: Path,
    input_pdb: Path,
    output_dir: Path,
    scorefile_path: Path,
) -> tuple[str, ...]:
    """Build the deterministic restrained Rosetta relax smoke command."""
    return (
        str(rosetta_executable),
        "-database",
        str(rosetta_database),
        "-in:file:fullatom",
        "-s",
        str(input_pdb),
        "-out:path:all",
        str(output_dir),
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
        "-relax:constrain_relax_to_start_coords",
        "-relax:ramp_constraints",
        "false",
        "-score:set_weights",
        "fa_elec",
        str(ROSETTA_FA_ELEC_WEIGHT),
    )


def parse_rosetta_relax_scorefile(path: Path) -> dict[str, float | str]:
    """Parse a single-pose Rosetta relax scorefile."""
    return parse_rosetta_scorefile(path)


def resolve_relaxed_pdb_path(
    *,
    relax_output_dir: Path,
    post_relax_scores: Mapping[str, float | str],
) -> Path:
    """Resolve the single relaxed PDB produced by the smoke relax run."""
    require_path(relax_output_dir)
    description = post_relax_scores.get("description")
    if isinstance(description, str):
        expected_path = relax_output_dir / f"{description}.pdb"
        if expected_path.exists():
            return expected_path
    pdb_paths = sorted(relax_output_dir.glob("*.pdb"))
    if len(pdb_paths) == 1:
        return pdb_paths[0]
    if not pdb_paths:
        raise FileNotFoundError(f"No relaxed PDB found in {relax_output_dir}")
    raise ValueError(f"Expected exactly one relaxed PDB in {relax_output_dir}, found {len(pdb_paths)}")


def rank_rosetta_relax_candidates(
    results: Sequence[RosettaRelaxSmokeResult],
) -> tuple[RosettaRelaxSmokeResult, ...]:
    """Rank successful relax results by post-relax total_score, then delta_total_score."""
    successful_results = [result for result in results if result.is_success]
    return tuple(
        sorted(
            successful_results,
            key=lambda result: (
                _require_numeric_score(result.post_relax_total_score, result.candidate.candidate_id, "post_relax_total_score"),
                _require_numeric_score(result.delta_total_score, result.candidate.candidate_id, "delta_total_score"),
                result.candidate.candidate_id,
            ),
        )
    )


def build_rosetta_relax_summary_rows(
    results: Sequence[RosettaRelaxSmokeResult],
) -> tuple[dict[str, object], ...]:
    """Build deterministic Phase 5A2 summary rows for all selected candidates."""
    pre_relax_score_terms = _ordered_score_terms(
        result.candidate.pre_relax_scores for result in results
    )
    post_relax_score_terms = _ordered_score_terms(
        result.post_relax_scores for result in results if result.post_relax_scores is not None
    )
    rows: list[dict[str, object]] = []
    for result in sorted(results, key=lambda item: item.candidate.rosetta_rank):
        candidate = result.candidate
        row: dict[str, object] = {
            "rosetta_rank": candidate.rosetta_rank,
            "shortlist_rank": candidate.shortlist_rank,
            "candidate_id": candidate.candidate_id,
            "campaign_id": candidate.campaign_id,
            "backbone_id": candidate.backbone_id,
            "preserved_metal_identity": candidate.preserved_metal_identity,
            "representative_design_id": candidate.representative_design_id,
            "representative_ligand_confidence": candidate.representative_ligand_confidence,
            "representative_overall_confidence": candidate.representative_overall_confidence,
            "representative_packed_pdb": _display_path(candidate.representative_packed_pdb),
            "pre_relax_scorefile_source": _display_path(candidate.pre_relax_scorefile_source),
            "copied_pre_relax_scorefile": _display_path(result.pre_relax_scorefile_copy),
            "relax_scorefile_path": _display_path(result.relax_scorefile_path),
            "relaxed_pdb_path": _display_path(result.relaxed_pdb_path),
            "pre_relax_total_score": result.pre_relax_total_score,
            "post_relax_total_score": result.post_relax_total_score,
            "delta_total_score": result.delta_total_score,
            "success_status": result.success_status,
            "failure_status": result.failure_status,
            "failure_reason": result.failure_reason,
            "pre_relax_description": candidate.pre_relax_scores.get("description", ""),
            "post_relax_description": (
                result.post_relax_scores.get("description", "") if result.post_relax_scores is not None else ""
            ),
        }
        for term in pre_relax_score_terms:
            row[f"pre_relax_{term}"] = candidate.pre_relax_scores.get(term)
        for term in post_relax_score_terms:
            row[f"post_relax_{term}"] = (
                result.post_relax_scores.get(term) if result.post_relax_scores is not None else None
            )
        rows.append(row)
    return tuple(rows)


def build_rosetta_relax_ranking_rows(
    results: Sequence[RosettaRelaxSmokeResult],
) -> tuple[dict[str, object], ...]:
    """Build successful post-relax ranking rows for MD triage."""
    summary_by_candidate = {
        str(row["candidate_id"]): row
        for row in build_rosetta_relax_summary_rows(results)
    }
    ranking_rows: list[dict[str, object]] = []
    for rank, result in enumerate(rank_rosetta_relax_candidates(results), start=1):
        ranking_row = {"rosetta_relax_rank": rank}
        ranking_row.update(summary_by_candidate[result.candidate.candidate_id])
        ranking_rows.append(ranking_row)
    return tuple(ranking_rows)


def render_rosetta_relax_smoke_markdown(
    *,
    rosetta_bin_dir: Path,
    rosetta_database: Path,
    results: Sequence[RosettaRelaxSmokeResult],
) -> str:
    """Render the Phase 5A2 restrained relax smoke report."""
    ranked_results = rank_rosetta_relax_candidates(results)
    selected_ids = ", ".join(result.candidate.candidate_id for result in sorted(results, key=lambda item: item.candidate.rosetta_rank))
    lines = [
        "# Rosetta Relax Smoke",
        "",
        "Deterministic Phase 5A2 restrained Rosetta `relax` smoke refinement for the top Phase 5A1 score-ranked shortlist candidates.",
        "",
        f"- rosetta_bin_dir: `{rosetta_bin_dir}`",
        f"- rosetta_database: `{rosetta_database}`",
        f"- relax application: `{ROSETTA_RELAX_APPLICATION_NAME}`",
        f"- selected candidates: `{selected_ids}`",
        f"- selected top-N: `{TOP_RELAX_CANDIDATE_COUNT}`",
        f"- seed: `{ROSETTA_SEED}`",
        f"- nstruct: `{ROSETTA_NSTRUCT}`",
        f"- fa_elec weight override: `{ROSETTA_FA_ELEC_WEIGHT}`",
        "- restraint mode: `-relax:constrain_relax_to_start_coords` with constraint ramping disabled",
        "- single-threading: `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `NUMEXPR_NUM_THREADS=1`, `VECLIB_MAXIMUM_THREADS=1`",
        "- preserved input handling: representative packed PDBs reused from Phase 5A1 with waters retained, auto metal setup enabled, headers preserved, and PDB renumbering disabled",
        "- excluded in this phase: `fixbb`, `MD`, `QM`, and quantum steps",
        "",
    ]
    if ranked_results:
        lines.extend(
            [
                "| rosetta_relax_rank | candidate_id | campaign_id | backbone_id | metal | pre_relax_total_score | post_relax_total_score | delta_total_score |",
                "| ---: | --- | --- | --- | --- | ---: | ---: | ---: |",
            ]
        )
        for rank, result in enumerate(ranked_results, start=1):
            candidate = result.candidate
            lines.append(
                "| "
                f"{rank} | "
                f"{candidate.candidate_id} | "
                f"{candidate.campaign_id} | "
                f"{candidate.backbone_id} | "
                f"{candidate.preserved_metal_identity} | "
                f"{result.pre_relax_total_score:.3f} | "
                f"{_require_numeric_score(result.post_relax_total_score, candidate.candidate_id, 'post_relax_total_score'):.3f} | "
                f"{_require_numeric_score(result.delta_total_score, candidate.candidate_id, 'delta_total_score'):.3f} |"
            )
    else:
        lines.append("No successful relax refinements were available for post-relax ranking.")

    failed_results = [result for result in results if not result.is_success]
    if failed_results:
        lines.extend(
            [
                "",
                "| candidate_id | failure_status | failure_reason |",
                "| --- | --- | --- |",
            ]
        )
        for result in sorted(failed_results, key=lambda item: item.candidate.rosetta_rank):
            lines.append(
                "| "
                f"{result.candidate.candidate_id} | "
                f"{result.failure_status or 'failed'} | "
                f"{result.failure_reason or 'Unknown failure'} |"
            )
    lines.append("")
    return "\n".join(lines)


def run_rosetta_relax_smoke(
    *,
    rosetta_bin_dir: Path,
    rosetta_database: Path,
    ranking_path: Path,
    output_root: Path,
    summary_path: Path,
    ranking_output_path: Path,
    report_path: Path,
) -> tuple[tuple[dict[str, object], ...], tuple[dict[str, object], ...]]:
    """Run the deterministic Phase 5A2 restrained relax smoke workflow."""
    validated_bin_dir, validated_database = validate_rosetta_paths(rosetta_bin_dir, rosetta_database)
    relax_executable = resolve_rosetta_executable(validated_bin_dir, ROSETTA_RELAX_APPLICATION_NAME)
    candidates = select_top_rosetta_relax_candidates(ranking_path=ranking_path, output_root=output_root)
    if not candidates:
        raise ValueError(f"No Rosetta relax smoke candidates discovered from {ranking_path}")

    _reset_output_directory(output_root)
    results: list[RosettaRelaxSmokeResult] = []
    for candidate in candidates:
        LOGGER.info("Running Rosetta relax smoke refinement for %s", candidate.candidate_id)
        candidate.output_dir.mkdir(parents=True, exist_ok=True)
        score_dir = candidate.output_dir / "score"
        relax_dir = candidate.output_dir / "relax"
        score_dir.mkdir(parents=True, exist_ok=True)
        relax_dir.mkdir(parents=True, exist_ok=True)

        pre_relax_scorefile_copy = score_dir / "pre_relax_score.sc"
        shutil.copy2(candidate.pre_relax_scorefile_source, pre_relax_scorefile_copy)
        write_csv_rows(score_dir / "selected_candidate_snapshot.csv", [_candidate_snapshot_row(candidate)])

        relax_scorefile_path = relax_dir / "relax.sc"
        relax_command = build_rosetta_relax_command(
            rosetta_executable=relax_executable,
            rosetta_database=validated_database,
            input_pdb=candidate.representative_packed_pdb,
            output_dir=relax_dir,
            scorefile_path=relax_scorefile_path,
        )
        atomic_write_text(relax_dir / "relax_command.txt", shlex.join(relax_command) + "\n")
        try:
            _run_logged_command(
                command=relax_command,
                cwd=relax_dir,
                stdout_path=relax_dir / "relax_stdout.txt",
                stderr_path=relax_dir / "relax_stderr.txt",
            )
            post_relax_scores = parse_rosetta_relax_scorefile(relax_scorefile_path)
            relaxed_pdb_path = resolve_relaxed_pdb_path(
                relax_output_dir=relax_dir,
                post_relax_scores=post_relax_scores,
            )
            results.append(
                RosettaRelaxSmokeResult(
                    candidate=candidate,
                    pre_relax_scorefile_copy=pre_relax_scorefile_copy,
                    relax_scorefile_path=relax_scorefile_path,
                    relaxed_pdb_path=relaxed_pdb_path,
                    post_relax_scores=post_relax_scores,
                    success_status="success",
                    failure_status="",
                    failure_reason="",
                )
            )
        except RuntimeError as exc:
            results.append(
                RosettaRelaxSmokeResult(
                    candidate=candidate,
                    pre_relax_scorefile_copy=pre_relax_scorefile_copy,
                    relax_scorefile_path=relax_scorefile_path,
                    relaxed_pdb_path=None,
                    post_relax_scores=None,
                    success_status="failed",
                    failure_status="relax_failed",
                    failure_reason=str(exc),
                )
            )
        except (FileNotFoundError, ValueError) as exc:
            results.append(
                RosettaRelaxSmokeResult(
                    candidate=candidate,
                    pre_relax_scorefile_copy=pre_relax_scorefile_copy,
                    relax_scorefile_path=relax_scorefile_path,
                    relaxed_pdb_path=None,
                    post_relax_scores=None,
                    success_status="failed",
                    failure_status="postprocess_failed",
                    failure_reason=str(exc),
                )
            )

    summary_rows = build_rosetta_relax_summary_rows(results)
    ranking_rows = build_rosetta_relax_ranking_rows(results)
    _write_csv_rows_with_fieldnames(
        summary_path,
        summary_rows,
        _summary_fieldnames(results),
    )
    _write_csv_rows_with_fieldnames(
        ranking_output_path,
        ranking_rows,
        _ranking_fieldnames(results),
    )
    atomic_write_text(
        report_path,
        render_rosetta_relax_smoke_markdown(
            rosetta_bin_dir=validated_bin_dir,
            rosetta_database=validated_database,
            results=results,
        ),
    )
    return summary_rows, ranking_rows


def _ordered_score_terms(score_mappings: Iterable[Mapping[str, float | str]]) -> tuple[str, ...]:
    ordered_terms: list[str] = []
    seen_terms: set[str] = set()
    for score_mapping in score_mappings:
        for key, value in score_mapping.items():
            if key in {"description", "total_score"} or key in seen_terms:
                continue
            if isinstance(value, float):
                seen_terms.add(key)
                ordered_terms.append(key)
    return tuple(ordered_terms)


def _summary_fieldnames(results: Sequence[RosettaRelaxSmokeResult]) -> tuple[str, ...]:
    fieldnames: list[str] = [
        "rosetta_rank",
        "shortlist_rank",
        "candidate_id",
        "campaign_id",
        "backbone_id",
        "preserved_metal_identity",
        "representative_design_id",
        "representative_ligand_confidence",
        "representative_overall_confidence",
        "representative_packed_pdb",
        "pre_relax_scorefile_source",
        "copied_pre_relax_scorefile",
        "relax_scorefile_path",
        "relaxed_pdb_path",
        "pre_relax_total_score",
        "post_relax_total_score",
        "delta_total_score",
        "success_status",
        "failure_status",
        "failure_reason",
        "pre_relax_description",
        "post_relax_description",
    ]
    for term in _ordered_score_terms(result.candidate.pre_relax_scores for result in results):
        fieldnames.append(f"pre_relax_{term}")
    for term in _ordered_score_terms(result.post_relax_scores for result in results if result.post_relax_scores is not None):
        fieldnames.append(f"post_relax_{term}")
    return tuple(fieldnames)


def _ranking_fieldnames(results: Sequence[RosettaRelaxSmokeResult]) -> tuple[str, ...]:
    return ("rosetta_relax_rank",) + _summary_fieldnames(results)


def _candidate_snapshot_row(candidate: RosettaRelaxSmokeCandidate) -> dict[str, object]:
    return {
        "rosetta_rank": candidate.rosetta_rank,
        "shortlist_rank": candidate.shortlist_rank,
        "candidate_id": candidate.candidate_id,
        "campaign_id": candidate.campaign_id,
        "backbone_id": candidate.backbone_id,
        "preserved_metal_identity": candidate.preserved_metal_identity,
        "representative_design_id": candidate.representative_design_id,
        "representative_ligand_confidence": candidate.representative_ligand_confidence,
        "representative_overall_confidence": candidate.representative_overall_confidence,
        "representative_packed_pdb": _display_path(candidate.representative_packed_pdb),
        "pre_relax_scorefile_source": _display_path(candidate.pre_relax_scorefile_source),
        "pre_relax_total_score": candidate.pre_relax_scores.get("total_score"),
    }


def _require_numeric_score(value: float | None, candidate_id: str, field_name: str) -> float:
    if value is None:
        raise ValueError(f"Missing numeric {field_name} for {candidate_id}")
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
    environment = os.environ.copy()
    environment.update(SINGLE_THREADED_ENVIRONMENT)
    completed = subprocess.run(
        list(command),
        cwd=cwd,
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )
    atomic_write_text(stdout_path, completed.stdout)
    atomic_write_text(stderr_path, completed.stderr)
    if completed.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {completed.returncode}: {shlex.join(command)}")


def _write_csv_rows_with_fieldnames(
    path: Path,
    rows: Sequence[Mapping[str, object]],
    fieldnames: Sequence[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("w", encoding="utf-8", newline="", delete=False, dir=path.parent) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        if rows:
            writer.writerows(rows)
        temp_path = Path(handle.name)
    temp_path.replace(path)
