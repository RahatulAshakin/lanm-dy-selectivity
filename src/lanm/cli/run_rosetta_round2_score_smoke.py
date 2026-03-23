"""Run the deterministic Phase 7F Rosetta round-2 score-only triage workflow."""

from __future__ import annotations

import argparse
from pathlib import Path

from lanm.analysis.rosetta_round2_score_smoke import run_rosetta_round2_score_smoke
from lanm.logging_utils import configure_logging
from lanm.paths import (
    LIGANDMPNN_ROUND2_SEQUENCE_CATALOG_PATH,
    LIGANDMPNN_ROUND2_SHORTLIST_PATH,
    RESULTS_ROSETTA_ROUND2_SCORE_SMOKE_DIR,
    ROSETTA_ROUND2_CANDIDATE_RANKING_PATH,
    ROSETTA_ROUND2_SCORE_SMOKE_REPORT_PATH,
    ROSETTA_ROUND2_SCORE_SUMMARY_PATH,
    ROUND2_MD_RESCREEN_PANEL_PATH,
    ensure_runtime_directories,
)


def _build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the deterministic Phase 7F Rosetta round-2 score-only triage workflow.",
    )
    parser.add_argument(
        "--rosetta-bin-dir",
        required=True,
        type=Path,
        help="Absolute path to the local Rosetta main/source/bin directory.",
    )
    parser.add_argument(
        "--rosetta-database",
        required=True,
        type=Path,
        help="Absolute path to the local Rosetta main/database directory.",
    )
    return parser


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    args = _build_argument_parser().parse_args()
    try:
        run_rosetta_round2_score_smoke(
            rosetta_bin_dir=args.rosetta_bin_dir,
            rosetta_database=args.rosetta_database,
            shortlist_path=LIGANDMPNN_ROUND2_SHORTLIST_PATH,
            sequence_catalog_path=LIGANDMPNN_ROUND2_SEQUENCE_CATALOG_PATH,
            output_root=RESULTS_ROSETTA_ROUND2_SCORE_SMOKE_DIR,
            summary_path=ROSETTA_ROUND2_SCORE_SUMMARY_PATH,
            ranking_path=ROSETTA_ROUND2_CANDIDATE_RANKING_PATH,
            panel_path=ROUND2_MD_RESCREEN_PANEL_PATH,
            report_path=ROSETTA_ROUND2_SCORE_SMOKE_REPORT_PATH,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
