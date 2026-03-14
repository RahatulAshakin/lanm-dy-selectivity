"""Run the deterministic Phase 5A1 Rosetta score_jd2 baseline workflow."""

from __future__ import annotations

import argparse
from pathlib import Path

from lanm.analysis.rosetta_score_smoke import run_rosetta_score_smoke
from lanm.logging_utils import configure_logging
from lanm.paths import (
    LIGANDMPNN_SEQUENCE_CATALOG_PATH,
    LIGANDMPNN_SHORTLIST_PATH,
    RESULTS_ROSETTA_SCORE_SMOKE_DIR,
    ROSETTA_SCORE_CANDIDATE_RANKING_PATH,
    ROSETTA_SCORE_SMOKE_REPORT_PATH,
    ROSETTA_SCORE_SMOKE_SUMMARY_PATH,
    ensure_runtime_directories,
)


def _build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the deterministic Phase 5A1 Rosetta score_jd2 baseline workflow.",
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
        run_rosetta_score_smoke(
            rosetta_bin_dir=args.rosetta_bin_dir,
            rosetta_database=args.rosetta_database,
            shortlist_path=LIGANDMPNN_SHORTLIST_PATH,
            sequence_catalog_path=LIGANDMPNN_SEQUENCE_CATALOG_PATH,
            output_root=RESULTS_ROSETTA_SCORE_SMOKE_DIR,
            summary_path=ROSETTA_SCORE_SMOKE_SUMMARY_PATH,
            ranking_path=ROSETTA_SCORE_CANDIDATE_RANKING_PATH,
            report_path=ROSETTA_SCORE_SMOKE_REPORT_PATH,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
