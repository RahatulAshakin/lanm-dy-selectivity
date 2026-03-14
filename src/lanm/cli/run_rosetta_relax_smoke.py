"""Run the deterministic Phase 5A2 Rosetta relax smoke workflow."""

from __future__ import annotations

import argparse
from pathlib import Path

from lanm.analysis.rosetta_relax_smoke import run_rosetta_relax_smoke
from lanm.logging_utils import configure_logging
from lanm.paths import (
    RESULTS_ROSETTA_RELAX_SMOKE_DIR,
    ROSETTA_RELAX_CANDIDATE_RANKING_PATH,
    ROSETTA_RELAX_SMOKE_REPORT_PATH,
    ROSETTA_RELAX_SMOKE_SUMMARY_PATH,
    ROSETTA_SCORE_CANDIDATE_RANKING_PATH,
    ensure_runtime_directories,
)


def _build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the deterministic Phase 5A2 Rosetta relax smoke workflow.",
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
        run_rosetta_relax_smoke(
            rosetta_bin_dir=args.rosetta_bin_dir,
            rosetta_database=args.rosetta_database,
            ranking_path=ROSETTA_SCORE_CANDIDATE_RANKING_PATH,
            output_root=RESULTS_ROSETTA_RELAX_SMOKE_DIR,
            summary_path=ROSETTA_RELAX_SMOKE_SUMMARY_PATH,
            ranking_output_path=ROSETTA_RELAX_CANDIDATE_RANKING_PATH,
            report_path=ROSETTA_RELAX_SMOKE_REPORT_PATH,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
