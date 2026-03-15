"""Run the deterministic Phase 7D LigandMPNN round-2 smoke workflow."""

from __future__ import annotations

import argparse
from pathlib import Path

from lanm.analysis.ligandmpnn_round2_smoke import run_ligandmpnn_round2_smoke
from lanm.logging_utils import configure_logging
from lanm.paths import (
    LIGANDMPNN_ROUND2_INPUT_MANIFEST_PATH,
    LIGANDMPNN_ROUND2_SEQUENCE_CATALOG_PATH,
    LIGANDMPNN_ROUND2_SMOKE_REPORT_PATH,
    LIGANDMPNN_ROUND2_SMOKE_SUMMARY_PATH,
    RESULTS_LIGANDMPNN_ROUND2_SMOKE_DIR,
    ensure_runtime_directories,
)


def _build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the deterministic Phase 7D LigandMPNN round-2 smoke workflow.",
    )
    parser.add_argument(
        "--ligandmpnn-root",
        required=True,
        type=Path,
        help="Absolute path to the local LigandMPNN checkout.",
    )
    return parser


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    args = _build_argument_parser().parse_args()
    try:
        run_ligandmpnn_round2_smoke(
            ligandmpnn_root=args.ligandmpnn_root,
            manifest_path=LIGANDMPNN_ROUND2_INPUT_MANIFEST_PATH,
            output_root=RESULTS_LIGANDMPNN_ROUND2_SMOKE_DIR,
            summary_path=LIGANDMPNN_ROUND2_SMOKE_SUMMARY_PATH,
            sequence_catalog_path=LIGANDMPNN_ROUND2_SEQUENCE_CATALOG_PATH,
            report_path=LIGANDMPNN_ROUND2_SMOKE_REPORT_PATH,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
