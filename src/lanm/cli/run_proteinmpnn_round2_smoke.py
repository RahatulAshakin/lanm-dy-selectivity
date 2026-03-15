"""Run the deterministic Phase 7B ProteinMPNN round-2 smoke workflow."""

from __future__ import annotations

import argparse
from pathlib import Path

from lanm.analysis.proteinmpnn_round2_smoke import run_proteinmpnn_round2_smoke
from lanm.logging_utils import configure_logging
from lanm.paths import (
    PROTEINMPNN_ROUND2_CAMPAIGNS_DIR,
    PROTEINMPNN_ROUND2_SEQUENCE_CATALOG_PATH,
    PROTEINMPNN_ROUND2_SMOKE_REPORT_PATH,
    PROTEINMPNN_ROUND2_SUMMARY_PATH,
    REDESIGN_ROUND2_CONFIG_PATH,
    REDESIGN_ROUND2_POSITIONS_PATH,
    REDESIGN_ROUND2_SEEDS_PATH,
    RESULTS_PROTEINMPNN_ROUND2_SMOKE_DIR,
    ensure_runtime_directories,
)


def _build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the deterministic Phase 7B ProteinMPNN round-2 smoke workflow.",
    )
    parser.add_argument(
        "--proteinmpnn-root",
        required=True,
        type=Path,
        help="Absolute path to the local ProteinMPNN checkout.",
    )
    return parser


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    args = _build_argument_parser().parse_args()
    try:
        run_proteinmpnn_round2_smoke(
            proteinmpnn_root=args.proteinmpnn_root,
            redesign_round2_config_path=REDESIGN_ROUND2_CONFIG_PATH,
            round2_seed_table_path=REDESIGN_ROUND2_SEEDS_PATH,
            round2_position_table_path=REDESIGN_ROUND2_POSITIONS_PATH,
            campaigns_dir=PROTEINMPNN_ROUND2_CAMPAIGNS_DIR,
            output_root=RESULTS_PROTEINMPNN_ROUND2_SMOKE_DIR,
            summary_path=PROTEINMPNN_ROUND2_SUMMARY_PATH,
            sequence_catalog_path=PROTEINMPNN_ROUND2_SEQUENCE_CATALOG_PATH,
            report_path=PROTEINMPNN_ROUND2_SMOKE_REPORT_PATH,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
