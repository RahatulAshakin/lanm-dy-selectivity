"""Run the deterministic Phase 3B ProteinMPNN smoke workflow."""

from __future__ import annotations

import argparse
from pathlib import Path

from lanm.analysis.proteinmpnn_smoke import run_proteinmpnn_smoke
from lanm.logging_utils import configure_logging
from lanm.paths import (
    DESIGN_BACKBONE_MANIFEST_PATH,
    DESIGN_CAMPAIGN_MANIFEST_PATH,
    PROTEINMPNN_CAMPAIGNS_DIR,
    PROTEINMPNN_SEQUENCE_CATALOG_PATH,
    PROTEINMPNN_SMOKE_REPORT_PATH,
    PROTEINMPNN_SMOKE_SUMMARY_PATH,
    RESULTS_PROTEINMPNN_SMOKE_DIR,
    ensure_runtime_directories,
)


def _build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the deterministic Phase 3B ProteinMPNN smoke workflow.",
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
        run_proteinmpnn_smoke(
            proteinmpnn_root=args.proteinmpnn_root,
            campaign_manifest_path=DESIGN_CAMPAIGN_MANIFEST_PATH,
            backbone_manifest_path=DESIGN_BACKBONE_MANIFEST_PATH,
            campaigns_dir=PROTEINMPNN_CAMPAIGNS_DIR,
            output_root=RESULTS_PROTEINMPNN_SMOKE_DIR,
            summary_path=PROTEINMPNN_SMOKE_SUMMARY_PATH,
            sequence_catalog_path=PROTEINMPNN_SEQUENCE_CATALOG_PATH,
            report_path=PROTEINMPNN_SMOKE_REPORT_PATH,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
