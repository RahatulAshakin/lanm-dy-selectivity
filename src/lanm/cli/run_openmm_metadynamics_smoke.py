"""Run Phase 6C1 reduced-panel well-tempered OpenMM metadynamics triage."""

from __future__ import annotations

import argparse

from lanm.analysis.openmm_metadynamics_smoke import run_openmm_metadynamics_smoke
from lanm.logging_utils import configure_logging
from lanm.paths import (
    OPENMM_METADYNAMICS_PANEL_STATUS_PATH,
    OPENMM_METADYNAMICS_REPORT_PATH,
    OPENMM_METADYNAMICS_SUMMARY_PATH,
    OPENMM_SCREENING_PANEL_STATUS_PATH,
    OPENMM_SYSTEM_BUILD_SUMMARY_PATH,
    RESULTS_OPENMM_METADYNAMICS_SMOKE_DIR,
    ensure_runtime_directories,
)


def _build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run Phase 6C1 reduced-panel well-tempered OpenMM metadynamics triage.",
    )
    parser.add_argument(
        "--max-workers",
        default=1,
        type=int,
        help="Maximum number of Phase 6C1 systems to process at once.",
    )
    return parser


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    args = _build_argument_parser().parse_args()
    try:
        run_openmm_metadynamics_smoke(
            screening_panel_status_path=OPENMM_SCREENING_PANEL_STATUS_PATH,
            system_build_summary_path=OPENMM_SYSTEM_BUILD_SUMMARY_PATH,
            output_root=RESULTS_OPENMM_METADYNAMICS_SMOKE_DIR,
            system_summary_path=OPENMM_METADYNAMICS_SUMMARY_PATH,
            panel_status_path=OPENMM_METADYNAMICS_PANEL_STATUS_PATH,
            report_path=OPENMM_METADYNAMICS_REPORT_PATH,
            max_workers=args.max_workers,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
