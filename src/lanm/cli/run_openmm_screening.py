"""Run Phase 6B2 deterministic OpenMM screening across the full MD panel."""

from __future__ import annotations

from lanm.analysis.openmm_screening import run_openmm_screening
from lanm.logging_utils import configure_logging
from lanm.paths import (
    OPENMM_EQUILIBRATION_SMOKE_SUMMARY_PATH,
    OPENMM_SCREENING_PANEL_STATUS_PATH,
    OPENMM_SCREENING_REPLICATES_PATH,
    OPENMM_SCREENING_REPORT_PATH,
    OPENMM_SCREENING_SUMMARY_PATH,
    OPENMM_SYSTEM_BUILD_SUMMARY_PATH,
    RESULTS_OPENMM_SCREENING_DIR,
    ensure_runtime_directories,
)


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    try:
        run_openmm_screening(
            system_build_summary_path=OPENMM_SYSTEM_BUILD_SUMMARY_PATH,
            equilibration_summary_path=OPENMM_EQUILIBRATION_SMOKE_SUMMARY_PATH,
            output_root=RESULTS_OPENMM_SCREENING_DIR,
            replicate_summary_path=OPENMM_SCREENING_REPLICATES_PATH,
            system_summary_path=OPENMM_SCREENING_SUMMARY_PATH,
            panel_status_path=OPENMM_SCREENING_PANEL_STATUS_PATH,
            report_path=OPENMM_SCREENING_REPORT_PATH,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
