"""Run Phase 6B1 deterministic OpenMM equilibration smoke tests for baseline systems."""

from __future__ import annotations

from lanm.analysis.openmm_equilibration_smoke import run_openmm_equilibration_smoke
from lanm.logging_utils import configure_logging
from lanm.paths import (
    OPENMM_EQUILIBRATION_SMOKE_REPORT_PATH,
    OPENMM_EQUILIBRATION_SMOKE_SUMMARY_PATH,
    OPENMM_SYSTEM_BUILD_SUMMARY_PATH,
    RESULTS_OPENMM_EQUILIBRATION_SMOKE_DIR,
    ensure_runtime_directories,
)


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    try:
        run_openmm_equilibration_smoke(
            system_build_summary_path=OPENMM_SYSTEM_BUILD_SUMMARY_PATH,
            output_root=RESULTS_OPENMM_EQUILIBRATION_SMOKE_DIR,
            summary_path=OPENMM_EQUILIBRATION_SMOKE_SUMMARY_PATH,
            report_path=OPENMM_EQUILIBRATION_SMOKE_REPORT_PATH,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
