"""Run Phase 7G1 reduced OpenMM MD rescreening for the round-2 finalists."""

from __future__ import annotations

from lanm.analysis.round2_openmm_screening import run_round2_openmm_screening
from lanm.logging_utils import configure_logging
from lanm.paths import (
    MD_PROTOCOL_CONFIG_PATH,
    METAL_PARAMETER_VALUES_PATH,
    RESULTS_ROUND2_OPENMM_SCREENING_DIR,
    ROUND2_MD_RESCREEN_PANEL_PATH,
    ROUND2_OPENMM_SCREENING_PANEL_STATUS_PATH,
    ROUND2_OPENMM_SCREENING_REPORT_PATH,
    ROUND2_OPENMM_SCREENING_SUMMARY_PATH,
    ensure_runtime_directories,
)


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    try:
        run_round2_openmm_screening(
            panel_path=ROUND2_MD_RESCREEN_PANEL_PATH,
            md_protocol_path=MD_PROTOCOL_CONFIG_PATH,
            metal_parameter_values_path=METAL_PARAMETER_VALUES_PATH,
            output_root=RESULTS_ROUND2_OPENMM_SCREENING_DIR,
            summary_path=ROUND2_OPENMM_SCREENING_SUMMARY_PATH,
            panel_status_path=ROUND2_OPENMM_SCREENING_PANEL_STATUS_PATH,
            report_path=ROUND2_OPENMM_SCREENING_REPORT_PATH,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
