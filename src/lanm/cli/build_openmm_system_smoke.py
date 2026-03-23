"""Run the deterministic Phase 6A3 representative OpenMM system-build smoke test."""

from __future__ import annotations

from lanm.analysis.openmm_system_smoke import run_openmm_system_smoke
from lanm.logging_utils import configure_logging
from lanm.paths import (
    MD_PROTOCOL_CONFIG_PATH,
    MD_SYSTEM_BUILD_MANIFEST_PATH,
    MD_VALIDATION_PANEL_PATH,
    METAL_BUILD_READINESS_PATH,
    METAL_PARAMETER_VALUES_PATH,
    OPENMM_SYSTEM_SMOKE_REPORT_PATH,
    OPENMM_SYSTEM_SMOKE_SUMMARY_PATH,
    RESULTS_OPENMM_SYSTEM_SMOKE_DIR,
    ensure_runtime_directories,
)


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    try:
        run_openmm_system_smoke(
            md_validation_panel_path=MD_VALIDATION_PANEL_PATH,
            md_system_build_manifest_path=MD_SYSTEM_BUILD_MANIFEST_PATH,
            metal_build_readiness_path=METAL_BUILD_READINESS_PATH,
            md_protocol_path=MD_PROTOCOL_CONFIG_PATH,
            metal_parameter_values_path=METAL_PARAMETER_VALUES_PATH,
            output_root=RESULTS_OPENMM_SYSTEM_SMOKE_DIR,
            summary_path=OPENMM_SYSTEM_SMOKE_SUMMARY_PATH,
            report_path=OPENMM_SYSTEM_SMOKE_REPORT_PATH,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
