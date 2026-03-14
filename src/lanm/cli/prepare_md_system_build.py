"""Prepare deterministic Phase 6A2 OpenMM system-build scaffolding."""

from __future__ import annotations

from lanm.analysis.md_system_build import prepare_md_system_build
from lanm.logging_utils import configure_logging
from lanm.paths import (
    MD_INPUT_MANIFEST_PATH,
    MD_PROTOCOL_CONFIG_PATH,
    MD_SYSTEM_BUILDING_PLAN_REPORT_PATH,
    MD_SYSTEM_BUILD_MANIFEST_PATH,
    METAL_MODEL_REGISTRY_PATH,
    RESULTS_MD_SYSTEM_BUILD_DIR,
    ensure_runtime_directories,
)


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    try:
        prepare_md_system_build(
            md_input_manifest_path=MD_INPUT_MANIFEST_PATH,
            md_protocol_path=MD_PROTOCOL_CONFIG_PATH,
            metal_model_registry_path=METAL_MODEL_REGISTRY_PATH,
            md_system_build_root=RESULTS_MD_SYSTEM_BUILD_DIR,
            md_system_build_manifest_path=MD_SYSTEM_BUILD_MANIFEST_PATH,
            md_system_build_report_path=MD_SYSTEM_BUILDING_PLAN_REPORT_PATH,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
