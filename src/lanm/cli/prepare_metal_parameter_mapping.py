"""Prepare deterministic Phase 6A2c metal-parameter readiness outputs."""

from __future__ import annotations

from lanm.analysis.metal_parameter_mapping import prepare_metal_parameter_mapping
from lanm.logging_utils import configure_logging
from lanm.paths import (
    MD_SYSTEM_BUILD_MANIFEST_PATH,
    METAL_BUILD_READINESS_PATH,
    METAL_MODEL_REGISTRY_PATH,
    METAL_PARAMETER_CONFIG_PATH,
    METAL_PARAMETER_MAPPING_PATH,
    METAL_PARAMETER_STRATEGY_REPORT_PATH,
    METAL_PARAMETER_VALUES_PATH,
    ensure_runtime_directories,
)


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    try:
        prepare_metal_parameter_mapping(
            metal_model_registry_path=METAL_MODEL_REGISTRY_PATH,
            md_system_build_manifest_path=MD_SYSTEM_BUILD_MANIFEST_PATH,
            metal_parameter_config_path=METAL_PARAMETER_CONFIG_PATH,
            metal_parameter_values_path=METAL_PARAMETER_VALUES_PATH,
            metal_parameter_mapping_path=METAL_PARAMETER_MAPPING_PATH,
            metal_build_readiness_path=METAL_BUILD_READINESS_PATH,
            metal_parameter_strategy_report_path=METAL_PARAMETER_STRATEGY_REPORT_PATH,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
