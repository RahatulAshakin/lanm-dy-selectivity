"""Prepare deterministic Phase 6A1 OpenMM-ready multi-metal MD inputs."""

from __future__ import annotations

from lanm.analysis.md_input_preparation import prepare_md_inputs
from lanm.logging_utils import configure_logging
from lanm.paths import (
    DESIGN_BACKBONE_MANIFEST_PATH,
    MD_INPUT_MANIFEST_PATH,
    MD_INPUT_PREPARATION_REPORT_PATH,
    MD_PANEL_CONFIG_PATH,
    MD_PROTOCOL_CONFIG_PATH,
    MD_VALIDATION_PANEL_PATH,
    RESULTS_MD_INPUTS_DIR,
    ensure_runtime_directories,
)


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    try:
        prepare_md_inputs(
            panel_path=MD_VALIDATION_PANEL_PATH,
            md_panel_config_path=MD_PANEL_CONFIG_PATH,
            design_backbone_manifest_path=DESIGN_BACKBONE_MANIFEST_PATH,
            manifest_path=MD_INPUT_MANIFEST_PATH,
            report_path=MD_INPUT_PREPARATION_REPORT_PATH,
            protocol_config_path=MD_PROTOCOL_CONFIG_PATH,
            md_inputs_root=RESULTS_MD_INPUTS_DIR,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
