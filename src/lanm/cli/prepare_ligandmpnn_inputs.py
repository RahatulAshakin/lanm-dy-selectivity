"""Prepare deterministic Phase 4A LigandMPNN inputs for shortlisted candidates."""

from __future__ import annotations

from lanm.analysis.ligandmpnn_inputs import prepare_ligandmpnn_inputs
from lanm.logging_utils import configure_logging
from lanm.paths import (
    DESIGN_BACKBONE_MANIFEST_PATH,
    DESIGN_CAMPAIGN_MANIFEST_PATH,
    DESIGN_CAMPAIGN_POSITIONS_PATH,
    LIGANDMPNN_INPUTS_CONFIG_PATH,
    LIGANDMPNN_INPUTS_DIR,
    LIGANDMPNN_INPUTS_REPORT_PATH,
    LIGANDMPNN_INPUT_MANIFEST_PATH,
    LIGANDMPNN_REDESIGN_POSITIONS_PATH,
    PROTEINMPNN_SHORTLIST_PATH,
    ensure_runtime_directories,
)


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    try:
        prepare_ligandmpnn_inputs(
            shortlist_path=PROTEINMPNN_SHORTLIST_PATH,
            campaign_manifest_path=DESIGN_CAMPAIGN_MANIFEST_PATH,
            backbone_manifest_path=DESIGN_BACKBONE_MANIFEST_PATH,
            design_campaign_positions_path=DESIGN_CAMPAIGN_POSITIONS_PATH,
            manifest_path=LIGANDMPNN_INPUT_MANIFEST_PATH,
            redesign_positions_path=LIGANDMPNN_REDESIGN_POSITIONS_PATH,
            report_path=LIGANDMPNN_INPUTS_REPORT_PATH,
            input_root=LIGANDMPNN_INPUTS_DIR,
            config_path=LIGANDMPNN_INPUTS_CONFIG_PATH,
        )
    except FileNotFoundError as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
