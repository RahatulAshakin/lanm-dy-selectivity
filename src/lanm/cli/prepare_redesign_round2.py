"""Prepare deterministic Phase 7A metadynamics-informed round-2 redesign inputs."""

from __future__ import annotations

from lanm.analysis.redesign_round2 import (
    build_redesign_round2_artifacts,
    write_redesign_round2_outputs,
)
from lanm.logging_utils import configure_logging
from lanm.paths import (
    DESIGN_CAMPAIGN_MANIFEST_PATH,
    DESIGN_CAMPAIGN_POSITIONS_PATH,
    DESIGN_MASK_CANDIDATES_PATH,
    DESIGN_MASKS_PATH,
    LIGANDMPNN_SHORTLIST_PATH,
    MD_VALIDATION_PANEL_PATH,
    OPENMM_METADYNAMICS_PANEL_STATUS_PATH,
    OPENMM_METADYNAMICS_SUMMARY_PATH,
    PROTEINMPNN_ROUND2_INPUTS_DIR,
    REDESIGN_ROUND2_CONFIG_PATH,
    REDESIGN_ROUND2_POSITIONS_PATH,
    REDESIGN_ROUND2_REPORT_PATH,
    REDESIGN_ROUND2_SEEDS_PATH,
    ensure_runtime_directories,
)


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    try:
        artifacts = build_redesign_round2_artifacts(
            openmm_metadynamics_summary_path=OPENMM_METADYNAMICS_SUMMARY_PATH,
            openmm_metadynamics_panel_status_path=OPENMM_METADYNAMICS_PANEL_STATUS_PATH,
            design_mask_candidates_path=DESIGN_MASK_CANDIDATES_PATH,
            design_masks_path=DESIGN_MASKS_PATH,
            ligandmpnn_shortlist_path=LIGANDMPNN_SHORTLIST_PATH,
            md_validation_panel_path=MD_VALIDATION_PANEL_PATH,
            design_campaign_manifest_path=DESIGN_CAMPAIGN_MANIFEST_PATH,
            design_campaign_positions_path=DESIGN_CAMPAIGN_POSITIONS_PATH,
            proteinmpnn_round2_root=PROTEINMPNN_ROUND2_INPUTS_DIR,
        )
        write_redesign_round2_outputs(
            artifacts=artifacts,
            config_path=REDESIGN_ROUND2_CONFIG_PATH,
            seed_table_path=REDESIGN_ROUND2_SEEDS_PATH,
            position_table_path=REDESIGN_ROUND2_POSITIONS_PATH,
            report_path=REDESIGN_ROUND2_REPORT_PATH,
            proteinmpnn_round2_root=PROTEINMPNN_ROUND2_INPUTS_DIR,
        )
    except FileNotFoundError as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
