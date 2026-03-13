"""Prepare deterministic Phase 3A design campaign inputs for ProteinMPNN."""

from __future__ import annotations

from lanm.analysis.design_campaigns import (
    build_design_campaign_artifacts,
    write_design_campaign_outputs,
)
from lanm.logging_utils import configure_logging
from lanm.paths import (
    CROSS_TEMPLATE_RESIDUE_ALIGNMENT_PATH,
    DESIGN_BACKBONE_MANIFEST_PATH,
    DESIGN_CAMPAIGN_MANIFEST_PATH,
    DESIGN_CAMPAIGN_POSITIONS_PATH,
    DESIGN_CAMPAIGNS_PATH,
    DESIGN_CAMPAIGNS_REPORT_PATH,
    DESIGN_MASK_CANDIDATES_PATH,
    DESIGN_MASKS_PATH,
    PROTEINMPNN_INPUTS_DIR,
    TEMPLATE_CHAIN_SUMMARY_PATH,
    ensure_runtime_directories,
)


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    try:
        artifacts = build_design_campaign_artifacts(
            design_masks_path=DESIGN_MASKS_PATH,
            design_mask_candidates_path=DESIGN_MASK_CANDIDATES_PATH,
            cross_template_alignment_path=CROSS_TEMPLATE_RESIDUE_ALIGNMENT_PATH,
            template_chain_summary_path=TEMPLATE_CHAIN_SUMMARY_PATH,
        )
        write_design_campaign_outputs(
            artifacts=artifacts,
            report_path=DESIGN_CAMPAIGNS_REPORT_PATH,
            campaign_positions_path=DESIGN_CAMPAIGN_POSITIONS_PATH,
            backbone_manifest_path=DESIGN_BACKBONE_MANIFEST_PATH,
            campaign_manifest_path=DESIGN_CAMPAIGN_MANIFEST_PATH,
            config_path=DESIGN_CAMPAIGNS_PATH,
            proteinmpnn_root=PROTEINMPNN_INPUTS_DIR,
        )
    except FileNotFoundError as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
