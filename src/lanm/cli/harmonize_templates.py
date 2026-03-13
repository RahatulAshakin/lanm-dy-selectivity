"""Generate deterministic Phase 2A/2B/2C template harmonization outputs."""

from __future__ import annotations

from lanm.analysis.template_harmonization import (
    build_template_harmonization_artifacts,
    write_template_harmonization_outputs,
)
from lanm.logging_utils import configure_logging
from lanm.paths import (
    CROSS_TEMPLATE_RESIDUE_ALIGNMENT_PATH,
    DESIGN_MASK_CANDIDATES_PATH,
    DESIGN_MASKS_PATH,
    LANMODULIN_SEQUENCES_PATH,
    LOCAL_STRUCTURE_MANIFEST_PATH,
    RESIDUE_ROLE_MAP_PATH,
    TEMPLATE_CHAIN_SUMMARY_PATH,
    TEMPLATE_HARMONIZATION_FIGURE_PATH,
    TEMPLATE_HARMONIZATION_REPORT_PATH,
    TEMPLATE_SITE_SUMMARY_PATH,
    ensure_runtime_directories,
)


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    try:
        artifacts = build_template_harmonization_artifacts(
            sequences_path=LANMODULIN_SEQUENCES_PATH,
            manifest_path=LOCAL_STRUCTURE_MANIFEST_PATH,
        )
        write_template_harmonization_outputs(
            artifacts=artifacts,
            report_path=TEMPLATE_HARMONIZATION_REPORT_PATH,
            chain_summary_path=TEMPLATE_CHAIN_SUMMARY_PATH,
            site_summary_path=TEMPLATE_SITE_SUMMARY_PATH,
            cross_template_alignment_path=CROSS_TEMPLATE_RESIDUE_ALIGNMENT_PATH,
            residue_role_map_path=RESIDUE_ROLE_MAP_PATH,
            design_mask_candidates_path=DESIGN_MASK_CANDIDATES_PATH,
            design_masks_path=DESIGN_MASKS_PATH,
            figure_path=TEMPLATE_HARMONIZATION_FIGURE_PATH,
        )
    except FileNotFoundError as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
