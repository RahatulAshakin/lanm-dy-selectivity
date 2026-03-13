"""Select deterministic Phase 3C ProteinMPNN candidates from existing smoke outputs."""

from __future__ import annotations

import argparse

from lanm.analysis.proteinmpnn_candidates import select_proteinmpnn_candidates
from lanm.logging_utils import configure_logging
from lanm.paths import (
    DESIGN_BACKBONE_MANIFEST_PATH,
    DESIGN_CAMPAIGN_MANIFEST_PATH,
    DESIGN_CAMPAIGN_POSITIONS_PATH,
    PROTEINMPNN_CAMPAIGNS_DIR,
    PROTEINMPNN_SHORTLIST_FASTA_PATH,
    PROTEINMPNN_SHORTLIST_PATH,
    PROTEINMPNN_SHORTLIST_REPORT_PATH,
    PROTEINMPNN_UNIQUE_SEQUENCES_PATH,
    RESULTS_PROTEINMPNN_SMOKE_DIR,
    ensure_runtime_directories,
)


def _build_argument_parser() -> argparse.ArgumentParser:
    return argparse.ArgumentParser(
        description="Select deterministic Phase 3C ProteinMPNN candidates from existing smoke outputs.",
    )


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    _build_argument_parser().parse_args()
    try:
        select_proteinmpnn_candidates(
            campaign_manifest_path=DESIGN_CAMPAIGN_MANIFEST_PATH,
            backbone_manifest_path=DESIGN_BACKBONE_MANIFEST_PATH,
            campaigns_dir=PROTEINMPNN_CAMPAIGNS_DIR,
            smoke_output_root=RESULTS_PROTEINMPNN_SMOKE_DIR,
            design_campaign_positions_path=DESIGN_CAMPAIGN_POSITIONS_PATH,
            unique_sequences_path=PROTEINMPNN_UNIQUE_SEQUENCES_PATH,
            shortlist_path=PROTEINMPNN_SHORTLIST_PATH,
            report_path=PROTEINMPNN_SHORTLIST_REPORT_PATH,
            shortlist_fasta_path=PROTEINMPNN_SHORTLIST_FASTA_PATH,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
