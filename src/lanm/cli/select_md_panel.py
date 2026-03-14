"""Select the deterministic Phase 5B integrated ranking and MD validation panel."""

from __future__ import annotations

import argparse

from lanm.analysis.md_panel_selection import run_md_panel_selection
from lanm.logging_utils import configure_logging
from lanm.paths import (
    INTEGRATED_CANDIDATE_RANKING_PATH,
    LIGANDMPNN_SHORTLIST_PATH,
    LIGANDMPNN_SMOKE_SUMMARY_PATH,
    MD_PANEL_CONFIG_PATH,
    MD_PANEL_SELECTION_REPORT_PATH,
    MD_VALIDATION_PANEL_PATH,
    PROTEINMPNN_SHORTLIST_PATH,
    ROSETTA_SCORE_CANDIDATE_RANKING_PATH,
    ensure_runtime_directories,
)


def _build_argument_parser() -> argparse.ArgumentParser:
    return argparse.ArgumentParser(
        description="Select the deterministic Phase 5B integrated ranking and MD validation panel.",
    )


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    _build_argument_parser().parse_args()
    try:
        run_md_panel_selection(
            proteinmpnn_shortlist_path=PROTEINMPNN_SHORTLIST_PATH,
            ligandmpnn_shortlist_path=LIGANDMPNN_SHORTLIST_PATH,
            ligandmpnn_smoke_summary_path=LIGANDMPNN_SMOKE_SUMMARY_PATH,
            rosetta_score_ranking_path=ROSETTA_SCORE_CANDIDATE_RANKING_PATH,
            integrated_ranking_path=INTEGRATED_CANDIDATE_RANKING_PATH,
            panel_path=MD_VALIDATION_PANEL_PATH,
            report_path=MD_PANEL_SELECTION_REPORT_PATH,
            md_panel_config_path=MD_PANEL_CONFIG_PATH,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
