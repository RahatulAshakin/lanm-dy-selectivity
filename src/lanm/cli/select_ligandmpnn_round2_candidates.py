"""Select deterministic Phase 7E LigandMPNN round-2 candidates for Rosetta triage."""

from __future__ import annotations

import argparse

from lanm.analysis.ligandmpnn_round2_candidates import select_ligandmpnn_round2_candidates
from lanm.logging_utils import configure_logging
from lanm.paths import (
    LIGANDMPNN_ROUND2_INPUT_MANIFEST_PATH,
    LIGANDMPNN_ROUND2_REDESIGN_POSITIONS_PATH,
    LIGANDMPNN_ROUND2_SHORTLIST_FASTA_PATH,
    LIGANDMPNN_ROUND2_SHORTLIST_PATH,
    LIGANDMPNN_ROUND2_SHORTLIST_REPORT_PATH,
    LIGANDMPNN_ROUND2_UNIQUE_SEQUENCES_PATH,
    PROTEINMPNN_ROUND2_SHORTLIST_PATH,
    RESULTS_LIGANDMPNN_ROUND2_SMOKE_DIR,
    ROSETTA_ROUND2_SHORTLIST_PATH,
    ensure_runtime_directories,
)


def _build_argument_parser() -> argparse.ArgumentParser:
    return argparse.ArgumentParser(
        description="Select deterministic Phase 7E LigandMPNN round-2 candidates for Rosetta triage.",
    )


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    _build_argument_parser().parse_args()
    try:
        select_ligandmpnn_round2_candidates(
            proteinmpnn_round2_shortlist_path=PROTEINMPNN_ROUND2_SHORTLIST_PATH,
            manifest_path=LIGANDMPNN_ROUND2_INPUT_MANIFEST_PATH,
            redesign_positions_path=LIGANDMPNN_ROUND2_REDESIGN_POSITIONS_PATH,
            smoke_output_root=RESULTS_LIGANDMPNN_ROUND2_SMOKE_DIR,
            unique_sequences_path=LIGANDMPNN_ROUND2_UNIQUE_SEQUENCES_PATH,
            shortlist_path=LIGANDMPNN_ROUND2_SHORTLIST_PATH,
            report_path=LIGANDMPNN_ROUND2_SHORTLIST_REPORT_PATH,
            shortlist_fasta_path=LIGANDMPNN_ROUND2_SHORTLIST_FASTA_PATH,
            rosetta_shortlist_path=ROSETTA_ROUND2_SHORTLIST_PATH,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
