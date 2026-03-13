"""Select deterministic Phase 4C LigandMPNN candidates from existing smoke outputs."""

from __future__ import annotations

import argparse

from lanm.analysis.ligandmpnn_candidates import select_ligandmpnn_candidates
from lanm.logging_utils import configure_logging
from lanm.paths import (
    LIGANDMPNN_INPUT_MANIFEST_PATH,
    LIGANDMPNN_REDESIGN_POSITIONS_PATH,
    LIGANDMPNN_SHORTLIST_FASTA_PATH,
    LIGANDMPNN_SHORTLIST_PATH,
    LIGANDMPNN_SHORTLIST_REPORT_PATH,
    LIGANDMPNN_UNIQUE_SEQUENCES_PATH,
    PROTEINMPNN_SHORTLIST_PATH,
    RESULTS_LIGANDMPNN_SMOKE_DIR,
    ROSETTA_SHORTLIST_PATH,
    ensure_runtime_directories,
)


def _build_argument_parser() -> argparse.ArgumentParser:
    return argparse.ArgumentParser(
        description="Select deterministic Phase 4C LigandMPNN candidates from existing smoke outputs.",
    )


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    _build_argument_parser().parse_args()
    try:
        select_ligandmpnn_candidates(
            proteinmpnn_shortlist_path=PROTEINMPNN_SHORTLIST_PATH,
            manifest_path=LIGANDMPNN_INPUT_MANIFEST_PATH,
            redesign_positions_path=LIGANDMPNN_REDESIGN_POSITIONS_PATH,
            smoke_output_root=RESULTS_LIGANDMPNN_SMOKE_DIR,
            unique_sequences_path=LIGANDMPNN_UNIQUE_SEQUENCES_PATH,
            shortlist_path=LIGANDMPNN_SHORTLIST_PATH,
            report_path=LIGANDMPNN_SHORTLIST_REPORT_PATH,
            shortlist_fasta_path=LIGANDMPNN_SHORTLIST_FASTA_PATH,
            rosetta_shortlist_path=ROSETTA_SHORTLIST_PATH,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
