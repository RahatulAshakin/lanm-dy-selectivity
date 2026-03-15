"""Prepare deterministic Phase 7C LigandMPNN round-2 inputs."""

from __future__ import annotations

from lanm.analysis.ligandmpnn_round2_inputs import prepare_ligandmpnn_round2_inputs
from lanm.logging_utils import configure_logging
from lanm.paths import (
    LIGANDMPNN_ROUND2_INPUTS_CONFIG_PATH,
    LIGANDMPNN_ROUND2_INPUTS_DIR,
    LIGANDMPNN_ROUND2_INPUTS_REPORT_PATH,
    LIGANDMPNN_ROUND2_INPUT_MANIFEST_PATH,
    LIGANDMPNN_ROUND2_REDESIGN_POSITIONS_PATH,
    PROTEINMPNN_ROUND2_CAMPAIGNS_DIR,
    PROTEINMPNN_ROUND2_SEQUENCE_CATALOG_PATH,
    PROTEINMPNN_ROUND2_SHORTLIST_PATH,
    PROTEINMPNN_ROUND2_SHORTLIST_REPORT_PATH,
    PROTEINMPNN_ROUND2_SUMMARY_PATH,
    PROTEINMPNN_ROUND2_UNIQUE_SEQUENCES_PATH,
    REDESIGN_ROUND2_CONFIG_PATH,
    REDESIGN_ROUND2_POSITIONS_PATH,
    REDESIGN_ROUND2_SEEDS_PATH,
    RESULTS_PROTEINMPNN_ROUND2_SMOKE_DIR,
    ensure_runtime_directories,
)


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    try:
        prepare_ligandmpnn_round2_inputs(
            redesign_round2_config_path=REDESIGN_ROUND2_CONFIG_PATH,
            round2_seed_table_path=REDESIGN_ROUND2_SEEDS_PATH,
            round2_position_table_path=REDESIGN_ROUND2_POSITIONS_PATH,
            round2_campaigns_dir=PROTEINMPNN_ROUND2_CAMPAIGNS_DIR,
            round2_smoke_output_root=RESULTS_PROTEINMPNN_ROUND2_SMOKE_DIR,
            round2_summary_path=PROTEINMPNN_ROUND2_SUMMARY_PATH,
            round2_sequence_catalog_path=PROTEINMPNN_ROUND2_SEQUENCE_CATALOG_PATH,
            unique_sequences_path=PROTEINMPNN_ROUND2_UNIQUE_SEQUENCES_PATH,
            shortlist_path=PROTEINMPNN_ROUND2_SHORTLIST_PATH,
            shortlist_report_path=PROTEINMPNN_ROUND2_SHORTLIST_REPORT_PATH,
            ligandmpnn_round2_root=LIGANDMPNN_ROUND2_INPUTS_DIR,
            manifest_path=LIGANDMPNN_ROUND2_INPUT_MANIFEST_PATH,
            redesign_positions_path=LIGANDMPNN_ROUND2_REDESIGN_POSITIONS_PATH,
            inputs_report_path=LIGANDMPNN_ROUND2_INPUTS_REPORT_PATH,
            config_path=LIGANDMPNN_ROUND2_INPUTS_CONFIG_PATH,
        )
    except FileNotFoundError as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
