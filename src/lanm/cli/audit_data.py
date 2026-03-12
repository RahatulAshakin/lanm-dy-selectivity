"""Normalize incoming data and write the Phase 1 dataset audit."""

from __future__ import annotations

from lanm.configuration import load_project_config
from lanm.data.audit import render_dataset_audit_report, summarize_local_bundle
from lanm.data.docx_extract import extract_project_brief
from lanm.data.normalize import normalize_incoming_bundle
from lanm.filesystem import atomic_write_text, write_csv_rows
from lanm.logging_utils import configure_logging
from lanm.paths import (
    DATASET_AUDIT_PATH,
    DATASET_INVENTORY_PATH,
    ensure_runtime_directories,
)


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    try:
        normalization = normalize_incoming_bundle()
        brief = extract_project_brief()
        config = load_project_config()
        inventory = summarize_local_bundle()
        write_csv_rows(DATASET_INVENTORY_PATH, inventory)
        atomic_write_text(
            DATASET_AUDIT_PATH,
            render_dataset_audit_report(config, normalization, brief, inventory),
        )
    except FileNotFoundError as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()

