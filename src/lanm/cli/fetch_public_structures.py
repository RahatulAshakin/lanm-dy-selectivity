"""Fetch missing public CIF structures referenced by the local manifest."""

from __future__ import annotations

from lanm.data.fetch import fetch_public_structures
from lanm.data.normalize import normalize_incoming_bundle
from lanm.filesystem import write_csv_rows
from lanm.logging_utils import configure_logging
from lanm.paths import PUBLIC_STRUCTURE_FETCH_LOG_PATH, ensure_runtime_directories


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    try:
        normalize_incoming_bundle()
        records = fetch_public_structures()
        write_csv_rows(PUBLIC_STRUCTURE_FETCH_LOG_PATH, records)
    except FileNotFoundError as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()

