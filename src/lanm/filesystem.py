"""Filesystem helpers for deterministic writes."""

from __future__ import annotations

import csv
import json
import shutil
from dataclasses import asdict, is_dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Iterable

try:
    import yaml  # type: ignore
except ModuleNotFoundError:  # pragma: no cover - optional dependency
    yaml = None


def copy_if_changed(source: Path, destination: Path) -> bool:
    """Copy a file only when content differs or the destination is missing."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and source.read_bytes() == destination.read_bytes():
        return False
    shutil.copy2(source, destination)
    return True


def atomic_write_text(path: Path, text: str) -> None:
    """Write text atomically for idempotent CLI outputs."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("w", encoding="utf-8", delete=False, dir=path.parent) as handle:
        handle.write(text)
        temp_path = Path(handle.name)
    temp_path.replace(path)


def atomic_write_bytes(path: Path, payload: bytes) -> None:
    """Write bytes atomically for idempotent binary outputs."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("wb", delete=False, dir=path.parent) as handle:
        handle.write(payload)
        temp_path = Path(handle.name)
    temp_path.replace(path)


def write_csv_rows(path: Path, rows: Iterable[object]) -> None:
    """Write dataclass or dict rows to CSV deterministically."""
    row_list: list[dict[str, object]] = []
    for row in rows:
        if is_dataclass(row):
            row_list.append(asdict(row))
        elif isinstance(row, dict):
            row_list.append(row)
        else:
            raise TypeError(f"Unsupported CSV row type: {type(row)!r}")
    fieldnames = list(row_list[0].keys()) if row_list else []
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("w", encoding="utf-8", newline="", delete=False, dir=path.parent) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if fieldnames:
            writer.writeheader()
            writer.writerows(row_list)
        temp_path = Path(handle.name)
    temp_path.replace(path)


def write_json(path: Path, payload: object) -> None:
    """Write JSON atomically."""
    atomic_write_text(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def write_jsonl(path: Path, rows: Iterable[object]) -> None:
    """Write newline-delimited JSON rows atomically."""
    payload_lines: list[str] = []
    for row in rows:
        if is_dataclass(row):
            serializable = asdict(row)
        else:
            serializable = row
        payload_lines.append(json.dumps(serializable, sort_keys=True))
    atomic_write_text(path, "\n".join(payload_lines) + ("\n" if payload_lines else ""))


def write_yaml(path: Path, payload: object) -> None:
    """Write YAML atomically, falling back to JSON-compatible YAML."""
    if yaml is not None:
        atomic_write_text(path, yaml.safe_dump(payload, sort_keys=False))
        return
    atomic_write_text(path, json.dumps(payload, indent=2) + "\n")
