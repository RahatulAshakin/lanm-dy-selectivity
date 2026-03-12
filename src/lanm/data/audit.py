"""Dataset audit helpers for the normalized local bundle."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from typing import Iterable

from lanm.configuration import ProjectConfig
from lanm.models import BriefExtractionResult, FileSummary, NormalizationResult
from lanm.paths import LOCAL_BUNDLE_DIR, REPO_ROOT

TEXTUAL_SUFFIXES = {".txt", ".md", ".cif"}
BINARY_SUFFIXES = {".docx", ".zip"}


def _infer_file_type(path: Path) -> str:
    suffix = path.suffix.lower()
    return {
        ".csv": "csv",
        ".json": "json",
        ".txt": "text",
        ".docx": "docx",
        ".md": "markdown",
        ".cif": "cif",
    }.get(suffix, suffix.lstrip(".") or "unknown")


def _display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def _should_ignore_bundle_file(path: Path) -> bool:
    return "zone.identifier" in path.name.lower()


def _summarize_csv(path: Path) -> FileSummary:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    counter: Counter[str] = Counter()
    chains: set[str] = set()
    residue_numbers_by_chain: dict[str, list[int]] = {}
    hetatm_count = 0
    for row in rows:
        for fieldname in fieldnames:
            if not str(row.get(fieldname, "")).strip():
                counter[fieldname] += 1
        chain_id = str(row.get("chain_id", "")).strip()
        residue_value = str(row.get("residue_seq", "")).strip()
        if chain_id:
            chains.add(chain_id)
        if chain_id and residue_value:
            try:
                residue_numbers_by_chain.setdefault(chain_id, []).append(int(residue_value))
            except ValueError:
                pass
        if str(row.get("record_type", "")).strip().upper() == "HETATM":
            hetatm_count += 1
    residue_ranges: list[str] = []
    for chain_id in sorted(residue_numbers_by_chain):
        values = residue_numbers_by_chain[chain_id]
        residue_ranges.append(f"{chain_id}:{min(values)}-{max(values)}")
    missing_fields = ";".join(
        f"{field}:{counter[field]}"
        for field in fieldnames
        if counter[field] > 0
    )
    return FileSummary(
        relative_path=_display_path(path),
        file_type="csv",
        file_size_bytes=path.stat().st_size,
        row_count=len(rows),
        column_names=";".join(fieldnames),
        chain_ids=";".join(sorted(chains)),
        residue_ranges=";".join(residue_ranges),
        hetatm_count=hetatm_count,
        missing_fields=missing_fields,
    )


def _summarize_json(path: Path) -> FileSummary:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(payload, list):
        row_count = len(payload)
    elif isinstance(payload, dict):
        row_count = len(payload)
    else:
        row_count = 1
    return FileSummary(
        relative_path=_display_path(path),
        file_type="json",
        file_size_bytes=path.stat().st_size,
        row_count=row_count,
        column_names="",
        chain_ids="",
        residue_ranges="",
        hetatm_count=0,
        missing_fields="",
    )


def _summarize_binary(path: Path) -> FileSummary:
    return FileSummary(
        relative_path=_display_path(path),
        file_type=_infer_file_type(path),
        file_size_bytes=path.stat().st_size,
        row_count=0,
        column_names="",
        chain_ids="",
        residue_ranges="",
        hetatm_count=0,
        missing_fields="",
    )


def _summarize_textual(path: Path) -> FileSummary:
    try:
        line_count = len(path.read_text(encoding="utf-8").splitlines())
    except UnicodeDecodeError:
        return _summarize_binary(path)
    return FileSummary(
        relative_path=_display_path(path),
        file_type=_infer_file_type(path),
        file_size_bytes=path.stat().st_size,
        row_count=line_count,
        column_names="",
        chain_ids="",
        residue_ranges="",
        hetatm_count=0,
        missing_fields="",
    )


def summarize_local_bundle(bundle_dir: Path = LOCAL_BUNDLE_DIR) -> list[FileSummary]:
    summaries: list[FileSummary] = []
    for path in sorted(
        item for item in bundle_dir.iterdir()
        if item.is_file() and not _should_ignore_bundle_file(item)
    ):
        suffix = path.suffix.lower()
        if suffix == ".csv":
            summaries.append(_summarize_csv(path))
        elif suffix == ".json":
            summaries.append(_summarize_json(path))
        elif suffix in BINARY_SUFFIXES:
            summaries.append(_summarize_binary(path))
        elif suffix in TEXTUAL_SUFFIXES:
            summaries.append(_summarize_textual(path))
        else:
            summaries.append(_summarize_textual(path))
    return summaries


def render_dataset_audit_report(
    config: ProjectConfig,
    normalization: NormalizationResult,
    brief: BriefExtractionResult,
    inventory: Iterable[FileSummary],
) -> str:
    inventory_rows = list(inventory)
    total_bytes = sum(row.file_size_bytes for row in inventory_rows)
    report_lines = [
        "# Dataset Audit",
        "",
        "## Configuration",
        f"- project_name: `{config.project_name}`",
        f"- target_metal: `{config.target_metal}`",
        f"- competitors: `{', '.join(config.competitors)}`",
        f"- first_shell_cutoff_A: `{config.first_shell_cutoff_A}`",
        f"- second_sphere_cutoff_A: `{config.second_sphere_cutoff_A}`",
        "",
        "## Normalization",
        f"- normalized_file_count: `{len(normalization.normalized_paths)}`",
        f"- normalized_paths: `{', '.join(normalization.normalized_paths)}`",
        "",
        "## Project Brief Extraction",
        f"- status: `{brief.status}`",
        f"- source_path: `{brief.source_path or 'not found'}`",
        f"- output_path: `{brief.output_path}`",
        f"- note: {brief.note}",
        "",
        "## Inventory Summary",
        f"- file_count: `{len(inventory_rows)}`",
        f"- total_size_bytes: `{total_bytes}`",
        "",
        "| relative_path | file_type | rows | size_bytes | chains | residue_ranges | hetatm_count |",
        "| --- | --- | ---: | ---: | --- | --- | ---: |",
    ]
    for row in inventory_rows:
        report_lines.append(
            f"| {row.relative_path} | {row.file_type} | {row.row_count} | "
            f"{row.file_size_bytes} | {row.chain_ids or '-'} | "
            f"{row.residue_ranges or '-'} | {row.hetatm_count} |"
        )
    report_lines.append("")
    return "\n".join(report_lines)
