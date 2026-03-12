"""Deterministic normalization for the local Phase 1 data bundle."""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

from lanm.filesystem import copy_if_changed
from lanm.models import NormalizationResult
from lanm.paths import INCOMING_DIR, LOCAL_BUNDLE_DIR, REPO_ROOT


def require_path(path: Path) -> Path:
    """Return an existing path or raise with the exact missing path string."""
    if not path.exists():
        raise FileNotFoundError(str(path))
    return path


def snake_case_filename(filename: str) -> str:
    """Normalize an arbitrary filename to a deterministic snake_case form."""
    path = Path(filename)
    stem = re.sub(r"[^0-9A-Za-z]+", "_", path.stem).strip("_").lower()
    suffix = path.suffix.lower()
    return f"{stem}{suffix}"


def should_skip_normalization(path: Path) -> bool:
    """Ignore Windows ADS sidecar names that should never enter the bundle."""
    return "zone.identifier" in path.name.lower()


def _display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def _resolve_destination_names(sources: list[Path]) -> dict[Path, str]:
    counts: defaultdict[str, int] = defaultdict(int)
    resolved: dict[Path, str] = {}
    for source in sorted(sources, key=lambda item: item.name.lower()):
        normalized = snake_case_filename(source.name)
        counts[normalized] += 1
        if counts[normalized] == 1:
            resolved[source] = normalized
            continue
        path = Path(normalized)
        resolved[source] = f"{path.stem}_{counts[normalized]}{path.suffix}"
    return resolved


def normalize_incoming_bundle(
    incoming_dir: Path = INCOMING_DIR,
    local_bundle_dir: Path = LOCAL_BUNDLE_DIR,
) -> NormalizationResult:
    """Copy incoming files into the raw bundle with stable filenames."""
    require_path(incoming_dir)
    if local_bundle_dir.exists():
        for path in sorted(
            (path for path in local_bundle_dir.iterdir() if path.is_file() and should_skip_normalization(path)),
            key=lambda item: item.name.lower(),
        ):
            path.unlink()
    source_files = sorted(
        path
        for path in incoming_dir.iterdir()
        if path.is_file() and not should_skip_normalization(path)
    )
    skipped_sidecars = tuple(
        _display_path(path)
        for path in sorted(
            (path for path in incoming_dir.iterdir() if path.is_file() and should_skip_normalization(path)),
            key=lambda item: item.name.lower(),
        )
    )
    destination_names = _resolve_destination_names(source_files)
    normalized_paths: list[str] = []
    for source in source_files:
        destination = local_bundle_dir / destination_names[source]
        copy_if_changed(source, destination)
        normalized_paths.append(_display_path(destination))
    return NormalizationResult(
        normalized_paths=tuple(normalized_paths),
        skipped_sidecars=skipped_sidecars,
    )
