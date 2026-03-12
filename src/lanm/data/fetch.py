"""Manifest loading and public structure fetch helpers."""

from __future__ import annotations

import csv
from pathlib import Path

import requests

from lanm.filesystem import atomic_write_bytes
from lanm.models import FetchRecord, StructureManifestEntry
from lanm.paths import LOCAL_BUNDLE_DIR, LOCAL_STRUCTURE_MANIFEST_PATH, PUBLIC_STRUCTURES_DIR, REPO_ROOT


def require_path(path: Path) -> Path:
    """Return an existing path or raise with the exact missing path string."""
    if not path.exists():
        raise FileNotFoundError(str(path))
    return path


def load_structure_manifest(
    manifest_path: Path = LOCAL_STRUCTURE_MANIFEST_PATH,
) -> list[StructureManifestEntry]:
    """Load the normalized structure manifest."""
    require_path(manifest_path)
    with manifest_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        entries = [
            StructureManifestEntry(
                pdb_id=str(row["pdb_id"]).upper(),
                title=str(row["title"]),
                bound_metal=str(row["bound_metal"]),
                organism=str(row["organism"]),
                method=str(row["method"]),
                released=str(row["released"]),
                rcsb_url=str(row["rcsb_url"]),
                cif_download_url=str(row["cif_download_url"]),
            )
            for row in reader
        ]
    return sorted(entries, key=lambda item: item.pdb_id)


def _has_local_structure_assets(bundle_dir: Path, pdb_id: str) -> bool:
    stem = pdb_id.lower()
    candidates = (
        bundle_dir / f"{stem}_atoms.csv",
        bundle_dir / f"{stem}.txt",
        bundle_dir / f"{stem}.cif",
        bundle_dir / f"{stem}.pdb",
    )
    return any(path.exists() for path in candidates)


def fetch_public_structures(
    manifest_path: Path = LOCAL_STRUCTURE_MANIFEST_PATH,
    local_bundle_dir: Path = LOCAL_BUNDLE_DIR,
    public_structures_dir: Path = PUBLIC_STRUCTURES_DIR,
    timeout_seconds: int = 60,
) -> list[FetchRecord]:
    """Fetch CIFs for manifest entries that are not already present locally."""
    records: list[FetchRecord] = []
    for entry in load_structure_manifest(manifest_path):
        destination = public_structures_dir / f"{entry.pdb_id}.cif"
        if destination.exists():
            records.append(
                FetchRecord(
                    pdb_id=entry.pdb_id,
                    bound_metal=entry.bound_metal,
                    source_url=entry.cif_download_url,
                    local_path=str(destination.relative_to(REPO_ROOT)),
                    status="present",
                    size_bytes=destination.stat().st_size,
                    http_status=200,
                    message="File already present.",
                )
            )
            continue
        if _has_local_structure_assets(local_bundle_dir, entry.pdb_id):
            records.append(
                FetchRecord(
                    pdb_id=entry.pdb_id,
                    bound_metal=entry.bound_metal,
                    source_url=entry.cif_download_url,
                    local_path=str(destination.relative_to(REPO_ROOT)),
                    status="available_locally",
                    size_bytes=0,
                    http_status=0,
                    message="Local structure assets already cover this entry.",
                )
            )
            continue
        response = requests.get(entry.cif_download_url, timeout=timeout_seconds)
        if response.status_code != 200:
            raise RuntimeError(f"{entry.pdb_id} {response.status_code} {entry.cif_download_url}")
        atomic_write_bytes(destination, response.content)
        records.append(
            FetchRecord(
                pdb_id=entry.pdb_id,
                bound_metal=entry.bound_metal,
                source_url=entry.cif_download_url,
                local_path=str(destination.relative_to(REPO_ROOT)),
                status="downloaded",
                size_bytes=destination.stat().st_size,
                http_status=response.status_code,
                message="Downloaded from RCSB.",
            )
        )
    return records

