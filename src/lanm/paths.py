"""Repository path helpers."""

from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = REPO_ROOT / "config" / "project.yaml"
INCOMING_DIR = REPO_ROOT / "data" / "incoming"
LOCAL_BUNDLE_DIR = REPO_ROOT / "data" / "raw" / "local_bundle"
PUBLIC_DATA_DIR = REPO_ROOT / "data" / "raw" / "public"
PUBLIC_STRUCTURES_DIR = PUBLIC_DATA_DIR / "structures"
RESULTS_DIR = REPO_ROOT / "results"
RESULTS_TABLES_DIR = RESULTS_DIR / "tables"
RESULTS_REPORTS_DIR = RESULTS_DIR / "reports"
RESULTS_FIGURES_DIR = RESULTS_DIR / "figures"
DOCS_DIR = REPO_ROOT / "docs"
PROJECT_BRIEF_PATH = DOCS_DIR / "project_brief_extracted.md"


def ensure_runtime_directories() -> None:
    """Create the project directories used by the phase-1 CLIs."""
    for path in (
        LOCAL_BUNDLE_DIR,
        PUBLIC_STRUCTURES_DIR,
        RESULTS_TABLES_DIR,
        RESULTS_REPORTS_DIR,
        RESULTS_FIGURES_DIR,
        DOCS_DIR,
    ):
        path.mkdir(parents=True, exist_ok=True)
