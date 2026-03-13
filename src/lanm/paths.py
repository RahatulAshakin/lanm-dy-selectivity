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
RESULTS_DESIGN_INPUTS_DIR = RESULTS_DIR / "design_inputs"
PROTEINMPNN_INPUTS_DIR = RESULTS_DESIGN_INPUTS_DIR / "proteinmpnn"
PROTEINMPNN_BACKBONES_DIR = PROTEINMPNN_INPUTS_DIR / "backbones"
PROTEINMPNN_CAMPAIGNS_DIR = PROTEINMPNN_INPUTS_DIR / "campaigns"
DOCS_DIR = REPO_ROOT / "docs"
PROJECT_BRIEF_PATH = DOCS_DIR / "project_brief_extracted.md"
LOCAL_STRUCTURE_MANIFEST_PATH = LOCAL_BUNDLE_DIR / "lanmodulin_lanthanide_structures.csv"
LANMODULIN_SEQUENCES_PATH = LOCAL_BUNDLE_DIR / "lanmodulin_sequences.csv"
DATASET_INVENTORY_PATH = RESULTS_TABLES_DIR / "dataset_inventory.csv"
DATASET_AUDIT_PATH = RESULTS_REPORTS_DIR / "dataset_audit.md"
PUBLIC_STRUCTURE_FETCH_LOG_PATH = RESULTS_TABLES_DIR / "public_structure_fetch_log.csv"
METAL_SITE_SUMMARY_PATH = RESULTS_TABLES_DIR / "metal_site_summary.csv"
SHELL_ANNOTATION_PATH = RESULTS_TABLES_DIR / "shell_annotation.csv"
SITE_OVERVIEW_FIGURE_PATH = RESULTS_FIGURES_DIR / "template_site_overview.png"
TEMPLATE_HARMONIZATION_REPORT_PATH = RESULTS_REPORTS_DIR / "template_harmonization.md"
TEMPLATE_CHAIN_SUMMARY_PATH = RESULTS_TABLES_DIR / "template_chain_summary.csv"
TEMPLATE_SITE_SUMMARY_PATH = RESULTS_TABLES_DIR / "template_site_summary.csv"
CROSS_TEMPLATE_RESIDUE_ALIGNMENT_PATH = RESULTS_TABLES_DIR / "cross_template_residue_alignment.csv"
RESIDUE_ROLE_MAP_PATH = RESULTS_TABLES_DIR / "residue_role_map.csv"
DESIGN_MASK_CANDIDATES_PATH = RESULTS_TABLES_DIR / "design_mask_candidates.csv"
TEMPLATE_HARMONIZATION_FIGURE_PATH = RESULTS_FIGURES_DIR / "template_harmonization_overview.png"
DESIGN_MASKS_PATH = REPO_ROOT / "config" / "design_masks.yaml"
DESIGN_CAMPAIGNS_PATH = REPO_ROOT / "config" / "design_campaigns.yaml"
DESIGN_CAMPAIGNS_REPORT_PATH = RESULTS_REPORTS_DIR / "design_campaigns.md"
DESIGN_CAMPAIGN_POSITIONS_PATH = RESULTS_TABLES_DIR / "design_campaign_positions.csv"
DESIGN_BACKBONE_MANIFEST_PATH = RESULTS_TABLES_DIR / "design_backbone_manifest.csv"
DESIGN_CAMPAIGN_MANIFEST_PATH = RESULTS_TABLES_DIR / "design_campaign_manifest.csv"


def ensure_runtime_directories() -> None:
    """Create the project directories used by the phase-1 CLIs."""
    for path in (
        LOCAL_BUNDLE_DIR,
        PUBLIC_STRUCTURES_DIR,
        RESULTS_TABLES_DIR,
        RESULTS_REPORTS_DIR,
        RESULTS_FIGURES_DIR,
        RESULTS_DESIGN_INPUTS_DIR,
        PROTEINMPNN_INPUTS_DIR,
        PROTEINMPNN_BACKBONES_DIR,
        PROTEINMPNN_CAMPAIGNS_DIR,
        DOCS_DIR,
    ):
        path.mkdir(parents=True, exist_ok=True)
