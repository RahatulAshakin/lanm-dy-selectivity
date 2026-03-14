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
RESULTS_FASTA_DIR = RESULTS_DIR / "fasta"
RESULTS_PROTEINMPNN_SMOKE_DIR = RESULTS_DIR / "proteinmpnn_smoke"
RESULTS_LIGANDMPNN_SMOKE_DIR = RESULTS_DIR / "ligandmpnn_smoke"
RESULTS_ROSETTA_SMOKE_DIR = RESULTS_DIR / "rosetta_smoke"
RESULTS_ROSETTA_SCORE_SMOKE_DIR = RESULTS_DIR / "rosetta_score_smoke"
RESULTS_ROSETTA_RELAX_SMOKE_DIR = RESULTS_DIR / "rosetta_relax_smoke"
RESULTS_DESIGN_INPUTS_DIR = RESULTS_DIR / "design_inputs"
RESULTS_MD_INPUTS_DIR = RESULTS_DIR / "md_inputs"
RESULTS_MD_SYSTEM_BUILD_DIR = RESULTS_DIR / "md_system_build"
RESULTS_OPENMM_SYSTEM_BUILD_DIR = RESULTS_DIR / "openmm_system_build"
RESULTS_OPENMM_SYSTEM_SMOKE_DIR = RESULTS_DIR / "openmm_system_smoke"
RESULTS_OPENMM_EQUILIBRATION_SMOKE_DIR = RESULTS_DIR / "openmm_equilibration_smoke"
PROTEINMPNN_INPUTS_DIR = RESULTS_DESIGN_INPUTS_DIR / "proteinmpnn"
LIGANDMPNN_INPUTS_DIR = RESULTS_DESIGN_INPUTS_DIR / "ligandmpnn"
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
PROTEINMPNN_SMOKE_SUMMARY_PATH = RESULTS_TABLES_DIR / "proteinmpnn_smoke_summary.csv"
PROTEINMPNN_SEQUENCE_CATALOG_PATH = RESULTS_TABLES_DIR / "proteinmpnn_sequence_catalog.csv"
PROTEINMPNN_SMOKE_REPORT_PATH = RESULTS_REPORTS_DIR / "proteinmpnn_smoke.md"
PROTEINMPNN_UNIQUE_SEQUENCES_PATH = RESULTS_TABLES_DIR / "proteinmpnn_unique_sequences.csv"
PROTEINMPNN_SHORTLIST_PATH = RESULTS_TABLES_DIR / "proteinmpnn_shortlist.csv"
PROTEINMPNN_SHORTLIST_REPORT_PATH = RESULTS_REPORTS_DIR / "proteinmpnn_shortlist.md"
PROTEINMPNN_SHORTLIST_FASTA_PATH = RESULTS_FASTA_DIR / "proteinmpnn_shortlist.fa"
LIGANDMPNN_INPUT_MANIFEST_PATH = RESULTS_TABLES_DIR / "ligandmpnn_input_manifest.csv"
LIGANDMPNN_REDESIGN_POSITIONS_PATH = RESULTS_TABLES_DIR / "ligandmpnn_redesign_positions.csv"
LIGANDMPNN_INPUTS_REPORT_PATH = RESULTS_REPORTS_DIR / "ligandmpnn_inputs.md"
LIGANDMPNN_INPUTS_CONFIG_PATH = REPO_ROOT / "config" / "ligandmpnn_inputs.yaml"
LIGANDMPNN_SMOKE_SUMMARY_PATH = RESULTS_TABLES_DIR / "ligandmpnn_smoke_summary.csv"
LIGANDMPNN_SEQUENCE_CATALOG_PATH = RESULTS_TABLES_DIR / "ligandmpnn_sequence_catalog.csv"
LIGANDMPNN_SMOKE_REPORT_PATH = RESULTS_REPORTS_DIR / "ligandmpnn_smoke.md"
LIGANDMPNN_UNIQUE_SEQUENCES_PATH = RESULTS_TABLES_DIR / "ligandmpnn_unique_sequences.csv"
LIGANDMPNN_SHORTLIST_PATH = RESULTS_TABLES_DIR / "ligandmpnn_shortlist.csv"
LIGANDMPNN_SHORTLIST_REPORT_PATH = RESULTS_REPORTS_DIR / "ligandmpnn_shortlist.md"
LIGANDMPNN_SHORTLIST_FASTA_PATH = RESULTS_FASTA_DIR / "ligandmpnn_shortlist.fa"
ROSETTA_SHORTLIST_PATH = REPO_ROOT / "config" / "rosetta_shortlist.yaml"
ROSETTA_SMOKE_SUMMARY_PATH = RESULTS_TABLES_DIR / "rosetta_smoke_summary.csv"
ROSETTA_CANDIDATE_RANKING_PATH = RESULTS_TABLES_DIR / "rosetta_candidate_ranking.csv"
ROSETTA_SMOKE_REPORT_PATH = RESULTS_REPORTS_DIR / "rosetta_smoke.md"
ROSETTA_SCORE_SMOKE_SUMMARY_PATH = RESULTS_TABLES_DIR / "rosetta_score_smoke_summary.csv"
ROSETTA_SCORE_CANDIDATE_RANKING_PATH = RESULTS_TABLES_DIR / "rosetta_score_candidate_ranking.csv"
ROSETTA_SCORE_SMOKE_REPORT_PATH = RESULTS_REPORTS_DIR / "rosetta_score_smoke.md"
ROSETTA_RELAX_SMOKE_SUMMARY_PATH = RESULTS_TABLES_DIR / "rosetta_relax_smoke_summary.csv"
ROSETTA_RELAX_CANDIDATE_RANKING_PATH = RESULTS_TABLES_DIR / "rosetta_relax_candidate_ranking.csv"
ROSETTA_RELAX_SMOKE_REPORT_PATH = RESULTS_REPORTS_DIR / "rosetta_relax_smoke.md"
INTEGRATED_CANDIDATE_RANKING_PATH = RESULTS_TABLES_DIR / "integrated_candidate_ranking.csv"
MD_VALIDATION_PANEL_PATH = RESULTS_TABLES_DIR / "md_validation_panel.csv"
MD_PANEL_SELECTION_REPORT_PATH = RESULTS_REPORTS_DIR / "md_panel_selection.md"
MD_PANEL_CONFIG_PATH = REPO_ROOT / "config" / "md_panel.yaml"
MD_INPUT_MANIFEST_PATH = RESULTS_TABLES_DIR / "md_input_manifest.csv"
MD_INPUT_PREPARATION_REPORT_PATH = RESULTS_REPORTS_DIR / "md_input_preparation.md"
MD_PROTOCOL_CONFIG_PATH = REPO_ROOT / "config" / "md_protocol.yaml"
METAL_MODEL_REGISTRY_PATH = REPO_ROOT / "config" / "metal_models.yaml"
METAL_PARAMETER_CONFIG_PATH = REPO_ROOT / "config" / "metal_parameters.yaml"
METAL_PARAMETER_VALUES_PATH = REPO_ROOT / "config" / "metal_parameter_values.yaml"
MD_SYSTEM_BUILD_MANIFEST_PATH = RESULTS_TABLES_DIR / "md_system_build_manifest.csv"
MD_SYSTEM_BUILDING_PLAN_REPORT_PATH = RESULTS_REPORTS_DIR / "md_system_building_plan.md"
METAL_PARAMETER_MAPPING_PATH = RESULTS_TABLES_DIR / "metal_parameter_mapping.csv"
METAL_BUILD_READINESS_PATH = RESULTS_TABLES_DIR / "metal_build_readiness.csv"
METAL_PARAMETER_STRATEGY_REPORT_PATH = RESULTS_REPORTS_DIR / "metal_parameter_strategy.md"
OPENMM_SYSTEM_SMOKE_SUMMARY_PATH = RESULTS_TABLES_DIR / "openmm_system_smoke_summary.csv"
OPENMM_SYSTEM_SMOKE_REPORT_PATH = RESULTS_REPORTS_DIR / "openmm_system_smoke.md"
OPENMM_SYSTEM_BUILD_SUMMARY_PATH = RESULTS_TABLES_DIR / "openmm_system_build_summary.csv"
OPENMM_SYSTEM_BUILD_REPORT_PATH = RESULTS_REPORTS_DIR / "openmm_system_building.md"
OPENMM_EQUILIBRATION_SMOKE_SUMMARY_PATH = RESULTS_TABLES_DIR / "openmm_equilibration_smoke_summary.csv"
OPENMM_EQUILIBRATION_SMOKE_REPORT_PATH = RESULTS_REPORTS_DIR / "openmm_equilibration_smoke.md"


def ensure_runtime_directories() -> None:
    """Create the project directories used by the phase-1 CLIs."""
    for path in (
        LOCAL_BUNDLE_DIR,
        PUBLIC_STRUCTURES_DIR,
        RESULTS_TABLES_DIR,
        RESULTS_REPORTS_DIR,
        RESULTS_FIGURES_DIR,
        RESULTS_FASTA_DIR,
        RESULTS_PROTEINMPNN_SMOKE_DIR,
        RESULTS_LIGANDMPNN_SMOKE_DIR,
        RESULTS_ROSETTA_SMOKE_DIR,
        RESULTS_ROSETTA_SCORE_SMOKE_DIR,
        RESULTS_ROSETTA_RELAX_SMOKE_DIR,
        RESULTS_DESIGN_INPUTS_DIR,
        RESULTS_MD_INPUTS_DIR,
        RESULTS_MD_SYSTEM_BUILD_DIR,
        RESULTS_OPENMM_SYSTEM_BUILD_DIR,
        RESULTS_OPENMM_SYSTEM_SMOKE_DIR,
        RESULTS_OPENMM_EQUILIBRATION_SMOKE_DIR,
        PROTEINMPNN_INPUTS_DIR,
        LIGANDMPNN_INPUTS_DIR,
        PROTEINMPNN_BACKBONES_DIR,
        PROTEINMPNN_CAMPAIGNS_DIR,
        DOCS_DIR,
    ):
        path.mkdir(parents=True, exist_ok=True)
