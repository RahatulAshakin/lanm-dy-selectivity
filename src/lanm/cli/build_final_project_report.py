"""Build the deterministic final reporting-only project package."""

from __future__ import annotations

from lanm.analysis.final_project_report import (
    FinalProjectReportInputs,
    FinalProjectReportOutputs,
    build_final_project_package,
)
from lanm.logging_utils import configure_logging
from lanm.paths import (
    DESIGN_MASK_CANDIDATES_PATH,
    FINAL_FIGURE_1_TEMPLATE_HARMONIZATION_PATH,
    FINAL_FIGURE_2_DESIGN_FUNNEL_PATH,
    FINAL_FIGURE_3_MULTIMETAL_VALIDATION_PATH,
    FINAL_FIGURE_4_QUANTUM_POCKET_SUMMARY_PATH,
    FINAL_PROJECT_REPORT_PATH,
    FINAL_TABLE_1_TEMPLATE_SUMMARY_PATH,
    FINAL_TABLE_2_CANDIDATE_PROGRESSION_PATH,
    FINAL_TABLE_3_SPECIFIC_POCKET_THERMODYNAMICS_PATH,
    INTEGRATED_CANDIDATE_RANKING_PATH,
    LIGANDMPNN_ROUND2_SHORTLIST_PATH,
    LIGANDMPNN_ROUND2_SMOKE_SUMMARY_PATH,
    LIGANDMPNN_SHORTLIST_PATH,
    PAPER_LONG_PATH,
    PAPER_SHORT_PATH,
    PROJECT_DETAILED_REPORT_PATH,
    PROJECT_README_PATH,
    PROJECT_SUMMARY_REPORT_PATH,
    PROTEINMPNN_ROUND2_SHORTLIST_PATH,
    PROTEINMPNN_SHORTLIST_PATH,
    RESIDUE_ROLE_MAP_PATH,
    ROSETTA_ROUND2_CANDIDATE_RANKING_PATH,
    ROUND2_MD_RESCREEN_PANEL_PATH,
    ROUND2_OPENMM_SCREENING_PANEL_STATUS_PATH,
    ROUND2_OPENMM_SCREENING_SUMMARY_PATH,
    TEMPLATE_CHAIN_SUMMARY_PATH,
    MD_VALIDATION_PANEL_PATH,
    ensure_runtime_directories,
)


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    try:
        build_final_project_package(
            inputs=FinalProjectReportInputs(
                template_chain_summary_path=TEMPLATE_CHAIN_SUMMARY_PATH,
                residue_role_map_path=RESIDUE_ROLE_MAP_PATH,
                design_mask_candidates_path=DESIGN_MASK_CANDIDATES_PATH,
                proteinmpnn_shortlist_path=PROTEINMPNN_SHORTLIST_PATH,
                ligandmpnn_shortlist_path=LIGANDMPNN_SHORTLIST_PATH,
                integrated_candidate_ranking_path=INTEGRATED_CANDIDATE_RANKING_PATH,
                md_validation_panel_path=MD_VALIDATION_PANEL_PATH,
                proteinmpnn_round2_shortlist_path=PROTEINMPNN_ROUND2_SHORTLIST_PATH,
                ligandmpnn_round2_smoke_summary_path=LIGANDMPNN_ROUND2_SMOKE_SUMMARY_PATH,
                ligandmpnn_round2_shortlist_path=LIGANDMPNN_ROUND2_SHORTLIST_PATH,
                rosetta_round2_candidate_ranking_path=ROSETTA_ROUND2_CANDIDATE_RANKING_PATH,
                round2_md_rescreen_panel_path=ROUND2_MD_RESCREEN_PANEL_PATH,
                round2_openmm_screening_summary_path=ROUND2_OPENMM_SCREENING_SUMMARY_PATH,
                round2_openmm_screening_panel_status_path=ROUND2_OPENMM_SCREENING_PANEL_STATUS_PATH,
            ),
            outputs=FinalProjectReportOutputs(
                final_project_report_path=FINAL_PROJECT_REPORT_PATH,
                project_summary_report_path=PROJECT_SUMMARY_REPORT_PATH,
                project_detailed_report_path=PROJECT_DETAILED_REPORT_PATH,
                paper_short_path=PAPER_SHORT_PATH,
                paper_long_path=PAPER_LONG_PATH,
                readme_path=PROJECT_README_PATH,
                final_table_1_path=FINAL_TABLE_1_TEMPLATE_SUMMARY_PATH,
                final_table_2_path=FINAL_TABLE_2_CANDIDATE_PROGRESSION_PATH,
                final_table_3_path=FINAL_TABLE_3_SPECIFIC_POCKET_THERMODYNAMICS_PATH,
                final_figure_1_path=FINAL_FIGURE_1_TEMPLATE_HARMONIZATION_PATH,
                final_figure_2_path=FINAL_FIGURE_2_DESIGN_FUNNEL_PATH,
                final_figure_3_path=FINAL_FIGURE_3_MULTIMETAL_VALIDATION_PATH,
                final_figure_4_path=FINAL_FIGURE_4_QUANTUM_POCKET_SUMMARY_PATH,
            ),
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
