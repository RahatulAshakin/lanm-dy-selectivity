from __future__ import annotations

import csv
import math
from pathlib import Path

from lanm.analysis.final_project_report import (
    FINALIST_CANDIDATE_IDS,
    FinalProjectReportInputs,
    FinalProjectReportOutputs,
    build_final_project_package,
    build_specific_pocket_thermodynamics_table,
)


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if fieldnames:
            writer.writeheader()
            writer.writerows(rows)


def _build_round2_summary_rows() -> list[dict[str, str]]:
    metals = {
        "Dy": ("stable_bound", "2.215379"),
        "Nd": ("stable_bound", "2.316352"),
        "Y": ("stable_bound", "2.194464"),
        "Al": ("persistent_capture_flag", "1.725170"),
        "Fe": ("persistent_capture_flag", "1.926443"),
    }
    overrides = {
        "hans_interface_ss_plus_if_u04_r2u01": {
            "Dy": ("stable_bound", "2.203839"),
            "Nd": ("stable_bound", "2.310945"),
            "Y": ("stable_bound", "2.190009"),
            "Al": ("persistent_capture_flag", "1.720473"),
            "Fe": ("persistent_capture_flag", "1.931271"),
        },
        "am1_mex_ss_only_u02_r2u02": {
            "Dy": ("stable_bound", "2.201751"),
            "Nd": ("stable_bound", "2.301707"),
            "Y": ("stable_bound", "2.185430"),
            "Al": ("persistent_capture_flag", "1.714569"),
            "Fe": ("persistent_capture_flag", "1.926456"),
        },
    }
    topology_by_candidate = {
        "hans_pocket_ss_only_u02_r2u01": "hans_monomer",
        "hans_interface_ss_plus_if_u04_r2u01": "hans_interface_multichain",
        "am1_mex_ss_only_u02_r2u02": "am1_monomer",
    }
    rows: list[dict[str, str]] = []
    for candidate_id in FINALIST_CANDIDATE_IDS:
        candidate_metals = overrides.get(candidate_id, metals)
        for metal, (status, distance) in candidate_metals.items():
            rows.append(
                {
                    "candidate_id": candidate_id,
                    "target_metal": metal,
                    "topology_class": topology_by_candidate[candidate_id],
                    "replicate_success_count": "3",
                    "median_mean_min_metal_oxygen_distance_A": distance,
                    "median_inner_sphere_occupancy_fraction": "1.0",
                    "screening_status": status,
                    "failure_reason": "",
                }
            )
    return rows


def _build_panel_status_rows() -> list[dict[str, str]]:
    return [
        {
            "candidate_id": "am1_mex_ss_only_u02_r2u02",
            "dy_status": "stable_bound",
            "nd_status": "stable_bound",
            "y_status": "stable_bound",
            "al_status": "persistent_capture_flag",
            "fe_status": "persistent_capture_flag",
            "candidate_keep_for_metadynamics": "False",
        },
        {
            "candidate_id": "hans_interface_ss_plus_if_u04_r2u01",
            "dy_status": "stable_bound",
            "nd_status": "stable_bound",
            "y_status": "stable_bound",
            "al_status": "persistent_capture_flag",
            "fe_status": "persistent_capture_flag",
            "candidate_keep_for_metadynamics": "False",
        },
        {
            "candidate_id": "hans_pocket_ss_only_u02_r2u01",
            "dy_status": "stable_bound",
            "nd_status": "stable_bound",
            "y_status": "stable_bound",
            "al_status": "persistent_capture_flag",
            "fe_status": "persistent_capture_flag",
            "candidate_keep_for_metadynamics": "False",
        },
    ]


def _build_finalist_rows() -> list[dict[str, str]]:
    return [
        {
            "panel_rank": "1",
            "panel_member_id": "am1_mex_ss_only_u02_r2u02",
            "panel_role": "top_am1_mex_round2_candidate",
            "candidate_id": "am1_mex_ss_only_u02_r2u02",
            "campaign_id": "am1_mex_ss_only",
            "backbone_id": "am1_mex_8fns_chain_a",
            "topology_class": "am1_monomer",
            "preserved_metal_identity": "ND",
            "rosetta_score_rank_within_topology": "1",
            "topology_candidate_count": "1",
            "representative_design_id": "1",
            "representative_sequence_id": "am1_mex_ss_only_u02_r2u02_design_01",
            "representative_ligand_confidence": "0.2213",
            "representative_overall_confidence": "0.3023",
            "total_score": "218.234",
            "representative_packed_pdb": "results/mock/am1_mex_ss_only_u02_r2u02.pdb",
            "selection_reason": "mock",
        },
        {
            "panel_rank": "2",
            "panel_member_id": "hans_pocket_ss_only_u02_r2u01",
            "panel_role": "top_hans_pocket_round2_candidate",
            "candidate_id": "hans_pocket_ss_only_u02_r2u01",
            "campaign_id": "hans_pocket_ss_only",
            "backbone_id": "hans_pocket_8fnr_chain_a",
            "topology_class": "hans_monomer",
            "preserved_metal_identity": "DY",
            "rosetta_score_rank_within_topology": "1",
            "topology_candidate_count": "1",
            "representative_design_id": "1",
            "representative_sequence_id": "hans_pocket_ss_only_u02_r2u01_design_01",
            "representative_ligand_confidence": "1.0",
            "representative_overall_confidence": "0.3399",
            "total_score": "144.07",
            "representative_packed_pdb": "results/mock/hans_pocket_ss_only_u02_r2u01.pdb",
            "selection_reason": "mock",
        },
        {
            "panel_rank": "3",
            "panel_member_id": "hans_interface_ss_plus_if_u04_r2u01",
            "panel_role": "best_hans_interface_aware_round2_candidate",
            "candidate_id": "hans_interface_ss_plus_if_u04_r2u01",
            "campaign_id": "hans_interface_ss_plus_if",
            "backbone_id": "hans_interface_8fnr_a_b_c_d",
            "topology_class": "hans_interface_multichain",
            "preserved_metal_identity": "DY",
            "rosetta_score_rank_within_topology": "1",
            "topology_candidate_count": "2",
            "representative_design_id": "2",
            "representative_sequence_id": "hans_interface_ss_plus_if_u04_r2u01_design_02",
            "representative_ligand_confidence": "0.3454",
            "representative_overall_confidence": "0.3303",
            "total_score": "530.712",
            "representative_packed_pdb": "results/mock/hans_interface_ss_plus_if_u04_r2u01.pdb",
            "selection_reason": "mock",
        },
    ]


def test_build_specific_pocket_thermodynamics_table_ranks_expected_main_pocket() -> None:
    pocket_rows = build_specific_pocket_thermodynamics_table(
        tuple(_build_finalist_rows()),
        tuple(_build_panel_status_rows()),
    )

    assert [row["candidate_id"] for row in pocket_rows] == [
        "hans_pocket_ss_only_u02_r2u01",
        "hans_interface_ss_plus_if_u04_r2u01",
        "am1_mex_ss_only_u02_r2u02",
    ]
    assert pocket_rows[0]["qm_priority_note"].startswith("Main Quantum Pocket;")
    assert math.isclose(float(pocket_rows[0]["pocket_confidence"]), 0.80197, rel_tol=0.0, abs_tol=1e-6)
    assert math.isclose(float(pocket_rows[0]["deltaG_pocket_proxy_kcal_mol"]), 0.130645, rel_tol=0.0, abs_tol=1e-6)


def test_build_final_project_package_writes_all_required_outputs(tmp_path: Path) -> None:
    inputs_dir = tmp_path / "inputs"
    outputs_dir = tmp_path / "outputs"

    _write_csv(
        inputs_dir / "template_chain_summary.csv",
        [
            {"template_id": "6MI5", "source_type": "cif", "chain_id": "X", "metal_identity": "Y", "metal_site_count": "3"},
            {"template_id": "8FNS", "source_type": "csv_atom_table", "chain_id": "A", "metal_identity": "ND", "metal_site_count": "4"},
            {"template_id": "8DQ2", "source_type": "csv_atom_table", "chain_id": "A", "metal_identity": "LA", "metal_site_count": "3"},
            {"template_id": "8FNR", "source_type": "cif", "chain_id": "A", "metal_identity": "DY", "metal_site_count": "4"},
        ],
    )
    _write_csv(
        inputs_dir / "residue_role_map.csv",
        [
            {"template_id": "6MI5", "role": "first_shell"},
            {"template_id": "6MI5", "role": "second_sphere"},
            {"template_id": "8FNS", "role": "first_shell"},
            {"template_id": "8DQ2", "role": "interchain_contact"},
            {"template_id": "8FNR", "role": "first_shell"},
            {"template_id": "8FNR", "role": "solvent_contact"},
        ],
    )
    _write_csv(
        inputs_dir / "design_mask_candidates.csv",
        [
            {"canonical_family_position": "1", "mutable_second_sphere": "False", "mutable_interface": "False"},
            {"canonical_family_position": "18", "mutable_second_sphere": "True", "mutable_interface": "False"},
            {"canonical_family_position": "27", "mutable_second_sphere": "False", "mutable_interface": "True"},
            {"canonical_family_position": "38", "mutable_second_sphere": "True", "mutable_interface": "True"},
        ],
    )
    _write_csv(
        inputs_dir / "proteinmpnn_shortlist.csv",
        [
            {"shortlist_rank": "1", "candidate_id": "am1_mex_ss_only_u02", "campaign_id": "am1_mex_ss_only", "backbone_id": "am1_mex_8fns_chain_a", "mutation_string": "A38:D>N"},
            {"shortlist_rank": "2", "candidate_id": "hans_pocket_ss_only_u02", "campaign_id": "hans_pocket_ss_only", "backbone_id": "hans_pocket_8fnr_chain_a", "mutation_string": "A44:A>D"},
            {"shortlist_rank": "3", "candidate_id": "hans_interface_ss_plus_if_u04", "campaign_id": "hans_interface_ss_plus_if", "backbone_id": "hans_interface_8fnr_a_b_c_d", "mutation_string": "A35:S>T"},
            {"shortlist_rank": "4", "candidate_id": "hans_interface_ss_plus_if_u02", "campaign_id": "hans_interface_ss_plus_if", "backbone_id": "hans_interface_8fnr_a_b_c_d", "mutation_string": "A35:S>T,A44:A>D"},
        ],
    )
    _write_csv(
        inputs_dir / "ligandmpnn_shortlist.csv",
        [
            {"candidate_id": "am1_mex_ss_only_u02", "mutation_string": "T12S"},
            {"candidate_id": "hans_pocket_ss_only_u02", "mutation_string": "T16D"},
            {"candidate_id": "hans_interface_ss_plus_if_u04", "mutation_string": "D21A"},
            {"candidate_id": "hans_interface_ss_plus_if_u02", "mutation_string": "D21P"},
        ],
    )
    _write_csv(
        inputs_dir / "integrated_candidate_ranking.csv",
        [
            {"candidate_id": "am1_mex_ss_only_u02"},
            {"candidate_id": "hans_pocket_ss_only_u02"},
            {"candidate_id": "hans_interface_ss_plus_if_u04"},
            {"candidate_id": "hans_interface_ss_plus_if_u02"},
        ],
    )
    _write_csv(
        inputs_dir / "md_validation_panel.csv",
        [
            {"candidate_id": "am1_mex_ss_only_u02", "panel_member_type": "designed_candidate"},
            {"candidate_id": "hans_pocket_ss_only_u02", "panel_member_type": "designed_candidate"},
            {"candidate_id": "hans_interface_ss_plus_if_u04", "panel_member_type": "designed_candidate"},
            {"candidate_id": "hans_interface_ss_plus_if_u02", "panel_member_type": "designed_candidate"},
        ],
    )
    _write_csv(
        inputs_dir / "proteinmpnn_round2_shortlist.csv",
        [
            {"shortlist_rank": "1", "candidate_id": "am1_mex_ss_only_u02_r2u02", "seed_candidate_id": "am1_mex_ss_only_u02", "backbone_id": "am1_mex_8fns_chain_a", "topology_class": "am1_monomer", "mutation_string": "A38:D>N,A56:D>K"},
            {"shortlist_rank": "2", "candidate_id": "hans_pocket_ss_only_u02_r2u01", "seed_candidate_id": "hans_pocket_ss_only_u02", "backbone_id": "hans_pocket_8fnr_chain_a", "topology_class": "hans_monomer", "mutation_string": "A44:A>D,A55:T>K"},
            {"shortlist_rank": "3", "candidate_id": "hans_interface_ss_plus_if_u04_r2u01", "seed_candidate_id": "hans_interface_ss_plus_if_u04", "backbone_id": "hans_interface_8fnr_a_b_c_d", "topology_class": "hans_interface_multichain", "mutation_string": "A35:S>T,A44:A>D"},
            {"shortlist_rank": "4", "candidate_id": "hans_interface_ss_plus_if_u04_r2u02", "seed_candidate_id": "hans_interface_ss_plus_if_u04", "backbone_id": "hans_interface_8fnr_a_b_c_d", "topology_class": "hans_interface_multichain", "mutation_string": "A35:S>T"},
        ],
    )
    _write_csv(
        inputs_dir / "ligandmpnn_round2_smoke_summary.csv",
        [
            {"candidate_id": "am1_mex_ss_only_u02_r2u02"},
            {"candidate_id": "hans_pocket_ss_only_u02_r2u01"},
            {"candidate_id": "hans_interface_ss_plus_if_u04_r2u01"},
            {"candidate_id": "hans_interface_ss_plus_if_u04_r2u02"},
        ],
    )
    _write_csv(
        inputs_dir / "ligandmpnn_round2_shortlist.csv",
        [
            {"candidate_id": "am1_mex_ss_only_u02_r2u02", "mutation_string": "A38:D>N,A56:D>K,A93:K>E"},
            {"candidate_id": "hans_pocket_ss_only_u02_r2u01", "mutation_string": "A44:A>D,A55:T>W,A92:M>L,A97:K>A"},
            {"candidate_id": "hans_interface_ss_plus_if_u04_r2u01", "mutation_string": "A35:S>T,A44:A>D,A55:W>V,A115:E>D"},
            {"candidate_id": "hans_interface_ss_plus_if_u04_r2u02", "mutation_string": "A100:R>K"},
        ],
    )
    _write_csv(
        inputs_dir / "rosetta_round2_candidate_ranking.csv",
        [
            {"candidate_id": "am1_mex_ss_only_u02_r2u02"},
            {"candidate_id": "hans_pocket_ss_only_u02_r2u01"},
            {"candidate_id": "hans_interface_ss_plus_if_u04_r2u01"},
            {"candidate_id": "hans_interface_ss_plus_if_u04_r2u02"},
        ],
    )
    _write_csv(inputs_dir / "round2_md_rescreen_panel.csv", _build_finalist_rows())
    _write_csv(inputs_dir / "round2_openmm_screening_summary.csv", _build_round2_summary_rows())
    _write_csv(inputs_dir / "round2_openmm_screening_panel_status.csv", _build_panel_status_rows())

    result = build_final_project_package(
        inputs=FinalProjectReportInputs(
            template_chain_summary_path=inputs_dir / "template_chain_summary.csv",
            residue_role_map_path=inputs_dir / "residue_role_map.csv",
            design_mask_candidates_path=inputs_dir / "design_mask_candidates.csv",
            proteinmpnn_shortlist_path=inputs_dir / "proteinmpnn_shortlist.csv",
            ligandmpnn_shortlist_path=inputs_dir / "ligandmpnn_shortlist.csv",
            integrated_candidate_ranking_path=inputs_dir / "integrated_candidate_ranking.csv",
            md_validation_panel_path=inputs_dir / "md_validation_panel.csv",
            proteinmpnn_round2_shortlist_path=inputs_dir / "proteinmpnn_round2_shortlist.csv",
            ligandmpnn_round2_smoke_summary_path=inputs_dir / "ligandmpnn_round2_smoke_summary.csv",
            ligandmpnn_round2_shortlist_path=inputs_dir / "ligandmpnn_round2_shortlist.csv",
            rosetta_round2_candidate_ranking_path=inputs_dir / "rosetta_round2_candidate_ranking.csv",
            round2_md_rescreen_panel_path=inputs_dir / "round2_md_rescreen_panel.csv",
            round2_openmm_screening_summary_path=inputs_dir / "round2_openmm_screening_summary.csv",
            round2_openmm_screening_panel_status_path=inputs_dir / "round2_openmm_screening_panel_status.csv",
        ),
        outputs=FinalProjectReportOutputs(
            final_project_report_path=outputs_dir / "results" / "reports" / "final_project_report.md",
            project_summary_report_path=outputs_dir / "results" / "reports" / "project_summary_report.md",
            project_detailed_report_path=outputs_dir / "results" / "reports" / "project_detailed_report.md",
            paper_short_path=outputs_dir / "docs" / "paper_short.md",
            paper_long_path=outputs_dir / "docs" / "paper_long.md",
            readme_path=outputs_dir / "README.md",
            final_table_1_path=outputs_dir / "results" / "tables" / "final_table_1_template_summary.csv",
            final_table_2_path=outputs_dir / "results" / "tables" / "final_table_2_candidate_progression.csv",
            final_table_3_path=outputs_dir / "results" / "tables" / "final_table_3_specific_pocket_thermodynamics.csv",
            final_figure_1_path=outputs_dir / "results" / "figures" / "final_figure_1_template_harmonization.png",
            final_figure_2_path=outputs_dir / "results" / "figures" / "final_figure_2_design_funnel.png",
            final_figure_3_path=outputs_dir / "results" / "figures" / "final_figure_3_multimetal_validation.png",
            final_figure_4_path=outputs_dir / "results" / "figures" / "final_figure_4_quantum_pocket_summary.png",
        ),
    )

    assert result.main_quantum_pocket_candidate_id == "hans_pocket_ss_only_u02_r2u01"
    assert all(path.exists() for path in result.output_paths)
    assert "Rare-earth retention was preserved across tested candidates." in (outputs_dir / "results" / "reports" / "final_project_report.md").read_text(encoding="utf-8")
    assert "No candidate advanced to full QM/QCT in the validated pipeline." in (outputs_dir / "docs" / "paper_short.md").read_text(encoding="utf-8")
    assert "The Main Quantum Pocket ranking is a deterministic next-step prioritization based on completed workflow outputs." in (outputs_dir / "docs" / "paper_long.md").read_text(encoding="utf-8")
    assert (outputs_dir / "results" / "tables" / "final_table_3_specific_pocket_thermodynamics.csv").read_text(encoding="utf-8").splitlines()[1].startswith("1,hans_pocket_ss_only_u02_r2u01")
