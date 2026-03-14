from pathlib import Path

from lanm.analysis.md_panel_selection import (
    IntegratedCandidateRow,
    LigandCandidateContext,
    RosettaScoreRankingRow,
    build_integrated_candidate_rows,
    normalize_rosetta_scores_by_backbone,
    resolve_ligand_candidate_contexts,
    select_md_validation_panel,
)
from lanm.models import LigandMPNNSmokeSummaryRow, ProteinMPNNShortlistRow


def _make_protein_row(
    *,
    shortlist_rank: int,
    candidate_id: str,
    campaign_id: str,
    backbone_id: str,
) -> ProteinMPNNShortlistRow:
    return ProteinMPNNShortlistRow(
        shortlist_rank=shortlist_rank,
        candidate_id=candidate_id,
        campaign_id=campaign_id,
        backbone_id=backbone_id,
        design_set_name="test_campaign",
        designed_sequence="AAAA",
        representative_sequence_id=f"{candidate_id}_seq",
        occurrence_count=1,
        campaign_rank=1,
        temperature=0.1,
        best_score=0.1,
        best_global_score=0.2,
        best_seq_recovery=0.3,
        mutation_count=1,
        mutation_string="A1V",
        canonical_family_positions_mutated="1",
        am1_mature_positions_mutated="1",
        includes_second_sphere_position=True,
        includes_interface_position=False,
        retention_reason="test",
    )


def _make_smoke_row(
    *,
    shortlist_rank: int,
    candidate_id: str,
    campaign_id: str,
    backbone_id: str,
    preserved_metal_identity: str,
    mean_ligand_confidence: float,
    max_ligand_confidence: float,
    mean_overall_confidence: float,
    max_overall_confidence: float,
    sequence_length: int,
) -> LigandMPNNSmokeSummaryRow:
    return LigandMPNNSmokeSummaryRow(
        shortlist_rank=shortlist_rank,
        candidate_id=candidate_id,
        campaign_id=campaign_id,
        backbone_id=backbone_id,
        source_structure_id="8FNS",
        designed_chains="A",
        preserved_metal_identity=preserved_metal_identity,
        redesigned_residues="A35",
        generated_sequence_count=2,
        unique_sequence_count=1,
        unique_sequence_fraction=0.5,
        sequence_length=sequence_length,
        mean_overall_confidence=mean_overall_confidence,
        min_overall_confidence=mean_overall_confidence - 0.01,
        max_overall_confidence=max_overall_confidence,
        mean_ligand_confidence=mean_ligand_confidence,
        min_ligand_confidence=mean_ligand_confidence - 0.01,
        max_ligand_confidence=max_ligand_confidence,
        mean_pairwise_identity=1.0,
        min_pairwise_identity=1.0,
        max_pairwise_identity=1.0,
        temperature=0.1,
        seed=37,
        output_fasta_path=f"results/ligandmpnn_smoke/{candidate_id}/seqs/{candidate_id}.fa",
    )


def _make_rosetta_row(
    *,
    rosetta_rank: int,
    shortlist_rank: int,
    candidate_id: str,
    campaign_id: str,
    backbone_id: str,
    preserved_metal_identity: str,
    ligand_confidence: float,
    overall_confidence: float,
    total_score: float,
) -> RosettaScoreRankingRow:
    return RosettaScoreRankingRow(
        rosetta_rank=rosetta_rank,
        shortlist_rank=shortlist_rank,
        candidate_id=candidate_id,
        campaign_id=campaign_id,
        backbone_id=backbone_id,
        preserved_metal_identity=preserved_metal_identity,
        representative_design_id=1,
        representative_ligand_confidence=ligand_confidence,
        representative_overall_confidence=overall_confidence,
        representative_packed_pdb=f"results/ligandmpnn_smoke/{candidate_id}/packed/{candidate_id}_packed_1_1.pdb",
        scorefile_path=f"results/rosetta_score_smoke/{candidate_id}/score.sc",
        total_score=total_score,
    )


def _make_ligand_context(
    *,
    shortlist_rank: int,
    candidate_id: str,
    campaign_id: str,
    backbone_id: str,
    preserved_metal_identity: str,
    ligand_confidence: float,
    overall_confidence: float,
    source: str = "ligandmpnn_shortlist.csv",
) -> LigandCandidateContext:
    return LigandCandidateContext(
        shortlist_rank=shortlist_rank,
        candidate_id=candidate_id,
        campaign_id=campaign_id,
        backbone_id=backbone_id,
        preserved_metal_identity=preserved_metal_identity,
        representative_design_id=1,
        representative_ligand_confidence=ligand_confidence,
        representative_overall_confidence=overall_confidence,
        representative_packed_pdb=f"results/ligandmpnn_smoke/{candidate_id}/packed/{candidate_id}_packed_1_1.pdb",
        source=source,
    )


def _make_integrated_row(
    *,
    integrated_rank: int,
    candidate_id: str,
    campaign_id: str,
    backbone_id: str,
    topology_class: str,
    proteinmpnn_rank: int,
    ligand_confidence: float,
    rosetta_rank_within_topology: int,
    rosetta_total_score: float,
) -> IntegratedCandidateRow:
    return IntegratedCandidateRow(
        integrated_rank=integrated_rank,
        candidate_id=candidate_id,
        campaign_id=campaign_id,
        backbone_id=backbone_id,
        topology_class=topology_class,
        preserved_metal_identity="ND" if topology_class == "am1_monomer" else "DY",
        proteinmpnn_rank=proteinmpnn_rank,
        ligandmpnn_shortlist_rank=integrated_rank,
        ligandmpnn_representative_design_id=1,
        ligandmpnn_representative_ligand_confidence=ligand_confidence,
        ligandmpnn_representative_overall_confidence=ligand_confidence + 0.01,
        ligandmpnn_mean_ligand_confidence=ligand_confidence - 0.01,
        ligandmpnn_max_ligand_confidence=ligand_confidence,
        ligandmpnn_mean_overall_confidence=ligand_confidence,
        ligandmpnn_max_overall_confidence=ligand_confidence + 0.01,
        rosetta_global_rank=integrated_rank,
        rosetta_score_rank_within_topology=rosetta_rank_within_topology,
        rosetta_total_score=rosetta_total_score,
        rosetta_total_score_delta_from_backbone_best=0.0 if rosetta_rank_within_topology == 1 else 1.0,
        rosetta_total_score_normalized_within_backbone=0.0 if rosetta_rank_within_topology == 1 else 1.0,
        integrated_priority_score=float(integrated_rank) / 10.0,
        direct_cross_topology_comparison_valid=False,
        cross_topology_comparison_note="test note",
        representative_packed_pdb=f"results/ligandmpnn_smoke/{candidate_id}/packed/{candidate_id}_packed_1_1.pdb",
        ligand_context_source="ligandmpnn_shortlist.csv",
    )


def test_normalize_rosetta_scores_by_backbone_separates_backbone_scales() -> None:
    rows = (
        _make_rosetta_row(
            rosetta_rank=1,
            shortlist_rank=1,
            candidate_id="am1_best",
            campaign_id="am1_mex_ss_only",
            backbone_id="am1_mex_8fns_chain_a",
            preserved_metal_identity="ND",
            ligand_confidence=0.30,
            overall_confidence=0.31,
            total_score=100.0,
        ),
        _make_rosetta_row(
            rosetta_rank=2,
            shortlist_rank=2,
            candidate_id="am1_worse",
            campaign_id="am1_mex_ss_only",
            backbone_id="am1_mex_8fns_chain_a",
            preserved_metal_identity="ND",
            ligand_confidence=0.29,
            overall_confidence=0.30,
            total_score=120.0,
        ),
        _make_rosetta_row(
            rosetta_rank=3,
            shortlist_rank=3,
            candidate_id="interface_best",
            campaign_id="hans_interface_ss_plus_if",
            backbone_id="hans_interface_8fnr_a_b_c_d",
            preserved_metal_identity="DY",
            ligand_confidence=0.45,
            overall_confidence=0.46,
            total_score=500.0,
        ),
        _make_rosetta_row(
            rosetta_rank=4,
            shortlist_rank=4,
            candidate_id="interface_worse",
            campaign_id="hans_interface_ss_plus_if",
            backbone_id="hans_interface_8fnr_a_b_c_d",
            preserved_metal_identity="DY",
            ligand_confidence=0.20,
            overall_confidence=0.21,
            total_score=800.0,
        ),
    )

    normalized = normalize_rosetta_scores_by_backbone(rows)

    assert normalized["am1_best"]["rank_within_topology"] == 1
    assert normalized["am1_best"]["normalized_within_backbone"] == 0.0
    assert normalized["am1_worse"]["normalized_within_backbone"] == 1.0
    assert normalized["interface_best"]["rank_within_topology"] == 1
    assert normalized["interface_best"]["normalized_within_backbone"] == 0.0
    assert normalized["interface_worse"]["normalized_within_backbone"] == 1.0


def test_resolve_ligand_candidate_contexts_falls_back_to_rosetta_metadata(tmp_path: Path) -> None:
    rosetta_rows = (
        _make_rosetta_row(
            rosetta_rank=1,
            shortlist_rank=4,
            candidate_id="am1_best",
            campaign_id="am1_mex_ss_only",
            backbone_id="am1_mex_8fns_chain_a",
            preserved_metal_identity="ND",
            ligand_confidence=0.35,
            overall_confidence=0.36,
            total_score=111.0,
        ),
    )

    contexts = resolve_ligand_candidate_contexts(
        ligandmpnn_shortlist_path=tmp_path / "missing_ligandmpnn_shortlist.csv",
        rosetta_rows=rosetta_rows,
    )

    assert len(contexts) == 1
    assert contexts[0].source == "inferred_from_ligandmpnn_smoke_and_rosetta_score"
    assert contexts[0].representative_design_id == 1
    assert contexts[0].representative_ligand_confidence == 0.35
    assert contexts[0].shortlist_rank == 4


def test_build_integrated_candidate_rows_uses_within_backbone_rosetta_not_raw_cross_topology() -> None:
    protein_rows = (
        _make_protein_row(
            shortlist_rank=1,
            candidate_id="am1_best",
            campaign_id="am1_mex_ss_only",
            backbone_id="am1_mex_8fns_chain_a",
        ),
        _make_protein_row(
            shortlist_rank=2,
            candidate_id="am1_other",
            campaign_id="am1_mex_ss_only",
            backbone_id="am1_mex_8fns_chain_a",
        ),
        _make_protein_row(
            shortlist_rank=4,
            candidate_id="interface_best",
            campaign_id="hans_interface_ss_plus_if",
            backbone_id="hans_interface_8fnr_a_b_c_d",
        ),
        _make_protein_row(
            shortlist_rank=5,
            candidate_id="interface_other",
            campaign_id="hans_interface_ss_plus_if",
            backbone_id="hans_interface_8fnr_a_b_c_d",
        ),
    )
    ligand_rows = (
        _make_ligand_context(
            shortlist_rank=1,
            candidate_id="am1_best",
            campaign_id="am1_mex_ss_only",
            backbone_id="am1_mex_8fns_chain_a",
            preserved_metal_identity="ND",
            ligand_confidence=0.30,
            overall_confidence=0.31,
        ),
        _make_ligand_context(
            shortlist_rank=2,
            candidate_id="am1_other",
            campaign_id="am1_mex_ss_only",
            backbone_id="am1_mex_8fns_chain_a",
            preserved_metal_identity="ND",
            ligand_confidence=0.28,
            overall_confidence=0.29,
        ),
        _make_ligand_context(
            shortlist_rank=3,
            candidate_id="interface_best",
            campaign_id="hans_interface_ss_plus_if",
            backbone_id="hans_interface_8fnr_a_b_c_d",
            preserved_metal_identity="DY",
            ligand_confidence=0.45,
            overall_confidence=0.46,
        ),
        _make_ligand_context(
            shortlist_rank=4,
            candidate_id="interface_other",
            campaign_id="hans_interface_ss_plus_if",
            backbone_id="hans_interface_8fnr_a_b_c_d",
            preserved_metal_identity="DY",
            ligand_confidence=0.20,
            overall_confidence=0.21,
        ),
    )
    smoke_rows = (
        _make_smoke_row(
            shortlist_rank=1,
            candidate_id="am1_best",
            campaign_id="am1_mex_ss_only",
            backbone_id="am1_mex_8fns_chain_a",
            preserved_metal_identity="ND",
            mean_ligand_confidence=0.29,
            max_ligand_confidence=0.30,
            mean_overall_confidence=0.30,
            max_overall_confidence=0.31,
            sequence_length=105,
        ),
        _make_smoke_row(
            shortlist_rank=2,
            candidate_id="am1_other",
            campaign_id="am1_mex_ss_only",
            backbone_id="am1_mex_8fns_chain_a",
            preserved_metal_identity="ND",
            mean_ligand_confidence=0.27,
            max_ligand_confidence=0.28,
            mean_overall_confidence=0.28,
            max_overall_confidence=0.29,
            sequence_length=105,
        ),
        _make_smoke_row(
            shortlist_rank=4,
            candidate_id="interface_best",
            campaign_id="hans_interface_ss_plus_if",
            backbone_id="hans_interface_8fnr_a_b_c_d",
            preserved_metal_identity="DY",
            mean_ligand_confidence=0.44,
            max_ligand_confidence=0.45,
            mean_overall_confidence=0.45,
            max_overall_confidence=0.46,
            sequence_length=435,
        ),
        _make_smoke_row(
            shortlist_rank=5,
            candidate_id="interface_other",
            campaign_id="hans_interface_ss_plus_if",
            backbone_id="hans_interface_8fnr_a_b_c_d",
            preserved_metal_identity="DY",
            mean_ligand_confidence=0.19,
            max_ligand_confidence=0.20,
            mean_overall_confidence=0.20,
            max_overall_confidence=0.21,
            sequence_length=435,
        ),
    )
    rosetta_rows = (
        _make_rosetta_row(
            rosetta_rank=1,
            shortlist_rank=1,
            candidate_id="am1_best",
            campaign_id="am1_mex_ss_only",
            backbone_id="am1_mex_8fns_chain_a",
            preserved_metal_identity="ND",
            ligand_confidence=0.30,
            overall_confidence=0.31,
            total_score=100.0,
        ),
        _make_rosetta_row(
            rosetta_rank=2,
            shortlist_rank=2,
            candidate_id="am1_other",
            campaign_id="am1_mex_ss_only",
            backbone_id="am1_mex_8fns_chain_a",
            preserved_metal_identity="ND",
            ligand_confidence=0.28,
            overall_confidence=0.29,
            total_score=120.0,
        ),
        _make_rosetta_row(
            rosetta_rank=3,
            shortlist_rank=3,
            candidate_id="interface_best",
            campaign_id="hans_interface_ss_plus_if",
            backbone_id="hans_interface_8fnr_a_b_c_d",
            preserved_metal_identity="DY",
            ligand_confidence=0.45,
            overall_confidence=0.46,
            total_score=500.0,
        ),
        _make_rosetta_row(
            rosetta_rank=4,
            shortlist_rank=4,
            candidate_id="interface_other",
            campaign_id="hans_interface_ss_plus_if",
            backbone_id="hans_interface_8fnr_a_b_c_d",
            preserved_metal_identity="DY",
            ligand_confidence=0.20,
            overall_confidence=0.21,
            total_score=800.0,
        ),
    )

    integrated = build_integrated_candidate_rows(
        protein_rows=protein_rows,
        ligand_context_rows=ligand_rows,
        smoke_rows=smoke_rows,
        rosetta_rows=rosetta_rows,
    )

    assert [row.candidate_id for row in integrated] == [
        "am1_best",
        "interface_best",
        "am1_other",
        "interface_other",
    ]
    assert integrated[1].rosetta_total_score == 500.0
    assert integrated[1].rosetta_score_rank_within_topology == 1
    assert integrated[1].direct_cross_topology_comparison_valid is False


def test_select_md_validation_panel_enforces_constraints_and_reference_baselines() -> None:
    integrated_rows = (
        _make_integrated_row(
            integrated_rank=1,
            candidate_id="pocket_top",
            campaign_id="hans_pocket_ss_only",
            backbone_id="hans_pocket_8fnr_chain_a",
            topology_class="hans_monomer",
            proteinmpnn_rank=3,
            ligand_confidence=0.43,
            rosetta_rank_within_topology=1,
            rosetta_total_score=150.0,
        ),
        _make_integrated_row(
            integrated_rank=2,
            candidate_id="am1_top",
            campaign_id="am1_mex_ss_only",
            backbone_id="am1_mex_8fns_chain_a",
            topology_class="am1_monomer",
            proteinmpnn_rank=1,
            ligand_confidence=0.35,
            rosetta_rank_within_topology=1,
            rosetta_total_score=111.0,
        ),
        _make_integrated_row(
            integrated_rank=3,
            candidate_id="interface_top",
            campaign_id="hans_interface_ss_plus_if",
            backbone_id="hans_interface_8fnr_a_b_c_d",
            topology_class="hans_interface_multichain",
            proteinmpnn_rank=2,
            ligand_confidence=0.39,
            rosetta_rank_within_topology=1,
            rosetta_total_score=531.0,
        ),
        _make_integrated_row(
            integrated_rank=4,
            candidate_id="pocket_other",
            campaign_id="hans_pocket_ss_only",
            backbone_id="hans_pocket_8fnr_chain_a",
            topology_class="hans_monomer",
            proteinmpnn_rank=5,
            ligand_confidence=0.42,
            rosetta_rank_within_topology=2,
            rosetta_total_score=156.0,
        ),
        _make_integrated_row(
            integrated_rank=5,
            candidate_id="interface_other",
            campaign_id="hans_interface_if_only",
            backbone_id="hans_interface_8fnr_a_b_c_d",
            topology_class="hans_interface_multichain",
            proteinmpnn_rank=6,
            ligand_confidence=0.24,
            rosetta_rank_within_topology=2,
            rosetta_total_score=604.0,
        ),
        _make_integrated_row(
            integrated_rank=6,
            candidate_id="am1_other",
            campaign_id="am1_mex_ss_only",
            backbone_id="am1_mex_8fns_chain_a",
            topology_class="am1_monomer",
            proteinmpnn_rank=7,
            ligand_confidence=0.33,
            rosetta_rank_within_topology=2,
            rosetta_total_score=112.0,
        ),
    )

    panel = select_md_validation_panel(integrated_rows)

    assert len(panel) == 6
    assert {row.panel_role for row in panel} == {
        "am1_mex_designed_candidate",
        "am1_mex_wild_type_reference",
        "hans_pocket_designed_candidate",
        "hans_interface_aware_candidate",
        "hans_wild_type_reference",
        "best_remaining_designed_candidate",
    }
    assert [row.candidate_id for row in panel if row.panel_role == "am1_mex_designed_candidate"] == ["am1_top"]
    assert [row.candidate_id for row in panel if row.panel_role == "hans_pocket_designed_candidate"] == ["pocket_top"]
    assert [row.candidate_id for row in panel if row.panel_role == "hans_interface_aware_candidate"] == ["interface_top"]
    assert [row.candidate_id for row in panel if row.panel_role == "best_remaining_designed_candidate"] == ["pocket_other"]
    assert [row.backbone_id for row in panel if row.panel_role == "hans_wild_type_reference"] == [
        "hans_interface_8fnr_a_b_c_d"
    ]
