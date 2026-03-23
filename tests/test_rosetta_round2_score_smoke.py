from pathlib import Path

import pytest

from lanm.analysis.rosetta_round2_score_smoke import (
    build_round2_rosetta_score_command,
    choose_round2_representative_packed_pdb,
    parse_round2_rosetta_scorefile,
    select_round2_md_rescreen_panel,
)
from lanm.models import LigandMPNNRound2SequenceCatalogRow


def _make_sequence_row(
    *,
    candidate_id: str = "example_candidate",
    topology_class: str = "am1_monomer",
    design_id: int,
    ligand_confidence: float,
    overall_confidence: float,
    packed_pdb_path: str,
) -> LigandMPNNRound2SequenceCatalogRow:
    return LigandMPNNRound2SequenceCatalogRow(
        shortlist_rank=1,
        candidate_id=candidate_id,
        seed_rank=1,
        seed_candidate_id="seed_candidate",
        campaign_id="example_campaign",
        backbone_id="example_backbone",
        topology_class=topology_class,
        designed_chains="A",
        fixed_context_chains="",
        preserved_metal_identity="ND",
        redesigned_residues="A35,A44",
        design_id=design_id,
        temperature=0.1,
        seed=37,
        overall_confidence=overall_confidence,
        ligand_confidence=ligand_confidence,
        seq_recovery=0.5,
        sequence_id=f"{candidate_id}_design_{design_id:02d}",
        designed_sequence="AAAA",
        sequence_length=4,
        output_fasta_path=f"results/ligandmpnn_round2_smoke/{candidate_id}/seqs/{candidate_id}.fa",
        backbone_pdb_path=f"results/ligandmpnn_round2_smoke/{candidate_id}/backbones/{candidate_id}_{design_id}.pdb",
        packed_pdb_path=packed_pdb_path,
    )


def _make_ranking_row(
    *,
    shortlist_rank: int,
    candidate_id: str,
    campaign_id: str,
    topology_class: str,
    rosetta_score_rank_within_topology: int,
    representative_ligand_confidence: float,
    representative_overall_confidence: float = 0.5,
    total_score: float = 100.0,
) -> dict[str, object]:
    return {
        "rosetta_score_rank_within_topology": rosetta_score_rank_within_topology,
        "topology_candidate_count": 2 if topology_class == "hans_interface_multichain" else 1,
        "topology_best_total_score": total_score,
        "total_score_delta_from_topology_best": 0.0,
        "cross_topology_total_score_comparison_valid": False,
        "cross_topology_comparison_note": "within-topology only",
        "shortlist_rank": shortlist_rank,
        "shortlist_sequence_id": f"{candidate_id}_design_01",
        "candidate_id": candidate_id,
        "seed_rank": 1,
        "seed_candidate_id": "seed_candidate",
        "campaign_id": campaign_id,
        "backbone_id": f"{candidate_id}_backbone",
        "topology_class": topology_class,
        "preserved_metal_identity": "DY",
        "shortlist_design_id": 1,
        "shortlist_retention_reason": "retained",
        "representative_design_id": 1,
        "representative_sequence_id": f"{candidate_id}_design_01",
        "representative_ligand_confidence": representative_ligand_confidence,
        "representative_overall_confidence": representative_overall_confidence,
        "representative_packed_pdb": f"results/rosetta_round2_score_smoke/{candidate_id}/{candidate_id}.pdb",
        "scorefile_path": f"results/rosetta_round2_score_smoke/{candidate_id}/score.sc",
        "total_score": total_score,
        "score_description": f"{candidate_id}_0001",
    }


def test_choose_round2_representative_packed_pdb_prefers_ligand_then_overall_then_smallest_design_id() -> None:
    rows = (
        _make_sequence_row(
            design_id=3,
            ligand_confidence=0.75,
            overall_confidence=0.95,
            packed_pdb_path="results/ligandmpnn_round2_smoke/example_candidate/packed/example_candidate_packed_3_1.pdb",
        ),
        _make_sequence_row(
            design_id=2,
            ligand_confidence=0.8,
            overall_confidence=0.7,
            packed_pdb_path="results/ligandmpnn_round2_smoke/example_candidate/packed/example_candidate_packed_2_1.pdb",
        ),
        _make_sequence_row(
            design_id=1,
            ligand_confidence=0.8,
            overall_confidence=0.7,
            packed_pdb_path="results/ligandmpnn_round2_smoke/example_candidate/packed/example_candidate_packed_1_1.pdb",
        ),
    )

    representative = choose_round2_representative_packed_pdb(rows)

    assert representative.design_id == 1
    assert representative.packed_pdb_path.endswith("example_candidate_packed_1_1.pdb")


def test_build_round2_rosetta_score_command() -> None:
    command = build_round2_rosetta_score_command(
        rosetta_executable=Path("/opt/rosetta/source/bin/score_jd2.default.linuxgccrelease"),
        rosetta_database=Path("/opt/rosetta/database"),
        input_pdb=Path("/tmp/example_input.pdb"),
        scorefile_path=Path("/tmp/score.sc"),
    )

    assert command == (
        "/opt/rosetta/source/bin/score_jd2.default.linuxgccrelease",
        "-database",
        "/opt/rosetta/database",
        "-in:file:fullatom",
        "-s",
        "/tmp/example_input.pdb",
        "-out:file:scorefile",
        "/tmp/score.sc",
        "-out:overwrite",
        "-nstruct",
        "1",
        "-ignore_waters",
        "false",
        "-in:auto_setup_metals",
        "true",
        "-run:preserve_header",
        "true",
        "-out:file:renumber_pdb",
        "false",
        "-run:constant_seed",
        "true",
        "-run:jran",
        "37",
        "-score:set_weights",
        "fa_elec",
        "0.0",
    )


def test_parse_round2_rosetta_scorefile_reads_single_score_row(tmp_path: Path) -> None:
    scorefile_path = tmp_path / "score.sc"
    scorefile_path.write_text(
        "\n".join(
            [
                "SEQUENCE:",
                "SCORE: total_score fa_atr fa_rep fa_sol description",
                "SCORE: -123.450 -456.700 89.100 12.340 example_input_0001",
                "",
            ]
        ),
        encoding="utf-8",
    )

    parsed = parse_round2_rosetta_scorefile(scorefile_path)

    assert parsed == {
        "total_score": pytest.approx(-123.45),
        "fa_atr": pytest.approx(-456.7),
        "fa_rep": pytest.approx(89.1),
        "fa_sol": pytest.approx(12.34),
        "description": "example_input_0001",
    }


def test_select_round2_md_rescreen_panel_keeps_required_categories_and_prefers_interface_rank() -> None:
    ranking_rows = (
        _make_ranking_row(
            shortlist_rank=4,
            candidate_id="am1_mex_ss_only_u02_r2u02",
            campaign_id="am1_mex_ss_only",
            topology_class="am1_monomer",
            rosetta_score_rank_within_topology=1,
            representative_ligand_confidence=0.2096,
            total_score=111.2,
        ),
        _make_ranking_row(
            shortlist_rank=1,
            candidate_id="hans_pocket_ss_only_u02_r2u01",
            campaign_id="hans_pocket_ss_only",
            topology_class="hans_monomer",
            rosetta_score_rank_within_topology=1,
            representative_ligand_confidence=1.0,
            total_score=150.7,
        ),
        _make_ranking_row(
            shortlist_rank=2,
            candidate_id="hans_interface_ss_plus_if_u04_r2u01",
            campaign_id="hans_interface_ss_plus_if",
            topology_class="hans_interface_multichain",
            rosetta_score_rank_within_topology=2,
            representative_ligand_confidence=0.3454,
            total_score=540.0,
        ),
        _make_ranking_row(
            shortlist_rank=3,
            candidate_id="hans_interface_ss_plus_if_u04_r2u02",
            campaign_id="hans_interface_ss_plus_if",
            topology_class="hans_interface_multichain",
            rosetta_score_rank_within_topology=1,
            representative_ligand_confidence=0.3300,
            total_score=530.0,
        ),
    )

    panel_rows = select_round2_md_rescreen_panel(ranking_rows)

    assert [row.candidate_id for row in panel_rows] == [
        "am1_mex_ss_only_u02_r2u02",
        "hans_pocket_ss_only_u02_r2u01",
        "hans_interface_ss_plus_if_u04_r2u02",
    ]
    assert panel_rows[2].panel_role == "best_hans_interface_aware_round2_candidate"
    assert "rank first and ligand_confidence second" in panel_rows[2].selection_reason


def test_select_round2_md_rescreen_panel_uses_ligand_confidence_when_interface_ranks_tie() -> None:
    ranking_rows = (
        _make_ranking_row(
            shortlist_rank=2,
            candidate_id="hans_interface_ss_plus_if_u04_r2u01",
            campaign_id="hans_interface_ss_plus_if",
            topology_class="hans_interface_multichain",
            rosetta_score_rank_within_topology=1,
            representative_ligand_confidence=0.3454,
            total_score=530.0,
        ),
        _make_ranking_row(
            shortlist_rank=3,
            candidate_id="hans_interface_ss_plus_if_u04_r2u02",
            campaign_id="hans_interface_ss_plus_if",
            topology_class="hans_interface_multichain",
            rosetta_score_rank_within_topology=1,
            representative_ligand_confidence=0.3300,
            total_score=530.0,
        ),
    )

    panel_rows = select_round2_md_rescreen_panel(ranking_rows)

    assert [row.candidate_id for row in panel_rows] == ["hans_interface_ss_plus_if_u04_r2u01"]
