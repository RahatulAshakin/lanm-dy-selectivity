from pathlib import Path

import pytest

from lanm.analysis.rosetta_score_smoke import (
    RosettaScoreSmokeCandidate,
    RosettaScoreSmokeResult,
    build_rosetta_score_command,
    choose_representative_packed_pdb,
    parse_rosetta_scorefile,
    rank_rosetta_score_candidates,
)
from lanm.models import LigandMPNNSequenceCatalogRow


def _make_sequence_row(
    *,
    design_id: int,
    ligand_confidence: float,
    overall_confidence: float,
    packed_pdb_path: str,
) -> LigandMPNNSequenceCatalogRow:
    return LigandMPNNSequenceCatalogRow(
        shortlist_rank=1,
        candidate_id="example_candidate",
        campaign_id="example_campaign",
        backbone_id="example_backbone",
        source_structure_id="8FNS",
        designed_chains="A",
        preserved_metal_identity="ND",
        redesigned_residues="A35,A44",
        design_id=design_id,
        temperature=0.1,
        seed=37,
        overall_confidence=overall_confidence,
        ligand_confidence=ligand_confidence,
        seq_recovery=0.5,
        sequence_id=f"example_candidate_design_{design_id:02d}",
        designed_sequence="AAAA",
        sequence_length=4,
        output_fasta_path="results/ligandmpnn_smoke/example_candidate/seqs/example_candidate.fa",
        backbone_pdb_path=f"results/ligandmpnn_smoke/example_candidate/backbones/example_candidate_{design_id}.pdb",
        packed_pdb_path=packed_pdb_path,
    )


def _make_result(*, candidate_id: str, shortlist_rank: int, total_score: float) -> RosettaScoreSmokeResult:
    candidate = RosettaScoreSmokeCandidate(
        shortlist_rank=shortlist_rank,
        candidate_id=candidate_id,
        campaign_id="example_campaign",
        backbone_id="example_backbone",
        preserved_metal_identity="ND",
        representative_design_id=1,
        representative_ligand_confidence=0.8,
        representative_overall_confidence=0.9,
        representative_packed_pdb=Path(f"/tmp/{candidate_id}.pdb"),
        output_dir=Path(f"/tmp/{candidate_id}"),
    )
    return RosettaScoreSmokeResult(
        candidate=candidate,
        scores={"total_score": total_score, "fa_atr": -10.0, "description": f"{candidate_id}_0001"},
        scorefile_path=Path(f"/tmp/{candidate_id}.sc"),
    )


def test_choose_representative_packed_pdb_prefers_ligand_then_overall_then_smallest_design_id() -> None:
    rows = (
        _make_sequence_row(
            design_id=3,
            ligand_confidence=0.75,
            overall_confidence=0.95,
            packed_pdb_path="results/ligandmpnn_smoke/example_candidate/packed/example_candidate_packed_3_1.pdb",
        ),
        _make_sequence_row(
            design_id=2,
            ligand_confidence=0.8,
            overall_confidence=0.7,
            packed_pdb_path="results/ligandmpnn_smoke/example_candidate/packed/example_candidate_packed_2_1.pdb",
        ),
        _make_sequence_row(
            design_id=1,
            ligand_confidence=0.8,
            overall_confidence=0.7,
            packed_pdb_path="results/ligandmpnn_smoke/example_candidate/packed/example_candidate_packed_1_1.pdb",
        ),
    )

    representative = choose_representative_packed_pdb(rows)

    assert representative.design_id == 1
    assert representative.packed_pdb_path.endswith("example_candidate_packed_1_1.pdb")


def test_build_rosetta_score_command() -> None:
    command = build_rosetta_score_command(
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


def test_parse_rosetta_scorefile_reads_single_score_row(tmp_path: Path) -> None:
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

    parsed = parse_rosetta_scorefile(scorefile_path)

    assert parsed == {
        "total_score": pytest.approx(-123.45),
        "fa_atr": pytest.approx(-456.7),
        "fa_rep": pytest.approx(89.1),
        "fa_sol": pytest.approx(12.34),
        "description": "example_input_0001",
    }


def test_rank_rosetta_score_candidates_uses_total_score_then_candidate_id() -> None:
    higher_score = _make_result(candidate_id="gamma_candidate", shortlist_rank=3, total_score=-19.0)
    tied_beta = _make_result(candidate_id="beta_candidate", shortlist_rank=1, total_score=-20.0)
    tied_alpha = _make_result(candidate_id="alpha_candidate", shortlist_rank=2, total_score=-20.0)

    ranked = rank_rosetta_score_candidates((tied_beta, higher_score, tied_alpha))

    assert [result.candidate.candidate_id for result in ranked] == [
        "alpha_candidate",
        "beta_candidate",
        "gamma_candidate",
    ]
