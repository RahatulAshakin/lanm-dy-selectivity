from pathlib import Path

import pytest

from lanm.analysis.rosetta_relax_smoke import (
    RosettaRelaxSmokeCandidate,
    RosettaRelaxSmokeResult,
    build_rosetta_relax_command,
    parse_rosetta_relax_scorefile,
    rank_rosetta_relax_candidates,
    select_top_rosetta_relax_candidates,
)


def test_select_top_rosetta_relax_candidates_reads_top_three_from_ranking_csv(tmp_path: Path) -> None:
    ranking_path = tmp_path / "rosetta_score_candidate_ranking.csv"
    output_root = tmp_path / "relax_outputs"

    rows = [
        _ranking_row(tmp_path=tmp_path, rosetta_rank=4, candidate_id="delta_candidate", total_score=40.0),
        _ranking_row(tmp_path=tmp_path, rosetta_rank=2, candidate_id="beta_candidate", total_score=20.0),
        _ranking_row(tmp_path=tmp_path, rosetta_rank=1, candidate_id="alpha_candidate", total_score=10.0),
        _ranking_row(tmp_path=tmp_path, rosetta_rank=3, candidate_id="gamma_candidate", total_score=30.0),
    ]
    ranking_path.write_text(
        "\n".join(
            [
                ",".join(rows[0].keys()),
                *(",".join(str(row[key]) for key in row.keys()) for row in rows),
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    candidates = select_top_rosetta_relax_candidates(
        ranking_path=ranking_path,
        output_root=output_root,
    )

    assert [candidate.candidate_id for candidate in candidates] == [
        "alpha_candidate",
        "beta_candidate",
        "gamma_candidate",
    ]
    assert [candidate.rosetta_rank for candidate in candidates] == [1, 2, 3]
    assert candidates[0].pre_relax_scores["total_score"] == pytest.approx(10.0)
    assert candidates[0].output_dir == output_root / "alpha_candidate"


def test_build_rosetta_relax_command() -> None:
    command = build_rosetta_relax_command(
        rosetta_executable=Path("/opt/rosetta/source/bin/relax.default.linuxgccrelease"),
        rosetta_database=Path("/opt/rosetta/database"),
        input_pdb=Path("/tmp/example_input.pdb"),
        output_dir=Path("/tmp/relax"),
        scorefile_path=Path("/tmp/relax/relax.sc"),
    )

    assert command == (
        "/opt/rosetta/source/bin/relax.default.linuxgccrelease",
        "-database",
        "/opt/rosetta/database",
        "-in:file:fullatom",
        "-s",
        "/tmp/example_input.pdb",
        "-out:path:all",
        "/tmp/relax",
        "-out:file:scorefile",
        "/tmp/relax/relax.sc",
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
        "-relax:constrain_relax_to_start_coords",
        "-relax:ramp_constraints",
        "false",
        "-score:set_weights",
        "fa_elec",
        "0.0",
    )
    assert "-in:file:keep_input_protonation_state" not in command


def test_parse_rosetta_relax_scorefile_reads_single_score_row(tmp_path: Path) -> None:
    scorefile_path = tmp_path / "relax.sc"
    scorefile_path.write_text(
        "\n".join(
            [
                "SEQUENCE:",
                "SCORE: total_score fa_atr fa_rep coordinate_constraint description",
                "SCORE: -98.765 -123.400 45.600 8.900 relaxed_input_0001",
                "",
            ]
        ),
        encoding="utf-8",
    )

    parsed = parse_rosetta_relax_scorefile(scorefile_path)

    assert parsed == {
        "total_score": pytest.approx(-98.765),
        "fa_atr": pytest.approx(-123.4),
        "fa_rep": pytest.approx(45.6),
        "coordinate_constraint": pytest.approx(8.9),
        "description": "relaxed_input_0001",
    }


def test_rank_rosetta_relax_candidates_skips_failures_and_breaks_ties_with_delta_total_score() -> None:
    alpha_candidate = _make_candidate(candidate_id="alpha_candidate", rosetta_rank=1, pre_relax_total_score=14.0)
    beta_candidate = _make_candidate(candidate_id="beta_candidate", rosetta_rank=2, pre_relax_total_score=12.0)
    failed_candidate = _make_candidate(candidate_id="failed_candidate", rosetta_rank=3, pre_relax_total_score=11.0)

    alpha_result = RosettaRelaxSmokeResult(
        candidate=alpha_candidate,
        pre_relax_scorefile_copy=Path("/tmp/alpha/pre_relax_score.sc"),
        relax_scorefile_path=Path("/tmp/alpha/relax/relax.sc"),
        relaxed_pdb_path=Path("/tmp/alpha/relax/alpha_candidate_0001.pdb"),
        post_relax_scores={"total_score": 5.0, "fa_atr": -10.0, "description": "alpha_candidate_0001"},
        success_status="success",
        failure_status="",
        failure_reason="",
    )
    beta_result = RosettaRelaxSmokeResult(
        candidate=beta_candidate,
        pre_relax_scorefile_copy=Path("/tmp/beta/pre_relax_score.sc"),
        relax_scorefile_path=Path("/tmp/beta/relax/relax.sc"),
        relaxed_pdb_path=Path("/tmp/beta/relax/beta_candidate_0001.pdb"),
        post_relax_scores={"total_score": 5.0, "fa_atr": -9.0, "description": "beta_candidate_0001"},
        success_status="success",
        failure_status="",
        failure_reason="",
    )
    failed_result = RosettaRelaxSmokeResult(
        candidate=failed_candidate,
        pre_relax_scorefile_copy=Path("/tmp/failed/pre_relax_score.sc"),
        relax_scorefile_path=Path("/tmp/failed/relax/relax.sc"),
        relaxed_pdb_path=None,
        post_relax_scores=None,
        success_status="failed",
        failure_status="relax_failed",
        failure_reason="Command failed",
    )

    ranked = rank_rosetta_relax_candidates((beta_result, failed_result, alpha_result))

    assert [result.candidate.candidate_id for result in ranked] == [
        "alpha_candidate",
        "beta_candidate",
    ]
    assert ranked[0].delta_total_score == pytest.approx(-9.0)
    assert ranked[1].delta_total_score == pytest.approx(-7.0)


def _ranking_row(*, tmp_path: Path, rosetta_rank: int, candidate_id: str, total_score: float) -> dict[str, object]:
    packed_pdb_path = tmp_path / f"{candidate_id}.pdb"
    packed_pdb_path.write_text("ATOM\n", encoding="utf-8")
    scorefile_path = tmp_path / f"{candidate_id}.sc"
    scorefile_path.write_text("SCORE:\n", encoding="utf-8")
    return {
        "rosetta_rank": rosetta_rank,
        "shortlist_rank": rosetta_rank,
        "candidate_id": candidate_id,
        "campaign_id": "example_campaign",
        "backbone_id": "example_backbone",
        "preserved_metal_identity": "ND",
        "representative_design_id": 1,
        "representative_ligand_confidence": 0.8,
        "representative_overall_confidence": 0.9,
        "representative_packed_pdb": packed_pdb_path,
        "scorefile_path": scorefile_path,
        "total_score": total_score,
        "score_description": f"{candidate_id}_0001",
        "fa_atr": -10.0 * rosetta_rank,
        "fa_rep": 5.0 * rosetta_rank,
    }


def _make_candidate(*, candidate_id: str, rosetta_rank: int, pre_relax_total_score: float) -> RosettaRelaxSmokeCandidate:
    return RosettaRelaxSmokeCandidate(
        rosetta_rank=rosetta_rank,
        shortlist_rank=rosetta_rank,
        candidate_id=candidate_id,
        campaign_id="example_campaign",
        backbone_id="example_backbone",
        preserved_metal_identity="ND",
        representative_design_id=1,
        representative_ligand_confidence=0.8,
        representative_overall_confidence=0.9,
        representative_packed_pdb=Path(f"/tmp/{candidate_id}.pdb"),
        pre_relax_scorefile_source=Path(f"/tmp/{candidate_id}.sc"),
        pre_relax_scores={
            "total_score": pre_relax_total_score,
            "fa_atr": -10.0,
            "description": f"{candidate_id}_0001",
        },
        output_dir=Path(f"/tmp/{candidate_id}"),
    )
