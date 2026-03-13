from pathlib import Path

from lanm.analysis.proteinmpnn_candidates import (
    build_balanced_shortlist,
    deduplicate_campaign_sequences,
    extract_mutation_summary,
)
from lanm.analysis.proteinmpnn_smoke import (
    ProteinMPNNGeneratedRecord,
    ProteinMPNNNativeRecord,
    ProteinMPNNParsedOutput,
    ProteinMPNNSmokeCampaign,
)
from lanm.models import DesignCampaignPositionRow, ProteinMPNNUniqueSequenceRow


def _make_campaign(tmp_path: Path, campaign_id: str = "example_campaign") -> ProteinMPNNSmokeCampaign:
    return ProteinMPNNSmokeCampaign(
        campaign_id=campaign_id,
        backbone_id="example_backbone",
        backbone_structure_id="8FNS",
        design_set_name="campaign_ss_only",
        designed_chains=("A",),
        fixed_context_chains=(),
        pdb_path=tmp_path / f"{campaign_id}.pdb",
        chain_assignment_path=tmp_path / "chain_id.jsonl",
        fixed_positions_path=tmp_path / "fixed_positions.jsonl",
        output_dir=tmp_path / "outputs" / campaign_id,
    )


def _make_position_row(
    *,
    campaign_id: str = "example_campaign",
    sequence_index: int,
    residue_name: str,
    canonical_family_position: int | None,
    am1_mature_position: int | None,
    mutable_second_sphere: bool = False,
    mutable_interface: bool = False,
) -> DesignCampaignPositionRow:
    return DesignCampaignPositionRow(
        campaign_id=campaign_id,
        design_set_name="campaign_ss_only",
        backbone_id="example_backbone",
        structure_id="8FNS",
        chain_id="A",
        sequence_index=sequence_index,
        residue_seq=20 + sequence_index,
        insertion_code="",
        residue_name=residue_name,
        canonical_family_position=canonical_family_position,
        am1_mature_position=am1_mature_position,
        fixed_first_shell=False,
        protected_positions=False,
        mutable_second_sphere=mutable_second_sphere,
        mutable_interface=mutable_interface,
        hard_fixed=False,
        designable=True,
        position_state="designable",
        position_reason="campaign_designable",
        rationale="",
    )


def _make_unique_row(
    *,
    campaign_id: str,
    campaign_rank: int,
    designed_sequence: str,
    best_score: float,
    best_global_score: float,
    best_seq_recovery: float = 0.7,
) -> ProteinMPNNUniqueSequenceRow:
    return ProteinMPNNUniqueSequenceRow(
        candidate_id=f"{campaign_id}_u{campaign_rank:02d}",
        campaign_id=campaign_id,
        backbone_id=f"{campaign_id}_backbone",
        design_set_name="campaign_ss_only",
        designed_sequence=designed_sequence,
        representative_sequence_id=f"{campaign_id}_T0.1_sample_{campaign_rank:02d}",
        occurrence_count=1,
        campaign_rank=campaign_rank,
        temperature=0.1,
        best_score=best_score,
        best_global_score=best_global_score,
        best_seq_recovery=best_seq_recovery,
        mutation_count=campaign_rank,
        mutation_string=f"A{campaign_rank}V",
        canonical_family_positions_mutated=str(10 + campaign_rank),
        am1_mature_positions_mutated=str(20 + campaign_rank),
        includes_second_sphere_position=False,
        includes_interface_position=False,
    )


def test_deduplicate_campaign_sequences_collapses_duplicate_outputs(tmp_path: Path) -> None:
    campaign = _make_campaign(tmp_path)
    position_rows = (
        _make_position_row(sequence_index=1, residue_name="ALA", canonical_family_position=11, am1_mature_position=10),
        _make_position_row(
            sequence_index=2,
            residue_name="CYS",
            canonical_family_position=18,
            am1_mature_position=16,
            mutable_second_sphere=True,
        ),
        _make_position_row(sequence_index=3, residue_name="ASP", canonical_family_position=22, am1_mature_position=20),
        _make_position_row(
            sequence_index=4,
            residue_name="GLU",
            canonical_family_position=27,
            am1_mature_position=25,
            mutable_interface=True,
        ),
    )
    parsed_output = ProteinMPNNParsedOutput(
        native_record=ProteinMPNNNativeRecord(
            name="example_campaign",
            score=1.0,
            global_score=1.1,
            fixed_chains=(),
            designed_chains=("A",),
            model_name="v_48_020",
            git_hash="abc123",
            seed=37,
            sequence="ACDE",
        ),
        generated_records=(
            ProteinMPNNGeneratedRecord(
                temperature=0.15,
                sample_number=1,
                score=0.60,
                global_score=1.00,
                seq_recovery=0.75,
                sequence="AGDE",
            ),
            ProteinMPNNGeneratedRecord(
                temperature=0.10,
                sample_number=2,
                score=0.50,
                global_score=1.10,
                seq_recovery=0.50,
                sequence="AGDE",
            ),
            ProteinMPNNGeneratedRecord(
                temperature=0.10,
                sample_number=3,
                score=0.55,
                global_score=0.90,
                seq_recovery=0.75,
                sequence="ACDF",
            ),
        ),
    )

    result = deduplicate_campaign_sequences(
        campaign=campaign,
        parsed_output=parsed_output,
        position_rows=position_rows,
    )

    assert result.raw_sequence_count == 3
    assert len(result.unique_rows) == 2

    top_row = result.unique_rows[0]
    assert top_row.candidate_id == "example_campaign_u01"
    assert top_row.designed_sequence == "AGDE"
    assert top_row.representative_sequence_id == "example_campaign_T0.1_sample_02"
    assert top_row.occurrence_count == 2
    assert top_row.temperature == 0.1
    assert top_row.best_score == 0.5
    assert top_row.best_global_score == 1.0
    assert top_row.best_seq_recovery == 0.75
    assert top_row.mutation_count == 1
    assert top_row.mutation_string == "C2G"
    assert top_row.canonical_family_positions_mutated == "18"

    second_row = result.unique_rows[1]
    assert second_row.candidate_id == "example_campaign_u02"
    assert second_row.designed_sequence == "ACDF"
    assert second_row.mutation_string == "E4F"
    assert second_row.includes_interface_position is True


def test_extract_mutation_summary_formats_mutations_and_position_sets() -> None:
    position_rows = (
        _make_position_row(sequence_index=1, residue_name="ALA", canonical_family_position=13, am1_mature_position=12),
        _make_position_row(
            sequence_index=2,
            residue_name="CYS",
            canonical_family_position=18,
            am1_mature_position=16,
            mutable_second_sphere=True,
            mutable_interface=True,
        ),
        _make_position_row(
            sequence_index=3,
            residue_name="ASP",
            canonical_family_position=22,
            am1_mature_position=20,
            mutable_second_sphere=True,
        ),
        _make_position_row(sequence_index=4, residue_name="GLU", canonical_family_position=24, am1_mature_position=22),
    )

    summary = extract_mutation_summary(
        native_sequence="ACDE",
        designed_sequence="AFGE",
        designed_chains=("A",),
        position_rows=position_rows,
    )

    assert summary.mutation_count == 2
    assert summary.mutation_string == "C2F,D3G"
    assert summary.canonical_family_positions_mutated == "18,22"
    assert summary.am1_mature_positions_mutated == "16,20"
    assert summary.includes_second_sphere_position is True
    assert summary.includes_interface_position is True


def test_build_balanced_shortlist_honors_campaign_minima_and_duplicate_exclusion() -> None:
    unique_rows = (
        _make_unique_row(
            campaign_id="hans_interface_ss_plus_if",
            campaign_rank=1,
            designed_sequence="SEQ_P1",
            best_score=0.10,
            best_global_score=0.90,
        ),
        _make_unique_row(
            campaign_id="hans_interface_ss_plus_if",
            campaign_rank=2,
            designed_sequence="SEQ_DUP",
            best_score=0.20,
            best_global_score=1.00,
        ),
        _make_unique_row(
            campaign_id="hans_interface_ss_plus_if",
            campaign_rank=3,
            designed_sequence="SEQ_P3",
            best_score=0.30,
            best_global_score=1.10,
        ),
        _make_unique_row(
            campaign_id="am1_mex_ss_only",
            campaign_rank=1,
            designed_sequence="SEQ_DUP",
            best_score=0.05,
            best_global_score=0.80,
        ),
        _make_unique_row(
            campaign_id="am1_mex_ss_only",
            campaign_rank=2,
            designed_sequence="SEQ_A2",
            best_score=0.25,
            best_global_score=0.95,
        ),
        _make_unique_row(
            campaign_id="hans_interface_if_only",
            campaign_rank=1,
            designed_sequence="SEQ_I1",
            best_score=0.15,
            best_global_score=0.92,
        ),
        _make_unique_row(
            campaign_id="hans_interface_if_only",
            campaign_rank=2,
            designed_sequence="SEQ_I2",
            best_score=0.45,
            best_global_score=1.20,
        ),
        _make_unique_row(
            campaign_id="hans_pocket_ss_only",
            campaign_rank=1,
            designed_sequence="SEQ_H1",
            best_score=0.12,
            best_global_score=0.91,
        ),
        _make_unique_row(
            campaign_id="hans_pocket_ss_only",
            campaign_rank=2,
            designed_sequence="SEQ_H2",
            best_score=0.32,
            best_global_score=1.05,
        ),
    )

    result = build_balanced_shortlist(
        unique_rows,
        per_campaign_limit=2,
        total_limit=7,
    )

    assert [row.candidate_id for row in result.shortlist_rows] == [
        "hans_interface_ss_plus_if_u01",
        "hans_interface_ss_plus_if_u02",
        "am1_mex_ss_only_u02",
        "hans_interface_if_only_u01",
        "hans_pocket_ss_only_u01",
        "hans_interface_if_only_u02",
        "hans_pocket_ss_only_u02",
    ]
    assert all(row.designed_sequence != "SEQ_DUP" or row.campaign_id == "hans_interface_ss_plus_if" for row in result.shortlist_rows)
    assert result.skipped_cross_campaign_duplicates == 1
    assert result.shortlist_rows[0].retention_reason == "required interface campaign coverage"
    assert result.shortlist_rows[2].retention_reason == "required campaign coverage"
    assert result.shortlist_rows[-1].retention_reason == "balanced round-robin fill (round 1)"
    assert all(row.candidate_id != "hans_interface_ss_plus_if_u03" for row in result.shortlist_rows)
