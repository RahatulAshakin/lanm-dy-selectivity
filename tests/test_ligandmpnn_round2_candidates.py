from pathlib import Path

from lanm.analysis.ligandmpnn_round2_candidates import (
    LigandMPNNRound2GeneratedDesign,
    LigandMPNNRound2SelectionCandidate,
    build_ligandmpnn_round2_shortlist,
    collapse_global_unique_sequences,
    deduplicate_candidate_sequences,
    parse_ligandmpnn_round2_candidate_outputs,
)
from lanm.models import LigandMPNNRound2RedesignPositionRow, LigandMPNNRound2UniqueSequenceRow


def _make_redesign_row(
    *,
    candidate_id: str,
    chain_id: str,
    sequence_index: int,
    residue_seq: int,
    seed_amino_acid: str,
    designed_amino_acid: str,
) -> LigandMPNNRound2RedesignPositionRow:
    return LigandMPNNRound2RedesignPositionRow(
        shortlist_rank=1,
        candidate_id=candidate_id,
        seed_rank=1,
        seed_candidate_id="seed_alpha",
        campaign_id="example_campaign",
        backbone_id="example_backbone",
        chain_id=chain_id,
        sequence_index=sequence_index,
        residue_seq=residue_seq,
        insertion_code="",
        ligandmpnn_residue_id=f"{chain_id}{residue_seq}",
        residue_name="ALA",
        seed_amino_acid=seed_amino_acid,
        designed_amino_acid=designed_amino_acid,
        mutation_token=f"{chain_id}{residue_seq}:{seed_amino_acid}>{designed_amino_acid}",
        canonical_family_position=18 + sequence_index,
        am1_mature_position=16 + sequence_index,
    )


def _make_candidate(
    tmp_path: Path,
    *,
    candidate_id: str,
    shortlist_rank: int,
    campaign_id: str = "example_campaign",
) -> LigandMPNNRound2SelectionCandidate:
    input_pdb_path = tmp_path / f"{candidate_id}.pdb"
    input_pdb_path.write_text("ATOM\nEND\n", encoding="utf-8")
    return LigandMPNNRound2SelectionCandidate(
        shortlist_rank=shortlist_rank,
        candidate_id=candidate_id,
        seed_rank=1,
        seed_candidate_id="seed_alpha",
        campaign_id=campaign_id,
        backbone_id="example_backbone",
        topology_class="example_topology",
        designed_chains=("A",),
        fixed_context_chains=(),
        output_chain_order=("A",),
        preserved_metal_identity="DY",
        redesign_rows=(
            _make_redesign_row(
                candidate_id=candidate_id,
                chain_id="A",
                sequence_index=2,
                residue_seq=20,
                seed_amino_acid="C",
                designed_amino_acid="Q",
            ),
            _make_redesign_row(
                candidate_id=candidate_id,
                chain_id="A",
                sequence_index=4,
                residue_seq=22,
                seed_amino_acid="E",
                designed_amino_acid="S",
            ),
        ),
        input_pdb_path=input_pdb_path,
        output_dir=tmp_path / "outputs" / candidate_id,
    )


def _make_generated_design(
    *,
    candidate_id: str,
    shortlist_rank: int,
    campaign_id: str,
    design_id: int,
    sequence_id: str,
    designed_chain_sequence: str,
    ligand_confidence: float,
    overall_confidence: float,
    seq_recovery: float,
    mutation_string: str,
) -> LigandMPNNRound2GeneratedDesign:
    mutation_count = 0 if not mutation_string else len(mutation_string.split(","))
    redesigned_residue_identifiers = "" if not mutation_string else ",".join(
        token.split(":", 1)[0] for token in mutation_string.split(",")
    )
    return LigandMPNNRound2GeneratedDesign(
        sequence_id=sequence_id,
        candidate_id=candidate_id,
        seed_rank=1,
        seed_candidate_id="seed_alpha",
        campaign_id=campaign_id,
        backbone_id="example_backbone",
        topology_class="example_topology",
        preserved_metal_identity="DY",
        candidate_shortlist_rank=shortlist_rank,
        design_id=design_id,
        overall_confidence=overall_confidence,
        ligand_confidence=ligand_confidence,
        seq_recovery=seq_recovery,
        mutation_count=mutation_count,
        mutation_string=mutation_string,
        redesigned_residue_identifiers=redesigned_residue_identifiers,
        designed_chain_sequence=designed_chain_sequence,
        input_pdb_path=f"inputs/{candidate_id}.pdb",
        output_fasta_path=f"seqs/{candidate_id}.fa",
        backbone_pdb_path=f"backbones/{sequence_id}.pdb",
        packed_pdb_path=f"packed/{sequence_id}.pdb",
    )


def _make_unique_row(
    *,
    sequence_id: str,
    candidate_id: str,
    shortlist_rank: int,
    campaign_id: str,
    mean_ligand_confidence: float,
    mean_overall_confidence: float,
    best_seq_recovery: float,
) -> LigandMPNNRound2UniqueSequenceRow:
    return LigandMPNNRound2UniqueSequenceRow(
        unique_sequence_rank=0,
        sequence_id=sequence_id,
        candidate_id=candidate_id,
        seed_rank=shortlist_rank,
        seed_candidate_id=f"seed_{candidate_id}",
        campaign_id=campaign_id,
        backbone_id=f"{campaign_id}_backbone",
        topology_class=f"{campaign_id}_topology",
        preserved_metal_identity="DY",
        candidate_shortlist_rank=shortlist_rank,
        candidate_unique_rank=1,
        design_id=1,
        total_occurrence_count=1,
        source_candidate_count=1,
        source_candidate_ids=candidate_id,
        mean_ligand_confidence=mean_ligand_confidence,
        mean_overall_confidence=mean_overall_confidence,
        best_seq_recovery=best_seq_recovery,
        mutation_count=1,
        mutation_string="A20:Q>E",
        redesigned_residue_identifiers="A20",
        designed_chain_sequence=sequence_id,
        input_pdb_path=f"inputs/{candidate_id}.pdb",
        output_fasta_path=f"seqs/{candidate_id}.fa",
        backbone_pdb_path=f"backbones/{sequence_id}.pdb",
        packed_pdb_path=f"packed/{sequence_id}.pdb",
    )


def test_parse_ligandmpnn_round2_candidate_outputs_extracts_multichain_designed_sequence(tmp_path: Path) -> None:
    candidate = LigandMPNNRound2SelectionCandidate(
        shortlist_rank=3,
        candidate_id="example_candidate",
        seed_rank=2,
        seed_candidate_id="seed_alpha",
        campaign_id="hans_interface_ss_plus_if",
        backbone_id="example_backbone",
        topology_class="example_interface",
        designed_chains=("A", "C"),
        fixed_context_chains=("B",),
        output_chain_order=("B", "A", "C"),
        preserved_metal_identity="DY",
        redesign_rows=(
            _make_redesign_row(
                candidate_id="example_candidate",
                chain_id="A",
                sequence_index=2,
                residue_seq=20,
                seed_amino_acid="Q",
                designed_amino_acid="Q",
            ),
            _make_redesign_row(
                candidate_id="example_candidate",
                chain_id="C",
                sequence_index=1,
                residue_seq=40,
                seed_amino_acid="C",
                designed_amino_acid="C",
            ),
        ),
        input_pdb_path=tmp_path / "example_candidate.pdb",
        output_dir=tmp_path / "outputs" / "example_candidate",
    )
    candidate.input_pdb_path.write_text("ATOM\nEND\n", encoding="utf-8")

    fasta_path = candidate.output_fasta_path
    fasta_path.parent.mkdir(parents=True)
    fasta_path.write_text(
        "\n".join(
            [
                ">example_candidate, T=0.1, seed=37, num_res=2, num_ligand_res=2, use_ligand_context=True, ligand_cutoff_distance=8.0, batch_size=2, number_of_batches=1, model_path=./model_params/ligandmpnn_v_32_010_25.pt",
                "BBBB:AQ:CC",
                ">example_candidate, id=1, T=0.1, seed=37, overall_confidence=0.8000, ligand_confidence=0.7000, seq_rec=0.5000",
                "BBBB:AR:DC",
                "",
            ]
        ),
        encoding="utf-8",
    )

    backbone_path = candidate.output_dir / "backbones" / f"{candidate.candidate_id}_1.pdb"
    packed_path = candidate.output_dir / "packed" / f"{candidate.candidate_id}_packed_1_1.pdb"
    backbone_path.parent.mkdir(parents=True, exist_ok=True)
    packed_path.parent.mkdir(parents=True, exist_ok=True)
    backbone_path.write_text("ATOM\nEND\n", encoding="utf-8")
    packed_path.write_text("ATOM\nEND\n", encoding="utf-8")

    designs = parse_ligandmpnn_round2_candidate_outputs(candidate)

    assert [design.sequence_id for design in designs] == ["example_candidate_design_01"]
    assert designs[0].designed_chain_sequence == "AR/DC"
    assert designs[0].mutation_count == 2
    assert designs[0].mutation_string == "A20:Q>R,C40:C>D"
    assert designs[0].redesigned_residue_identifiers == "A20,C40"


def test_collapse_global_unique_sequences_deduplicates_with_weighted_means(tmp_path: Path) -> None:
    candidate_one = _make_candidate(
        tmp_path,
        candidate_id="candidate_one",
        shortlist_rank=1,
        campaign_id="am1_mex_ss_only",
    )
    candidate_two = _make_candidate(
        tmp_path,
        candidate_id="candidate_two",
        shortlist_rank=2,
        campaign_id="hans_interface_ss_plus_if",
    )

    result_one = deduplicate_candidate_sequences(
        candidate_one,
        (
            _make_generated_design(
                candidate_id="candidate_one",
                shortlist_rank=1,
                campaign_id="am1_mex_ss_only",
                design_id=1,
                sequence_id="candidate_one_design_01",
                designed_chain_sequence="AAAA",
                ligand_confidence=0.90,
                overall_confidence=0.40,
                seq_recovery=0.20,
                mutation_string="A20:Q>E",
            ),
            _make_generated_design(
                candidate_id="candidate_one",
                shortlist_rank=1,
                campaign_id="am1_mex_ss_only",
                design_id=2,
                sequence_id="candidate_one_design_02",
                designed_chain_sequence="AAAA",
                ligand_confidence=0.60,
                overall_confidence=0.80,
                seq_recovery=0.50,
                mutation_string="A20:Q>E",
            ),
            _make_generated_design(
                candidate_id="candidate_one",
                shortlist_rank=1,
                campaign_id="am1_mex_ss_only",
                design_id=3,
                sequence_id="candidate_one_design_03",
                designed_chain_sequence="BBBB",
                ligand_confidence=0.30,
                overall_confidence=0.20,
                seq_recovery=0.10,
                mutation_string="A22:S>T",
            ),
        ),
    )
    result_two = deduplicate_candidate_sequences(
        candidate_two,
        (
            _make_generated_design(
                candidate_id="candidate_two",
                shortlist_rank=2,
                campaign_id="hans_interface_ss_plus_if",
                design_id=1,
                sequence_id="candidate_two_design_01",
                designed_chain_sequence="AAAA",
                ligand_confidence=0.30,
                overall_confidence=0.60,
                seq_recovery=0.90,
                mutation_string="A20:Q>E",
            ),
        ),
    )

    unique_rows = collapse_global_unique_sequences((result_one, result_two))

    assert len(unique_rows) == 2
    top_row = unique_rows[0]
    assert top_row.candidate_id == "candidate_one"
    assert top_row.total_occurrence_count == 3
    assert top_row.source_candidate_count == 2
    assert top_row.source_candidate_ids == "candidate_one,candidate_two"
    assert top_row.mean_ligand_confidence == 0.6
    assert round(top_row.mean_overall_confidence, 4) == 0.6
    assert top_row.best_seq_recovery == 0.9
    assert top_row.mutation_string == "A20:Q>E"


def test_build_ligandmpnn_round2_shortlist_ranks_selected_rows_globally() -> None:
    unique_rows = (
        _make_unique_row(
            sequence_id="am1_sequence",
            candidate_id="am1_candidate",
            shortlist_rank=1,
            campaign_id="am1_mex_ss_only",
            mean_ligand_confidence=0.20,
            mean_overall_confidence=0.40,
            best_seq_recovery=0.70,
        ),
        _make_unique_row(
            sequence_id="pocket_sequence",
            candidate_id="pocket_candidate",
            shortlist_rank=2,
            campaign_id="hans_pocket_ss_only",
            mean_ligand_confidence=0.95,
            mean_overall_confidence=0.50,
            best_seq_recovery=0.10,
        ),
        _make_unique_row(
            sequence_id="interface_sequence",
            candidate_id="interface_candidate",
            shortlist_rank=3,
            campaign_id="hans_interface_ss_plus_if",
            mean_ligand_confidence=0.55,
            mean_overall_confidence=0.60,
            best_seq_recovery=0.90,
        ),
        _make_unique_row(
            sequence_id="wildcard_sequence",
            candidate_id="wildcard_candidate",
            shortlist_rank=4,
            campaign_id="hans_interface_if_only",
            mean_ligand_confidence=0.80,
            mean_overall_confidence=0.65,
            best_seq_recovery=0.80,
        ),
        _make_unique_row(
            sequence_id="interface_backup",
            candidate_id="interface_backup_candidate",
            shortlist_rank=5,
            campaign_id="hans_interface_ss_plus_if",
            mean_ligand_confidence=0.50,
            mean_overall_confidence=0.70,
            best_seq_recovery=0.85,
        ),
    )

    shortlist = build_ligandmpnn_round2_shortlist(unique_rows)

    assert [row.sequence_id for row in shortlist.shortlist_rows] == [
        "pocket_sequence",
        "wildcard_sequence",
        "interface_sequence",
        "am1_sequence",
    ]
    assert [row.retention_reason for row in shortlist.shortlist_rows] == [
        "top Hans pocket round-2 candidate",
        "best remaining unique candidate",
        "top Hans interface-aware round-2 candidate",
        "top AM1/Mex round-2 candidate",
    ]
