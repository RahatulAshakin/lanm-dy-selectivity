from pathlib import Path

from lanm.analysis.ligandmpnn_candidates import (
    CandidateUniqueSequenceResult,
    LigandMPNNGeneratedDesign,
    LigandMPNNSelectionCandidate,
    build_ligandmpnn_shortlist,
    collapse_global_unique_sequences,
    deduplicate_candidate_sequences,
    extract_designed_chain_sequence,
    parse_ligandmpnn_candidate_outputs,
)
from lanm.models import LigandMPNNRedesignPositionRow, LigandMPNNUniqueSequenceRow


def _make_redesign_row(
    *,
    candidate_id: str,
    sequence_index: int,
    residue_seq: int,
    native_amino_acid: str,
    designed_amino_acid: str,
) -> LigandMPNNRedesignPositionRow:
    return LigandMPNNRedesignPositionRow(
        shortlist_rank=1,
        candidate_id=candidate_id,
        campaign_id="example_campaign",
        backbone_id="example_backbone",
        chain_id="A",
        sequence_index=sequence_index,
        residue_seq=residue_seq,
        insertion_code="",
        ligandmpnn_residue_id=f"A{residue_seq}",
        residue_name="ALA",
        native_amino_acid=native_amino_acid,
        designed_amino_acid=designed_amino_acid,
        mutation_token=f"{native_amino_acid}{sequence_index}{designed_amino_acid}",
        canonical_family_position=None,
        am1_mature_position=None,
    )


def _make_candidate(tmp_path: Path, *, candidate_id: str) -> LigandMPNNSelectionCandidate:
    input_pdb_path = tmp_path / f"{candidate_id}.pdb"
    input_pdb_path.write_text("ATOM\nEND\n", encoding="utf-8")
    return LigandMPNNSelectionCandidate(
        shortlist_rank=1,
        candidate_id=candidate_id,
        campaign_id="example_campaign",
        backbone_id="example_backbone",
        designed_chains=("A",),
        preserved_chains=("A",),
        preserved_metal_identity="DY",
        input_sequence="AQRS",
        redesign_rows=(
            _make_redesign_row(
                candidate_id=candidate_id,
                sequence_index=2,
                residue_seq=20,
                native_amino_acid="C",
                designed_amino_acid="Q",
            ),
            _make_redesign_row(
                candidate_id=candidate_id,
                sequence_index=4,
                residue_seq=22,
                native_amino_acid="E",
                designed_amino_acid="S",
            ),
        ),
        input_pdb_path=input_pdb_path,
        output_dir=tmp_path / "outputs" / candidate_id,
    )


def _make_generated_design(
    *,
    candidate_id: str,
    campaign_id: str,
    candidate_shortlist_rank: int,
    design_id: int,
    sequence_id: str,
    designed_chain_sequence: str,
    ligand_confidence: float,
    overall_confidence: float,
    seq_recovery: float,
    mutation_count: int,
) -> LigandMPNNGeneratedDesign:
    return LigandMPNNGeneratedDesign(
        sequence_id=sequence_id,
        candidate_id=candidate_id,
        campaign_id=campaign_id,
        backbone_id="example_backbone",
        preserved_metal_identity="DY",
        candidate_shortlist_rank=candidate_shortlist_rank,
        design_id=design_id,
        overall_confidence=overall_confidence,
        ligand_confidence=ligand_confidence,
        seq_recovery=seq_recovery,
        mutation_count=mutation_count,
        mutation_string="" if mutation_count == 0 else f"A1V_{sequence_id}",
        redesigned_residue_identifiers="" if mutation_count == 0 else "A20",
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
    campaign_id: str,
    candidate_shortlist_rank: int,
    candidate_unique_rank: int,
    ligand_confidence: float,
    overall_confidence: float,
    seq_recovery: float,
    candidate_prefilter_selected: bool = True,
) -> LigandMPNNUniqueSequenceRow:
    return LigandMPNNUniqueSequenceRow(
        unique_sequence_rank=0,
        sequence_id=sequence_id,
        candidate_id=candidate_id,
        campaign_id=campaign_id,
        backbone_id=f"{campaign_id}_backbone",
        preserved_metal_identity="DY",
        candidate_shortlist_rank=candidate_shortlist_rank,
        candidate_unique_rank=candidate_unique_rank,
        candidate_prefilter_selected=candidate_prefilter_selected,
        design_id=candidate_unique_rank,
        total_occurrence_count=1,
        source_candidate_count=1,
        source_candidate_ids=candidate_id,
        overall_confidence=overall_confidence,
        ligand_confidence=ligand_confidence,
        seq_recovery=seq_recovery,
        mutation_count=candidate_unique_rank,
        mutation_string=f"A{candidate_unique_rank}V",
        redesigned_residue_identifiers="A20",
        designed_chain_sequence=sequence_id,
        input_pdb_path=f"inputs/{candidate_id}.pdb",
        output_fasta_path=f"seqs/{candidate_id}.fa",
        backbone_pdb_path=f"backbones/{sequence_id}.pdb",
        packed_pdb_path=f"packed/{sequence_id}.pdb",
    )


def test_parse_ligandmpnn_candidate_outputs_extracts_designed_chain_sequences(tmp_path: Path) -> None:
    candidate = LigandMPNNSelectionCandidate(
        shortlist_rank=1,
        candidate_id="example_candidate",
        campaign_id="example_campaign",
        backbone_id="example_backbone",
        designed_chains=("A",),
        preserved_chains=("A", "B", "C"),
        preserved_metal_identity="DY",
        input_sequence="AQRS",
        redesign_rows=(
            _make_redesign_row(
                candidate_id="example_candidate",
                sequence_index=2,
                residue_seq=20,
                native_amino_acid="C",
                designed_amino_acid="Q",
            ),
            _make_redesign_row(
                candidate_id="example_candidate",
                sequence_index=4,
                residue_seq=22,
                native_amino_acid="E",
                designed_amino_acid="S",
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
                "ACRE:BBBB:CCCC",
                ">example_candidate, id=1, T=0.1, seed=37, overall_confidence=0.8000, ligand_confidence=0.7000, seq_rec=0.5000",
                "AERS:BBBB:CCCC",
                ">example_candidate, id=2, T=0.1, seed=37, overall_confidence=0.7500, ligand_confidence=0.6900, seq_rec=0.4000",
                "AQRT:BBBB:CCCC",
                "",
            ]
        ),
        encoding="utf-8",
    )

    for design_id in (1, 2):
        backbone_path = candidate.output_dir / "backbones" / f"{candidate.candidate_id}_{design_id}.pdb"
        packed_path = candidate.output_dir / "packed" / f"{candidate.candidate_id}_packed_{design_id}_1.pdb"
        backbone_path.parent.mkdir(parents=True, exist_ok=True)
        packed_path.parent.mkdir(parents=True, exist_ok=True)
        backbone_path.write_text("ATOM\nEND\n", encoding="utf-8")
        packed_path.write_text("ATOM\nEND\n", encoding="utf-8")

    designs = parse_ligandmpnn_candidate_outputs(candidate)

    assert [design.sequence_id for design in designs] == [
        "example_candidate_design_01",
        "example_candidate_design_02",
    ]
    assert designs[0].designed_chain_sequence == "AERS"
    assert designs[0].mutation_count == 1
    assert designs[0].mutation_string == "Q2E"
    assert designs[0].redesigned_residue_identifiers == "A20"
    assert designs[1].mutation_string == "S4T"
    assert designs[1].packed_pdb_path.endswith("example_candidate_packed_2_1.pdb")


def test_extract_designed_chain_sequence_handles_multichain_outputs() -> None:
    designed_sequence = extract_designed_chain_sequence(
        "AAAA:BBBB:CCCC",
        preserved_chains=("A", "B", "C"),
        designed_chains=("A", "C"),
    )

    assert designed_sequence == "AAAA/CCCC"


def test_collapse_global_unique_sequences_deduplicates_within_and_across_candidates(tmp_path: Path) -> None:
    candidate_one = _make_candidate(tmp_path, candidate_id="candidate_one")
    candidate_two = _make_candidate(tmp_path, candidate_id="candidate_two")

    result_one = deduplicate_candidate_sequences(
        candidate_one,
        (
            _make_generated_design(
                candidate_id="candidate_one",
                campaign_id="campaign_alpha",
                candidate_shortlist_rank=1,
                design_id=1,
                sequence_id="candidate_one_design_01",
                designed_chain_sequence="AAAA",
                ligand_confidence=0.80,
                overall_confidence=0.70,
                seq_recovery=0.40,
                mutation_count=1,
            ),
            _make_generated_design(
                candidate_id="candidate_one",
                campaign_id="campaign_alpha",
                candidate_shortlist_rank=1,
                design_id=2,
                sequence_id="candidate_one_design_02",
                designed_chain_sequence="AAAA",
                ligand_confidence=0.80,
                overall_confidence=0.90,
                seq_recovery=0.45,
                mutation_count=1,
            ),
            _make_generated_design(
                candidate_id="candidate_one",
                campaign_id="campaign_alpha",
                candidate_shortlist_rank=1,
                design_id=3,
                sequence_id="candidate_one_design_03",
                designed_chain_sequence="AAAT",
                ligand_confidence=0.70,
                overall_confidence=0.80,
                seq_recovery=0.35,
                mutation_count=2,
            ),
        ),
    )
    result_two = deduplicate_candidate_sequences(
        candidate_two,
        (
            _make_generated_design(
                candidate_id="candidate_two",
                campaign_id="campaign_beta",
                candidate_shortlist_rank=2,
                design_id=1,
                sequence_id="candidate_two_design_01",
                designed_chain_sequence="GGGG",
                ligand_confidence=1.00,
                overall_confidence=0.95,
                seq_recovery=0.55,
                mutation_count=1,
            ),
            _make_generated_design(
                candidate_id="candidate_two",
                campaign_id="campaign_beta",
                candidate_shortlist_rank=2,
                design_id=2,
                sequence_id="candidate_two_design_02",
                designed_chain_sequence="GGGT",
                ligand_confidence=0.95,
                overall_confidence=0.85,
                seq_recovery=0.50,
                mutation_count=1,
            ),
            _make_generated_design(
                candidate_id="candidate_two",
                campaign_id="campaign_beta",
                candidate_shortlist_rank=2,
                design_id=3,
                sequence_id="candidate_two_design_03",
                designed_chain_sequence="AAAA",
                ligand_confidence=0.90,
                overall_confidence=0.88,
                seq_recovery=0.45,
                mutation_count=2,
            ),
        ),
    )

    unique_rows = collapse_global_unique_sequences((result_one, result_two))

    assert result_one.raw_sequence_count == 3
    assert len(result_one.unique_rows) == 2
    assert [row.candidate_unique_rank for row in result_one.unique_rows] == [1, 2]
    assert result_one.unique_rows[0].sequence_id == "candidate_one_design_02"

    assert len(unique_rows) == 4
    duplicate_row = next(row for row in unique_rows if row.designed_chain_sequence == "AAAA")
    assert duplicate_row.candidate_id == "candidate_one"
    assert duplicate_row.sequence_id == "candidate_one_design_02"
    assert duplicate_row.total_occurrence_count == 3
    assert duplicate_row.source_candidate_count == 2
    assert duplicate_row.source_candidate_ids == "candidate_one,candidate_two"
    assert duplicate_row.candidate_prefilter_selected is True


def test_build_ligandmpnn_shortlist_honors_required_campaign_coverage() -> None:
    unique_rows = (
        _make_unique_row(
            sequence_id="IFACE_01",
            candidate_id="iface_candidate_01",
            campaign_id="hans_interface_ss_plus_if",
            candidate_shortlist_rank=1,
            candidate_unique_rank=1,
            ligand_confidence=0.95,
            overall_confidence=0.90,
            seq_recovery=0.60,
        ),
        _make_unique_row(
            sequence_id="IFACE_02",
            candidate_id="iface_candidate_02",
            campaign_id="hans_interface_ss_plus_if",
            candidate_shortlist_rank=2,
            candidate_unique_rank=1,
            ligand_confidence=0.93,
            overall_confidence=0.89,
            seq_recovery=0.58,
        ),
        _make_unique_row(
            sequence_id="IFACE_03",
            candidate_id="iface_candidate_03",
            campaign_id="hans_interface_ss_plus_if",
            candidate_shortlist_rank=3,
            candidate_unique_rank=2,
            ligand_confidence=0.88,
            overall_confidence=0.84,
            seq_recovery=0.52,
        ),
        _make_unique_row(
            sequence_id="POCKET_01",
            candidate_id="pocket_candidate_01",
            campaign_id="hans_pocket_ss_only",
            candidate_shortlist_rank=4,
            candidate_unique_rank=1,
            ligand_confidence=0.91,
            overall_confidence=0.87,
            seq_recovery=0.57,
        ),
        _make_unique_row(
            sequence_id="AM1_01",
            candidate_id="am1_candidate_01",
            campaign_id="am1_mex_ss_only",
            candidate_shortlist_rank=5,
            candidate_unique_rank=1,
            ligand_confidence=0.89,
            overall_confidence=0.86,
            seq_recovery=0.54,
        ),
        _make_unique_row(
            sequence_id="AM1_02",
            candidate_id="am1_candidate_02",
            campaign_id="am1_mex_ss_only",
            candidate_shortlist_rank=6,
            candidate_unique_rank=2,
            ligand_confidence=0.84,
            overall_confidence=0.82,
            seq_recovery=0.50,
        ),
        _make_unique_row(
            sequence_id="IFONLY_01",
            candidate_id="ifonly_candidate_01",
            campaign_id="hans_interface_if_only",
            candidate_shortlist_rank=7,
            candidate_unique_rank=1,
            ligand_confidence=0.85,
            overall_confidence=0.81,
            seq_recovery=0.49,
        ),
        _make_unique_row(
            sequence_id="IFONLY_02",
            candidate_id="ifonly_candidate_02",
            campaign_id="hans_interface_if_only",
            candidate_shortlist_rank=8,
            candidate_unique_rank=3,
            ligand_confidence=0.99,
            overall_confidence=0.98,
            seq_recovery=0.70,
            candidate_prefilter_selected=False,
        ),
    )

    shortlist = build_ligandmpnn_shortlist(unique_rows, total_limit=7)

    assert [row.sequence_id for row in shortlist.shortlist_rows] == [
        "IFACE_01",
        "IFACE_02",
        "POCKET_01",
        "AM1_01",
        "IFACE_03",
        "AM1_02",
        "IFONLY_01",
    ]
    assert shortlist.shortlist_rows[0].retention_reason == "required interface campaign coverage"
    assert shortlist.shortlist_rows[2].retention_reason == "required pocket campaign coverage"
    assert shortlist.shortlist_rows[3].retention_reason == "required AM1/Mex coverage"
    assert shortlist.shortlist_rows[-1].retention_reason == "balanced round-robin fill (round 1)"
    assert all(row.sequence_id != "IFONLY_02" for row in shortlist.shortlist_rows)
