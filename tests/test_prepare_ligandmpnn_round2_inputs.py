import csv
from pathlib import Path

from lanm.analysis.ligandmpnn_round2_inputs import (
    build_round2_shortlist,
    deduplicate_round2_campaign_sequences,
    export_ligandmpnn_round2_inputs,
    extract_round2_mutation_summary,
)
from lanm.analysis.proteinmpnn_round2_smoke import (
    ProteinMPNNRound2Campaign,
    Round2DesignedResidue,
    Round2RedesignPosition,
    summarize_round2_proteinmpnn_output,
)
from lanm.models import ProteinMPNNRound2ShortlistRow, ProteinMPNNRound2UniqueSequenceRow

ONE_TO_THREE = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "K": "LYS",
    "T": "THR",
}


def _single_chain_campaign(
    tmp_path: Path,
    *,
    candidate_id: str,
    seed_rank: int = 1,
    campaign_id: str = "alpha_campaign",
    backbone_id: str = "alpha_backbone",
    topology_class: str = "alpha_topology",
    design_set_name: str = "campaign_ss_only",
    seed_sequence: str,
    residue_start: int = 101,
    insertion_codes: tuple[str, ...] | None = None,
    redesignable_sequence_indices: tuple[int, ...],
    canonical_positions: dict[int, int],
    am1_positions: dict[int, int],
    fixed_context_chains: tuple[str, ...] = (),
    pdb_path: Path | None = None,
) -> ProteinMPNNRound2Campaign:
    insertion_codes = insertion_codes or ("",) * len(seed_sequence)
    designed_residues: list[Round2DesignedResidue] = []
    redesignable_positions: list[Round2RedesignPosition] = []
    for sequence_index, amino_acid in enumerate(seed_sequence, start=1):
        insertion_code = insertion_codes[sequence_index - 1]
        residue_seq = residue_start + sequence_index - 1
        residue_id = f"A{residue_seq}{insertion_code}".strip()
        designed_residues.append(
            Round2DesignedResidue(
                chain_id="A",
                sequence_index=sequence_index,
                residue_seq=residue_seq,
                insertion_code=insertion_code,
                residue_name=ONE_TO_THREE[amino_acid],
                seed_amino_acid=amino_acid,
            )
        )
        if sequence_index in redesignable_sequence_indices:
            redesignable_positions.append(
                Round2RedesignPosition(
                    chain_id="A",
                    sequence_index=sequence_index,
                    residue_id=residue_id,
                    seed_amino_acid=amino_acid,
                    canonical_family_position=canonical_positions[sequence_index],
                    am1_mature_position=am1_positions[sequence_index],
                )
            )

    source_pdb = pdb_path or (tmp_path / f"{candidate_id}.pdb")
    return ProteinMPNNRound2Campaign(
        seed_rank=seed_rank,
        candidate_id=candidate_id,
        campaign_id=campaign_id,
        backbone_id=backbone_id,
        topology_class=topology_class,
        design_set_name=design_set_name,
        designed_chains=("A",),
        fixed_context_chains=fixed_context_chains,
        pdb_path=source_pdb,
        chain_assignment_path=tmp_path / f"{candidate_id}_chain_id.jsonl",
        fixed_positions_path=tmp_path / f"{candidate_id}_fixed_positions.jsonl",
        output_dir=tmp_path / "proteinmpnn_round2_smoke" / candidate_id,
        seed_scaffold_sequence=seed_sequence,
        designed_residues=tuple(designed_residues),
        redesignable_positions=tuple(redesignable_positions),
    )


def _write_round2_fasta(
    path: Path,
    *,
    campaign: ProteinMPNNRound2Campaign,
    generated_rows: tuple[tuple[float, int, float, float, float, str], ...],
) -> None:
    lines = [
        (
            f">{campaign.candidate_id}, "
            f"score=1.0000, "
            f"global_score=1.1000, "
            f"fixed_chains={list(campaign.fixed_context_chains)!r}, "
            f"designed_chains={list(campaign.designed_chains)!r}, "
            f"model_name=v_test, "
            f"git_hash=testhash, "
            f"seed=37"
        ),
        campaign.seed_scaffold_sequence,
    ]
    for temperature, sample_number, score, global_score, seq_recovery, sequence in generated_rows:
        lines.append(
            f">T={temperature:g}, sample={sample_number}, score={score:.4f}, "
            f"global_score={global_score:.4f}, seq_recovery={seq_recovery:.4f}"
        )
        lines.append(sequence)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _unique_row(
    *,
    candidate_id: str,
    seed_rank: int,
    seed_candidate_id: str,
    designed_sequence: str,
    campaign_rank: int,
    best_score: float,
) -> ProteinMPNNRound2UniqueSequenceRow:
    return ProteinMPNNRound2UniqueSequenceRow(
        candidate_id=candidate_id,
        seed_rank=seed_rank,
        seed_candidate_id=seed_candidate_id,
        campaign_id=f"{seed_candidate_id}_campaign",
        backbone_id=f"{seed_candidate_id}_backbone",
        topology_class="alpha",
        design_set_name="campaign_ss_only",
        designed_chains="A",
        designed_sequence=designed_sequence,
        representative_sequence_id=f"{candidate_id}_seq",
        occurrence_count=1,
        campaign_rank=campaign_rank,
        temperature=0.1,
        best_score=best_score,
        best_global_score=best_score + 0.1,
        best_seq_recovery=0.8,
        mutation_count=1,
        mutation_string="A101:A>G",
        redesigned_residue_ids="A101",
        redesigned_canonical_positions="18",
        redesigned_am1_positions="16",
    )


def test_deduplicate_round2_campaign_sequences_collapses_identical_sequences(tmp_path: Path) -> None:
    campaign = _single_chain_campaign(
        tmp_path,
        candidate_id="alpha_seed",
        seed_sequence="ACDE",
        redesignable_sequence_indices=(2, 4),
        canonical_positions={2: 18, 4: 22},
        am1_positions={2: 16, 4: 20},
    )
    fasta_path = campaign.output_dir / "seqs" / f"{campaign.candidate_id}.fa"
    _write_round2_fasta(
        fasta_path,
        campaign=campaign,
        generated_rows=(
            (0.1, 1, 0.5000, 0.7000, 0.7000, "AFDE"),
            (0.1, 2, 0.4000, 0.6500, 0.8500, "AFDE"),
            (0.1, 3, 0.6000, 0.8000, 0.7500, "ACDK"),
        ),
    )

    _, catalog_rows = summarize_round2_proteinmpnn_output(
        campaign=campaign,
        fasta_path=fasta_path,
    )
    result = deduplicate_round2_campaign_sequences(
        campaign=campaign,
        sequence_catalog_rows=catalog_rows,
    )

    assert result.raw_sequence_count == 3
    assert len(result.unique_rows) == 2

    top_row = result.unique_rows[0]
    assert top_row.candidate_id == "alpha_seed_r2u01"
    assert top_row.designed_sequence == "AFDE"
    assert top_row.occurrence_count == 2
    assert top_row.representative_sequence_id == "alpha_seed_T0.1_sample_02"
    assert top_row.best_score == 0.4
    assert top_row.best_global_score == 0.65
    assert top_row.best_seq_recovery == 0.85
    assert top_row.mutation_string == "A102:C>F"
    assert top_row.redesigned_residue_ids == "A102"


def test_extract_round2_mutation_summary_preserves_insertion_codes(tmp_path: Path) -> None:
    campaign = _single_chain_campaign(
        tmp_path,
        candidate_id="beta_seed",
        seed_sequence="AC",
        residue_start=10,
        insertion_codes=("", "A"),
        redesignable_sequence_indices=(2,),
        canonical_positions={2: 24},
        am1_positions={2: 22},
    )

    summary = extract_round2_mutation_summary(
        campaign=campaign,
        designed_sequence="AT",
    )

    assert summary.mutation_count == 1
    assert summary.mutation_string == "A11A:C>T"
    assert summary.redesigned_residue_ids == "A11A"
    assert summary.redesigned_canonical_positions == "24"
    assert summary.redesigned_am1_positions == "22"


def test_build_round2_shortlist_balances_seed_coverage_and_skips_duplicates() -> None:
    unique_rows = (
        _unique_row(
            candidate_id="seed1_r2u01",
            seed_rank=1,
            seed_candidate_id="seed1",
            designed_sequence="AAAA",
            campaign_rank=1,
            best_score=0.10,
        ),
        _unique_row(
            candidate_id="seed1_r2u02",
            seed_rank=1,
            seed_candidate_id="seed1",
            designed_sequence="BBBB",
            campaign_rank=2,
            best_score=0.20,
        ),
        _unique_row(
            candidate_id="seed2_r2u01",
            seed_rank=2,
            seed_candidate_id="seed2",
            designed_sequence="AAAA",
            campaign_rank=1,
            best_score=0.05,
        ),
        _unique_row(
            candidate_id="seed2_r2u02",
            seed_rank=2,
            seed_candidate_id="seed2",
            designed_sequence="CCCC",
            campaign_rank=2,
            best_score=0.30,
        ),
        _unique_row(
            candidate_id="seed3_r2u01",
            seed_rank=3,
            seed_candidate_id="seed3",
            designed_sequence="DDDD",
            campaign_rank=1,
            best_score=0.15,
        ),
        _unique_row(
            candidate_id="seed3_r2u02",
            seed_rank=3,
            seed_candidate_id="seed3",
            designed_sequence="BBBB",
            campaign_rank=2,
            best_score=0.25,
        ),
    )

    shortlist = build_round2_shortlist(unique_rows)

    assert shortlist.skipped_cross_campaign_duplicates == 2
    assert [row.candidate_id for row in shortlist.shortlist_rows] == [
        "seed1_r2u01",
        "seed2_r2u02",
        "seed3_r2u01",
        "seed1_r2u02",
    ]
    assert {row.seed_candidate_id for row in shortlist.shortlist_rows[:3]} == {"seed1", "seed2", "seed3"}
    assert len({row.designed_sequence for row in shortlist.shortlist_rows}) == len(shortlist.shortlist_rows)


def test_export_ligandmpnn_round2_inputs_copies_source_pdb_and_exports_redesign_positions(
    tmp_path: Path,
) -> None:
    source_pdb = tmp_path / "alpha_seed_source.pdb"
    source_pdb.write_text(
        "\n".join(
            (
                "REMARK phase-7c test structure",
                "ATOM      1  N   ALA A  10      11.104  13.207   7.980  1.00 20.00           N  ",
                "ATOM      2  CA  ALA A  10      12.000  13.500   8.500  1.00 20.00           C  ",
                "ATOM      3  N   CYS A  11A     13.104  14.207   9.980  1.00 20.00           N  ",
                "ATOM      4  CA  CYS A  11A     14.000  14.500  10.500  1.00 20.00           C  ",
                "TER",
                "HETATM    5 FE   FE  A 201      15.000  12.000   9.000  1.00 10.00          FE  ",
                "HETATM    6  O   HOH A 301      16.000  11.000   8.000  1.00 30.00           O  ",
                "END",
            )
        )
        + "\n",
        encoding="utf-8",
    )
    campaign = _single_chain_campaign(
        tmp_path,
        candidate_id="alpha_seed",
        seed_sequence="AC",
        residue_start=10,
        insertion_codes=("", "A"),
        redesignable_sequence_indices=(2,),
        canonical_positions={2: 24},
        am1_positions={2: 22},
        pdb_path=source_pdb,
    )
    shortlist_row = ProteinMPNNRound2ShortlistRow(
        shortlist_rank=1,
        candidate_id="alpha_seed_r2u01",
        seed_rank=campaign.seed_rank,
        seed_candidate_id=campaign.candidate_id,
        campaign_id=campaign.campaign_id,
        backbone_id=campaign.backbone_id,
        topology_class=campaign.topology_class,
        design_set_name=campaign.design_set_name,
        designed_chains=campaign.designed_chains_text,
        designed_sequence="AT",
        representative_sequence_id="alpha_seed_T0.1_sample_01",
        occurrence_count=1,
        campaign_rank=1,
        temperature=0.1,
        best_score=0.5,
        best_global_score=0.7,
        best_seq_recovery=0.8,
        mutation_count=1,
        mutation_string="A11A:C>T",
        redesigned_residue_ids="A11A",
        redesigned_canonical_positions="24",
        redesigned_am1_positions="22",
        retention_reason="required seed scaffold coverage",
    )

    manifest_path = tmp_path / "ligandmpnn_round2_input_manifest.csv"
    redesign_positions_path = tmp_path / "ligandmpnn_round2_redesign_positions.csv"
    report_path = tmp_path / "ligandmpnn_round2_inputs.md"
    input_root = tmp_path / "ligandmpnn_round2"
    config_path = tmp_path / "ligandmpnn_round2_inputs.yaml"

    manifest_rows, redesign_rows = export_ligandmpnn_round2_inputs(
        shortlist_rows=(shortlist_row,),
        campaigns_by_seed={campaign.candidate_id: campaign},
        manifest_path=manifest_path,
        redesign_positions_path=redesign_positions_path,
        report_path=report_path,
        input_root=input_root,
        config_path=config_path,
    )

    assert len(manifest_rows) == 1
    assert len(redesign_rows) == 1
    assert manifest_rows[0].redesigned_residue_ids == "A11A"
    assert redesign_rows[0].ligandmpnn_residue_id == "A11A"
    assert redesign_rows[0].mutation_token == "A11A:C>T"

    exported_pdb = input_root / "alpha_seed_r2u01" / "alpha_seed_r2u01.pdb"
    exported_residues = input_root / "alpha_seed_r2u01" / "redesigned_residues.txt"
    assert exported_pdb.read_text(encoding="utf-8") == source_pdb.read_text(encoding="utf-8")
    assert exported_residues.read_text(encoding="utf-8") == "A11A\n"

    with redesign_positions_path.open("r", encoding="utf-8", newline="") as handle:
        exported_rows = list(csv.DictReader(handle))
    assert exported_rows[0]["ligandmpnn_residue_id"] == "A11A"
    assert exported_rows[0]["mutation_token"] == "A11A:C>T"

    assert "alpha_seed_r2u01" in report_path.read_text(encoding="utf-8")
    assert "alpha_seed_r2u01" in config_path.read_text(encoding="utf-8")
