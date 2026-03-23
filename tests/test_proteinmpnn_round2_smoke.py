import csv
import json
import textwrap
from pathlib import Path

from lanm.analysis.proteinmpnn_round2_smoke import (
    ProteinMPNNRound2Campaign,
    Round2DesignedResidue,
    Round2RedesignPosition,
    build_round2_parse_command,
    build_round2_run_command,
    discover_round2_proteinmpnn_campaigns,
    render_proteinmpnn_round2_smoke_markdown,
    summarize_round2_proteinmpnn_output,
)
from lanm.models import AtomRecord
from lanm.structure.pdb import render_protein_pdb

ONE_TO_THREE = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "G": "GLY",
    "K": "LYS",
    "Q": "GLN",
    "T": "THR",
    "V": "VAL",
}


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _write_jsonl_object(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload) + "\n", encoding="utf-8")


def _write_round2_config(path: Path, payload: str) -> None:
    path.write_text(textwrap.dedent(payload).strip() + "\n", encoding="utf-8")


def _write_test_pdb(path: Path, chain_sequences: dict[str, str], chain_starts: dict[str, int]) -> None:
    atoms: list[AtomRecord] = []
    atom_serial = 1
    for chain_id, sequence in chain_sequences.items():
        start = chain_starts[chain_id]
        for offset, amino_acid in enumerate(sequence):
            atoms.append(
                AtomRecord(
                    structure_id=path.stem.upper(),
                    record_type="ATOM",
                    atom_serial=atom_serial,
                    atom_name="CA",
                    alt_loc="",
                    residue_name=ONE_TO_THREE[amino_acid],
                    chain_id=chain_id,
                    residue_seq=start + offset,
                    insertion_code="",
                    x=float(atom_serial),
                    y=0.0,
                    z=0.0,
                    occupancy=1.0,
                    b_factor=1.0,
                    element="C",
                    charge="",
                )
            )
            atom_serial += 1
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        render_protein_pdb(atoms, selected_chain_ids=tuple(chain_sequences.keys())),
        encoding="utf-8",
    )


def test_discover_round2_proteinmpnn_campaigns_reads_phase_7a_exports(tmp_path: Path) -> None:
    campaigns_dir = tmp_path / "results" / "design_inputs" / "proteinmpnn_round2" / "campaigns"
    output_root = tmp_path / "results" / "proteinmpnn_round2_smoke"
    alpha_dir = campaigns_dir / "alpha_seed"
    beta_dir = campaigns_dir / "beta_seed"
    alpha_pdb = alpha_dir / "alpha_seed.pdb"
    beta_pdb = beta_dir / "beta_seed.pdb"
    _write_test_pdb(alpha_pdb, {"A": "ACDE"}, {"A": 101})
    _write_test_pdb(beta_pdb, {"A": "AK", "B": "GG", "C": "GT"}, {"A": 21, "B": 51, "C": 81})

    _write_jsonl_object(alpha_dir / "chain_id.jsonl", {"alpha_seed": [["A"], []]})
    _write_jsonl_object(alpha_dir / "fixed_positions.jsonl", {"alpha_seed": {"A": [1, 3]}})
    _write_jsonl_object(beta_dir / "chain_id.jsonl", {"beta_seed": [["A", "C"], ["B"]]})
    _write_jsonl_object(beta_dir / "fixed_positions.jsonl", {"beta_seed": {"A": [2], "C": [1]}})

    config_path = tmp_path / "config" / "redesign_round2.yaml"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    _write_round2_config(
        config_path,
        f"""
        version: 1
        phase: 7A
        selection_policy:
          seed_candidate_ids:
            - alpha_seed
            - beta_seed
        seeds:
          alpha_seed:
            design_set_name: campaign_ss_only
            designed_chains: [A]
            fixed_context_chains: []
            exported_pdb_path: {alpha_pdb}
            chain_assignment_path: {alpha_dir / 'chain_id.jsonl'}
            fixed_positions_path: {alpha_dir / 'fixed_positions.jsonl'}
          beta_seed:
            design_set_name: campaign_if_only
            designed_chains: [A, C]
            fixed_context_chains: [B]
            exported_pdb_path: {beta_pdb}
            chain_assignment_path: {beta_dir / 'chain_id.jsonl'}
            fixed_positions_path: {beta_dir / 'fixed_positions.jsonl'}
        """,
    )

    seed_table_path = tmp_path / "results" / "tables" / "redesign_round2_seeds.csv"
    seed_table_path.parent.mkdir(parents=True, exist_ok=True)
    _write_csv(
        seed_table_path,
        [
            {
                "seed_rank": 1,
                "candidate_id": "alpha_seed",
                "campaign_id": "alpha_campaign",
                "backbone_id": "alpha_backbone",
                "topology_class": "alpha_topology",
                "designed_chains": "A",
                "fixed_context_chains": "",
                "exported_pdb_path": str(alpha_pdb),
                "chain_assignment_path": str(alpha_dir / "chain_id.jsonl"),
                "fixed_positions_path": str(alpha_dir / "fixed_positions.jsonl"),
                "redesignable_residue_ids": "A102,A104",
                "redesignable_canonical_positions": "10,12",
                "redesignable_am1_positions": "8,10",
            },
            {
                "seed_rank": 2,
                "candidate_id": "beta_seed",
                "campaign_id": "beta_campaign",
                "backbone_id": "beta_backbone",
                "topology_class": "beta_topology",
                "designed_chains": "A,C",
                "fixed_context_chains": "B",
                "exported_pdb_path": str(beta_pdb),
                "chain_assignment_path": str(beta_dir / "chain_id.jsonl"),
                "fixed_positions_path": str(beta_dir / "fixed_positions.jsonl"),
                "redesignable_residue_ids": "A21,C82",
                "redesignable_canonical_positions": "18,44",
                "redesignable_am1_positions": "16,40",
            },
        ],
    )

    position_table_path = tmp_path / "results" / "tables" / "redesign_round2_positions.csv"
    _write_csv(
        position_table_path,
        [
            {
                "candidate_id": "alpha_seed",
                "chain_id": "A",
                "sequence_index": 2,
                "residue_id": "A102",
                "seed_amino_acid": "C",
                "canonical_family_position": 10,
                "am1_mature_position": 8,
            },
            {
                "candidate_id": "alpha_seed",
                "chain_id": "A",
                "sequence_index": 4,
                "residue_id": "A104",
                "seed_amino_acid": "E",
                "canonical_family_position": 12,
                "am1_mature_position": 10,
            },
            {
                "candidate_id": "beta_seed",
                "chain_id": "A",
                "sequence_index": 1,
                "residue_id": "A21",
                "seed_amino_acid": "A",
                "canonical_family_position": 18,
                "am1_mature_position": 16,
            },
            {
                "candidate_id": "beta_seed",
                "chain_id": "C",
                "sequence_index": 2,
                "residue_id": "C82",
                "seed_amino_acid": "T",
                "canonical_family_position": 44,
                "am1_mature_position": 40,
            },
        ],
    )

    campaigns = discover_round2_proteinmpnn_campaigns(
        redesign_round2_config_path=config_path,
        round2_seed_table_path=seed_table_path,
        round2_position_table_path=position_table_path,
        campaigns_dir=campaigns_dir,
        output_root=output_root,
    )

    assert [campaign.candidate_id for campaign in campaigns] == ["alpha_seed", "beta_seed"]
    assert campaigns[0].seed_scaffold_sequence == "ACDE"
    assert campaigns[0].redesignable_residue_ids_text == "A102,A104"
    assert campaigns[1].designed_chains == ("A", "C")
    assert campaigns[1].fixed_context_chains == ("B",)
    assert campaigns[1].seed_scaffold_sequence == "AK/GT"
    assert campaigns[1].redesignable_canonical_positions_text == "18,44"
    assert campaigns[1].output_dir == output_root / "beta_seed"


def test_build_round2_proteinmpnn_commands() -> None:
    campaign = ProteinMPNNRound2Campaign(
        seed_rank=1,
        candidate_id="alpha_seed",
        campaign_id="alpha_campaign",
        backbone_id="alpha_backbone",
        topology_class="alpha_topology",
        design_set_name="campaign_ss_only",
        designed_chains=("A",),
        fixed_context_chains=(),
        pdb_path=Path("/tmp/alpha_seed/alpha_seed.pdb"),
        chain_assignment_path=Path("/tmp/alpha_seed/chain_id.jsonl"),
        fixed_positions_path=Path("/tmp/alpha_seed/fixed_positions.jsonl"),
        output_dir=Path("/tmp/output/alpha_seed"),
        seed_scaffold_sequence="ACDE",
        designed_residues=(),
        redesignable_positions=(),
    )

    parse_command = build_round2_parse_command(
        campaign=campaign,
        proteinmpnn_root=Path("/opt/ProteinMPNN"),
        python_executable="/usr/bin/python3",
    )
    assert parse_command == (
        "/usr/bin/python3",
        "/opt/ProteinMPNN/helper_scripts/parse_multiple_chains.py",
        "--input_path",
        "/tmp/alpha_seed",
        "--output_path",
        "/tmp/output/alpha_seed/parsed_pdbs.jsonl",
    )

    run_command = build_round2_run_command(
        campaign=campaign,
        proteinmpnn_root=Path("/opt/ProteinMPNN"),
        parsed_jsonl_path=Path("/tmp/output/alpha_seed/parsed_pdbs.jsonl"),
        python_executable="/usr/bin/python3",
    )
    assert run_command == (
        "/usr/bin/python3",
        "/opt/ProteinMPNN/protein_mpnn_run.py",
        "--jsonl_path",
        "/tmp/output/alpha_seed/parsed_pdbs.jsonl",
        "--chain_id_jsonl",
        "/tmp/alpha_seed/chain_id.jsonl",
        "--fixed_positions_jsonl",
        "/tmp/alpha_seed/fixed_positions.jsonl",
        "--out_folder",
        "/tmp/output/alpha_seed",
        "--num_seq_per_target",
        "20",
        "--sampling_temp",
        "0.1 0.15",
        "--seed",
        "37",
        "--batch_size",
        "1",
    )


def test_summarize_round2_proteinmpnn_output_reports_unique_sequences_and_mutation_counts(tmp_path: Path) -> None:
    fasta_path = tmp_path / "seqs" / "alpha_seed.fa"
    fasta_path.parent.mkdir(parents=True)
    fasta_path.write_text(
        "\n".join(
            [
                ">alpha_seed, score=1.2000, global_score=1.5000, fixed_chains=[], designed_chains=['A'], model_name=v_48_020, git_hash=abc123, seed=37",
                "ACDE",
                ">T=0.1, sample=1, score=0.7000, global_score=0.9000, seq_recovery=1.0000",
                "ACDE",
                ">T=0.1, sample=2, score=0.7100, global_score=0.9100, seq_recovery=0.7500",
                "AVDE",
                ">T=0.15, sample=1, score=0.7200, global_score=0.9200, seq_recovery=0.5000",
                "AVDQ",
                "",
            ]
        ),
        encoding="utf-8",
    )

    campaign = ProteinMPNNRound2Campaign(
        seed_rank=1,
        candidate_id="alpha_seed",
        campaign_id="alpha_campaign",
        backbone_id="alpha_backbone",
        topology_class="alpha_topology",
        design_set_name="campaign_ss_only",
        designed_chains=("A",),
        fixed_context_chains=(),
        pdb_path=tmp_path / "alpha_seed.pdb",
        chain_assignment_path=tmp_path / "chain_id.jsonl",
        fixed_positions_path=tmp_path / "fixed_positions.jsonl",
        output_dir=tmp_path / "outputs" / "alpha_seed",
        seed_scaffold_sequence="ACDE",
        designed_residues=(
            Round2DesignedResidue("A", 1, 101, "", "ALA", "A"),
            Round2DesignedResidue("A", 2, 102, "", "CYS", "C"),
            Round2DesignedResidue("A", 3, 103, "", "ASP", "D"),
            Round2DesignedResidue("A", 4, 104, "", "GLU", "E"),
        ),
        redesignable_positions=(
            Round2RedesignPosition("A", 2, "A102", "C", 10, 8),
            Round2RedesignPosition("A", 4, "A104", "E", 12, 10),
        ),
    )

    summary_row, catalog_rows = summarize_round2_proteinmpnn_output(
        campaign=campaign,
        fasta_path=fasta_path,
    )

    assert summary_row.seed_scaffold_sequence == "ACDE"
    assert summary_row.generated_sequence_count == 3
    assert summary_row.unique_sequence_count == 3
    assert summary_row.mutation_count_distribution == "0:1,1:1,2:1"
    assert summary_row.redesignable_residue_ids == "A102,A104"

    assert [row.sequence_id for row in catalog_rows] == [
        "alpha_seed_T0.1_sample_01",
        "alpha_seed_T0.1_sample_02",
        "alpha_seed_T0.15_sample_01",
    ]
    assert [row.mutation_count_vs_seed for row in catalog_rows] == [0, 1, 2]
    assert catalog_rows[1].mutation_string_vs_seed == "A102:C>V"
    assert catalog_rows[2].mutated_residue_ids == "A102,A104"

    markdown = render_proteinmpnn_round2_smoke_markdown(
        proteinmpnn_root=Path("/home/ashak/apps/ProteinMPNN"),
        summary_rows=(summary_row,),
    )
    assert "alpha_seed" in markdown
    assert "A102,A104" in markdown
    assert "0:1,1:1,2:1" in markdown
