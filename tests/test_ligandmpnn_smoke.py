import csv
from pathlib import Path

import pytest

from lanm.analysis.ligandmpnn_smoke import (
    LigandMPNNSmokeCandidate,
    build_ligandmpnn_run_command,
    discover_ligandmpnn_smoke_candidates,
    parse_ligandmpnn_fasta,
    render_ligandmpnn_smoke_markdown,
    summarize_ligandmpnn_output,
)


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def test_discover_ligandmpnn_smoke_candidates_reads_manifest_and_inputs(tmp_path: Path) -> None:
    inputs_dir = tmp_path / "inputs"
    output_root = tmp_path / "outputs"

    alpha_dir = inputs_dir / "alpha_candidate"
    alpha_dir.mkdir(parents=True)
    (alpha_dir / "alpha_candidate.pdb").write_text("ATOM\nEND\n", encoding="utf-8")
    (alpha_dir / "redesigned_residues.txt").write_text("A35 C8\n", encoding="utf-8")

    beta_dir = inputs_dir / "beta_candidate"
    beta_dir.mkdir(parents=True)
    (beta_dir / "beta_candidate.pdb").write_text("ATOM\nEND\n", encoding="utf-8")
    (beta_dir / "redesigned_residues.txt").write_text("A10\n", encoding="utf-8")

    manifest_path = tmp_path / "ligandmpnn_input_manifest.csv"
    _write_csv(
        manifest_path,
        [
            {
                "shortlist_rank": "2",
                "candidate_id": "beta_candidate",
                "campaign_id": "beta_campaign",
                "backbone_id": "beta_backbone",
                "source_structure_id": "8FNR",
                "source_kind": "cif",
                "source_path": "source_b",
                "preserved_chains": "A,B",
                "designed_chains": "A",
                "fixed_context_chains": "B",
                "preserved_metal_identity": "DY",
                "preserved_metal_site_count": "4",
                "preserved_solvent_residue_count": "12",
                "redesigned_residue_count": "1",
                "pdb_path": str(beta_dir / "beta_candidate.pdb"),
                "redesigned_residues_path": str(beta_dir / "redesigned_residues.txt"),
            },
            {
                "shortlist_rank": "1",
                "candidate_id": "alpha_candidate",
                "campaign_id": "alpha_campaign",
                "backbone_id": "alpha_backbone",
                "source_structure_id": "8FNS",
                "source_kind": "csv_atom_table",
                "source_path": "source_a",
                "preserved_chains": "A,C",
                "designed_chains": "A,C",
                "fixed_context_chains": "B",
                "preserved_metal_identity": "ND",
                "preserved_metal_site_count": "2",
                "preserved_solvent_residue_count": "5",
                "redesigned_residue_count": "2",
                "pdb_path": str(alpha_dir / "alpha_candidate.pdb"),
                "redesigned_residues_path": str(alpha_dir / "redesigned_residues.txt"),
            },
        ],
    )

    candidates = discover_ligandmpnn_smoke_candidates(
        manifest_path=manifest_path,
        output_root=output_root,
    )

    assert [candidate.candidate_id for candidate in candidates] == ["alpha_candidate", "beta_candidate"]
    assert candidates[0].designed_chains == ("A", "C")
    assert candidates[0].redesigned_residue_ids == ("A35", "C8")
    assert candidates[0].preserved_metal_identity == "ND"
    assert candidates[0].output_dir == output_root / "alpha_candidate"
    assert candidates[1].designed_chains_text == "A"
    assert candidates[1].redesigned_residues_text == "A10"


def test_build_ligandmpnn_run_command() -> None:
    candidate = LigandMPNNSmokeCandidate(
        shortlist_rank=1,
        candidate_id="example_candidate",
        campaign_id="example_campaign",
        backbone_id="example_backbone",
        source_structure_id="8FNS",
        designed_chains=("A", "C"),
        preserved_metal_identity="ND",
        redesigned_residue_ids=("A35", "C8"),
        pdb_path=Path("/tmp/example_candidate.pdb"),
        redesigned_residues_path=Path("/tmp/redesigned_residues.txt"),
        output_dir=Path("/tmp/results/example_candidate"),
    )

    command = build_ligandmpnn_run_command(
        ligandmpnn_root=Path("/opt/LigandMPNN"),
        candidate=candidate,
        python_executable="/usr/bin/python3",
    )

    assert command == (
        "/usr/bin/python3",
        "/opt/LigandMPNN/run.py",
        "--model_type",
        "ligand_mpnn",
        "--seed",
        "37",
        "--pdb_path",
        "/tmp/example_candidate.pdb",
        "--out_folder",
        "/tmp/results/example_candidate",
        "--chains_to_design",
        "A,C",
        "--redesigned_residues",
        "A35 C8",
        "--ligand_mpnn_use_atom_context",
        "1",
        "--ligand_mpnn_use_side_chain_context",
        "1",
        "--pack_side_chains",
        "1",
        "--number_of_packs_per_design",
        "1",
        "--pack_with_ligand_context",
        "1",
        "--batch_size",
        "2",
        "--number_of_batches",
        "1",
        "--temperature",
        "0.1",
        "--verbose",
        "0",
    )


def test_parse_and_summarize_ligandmpnn_output(tmp_path: Path) -> None:
    fasta_path = tmp_path / "results" / "example_candidate" / "seqs" / "example_candidate.fa"
    fasta_path.parent.mkdir(parents=True)
    fasta_path.write_text(
        "\n".join(
            [
                ">example_candidate, T=0.1, seed=37, num_res=2, num_ligand_res=2, use_ligand_context=True, ligand_cutoff_distance=8.0, batch_size=1, number_of_batches=2, model_path=./model_params/ligandmpnn_v_32_010_25.pt",
                "AAAA",
                ">example_candidate, id=1, T=0.1, seed=37, overall_confidence=0.8000, ligand_confidence=0.7000, seq_rec=1.0000",
                "AAAA",
                ">example_candidate, id=2, T=0.1, seed=37, overall_confidence=0.9000, ligand_confidence=0.6000, seq_rec=0.7500",
                "AAAT",
                ">example_candidate, id=3, T=0.1, seed=37, overall_confidence=0.7000, ligand_confidence=0.8000, seq_rec=0.5000",
                "AATT",
                "",
            ]
        ),
        encoding="utf-8",
    )

    candidate = LigandMPNNSmokeCandidate(
        shortlist_rank=1,
        candidate_id="example_candidate",
        campaign_id="example_campaign",
        backbone_id="am1_mex_8fns_chain_a",
        source_structure_id="8FNS",
        designed_chains=("A",),
        preserved_metal_identity="ND",
        redesigned_residue_ids=("A35", "A44"),
        pdb_path=tmp_path / "example_candidate.pdb",
        redesigned_residues_path=tmp_path / "redesigned_residues.txt",
        output_dir=tmp_path / "results" / "example_candidate",
    )

    parsed_output = parse_ligandmpnn_fasta(fasta_path)
    summary_row, catalog_rows = summarize_ligandmpnn_output(
        candidate=candidate,
        parsed_output=parsed_output,
        fasta_path=fasta_path,
    )

    assert parsed_output.native_record.name == "example_candidate"
    assert parsed_output.native_record.use_ligand_context is True
    assert len(parsed_output.generated_records) == 3

    assert summary_row.generated_sequence_count == 3
    assert summary_row.unique_sequence_count == 3
    assert summary_row.unique_sequence_fraction == 1.0
    assert summary_row.sequence_length == 4
    assert summary_row.mean_overall_confidence == pytest.approx(0.8)
    assert summary_row.mean_ligand_confidence == pytest.approx(0.7)
    assert summary_row.mean_pairwise_identity == pytest.approx(2.0 / 3.0)
    assert summary_row.min_pairwise_identity == pytest.approx(0.5)
    assert summary_row.max_pairwise_identity == pytest.approx(0.75)
    assert summary_row.redesigned_residues == "A35,A44"
    assert summary_row.output_fasta_path.endswith("example_candidate.fa")

    assert [row.sequence_id for row in catalog_rows] == [
        "example_candidate_design_01",
        "example_candidate_design_02",
        "example_candidate_design_03",
    ]
    assert catalog_rows[0].designed_sequence == "AAAA"
    assert catalog_rows[1].overall_confidence == pytest.approx(0.9)
    assert catalog_rows[2].packed_pdb_path.endswith("example_candidate_packed_3_1.pdb")

    markdown = render_ligandmpnn_smoke_markdown(
        ligandmpnn_root=Path("/home/ashak/apps/LigandMPNN"),
        summary_rows=(summary_row,),
    )
    assert "example_candidate" in markdown
    assert "am1_mex_8fns_chain_a" in markdown
    assert "ND" in markdown
    assert "A35,A44" in markdown
    assert "3/3 unique" in markdown
