import csv
from pathlib import Path

import pytest

from lanm.analysis.proteinmpnn_smoke import (
    ProteinMPNNSmokeCampaign,
    build_parse_multiple_chains_command,
    build_proteinmpnn_run_command,
    discover_proteinmpnn_smoke_campaigns,
    parse_proteinmpnn_fasta,
    render_proteinmpnn_smoke_markdown,
    summarize_proteinmpnn_output,
)


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def test_build_proteinmpnn_smoke_commands() -> None:
    proteinmpnn_root = Path("/opt/ProteinMPNN")
    parsed_jsonl_path = Path("/tmp/run/parsed_pdbs.jsonl")
    chain_assignment_path = Path("/tmp/run/chain_id.jsonl")
    fixed_positions_path = Path("/tmp/run/fixed_positions.jsonl")
    output_dir = Path("/tmp/run")

    parse_command = build_parse_multiple_chains_command(
        proteinmpnn_root=proteinmpnn_root,
        input_path=Path("/tmp/input_campaign"),
        output_path=parsed_jsonl_path,
        python_executable="/usr/bin/python3",
    )
    assert parse_command == (
        "/usr/bin/python3",
        "/opt/ProteinMPNN/helper_scripts/parse_multiple_chains.py",
        "--input_path",
        "/tmp/input_campaign",
        "--output_path",
        "/tmp/run/parsed_pdbs.jsonl",
    )

    run_command = build_proteinmpnn_run_command(
        proteinmpnn_root=proteinmpnn_root,
        parsed_jsonl_path=parsed_jsonl_path,
        chain_assignment_path=chain_assignment_path,
        fixed_positions_path=fixed_positions_path,
        output_dir=output_dir,
        python_executable="/usr/bin/python3",
    )
    assert run_command == (
        "/usr/bin/python3",
        "/opt/ProteinMPNN/protein_mpnn_run.py",
        "--jsonl_path",
        "/tmp/run/parsed_pdbs.jsonl",
        "--chain_id_jsonl",
        "/tmp/run/chain_id.jsonl",
        "--fixed_positions_jsonl",
        "/tmp/run/fixed_positions.jsonl",
        "--out_folder",
        "/tmp/run",
        "--num_seq_per_target",
        "20",
        "--sampling_temp",
        "0.1 0.15",
        "--seed",
        "37",
        "--batch_size",
        "1",
    )


def test_discover_proteinmpnn_smoke_campaigns_reads_manifest_and_directories(tmp_path: Path) -> None:
    campaigns_dir = tmp_path / "campaigns"
    output_root = tmp_path / "outputs"
    for campaign_id in ("beta_campaign", "alpha_campaign"):
        campaign_dir = campaigns_dir / campaign_id
        campaign_dir.mkdir(parents=True)
        (campaign_dir / f"{campaign_id}.pdb").write_text("ATOM\n", encoding="utf-8")
        (campaign_dir / "chain_id.jsonl").write_text("{}", encoding="utf-8")
        (campaign_dir / "fixed_positions.jsonl").write_text("{}", encoding="utf-8")

    campaign_manifest_path = tmp_path / "design_campaign_manifest.csv"
    _write_csv(
        campaign_manifest_path,
        [
            {
                "campaign_id": "alpha_campaign",
                "backbone_id": "alpha_backbone",
                "design_set_name": "campaign_ss_only",
                "designed_chains": "A",
                "fixed_context_chains": "",
                "exported_pdb_path": str(campaigns_dir / "alpha_campaign" / "alpha_campaign.pdb"),
                "chain_assignment_path": str(campaigns_dir / "alpha_campaign" / "chain_id.jsonl"),
                "fixed_positions_path": str(campaigns_dir / "alpha_campaign" / "fixed_positions.jsonl"),
                "designable_residue_count": "16",
                "fixed_residue_count": "89",
                "designable_canonical_positions": "13,18",
                "designable_am1_positions": "12,16",
            },
            {
                "campaign_id": "beta_campaign",
                "backbone_id": "beta_backbone",
                "design_set_name": "campaign_if_only",
                "designed_chains": "A,C",
                "fixed_context_chains": "B",
                "exported_pdb_path": str(campaigns_dir / "beta_campaign" / "beta_campaign.pdb"),
                "chain_assignment_path": str(campaigns_dir / "beta_campaign" / "chain_id.jsonl"),
                "fixed_positions_path": str(campaigns_dir / "beta_campaign" / "fixed_positions.jsonl"),
                "designable_residue_count": "10",
                "fixed_residue_count": "100",
                "designable_canonical_positions": "18,27",
                "designable_am1_positions": "16,25",
            },
            {
                "campaign_id": "manifest_only",
                "backbone_id": "alpha_backbone",
                "design_set_name": "campaign_ss_only",
                "designed_chains": "A",
                "fixed_context_chains": "",
                "exported_pdb_path": str(campaigns_dir / "manifest_only" / "manifest_only.pdb"),
                "chain_assignment_path": str(campaigns_dir / "manifest_only" / "chain_id.jsonl"),
                "fixed_positions_path": str(campaigns_dir / "manifest_only" / "fixed_positions.jsonl"),
                "designable_residue_count": "16",
                "fixed_residue_count": "89",
                "designable_canonical_positions": "13,18",
                "designable_am1_positions": "12,16",
            },
        ],
    )

    backbone_manifest_path = tmp_path / "design_backbone_manifest.csv"
    _write_csv(
        backbone_manifest_path,
        [
            {
                "backbone_id": "alpha_backbone",
                "structure_id": "8FNS",
                "source_kind": "csv_atom_table",
                "source_path": "source_a",
                "experimental_method": "X-RAY DIFFRACTION",
                "selected_chains": "A",
                "design_chains": "A",
                "fixed_context_chains": "",
                "selected_residue_count": "105",
                "design_chain_residue_count": "105",
                "design_chain_metal_site_count": "4",
                "representative_note": "note",
                "backbone_pdb_path": "backbone_a.pdb",
            },
            {
                "backbone_id": "beta_backbone",
                "structure_id": "8FNR",
                "source_kind": "cif",
                "source_path": "source_b",
                "experimental_method": "X-RAY DIFFRACTION",
                "selected_chains": "A,B,C",
                "design_chains": "A,C",
                "fixed_context_chains": "B",
                "selected_residue_count": "300",
                "design_chain_residue_count": "200",
                "design_chain_metal_site_count": "4",
                "representative_note": "note",
                "backbone_pdb_path": "backbone_b.pdb",
            },
        ],
    )

    campaigns = discover_proteinmpnn_smoke_campaigns(
        campaign_manifest_path=campaign_manifest_path,
        backbone_manifest_path=backbone_manifest_path,
        campaigns_dir=campaigns_dir,
        output_root=output_root,
    )

    assert [campaign.campaign_id for campaign in campaigns] == ["alpha_campaign", "beta_campaign"]
    assert campaigns[0].backbone_id == "alpha_backbone"
    assert campaigns[0].backbone_structure_id == "8FNS"
    assert campaigns[1].designed_chains == ("A", "C")
    assert campaigns[1].fixed_context_chains == ("B",)
    assert campaigns[1].output_dir == output_root / "beta_campaign"


def test_parse_and_summarize_proteinmpnn_output(tmp_path: Path) -> None:
    fasta_path = tmp_path / "seqs" / "example_campaign.fa"
    fasta_path.parent.mkdir(parents=True)
    fasta_path.write_text(
        "\n".join(
            [
                ">example_campaign, score=1.2000, global_score=1.5000, fixed_chains=['B'], designed_chains=['A'], model_name=v_48_020, git_hash=abc123, seed=37",
                "AAAA",
                ">T=0.1, sample=1, score=0.7000, global_score=0.9000, seq_recovery=0.5000",
                "AAAA",
                ">T=0.1, sample=2, score=0.7100, global_score=0.9100, seq_recovery=0.7500",
                "AAAT",
                ">T=0.15, sample=1, score=0.7200, global_score=0.9200, seq_recovery=0.2500",
                "AATT",
                "",
            ]
        ),
        encoding="utf-8",
    )

    campaign = ProteinMPNNSmokeCampaign(
        campaign_id="example_campaign",
        backbone_id="am1_mex_8fns_chain_a",
        backbone_structure_id="8FNS",
        design_set_name="campaign_ss_only",
        designed_chains=("A",),
        fixed_context_chains=("B",),
        pdb_path=tmp_path / "example_campaign.pdb",
        chain_assignment_path=tmp_path / "chain_id.jsonl",
        fixed_positions_path=tmp_path / "fixed_positions.jsonl",
        output_dir=tmp_path / "outputs" / "example_campaign",
    )

    parsed_output = parse_proteinmpnn_fasta(fasta_path)
    summary_row, catalog_rows = summarize_proteinmpnn_output(
        campaign=campaign,
        parsed_output=parsed_output,
        fasta_path=fasta_path,
    )

    assert parsed_output.native_record.name == "example_campaign"
    assert parsed_output.native_record.fixed_chains == ("B",)
    assert len(parsed_output.generated_records) == 3

    assert summary_row.generated_sequence_count == 3
    assert summary_row.unique_sequence_count == 3
    assert summary_row.unique_sequence_fraction == 1.0
    assert summary_row.sequence_length == 4
    assert summary_row.native_score == 1.2
    assert summary_row.mean_pairwise_identity == pytest.approx(2.0 / 3.0)
    assert summary_row.min_pairwise_identity == pytest.approx(0.5)
    assert summary_row.max_pairwise_identity == pytest.approx(0.75)
    assert summary_row.output_fasta_path.endswith("example_campaign.fa")

    assert [row.sequence_id for row in catalog_rows] == [
        "example_campaign_T0.1_sample_01",
        "example_campaign_T0.1_sample_02",
        "example_campaign_T0.15_sample_01",
    ]
    assert catalog_rows[0].designed_sequence == "AAAA"
    assert catalog_rows[1].score == pytest.approx(0.71)
    assert catalog_rows[2].seq_recovery == pytest.approx(0.25)

    markdown = render_proteinmpnn_smoke_markdown(
        proteinmpnn_root=Path("/home/ashak/apps/ProteinMPNN"),
        summary_rows=(summary_row,),
    )
    assert "example_campaign" in markdown
    assert "3/3 unique" in markdown
