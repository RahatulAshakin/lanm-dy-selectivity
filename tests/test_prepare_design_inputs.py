import csv
import json
from pathlib import Path

from lanm.analysis.design_campaigns import (
    build_design_campaign_artifacts,
    write_design_campaign_outputs,
)
from lanm.paths import (
    CROSS_TEMPLATE_RESIDUE_ALIGNMENT_PATH,
    DESIGN_MASK_CANDIDATES_PATH,
    DESIGN_MASKS_PATH,
    TEMPLATE_CHAIN_SUMMARY_PATH,
)


def _build_artifacts():
    return build_design_campaign_artifacts(
        design_masks_path=DESIGN_MASKS_PATH,
        design_mask_candidates_path=DESIGN_MASK_CANDIDATES_PATH,
        cross_template_alignment_path=CROSS_TEMPLATE_RESIDUE_ALIGNMENT_PATH,
        template_chain_summary_path=TEMPLATE_CHAIN_SUMMARY_PATH,
    )


def test_build_design_campaign_artifacts_derives_expected_campaigns() -> None:
    artifacts = _build_artifacts()
    design_set_map = {
        design_set.name: design_set
        for design_set in artifacts.design_sets
    }
    backbone_map = {
        backbone.backbone_id: backbone
        for backbone in artifacts.backbones
    }
    campaign_map = {
        campaign.campaign_id: campaign
        for campaign in artifacts.campaigns
    }

    assert list(design_set_map) == [
        "hard_fixed",
        "campaign_ss_only",
        "campaign_ss_plus_if",
        "campaign_if_only",
    ]
    assert design_set_map["hard_fixed"].canonical_positions == (
        1, 11, 14, 15, 16, 17, 19, 20, 21, 23, 25, 26, 28, 30, 31, 33, 41,
        42, 43, 45, 46, 47, 49, 52, 55, 66, 67, 68, 70, 71, 72, 74, 77, 90,
        91, 92, 94, 95, 96, 101, 105, 107, 110,
    )
    assert design_set_map["campaign_ss_only"].canonical_positions == (
        13, 18, 22, 24, 40, 44, 48, 65, 69, 73, 83, 89, 93, 97, 98, 100,
    )
    assert design_set_map["campaign_ss_plus_if"].canonical_positions == (
        13, 18, 22, 24, 27, 38, 39, 40, 44, 48, 65, 69, 73, 75, 76, 79, 80,
        83, 89, 93, 97, 98, 100,
    )
    assert design_set_map["campaign_if_only"].canonical_positions == (
        18, 27, 38, 39, 44, 75, 76, 79, 80, 83,
    )

    assert list(backbone_map) == [
        "am1_mex_8fns_chain_a",
        "hans_pocket_8fnr_chain_a",
        "hans_interface_8fnr_a_b_c_d",
    ]
    assert backbone_map["am1_mex_8fns_chain_a"].design_chains == ("A",)
    assert backbone_map["hans_pocket_8fnr_chain_a"].design_chains == ("A",)
    assert backbone_map["hans_pocket_8fnr_chain_a"].design_chain_metal_site_count == 4
    assert backbone_map["hans_interface_8fnr_a_b_c_d"].selected_chains == ("A", "B", "C", "D")
    assert backbone_map["hans_interface_8fnr_a_b_c_d"].fixed_context_chains == ("B", "C", "D")

    assert list(campaign_map) == [
        "am1_mex_ss_only",
        "hans_pocket_ss_only",
        "hans_interface_ss_plus_if",
        "hans_interface_if_only",
    ]
    assert campaign_map["am1_mex_ss_only"].designable_residue_count == 16
    assert campaign_map["am1_mex_ss_only"].fixed_residue_count == 89
    assert campaign_map["hans_pocket_ss_only"].designable_residue_count == 16
    assert campaign_map["hans_interface_ss_plus_if"].designable_residue_count == 23
    assert campaign_map["hans_interface_if_only"].designable_residue_count == 10
    assert len(artifacts.campaign_position_rows) == 435

    designable_row = next(
        row
        for row in artifacts.campaign_position_rows
        if row.campaign_id == "am1_mex_ss_only"
        and row.chain_id == "A"
        and row.sequence_index == 6
    )
    assert designable_row.residue_seq == 34
    assert designable_row.canonical_family_position == 13
    assert designable_row.designable is True
    assert designable_row.position_state == "designable"
    assert designable_row.mutable_second_sphere is True

    hard_fixed_row = next(
        row
        for row in artifacts.campaign_position_rows
        if row.campaign_id == "am1_mex_ss_only"
        and row.chain_id == "A"
        and row.sequence_index == 7
    )
    assert hard_fixed_row.residue_seq == 35
    assert hard_fixed_row.canonical_family_position == 14
    assert hard_fixed_row.designable is False
    assert hard_fixed_row.hard_fixed is True
    assert hard_fixed_row.position_reason == "hard_fixed"


def test_write_design_campaign_outputs_exports_fixed_positions_and_metadata(tmp_path: Path) -> None:
    artifacts = _build_artifacts()
    report_path = tmp_path / "design_campaigns.md"
    campaign_positions_path = tmp_path / "design_campaign_positions.csv"
    backbone_manifest_path = tmp_path / "design_backbone_manifest.csv"
    campaign_manifest_path = tmp_path / "design_campaign_manifest.csv"
    config_path = tmp_path / "design_campaigns.yaml"
    proteinmpnn_root = tmp_path / "proteinmpnn"

    write_design_campaign_outputs(
        artifacts=artifacts,
        report_path=report_path,
        campaign_positions_path=campaign_positions_path,
        backbone_manifest_path=backbone_manifest_path,
        campaign_manifest_path=campaign_manifest_path,
        config_path=config_path,
        proteinmpnn_root=proteinmpnn_root,
    )

    assert report_path.exists()
    assert campaign_positions_path.exists()
    assert backbone_manifest_path.exists()
    assert campaign_manifest_path.exists()
    assert config_path.exists()
    assert (proteinmpnn_root / "backbones" / "am1_mex_8fns_chain_a.pdb").exists()
    assert (proteinmpnn_root / "backbones" / "hans_interface_8fnr_a_b_c_d.pdb").exists()
    assert (proteinmpnn_root / "campaigns" / "am1_mex_ss_only" / "am1_mex_ss_only.pdb").exists()
    assert (proteinmpnn_root / "campaigns" / "hans_interface_ss_plus_if" / "chain_id.jsonl").exists()
    assert (proteinmpnn_root / "campaigns" / "hans_interface_if_only" / "fixed_positions.jsonl").exists()

    chain_assignment = json.loads(
        (proteinmpnn_root / "campaigns" / "hans_interface_ss_plus_if" / "chain_id.jsonl").read_text(encoding="utf-8")
    )
    assert chain_assignment == {"hans_interface_ss_plus_if": [["A"], ["B", "C", "D"]]}

    fixed_positions = json.loads(
        (proteinmpnn_root / "campaigns" / "am1_mex_ss_only" / "fixed_positions.jsonl").read_text(encoding="utf-8")
    )
    am1_fixed = fixed_positions["am1_mex_ss_only"]["A"]
    assert len(am1_fixed) == 89
    assert 6 not in am1_fixed
    assert 7 in am1_fixed
    assert 10 not in am1_fixed
    assert 11 in am1_fixed

    interface_pdb_path = proteinmpnn_root / "campaigns" / "hans_interface_ss_plus_if" / "hans_interface_ss_plus_if.pdb"
    atom_chain_ids = {
        line[21].strip()
        for line in interface_pdb_path.read_text(encoding="utf-8").splitlines()
        if line.startswith("ATOM")
    }
    assert atom_chain_ids == {"A", "B", "C", "D"}

    with campaign_positions_path.open("r", encoding="utf-8", newline="") as handle:
        position_rows = list(csv.DictReader(handle))
    row = next(
        item
        for item in position_rows
        if item["campaign_id"] == "hans_interface_if_only"
        and item["canonical_family_position"] == "76"
    )
    assert row["designable"] == "True"
    assert row["position_state"] == "designable"
    assert row["chain_id"] == "A"

    row = next(
        item
        for item in position_rows
        if item["campaign_id"] == "hans_interface_if_only"
        and item["canonical_family_position"] == "73"
    )
    assert row["designable"] == "False"
    assert row["position_reason"] == "not_in_campaign_design_set"

    assert "Phase 3A deterministic ProteinMPNN input preparation" in report_path.read_text(encoding="utf-8")
    assert "campaign_ss_plus_if" in config_path.read_text(encoding="utf-8")
    assert "hans_interface_ss_plus_if" in campaign_manifest_path.read_text(encoding="utf-8")
