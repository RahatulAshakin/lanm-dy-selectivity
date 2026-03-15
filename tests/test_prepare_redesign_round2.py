import csv
import json
from pathlib import Path

from lanm.analysis.md_panel_selection import load_ligandmpnn_shortlist_rows
from lanm.analysis.proteinmpnn_candidates import load_design_campaign_position_rows
from lanm.analysis.redesign_round2 import (
    ROUND2_SEED_IDS,
    MDValidationPanelRecord,
    Round2Artifacts,
    Round2PositionRow,
    Round2SeedContext,
    Round2SeedRow,
    SeedResolvedPosition,
    build_redesign_round2_artifacts,
    discover_round2_seed_contexts,
    load_design_campaign_manifest_rows,
    load_md_validation_panel_rows,
    load_openmm_metadynamics_panel_status_rows,
    load_openmm_metadynamics_summary_rows,
    select_round2_redesign_positions,
    write_redesign_round2_outputs,
)
from lanm.models import DesignMaskCandidate
from lanm.paths import (
    DESIGN_CAMPAIGN_MANIFEST_PATH,
    DESIGN_CAMPAIGN_POSITIONS_PATH,
    DESIGN_MASK_CANDIDATES_PATH,
    DESIGN_MASKS_PATH,
    LIGANDMPNN_SHORTLIST_PATH,
    MD_VALIDATION_PANEL_PATH,
    OPENMM_METADYNAMICS_PANEL_STATUS_PATH,
    OPENMM_METADYNAMICS_SUMMARY_PATH,
)


def _resolved_position(
    *,
    sequence_index: int,
    residue_seq: int,
    residue_name: str,
    seed_amino_acid: str,
    canonical_family_position: int | None,
    am1_mature_position: int | None,
    fixed_first_shell: bool = False,
    protected_positions: bool = False,
    mutable_second_sphere: bool = False,
    mutable_interface: bool = False,
) -> SeedResolvedPosition:
    return SeedResolvedPosition(
        chain_id="A",
        sequence_index=sequence_index,
        residue_seq=residue_seq,
        insertion_code="",
        residue_name=residue_name,
        seed_amino_acid=seed_amino_acid,
        canonical_family_position=canonical_family_position,
        am1_mature_position=am1_mature_position,
        fixed_first_shell=fixed_first_shell,
        protected_positions=protected_positions,
        mutable_second_sphere=mutable_second_sphere,
        mutable_interface=mutable_interface,
    )


def _mask_row(
    *,
    canonical_family_position: int,
    am1_mature_position: int | None,
    mutable_second_sphere: bool = False,
    mutable_interface: bool = False,
    fixed_first_shell: bool = False,
    protected_positions: bool = False,
) -> DesignMaskCandidate:
    return DesignMaskCandidate(
        canonical_family_position=canonical_family_position,
        am1_mature_position=am1_mature_position,
        am1_reference_residue="ALA",
        observed_residue_identities="ALA",
        template_coverage_count=1,
        first_shell_observation_count=1 if fixed_first_shell else 0,
        second_sphere_observation_count=1 if mutable_second_sphere else 0,
        hans_interface_observation_count=1 if mutable_interface else 0,
        interface_neighborhood=mutable_interface,
        fixed_first_shell=fixed_first_shell,
        mutable_second_sphere=mutable_second_sphere,
        mutable_interface=mutable_interface,
        protected_positions=protected_positions,
        protection_reasons="",
        rationale="synthetic",
    )


def test_discover_round2_seed_contexts_matches_phase_7a_seed_snapshot(tmp_path: Path) -> None:
    artifacts = build_redesign_round2_artifacts(
        openmm_metadynamics_summary_path=OPENMM_METADYNAMICS_SUMMARY_PATH,
        openmm_metadynamics_panel_status_path=OPENMM_METADYNAMICS_PANEL_STATUS_PATH,
        design_mask_candidates_path=DESIGN_MASK_CANDIDATES_PATH,
        design_masks_path=DESIGN_MASKS_PATH,
        ligandmpnn_shortlist_path=LIGANDMPNN_SHORTLIST_PATH,
        md_validation_panel_path=MD_VALIDATION_PANEL_PATH,
        design_campaign_manifest_path=DESIGN_CAMPAIGN_MANIFEST_PATH,
        design_campaign_positions_path=DESIGN_CAMPAIGN_POSITIONS_PATH,
        proteinmpnn_round2_root=tmp_path / "proteinmpnn_round2",
    )

    assert [seed.candidate_id for seed in artifacts.seed_contexts] == list(ROUND2_SEED_IDS)
    assert [seed.panel_role for seed in artifacts.seed_contexts] == [
        "am1_mex_designed_candidate",
        "hans_pocket_designed_candidate",
        "hans_interface_aware_candidate",
    ]
    assert all(seed.dy_status == "retained_bound" for seed in artifacts.seed_contexts)
    assert all(seed.al_status == "persistent_capture" for seed in artifacts.seed_contexts)
    assert all(seed.fe_status == "persistent_capture" for seed in artifacts.seed_contexts)
    assert all(seed.candidate_keep_for_qm is False for seed in artifacts.seed_contexts)
    assert artifacts.seed_contexts[0].seed_mutation_canonical_positions == (18, 65, 69, 100)
    assert artifacts.seed_contexts[1].seed_mutation_canonical_positions == (18, 83, 93, 98)
    assert artifacts.seed_contexts[2].seed_mutation_canonical_positions == (18, 38, 93, 98)
    assert "Why QM Is Not Started Yet" in artifacts.report_markdown
    assert "Why These Three Seed Scaffolds" in artifacts.report_markdown
    assert "Why The Scope Is Intentionally Narrow" in artifacts.report_markdown


def test_select_round2_redesign_positions_applies_required_filters() -> None:
    seed_context = Round2SeedContext(
        seed_rank=1,
        candidate_id="seed_alpha",
        campaign_id="seed_campaign",
        backbone_id="seed_backbone",
        design_set_name="campaign_ss_plus_if",
        topology_class="am1_monomer",
        designed_chains=("A",),
        fixed_context_chains=(),
        starting_structure_path=Path("/tmp/seed_alpha.pdb"),
        designed_chain_sequence="AAAAAA",
        panel_rank=1,
        panel_role="am1_mex_designed_candidate",
        selection_reason="synthetic",
        dy_status="retained_bound",
        nd_status="retained_bound",
        y_status="retained_bound",
        al_status="persistent_capture",
        fe_status="persistent_capture",
        candidate_keep_for_qm=False,
        seed_mutation_count=1,
        seed_mutation_tokens=("K10D",),
        seed_mutation_residue_ids=("A100",),
        seed_mutation_canonical_positions=(18,),
        seed_mutation_am1_positions=(16,),
        resolved_positions=(
            _resolved_position(
                sequence_index=1,
                residue_seq=100,
                residue_name="LYS",
                seed_amino_acid="D",
                canonical_family_position=18,
                am1_mature_position=16,
                mutable_second_sphere=True,
            ),
            _resolved_position(
                sequence_index=2,
                residue_seq=101,
                residue_name="ASP",
                seed_amino_acid="D",
                canonical_family_position=19,
                am1_mature_position=17,
                mutable_second_sphere=True,
                fixed_first_shell=True,
            ),
            _resolved_position(
                sequence_index=3,
                residue_seq=102,
                residue_name="GLU",
                seed_amino_acid="E",
                canonical_family_position=25,
                am1_mature_position=23,
                mutable_interface=True,
            ),
            _resolved_position(
                sequence_index=4,
                residue_seq=103,
                residue_name="LEU",
                seed_amino_acid="L",
                canonical_family_position=40,
                am1_mature_position=36,
                mutable_second_sphere=True,
            ),
            _resolved_position(
                sequence_index=5,
                residue_seq=104,
                residue_name="ASN",
                seed_amino_acid="N",
                canonical_family_position=41,
                am1_mature_position=None,
                mutable_second_sphere=True,
            ),
            _resolved_position(
                sequence_index=6,
                residue_seq=105,
                residue_name="THR",
                seed_amino_acid="T",
                canonical_family_position=17,
                am1_mature_position=15,
                mutable_second_sphere=True,
                protected_positions=True,
            ),
        ),
    )
    design_mask_rows = (
        _mask_row(canonical_family_position=18, am1_mature_position=16, mutable_second_sphere=True),
        _mask_row(canonical_family_position=38, am1_mature_position=34, mutable_interface=True),
    )

    selected_rows = select_round2_redesign_positions(
        seed_context=seed_context,
        design_mask_rows=design_mask_rows,
        interface_window_size=3,
    )

    assert [row.canonical_family_position for row in selected_rows] == [18, 40]
    assert [row.selection_reason for row in selected_rows] == [
        "seed_mutation_window",
        "known_interface_window",
    ]
    assert selected_rows[0].is_seed_mutation_site is True
    assert selected_rows[1].near_known_interface_site is True


def test_write_redesign_round2_outputs_exports_campaign_layout(tmp_path: Path) -> None:
    source_pdb = tmp_path / "source_seed.pdb"
    source_pdb.write_text("ATOM      1  CA  ALA A   1       0.0   0.0   0.0  1.00  1.00           C\nEND\n", encoding="utf-8")
    proteinmpnn_round2_root = tmp_path / "proteinmpnn_round2"
    seed_context = Round2SeedContext(
        seed_rank=1,
        candidate_id="seed_alpha",
        campaign_id="seed_campaign",
        backbone_id="seed_backbone",
        design_set_name="campaign_ss_plus_if",
        topology_class="hans_interface_multichain",
        designed_chains=("A",),
        fixed_context_chains=("B",),
        starting_structure_path=source_pdb,
        designed_chain_sequence="ACD",
        panel_rank=4,
        panel_role="hans_interface_aware_candidate",
        selection_reason="synthetic seed",
        dy_status="retained_bound",
        nd_status="retained_bound",
        y_status="retained_bound",
        al_status="persistent_capture",
        fe_status="persistent_capture",
        candidate_keep_for_qm=False,
        seed_mutation_count=1,
        seed_mutation_tokens=("A2C",),
        seed_mutation_residue_ids=("A2",),
        seed_mutation_canonical_positions=(40,),
        seed_mutation_am1_positions=(36,),
        resolved_positions=(
            _resolved_position(
                sequence_index=1,
                residue_seq=1,
                residue_name="ALA",
                seed_amino_acid="A",
                canonical_family_position=38,
                am1_mature_position=34,
                mutable_interface=True,
            ),
            _resolved_position(
                sequence_index=2,
                residue_seq=2,
                residue_name="ALA",
                seed_amino_acid="C",
                canonical_family_position=40,
                am1_mature_position=36,
                mutable_second_sphere=True,
            ),
            _resolved_position(
                sequence_index=3,
                residue_seq=3,
                residue_name="ASP",
                seed_amino_acid="D",
                canonical_family_position=44,
                am1_mature_position=40,
                mutable_second_sphere=True,
            ),
        ),
    )
    position_row = Round2PositionRow(
        seed_rank=1,
        candidate_id="seed_alpha",
        campaign_id="seed_campaign",
        backbone_id="seed_backbone",
        topology_class="hans_interface_multichain",
        chain_id="A",
        sequence_index=2,
        residue_seq=2,
        insertion_code="",
        residue_id="A2",
        residue_name="ALA",
        native_amino_acid="A",
        seed_amino_acid="C",
        canonical_family_position=40,
        am1_mature_position=36,
        mutable_second_sphere=True,
        mutable_interface=False,
        near_seed_mutation=True,
        near_known_interface_site=True,
        is_seed_mutation_site=True,
        selection_reason="seed_mutation_window,known_interface_window",
    )
    seed_row = Round2SeedRow(
        seed_rank=1,
        candidate_id="seed_alpha",
        campaign_id="seed_campaign",
        backbone_id="seed_backbone",
        topology_class="hans_interface_multichain",
        panel_rank=4,
        panel_role="hans_interface_aware_candidate",
        designed_chains="A",
        fixed_context_chains="B",
        dy_status="retained_bound",
        nd_status="retained_bound",
        y_status="retained_bound",
        al_status="persistent_capture",
        fe_status="persistent_capture",
        candidate_keep_for_qm=False,
        starting_structure_path=str(source_pdb),
        exported_backbone_path=str(proteinmpnn_round2_root / "backbones" / "seed_alpha.pdb"),
        exported_pdb_path=str(proteinmpnn_round2_root / "campaigns" / "seed_alpha" / "seed_alpha.pdb"),
        chain_assignment_path=str(proteinmpnn_round2_root / "campaigns" / "seed_alpha" / "chain_id.jsonl"),
        fixed_positions_path=str(proteinmpnn_round2_root / "campaigns" / "seed_alpha" / "fixed_positions.jsonl"),
        seed_mutation_count=1,
        seed_mutation_tokens="A2C",
        seed_mutation_residue_ids="A2",
        seed_mutation_canonical_positions="40",
        seed_mutation_am1_positions="36",
        redesignable_position_count=1,
        redesignable_residue_ids="A2",
        redesignable_canonical_positions="40",
        redesignable_am1_positions="36",
        selection_reason="synthetic seed",
    )
    artifacts = Round2Artifacts(
        interface_window_size=3,
        seed_contexts=(seed_context,),
        seed_rows=(seed_row,),
        position_rows=(position_row,),
        config_payload={"version": 1, "phase": "7A"},
        report_markdown="# synthetic\n",
    )
    config_path = tmp_path / "redesign_round2.yaml"
    seed_table_path = tmp_path / "redesign_round2_seeds.csv"
    position_table_path = tmp_path / "redesign_round2_positions.csv"
    report_path = tmp_path / "redesign_round2_plan.md"

    write_redesign_round2_outputs(
        artifacts=artifacts,
        config_path=config_path,
        seed_table_path=seed_table_path,
        position_table_path=position_table_path,
        report_path=report_path,
        proteinmpnn_round2_root=proteinmpnn_round2_root,
    )

    assert config_path.exists()
    assert seed_table_path.exists()
    assert position_table_path.exists()
    assert report_path.exists()
    exported_backbone = proteinmpnn_round2_root / "backbones" / "seed_alpha.pdb"
    exported_pdb = proteinmpnn_round2_root / "campaigns" / "seed_alpha" / "seed_alpha.pdb"
    chain_assignment_path = proteinmpnn_round2_root / "campaigns" / "seed_alpha" / "chain_id.jsonl"
    fixed_positions_path = proteinmpnn_round2_root / "campaigns" / "seed_alpha" / "fixed_positions.jsonl"

    assert exported_backbone.read_text(encoding="utf-8") == source_pdb.read_text(encoding="utf-8")
    assert exported_pdb.read_text(encoding="utf-8") == source_pdb.read_text(encoding="utf-8")
    assert json.loads(chain_assignment_path.read_text(encoding="utf-8")) == {
        "seed_alpha": [["A"], ["B"]]
    }
    assert json.loads(fixed_positions_path.read_text(encoding="utf-8")) == {
        "seed_alpha": {"A": [1, 3]}
    }

    with seed_table_path.open("r", encoding="utf-8", newline="") as handle:
        seed_rows = list(csv.DictReader(handle))
    with position_table_path.open("r", encoding="utf-8", newline="") as handle:
        position_rows = list(csv.DictReader(handle))
    assert seed_rows[0]["candidate_id"] == "seed_alpha"
    assert position_rows[0]["residue_id"] == "A2"
