import csv
from pathlib import Path

import yaml

from lanm.analysis.md_input_preparation import (
    TARGET_METALS,
    discover_md_validation_panel,
    load_design_backbone_sources,
    prepare_md_inputs,
    render_target_metal_starting_structure,
    resolve_panel_member_structure_template,
)
from lanm.paths import DESIGN_BACKBONE_MANIFEST_PATH, MD_PANEL_CONFIG_PATH, MD_VALIDATION_PANEL_PATH


def _write_text(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _hetatm_lines(pdb_text: str) -> list[str]:
    return [line for line in pdb_text.splitlines() if line.startswith("HETATM")]


def test_discover_md_validation_panel_matches_phase_5b_snapshot() -> None:
    panel_rows = discover_md_validation_panel(
        panel_path=MD_VALIDATION_PANEL_PATH,
        md_panel_config_path=MD_PANEL_CONFIG_PATH,
    )

    assert len(panel_rows) == 6
    assert [row.panel_rank for row in panel_rows] == [1, 2, 3, 4, 5, 6]
    assert [row.panel_member_id for row in panel_rows] == [
        "am1_mex_ss_only_u02",
        "am1_mex_wt_reference",
        "hans_pocket_ss_only_u02",
        "hans_interface_ss_plus_if_u04",
        "hans_interface_wt_reference",
        "hans_interface_ss_plus_if_u02",
    ]


def test_render_target_metal_starting_structure_preserves_multichain_numbering_and_insertion_codes(
    tmp_path: Path,
) -> None:
    designed_pdb_path = _write_text(
        tmp_path / "inputs" / "designed_candidate.pdb",
        "\n".join(
            [
                "ATOM      1  N   GLY A  10A      1.000   2.000   3.000  1.00 20.00           N  ",
                "ATOM      2  CA  GLY A  10A      1.500   2.500   3.500  1.00 20.00           C  ",
                "ATOM      3  N   ALA B  20       4.000   5.000   6.000  1.00 20.00           N  ",
                "ATOM      4  CA  ALA B  20       4.500   5.500   6.500  1.00 20.00           C  ",
                "HETATM    5  DY  DY  A 201       1.200   2.200   3.200  1.00  0.00          DY  ",
                "HETATM    6  O   HOH A 301       1.350   2.350   3.350  1.00 10.00           O  ",
                "END",
                "",
            ]
        ),
    )
    panel_csv = _write_text(
        tmp_path / "md_validation_panel.csv",
        "\n".join(
            [
                "panel_rank,panel_member_id,panel_member_type,panel_role,candidate_id,campaign_id,backbone_id,topology_class,preserved_metal_identity,proteinmpnn_rank,ligandmpnn_representative_ligand_confidence,ligandmpnn_representative_overall_confidence,ligandmpnn_mean_ligand_confidence,rosetta_score_rank_within_topology,rosetta_total_score,starting_structure_path,selection_reason",
                f"1,synthetic_interface,designed_candidate,hans_interface_aware_candidate,synthetic_interface,synthetic_campaign,synthetic_backbone,hans_interface_multichain,DY,1,0.5,0.6,0.4,1,10.0,{designed_pdb_path},synthetic test row",
                "",
            ]
        ),
    )
    md_panel_yaml_path = tmp_path / "md_panel.yaml"
    md_panel_yaml_path.write_text(
        yaml.safe_dump(
            {
                "panel": [
                    {
                        "panel_rank": 1,
                        "panel_member_id": "synthetic_interface",
                        "panel_member_type": "designed_candidate",
                        "panel_role": "hans_interface_aware_candidate",
                        "candidate_id": "synthetic_interface",
                        "backbone_id": "synthetic_backbone",
                        "topology_class": "hans_interface_multichain",
                        "starting_structure_path": str(designed_pdb_path),
                    }
                ]
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    backbone_manifest = _write_text(
        tmp_path / "design_backbone_manifest.csv",
        "\n".join(
            [
                "backbone_id,structure_id,source_kind,source_path,experimental_method,selected_chains,design_chains,fixed_context_chains,selected_residue_count,design_chain_residue_count,design_chain_metal_site_count,representative_note,backbone_pdb_path",
                f"synthetic_backbone,SYN,csv_atom_table,{tmp_path / 'unused.csv'},X-RAY DIFFRACTION,\"A,B\",A,B,2,1,1,synthetic,{designed_pdb_path}",
                "",
            ]
        ),
    )

    panel_member = discover_md_validation_panel(
        panel_path=panel_csv,
        md_panel_config_path=md_panel_yaml_path,
    )[0]
    template = resolve_panel_member_structure_template(
        panel_member=panel_member,
        design_backbone_sources=load_design_backbone_sources(backbone_manifest),
        nearby_solvent_cutoff_A=6.0,
    )
    rendered = render_target_metal_starting_structure(template, target_metal="Fe")

    assert rendered.metal_site_count == 1
    assert rendered.preserved_solvent_residue_count == 1
    hetatm_lines = _hetatm_lines(rendered.pdb_text)
    assert len(hetatm_lines) == 2
    assert any(
        line[12:16].strip() == "FE"
        and line[17:20].strip() == "FE"
        and line[21].strip() == "A"
        and line[22:26].strip() == "201"
        for line in hetatm_lines
    )
    assert not any(line[12:16].strip() == "DY" or line[17:20].strip() == "DY" for line in hetatm_lines)
    assert any(
        line[17:20].strip() == "HOH"
        and line[21].strip() == "A"
        and line[22:26].strip() == "301"
        for line in hetatm_lines
    )
    assert "GLY A  10A" in rendered.pdb_text
    assert "ALA B  20 " in rendered.pdb_text
    atom_chain_ids = {
        line[21].strip()
        for line in rendered.pdb_text.splitlines()
        if line.startswith("ATOM")
    }
    assert atom_chain_ids == {"A", "B"}


def test_prepare_md_inputs_writes_manifest_metadata_and_protocol(tmp_path: Path) -> None:
    manifest_path = tmp_path / "results" / "tables" / "md_input_manifest.csv"
    report_path = tmp_path / "results" / "reports" / "md_input_preparation.md"
    protocol_path = tmp_path / "config" / "md_protocol.yaml"
    md_inputs_root = tmp_path / "results" / "md_inputs"

    manifest_rows = prepare_md_inputs(
        panel_path=MD_VALIDATION_PANEL_PATH,
        md_panel_config_path=MD_PANEL_CONFIG_PATH,
        design_backbone_manifest_path=DESIGN_BACKBONE_MANIFEST_PATH,
        manifest_path=manifest_path,
        report_path=report_path,
        protocol_config_path=protocol_path,
        md_inputs_root=md_inputs_root,
    )

    assert len(manifest_rows) == 30
    assert manifest_path.exists()
    assert report_path.exists()
    assert protocol_path.exists()
    assert {row.target_metal for row in manifest_rows} == set(TARGET_METALS)

    with manifest_path.open("r", encoding="utf-8", newline="") as handle:
        exported_rows = list(csv.DictReader(handle))
    assert len(exported_rows) == 30

    dy_row = next(
        row
        for row in exported_rows
        if row["panel_member_id"] == "am1_mex_ss_only_u02"
        and row["target_metal"] == "Dy"
    )
    metadata_path = Path(dy_row["system_metadata_path"])
    starting_structure_path = Path(dy_row["starting_structure_path"])

    assert metadata_path.exists()
    assert starting_structure_path.exists()

    metadata = yaml.safe_load(metadata_path.read_text(encoding="utf-8"))
    assert metadata["panel_member_id"] == "am1_mex_ss_only_u02"
    assert metadata["candidate_id"] == "am1_mex_ss_only_u02"
    assert metadata["target_metal"] == "Dy"
    assert metadata["replicate_count"] == 3
    assert metadata["target_temperature_K"] == 298
    assert metadata["suggested_collective_variables"] == [
        "coordination_number",
        "mean_metal_oxygen_distance",
    ]

    protocol_payload = yaml.safe_load(protocol_path.read_text(encoding="utf-8"))
    assert protocol_payload["phase"] == "6A1"
    assert protocol_payload["target_metals"] == list(TARGET_METALS)
    assert len(protocol_payload["panel_members"]) == 6
    assert "OpenMM-ready" in report_path.read_text(encoding="utf-8")


def test_resolve_panel_member_structure_template_restores_wild_type_metal_context() -> None:
    panel_rows = discover_md_validation_panel(
        panel_path=MD_VALIDATION_PANEL_PATH,
        md_panel_config_path=MD_PANEL_CONFIG_PATH,
    )
    wild_type_row = next(row for row in panel_rows if row.panel_member_id == "am1_mex_wt_reference")

    template = resolve_panel_member_structure_template(
        panel_member=wild_type_row,
        design_backbone_sources=load_design_backbone_sources(DESIGN_BACKBONE_MANIFEST_PATH),
        nearby_solvent_cutoff_A=6.0,
    )
    rendered = render_target_metal_starting_structure(template, target_metal="Dy")

    assert template.source_structure == "results/design_inputs/proteinmpnn/backbones/am1_mex_8fns_chain_a.pdb"
    assert template.metal_site_template_source == "data/raw/local_bundle/8fns_atoms.csv"
    assert len(template.metal_site_atoms) == 4
    assert len({atom.residue_key for atom in template.nearby_solvent_atoms}) == 25
    assert rendered.metal_site_count == 4
    assert rendered.preserved_solvent_residue_count == 25
    hetatm_lines = _hetatm_lines(rendered.pdb_text)
    assert len(hetatm_lines) == 29
    assert sum(1 for line in hetatm_lines if line[17:20].strip() == "DY") == 4
    assert any(
        line[12:16].strip() == "DY"
        and line[17:20].strip() == "DY"
        and line[21].strip() == "A"
        and line[22:26].strip() == "201"
        for line in hetatm_lines
    )
