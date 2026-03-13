import csv
from pathlib import Path

from lanm.analysis.ligandmpnn_inputs import (
    build_ligandmpnn_pdb_export,
    discover_shortlisted_ligandmpnn_candidates,
    prepare_ligandmpnn_inputs,
)
from lanm.paths import (
    DESIGN_BACKBONE_MANIFEST_PATH,
    DESIGN_CAMPAIGN_MANIFEST_PATH,
    DESIGN_CAMPAIGN_POSITIONS_PATH,
    LOCAL_BUNDLE_DIR,
    PROTEINMPNN_SHORTLIST_PATH,
)
from lanm.structure.atoms import read_atom_records


def _discover_candidates():
    return discover_shortlisted_ligandmpnn_candidates(
        shortlist_path=PROTEINMPNN_SHORTLIST_PATH,
        campaign_manifest_path=DESIGN_CAMPAIGN_MANIFEST_PATH,
        backbone_manifest_path=DESIGN_BACKBONE_MANIFEST_PATH,
    )


def test_discover_shortlisted_ligandmpnn_candidates_resolves_expected_context() -> None:
    candidates = _discover_candidates()

    assert len(candidates) == 12
    assert candidates[0].candidate_id == "hans_interface_ss_plus_if_u01"
    assert candidates[0].source_structure_id == "8FNR"
    assert candidates[0].preserved_chains == ("A", "B", "C", "D")
    assert candidates[0].designed_chains == ("A",)
    assert candidates[0].fixed_context_chains == ("B", "C", "D")
    assert candidates[0].preserved_metal_identity == "DY"

    am1_candidate = next(candidate for candidate in candidates if candidate.backbone_id == "am1_mex_8fns_chain_a")
    assert am1_candidate.source_structure_id == "8FNS"
    assert am1_candidate.preserved_chains == ("A",)
    assert am1_candidate.preserved_metal_identity == "ND"
    assert am1_candidate.source_path == LOCAL_BUNDLE_DIR / "8fns_atoms.csv"

    pocket_candidate = next(candidate for candidate in candidates if candidate.backbone_id == "hans_pocket_8fnr_chain_a")
    assert pocket_candidate.source_structure_id == "8FNR"
    assert pocket_candidate.preserved_chains == ("A",)
    assert pocket_candidate.fixed_context_chains == ()
    assert pocket_candidate.preserved_metal_identity == "DY"


def test_prepare_ligandmpnn_inputs_exports_redesign_positions(tmp_path: Path) -> None:
    manifest_path = tmp_path / "ligandmpnn_input_manifest.csv"
    redesign_positions_path = tmp_path / "ligandmpnn_redesign_positions.csv"
    report_path = tmp_path / "ligandmpnn_inputs.md"
    input_root = tmp_path / "ligandmpnn"
    config_path = tmp_path / "ligandmpnn_inputs.yaml"

    manifest_rows, redesign_rows = prepare_ligandmpnn_inputs(
        shortlist_path=PROTEINMPNN_SHORTLIST_PATH,
        campaign_manifest_path=DESIGN_CAMPAIGN_MANIFEST_PATH,
        backbone_manifest_path=DESIGN_BACKBONE_MANIFEST_PATH,
        design_campaign_positions_path=DESIGN_CAMPAIGN_POSITIONS_PATH,
        manifest_path=manifest_path,
        redesign_positions_path=redesign_positions_path,
        report_path=report_path,
        input_root=input_root,
        config_path=config_path,
    )

    assert len(manifest_rows) == 12
    assert len(redesign_rows) == sum(row.redesigned_residue_count for row in manifest_rows)

    redesigned_residues_text = (
        input_root / "am1_mex_ss_only_u01" / "redesigned_residues.txt"
    ).read_text(encoding="utf-8")
    assert redesigned_residues_text == "A38 A83 A87 A111 A118\n"

    with redesign_positions_path.open("r", encoding="utf-8", newline="") as handle:
        exported_rows = list(csv.DictReader(handle))
    exported_ids = [
        row["ligandmpnn_residue_id"]
        for row in exported_rows
        if row["candidate_id"] == "am1_mex_ss_only_u01"
    ]
    assert exported_ids == ["A38", "A83", "A87", "A111", "A118"]

    mutation_tokens = [
        row["mutation_token"]
        for row in exported_rows
        if row["candidate_id"] == "hans_interface_ss_plus_if_u02"
    ]
    assert mutation_tokens == ["K12T", "D16T", "A21D", "A59Y", "R77K", "K87G", "T92D"]

    assert "hans_interface_ss_plus_if_u01" in report_path.read_text(encoding="utf-8")
    assert "am1_mex_ss_only_u01" in config_path.read_text(encoding="utf-8")


def test_build_ligandmpnn_pdb_export_preserves_metals_and_nearby_solvent() -> None:
    atoms = tuple(read_atom_records(LOCAL_BUNDLE_DIR / "8fns_atoms.csv", structure_id="8FNS"))

    export = build_ligandmpnn_pdb_export(
        atoms=atoms,
        preserved_chains=("A",),
        preserved_metal_identity="ND",
        nearby_solvent_cutoff_A=6.0,
    )

    assert export.preserved_metal_site_count == 4
    assert export.preserved_solvent_residue_count == 25

    lines = export.pdb_text.splitlines()
    hetatm_lines = [line for line in lines if line.startswith("HETATM")]

    assert lines[-1] == "END"
    assert any(line.startswith("TER") for line in lines)
    assert sum(1 for line in hetatm_lines if line[17:20].strip() == "ND") == 4
    assert sum(1 for line in hetatm_lines if line[17:20].strip() == "HOH") == 25
    assert not any(line[17:20].strip() == "HOH" and line[22:26].strip() == "301" for line in hetatm_lines)
    assert any(line[17:20].strip() == "HOH" and line[22:26].strip() == "326" for line in hetatm_lines)
    assert all(line[21].strip() == "A" for line in hetatm_lines)
