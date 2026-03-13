from pathlib import Path

import pytest

from lanm.analysis.template_harmonization import (
    build_template_harmonization_artifacts,
    identify_am1_mature_sequence_reference,
    load_sequence_rows,
    write_template_harmonization_outputs,
)
from lanm.paths import LOCAL_STRUCTURE_MANIFEST_PATH, REPO_ROOT
from lanm.structure.templates import parse_cif_atom_records, read_cif_experimental_method


def test_parse_cif_atom_records_uses_first_model_only_for_6mi5() -> None:
    path = REPO_ROOT / "data" / "raw" / "public" / "structures" / "6MI5.cif"

    atoms = parse_cif_atom_records(path)
    metal_atoms = [atom for atom in atoms if atom.record_type == "HETATM" and atom.element == "Y"]
    polymer_residues = sorted({atom.residue_seq for atom in atoms if atom.record_type == "ATOM"})

    assert read_cif_experimental_method(path) == "SOLUTION NMR"
    assert len(metal_atoms) == 3
    assert {atom.chain_id for atom in metal_atoms} == {"X"}
    assert polymer_residues[0] == 23
    assert polymer_residues[-1] == 139


def test_build_template_harmonization_artifacts_matches_expected_template_counts() -> None:
    artifacts = build_template_harmonization_artifacts(
        sequences_path=REPO_ROOT / "data" / "raw" / "local_bundle" / "lanmodulin_sequences.csv",
        manifest_path=LOCAL_STRUCTURE_MANIFEST_PATH,
    )
    summary_map = {summary.template_id: summary for summary in artifacts.template_summaries}

    assert artifacts.sequence_record_count == 3
    assert [summary.template_id for summary in artifacts.template_summaries] == [
        "6MI5",
        "8FNS",
        "8DQ2",
        "8FNR",
    ]
    assert len(artifacts.chain_rows) == 10
    assert len(artifacts.site_rows) == 33
    assert len(artifacts.cross_template_rows) == 1097
    assert any(row.role == "interchain_contact" for row in artifacts.residue_role_rows)
    assert any(row.role == "solvent_contact" for row in artifacts.residue_role_rows)
    assert summary_map["6MI5"].source_type == "cif"
    assert summary_map["6MI5"].chain_ids == ("X",)
    assert summary_map["6MI5"].metal_site_count == 3
    assert summary_map["8FNS"].residue_ranges == ("A:29-133",)
    assert summary_map["8DQ2"].metal_site_count == 12
    assert summary_map["8FNR"].experimental_method == "X-RAY DIFFRACTION"

    chain_d_row = next(
        row
        for row in artifacts.chain_rows
        if row.template_id == "8FNR" and row.chain_id == "D"
    )
    assert chain_d_row.residue_count == 105
    assert chain_d_row.metal_site_count == 4

    am1_start = next(
        row for row in artifacts.cross_template_rows
        if row.template_id == "6MI5" and row.chain_id == "X" and row.template_residue_seq == 23
    )
    assert am1_start.canonical_family_position == 1
    assert am1_start.am1_mature_position == 1
    assert am1_start.alignment_status == "aligned_match"

    am1_his_tag = next(
        row for row in artifacts.cross_template_rows
        if row.template_id == "6MI5" and row.chain_id == "X" and row.template_residue_seq == 139
    )
    assert am1_his_tag.canonical_family_position is None
    assert am1_his_tag.am1_mature_position is None
    assert am1_his_tag.alignment_status == "outside_am1_mature_reference"

    fns_start = next(
        row for row in artifacts.cross_template_rows
        if row.template_id == "8FNS" and row.chain_id == "A" and row.template_residue_seq == 29
    )
    assert fns_start.am1_mature_position == 7
    assert fns_start.alignment_status == "aligned_match"

    first_shell_row = next(
        row for row in artifacts.residue_role_rows
        if row.template_id == "8FNS"
        and row.chain_id == "A"
        and row.site_index == 1
        and row.role == "first_shell"
        and row.residue_seq == 35
    )
    assert first_shell_row.site_label == "A:201"
    assert first_shell_row.canonical_family_position is not None
    assert first_shell_row.am1_mature_position == 13
    assert first_shell_row.alignment_status == "aligned_match"

    solvent_row = next(
        row for row in artifacts.residue_role_rows
        if row.template_id == "8FNS"
        and row.chain_id == "A"
        and row.site_index == 1
        and row.role == "solvent_contact"
        and row.residue_name == "HOH"
        and row.residue_seq == 326
    )
    assert solvent_row.alignment_status == "non_polymer"

    interchain_row = next(
        row for row in artifacts.residue_role_rows
        if row.template_id == "8DQ2"
        and row.role == "interchain_contact"
    )
    assert interchain_row.site_index is None
    assert interchain_row.site_label == ""
    assert interchain_row.note.startswith("partner_chains=")

    assert "## AM1/Mex family" in artifacts.report_markdown
    assert "## Hans family" in artifacts.report_markdown
    assert "## Alignment and residue roles" in artifacts.report_markdown
    assert "`8DQ2` and `8FNR` both resolve chains A, B, C, D" in artifacts.report_markdown


def test_write_template_harmonization_outputs_creates_five_files(tmp_path: Path) -> None:
    artifacts = build_template_harmonization_artifacts(
        sequences_path=REPO_ROOT / "data" / "raw" / "local_bundle" / "lanmodulin_sequences.csv",
        manifest_path=LOCAL_STRUCTURE_MANIFEST_PATH,
    )
    report_path = tmp_path / "template_harmonization.md"
    chain_summary_path = tmp_path / "template_chain_summary.csv"
    site_summary_path = tmp_path / "template_site_summary.csv"
    cross_template_alignment_path = tmp_path / "cross_template_residue_alignment.csv"
    residue_role_map_path = tmp_path / "residue_role_map.csv"

    write_template_harmonization_outputs(
        artifacts=artifacts,
        report_path=report_path,
        chain_summary_path=chain_summary_path,
        site_summary_path=site_summary_path,
        cross_template_alignment_path=cross_template_alignment_path,
        residue_role_map_path=residue_role_map_path,
    )

    assert report_path.exists()
    assert chain_summary_path.exists()
    assert site_summary_path.exists()
    assert cross_template_alignment_path.exists()
    assert residue_role_map_path.exists()
    assert "template_id,source_type,chain_id" in chain_summary_path.read_text(encoding="utf-8")
    assert "template_id,chain_id,template_residue_seq" in cross_template_alignment_path.read_text(encoding="utf-8")
    assert "template_id,chain_id,site_index,site_label,role" in residue_role_map_path.read_text(encoding="utf-8")


def test_identify_am1_mature_sequence_reference_returns_confirmed_reference() -> None:
    am1_reference = identify_am1_mature_sequence_reference(
        load_sequence_rows(REPO_ROOT / "data" / "raw" / "local_bundle" / "lanmodulin_sequences.csv")
    )

    assert len(am1_reference.sequence) == 111
    assert am1_reference.sequence.startswith("PTTTTKVDIAAF")
    assert "23-133" in am1_reference.note


def test_load_sequence_rows_raises_exact_missing_path(tmp_path: Path) -> None:
    missing_path = tmp_path / "missing_sequences.csv"

    with pytest.raises(FileNotFoundError, match=f"^{missing_path}$"):
        load_sequence_rows(missing_path)
