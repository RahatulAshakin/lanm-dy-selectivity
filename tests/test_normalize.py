from pathlib import Path

from lanm.data.normalize import normalize_incoming_bundle, should_skip_normalization, snake_case_filename


def test_snake_case_filename_is_deterministic() -> None:
    assert snake_case_filename("Project description .docx") == "project_description.docx"
    assert snake_case_filename("8FNS_atoms.csv") == "8fns_atoms.csv"
    assert snake_case_filename("lanthanide repository-download links.TXT") == "lanthanide_repository_download_links.txt"


def test_should_skip_zone_identifier_sidecars() -> None:
    assert should_skip_normalization(Path("foo.csv:Zone.Identifier"))
    assert not should_skip_normalization(Path("foo.csv"))


def test_normalize_incoming_bundle_skips_zone_identifier(tmp_path: Path) -> None:
    incoming_dir = tmp_path / "incoming"
    bundle_dir = tmp_path / "bundle"
    incoming_dir.mkdir()
    bundle_dir.mkdir()
    (incoming_dir / "useful.txt").write_text("ok\n", encoding="utf-8")
    (incoming_dir / "useful.txt:Zone.Identifier").write_text("ignored\n", encoding="utf-8")
    (bundle_dir / "useful_txt_zone.identifier").write_text("stale\n", encoding="utf-8")

    result = normalize_incoming_bundle(incoming_dir=incoming_dir, local_bundle_dir=bundle_dir)

    assert (bundle_dir / "useful.txt").exists()
    assert not (bundle_dir / "useful_txt_zone.identifier").exists()
    assert not (bundle_dir / "useful.txt_zone_identifier").exists()
    assert len(result.normalized_paths) == 1
    assert len(result.skipped_sidecars) == 1
