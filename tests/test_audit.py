from pathlib import Path

from lanm.data.audit import summarize_local_bundle


def test_summarize_local_bundle_handles_docx_like_binary_without_crashing(tmp_path: Path) -> None:
    bundle_dir = tmp_path / "bundle"
    bundle_dir.mkdir()
    docx_path = bundle_dir / "project_description.docx"
    docx_path.write_bytes(b"PK\x03\x04\x14\x00\x00\x00binary-docx")

    summaries = summarize_local_bundle(bundle_dir)

    assert len(summaries) == 1
    assert summaries[0].file_type == "docx"
    assert summaries[0].row_count == 0
    assert summaries[0].file_size_bytes == docx_path.stat().st_size
