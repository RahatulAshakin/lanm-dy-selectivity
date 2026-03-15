from __future__ import annotations

from pathlib import Path

from lanm.md import openmm_forcefields
from lanm.md.openmm_forcefields import (
    AMBER19_ALL_COMPONENTS,
    AMBER19_ALL_XML,
    resolve_openmm_forcefield_files,
)


def test_resolve_openmm_forcefield_files_resolves_bundled_resource_names() -> None:
    resolved = resolve_openmm_forcefield_files((AMBER19_ALL_XML, "amber19/opc3.xml"))

    assert all(Path(path).is_absolute() for path in resolved)
    assert all(Path(path).is_file() for path in resolved)
    assert Path(resolved[-1]).as_posix().endswith("amber19/opc3.xml")
    protein_forcefield_paths = resolved[:-1]
    if len(protein_forcefield_paths) == 1:
        assert Path(protein_forcefield_paths[0]).name == AMBER19_ALL_XML
    else:
        assert len(protein_forcefield_paths) == len(AMBER19_ALL_COMPONENTS)
        assert [
            Path(path).as_posix().endswith(component)
            for path, component in zip(protein_forcefield_paths, AMBER19_ALL_COMPONENTS, strict=True)
        ] == [True] * len(AMBER19_ALL_COMPONENTS)


def test_resolve_openmm_forcefield_files_preserves_absolute_paths() -> None:
    bundled_paths = resolve_openmm_forcefield_files((AMBER19_ALL_XML, "amber19/opc3.xml"))

    resolved = resolve_openmm_forcefield_files(bundled_paths)

    assert resolved == bundled_paths


def test_resolve_openmm_forcefield_files_expands_amber19_all_shortcut_when_missing(
    tmp_path: Path,
    monkeypatch,
) -> None:
    data_directory = tmp_path / "openmm-data"
    for relative_path in (*AMBER19_ALL_COMPONENTS, "amber19/opc3.xml"):
        candidate = data_directory / relative_path
        candidate.parent.mkdir(parents=True, exist_ok=True)
        candidate.write_text("<ForceField/>", encoding="utf-8")

    monkeypatch.setattr(
        openmm_forcefields,
        "_iter_openmm_forcefield_data_directories",
        lambda: (data_directory,),
    )

    resolved = resolve_openmm_forcefield_files((AMBER19_ALL_XML, "amber19/opc3.xml"))

    assert resolved == tuple(
        str((data_directory / relative_path).resolve())
        for relative_path in (*AMBER19_ALL_COMPONENTS, "amber19/opc3.xml")
    )
