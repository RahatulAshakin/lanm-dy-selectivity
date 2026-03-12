"""Annotate metal-centered coordination shells from local atom tables."""

from __future__ import annotations

from pathlib import Path

from lanm.configuration import load_project_config
from lanm.data.fetch import load_structure_manifest
from lanm.data.normalize import normalize_incoming_bundle
from lanm.filesystem import write_csv_rows
from lanm.logging_utils import configure_logging
from lanm.models import AtomRecord, MetalSiteSummary, ShellAnnotation
from lanm.paths import (
    LOCAL_BUNDLE_DIR,
    METAL_SITE_SUMMARY_PATH,
    SHELL_ANNOTATION_PATH,
    SITE_OVERVIEW_FIGURE_PATH,
    ensure_runtime_directories,
)
from lanm.structure.atoms import read_atom_records
from lanm.structure.geometry import annotate_metal_site
from lanm.viz.plots import render_site_overview

ATOM_TABLES = ("8fns_atoms.csv", "8dq2_atoms.csv")


def _manifest_bound_metal_map() -> dict[str, str]:
    return {entry.pdb_id: entry.bound_metal_symbol for entry in load_structure_manifest()}


def _require_atom_tables(bundle_dir: Path) -> list[Path]:
    resolved: list[Path] = []
    for name in ATOM_TABLES:
        path = bundle_dir / name
        if not path.exists():
            raise FileNotFoundError(str(path))
        resolved.append(path)
    return resolved


def _select_bound_metal_atoms(atoms: list[AtomRecord], expected_element: str) -> list[AtomRecord]:
    return sorted(
        [
            atom for atom in atoms
            if atom.record_type == "HETATM" and atom.element_upper == expected_element
        ],
        key=lambda atom: (atom.structure_id, atom.chain_id, atom.residue_seq, atom.atom_serial),
    )


def main() -> None:
    configure_logging()
    ensure_runtime_directories()
    try:
        normalize_incoming_bundle()
        config = load_project_config()
        manifest_metals = _manifest_bound_metal_map()
        summaries: list[MetalSiteSummary] = []
        annotations: list[ShellAnnotation] = []
        for atom_table_path in _require_atom_tables(LOCAL_BUNDLE_DIR):
            atoms = read_atom_records(atom_table_path)
            structure_id = atoms[0].structure_id if atoms else atom_table_path.stem.split("_", 1)[0].upper()
            expected_element = manifest_metals.get(structure_id)
            if expected_element is None:
                raise FileNotFoundError(str(atom_table_path))
            for metal_atom in _select_bound_metal_atoms(atoms, expected_element):
                summary, shell_rows = annotate_metal_site(
                    metal_atom=metal_atom,
                    atoms=atoms,
                    first_shell_cutoff_A=config.first_shell_cutoff_A,
                    second_sphere_cutoff_A=config.second_sphere_cutoff_A,
                )
                summaries.append(summary)
                annotations.extend(shell_rows)
        write_csv_rows(METAL_SITE_SUMMARY_PATH, summaries)
        write_csv_rows(SHELL_ANNOTATION_PATH, annotations)
        render_site_overview(summaries, config, SITE_OVERVIEW_FIGURE_PATH)
    except FileNotFoundError as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
