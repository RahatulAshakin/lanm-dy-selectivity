"""Helpers for resolving OpenMM force-field XML resources."""

from __future__ import annotations

from pathlib import Path

AMBER19_ALL_XML = "amber19-all.xml"
AMBER19_ALL_COMPONENTS = (
    "amber19/protein.ff19SB.xml",
    "amber19/DNA.OL21.xml",
    "amber14/RNA.OL3.xml",
    "amber19/lipid21.xml",
)


def _iter_openmm_forcefield_data_directories() -> tuple[Path, ...]:
    try:
        from openmm.app import forcefield as forcefield_module
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for force-field resolution") from exc

    get_data_directories = getattr(forcefield_module, "_getDataDirectories", None)
    if callable(get_data_directories):
        return tuple(Path(path) for path in get_data_directories())
    return (Path(forcefield_module.__file__).resolve().parent / "data",)


def _resolve_bundled_forcefield_file(
    file_name: str | Path,
    *,
    data_directories: tuple[Path, ...],
) -> Path | None:
    candidate = Path(file_name)
    for data_directory in data_directories:
        bundled_candidate = data_directory / candidate
        if bundled_candidate.is_file():
            return bundled_candidate.resolve()
    return None


def _resolve_required_bundled_forcefield_file(
    file_name: str,
    *,
    data_directories: tuple[Path, ...],
) -> Path:
    resolved_path = _resolve_bundled_forcefield_file(file_name, data_directories=data_directories)
    if resolved_path is not None:
        return resolved_path
    data_directory_text = ", ".join(str(path) for path in data_directories)
    raise FileNotFoundError(
        f'Could not resolve OpenMM bundled force-field XML "{file_name}" from {data_directory_text}'
    )


def resolve_openmm_forcefield_files(forcefield_files: tuple[str, ...]) -> tuple[str, ...]:
    """Resolve built-in OpenMM force-field resource names to filesystem paths."""
    data_directories = _iter_openmm_forcefield_data_directories()
    resolved_files: list[str] = []
    for file_name in forcefield_files:
        candidate = Path(file_name).expanduser()
        if candidate.is_file():
            resolved_files.append(str(candidate.resolve()))
            continue
        bundled_candidate = _resolve_bundled_forcefield_file(candidate, data_directories=data_directories)
        if bundled_candidate is not None:
            resolved_files.append(str(bundled_candidate))
            continue
        if file_name == AMBER19_ALL_XML:
            resolved_files.extend(
                str(_resolve_required_bundled_forcefield_file(component, data_directories=data_directories))
                for component in AMBER19_ALL_COMPONENTS
            )
            continue
        raise FileNotFoundError(
            f'Could not resolve OpenMM force-field XML "{file_name}" as a filesystem path or bundled resource'
        )
    return tuple(resolved_files)
