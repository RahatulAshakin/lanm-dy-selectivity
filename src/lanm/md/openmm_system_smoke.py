"""OpenMM helpers for Phase 6A3 system-build smoke tests."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from io import StringIO
from typing import Iterator

from lanm.md.openmm_forcefields import resolve_openmm_forcefield_files


@dataclass(frozen=True, slots=True)
class OpenMMSerializedSystem:
    prepared_structure_pdb_text: str
    system_xml: str
    integrator_xml: str
    atom_count: int
    residue_count: int


@contextmanager
def _patch_forcefield_create_system_for_truncated_chains(forcefield: object) -> Iterator[object]:
    """Temporarily force OpenMM to ignore external bonds for clipped panel chains."""
    try:
        from openmm import app, unit
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for Phase 6A3 system-build smoke tests") from exc

    original_create_system = forcefield.createSystem

    def patched_create_system(
        topology: object,
        nonbondedMethod: object = app.NoCutoff,
        nonbondedCutoff: object = 1.0 * unit.nanometer,
        constraints: object | None = None,
        rigidWater: object | None = None,
        removeCMMotion: bool = True,
        hydrogenMass: object | None = None,
        residueTemplates: dict[object, object] = {},
        ignoreExternalBonds: bool = False,
        switchDistance: object | None = None,
        flexibleConstraints: bool = False,
        drudeMass: object = 0.4 * unit.daltons,
        **args: object,
    ) -> object:
        return original_create_system(
            topology,
            nonbondedMethod=nonbondedMethod,
            nonbondedCutoff=nonbondedCutoff,
            constraints=constraints,
            rigidWater=rigidWater,
            removeCMMotion=removeCMMotion,
            hydrogenMass=hydrogenMass,
            residueTemplates=residueTemplates,
            ignoreExternalBonds=True,
            switchDistance=switchDistance,
            flexibleConstraints=flexibleConstraints,
            drudeMass=drudeMass,
            **args,
        )

    forcefield.createSystem = patched_create_system  # type: ignore[attr-defined]
    try:
        yield original_create_system
    finally:
        forcefield.createSystem = original_create_system  # type: ignore[attr-defined]


def build_openmm_serialized_system(
    *,
    prepared_structure_pdb_text: str,
    forcefield_files: tuple[str, ...],
    target_temperature_K: int,
) -> OpenMMSerializedSystem:
    """Hydrogenate a prepared structure and serialize the OpenMM System and integrator."""
    try:
        from openmm import LangevinMiddleIntegrator, XmlSerializer, app, unit
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for Phase 6A3 system-build smoke tests") from exc

    pdb = app.PDBFile(StringIO(prepared_structure_pdb_text))
    resolved_forcefield_files = resolve_openmm_forcefield_files(forcefield_files)
    forcefield = app.ForceField(*resolved_forcefield_files)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    with _patch_forcefield_create_system_for_truncated_chains(forcefield) as original_create_system:
        modeller.addHydrogens(forcefield)
        system = original_create_system(
            modeller.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=app.HBonds,
            ignoreExternalBonds=True,
        )
    integrator = LangevinMiddleIntegrator(
        target_temperature_K * unit.kelvin,
        1.0 / unit.picosecond,
        0.004 * unit.picoseconds,
    )
    prepared_structure_handle = StringIO()
    app.PDBFile.writeFile(
        modeller.topology,
        modeller.positions,
        prepared_structure_handle,
        keepIds=True,
    )
    return OpenMMSerializedSystem(
        prepared_structure_pdb_text=prepared_structure_handle.getvalue(),
        system_xml=XmlSerializer.serialize(system),
        integrator_xml=XmlSerializer.serialize(integrator),
        atom_count=modeller.topology.getNumAtoms(),
        residue_count=modeller.topology.getNumResidues(),
    )
