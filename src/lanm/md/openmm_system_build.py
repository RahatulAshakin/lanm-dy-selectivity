"""Configurable OpenMM serialization helpers for Phase 6A4 system builds."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from io import StringIO
from typing import Iterator


DEFAULT_NONBONDED_METHOD = "NoCutoff"
DEFAULT_CUTOFF_NM = 1.0
DEFAULT_FRICTION_COEFF_PS = 1.0
DEFAULT_TIMESTEP_FS = 2.0
DEFAULT_HYDROGEN_MASS_REPARTITIONING = False
DEFAULT_HYDROGEN_MASS_DA = 3.024
SUPPORTED_NONBONDED_METHODS = (
    "NoCutoff",
    "CutoffNonPeriodic",
    "CutoffPeriodic",
    "PME",
)
EXPLICIT_SINGLE_ATOM_RESIDUE_TEMPLATES = {
    "FE": "FE",
    "FE2": "FE2",
}


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
        raise RuntimeError("OpenMM is required for Phase 6A4 system builds") from exc

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


def _normalize_nonbonded_method(method_name: str) -> str:
    normalized = method_name.strip()
    if normalized not in SUPPORTED_NONBONDED_METHODS:
        raise ValueError(
            f"Unsupported nonbonded_method {method_name!r}; expected one of {SUPPORTED_NONBONDED_METHODS}"
        )
    return normalized


def _sanitize_prepared_structure_pdb_text(pdb_text: str) -> str:
    """Strip OpenMM's dated REMARK header so the artifact remains deterministic."""
    lines = [
        line
        for line in pdb_text.splitlines()
        if not line.startswith("REMARK   1 CREATED WITH OPENMM")
    ]
    return "\n".join(lines) + ("\n" if lines else "")


def _build_explicit_residue_templates(topology: object) -> dict[object, str]:
    """Resolve deterministic residue-template choices for ambiguous single-atom ions."""
    residue_templates: dict[object, str] = {}
    for residue in topology.residues():
        atoms = tuple(residue.atoms())
        if len(atoms) != 1:
            continue
        template_name = EXPLICIT_SINGLE_ATOM_RESIDUE_TEMPLATES.get(str(residue.name).strip())
        if template_name is None:
            continue
        residue_templates[residue] = template_name
    return residue_templates


def build_openmm_serialized_system(
    *,
    prepared_structure_pdb_text: str,
    forcefield_files: tuple[str, ...],
    target_temperature_K: int,
    friction_coeff_ps: float = DEFAULT_FRICTION_COEFF_PS,
    timestep_fs: float = DEFAULT_TIMESTEP_FS,
    nonbonded_method: str = DEFAULT_NONBONDED_METHOD,
    cutoff_nm: float = DEFAULT_CUTOFF_NM,
    hydrogen_mass_repartitioning: bool = DEFAULT_HYDROGEN_MASS_REPARTITIONING,
    hydrogen_mass_da: float = DEFAULT_HYDROGEN_MASS_DA,
) -> OpenMMSerializedSystem:
    """Hydrogenate a prepared structure and serialize the OpenMM System and integrator."""
    try:
        from openmm import LangevinMiddleIntegrator, XmlSerializer, app, unit
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for Phase 6A4 system builds") from exc

    normalized_nonbonded_method = _normalize_nonbonded_method(nonbonded_method)
    if target_temperature_K <= 0:
        raise ValueError("target_temperature_K must be positive")
    if friction_coeff_ps <= 0.0:
        raise ValueError("friction_coeff_ps must be positive")
    if timestep_fs <= 0.0:
        raise ValueError("timestep_fs must be positive")
    if cutoff_nm <= 0.0:
        raise ValueError("cutoff_nm must be positive")
    if hydrogen_mass_repartitioning and hydrogen_mass_da <= 0.0:
        raise ValueError("hydrogen_mass_da must be positive when hydrogen_mass_repartitioning is enabled")

    nonbonded_method_value = getattr(app, normalized_nonbonded_method)
    pdb = app.PDBFile(StringIO(prepared_structure_pdb_text))
    forcefield = app.ForceField(*forcefield_files)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    hydrogen_mass = hydrogen_mass_da * unit.amu if hydrogen_mass_repartitioning else None
    input_residue_templates = _build_explicit_residue_templates(pdb.topology)
    with _patch_forcefield_create_system_for_truncated_chains(forcefield) as original_create_system:
        modeller.addHydrogens(forcefield, residueTemplates=input_residue_templates)
        modeller_residue_templates = _build_explicit_residue_templates(modeller.topology)
        system = original_create_system(
            modeller.topology,
            nonbondedMethod=nonbonded_method_value,
            nonbondedCutoff=cutoff_nm * unit.nanometer,
            constraints=app.HBonds,
            hydrogenMass=hydrogen_mass,
            residueTemplates=modeller_residue_templates,
            ignoreExternalBonds=True,
        )
    integrator = LangevinMiddleIntegrator(
        target_temperature_K * unit.kelvin,
        friction_coeff_ps / unit.picosecond,
        timestep_fs * unit.femtoseconds,
    )
    prepared_structure_handle = StringIO()
    app.PDBFile.writeFile(
        modeller.topology,
        modeller.positions,
        prepared_structure_handle,
        keepIds=True,
    )
    return OpenMMSerializedSystem(
        prepared_structure_pdb_text=_sanitize_prepared_structure_pdb_text(prepared_structure_handle.getvalue()),
        system_xml=XmlSerializer.serialize(system),
        integrator_xml=XmlSerializer.serialize(integrator),
        atom_count=modeller.topology.getNumAtoms(),
        residue_count=modeller.topology.getNumResidues(),
    )
