"""Deterministic well-tempered metadynamics helpers for Phase 6C1."""

from __future__ import annotations

import math
from dataclasses import dataclass
from io import StringIO
from typing import Sequence
from xml.etree import ElementTree

import numpy as np

from lanm.filesystem import atomic_write_text, write_csv_rows
from lanm.md.openmm_equilibration_smoke import (
    _center_positions_in_box,
    _orthorhombic_box_lengths_nm,
    _orthorhombic_box_vectors,
    _render_pdb_text,
)

DEFAULT_METADYNAMICS_DURATION_PS = 500.0
DEFAULT_METADYNAMICS_REPORT_STRIDE_STEPS = 500
DEFAULT_METADYNAMICS_DEPOSITION_FREQUENCY_STEPS = 250
DEFAULT_METADYNAMICS_BIAS_FACTOR = 8.0
DEFAULT_METADYNAMICS_GAUSSIAN_HEIGHT_KJ_PER_MOL = 1.0
DEFAULT_DONOR_CUTOFF_A = 3.2
DEFAULT_COORDINATION_SWITCH_DISTANCE_A = 3.2
DEFAULT_COORDINATION_SWITCH_POWER = 12
DEFAULT_COORDINATION_NUMBER_MIN = 0.0
DEFAULT_COORDINATION_NUMBER_MAX = 12.0
DEFAULT_COORDINATION_NUMBER_BIAS_WIDTH = 0.25
DEFAULT_MEAN_DISTANCE_MIN_A = 1.8
DEFAULT_MEAN_DISTANCE_MAX_A = 6.5
DEFAULT_MEAN_DISTANCE_BIAS_WIDTH_A = 0.2
DEFAULT_ESCAPE_COORDINATION_THRESHOLD = 5.0
DEFAULT_ESCAPE_DISTANCE_THRESHOLD_A = 3.2
DEFAULT_NPT_BOX_PADDING_NM = 2.0
DEFAULT_METADYNAMICS_CPU_THREAD_COUNT = 3


@dataclass(frozen=True, slots=True)
class CollectiveVariableAtomSelection:
    metal_atom_indices: tuple[int, ...]
    pocket_oxygen_atom_indices: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class CollectiveVariableDescriptor:
    name: str
    minimum_value: float
    maximum_value: float
    bias_width: float
    grid_width: int
    units: str


@dataclass(slots=True)
class OpenMMMetadynamicsCollectiveVariableSetup:
    atom_selection: CollectiveVariableAtomSelection
    coordination_number_variable: object
    mean_metal_oxygen_distance_variable: object
    coordination_number_descriptor: CollectiveVariableDescriptor
    mean_metal_oxygen_distance_descriptor: CollectiveVariableDescriptor


@dataclass(frozen=True, slots=True)
class MetadynamicsCVTimeSeriesRow:
    step: int
    time_ps: float
    coordination_number: float
    mean_metal_oxygen_distance_A: float
    bias_potential_kj_per_mol: float
    escape_condition_met: str


@dataclass(frozen=True, slots=True)
class OpenMMMetadynamicsArtifacts:
    cv_timeseries_rows: tuple[MetadynamicsCVTimeSeriesRow, ...]
    bias_state_xml: str
    final_state_xml: str
    final_structure_pdb_text: str
    platform_name: str
    platform_properties: str
    metal_atom_count: int
    pocket_oxygen_atom_count: int
    minimum_coordination_number: float
    maximum_mean_metal_oxygen_distance_A: float
    escape_event_detected: bool
    final_coordination_number: float
    final_mean_metal_oxygen_distance_A: float


def _choose_metadynamics_cuda_platform() -> tuple[object, dict[str, str], str, str]:
    try:
        from openmm import Platform
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for Phase 6C1 metadynamics") from exc

    properties = {
        "Precision": "mixed",
        "DeviceIndex": "0",
        "DeterministicForces": "true",
        "UseCpuPme": "false",
    }
    try:
        platform = Platform.getPlatformByName("CUDA")
    except Exception as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("CUDA OpenMM platform is required for Phase 6C1 metadynamics") from exc
    properties_display = ", ".join(f"{key}={value}" for key, value in properties.items())
    return platform, properties, "CUDA", properties_display


def _box_lengths_from_vectors_nm(box_vectors: object) -> tuple[float, float, float]:
    try:
        from openmm import unit
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for Phase 6C1 metadynamics") from exc

    return (
        abs(float(box_vectors[0][0].value_in_unit(unit.nanometer))),
        abs(float(box_vectors[1][1].value_in_unit(unit.nanometer))),
        abs(float(box_vectors[2][2].value_in_unit(unit.nanometer))),
    )


def _box_volume_from_vectors_nm3(box_vectors: object) -> float:
    lengths = _box_lengths_from_vectors_nm(box_vectors)
    return float(lengths[0] * lengths[1] * lengths[2])


def _distance_angstroms(coord_a: Sequence[float], coord_b: Sequence[float]) -> float:
    dx = float(coord_a[0] - coord_b[0])
    dy = float(coord_a[1] - coord_b[1])
    dz = float(coord_a[2] - coord_b[2])
    return math.sqrt((dx * dx) + (dy * dy) + (dz * dz))


def discover_local_pocket_cv_atoms(
    *,
    equilibrated_structure_pdb_text: str,
    target_metal: str,
    donor_cutoff_A: float = DEFAULT_DONOR_CUTOFF_A,
) -> CollectiveVariableAtomSelection:
    """Select target-metal atoms and their starting first-shell pocket oxygens."""
    try:
        from openmm import app, unit
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for Phase 6C1 metadynamics") from exc

    if donor_cutoff_A <= 0.0:
        raise ValueError("donor_cutoff_A must be positive")

    pdb = app.PDBFile(StringIO(equilibrated_structure_pdb_text))
    coordinates_A = pdb.positions.value_in_unit(unit.angstrom)
    normalized_target_metal = target_metal.strip().upper()

    metal_atom_indices: list[int] = []
    for atom in pdb.topology.atoms():
        symbol = atom.element.symbol.upper() if atom.element is not None else ""
        if symbol == normalized_target_metal:
            metal_atom_indices.append(atom.index)
    if not metal_atom_indices:
        raise ValueError(f"Could not find any {target_metal} atoms in the metadynamics topology")

    pocket_oxygen_atom_indices: list[int] = []
    for atom in pdb.topology.atoms():
        symbol = atom.element.symbol.upper() if atom.element is not None else ""
        if symbol != "O":
            continue
        oxygen_position = coordinates_A[atom.index]
        if any(
            _distance_angstroms(oxygen_position, coordinates_A[metal_index]) <= donor_cutoff_A
            for metal_index in metal_atom_indices
        ):
            pocket_oxygen_atom_indices.append(atom.index)
    if not pocket_oxygen_atom_indices:
        raise ValueError(
            f"Could not find any pocket oxygen atoms within {donor_cutoff_A:.1f} A of {target_metal}"
        )

    return CollectiveVariableAtomSelection(
        metal_atom_indices=tuple(metal_atom_indices),
        pocket_oxygen_atom_indices=tuple(sorted(set(pocket_oxygen_atom_indices))),
    )


def build_metadynamics_collective_variable_setup(
    *,
    atom_selection: CollectiveVariableAtomSelection,
    coordination_switch_distance_A: float = DEFAULT_COORDINATION_SWITCH_DISTANCE_A,
    coordination_switch_power: int = DEFAULT_COORDINATION_SWITCH_POWER,
    coordination_number_min: float = DEFAULT_COORDINATION_NUMBER_MIN,
    coordination_number_max: float = DEFAULT_COORDINATION_NUMBER_MAX,
    coordination_number_bias_width: float = DEFAULT_COORDINATION_NUMBER_BIAS_WIDTH,
    mean_distance_min_A: float = DEFAULT_MEAN_DISTANCE_MIN_A,
    mean_distance_max_A: float = DEFAULT_MEAN_DISTANCE_MAX_A,
    mean_distance_bias_width_A: float = DEFAULT_MEAN_DISTANCE_BIAS_WIDTH_A,
) -> OpenMMMetadynamicsCollectiveVariableSetup:
    """Build deterministic coordination-number and mean-distance bias variables."""
    try:
        from openmm import CustomBondForce, CustomCVForce
        from openmm.app.metadynamics import BiasVariable
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for Phase 6C1 metadynamics") from exc

    if coordination_switch_distance_A <= 0.0:
        raise ValueError("coordination_switch_distance_A must be positive")
    if coordination_switch_power <= 0:
        raise ValueError("coordination_switch_power must be positive")
    if coordination_number_min >= coordination_number_max:
        raise ValueError("coordination_number bounds must be increasing")
    if coordination_number_bias_width <= 0.0:
        raise ValueError("coordination_number_bias_width must be positive")
    if mean_distance_min_A >= mean_distance_max_A:
        raise ValueError("mean_distance bounds must be increasing")
    if mean_distance_bias_width_A <= 0.0:
        raise ValueError("mean_distance_bias_width_A must be positive")
    if not atom_selection.metal_atom_indices:
        raise ValueError("Expected at least one target-metal atom for collective-variable setup")
    if not atom_selection.pocket_oxygen_atom_indices:
        raise ValueError("Expected at least one pocket oxygen atom for collective-variable setup")

    switch_distance_nm = coordination_switch_distance_A / 10.0
    mean_distance_min_nm = mean_distance_min_A / 10.0
    mean_distance_max_nm = mean_distance_max_A / 10.0
    mean_distance_bias_width_nm = mean_distance_bias_width_A / 10.0
    metal_count = len(atom_selection.metal_atom_indices)

    coordination_expression = f"1/(1 + (r/{switch_distance_nm:.8f})^{coordination_switch_power})"
    weighted_distance_expression = (
        f"r/(1 + (r/{switch_distance_nm:.8f})^{coordination_switch_power})"
    )
    coordination_cv_force = CustomBondForce(
        f"{1.0 / metal_count:.8f}/(1 + (r/{switch_distance_nm:.8f})^{coordination_switch_power})"
    )
    coordination_sum_force = CustomBondForce(coordination_expression)
    weighted_distance_force = CustomBondForce(weighted_distance_expression)
    for metal_index in atom_selection.metal_atom_indices:
        for oxygen_index in atom_selection.pocket_oxygen_atom_indices:
            coordination_cv_force.addBond(metal_index, oxygen_index)
            coordination_sum_force.addBond(metal_index, oxygen_index)
            weighted_distance_force.addBond(metal_index, oxygen_index)

    mean_distance_cv_force = CustomCVForce("weighted_distance_sum/(coordination_sum + 1.0e-6)")
    mean_distance_cv_force.addCollectiveVariable("weighted_distance_sum", weighted_distance_force)
    mean_distance_cv_force.addCollectiveVariable("coordination_sum", coordination_sum_force)

    coordination_variable = BiasVariable(
        coordination_cv_force,
        coordination_number_min,
        coordination_number_max,
        coordination_number_bias_width,
        periodic=False,
    )
    mean_distance_variable = BiasVariable(
        mean_distance_cv_force,
        mean_distance_min_nm,
        mean_distance_max_nm,
        mean_distance_bias_width_nm,
        periodic=False,
    )

    return OpenMMMetadynamicsCollectiveVariableSetup(
        atom_selection=atom_selection,
        coordination_number_variable=coordination_variable,
        mean_metal_oxygen_distance_variable=mean_distance_variable,
        coordination_number_descriptor=CollectiveVariableDescriptor(
            name="coordination_number",
            minimum_value=coordination_number_min,
            maximum_value=coordination_number_max,
            bias_width=coordination_number_bias_width,
            grid_width=coordination_variable.gridWidth,
            units="dimensionless",
        ),
        mean_metal_oxygen_distance_descriptor=CollectiveVariableDescriptor(
            name="mean_metal_oxygen_distance_A",
            minimum_value=mean_distance_min_A,
            maximum_value=mean_distance_max_A,
            bias_width=mean_distance_bias_width_A,
            grid_width=mean_distance_variable.gridWidth,
            units="angstrom",
        ),
    )


def _prepare_periodic_metadynamics_system(
    *,
    system_xml: str,
    starting_positions: object,
    starting_box_vectors: object,
    cutoff_nm: float,
    box_padding_nm: float,
) -> tuple[object, object, object]:
    try:
        from openmm import XmlSerializer, unit
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for Phase 6C1 metadynamics") from exc

    periodic_system = XmlSerializer.deserialize(system_xml)
    for force_index in range(periodic_system.getNumForces()):
        force = periodic_system.getForce(force_index)
        if force.__class__.__name__ != "NonbondedForce":
            continue
        if force.getNonbondedMethod() in (force.NoCutoff, force.CutoffNonPeriodic):
            force.setNonbondedMethod(force.CutoffPeriodic)
        force.setCutoffDistance(cutoff_nm * unit.nanometer)

    periodic_box_vectors = starting_box_vectors
    periodic_positions = starting_positions
    if _box_volume_from_vectors_nm3(starting_box_vectors) <= 0.0:
        box_lengths_nm = _orthorhombic_box_lengths_nm(
            starting_positions,
            cutoff_nm=cutoff_nm,
            padding_nm=box_padding_nm,
        )
        periodic_box_vectors = _orthorhombic_box_vectors(box_lengths_nm)
        periodic_positions = _center_positions_in_box(starting_positions, box_lengths_nm)

    periodic_system.setDefaultPeriodicBoxVectors(*periodic_box_vectors)
    return periodic_system, periodic_positions, periodic_box_vectors


def _detect_escape_event(
    *,
    coordination_number: float,
    mean_metal_oxygen_distance_A: float,
    coordination_threshold: float,
    distance_threshold_A: float,
) -> bool:
    return (
        coordination_number <= coordination_threshold
        and mean_metal_oxygen_distance_A >= distance_threshold_A
    )


def _sample_collective_variables(
    *,
    simulation: object,
    metadynamics: object,
    bias_force_group: int,
    current_step: int,
    timestep_fs: float,
    escape_coordination_threshold: float,
    escape_distance_threshold_A: float,
) -> MetadynamicsCVTimeSeriesRow:
    try:
        from openmm import unit
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for Phase 6C1 metadynamics") from exc

    coordination_number, mean_distance_nm = tuple(metadynamics.getCollectiveVariables(simulation))
    mean_distance_A = float(mean_distance_nm) * 10.0
    bias_energy_kj_per_mol = float(
        simulation.context.getState(energy=True, groups=1 << bias_force_group)
        .getPotentialEnergy()
        .value_in_unit(unit.kilojoules_per_mole)
    )
    escape_condition_met = _detect_escape_event(
        coordination_number=float(coordination_number),
        mean_metal_oxygen_distance_A=mean_distance_A,
        coordination_threshold=escape_coordination_threshold,
        distance_threshold_A=escape_distance_threshold_A,
    )
    return MetadynamicsCVTimeSeriesRow(
        step=current_step,
        time_ps=(current_step * timestep_fs) / 1000.0,
        coordination_number=float(coordination_number),
        mean_metal_oxygen_distance_A=mean_distance_A,
        bias_potential_kj_per_mol=bias_energy_kj_per_mol,
        escape_condition_met="TRUE" if escape_condition_met else "FALSE",
    )


def _serialize_bias_state_xml(
    *,
    cv_setup: OpenMMMetadynamicsCollectiveVariableSetup,
    metadynamics: object,
    final_row: MetadynamicsCVTimeSeriesRow,
    platform_name: str,
    platform_properties: str,
) -> str:
    total_bias_kj_per_mol = np.asarray(metadynamics._totalBias, dtype=float)

    root = ElementTree.Element("metadynamics_bias_state")
    ElementTree.SubElement(
        root,
        "platform",
        name=platform_name,
        properties=platform_properties,
    )
    ElementTree.SubElement(
        root,
        "atom_selection",
        metal_atom_count=str(len(cv_setup.atom_selection.metal_atom_indices)),
        pocket_oxygen_atom_count=str(len(cv_setup.atom_selection.pocket_oxygen_atom_indices)),
    )
    ElementTree.SubElement(
        root,
        "variable",
        name=cv_setup.coordination_number_descriptor.name,
        minimum=f"{cv_setup.coordination_number_descriptor.minimum_value:.6f}",
        maximum=f"{cv_setup.coordination_number_descriptor.maximum_value:.6f}",
        bias_width=f"{cv_setup.coordination_number_descriptor.bias_width:.6f}",
        grid_width=str(cv_setup.coordination_number_descriptor.grid_width),
        units=cv_setup.coordination_number_descriptor.units,
    )
    ElementTree.SubElement(
        root,
        "variable",
        name=cv_setup.mean_metal_oxygen_distance_descriptor.name,
        minimum=f"{cv_setup.mean_metal_oxygen_distance_descriptor.minimum_value:.6f}",
        maximum=f"{cv_setup.mean_metal_oxygen_distance_descriptor.maximum_value:.6f}",
        bias_width=f"{cv_setup.mean_metal_oxygen_distance_descriptor.bias_width:.6f}",
        grid_width=str(cv_setup.mean_metal_oxygen_distance_descriptor.grid_width),
        units=cv_setup.mean_metal_oxygen_distance_descriptor.units,
    )
    ElementTree.SubElement(
        root,
        "final_sample",
        coordination_number=f"{final_row.coordination_number:.6f}",
        mean_metal_oxygen_distance_A=f"{final_row.mean_metal_oxygen_distance_A:.6f}",
        bias_potential_kj_per_mol=f"{final_row.bias_potential_kj_per_mol:.6f}",
        step=str(final_row.step),
        time_ps=f"{final_row.time_ps:.6f}",
    )
    grid_element = ElementTree.SubElement(
        root,
        "total_bias_grid_kj_per_mol",
        axis_order="mean_metal_oxygen_distance_A,coordination_number",
        shape="x".join(str(length) for length in total_bias_kj_per_mol.shape),
    )
    grid_element.text = " ".join(f"{value:.8f}" for value in total_bias_kj_per_mol.ravel(order="C"))
    ElementTree.indent(root, space="  ")
    return ElementTree.tostring(root, encoding="unicode") + "\n"


def run_openmm_metadynamics_smoke_system(
    *,
    equilibrated_structure_pdb_text: str,
    system_xml: str,
    equilibrated_state_xml: str,
    target_metal: str,
    target_temperature_K: int,
    friction_coeff_ps: float,
    timestep_fs: float,
    cutoff_nm: float,
    duration_ps: float = DEFAULT_METADYNAMICS_DURATION_PS,
    report_stride_steps: int = DEFAULT_METADYNAMICS_REPORT_STRIDE_STEPS,
    bias_factor: float = DEFAULT_METADYNAMICS_BIAS_FACTOR,
    gaussian_height_kj_per_mol: float = DEFAULT_METADYNAMICS_GAUSSIAN_HEIGHT_KJ_PER_MOL,
    deposition_frequency_steps: int = DEFAULT_METADYNAMICS_DEPOSITION_FREQUENCY_STEPS,
    donor_cutoff_A: float = DEFAULT_DONOR_CUTOFF_A,
    coordination_switch_distance_A: float = DEFAULT_COORDINATION_SWITCH_DISTANCE_A,
    coordination_switch_power: int = DEFAULT_COORDINATION_SWITCH_POWER,
    coordination_number_min: float = DEFAULT_COORDINATION_NUMBER_MIN,
    coordination_number_max: float = DEFAULT_COORDINATION_NUMBER_MAX,
    coordination_number_bias_width: float = DEFAULT_COORDINATION_NUMBER_BIAS_WIDTH,
    mean_distance_min_A: float = DEFAULT_MEAN_DISTANCE_MIN_A,
    mean_distance_max_A: float = DEFAULT_MEAN_DISTANCE_MAX_A,
    mean_distance_bias_width_A: float = DEFAULT_MEAN_DISTANCE_BIAS_WIDTH_A,
    escape_coordination_threshold: float = DEFAULT_ESCAPE_COORDINATION_THRESHOLD,
    escape_distance_threshold_A: float = DEFAULT_ESCAPE_DISTANCE_THRESHOLD_A,
    seed: int = 101,
    cv_timeseries_path: object | None = None,
    bias_state_path: object | None = None,
    final_state_path: object | None = None,
    final_structure_path: object | None = None,
    box_padding_nm: float = DEFAULT_NPT_BOX_PADDING_NM,
    cpu_thread_count: int = DEFAULT_METADYNAMICS_CPU_THREAD_COUNT,
) -> OpenMMMetadynamicsArtifacts:
    """Run one deterministic well-tempered metadynamics triage simulation."""
    try:
        from openmm import LangevinMiddleIntegrator, XmlSerializer, app, unit
        from openmm.app.metadynamics import Metadynamics
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for Phase 6C1 metadynamics") from exc

    if target_temperature_K <= 0:
        raise ValueError("target_temperature_K must be positive")
    if friction_coeff_ps <= 0.0:
        raise ValueError("friction_coeff_ps must be positive")
    if timestep_fs <= 0.0:
        raise ValueError("timestep_fs must be positive")
    if cutoff_nm <= 0.0:
        raise ValueError("cutoff_nm must be positive")
    if duration_ps <= 0.0:
        raise ValueError("duration_ps must be positive")
    if report_stride_steps <= 0:
        raise ValueError("report_stride_steps must be positive")
    if bias_factor <= 1.0:
        raise ValueError("bias_factor must be greater than 1")
    if gaussian_height_kj_per_mol <= 0.0:
        raise ValueError("gaussian_height_kj_per_mol must be positive")
    if deposition_frequency_steps <= 0:
        raise ValueError("deposition_frequency_steps must be positive")
    if seed < 0:
        raise ValueError("seed must be non-negative")
    if (
        cv_timeseries_path is None
        or bias_state_path is None
        or final_state_path is None
        or final_structure_path is None
    ):
        raise ValueError("All Phase 6C1 artifact output paths must be provided")

    step_count = int(round((duration_ps * 1000.0) / timestep_fs))
    if step_count <= 0:
        raise ValueError("Phase 6C1 run must contain at least one integration step")

    pdb = app.PDBFile(StringIO(equilibrated_structure_pdb_text))
    atom_selection = discover_local_pocket_cv_atoms(
        equilibrated_structure_pdb_text=equilibrated_structure_pdb_text,
        target_metal=target_metal,
        donor_cutoff_A=donor_cutoff_A,
    )
    cv_setup = build_metadynamics_collective_variable_setup(
        atom_selection=atom_selection,
        coordination_switch_distance_A=coordination_switch_distance_A,
        coordination_switch_power=coordination_switch_power,
        coordination_number_min=coordination_number_min,
        coordination_number_max=coordination_number_max,
        coordination_number_bias_width=coordination_number_bias_width,
        mean_distance_min_A=mean_distance_min_A,
        mean_distance_max_A=mean_distance_max_A,
        mean_distance_bias_width_A=mean_distance_bias_width_A,
    )
    equilibrated_state = XmlSerializer.deserialize(equilibrated_state_xml)
    metadynamics_system, starting_positions, starting_box_vectors = _prepare_periodic_metadynamics_system(
        system_xml=system_xml,
        starting_positions=equilibrated_state.getPositions(),
        starting_box_vectors=equilibrated_state.getPeriodicBoxVectors(),
        cutoff_nm=cutoff_nm,
        box_padding_nm=box_padding_nm,
    )
    metadynamics = Metadynamics(
        metadynamics_system,
        [
            cv_setup.coordination_number_variable,
            cv_setup.mean_metal_oxygen_distance_variable,
        ],
        target_temperature_K * unit.kelvin,
        bias_factor,
        gaussian_height_kj_per_mol * unit.kilojoules_per_mole,
        deposition_frequency_steps,
    )
    integrator = LangevinMiddleIntegrator(
        target_temperature_K * unit.kelvin,
        friction_coeff_ps / unit.picosecond,
        timestep_fs * unit.femtoseconds,
    )
    integrator.setRandomNumberSeed(seed)
    platform, platform_properties, platform_name, platform_properties_display = (
        _choose_metadynamics_cuda_platform()
    )
    simulation = app.Simulation(
        pdb.topology,
        metadynamics_system,
        integrator,
        platform,
        platform_properties,
    )
    simulation.context.setPeriodicBoxVectors(*starting_box_vectors)
    simulation.context.setPositions(starting_positions)
    starting_velocities = equilibrated_state.getVelocities()
    if starting_velocities is None:
        simulation.context.setVelocitiesToTemperature(target_temperature_K * unit.kelvin, seed)
    else:
        simulation.context.setVelocities(starting_velocities)

    bias_force_group = metadynamics._force.getForceGroup()
    time_series_rows: list[MetadynamicsCVTimeSeriesRow] = [
        _sample_collective_variables(
            simulation=simulation,
            metadynamics=metadynamics,
            bias_force_group=bias_force_group,
            current_step=0,
            timestep_fs=timestep_fs,
            escape_coordination_threshold=escape_coordination_threshold,
            escape_distance_threshold_A=escape_distance_threshold_A,
        )
    ]
    integrated_steps = 0
    while integrated_steps < step_count:
        chunk_steps = min(report_stride_steps, step_count - integrated_steps)
        metadynamics.step(simulation, chunk_steps)
        integrated_steps += chunk_steps
        time_series_rows.append(
            _sample_collective_variables(
                simulation=simulation,
                metadynamics=metadynamics,
                bias_force_group=bias_force_group,
                current_step=integrated_steps,
                timestep_fs=timestep_fs,
                escape_coordination_threshold=escape_coordination_threshold,
                escape_distance_threshold_A=escape_distance_threshold_A,
            )
        )

    final_state = simulation.context.getState(
        getPositions=True,
        getVelocities=True,
        getEnergy=True,
        enforcePeriodicBox=True,
    )
    final_row = time_series_rows[-1]
    bias_state_xml = _serialize_bias_state_xml(
        cv_setup=cv_setup,
        metadynamics=metadynamics,
        final_row=final_row,
        platform_name=platform_name,
        platform_properties=platform_properties_display,
    )
    final_state_xml = XmlSerializer.serialize(final_state)
    final_structure_pdb_text = _render_pdb_text(pdb.topology, final_state.getPositions())

    write_csv_rows(cv_timeseries_path, tuple(time_series_rows))
    atomic_write_text(bias_state_path, bias_state_xml)
    atomic_write_text(final_state_path, final_state_xml)
    atomic_write_text(final_structure_path, final_structure_pdb_text)

    return OpenMMMetadynamicsArtifacts(
        cv_timeseries_rows=tuple(time_series_rows),
        bias_state_xml=bias_state_xml,
        final_state_xml=final_state_xml,
        final_structure_pdb_text=final_structure_pdb_text,
        platform_name=platform_name,
        platform_properties=platform_properties_display,
        metal_atom_count=len(atom_selection.metal_atom_indices),
        pocket_oxygen_atom_count=len(atom_selection.pocket_oxygen_atom_indices),
        minimum_coordination_number=min(row.coordination_number for row in time_series_rows),
        maximum_mean_metal_oxygen_distance_A=max(
            row.mean_metal_oxygen_distance_A for row in time_series_rows
        ),
        escape_event_detected=any(row.escape_condition_met == "TRUE" for row in time_series_rows),
        final_coordination_number=final_row.coordination_number,
        final_mean_metal_oxygen_distance_A=final_row.mean_metal_oxygen_distance_A,
    )
