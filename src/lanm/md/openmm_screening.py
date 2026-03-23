"""Deterministic Phase 6B2 OpenMM screening helpers."""

from __future__ import annotations

import math
from dataclasses import dataclass
from io import StringIO
from pathlib import Path

from lanm.filesystem import atomic_write_text
from lanm.md.openmm_equilibration_smoke import (
    DEFAULT_BAROSTAT_FREQUENCY_STEPS,
    DEFAULT_CPU_THREAD_COUNT,
    DEFAULT_NPT_BOX_PADDING_NM,
    _center_positions_in_box,
    _choose_openmm_platform,
    _orthorhombic_box_lengths_nm,
    _orthorhombic_box_vectors,
    _render_pdb_text,
)

DEFAULT_SCREENING_DURATION_PS = 100.0
DEFAULT_PRESSURE_ATM = 1.0
DEFAULT_REPORT_STRIDE_STEPS = 500
DEFAULT_INNER_SPHERE_CUTOFF_A = 3.2


@dataclass(frozen=True, slots=True)
class OpenMMScreeningArtifacts:
    platform_name: str
    platform_properties: str
    metal_atom_count: int
    oxygen_atom_count: int
    frame_count: int
    mean_min_metal_oxygen_distance_A: float
    inner_sphere_occupancy_fraction: float


def _box_lengths_from_vectors_nm(box_vectors: object) -> tuple[float, float, float]:
    try:
        from openmm import unit
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for Phase 6B2 screening") from exc

    return (
        abs(float(box_vectors[0][0].value_in_unit(unit.nanometer))),
        abs(float(box_vectors[1][1].value_in_unit(unit.nanometer))),
        abs(float(box_vectors[2][2].value_in_unit(unit.nanometer))),
    )


def _box_volume_from_vectors_nm3(box_vectors: object) -> float:
    lengths = _box_lengths_from_vectors_nm(box_vectors)
    return float(lengths[0] * lengths[1] * lengths[2])


def _discover_metal_and_oxygen_atom_indices(
    *,
    topology: object,
    target_metal: str,
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    metal_atom_indices: list[int] = []
    oxygen_atom_indices: list[int] = []
    for atom in topology.atoms():
        symbol = atom.element.symbol if atom.element is not None else ""
        if symbol == target_metal:
            metal_atom_indices.append(atom.index)
        if symbol == "O":
            oxygen_atom_indices.append(atom.index)
    if not metal_atom_indices:
        raise ValueError(f"Could not find any {target_metal} atoms in the screening topology")
    if not oxygen_atom_indices:
        raise ValueError("Could not find any oxygen atoms in the screening topology")
    return tuple(metal_atom_indices), tuple(oxygen_atom_indices)


class MetalOxygenScreeningReporter:
    """Collect coarse metal-oxygen screening metrics during production."""

    def __init__(
        self,
        *,
        report_interval: int,
        metal_atom_indices: tuple[int, ...],
        oxygen_atom_indices: tuple[int, ...],
        inner_sphere_cutoff_A: float,
    ) -> None:
        if report_interval <= 0:
            raise ValueError("report_interval must be positive")
        if inner_sphere_cutoff_A <= 0.0:
            raise ValueError("inner_sphere_cutoff_A must be positive")
        self._report_interval = report_interval
        self._metal_atom_indices = metal_atom_indices
        self._oxygen_atom_indices = oxygen_atom_indices
        self._inner_sphere_cutoff_A = inner_sphere_cutoff_A
        self._frame_mean_min_distances_A: list[float] = []
        self._frame_inner_sphere_occupancy_fractions: list[float] = []

    def describeNextReport(self, simulation: object) -> tuple[int, bool, bool, bool, bool, bool]:
        steps = self._report_interval - (simulation.currentStep % self._report_interval)
        return (steps, True, False, False, False, True)

    def report(self, simulation: object, state: object) -> None:
        del simulation
        try:
            from openmm import unit
        except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
            raise RuntimeError("OpenMM is required for Phase 6B2 screening") from exc

        coordinates = state.getPositions().value_in_unit(unit.nanometer)
        box_lengths_nm = _box_lengths_from_vectors_nm(state.getPeriodicBoxVectors())

        oxygen_coordinates = [coordinates[index] for index in self._oxygen_atom_indices]
        frame_min_distances_A: list[float] = []
        for metal_index in self._metal_atom_indices:
            metal_coordinate = coordinates[metal_index]
            minimum_distance_nm = math.inf
            for oxygen_coordinate in oxygen_coordinates:
                dx = float(oxygen_coordinate[0] - metal_coordinate[0])
                dy = float(oxygen_coordinate[1] - metal_coordinate[1])
                dz = float(oxygen_coordinate[2] - metal_coordinate[2])
                if box_lengths_nm[0] > 0.0:
                    dx -= box_lengths_nm[0] * round(dx / box_lengths_nm[0])
                if box_lengths_nm[1] > 0.0:
                    dy -= box_lengths_nm[1] * round(dy / box_lengths_nm[1])
                if box_lengths_nm[2] > 0.0:
                    dz -= box_lengths_nm[2] * round(dz / box_lengths_nm[2])
                distance_nm = math.sqrt((dx * dx) + (dy * dy) + (dz * dz))
                if distance_nm < minimum_distance_nm:
                    minimum_distance_nm = distance_nm
            frame_min_distances_A.append(minimum_distance_nm * 10.0)

        frame_mean_min_distance_A = sum(frame_min_distances_A) / len(frame_min_distances_A)
        frame_inner_sphere_fraction = (
            sum(distance_A <= self._inner_sphere_cutoff_A for distance_A in frame_min_distances_A)
            / len(frame_min_distances_A)
        )
        self._frame_mean_min_distances_A.append(frame_mean_min_distance_A)
        self._frame_inner_sphere_occupancy_fractions.append(frame_inner_sphere_fraction)

    @property
    def frame_count(self) -> int:
        return len(self._frame_mean_min_distances_A)

    @property
    def mean_min_metal_oxygen_distance_A(self) -> float:
        if not self._frame_mean_min_distances_A:
            raise ValueError("No screening frames were collected")
        return sum(self._frame_mean_min_distances_A) / len(self._frame_mean_min_distances_A)

    @property
    def inner_sphere_occupancy_fraction(self) -> float:
        if not self._frame_inner_sphere_occupancy_fractions:
            raise ValueError("No screening frames were collected")
        return sum(self._frame_inner_sphere_occupancy_fractions) / len(
            self._frame_inner_sphere_occupancy_fractions
        )


def _prepare_periodic_screening_system(
    *,
    system_xml: str,
    starting_positions: object,
    starting_box_vectors: object,
    target_temperature_K: int,
    pressure_atm: float,
    cutoff_nm: float,
    box_padding_nm: float,
    barostat_frequency_steps: int,
    seed: int,
) -> tuple[object, object, object]:
    try:
        from openmm import MonteCarloBarostat, XmlSerializer, unit
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for Phase 6B2 screening") from exc

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
    barostat = MonteCarloBarostat(
        pressure_atm * unit.atmosphere,
        target_temperature_K * unit.kelvin,
        barostat_frequency_steps,
    )
    barostat.setRandomNumberSeed(seed)
    periodic_system.addForce(barostat)
    return periodic_system, periodic_positions, periodic_box_vectors


def run_openmm_screening_replicate(
    *,
    prepared_structure_pdb_text: str,
    system_xml: str,
    equilibrated_state_xml: str,
    target_metal: str,
    target_temperature_K: int,
    friction_coeff_ps: float,
    timestep_fs: float,
    cutoff_nm: float,
    production_duration_ps: float = DEFAULT_SCREENING_DURATION_PS,
    pressure_atm: float = DEFAULT_PRESSURE_ATM,
    report_stride_steps: int = DEFAULT_REPORT_STRIDE_STEPS,
    inner_sphere_cutoff_A: float = DEFAULT_INNER_SPHERE_CUTOFF_A,
    seed: int = 101,
    state_data_path: Path | None = None,
    final_state_path: Path | None = None,
    final_structure_path: Path | None = None,
    screening_log_path: Path | None = None,
    trajectory_path: Path | None = None,
    barostat_frequency_steps: int = DEFAULT_BAROSTAT_FREQUENCY_STEPS,
    box_padding_nm: float = DEFAULT_NPT_BOX_PADDING_NM,
    cpu_thread_count: int = DEFAULT_CPU_THREAD_COUNT,
) -> OpenMMScreeningArtifacts:
    """Run one deterministic Phase 6B2 NPT screening replicate."""
    try:
        from openmm import LangevinMiddleIntegrator, XmlSerializer, app, unit
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for Phase 6B2 screening") from exc

    if target_temperature_K <= 0:
        raise ValueError("target_temperature_K must be positive")
    if friction_coeff_ps <= 0.0:
        raise ValueError("friction_coeff_ps must be positive")
    if timestep_fs <= 0.0:
        raise ValueError("timestep_fs must be positive")
    if cutoff_nm <= 0.0:
        raise ValueError("cutoff_nm must be positive")
    if production_duration_ps <= 0.0:
        raise ValueError("production_duration_ps must be positive")
    if pressure_atm <= 0.0:
        raise ValueError("pressure_atm must be positive")
    if report_stride_steps <= 0:
        raise ValueError("report_stride_steps must be positive")
    if seed < 0:
        raise ValueError("seed must be non-negative")

    if state_data_path is None or final_state_path is None or final_structure_path is None or trajectory_path is None:
        raise ValueError("All Phase 6B2 artifact output paths must be provided")
    del screening_log_path

    for artifact_path in (state_data_path, final_state_path, final_structure_path, trajectory_path):
        artifact_path.parent.mkdir(parents=True, exist_ok=True)

    step_count = int(round((production_duration_ps * 1000.0) / timestep_fs))
    if step_count <= 0:
        raise ValueError("Phase 6B2 replicate must contain at least one production step")

    pdb = app.PDBFile(StringIO(prepared_structure_pdb_text))
    metal_atom_indices, oxygen_atom_indices = _discover_metal_and_oxygen_atom_indices(
        topology=pdb.topology,
        target_metal=target_metal,
    )
    equilibrated_state = XmlSerializer.deserialize(equilibrated_state_xml)
    screening_system, starting_positions, starting_box_vectors = _prepare_periodic_screening_system(
        system_xml=system_xml,
        starting_positions=equilibrated_state.getPositions(),
        starting_box_vectors=equilibrated_state.getPeriodicBoxVectors(),
        target_temperature_K=target_temperature_K,
        pressure_atm=pressure_atm,
        cutoff_nm=cutoff_nm,
        box_padding_nm=box_padding_nm,
        barostat_frequency_steps=barostat_frequency_steps,
        seed=seed,
    )
    integrator = LangevinMiddleIntegrator(
        target_temperature_K * unit.kelvin,
        friction_coeff_ps / unit.picosecond,
        timestep_fs * unit.femtoseconds,
    )
    integrator.setRandomNumberSeed(seed)
    platform, platform_properties, platform_name, platform_properties_display = _choose_openmm_platform(
        cpu_thread_count
    )
    simulation = app.Simulation(
        pdb.topology,
        screening_system,
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

    state_data_reporter = app.StateDataReporter(
        str(state_data_path),
        report_stride_steps,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=True,
        separator=",",
    )
    trajectory_reporter = app.DCDReporter(
        str(trajectory_path),
        report_stride_steps,
        enforcePeriodicBox=True,
    )
    screening_reporter = MetalOxygenScreeningReporter(
        report_interval=report_stride_steps,
        metal_atom_indices=metal_atom_indices,
        oxygen_atom_indices=oxygen_atom_indices,
        inner_sphere_cutoff_A=inner_sphere_cutoff_A,
    )
    simulation.reporters.append(state_data_reporter)
    simulation.reporters.append(trajectory_reporter)
    simulation.reporters.append(screening_reporter)
    simulation.step(step_count)

    final_state = simulation.context.getState(
        getPositions=True,
        getVelocities=True,
        getEnergy=True,
        enforcePeriodicBox=True,
    )
    atomic_write_text(final_state_path, XmlSerializer.serialize(final_state))
    atomic_write_text(
        final_structure_path,
        _render_pdb_text(pdb.topology, final_state.getPositions()),
    )
    return OpenMMScreeningArtifacts(
        platform_name=platform_name,
        platform_properties=platform_properties_display,
        metal_atom_count=len(metal_atom_indices),
        oxygen_atom_count=len(oxygen_atom_indices),
        frame_count=screening_reporter.frame_count,
        mean_min_metal_oxygen_distance_A=screening_reporter.mean_min_metal_oxygen_distance_A,
        inner_sphere_occupancy_fraction=screening_reporter.inner_sphere_occupancy_fraction,
    )
