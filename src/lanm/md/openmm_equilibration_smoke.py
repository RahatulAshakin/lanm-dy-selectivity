"""Deterministic OpenMM minimization/NVT/NPT smoke helpers for Phase 6B1."""

from __future__ import annotations

import os
from dataclasses import dataclass
from io import StringIO
from typing import Sequence

DEFAULT_NVT_DURATION_PS = 20.0
DEFAULT_NPT_DURATION_PS = 20.0
DEFAULT_PRESSURE_ATM = 1.0
DEFAULT_BAROSTAT_FREQUENCY_STEPS = 25
DEFAULT_NPT_BOX_PADDING_NM = 2.0
DEFAULT_MINIMIZATION_MAX_ITERATIONS_SCHEDULE = (200, 1000, 0)
DEFAULT_CPU_THREAD_COUNT = max(1, min(16, os.cpu_count() or 1))


@dataclass(frozen=True, slots=True)
class EquilibrationLogRow:
    stage_name: str
    status: str
    attempt_index: int
    minimization_max_iterations: int
    step_count: int
    duration_ps: float
    potential_energy_kj_per_mol: float | None
    temperature_K: float | None
    box_volume_nm3: float | None
    notes: str


@dataclass(frozen=True, slots=True)
class OpenMMEquilibrationSmokeArtifacts:
    minimized_structure_pdb_text: str
    nvt_final_pdb_text: str
    npt_final_pdb_text: str
    final_state_xml: str
    log_rows: tuple[EquilibrationLogRow, ...]
    minimized_potential_energy: float
    final_nvt_potential_energy: float
    final_npt_potential_energy: float
    final_temperature_K: float
    final_box_volume_nm3: float | None
    platform_name: str
    platform_properties: str
    minimization_max_iterations_used: int


def _sanitize_pdb_text(pdb_text: str) -> str:
    lines = [
        line
        for line in pdb_text.splitlines()
        if not line.startswith("REMARK   1 CREATED WITH OPENMM")
    ]
    return "\n".join(lines) + ("\n" if lines else "")


def _choose_openmm_platform(cpu_thread_count: int) -> tuple[object, dict[str, str], str, str]:
    try:
        from openmm import Platform
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for Phase 6B1 equilibration smoke tests") from exc

    try:
        platform = Platform.getPlatformByName("CPU")
        threads = str(max(1, cpu_thread_count))
        properties = {"Threads": threads}
        return platform, properties, "CPU", f"Threads={threads}"
    except Exception:
        try:
            platform = Platform.getPlatformByName("Reference")
            return platform, {}, "Reference", "-"
        except Exception as exc:  # pragma: no cover - exercised through runtime CLI
            raise RuntimeError("No usable OpenMM platform found for Phase 6B1") from exc


def _render_pdb_text(topology: object, positions: object) -> str:
    try:
        from openmm import app
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for Phase 6B1 equilibration smoke tests") from exc

    handle = StringIO()
    app.PDBFile.writeFile(topology, positions, handle, keepIds=True)
    return _sanitize_pdb_text(handle.getvalue())


def _compute_degrees_of_freedom(system: object) -> int:
    try:
        from openmm import CMMotionRemover, unit
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for Phase 6B1 equilibration smoke tests") from exc

    dof = 0
    zero_mass = 0.0 * unit.dalton
    for particle_index in range(system.getNumParticles()):
        if system.getParticleMass(particle_index) > zero_mass:
            dof += 3
    for constraint_index in range(system.getNumConstraints()):
        particle_a, particle_b, _ = system.getConstraintParameters(constraint_index)
        if (
            system.getParticleMass(particle_a) > zero_mass
            and system.getParticleMass(particle_b) > zero_mass
        ):
            dof -= 1
    for force_index in range(system.getNumForces()):
        if isinstance(system.getForce(force_index), CMMotionRemover):
            dof -= 3
            break
    if dof <= 0:
        raise ValueError("OpenMM system has non-positive degrees of freedom")
    return dof


def _compute_temperature_K(state: object, degrees_of_freedom: int) -> float:
    try:
        from openmm import unit
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for Phase 6B1 equilibration smoke tests") from exc

    return float(
        (2.0 * state.getKineticEnergy() / (degrees_of_freedom * unit.MOLAR_GAS_CONSTANT_R)).value_in_unit(
            unit.kelvin
        )
    )


def _compute_box_volume_nm3(state: object) -> float | None:
    try:
        from openmm import unit
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for Phase 6B1 equilibration smoke tests") from exc

    box_vectors = state.getPeriodicBoxVectors()
    a = tuple(float(value) for value in box_vectors[0].value_in_unit(unit.nanometer))
    b = tuple(float(value) for value in box_vectors[1].value_in_unit(unit.nanometer))
    c = tuple(float(value) for value in box_vectors[2].value_in_unit(unit.nanometer))
    volume = (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )
    absolute_volume = abs(volume)
    if absolute_volume == 0.0:
        return None
    return absolute_volume


def _orthorhombic_box_lengths_nm(
    positions: object,
    *,
    cutoff_nm: float,
    padding_nm: float,
) -> tuple[float, float, float]:
    try:
        from openmm import unit
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for Phase 6B1 equilibration smoke tests") from exc

    coords = positions.value_in_unit(unit.nanometer)
    minima = [min(coord[axis] for coord in coords) for axis in range(3)]
    maxima = [max(coord[axis] for coord in coords) for axis in range(3)]
    minimum_length_nm = max(2.0 * cutoff_nm + 0.2, 2.2)
    lengths = [
        max(float(maxima[axis] - minima[axis]) + padding_nm, minimum_length_nm)
        for axis in range(3)
    ]
    return (lengths[0], lengths[1], lengths[2])


def _center_positions_in_box(positions: object, box_lengths_nm: Sequence[float]) -> object:
    try:
        from openmm import Vec3, unit
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for Phase 6B1 equilibration smoke tests") from exc

    coords = positions.value_in_unit(unit.nanometer)
    minima = [min(coord[axis] for coord in coords) for axis in range(3)]
    maxima = [max(coord[axis] for coord in coords) for axis in range(3)]
    centers = [(minima[axis] + maxima[axis]) / 2.0 for axis in range(3)]
    target_centers = [float(length) / 2.0 for length in box_lengths_nm]
    shifted = [
        Vec3(
            float(coord[0] - centers[0] + target_centers[0]),
            float(coord[1] - centers[1] + target_centers[1]),
            float(coord[2] - centers[2] + target_centers[2]),
        )
        for coord in coords
    ]
    return shifted * unit.nanometer


def _orthorhombic_box_vectors(box_lengths_nm: Sequence[float]) -> tuple[object, object, object]:
    try:
        from openmm import Vec3, unit
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for Phase 6B1 equilibration smoke tests") from exc

    return (
        Vec3(float(box_lengths_nm[0]), 0.0, 0.0) * unit.nanometer,
        Vec3(0.0, float(box_lengths_nm[1]), 0.0) * unit.nanometer,
        Vec3(0.0, 0.0, float(box_lengths_nm[2])) * unit.nanometer,
    )


def _prepare_npt_system(
    *,
    system_xml: str,
    starting_positions: object,
    target_temperature_K: int,
    pressure_atm: float,
    cutoff_nm: float,
    box_padding_nm: float,
    barostat_frequency_steps: int,
) -> tuple[object, object, tuple[float, float, float], str]:
    try:
        from openmm import MonteCarloBarostat, XmlSerializer, unit
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for Phase 6B1 equilibration smoke tests") from exc

    npt_system = XmlSerializer.deserialize(system_xml)
    adjusted_nonbonded_force = False
    for force_index in range(npt_system.getNumForces()):
        force = npt_system.getForce(force_index)
        if force.__class__.__name__ != "NonbondedForce":
            continue
        if force.getNonbondedMethod() in (force.NoCutoff, force.CutoffNonPeriodic):
            force.setNonbondedMethod(force.CutoffPeriodic)
            adjusted_nonbonded_force = True
        force.setCutoffDistance(cutoff_nm * unit.nanometer)

    box_lengths_nm = _orthorhombic_box_lengths_nm(
        starting_positions,
        cutoff_nm=cutoff_nm,
        padding_nm=box_padding_nm,
    )
    box_vectors = _orthorhombic_box_vectors(box_lengths_nm)
    npt_system.setDefaultPeriodicBoxVectors(*box_vectors)
    barostat = MonteCarloBarostat(
        pressure_atm * unit.atmosphere,
        target_temperature_K * unit.kelvin,
        barostat_frequency_steps,
    )
    barostat.setRandomNumberSeed(0)
    npt_system.addForce(barostat)
    centered_positions = _center_positions_in_box(starting_positions, box_lengths_nm)
    notes = (
        "periodicized_npt_box_nm="
        f"{box_lengths_nm[0]:.3f}x{box_lengths_nm[1]:.3f}x{box_lengths_nm[2]:.3f}; "
        f"nonbonded_periodicized={adjusted_nonbonded_force}"
    )
    return npt_system, centered_positions, box_lengths_nm, notes


def run_openmm_equilibration_smoke_protocol(
    *,
    prepared_structure_pdb_text: str,
    system_xml: str,
    target_temperature_K: int,
    friction_coeff_ps: float,
    timestep_fs: float,
    cutoff_nm: float,
    nvt_duration_ps: float = DEFAULT_NVT_DURATION_PS,
    npt_duration_ps: float = DEFAULT_NPT_DURATION_PS,
    pressure_atm: float = DEFAULT_PRESSURE_ATM,
    barostat_frequency_steps: int = DEFAULT_BAROSTAT_FREQUENCY_STEPS,
    box_padding_nm: float = DEFAULT_NPT_BOX_PADDING_NM,
    minimization_max_iterations_schedule: Sequence[int] = DEFAULT_MINIMIZATION_MAX_ITERATIONS_SCHEDULE,
    cpu_thread_count: int = DEFAULT_CPU_THREAD_COUNT,
) -> OpenMMEquilibrationSmokeArtifacts:
    """Run deterministic minimization, NVT, and NPT smoke stages for one system."""
    try:
        from openmm import LangevinMiddleIntegrator, XmlSerializer, app, unit
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised through runtime CLI
        raise RuntimeError("OpenMM is required for Phase 6B1 equilibration smoke tests") from exc

    if target_temperature_K <= 0:
        raise ValueError("target_temperature_K must be positive")
    if friction_coeff_ps <= 0.0:
        raise ValueError("friction_coeff_ps must be positive")
    if timestep_fs <= 0.0:
        raise ValueError("timestep_fs must be positive")
    if cutoff_nm <= 0.0:
        raise ValueError("cutoff_nm must be positive")
    if nvt_duration_ps <= 0.0 or npt_duration_ps <= 0.0:
        raise ValueError("nvt_duration_ps and npt_duration_ps must be positive")
    if pressure_atm <= 0.0:
        raise ValueError("pressure_atm must be positive")
    if box_padding_nm <= 0.0:
        raise ValueError("box_padding_nm must be positive")
    if barostat_frequency_steps <= 0:
        raise ValueError("barostat_frequency_steps must be positive")

    nvt_step_count = int(round((nvt_duration_ps * 1000.0) / timestep_fs))
    npt_step_count = int(round((npt_duration_ps * 1000.0) / timestep_fs))
    if nvt_step_count <= 0 or npt_step_count <= 0:
        raise ValueError("NVT and NPT smoke stages must contain at least one step")

    platform, platform_properties, platform_name, platform_properties_display = _choose_openmm_platform(
        cpu_thread_count
    )
    pdb = app.PDBFile(StringIO(prepared_structure_pdb_text))
    last_exception: Exception | None = None

    for attempt_index, minimization_max_iterations in enumerate(
        tuple(minimization_max_iterations_schedule),
        start=1,
    ):
        try:
            nvt_system = XmlSerializer.deserialize(system_xml)
            nvt_integrator = LangevinMiddleIntegrator(
                target_temperature_K * unit.kelvin,
                friction_coeff_ps / unit.picosecond,
                timestep_fs * unit.femtoseconds,
            )
            nvt_integrator.setRandomNumberSeed(0)
            nvt_simulation = app.Simulation(
                pdb.topology,
                nvt_system,
                nvt_integrator,
                platform,
                platform_properties,
            )
            nvt_simulation.context.setPositions(pdb.positions)
            if minimization_max_iterations > 0:
                nvt_simulation.minimizeEnergy(maxIterations=minimization_max_iterations)
            else:
                nvt_simulation.minimizeEnergy()
            minimized_state = nvt_simulation.context.getState(getPositions=True, getEnergy=True)
            minimized_pdb_text = _render_pdb_text(pdb.topology, minimized_state.getPositions())
            minimized_potential_energy = float(
                minimized_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            )
            nvt_degrees_of_freedom = _compute_degrees_of_freedom(nvt_system)
            nvt_simulation.context.setVelocitiesToTemperature(target_temperature_K * unit.kelvin, 0)
            nvt_simulation.step(nvt_step_count)
            nvt_state = nvt_simulation.context.getState(
                getPositions=True,
                getVelocities=True,
                getEnergy=True,
            )
            nvt_pdb_text = _render_pdb_text(pdb.topology, nvt_state.getPositions())
            final_nvt_potential_energy = float(
                nvt_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            )
            nvt_temperature_K = _compute_temperature_K(nvt_state, nvt_degrees_of_freedom)

            npt_system, npt_starting_positions, _, npt_notes = _prepare_npt_system(
                system_xml=system_xml,
                starting_positions=nvt_state.getPositions(),
                target_temperature_K=target_temperature_K,
                pressure_atm=pressure_atm,
                cutoff_nm=cutoff_nm,
                box_padding_nm=box_padding_nm,
                barostat_frequency_steps=barostat_frequency_steps,
            )
            npt_degrees_of_freedom = _compute_degrees_of_freedom(npt_system)
            npt_integrator = LangevinMiddleIntegrator(
                target_temperature_K * unit.kelvin,
                friction_coeff_ps / unit.picosecond,
                timestep_fs * unit.femtoseconds,
            )
            npt_integrator.setRandomNumberSeed(0)
            npt_simulation = app.Simulation(
                pdb.topology,
                npt_system,
                npt_integrator,
                platform,
                platform_properties,
            )
            npt_box_vectors = _orthorhombic_box_vectors(
                _orthorhombic_box_lengths_nm(
                    nvt_state.getPositions(),
                    cutoff_nm=cutoff_nm,
                    padding_nm=box_padding_nm,
                )
            )
            npt_simulation.context.setPeriodicBoxVectors(*npt_box_vectors)
            npt_simulation.context.setPositions(npt_starting_positions)
            npt_simulation.context.setVelocities(nvt_state.getVelocities())
            npt_simulation.step(npt_step_count)
            npt_state = npt_simulation.context.getState(
                getPositions=True,
                getVelocities=True,
                getEnergy=True,
                enforcePeriodicBox=True,
            )
            npt_pdb_text = _render_pdb_text(pdb.topology, npt_state.getPositions())
            final_npt_potential_energy = float(
                npt_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            )
            final_temperature_K = _compute_temperature_K(npt_state, npt_degrees_of_freedom)
            final_box_volume_nm3 = _compute_box_volume_nm3(npt_state)
            final_state = npt_simulation.context.getState(
                getPositions=True,
                getVelocities=True,
                getEnergy=True,
                enforcePeriodicBox=True,
            )
            final_state_xml = XmlSerializer.serialize(final_state)
            attempt_notes = (
                f"platform={platform_name}; properties={platform_properties_display}; "
                f"temperature_target_K={target_temperature_K}; pressure_atm={pressure_atm}; "
                f"{npt_notes}"
            )
            log_rows = (
                EquilibrationLogRow(
                    stage_name="minimization",
                    status="success",
                    attempt_index=attempt_index,
                    minimization_max_iterations=minimization_max_iterations,
                    step_count=0,
                    duration_ps=0.0,
                    potential_energy_kj_per_mol=minimized_potential_energy,
                    temperature_K=None,
                    box_volume_nm3=None,
                    notes=attempt_notes,
                ),
                EquilibrationLogRow(
                    stage_name="nvt",
                    status="success",
                    attempt_index=attempt_index,
                    minimization_max_iterations=minimization_max_iterations,
                    step_count=nvt_step_count,
                    duration_ps=nvt_duration_ps,
                    potential_energy_kj_per_mol=final_nvt_potential_energy,
                    temperature_K=nvt_temperature_K,
                    box_volume_nm3=None,
                    notes=attempt_notes,
                ),
                EquilibrationLogRow(
                    stage_name="npt",
                    status="success",
                    attempt_index=attempt_index,
                    minimization_max_iterations=minimization_max_iterations,
                    step_count=npt_step_count,
                    duration_ps=npt_duration_ps,
                    potential_energy_kj_per_mol=final_npt_potential_energy,
                    temperature_K=final_temperature_K,
                    box_volume_nm3=final_box_volume_nm3,
                    notes=attempt_notes,
                ),
            )
            return OpenMMEquilibrationSmokeArtifacts(
                minimized_structure_pdb_text=minimized_pdb_text,
                nvt_final_pdb_text=nvt_pdb_text,
                npt_final_pdb_text=npt_pdb_text,
                final_state_xml=final_state_xml,
                log_rows=log_rows,
                minimized_potential_energy=minimized_potential_energy,
                final_nvt_potential_energy=final_nvt_potential_energy,
                final_npt_potential_energy=final_npt_potential_energy,
                final_temperature_K=final_temperature_K,
                final_box_volume_nm3=final_box_volume_nm3,
                platform_name=platform_name,
                platform_properties=platform_properties_display,
                minimization_max_iterations_used=minimization_max_iterations,
            )
        except Exception as exc:
            last_exception = exc

    if last_exception is None:  # pragma: no cover - defensive
        raise RuntimeError("Phase 6B1 equilibration smoke terminated without producing a result")
    raise last_exception
