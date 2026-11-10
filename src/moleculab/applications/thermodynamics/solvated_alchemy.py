# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Alchemical decoupling of a neutral ligand in an explicit water box.

This module extends the molecular free-energy path to a multi-component
System containing rigid water plus a fixed neutral ligand.  Water-water
interactions use the water box's existing nonbonded model (including Ewald
electrostatics when configured).  Ligand-water interactions are soft-core
Lennard-Jones potentials applied only to the ligand-water pairs, so the ligand
can be decoupled without scaling the solvent force field.

The ligand is held fixed at ``ligand_position`` while the water is sampled by
thermostatted Velocity Verlet.  FEP, TI, and MBAR are computed from the
ligand-water interaction energies only; the water-water energy is common to
all alchemical states and cancels from the free-energy differences.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from ...algorithms.integrators import velocity_verlet_step_system
from ...composition.system import System
from ...foundation.particles import Particle
from ...foundation.units import MOLAR_GAS_CONSTANT
from ...foundation.validation import as_positive_float, as_vector3
from ...models.nonbonded import LennardJonesParameters
from ...protocols.thermostats import berendsen_thermostat
from .fep import free_energy_perturbation
from .lj_alchemy import MultistageFEPResult
from .mbar import MBARResult
from .mbar import mbar as mbar_solver
from .mc_alchemy import SoftCoreLennardJonesPair
from .ti import TIResult, integrate_thermodynamic_derivative


@dataclass
class SolvatedLigandAlchemy:
    """Soft-core decoupling of a fixed neutral ligand in a water box."""

    water_system: System
    ligand_epsilon: float = 1.0
    ligand_sigma: float = 0.3
    ligand_mass: float = 1.0
    alpha: float = 0.5
    temperature: float = 300.0
    ligand_position: np.ndarray | None = None
    dt: float = 5e-5
    coupling_time: float = 0.05
    equilibration_steps: int = 200
    sample_every: int = 10

    def __post_init__(self) -> None:
        if not isinstance(self.water_system, System):
            raise TypeError("water_system must be a System.")
        if self.water_system.box is None:
            raise ValueError("water_system must have a periodic Box.")
        if self.water_system.nonbonded_model is None:
            raise ValueError("water_system must have a nonbonded_model.")
        object.__setattr__(
            self,
            "ligand_epsilon",
            as_positive_float(self.ligand_epsilon, "ligand_epsilon"),
        )
        object.__setattr__(
            self,
            "ligand_sigma",
            as_positive_float(self.ligand_sigma, "ligand_sigma"),
        )
        object.__setattr__(
            self,
            "ligand_mass",
            as_positive_float(self.ligand_mass, "ligand_mass"),
        )
        object.__setattr__(self, "alpha", as_positive_float(self.alpha, "alpha"))
        object.__setattr__(
            self,
            "temperature",
            as_positive_float(self.temperature, "temperature"),
        )
        object.__setattr__(self, "dt", as_positive_float(self.dt, "dt"))
        object.__setattr__(
            self,
            "coupling_time",
            as_positive_float(self.coupling_time, "coupling_time"),
        )
        for name in ("equilibration_steps", "sample_every"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
                raise ValueError(f"{name} must be an integer.")
            if int(value) < 1:
                raise ValueError(f"{name} must be at least 1.")
            object.__setattr__(self, name, int(value))
        if self.ligand_position is None:
            position = self.water_system.box.lengths / 2.0
        else:
            position = as_vector3(self.ligand_position, "ligand_position")
        object.__setattr__(self, "ligand_position", position.copy())

    @property
    def n_water(self) -> int:
        return self.water_system.n_particles

    @property
    def ligand_index(self) -> int:
        return self.n_water

    @property
    def rt(self) -> float:
        return MOLAR_GAS_CONSTANT * self.temperature

    def pair_potential(self, lambda_value: float) -> SoftCoreLennardJonesPair:
        return SoftCoreLennardJonesPair(
            epsilon=self.ligand_epsilon,
            sigma=self.ligand_sigma,
            alpha=self.alpha,
            lambda_value=lambda_value,
        )

    def make_system(self, lambda_value: float) -> System:
        """Build the water + fixed ligand System at one lambda value."""

        base_model = self.water_system.nonbonded_model
        parameters = dict(base_model.parameters)
        parameters["L"] = LennardJonesParameters(epsilon=0.0, sigma=0.1)
        ligand_water_pairs = {
            (min(self.ligand_index, index), max(self.ligand_index, index))
            for index in range(self.n_water)
        }
        model = replace(
            base_model,
            parameters=parameters,
            exclusions=frozenset(base_model.exclusions | ligand_water_pairs),
        )
        particles = [
            particle.copy() for particle in self.water_system.particles
        ]
        particles.append(
            Particle(
                mass=self.ligand_mass,
                position=self.ligand_position,
                atom_type="L",
                charge=0.0,
            )
        )
        pair = self.pair_potential(lambda_value)
        indexed = tuple(
            (self.ligand_index, index, pair) for index in range(self.n_water)
        )
        return System(
            particles=particles,
            pair_potentials=self.water_system.pair_potentials,
            topology=self.water_system.topology,
            box=self.water_system.box,
            nonbonded_model=model,
            constraints=self.water_system.constraints,
            indexed_pair_potentials=indexed,
        )

    def _ligand_water_energies(
        self,
        lambda_value: float,
        water_positions: np.ndarray,
    ) -> np.ndarray:
        positions = np.asarray(water_positions, dtype=float)
        if positions.ndim != 3 or positions.shape[1:] != (self.n_water, 3):
            raise ValueError(
                "water_positions must have shape "
                f"(n_samples, {self.n_water}, 3), got {positions.shape}."
            )
        distances = np.linalg.norm(
            positions - self.ligand_position[np.newaxis, np.newaxis, :],
            axis=2,
        )
        energies = self.pair_potential(lambda_value).energy_from_radius(distances)
        return np.sum(energies, axis=1)

    def _ligand_water_derivatives(
        self,
        lambda_value: float,
        water_positions: np.ndarray,
    ) -> np.ndarray:
        positions = np.asarray(water_positions, dtype=float)
        if positions.ndim != 3 or positions.shape[1:] != (self.n_water, 3):
            raise ValueError(
                "water_positions must have shape "
                f"(n_samples, {self.n_water}, 3), got {positions.shape}."
            )
        distances = np.linalg.norm(
            positions - self.ligand_position[np.newaxis, np.newaxis, :],
            axis=2,
        )
        values = self.pair_potential(lambda_value).derivative_from_radius(
            distances
        )
        return np.sum(values, axis=1)

    def sample_water_positions(
        self,
        lambda_value: float,
        n_samples: int,
        seed: int | None = None,
    ) -> np.ndarray:
        """Run thermostatted MD and return saved water configurations."""

        if isinstance(n_samples, bool) or not isinstance(
            n_samples, (int, np.integer)
        ):
            raise ValueError("n_samples must be an integer.")
        if int(n_samples) < 1:
            raise ValueError("n_samples must be at least 1.")
        n_samples = int(n_samples)

        system = self.make_system(lambda_value)
        rng = np.random.default_rng(seed)
        for particle in system.particles[: self.n_water]:
            particle.velocity = rng.normal(
                scale=np.sqrt(self.rt / particle.mass),
                size=3,
            )
        system.particles[self.ligand_index].velocity = np.zeros(3)

        total_steps = self.equilibration_steps + n_samples * self.sample_every
        samples: list[np.ndarray] = []
        for step in range(total_steps):
            system = velocity_verlet_step_system(system, self.dt)
            system.particles[self.ligand_index].position = (
                self.ligand_position.copy()
            )
            system.particles[self.ligand_index].velocity = np.zeros(3)
            system = berendsen_thermostat(
                system,
                temperature=self.temperature,
                coupling_time=self.coupling_time,
                dt=self.dt,
            )
            system.particles[self.ligand_index].position = (
                self.ligand_position.copy()
            )
            system.particles[self.ligand_index].velocity = np.zeros(3)
            if (
                step >= self.equilibration_steps
                and (step - self.equilibration_steps + 1) % self.sample_every == 0
            ):
                samples.append(system.positions[: self.n_water].copy())

        return np.asarray(samples, dtype=float)

    def fep(
        self,
        n_windows: int = 7,
        n_samples: int = 100,
        seed: int | None = None,
    ) -> MultistageFEPResult:
        if int(n_windows) < 2:
            raise ValueError("n_windows must be at least 2.")
        lambda_values = np.linspace(0.0, 1.0, int(n_windows))
        stage_delta_g = np.empty(len(lambda_values) - 1, dtype=float)
        stage_errors = np.empty(len(lambda_values) - 1, dtype=float)

        for index in range(len(lambda_values) - 1):
            positions = self.sample_water_positions(
                lambda_values[index],
                n_samples,
                seed=None if seed is None else seed + index,
            )
            work = self._ligand_water_energies(
                lambda_values[index + 1], positions
            ) - self._ligand_water_energies(
                lambda_values[index], positions
            )
            result = free_energy_perturbation(work, self.temperature)
            stage_delta_g[index] = result.delta_g
            stage_errors[index] = result.standard_error

        return MultistageFEPResult(
            delta_g=float(np.sum(stage_delta_g)),
            standard_error=float(np.sqrt(np.sum(stage_errors**2))),
            lambda_values=lambda_values,
            stage_delta_g=stage_delta_g,
            stage_standard_errors=stage_errors,
        )

    def ti(
        self,
        n_windows: int = 7,
        n_samples: int = 100,
        seed: int | None = None,
        method: str = "simpson",
    ) -> TIResult:
        if int(n_windows) < 2:
            raise ValueError("n_windows must be at least 2.")
        if method == "simpson" and int(n_windows) % 2 == 0:
            raise ValueError("Simpson's rule requires an odd n_windows.")
        lambda_values = np.linspace(0.0, 1.0, int(n_windows))
        derivatives = np.empty(len(lambda_values), dtype=float)
        errors = np.empty(len(lambda_values), dtype=float)

        for index, lambda_value in enumerate(lambda_values):
            positions = self.sample_water_positions(
                lambda_value,
                n_samples,
                seed=None if seed is None else seed + index,
            )
            values = self._ligand_water_derivatives(lambda_value, positions)
            derivatives[index] = float(np.mean(values))
            errors[index] = float(
                np.std(values, ddof=1) / np.sqrt(len(values))
            )

        return integrate_thermodynamic_derivative(
            lambda_values,
            derivatives,
            errors,
            method=method,
        )

    def mbar(
        self,
        n_windows: int = 7,
        n_samples: int = 100,
        seed: int | None = None,
    ) -> MBARResult:
        if int(n_windows) < 2:
            raise ValueError("n_windows must be at least 2.")
        lambda_values = np.linspace(0.0, 1.0, int(n_windows))
        all_positions = []
        sample_states = []
        for index, lambda_value in enumerate(lambda_values):
            positions = self.sample_water_positions(
                lambda_value,
                n_samples,
                seed=None if seed is None else seed + index,
            )
            all_positions.append(positions)
            sample_states.extend([index] * int(n_samples))

        positions = np.concatenate(all_positions, axis=0)
        reduced_potentials = np.empty(
            (positions.shape[0], len(lambda_values)),
            dtype=float,
        )
        for index, lambda_value in enumerate(lambda_values):
            reduced_potentials[:, index] = (
                self._ligand_water_energies(lambda_value, positions) / self.rt
            )

        return mbar_solver(
            reduced_potentials,
            np.asarray(sample_states, dtype=int),
            self.temperature,
        )
