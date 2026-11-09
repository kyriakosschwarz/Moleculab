"""Multi-particle Monte Carlo free-energy path for a bonded LJ pair.

This module samples a real two-particle :class:`~moleculab.system.System` by
Metropolis Monte Carlo.  The particles are connected by a harmonic bond and
also interact through a soft-core Lennard-Jones pair potential.  The
alchemical coupling is

.. math::

    U_\\lambda = U_{\\mathrm{bond}} + U_{\\mathrm{LJ}}(\\lambda),

with the same soft-core form used by :mod:`moleculab.lj_alchemy`.

Unlike the radial toy model, sampling here moves actual particles in a
multi-particle System and accepts or rejects moves using
``System.potential_energy()``.  The reference is still an independent
high-accuracy radial integral because the bond and pair potentials depend only
on the relative distance.

FEP, TI, and MBAR are then computed from the sampled relative coordinates.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._validation import as_finite_float, as_positive_float
from .bonded import HarmonicBondPotential
from .fep import free_energy_perturbation
from .lj_alchemy import MultistageFEPResult
from .mbar import MBARResult, mbar as mbar_solver
from .particles import Particle
from .system import System
from .ti import TIResult, integrate_thermodynamic_derivative
from .topology import Bond, Topology
from .units import BOLTZMANN_CONSTANT


@dataclass(frozen=True)
class SoftCoreLennardJonesPair:
    """A soft-core Lennard-Jones pair potential."""

    epsilon: float
    sigma: float
    alpha: float
    lambda_value: float

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "epsilon", as_positive_float(self.epsilon, "epsilon")
        )
        object.__setattr__(self, "sigma", as_positive_float(self.sigma, "sigma"))
        object.__setattr__(self, "alpha", as_positive_float(self.alpha, "alpha"))
        lambda_value = as_finite_float(self.lambda_value, "lambda_value")
        if not 0.0 <= lambda_value <= 1.0:
            raise ValueError("lambda_value must be between 0 and 1.")
        object.__setattr__(self, "lambda_value", lambda_value)

    def energy_from_radius(self, radius: np.ndarray | float) -> np.ndarray | float:
        """Return the pair energy in kJ/mol for a radial distance."""

        r = np.asarray(radius, dtype=float)
        if np.any(r <= 0.0) or not np.all(np.isfinite(r)):
            raise ValueError("radius must contain finite positive values.")
        s = self.alpha * (1.0 - self.lambda_value) + (r / self.sigma) ** 6
        return 4.0 * self.epsilon * self.lambda_value * (1.0 / s**2 - 1.0 / s)

    def derivative_from_radius(self, radius: np.ndarray | float) -> np.ndarray | float:
        """Return dU/dlambda from a radial distance."""

        r = np.asarray(radius, dtype=float)
        if np.any(r <= 0.0) or not np.all(np.isfinite(r)):
            raise ValueError("radius must contain finite positive values.")
        s = self.alpha * (1.0 - self.lambda_value) + (r / self.sigma) ** 6
        return 4.0 * self.epsilon * (
            (1.0 / s**2 - 1.0 / s)
            + self.lambda_value
            * (2.0 * self.alpha / s**3 - self.alpha / s**2)
        )

    def energy(self, displacement: np.ndarray) -> float:
        radius = float(np.linalg.norm(displacement))
        return float(self.energy_from_radius(radius))

    def force(self, displacement: np.ndarray) -> np.ndarray:
        vector = np.asarray(displacement, dtype=float)
        radius = float(np.linalg.norm(vector))
        if radius == 0.0:
            raise ValueError("Pair potential is undefined at zero separation.")
        s = self.alpha * (1.0 - self.lambda_value) + (radius / self.sigma) ** 6
        ds_dr = 6.0 * radius**5 / self.sigma**6
        dU_dr = (
            4.0
            * self.epsilon
            * self.lambda_value
            * (-2.0 / s**3 + 1.0 / s**2)
            * ds_dr
        )
        return -dU_dr * vector / radius


@dataclass(frozen=True)
class BondedLennardJonesMC:
    """Monte Carlo sampler and free-energy estimators for a bonded LJ pair.

    The two particles have equal mass and are connected by a harmonic bond.
    A soft-core Lennard-Jones pair potential is coupled through ``lambda``.
    """

    epsilon: float = 1.0
    sigma: float = 0.3
    bond_stiffness: float = 1000.0
    bond_length: float = 0.4
    temperature: float = 300.0
    alpha: float = 0.5
    mass: float = 1.0
    max_displacement: float = 0.015
    equilibration_steps: int = 1000
    sample_every: int = 5
    reference_r_max: float = 2.0
    reference_grid: int = 100_000

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "epsilon", as_positive_float(self.epsilon, "epsilon")
        )
        object.__setattr__(self, "sigma", as_positive_float(self.sigma, "sigma"))
        object.__setattr__(
            self,
            "bond_stiffness",
            as_positive_float(self.bond_stiffness, "bond_stiffness"),
        )
        object.__setattr__(
            self,
            "bond_length",
            as_positive_float(self.bond_length, "bond_length"),
        )
        object.__setattr__(
            self, "temperature", as_positive_float(self.temperature, "temperature")
        )
        object.__setattr__(self, "alpha", as_positive_float(self.alpha, "alpha"))
        object.__setattr__(self, "mass", as_positive_float(self.mass, "mass"))
        object.__setattr__(
            self,
            "max_displacement",
            as_positive_float(self.max_displacement, "max_displacement"),
        )
        for name in ("equilibration_steps", "sample_every", "reference_grid"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
                raise ValueError(f"{name} must be an integer.")
            if int(value) < 1:
                raise ValueError(f"{name} must be at least 1.")
            object.__setattr__(self, name, int(value))
        object.__setattr__(
            self,
            "reference_r_max",
            as_positive_float(self.reference_r_max, "reference_r_max"),
        )

    @property
    def rt(self) -> float:
        return BOLTZMANN_CONSTANT * self.temperature

    def pair_potential(self, lambda_value: float) -> SoftCoreLennardJonesPair:
        return SoftCoreLennardJonesPair(
            epsilon=self.epsilon,
            sigma=self.sigma,
            alpha=self.alpha,
            lambda_value=lambda_value,
        )

    def make_system(self, lambda_value: float) -> System:
        """Build the two-particle System at one alchemical state."""

        half = 0.5 * self.bond_length
        return System(
            particles=[
                Particle(self.mass, [-half, 0.0, 0.0]),
                Particle(self.mass, [half, 0.0, 0.0]),
            ],
            pair_potentials=[self.pair_potential(lambda_value)],
            topology=Topology(
                bonds=[
                    Bond(
                        0,
                        1,
                        HarmonicBondPotential(
                            stiffness=self.bond_stiffness,
                            length0=self.bond_length,
                        ),
                    )
                ]
            ),
        )

    def _energy_at_radius(self, lambda_value: float, radius: np.ndarray) -> np.ndarray:
        bond = 0.5 * self.bond_stiffness * (radius - self.bond_length) ** 2
        return bond + self.pair_potential(lambda_value).energy_from_radius(radius)

    def sample_positions(
        self,
        lambda_value: float,
        n_samples: int,
        seed: int | None = None,
    ) -> np.ndarray:
        """Run Metropolis MC and return saved two-particle configurations."""

        if isinstance(n_samples, bool) or not isinstance(
            n_samples, (int, np.integer)
        ):
            raise ValueError("n_samples must be an integer.")
        if int(n_samples) < 1:
            raise ValueError("n_samples must be at least 1.")
        n_samples = int(n_samples)

        system = self.make_system(lambda_value)
        rng = np.random.default_rng(seed)
        current_energy = system.potential_energy()
        samples: list[np.ndarray] = []
        total_steps = self.equilibration_steps + n_samples * self.sample_every

        for step in range(total_steps):
            old_first = system.particles[0].position.copy()
            old_second = system.particles[1].position.copy()
            delta = rng.normal(scale=self.max_displacement, size=3)
            system.particles[0].position += delta
            system.particles[1].position -= delta
            proposed_energy = system.potential_energy()

            if proposed_energy <= current_energy or rng.random() < np.exp(
                -(proposed_energy - current_energy) / self.rt
            ):
                current_energy = proposed_energy
            else:
                system.particles[0].position = old_first
                system.particles[1].position = old_second

            if (
                step >= self.equilibration_steps
                and (step - self.equilibration_steps + 1) % self.sample_every == 0
            ):
                samples.append(system.positions.copy())

        return np.asarray(samples, dtype=float)

    def _sample_radii(
        self,
        lambda_value: float,
        n_samples: int,
        seed: int | None = None,
    ) -> np.ndarray:
        positions = self.sample_positions(lambda_value, n_samples, seed=seed)
        return np.linalg.norm(positions[:, 0] - positions[:, 1], axis=1)

    def reference_delta_g(self) -> float:
        """Return the endpoint free-energy difference by radial quadrature."""

        radius = np.linspace(0.0, self.reference_r_max, self.reference_grid)
        radius[0] = 1e-12

        def log_partition(lambda_value: float) -> float:
            log_probability = (
                2.0 * np.log(radius)
                - self._energy_at_radius(lambda_value, radius) / self.rt
            )
            maximum = float(np.max(log_probability))
            integral = np.trapezoid(
                np.exp(log_probability - maximum),
                radius,
            )
            return maximum + float(np.log(integral))

        return float(-self.rt * (log_partition(1.0) - log_partition(0.0)))

    def fep(
        self,
        n_windows: int = 11,
        n_samples: int = 2000,
        seed: int | None = None,
    ) -> MultistageFEPResult:
        """Run multistage FEP from actual Monte Carlo samples."""

        if isinstance(n_windows, bool) or not isinstance(
            n_windows, (int, np.integer)
        ):
            raise ValueError("n_windows must be an integer.")
        if int(n_windows) < 2:
            raise ValueError("n_windows must be at least 2.")

        lambda_values = np.linspace(0.0, 1.0, int(n_windows))
        stage_delta_g = np.empty(len(lambda_values) - 1, dtype=float)
        stage_errors = np.empty(len(lambda_values) - 1, dtype=float)

        for index in range(len(lambda_values) - 1):
            radii = self._sample_radii(
                lambda_values[index],
                n_samples,
                seed=None if seed is None else seed + index,
            )
            work = (
                self.pair_potential(lambda_values[index + 1]).energy_from_radius(radii)
                - self.pair_potential(lambda_values[index]).energy_from_radius(radii)
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
        n_windows: int = 21,
        n_samples: int = 2000,
        seed: int | None = None,
        method: str = "simpson",
    ) -> TIResult:
        """Run thermodynamic integration from actual Monte Carlo samples."""

        if isinstance(n_windows, bool) or not isinstance(
            n_windows, (int, np.integer)
        ):
            raise ValueError("n_windows must be an integer.")
        if int(n_windows) < 2:
            raise ValueError("n_windows must be at least 2.")
        if method == "simpson" and int(n_windows) % 2 == 0:
            raise ValueError("Simpson's rule requires an odd n_windows.")

        lambda_values = np.linspace(0.0, 1.0, int(n_windows))
        derivatives = np.empty(len(lambda_values), dtype=float)
        errors = np.empty(len(lambda_values), dtype=float)

        for index, lambda_value in enumerate(lambda_values):
            radii = self._sample_radii(
                lambda_value,
                n_samples,
                seed=None if seed is None else seed + index,
            )
            values = self.pair_potential(lambda_value).derivative_from_radius(radii)
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
        n_windows: int = 11,
        n_samples: int = 2000,
        seed: int | None = None,
    ) -> MBARResult:
        """Run MBAR from actual Monte Carlo samples."""

        if isinstance(n_windows, bool) or not isinstance(
            n_windows, (int, np.integer)
        ):
            raise ValueError("n_windows must be an integer.")
        if int(n_windows) < 2:
            raise ValueError("n_windows must be at least 2.")

        lambda_values = np.linspace(0.0, 1.0, int(n_windows))
        all_radii = []
        sample_states = []
        for index, lambda_value in enumerate(lambda_values):
            radii = self._sample_radii(
                lambda_value,
                n_samples,
                seed=None if seed is None else seed + index,
            )
            all_radii.append(radii)
            sample_states.extend([index] * int(n_samples))

        radii = np.concatenate(all_radii)
        reduced_potentials = np.empty(
            (len(radii), len(lambda_values)),
            dtype=float,
        )
        for index, lambda_value in enumerate(lambda_values):
            reduced_potentials[:, index] = (
                self._energy_at_radius(lambda_value, radii) / self.rt
            )

        return mbar_solver(
            reduced_potentials,
            np.asarray(sample_states, dtype=int),
            self.temperature,
        )
