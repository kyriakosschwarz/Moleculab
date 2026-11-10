"""Alchemical Lennard-Jones pair in a harmonic trap.

This module extends the free-energy estimators from the harmonic oscillator to
a simple molecular interaction with a repulsive core.  A mobile particle moves
relative to a fixed particle at the origin under

.. math::

    U_\\lambda(r) = \\frac{1}{2}k r^2
    + 4\\epsilon\\lambda
    \\left[
    \\frac{1}{s(r,\\lambda)^2}
    - \\frac{1}{s(r,\\lambda)}
    \\right],

where

.. math::

    s(r,\\lambda) = \\alpha(1-\\lambda) + (r/\\sigma)^6.

The first term is a harmonic trap.  The second is a **soft-core**
Lennard-Jones interaction: at :math:`\\lambda=0` it vanishes, at
:math:`\\lambda=1` it is the ordinary 12-6 Lennard-Jones potential, and the
soft core removes the singular endpoint that would make linear coupling have
infinite TI/MBAR variance.

The relative-coordinate radial distribution is

.. math::

    p_\\lambda(r) \\propto r^2 e^{-U_\\lambda(r)/(RT)}.

Samples are drawn exactly from this distribution by inverse-CDF interpolation
on a radial grid.  The reference free-energy difference is obtained by
high-accuracy numerical integration of the radial partition function.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..foundation.units import BOLTZMANN_CONSTANT
from ..foundation.validation import as_finite_float, as_positive_float
from .fep import free_energy_perturbation
from .mbar import MBARResult
from .mbar import mbar as mbar_solver
from .ti import TIResult, integrate_thermodynamic_derivative


@dataclass(frozen=True)
class MultistageFEPResult:
    """Result of multistage FEP along the alchemical path."""

    delta_g: float
    standard_error: float
    lambda_values: np.ndarray
    stage_delta_g: np.ndarray
    stage_standard_errors: np.ndarray


@dataclass(frozen=True)
class SoftCoreLennardJones:
    """A soft-core LJ pair coupled to a harmonic trap.

    Parameters
    ----------
    epsilon:
        Lennard-Jones well depth in kJ/mol.
    sigma:
        Lennard-Jones diameter in nm.
    trap_stiffness:
        Harmonic trap stiffness in kJ/(mol nm^2).
    temperature:
        Temperature in K.
    alpha:
        Soft-core parameter.  Smaller values soften the core more strongly.
    r_max:
        Radial grid maximum in nm.
    n_grid:
        Number of radial grid points used for sampling and the reference.
    """

    epsilon: float = 1.0
    sigma: float = 0.3
    trap_stiffness: float = 20.0
    temperature: float = 300.0
    alpha: float = 0.5
    r_max: float = 2.0
    n_grid: int = 100_000

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "epsilon", as_positive_float(self.epsilon, "epsilon")
        )
        object.__setattr__(self, "sigma", as_positive_float(self.sigma, "sigma"))
        object.__setattr__(
            self,
            "trap_stiffness",
            as_positive_float(self.trap_stiffness, "trap_stiffness"),
        )
        object.__setattr__(
            self,
            "temperature",
            as_positive_float(self.temperature, "temperature"),
        )
        object.__setattr__(self, "alpha", as_positive_float(self.alpha, "alpha"))
        object.__setattr__(self, "r_max", as_positive_float(self.r_max, "r_max"))
        if isinstance(self.n_grid, bool) or not isinstance(
            self.n_grid, (int, np.integer)
        ):
            raise ValueError("n_grid must be an integer.")
        if int(self.n_grid) < 1000:
            raise ValueError("n_grid must be at least 1000.")
        object.__setattr__(self, "n_grid", int(self.n_grid))

    @property
    def rt(self) -> float:
        """Molar thermal energy RT in kJ/mol."""

        return BOLTZMANN_CONSTANT * self.temperature

    def soft_core_energy(self, lambda_value: float, radius: np.ndarray) -> np.ndarray:
        """Return the soft-core LJ energy (excluding the trap) in kJ/mol."""

        lam = as_finite_float(lambda_value, "lambda_value")
        if not 0.0 <= lam <= 1.0:
            raise ValueError("lambda_value must be between 0 and 1.")
        r = np.asarray(radius, dtype=float)
        if np.any(r <= 0.0) or not np.all(np.isfinite(r)):
            raise ValueError("radius must contain finite positive values.")
        s = self.alpha * (1.0 - lam) + (r / self.sigma) ** 6
        return 4.0 * self.epsilon * lam * (1.0 / s**2 - 1.0 / s)

    def soft_core_derivative(
        self,
        lambda_value: float,
        radius: np.ndarray,
    ) -> np.ndarray:
        """Return dU_softcore/dlambda in kJ/mol."""

        lam = as_finite_float(lambda_value, "lambda_value")
        if not 0.0 <= lam <= 1.0:
            raise ValueError("lambda_value must be between 0 and 1.")
        r = np.asarray(radius, dtype=float)
        if np.any(r <= 0.0) or not np.all(np.isfinite(r)):
            raise ValueError("radius must contain finite positive values.")
        s = self.alpha * (1.0 - lam) + (r / self.sigma) ** 6
        return 4.0 * self.epsilon * (
            (1.0 / s**2 - 1.0 / s)
            + lam * (2.0 * self.alpha / s**3 - self.alpha / s**2)
        )

    def energy(self, lambda_value: float, radius: np.ndarray) -> np.ndarray:
        """Return total trap + soft-core energy in kJ/mol."""

        r = np.asarray(radius, dtype=float)
        return 0.5 * self.trap_stiffness * r * r + self.soft_core_energy(
            lambda_value, r
        )

    def _radial_grid(self) -> np.ndarray:
        grid = np.linspace(0.0, self.r_max, self.n_grid)
        grid[0] = 1e-12
        return grid

    def _radial_cdf(self, lambda_value: float) -> tuple[np.ndarray, np.ndarray]:
        radius = self._radial_grid()
        log_probability = (
            2.0 * np.log(radius)
            - self.energy(lambda_value, radius) / self.rt
        )
        probability = np.exp(log_probability - np.max(log_probability))
        cumulative = np.concatenate(
            [
                [0.0],
                np.cumsum(
                    0.5 * (probability[1:] + probability[:-1]) * np.diff(radius)
                ),
            ]
        )
        return radius, cumulative / cumulative[-1]

    def sample_radii(
        self,
        lambda_value: float,
        n_samples: int,
        seed: int | None = None,
    ) -> np.ndarray:
        """Draw radii from the lambda state's radial Boltzmann distribution."""

        if isinstance(n_samples, bool) or not isinstance(
            n_samples, (int, np.integer)
        ):
            raise ValueError("n_samples must be an integer.")
        if int(n_samples) < 1:
            raise ValueError("n_samples must be at least 1.")
        rng = np.random.default_rng(seed)
        radius, cumulative = self._radial_cdf(lambda_value)
        return np.interp(rng.random(int(n_samples)), cumulative, radius)

    def sample_positions(
        self,
        lambda_value: float,
        n_samples: int,
        seed: int | None = None,
    ) -> np.ndarray:
        """Draw 3D relative positions from the lambda state."""

        radii = self.sample_radii(lambda_value, n_samples, seed=seed)
        rng = np.random.default_rng(None if seed is None else seed + 1)
        directions = rng.normal(size=(int(n_samples), 3))
        directions /= np.linalg.norm(directions, axis=1, keepdims=True)
        return directions * radii[:, np.newaxis]

    def reference_delta_g(self) -> float:
        """Return the endpoint free-energy difference by radial quadrature."""

        radius = self._radial_grid()

        def log_partition(lambda_value: float) -> float:
            log_probability = (
                2.0 * np.log(radius)
                - self.energy(lambda_value, radius) / self.rt
            )
            maximum = float(np.max(log_probability))
            integral = np.trapezoid(
                np.exp(log_probability - maximum),
                radius,
            )
            return maximum + float(np.log(integral))

        return float(
            -self.rt * (log_partition(1.0) - log_partition(0.0))
        )

    def fep(
        self,
        n_windows: int = 11,
        n_samples: int = 50_000,
        seed: int | None = None,
    ) -> MultistageFEPResult:
        """Run multistage FEP between lambda=0 and lambda=1."""

        if isinstance(n_windows, bool) or not isinstance(
            n_windows, (int, np.integer)
        ):
            raise ValueError("n_windows must be an integer.")
        if int(n_windows) < 2:
            raise ValueError("n_windows must be at least 2.")
        if isinstance(n_samples, bool) or not isinstance(
            n_samples, (int, np.integer)
        ):
            raise ValueError("n_samples must be an integer.")
        if int(n_samples) < 2:
            raise ValueError("n_samples must be at least 2.")

        lambda_values = np.linspace(0.0, 1.0, int(n_windows))
        stage_delta_g = np.empty(len(lambda_values) - 1, dtype=float)
        stage_errors = np.empty(len(lambda_values) - 1, dtype=float)

        for index in range(len(lambda_values) - 1):
            radii = self.sample_radii(
                lambda_values[index],
                int(n_samples),
                seed=None if seed is None else seed + index,
            )
            work = (
                self.energy(lambda_values[index + 1], radii)
                - self.energy(lambda_values[index], radii)
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
        n_samples: int = 50_000,
        seed: int | None = None,
        method: str = "simpson",
    ) -> TIResult:
        """Run thermodynamic integration between lambda=0 and lambda=1."""

        if isinstance(n_windows, bool) or not isinstance(
            n_windows, (int, np.integer)
        ):
            raise ValueError("n_windows must be an integer.")
        if int(n_windows) < 2:
            raise ValueError("n_windows must be at least 2.")
        if method == "simpson" and int(n_windows) % 2 == 0:
            raise ValueError("Simpson's rule requires an odd n_windows.")
        if isinstance(n_samples, bool) or not isinstance(
            n_samples, (int, np.integer)
        ):
            raise ValueError("n_samples must be an integer.")
        if int(n_samples) < 2:
            raise ValueError("n_samples must be at least 2.")

        lambda_values = np.linspace(0.0, 1.0, int(n_windows))
        derivatives = np.empty(len(lambda_values), dtype=float)
        errors = np.empty(len(lambda_values), dtype=float)

        for index, lambda_value in enumerate(lambda_values):
            radii = self.sample_radii(
                lambda_value,
                int(n_samples),
                seed=None if seed is None else seed + index,
            )
            values = self.soft_core_derivative(lambda_value, radii)
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
        n_samples: int = 20_000,
        seed: int | None = None,
    ) -> MBARResult:
        """Run MBAR over the soft-core alchemical windows."""

        if isinstance(n_windows, bool) or not isinstance(
            n_windows, (int, np.integer)
        ):
            raise ValueError("n_windows must be an integer.")
        if int(n_windows) < 2:
            raise ValueError("n_windows must be at least 2.")
        if isinstance(n_samples, bool) or not isinstance(
            n_samples, (int, np.integer)
        ):
            raise ValueError("n_samples must be an integer.")
        if int(n_samples) < 1:
            raise ValueError("n_samples must be at least 1.")

        lambda_values = np.linspace(0.0, 1.0, int(n_windows))
        radii = []
        sample_states = []
        for index, lambda_value in enumerate(lambda_values):
            radii.append(
                self.sample_radii(
                    lambda_value,
                    int(n_samples),
                    seed=None if seed is None else seed + index,
                )
            )
            sample_states.extend([index] * int(n_samples))

        all_radii = np.concatenate(radii)
        reduced_potentials = np.empty(
            (len(all_radii), len(lambda_values)),
            dtype=float,
        )
        for index, lambda_value in enumerate(lambda_values):
            reduced_potentials[:, index] = (
                self.energy(lambda_value, all_radii) / self.rt
            )

        return mbar_solver(
            reduced_potentials,
            np.asarray(sample_states, dtype=int),
            self.temperature,
        )
