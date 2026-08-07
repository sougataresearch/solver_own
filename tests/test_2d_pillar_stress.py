"""Phase 4b — near-degenerate / ill-conditioned stress tests for
`solve_layer_eigenmodes_patterned`.

Per `phases.md` Phase 4b and `eigenmodes.py`'s docstring: a deliberate
stress sweep (index contrast from `3.48` to a lossy-metal-like `-20+2j`,
`num_orders` up to 225, near-touching pillars, sub-percent sliver
rectangles, near-degenerate nested shapes) found no catastrophic failure —
condition numbers grew but stayed in the tens-to-low-hundreds, with energy
conservation and the independent `RigorousCoupledWaveAnalysis.jl`
eigenvalue oracle (`tests/oracles/rcwa_2djl_eigenvalues.py`) both holding
to ~1e-10 throughout. This file freezes representative cases from that
sweep as permanent regression tests (an honest empirical finding, not a
promise no pathological case exists elsewhere — see the docstring's
caveat), plus tests for the condition-number `WARNING` logging itself.
"""

from __future__ import annotations

import logging

import numpy as np
import pytest
from oracles.rcwa_2djl_eigenvalues import eigenoperator_eigenvalues

from sougata_solver import eigenmodes
from sougata_solver.eigenmodes import solve_layer_eigenmodes_patterned
from sougata_solver.fourier_basis import truncate_fourier_orders
from sougata_solver.fourier_factorization import toeplitz_matrix
from sougata_solver.geometry import Circle, Lattice, Pattern, Rectangle
from sougata_solver.materials import Material


pytestmark = pytest.mark.oracle  # Category 17 target 17.1: system-tier test, cross-checked against a named external oracle

PERIOD = 0.7
AIR = Material("air", 1.0)


def _kx_ky(lattice: Lattice, g: np.ndarray, omega: float) -> tuple[np.ndarray, np.ndarray]:
    lk = lattice.reciprocal_vectors()
    kx = 2 * np.pi * (g[:, 0] * lk[0, 0] + g[:, 1] * lk[1, 0])
    ky = 2 * np.pi * (g[:, 0] * lk[0, 1] + g[:, 1] * lk[1, 1])
    return kx, ky


def _check_case(pattern: Pattern, num_orders: int, cond_ceiling: float) -> None:
    """Shared assertion body: energy conservation via the RCWA.jl
    eigenvalue oracle agreement (a proxy for R/T correctness at the
    eigenoperator level, cheaper than a full `Simulation.solve()` per
    case) plus an explicit condition-number sanity ceiling, so a future
    regression that pushes conditioning far past what was observed this
    session shows up as a test failure, not just a quiet WARNING."""
    wavelength = 1.0
    omega = 2 * np.pi / wavelength
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    g = truncate_fourier_orders(lattice, num_orders, "circular")
    kx, ky = _kx_ky(lattice, g, omega)
    epsilon_hat = toeplitz_matrix(pattern, lattice, g, wavelength, inverse=False)

    modes = solve_layer_eigenmodes_patterned(omega, kx, ky, epsilon_hat)
    ours = np.sort((modes.q**2).real)
    oracle = np.sort(eigenoperator_eigenvalues(omega, kx, ky, epsilon_hat).real)

    assert ours == pytest.approx(oracle, abs=1e-6)
    assert np.linalg.cond(epsilon_hat) < cond_ceiling
    assert np.linalg.cond(modes.phi) < cond_ceiling


@pytest.mark.parametrize(
    "n_index,radius_frac,num_orders",
    [
        (3.48, 0.18, 25),
        (10.0, 0.18, 25),
        (-20 + 2j, 0.3, 25),  # lossy-metal-like
        (3.48, 0.49, 25),  # near-touching pillar
    ],
)
def test_high_contrast_pillar_stress_cases(n_index, radius_frac, num_orders):
    si = Material("stress", n_index**2 if not isinstance(n_index, complex) else n_index)
    pattern = Pattern(background=AIR, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=radius_frac * PERIOD, material=si)])
    _check_case(pattern, num_orders, cond_ceiling=1e3)


@pytest.mark.slow
@pytest.mark.parametrize("num_orders", [81, 121, 225])
def test_high_num_orders_stress(num_orders):
    si = Material("si", 3.48**2)
    pattern = Pattern(background=AIR, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.18 * PERIOD, material=si)])
    _check_case(pattern, num_orders, cond_ceiling=1e3)


def test_sliver_rectangle_stress():
    si = Material("si", 3.48**2)
    pattern = Pattern(background=AIR, shapes=[Rectangle(center=(PERIOD / 2, PERIOD / 2), halfwidth=(0.005, 0.29), material=si)])
    _check_case(pattern, num_orders=49, cond_ceiling=1e3)


def test_near_degenerate_nested_circles_stress():
    """Two nested circles with an almost-equal radius (`1e-4`-scale
    difference) -- the containment-tree subtraction rule
    (`geometry.py::Pattern.containment_tree`) computes a near-cancelling
    `dval` for the inner shape, the adversarial case for that code path."""
    si = Material("si", 3.48**2)
    glass = Material("glass", 2.0**2)
    pattern = Pattern(
        background=AIR,
        shapes=[
            Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.3 * PERIOD, material=si),
            Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.3 * PERIOD - 1e-4, material=glass),
        ],
    )
    _check_case(pattern, num_orders=81, cond_ceiling=1e3)


# ---------------------------------------------------------------------------
# Condition-number WARNING logging
# ---------------------------------------------------------------------------


def test_no_warning_for_well_conditioned_case(caplog):
    si = Material("si", 3.48**2)
    pattern = Pattern(background=AIR, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.18 * PERIOD, material=si)])
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    g = truncate_fourier_orders(lattice, 25, "circular")
    omega = 2 * np.pi
    kx, ky = _kx_ky(lattice, g, omega)
    epsilon_hat = toeplitz_matrix(pattern, lattice, g, 1.0, inverse=False)

    with caplog.at_level(logging.WARNING, logger="sougata_solver.eigenmodes"):
        solve_layer_eigenmodes_patterned(omega, kx, ky, epsilon_hat)
    assert caplog.records == []


def test_warning_fires_when_threshold_is_exceeded(caplog, monkeypatch):
    """Forces the WARNING path by lowering `ILL_CONDITIONED_THRESHOLD`
    below an ordinary, otherwise-unremarkable case's actual condition
    number -- confirms the logging mechanism itself works (message,
    logger name, level), independent of whether a real geometry happens
    to be ill-conditioned enough to trigger it naturally."""
    monkeypatch.setattr(eigenmodes, "ILL_CONDITIONED_THRESHOLD", 1.0)

    si = Material("si", 3.48**2)
    pattern = Pattern(background=AIR, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.18 * PERIOD, material=si)])
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    g = truncate_fourier_orders(lattice, 25, "circular")
    omega = 2 * np.pi
    kx, ky = _kx_ky(lattice, g, omega)
    epsilon_hat = toeplitz_matrix(pattern, lattice, g, 1.0, inverse=False)

    with caplog.at_level(logging.WARNING, logger="sougata_solver.eigenmodes"):
        solve_layer_eigenmodes_patterned(omega, kx, ky, epsilon_hat)

    messages = [r.getMessage() for r in caplog.records]
    assert any("epsilon_hat is ill-conditioned" in m for m in messages)
    assert any("phi is ill-conditioned" in m for m in messages)
