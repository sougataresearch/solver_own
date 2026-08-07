"""Phase 4a tests for 2D-periodic patterned layers.

Tiers enforced here:
- unit/reduction: ky=0 reduction to the already-validated 1D solver
  (TE-like block only -- see below for why the TM-like block correctly
  does NOT reduce), and fully-patterned degenerates to uniform.
- independent oracle (eigenvalue-level): `solve_layer_eigenmodes_patterned`'s
  eigenoperator cross-checked against a structurally different formula
  hand-transcribed from `RigorousCoupledWaveAnalysis.jl`
  (`tests/oracles/rcwa_2djl_eigenvalues.py`) -- this is the test that
  would have caught the Epsilon2 bug described below (the old "ky=0
  reduces to 1D" test could not, since both solvers shared the same wrong
  formula and so trivially agreed).
- physical-invariant: energy conservation for moderate-contrast pillar cases.
- no full-R/T external-oracle comparison test yet -- see
  `tests/oracles/rcwa_2d_pillar.py` for what's still open (Phase 4b).
"""

from __future__ import annotations

import numpy as np
import pytest
from oracles.rcwa_2djl_eigenvalues import eigenoperator_eigenvalues

from sougata_solver.eigenmodes import solve_layer_eigenmodes_1d, solve_layer_eigenmodes_patterned
from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.fourier_basis import truncate_fourier_orders
from sougata_solver.fourier_factorization import toeplitz_matrix
from sougata_solver.geometry import Circle, Lattice, Lattice1D, Pattern, Rectangle, Slab
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation


pytestmark = pytest.mark.oracle  # Category 17 target 17.1: system-tier test, cross-checked against a named external oracle


PERIOD = 0.7
AIR = Material("air", 1.0)
SI = Material("si", 3.48**2)
LAYER_THICKNESS = 0.46


def _pillar_pattern() -> Pattern:
    return Pattern(
        background=AIR,
        shapes=[
            Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.18, material=SI),
            Rectangle(center=(0.15, 0.15), halfwidth=(0.08, 0.12), material=SI),
        ],
    )


def _pillar_simulation(num_orders: int = 5) -> Simulation:
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    layer = Layer("pillar", LAYER_THICKNESS, pattern=_pillar_pattern())
    return Simulation(lattice, [layer], num_orders=num_orders, incidence=AIR, transmission=AIR)


# ---------------------------------------------------------------------------
# Reduction: patterned -> uniform when shape == background
# ---------------------------------------------------------------------------


def test_2d_patterned_layer_reduces_to_uniform_when_shapes_match_background():
    pattern = Pattern(
        background=AIR,
        shapes=[
            Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.18, material=AIR),
            Rectangle(center=(0.15, 0.15), halfwidth=(0.08, 0.12), material=AIR),
        ],
    )
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    layer = Layer("patch", LAYER_THICKNESS, pattern=pattern)
    sim = Simulation(lattice, [layer], num_orders=7, incidence=AIR, transmission=AIR)
    excitation = PlaneWaveExcitation(wavelength=1.0, theta=0.0, phi=0.0, s_amplitude=0.7, p_amplitude=0.3)
    result = sim.solve(excitation)
    assert result.reflectance() == pytest.approx(0.0, abs=1e-10)
    assert result.transmittance() == pytest.approx(1.0, abs=1e-8)


# ---------------------------------------------------------------------------
# ky=0 cross-check: 2D-general vs. 1D solver
#
# These do NOT produce the same result in general, and that's correct, not
# a bug: solve_layer_eigenmodes_patterned transcribes S4's true-2D,
# no-polarization-basis closed-form path (fmm_closed.cpp:133-139), which
# uses ordinary Laurent's rule throughout (Epsilon2 = block_diag(epsilon_hat,
# epsilon_hat), kp built from inv(epsilon_hat)) -- S4 only applies Li's
# (1996) inverse-rule correction inside the *separate* 1D branch
# (fmm_closed.cpp:110-132, 0==Lr[2]&&0==Lr[3]), which solve_layer_eigenmodes_1d
# transcribes. An earlier draft of solve_layer_eigenmodes_patterned
# incorrectly copied the 1D branch's Epsilon2 construction into the general
# 2D solver (making the two solvers trivially/circularly agree at ky=0);
# this was caught by actually reading fmm_closed.cpp's true-2D branch and
# fixed -- see eigenmodes.py's docstring for the full citation.
# ---------------------------------------------------------------------------


def test_2d_patterned_ky_zero_te_block_matches_1d():
    """The TE-like (top-left) block never depends on which
    Fourier-factorization rule is used for `einv` -- `kp[:n,:n]` is
    `omega^2*I` regardless of `einv`'s value whenever `ky=0` (its `kappa`
    entries all carry a `ky` factor), and both solvers use the same
    direct-rule `epsilon_hat` for `Epsilon2`'s top-left block -- so this
    part of the two solvers' output is expected, and verified, to agree."""
    num_ord = 5
    n = 2 * num_ord + 1
    omega = 2 * np.pi * 3.48 / 1.0
    kx0 = 0.21
    kx = kx0 + 2 * np.pi * np.arange(-num_ord, num_ord + 1) / PERIOD
    ky = np.zeros(n)
    g = np.column_stack([np.arange(-num_ord, num_ord + 1), np.zeros(n)])

    slab = Slab(center_x=PERIOD / 2, halfwidth=0.18, material=SI)
    lattice = Lattice1D(PERIOD)
    pattern = Pattern(background=AIR, shapes=[slab])
    epsilon_hat = toeplitz_matrix(pattern, lattice, g, wavelength=1.0, inverse=False)
    epsilon_inv_hat = toeplitz_matrix(pattern, lattice, g, wavelength=1.0, inverse=True)

    modes_2d = solve_layer_eigenmodes_patterned(omega, kx, ky, epsilon_hat)
    modes_1d = solve_layer_eigenmodes_1d(omega, kx, ky, epsilon_hat, epsilon_inv_hat)

    assert modes_2d.q.shape == modes_1d.q.shape
    te_2d = sorted((modes_2d.q[:n] ** 2).real.tolist())
    te_1d = sorted((modes_1d.q[:n] ** 2).real.tolist())
    assert te_2d == pytest.approx(te_1d, abs=1e-10)


def test_2d_patterned_ky_zero_tm_block_differs_from_1d():
    """The TM-like (bottom-right) block is exactly where the two
    Fourier-factorization rules diverge (Laurent's rule vs. Li's 1996
    inverse rule) -- for a genuinely discontinuous pattern the two should
    give visibly different eigenvalues, not agree. A future regression that
    reintroduces the ky=0-reduces-to-1D bug would make this assertion fail
    (the two would go back to agreeing)."""
    num_ord = 5
    n = 2 * num_ord + 1
    omega = 2 * np.pi * 3.48 / 1.0
    kx0 = 0.21
    kx = kx0 + 2 * np.pi * np.arange(-num_ord, num_ord + 1) / PERIOD
    ky = np.zeros(n)
    g = np.column_stack([np.arange(-num_ord, num_ord + 1), np.zeros(n)])

    slab = Slab(center_x=PERIOD / 2, halfwidth=0.18, material=SI)
    lattice = Lattice1D(PERIOD)
    pattern = Pattern(background=AIR, shapes=[slab])
    epsilon_hat = toeplitz_matrix(pattern, lattice, g, wavelength=1.0, inverse=False)
    epsilon_inv_hat = toeplitz_matrix(pattern, lattice, g, wavelength=1.0, inverse=True)

    modes_2d = solve_layer_eigenmodes_patterned(omega, kx, ky, epsilon_hat)
    modes_1d = solve_layer_eigenmodes_1d(omega, kx, ky, epsilon_hat, epsilon_inv_hat)

    tm_2d = sorted((modes_2d.q[n:] ** 2).real.tolist())
    tm_1d = sorted((modes_1d.q[n:] ** 2).real.tolist())
    assert max(abs(a - b) for a, b in zip(tm_2d, tm_1d)) > 1e-3


# ---------------------------------------------------------------------------
# Independent eigenoperator oracle (RigorousCoupledWaveAnalysis.jl)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("num_orders", [5, 7, 9])
@pytest.mark.parametrize("theta_deg", [0.0, 20.0])
def test_2d_patterned_eigenvalues_match_rcwa_jl_oracle(num_orders, theta_deg):
    """Cross-check `solve_layer_eigenmodes_patterned`'s `q^2` eigenvalues
    against `tests/oracles/rcwa_2djl_eigenvalues.py` (hand-transcribed from
    `RigorousCoupledWaveAnalysis.jl`, a structurally different eigenoperator
    derivation than S4's `Epsilon2 @ kp - coupling` route). Both consume
    the *same*, already-Phase-2-validated `epsilon_hat`, isolating the
    eigenoperator-construction step -- this is the test that would have
    caught the earlier `Epsilon2` bug (see module docstring)."""
    wavelength = 1.0
    omega = 2 * np.pi / wavelength
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    g = truncate_fourier_orders(lattice, num_orders, "circular")
    lk = lattice.reciprocal_vectors()
    kx0 = omega * np.sin(np.radians(theta_deg))
    kx = kx0 + 2 * np.pi * (g[:, 0] * lk[0, 0] + g[:, 1] * lk[1, 0])
    ky = 2 * np.pi * (g[:, 0] * lk[0, 1] + g[:, 1] * lk[1, 1])

    epsilon_hat = toeplitz_matrix(_pillar_pattern(), lattice, g, wavelength, inverse=False)
    modes = solve_layer_eigenmodes_patterned(omega, kx, ky, epsilon_hat)
    ours = np.sort((modes.q**2).real)

    oracle = np.sort(eigenoperator_eigenvalues(omega, kx, ky, epsilon_hat).real)

    assert ours == pytest.approx(oracle, abs=1e-6)


# ---------------------------------------------------------------------------
# Physical invariant: energy conservation for moderate-contrast pillar cases
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("theta_deg", [0.0, 30.0])
@pytest.mark.parametrize("s_amp,p_amp", [(1.0, 0.0), (0.0, 1.0), (0.6, 0.8)])
def test_energy_conservation_moderate_contrast_pillar(theta_deg, s_amp, p_amp):
    sim = _pillar_simulation(num_orders=7)
    excitation = PlaneWaveExcitation(
        wavelength=1.0,
        theta=np.radians(theta_deg),
        phi=0.0,
        s_amplitude=s_amp,
        p_amplitude=p_amp,
    )
    result = sim.solve(excitation)
    de = result.diffraction_efficiencies()
    total = sum(de_r + de_t for de_r, de_t in de.values())
    assert total == pytest.approx(1.0, abs=1e-8)
    assert result.reflectance() + result.transmittance() == pytest.approx(1.0, abs=1e-8)


# ---------------------------------------------------------------------------
# Moderate-orders sanity: runnable and finite R/T for a moderate case
# ---------------------------------------------------------------------------


def test_2d_pillar_runs_end_to_end():
    sim = _pillar_simulation(num_orders=5)
    excitation = PlaneWaveExcitation(wavelength=1.0, theta=0.0, phi=0.0, s_amplitude=1.0, p_amplitude=0.0)
    result = sim.solve(excitation)
    assert 0.0 <= result.reflectance() <= 1.0
    assert 0.0 <= result.transmittance() <= 1.0
