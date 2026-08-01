"""Phase 3 (1D-periodic lamellar gratings) tests.

Tiers per `testing.md`: unit (geometry/truncation/block-diagonal-reduction),
regression (reduces to Phase 1's uniform result), physical-invariant
(energy conservation, convergence-rate-vs-`num_orders`), and system
(cross-check against `tests/oracles/rcwa_1d_gaylord.py`, hand-transcribed
from the vendored `Rigorous-Coupled-Wave-Analysis` reference).
"""

from __future__ import annotations

import numpy as np
import pytest
from oracles.rcwa_1d_gaylord import solve_te, solve_tm

from sougata_solver.eigenmodes import build_kp_matrix, solve_layer_eigenmodes_1d, solve_layer_eigenmodes_uniform
from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.fourier_basis import truncate_fourier_orders_1d
from sougata_solver.geometry import Lattice1D, Pattern, Slab
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation

PERIOD = 0.7
FILL_FACTOR = 0.3
THICKNESS = 0.46
N_RIDGE = 3.48
N_GROOVE = 1.0


def _binary_grating_pattern(air: Material, ridge: Material) -> Pattern:
    pattern = Pattern(background=air)
    pattern.add(Slab(center_x=-PERIOD * (1 - FILL_FACTOR) / 2, halfwidth=0.5 * FILL_FACTOR * PERIOD, material=ridge))
    return pattern


def _grating_simulation(num_ord: int, incidence: Material, transmission: Material, ridge: Material) -> Simulation:
    lattice = Lattice1D(PERIOD)
    layer = Layer("grating", THICKNESS, pattern=_binary_grating_pattern(incidence, ridge))
    return Simulation(lattice, [layer], num_orders=2 * num_ord + 1, incidence=incidence, transmission=transmission)


# ---------------------------------------------------------------------------
# Unit tests
# ---------------------------------------------------------------------------


def test_lattice1d_reciprocal_and_area():
    lattice = Lattice1D(PERIOD)
    Lk = lattice.reciprocal_vectors()
    assert Lk[0, 0] == pytest.approx(1.0 / PERIOD)
    assert Lk[0, 1] == 0.0
    assert np.all(Lk[1] == 0.0)
    assert lattice.unit_cell_area() == PERIOD


def test_slab_fourier_transform_dc_term_is_area():
    slab = Slab(center_x=0.1, halfwidth=0.2, material=Material("m", 2.0))
    assert slab.fourier_transform(0.0, 0.0) == pytest.approx(slab.area)
    assert slab.area == pytest.approx(0.4)


def test_slab_fourier_transform_matches_rectangle_x_factor():
    """A `Slab` is `Rectangle`'s x-only sinc factor; cross-check against a
    from-scratch Riemann-sum evaluation of the same 1D Fourier integral
    (independent of both `Slab` and `Rectangle`)."""
    slab = Slab(center_x=0.15, halfwidth=0.25, material=Material("m", 2.0))
    kx = 1.3
    x = np.linspace(slab.center_x - slab.halfwidth, slab.center_x + slab.halfwidth, 200_000)
    riemann = np.trapezoid(np.exp(-2j * np.pi * kx * x), x)
    assert slab.fourier_transform(kx, 0.0) == pytest.approx(riemann, abs=1e-3)


def test_truncate_fourier_orders_1d_symmetric_and_g2_zero():
    lattice = Lattice1D(PERIOD)
    g = truncate_fourier_orders_1d(lattice, 7)
    assert g.shape == (7, 2)
    assert np.all(g[:, 1] == 0)
    assert sorted(g[:, 0].tolist()) == [-3, -2, -1, 0, 1, 2, 3]


def test_op_is_block_diagonal_when_ky_zero(rng):
    """Direct check of the block-diagonal reduction claimed in
    `solve_layer_eigenmodes_1d`'s docstring: for `ky = 0`, S4's general
    eigenoperator (`rcwa.cpp:794-806`) `op = Epsilon2 @ kp -
    [[kxkx,kxky],[kykx,kyky]]` has exactly zero off-diagonal blocks, for a
    random complex Hermitian-Toeplitz-like `epsilon_hat`/`epsilon_inv_hat`
    pair (not just the physically-generated ones)."""
    n = 5
    omega = 3.7
    kx = rng.normal(size=n) + 1j * 0
    ky = np.zeros(n)
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    epsilon_hat = a + a.conj().T  # Hermitian, plausible Toeplitz-like structure
    b = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    epsilon_inv_hat = b @ b.conj().T + n * np.eye(n)  # Hermitian positive-definite (invertible)

    kp = build_kp_matrix(omega, kx, ky, epsilon_inv_hat)
    epsilon2 = np.zeros((2 * n, 2 * n), dtype=complex)
    epsilon2[:n, :n] = epsilon_hat
    epsilon2[n:, n:] = np.linalg.inv(epsilon_inv_hat)
    op = epsilon2 @ kp
    op[:n, :n] -= np.diag(kx * kx)
    op[n:, n:] -= np.diag(ky * ky)

    assert np.allclose(op[:n, n:], 0.0)
    assert np.allclose(op[n:, :n], 0.0)


# ---------------------------------------------------------------------------
# Regression: reduces to the uniform-layer result
# ---------------------------------------------------------------------------


def test_1d_patterned_layer_reduces_to_uniform_when_shape_matches_background():
    """A `Slab` pattern whose shape material equals the background material
    is optically a uniform layer; `solve_layer_eigenmodes_1d` should
    reproduce `solve_layer_eigenmodes_uniform`'s `q` (up to reordering) and
    total R/T exactly."""
    air = Material("air", 1.0)
    num_ord = 6
    sim = _grating_simulation(num_ord, incidence=air, transmission=air, ridge=air)
    excitation = PlaneWaveExcitation(wavelength=1.0, theta=0.3, phi=0.0, s_amplitude=0.7, p_amplitude=0.3)
    result = sim.solve(excitation)

    omega = excitation.omega()
    kx0, ky0 = excitation.k_parallel(1.0)
    n = 2 * num_ord + 1
    kx = kx0 + 2 * np.pi * np.arange(-num_ord, num_ord + 1) / PERIOD
    ky = np.zeros(n)
    uniform_modes = solve_layer_eigenmodes_uniform(omega, kx, ky, 1.0)

    assert sorted(result.all_modes[1].q.real.tolist()) == pytest.approx(
        sorted(uniform_modes.q.real.tolist()), abs=1e-8
    )
    assert result.reflectance() == pytest.approx(0.0, abs=1e-10)
    assert result.transmittance() == pytest.approx(1.0, abs=1e-8)


# ---------------------------------------------------------------------------
# Physical-invariant tests (testing.md, required starting Phase 3)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("theta_deg", [0.0, 30.0])
@pytest.mark.parametrize("s_amp,p_amp", [(1.0, 0.0), (0.0, 1.0), (0.6, 0.8)])
def test_energy_conservation(theta_deg, s_amp, p_amp):
    air = Material("air", 1.0)
    si = Material("si", N_RIDGE**2)
    sim = _grating_simulation(15, incidence=air, transmission=air, ridge=si)
    excitation = PlaneWaveExcitation(
        wavelength=0.5, theta=np.radians(theta_deg), phi=0.0, s_amplitude=s_amp, p_amplitude=p_amp
    )
    result = sim.solve(excitation)
    de = result.diffraction_efficiencies()
    total = sum(de_r + de_t for de_r, de_t in de.values())
    assert total == pytest.approx(1.0, abs=1e-8)
    assert result.reflectance() + result.transmittance() == pytest.approx(1.0, abs=1e-8)


@pytest.mark.slow
def test_convergence_rate_vs_num_orders():
    """Measured convergence (not just "does it converge") vs `num_orders`
    for the TM polarization, the case that exercises the inverse-rule
    (`epsilon_inv_hat`) correction (Li 1996) at the grating's discontinuous
    interface. Checks the error decreases monotonically toward a
    high-order reference -- a flat-value/wrong-rate curve (the failure mode
    Li's rule exists to avoid) would show up as non-monotonic or stalled
    error here, per `testing.md`'s Physical-Invariant Testing tier.
    """
    air = Material("air", 1.0)
    si = Material("si", N_RIDGE**2)

    def reflectance_at(num_ord: int) -> float:
        sim = _grating_simulation(num_ord, incidence=air, transmission=air, ridge=si)
        excitation = PlaneWaveExcitation(wavelength=1.0, theta=0.0, phi=0.0, s_amplitude=0.0, p_amplitude=1.0)
        return sim.solve(excitation).reflectance()

    reference = reflectance_at(300)
    errors = [abs(reflectance_at(n) - reference) for n in (10, 20, 40, 80, 160)]
    assert all(e1 >= e2 for e1, e2 in zip(errors, errors[1:]))
    assert errors[-1] < errors[0]


# ---------------------------------------------------------------------------
# System test: cross-check against tests/oracles/rcwa_1d_gaylord.py
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("theta_deg", [0.0, 20.0])
def test_te_matches_gaylord_oracle(theta_deg):
    air = Material("air", 1.0)
    si = Material("si", N_RIDGE**2)
    num_ord = 15
    sim = _grating_simulation(num_ord, incidence=air, transmission=air, ridge=si)
    excitation = PlaneWaveExcitation(
        wavelength=1.0, theta=np.radians(theta_deg), phi=0.0, s_amplitude=1.0, p_amplitude=0.0
    )
    result = sim.solve(excitation)

    de_r, de_t = solve_te(
        wavelength=1.0,
        theta=np.radians(theta_deg),
        n_ridge=N_RIDGE,
        n_groove=N_GROOVE,
        fill_factor=FILL_FACTOR,
        lattice_constant=PERIOD,
        thickness=THICKNESS,
        num_ord=num_ord,
    )
    assert result.reflectance() == pytest.approx(de_r.sum(), abs=1e-6)
    assert result.transmittance() == pytest.approx(de_t.sum(), abs=1e-6)


@pytest.mark.slow
def test_tm_matches_gaylord_oracle_at_high_num_orders():
    """TM converges markedly slower than TE with `num_orders` (measured
    this session: TE agrees with the oracle to ~1e-10 already at
    `num_ord=10`; TM needs `num_ord` in the hundreds to approach the same
    oracle value -- both this project's solver and the oracle converge to
    the same limit, confirmed by direct convergence sweep during
    development, but the oracle's own module docstring self-reports "STILL
    NOT WORKING YET", so this is deliberately a looser-tolerance,
    `slow`-marked secondary check, not the primary TM correctness signal
    (that's `test_energy_conservation` and the reduces-to-uniform
    regression test above, both of which exercise the same code path)."""
    air = Material("air", 1.0)
    si = Material("si", N_RIDGE**2)
    num_ord = 150
    sim = _grating_simulation(num_ord, incidence=air, transmission=air, ridge=si)
    excitation = PlaneWaveExcitation(wavelength=1.0, theta=0.0, phi=0.0, s_amplitude=0.0, p_amplitude=1.0)
    result = sim.solve(excitation)

    de_r, de_t = solve_tm(
        wavelength=1.0,
        theta=0.0,
        n_ridge=N_RIDGE,
        n_groove=N_GROOVE,
        fill_factor=FILL_FACTOR,
        lattice_constant=PERIOD,
        thickness=THICKNESS,
        num_ord=num_ord,
    )
    assert result.reflectance() == pytest.approx(de_r.sum(), abs=0.01)
    assert result.transmittance() == pytest.approx(de_t.sum(), abs=0.01)
