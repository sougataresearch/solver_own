"""Category 9 targets 9.2-9.7 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`, Phase
7): real-space field reconstruction. See `CONVENTIONS.md`'s "Real-space
field reconstruction" section for the transcribed formulas
(`fields.modal_field_components`/`reconstruct_field_at_points`/
`propagate_amplitudes`) and `decisions.md` ADR-015 for
`smatrix.interior_amplitudes`'s independent-derivation status and the
validation this file provides for it.

Every numeric claim in this file's docstrings was confirmed by direct
interactive computation before being encoded as a test (per `rules.md`'s
"never fabricate a benchmark" rule) -- see `memory.md`'s Phase 7 entry for
the session account.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.eigenmodes import solve_layer_eigenmodes_uniform
from sougata_solver.fields import (
    modal_field_components,
    propagate_amplitudes,
    reconstruct_field_at_points,
    save_field_grid_npz,
    z_poynting_flux,
)
from sougata_solver.fourier_factorization import toeplitz_matrix
from sougata_solver.geometry import Circle, Lattice, Lattice1D, Pattern, Slab
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation
from sougata_solver.smatrix import SMatrixStack, interior_amplitudes

WAVELENGTH = 0.6e-6


# ---------------------------------------------------------------------------
# 9.2 Uniform-layer reconstruction vs. the analytic plane wave
# ---------------------------------------------------------------------------


def test_uniform_layer_field_at_z0_matches_incident_field_xy():
    """`modal_field_components`' Ex/Ey at the incidence half-space's own
    z=0 reference plane must reproduce `PlaneWaveExcitation.incident_field_xy()`
    exactly -- the two are supposed to describe the same physical quantity
    by construction (`incident_mode_amplitude` was built to invert this
    exact relation)."""
    theta, phi = math.radians(20.0), math.radians(10.0)
    excitation = PlaneWaveExcitation(WAVELENGTH, theta, phi, s_amplitude=0.6, p_amplitude=0.8)
    omega = excitation.omega()
    kx0, ky0 = excitation.k_parallel(1.0)
    kx, ky = np.array([kx0]), np.array([ky0])
    modes = solve_layer_eigenmodes_uniform(omega, kx, ky, eps=1.0)
    a0 = excitation.incident_mode_amplitude(modes, num_orders=1, zeroth_order_index=0)
    b0 = np.zeros_like(a0)

    ex, ey, _ez, _hx, _hy, _hz = modal_field_components(omega, kx, ky, modes.q, modes.kp, modes.phi, 1.0, a0, b0)
    ex_expected, ey_expected = excitation.incident_field_xy()
    assert ex[0] == pytest.approx(ex_expected, abs=1e-10)
    assert ey[0] == pytest.approx(ey_expected, abs=1e-10)


def test_uniform_layer_reconstruction_matches_analytic_plane_wave():
    """Full pipeline (`propagate_amplitudes` + `reconstruct_field_at_points`)
    at an arbitrary `(x, y, z)`, compared against the closed-form plane
    wave `E(x,y,z) = E0 * exp(i*(kx*x + ky*y + q*z))` -- confirmed to match
    to `~1e-10` (double-precision floating point), not just "close"."""
    theta, phi = math.radians(20.0), math.radians(10.0)
    excitation = PlaneWaveExcitation(WAVELENGTH, theta, phi, s_amplitude=0.6, p_amplitude=0.8)
    omega = excitation.omega()
    kx0, ky0 = excitation.k_parallel(1.0)
    kx, ky = np.array([kx0]), np.array([ky0])
    modes = solve_layer_eigenmodes_uniform(omega, kx, ky, eps=1.0)
    a0 = excitation.incident_mode_amplitude(modes, num_orders=1, zeroth_order_index=0)
    b0 = np.zeros_like(a0)
    ex0, ey0 = excitation.incident_field_xy()

    x, y, z = 0.3e-6, -0.2e-6, 0.15e-6
    a_z, b_z = propagate_amplitudes(modes.q, z, a0, b0)
    ex, ey, _ez, _hx, _hy, _hz = modal_field_components(omega, kx, ky, modes.q, modes.kp, modes.phi, 1.0, a_z, b_z)
    ex_recon = reconstruct_field_at_points(kx, ky, np.array(x), np.array(y), ex)
    ey_recon = reconstruct_field_at_points(kx, ky, np.array(x), np.array(y), ey)

    phase = np.exp(1j * (kx0 * x + ky0 * y + modes.q[0] * z))
    assert ex_recon == pytest.approx(ex0 * phase, abs=1e-10)
    assert ey_recon == pytest.approx(ey0 * phase, abs=1e-10)


def test_uniform_layer_field_is_transverse():
    """Independent physics check (not just formula reproduction): a plane
    wave in a source-free isotropic medium must satisfy `k . E = 0` and
    `k . H = 0` -- confirmed here to `~1e-9` relative to field magnitudes
    of order 1 (i.e. to full double precision given the ~1e7 rad/m scale
    of kx/ky/q)."""
    theta, phi = math.radians(20.0), math.radians(10.0)
    excitation = PlaneWaveExcitation(WAVELENGTH, theta, phi, s_amplitude=0.6, p_amplitude=0.8)
    omega = excitation.omega()
    kx0, ky0 = excitation.k_parallel(1.0)
    kx, ky = np.array([kx0]), np.array([ky0])
    modes = solve_layer_eigenmodes_uniform(omega, kx, ky, eps=1.0)
    a0 = excitation.incident_mode_amplitude(modes, num_orders=1, zeroth_order_index=0)
    b0 = np.zeros_like(a0)

    ex, ey, ez, hx, hy, hz = modal_field_components(omega, kx, ky, modes.q, modes.kp, modes.phi, 1.0, a0, b0)
    k_dot_e = kx0 * ex[0] + ky0 * ey[0] + modes.q[0] * ez[0]
    k_dot_h = kx0 * hx[0] + ky0 * hy[0] + modes.q[0] * hz[0]
    assert abs(k_dot_e) < 1e-6
    assert abs(k_dot_h) < 1e-6


# ---------------------------------------------------------------------------
# 9.3 Interior amplitudes
# ---------------------------------------------------------------------------


def _three_layer_stack(num_orders=1):
    air = Material("air", 1.0)
    glass = Material("glass", 1.5**2)
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))
    l1 = Layer("L1", 0.2e-6, material=Material("L1", 2.0**2))
    l2 = Layer("L2", 0.15e-6, material=Material("L2", 3.0**2))
    sim = Simulation(lattice, [l1, l2], num_orders=num_orders, incidence=air, transmission=glass)
    return sim, l1, l2


def test_interior_amplitudes_at_full_stack_reproduce_transmitted_amplitude():
    """Zero-free-parameter consistency check (`decisions.md` ADR-015):
    `interior_amplitudes` called with the *full* S-matrix (the partial
    stack up to the very last layer) must reproduce `a_transmitted` and
    `b~0` (nothing incident from the transmission side) exactly -- both
    quantities are already independently computed by `Simulation.solve`."""
    sim, l1, l2 = _three_layer_stack()
    excitation = PlaneWaveExcitation(WAVELENGTH, math.radians(25.0), math.radians(15.0), s_amplitude=0.7, p_amplitude=0.5)
    result = sim.solve(excitation)

    n2 = 2 * sim.num_orders
    thicknesses = [layer.thickness for layer in sim.layer_stack]
    stack = SMatrixStack(thicknesses, result.all_modes)

    a_full, b_full = interior_amplitudes(stack.full_smatrix(), n2, result.a0, result.b_reflected)
    assert a_full == pytest.approx(result.a_transmitted, abs=1e-12)
    assert np.max(np.abs(b_full)) < 1e-12


def test_interior_amplitudes_uniform_stack_energy_conservation():
    """Recovered interior amplitudes describe a genuinely energy-
    conserving local field: forward minus backward modal power at an
    interior interface must be conserved along a lossless stack and equal
    the already-validated transmitted power."""
    sim, l1, l2 = _three_layer_stack()
    excitation = PlaneWaveExcitation(WAVELENGTH, math.radians(25.0), math.radians(15.0), s_amplitude=0.7, p_amplitude=0.5)
    result = sim.solve(excitation)
    n2 = 2 * sim.num_orders
    thicknesses = [layer.thickness for layer in sim.layer_stack]
    stack = SMatrixStack(thicknesses, result.all_modes)

    a1_top, b1_top = interior_amplitudes(stack.partial_smatrix_up_to(1), n2, result.a0, result.b_reflected)
    modes1 = result.all_modes[1]
    fwd, bwd = z_poynting_flux(excitation.omega(), modes1.q, modes1.kp, modes1.phi, a1_top, b1_top)
    zeros = np.zeros_like(result.a0)
    incident_power, _ = z_poynting_flux(
        excitation.omega(), result.all_modes[0].q, result.all_modes[0].kp, result.all_modes[0].phi, result.a0, zeros
    )
    net_power_fraction = (fwd + bwd).real / incident_power.real
    assert net_power_fraction == pytest.approx(result.transmittance(), abs=1e-8)


# ---------------------------------------------------------------------------
# 9.5 Field continuity (no surface current at a dielectric interface)
# ---------------------------------------------------------------------------


def test_tangential_field_continuity_across_interior_interface():
    """Reconstruct Ex/Ey/Hx/Hy at the *same* physical point from both
    sides of the L1/L2 interface (layer 1's bottom vs. layer 2's top,
    using each layer's own kp/phi/material) -- must match exactly, a
    genuine boundary-condition check (no surface current at an ordinary
    dielectric interface), confirmed to `~1e-10`."""
    sim, l1, l2 = _three_layer_stack()
    excitation = PlaneWaveExcitation(WAVELENGTH, math.radians(25.0), math.radians(15.0), s_amplitude=0.7, p_amplitude=0.5)
    result = sim.solve(excitation)
    omega = excitation.omega()
    n2 = 2 * sim.num_orders
    thicknesses = [layer.thickness for layer in sim.layer_stack]
    stack = SMatrixStack(thicknesses, result.all_modes)

    a1_top, b1_top = interior_amplitudes(stack.partial_smatrix_up_to(1), n2, result.a0, result.b_reflected)
    a2_top, b2_top = interior_amplitudes(stack.partial_smatrix_up_to(2), n2, result.a0, result.b_reflected)

    modes1, modes2 = result.all_modes[1], result.all_modes[2]
    a1_bot, b1_bot = propagate_amplitudes(modes1.q, l1.thickness, a1_top, b1_top)

    kx0, ky0 = excitation.k_parallel(1.0)
    kx, ky = np.array([kx0]), np.array([ky0])
    eps1 = complex(Material("L1", 2.0**2).epsilon_tensor(WAVELENGTH)[0, 0])
    eps2 = complex(Material("L2", 3.0**2).epsilon_tensor(WAVELENGTH)[0, 0])

    field_below = modal_field_components(omega, kx, ky, modes1.q, modes1.kp, modes1.phi, 1.0 / eps1, a1_bot, b1_bot)
    field_above = modal_field_components(omega, kx, ky, modes2.q, modes2.kp, modes2.phi, 1.0 / eps2, a2_top, b2_top)

    ex1, ey1, _ez1, hx1, hy1, _hz1 = field_below
    ex2, ey2, _ez2, hx2, hy2, _hz2 = field_above
    assert ex1 == pytest.approx(ex2, abs=1e-10)
    assert ey1 == pytest.approx(ey2, abs=1e-10)
    assert hx1 == pytest.approx(hx2, abs=1e-10)
    assert hy1 == pytest.approx(hy2, abs=1e-10)


# ---------------------------------------------------------------------------
# 9.4 1D inverse Fourier sum (periodicity) and 9.6 field-flux (2D pillar)
# ---------------------------------------------------------------------------


def _pillar_result(num_orders=9):
    period = 0.7e-6
    air = Material("air", 1.0)
    si = Material("si", 3.48**2)
    lattice = Lattice((period, 0.0), (0.0, period))
    pattern = Pattern(background=air, shapes=[Circle(center=(period / 2, period / 2), radius=0.18e-6, material=si)])
    layer = Layer("pillar", 0.3e-6, pattern=pattern)
    sim = Simulation(lattice, [layer], num_orders=num_orders, incidence=air, transmission=air)
    excitation = PlaneWaveExcitation(WAVELENGTH, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    result = sim.solve(excitation)
    return sim, lattice, pattern, layer, excitation, result, period


def test_1d_grating_field_reconstruction_is_periodic():
    """9.4: reconstruct Ex along x at fixed z within a 1D lamellar grating
    and confirm `E(x + period) == E(x)` -- not automatic from the code
    structure alone (a wrong reciprocal-lattice `kx` would break this),
    confirmed directly over more than one period."""
    period = 0.7e-6
    fill_factor = 0.3
    air = Material("air", 1.0)
    si = Material("si", 3.48**2)
    lattice = Lattice1D(period)
    pattern = Pattern(background=air)
    pattern.add(Slab(center_x=0.0, halfwidth=0.5 * fill_factor * period, material=si))
    layer = Layer("grating", 0.46e-6, pattern=pattern)
    sim = Simulation(lattice, [layer], num_orders=9, incidence=air, transmission=air)
    excitation = PlaneWaveExcitation(WAVELENGTH, math.radians(10.0), 0.0, s_amplitude=1.0, p_amplitude=0.0)
    result = sim.solve(excitation)

    n2 = 2 * sim.num_orders
    thicknesses = [ly.thickness for ly in sim.layer_stack]
    stack = SMatrixStack(thicknesses, result.all_modes)
    a_top, b_top = interior_amplitudes(stack.partial_smatrix_up_to(1), n2, result.a0, result.b_reflected)
    modes = result.all_modes[1]
    a_z, b_z = propagate_amplitudes(modes.q, 0.2e-6, a_top, b_top)

    lk = lattice.reciprocal_vectors()
    kx0, _ky0 = excitation.k_parallel(1.0)
    kx = kx0 + 2 * np.pi * result.g[:, 0] * lk[0, 0]
    ky = np.zeros_like(kx)

    epsilon_hat = toeplitz_matrix(pattern, lattice, result.g, WAVELENGTH, inverse=False)
    epsilon_inv_hat = toeplitz_matrix(pattern, lattice, result.g, WAVELENGTH, inverse=True)
    ex, _ey, _ez, _hx, _hy, _hz = modal_field_components(
        excitation.omega(), kx, ky, modes.q, modes.kp, modes.phi, epsilon_inv_hat, a_z, b_z
    )

    x = np.linspace(-0.3, 2.3, 53) * period  # spans more than 2 periods
    y = np.zeros_like(x)
    field = reconstruct_field_at_points(kx, ky, x, y, ex)
    field_shifted = reconstruct_field_at_points(kx, ky, x + period, y, ex)
    assert field_shifted == pytest.approx(field, abs=1e-9)


def test_2d_pillar_field_reconstruction_flux_matches_transmittance():
    """9.6/9.7 combined: reconstruct the full 6-component field over a 2D
    grid spanning one unit cell at mid-layer depth, integrate the
    real-space Poynting flux (`Sz = Re(Ex*conj(Hy) - Ey*conj(Hx))`, note
    the missing `0.5` -- see `CONVENTIONS.md`), and confirm it matches
    `SimulationResult.transmittance()` -- the category's own exit
    criterion ("field-derived flux agrees with solver R/T"), and this
    project's first genuinely 2D field-reconstruction test."""
    sim, lattice, pattern, layer, excitation, result, period = _pillar_result()
    omega = excitation.omega()
    n2 = 2 * sim.num_orders
    thicknesses = [ly.thickness for ly in sim.layer_stack]
    stack = SMatrixStack(thicknesses, result.all_modes)
    a_top, b_top = interior_amplitudes(stack.partial_smatrix_up_to(1), n2, result.a0, result.b_reflected)
    modes = result.all_modes[1]
    a_z, b_z = propagate_amplitudes(modes.q, layer.thickness / 2, a_top, b_top)

    lk = lattice.reciprocal_vectors()
    kx = 2 * np.pi * (result.g[:, 0] * lk[0, 0] + result.g[:, 1] * lk[1, 0])
    ky = 2 * np.pi * (result.g[:, 0] * lk[0, 1] + result.g[:, 1] * lk[1, 1])

    epsilon_hat = toeplitz_matrix(pattern, lattice, result.g, WAVELENGTH, inverse=False)
    einv = np.linalg.inv(epsilon_hat)
    ex, ey, _ez, hx, hy, _hz = modal_field_components(omega, kx, ky, modes.q, modes.kp, modes.phi, einv, a_z, b_z)

    n_grid = 60
    xs = (np.arange(n_grid) + 0.5) / n_grid * period
    ys = (np.arange(n_grid) + 0.5) / n_grid * period
    x_grid, y_grid = np.meshgrid(xs, ys, indexing="ij")
    ex_r = reconstruct_field_at_points(kx, ky, x_grid, y_grid, ex)
    ey_r = reconstruct_field_at_points(kx, ky, x_grid, y_grid, ey)
    hx_r = reconstruct_field_at_points(kx, ky, x_grid, y_grid, hx)
    hy_r = reconstruct_field_at_points(kx, ky, x_grid, y_grid, hy)

    sz = np.real(ex_r * np.conj(hy_r) - ey_r * np.conj(hx_r))
    net_flux_fraction = np.mean(sz)  # normalized: incident power is 1 (unit s_amplitude, normal incidence)
    assert net_flux_fraction == pytest.approx(result.transmittance(), abs=1e-6)


# ---------------------------------------------------------------------------
# 9.8 Field export (NumPy .npz)
# ---------------------------------------------------------------------------


def test_save_field_grid_npz_roundtrip(tmp_path):
    x = np.linspace(0.0, 0.7e-6, 5)
    y = np.linspace(0.0, 0.7e-6, 5)
    x_grid, y_grid = np.meshgrid(x, y, indexing="ij")
    ex = (x_grid + 1j * y_grid).astype(complex)
    ey = np.zeros_like(ex)
    ez = np.zeros_like(ex)
    hx = np.zeros_like(ex)
    hy = np.ones_like(ex)
    hz = np.zeros_like(ex)

    path = save_field_grid_npz(
        str(tmp_path / "field_grid"), x_grid, y_grid, 0.15e-6, ex, ey, ez, hx, hy, hz,
        wavelength=0.6e-6, theta_deg=10.0,
    )
    assert path.exists()
    assert path.name == "field_grid.npz"

    loaded = np.load(path)
    assert np.array_equal(loaded["x"], x_grid)
    assert np.allclose(loaded["Ex"], ex)
    assert np.allclose(loaded["Hy"], hy)
    assert float(loaded["wavelength"]) == pytest.approx(0.6e-6)
    assert float(loaded["theta_deg"]) == pytest.approx(10.0)


def test_save_field_grid_npz_does_not_double_extension(tmp_path):
    path = save_field_grid_npz(
        str(tmp_path / "already_named.npz"),
        np.zeros(2), np.zeros(2), 0.0,
        np.zeros(2, dtype=complex), np.zeros(2, dtype=complex), np.zeros(2, dtype=complex),
        np.zeros(2, dtype=complex), np.zeros(2, dtype=complex), np.zeros(2, dtype=complex),
    )
    assert path.name == "already_named.npz"
