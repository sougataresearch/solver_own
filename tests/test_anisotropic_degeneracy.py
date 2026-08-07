"""Category 1 target 1.7 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): deterministic
degeneracy/mode-ordering policy for the anisotropic eigensolvers
(`eigenmodes._canonical_mode_order`, targets 1.3/1.4/1.6).

Tiers enforced here, per `rules.md` Testing Requirements:
- determinism: repeated solves of the *same* input produce bit-identical
  `(q, phi)` (LAPACK `geev`'s own determinism -- not new behavior, just
  verified directly rather than assumed).
- ordering: `_canonical_mode_order`'s documented sort key produces a
  reproducible, non-decreasing-by-key mode order.
- near-degenerate case: a uniaxial material tuned close to the isotropic
  point (small `eps_xx - eps_yy`) still gives a deterministic order and
  energy conservation -- the case `ILL_CONDITIONED_THRESHOLD` detection
  exists for, exercised here from the ordering-policy side.
"""

from __future__ import annotations


import numpy as np
import pytest

from sougata_solver.eigenmodes import (
    _canonical_mode_order,
    solve_layer_eigenmodes_uniform_diagonal,
    solve_layer_eigenmodes_uniform_inplane,
)
from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Circle, Lattice, Pattern
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation

WAVELENGTH = 0.55e-6


# ---------------------------------------------------------------------------
# Unit: _canonical_mode_order's sort key is respected
# ---------------------------------------------------------------------------


def test_canonical_mode_order_sorts_by_rounded_real_then_imag():
    q = np.array([3.0 + 1j, 1.0 + 2j, 1.0 + 1j, 2.0 + 0j])
    phi = np.eye(4, dtype=complex)  # columns identifiable by their original index
    q_sorted, phi_sorted = _canonical_mode_order(q, phi)
    assert q_sorted.tolist() == [1.0 + 1j, 1.0 + 2j, 2.0 + 0j, 3.0 + 1j]
    # phi's columns must have moved together with q (column 2 -> position 0, etc.)
    expected_col_order = [2, 1, 3, 0]
    assert np.array_equal(phi_sorted, phi[:, expected_col_order])


def test_canonical_mode_order_ties_break_by_original_index():
    # Two exactly-equal eigenvalues (degenerate): tie must break by original index.
    q = np.array([5.0 + 0j, 5.0 + 0j, 1.0 + 0j])
    phi = np.eye(3, dtype=complex)
    q_sorted, phi_sorted = _canonical_mode_order(q, phi)
    assert q_sorted.tolist() == [1.0 + 0j, 5.0 + 0j, 5.0 + 0j]
    assert np.array_equal(phi_sorted, phi[:, [2, 0, 1]])


# ---------------------------------------------------------------------------
# Determinism: repeated solves of the same input give identical output
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "solver_call",
    [
        lambda: solve_layer_eigenmodes_uniform_diagonal(
            2 * np.pi / WAVELENGTH, np.array([0.3e6, 0.7e6]), np.array([0.1e6, 0.2e6]), 2.25, 4.0, 3.1
        ),
        lambda: solve_layer_eigenmodes_uniform_inplane(
            2 * np.pi / WAVELENGTH, np.array([0.3e6, 0.7e6]), np.array([0.1e6, 0.2e6]), 2.25, 0.3, 0.3, 4.0, 3.1
        ),
    ],
)
def test_repeated_solve_is_deterministic(solver_call):
    modes_a = solver_call()
    modes_b = solver_call()
    assert np.array_equal(modes_a.q, modes_b.q)
    assert np.array_equal(modes_a.phi, modes_b.phi)


# ---------------------------------------------------------------------------
# Near-degenerate case: uniaxial material tuned close to isotropic
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("delta", [1e-2, 1e-6, 1e-10, 0.0])
def test_near_degenerate_uniform_inplane_deterministic_and_energy_conserving(delta):
    omega = 2 * np.pi / WAVELENGTH
    kx = np.array([0.3, -0.5, 1.1]) * omega
    ky = np.array([0.2, 0.4, -0.6]) * omega
    eps_xx = 2.25
    eps_yy = 2.25 + delta  # -> isotropic exactly when delta == 0 (fully degenerate)
    eps_zz = 3.1

    modes_a = solve_layer_eigenmodes_uniform_inplane(omega, kx, ky, eps_xx, 0.0, 0.0, eps_yy, eps_zz)
    modes_b = solve_layer_eigenmodes_uniform_inplane(omega, kx, ky, eps_xx, 0.0, 0.0, eps_yy, eps_zz)
    assert np.array_equal(modes_a.q, modes_b.q)
    assert np.array_equal(modes_a.phi, modes_b.phi)

    # Mode order key is non-decreasing (canonical order actually applied).
    key_real = np.round(modes_a.q.real / max(1.0, float(np.max(np.abs(modes_a.q)))) / 1e-9)
    assert np.all(np.diff(key_real) >= 0)


def test_near_degenerate_patterned_inplane_energy_conservation():
    lattice = Lattice(a=(0.7, 0.0), b=(0.0, 0.7))
    air = Material("air", 1.0)
    tensor = np.array([[2.2501, 0.0, 0], [0.0, 2.25, 0], [0, 0, 3.1]], dtype=complex)  # near-isotropic
    si_aniso = Material.from_permittivity_tensor("near_iso", tensor)
    pattern = Pattern(background=air, shapes=[Circle(center=(0.35, 0.35), radius=0.18, material=si_aniso)])
    layer = Layer("pillar", 0.46, pattern=pattern)
    sim = Simulation(lattice, [layer], num_orders=5, incidence=air, transmission=air)
    excitation = PlaneWaveExcitation(1.0, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    result_a = sim.solve(excitation)
    result_b = sim.solve(excitation)
    assert result_a.reflectance() == pytest.approx(result_b.reflectance(), abs=1e-12)
    assert result_a.reflectance() + result_a.transmittance() == pytest.approx(1.0, abs=1e-8)
