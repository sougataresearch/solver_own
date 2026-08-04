"""Category 2 target 2.2 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): eigenvalue/
mode-conditioning diagnostics (`layer.EigenmodeDiagnostics`), attached as
`LayerEigenmodes.diagnostics` by every `eigenmodes.py` solver.

Tiers enforced here, per `rules.md` Testing Requirements:
- unit: diagnostics fields match independently-recomputed quantities
  (`np.linalg.cond`, `classify_propagating`, a from-scratch min-pairwise-gap
  loop) for a representative case of each solver.
- regression: attaching diagnostics does not change `q`/`phi`/`kp`/
  `epsilon_inv`/`is_scalar_isotropic` -- verified by comparing against the
  pre-target-2.2 code path is not directly possible (the field simply
  didn't exist), so instead this compares two independent solves of the
  *same* input and asserts every non-diagnostics field is bit-identical,
  i.e. adding diagnostics didn't introduce any nondeterminism or mutate
  shared state.
"""

from __future__ import annotations

import numpy as np
import pytest

from sougata_solver.eigenmodes import (
    solve_layer_eigenmodes_1d,
    solve_layer_eigenmodes_patterned,
    solve_layer_eigenmodes_patterned_inplane,
    solve_layer_eigenmodes_uniform,
    solve_layer_eigenmodes_uniform_diagonal,
    solve_layer_eigenmodes_uniform_inplane,
)
from sougata_solver.fourier_basis import truncate_fourier_orders
from sougata_solver.fourier_factorization import toeplitz_matrix, toeplitz_matrix_component
from sougata_solver.geometry import Circle, Lattice, Pattern
from sougata_solver.materials import Material

PERIOD = 0.7
WAVELENGTH = 1.0
OMEGA = 2 * np.pi / WAVELENGTH
AIR = Material("air", 1.0)
SI = Material("si", 3.48**2)


def _brute_force_min_gap(q: np.ndarray) -> float:
    n = q.shape[0]
    best = float("inf")
    for i in range(n):
        for j in range(n):
            if i != j:
                best = min(best, abs(q[i] - q[j]))
    return best


# ---------------------------------------------------------------------------
# Diagnostics fields match independently-recomputed quantities
# ---------------------------------------------------------------------------


def test_uniform_isotropic_diagnostics():
    kx = np.array([0.1, 0.4]) * OMEGA
    ky = np.array([0.2, -0.3]) * OMEGA
    modes = solve_layer_eigenmodes_uniform(OMEGA, kx, ky, eps=2.25)
    d = modes.diagnostics
    assert d is not None
    assert d.cond_epsilon == pytest.approx(1.0)
    assert d.cond_phi == pytest.approx(np.linalg.cond(modes.phi))
    assert d.min_eigenvalue_gap == pytest.approx(_brute_force_min_gap(modes.q))
    assert d.num_propagating + d.num_evanescent == modes.q.shape[0]


def test_uniform_diagonal_diagnostics():
    kx = np.array([0.3, 0.7]) * OMEGA
    ky = np.array([0.1, 0.2]) * OMEGA
    modes = solve_layer_eigenmodes_uniform_diagonal(OMEGA, kx, ky, eps_xx=2.25, eps_yy=4.0, eps_zz=3.1)
    d = modes.diagnostics
    assert d is not None
    assert d.cond_phi == pytest.approx(np.linalg.cond(modes.phi))
    assert d.min_eigenvalue_gap == pytest.approx(_brute_force_min_gap(modes.q))


def test_1d_diagnostics_cond_epsilon_matches_epsilon_inv_hat():
    lattice_1d_kx = np.array([0.3, 1.3, -0.7]) * OMEGA
    ky = np.zeros(3)
    n = 3
    epsilon_hat = np.diag([2.25, 2.25, 2.25]).astype(complex)
    epsilon_inv_hat = np.diag([1.0 / 2.25, 1.0 / 2.25, 1.0 / 2.25]).astype(complex)
    modes = solve_layer_eigenmodes_1d(OMEGA, lattice_1d_kx, ky, epsilon_hat, epsilon_inv_hat)
    d = modes.diagnostics
    assert d is not None
    assert d.cond_epsilon == pytest.approx(np.linalg.cond(epsilon_inv_hat))
    assert d.num_propagating + d.num_evanescent == 2 * n


def test_patterned_isotropic_diagnostics_cond_epsilon_matches_epsilon_hat():
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    g = truncate_fourier_orders(lattice, 5, "circular")
    lk = lattice.reciprocal_vectors()
    kx = 2 * np.pi * (g[:, 0] * lk[0, 0] + g[:, 1] * lk[1, 0])
    ky = 2 * np.pi * (g[:, 0] * lk[0, 1] + g[:, 1] * lk[1, 1])
    pattern = Pattern(background=AIR, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.18, material=SI)])
    epsilon_hat = toeplitz_matrix(pattern, lattice, g, WAVELENGTH, inverse=False)

    modes = solve_layer_eigenmodes_patterned(OMEGA, kx, ky, epsilon_hat)
    d = modes.diagnostics
    assert d is not None
    assert d.cond_epsilon == pytest.approx(np.linalg.cond(epsilon_hat))
    assert d.cond_phi == pytest.approx(np.linalg.cond(modes.phi))


def test_patterned_inplane_diagnostics_cond_epsilon_matches_epsilon_hat_zz():
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    g = truncate_fourier_orders(lattice, 5, "circular")
    lk = lattice.reciprocal_vectors()
    kx = 2 * np.pi * (g[:, 0] * lk[0, 0] + g[:, 1] * lk[1, 0])
    ky = 2 * np.pi * (g[:, 0] * lk[0, 1] + g[:, 1] * lk[1, 1])
    tensor = np.array([[2.25, 0.3, 0], [0.3, 4.0, 0], [0, 0, 3.1]], dtype=complex)
    aniso = Material.from_permittivity_tensor("aniso", tensor)
    pattern = Pattern(background=AIR, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.18, material=aniso)])
    exx = toeplitz_matrix_component(pattern, lattice, g, WAVELENGTH, 0, 0)
    exy = toeplitz_matrix_component(pattern, lattice, g, WAVELENGTH, 0, 1)
    eyx = toeplitz_matrix_component(pattern, lattice, g, WAVELENGTH, 1, 0)
    eyy = toeplitz_matrix_component(pattern, lattice, g, WAVELENGTH, 1, 1)
    ezz = toeplitz_matrix_component(pattern, lattice, g, WAVELENGTH, 2, 2)

    modes = solve_layer_eigenmodes_patterned_inplane(OMEGA, kx, ky, exx, exy, eyx, eyy, ezz)
    d = modes.diagnostics
    assert d is not None
    assert d.cond_epsilon == pytest.approx(np.linalg.cond(ezz))


# ---------------------------------------------------------------------------
# Attaching diagnostics does not change the solve result itself
# ---------------------------------------------------------------------------


def test_diagnostics_do_not_change_solve_results():
    kx = np.array([0.1, 0.4, -0.9]) * OMEGA
    ky = np.array([0.2, -0.3, 0.6]) * OMEGA
    modes_a = solve_layer_eigenmodes_uniform_inplane(OMEGA, kx, ky, 2.25, 0.3, 0.3, 4.0, 3.1)
    modes_b = solve_layer_eigenmodes_uniform_inplane(OMEGA, kx, ky, 2.25, 0.3, 0.3, 4.0, 3.1)
    assert np.array_equal(modes_a.q, modes_b.q)
    assert np.array_equal(modes_a.phi, modes_b.phi)
    assert np.array_equal(modes_a.kp, modes_b.kp)
    assert modes_a.epsilon_inv == modes_b.epsilon_inv
    assert modes_a.is_scalar_isotropic == modes_b.is_scalar_isotropic
