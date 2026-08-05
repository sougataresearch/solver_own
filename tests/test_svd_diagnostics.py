"""Category 12 target 12.4 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
`eigenmodes.svd_diagnostics` -- an opt-in (never automatically called
during a solve), more detailed singular-value diagnostic than the
condition-number-only `EigenmodeDiagnostics.cond_phi` (Category 2 target
2.2) that every dense anisotropic/patterned solver already attaches.
"""

from __future__ import annotations

import numpy as np
import pytest

from sougata_solver.eigenmodes import SVDDiagnostics, solve_layer_eigenmodes_patterned, svd_diagnostics
from sougata_solver.fourier_basis import truncate_fourier_orders
from sougata_solver.fourier_factorization import toeplitz_matrix
from sougata_solver.geometry import Circle, Lattice, Pattern
from sougata_solver.materials import Material

PERIOD = 0.7
AIR = Material("air", 1.0)


def test_svd_diagnostics_identity_matrix():
    result = svd_diagnostics(np.eye(4, dtype=complex))
    assert isinstance(result, SVDDiagnostics)
    np.testing.assert_allclose(result.singular_values, np.ones(4))
    assert result.condition_number == pytest.approx(1.0)
    assert result.num_small_singular_values == 0


def test_svd_diagnostics_singular_values_are_descending():
    rng = np.random.default_rng(0)
    a = rng.normal(size=(6, 6)) + 1j * rng.normal(size=(6, 6))
    result = svd_diagnostics(a)
    assert np.all(result.singular_values[:-1] >= result.singular_values[1:])


def test_svd_diagnostics_condition_number_matches_numpy_cond():
    rng = np.random.default_rng(1)
    a = rng.normal(size=(5, 5)) + 1j * rng.normal(size=(5, 5))
    result = svd_diagnostics(a)
    assert result.condition_number == pytest.approx(np.linalg.cond(a), rel=1e-8)


def test_svd_diagnostics_detects_near_rank_deficient_matrix():
    """A matrix with one near-null column must show exactly one small
    singular value, not just a large condition number -- the specific
    added value over `cond_phi` alone (target 12.4's own "how many modes"
    framing)."""
    rng = np.random.default_rng(2)
    a = rng.normal(size=(5, 5)) + 1j * rng.normal(size=(5, 5))
    a[:, -1] = a[:, 0] * 1e-10  # nearly linearly dependent column
    result = svd_diagnostics(a, relative_threshold=1e-6)
    assert result.num_small_singular_values == 1
    assert result.condition_number > 1e6


def test_svd_diagnostics_relative_threshold_is_configurable():
    rng = np.random.default_rng(3)
    a = rng.normal(size=(5, 5)) + 1j * rng.normal(size=(5, 5))
    a[:, -1] = a[:, 0] * 1e-3
    loose = svd_diagnostics(a, relative_threshold=1e-2)
    tight = svd_diagnostics(a, relative_threshold=1e-6)
    assert loose.num_small_singular_values >= tight.num_small_singular_values


def test_svd_diagnostics_on_ill_conditioned_pillar_eigenvector_matrix():
    """Reuses Phase 4b's near-touching-pillar stress fixture
    (`tests/test_2d_pillar_stress.py`), the most ill-conditioned `phi`
    already characterized in this project -- confirms `svd_diagnostics`
    runs on a genuine solver-produced eigenvector matrix, not just
    synthetic random matrices, and that its `condition_number` agrees
    with the already-established `cond_phi` diagnostic on the same
    matrix."""
    wavelength = 1.0
    omega = 2 * np.pi / wavelength
    si = Material("stress", 3.48**2)
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.49 * PERIOD, material=si)])
    g = truncate_fourier_orders(lattice, 25, "circular")
    lk = lattice.reciprocal_vectors()
    kx = 2 * np.pi * (g[:, 0] * lk[0, 0] + g[:, 1] * lk[1, 0])
    ky = 2 * np.pi * (g[:, 0] * lk[0, 1] + g[:, 1] * lk[1, 1])
    epsilon_hat = toeplitz_matrix(pattern, lattice, g, wavelength, inverse=False)

    modes = solve_layer_eigenmodes_patterned(omega, kx, ky, epsilon_hat)
    assert modes.diagnostics is not None

    result = svd_diagnostics(modes.phi)
    assert result.condition_number == pytest.approx(modes.diagnostics.cond_phi, rel=1e-6)
    assert result.singular_values.shape == (modes.phi.shape[0],)
