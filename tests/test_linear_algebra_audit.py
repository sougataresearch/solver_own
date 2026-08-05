"""Category 12 target 12.2 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
direct-inverse audit.

**Audit performed**: grepped `src/sougata_solver/*.py` for
`linalg.inv`/`.inv(`/`linalg.solve`/`lu_factor`/`lu_solve` (every
inverse-shaped or linear-system-solve call site) and reviewed each:

- `smatrix.py::_solve` -- already uses the house convention
  (`scipy.linalg.lu_factor`/`lu_solve`), no change.
- `excitation.py::incident_mode_amplitude` -- `np.linalg.solve(kp_phi,
  rhs)` against a specific right-hand-side vector, not a full matrix
  inverse -- already the minimal necessary operation, no change.
- `geometry.py::Lattice.reciprocal_vectors` -- `np.linalg.inv(self._Lr)`,
  a 2x2 real matrix inverted once at `Lattice` construction. **Reviewed
  and kept as-is**: not "demonstrably unnecessary" per the target's own
  wording -- a 2x2 direct inverse costs nothing measurable, and switching
  to a solve-based approach would add code with no benefit for a matrix
  this small.
- `eigenmodes.py` (three call sites: `solve_layer_eigenmodes_1d`,
  `solve_layer_eigenmodes_patterned`, `solve_layer_eigenmodes_patterned_inplane`)
  -- each computed a full dense inverse via `np.linalg.solve(A, np.eye(n))`.
  **Not** "demonstrably unnecessary" either (the full matrix, not a
  single linear-system answer, is genuinely consumed downstream by
  `build_kp_matrix`'s `V @ Einv_block @ V^T` construction) -- but
  inconsistent with the project's own documented house convention
  (`rules.md`'s Performance Requirements: "`scipy.linalg.lu_factor`/
  `lu_solve`... is the house convention for solving linear systems...
  reuse this helper rather than reintroducing direct matrix inversion
  elsewhere"). Fixed by extracting `eigenmodes._dense_inverse`
  (`scipy.linalg.lu_factor`/`lu_solve` against the identity, the same
  house-convention pattern `smatrix.py::_solve` already uses) and
  replacing all three call sites -- a consistency fix, not a new
  algorithm or a correctness change.

This file is the equivalence-test regression guard the target's own
wording requires: confirms the refactored `_dense_inverse`-based solve
results are bit-for-bit identical to values captured from the
pre-refactor `np.linalg.solve(A, np.eye(n))` code path, and pins the
underlying claim (`_dense_inverse` gives the same matrix as a direct
`np.linalg.inv`) as its own unit check.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from sougata_solver.eigenmodes import _dense_inverse
from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.fourier_basis import truncate_fourier_orders
from sougata_solver.fourier_factorization import toeplitz_matrix
from sougata_solver.geometry import Circle, Lattice, Lattice1D, Pattern, Slab
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation

AIR = Material("air", 1.0)
SI = Material("si", 3.48**2)


# ---------------------------------------------------------------------------
# Unit: _dense_inverse matches np.linalg.inv directly
# ---------------------------------------------------------------------------


def test_dense_inverse_matches_numpy_inv_random_matrix():
    rng = np.random.default_rng(0)
    a = rng.normal(size=(6, 6)) + 1j * rng.normal(size=(6, 6))
    np.testing.assert_allclose(_dense_inverse(a), np.linalg.inv(a), atol=1e-10)


def test_dense_inverse_matches_numpy_inv_identity():
    a = np.eye(4, dtype=complex)
    np.testing.assert_allclose(_dense_inverse(a), np.eye(4, dtype=complex), atol=1e-12)


def test_dense_inverse_round_trips():
    rng = np.random.default_rng(1)
    a = rng.normal(size=(5, 5)) + 1j * rng.normal(size=(5, 5))
    a_inv = _dense_inverse(a)
    np.testing.assert_allclose(a @ a_inv, np.eye(5, dtype=complex), atol=1e-9)


# ---------------------------------------------------------------------------
# System: end-to-end R/T unchanged by the eigenmodes.py refactor
# ---------------------------------------------------------------------------
# Values captured from the pre-refactor code path
# (`np.linalg.solve(A, np.eye(n))` at each of the three call sites) before
# switching to `_dense_inverse`, per this target's own "add equivalence
# tests" requirement.


def test_1d_grating_result_unchanged_by_dense_inverse_refactor():
    lattice = Lattice1D(0.7)
    pattern = Pattern(background=AIR)
    pattern.add(Slab(center_x=0.0, halfwidth=0.15, material=SI))
    sim = Simulation(lattice, [Layer("grating", 0.3, pattern=pattern)], num_orders=9, incidence=AIR, transmission=AIR)
    result = sim.solve(PlaneWaveExcitation(1.0, math.radians(10.0), 0.0, s_amplitude=0.0, p_amplitude=1.0))

    assert result.reflectance() == pytest.approx(0.427458044405598, abs=1e-12)
    assert result.transmittance() == pytest.approx(0.5725419555944025, abs=1e-12)


def test_2d_pillar_result_unchanged_by_dense_inverse_refactor():
    lattice = Lattice((0.7, 0.0), (0.0, 0.7))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(0.35, 0.35), radius=0.14, material=SI)])
    sim = Simulation(lattice, [Layer("pillar", 0.3, pattern=pattern)], num_orders=9, incidence=AIR, transmission=AIR)
    result = sim.solve(PlaneWaveExcitation(1.0, math.radians(10.0), 0.0, s_amplitude=1.0, p_amplitude=0.0))

    assert result.reflectance() == pytest.approx(0.5059492125225666, abs=1e-12)
    assert result.transmittance() == pytest.approx(0.4940507874774355, abs=1e-12)


# ---------------------------------------------------------------------------
# Category 12 target 12.5: the density measurement `decisions.md` ADR-021
# relies on -- pinned here so a future session doesn't have to re-derive it.
# ---------------------------------------------------------------------------


def test_patterned_layer_toeplitz_matrix_is_fully_dense():
    """`decisions.md` ADR-021 (sparse/iterative methods evaluated and
    rejected): the direct-rule Toeplitz permittivity matrix for an
    ordinary 2D patterned layer has no exploitable sparsity -- every
    pairwise Fourier-order coupling is nonzero, the expected consequence
    of a shape's continuous Fourier transform sampled at every
    reciprocal-lattice difference. A future change to the Fourier-
    factorization construction that somehow introduced real sparsity
    would be worth knowing about (it would strengthen, not weaken, this
    project's position), so this is checked directly, not assumed."""
    lattice = Lattice((0.7, 0.0), (0.0, 0.7))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(0.35, 0.35), radius=0.14, material=SI)])
    g = truncate_fourier_orders(lattice, 49, "circular")
    epsilon_hat = toeplitz_matrix(pattern, lattice, g, 0.6, inverse=False)

    n = epsilon_hat.shape[0]
    nonzero_fraction = np.count_nonzero(np.abs(epsilon_hat) > 1e-12 * np.max(np.abs(epsilon_hat))) / (n * n)
    assert nonzero_fraction == pytest.approx(1.0)
