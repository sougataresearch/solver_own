"""Category 2 target 2.4 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): configurable
small-eigenvalue-gap `WARNING` (`eigenmodes.DEGENERATE_GAP_THRESHOLD`,
`eigenmodes._warn_on_small_eigenvalue_gap`), applied to the three
anisotropic dense eigensolvers that already carry the target 1.7 canonical-
ordering policy -- **not** `solve_layer_eigenmodes_patterned` (Phase 4a,
isotropic 2D), which was found during this target's own testing to have
routine, harmless `C4v`-symmetry degeneracy for an ordinary well-conditioned
case; see that function's docstring for the full account of why it's
excluded.

Same "detect, don't silently correct" pattern, and the same test structure
(`caplog`/`monkeypatch`), as `tests/test_2d_pillar_stress.py`'s
`ILL_CONDITIONED_THRESHOLD` tests -- both warning-fires and no-warning
paths are required per `rules.md` Testing Requirements.
"""

from __future__ import annotations

import logging

import numpy as np

from sougata_solver import eigenmodes
from sougata_solver.eigenmodes import (
    solve_layer_eigenmodes_patterned_inplane,
    solve_layer_eigenmodes_uniform_diagonal,
    solve_layer_eigenmodes_uniform_inplane,
)
from sougata_solver.fourier_basis import truncate_fourier_orders
from sougata_solver.fourier_factorization import toeplitz_matrix_component
from sougata_solver.geometry import Circle, Lattice, Pattern
from sougata_solver.materials import Material

PERIOD = 0.7
WAVELENGTH = 1.0
OMEGA = 2 * np.pi / WAVELENGTH
AIR = Material("air", 1.0)


# ---------------------------------------------------------------------------
# No warning for an ordinary, clearly-non-degenerate case
# ---------------------------------------------------------------------------


def test_no_warning_for_well_separated_eigenvalues(caplog):
    kx = np.array([0.1, 0.4, -0.9]) * OMEGA
    ky = np.array([0.2, -0.3, 0.6]) * OMEGA
    with caplog.at_level(logging.WARNING, logger="sougata_solver.eigenmodes"):
        solve_layer_eigenmodes_uniform_diagonal(OMEGA, kx, ky, eps_xx=2.25, eps_yy=4.0, eps_zz=3.1)
    assert [r for r in caplog.records if "near-degenerate" in r.getMessage()] == []


def test_no_warning_for_well_separated_eigenvalues_inplane(caplog):
    kx = np.array([0.1, 0.4, -0.9]) * OMEGA
    ky = np.array([0.2, -0.3, 0.6]) * OMEGA
    with caplog.at_level(logging.WARNING, logger="sougata_solver.eigenmodes"):
        solve_layer_eigenmodes_uniform_inplane(OMEGA, kx, ky, 2.25, 0.3, 0.3, 4.0, 3.1)
    assert [r for r in caplog.records if "near-degenerate" in r.getMessage()] == []


# ---------------------------------------------------------------------------
# Warning fires for a deliberately near-isotropic (near-degenerate) case
# ---------------------------------------------------------------------------


def test_warning_fires_for_near_isotropic_uniform_diagonal(caplog):
    """`eps_xx` and `eps_yy` almost equal -> the two normal-incidence
    branches (`q_x^2=eps_xx*omega^2`, `q_y^2=eps_yy*omega^2`, per
    `solve_layer_eigenmodes_uniform_diagonal`'s own closed-form derivation)
    sit very close together, a genuine near-degeneracy, not a forced
    threshold like the `monkeypatch` test below."""
    kx = np.array([0.0])
    ky = np.array([0.0])
    with caplog.at_level(logging.WARNING, logger="sougata_solver.eigenmodes"):
        solve_layer_eigenmodes_uniform_diagonal(OMEGA, kx, ky, eps_xx=2.25, eps_yy=2.25 + 1e-10, eps_zz=3.1)
    messages = [r.getMessage() for r in caplog.records]
    assert any("near-degenerate eigenvalues detected" in m for m in messages)


def test_warning_fires_when_threshold_is_raised_above_an_ordinary_gap(caplog, monkeypatch):
    """Forces the WARNING path by raising `DEGENERATE_GAP_THRESHOLD` far
    above an ordinary, otherwise well-separated case's actual gap --
    confirms the logging mechanism itself works (message, logger name,
    level), independent of whether a real case happens to be
    near-degenerate enough to trigger it naturally."""
    monkeypatch.setattr(eigenmodes, "DEGENERATE_GAP_THRESHOLD", 1e9)

    kx = np.array([0.1, 0.4, -0.9]) * OMEGA
    ky = np.array([0.2, -0.3, 0.6]) * OMEGA
    with caplog.at_level(logging.WARNING, logger="sougata_solver.eigenmodes"):
        solve_layer_eigenmodes_uniform_inplane(OMEGA, kx, ky, 2.25, 0.3, 0.3, 4.0, 3.1)
    messages = [r.getMessage() for r in caplog.records]
    assert any("near-degenerate eigenvalues detected" in m for m in messages)


def test_patterned_inplane_warning_fires_for_near_isotropic_case(caplog):
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    g = truncate_fourier_orders(lattice, 5, "circular")
    lk = lattice.reciprocal_vectors()
    kx = 2 * np.pi * (g[:, 0] * lk[0, 0] + g[:, 1] * lk[1, 0])
    ky = 2 * np.pi * (g[:, 0] * lk[0, 1] + g[:, 1] * lk[1, 1])
    tensor = np.array([[2.25 + 1e-10, 0.0, 0], [0.0, 2.25, 0], [0, 0, 3.1]], dtype=complex)  # near-isotropic
    aniso = Material.from_permittivity_tensor("near_iso", tensor)
    pattern = Pattern(background=AIR, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.18, material=aniso)])
    exx = toeplitz_matrix_component(pattern, lattice, g, WAVELENGTH, 0, 0)
    exy = toeplitz_matrix_component(pattern, lattice, g, WAVELENGTH, 0, 1)
    eyx = toeplitz_matrix_component(pattern, lattice, g, WAVELENGTH, 1, 0)
    eyy = toeplitz_matrix_component(pattern, lattice, g, WAVELENGTH, 1, 1)
    ezz = toeplitz_matrix_component(pattern, lattice, g, WAVELENGTH, 2, 2)

    with caplog.at_level(logging.WARNING, logger="sougata_solver.eigenmodes"):
        solve_layer_eigenmodes_patterned_inplane(OMEGA, kx, ky, exx, exy, eyx, eyy, ezz)
    messages = [r.getMessage() for r in caplog.records]
    assert any("near-degenerate eigenvalues detected" in m for m in messages)
