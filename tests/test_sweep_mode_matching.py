"""Category 2 target 2.3 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): deterministic
mode ordering across a small wavelength sweep, for the three anisotropic
dense eigensolvers that already apply `eigenmodes._canonical_mode_order`
(target 1.7): `solve_layer_eigenmodes_uniform_diagonal`,
`solve_layer_eigenmodes_uniform_inplane`, `solve_layer_eigenmodes_patterned_inplane`.

`solve_layer_eigenmodes_patterned` (Phase 4a, isotropic 2D) is deliberately
**not** included -- extending canonical ordering to it was tried and
reverted, since it broke `tests/test_2d_pillar.py`'s block-structure
regression tests; see that function's own docstring ("Category 2 target 2.3
-- evaluated and deliberately NOT applied here") for the full account.

What "does not arbitrarily permute" means, tested here: for a small,
non-degenerate wavelength step, the canonically-ordered `q` at each sweep
point, compared index-by-index to the previous point, changes by an amount
consistent with the physical eigenvalue's own continuous variation (small
relative to the smallest inter-mode gap), not a discontinuous jump that
would indicate two modes swapped positions in the ordering. This is
explicitly *not* a claim of continuity through an eigenvalue crossing/
degeneracy (`_canonical_mode_order`'s docstring already states that caveat)
-- the sweep below is deliberately chosen non-degenerate throughout (no
`eps_xx == eps_yy` crossing) so no crossing is expected.
"""

from __future__ import annotations

import numpy as np

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
AIR = Material("air", 1.0)

WAVELENGTHS = np.linspace(1.00, 1.02, 6)  # small, close-together sweep points


def _kx_ky(omega: float) -> tuple[np.ndarray, np.ndarray]:
    kx = np.array([0.1, 0.4, -0.6]) * omega
    ky = np.array([0.2, -0.3, 0.5]) * omega
    return kx, ky


def _assert_no_arbitrary_permutation(q_by_step: list[np.ndarray]) -> None:
    """A permutation-free sweep has each mode's trajectory changing
    smoothly step-to-step; a step where the canonical order swapped two
    modes shows up as a jump comparable to the *gap between* those modes
    rather than the sweep's own small parameter step. Checked here by:
    for every step, the max index-wise change is much smaller than the
    minimum inter-mode gap at that step (a permutation would instead
    produce a jump on the order of that gap)."""
    for i in range(1, len(q_by_step)):
        prev, curr = q_by_step[i - 1], q_by_step[i]
        step_change = np.max(np.abs(curr - prev))
        diffs = np.abs(curr[:, None] - curr[None, :])
        np.fill_diagonal(diffs, np.inf)
        min_gap = np.min(diffs)
        assert step_change < 0.5 * min_gap, (
            f"sweep step {i}: index-wise q change ({step_change}) is not small "
            f"relative to the inter-mode gap ({min_gap}) -- looks like an "
            "arbitrary mode-order permutation, not smooth variation"
        )


# ---------------------------------------------------------------------------
# Uniform diagonal / in-plane solvers
# ---------------------------------------------------------------------------


def test_uniform_diagonal_mode_order_does_not_arbitrarily_permute_across_sweep():
    q_by_step = []
    for wavelength in WAVELENGTHS:
        omega = 2 * np.pi / wavelength
        kx, ky = _kx_ky(omega)
        modes = solve_layer_eigenmodes_uniform_diagonal(omega, kx, ky, eps_xx=2.25, eps_yy=4.0, eps_zz=3.1)
        q_by_step.append(modes.q)
    _assert_no_arbitrary_permutation(q_by_step)


def test_uniform_inplane_mode_order_does_not_arbitrarily_permute_across_sweep():
    q_by_step = []
    for wavelength in WAVELENGTHS:
        omega = 2 * np.pi / wavelength
        kx, ky = _kx_ky(omega)
        modes = solve_layer_eigenmodes_uniform_inplane(omega, kx, ky, 2.25, 0.3, 0.3, 4.0, 3.1)
        q_by_step.append(modes.q)
    _assert_no_arbitrary_permutation(q_by_step)


def test_patterned_inplane_mode_order_does_not_arbitrarily_permute_across_sweep():
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    g = truncate_fourier_orders(lattice, 5, "circular")
    lk = lattice.reciprocal_vectors()
    tensor = np.array([[2.25, 0.3, 0], [0.3, 4.0, 0], [0, 0, 3.1]], dtype=complex)
    aniso = Material.from_permittivity_tensor("aniso", tensor)
    pattern = Pattern(background=AIR, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.18, material=aniso)])

    q_by_step = []
    for wavelength in WAVELENGTHS:
        omega = 2 * np.pi / wavelength
        kx = 2 * np.pi * (g[:, 0] * lk[0, 0] + g[:, 1] * lk[1, 0])
        ky = 2 * np.pi * (g[:, 0] * lk[0, 1] + g[:, 1] * lk[1, 1])
        exx = toeplitz_matrix_component(pattern, lattice, g, wavelength, 0, 0)
        exy = toeplitz_matrix_component(pattern, lattice, g, wavelength, 0, 1)
        eyx = toeplitz_matrix_component(pattern, lattice, g, wavelength, 1, 0)
        eyy = toeplitz_matrix_component(pattern, lattice, g, wavelength, 1, 1)
        ezz = toeplitz_matrix_component(pattern, lattice, g, wavelength, 2, 2)
        modes = solve_layer_eigenmodes_patterned_inplane(omega, kx, ky, exx, exy, eyx, eyy, ezz)
        q_by_step.append(modes.q)
    _assert_no_arbitrary_permutation(q_by_step)


# ---------------------------------------------------------------------------
# Repeated identical-input solves within a "sweep" (step size 0) stay
# bit-identical -- the trivial case of "does not arbitrarily permute".
# ---------------------------------------------------------------------------


def test_repeated_identical_step_is_bit_identical():
    omega = 2 * np.pi / 1.0
    kx, ky = _kx_ky(omega)
    modes_a = solve_layer_eigenmodes_uniform_diagonal(omega, kx, ky, eps_xx=2.25, eps_yy=4.0, eps_zz=3.1)
    modes_b = solve_layer_eigenmodes_uniform_diagonal(omega, kx, ky, eps_xx=2.25, eps_yy=4.0, eps_zz=3.1)
    assert np.array_equal(modes_a.q, modes_b.q)
