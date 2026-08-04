"""Category 3 target 3.1 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): pins
`design.md`'s Fourier-Factorization Rule Inventory table against actual
solver behavior, so a future refactor that silently changes which matrix a
solver factorizes/inverts shows up as a test failure here, not just as a
stale table.

Two complementary check styles:
- black-box: `LayerEigenmodes.epsilon_inv` (already a public field) directly
  reveals which matrix each solver treats as "the" permittivity-inverse for
  Fourier-factorization purposes, without reimplementing any solver
  internals.
- white-box (1D only, where `epsilon_inv` alone can't distinguish the
  TE/TM blocks' different `Epsilon2` source): reproduce the two candidate
  `Epsilon2` block formulas directly from `build_kp_matrix` + the raw
  Toeplitz inputs (the same technique `test_1d_grating.py`'s
  `test_op_is_block_diagonal_when_ky_zero` already uses for a different
  invariant) and confirm the TE block matches the direct-rule formula, not
  the inverse-rule one it could have been confused with.
"""

from __future__ import annotations

import numpy as np
import pytest

from sougata_solver.eigenmodes import (
    build_kp_matrix,
    solve_layer_eigenmodes_1d,
    solve_layer_eigenmodes_patterned,
    solve_layer_eigenmodes_patterned_inplane,
    solve_layer_eigenmodes_uniform,
    solve_layer_eigenmodes_uniform_diagonal,
)


def _hermitian_positive_definite(rng, n: int) -> np.ndarray:
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return a @ a.conj().T + n * np.eye(n)


# ---------------------------------------------------------------------------
# Uniform solvers: exact, no Fourier factorization at all
# ---------------------------------------------------------------------------


def test_uniform_isotropic_has_no_fourier_factorized_matrix():
    modes = solve_layer_eigenmodes_uniform(2 * np.pi, np.array([0.3]), np.array([0.1]), eps=2.25)
    assert modes.epsilon_inv is None


def test_uniform_diagonal_epsilon_inv_is_exact_scalar_not_a_toeplitz():
    modes = solve_layer_eigenmodes_uniform_diagonal(
        2 * np.pi, np.array([0.3]), np.array([0.1]), eps_xx=2.25, eps_yy=4.0, eps_zz=3.1
    )
    assert modes.epsilon_inv == pytest.approx(1.0 / 3.1)


# ---------------------------------------------------------------------------
# 1D solver: TE block uses the direct-rule Toeplitz, TM block the
# (numerically re-inverted) inverse-rule Toeplitz -- the one place in the
# project epsilon_inv_hat is actually consumed as such.
# ---------------------------------------------------------------------------


def test_1d_solver_stores_the_inverse_rule_toeplitz_itself_as_epsilon_inv(rng):
    n = 3
    epsilon_hat = _hermitian_positive_definite(rng, n)
    epsilon_inv_hat = _hermitian_positive_definite(rng, n)
    modes = solve_layer_eigenmodes_1d(
        2 * np.pi, rng.normal(size=n), np.zeros(n), epsilon_hat, epsilon_inv_hat
    )
    # Returned unchanged -- not further inverted -- unlike the 2D solvers below.
    assert np.array_equal(modes.epsilon_inv, epsilon_inv_hat)


def test_1d_te_block_uses_direct_rule_tm_block_uses_inverse_rule(rng):
    n = 4
    omega = 3.7
    kx = rng.normal(size=n)
    ky = np.zeros(n)
    epsilon_hat = _hermitian_positive_definite(rng, n)
    epsilon_inv_hat = _hermitian_positive_definite(rng, n)

    modes = solve_layer_eigenmodes_1d(omega, kx, ky, epsilon_hat, epsilon_inv_hat)
    kp = build_kp_matrix(omega, kx, ky, epsilon_inv_hat)

    op_te_direct_rule = epsilon_hat @ kp[:n, :n] - np.diag(kx * kx)
    op_tm_inverse_rule = np.linalg.solve(epsilon_inv_hat, np.eye(n, dtype=complex)) @ kp[n:, n:]

    q_te_expected = np.sort_complex(np.linalg.eigvals(op_te_direct_rule))
    q_tm_expected = np.sort_complex(np.linalg.eigvals(op_tm_inverse_rule))
    q_te_actual = np.sort_complex(modes.q[:n] ** 2)
    q_tm_actual = np.sort_complex(modes.q[n:] ** 2)

    assert q_te_actual == pytest.approx(q_te_expected, abs=1e-8)
    assert q_tm_actual == pytest.approx(q_tm_expected, abs=1e-8)

    # The test must actually discriminate between the two rules, not just
    # reproduce solve_layer_eigenmodes_1d's own formula tautologically:
    # confirm the TE block would give a *different* answer under the wrong
    # (inverse) rule, for this random input.
    op_te_wrong_rule = epsilon_inv_hat @ kp[:n, :n] - np.diag(kx * kx)
    q_te_wrong = np.sort_complex(np.linalg.eigvals(op_te_wrong_rule))
    assert not np.allclose(q_te_actual, q_te_wrong)


# ---------------------------------------------------------------------------
# 2D solvers: both use a numerical matrix-inverse of the *direct*-rule
# Toeplitz, never a separately-factorized inverse-rule Toeplitz.
# ---------------------------------------------------------------------------


def test_2d_isotropic_epsilon_inv_is_numerical_inverse_of_direct_rule(rng):
    n = 3
    epsilon_hat = _hermitian_positive_definite(rng, n)
    kx = rng.normal(size=n)
    ky = rng.normal(size=n)
    modes = solve_layer_eigenmodes_patterned(2 * np.pi, kx, ky, epsilon_hat)
    assert modes.epsilon_inv == pytest.approx(np.linalg.inv(epsilon_hat))


def test_2d_anisotropic_epsilon_inv_is_numerical_inverse_of_direct_rule_ezz(rng):
    n = 3
    epsilon_hat_xx = _hermitian_positive_definite(rng, n)
    epsilon_hat_yy = _hermitian_positive_definite(rng, n)
    epsilon_hat_zz = _hermitian_positive_definite(rng, n)
    zero = np.zeros((n, n), dtype=complex)
    kx = rng.normal(size=n)
    ky = rng.normal(size=n)
    modes = solve_layer_eigenmodes_patterned_inplane(
        2 * np.pi, kx, ky, epsilon_hat_xx, zero, zero, epsilon_hat_yy, epsilon_hat_zz
    )
    assert modes.epsilon_inv == pytest.approx(np.linalg.inv(epsilon_hat_zz))
