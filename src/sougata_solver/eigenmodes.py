"""Per-layer eigenmode solver.

Formulas verified directly against `S4/S4/rcwa.cpp` (not paraphrased):
`MakeKPMatrix` (lines 242-293), `SolveLayerEigensystem_uniform`
(lines 422-502), and the general eigen-operator construction in
`SolveLayerEigensystem` (lines 794-827).
"""

from __future__ import annotations

import numpy as np

from sougata_solver.layer import LayerEigenmodes


def build_kp_matrix(omega: complex, kx: np.ndarray, ky: np.ndarray, epsilon_inv) -> np.ndarray:
    """Build the k-parallel operator:

        kp = omega^2 * I_2n - V @ Einv_block @ V^T

    where `V = [[diag(ky)], [-diag(kx)]]` (2n x n). `epsilon_inv` is either
    a scalar (isotropic uniform layer, `Einv_block = epsilon_inv * I_n`) or
    an `(n, n)` Fourier-space matrix (patterned/anisotropic layer).

    Source: `rcwa.cpp::MakeKPMatrix`, lines 242-293.
    """
    kx = np.asarray(kx, dtype=complex)
    ky = np.asarray(ky, dtype=complex)
    n = kx.shape[0]
    n2 = 2 * n
    idx = np.arange(n)

    kappa = np.zeros((n2, n2), dtype=complex)
    if np.ndim(epsilon_inv) == 0:
        epsinv = complex(epsilon_inv)
        kappa[idx, idx] = ky * epsinv * ky
        kappa[idx, idx + n] = -ky * epsinv * kx
        kappa[idx + n, idx] = -kx * epsinv * ky
        kappa[idx + n, idx + n] = kx * epsinv * kx
    else:
        einv = np.asarray(epsilon_inv, dtype=complex)
        ky_diag = np.diag(ky)
        kx_diag = np.diag(kx)
        kappa[:n, :n] = ky_diag @ einv @ ky_diag
        kappa[:n, n:] = -ky_diag @ einv @ kx_diag
        kappa[n:, :n] = -kx_diag @ einv @ ky_diag
        kappa[n:, n:] = kx_diag @ einv @ kx_diag

    kp = -kappa
    kp[idx, idx] += omega**2
    kp[idx + n, idx + n] += omega**2
    return kp


def _select_q_branch(q_sq: np.ndarray, tol: float = 4 * np.finfo(float).eps) -> np.ndarray:
    """Select the outgoing/decaying branch of `q = sqrt(q_sq)`, matching
    S4's real-frequency convention (`rcwa.cpp:455-467` / `847-861`):
    near-real `q_sq` is handled specially to keep purely-evanescent modes
    exactly on the positive-imaginary axis; otherwise the principal branch
    is flipped to have `Im(q) >= 0`.
    """
    q_sq = np.asarray(q_sq, dtype=complex)
    q = np.empty_like(q_sq)

    near_real = np.abs(q_sq.imag) <= tol * np.abs(q_sq.real)
    real_part = q_sq.real

    pos = near_real & (real_part >= 0)
    neg = near_real & (real_part < 0)
    q[pos] = np.sqrt(real_part[pos])
    q[neg] = 1j * np.sqrt(-real_part[neg])

    other = ~near_real
    q_other = np.sqrt(q_sq[other])
    flip = q_other.imag < 0
    q_other[flip] = -q_other[flip]
    q[other] = q_other
    return q


def solve_layer_eigenmodes_uniform(omega: complex, kx: np.ndarray, ky: np.ndarray, eps: complex) -> LayerEigenmodes:
    """Closed-form eigenmode solve for a uniform isotropic layer.

    `q[i] = q[i+n] = branch_select(eps*omega^2 - kx[i]^2 - ky[i]^2)`,
    `phi = I_2n` (the eigenbasis coincides with the plane-wave field basis
    for a homogeneous isotropic medium). Verified algebraically consistent
    with the general eigen-operator (`op = Epsilon2 @ kp - U@U^T` reduces
    to `(eps*omega^2 - kx^2 - ky^2) * I_2n` per order when `Epsilon2 = eps*I`).

    Source: `rcwa.cpp::SolveLayerEigensystem_uniform`, lines 422-502.
    """
    kx = np.asarray(kx, dtype=complex)
    ky = np.asarray(ky, dtype=complex)
    n = kx.shape[0]
    n2 = 2 * n

    q_sq = eps * omega**2 - kx**2 - ky**2
    q_half = _select_q_branch(q_sq)
    q = np.concatenate([q_half, q_half])

    phi = np.eye(n2, dtype=complex)
    kp = build_kp_matrix(omega, kx, ky, 1.0 / eps)

    return LayerEigenmodes(q=q, phi=phi, kp=kp, epsilon_inv=None, is_scalar_isotropic=True)
