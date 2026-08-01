"""Independent 1D-grating RCWA oracle for Phase 3's system test.

Hand-transcribed (not imported, per `rules.md` AI Coding Rule 7) from
`REFERENCE/Rigorous-Coupled-Wave-Analysis/RCWA_1D_examples/`, which cites
Moharam, Grann, Pommet & Gaylord, "Formulation for stable and efficient
implementation of the rigorous coupled-wave analysis of binary gratings,"
J. Opt. Soc. Am. A 12(5), 1995, in its own header comment (`references.md`'s
planned Phase 3 oracle). Per `references.md`'s explicit instruction to "run
it directly rather than hand-copying paper tables" -- confirmed this
session that neither source file hard-codes any paper table numbers, they
only compute and plot a spectral sweep -- and following the precedent set
by `tests/oracles/empy_tmm.py` (transcribe, fix bugs found, don't reproduce
them), this module transcribes the physics only, with plotting/CLI code
(both source files interleave `matplotlib` calls with no `__main__` guard)
stripped out.

**TE** (`solve_te`): transcribed from `1D_Grating_Gaylord_TE.py:139-257`.
Uses Rumpf/Moharam's TE formulation (`A = KX2 - E`, eigenvalue problem for
the symmetric operator, boundary-matching via the `a`/`b` auxiliary
decomposition at lines 231-247) -- a different normalization/convention
from this project's own S4-derived `(q, phi, kp)` (see
`eigenmodes.solve_layer_eigenmodes_1d`'s docstring for why that
convention, not this one, is what `sougata_solver` implements) -- so this
is used purely as an independent cross-check oracle, not transcribed into
`src/`.

**TM** (`solve_tm`): transcribed from `1D_Grating_Gaylord_TM.py:169-246`
(the `A = inv(E_conv_inv) @ (KX @ solve(E, KX) - I)` formulation, general
`eig` since `A` is not Hermitian). One real, source-file-acknowledged
caveat: `1D_Grating_Gaylord_TM.py`'s own module docstring says "STILL NOT
WORKING YET" (line 14), and its hard-coded example additionally sets
`n_groove = 3.48` (same index as the ridge, i.e. not a real grating) rather
than `1.0` like the TE file -- this transcription fixes the latter (uses
the same binary-grating geometry as `solve_te`, ridge/groove passed as
parameters instead of hard-coded), but the former (the file's own
"not working yet" admission) is a real, unresolved caveat this project did
not independently re-derive or prove -- `test_1d_grating.py`'s TM
cross-check against this oracle is therefore run with a looser tolerance
and documented as secondary evidence, not the primary TM correctness
signal (that role is filled by the energy-conservation invariant and the
reduces-to-uniform-layer regression test, both of which exercise the exact
same `solve_layer_eigenmodes_1d` code path as TE).

Sign/amplitude convention note: both source files use `exp(-i*k*r)` for the
forward-propagating wave (module docstring, both files, line ~20/16) --
opposite time convention from some RCWA texts, but irrelevant here since
only diffraction-efficiency (power) ratios are compared, not raw fields or
phases, matching how `tests/oracles/empy_tmm.py` sidesteps the same kind of
sign-convention mismatch for Phase 1.
"""

from __future__ import annotations

import numpy as np


def _grating_fourier_harmonics(order: int, fill_factor: float, eps_ridge: complex, eps_groove: complex) -> complex:
    """`1D_Grating_Gaylord_TE.py:26-38`, generalized to accept permittivity
    directly (the source computes `n_ridge**2`/`n_groove**2` inline;
    passing `eps` instead avoids re-deriving that trivial step)."""
    if order == 0:
        return eps_ridge * fill_factor + eps_groove * (1 - fill_factor)
    return (eps_ridge - eps_groove) * np.sin(np.pi * order * fill_factor) / (np.pi * order)


def _toeplitz_convolution_matrix(num_ord: int, fill_factor: float, eps_ridge: complex, eps_groove: complex) -> np.ndarray:
    """`1D_Grating_Gaylord_TE.py:139-151`'s convolution-matrix assembly
    (`E[prow, pcol] = fourier_array[p0 + pfft]`), using the closed-form
    harmonics directly (`grating_fourier_harmonics`) instead of the
    source's FFT-of-rasterized-profile route (`grating_fft`,
    `fft_fourier_array`) -- both compute the same Fourier coefficients of a
    binary step profile; the closed-form route avoids the source's
    rasterization-resolution (`Nx`) as an extra free parameter."""
    n = 2 * num_ord + 1
    p_index = np.arange(-num_ord, num_ord + 1)
    matrix = np.zeros((n, n), dtype=complex)
    for prow in range(n):
        for pcol in range(n):
            dg = int(p_index[prow] - p_index[pcol])
            matrix[prow, pcol] = _grating_fourier_harmonics(dg, fill_factor, eps_ridge, eps_groove)
    return matrix


def solve_te(
    wavelength: float,
    theta: float,
    n_ridge: float,
    n_groove: float,
    fill_factor: float,
    lattice_constant: float,
    thickness: float,
    num_ord: int,
    n1: float = 1.0,
    n2: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-order TE diffraction efficiencies `(DE_r, DE_t)`, each length
    `2*num_ord+1`, ordered `[-num_ord, ..., num_ord]`.

    Transcribed from `1D_Grating_Gaylord_TE.py:160-257` (the per-wavelength
    loop body), unmodified in structure (variable names kept close to the
    source for line-by-line auditability).
    """
    indices = np.arange(-num_ord, num_ord + 1)
    eps_ridge = complex(n_ridge**2)
    eps_groove = complex(n_groove**2)
    identity = np.eye(2 * num_ord + 1, dtype=complex)
    e_conv = _toeplitz_convolution_matrix(num_ord, fill_factor, eps_ridge, eps_groove)

    j = 1j
    lam0 = wavelength
    k0 = 2 * np.pi / lam0

    kx_array = k0 * (n1 * np.sin(theta) + indices * (lam0 / lattice_constant))
    k_xi = kx_array
    kx2 = np.diag((k_xi / k0) ** 2)

    a_op = kx2 - e_conv
    eigenvals, w = np.linalg.eig(a_op)
    eigenvals = eigenvals.astype(complex)
    q_diag = np.sqrt(eigenvals)
    q_mat = np.diag(q_diag)
    v_mat = w @ q_mat
    x_mat = np.diag(np.exp(-k0 * q_diag * thickness))

    k_i = (k0**2 * (n1**2 - (k_xi / k0) ** 2)).astype(complex)
    k_ii = (k0**2 * (n2**2 - (k_xi / k0) ** 2)).astype(complex)
    k_i = np.sqrt(k_i)
    k_ii = np.sqrt(k_ii)
    y_i = np.diag(k_i / k0)
    y_ii = np.diag(k_ii / k0)

    delta_i0 = np.zeros((len(kx_array), 1), dtype=complex)
    delta_i0[num_ord] = 1
    n_delta_i0 = delta_i0 * j * n1 * np.cos(theta)

    w_inv = np.linalg.inv(w)
    v_inv = np.linalg.inv(v_mat)
    a_aux = 0.5 * (w_inv + j * v_inv @ y_ii)
    b_aux = 0.5 * (w_inv - j * v_inv @ y_ii)
    fbi_x = np.linalg.solve(b_aux, identity) @ x_mat

    term = x_mat @ a_aux @ fbi_x
    f = w @ (identity + term)
    g = v_mat @ (-identity + term)
    t = np.linalg.solve(j * y_i @ f + g, j * y_i @ delta_i0 + n_delta_i0)
    r = f @ t - delta_i0
    t = fbi_x @ t

    de_r = (r * np.conj(r) * np.real(k_i)[:, None] / (k0 * n1 * np.cos(theta))).real.ravel()
    de_t = (t * np.conj(t) * np.real(k_ii)[:, None] / (k0 * n1 * np.cos(theta))).real.ravel()
    return de_r, de_t


def solve_tm(
    wavelength: float,
    theta: float,
    n_ridge: float,
    n_groove: float,
    fill_factor: float,
    lattice_constant: float,
    thickness: float,
    num_ord: int,
    n1: float = 1.0,
    n2: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-order TM diffraction efficiencies `(DE_r, DE_t)`. See module
    docstring: `1D_Grating_Gaylord_TM.py` self-documents as "STILL NOT
    WORKING YET" -- this transcription (from lines 169-246) fixes the
    groove-material transcription artifact but not that underlying,
    unresolved caveat; treat as secondary evidence only.
    """
    indices = np.arange(-num_ord, num_ord + 1)
    eps_ridge = complex(n_ridge**2)
    eps_groove = complex(n_groove**2)
    identity = np.eye(2 * num_ord + 1, dtype=complex)
    e_conv = _toeplitz_convolution_matrix(num_ord, fill_factor, eps_ridge, eps_groove)
    e_conv_inv = _toeplitz_convolution_matrix(num_ord, fill_factor, 1.0 / eps_ridge, 1.0 / eps_groove)

    j = 1j
    lam0 = wavelength
    k0 = 2 * np.pi / lam0

    kx_array = k0 * (n1 * np.sin(theta) + indices * (lam0 / lattice_constant))
    k_xi = kx_array
    kx_diag = np.diag(k_xi / k0)

    a_op = np.linalg.inv(e_conv_inv) @ (kx_diag @ np.linalg.solve(e_conv, kx_diag) - identity)

    eigenvals, w = np.linalg.eig(a_op)
    eigenvals = eigenvals.astype(complex)
    q_diag = np.sqrt(eigenvals)
    q_mat = np.diag(q_diag)
    v_mat = e_conv_inv @ (w @ q_mat)
    x_mat = np.diag(np.exp(-k0 * q_diag * thickness))

    k_i = (k0**2 * (n1**2 - (k_xi / k0) ** 2)).astype(complex)
    k_ii = (k0**2 * (n2**2 - (k_xi / k0) ** 2)).astype(complex)
    k_i = np.sqrt(k_i)
    k_ii = np.sqrt(k_ii)
    z_i = np.diag(k_i / (n1**2 * k0))
    z_ii = np.diag(k_ii / (n2**2 * k0))

    delta_i0 = np.zeros((len(kx_array), 1), dtype=complex)
    delta_i0[num_ord] = 1
    n_delta_i0 = delta_i0 * j * np.cos(theta) / n1

    n2q = 2 * num_ord + 1
    o_mat = np.block([[w, w], [v_mat, -v_mat]])
    fg = np.concatenate([identity, j * z_ii], axis=0)
    ab = np.linalg.solve(o_mat, fg)
    a_aux = ab[:n2q, :]
    b_aux = ab[n2q:, :]

    term = x_mat @ a_aux @ np.linalg.solve(b_aux, x_mat)
    f = w @ (identity + term)
    g = v_mat @ (-identity + term)
    t = np.linalg.solve(j * z_i @ f + g, j * z_i @ delta_i0 + n_delta_i0)
    r = f @ t - delta_i0
    t = np.linalg.solve(b_aux, x_mat) @ t

    de_r = (r * np.conj(r) * np.real(k_i)[:, None] / (k0 * n1 * np.cos(theta))).real.ravel()
    de_t = (t * np.conj(t) * np.real(k_ii)[:, None] / n2**2 / (k0 * np.cos(theta) / n1)).real.ravel()
    return de_r, de_t
