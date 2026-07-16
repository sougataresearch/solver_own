"""Interface matrices and S-matrix (scattering matrix) cascading.

Interface-matrix derivation verified directly against `S4/S4/rcwa.cpp`
(`GetSMatrix`, lines 936-1096, particularly the `P`/`Q`/`Ta`/`Tb` comment
block at 944-957). The transfer-to-scattering conversion and star-product
recursion are transcribed directly from `S4/S4r/StarProduct.hpp`
(`T2Sblocks` lines 51-65, `StarProduct` lines 83-110), which is a standard
Redheffer star product for optical scattering matrices.

All S-matrices here are stored as a single `(2*n2, 2*n2)` block array
`[[S00, S01], [S10, S11]]` with a *fixed* block size `n2 = 2*num_orders`
throughout a simulation (every layer shares the same Fourier order count),
satisfying the convention::

    [a_right; b_left] = S @ [a_left; b_right]

i.e. `a` is the forward-going amplitude vector, `b` the backward-going one,
and `S` maps the two *incoming* amplitudes (from the left and right) to the
two *outgoing* ones.
"""

from __future__ import annotations

import numpy as np
import scipy.linalg as sla

from sougata_solver.layer import LayerEigenmodes


def _solve(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Robust `inv(a) @ b` via LU factorization (never form `inv(a)` alone)."""
    lu, piv = sla.lu_factor(a)
    return sla.lu_solve((lu, piv), b)


def _is_trivial_interface(modes_l: LayerEigenmodes, modes_lp1: LayerEigenmodes) -> bool:
    return (
        modes_l is modes_lp1
        or (
            np.array_equal(modes_l.q, modes_lp1.q)
            and np.array_equal(modes_l.phi, modes_lp1.phi)
            and np.array_equal(modes_l.kp, modes_lp1.kp)
        )
    )


def interface_smatrix(modes_l: LayerEigenmodes, modes_lp1: LayerEigenmodes) -> np.ndarray:
    """Scattering matrix of the interface between layer `l` and `l+1`.

    Transfer-style relation (from field continuity, `rcwa.cpp:944-957`)::

        A = phi,  B = kp @ phi @ diag(1/q)
        P = inv(A_l) @ A_{l+1}
        Q = inv(B_l) @ B_{l+1}
        Ta = 0.5*(P+Q), Tb = 0.5*(P-Q)
        [a_l; b_l] = [[Ta, Tb], [Tb, Ta]] @ [a_{l+1}; b_{l+1}]

    Converted to scattering form via `T2Sblocks` (`StarProduct.hpp:51-65`)::

        S00 = inv(Ta)
        S10 = Tb @ S00
        S01 = -S00 @ Tb
        S11 = Ta + Tb @ S01
    """
    n2 = modes_l.q.shape[0]

    if _is_trivial_interface(modes_l, modes_lp1):
        s = np.zeros((2 * n2, 2 * n2), dtype=complex)
        s[:n2, :n2] = np.eye(n2)
        s[n2:, n2:] = np.eye(n2)
        return s

    a_l, a_lp1 = modes_l.phi, modes_lp1.phi
    b_l = modes_l.kp @ modes_l.phi / modes_l.q[None, :]
    b_lp1 = modes_lp1.kp @ modes_lp1.phi / modes_lp1.q[None, :]

    p = _solve(a_l, a_lp1)
    q = _solve(b_l, b_lp1)
    ta = 0.5 * (p + q)
    tb = 0.5 * (p - q)

    identity = np.eye(n2, dtype=complex)
    s00 = _solve(ta, identity)
    s10 = tb @ s00
    s01 = -s00 @ tb
    s11 = ta + tb @ s01

    s = np.zeros((2 * n2, 2 * n2), dtype=complex)
    s[:n2, :n2] = s00
    s[:n2, n2:] = s01
    s[n2:, :n2] = s10
    s[n2:, n2:] = s11
    return s


def propagation_smatrix(q: np.ndarray, thickness: float) -> np.ndarray:
    """Scattering matrix for free propagation through a layer of finite
    `thickness`: a pure phase delay, no forward/backward mode mixing.

        a_out = diag(exp(i*q*thickness)) @ a_in
        b_out = diag(exp(i*q*thickness)) @ b_in
        => S00 = S11 = diag(exp(i*q*thickness)), S01 = S10 = 0

    (`q` was branch-selected with `Im(q) >= 0`, so this phase factor decays
    for a forward wave traversing a lossy/evanescent layer of thickness > 0.)
    """
    n2 = q.shape[0]
    phase = np.exp(1j * q * thickness)
    s = np.zeros((2 * n2, 2 * n2), dtype=complex)
    idx = np.arange(n2)
    s[idx, idx] = phase
    s[idx + n2, idx + n2] = phase
    return s


def star_product(n2: int, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Cascade two `(2*n2, 2*n2)` scattering matrices `A` (left interface)
    and `B` (right interface) into their combined scattering matrix,
    following the Redheffer star product (`StarProduct.hpp:83-110`), with a
    fixed block size `n2` throughout (every layer has the same order count)::

        C00 = B00 @ inv(I - A01@B10) @ A00
        C01 = B00 @ inv(I - A01@B10) @ A01@B11 + B01
        C10 = A10 + A11 @ inv(I - B10@A01) @ B10@A00
        C11 =        A11 @ inv(I - B10@A01) @ B11
    """
    a00, a01 = a[:n2, :n2], a[:n2, n2:]
    a10, a11 = a[n2:, :n2], a[n2:, n2:]
    b00, b01 = b[:n2, :n2], b[:n2, n2:]
    b10, b11 = b[n2:, :n2], b[n2:, n2:]
    identity = np.eye(n2, dtype=complex)

    t1 = identity - a01 @ b10
    c00 = b00 @ _solve(t1, a00)
    c01 = b00 @ _solve(t1, a01 @ b11) + b01

    t2 = identity - b10 @ a01
    c10 = a10 + a11 @ _solve(t2, b10 @ a00)
    c11 = a11 @ _solve(t2, b11)

    c = np.zeros((2 * n2, 2 * n2), dtype=complex)
    c[:n2, :n2] = c00
    c[:n2, n2:] = c01
    c[n2:, :n2] = c10
    c[n2:, n2:] = c11
    return c


class SMatrixStack:
    """Builds and caches the full S-matrix cascade for a layer stack at one
    wavelength. Layer `0` and `len(modes)-1` are the semi-infinite
    incidence/transmission half-spaces (no propagation phase); every layer
    in between contributes a propagation S-matrix over its own thickness.
    """

    def __init__(self, thicknesses: list[float], all_modes: list[LayerEigenmodes]):
        if len(thicknesses) != len(all_modes):
            raise ValueError("thicknesses and all_modes must have the same length")
        n2 = all_modes[0].q.shape[0]
        self.n2 = n2

        cumulative = interface_smatrix(all_modes[0], all_modes[1])
        self._partial = [cumulative]
        for i in range(1, len(all_modes) - 1):
            prop = propagation_smatrix(all_modes[i].q, thicknesses[i])
            cumulative = star_product(n2, cumulative, prop)
            iface = interface_smatrix(all_modes[i], all_modes[i + 1])
            cumulative = star_product(n2, cumulative, iface)
            self._partial.append(cumulative)

    def full_smatrix(self) -> np.ndarray:
        return self._partial[-1]

    def partial_smatrix_up_to(self, layer_index: int) -> np.ndarray:
        """S-matrix of layers `[0, layer_index]` only (0-indexed into the
        original `all_modes` list), used for field reconstruction at
        intermediate depths without re-cascading from scratch."""
        return self._partial[layer_index - 1]
