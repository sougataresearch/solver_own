"""Poynting flux and (later) real-space field reconstruction.

`z_poynting_flux` is transcribed directly from `S4/S4/rcwa.cpp::GetZPoyntingFlux`
(lines 1846-1897), not re-derived, since a from-scratch re-derivation of the
sign/normalization conventions embedded in the `kp`/`phi` operators risked
introducing exactly the kind of subtle error this module needs to avoid.
Correctness is then checked empirically in Phase 1 against an independent
analytic Fresnel/TMM oracle.
"""

from __future__ import annotations

import numpy as np


def z_poynting_flux(
    omega: complex,
    q: np.ndarray,
    kp: np.ndarray,
    phi: np.ndarray,
    avec: np.ndarray,
    bvec: np.ndarray,
) -> tuple[complex, complex]:
    """Time-averaged z-Poynting flux carried by forward amplitudes `avec`
    and backward amplitudes `bvec` at one reference plane (summed over all
    Fourier orders). Returns `(forward, backward)`, both nominally real for
    propagating, lossless configurations; take `.real` when reporting power.

    Source: `rcwa.cpp::GetZPoyntingFlux`, lines 1846-1897.
    """
    a2 = avec / (omega * q)
    b2 = bvec / (omega * q)
    a3 = phi @ a2
    b3 = phi @ b2
    ka = kp @ a3
    kb = kp @ b3
    alpha = phi.conj().T @ ka
    beta = phi.conj().T @ kb

    forward = np.vdot(avec, alpha).real
    backward = -np.vdot(bvec, beta).real
    diff = 0.5 * (np.vdot(bvec, alpha) - np.vdot(beta, avec))
    forward = forward + diff
    backward = backward + np.conj(diff)
    return forward, backward


def tangential_e_field(
    omega: complex,
    q: np.ndarray,
    kp: np.ndarray,
    phi: np.ndarray,
    avec: np.ndarray,
    bvec: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Tangential `(Ex, Ey)`, each length `n = num_orders`, from mode
    amplitudes `avec` (forward), `bvec` (backward).

    Transcribed from the E-field half of `rcwa.cpp::GetInPlaneFieldVector`
    (lines 1959-1995) -- note this is *not* the naive `E = phi @ (a+b)`;
    that quantity is actually `H`. `E` uses `(a-b)` with an index swap and
    a sign flip::

        u = kp @ phi @ (avec - bvec) / (omega * q)
        Ex = u[n:2n]
        Ey = -u[0:n]
    """
    n2 = avec.shape[0]
    n = n2 // 2
    u = kp @ (phi @ ((avec - bvec) / (omega * q)))
    ex = u[n:]
    ey = -u[:n]
    return ex, ey
