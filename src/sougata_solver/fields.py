"""Poynting flux and real-space field reconstruction.

`z_poynting_flux` is transcribed directly from `S4/S4/rcwa.cpp::GetZPoyntingFlux`
(lines 1846-1897), not re-derived, since a from-scratch re-derivation of the
sign/normalization conventions embedded in the `kp`/`phi` operators risked
introducing exactly the kind of subtle error this module needs to avoid.
Correctness is then checked empirically in Phase 1 against an independent
analytic Fresnel/TMM oracle.

Category 9 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`, Phase 7) adds real-space
E/H reconstruction (`modal_field_components`, `propagate_amplitudes`,
`reconstruct_field_at_points`) and NumPy field-grid export
(`save_field_grid_npz`) -- see `CONVENTIONS.md`'s "Real-space field
reconstruction" section for the full formula citations and conventions
(including a genuine finding: `z_poynting_flux` is missing the textbook
`0.5` time-average factor, harmless for R/T ratios but relevant for an
absolute real-space flux integral).
"""

from __future__ import annotations

from pathlib import Path

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


def modal_field_components(
    omega: complex,
    kx: np.ndarray,
    ky: np.ndarray,
    q: np.ndarray,
    kp: np.ndarray,
    phi: np.ndarray,
    epsilon_inv,
    avec: np.ndarray,
    bvec: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Category 9 target 9.1/9.2 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`, Phase
    7): full 6-component per-Fourier-order field `(Ex, Ey, Ez, Hx, Hy, Hz)`
    at one z-reference plane, from mode amplitudes `avec` (forward), `bvec`
    (backward). Each returned array has length `n = num_orders`.

    Transverse components transcribed from
    `S4/S4/rcwa.cpp::GetInPlaneFieldVector` (lines 1959-1995) -- `Ex`/`Ey`
    reuse `tangential_e_field` (already transcribed from the same source);
    `Hx`/`Hy` are the *other* half of that same source function, the
    quantity `tangential_e_field`'s own docstring already warns looks like
    (but is not) `E`::

        Hx, Hy = split(phi @ (avec + bvec))

    Longitudinal components transcribed from the same function's caller,
    `GetFieldAtPoint` (lines 1997-2074), which derives them from the
    source-free Maxwell curl equations restricted to their z-component
    (`del -> i*(kx, ky, d/dz)` per Fourier order)::

        Ez = epsilon_inv @ (ky*Hx - kx*Hy) / omega
        Hz = (kx*Ey - ky*Ex) / omega

    `epsilon_inv` is either a scalar (uniform isotropic layer) or an
    `(n, n)` Fourier-space matrix (patterned/anisotropic layer) -- same
    scalar-or-matrix dispatch as `eigenmodes.build_kp_matrix`.
    """
    ex, ey = tangential_e_field(omega, q, kp, phi, avec, bvec)
    h = phi @ (avec + bvec)
    n = h.shape[0] // 2
    hx, hy = h[:n], h[n:]

    if np.ndim(epsilon_inv) == 0:
        ez = complex(epsilon_inv) * (ky * hx - kx * hy) / omega
    else:
        ez = (epsilon_inv @ (ky * hx - kx * hy)) / omega
    hz = (kx * ey - ky * ex) / omega

    return ex, ey, ez, hx, hy, hz


def propagate_amplitudes(q: np.ndarray, z: float, a_top: np.ndarray, b_top: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Category 9 target 9.1 (Phase 7): forward/backward mode amplitudes at
    depth `z` (`0 <= z <= thickness`) within a layer, given the amplitudes
    `(a_top, b_top)` at that layer's top reference plane (`z=0` local to
    the layer -- exactly what `smatrix.interior_amplitudes` recovers).

    Independently derived (see `CONVENTIONS.md`'s "Real-space field
    reconstruction" section and `decisions.md` ADR-015), algebraically
    consistent with `smatrix.propagation_smatrix`'s already-established
    `exp(+i*q*thickness)` convention (confirmed directly, not assumed --
    `tests/test_field_reconstruction.py` checks that propagating a
    layer's top amplitudes to `z=thickness` reproduces the amplitudes the
    existing `propagation_smatrix`-based cascade already computes for the
    next interface)::

        a(z) = a_top * exp(+i*q*z)
        b(z) = b_top * exp(-i*q*z)
    """
    phase_fwd = np.exp(1j * q * z)
    phase_bwd = np.exp(-1j * q * z)
    return a_top * phase_fwd, b_top * phase_bwd


def reconstruct_field_at_points(
    kx: np.ndarray,
    ky: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    field_component: np.ndarray,
) -> np.ndarray:
    """Category 9 targets 9.2/9.4/9.7 (Phase 7): inverse-Fourier-sum one
    real-space field component at point(s) `(x, y)` from its per-order
    Fourier coefficients `field_component` (length `n = num_orders`, one of
    the six arrays `modal_field_components` returns).

    Dimension-agnostic (serves a single point (target 9.2), a 1D line
    (target 9.4), or a 2D grid (target 9.7) alike -- `x`/`y` broadcast
    against the `n`-length Fourier-order axis via a trailing `[..., None]`,
    so any `x`/`y` shape works, not three separate reconstruction
    functions): transcribed from `S4/S4/rcwa.cpp::GetFieldAtPoint`'s phase
    sum (lines 2044-2050)::

        F(x, y) = sum_i F_i * exp(i*(kx_i*x + ky_i*y))

    `kx`, `ky` are this project's already-angular per-order wavevectors
    (the same arrays `eigenmodes.py`/`simulation.py` use) -- see
    `CONVENTIONS.md` for why this is a *different* convention from
    `geometry.py`'s "cycles per unit length" `kx`/`ky` used only for shape
    Fourier transforms.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    phase = np.exp(1j * (kx * x[..., None] + ky * y[..., None]))
    return phase @ field_component


def save_field_grid_npz(
    path: str,
    x: np.ndarray,
    y: np.ndarray,
    z: float,
    ex: np.ndarray,
    ey: np.ndarray,
    ez: np.ndarray,
    hx: np.ndarray,
    hy: np.ndarray,
    hz: np.ndarray,
    **metadata: object,
) -> Path:
    """Category 9 target 9.8 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`, Phase
    7): save a reconstructed real-space field grid to a `.npz` file via
    `numpy.savez` -- lossless, no schema-design risk (unlike CSV, which
    would need a flattening/column convention decided first, or HDF5,
    which is a new dependency), so this is the safe "NumPy first" half of
    this target's own wording. CSV/HDF5 export is deliberately **not**
    added here -- left open pending the schema design that target's
    wording calls for, not silently done ad hoc.

    `x`/`y` are the real-space coordinate arrays the grid was evaluated at
    (any shape `reconstruct_field_at_points` accepts), `z` the single
    depth the grid was reconstructed at, `ex..hz` the six field-component
    arrays (same shape as `x`/`y`). Any additional `metadata` keyword
    arguments (wavelength, angle, material names, ...) are stored
    alongside the arrays, so a saved file is self-describing without a
    companion `run_metadata.txt` (unlike `output_paths.write_run_metadata`,
    which is for scalar run parameters, not array data).
    """
    path_obj = Path(path)
    np.savez(
        path_obj,
        x=x,
        y=y,
        z=np.asarray(z),
        Ex=ex,
        Ey=ey,
        Ez=ez,
        Hx=hx,
        Hy=hy,
        Hz=hz,
        **{k: np.asarray(v) for k, v in metadata.items()},
    )
    # np.savez appends ".npz" unless the given path already ends with it.
    return path_obj if str(path_obj).endswith(".npz") else Path(str(path_obj) + ".npz")
