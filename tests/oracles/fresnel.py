"""Independent N-layer Fresnel/TMM reference (characteristic-matrix method).

Written from scratch from the standard Born & Wolf / Macleod thin-film
formulation, not derived from EMpy or sougata_solver, so it serves as an
independent oracle for `test_analytic_fresnel.py`.
"""

from __future__ import annotations

import cmath
import math


def _kz(n: complex, k0: float, kx: float) -> complex:
    """z-wavenumber with the outgoing/decaying branch (`Im(kz) >= 0`)."""
    val = (n * k0) ** 2 - kx**2
    kz = cmath.sqrt(val)
    if kz.imag < 0:
        kz = -kz
    return kz


def _admittance(n: complex, kz: complex, k0: float, polarization: str) -> complex:
    costheta = kz / (n * k0)
    if polarization == "s":
        return n * costheta
    elif polarization == "p":
        return n / costheta
    raise ValueError(f"polarization must be 's' or 'p', got {polarization!r}")


def multilayer_rt(
    wavelength: float,
    theta0: float,
    polarization: str,
    n_incidence: complex,
    layers: list[tuple[complex, float]],
    n_substrate: complex,
) -> tuple[float, float]:
    """Reflectance and transmittance of an N-layer isotropic stack.

    `layers` is an ordered list of `(refractive_index, thickness)` for the
    internal layers, from the incidence side to the substrate side.
    """
    k0 = 2.0 * math.pi / wavelength
    n_incidence = complex(n_incidence)
    n_substrate = complex(n_substrate)
    kx = (n_incidence * k0 * math.sin(theta0)).real

    kz0 = _kz(n_incidence, k0, kx)
    eta0 = _admittance(n_incidence, kz0, k0, polarization)

    kz_sub = _kz(n_substrate, k0, kx)
    eta_sub = _admittance(n_substrate, kz_sub, k0, polarization)

    m = [[1.0 + 0j, 0j], [0j, 1.0 + 0j]]
    for n_j, d_j in layers:
        n_j = complex(n_j)
        kz_j = _kz(n_j, k0, kx)
        eta_j = _admittance(n_j, kz_j, k0, polarization)
        delta = kz_j * d_j
        cos_d = cmath.cos(delta)
        sin_d = cmath.sin(delta)
        # Sign verified by direct derivation from E(z)=E+*exp(i*kz*z)+E-*exp(-i*kz*z)
        # (matching sougata_solver's forward-wave convention): solving for the front-face
        # (E1,H1) in terms of the back-face (E2,H2) gives -i*sin(delta), not +i.
        # This only matters when delta is complex (absorbing layers) -- with the
        # wrong sign, lossless cases still pass (real delta hides the error) while
        # every absorbing-layer case silently violates R+T<=1.
        mj = [[cos_d, -1j * sin_d / eta_j], [-1j * eta_j * sin_d, cos_d]]
        m = [
            [m[0][0] * mj[0][0] + m[0][1] * mj[1][0], m[0][0] * mj[0][1] + m[0][1] * mj[1][1]],
            [m[1][0] * mj[0][0] + m[1][1] * mj[1][0], m[1][0] * mj[0][1] + m[1][1] * mj[1][1]],
        ]

    b = m[0][0] * 1.0 + m[0][1] * eta_sub
    c = m[1][0] * 1.0 + m[1][1] * eta_sub

    r = (eta0 * b - c) / (eta0 * b + c)
    t = 2.0 * eta0 / (eta0 * b + c)

    reflectance = abs(r) ** 2
    transmittance = (eta_sub.real / eta0.real) * abs(t) ** 2
    return reflectance, transmittance
