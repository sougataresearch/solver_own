"""Second independent N-layer TMM oracle, transcribed from the vendored
`EMpy` reference implementation (`EMpy/EMpy/transfer_matrix.py`,
`IsotropicTransferMatrix.solve`, lines 52-134, plus the `snell` helper at
`EMpy/EMpy/utils.py:1006-1021`) -- not imported (`EMpy` is reference-only
per `decisions.md`, never a runtime dependency), transcribed by hand.

Chosen over the other vendored TMM sources (`Rigorous-Coupled-Wave-
Analysis/TMM_functions`, which mixes plotting/prints into the physics and
has an anisotropic function that references undefined free variables) as
the best-quality standalone isotropic-TMM code available: a clean classic
Abeles dynamical-matrix (characteristic-matrix) formulation, no framework
dependencies, one function per concern.

This is deliberately a *second, differently-sourced* oracle alongside
`fresnel.py` (written from scratch from Born & Wolf/Macleod) -- agreement
between two independently-derived implementations is stronger evidence than
either alone.

Three intentional deviations from the EMpy source, all found by running this
transcription against real (absorbing) materials during development -- i.e.
this cross-check caught real bugs in the source it's transcribed from:

1. The transmittance prefactor there is
   `abs((n_last*cos(theta_last)) / (n0*cos(theta0)))` (EMpy
   `transfer_matrix.py:121,126`). The physically correct quantity is the
   *real part* of that ratio (it's a Poynting-flux ratio), not its modulus
   -- see `fresnel.py`'s own note on this exact issue.
2. EMpy computes each layer's angle via a `snell()` `arcsin` chain
   (`utils.py:1006-1021`) and takes `cos()` of it to build the propagation
   phase. For a complex (absorbing) index, `cmath.asin` does not reliably
   return the branch corresponding to the physically decaying wave. Fixed
   here by computing each layer's `kz` directly as a branch-selected
   complex square root (`Im(kz) >= 0`, the same convention `fresnel.py` and
   sougata_solver's own eigenmode solver use) instead of chaining
   `arcsin`/`cos`.
3. The dynamical-matrix product `D @ P @ inv(D)` has a sign ambiguity in
   which of the two `exp(+-i*phi)` propagation eigenvalues pairs with which
   of `D`'s two columns -- get it backwards and lossless cases still pass
   (real `phi` hides the sign error) while every absorbing-layer case
   produces unphysical R>1. This is the *exact* subtlety `fresnel.py`'s own
   comment documents ("front-face in terms of back-face gives -i*sin(delta),
   not +i"). Verified by direct numerical comparison against `fresnel.py`
   for single- and multi-layer absorbing cases: `prop =
   diag(exp(-i*phi), exp(+i*phi))`, not `diag(exp(+i*phi), exp(-i*phi))`,
   is the convention consistent with front-face-in-terms-of-back-face.
"""

from __future__ import annotations

import cmath
import math


def _kz_branch(n: complex, k0: float, kx: float) -> complex:
    """z-wavenumber with the outgoing/decaying branch (`Im(kz) >= 0`) --
    see point 2 in the module docstring."""
    kz = cmath.sqrt((n * k0) ** 2 - kx**2)
    if kz.imag < 0:
        kz = -kz
    return kz


def isotropic_multilayer_rt(
    wavelength: float,
    theta_inc: float,
    n_incidence: complex,
    layers: list[tuple[complex, float]],
    n_substrate: complex,
    polarization: str,
) -> tuple[float, float]:
    """Reflectance and transmittance of an N-layer isotropic stack at one
    wavelength/angle, `polarization` in `{"s", "p"}`.

    `layers` is an ordered list of `(refractive_index, thickness)` for the
    internal layers, incidence side to substrate side -- same calling
    convention as `fresnel.multilayer_rt`, for direct comparison.
    """
    n = [complex(n_incidence)] + [complex(nj) for nj, _ in layers] + [complex(n_substrate)]
    k0 = 2.0 * math.pi / wavelength
    # kx is real and conserved across all layers (phase matching / Snell's
    # law in k-space); each layer's kz is then a branch-selected complex
    # sqrt -- see point 2 in the module docstring.
    kx = (n[0] * k0 * math.sin(theta_inc)).real
    kz = [_kz_branch(nj, k0, kx) for nj in n]
    cos_theta = [kzj / (nj * k0) for kzj, nj in zip(kz, n)]

    def dynamical_matrix(nj: complex, cj: complex) -> list[list[complex]]:
        if polarization == "s":
            return [[1.0, 1.0], [nj * cj, -nj * cj]]
        elif polarization == "p":
            return [[cj, cj], [nj, -nj]]
        raise ValueError(f"polarization must be 's' or 'p', got {polarization!r}")

    def mat2(a, b):
        return [
            [a[0][0] * b[0][0] + a[0][1] * b[1][0], a[0][0] * b[0][1] + a[0][1] * b[1][1]],
            [a[1][0] * b[0][0] + a[1][1] * b[1][0], a[1][0] * b[0][1] + a[1][1] * b[1][1]],
        ]

    def inv2(m):
        det = m[0][0] * m[1][1] - m[0][1] * m[1][0]
        return [[m[1][1] / det, -m[0][1] / det], [-m[1][0] / det, m[0][0] / det]]

    d0 = dynamical_matrix(n[0], cos_theta[0])
    m = inv2(d0)

    for i, (nj, dj) in enumerate(layers, start=1):
        nj = complex(nj)
        dj_mat = dynamical_matrix(nj, cos_theta[i])
        phi = kz[i] * dj
        # -i/+i, not +i/-i -- see point 3 in the module docstring.
        prop = [[cmath.exp(-1j * phi), 0j], [0j, cmath.exp(1j * phi)]]
        m = mat2(m, mat2(dj_mat, mat2(prop, inv2(dj_mat))))

    d_last = dynamical_matrix(n[-1], cos_theta[-1])
    m = mat2(m, d_last)

    r = m[1][0] / m[0][0]
    t = 1.0 / m[0][0]

    reflectance = abs(r) ** 2
    # Re(...) ratio, not abs(...) -- see module docstring.
    prefactor = ((n[-1] * cos_theta[-1]) / (n[0] * cos_theta[0])).real
    transmittance = prefactor * abs(t) ** 2
    return reflectance, transmittance
