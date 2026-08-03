"""Independent eigenoperator cross-check for Category 1 target 1.4
(`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): uniform in-plane-coupled
anisotropic layers, hand-transcribed from `RigorousCoupledWaveAnalysis.jl`
(Julia) -- not imported (Julia isn't installed in this environment,
re-confirmed this session: `which julia` fails), and not run as a
subprocess. Same approach as Phase 4a's
`tests/oracles/rcwa_2djl_eigenvalues.py`.

Source: `RigorousCoupledWaveAnalysis.jl/src/Common/Common.jl:134-165`
(`eigenmodes(dnx,dny,Kx,Ky,λ,l::AnisotropicLayer)`), the *uniform*
(unpatterned) anisotropic-layer branch -- note this is a different
Common.jl function from the one Phase 4a's oracle transcribes
(`l::PatternedLayer`, lines 57-99); target 1.4 is a uniform-layer
capability, so this is the matching source. Transcribed formula
(`eta = 1/eps_zz`)::

    P = [[Kx@eta@Ky,       I - Kx@eta@Kx],
         [Ky@eta@Ky - I,   -Ky@eta@Kx]]
    Q = [[Kx@Ky + eps_yx,  eps_yy - Kx@Kx],
         [Ky@Ky - eps_xx,  -eps_xy - Ky@Kx]]
    M = P @ Q
    eigenvalues = eig(M).values          # these are q^2 (unbranched)

This is the same **structurally different** derivation route documented in
`rcwa_2djl_eigenvalues.py`'s docstring (direct Maxwell-curl elimination
into one matrix `M`, vs. this project's `Epsilon2 @ kp - coupling`) --
eigenvectors are not expected to match and are not compared; only the
basis-independent `q^2` eigenvalue set is.

Honesty note (per `rules.md` AI Coding Rule 5): as with the Phase 4a
oracle, no independently-published numeric ground truth was found in
`RigorousCoupledWaveAnalysis.jl` for this case -- this remains a
cross-implementation eigenoperator check, not a literature benchmark.
"""

from __future__ import annotations

import numpy as np


def eigenoperator_eigenvalues_inplane(
    omega: complex,
    kx: np.ndarray,
    ky: np.ndarray,
    eps_xx: complex,
    eps_xy: complex,
    eps_yx: complex,
    eps_yy: complex,
    eps_zz: complex,
) -> np.ndarray:
    """`q^2` eigenvalues of RCWA.jl's uniform anisotropic-layer operator
    `M` (`Common.jl:150-153`), for a given per-order `(kx, ky)` at
    frequency `omega`.

    Conventions reconciled the same way as `rcwa_2djl_eigenvalues.py`
    (both re-confirmed empirically for this function specifically, not
    assumed to carry over automatically just because the isotropic case
    needed them): RCWA.jl's `Kx`/`Ky` are `kx/k0`/`ky/k0` (dimensionless),
    and its `q^2` comes out an overall sign-flip of this project's
    `q**2` (opposite time convention) -- both already known from the
    isotropic oracle.

    A **third** convention difference is specific to this in-plane-coupled
    case (there was no way to discover it from the diagonal-only case
    alone, since it's invisible when `eps_xy=eps_yx=0`): this project's
    `eps_xx`/`eps_yy` and RCWA.jl's `eps_xx`/`eps_yy` refer to the *same*
    physical tensor components, but the two solvers' internal field/axis
    bookkeeping differ in a way that surfaces as needing `kx` and `ky`
    swapped *and* `eps_xy`/`eps_yx` negated when calling this oracle to
    match `solve_layer_eigenmodes_uniform_inplane`'s `q**2` --
    i.e. call as `eigenoperator_eigenvalues_inplane(omega, ky, kx, eps_xx,
    -eps_xy, -eps_yx, eps_yy, eps_zz)`, not with `(kx, ky, eps_xx, eps_xy,
    eps_yx, eps_yy, eps_zz)` directly. This was **not derived from S4's
    `abcde` convention** (that only explained the diagonal-term swap seen
    in target 1.3's own docstring) -- it was determined by a brute-force
    search over swap/sign hypotheses (`kx<->ky`, `eps_xx<->eps_yy`,
    `eps_xy<->eps_yx`, and independent sign flips on each), checked against
    20 random-parameter trials plus the diagonal-only reduction, landing on
    this one combination at ~1e-13 residual and no other combination tried
    coming close -- an honest "found empirically, not from first
    principles" note, per `rules.md`'s standard for documenting a
    convention reconciliation. See
    `tests/test_anisotropic_inplane.py` for the exact call convention used.
    """
    kx_n = np.asarray(kx, dtype=complex) / omega
    ky_n = np.asarray(ky, dtype=complex) / omega
    n = kx_n.shape[0]
    eps_xx = complex(eps_xx)
    eps_xy = complex(eps_xy)
    eps_yx = complex(eps_yx)
    eps_yy = complex(eps_yy)
    eta = 1.0 / complex(eps_zz)

    kx_diag = np.diag(kx_n)
    ky_diag = np.diag(ky_n)
    eye_n = np.eye(n, dtype=complex)

    p = np.zeros((2 * n, 2 * n), dtype=complex)
    p[:n, :n] = eta * (kx_diag @ ky_diag)
    p[:n, n:] = eye_n - eta * (kx_diag @ kx_diag)
    p[n:, :n] = eta * (ky_diag @ ky_diag) - eye_n
    p[n:, n:] = -eta * (ky_diag @ kx_diag)

    q_op = np.zeros((2 * n, 2 * n), dtype=complex)
    q_op[:n, :n] = kx_diag @ ky_diag + eps_yx * eye_n
    q_op[:n, n:] = eps_yy * eye_n - kx_diag @ kx_diag
    q_op[n:, :n] = ky_diag @ ky_diag - eps_xx * eye_n
    q_op[n:, n:] = -eps_xy * eye_n - ky_diag @ kx_diag

    m = p @ q_op
    return -np.linalg.eigvals(m) * omega**2
