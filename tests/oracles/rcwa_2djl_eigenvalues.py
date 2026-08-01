"""Independent eigenoperator cross-check for Phase 4a's 2D patterned-layer
solver, hand-transcribed from `RigorousCoupledWaveAnalysis.jl` (Julia) --
not imported (Julia isn't installed in this environment; `which julia`
fails), and not run as a subprocess.

This is deliberately narrow in scope, and that scope is the point: it
checks only the **eigenoperator construction**, not a full R/T pipeline,
by feeding it this project's own already-Phase-2-validated `epsilon_hat`
Toeplitz matrix rather than re-deriving a shape's Fourier transform too.
That isolates exactly the step that broke in the bug this oracle exists to
guard against (`eigenmodes.solve_layer_eigenmodes_patterned`'s `Epsilon2`
construction was copied from the wrong branch of S4's source -- see that
function's docstring for the full account).

Source: `RigorousCoupledWaveAnalysis.jl/src/Common/Common.jl:57-99`
(`eigenmodes(dnx,dny,Kx,Ky,λ,l::PatternedLayer)`), isotropic branch only
(`εxy=εyx=0`, `εxx=εyy=εzz`). Transcribed formula::

    A = inv(E) @ (Kx @ E)
    B = inv(E) @ (Ky @ E)
    M = [[Ky^2 - E + Kx@A,  -Ky@Kx + Kx@B],
         [-Kx@Ky + Ky@A,     Kx^2 - E + Ky@B]]
    eigenvalues = eig(M).values          # these are q^2 (unbranched)

This is a **completely different derivation route** from S4's
`Epsilon2 @ kp - coupling` construction (`eigenmodes.py`'s
`solve_layer_eigenmodes_patterned`): RCWA.jl eliminates directly from
Maxwell's curl equations into a single matrix `M`, in a different
field-basis/gauge (its `W` is explicitly "transform towards Electric
fields," vs. this project's `phi` which is tied to S4's specific H-field-
adjacent convention -- see `eigenmodes.py`'s and `fields.py`'s docstrings).
Eigenvectors/normalization are therefore not expected to match and are not
compared here. What *is* expected to match, if both formulas are correct,
is the set of eigenvalues (`q^2`, the physical propagation constants) --
these are basis-independent. This is analogous to how Phase 2 used two
independently-coded FFT-of-rasterized-mask implementations
(`Rigorous-Coupled-Wave-Analysis`'s `convmat2D.py`,
`RigorousCoupledWaveAnalysis.jl`'s `ft2d.jl::real2recip`) to cross-check
the *same* Toeplitz-matrix values from a different code path -- here the
cross-check is one level up, on the eigenoperator built from that matrix.

Honesty note (per `rules.md` AI Coding Rule 5): no independently-published
numeric ground truth was found in `RigorousCoupledWaveAnalysis.jl` for a
2D-patterned case -- its own `test/runtests.jl` (lines 70-111) uses
`rand()` parameters and only checks internal self-consistency (its own ETM
engine vs. its own SRCWA engine agree, R+T=1). So this remains a
cross-implementation eigenoperator check, not a literature benchmark --
the same honesty bar already applied to Phase 3's TM oracle caveat
(`tests/oracles/rcwa_1d_gaylord.py`). A true external, independently-
published 2D R/T oracle is still open work for Phase 4b (see
`tests/oracles/rcwa_2d_pillar.py`).
"""

from __future__ import annotations

import numpy as np


def eigenoperator_eigenvalues(omega: complex, kx: np.ndarray, ky: np.ndarray, epsilon_hat: np.ndarray) -> np.ndarray:
    """`q^2` eigenvalues of RCWA.jl's isotropic patterned-layer operator
    `M` (`Common.jl:93`), built from this project's own direct-rule
    Toeplitz matrix `epsilon_hat` (Phase 2's `toeplitz_matrix(...,
    inverse=False)`), for a given per-order `(kx, ky)` at frequency `omega`.

    Two conventions are reconciled here so the returned values are
    directly comparable to `LayerEigenmodes.q**2` from
    `solve_layer_eigenmodes_patterned`, both confirmed empirically (not
    just assumed) while building this oracle:

    1. **Normalization.** `RigorousCoupledWaveAnalysis.jl/src/Common/grids.jl:78-84`
       ("all k vectors are generally normalized to k0 here") builds `Kx`,
       `Ky` as `kx/k0`, `ky/k0` (dimensionless) -- `M` itself has no
       `omega`/`k0` term anywhere, unlike this project's `build_kp_matrix`,
       which is only consistent if its `Kx`/`Ky` inputs are already
       normalized. So `kx/omega`, `ky/omega` are fed into `M` here, and the
       raw eigenvalues are scaled by `omega**2` to undo the normalization.
    2. **Sign.** Even after normalization, `M`'s eigenvalues came out as
       exactly `-1` times this project's `q**2` (verified numerically to
       ~1e-12 across several `num_orders`/angle/pattern combinations before
       this function was finalized) -- an overall opposite time-convention
       choice (`exp(+iwt)` vs. `exp(-iwt)`, the same class of sign
       footnote already documented in `tests/oracles/rcwa_1d_gaylord.py`
       and `tests/oracles/empy_tmm.py` for the other vendored oracles), not
       a formula discrepancy. Negated here so the caller can compare
       directly against `modes.q**2` with no sign bookkeeping of their own.
    """
    kx = np.asarray(kx, dtype=complex) / omega
    ky = np.asarray(ky, dtype=complex) / omega
    e = np.asarray(epsilon_hat, dtype=complex)
    n = kx.shape[0]

    kx_diag = np.diag(kx)
    ky_diag = np.diag(ky)

    a = np.linalg.solve(e, kx_diag @ e)
    b = np.linalg.solve(e, ky_diag @ e)

    m = np.zeros((2 * n, 2 * n), dtype=complex)
    m[:n, :n] = np.diag(ky * ky) - e + kx_diag @ a
    m[:n, n:] = -ky_diag @ kx_diag + kx_diag @ b
    m[n:, :n] = -kx_diag @ ky_diag + ky_diag @ a
    m[n:, n:] = np.diag(kx * kx) - e + ky_diag @ b

    return -np.linalg.eigvals(m) * omega**2
