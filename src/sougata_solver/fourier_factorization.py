"""Toeplitz permittivity matrices for patterned layers.

Formulas transcribed directly from S4 (not re-derived), for two pieces:

1. The per-G-vector Fourier coefficient of a `Pattern`
   (`pattern_epsilon_hat`): `S4/S4/pattern/pattern.c::pattern_get_fourier_transform`
   (lines 889-1029). The DC (`G=0`) term is the background material's value;
   every shape thereafter contributes
   `(value_shape - value_parent) * shape.fourier_transform(k) / unit_cell_area`,
   where `value_parent` is the background's value if the shape has no parent
   (S4's subtraction rule for nested/composite shapes, mirrored here by
   `Pattern.containment_tree`). This matches the sign/phase convention
   already established in `geometry.py` (see its module docstring), so no
   extra sign-flip is introduced here.

2. Which Fourier coefficient goes into which matrix entry
   (`toeplitz_matrix`): `S4/S4/fmm/fmm_closed.cpp::FMMGetEpsilon_ClosedForm`
   (lines 77-106 for the direct matrix, lines 111-127 for the inverse-rule
   matrix before any further Li-factorization processing) — entry `(i, j)`
   is `hat{eps}(G_i - G_j)` (or `hat{1/eps}(G_i - G_j)` for the inverse
   rule), evaluated with `pattern_epsilon_hat`.

Cross-verified independently against two other implementations of the same
Toeplitz-matrix structure, both vendored under `REFERENCE/`:
`Rigorous-Coupled-Wave-Analysis/convolution_matrices/convmat2D.py` (Python)
and `RigorousCoupledWaveAnalysis.jl/src/Common/ft2d.jl::real2recip` (Julia).
Both independently build `matrix[i, j] = hat{eps}(G_i - G_j)` via
raster+FFT rather than S4's/this module's analytic shape transforms
(the raster+FFT alternative was considered and explicitly rejected for
this project, see `decisions.md`) — but the underlying Toeplitz structure
matches, which is exactly what `tests/test_fourier_factorization.py`'s
`test_toeplitz_matrix_matches_rcwa_fft_convmat_reference_*` tests check
(a from-scratch FFT-of-rasterized-pattern reference, not calling into this
module), per `rules.md`'s explicit sanctioning of an "S4/EMpy/RCWA.jl
cross-check" as an oracle.

Phase 2 scope only: scalar isotropic materials. `pattern_epsilon_hat` raises
`NotImplementedError` for anisotropic materials in a pattern; full tensor
factorization is Phase 6, per `phases.md`.
"""

from __future__ import annotations

import numpy as np

from sougata_solver.geometry import Lattice, Pattern
from sougata_solver.materials import Material


def _scalar_value(material: Material, wavelength: float, inverse: bool) -> complex:
    if not material.is_isotropic:
        raise NotImplementedError(
            f"material {material.name!r} is anisotropic; scalar Fourier "
            "factorization is Phase 2 scope only (see phases.md Phase 6 "
            "for tensor permittivity factorization)"
        )
    eps = complex(material.epsilon_tensor(wavelength)[0, 0])
    return 1.0 / eps if inverse else eps


def pattern_epsilon_hat(
    pattern: Pattern,
    lattice: Lattice,
    g1: int,
    g2: int,
    wavelength: float,
    inverse: bool = False,
) -> complex:
    """Fourier coefficient `hat{eps}(G)` (or `hat{1/eps}(G)` if
    `inverse=True`) of the scalar permittivity pattern at reciprocal
    lattice index `(g1, g2)`. See module docstring for the source formula.
    """
    Lk = lattice.reciprocal_vectors()
    k = g1 * Lk[0] + g2 * Lk[1]
    kx, ky = float(k[0]), float(k[1])
    is_dc = g1 == 0 and g2 == 0

    total = _scalar_value(pattern.background, wavelength, inverse) if is_dc else complex(0.0)

    area = lattice.unit_cell_area()
    parents = pattern.containment_tree()
    for i, shape in enumerate(pattern.shapes):
        parent_material = pattern.background if parents[i] is None else pattern.shapes[parents[i]].material
        dval = _scalar_value(shape.material, wavelength, inverse) - _scalar_value(
            parent_material, wavelength, inverse
        )
        total += dval * complex(shape.fourier_transform(kx, ky)) / area

    return total


def toeplitz_matrix(
    pattern: Pattern,
    lattice: Lattice,
    g_vectors: np.ndarray,
    wavelength: float,
    inverse: bool = False,
) -> np.ndarray:
    """Build the `(n, n)` Toeplitz permittivity matrix with entry
    `[i, j] = hat{eps}(G_i - G_j)` (or the inverse-rule variant if
    `inverse=True`), for the truncated G-vector set `g_vectors`
    (`(n, 2)` int array, e.g. from `fourier_basis.truncate_fourier_orders`).
    """
    n = len(g_vectors)
    matrix = np.zeros((n, n), dtype=complex)
    cache: dict[tuple[int, int], complex] = {}
    for i in range(n):
        for j in range(n):
            dg = (int(g_vectors[i, 0] - g_vectors[j, 0]), int(g_vectors[i, 1] - g_vectors[j, 1]))
            if dg not in cache:
                cache[dg] = pattern_epsilon_hat(pattern, lattice, dg[0], dg[1], wavelength, inverse=inverse)
            matrix[i, j] = cache[dg]
    return matrix
