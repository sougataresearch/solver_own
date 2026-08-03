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

Phase 2 scope: scalar isotropic materials (`pattern_epsilon_hat`/
`toeplitz_matrix`, unchanged since Phase 2). Category 1 target 1.6
(`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`) adds per-tensor-component Toeplitz
construction (`pattern_epsilon_hat_component`/`toeplitz_matrix_component`)
for patterned layers containing diagonal or in-plane-coupled anisotropic
materials, transcribed from `S4/S4/fmm/fmm_closed.cpp`'s `have_tensor`
branch (lines 165-256) -- see that function's docstring for the full
citation. Longitudinal coupling (eps_xz/eps_yz/eps_zx/eps_zy) remains
out of scope pending target 1.5.
"""

from __future__ import annotations

from typing import Callable

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


def _pattern_fourier_sum(
    pattern: Pattern,
    lattice: Lattice,
    g1: int,
    g2: int,
    value_fn: Callable[[Material], complex],
) -> complex:
    """Shared subtraction-rule accumulation used by both
    `pattern_epsilon_hat` (scalar `value_fn`) and
    `pattern_epsilon_hat_component` (tensor-component `value_fn`) -- the
    per-shape/parent bookkeeping (`Pattern.containment_tree`) is identical
    either way; only which scalar is extracted from a `Material` differs.
    See the module docstring for the source citation.
    """
    Lk = lattice.reciprocal_vectors()
    k = g1 * Lk[0] + g2 * Lk[1]
    kx, ky = float(k[0]), float(k[1])
    is_dc = g1 == 0 and g2 == 0

    total = value_fn(pattern.background) if is_dc else complex(0.0)

    area = lattice.unit_cell_area()
    parents = pattern.containment_tree()
    for i, shape in enumerate(pattern.shapes):
        parent_material = pattern.background if parents[i] is None else pattern.shapes[parents[i]].material
        dval = value_fn(shape.material) - value_fn(parent_material)
        total += dval * complex(shape.fourier_transform(kx, ky)) / area

    return total


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
    return _pattern_fourier_sum(pattern, lattice, g1, g2, lambda m: _scalar_value(m, wavelength, inverse))


def pattern_epsilon_hat_component(
    pattern: Pattern,
    lattice: Lattice,
    g1: int,
    g2: int,
    wavelength: float,
    row: int,
    col: int,
) -> complex:
    """Fourier coefficient `hat{eps_rc}(G)` of tensor component `(row, col)`
    of the (possibly anisotropic) permittivity pattern -- the direct-rule
    (Laurent's-rule) accumulation only, matching `fmm_closed.cpp`'s
    `have_tensor` branch (lines 214-256), which builds each of `Epsilon2`'s
    four in-plane quadrants (`eps_xx, eps_xy, eps_yx, eps_yy`) and the
    `eps_zz` matrix this same way -- no separate inverse-rule variant exists
    for tensor components (Li's inverse rule is 1D-scalar-only in this
    project, per `solve_layer_eigenmodes_1d`'s docstring; the anisotropic
    patterned solver instead numerically inverts the direct-rule `eps_zz`
    Toeplitz, see `eigenmodes.solve_layer_eigenmodes_patterned_inplane`).
    """
    return _pattern_fourier_sum(pattern, lattice, g1, g2, lambda m: complex(m.epsilon_tensor(wavelength)[row, col]))


def _toeplitz(n: int, g_vectors: np.ndarray, coefficient: Callable[[int, int], complex]) -> np.ndarray:
    matrix = np.zeros((n, n), dtype=complex)
    cache: dict[tuple[int, int], complex] = {}
    for i in range(n):
        for j in range(n):
            dg = (int(g_vectors[i, 0] - g_vectors[j, 0]), int(g_vectors[i, 1] - g_vectors[j, 1]))
            if dg not in cache:
                cache[dg] = coefficient(dg[0], dg[1])
            matrix[i, j] = cache[dg]
    return matrix


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
    return _toeplitz(
        n, g_vectors, lambda g1, g2: pattern_epsilon_hat(pattern, lattice, g1, g2, wavelength, inverse=inverse)
    )


def toeplitz_matrix_component(
    pattern: Pattern,
    lattice: Lattice,
    g_vectors: np.ndarray,
    wavelength: float,
    row: int,
    col: int,
) -> np.ndarray:
    """Build the `(n, n)` direct-rule Toeplitz matrix of permittivity
    tensor component `(row, col)`, entry `[i, j] = hat{eps_rc}(G_i - G_j)`
    -- Category 1 target 1.6's per-component generalization of
    `toeplitz_matrix`. See `pattern_epsilon_hat_component`'s docstring for
    the source citation.
    """
    n = len(g_vectors)
    return _toeplitz(
        n, g_vectors, lambda g1, g2: pattern_epsilon_hat_component(pattern, lattice, g1, g2, wavelength, row, col)
    )
