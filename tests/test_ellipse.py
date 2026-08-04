"""Category 4 target 4.3 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): `Ellipse`
primitive. Tiers per `rules.md` Testing Requirements: DC/area unit checks,
a nonzero-Fourier-coefficient cross-check against a from-scratch rasterized
reference (same style as `test_fourier_factorization.py`'s `Circle`/
`Rectangle` checks), reduces-to-`Circle` regression (both the Fourier
transform and the independently-derived `contains`/`signed_distance_normal`
methods), and one end-to-end RCWA example
(`structures/via/elliptical_pillar.py`, per the Category 4 exit criteria's
"one end-to-end RCWA example" requirement).
"""

from __future__ import annotations

import numpy as np
import pytest

from sougata_solver.geometry import Circle, Ellipse, Lattice, Pattern
from sougata_solver.materials import Material

CORE = Material("core", 4.0)
BG = Material("bg", 1.0)


# ---------------------------------------------------------------------------
# DC / area
# ---------------------------------------------------------------------------


def test_ellipse_area():
    e = Ellipse(center=(0.0, 0.0), halfwidth=(0.3, 0.15), material=CORE)
    assert e.area == pytest.approx(np.pi * 0.3 * 0.15)


def test_ellipse_dc_value_equals_area():
    e = Ellipse(center=(0.1, -0.2), halfwidth=(0.3, 0.15), material=CORE)
    assert e.fourier_transform(0.0, 0.0) == pytest.approx(e.area)


# ---------------------------------------------------------------------------
# Nonzero Fourier coefficient vs. a from-scratch rasterized reference
# ---------------------------------------------------------------------------

_LX, _LY = 1.3, 0.9
_N_GRID = 900


def _rasterized_ellipse_fourier_transform(e: Ellipse, kx: float, ky: float) -> complex:
    """Independent Riemann-sum evaluation of the same continuous Fourier
    integral (does not call `Ellipse.fourier_transform`), same style as
    `test_fourier_factorization.py`'s `_rasterized_coefficient`."""
    x = (np.arange(_N_GRID) / _N_GRID - 0.5) * _LX
    y = (np.arange(_N_GRID) / _N_GRID - 0.5) * _LY
    X, Y = np.meshgrid(x, y, indexing="ij")
    dA = (_LX / _N_GRID) * (_LY / _N_GRID)

    ca, sa = np.cos(e.angle), np.sin(e.angle)
    lx = ca * (X - e.center[0]) + sa * (Y - e.center[1])
    ly = -sa * (X - e.center[0]) + ca * (Y - e.center[1])
    hx, hy = e.halfwidth
    mask = (lx / hx) ** 2 + (ly / hy) ** 2 <= 1.0

    phase = np.exp(-2j * np.pi * (kx * X + ky * Y))
    return complex(np.sum(np.where(mask, phase, 0.0)) * dA)


@pytest.mark.parametrize("angle", [0.0, 0.4, -0.9])
@pytest.mark.parametrize("kx,ky", [(1.0, 0.0), (0.0, 1.3), (1.5, -0.7)])
def test_ellipse_fourier_transform_matches_rasterized_reference(kx, ky, angle):
    e = Ellipse(center=(0.05, -0.03), halfwidth=(0.3, 0.15), material=CORE, angle=angle)
    analytic = e.fourier_transform(kx, ky)
    reference = _rasterized_ellipse_fourier_transform(e, kx, ky)
    assert analytic == pytest.approx(reference, abs=6e-3)


# ---------------------------------------------------------------------------
# Reduction to Circle when hx == hy
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("kx,ky", [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (0.7, -1.1)])
def test_ellipse_reduces_to_circle_fourier_transform(kx, ky):
    ellipse = Ellipse(center=(0.1, 0.2), halfwidth=(0.25, 0.25), material=CORE)
    circle = Circle(center=(0.1, 0.2), radius=0.25, material=CORE)
    assert ellipse.fourier_transform(kx, ky) == pytest.approx(circle.fourier_transform(kx, ky), abs=1e-10)


def test_ellipse_reduces_to_circle_contains():
    ellipse = Ellipse(center=(0.1, 0.2), halfwidth=(0.25, 0.25), material=CORE)
    circle = Circle(center=(0.1, 0.2), radius=0.25, material=CORE)
    for x, y in [(0.1, 0.2), (0.3, 0.2), (0.34, 0.2), (0.36, 0.2), (0.5, 0.5)]:
        assert ellipse.contains(x, y) == circle.contains(x, y)


def test_ellipse_reduces_to_circle_normal():
    ellipse = Ellipse(center=(0.0, 0.0), halfwidth=(0.25, 0.25), material=CORE)
    circle = Circle(center=(0.0, 0.0), radius=0.25, material=CORE)
    for x, y in [(0.25, 0.0), (0.0, 0.25), (0.18, 0.18)]:
        assert ellipse.signed_distance_normal(x, y) == pytest.approx(circle.signed_distance_normal(x, y), abs=1e-10)


# ---------------------------------------------------------------------------
# Rotation sanity
# ---------------------------------------------------------------------------


def test_rotated_ellipse_contains_matches_unrotated_after_inverse_rotation():
    e_flat = Ellipse(center=(0.0, 0.0), halfwidth=(0.3, 0.1), material=CORE, angle=0.0)
    e_rot = Ellipse(center=(0.0, 0.0), halfwidth=(0.3, 0.1), material=CORE, angle=np.pi / 2)
    # A 90-degree rotation swaps the roles of the semi-axes for containment purposes.
    assert e_flat.contains(0.25, 0.0)
    assert e_rot.contains(0.0, 0.25)
    assert not e_rot.contains(0.25, 0.0)


# ---------------------------------------------------------------------------
# Whole-Pattern integration: Toeplitz entry sanity (nonzero, finite)
# ---------------------------------------------------------------------------


def test_ellipse_pattern_toeplitz_entries_are_finite():
    from sougata_solver.fourier_basis import truncate_fourier_orders
    from sougata_solver.fourier_factorization import toeplitz_matrix

    lattice = Lattice(a=(0.7, 0.0), b=(0.0, 0.7))
    pattern = Pattern(background=BG, shapes=[Ellipse(center=(0.35, 0.35), halfwidth=(0.2, 0.1), material=CORE)])
    g = truncate_fourier_orders(lattice, num_orders=9)
    matrix = toeplitz_matrix(pattern, lattice, g, wavelength=1.0)
    assert np.all(np.isfinite(matrix))
