"""Category 4 targets 4.4/4.5 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
`Polygon` primitive. Tiers per `rules.md` Testing Requirements: reduction
to `Rectangle` (a square expressed as a `Polygon` must match exactly --
a strong regression check since `Rectangle` is already oracle-adjacent
tested), a from-scratch rasterized-reference cross-check for a genuinely
non-rectangular case (triangle, and a non-convex L-shape, per target 4.5's
"validate ... against numerical integration" wording), input validation
(too few vertices, CW winding), and one end-to-end RCWA example
(`structures/via/triangular_pillar.py`).
"""

from __future__ import annotations

import numpy as np
import pytest

from sougata_solver.geometry import Circle, Lattice, Pattern, Polygon, Rectangle
from sougata_solver.materials import Material

CORE = Material("core", 4.0)
BG = Material("bg", 1.0)


# ---------------------------------------------------------------------------
# Reduction to Rectangle (a square polygon)
# ---------------------------------------------------------------------------

_SQUARE_VERTICES = ((-0.2, -0.15), (0.2, -0.15), (0.2, 0.15), (-0.2, 0.15))


def test_square_polygon_area_matches_rectangle():
    poly = Polygon(center=(0.1, -0.2), vertices=_SQUARE_VERTICES, material=CORE)
    rect = Rectangle(center=(0.1, -0.2), halfwidth=(0.2, 0.15), material=CORE)
    assert poly.area == pytest.approx(rect.area)


@pytest.mark.parametrize("kx,ky", [(0.0, 0.0), (1.0, 0.0), (0.0, 1.3), (0.7, -1.1), (-2.2, 0.9)])
def test_square_polygon_fourier_transform_matches_rectangle(kx, ky):
    poly = Polygon(center=(0.1, -0.2), vertices=_SQUARE_VERTICES, material=CORE)
    rect = Rectangle(center=(0.1, -0.2), halfwidth=(0.2, 0.15), material=CORE)
    assert poly.fourier_transform(kx, ky) == pytest.approx(rect.fourier_transform(kx, ky), abs=1e-10)


def test_square_polygon_contains_matches_rectangle():
    poly = Polygon(center=(0.1, -0.2), vertices=_SQUARE_VERTICES, material=CORE)
    rect = Rectangle(center=(0.1, -0.2), halfwidth=(0.2, 0.15), material=CORE)
    for x, y in [(0.1, -0.2), (0.29, -0.2), (0.31, -0.2), (0.1, -0.05), (0.1, -0.04)]:
        assert poly.contains(x, y) == rect.contains(x, y)


def test_square_polygon_normal_matches_rectangle_away_from_corners():
    poly = Polygon(center=(0.1, -0.2), vertices=_SQUARE_VERTICES, material=CORE)
    rect = Rectangle(center=(0.1, -0.2), halfwidth=(0.2, 0.15), material=CORE)
    for x, y in [(0.3, -0.2), (0.1, -0.05), (-0.1, -0.2), (0.1, -0.35)]:
        assert poly.signed_distance_normal(x, y) == pytest.approx(rect.signed_distance_normal(x, y), abs=1e-10)


def test_rotated_square_polygon_matches_rotated_rectangle():
    poly = Polygon(center=(0.0, 0.0), vertices=_SQUARE_VERTICES, material=CORE, angle=0.4)
    rect = Rectangle(center=(0.0, 0.0), halfwidth=(0.2, 0.15), material=CORE, angle=0.4)
    assert poly.fourier_transform(1.3, -0.6) == pytest.approx(rect.fourier_transform(1.3, -0.6), abs=1e-10)


# ---------------------------------------------------------------------------
# From-scratch rasterized reference: triangle and non-convex L-shape
# ---------------------------------------------------------------------------

_LX, _LY = 1.3, 0.9
_N_GRID = 900


def _point_in_polygon_numpy(X: np.ndarray, Y: np.ndarray, vertices) -> np.ndarray:
    """Independent (vectorized) PNPoly re-implementation for the raster
    reference -- deliberately not calling `Polygon.contains`."""
    n = len(vertices)
    inside = np.zeros(X.shape, dtype=bool)
    vjx, vjy = vertices[-1]
    for i in range(n):
        vix, viy = vertices[i]
        cond = (viy > Y) != (vjy > Y)
        with np.errstate(divide="ignore", invalid="ignore"):
            x_intersect = (vjx - vix) * (Y - viy) / (vjy - viy) + vix
        cross = cond & (X < x_intersect)
        inside = np.where(cross, ~inside, inside)
        vjx, vjy = vix, viy
    return inside


def _rasterized_polygon_fourier_transform(poly: Polygon, kx: float, ky: float) -> complex:
    x = (np.arange(_N_GRID) / _N_GRID - 0.5) * _LX
    y = (np.arange(_N_GRID) / _N_GRID - 0.5) * _LY
    X, Y = np.meshgrid(x, y, indexing="ij")
    dA = (_LX / _N_GRID) * (_LY / _N_GRID)

    ca, sa = np.cos(poly.angle), np.sin(poly.angle)
    lx = ca * (X - poly.center[0]) + sa * (Y - poly.center[1])
    ly = -sa * (X - poly.center[0]) + ca * (Y - poly.center[1])
    mask = _point_in_polygon_numpy(lx, ly, poly.vertices)

    phase = np.exp(-2j * np.pi * (kx * X + ky * Y))
    return complex(np.sum(np.where(mask, phase, 0.0)) * dA)


_TRIANGLE_VERTICES = ((0.0, 0.2), (-0.18, -0.12), (0.18, -0.12))
# Non-convex L-shape (6 vertices, CCW).
_L_SHAPE_VERTICES = ((-0.2, -0.2), (0.2, -0.2), (0.2, 0.0), (0.0, 0.0), (0.0, 0.2), (-0.2, 0.2))


@pytest.mark.parametrize("vertices", [_TRIANGLE_VERTICES, _L_SHAPE_VERTICES])
def test_polygon_area_matches_rasterized_reference(vertices):
    poly = Polygon(center=(0.0, 0.0), vertices=vertices, material=CORE)
    x = (np.arange(_N_GRID) / _N_GRID - 0.5) * _LX
    y = (np.arange(_N_GRID) / _N_GRID - 0.5) * _LY
    X, Y = np.meshgrid(x, y, indexing="ij")
    dA = (_LX / _N_GRID) * (_LY / _N_GRID)
    mask = _point_in_polygon_numpy(X, Y, vertices)
    rasterized_area = float(np.sum(mask)) * dA
    assert poly.area == pytest.approx(rasterized_area, rel=5e-3)


@pytest.mark.parametrize("vertices", [_TRIANGLE_VERTICES, _L_SHAPE_VERTICES])
@pytest.mark.parametrize("kx,ky", [(1.0, 0.0), (0.0, 1.3), (1.5, -0.7), (2.0, 2.0)])
def test_polygon_fourier_transform_matches_rasterized_reference(kx, ky, vertices):
    poly = Polygon(center=(0.05, -0.03), vertices=vertices, material=CORE)
    analytic = poly.fourier_transform(kx, ky)
    reference = _rasterized_polygon_fourier_transform(poly, kx, ky)
    assert analytic == pytest.approx(reference, abs=8e-3)


def test_polygon_contains_matches_pnpoly_reference():
    poly = Polygon(center=(0.0, 0.0), vertices=_L_SHAPE_VERTICES, material=CORE)
    # Inside the "notch" cut out of the L-shape -- must be excluded.
    assert not poly.contains(0.1, 0.1)
    # Inside the solid part of the L.
    assert poly.contains(-0.1, 0.1)
    assert poly.contains(0.1, -0.1)
    assert not poly.contains(0.3, 0.3)


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def test_polygon_rejects_fewer_than_three_vertices():
    with pytest.raises(ValueError, match="at least 3"):
        Polygon(center=(0.0, 0.0), vertices=((0.0, 0.0), (1.0, 0.0)), material=CORE)


def test_polygon_rejects_clockwise_winding():
    cw_square = tuple(reversed(_SQUARE_VERTICES))
    with pytest.raises(ValueError, match="CCW"):
        Polygon(center=(0.0, 0.0), vertices=cw_square, material=CORE)


def test_polygon_rejects_non_finite_vertex():
    with pytest.raises(ValueError, match="finite"):
        Polygon(center=(0.0, 0.0), vertices=((0.0, 0.0), (1.0, float("nan")), (0.0, 1.0)), material=CORE)


# ---------------------------------------------------------------------------
# Whole-Pattern integration: Toeplitz entry sanity (finite)
# ---------------------------------------------------------------------------


def test_polygon_pattern_toeplitz_entries_are_finite():
    from sougata_solver.fourier_basis import truncate_fourier_orders
    from sougata_solver.fourier_factorization import toeplitz_matrix

    lattice = Lattice(a=(0.7, 0.0), b=(0.0, 0.7))
    pattern = Pattern(
        background=BG, shapes=[Polygon(center=(0.35, 0.35), vertices=_TRIANGLE_VERTICES, material=CORE)]
    )
    g = truncate_fourier_orders(lattice, num_orders=9)
    matrix = toeplitz_matrix(pattern, lattice, g, wavelength=1.0)
    assert np.all(np.isfinite(matrix))
