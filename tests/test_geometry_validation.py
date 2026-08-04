"""Category 4 target 4.1 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): validate at
construction (`design.md`'s Error Handling conventions, already followed by
`Layer.__post_init__`) that `Lattice`/`Lattice1D`/`Circle`/`Rectangle`/`Slab`
reject non-finite dimensions, degenerate lattice vectors, and invalid
(non-positive) shape sizes -- rather than letting a NaN/inf propagate into
a downstream Fourier-factorization sum or eigensolve as a silently wrong
number, or a degenerate lattice surface as a cryptic `LinAlgError` from
`Lattice.reciprocal_vectors()`'s `np.linalg.inv` call.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from sougata_solver.geometry import Circle, Lattice, Lattice1D, Rectangle, Slab
from sougata_solver.materials import Material

MAT = Material("m", 2.0)


# ---------------------------------------------------------------------------
# Lattice
# ---------------------------------------------------------------------------


def test_lattice_accepts_ordinary_basis():
    lattice = Lattice(a=(0.7, 0.0), b=(0.0, 0.7))
    assert lattice.unit_cell_area() == pytest.approx(0.49)


@pytest.mark.parametrize("a", [(math.nan, 0.0), (math.inf, 0.0), (0.7, math.nan)])
def test_lattice_rejects_non_finite_basis_vector(a):
    with pytest.raises(ValueError, match="finite"):
        Lattice(a=a, b=(0.0, 0.7))


def test_lattice_rejects_collinear_basis_vectors():
    with pytest.raises(ValueError, match="degenerate"):
        Lattice(a=(1.0, 0.0), b=(2.0, 0.0))


def test_lattice_rejects_zero_basis_vector():
    with pytest.raises(ValueError, match="degenerate"):
        Lattice(a=(0.0, 0.0), b=(0.0, 0.7))


# ---------------------------------------------------------------------------
# Lattice1D
# ---------------------------------------------------------------------------


def test_lattice1d_accepts_positive_period():
    assert Lattice1D(0.7).period == pytest.approx(0.7)


@pytest.mark.parametrize("period", [0.0, -0.5, math.nan, math.inf])
def test_lattice1d_rejects_invalid_period(period):
    with pytest.raises(ValueError):
        Lattice1D(period)


# ---------------------------------------------------------------------------
# Circle
# ---------------------------------------------------------------------------


def test_circle_accepts_ordinary_params():
    Circle(center=(0.1, -0.2), radius=0.3, material=MAT)


@pytest.mark.parametrize("radius", [0.0, -0.1, math.nan, math.inf])
def test_circle_rejects_invalid_radius(radius):
    with pytest.raises(ValueError):
        Circle(center=(0.0, 0.0), radius=radius, material=MAT)


def test_circle_rejects_non_finite_center():
    with pytest.raises(ValueError, match="finite"):
        Circle(center=(math.nan, 0.0), radius=0.2, material=MAT)


# ---------------------------------------------------------------------------
# Rectangle
# ---------------------------------------------------------------------------


def test_rectangle_accepts_ordinary_params():
    Rectangle(center=(0.0, 0.0), halfwidth=(0.2, 0.1), material=MAT, angle=0.3)


@pytest.mark.parametrize("halfwidth", [(0.0, 0.1), (0.1, -0.2), (math.nan, 0.1), (math.inf, 0.1)])
def test_rectangle_rejects_invalid_halfwidth(halfwidth):
    with pytest.raises(ValueError):
        Rectangle(center=(0.0, 0.0), halfwidth=halfwidth, material=MAT)


def test_rectangle_rejects_non_finite_angle():
    with pytest.raises(ValueError, match="finite"):
        Rectangle(center=(0.0, 0.0), halfwidth=(0.1, 0.1), material=MAT, angle=math.inf)


# ---------------------------------------------------------------------------
# Slab
# ---------------------------------------------------------------------------


def test_slab_accepts_ordinary_params():
    Slab(center_x=0.1, halfwidth=0.2, material=MAT)


@pytest.mark.parametrize("halfwidth", [0.0, -0.1, math.nan, math.inf])
def test_slab_rejects_invalid_halfwidth(halfwidth):
    with pytest.raises(ValueError):
        Slab(center_x=0.0, halfwidth=halfwidth, material=MAT)


def test_slab_rejects_non_finite_center_x():
    with pytest.raises(ValueError, match="finite"):
        Slab(center_x=math.nan, halfwidth=0.2, material=MAT)
