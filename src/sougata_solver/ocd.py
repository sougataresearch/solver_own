"""Semiconductor OCD (optical critical dimension) scatterometry helpers.

`COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 11. This module adds no new
physics: `OCDTrapezoidParams`/`trapezoid_trench_layers` are a validated,
CD-first (linewidth, not halfwidth) input parametrization built as a thin
wrapper around Phase 5's already-oracle-adjacent `staircase.staircase_slab_layers`
(the same staircase discretization Phase 5 already validated via a
convergence-vs-`num_slices` study); `rounded_rectangle_polygon` builds an
ordinary `geometry.Polygon` (already analytically Fourier-transformed and
validated, Category 4 targets 4.4/4.5) from arc-sampled vertices -- no new
Fourier-factorization or eigenmode formula anywhere in this module.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from sougata_solver.geometry import Polygon
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.staircase import staircase_slab_layers


@dataclass
class OCDTrapezoidParams:
    """Category 11 target 11.1: validated CD/period/height/material
    parameters for a trapezoidal trench cross-section -- the standard
    scatterometry parametrization (critical dimension = linewidth, not
    the `halfwidth` `staircase_slab_layers` itself takes).

    `sidewall_angle_deg` is deliberately a computed property, not a
    separately stored field: for a linear taper, the sidewall angle is
    fully determined by `top_cd`/`bottom_cd`/`height` (three numbers
    already fix the trapezoid), so storing a fourth, independently
    settable "sidewall angle" field would let it silently drift out of
    sync with the other three -- a validation liability with no
    corresponding benefit, per `rules.md`'s "don't add features not
    needed" guidance.
    """

    top_cd: float
    bottom_cd: float
    period: float
    height: float
    shape_material: Material
    background_material: Material

    def __post_init__(self) -> None:
        for label, value in (
            ("top_cd", self.top_cd),
            ("bottom_cd", self.bottom_cd),
            ("period", self.period),
            ("height", self.height),
        ):
            if not math.isfinite(value):
                raise ValueError(f"OCDTrapezoidParams.{label} must be finite, got {value!r}")
            if not (value > 0):
                raise ValueError(f"OCDTrapezoidParams.{label} must be > 0, got {value!r}")
        if self.top_cd > self.period:
            raise ValueError(f"top_cd ({self.top_cd!r}) must be <= period ({self.period!r})")
        if self.bottom_cd > self.period:
            raise ValueError(f"bottom_cd ({self.bottom_cd!r}) must be <= period ({self.period!r})")

    @property
    def sidewall_angle_deg(self) -> float:
        """Angle of the sidewall from vertical (`0` = untapered/vertical),
        derived from `top_cd`/`bottom_cd`/`height` -- see class docstring
        for why this is a property, not a stored field."""
        half_cd_diff = abs(self.top_cd - self.bottom_cd) / 2.0
        if half_cd_diff == 0.0:
            return 0.0
        return math.degrees(math.atan(half_cd_diff / self.height))


def trapezoid_trench_layers(
    params: OCDTrapezoidParams, num_slices: int, name_prefix: str = "trench"
) -> list[Layer]:
    """Category 11 target 11.2: staircase trapezoidal trench from
    `OCDTrapezoidParams` -- a thin CD-to-halfwidth wrapper
    (`halfwidth = CD / 2`) around `staircase.staircase_slab_layers`
    (Phase 5, already validated by a convergence-vs-`num_slices` study,
    `tests/test_staircase.py`), not a new discretization scheme. Zero
    taper (`top_cd == bottom_cd`) reduces exactly to that module's own
    zero-taper regression case.
    """
    return staircase_slab_layers(
        center_x=0.0,
        top_halfwidth=params.top_cd / 2.0,
        bottom_halfwidth=params.bottom_cd / 2.0,
        thickness=params.height,
        num_slices=num_slices,
        shape_material=params.shape_material,
        background_material=params.background_material,
        name_prefix=name_prefix,
    )


def rounded_rectangle_polygon(
    center: tuple[float, float],
    halfwidth: tuple[float, float],
    corner_radius: float,
    material: Material,
    num_arc_points: int = 8,
    angle: float = 0.0,
) -> Polygon:
    """Category 11 targets 11.3 (design)/11.4 (implementation): a
    rounded-rectangle cross-section (the standard OCD corner-rounding
    approximation for vias/pillars/pads), built as a `geometry.Polygon`
    with each of the four corners replaced by an arc of `num_arc_points`
    vertices sampling a quarter circle of radius `corner_radius` -- not a
    new geometry primitive or a new Fourier-transform formula: `Polygon`'s
    existing analytic Fourier transform (Category 4 target 4.5, already
    validated) handles the resulting vertex list unchanged.

    **Design (target 11.3)**: the periodic geometry approximation chosen
    is "polygon with arc-sampled corners" (not, e.g., a superellipse or a
    smoothed indicator function) because it reuses `Polygon` exactly as-is
    -- no new Shape subclass, no new Fourier-transform derivation.
    `num_arc_points` (`>= 2`, the two arc endpoints) is the convergence
    parameter: as it increases, the polygon's boundary approaches the true
    quarter-circle arc, and its area approaches the closed-form rounded-
    rectangle area `4*hx*hy - (4 - pi)*r^2` (a full rectangle minus the
    four corner squares' excess over the quarter circles removed from
    them) -- confirmed directly, not assumed
    (`tests/test_ocd.py::test_rounded_rectangle_area_converges_to_closed_form`).
    `corner_radius=0` reduces exactly to `Rectangle`'s four corners
    regardless of `num_arc_points` (every "arc" collapses to the single
    corner point).
    """
    hx, hy = halfwidth
    if not (0.0 <= corner_radius <= min(hx, hy)):
        raise ValueError(f"corner_radius must be in [0, min(halfwidth)={min(hx, hy)!r}], got {corner_radius!r}")
    if num_arc_points < 2:
        raise ValueError(f"num_arc_points must be >= 2, got {num_arc_points}")

    r = corner_radius
    corners = [
        (hx - r, hy - r, 0.0),      # top-right corner, arc from 0 to 90 deg
        (-(hx - r), hy - r, 90.0),  # top-left, 90 to 180 deg
        (-(hx - r), -(hy - r), 180.0),  # bottom-left, 180 to 270 deg
        (hx - r, -(hy - r), 270.0),  # bottom-right, 270 to 360 deg
    ]
    vertices: list[tuple[float, float]] = []
    for cx, cy, start_deg in corners:
        thetas = np.radians(start_deg + np.linspace(0.0, 90.0, num_arc_points))
        for theta in thetas:
            vertices.append((cx + r * math.cos(theta), cy + r * math.sin(theta)))

    return Polygon(center=center, vertices=tuple(vertices), material=material, angle=angle)
