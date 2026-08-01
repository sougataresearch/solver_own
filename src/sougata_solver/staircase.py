"""Staircase (z-discretized) layer-stack generator for tapered sidewalls.

Phase 5. Approximates a linearly-tapered via/trench sidewall as
`num_slices` uniform-in-z layers, each a `Circle`/`Rectangle`/`Slab`
pattern whose size is linearly interpolated between a `top` and `bottom`
value -- the standard staircase/multi-slice approximation used throughout
RCWA literature for slanted-sidewall gratings, per `decisions.md` ADR-004
("Tapered sidewalls via staircase discretization, not new Fourier math").

This module is **not** transcribed from any vendored `REFERENCE/` repo: per
the `phase-reference-picker` skill's procedure, a grep across every
RCWA-family repo (`S4`, `EMpy`, `RigorousCoupledWaveAnalysis.jl`,
`Rigorous-Coupled-Wave-Analysis`) for "stair"/"taper" found no matching
staircase-generator code to cite (only unrelated hits in `meep`/`gprMax`
docs, a different numerical method per that skill's own guidance) -- so
this is independently derived, per `rules.md` AI Coding Rule 1, and flagged
for the extra test scrutiny that entails. Its correctness claim rests on
the convergence-vs-`num_slices` study in `tests/test_staircase.py`, not an
external oracle: each individual slice reuses Phase 3/4a's
already-oracle-validated per-layer eigenmode solve unchanged, so the only
new risk here is the discretization itself, not a new physics formula
(matches `phases.md` Phase 5's "no new Fourier/eigenmode math" scoping).

Interpolation convention: slice `i` (0-indexed, `i=0` nearest the `top`
size, `i=num_slices-1` nearest `bottom`) uses the linearly-interpolated
size at its z-midpoint, `frac = (i + 0.5) / num_slices` -- a symmetric
choice with no a-priori reason to bias toward either slice edge. All
`num_slices` slices have equal thickness `thickness / num_slices`.
"""

from __future__ import annotations

import numpy as np

from sougata_solver.geometry import Circle, Pattern, Rectangle, Slab
from sougata_solver.layer import Layer
from sougata_solver.materials import Material


def _slice_fractions(num_slices: int) -> np.ndarray:
    if num_slices < 1:
        raise ValueError("num_slices must be >= 1")
    return (np.arange(num_slices) + 0.5) / num_slices


def staircase_circle_layers(
    center: tuple[float, float],
    top_radius: float,
    bottom_radius: float,
    thickness: float,
    num_slices: int,
    shape_material: Material,
    background_material: Material,
    name_prefix: str = "via",
) -> list[Layer]:
    """`num_slices` uniform-in-z `Circle` layers approximating a linearly
    tapered via: radius `top_radius` at z=0 tapering to `bottom_radius` at
    z=thickness. See module docstring for the interpolation convention."""
    fracs = _slice_fractions(num_slices)
    radii = top_radius + (bottom_radius - top_radius) * fracs
    slice_thickness = thickness / num_slices
    layers = []
    for i, radius in enumerate(radii):
        pattern = Pattern(
            background=background_material,
            shapes=[Circle(center=center, radius=float(radius), material=shape_material)],
        )
        layers.append(Layer(f"{name_prefix}_{i}", slice_thickness, pattern=pattern))
    return layers


def staircase_rectangle_layers(
    center: tuple[float, float],
    top_halfwidth: tuple[float, float],
    bottom_halfwidth: tuple[float, float],
    thickness: float,
    num_slices: int,
    shape_material: Material,
    background_material: Material,
    angle: float = 0.0,
    name_prefix: str = "via",
) -> list[Layer]:
    """Same as `staircase_circle_layers`, for a `Rectangle` via whose
    `(hx, hy)` halfwidths taper independently between `top_halfwidth` and
    `bottom_halfwidth`."""
    fracs = _slice_fractions(num_slices)
    top_hw = np.asarray(top_halfwidth, dtype=float)
    bottom_hw = np.asarray(bottom_halfwidth, dtype=float)
    slice_thickness = thickness / num_slices
    layers = []
    for i, frac in enumerate(fracs):
        hw = top_hw + (bottom_hw - top_hw) * frac
        pattern = Pattern(
            background=background_material,
            shapes=[
                Rectangle(
                    center=center,
                    halfwidth=(float(hw[0]), float(hw[1])),
                    material=shape_material,
                    angle=angle,
                )
            ],
        )
        layers.append(Layer(f"{name_prefix}_{i}", slice_thickness, pattern=pattern))
    return layers


def staircase_slab_layers(
    center_x: float,
    top_halfwidth: float,
    bottom_halfwidth: float,
    thickness: float,
    num_slices: int,
    shape_material: Material,
    background_material: Material,
    name_prefix: str = "trench",
) -> list[Layer]:
    """Same as `staircase_circle_layers`, for a `Slab` (1D trench) whose
    halfwidth tapers between `top_halfwidth` and `bottom_halfwidth`, center
    fixed at `center_x`."""
    fracs = _slice_fractions(num_slices)
    slice_thickness = thickness / num_slices
    layers = []
    for i, frac in enumerate(fracs):
        halfwidth = top_halfwidth + (bottom_halfwidth - top_halfwidth) * frac
        pattern = Pattern(
            background=background_material,
            shapes=[Slab(center_x=center_x, halfwidth=float(halfwidth), material=shape_material)],
        )
        layers.append(Layer(f"{name_prefix}_{i}", slice_thickness, pattern=pattern))
    return layers
