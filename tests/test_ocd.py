"""Category 11 targets 11.1-11.4 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
`ocd.OCDTrapezoidParams`, `ocd.trapezoid_trench_layers`, and
`ocd.rounded_rectangle_polygon`. No new physics anywhere in this module
(see `ocd.py`'s module docstring) -- these tests check construction-time
validation, reduction to already-validated cases, and (for the corner-
rounding geometry) convergence to a closed-form analytic area, not a new
oracle-comparison (none is needed; nothing new is being solved).
"""

from __future__ import annotations

import math

import pytest

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice1D
from sougata_solver.materials import Material
from sougata_solver.ocd import OCDTrapezoidParams, rounded_rectangle_polygon, trapezoid_trench_layers
from sougata_solver.simulation import Simulation
from sougata_solver.staircase import staircase_slab_layers

AIR = Material("air", 1.0)
SI = Material("si", 3.48**2)


# ---------------------------------------------------------------------------
# 11.1 OCDTrapezoidParams validation
# ---------------------------------------------------------------------------


def test_ocd_trapezoid_params_valid_construction():
    params = OCDTrapezoidParams(top_cd=0.2, bottom_cd=0.15, period=0.7, height=0.3, shape_material=SI, background_material=AIR)
    assert params.top_cd == 0.2
    assert params.bottom_cd == 0.15


@pytest.mark.parametrize("field", ["top_cd", "bottom_cd", "period", "height"])
@pytest.mark.parametrize("bad_value", [0.0, -0.1, float("nan"), float("inf")])
def test_ocd_trapezoid_params_rejects_invalid_values(field, bad_value):
    kwargs = dict(top_cd=0.2, bottom_cd=0.15, period=0.7, height=0.3, shape_material=SI, background_material=AIR)
    kwargs[field] = bad_value
    with pytest.raises(ValueError):
        OCDTrapezoidParams(**kwargs)


def test_ocd_trapezoid_params_rejects_cd_wider_than_period():
    with pytest.raises(ValueError, match="period"):
        OCDTrapezoidParams(top_cd=0.8, bottom_cd=0.15, period=0.7, height=0.3, shape_material=SI, background_material=AIR)


def test_sidewall_angle_zero_for_untapered_trapezoid():
    params = OCDTrapezoidParams(top_cd=0.2, bottom_cd=0.2, period=0.7, height=0.3, shape_material=SI, background_material=AIR)
    assert params.sidewall_angle_deg == 0.0


def test_sidewall_angle_matches_hand_computed_geometry():
    top_cd, bottom_cd, height = 0.3, 0.1, 0.3
    params = OCDTrapezoidParams(top_cd=top_cd, bottom_cd=bottom_cd, period=0.7, height=height, shape_material=SI, background_material=AIR)
    expected = math.degrees(math.atan(((top_cd - bottom_cd) / 2) / height))
    assert params.sidewall_angle_deg == pytest.approx(expected)


# ---------------------------------------------------------------------------
# 11.2 trapezoid_trench_layers
# ---------------------------------------------------------------------------


def test_trapezoid_trench_layers_zero_taper_matches_staircase_slab_layers_directly():
    """Zero taper (`top_cd == bottom_cd`) must reduce exactly to calling
    `staircase_slab_layers` directly with `halfwidth = CD/2` -- confirming
    the CD-to-halfwidth wrapper introduces no discrepancy."""
    params = OCDTrapezoidParams(top_cd=0.3, bottom_cd=0.3, period=0.7, height=0.4, shape_material=SI, background_material=AIR)
    ocd_layers = trapezoid_trench_layers(params, num_slices=4)
    direct_layers = staircase_slab_layers(
        center_x=0.0, top_halfwidth=0.15, bottom_halfwidth=0.15, thickness=0.4,
        num_slices=4, shape_material=SI, background_material=AIR,
    )
    for a, b in zip(ocd_layers, direct_layers):
        assert a.thickness == pytest.approx(b.thickness)
        assert a.pattern.shapes[0].halfwidth == pytest.approx(b.pattern.shapes[0].halfwidth)


def test_trapezoid_trench_layers_end_to_end_energy_conservation():
    params = OCDTrapezoidParams(top_cd=0.3, bottom_cd=0.15, period=0.7, height=0.4, shape_material=SI, background_material=AIR)
    layers = trapezoid_trench_layers(params, num_slices=8)
    lattice = Lattice1D(params.period)
    sim = Simulation(lattice, layers, num_orders=9, incidence=AIR, transmission=AIR)
    excitation = PlaneWaveExcitation(1.0, math.radians(10.0), 0.0, s_amplitude=1.0, p_amplitude=0.0)
    result = sim.solve(excitation)
    assert result.reflectance() + result.transmittance() == pytest.approx(1.0, abs=1e-6)


# ---------------------------------------------------------------------------
# 11.3/11.4 rounded_rectangle_polygon
# ---------------------------------------------------------------------------


def test_rounded_rectangle_zero_radius_matches_plain_rectangle_area():
    poly = rounded_rectangle_polygon(center=(0.0, 0.0), halfwidth=(0.2, 0.15), corner_radius=0.0, material=SI, num_arc_points=8)
    assert poly.area == pytest.approx(4 * 0.2 * 0.15)


def test_rounded_rectangle_area_converges_to_closed_form():
    """Closed-form rounded-rectangle area: full rectangle minus the four
    corner squares' excess over the quarter circles removed from them,
    `4*hx*hy - (4-pi)*r^2`. As `num_arc_points` increases, the polygon's
    (exact, analytic-Fourier-transform-derived) area must converge to
    this value monotonically with shrinking error."""
    hx, hy, r = 0.2, 0.15, 0.05
    analytic_area = 4 * hx * hy - (4 - math.pi) * r**2
    orders = [2, 4, 8, 16, 32, 64]
    errors = []
    for n in orders:
        poly = rounded_rectangle_polygon(center=(0.0, 0.0), halfwidth=(hx, hy), corner_radius=r, material=SI, num_arc_points=n)
        errors.append(abs(poly.area - analytic_area))
    assert all(e1 >= e2 for e1, e2 in zip(errors, errors[1:]))
    assert errors[-1] < 1e-6


def test_rounded_rectangle_rejects_corner_radius_too_large():
    with pytest.raises(ValueError, match="corner_radius"):
        rounded_rectangle_polygon(center=(0.0, 0.0), halfwidth=(0.2, 0.15), corner_radius=0.2, material=SI, num_arc_points=8)


def test_rounded_rectangle_rejects_too_few_arc_points():
    with pytest.raises(ValueError, match="num_arc_points"):
        rounded_rectangle_polygon(center=(0.0, 0.0), halfwidth=(0.2, 0.15), corner_radius=0.05, material=SI, num_arc_points=1)


def test_rounded_rectangle_end_to_end_energy_conservation():
    from sougata_solver.geometry import Lattice, Pattern
    from sougata_solver.layer import Layer

    poly = rounded_rectangle_polygon(center=(0.35, 0.35), halfwidth=(0.15, 0.12), corner_radius=0.04, material=SI, num_arc_points=8)
    lattice = Lattice((0.7, 0.0), (0.0, 0.7))
    pattern = Pattern(background=AIR, shapes=[poly])
    sim = Simulation(lattice, [Layer("rounded_pad", 0.3, pattern=pattern)], num_orders=25, incidence=AIR, transmission=AIR)
    excitation = PlaneWaveExcitation(1.0, math.radians(10.0), 0.0, s_amplitude=1.0, p_amplitude=0.0)
    result = sim.solve(excitation)
    assert result.reflectance() + result.transmittance() == pytest.approx(1.0, abs=1e-6)
