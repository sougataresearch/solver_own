"""Category 4 target 4.2 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): unit-cell
bounds policy. Two halves, per `geometry.validate_pattern_fits_lattice`'s
docstring:

1. A shape whose footprint crosses a conceptual cell edge is already
   handled correctly by the existing analytic Fourier-transform machinery
   with no code change -- verified here against a from-scratch raster
   reference that explicitly tiles the shape's periodic images (not just
   asserted from the Poisson-summation argument in prose).
2. A shape that could overlap its own periodic images is explicitly
   rejected at `Simulation` construction time, per the conservative
   two-primitive-vector policy documented in
   `validate_pattern_fits_lattice`.
"""

from __future__ import annotations

import numpy as np
import pytest

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.fourier_factorization import pattern_epsilon_hat
from sougata_solver.geometry import Circle, Lattice, Pattern, Polygon, Rectangle, validate_pattern_fits_lattice
from sougata_solver.staircase import staircase_rectangle_layers
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation

AIR = Material("air", 1.0)
CORE = Material("core", 4.0)


# ---------------------------------------------------------------------------
# 1. Edge-crossing shapes are automatically, correctly periodized
# ---------------------------------------------------------------------------

_LX, _LY = 1.0, 0.8
_N_GRID = 900


def _periodic_tiling_raster_coefficient(pattern: Pattern, lattice: Lattice, g1: int, g2: int) -> complex:
    """Independent reference: rasterize the fundamental cell `[0, Lx) x
    [0, Ly)` and, for each grid point, check containment against every
    shape's periodic images in a 3x3 block of neighboring cells (m, n in
    {-1,0,1}) -- i.e. explicitly build the periodic tiling by real-space
    translation and wraparound, not by relying on the analytic-Fourier-
    transform-at-reciprocal-lattice-points argument under test."""
    x = (np.arange(_N_GRID) + 0.5) / _N_GRID * _LX
    y = (np.arange(_N_GRID) + 0.5) / _N_GRID * _LY
    X, Y = np.meshgrid(x, y, indexing="ij")

    eps_grid = np.full(X.shape, complex(pattern.background.epsilon_tensor(1.0)[0, 0]))
    for shape in pattern.shapes:
        assert isinstance(shape, Circle), "reference only supports Circle for this test"
        mask = np.zeros(X.shape, dtype=bool)
        for m in (-1, 0, 1):
            for n in (-1, 0, 1):
                dx = X - (shape.center[0] + m * _LX)
                dy = Y - (shape.center[1] + n * _LY)
                mask |= dx * dx + dy * dy <= shape.radius**2
        eps_grid[mask] = complex(shape.material.epsilon_tensor(1.0)[0, 0])

    Lk = lattice.reciprocal_vectors()
    k = g1 * Lk[0] + g2 * Lk[1]
    phase = np.exp(-2j * np.pi * (k[0] * X + k[1] * Y))
    return complex(np.sum(eps_grid * phase) / (_N_GRID * _N_GRID))


def test_shape_crossing_cell_edge_matches_periodic_tiling_reference():
    """Circle centered at `(0.02*Lx, 0.4*Ly)`, radius `0.1*Lx` -- its
    footprint spans `x in [-0.08, 0.12]*Lx`, i.e. it pokes past `x=0`
    into negative territory relative to the `[0, Lx)` fundamental-cell
    convention used by the reference raster above, a genuine edge-crossing
    case."""
    lattice = Lattice(a=(_LX, 0.0), b=(0.0, _LY))
    pattern = Pattern(background=AIR)
    pattern.add(Circle(center=(0.02 * _LX, 0.4 * _LY), radius=0.1 * _LX, material=CORE))

    for g1, g2 in [(0, 0), (1, 0), (0, 1), (1, 1), (-1, 2)]:
        analytic = pattern_epsilon_hat(pattern, lattice, g1, g2, wavelength=1.0)
        reference = _periodic_tiling_raster_coefficient(pattern, lattice, g1, g2)
        assert analytic == pytest.approx(reference, abs=6e-3), (g1, g2, analytic, reference)


# ---------------------------------------------------------------------------
# 2. Self-overlap-across-periodic-boundary rejection
# ---------------------------------------------------------------------------


def test_validate_pattern_fits_lattice_accepts_ordinary_pattern():
    lattice = Lattice(a=(0.7, 0.0), b=(0.0, 0.7))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(0.35, 0.35), radius=0.18, material=CORE)])
    validate_pattern_fits_lattice(pattern, lattice)  # must not raise


def test_validate_pattern_fits_lattice_rejects_self_overlapping_circle():
    lattice = Lattice(a=(0.7, 0.0), b=(0.0, 0.7))
    # diameter = 0.8 > period = 0.7 -- guaranteed to overlap its own periodic image
    pattern = Pattern(background=AIR, shapes=[Circle(center=(0.35, 0.35), radius=0.4, material=CORE)])
    with pytest.raises(ValueError, match="periodic images"):
        validate_pattern_fits_lattice(pattern, lattice)


def test_validate_pattern_fits_lattice_rejects_self_overlapping_rectangle():
    lattice = Lattice(a=(0.7, 0.0), b=(0.0, 0.7))
    pattern = Pattern(background=AIR, shapes=[Rectangle(center=(0.35, 0.35), halfwidth=(0.5, 0.5), material=CORE)])
    with pytest.raises(ValueError, match="periodic images"):
        validate_pattern_fits_lattice(pattern, lattice)


def test_simulation_construction_rejects_self_overlapping_pattern():
    lattice = Lattice(a=(0.7, 0.0), b=(0.0, 0.7))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(0.35, 0.35), radius=0.4, material=CORE)])
    layer = Layer("bad", 0.3, pattern=pattern)
    with pytest.raises(ValueError, match="periodic images"):
        Simulation(lattice, [layer], num_orders=5, incidence=AIR, transmission=AIR)


def test_simulation_construction_accepts_near_touching_but_non_overlapping_pattern():
    """`radius = 0.49 * period` (the same near-touching-pillar case Phase 4b's
    stress sweep already exercises, `tests/test_2d_pillar_stress.py`) has
    `2*radius = 0.98*period < period` -- must not be rejected."""
    period = 0.7
    lattice = Lattice(a=(period, 0.0), b=(0.0, period))
    pattern = Pattern(
        background=AIR, shapes=[Circle(center=(period / 2, period / 2), radius=0.49 * period, material=CORE)]
    )
    layer = Layer("near_touching", 0.3, pattern=pattern)
    sim = Simulation(lattice, [layer], num_orders=5, incidence=AIR, transmission=AIR)
    excitation = PlaneWaveExcitation(1.0, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    result = sim.solve(excitation)
    assert result.reflectance() + result.transmittance() == pytest.approx(1.0, abs=1e-6)


# ---------------------------------------------------------------------------
# 3. `Pattern.skip_bounds_check` (`decisions.md` ADR-035): a narrow, explicit
#    escape hatch for a shape the conservative bounding-radius test flags as
#    a false positive -- an elongated polygon whose y-extent exactly spans
#    the lattice period (never crosses it) but whose corner-to-center
#    distance exceeds the period, purely from its nonzero x-width at that
#    same y-extreme.
# ---------------------------------------------------------------------------

# Half-widths and y-extent taken from `structures/trench/y_tapered_polygon_trench.py`
# -- the actual Lumerical structure the project owner shared, confirmed (via
# a direct check of the RCWA region's General tab: propagation axis = z,
# backward) to be z-uniform with the taper genuinely in-plane along y. An
# earlier verbal claim that propagation axis = y led to a different,
# discarded reading of this same shape (`decisions.md` ADR-035's
# "Second correction") -- this is the confirmed-correct case.
_TRAP_PERIOD = 2.028
_TRAP_VERTICES = (
    (-0.2295, -_TRAP_PERIOD / 2),
    (0.2295, -_TRAP_PERIOD / 2),
    (0.2625, _TRAP_PERIOD / 2),
    (-0.2625, _TRAP_PERIOD / 2),
)


def _polygon_mask(vertices, cx_offset, cy_offset, x, y):
    """Point-in-polygon test (matplotlib-free ray casting) for a polygon
    translated by `(cx_offset, cy_offset)`, evaluated on grid `x, y`."""
    verts = [(vx + cx_offset, vy + cy_offset) for vx, vy in vertices]
    n = len(verts)
    inside = np.zeros(x.shape, dtype=bool)
    for i in range(n):
        x1, y1 = verts[i]
        x2, y2 = verts[(i + 1) % n]
        crosses = (y1 > y) != (y2 > y)
        with np.errstate(invalid="ignore", divide="ignore"):
            x_intersect = (x2 - x1) * (y - y1) / (y2 - y1 + 1e-300) + x1
        inside ^= crosses & (x < x_intersect)
    return inside


def test_y_tapered_trapezoid_bounding_radius_flags_false_positive():
    """Confirms this exact shape trips the conservative check (the reason
    `skip_bounds_check` exists) but does NOT actually self-overlap -- verified
    directly by rasterizing the shape and its +/-period periodic images and
    checking for pixel overlap, independent of the analytic bounding-radius
    argument under test."""
    lattice = Lattice(a=(_TRAP_PERIOD, 0.0), b=(0.0, _TRAP_PERIOD))
    pattern = Pattern(background=AIR, shapes=[Polygon(center=(0.0, 0.0), vertices=_TRAP_VERTICES, material=CORE)])
    with pytest.raises(ValueError, match="periodic images"):
        validate_pattern_fits_lattice(pattern, lattice)

    n_grid = 600
    x = (np.arange(n_grid) + 0.5) / n_grid * _TRAP_PERIOD - _TRAP_PERIOD / 2
    y = (np.arange(n_grid) + 0.5) / n_grid * (3 * _TRAP_PERIOD) - 1.5 * _TRAP_PERIOD
    X, Y = np.meshgrid(x, y, indexing="ij")

    base = _polygon_mask(_TRAP_VERTICES, 0.0, 0.0, X, Y)
    image_up = _polygon_mask(_TRAP_VERTICES, 0.0, _TRAP_PERIOD, X, Y)
    image_down = _polygon_mask(_TRAP_VERTICES, 0.0, -_TRAP_PERIOD, X, Y)
    assert not np.any(base & image_up), "shape overlaps its +period periodic image"
    assert not np.any(base & image_down), "shape overlaps its -period periodic image"


def test_skip_bounds_check_lets_verified_safe_pattern_construct_and_solve():
    """With `skip_bounds_check=True`, `Simulation` construction and solving
    both proceed for the exact shape the previous test independently
    confirmed does not truly self-overlap."""
    lattice = Lattice(a=(_TRAP_PERIOD, 0.0), b=(0.0, _TRAP_PERIOD))
    pattern = Pattern(
        background=CORE,
        shapes=[Polygon(center=(0.0, 0.0), vertices=_TRAP_VERTICES, material=AIR)],
        skip_bounds_check=True,
    )
    layer = Layer("y_tapered_trench", 0.5, pattern=pattern)
    sim = Simulation(lattice, [layer], num_orders=5, incidence=AIR, transmission=AIR)
    excitation = PlaneWaveExcitation(1.0, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    result = sim.solve(excitation)
    assert result.reflectance() + result.transmittance() == pytest.approx(1.0, abs=1e-6)


def test_skip_bounds_check_defaults_to_false_and_does_not_weaken_other_patterns():
    """`skip_bounds_check` is per-`Pattern`, opt-in -- a genuinely
    self-overlapping shape with the default `False` is still rejected."""
    lattice = Lattice(a=(0.7, 0.0), b=(0.0, 0.7))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(0.35, 0.35), radius=0.4, material=CORE)])
    assert pattern.skip_bounds_check is False
    with pytest.raises(ValueError, match="periodic images"):
        Simulation(lattice, [Layer("bad", 0.3, pattern=pattern)], num_orders=5, incidence=AIR, transmission=AIR)


# ---------------------------------------------------------------------------
# 4. `skip_bounds_check` on `staircase_rectangle_layers` output -- a synthetic
#    case (no longer tied to a live `structures/` script; it was built for a
#    y-is-depth reading of the Lumerical file that a direct General-tab
#    check later disproved -- `decisions.md` ADR-035's "Second correction")
#    kept as a distinct regression case: a depth-tapered rectangle whose
#    y-halfwidth (2.286 um, fixed) is much larger than the x-period
#    (2.028 um), on a rectangular (non-square) lattice -- the circular
#    bounding-radius check flags this even though the shape is comfortably
#    clear of overlap on BOTH axes independently.
# ---------------------------------------------------------------------------


def test_elongated_rectangle_on_narrow_axis_flags_false_positive_then_verified_safe():
    period_x, period_y = 2.028, 6.572
    half_hy = 4.572 / 2  # 2.286 -- exceeds period_x (2.028) alone
    lattice = Lattice(a=(period_x, 0.0), b=(0.0, period_y))
    pattern = Pattern(
        background=AIR, shapes=[Rectangle(center=(0.0, -2.286), halfwidth=(0.2625, half_hy), material=CORE)]
    )
    with pytest.raises(ValueError, match="periodic images"):
        validate_pattern_fits_lattice(pattern, lattice)

    # Independent per-axis check: no overlap along x (2*0.2625 < period_x)
    # or y (2*half_hy < period_y) -- confirms the rejection above is a false
    # positive from the circular bound, not a genuine overlap.
    assert 2 * 0.2625 < period_x
    assert 2 * half_hy < period_y


def test_skip_bounds_check_lets_offset_tapered_trench_construct_and_solve():
    """Synthetic (see module note above): a depth-tapered `Rectangle` via
    `staircase_rectangle_layers`, offset on a 2D lattice, minus dispersive-
    material file I/O -- exercises the same `skip_bounds_check` mechanism
    for an axis-aligned elongated shape, end to end."""
    period_x, period_y = 2.028e-6, 6.572e-6
    lattice = Lattice(a=(period_x, 0.0), b=(0.0, period_y))
    layers = staircase_rectangle_layers(
        center=(0.0, -2.286e-6),
        top_halfwidth=(0.5 * 0.525e-6, 4.572e-6 / 2),
        bottom_halfwidth=(0.5 * 0.459e-6, 4.572e-6 / 2),
        thickness=2.028e-6,
        num_slices=4,
        shape_material=AIR,
        background_material=CORE,
        name_prefix="trench",
    )
    for layer in layers:
        if layer.pattern is not None:
            layer.pattern.skip_bounds_check = True

    sim = Simulation(lattice, layers, num_orders=5, incidence=AIR, transmission=AIR)
    excitation = PlaneWaveExcitation(1.0, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    result = sim.solve(excitation)
    assert result.reflectance() + result.transmittance() == pytest.approx(1.0, abs=1e-6)
