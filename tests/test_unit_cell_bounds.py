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
from sougata_solver.geometry import Circle, Lattice, Pattern, Rectangle, validate_pattern_fits_lattice
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
