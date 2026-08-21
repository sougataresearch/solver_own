"""Phase 4c tests: analytic Normal Vector Method Fourier factorization.

Independently derived (not transcribed -- see `decisions.md` ADR-012), so
per `rules.md`'s Testing Requirements this needs its own numerical
cross-check (rasterize-and-sum, same technique
`test_fourier_factorization.py` uses for `epsilon_hat`) before any of it is
trusted, plus the physical-invariant checks (energy conservation,
convergence-rate improvement over the plain-Laurent baseline)
`phases.md` Phase 4c's deliverables call for. No second independent
oracle exists for this specific formula in any vendored repo (unlike
Phase 4a's `RigorousCoupledWaveAnalysis.jl` eigenvalue cross-check) -- see
`references.md`'s Reference Survey.
"""

from __future__ import annotations

import numpy as np
import pytest

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.fourier_factorization import pattern_normal_outer_hat
from sougata_solver.geometry import Circle, Lattice, Pattern, Rectangle
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation

# ---------------------------------------------------------------------------
# Rasterize-and-sum cross-check, mirroring test_fourier_factorization.py's
# _rasterized_coefficient technique but for (nx*nx, nx*ny, ny*ny) instead of
# eps -- an independent, from-scratch numerical evaluation that never calls
# into geometry.py's/fourier_factorization.py's analytic formulas.
# ---------------------------------------------------------------------------

_LX, _LY = 1.3, 0.9
_N_GRID = 900


def _rasterized_normal_outer(
    pattern: Pattern, lattice: Lattice, g1: int, g2: int
) -> tuple[complex, complex, complex]:
    x = (np.arange(_N_GRID) / _N_GRID - 0.5) * _LX
    y = (np.arange(_N_GRID) / _N_GRID - 0.5) * _LY
    X, Y = np.meshgrid(x, y, indexing="ij")

    pxx_grid = np.zeros(X.shape)
    pxy_grid = np.zeros(X.shape)
    pyy_grid = np.zeros(X.shape)
    for shape in pattern.shapes:
        if isinstance(shape, Circle):
            dx, dy = X - shape.center[0], Y - shape.center[1]
            r = np.hypot(dx, dy)
            mask = r <= shape.radius
            r_safe = np.where(r < 1e-14, 1.0, r)
            nx, ny = dx / r_safe, dy / r_safe
        elif isinstance(shape, Rectangle):
            assert shape.angle == 0.0, "rasterized reference only supports axis-aligned rectangles"
            hx, hy = shape.halfwidth
            lx, ly = X - shape.center[0], Y - shape.center[1]
            mask = (np.abs(lx) <= hx) & (np.abs(ly) <= hy)
            close_to_x = np.abs(hx - np.abs(lx)) <= np.abs(hy - np.abs(ly))
            lx_safe = np.where(lx == 0.0, 1.0, lx)
            ly_safe = np.where(ly == 0.0, 1.0, ly)
            nx = np.where(close_to_x, np.sign(lx_safe), 0.0)
            ny = np.where(close_to_x, 0.0, np.sign(ly_safe))
        else:
            raise NotImplementedError(type(shape))
        pxx_grid[mask] = (nx * nx)[mask]
        pxy_grid[mask] = (nx * ny)[mask]
        pyy_grid[mask] = (ny * ny)[mask]

    Lk = lattice.reciprocal_vectors()
    k = g1 * Lk[0] + g2 * Lk[1]
    phase = np.exp(-2j * np.pi * (k[0] * X + k[1] * Y))
    norm = _N_GRID * _N_GRID
    pxx = complex(np.sum(pxx_grid * phase) / norm)
    pxy = complex(np.sum(pxy_grid * phase) / norm)
    pyy = complex(np.sum(pyy_grid * phase) / norm)
    return pxx, pxy, pyy


def test_circle_normal_outer_matches_rasterized_reference():
    lattice = Lattice(a=(_LX, 0.0), b=(0.0, _LY))
    core = Material("core", 4.0)
    pattern = Pattern(background=Material("bg", 1.0))
    pattern.add(Circle(center=(0.05, -0.03), radius=0.3, material=core))

    for g1, g2 in [(0, 0), (1, 0), (0, 1), (1, 1), (-1, 2)]:
        analytic = pattern_normal_outer_hat(pattern, lattice, g1, g2)
        reference = _rasterized_normal_outer(pattern, lattice, g1, g2)
        for a, r in zip(analytic, reference):
            assert np.isclose(a, r, atol=5e-3), (g1, g2, analytic, reference)


def test_rectangle_normal_outer_matches_rasterized_reference():
    lattice = Lattice(a=(_LX, 0.0), b=(0.0, _LY))
    core = Material("core", 9.0)
    pattern = Pattern(background=Material("bg", 2.25))
    pattern.add(Rectangle(center=(-0.02, 0.04), halfwidth=(0.25, 0.18), material=core))

    for g1, g2 in [(0, 0), (1, 0), (0, 1), (1, 1), (2, -1)]:
        analytic = pattern_normal_outer_hat(pattern, lattice, g1, g2)
        reference = _rasterized_normal_outer(pattern, lattice, g1, g2)
        for a, r in zip(analytic, reference):
            assert np.isclose(a, r, atol=8e-3), (g1, g2, analytic, reference)


def test_rectangle_wide_and_tall_cases_match_rasterized_reference():
    """`hx > hy` and `hx < hy` exercise the trapezoid vs. triangle
    degeneration of the same 4-vertex-list formula (see
    `Rectangle.normal_outer_fourier_transform`'s docstring) -- both
    branches need their own cross-check, not just a near-square case."""
    lattice = Lattice(a=(_LX, 0.0), b=(0.0, _LY))
    core = Material("core", 6.25)
    for halfwidth in [(0.35, 0.12), (0.1, 0.3)]:
        pattern = Pattern(background=Material("bg", 1.0))
        pattern.add(Rectangle(center=(0.0, 0.0), halfwidth=halfwidth, material=core))
        for g1, g2 in [(1, 0), (0, 1), (1, 1)]:
            analytic = pattern_normal_outer_hat(pattern, lattice, g1, g2)
            reference = _rasterized_normal_outer(pattern, lattice, g1, g2)
            for a, r in zip(analytic, reference):
                assert np.isclose(a, r, atol=8e-3), (halfwidth, g1, g2, analytic, reference)


# ---------------------------------------------------------------------------
# Self-consistency: Pxx+Pyy must equal the shape's already-validated plain
# indicator-function Fourier transform, since nx^2+ny^2=1 everywhere. This
# holds independently of the rasterize-and-sum reference above, so it's a
# second, algebraically-motivated check, not a restatement of the first.
# ---------------------------------------------------------------------------


def test_circle_pxx_plus_pyy_equals_plain_fourier_transform():
    circle = Circle(center=(0.1, -0.05), radius=0.25, material=Material("core", 4.0))
    for kx, ky in [(0.0, 0.0), (0.7, 0.0), (0.0, 1.1), (0.6, -0.9)]:
        pxx, _pxy, pyy = circle.normal_outer_fourier_transform(kx, ky)
        plain = circle.fourier_transform(kx, ky)
        assert np.isclose(pxx + pyy, plain, atol=1e-9), (kx, ky, pxx + pyy, plain)


def test_rectangle_pxx_plus_pyy_equals_plain_fourier_transform():
    rect = Rectangle(center=(0.0, 0.0), halfwidth=(0.2, 0.35), material=Material("core", 4.0))
    for kx, ky in [(0.0, 0.0), (0.7, 0.0), (0.0, 1.1), (0.6, -0.9)]:
        pxx, _pxy, pyy = rect.normal_outer_fourier_transform(kx, ky)
        plain = rect.fourier_transform(kx, ky)
        assert np.isclose(pxx + pyy, plain, atol=1e-9), (kx, ky, pxx + pyy, plain)


def test_rectangle_rotated_pxx_plus_pyy_equals_plain_fourier_transform():
    rect = Rectangle(
        center=(0.05, -0.1), halfwidth=(0.2, 0.35), material=Material("core", 4.0), angle=np.radians(27.0)
    )
    for kx, ky in [(0.7, 0.0), (0.0, 1.1), (0.6, -0.9)]:
        pxx, _pxy, pyy = rect.normal_outer_fourier_transform(kx, ky)
        plain = rect.fourier_transform(kx, ky)
        assert np.isclose(pxx + pyy, plain, atol=1e-9), (kx, ky, pxx + pyy, plain)


def test_circle_dc_limit_matches_direct_angular_integral():
    """At k=0, `integral cos^2(theta) r dr dtheta` over a disk is `area/2`
    by direct integration (`integral_0^{2pi} cos^2(theta) dtheta = pi`); this
    is a closed-form value independent of both the rasterized reference and
    the Bessel-function derivation, so it's a third, distinct check."""
    circle = Circle(center=(0.0, 0.0), radius=0.3, material=Material("core", 4.0))
    pxx, pxy, pyy = circle.normal_outer_fourier_transform(0.0, 0.0)
    assert np.isclose(pxx, circle.area / 2, atol=1e-9)
    assert np.isclose(pyy, circle.area / 2, atol=1e-9)
    assert np.isclose(pxy, 0.0, atol=1e-9)


def test_rectangle_dc_limit_pxy_vanishes():
    rect = Rectangle(center=(0.0, 0.0), halfwidth=(0.2, 0.35), material=Material("core", 4.0))
    pxx, pxy, pyy = rect.normal_outer_fourier_transform(0.0, 0.0)
    assert np.isclose(pxx + pyy, rect.area, atol=1e-9)
    assert np.isclose(pxy, 0.0, atol=1e-9)


def test_nested_shapes_raise_not_implemented():
    lattice = Lattice(a=(_LX, 0.0), b=(0.0, _LY))
    outer = Circle(center=(0.0, 0.0), radius=0.3, material=Material("outer", 4.0))
    inner = Circle(center=(0.0, 0.0), radius=0.1, material=Material("inner", 9.0))
    pattern = Pattern(background=Material("bg", 1.0), shapes=[outer, inner])
    with pytest.raises(NotImplementedError):
        pattern_normal_outer_hat(pattern, lattice, 1, 0)


# ---------------------------------------------------------------------------
# Simulation-level: reduction, energy conservation, convergence-rate
# ---------------------------------------------------------------------------

PERIOD = 0.7
AIR = Material("air", 1.0)
SI = Material("si", 3.48**2)
LAYER_THICKNESS = 0.46


def test_nvm_reduces_to_uniform_when_shape_matches_background():
    """`Delta = epsilon_hat - inv(epsilon_inv_hat)` must vanish exactly when
    there's no real material discontinuity, regardless of `Pxx`/`Pxy`/`Pyy`
    still being geometrically nonzero (the normal-vector field only depends
    on the shape's shape, not its material) -- this is a real test of the
    `Delta` mechanism nulling the correction, not just a trivial P=0 case."""
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.18, material=AIR)])
    layer = Layer("patch", LAYER_THICKNESS, pattern=pattern)
    sim = Simulation(lattice, [layer], num_orders=7, incidence=AIR, transmission=AIR, fourier_rule="nvm")
    excitation = PlaneWaveExcitation(wavelength=1.0, theta=0.0, phi=0.0, s_amplitude=0.7, p_amplitude=0.3)
    result = sim.solve(excitation)
    assert result.reflectance() == pytest.approx(0.0, abs=1e-8)
    assert result.transmittance() == pytest.approx(1.0, abs=1e-6)


def _pillar_simulation(num_orders: int, fourier_rule: str) -> Simulation:
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.18, material=SI)])
    layer = Layer("pillar", LAYER_THICKNESS, pattern=pattern)
    return Simulation(
        lattice, [layer], num_orders=num_orders, incidence=AIR, transmission=AIR, fourier_rule=fourier_rule
    )


@pytest.mark.parametrize("theta_deg", [0.0, 30.0])
@pytest.mark.parametrize("s_amp,p_amp", [(1.0, 0.0), (0.0, 1.0), (0.6, 0.8)])
def test_energy_conservation_nvm_pillar(theta_deg, s_amp, p_amp):
    sim = _pillar_simulation(num_orders=7, fourier_rule="nvm")
    excitation = PlaneWaveExcitation(
        wavelength=1.0, theta=np.radians(theta_deg), phi=0.0, s_amplitude=s_amp, p_amplitude=p_amp
    )
    result = sim.solve(excitation)
    de = result.diffraction_efficiencies()
    total = sum(de_r + de_t for de_r, de_t in de.values())
    assert total == pytest.approx(1.0, abs=1e-8)
    assert result.reflectance() + result.transmittance() == pytest.approx(1.0, abs=1e-8)


def test_nvm_and_laurent_diverge_for_high_contrast_pattern():
    """Sanity guard for the Simulation-level wiring itself (`fourier_rule`
    actually selects a different code path): a regression that silently
    made `"nvm"` behave identically to `"laurent"` (e.g. a copy-paste of
    the wrong kwarg) would make this assertion fail."""
    excitation = PlaneWaveExcitation(wavelength=1.0, theta=0.0, phi=0.0, s_amplitude=1.0, p_amplitude=0.0)
    r_laurent = _pillar_simulation(num_orders=7, fourier_rule="laurent").solve(excitation).reflectance()
    r_nvm = _pillar_simulation(num_orders=7, fourier_rule="nvm").solve(excitation).reflectance()
    assert abs(r_laurent - r_nvm) > 1e-4, (r_laurent, r_nvm)


@pytest.mark.slow
def test_nvm_converges_faster_than_laurent_for_high_contrast_rectangle():
    """Following the structure of S4's own worked examples in
    `S4/examples/2d/Li_JOSA_14_2758_1997/` (a high-contrast dielectric
    rectangle, `num_orders` swept) -- not reproducing that paper's exact
    numeric table (not available in this environment, see
    `references.md`'s Reference Survey), but testing the qualitative
    convergence-rate claim Phase 4c's whole point rests on: the Normal
    Vector Method's reflectance should approach a shared high-order
    reference value faster (smaller error at a given, moderate
    `num_orders`) than plain Laurent's rule, for the same high-index-contrast
    pattern. Uses NVM's own high-`num_orders` result as the reference value,
    since NVM is expected to be the more accurate rule -- both rules must
    converge to the same true continuum answer as `num_orders -> infinity`,
    so this is a legitimate common target, not circular (Laurent's rule is
    being compared against a target derived from the *other* rule, not
    itself)."""
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    core = Material("core", 2.25)  # matches Li (1997)'s example contrast
    pattern = Pattern(background=AIR, shapes=[Rectangle(center=(0.0, 0.0), halfwidth=(0.25, 0.25), material=core)])
    layer = Layer("grating", LAYER_THICKNESS, pattern=pattern)
    excitation = PlaneWaveExcitation(wavelength=1.0, theta=0.0, phi=0.0, s_amplitude=1.0, p_amplitude=0.0)

    def reflectance(num_orders: int, rule: str) -> float:
        sim = Simulation(lattice, [layer], num_orders=num_orders, incidence=AIR, transmission=AIR, fourier_rule=rule)
        return sim.solve(excitation).reflectance()

    ref = reflectance(121, "nvm")
    orders_to_check = [9, 25, 49]
    errors_nvm = [abs(reflectance(n, "nvm") - ref) for n in orders_to_check]
    errors_laurent = [abs(reflectance(n, "laurent") - ref) for n in orders_to_check]

    for n, e_nvm, e_laurent in zip(orders_to_check, errors_nvm, errors_laurent):
        assert e_nvm <= e_laurent + 1e-6, (n, e_nvm, e_laurent)
    assert errors_nvm[-1] < errors_laurent[-1], (errors_nvm, errors_laurent)
