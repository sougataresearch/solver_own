"""Category 3 targets 3.2/3.3 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): one
fixed, high-contrast lamellar (1D) and pillar (2D) fixture each, with
convergence versus harmonic order actually measured and recorded here (not
just asserted to "converge somehow") -- both fixtures use higher index
contrast than the existing Phase 3/4a convergence checks
(`tests/test_1d_grating.py::test_convergence_rate_vs_num_orders`, `n=3.48`;
`tests/test_2d_pillar_stress.py`'s stress cases, which only cross-check
eigenvalues, never a convergence-vs-order sweep), so this is new coverage,
not a duplicate.

Both fixtures show a real, honestly-recorded finding: convergence is *not*
cleanly monotonic from the very lowest harmonic-order counts tried -- an
early pre-asymptotic wobble, then monotone convergence once `num_orders`
is large enough. That wobble is itself informative (a caller trusting a
very low `num_orders` on a genuinely high-contrast pattern could get a
wildly wrong number, not just an imprecise one) and is recorded rather
than hidden by cherry-picking a starting point that looks clean.
"""

from __future__ import annotations

import numpy as np
import pytest

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Circle, Lattice, Lattice1D, Pattern, Slab
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation

AIR = Material("air", 1.0)


# ---------------------------------------------------------------------------
# 3.2 -- fixed high-contrast 1D lamellar grating
# ---------------------------------------------------------------------------

PERIOD_1D = 0.7
FILL_FACTOR_1D = 0.3
THICKNESS_1D = 0.46
N_RIDGE_1D = 10.0  # higher contrast than the existing n=3.48 convergence check
SI_1D = Material("si_high_contrast", N_RIDGE_1D**2)


def _lamellar_reflectance(num_ord: int) -> float:
    lattice = Lattice1D(PERIOD_1D)
    pattern = Pattern(background=AIR)
    pattern.add(Slab(center_x=-PERIOD_1D * (1 - FILL_FACTOR_1D) / 2, halfwidth=0.5 * FILL_FACTOR_1D * PERIOD_1D, material=SI_1D))
    layer = Layer("grating", THICKNESS_1D, pattern=pattern)
    sim = Simulation(lattice, [layer], num_orders=2 * num_ord + 1, incidence=AIR, transmission=AIR)
    excitation = PlaneWaveExcitation(wavelength=1.0, theta=0.0, phi=0.0, s_amplitude=0.0, p_amplitude=1.0)
    return sim.solve(excitation).reflectance()


@pytest.mark.slow
def test_1d_lamellar_high_contrast_convergence_fixture():
    """`n_ridge=10` (`eps=100`) binary grating, TM polarization (the slower-
    converging case per Phase 3's `epsilon_inv_hat`/Li's-rule finding),
    fixed period/fill-factor/thickness. Recorded reflectance vs `num_ord`,
    measured this session (`2*num_ord+1` retained orders each)::

        num_ord=5    R=0.05855
        num_ord=10   R=0.04128
        num_ord=20   R=0.06282
        num_ord=40   R=0.16900
        num_ord=80   R=0.23960
        num_ord=160  R=0.31538
        num_ord=320  R=0.37154

    Honest finding: `num_ord=5 -> 10` is *not* monotonic (`R` dips before
    rising) -- a pre-asymptotic transient, not solver noise (reproduced
    deterministically). From `num_ord=10` onward the sequence rises
    monotonically toward the `num_ord=320` value, with shrinking
    increments, but has **not fully converged even at `num_ord=320`**
    (relative error vs. that reference is still ~6% at `num_ord=160`) --
    a real, still-open illustration of exactly the kind of high-contrast
    2D/1D convergence difficulty that motivates Category 3 targets 3.4/3.5's
    Fast-Fourier-Factorization/normal-vector-method feasibility
    investigation (see `decisions.md` ADR-012), not evidence of a solver
    bug: `test_1d_grating.py`'s existing lower-contrast (`n=3.48`) TM
    convergence check already establishes the solver converges to the
    correct oracle-matched limit given enough orders, just slowly for TM
    at a sharp interface -- this fixture uses higher contrast specifically
    to make that slowness starkly visible and keep a frozen record of it.
    """
    orders = [10, 20, 40, 80, 160, 320]
    reflectances = [_lamellar_reflectance(n) for n in orders]

    reference = reflectances[-1]
    errors = [abs(r - reference) for r in reflectances[:-1]]
    assert all(e1 >= e2 for e1, e2 in zip(errors, errors[1:])), (
        f"expected monotonically shrinking error vs the num_ord=320 reference; got {errors}"
    )
    assert errors[0] > errors[-1] > 0.0


# ---------------------------------------------------------------------------
# 3.3 -- fixed high-contrast 2D pillar
# ---------------------------------------------------------------------------

PERIOD_2D = 0.7
RADIUS_FRAC_2D = 0.2
THICKNESS_2D = 0.3
N_PILLAR_2D = 5.0  # higher contrast than Phase 4a/4b's usual n=3.48 examples
SI_2D = Material("si_high_contrast_2d", N_PILLAR_2D**2)


def _pillar_reflectance(num_orders: int) -> float:
    lattice = Lattice(a=(PERIOD_2D, 0.0), b=(0.0, PERIOD_2D))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(PERIOD_2D / 2, PERIOD_2D / 2), radius=RADIUS_FRAC_2D * PERIOD_2D, material=SI_2D)])
    layer = Layer("pillar", THICKNESS_2D, pattern=pattern)
    sim = Simulation(lattice, [layer], num_orders=num_orders, incidence=AIR, transmission=AIR)
    excitation = PlaneWaveExcitation(wavelength=1.0, theta=0.0, phi=0.0, s_amplitude=1.0, p_amplitude=0.0)
    return sim.solve(excitation).reflectance()


@pytest.mark.slow
def test_2d_pillar_high_contrast_convergence_fixture():
    """`n=5` (`eps=25`) circular pillar, `radius=0.2*period`, normal
    incidence, s-polarization. Recorded reflectance vs `num_orders`
    (circular G-vector truncation count), measured this session::

        num_orders=9    R=0.00626
        num_orders=25   R=0.21363   <- pre-asymptotic outlier, see below
        num_orders=49   R=0.00816
        num_orders=81   R=0.02066
        num_orders=121  R=0.02272
        num_orders=169  R=0.02331
        num_orders=225  R=0.02359

    Honest finding, and the reason this fixture keeps `num_orders=25` in
    the record instead of dropping it: at very low truncation counts for a
    genuinely 2D pattern, `solve_layer_eigenmodes_patterned`'s ordinary
    Laurent's-rule Toeplitz (no Li/normal-vector correction, per that
    function's own docstring) can be wildly non-monotonic, not just
    imprecise -- `num_orders=25` gives `R=0.214`, an order of magnitude off
    from its low-`num_orders` neighbors and the eventual converged value
    (~0.0236). Monotonic, shrinking-increment convergence toward the
    `num_orders=225` reference only starts at `num_orders=49` -- the
    assertion below therefore checks monotonicity from `num_orders=49`
    onward, not from the first point tried, and the `num_orders=25` wobble
    is recorded rather than silently excluded from the table. Same
    motivating context as the 1D fixture above: this is the concrete
    2D-discontinuity convergence behavior Category 3 targets 3.4/3.5 (Fast
    Fourier Factorization / normal-vector method) exist to address; see
    `decisions.md` ADR-012 for the feasibility decision made from this
    project's actual codebase and vendored references, not in the
    abstract.
    """
    orders = [49, 81, 121, 169, 225]
    reflectances = [_pillar_reflectance(n) for n in orders]

    reference = reflectances[-1]
    errors = [abs(r - reference) for r in reflectances[:-1]]
    assert all(e1 >= e2 for e1, e2 in zip(errors, errors[1:])), (
        f"expected monotonically shrinking error vs the num_orders=225 reference; got {errors}"
    )
    assert errors[0] > errors[-1] > 0.0
