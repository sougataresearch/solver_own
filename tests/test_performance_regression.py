"""Category 17 target 17.6 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): a
performance regression guard, added now that Category 12/13 have stable,
measured benchmark baselines (`profiling/baseline_profile.py`,
`profiling/benchmark_suite.py`) for this target's own gating condition
("add only after Category 13 has stable benchmark baselines").

**Design: relative scaling, never an absolute wall-clock threshold.**
`rules.md`'s Performance Requirements explicitly rule out hard-coded
absolute time assertions (wall-clock time is machine-dependent) --
`profiling/baseline_profile.py`'s own docstring makes the same point.
This test instead asserts a **ratio measured within a single run on a
single machine**: `time(num_orders=81) / time(num_orders=9)` for the
same 2D-pillar fixture `profiling/baseline_profile.py` already uses
(`profile_2d_pillar_eigensolve`'s exact `period`/`Circle` parameters,
reused verbatim, not a new fixture). A ratio is machine-independent in a
way an absolute time never is -- a faster or slower CPU shifts both
numerator and denominator together.

**Bound rationale**: Category 12's `design.md` "Linear-Algebra Baseline &
Factorization-Reuse Design" section measured this exact fixture's
eigensolve time growing ~160x from `num_orders=9` to `num_orders=81` (a
~9x increase in Fourier-order count) on the development machine. This
test asserts the ratio stays below `1000` -- roughly 6x headroom above
the measured value, generous enough to absorb ordinary machine-to-machine
and run-to-run variance (this is a `slow`-marked, occasionally-run guard,
not a tight per-commit gate) while still catching a genuine algorithmic
regression, e.g. an accidental `O(n^4)`-or-worse step that would blow the
ratio out to `10000`+, not a legitimate small optimization or slowdown
from an unrelated change. A **lower** bound (`ratio > 1`) is asserted
too: `num_orders=81` must never be *faster* than `num_orders=9` for the
same fixture -- if it ever were, that's not "a nice speedup," it's a sign
the timing harness itself is broken (e.g. accidentally measuring a cached
result Category 13's `_eigenmode_cache` short-circuited).
"""

from __future__ import annotations

import time

import numpy as np
import pytest

from sougata_solver.eigenmodes import solve_layer_eigenmodes_patterned
from sougata_solver.fourier_basis import truncate_fourier_orders
from sougata_solver.fourier_factorization import toeplitz_matrix
from sougata_solver.geometry import Circle, Lattice, Pattern
from sougata_solver.materials import Material

AIR = Material("air", 1.0)
SI = Material("si", 3.48**2)
WAVELENGTH = 0.6e-6
REPEATS = 3
RATIO_UPPER_BOUND = 1000.0  # ~6x headroom above the ~160x measured baseline (design.md)


def _pillar_eigensolve_time(num_orders: int) -> float:
    """Same fixture and call signature as
    `profiling/baseline_profile.py::profile_2d_pillar_eigensolve`."""
    period = 0.7e-6
    lattice = Lattice((period, 0.0), (0.0, period))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(period / 2, period / 2), radius=0.2 * period, material=SI)])
    g = truncate_fourier_orders(lattice, num_orders, "circular")
    n = len(g)
    epsilon_hat = toeplitz_matrix(pattern, lattice, g, WAVELENGTH, inverse=False)
    omega = 2 * np.pi / WAVELENGTH
    kx = np.zeros(n)
    ky = np.zeros(n)

    best = float("inf")
    for _ in range(REPEATS):
        t0 = time.perf_counter()
        solve_layer_eigenmodes_patterned(omega, kx, ky, epsilon_hat)
        best = min(best, time.perf_counter() - t0)
    return best


@pytest.mark.slow
def test_eigensolve_scaling_stays_within_measured_headroom():
    small_time = _pillar_eigensolve_time(9)
    large_time = _pillar_eigensolve_time(81)
    ratio = large_time / small_time

    assert ratio > 1.0, (
        f"num_orders=81 ({large_time:.4f}s) was not slower than num_orders=9 "
        f"({small_time:.4f}s) -- timing harness likely broken, not a real speedup"
    )
    assert ratio < RATIO_UPPER_BOUND, (
        f"eigensolve time scaling from num_orders=9 to 81 grew {ratio:.1f}x, "
        f"exceeding the {RATIO_UPPER_BOUND:.0f}x regression bound (~160x measured "
        "baseline in design.md, ~6x headroom) -- possible algorithmic regression"
    )
