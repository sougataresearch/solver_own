"""Category 2 target 2.5 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): one lossy,
high-contrast stress fixture through the **full** `Simulation.solve()`
pipeline (eigensolve -> S-matrix cascade -> R/T), not just the eigenoperator
-level cross-check Phase 4b already did
(`tests/test_2d_pillar_stress.py`) -- this closes a real gap, since Phase
4b's stress sweep never actually called `Simulation.solve()` end to end for
its lossy cases, only compared eigenvalues against the RCWA.jl oracle.

Two cases, per the target's "assert either valid conservation or a
documented numerical failure" wording:

- an isotropic lossy metal-like pillar, through the full pipeline;
- a lossy, in-plane-coupled anisotropic pillar (Category 1 target 1.6's
  `solve_layer_eigenmodes_patterned_inplane`), which has **never** been
  stress-tested with a high-contrast lossy material before (target 1.6's
  own test file only used lossless/Hermitian tensors, deliberately, to keep
  energy conservation meaningful there -- see `test_anisotropic_patterned.py`).

**Finding made while writing this fixture (not a solver bug -- a sign-
convention trap, worth recording so it isn't repeated).** The first attempt
at the isotropic case reused Phase 4b's `n = -20+2j` index verbatim
(`tests/test_2d_pillar_stress.py`, labeled "lossy-metal-like" there) and
got `R+T` up to ~17 -- badly non-passive. `n = -20+2j` squares to
`eps = 396 - 80j`, i.e. `Im(eps) < 0`. Per `CONVENTIONS.md`'s documented
phasor convention (`d/dt -> -i*omega`, the standard physics convention
where a passive/absorbing medium has `Im(eps) > 0`), that index is
actually a **gain** medium in this project's convention, not a lossy one --
`R+T > 1` was therefore the numerically-correct answer for the input given,
not a bug. Phase 4b's own test never caught this because it only compared
eigenvalues against an oracle and checked condition numbers, never called
`Simulation.solve()` for R/T -- exactly the coverage gap this target closes.
Phase 4b's docstring label is left uncorrected (per `rules.md` AI Coding
Rule 3, that file's existing passing tests are not touched), but the fixture
below uses a corrected-sign lossy metal (`Im(eps) > 0`) instead, confirmed
directly (not assumed) to give `R + T < 1` before trusting the passivity
assertion.

`testing.md`'s Physical-Invariant Testing tier defines the lossy energy
identity as `R + T + sum(DE) + A = 1`, with `A` computed from the imaginary
part of the layer permittivities (Poynting-flux divergence) -- that `A`
computation is **not implemented** in this project yet (layer-wise
absorption is `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 7 targets 7.5/
7.6, still open), so the full identity cannot be checked here. What *can*
be checked without it, and is checked below, is the necessary (not
sufficient) consequence of passivity for a lossy medium: `R >= 0`,
`T >= 0`, and `R + T <= 1` (a passive/lossy structure can only remove
power, never add it) -- a genuine, if weaker, correctness signal, honestly
scoped to what's actually implemented rather than a fabricated full-balance
check.
"""

from __future__ import annotations

import numpy as np
import pytest

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Circle, Lattice, Pattern
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation

PERIOD = 0.7
WAVELENGTH = 1.0
AIR = Material("air", 1.0)


def _assert_passive_and_finite(result) -> None:
    r, t = result.reflectance(), result.transmittance()
    assert np.isfinite(r) and np.isfinite(t), f"R/T not finite (R={r}, T={t}) -- see troubleshooting.md for known NaN causes"
    assert r >= -1e-8, f"R={r} is negative beyond numerical tolerance"
    assert t >= -1e-8, f"T={t} is negative beyond numerical tolerance"
    assert r + t <= 1.0 + 1e-6, f"R+T={r + t} exceeds 1 -- a lossy/passive medium cannot amplify"


def test_isotropic_lossy_metal_pillar_full_pipeline():
    """`eps = -396+80j`: same magnitude of index contrast as Phase 4b's
    eigenoperator-level stress sweep
    (`tests/test_2d_pillar_stress.py::test_high_contrast_pillar_stress_cases`'s
    `n = -20+2j` case), but with the sign of the imaginary part corrected
    for this project's `Im(eps) > 0 == lossy` convention (see the module
    docstring's finding) -- now run through the full `Simulation.solve()`
    -> R/T pipeline that sweep never exercised."""
    metal = Material("metal", -396 + 80j)
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.3 * PERIOD, material=metal)])
    layer = Layer("metal_pillar", 0.3, pattern=pattern)
    sim = Simulation(lattice, [layer], num_orders=25, incidence=AIR, transmission=AIR)
    excitation = PlaneWaveExcitation(WAVELENGTH, theta=0.0, phi=0.0, s_amplitude=1.0, p_amplitude=0.0)

    result = sim.solve(excitation)
    _assert_passive_and_finite(result)

    diagnostics = result.all_modes[1].diagnostics  # the patterned layer, index 1 (0 is incidence half-space)
    assert diagnostics is not None
    assert diagnostics.cond_phi < 1e4, "unexpectedly ill-conditioned for a case Phase 4b already characterized"


def test_anisotropic_lossy_high_contrast_pillar_full_pipeline():
    """First stress test (any tolerance level) of the target-1.6 patterned
    anisotropic eigensolver with a genuinely lossy, high-contrast,
    non-Hermitian (absorbing) tensor -- target 1.6's own test file
    deliberately stuck to lossless/Hermitian tensors so R+T=1 would hold;
    this fixture intentionally goes past that to see whether the solver
    stays well-behaved (finite, passive) once absorption is introduced."""
    tensor = np.array(
        [[-15 + 3j, 2 + 1j, 0], [2 + 1j, -12 + 4j, 0], [0, 0, -18 + 2j]], dtype=complex
    )
    lossy_aniso = Material.from_permittivity_tensor("lossy_aniso", tensor)
    lattice = Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))
    pattern = Pattern(
        background=AIR, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.25 * PERIOD, material=lossy_aniso)]
    )
    layer = Layer("lossy_aniso_pillar", 0.3, pattern=pattern)
    sim = Simulation(lattice, [layer], num_orders=13, incidence=AIR, transmission=AIR)
    excitation = PlaneWaveExcitation(WAVELENGTH, theta=10.0, phi=0.0, s_amplitude=1.0, p_amplitude=1.0)

    try:
        result = sim.solve(excitation)
    except np.linalg.LinAlgError as exc:
        pytest.xfail(f"documented numerical failure (LinAlgError) for this stress fixture: {exc}")
        return

    _assert_passive_and_finite(result)
