"""Category 7 targets 7.5/7.6 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
`SimulationResult.layer_absorption()` -- see `design.md`'s "Layer-Wise
Absorption Design" and `decisions.md` ADR-017 for the flux-divergence
formula (a composition of already-validated Category 9/Phase 7 pieces,
not a new physics formula) and why it was chosen over a volumetric
`Im(eps)*|E|^2` integral.

This closes the gap `tests/test_stress_regression.py`'s module docstring
explicitly flagged: "the full `R+T+A=1` lossy energy identity can't be
checked here [layer-wise absorption isn't implemented yet]" -- reuses that
same file's already-vetted `eps=-396+80j` lossy-metal fixture (sign-checked
against `CONVENTIONS.md`'s `d/dt -> -i*omega` passivity convention) as the
energy-balance validation target 7.5 requires before exposing an API.
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
SI = Material("si", 3.48**2)
METAL = Material("metal", -396 + 80j)  # same sign-corrected lossy fixture as test_stress_regression.py


def _lattice() -> Lattice:
    return Lattice(a=(PERIOD, 0.0), b=(0.0, PERIOD))


def _excitation() -> PlaneWaveExcitation:
    return PlaneWaveExcitation(WAVELENGTH, theta=0.0, phi=0.0, s_amplitude=1.0, p_amplitude=0.0)


def test_lossless_stack_has_zero_layer_absorption():
    """Energy conservation for a lossless stack: absorption must be zero
    for every layer, to numerical precision -- a free, oracle-independent
    correctness check on the flux-divergence formula itself."""
    pattern = Pattern(background=AIR, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.2 * PERIOD, material=SI)])
    layers = [Layer("l1", 0.3, pattern=pattern), Layer("l2", 0.2, material=SI)]
    sim = Simulation(_lattice(), layers, num_orders=9, incidence=AIR, transmission=AIR)
    result = sim.solve(_excitation())

    absorption = result.layer_absorption()
    assert len(absorption) == 2
    for a in absorption:
        assert abs(a) < 1e-8, f"lossless layer absorbed {a}, expected ~0"

    r, t = result.reflectance(), result.transmittance()
    assert r + t + sum(absorption) == pytest.approx(1.0, abs=1e-8)


def test_lossy_metal_pillar_satisfies_full_energy_balance_identity():
    """The identity `R + T + sum(layer_absorption()) == 1` finally closeable
    for a genuinely lossy structure, using the same lossy material as
    `test_stress_regression.py::test_isotropic_lossy_metal_pillar_full_pipeline`
    (that test already validated `R>=0`/`T>=0`/`R+T<=1`; this test goes
    further and checks the full identity, not just that weaker necessary
    condition).

    `thickness=0.05` (not that test's `0.3`) is deliberate, not an
    arbitrary shrink: `layer_absorption()` inherits `propagate_amplitudes`'s
    known numerical-stability envelope (`troubleshooting.md`'s "Interior
    field reconstruction can numerically blow up for thick, highly lossy,
    high-`num_orders` layers" entry) -- `max(Im(q))*thickness` must stay
    a modest double-digit-or-less exponent, or the deepest evanescent
    modes' backward-propagated amplitude (`b(z)=b_top*exp(-i*q*z)`,
    correct and already validated by Category 9's continuity tests) grows
    catastrophically and swamps the flux sum. Confirmed directly before
    picking this value: `max(Im(q))*thickness ~= 6.3` here (safe); the
    original `thickness=0.3`/`num_orders=25` combination reaches `~38`
    and gives a nonsensical `layer_absorption() ~= 573` -- not silently
    avoided, see the troubleshooting.md entry for the full account."""
    pattern = Pattern(background=AIR, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.3 * PERIOD, material=METAL)])
    layer = Layer("metal_pillar", 0.05, pattern=pattern)
    sim = Simulation(_lattice(), [layer], num_orders=25, incidence=AIR, transmission=AIR)
    result = sim.solve(_excitation())

    absorption = result.layer_absorption()
    assert len(absorption) == 1
    r, t = result.reflectance(), result.transmittance()

    assert absorption[0] > 0, "a lossy metal pillar must show positive absorption, not near-zero"
    assert r + t + sum(absorption) == pytest.approx(1.0, abs=1e-6)


def test_interior_amplitude_reconstruction_can_numerically_overflow_for_thick_lossy_layers():
    """Honest documented limitation, not a formula bug (`troubleshooting.md`):
    a thick, highly lossy, high-`num_orders` layer pushes
    `max(Im(q))*thickness` into a regime (~38 here) where
    `propagate_amplitudes`'s backward-wave exponential growth
    (`exp(-i*q*z)`, correct and validated for moderate cases by Category
    9's field-continuity tests) numerically overflows the deepest
    evanescent modes' amplitude, breaking the energy-balance identity this
    same fixture satisfies exactly at a smaller thickness (see the test
    above). Encoded as a regression guard on the *symptom*
    (`R+T+sum(A)` deviates far outside physical bounds) so a future change
    that silently "fixes" this by clipping/renormalizing values doesn't go
    unnoticed without a deliberate decision -- see `decisions.md` ADR-017's
    "known limitation" note."""
    pattern = Pattern(background=AIR, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.3 * PERIOD, material=METAL)])
    layer = Layer("metal_pillar", 0.3, pattern=pattern)
    sim = Simulation(_lattice(), [layer], num_orders=25, incidence=AIR, transmission=AIR)
    result = sim.solve(_excitation())

    absorption = result.layer_absorption()
    r, t = result.reflectance(), result.transmittance()
    # Confirms the failure mode is real and reproducible (not confirming
    # it's "correct" -- the opposite): the identity is badly violated.
    assert abs(r + t + sum(absorption) - 1.0) > 1.0


def test_lossy_multilayer_absorption_is_positive_per_lossy_layer_and_zero_for_lossless_layers():
    """A three-layer stack with one lossy patterned layer sandwiched
    between two lossless uniform layers: only the lossy layer should show
    non-negligible absorption; the lossless layers must stay at ~0."""
    lossless_top = Layer("lossless_top", 0.2, material=SI)
    lossy_pattern = Pattern(background=AIR, shapes=[Circle(center=(PERIOD / 2, PERIOD / 2), radius=0.25 * PERIOD, material=METAL)])
    lossy_middle = Layer("lossy_middle", 0.15, pattern=lossy_pattern)
    lossless_bottom = Layer("lossless_bottom", 0.2, material=SI)

    sim = Simulation(_lattice(), [lossless_top, lossy_middle, lossless_bottom], num_orders=9, incidence=AIR, transmission=AIR)
    result = sim.solve(_excitation())

    absorption = result.layer_absorption()
    assert len(absorption) == 3
    assert abs(absorption[0]) < 1e-6, f"lossless top layer absorbed {absorption[0]}, expected ~0"
    assert absorption[1] > 1e-4, "lossy middle layer should show real absorption"
    assert abs(absorption[2]) < 1e-6, f"lossless bottom layer absorbed {absorption[2]}, expected ~0"

    r, t = result.reflectance(), result.transmittance()
    assert r + t + sum(absorption) == pytest.approx(1.0, abs=1e-6)
