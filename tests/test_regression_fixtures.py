"""Category 17 target 17.4 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): a
compact, trusted, frozen reference spectrum compared against the current
solver's output on every run -- a snapshot regression guard, distinct
from (not a replacement for) the oracle-comparison tests in
`test_analytic_fresnel.py`/`test_thin_film_empy_cross_check.py`.

**Provenance**: `tests/regression_fixtures/thin_film_ar_coating_reference.npz`
was generated once, by hand, from the same quarter-wave MgF2-on-glass
anti-reflection-coating parameters as
`structures/thin_film/anti_reflection_coating.py` (air/MgF2 (`n=1.38`,
quarter-wave thick at `550 nm`)/glass (`n=1.5`), `theta=30 deg`,
s-polarized, `num_orders=1`), swept over 21 wavelength points from `400`
to `800 nm`. The generating script is reproduced in this file's
`_regenerate_fixture` function (not run automatically -- see its
docstring) so the fixture's exact provenance is never separated from the
test that consumes it.

**Tolerance rationale**: `abs=1e-10`, far tighter than any of this
project's oracle-comparison tolerances (`testing.md`'s Validation
Report tabulates those at `1e-6`-`1e-9`). This is deliberate: a
regression-snapshot test is checking for **bit-for-bit-scale
reproducibility of a fixed deterministic computation**, not agreement
with an independent physical model -- the uniform-multilayer solve path
exercised here is *itself* already independently oracle-validated
elsewhere (`test_analytic_fresnel.py`'s Fresnel-oracle comparison,
`test_thin_film_empy_cross_check.py`'s EMpy-oracle comparison, both
exercising the identical code path for a general uniform stack); this
fixture's job is only to catch an *unintended future change* to that
already-validated path (a refactor that accidentally perturbs a formula,
a dependency upgrade that changes floating-point rounding behavior
beyond what `1e-10` tolerates), not to re-establish physical correctness
from scratch.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation

FIXTURE_PATH = Path(__file__).parent / "regression_fixtures" / "thin_film_ar_coating_reference.npz"


def _build_simulation() -> Simulation:
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))
    air = Material("air", 1.0)
    glass = Material("glass", 1.5**2)
    mgf2_thickness = 0.55e-6 / (4 * 1.38)
    layers = [Layer("MgF2", mgf2_thickness, material=Material("MgF2", 1.38**2))]
    return Simulation(lattice, layers, num_orders=1, incidence=air, transmission=glass)


def _regenerate_fixture() -> None:
    """Not called by any test -- reproduced here purely as the
    fixture's documented provenance (target 17.4's "with provenance"
    requirement), so regenerating it later (should the reference
    intentionally need to change) doesn't require reverse-engineering the
    original parameters from the `.npz` file alone."""
    sim = _build_simulation()
    wavelengths = np.linspace(0.4e-6, 0.8e-6, 21)
    theta_deg = 30.0
    reflectance = np.zeros_like(wavelengths)
    transmittance = np.zeros_like(wavelengths)
    for i, wavelength in enumerate(wavelengths):
        excitation = PlaneWaveExcitation(float(wavelength), math.radians(theta_deg), 0.0, s_amplitude=1.0, p_amplitude=0.0)
        result = sim.solve(excitation)
        reflectance[i] = result.reflectance()
        transmittance[i] = result.transmittance()
    np.savez(FIXTURE_PATH, wavelengths_m=wavelengths, theta_deg=theta_deg, reflectance=reflectance, transmittance=transmittance)


def test_fixture_file_exists():
    assert FIXTURE_PATH.exists(), f"{FIXTURE_PATH} is missing -- see _regenerate_fixture for how to recreate it"


def test_current_solver_matches_frozen_reference_spectrum():
    with np.load(FIXTURE_PATH) as data:
        wavelengths = data["wavelengths_m"]
        theta_deg = float(data["theta_deg"])
        expected_reflectance = data["reflectance"]
        expected_transmittance = data["transmittance"]

    sim = _build_simulation()
    reflectance = np.zeros_like(wavelengths)
    transmittance = np.zeros_like(wavelengths)
    for i, wavelength in enumerate(wavelengths):
        excitation = PlaneWaveExcitation(float(wavelength), math.radians(theta_deg), 0.0, s_amplitude=1.0, p_amplitude=0.0)
        result = sim.solve(excitation)
        reflectance[i] = result.reflectance()
        transmittance[i] = result.transmittance()

    assert reflectance == pytest.approx(expected_reflectance, abs=1e-10)
    assert transmittance == pytest.approx(expected_transmittance, abs=1e-10)


def test_frozen_reference_itself_satisfies_energy_conservation():
    """A sanity check on the fixture data itself (not the live solver) --
    if the frozen numbers ever get corrupted or hand-edited, this catches
    it independently of whether the live solver still matches them."""
    with np.load(FIXTURE_PATH) as data:
        total = data["reflectance"] + data["transmittance"]
    assert total == pytest.approx(1.0, abs=1e-9)
