"""Category 10 targets 10.1-10.3, 10.6 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
`SimulationResult.complex_amplitudes()`, `.diffraction_angles()`,
`.energy_balance()`, and a frozen compact output-schema check across
uniform/1D/2D fixtures.

Target 10.4 (loss-accounting design) was already satisfied by Category 7's
`layer_absorption()` design (`decisions.md` ADR-017) -- no new test here.
Target 10.5 (per-order s/p conversion) is evaluated and explicitly
deferred, not implemented -- see `references.md`'s "Target 10.5 bounded
external-validation attempt" entry and `CONVENTIONS.md`'s Category 10
addendum.
"""

from __future__ import annotations

import math

import pytest
from oracles.fresnel import multilayer_complex_rt

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Circle, Lattice, Lattice1D, Pattern, Slab
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation


pytestmark = pytest.mark.oracle  # Category 17 target 17.1: system-tier test, cross-checked against a named external oracle

AIR = Material("air", 1.0)


# ---------------------------------------------------------------------------
# 10.1 complex_amplitudes vs. the Fresnel oracle (both polarizations)
# ---------------------------------------------------------------------------

WAVELENGTH_BARE = 0.55e-6
N1, N2 = 1.0, 1.5


def _bare_interface_result(theta_deg: float, s_amplitude: complex, p_amplitude: complex):
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))
    incidence = Material("incidence", N1**2)
    transmission = Material("transmission", N2**2)
    sim = Simulation(lattice, [], num_orders=1, incidence=incidence, transmission=transmission)
    theta = math.radians(theta_deg)
    excitation = PlaneWaveExcitation(WAVELENGTH_BARE, theta, 0.0, s_amplitude=s_amplitude, p_amplitude=p_amplitude)
    result = sim.solve(excitation)
    return result, excitation, theta


def test_complex_amplitudes_s_polarization_matches_fresnel_oracle():
    result, excitation, theta = _bare_interface_result(25.0, s_amplitude=1.0, p_amplitude=0.0)
    ex0, ey0 = excitation.incident_field_xy()
    amps = result.complex_amplitudes()[(0, 0)]

    r_oracle, t_oracle = multilayer_complex_rt(WAVELENGTH_BARE, theta, "s", N1, [], N2)
    assert amps["Ey_r"] / ey0 == pytest.approx(r_oracle, abs=1e-10)
    assert amps["Ey_t"] / ey0 == pytest.approx(t_oracle, abs=1e-10)
    assert abs(ex0) < 1e-12  # s-pol at phi=0 has no Ex component


def test_complex_amplitudes_p_polarization_matches_fresnel_oracle():
    """Same oracle, p-polarization -- confirmed to match exactly (see
    `simulation.py::complex_amplitudes`'s docstring for the non-obvious
    finding this test pins: the oracle's admittance-based `r_p` sign
    convention agrees with this solver's, even though a naively
    hand-written textbook `r_p` formula would not)."""
    result, excitation, theta = _bare_interface_result(25.0, s_amplitude=0.0, p_amplitude=1.0)
    ex0, ey0 = excitation.incident_field_xy()
    amps = result.complex_amplitudes()[(0, 0)]

    r_oracle, t_oracle = multilayer_complex_rt(WAVELENGTH_BARE, theta, "p", N1, [], N2)
    assert amps["Ex_r"] / ex0 == pytest.approx(r_oracle, abs=1e-10)
    assert amps["Ex_t"] / ex0 == pytest.approx(t_oracle, abs=1e-10)
    assert abs(ey0) < 1e-12  # p-pol at phi=0 has no Ey component


def test_complex_amplitudes_p_polarization_sign_differs_from_naive_hand_formula():
    """Documents (does not "fix") the convention-ambiguity finding: a
    naive `r_p = (n2*cos(ti) - n1*cos(tt)) / (n2*cos(ti) + n1*cos(tt))`
    formula disagrees in *sign* with both this solver and `fresnel.py`'s
    oracle (which agree with each other) -- a pre-existing ambiguity in
    how p-polarization's positive direction is chosen, not a bug."""
    theta = math.radians(25.0)
    costi = math.cos(theta)
    sint_t = N1 * math.sin(theta) / N2
    costt = math.sqrt(1 - sint_t**2)
    r_p_naive = (N2 * costi - N1 * costt) / (N2 * costi + N1 * costt)

    r_oracle, _t_oracle = multilayer_complex_rt(WAVELENGTH_BARE, theta, "p", N1, [], N2)
    assert r_oracle == pytest.approx(-r_p_naive, abs=1e-10)


# ---------------------------------------------------------------------------
# 10.2 diffraction_angles
# ---------------------------------------------------------------------------


def test_diffraction_angles_zeroth_order_matches_incidence_angle():
    lattice = Lattice((0.7e-6, 0.0), (0.0, 0.7e-6))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(0.35e-6, 0.35e-6), radius=0.14e-6, material=Material("si", 3.48**2))])
    sim = Simulation(lattice, [Layer("pillar", 0.3e-6, pattern=pattern)], num_orders=9, incidence=AIR, transmission=AIR)
    theta = math.radians(20.0)
    excitation = PlaneWaveExcitation(1.0e-6, theta, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    result = sim.solve(excitation)

    angles = result.diffraction_angles()[(0, 0)]
    assert angles["theta_r"] == pytest.approx(theta, abs=1e-9)
    # incidence == transmission medium here, so the zeroth transmitted
    # order propagates straight through at the same angle.
    assert angles["theta_t"] == pytest.approx(theta, abs=1e-9)
    assert angles["phi_r"] == pytest.approx(0.0, abs=1e-9)


def test_diffraction_angles_evanescent_order_reports_none_theta_but_defined_phi():
    """Reuses `tests/test_mode_classification.py`'s already-established
    Rayleigh-threshold fixture: above `lambda_threshold`, the `(1,0)`
    order is evanescent on both sides at this wavelength (`n_trans`'s
    threshold is longer than air's own, so 1.1x the glass threshold is
    also well past air's) -- `theta_r`/`theta_t` must both be `None` there
    (never a fabricated angle), while `phi_r`/`phi_t` (a purely geometric,
    q-independent quantity) stay defined. The zeroth order, in contrast,
    is always propagating (normal incidence) and must report a real
    `theta_r`."""
    period = 1e-6
    n_trans = 1.5
    lattice = Lattice((period, 0.0), (0.0, period))
    glass = Material("glass", n_trans**2)
    layer = Layer("filler", 0.1e-6, material=Material("filler", 2.0))
    sim = Simulation(lattice, [layer], num_orders=9, incidence=AIR, transmission=glass, truncation="circular")

    lambda_threshold = n_trans * period / 1
    excitation = PlaneWaveExcitation(lambda_threshold * 1.1, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    result = sim.solve(excitation)

    classification = result.order_classification()
    assert classification[(1, 0)]["incidence"] == "evanescent"
    assert classification[(1, 0)]["transmission"] == "evanescent"

    angles = result.diffraction_angles()
    assert angles[(1, 0)]["theta_r"] is None
    assert angles[(1, 0)]["theta_t"] is None
    assert angles[(1, 0)]["phi_r"] is not None
    assert angles[(1, 0)]["phi_t"] is not None
    assert angles[(0, 0)]["theta_r"] is not None


def test_diffraction_angles_matches_classical_1d_grating_equation():
    """`sin(theta_m) = sin(theta_inc) - m*wavelength/period` for a 1D
    grating along x -- an independent, textbook cross-check (not derived
    from `SimulationResult.diffraction_angles`'s own formula)."""
    period = 0.7e-6
    wavelength = 1.0e-6
    theta_inc = math.radians(10.0)
    lattice = Lattice1D(period)
    pattern = Pattern(background=AIR)
    pattern.add(Slab(center_x=0.0, halfwidth=0.15 * period, material=Material("si", 3.48**2)))
    sim = Simulation(lattice, [Layer("grating", 0.3e-6, pattern=pattern)], num_orders=5, incidence=AIR, transmission=AIR)
    excitation = PlaneWaveExcitation(wavelength, theta_inc, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    result = sim.solve(excitation)

    angles = result.diffraction_angles()
    classification = result.order_classification()
    for (g1, g2), entry in angles.items():
        if g2 != 0:
            continue
        if classification[(g1, g2)]["incidence"] != "propagating":
            continue
        expected_sin_theta = math.sin(theta_inc) - g1 * wavelength / period
        if abs(expected_sin_theta) > 1.0:
            continue
        expected_theta = math.asin(expected_sin_theta)
        assert abs(entry["theta_r"]) == pytest.approx(abs(expected_theta), abs=1e-6)


# ---------------------------------------------------------------------------
# 10.3 energy_balance
# ---------------------------------------------------------------------------


def test_energy_balance_lossless_residual_near_zero():
    lattice = Lattice((0.7e-6, 0.0), (0.0, 0.7e-6))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(0.35e-6, 0.35e-6), radius=0.14e-6, material=Material("si", 3.48**2))])
    sim = Simulation(lattice, [Layer("pillar", 0.3e-6, pattern=pattern)], num_orders=9, incidence=AIR, transmission=AIR)
    excitation = PlaneWaveExcitation(0.6e-6, math.radians(10.0), 0.0, s_amplitude=1.0, p_amplitude=0.0)
    result = sim.solve(excitation)

    balance = result.energy_balance()
    assert balance["incident"] == 1.0
    assert balance["reflected"] == pytest.approx(result.reflectance())
    assert balance["transmitted"] == pytest.approx(result.transmittance())
    assert abs(balance["residual"]) < 1e-8


def test_energy_balance_lossy_matches_layer_absorption():
    lattice = Lattice((0.7e-6, 0.0), (0.0, 0.7e-6))
    metal = Material("metal", -396 + 80j)
    pattern = Pattern(background=AIR, shapes=[Circle(center=(0.35e-6, 0.35e-6), radius=0.21e-6, material=metal)])
    sim = Simulation(lattice, [Layer("metal_pillar", 0.05e-6, pattern=pattern)], num_orders=25, incidence=AIR, transmission=AIR)
    excitation = PlaneWaveExcitation(1.0e-6, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    result = sim.solve(excitation)

    balance = result.energy_balance()
    assert balance["absorbed"] == pytest.approx(sum(result.layer_absorption()))
    assert abs(balance["residual"]) < 1e-6
    assert balance["absorbed"] > 0


# ---------------------------------------------------------------------------
# 10.6 Frozen compact output schema (uniform / 1D / 2D)
# ---------------------------------------------------------------------------


def _output_schema(result) -> dict:
    """The compact, cross-consistent public-output fixture target 10.6
    freezes: every public `SimulationResult` output method, keyed
    consistently, so a future accidental schema change (a renamed key, a
    dropped order) is caught here rather than discovered downstream."""
    return {
        "reflectance": result.reflectance(),
        "transmittance": result.transmittance(),
        "energy_balance": result.energy_balance(),
        "diffraction_efficiencies": result.diffraction_efficiencies(),
        "complex_amplitudes": result.complex_amplitudes(),
        "diffraction_angles": result.diffraction_angles(),
        "order_classification": result.order_classification(),
        "layer_absorption": result.layer_absorption(),
    }


def _assert_schema_is_self_consistent(schema: dict, num_orders: int, num_interior_layers: int) -> None:
    de_keys = set(schema["diffraction_efficiencies"].keys())
    amp_keys = set(schema["complex_amplitudes"].keys())
    ang_keys = set(schema["diffraction_angles"].keys())
    cls_keys = set(schema["order_classification"].keys())
    assert de_keys == amp_keys == ang_keys == cls_keys
    assert len(de_keys) == num_orders

    for key, (de_r, de_t) in schema["diffraction_efficiencies"].items():
        amps = schema["complex_amplitudes"][key]
        assert set(amps.keys()) == {"Ex_r", "Ey_r", "Ex_t", "Ey_t"}
        angs = schema["diffraction_angles"][key]
        assert set(angs.keys()) == {"theta_r", "phi_r", "theta_t", "phi_t"}

    assert sum(r for r, _t in schema["diffraction_efficiencies"].values()) == pytest.approx(schema["reflectance"], abs=1e-8)
    assert sum(t for _r, t in schema["diffraction_efficiencies"].values()) == pytest.approx(schema["transmittance"], abs=1e-8)
    assert schema["energy_balance"]["reflected"] == pytest.approx(schema["reflectance"])
    assert schema["energy_balance"]["transmitted"] == pytest.approx(schema["transmittance"])
    assert len(schema["layer_absorption"]) == num_interior_layers


def test_output_schema_uniform_thin_film():
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))
    sio2 = Material("sio2", 1.46**2)
    sim = Simulation(lattice, [Layer("sio2", 0.1e-6, material=sio2)], num_orders=1, incidence=AIR, transmission=Material("si", 3.48**2))
    excitation = PlaneWaveExcitation(0.6e-6, math.radians(10.0), 0.0, s_amplitude=1.0, p_amplitude=0.0)
    schema = _output_schema(sim.solve(excitation))
    _assert_schema_is_self_consistent(schema, num_orders=1, num_interior_layers=1)


def test_output_schema_1d_trench():
    lattice = Lattice1D(0.7e-6)
    pattern = Pattern(background=AIR)
    pattern.add(Slab(center_x=0.0, halfwidth=0.15e-6, material=Material("si", 3.48**2)))
    sim = Simulation(lattice, [Layer("grating", 0.3e-6, pattern=pattern)], num_orders=5, incidence=AIR, transmission=AIR)
    excitation = PlaneWaveExcitation(1.0e-6, math.radians(10.0), 0.0, s_amplitude=1.0, p_amplitude=0.0)
    schema = _output_schema(sim.solve(excitation))
    _assert_schema_is_self_consistent(schema, num_orders=5, num_interior_layers=1)


def test_output_schema_2d_pillar():
    lattice = Lattice((0.7e-6, 0.0), (0.0, 0.7e-6))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(0.35e-6, 0.35e-6), radius=0.14e-6, material=Material("si", 3.48**2))])
    sim = Simulation(lattice, [Layer("pillar", 0.3e-6, pattern=pattern)], num_orders=9, incidence=AIR, transmission=AIR)
    excitation = PlaneWaveExcitation(0.6e-6, math.radians(10.0), 0.0, s_amplitude=1.0, p_amplitude=0.0)
    schema = _output_schema(sim.solve(excitation))
    _assert_schema_is_self_consistent(schema, num_orders=9, num_interior_layers=1)
