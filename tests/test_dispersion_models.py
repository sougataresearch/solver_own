"""Category 5 targets 5.2-5.6 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
analytic dispersion models (`Material.from_sellmeier`/`from_cauchy`/
`from_lorentz`/`from_drude`/`from_drude_lorentz`). One section per target;
see each `Material` classmethod's own docstring for the exact source
citation.
"""

from __future__ import annotations

import numpy as np
import pytest

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice
from sougata_solver.layer import Layer
from sougata_solver.materials import RAKIC_ALUMINUM, RAKIC_GOLD, RAKIC_SILVER, RAKIC_TITANIUM, LorentzOscillator, Material
from sougata_solver.simulation import Simulation


def _n(material: Material, wavelength_m: float) -> complex:
    eps = complex(material.epsilon_tensor(wavelength_m)[0, 0])
    return np.sqrt(eps)


# ---------------------------------------------------------------------------
# 5.2 Sellmeier model
# ---------------------------------------------------------------------------

# BK7 SCHOTT catalog Sellmeier coefficients (um^2), confirmed via WebSearch
# this session against refractiveindex.info / SCHOTT data -- not from memory.
_BK7_SELLMEIER = dict(b1=1.03961212, c1=0.00600069867, b2=0.231792344, c2=0.0200179144, b3=1.01046945, c3=103.560653)


def test_sellmeier_bk7_matches_published_nd():
    """BK7's independently-published `n_d = 1.5168` at the Fraunhofer
    d-line (587.56 nm), confirmed via `WebSearch` this session."""
    bk7 = Material.from_sellmeier("BK7", **_BK7_SELLMEIER)
    n = _n(bk7, 587.56e-9)
    assert n.real == pytest.approx(1.5168, abs=2e-4)
    assert n.imag == pytest.approx(0.0, abs=1e-12)  # lossless, no absorption data


@pytest.mark.parametrize("wavelength_nm", [486.1, 587.56, 656.3])
def test_sellmeier_bk7_matches_direct_formula_evaluation(wavelength_nm):
    """Independent (not calling `Material`/`from_sellmeier`) direct
    evaluation of the Sellmeier formula, at three different wavelengths --
    guards against a transcription slip that only shows up away from the
    single already-checked d-line point."""
    bk7 = Material.from_sellmeier("BK7", **_BK7_SELLMEIER)
    w = wavelength_nm * 1e-3  # nm -> um
    w2 = w * w
    b = _BK7_SELLMEIER
    n2_expected = 1.0 + b["b1"] * w2 / (w2 - b["c1"]) + b["b2"] * w2 / (w2 - b["c2"]) + b["b3"] * w2 / (w2 - b["c3"])
    assert _n(bk7, wavelength_nm * 1e-9).real == pytest.approx(np.sqrt(n2_expected), abs=1e-10)


def test_sellmeier_is_dispersive_not_constant():
    bk7 = Material.from_sellmeier("BK7", **_BK7_SELLMEIER)
    n_blue = _n(bk7, 486.1e-9).real
    n_red = _n(bk7, 656.3e-9).real
    assert n_blue > n_red  # normal dispersion: n decreases with increasing wavelength


# ---------------------------------------------------------------------------
# 5.3 Cauchy model
# ---------------------------------------------------------------------------

# EMpy's own documented worked example (materials.py:65-72): a SiN Cauchy fit.
_SIN_CAUCHY = dict(a=1.887, b=0.01929, c=1.6662e-4)


def test_cauchy_sin_matches_documented_example_value():
    sin_material = Material.from_cauchy("SiN", **_SIN_CAUCHY)
    # Direct hand evaluation of A + B/w^2 + C/w^4 at w=1.0 um.
    expected = 1.887 + 0.01929 / 1.0**2 + 1.6662e-4 / 1.0**4
    assert _n(sin_material, 1.0e-6).real == pytest.approx(expected, abs=1e-12)


def test_cauchy_two_term_form_defaults_c_to_zero():
    two_term = Material.from_cauchy("two_term", a=1.887, b=0.01929)
    three_term_c_zero = Material.from_cauchy("three_term_c_zero", a=1.887, b=0.01929, c=0.0)
    assert _n(two_term, 1.0e-6) == pytest.approx(_n(three_term_c_zero, 1.0e-6))


def test_cauchy_is_lossless():
    sin_material = Material.from_cauchy("SiN", **_SIN_CAUCHY)
    assert _n(sin_material, 1.0e-6).imag == pytest.approx(0.0, abs=1e-12)


def _wavelength_nm_to_m(wavelength_nm: float) -> float:
    return wavelength_nm * 1e-9


def _ev_to_wavelength_m(energy_ev: float) -> float:
    return _wavelength_nm_to_m(1239.8 / energy_ev)


# ---------------------------------------------------------------------------
# 5.4 Lorentz model
# ---------------------------------------------------------------------------


def test_lorentz_matches_direct_formula_evaluation():
    eps_inf, strength, omega0_ev, gamma_ev = 1.0, 20.0, 2.0, 0.3
    material = Material.from_lorentz("osc", eps_inf, strength, omega0_ev, gamma_ev)
    for energy_ev in [0.5, 1.0, 1.8, 2.0, 2.2, 3.0, 5.0]:
        w = _ev_to_wavelength_m(energy_ev)
        expected = eps_inf + strength / (omega0_ev**2 - energy_ev**2 - 1j * gamma_ev * energy_ev)
        assert complex(material.epsilon_tensor(w)[0, 0]) == pytest.approx(expected, abs=1e-9)


def test_lorentz_is_lossy_not_gain_at_resonance():
    """Causality/sign-convention check (target 5.4's explicit requirement):
    under this project's `d/dt -> -i*omega` convention, a passive oscillator
    must have `Im(eps) > 0` at resonance -- exactly the sign Category 2
    target 2.5 found a naively-reused index violate. Also checks the exact
    hand-derived value `i*strength/(gamma*omega0)`, not just the sign."""
    eps_inf, strength, omega0_ev, gamma_ev = 1.0, 20.0, 2.0, 0.3
    material = Material.from_lorentz("osc", eps_inf, strength, omega0_ev, gamma_ev)
    eps_at_resonance = complex(material.epsilon_tensor(_ev_to_wavelength_m(omega0_ev))[0, 0])
    assert eps_at_resonance.imag > 0
    expected = eps_inf + 1j * strength / (gamma_ev * omega0_ev)
    assert eps_at_resonance == pytest.approx(expected, abs=1e-9)


def test_lorentz_is_passive_across_a_wavelength_sweep():
    """`Im(eps) > 0` (or exactly `0` far from resonance) everywhere, for a
    positive damping rate -- passivity should hold at every energy, not
    just exactly at resonance."""
    material = Material.from_lorentz("osc", eps_inf=1.0, strength=20.0, omega0_ev=2.0, gamma_ev=0.3)
    for energy_ev in np.linspace(0.2, 6.0, 30):
        eps = complex(material.epsilon_tensor(_ev_to_wavelength_m(energy_ev))[0, 0])
        assert eps.imag >= -1e-12


# ---------------------------------------------------------------------------
# 5.5 Drude model
# ---------------------------------------------------------------------------


def test_drude_matches_direct_formula_evaluation():
    eps_inf, omega_p_ev, gamma_ev = 1.0, 9.0, 0.05
    material = Material.from_drude("free_electron", eps_inf, omega_p_ev, gamma_ev)
    for energy_ev in [0.5, 1.0, 2.0, 5.0]:
        w = _ev_to_wavelength_m(energy_ev)
        expected = eps_inf - omega_p_ev**2 / (energy_ev * (energy_ev + 1j * gamma_ev))
        assert complex(material.epsilon_tensor(w)[0, 0]) == pytest.approx(expected, abs=1e-9)


def test_drude_is_metallic_below_plasma_energy():
    """`Re(eps) < 0` below the plasma energy is the defining qualitative
    signature of a Drude metal (total internal reflection / evanescent
    propagation) -- a physically meaningful sanity check, not just formula
    reproduction."""
    material = Material.from_drude("free_electron", eps_inf=1.0, omega_p_ev=9.0, gamma_ev=0.05)
    eps_below_plasma = complex(material.epsilon_tensor(_ev_to_wavelength_m(2.0))[0, 0])
    assert eps_below_plasma.real < 0


def test_drude_is_lossy():
    material = Material.from_drude("free_electron", eps_inf=1.0, omega_p_ev=9.0, gamma_ev=0.05)
    eps = complex(material.epsilon_tensor(_ev_to_wavelength_m(2.0))[0, 0])
    assert eps.imag > 0


# ---------------------------------------------------------------------------
# 5.6 Drude-Lorentz composition (Rakic Lorentz-Drude metal model)
# ---------------------------------------------------------------------------


def _rakic_reference_eps(omega_p_ev, f0, gamma0_ev, oscillators, energy_ev: float) -> complex:
    """Independent (not calling `Material.from_drude_lorentz`) re-evaluation
    of `rakic.jl`'s `LorentzDrude` formula, straight from its source lines
    -- the cross-check target 5.6 requires, not a tautological re-test of
    the same code path."""
    omega_p_effective = np.sqrt(f0) * omega_p_ev
    eps = 1.0 - omega_p_effective**2 / (energy_ev * (energy_ev + 1j * gamma0_ev))
    for osc in oscillators:
        eps += osc.f * omega_p_ev**2 / (osc.omega_ev**2 - energy_ev**2 - 1j * energy_ev * osc.gamma_ev)
    return eps


@pytest.mark.parametrize(
    "preset,name", [(RAKIC_GOLD, "Au"), (RAKIC_SILVER, "Ag"), (RAKIC_ALUMINUM, "Al"), (RAKIC_TITANIUM, "Ti")]
)
def test_rakic_metal_matches_independent_reference_evaluation(preset, name):
    material = Material.from_drude_lorentz(name, *preset)
    omega_p_ev, f0, gamma0_ev, oscillators = preset
    for energy_ev in [0.5, 1.0, 1.96, 2.5, 4.0]:  # 1.96 eV ~= 632.8 nm (HeNe)
        w = _ev_to_wavelength_m(energy_ev)
        expected = _rakic_reference_eps(omega_p_ev, f0, gamma0_ev, oscillators, energy_ev)
        assert complex(material.epsilon_tensor(w)[0, 0]) == pytest.approx(expected, abs=1e-9)


@pytest.mark.parametrize(
    "preset,name", [(RAKIC_GOLD, "Au"), (RAKIC_SILVER, "Ag"), (RAKIC_ALUMINUM, "Al"), (RAKIC_TITANIUM, "Ti")]
)
def test_rakic_metals_are_passive_and_metallic_across_visible_nir(preset, name):
    """Sanity check against the well-known qualitative behavior of real
    metals in the visible/near-IR (large negative Re(eps), positive
    Im(eps)) -- Rakic's own published coefficients (Appl. Opt. 37,
    5271-5283 (1998)) are the "published, tabulated reference" this target
    and target 5.5 validate against; an exact external digit-for-digit
    cross-check (e.g. Johnson & Christy 1972 tabulated n,k) was attempted
    via `WebSearch`/`WebFetch` this session but the actual data table was
    not fetchable in this environment (interactive/JS-rendered page, same
    class of limitation `references.md`'s target-1.5 bounded search already
    documented) -- not silently skipped, recorded here and in
    `references.md`."""
    material = Material.from_drude_lorentz(name, *preset)
    for wavelength_nm in [400, 500, 600, 700, 800, 1000]:
        eps = complex(material.epsilon_tensor(_wavelength_nm_to_m(wavelength_nm))[0, 0])
        assert eps.real < 0, (name, wavelength_nm, eps)
        assert eps.imag > 0, (name, wavelength_nm, eps)


def test_drude_lorentz_zero_oscillators_reduces_to_drude():
    omega_p_ev, gamma0_ev, eps_inf = 9.0, 0.05, 1.0
    drude_only = Material.from_drude_lorentz("drude_only", omega_p_ev, f0=1.0, gamma0_ev=gamma0_ev, oscillators=())
    equivalent_drude = Material.from_drude("equivalent", eps_inf, omega_p_ev, gamma0_ev)
    for energy_ev in [0.5, 1.0, 2.0, 5.0]:
        w = _ev_to_wavelength_m(energy_ev)
        assert complex(drude_only.epsilon_tensor(w)[0, 0]) == pytest.approx(
            complex(equivalent_drude.epsilon_tensor(w)[0, 0]), abs=1e-9
        )


def test_drude_lorentz_zero_strength_oscillator_contributes_nothing():
    omega_p_ev, f0, gamma0_ev = 9.0, 0.8, 0.05
    zero_strength_osc = LorentzOscillator(f=0.0, gamma_ev=0.3, omega_ev=2.0)
    with_zero_osc = Material.from_drude_lorentz("with_zero", omega_p_ev, f0, gamma0_ev, (zero_strength_osc,))
    without_osc = Material.from_drude_lorentz("without", omega_p_ev, f0, gamma0_ev, ())
    for energy_ev in [0.5, 1.0, 2.0, 5.0]:
        w = _ev_to_wavelength_m(energy_ev)
        assert complex(with_zero_osc.epsilon_tensor(w)[0, 0]) == pytest.approx(
            complex(without_osc.epsilon_tensor(w)[0, 0]), abs=1e-9
        )


def test_rakic_gold_full_pipeline_thin_film_is_passive():
    """End-to-end usability check (not just formula correctness): a Rakic
    gold Material solves through the full `Simulation` pipeline as an
    ordinary uniform layer. Since layer-wise absorption isn't implemented
    yet (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 7, open), only the
    weaker passivity check (`R+T<=1`) is available for a lossy material --
    same precedent as Category 2 target 2.5's stress-regression fixture."""
    au = Material.from_drude_lorentz("Au", *RAKIC_GOLD)
    air = Material("air", 1.0)
    lattice = Lattice(a=(0.5e-6, 0.0), b=(0.0, 0.5e-6))
    layer = Layer("gold_film", 50e-9, material=au)
    sim = Simulation(lattice, [layer], num_orders=1, incidence=air, transmission=air)
    excitation = PlaneWaveExcitation(wavelength=600e-9, theta=0.0, phi=0.0, s_amplitude=1.0, p_amplitude=0.0)
    result = sim.solve(excitation)
    r, t = result.reflectance(), result.transmittance()
    assert np.isfinite(r) and np.isfinite(t)
    assert r >= -1e-8 and t >= -1e-8
    assert r + t <= 1.0 + 1e-6
