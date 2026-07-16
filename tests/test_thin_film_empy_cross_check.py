"""Second Phase 1 validation gate: the actual SiO2-on-Si structure from
structures/thin_film/sio2_on_si_thin_film.py, cross-checked against a TMM
oracle transcribed from the vendored EMpy reference implementation
(tests/oracles/empy_tmm.py) -- independent both in source code (a real
published library, not sougata_solver-derived) and in derivation from
tests/oracles/fresnel.py (written from scratch). Two independently-sourced
oracles agreeing is stronger evidence than either alone.
"""

import math

import pytest
from oracles.empy_tmm import isotropic_multilayer_rt

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation

# Same placeholder indices structures/thin_film/sio2_on_si_thin_film.py falls
# back to when no NK_FILE CSVs are found.
SI_INDEX = 3.9 + 0.02j
SIO2_INDEX = 1.46
SIO2_THICKNESS = 50e-9
SI_THICKNESS = 12e-6


def _run(wavelength, theta, polarization):
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))
    air = Material("air", 1.0)
    si = Material("Si", SI_INDEX**2)
    sio2 = Material("SiO2", SIO2_INDEX**2)
    layers = [
        Layer("SiO2", SIO2_THICKNESS, material=sio2),
        Layer("Si", SI_THICKNESS, material=si),
    ]
    sim = Simulation(lattice, layers, num_orders=1, incidence=air, transmission=air)

    s_amp = 1.0 if polarization == "s" else 0.0
    p_amp = 1.0 if polarization == "p" else 0.0
    excitation = PlaneWaveExcitation(wavelength, theta, 0.0, s_amplitude=s_amp, p_amplitude=p_amp)
    result = sim.solve(excitation)
    r_sougata_solver = result.reflectance()
    t_sougata_solver = result.transmittance()

    r_oracle, t_oracle = isotropic_multilayer_rt(
        wavelength,
        theta,
        n_incidence=1.0,
        layers=[(SIO2_INDEX, SIO2_THICKNESS), (SI_INDEX, SI_THICKNESS)],
        n_substrate=1.0,
        polarization=polarization,
    )
    return r_sougata_solver, t_sougata_solver, r_oracle, t_oracle


@pytest.mark.parametrize("wavelength_nm", [400.0, 550.0, 700.0, 800.0])
@pytest.mark.parametrize("theta_deg", [0.0, 30.0, 65.0])
@pytest.mark.parametrize("polarization", ["s", "p"])
def test_sio2_on_si_matches_empy_tmm_oracle(wavelength_nm, theta_deg, polarization):
    wavelength = wavelength_nm * 1e-9
    theta = math.radians(theta_deg)
    r_sougata_solver, t_sougata_solver, r_oracle, t_oracle = _run(wavelength, theta, polarization)
    assert r_sougata_solver == pytest.approx(r_oracle, abs=1e-8)
    assert t_sougata_solver == pytest.approx(t_oracle, abs=1e-8)
