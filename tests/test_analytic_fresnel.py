"""Phase 1 validation gate: uniform multilayer stack vs. an independent
analytic Fresnel/TMM oracle (`tests/oracles/fresnel.py`, written from
scratch, not derived from EMpy or pyrcwa)."""

import math

import numpy as np
import pytest
from oracles.fresnel import multilayer_rt

from pyrcwa.excitation import PlaneWaveExcitation
from pyrcwa.geometry import Lattice
from pyrcwa.layer import Layer
from pyrcwa.materials import Material
from pyrcwa.simulation import Simulation

WAVELENGTH = 0.55e-6


def _run(n_incidence, layer_specs, n_substrate, theta, polarization):
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))
    mat_inc = Material("incidence", complex(n_incidence) ** 2)
    mat_sub = Material("substrate", complex(n_substrate) ** 2)
    layers = [
        Layer(f"layer{i}", thickness, material=Material(f"mat{i}", complex(n) ** 2))
        for i, (n, thickness) in enumerate(layer_specs)
    ]
    sim = Simulation(lattice, layers, num_orders=1, incidence=mat_inc, transmission=mat_sub)

    s_amp = 1.0 if polarization == "s" else 0.0
    p_amp = 1.0 if polarization == "p" else 0.0
    excitation = PlaneWaveExcitation(WAVELENGTH, theta, 0.0, s_amplitude=s_amp, p_amplitude=p_amp)
    result = sim.solve(excitation)
    r_pyrcwa = result.reflectance()
    t_pyrcwa = result.transmittance()

    r_oracle, t_oracle = multilayer_rt(WAVELENGTH, theta, polarization, n_incidence, layer_specs, n_substrate)
    return r_pyrcwa, t_pyrcwa, r_oracle, t_oracle


STRUCTURES = {
    "bare_interface": [],
    "single_film": [(2.0, 0.1e-6)],
    "ar_coating": [(1.46, 0.55e-6 / (4 * 1.46)), (2.35, 0.55e-6 / (4 * 2.35))],
}


@pytest.mark.parametrize("structure_name", list(STRUCTURES.keys()))
@pytest.mark.parametrize("theta_deg", [0.0, 20.0, 45.0])
@pytest.mark.parametrize("polarization", ["s", "p"])
def test_lossless_multilayer_matches_fresnel_oracle(structure_name, theta_deg, polarization):
    theta = math.radians(theta_deg)
    r_pyrcwa, t_pyrcwa, r_oracle, t_oracle = _run(
        n_incidence=1.0,
        layer_specs=STRUCTURES[structure_name],
        n_substrate=1.5,
        theta=theta,
        polarization=polarization,
    )
    assert r_pyrcwa == pytest.approx(r_oracle, abs=1e-8)
    assert t_pyrcwa == pytest.approx(t_oracle, abs=1e-8)
    assert (r_pyrcwa + t_pyrcwa) == pytest.approx(1.0, abs=1e-8)


@pytest.mark.parametrize("theta_deg", [0.0, 30.0])
@pytest.mark.parametrize("polarization", ["s", "p"])
def test_absorbing_layer_matches_fresnel_oracle(theta_deg, polarization):
    theta = math.radians(theta_deg)
    layer_specs = [(2.0 + 0.3j, 0.2e-6)]
    r_pyrcwa, t_pyrcwa, r_oracle, t_oracle = _run(
        n_incidence=1.0,
        layer_specs=layer_specs,
        n_substrate=1.5,
        theta=theta,
        polarization=polarization,
    )
    assert r_pyrcwa == pytest.approx(r_oracle, abs=1e-8)
    assert t_pyrcwa == pytest.approx(t_oracle, abs=1e-8)
    assert r_pyrcwa + t_pyrcwa < 1.0  # lossy: some power absorbed


def test_thick_absorptive_layer_is_numerically_stable():
    """S-matrix (not transfer-matrix) stacking must not blow up for a very
    thick, lossy layer where exp(+q*thickness) would overflow a naive
    transfer-matrix implementation."""
    theta = 0.0
    layer_specs = [(2.0 + 1.0j, 50e-6)]  # very thick & lossy vs. 0.55um wavelength
    r_pyrcwa, t_pyrcwa, r_oracle, t_oracle = _run(
        n_incidence=1.0,
        layer_specs=layer_specs,
        n_substrate=1.5,
        theta=theta,
        polarization="s",
    )
    assert np.isfinite(r_pyrcwa)
    assert np.isfinite(t_pyrcwa)
    assert r_pyrcwa == pytest.approx(r_oracle, abs=1e-8)
    assert t_pyrcwa == pytest.approx(t_oracle, abs=1e-6)
    assert t_pyrcwa < 1e-10  # essentially fully absorbed
