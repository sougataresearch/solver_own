"""Category 15 target 15.7: NumPy export/round-trip tests for
`export.export_sweep_npz`/`export.load_sweep_npz`."""

from __future__ import annotations

import math

import numpy as np
import pytest

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.export import export_sweep_npz, load_sweep_npz
from sougata_solver.geometry import Lattice
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation
from sougata_solver.sweep import sweep_wavelength


def _thin_film_sweep():
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))
    air = Material("air", 1.0)
    glass = Material("glass", 1.5**2)
    layers = [Layer("MgF2", 0.1e-6, material=Material("MgF2", 1.38**2))]
    sim = Simulation(lattice, layers, num_orders=1, incidence=air, transmission=glass)
    wavelengths = [0.5e-6, 0.55e-6, 0.6e-6, 0.65e-6]
    return sweep_wavelength(sim, wavelengths, theta=math.radians(10.0), phi=0.0, s_amplitude=1.0, p_amplitude=0.0)


def test_export_and_reload_round_trips_reflectance_and_transmittance(tmp_path):
    sweep = _thin_film_sweep()
    path = export_sweep_npz(sweep, tmp_path / "sweep.npz")
    assert path.exists()

    loaded = load_sweep_npz(path)
    np.testing.assert_allclose(loaded["parameter_values"], sweep.parameter_values)
    np.testing.assert_allclose(loaded["reflectance"], sweep.reflectance())
    np.testing.assert_allclose(loaded["transmittance"], sweep.transmittance())


def test_exported_metadata_contains_parameter_name_and_unit(tmp_path):
    sweep = _thin_film_sweep()
    path = export_sweep_npz(sweep, tmp_path / "sweep.npz")
    loaded = load_sweep_npz(path)
    assert loaded["metadata"]["parameter_name"] == sweep.parameter_name
    assert loaded["metadata"]["parameter_unit"] == sweep.parameter_unit


def test_exported_file_loads_without_allow_pickle(tmp_path):
    """Confirms the security rationale in export.py's docstring directly:
    a file written by export_sweep_npz must be loadable with NumPy's
    default (`allow_pickle=False`) -- this call already exercises that
    default via load_sweep_npz, but this test additionally calls
    np.load itself to make the guarantee explicit and independent of
    load_sweep_npz's own implementation."""
    sweep = _thin_film_sweep()
    path = export_sweep_npz(sweep, tmp_path / "sweep.npz")
    with np.load(path, allow_pickle=False) as data:
        assert set(data.files) == {"parameter_values", "reflectance", "transmittance", "metadata"}


def test_discrete_labeled_sweep_raises_not_silently_truncated(tmp_path):
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))
    air = Material("air", 1.0)
    glass = Material("glass", 1.5**2)
    layers = [Layer("MgF2", 0.1e-6, material=Material("MgF2", 1.38**2))]
    sim = Simulation(lattice, layers, num_orders=1, incidence=air, transmission=glass)
    excitation = PlaneWaveExcitation(0.55e-6, math.radians(10.0), 0.0, s_amplitude=1.0, p_amplitude=0.0)
    from sougata_solver.sweep import SweepResult

    result = sim.solve(excitation)
    labeled = SweepResult(
        parameter_name="polarization",
        parameter_unit="",
        parameter_values=[(1.0, 0.0), (0.0, 1.0)],
        results=[result, result],
    )
    with pytest.raises(ValueError, match="numeric parameter_values"):
        export_sweep_npz(labeled, tmp_path / "labeled.npz")
