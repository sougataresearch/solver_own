"""Category 15 targets 15.3/15.4 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
malformed-configuration validation tests (15.3), plus a test that a
configuration-file-driven run reproduces an existing example (15.4).

The reproduction target is `structures/thin_film/anti_reflection_coating.py`
(single quarter-wave MgF2-on-glass coating) -- chosen over
`structures/thin_film/sio2_on_si_thin_film.py` because that second file is
excluded from this session's git operations (pre-existing unrelated local
modification, see project history), not because of anything about the
physics itself.
"""

from __future__ import annotations

import json
import math

import pytest

from sougata_solver import config
from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation

# ---------------------------------------------------------------------------
# 15.4: a configuration file reproduces an existing thin-film example
# ---------------------------------------------------------------------------

_ANTI_REFLECTION_CONFIG = {
    "unit": "m",
    "lattice": {"type": "2d", "a": [1e-6, 0.0], "b": [0.0, 1e-6]},
    "materials": {
        "air": {"eps_re": 1.0},
        "glass": {"eps_re": 2.25},
        "MgF2": {"eps_re": 1.9044},
    },
    "layers": [
        {"name": "MgF2", "thickness": 0.55e-6 / (4 * 1.38), "material": "MgF2"},
    ],
    "incidence": "air",
    "transmission": "glass",
    "num_orders": 1,
    "excitation": {"wavelength": 0.55e-6, "theta_deg": 30.0, "phi_deg": 0.0, "s_amplitude_re": 1.0, "p_amplitude_re": 0.0},
}


def _reference_anti_reflection_result():
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))
    air = Material("air", 1.0)
    glass = Material("glass", 1.5**2)
    mgf2_thickness = 0.55e-6 / (4 * 1.38)
    layers = [Layer("MgF2", mgf2_thickness, material=Material("MgF2", 1.38**2))]
    sim = Simulation(lattice, layers, num_orders=1, incidence=air, transmission=glass)
    excitation = PlaneWaveExcitation(wavelength=0.55e-6, theta=math.radians(30.0), phi=0.0, s_amplitude=1.0, p_amplitude=0.0)
    return sim.solve(excitation)


def test_config_reproduces_anti_reflection_coating_example():
    sim, excitation = config.simulation_from_dict(_ANTI_REFLECTION_CONFIG)
    result = sim.solve(excitation)
    reference = _reference_anti_reflection_result()
    assert result.reflectance() == pytest.approx(reference.reflectance(), abs=1e-12)
    assert result.transmittance() == pytest.approx(reference.transmittance(), abs=1e-12)


def test_config_from_json_string_round_trip():
    text = json.dumps(_ANTI_REFLECTION_CONFIG)
    sim, excitation = config.simulation_from_json_string(text)
    result = sim.solve(excitation)
    r, t = result.reflectance(), result.transmittance()
    assert r + t == pytest.approx(1.0, abs=1e-9)


def test_config_from_json_file_round_trip(tmp_path):
    path = tmp_path / "anti_reflection.json"
    path.write_text(json.dumps(_ANTI_REFLECTION_CONFIG), encoding="utf-8")
    sim, excitation = config.simulation_from_json_file(str(path))
    result = sim.solve(excitation)
    r, t = result.reflectance(), result.transmittance()
    assert r + t == pytest.approx(1.0, abs=1e-9)


def test_config_supports_patterned_layer_via_geometry_io_schema():
    cfg = {
        "unit": "m",
        "lattice": {"type": "2d", "a": [0.7e-6, 0.0], "b": [0.0, 0.7e-6]},
        "materials": {"air": {"eps_re": 1.0}, "si": {"eps_re": 3.48**2}},
        "layers": [
            {
                "name": "pillar",
                "thickness": 0.3e-6,
                "pattern": {
                    "background": {"eps_re": 1.0},
                    "shapes": [{"type": "circle", "center": [0.35e-6, 0.35e-6], "radius": 0.14e-6, "material": {"eps_re": 3.48**2}}],
                },
            }
        ],
        "incidence": "air",
        "transmission": "air",
        "num_orders": 9,
        "excitation": {"wavelength": 0.6e-6, "theta_deg": 0.0, "phi_deg": 0.0, "s_amplitude_re": 1.0},
    }
    sim, excitation = config.simulation_from_dict(cfg)
    result = sim.solve(excitation)
    r, t = result.reflectance(), result.transmittance()
    assert 0.0 <= r <= 1.0
    assert r + t == pytest.approx(1.0, abs=1e-6)


# ---------------------------------------------------------------------------
# 15.3: malformed-configuration validation, all raising before any solve
# ---------------------------------------------------------------------------


def _valid_config():
    return json.loads(json.dumps(_ANTI_REFLECTION_CONFIG))


def test_missing_top_level_key_raises():
    cfg = _valid_config()
    del cfg["materials"]
    with pytest.raises(ValueError, match="missing required key"):
        config.simulation_from_dict(cfg)


def test_unknown_unit_raises():
    cfg = _valid_config()
    cfg["unit"] = "furlong"
    with pytest.raises(ValueError, match="config.unit"):
        config.simulation_from_dict(cfg)


def test_unknown_lattice_type_raises():
    cfg = _valid_config()
    cfg["lattice"] = {"type": "3d"}
    with pytest.raises(ValueError, match="lattice.type"):
        config.simulation_from_dict(cfg)


def test_material_missing_eps_re_raises():
    cfg = _valid_config()
    cfg["materials"]["air"] = {"eps_im": 0.0}
    with pytest.raises(ValueError, match="missing required key"):
        config.simulation_from_dict(cfg)


def test_layer_referencing_unknown_material_raises():
    cfg = _valid_config()
    cfg["layers"][0]["material"] = "does_not_exist"
    with pytest.raises(ValueError, match="unknown material"):
        config.simulation_from_dict(cfg)


def test_layer_with_both_material_and_pattern_raises():
    cfg = _valid_config()
    cfg["layers"][0]["pattern"] = {"background": {"eps_re": 1.0}, "shapes": []}
    with pytest.raises(ValueError, match="exactly one of"):
        config.simulation_from_dict(cfg)


def test_layer_with_neither_material_nor_pattern_raises():
    cfg = _valid_config()
    del cfg["layers"][0]["material"]
    with pytest.raises(ValueError, match="exactly one of"):
        config.simulation_from_dict(cfg)


def test_incidence_referencing_unknown_material_raises():
    cfg = _valid_config()
    cfg["incidence"] = "vacuum"
    with pytest.raises(ValueError, match="unknown material"):
        config.simulation_from_dict(cfg)


def test_negative_num_orders_raises():
    cfg = _valid_config()
    cfg["num_orders"] = -1
    with pytest.raises(ValueError, match="num_orders"):
        config.simulation_from_dict(cfg)


def test_non_integer_num_orders_raises():
    cfg = _valid_config()
    cfg["num_orders"] = 1.5
    with pytest.raises(ValueError, match="num_orders"):
        config.simulation_from_dict(cfg)


def test_unknown_truncation_raises():
    cfg = _valid_config()
    cfg["truncation"] = "hexagonal"
    with pytest.raises(ValueError, match="truncation"):
        config.simulation_from_dict(cfg)


def test_non_string_material_field_raises():
    cfg = _valid_config()
    cfg["materials"]["air"]["eps_re"] = "one"
    with pytest.raises(ValueError, match="expected a number"):
        config.simulation_from_dict(cfg)


def test_empty_layers_list_raises():
    cfg = _valid_config()
    cfg["layers"] = []
    with pytest.raises(ValueError, match="non-empty"):
        config.simulation_from_dict(cfg)


def test_malformed_json_string_raises_before_construction():
    with pytest.raises(json.JSONDecodeError):
        config.simulation_from_json_string("{not valid json")


def test_validation_never_reaches_a_numerical_solve(monkeypatch):
    """Confirms target 15.3's ordering requirement directly: a bad config
    raises from `simulation_from_dict` itself, never from inside
    `Simulation.solve` -- patch `solve` to fail loudly if ever reached."""

    def _must_not_be_called(self, excitation):
        raise AssertionError("Simulation.solve was reached from an invalid config")

    monkeypatch.setattr(Simulation, "solve", _must_not_be_called)
    cfg = _valid_config()
    cfg["num_orders"] = -5
    with pytest.raises(ValueError, match="num_orders"):
        config.simulation_from_dict(cfg)
