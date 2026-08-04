"""Category 4 target 4.6 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): parser and
validation tests for `geometry_io`'s minimal JSON pattern-import format.
Deliberately **not** integration-tested against `Simulation`/`Layer` here
-- per the target's own "add parser tests before solver integration"
wording, this file only checks that `pattern_from_dict`/
`pattern_from_json_string`/`pattern_from_json_file` produce the correct
`Pattern` (or raise a clear `ValueError`) from data; one test does feed the
result into `Layer`/`Simulation` to confirm the returned `Pattern` is a
completely ordinary, usable one via the *existing* public API -- not new
solver-side wiring.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Circle, Ellipse, Lattice, Polygon, Rectangle, Slab
from sougata_solver.geometry_io import pattern_from_dict, pattern_from_json_file, pattern_from_json_string
from sougata_solver.layer import Layer
from sougata_solver.simulation import Simulation


def _base_doc(shapes=()):
    return {
        "unit": "um",
        "background": {"eps_re": 1.0},
        "shapes": list(shapes),
    }


# ---------------------------------------------------------------------------
# Successful parsing, one per shape type
# ---------------------------------------------------------------------------


def test_pattern_from_dict_circle():
    doc = _base_doc([{"type": "circle", "center": [0.35, 0.35], "radius": 0.18, "material": {"eps_re": 12.11}}])
    pattern = pattern_from_dict(doc)
    assert len(pattern.shapes) == 1
    shape = pattern.shapes[0]
    assert isinstance(shape, Circle)
    assert shape.center == pytest.approx((0.35e-6, 0.35e-6))
    assert shape.radius == pytest.approx(0.18e-6)
    assert complex(shape.material.epsilon_tensor(1.0)[0, 0]) == pytest.approx(12.11 + 0j)


def test_pattern_from_dict_rectangle_with_angle():
    doc = _base_doc(
        [
            {
                "type": "rectangle",
                "center": [0.0, 0.0],
                "halfwidth": [0.2, 0.1],
                "angle_deg": 30.0,
                "material": {"eps_re": 4.0},
            }
        ]
    )
    shape = pattern_from_dict(doc).shapes[0]
    assert isinstance(shape, Rectangle)
    assert shape.halfwidth == pytest.approx((0.2e-6, 0.1e-6))
    assert shape.angle == pytest.approx(np.radians(30.0))


def test_pattern_from_dict_ellipse():
    doc = _base_doc([{"type": "ellipse", "center": [0.0, 0.0], "halfwidth": [0.3, 0.15], "material": {"eps_re": 4.0}}])
    shape = pattern_from_dict(doc).shapes[0]
    assert isinstance(shape, Ellipse)
    assert shape.halfwidth == pytest.approx((0.3e-6, 0.15e-6))


def test_pattern_from_dict_polygon():
    doc = _base_doc(
        [
            {
                "type": "polygon",
                "center": [0.0, 0.0],
                "vertices": [[-0.2, -0.15], [0.2, -0.15], [0.2, 0.15], [-0.2, 0.15]],
                "material": {"eps_re": 4.0, "eps_im": 0.1},
            }
        ]
    )
    shape = pattern_from_dict(doc).shapes[0]
    assert isinstance(shape, Polygon)
    assert shape.vertices[0] == pytest.approx((-0.2e-6, -0.15e-6))
    assert complex(shape.material.epsilon_tensor(1.0)[0, 0]) == pytest.approx(4.0 + 0.1j)


def test_pattern_from_dict_slab_and_default_unit_is_meters():
    doc = {
        "background": {"eps_re": 1.0},
        "shapes": [{"type": "slab", "center_x": 0.1, "halfwidth": 0.2, "material": {"eps_re": 12.11}}],
    }
    shape = pattern_from_dict(doc).shapes[0]
    assert isinstance(shape, Slab)
    assert shape.center_x == pytest.approx(0.1)  # unit defaults to "m", no scaling
    assert shape.halfwidth == pytest.approx(0.2)


def test_pattern_from_dict_empty_shapes_is_valid():
    pattern = pattern_from_dict(_base_doc())
    assert pattern.shapes == []


def test_pattern_from_dict_forwards_material_source():
    """Category 5 target 5.8: the optional `"source"` material key is
    passed through to `Material.source` unchanged."""
    doc = _base_doc(
        [
            {
                "type": "circle",
                "center": [0.0, 0.0],
                "radius": 0.1,
                "material": {"eps_re": 4.0, "source": "hand-entered test value"},
            }
        ]
    )
    shape = pattern_from_dict(doc).shapes[0]
    assert shape.material.source == "hand-entered test value"


def test_pattern_from_dict_rejects_non_string_material_source():
    doc = _base_doc([{"type": "circle", "center": [0.0, 0.0], "radius": 0.1, "material": {"eps_re": 4.0, "source": 42}}])
    with pytest.raises(ValueError, match="source"):
        pattern_from_dict(doc)


# ---------------------------------------------------------------------------
# JSON string / file entry points
# ---------------------------------------------------------------------------


def test_pattern_from_json_string():
    doc = _base_doc([{"type": "circle", "center": [0.0, 0.0], "radius": 0.1, "material": {"eps_re": 4.0}}])
    pattern = pattern_from_json_string(json.dumps(doc))
    assert isinstance(pattern.shapes[0], Circle)


def test_pattern_from_json_string_rejects_malformed_json():
    with pytest.raises(ValueError, match="invalid JSON"):
        pattern_from_json_string("{not valid json")


def test_pattern_from_json_file(tmp_path):
    doc = _base_doc([{"type": "circle", "center": [0.0, 0.0], "radius": 0.1, "material": {"eps_re": 4.0}}])
    path = tmp_path / "pattern.json"
    path.write_text(json.dumps(doc))
    pattern = pattern_from_json_file(str(path))
    assert isinstance(pattern.shapes[0], Circle)


# ---------------------------------------------------------------------------
# Validation: never eval/exec, always a clear ValueError on malformed input
# ---------------------------------------------------------------------------


def test_pattern_from_dict_requires_background():
    with pytest.raises(ValueError, match="background"):
        pattern_from_dict({"shapes": []})


def test_pattern_from_dict_rejects_unknown_unit():
    with pytest.raises(ValueError, match="unit"):
        pattern_from_dict({"unit": "furlongs", "background": {"eps_re": 1.0}, "shapes": []})


def test_pattern_from_dict_rejects_unknown_shape_type():
    doc = _base_doc([{"type": "hexagon", "material": {"eps_re": 1.0}}])
    with pytest.raises(ValueError, match="unknown shape type"):
        pattern_from_dict(doc)


def test_pattern_from_dict_rejects_missing_shape_field():
    doc = _base_doc([{"type": "circle", "center": [0.0, 0.0], "material": {"eps_re": 4.0}}])  # missing radius
    with pytest.raises(ValueError, match="radius"):
        pattern_from_dict(doc)


def test_pattern_from_dict_rejects_missing_material():
    doc = _base_doc([{"type": "circle", "center": [0.0, 0.0], "radius": 0.1}])
    with pytest.raises(ValueError, match="material"):
        pattern_from_dict(doc)


def test_pattern_from_dict_rejects_non_numeric_field():
    doc = _base_doc([{"type": "circle", "center": [0.0, 0.0], "radius": "big", "material": {"eps_re": 4.0}}])
    with pytest.raises(ValueError, match="expected a number"):
        pattern_from_dict(doc)


def test_pattern_from_dict_rejects_malformed_center():
    doc = _base_doc([{"type": "circle", "center": [0.0, 0.0, 0.0], "radius": 0.1, "material": {"eps_re": 4.0}}])
    with pytest.raises(ValueError, match="2-element"):
        pattern_from_dict(doc)


# ---------------------------------------------------------------------------
# The returned Pattern is a completely ordinary one, usable via the
# existing public API (no new solver-side wiring involved).
# ---------------------------------------------------------------------------


def test_imported_pattern_solves_end_to_end():
    doc = {
        "unit": "um",
        "background": {"eps_re": 1.0},
        "shapes": [{"type": "circle", "center": [0.35, 0.35], "radius": 0.18, "material": {"eps_re": 12.1104}}],
    }
    pattern = pattern_from_json_string(json.dumps(doc))
    lattice = Lattice(a=(0.7e-6, 0.0), b=(0.0, 0.7e-6))
    layer = Layer("imported", 0.46e-6, pattern=pattern)
    air = pattern.background
    sim = Simulation(lattice, [layer], num_orders=7, incidence=air, transmission=air)
    excitation = PlaneWaveExcitation(1.0e-6, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    result = sim.solve(excitation)
    assert result.reflectance() + result.transmittance() == pytest.approx(1.0, abs=1e-6)
