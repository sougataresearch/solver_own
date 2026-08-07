"""Category 15 targets 15.2-15.4 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`):
minimal JSON simulation-configuration schema, validation, and a runner.

**No new parsing dependency** (target 15.2's own wording) -- standard-
library `json` only, same discipline `geometry_io.py` already established
for Category 4 target 4.6 (never `eval`/`exec`/`pickle`, per `rules.md`
Security Rules). Patterned layers reuse `geometry_io.pattern_from_dict`
directly (its own "background"+"shapes" sub-schema embedded unchanged, not
duplicated) rather than re-inventing a second shape-JSON format.

**Validated before any numerical calculation** (target 15.3): every
function here only *constructs* `Material`/`Layer`/`Simulation`/
`PlaneWaveExcitation` objects -- it never calls `Simulation.solve()`.
Construction-time validation already raises immediately on malformed
input (missing keys, wrong types, an unknown material name) before a
caller can reach a numerical solve, the same "fail loud, fail early"
convention every other construction-time validator in this project
follows (Category 4 target 4.1, Category 7 target 7.1).

Schema (top level)::

    {
      "unit": "m" | "um" | "nm",              # optional, default "m"
      "lattice": {"type": "2d", "a": [ax,ay], "b": [bx,by]}
                | {"type": "1d", "period": p},
      "materials": {
        "<name>": {"eps_re": <num>, "eps_im": <num, optional, default 0>,
                    "source": <str, optional>}
        ...
      },
      "layers": [
        {"name": "<str>", "thickness": <num>, "material": "<material name>"}
        | {"name": "<str>", "thickness": <num>, "pattern": <geometry_io.pattern_from_dict schema>}
        ...
      ],
      "incidence": "<material name>",
      "transmission": "<material name>",
      "num_orders": <int>,
      "truncation": "circular" | "square",     # optional, default "circular"
      "excitation": {
        "wavelength": <num>,
        "theta_deg": <num>, "phi_deg": <num>,
        "s_amplitude_re": <num>, "s_amplitude_im": <num, optional, default 0>,
        "p_amplitude_re": <num>, "p_amplitude_im": <num, optional, default 0>
      }
    }

All lengths (`a`/`b`, `period`, layer `thickness`, `excitation.wavelength`)
are in `unit`, converted to meters before constructing anything -- the
same `unit`/`_UNIT_SCALE_M` convention `geometry_io.py` already uses.
"""

from __future__ import annotations

import json
import math
from typing import Any

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice, Lattice1D
from sougata_solver.geometry_io import pattern_from_dict
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation

_UNIT_SCALE_M = {"m": 1.0, "um": 1e-6, "nm": 1e-9}


def _require_keys(d: Any, keys: tuple[str, ...], context: str) -> None:
    if not isinstance(d, dict):
        raise ValueError(f"{context}: expected a JSON object, got {type(d).__name__}")
    missing = [k for k in keys if k not in d]
    if missing:
        raise ValueError(f"{context}: missing required key(s) {missing}")


def _require_number(value: Any, context: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{context}: expected a number, got {value!r}")
    return float(value)


def _material_from_dict(name: str, d: Any) -> Material:
    context = f"materials.{name}"
    _require_keys(d, ("eps_re",), context)
    eps_re = _require_number(d["eps_re"], f"{context}.eps_re")
    eps_im = _require_number(d.get("eps_im", 0.0), f"{context}.eps_im")
    source = d.get("source")
    if source is not None and not isinstance(source, str):
        raise ValueError(f"{context}.source: expected a string, got {source!r}")
    return Material(name, complex(eps_re, eps_im), source=source)


def _lattice_from_dict(d: Any, scale: float) -> Lattice | Lattice1D:
    _require_keys(d, ("type",), "lattice")
    kind = d["type"]
    if kind == "2d":
        _require_keys(d, ("a", "b"), "lattice")
        a = tuple(_require_number(v, "lattice.a[i]") * scale for v in d["a"])
        b = tuple(_require_number(v, "lattice.b[i]") * scale for v in d["b"])
        return Lattice(a, b)
    if kind == "1d":
        _require_keys(d, ("period",), "lattice")
        period = _require_number(d["period"], "lattice.period") * scale
        return Lattice1D(period)
    raise ValueError(f"lattice.type: expected '2d' or '1d', got {kind!r}")


def _lookup_material(materials: dict[str, Material], name: Any, context: str) -> Material:
    if not isinstance(name, str):
        raise ValueError(f"{context}: expected a material name (string), got {name!r}")
    if name not in materials:
        raise ValueError(f"{context}: unknown material {name!r} (declared materials: {sorted(materials)})")
    return materials[name]


def _layer_from_dict(d: Any, index: int, materials: dict[str, Material], scale: float) -> Layer:
    context = f"layers[{index}]"
    _require_keys(d, ("name", "thickness"), context)
    name = d["name"]
    if not isinstance(name, str):
        raise ValueError(f"{context}.name: expected a string, got {name!r}")
    thickness = _require_number(d["thickness"], f"{context}.thickness") * scale

    has_material = "material" in d
    has_pattern = "pattern" in d
    if has_material == has_pattern:
        raise ValueError(f"{context}: exactly one of 'material' or 'pattern' is required")
    if has_material:
        material = _lookup_material(materials, d["material"], f"{context}.material")
        return Layer(name, thickness, material=material)
    pattern = pattern_from_dict(d["pattern"])
    return Layer(name, thickness, pattern=pattern)


def _excitation_from_dict(d: Any, scale: float) -> PlaneWaveExcitation:
    _require_keys(d, ("wavelength", "theta_deg", "phi_deg"), "excitation")
    wavelength = _require_number(d["wavelength"], "excitation.wavelength") * scale
    theta = math.radians(_require_number(d["theta_deg"], "excitation.theta_deg"))
    phi = math.radians(_require_number(d["phi_deg"], "excitation.phi_deg"))
    s_re = _require_number(d.get("s_amplitude_re", 1.0), "excitation.s_amplitude_re")
    s_im = _require_number(d.get("s_amplitude_im", 0.0), "excitation.s_amplitude_im")
    p_re = _require_number(d.get("p_amplitude_re", 0.0), "excitation.p_amplitude_re")
    p_im = _require_number(d.get("p_amplitude_im", 0.0), "excitation.p_amplitude_im")
    return PlaneWaveExcitation(wavelength, theta, phi, complex(s_re, s_im), complex(p_re, p_im))


def simulation_from_dict(data: Any) -> tuple[Simulation, PlaneWaveExcitation]:
    """Category 15 targets 15.2/15.3: parse and validate a simulation
    configuration dict, returning `(Simulation, PlaneWaveExcitation)` --
    the caller still calls `.solve()` explicitly (this function never
    does), keeping "validated" and "numerically solved" as clearly
    separate steps, per target 15.3's own wording.
    """
    _require_keys(data, ("lattice", "materials", "layers", "incidence", "transmission", "num_orders", "excitation"), "config")

    unit = data.get("unit", "m")
    if unit not in _UNIT_SCALE_M:
        raise ValueError(f"config.unit: expected one of {sorted(_UNIT_SCALE_M)}, got {unit!r}")
    scale = _UNIT_SCALE_M[unit]

    lattice = _lattice_from_dict(data["lattice"], scale)

    materials_data = data["materials"]
    if not isinstance(materials_data, dict):
        raise ValueError(f"config.materials: expected a JSON object, got {type(materials_data).__name__}")
    materials = {name: _material_from_dict(name, d) for name, d in materials_data.items()}

    layers_data = data["layers"]
    if not isinstance(layers_data, list) or len(layers_data) == 0:
        raise ValueError("config.layers: expected a non-empty JSON array")
    layers = [_layer_from_dict(d, i, materials, scale) for i, d in enumerate(layers_data)]

    incidence = _lookup_material(materials, data["incidence"], "config.incidence")
    transmission = _lookup_material(materials, data["transmission"], "config.transmission")

    num_orders = data["num_orders"]
    if isinstance(num_orders, bool) or not isinstance(num_orders, int) or num_orders < 1:
        raise ValueError(f"config.num_orders: expected a positive integer, got {num_orders!r}")

    truncation = data.get("truncation", "circular")
    if truncation not in ("circular", "square"):
        raise ValueError(f"config.truncation: expected 'circular' or 'square', got {truncation!r}")

    excitation = _excitation_from_dict(data["excitation"], scale)

    simulation = Simulation(lattice, layers, num_orders, incidence, transmission, truncation)
    return simulation, excitation


def simulation_from_json_string(text: str) -> tuple[Simulation, PlaneWaveExcitation]:
    return simulation_from_dict(json.loads(text))


def simulation_from_json_file(path: str) -> tuple[Simulation, PlaneWaveExcitation]:
    with open(path, encoding="utf-8") as f:
        return simulation_from_dict(json.load(f))
