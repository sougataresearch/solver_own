"""Category 15 target 15.7 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): serialize
a `sweep.SweepResult` (parameter values, reflectance, transmittance) and its
metadata to a NumPy `.npz` archive.

Metadata is JSON-encoded into a plain string array (`np.array(json_text)`)
rather than stored as a pickled object array -- `np.load` on the resulting
file needs no `allow_pickle=True`, avoiding the same untrusted-deserialization
class of risk `rules.md`'s Security Rules already rule out for JSON
(`config.py`, `geometry_io.py`: never `eval`/`exec`/`pickle` on external
input). `parameter_values` must be numeric (wavelength/angle/thickness
sweeps) since `.npz` arrays need a homogeneous dtype; a discrete/labeled
sweep (e.g. the polarization sweep's `(s_amplitude, p_amplitude)` tuples)
is out of scope for this exporter -- `export_sweep_npz` raises rather than
silently truncating such data.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from sougata_solver.sweep import SweepResult


def export_sweep_npz(sweep: SweepResult, path: str | Path) -> Path:
    """Write `sweep` to `path` (a `.npz` archive) with arrays
    `parameter_values`, `reflectance`, `transmittance`, and a JSON-encoded
    `metadata` string array holding `parameter_name`, `parameter_unit`, and
    `sweep.metadata`. Returns `path` for chaining."""
    try:
        parameter_values = np.asarray(sweep.parameter_values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "export_sweep_npz requires numeric parameter_values (a discrete/"
            "labeled sweep, e.g. polarization tuples, is not supported)"
        ) from exc
    if parameter_values.ndim != 1:
        raise ValueError(
            "export_sweep_npz requires numeric parameter_values (a discrete/"
            "labeled sweep, e.g. polarization tuples, is not supported)"
        )

    metadata = {
        "parameter_name": sweep.parameter_name,
        "parameter_unit": sweep.parameter_unit,
        **sweep.metadata,
    }
    path = Path(path)
    np.savez(
        path,
        parameter_values=parameter_values,
        reflectance=sweep.reflectance(),
        transmittance=sweep.transmittance(),
        metadata=np.array(json.dumps(metadata)),
    )
    return path


def load_sweep_npz(path: str | Path) -> dict:
    """Read back a `.npz` archive written by `export_sweep_npz` as a plain
    dict (`parameter_values`, `reflectance`, `transmittance`, `metadata`) --
    no `SweepResult` reconstruction, since that would need the original
    `SimulationResult` objects this exporter deliberately does not retain.
    `allow_pickle` is never set (default `False`), matching the security
    rationale in this module's docstring."""
    with np.load(path) as data:
        return {
            "parameter_values": data["parameter_values"],
            "reflectance": data["reflectance"],
            "transmittance": data["transmittance"],
            "metadata": json.loads(str(data["metadata"])),
        }
