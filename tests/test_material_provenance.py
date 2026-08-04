"""Category 5 target 5.8 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): optional
source/citation metadata on `Material` and threading it into serialized
(`run_metadata.txt`) output. `source` is purely informational (never read
by any solver code) -- these tests check it's stored, forwarded by every
`from_*` classmethod, defaults to `None`, and can be written into
`output_paths.write_run_metadata`'s output without any change needed to
that already-generic function.
"""

from __future__ import annotations

import numpy as np
import pytest

from sougata_solver.materials import RAKIC_GOLD, Material
from sougata_solver.output_paths import write_run_metadata


def test_material_source_defaults_to_none():
    assert Material("m", 2.25).source is None


def test_material_source_stored_and_public():
    m = Material("m", 2.25, source="hand-picked test value")
    assert m.source == "hand-picked test value"


@pytest.mark.parametrize(
    "build",
    [
        lambda: Material.from_nk("m", 1.5, 0.0, source="s"),
        lambda: Material.from_sellmeier("m", 1.0, 0.006, 0.2, 0.02, 1.0, 100.0, source="s"),
        lambda: Material.from_cauchy("m", 1.5, 0.01, source="s"),
        lambda: Material.from_lorentz("m", 1.0, 20.0, 2.0, 0.3, source="s"),
        lambda: Material.from_drude("m", 1.0, 9.0, 0.05, source="s"),
        lambda: Material.from_drude_lorentz("m", *RAKIC_GOLD, source="s"),
        lambda: Material.from_permittivity_tensor("m", np.eye(3, dtype=complex) * 2.25, source="s"),
    ],
)
def test_every_from_classmethod_forwards_source(build):
    assert build().source == "s"


def test_from_nk_file_forwards_source(tmp_path):
    path = tmp_path / "nk.txt"
    path.write_text("0.5 1.5 0.0\n0.6 1.5 0.0\n0.7 1.5 0.0\n")
    material = Material.from_nk_file("m", str(path), wavelength_unit="um", source="local test file")
    assert material.source == "local test file"


def test_material_provenance_flows_into_run_metadata(tmp_path):
    """Worked example for the "serialized output" half of this target:
    `write_run_metadata` already accepts arbitrary `**params` (no change
    needed to `output_paths.py` itself), so `material.source` threads
    through as one more keyword argument. Uses `tmp_path`, not
    `output_paths.run_output_dir`, so this test doesn't leave a file behind
    in the project's real `outputs/` tree (that tree is for actual
    `structures/`-script runs, not test artifacts)."""
    gold = Material.from_drude_lorentz(
        "Au", *RAKIC_GOLD, source="Rakić et al., Appl. Opt. 37, 5271 (1998)"
    )
    metadata_path = write_run_metadata(
        tmp_path, __file__, pillar_material=gold.name, pillar_material_source=gold.source
    )
    text = metadata_path.read_text(encoding="utf-8")
    assert "pillar_material: Au" in text
    assert "pillar_material_source: Rakić et al., Appl. Opt. 37, 5271 (1998)" in text
