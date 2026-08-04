"""Category 2 target 2.1 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): tests
backing `design.md`'s Failure Contract section -- one test per documented
`ValueError`/`NotImplementedError`/`LinAlgError` condition, so that table
stays honest (a passing suite here means the documented contract still
matches actual code behavior, not just what was true when it was written).

Not a test of solver *correctness* -- purely "does this input trigger the
documented failure mode." `logging.WARNING` conditions (the fourth
category in the Failure Contract) are covered by
`tests/test_2d_pillar_stress.py` (`ILL_CONDITIONED_THRESHOLD`) and
`tests/test_degeneracy_warning.py` (target 2.4's `DEGENERATE_GAP_THRESHOLD`),
not duplicated here.
"""

from __future__ import annotations

import numpy as np
import pytest

from sougata_solver.eigenmodes import (
    solve_layer_eigenmodes_1d,
    solve_layer_eigenmodes_patterned,
    solve_layer_eigenmodes_patterned_inplane,
)
from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.fourier_basis import truncate_fourier_orders
from sougata_solver.fourier_factorization import _scalar_value
from sougata_solver.geometry import Circle, Lattice
from sougata_solver.layer import Layer, LayerEigenmodes
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation
from sougata_solver.smatrix import SMatrixStack

AIR = Material("air", 1.0)
SI = Material("si", 3.48**2)


# ---------------------------------------------------------------------------
# ValueError
# ---------------------------------------------------------------------------


def test_layer_requires_exactly_one_of_material_or_pattern():
    with pytest.raises(ValueError, match="exactly one"):
        Layer("bad", 1.0)
    with pytest.raises(ValueError, match="exactly one"):
        from sougata_solver.geometry import Pattern

        Layer("bad", 1.0, material=AIR, pattern=Pattern(background=AIR, shapes=[]))


def test_smatrix_stack_requires_matching_lengths():
    modes = LayerEigenmodes(
        q=np.array([1j, 1j]), phi=np.eye(2, dtype=complex), kp=np.eye(2, dtype=complex),
        epsilon_inv=None, is_scalar_isotropic=True,
    )
    with pytest.raises(ValueError, match="same length"):
        SMatrixStack(thicknesses=[1.0], all_modes=[modes, modes])


def test_staircase_rejects_zero_slices():
    from sougata_solver.staircase import staircase_circle_layers

    with pytest.raises(ValueError, match="num_slices"):
        staircase_circle_layers((0.5, 0.5), 0.2, 0.1, 1.0, num_slices=0, shape_material=SI, background_material=AIR)


def test_1d_eigensolver_rejects_nonzero_ky():
    omega = 2 * np.pi
    with pytest.raises(ValueError, match="ky == 0"):
        solve_layer_eigenmodes_1d(
            omega, np.array([0.0, 1.0]), np.array([0.0, 0.5]), np.eye(2, dtype=complex), np.eye(2, dtype=complex)
        )


def test_patterned_eigensolver_rejects_wrong_shape_epsilon_hat():
    omega = 2 * np.pi
    with pytest.raises(ValueError, match="\\(n, n\\)"):
        solve_layer_eigenmodes_patterned(omega, np.array([0.0, 1.0]), np.array([0.0, 0.0]), np.eye(3, dtype=complex))


def test_patterned_inplane_eigensolver_rejects_wrong_shape_component():
    omega = 2 * np.pi
    kx = np.array([0.0, 1.0])
    ky = np.array([0.0, 0.0])
    good = np.eye(2, dtype=complex)
    bad = np.eye(3, dtype=complex)
    with pytest.raises(ValueError, match="epsilon_hat_xy"):
        solve_layer_eigenmodes_patterned_inplane(omega, kx, ky, good, bad, good, good, good)


def test_material_from_nk_file_rejects_empty_file(tmp_path):
    path = tmp_path / "empty.txt"
    path.write_text("not numeric data\nalso not numeric\n")
    with pytest.raises(ValueError, match="No numeric data"):
        Material.from_nk_file("bad", str(path))


def test_material_from_nk_file_rejects_wrong_column_count(tmp_path):
    path = tmp_path / "bad_columns.txt"
    path.write_text("0.5 1.5 2.5 3.5\n0.6 1.6 2.6 3.6\n")
    with pytest.raises(ValueError, match="2 or 3 columns"):
        Material.from_nk_file("bad", str(path))


def test_material_rejects_non_scalar_non_3x3_tensor():
    with pytest.raises(ValueError, match="scalar or a 3x3"):
        Material("bad", np.eye(2))


# ---------------------------------------------------------------------------
# NotImplementedError (unimplemented-phase guards, always naming the target)
# ---------------------------------------------------------------------------


def test_simulation_rejects_longitudinal_coupling_uniform_layer():
    tensor = np.array([[2.25, 0, 0.3], [0, 2.25, 0], [0.3, 0, 3.1]], dtype=complex)
    aniso = Material.from_permittivity_tensor("longitudinal", tensor)
    lattice = Lattice(a=(0.7, 0.0), b=(0.0, 0.7))
    sim = Simulation(lattice, [Layer("l", 0.5, material=aniso)], num_orders=1, incidence=AIR, transmission=AIR)
    excitation = PlaneWaveExcitation(1.0, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    with pytest.raises(NotImplementedError, match="target 1.5"):
        sim.solve(excitation)


def test_simulation_rejects_longitudinal_coupling_patterned_layer():
    from sougata_solver.geometry import Pattern

    tensor = np.array([[2.25, 0, 0.3], [0, 2.25, 0], [0.3, 0, 3.1]], dtype=complex)
    aniso = Material.from_permittivity_tensor("longitudinal", tensor)
    lattice = Lattice(a=(0.7, 0.0), b=(0.0, 0.7))
    pattern = Pattern(background=AIR, shapes=[Circle(center=(0.35, 0.35), radius=0.18, material=aniso)])
    sim = Simulation(lattice, [Layer("l", 0.5, pattern=pattern)], num_orders=1, incidence=AIR, transmission=AIR)
    excitation = PlaneWaveExcitation(1.0, 0.0, 0.0, s_amplitude=1.0, p_amplitude=0.0)
    with pytest.raises(NotImplementedError, match="target 1.5"):
        sim.solve(excitation)


def test_scalar_fourier_factorization_rejects_anisotropic_material():
    aniso = Material.from_permittivity_tensor("aniso", np.diag([2.25, 4.0, 3.1]).astype(complex))
    with pytest.raises(NotImplementedError, match="Phase 2 scope"):
        _scalar_value(aniso, wavelength=1.0, inverse=False)


def test_truncate_fourier_orders_rejects_unknown_method():
    lattice = Lattice(a=(0.7, 0.0), b=(0.0, 0.7))
    with pytest.raises(NotImplementedError, match="not implemented"):
        truncate_fourier_orders(lattice, 5, method="square")


def test_refractiveindex_formula_file_rejects_unsupported_type(tmp_path):
    path = tmp_path / "bad_formula.yml"
    path.write_text("DATA:\n  - type: formula 1\n    coefficients: 1 2 3\n")
    with pytest.raises(NotImplementedError, match="formula 4"):
        Material.from_refractiveindex_formula_file("bad", str(path))


# ---------------------------------------------------------------------------
# numpy.linalg.LinAlgError (never caught -- propagates from a deliberately
# singular matrix)
# ---------------------------------------------------------------------------


def test_patterned_eigensolver_propagates_linalg_error_for_singular_epsilon_hat():
    omega = 2 * np.pi
    kx = np.array([0.5, 1.5])
    ky = np.array([0.3, -0.2])
    singular = np.array([[1.0, 1.0], [1.0, 1.0]], dtype=complex)  # rank 1, exactly singular
    with pytest.raises(np.linalg.LinAlgError):
        solve_layer_eigenmodes_patterned(omega, kx, ky, singular)


def test_patterned_inplane_eigensolver_propagates_linalg_error_for_singular_ezz():
    omega = 2 * np.pi
    kx = np.array([0.5, 1.5])
    ky = np.array([0.3, -0.2])
    good = np.eye(2, dtype=complex)
    singular = np.array([[1.0, 1.0], [1.0, 1.0]], dtype=complex)
    with pytest.raises(np.linalg.LinAlgError):
        solve_layer_eigenmodes_patterned_inplane(omega, kx, ky, good, good * 0, good * 0, good, singular)


def test_1d_eigensolver_propagates_linalg_error_for_singular_epsilon_inv_hat():
    omega = 2 * np.pi
    kx = np.array([0.5, 1.5])
    ky = np.array([0.0, 0.0])
    good = np.eye(2, dtype=complex)
    singular = np.array([[1.0, 1.0], [1.0, 1.0]], dtype=complex)
    with pytest.raises(np.linalg.LinAlgError):
        solve_layer_eigenmodes_1d(omega, kx, ky, good, singular)
