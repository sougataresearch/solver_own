"""Category 8 (`COMMERCIAL_RCWA_ATOMIC_TARGETS.md`): library-level
parameter sweeps built on top of `Simulation.solve()`.

Every `structures/*.py` wavelength-sweep script (e.g.
`structures/trench/trench_grating.py`) independently hand-wrote the same
"for value in values: build an excitation/layer edit; call `solve()`;
collect R/T" loop. This module promotes that pattern into a small set of
typed, reusable functions, per the category's own exit criterion: **each
sweep here is equivalent to repeated scalar `Simulation.solve()` calls and
introduces no new solver-formula risk** -- every function in this module
calls only already-validated `Simulation`/`SimulationResult` methods and
does no numerics of its own beyond bookkeeping (the one exception,
`find_convergence_index`, is a data-selection rule, not a physics formula,
and is validated against already-recorded convergence behavior --
see its own docstring).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Sequence

import numpy as np

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.simulation import Simulation, SimulationResult


@dataclass
class SweepResult:
    """Category 8 target 8.1: typed container for a one-parameter sweep.

    `parameter_values` holds whatever the swept parameter's natural
    representation is -- a `float` for wavelength/angle/thickness sweeps
    (unit given by `parameter_unit`, e.g. `"m"`/`"rad"`), or an opaque
    per-point label for a discrete sweep (target 8.4's polarization sweep
    uses `(s_amplitude, p_amplitude)` Jones-state tuples,
    `parameter_unit=""`). `metadata` is a plain dict of the *fixed*
    parameters that produced this sweep (num_orders, the excitation
    parameters held constant, ...) -- intended to be passed straight into
    `output_paths.write_run_metadata(output_dir, script_path,
    **sweep_result.metadata)` by a `structures/*.py` caller (this module
    itself never touches disk, per ADR-009/010's `structures/`-vs-library
    boundary). `extra` holds sweep-specific derived data that doesn't fit
    every sweep type (e.g. target 8.6's per-point conservation residual).
    """

    parameter_name: str
    parameter_unit: str
    parameter_values: list
    results: list[SimulationResult]
    metadata: dict = field(default_factory=dict)
    extra: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if len(self.parameter_values) != len(self.results):
            raise ValueError(
                f"parameter_values ({len(self.parameter_values)}) and results "
                f"({len(self.results)}) must have the same length"
            )
        if len(self.parameter_values) == 0:
            raise ValueError("a sweep must have at least one point")

    def reflectance(self) -> np.ndarray:
        return np.array([r.reflectance() for r in self.results])

    def transmittance(self) -> np.ndarray:
        return np.array([r.transmittance() for r in self.results])


def sweep_wavelength(
    simulation: Simulation,
    wavelengths: Sequence[float],
    theta: float,
    phi: float,
    s_amplitude: complex,
    p_amplitude: complex,
) -> SweepResult:
    """Category 8 target 8.2: promotes the loop every `structures/*.py`
    wavelength-sweep script hand-wrote (e.g.
    `structures/trench/trench_grating.py:71-82`) into a library function.
    """
    results = [
        simulation.solve(PlaneWaveExcitation(wavelength, theta, phi, s_amplitude, p_amplitude))
        for wavelength in wavelengths
    ]
    metadata = {
        "num_orders": simulation.num_orders,
        "truncation": simulation.truncation,
        "theta_rad": theta,
        "phi_rad": phi,
        "s_amplitude": s_amplitude,
        "p_amplitude": p_amplitude,
    }
    return SweepResult("wavelength", "m", list(wavelengths), results, metadata)


def sweep_theta(
    simulation: Simulation,
    wavelength: float,
    thetas: Sequence[float],
    phi: float,
    s_amplitude: complex,
    p_amplitude: complex,
) -> SweepResult:
    """Category 8 target 8.3: polar-angle sweep at fixed wavelength/phi.

    `toeplitz_matrix`/`toeplitz_matrix_component` depend only on
    `(pattern, wavelength)`, not on incidence angle (`design.md`'s
    "Layer/Toeplitz Caching Design", `decisions.md` ADR-016) -- reusing
    one `Simulation` instance across an entire theta sweep is exactly the
    scenario that cache was measured against, so this sweep benefits from
    it automatically with no extra code here.
    """
    results = [
        simulation.solve(PlaneWaveExcitation(wavelength, theta, phi, s_amplitude, p_amplitude))
        for theta in thetas
    ]
    metadata = {
        "num_orders": simulation.num_orders,
        "truncation": simulation.truncation,
        "wavelength_m": wavelength,
        "phi_rad": phi,
        "s_amplitude": s_amplitude,
        "p_amplitude": p_amplitude,
    }
    return SweepResult("theta", "rad", list(thetas), results, metadata)


def sweep_phi(
    simulation: Simulation,
    wavelength: float,
    theta: float,
    phis: Sequence[float],
    s_amplitude: complex,
    p_amplitude: complex,
) -> SweepResult:
    """Category 8 target 8.3: azimuthal-angle sweep at fixed wavelength/theta
    -- same fixed-wavelength Toeplitz-cache reuse as `sweep_theta`."""
    results = [
        simulation.solve(PlaneWaveExcitation(wavelength, theta, phi, s_amplitude, p_amplitude))
        for phi in phis
    ]
    metadata = {
        "num_orders": simulation.num_orders,
        "truncation": simulation.truncation,
        "wavelength_m": wavelength,
        "theta_rad": theta,
        "s_amplitude": s_amplitude,
        "p_amplitude": p_amplitude,
    }
    return SweepResult("phi", "rad", list(phis), results, metadata)


def sweep_polarization(
    simulation: Simulation,
    wavelength: float,
    theta: float,
    phi: float,
    jones_states: Sequence[tuple[complex, complex]],
) -> SweepResult:
    """Category 8 target 8.4: sweep over a finite, explicit list of Jones
    `(s_amplitude, p_amplitude)` states (e.g. TE/TM/linear/RCP/LCP,
    `CONVENTIONS.md`'s "Worked polarization examples" table) at fixed
    wavelength/angle -- deliberately a finite list, not a continuous
    parametrization, per the target's own wording."""
    if len(jones_states) == 0:
        raise ValueError("jones_states must be a non-empty finite list")
    results = [
        simulation.solve(PlaneWaveExcitation(wavelength, theta, phi, s, p)) for s, p in jones_states
    ]
    metadata = {
        "num_orders": simulation.num_orders,
        "truncation": simulation.truncation,
        "wavelength_m": wavelength,
        "theta_rad": theta,
        "phi_rad": phi,
    }
    return SweepResult("jones_state", "(s_amplitude, p_amplitude)", list(jones_states), results, metadata)


def sweep_thickness(
    simulation: Simulation,
    layer_name: str,
    thicknesses: Sequence[float],
    excitation: PlaneWaveExcitation,
) -> SweepResult:
    """Category 8 target 8.5: sweep one named layer's thickness, holding
    everything else in the stack and the excitation fixed.

    `thicknesses` are validated explicitly here (finite, `> 0`), not via
    `Layer.__post_init__` (Category 7 target 7.1's
    `layer._require_valid_thickness`) -- that validation only fires at
    construction time, and this sweep must mutate an already-constructed
    `Layer`'s `thickness` attribute in place between `solve()` calls
    (`Simulation.solve()` reads each layer's `thickness` fresh from
    `self.layer_stack` every call, so mutating it here is how
    `Simulation`'s otherwise-fixed-at-`__init__` layer stack can still be
    swept). The layer's original thickness is restored afterward (even on
    error), so this function leaves no surprising side effect on the
    `Simulation` instance for the caller.
    """
    target_layer = None
    for layer in simulation.layer_stack:
        if layer.name == layer_name:
            target_layer = layer
            break
    if target_layer is None:
        names = [layer.name for layer in simulation.layer_stack]
        raise ValueError(f"no layer named {layer_name!r} in this stack (have: {names})")

    for t in thicknesses:
        if not math.isfinite(t) or not (t > 0):
            raise ValueError(f"thickness {t!r} for layer {layer_name!r} must be finite and > 0")

    original_thickness = target_layer.thickness
    results = []
    try:
        for t in thicknesses:
            target_layer.thickness = t
            results.append(simulation.solve(excitation))
    finally:
        target_layer.thickness = original_thickness

    metadata = {
        "layer_name": layer_name,
        "num_orders": simulation.num_orders,
        "wavelength_m": excitation.wavelength,
        "theta_rad": excitation.theta,
        "phi_rad": excitation.phi,
    }
    return SweepResult(f"thickness[{layer_name}]", "m", list(thicknesses), results, metadata)


def harmonic_study(
    build_simulation: Callable[[int], Simulation],
    num_orders_values: Sequence[int],
    excitation: PlaneWaveExcitation,
) -> SweepResult:
    """Category 8 target 8.6: R/T and the `R+T+sum(layer_absorption())`
    conservation residual (Category 7 target 7.6) versus harmonic-order
    count. Makes **no** automatic stopping decision (that is target 8.8,
    gated on target 8.7's criterion) -- purely a data-collection pass.

    Takes a `Simulation`-*builder* callable, not a single `Simulation`
    instance, because `num_orders` is a construction-time parameter:
    `design.md`'s "Layer/Toeplitz Caching Design" (Category 7 target 7.3)
    explicitly documents and relies on the reciprocal-lattice truncation
    `g` being fixed for a given `Simulation` instance's entire lifetime --
    there is no supported way to change `num_orders` on one live instance,
    so each point in this study gets its own freshly-built `Simulation`.
    """
    if len(num_orders_values) == 0:
        raise ValueError("num_orders_values must be non-empty")
    results = []
    residuals = []
    for n in num_orders_values:
        sim = build_simulation(n)
        result = sim.solve(excitation)
        results.append(result)
        r, t = result.reflectance(), result.transmittance()
        absorbed = sum(result.layer_absorption())
        residuals.append(abs(1.0 - (r + t + absorbed)))
    metadata = {"wavelength_m": excitation.wavelength, "theta_rad": excitation.theta, "phi_rad": excitation.phi}
    return SweepResult(
        "num_orders", "", list(num_orders_values), results, metadata,
        extra={"conservation_residual": np.array(residuals)},
    )


def rayleigh_wood_wavelengths(
    period: float, num_orders: int, theta: float = 0.0, n_incidence: float = 1.0, n_check: float = 1.0
) -> list[float]:
    """Wavelengths where a diffraction order's out-of-plane wavenumber `q`
    is exactly zero (the Rayleigh/Wood's-anomaly condition) for a 1D
    lattice of the given `period`, within the `[-num_orders, num_orders]`
    truncated order set -- the exact singular points where
    `smatrix.py::interface_smatrix`'s `kp @ phi / q` divides by zero
    (`troubleshooting.md`'s documented, pre-existing Rayleigh-threshold
    limitation).

    Derived directly from this solver's own formulas, not a separately
    assumed grating equation: `kx0 = n_incidence*omega*sin(theta)`
    (`excitation.py::k_parallel`), `kx_m = kx0 + 2*pi*m/period`
    (`simulation.py:498`, for a `Lattice1D` whose reciprocal vector is
    `1/period`), and `q_sq = eps*omega**2 - kx**2` (`eigenmodes.py:332`,
    `ky=0` case). Setting `q_sq=0` for a lossless medium of index
    `n_check` (`eps = n_check**2`) and solving for wavelength gives::

        wavelength = period * (s*n_check - n_incidence*sin(theta)) / m

    for each nonzero integer order `m` in `[-num_orders, num_orders]` and
    `s` in `{+1, -1}`. `n_check` defaults to air (`1.0`) since that's the
    lossless medium (incidence half-space, or an air-filled groove) this
    singularity actually bites in practice -- an absorbing medium's `q_sq`
    is generically complex and does not hit exactly zero at a real
    wavelength, so it never needs checking here.

    Confirmed against an actual observed failure (a `PERIOD=2.032e-6`
    grating's divide-by-zero at exactly `508e-9`, `structures/trench/
    trench_ocd_sweep.py`): `rayleigh_wood_wavelengths(2.032e-6, 15)`
    includes `508e-9` (order `m=4`), matching what was found by hand.
    """
    if not (period > 0):
        raise ValueError(f"period must be > 0, got {period!r}")
    if num_orders < 1:
        raise ValueError(f"num_orders must be >= 1, got {num_orders!r}")

    sin_theta = math.sin(theta)
    wavelengths = []
    for m in list(range(-num_orders, 0)) + list(range(1, num_orders + 1)):
        for s in (1.0, -1.0):
            denom = m
            wl = period * (s * n_check - n_incidence * sin_theta) / denom
            if wl > 0:
                wavelengths.append(wl)
    return sorted(set(wavelengths))


def avoid_rayleigh_wood_anomalies(
    wavelengths: np.ndarray,
    period: float,
    num_orders: int,
    theta: float = 0.0,
    n_incidence: float = 1.0,
    n_check: float = 1.0,
    nudge: float = 0.5e-9,
) -> np.ndarray:
    """Return a copy of `wavelengths` with any grid point that lands
    (to float precision) on a `rayleigh_wood_wavelengths(...)` singularity
    nudged by `nudge` (meters) -- so a `structures/*.py` script never has
    to hand-derive/hardcode the singular wavelength for its own
    `period`/`num_orders`/`theta` the way `troubleshooting.md`'s original
    entries did. Only the offending point(s) move; every other point in
    `wavelengths`, including its endpoints, is returned unchanged.

    `nudge` defaults to 0.5nm, confirmed by direct probing (see the
    conversation this was added in) to already be comfortably clear of the
    singularity's divide-by-zero (0.1nm away was already clean).
    """
    wavelengths = np.array(wavelengths, dtype=float, copy=True)
    singular = rayleigh_wood_wavelengths(period, num_orders, theta, n_incidence, n_check)
    if len(wavelengths) < 2:
        return wavelengths
    step = np.median(np.abs(np.diff(wavelengths)))
    tol = max(step * 1e-6, 1e-15)
    for wl_singular in singular:
        idx = np.flatnonzero(np.abs(wavelengths - wl_singular) < tol)
        wavelengths[idx] += nudge
    return wavelengths


def find_convergence_index(values: Sequence[float], tolerance: float) -> int | None:
    """Category 8 target 8.7: a conservative convergence criterion.

    Returns the smallest index `i` such that **every later value** in
    `values` (not just the immediate next one) stays within `tolerance` of
    `values[i]` -- deliberately conservative, per this project's own
    already-recorded finding (`tests/test_fourier_convergence.py`,
    Category 3 targets 3.2/3.3) that a high-contrast pattern's
    convergence-vs-`num_orders` curve can be non-monotonic at low
    harmonic-order counts (a transient dip/spike a naive "within tolerance
    of the immediate next point" criterion would wrongly report as
    converged).

    An index requires confirmation from **at least one later point** to
    count as converged -- the last index in `values` is deliberately never
    itself eligible (it would otherwise trivially "converge" against zero
    remaining points to disagree with, a vacuous truth, not a genuine
    convergence confirmation). Consequently, `None` is returned both when
    fewer than two points are given and when the sequence never actually
    settles within the given data -- never claims convergence beyond what
    was actually measured.
    """
    values_arr = np.asarray(values, dtype=float)
    n = len(values_arr)
    for i in range(n - 1):
        if np.all(np.abs(values_arr[i:] - values_arr[i]) <= tolerance):
            return i
    return None


def auto_select_num_orders(
    build_simulation: Callable[[int], Simulation],
    num_orders_candidates: Sequence[int],
    excitation: PlaneWaveExcitation,
    tolerance: float,
    metric: str = "reflectance",
) -> tuple[int, SweepResult]:
    """Category 8 target 8.8: implemented only after target 8.7's
    criterion was validated to succeed on thin-film, trench, and pillar
    fixtures (`tests/test_harmonic_convergence.py`), per the target's own
    "implement only after 8.7 succeeds" wording.

    Runs `harmonic_study` over `num_orders_candidates` (must be given in
    ascending order), applies `find_convergence_index` to the chosen
    `metric` (`"reflectance"` or `"transmittance"`), and returns
    `(num_orders_candidates[index], sweep_result)`. Raises `ValueError` if
    no candidate converges within the given list -- never silently guesses
    or falls back to the largest candidate.
    """
    candidates = list(num_orders_candidates)
    if candidates != sorted(candidates):
        raise ValueError("num_orders_candidates must be given in ascending order")
    if metric not in ("reflectance", "transmittance"):
        raise ValueError(f"metric must be 'reflectance' or 'transmittance', got {metric!r}")

    sweep = harmonic_study(build_simulation, candidates, excitation)
    values = sweep.reflectance() if metric == "reflectance" else sweep.transmittance()
    index = find_convergence_index(values, tolerance)
    if index is None:
        raise ValueError(
            f"no num_orders in {candidates} converged to tolerance={tolerance} "
            f"for metric={metric!r} -- try more/larger candidates"
        )
    return candidates[index], sweep
