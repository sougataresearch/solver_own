"""Depth-tapered trench, staircase-discretized -- built from the project
owner's senior's FDTD reference file (`Trench_Result_0.3.fsp`), confirmed
directly (not assumed) via that file's own Lumerical dialogs:

    - `dimension: 2D` and boundary conditions `x=Periodic, y=PML` mean the
      structure is periodic along x (period 2.028 um) and uniform/invariant
      along whatever axis a 2D simulation doesn't include -- a standard 1D
      lamellar trench, matching `sougata_solver`'s `Lattice1D` +
      `staircase_slab_layers` exactly.
    - Source `injection axis = y-axis`, monitor `type = 2D Y-normal` confirm
      y is depth/propagation *in that FDTD file's own convention* (not
      `sougata_solver`'s, which is always z -- see `CONVENTIONS.md`).
    - `Etch_0`'s vertices (object center y=0.0303209, relative vertices at
      y=+-2.25568), read as y=depth: top edge (y=2.286, exactly flush with
      the `Trench` background rectangle's own y max) has half-width
      0.23009 um (TCD=0.46018 um); bottom edge (y=-2.2253591) has
      half-width 0.243176 um (BCD=0.486352 um) -- narrower at the surface,
      wider at depth (an inverse taper vs. the usual outward-flaring kind;
      built exactly as measured, not "corrected" toward a more common
      convention). Depth = 4.5113591 um.
    - `Trench` (the background rectangle)'s own y min (-3.286) is below the
      etch's bottom edge (-2.2253591) by 1.0606409 um -- a uniform,
      un-etched residual layer left beneath the taper before the
      semi-infinite half-space.

**Found and flagged, not silently reconciled**: a separate, related
Lumerical RCWA attempt at this same structure (`my trench.fsp`) was
confirmed (via its own General tab) to use `propagation axis = z`, and its
own `etch` object was found to be a single z-uniform polygon varying only
*in-plane* along y -- a genuinely different shape from this FDTD
reference's *depth* taper. This script builds the FDTD reference's
structure (the authoritative source, being the senior's file), not that
RCWA attempt's geometry.

RIDGE_MATERIAL / TRANSMISSION_MATERIAL are each independently selectable
below -- confirm these match your actual Lumerical material picks before
comparing results; nothing in this script can infer that mapping. The
residual layer and etch fill are assumed Si / air respectively (matching
the project owner's earlier confirmation for the related structure); the
half-space below the residual layer is assumed air (free-standing) --
neither assumption is visible in the FDTD screenshots shown so far.

`grating_number` (finite repeat count in the Lumerical FDTD-style geometry
view) has no equivalent parameter here: this solver's RCWA formalism is
already exactly, infinitely periodic via `Lattice1D`, so there's nothing to
truncate -- it's not silently dropped, it simply doesn't apply.

NOT modeled here: no corner rounding (R_top/R_bot), no bowed (non-linear)
sidewall.

Sweeps `num_slices` at a fixed wavelength/angle and prints R/T/A.

Run with:  python structures/trench/tapered_trench.py
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice1D
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.output_paths import run_output_dir, write_run_metadata
from sougata_solver.simulation import Simulation
from sougata_solver.staircase import staircase_slab_layers
from sougata_solver.sweep import avoid_rayleigh_wood_anomalies

# KLA material files have columns: wavelength [nm], n, k.
NK_DIR = Path(__file__).resolve().parents[3] / "NK_FILE"
MATERIAL_NK_PATHS = {
    "SiO2": NK_DIR / "sio2_KLA.txt",
    "SiO": NK_DIR / "sio_KLA.txt",
    "Si": NK_DIR / "si_KLA.txt",
    "Ni": NK_DIR / "ni_KLA.txt",
}
NK_WAVELENGTH_UNIT = "nm"

# ============================================================================
# EDIT: select the surrounding-slab / transmission-side materials
# independently. Use one of MATERIAL_NK_PATHS' keys, or "air" for
# vacuum/air (constant n=1). Confirm each against your actual Lumerical
# simulation before comparing -- none of this can be inferred from the
# FDTD dialogs alone.
# ============================================================================
SLAB_MATERIAL = "Si"           # the surrounding solid material (matches "Trench"/"Si_slab")
ETCH_MATERIAL = "air"          # the etched, tapered groove's fill
TRANSMISSION_MATERIAL = "air"  # semi-infinite half-space below the residual layer


# ============================================================================
# EDIT: the four CDs plus the residual layer, matching Trench_Result_0.3.fsp's
# `Etch_0`/`Trench` dialogs exactly (see module docstring for the derivation).
# ============================================================================
TCD = 0.46018e-6           # critical dimension at the surface (top), meters -- narrower
BCD = 0.486352e-6          # critical dimension at depth (bottom), meters -- wider
DEPTH = 4.5113591e-6       # etch depth (meters)
PERIOD = 2.028e-6          # lattice period (meters)
SPACING = PERIOD - TCD     # derived, for reference/sanity-check only
RESIDUAL_THICKNESS = 1.0606409e-6  # uniform, un-etched layer beneath the taper (meters)


# ============================================================================
# EDIT: incident light -- angle, polarization, order truncation.
# ============================================================================
INCIDENT_ANGLE_DEG = 0.0
AZIMUTHAL_ANGLE_DEG = 0.0
NUM_ORD = 15                # orders per side; total Fourier orders = 2*NUM_ORD+1
S_AMPLITUDE = 1.0
P_AMPLITUDE = 0.0

# Any grid point landing exactly on a Rayleigh/Wood's-anomaly wavelength for
# this PERIOD/NUM_ORD/angle (troubleshooting.md's documented q==0 divide-by-
# zero) is nudged automatically -- no manual per-range recomputation needed.
WAVELENGTHS = avoid_rayleigh_wood_anomalies(
    np.linspace(0.40e-6, 0.80e-6, 401), period=PERIOD, num_orders=NUM_ORD, theta=math.radians(INCIDENT_ANGLE_DEG)
)

# ============================================================================
# EDIT: slice-count sweep
# ============================================================================
SLICE_COUNTS = [1, 2, 4, 8, 16, 32, 64]

# Found via a targeted convergence check (4 sample wavelengths across
# SLICE_COUNTS): R stabilizes to ~1e-4 vs. the finest (64-slice) reference
# by 32 slices -- see `decisions.md` ADR-036. Any other script building on
# this same real-device geometry (e.g. `trench_ocd_sweep.py`) should import
# this constant directly rather than hand-copying a slice count, so a future
# re-run of the convergence check here can't silently leave that script
# using a stale value.
RECOMMENDED_NUM_SLICES = 32

OUTPUT_CSV = "output_trench_RT.csv"


def _material(name: str) -> Material:
    if name == "air":
        return Material("air", 1.0)
    if name not in MATERIAL_NK_PATHS:
        raise ValueError(f"material must be 'air' or one of {tuple(MATERIAL_NK_PATHS)}, got {name!r}")
    return Material.from_nk_file(name, str(MATERIAL_NK_PATHS[name]), NK_WAVELENGTH_UNIT)


def build_geometry(num_slices=None, period=None, tcd=None, bcd=None, depth=None, residual_thickness=None):
    """Returns (layers, lattice, incidence, transmission).

    `num_slices` defaults to `SLICE_COUNTS[-1]` (the finest/most
    representative geometry) -- this file's own `main()` sweeps over every
    value in `SLICE_COUNTS` itself, calling this with each one explicitly.
    """
    num_slices = num_slices if num_slices is not None else SLICE_COUNTS[-1]
    period = period if period is not None else PERIOD
    tcd = tcd if tcd is not None else TCD
    bcd = bcd if bcd is not None else BCD
    depth = depth if depth is not None else DEPTH
    residual_thickness = residual_thickness if residual_thickness is not None else RESIDUAL_THICKNESS

    slab = _material(SLAB_MATERIAL)
    etch = _material(ETCH_MATERIAL)
    transmission = _material(TRANSMISSION_MATERIAL)
    air = Material("air", 1.0)
    lattice = Lattice1D(period)
    layers = staircase_slab_layers(
        center_x=0.0,
        top_halfwidth=0.5 * tcd,
        bottom_halfwidth=0.5 * bcd,
        thickness=depth,
        num_slices=num_slices,
        shape_material=etch,
        background_material=slab,
    )
    if residual_thickness > 0:
        layers.append(Layer("residual_slab", residual_thickness, material=slab))
    return layers, lattice, air, transmission


def main() -> None:
    print(f"Period (pitch)     = {PERIOD * 1e6:.4f} um")
    print(f"TCD / BCD          = {TCD * 1e6:.4f} / {BCD * 1e6:.4f} um  (narrower at top, wider at depth)")
    print(f"Depth              = {DEPTH * 1e6:.4f} um  (+ {RESIDUAL_THICKNESS * 1e6:.4f} um residual slab)")
    print(f"Slab / Etch / Transmission = {SLAB_MATERIAL} / {ETCH_MATERIAL} / {TRANSMISSION_MATERIAL}")
    print()

    num_orders = 2 * NUM_ORD + 1

    reflectance = np.zeros((len(SLICE_COUNTS), len(WAVELENGTHS)))
    transmittance = np.zeros((len(SLICE_COUNTS), len(WAVELENGTHS)))
    absorptance = np.zeros((len(SLICE_COUNTS), len(WAVELENGTHS)))

    for si_idx, num_slices in enumerate(SLICE_COUNTS):
        layers, lattice, air, transmission = build_geometry(num_slices=num_slices)
        sim = Simulation(lattice, layers, num_orders=num_orders, incidence=air, transmission=transmission)

        for wl_idx, wavelength in enumerate(WAVELENGTHS):
            excitation = PlaneWaveExcitation(
                wavelength=wavelength,
                theta=math.radians(INCIDENT_ANGLE_DEG),
                phi=math.radians(AZIMUTHAL_ANGLE_DEG),
                s_amplitude=S_AMPLITUDE,
                p_amplitude=P_AMPLITUDE,
            )
            result = sim.solve(excitation)
            reflectance[si_idx, wl_idx] = result.reflectance()
            transmittance[si_idx, wl_idx] = result.transmittance()
            absorptance[si_idx, wl_idx] = 1.0 - reflectance[si_idx, wl_idx] - transmittance[si_idx, wl_idx]

        mid = len(WAVELENGTHS) // 2
        rta = reflectance[si_idx, mid] + transmittance[si_idx, mid] + absorptance[si_idx, mid]
        print(
            f"num_slices={num_slices:3d}  "
            f"R({WAVELENGTHS[mid]*1e9:.0f}nm)={reflectance[si_idx, mid]:.6f}  "
            f"T={transmittance[si_idx, mid]:.6f}  A={absorptance[si_idx, mid]:.6f}  R+T+A={rta:.6f}"
        )

    finest = SLICE_COUNTS[-1]
    out = run_output_dir("tapered_trench")
    csv_path = out / OUTPUT_CSV
    table = np.column_stack([WAVELENGTHS * 1e9, reflectance[-1], transmittance[-1], absorptance[-1]])
    np.savetxt(csv_path, table, delimiter=",", header="wavelength_nm,R,T,A", comments="")
    write_run_metadata(
        out,
        __file__,
        period_m=PERIOD,
        tcd_m=TCD,
        bcd_m=BCD,
        depth_m=DEPTH,
        spacing_m=SPACING,
        residual_thickness_m=RESIDUAL_THICKNESS,
        slab_material=SLAB_MATERIAL,
        etch_material=ETCH_MATERIAL,
        transmission_material=TRANSMISSION_MATERIAL,
        incident_angle_deg=INCIDENT_ANGLE_DEG,
        azimuthal_angle_deg=AZIMUTHAL_ANGLE_DEG,
        num_orders=num_orders,
        s_amplitude=S_AMPLITUDE,
        p_amplitude=P_AMPLITUDE,
        wavelength_range_m=(WAVELENGTHS[0], WAVELENGTHS[-1], len(WAVELENGTHS)),
        finest_num_slices=finest,
    )
    print(f"\nSaved finest-slice-count ({finest}) full spectrum: {csv_path}")
    print(f"Run metadata: {out / 'run_metadata.txt'}")


if __name__ == "__main__":
    main()
