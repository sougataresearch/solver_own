"""Single-taper ridge grating, staircase-discretized.

Matches the structure built by the Lumerical script actually in active use
(the "Edit structure group" script tab, keyed to the `tcd/bcd/depth/spacing/
grating_number/zSpan/yCompensation/zCompensation` user properties shown in
that dialog) -- a single trapezoid per period, NOT the two-segment
`ttcd`/`DoT` "5CD" version from `Code2_grating generate_5CD_hui.py` (that
script's own embedded structure-group script has an extra "upper" polygon;
this one doesn't, confirmed directly against the property panel, which
lists no `ttcd`/`DoT` property at all).

Stack (top to bottom, i.e. incidence-side to substrate-side):
    air                                (incidence, semi-infinite)
    ridge: tapers tcd (top) -> bcd     (single staircase-discretized
      (bottom), over thickness depth    segment, 1D-periodic along x)
    substrate                          (transmission, semi-infinite)

RIDGE_MATERIAL / SUBSTRATE_MATERIAL / GROOVE_MATERIAL are each independently
selectable below (not hardcoded) -- per the property panel, the ridge here
is SiO2 (Glass), but substrate and groove-fill aren't shown in that dialog
and must be set to whatever your actual Lumerical simulation actually uses.
Use whichever KLA n,k file matches the *same* database entry Lumerical's
material picker resolved to (e.g. Palik vs. KLA vs. Lumerical's own SiO2
model are different datasets) -- confirm this mapping manually; nothing in
this script can infer it for you.

`grating_number` (finite repeat count in the Lumerical FDTD-style geometry
view) has no equivalent parameter here: this solver's RCWA formalism is
already exactly, infinitely periodic via `Lattice1D`, so there's nothing to
truncate -- it's not silently dropped, it simply doesn't apply.

NOT modeled here (same caveats as the single-segment predecessor of this
script): no corner rounding (R_top/R_bot), no bowed (non-linear) sidewall.

Sweeps `num_slices` at a fixed wavelength/angle and prints R/T/A.

Run with:  python structures/trench/tapered_trench.py
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice1D
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
# EDIT: select ridge / substrate / groove-fill materials independently.  Use
# one of MATERIAL_NK_PATHS' keys, or "air" for vacuum/air (constant n=1).
# Confirm each against your actual Lumerical simulation before comparing --
# none of this can be inferred from the structure-group script alone.
# ============================================================================
RIDGE_MATERIAL = "Si"        # matches the property panel's "material" = SiO2 (Glass)
SUBSTRATE_MATERIAL = "Si"     # NOT shown in the property panel -- set to match your actual Lumerical substrate
GROOVE_MATERIAL = "air"       # gap fill between ridges -- NOT shown in the property panel either


# ============================================================================
# EDIT: the four CDs, matching the property panel's values exactly.
# ============================================================================
TCD = 1.383e-6      # critical dimension at the ridge top (incidence side), meters
BCD = 1.322e-6      # critical dimension at the ridge base (substrate side), meters
DEPTH = 4.981e-6          # ridge height (meters)
SPACING = 0.649e-6  # gap between adjacent ridges (meters)
PERIOD = TCD + SPACING     # = 3.0 um for the panel's values -- sanity-check this against your own numbers


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
    np.linspace(0.20e-6, 0.40e-6, 100), period=PERIOD, num_orders=NUM_ORD, theta=math.radians(INCIDENT_ANGLE_DEG)
)

# ============================================================================
# EDIT: slice-count sweep
# ============================================================================
SLICE_COUNTS = [1, 2, 4, 8, 16, 32, 64]

OUTPUT_CSV = "output_trench_RT.csv"


def _material(name: str) -> Material:
    if name == "air":
        return Material("air", 1.0)
    if name not in MATERIAL_NK_PATHS:
        raise ValueError(f"material must be 'air' or one of {tuple(MATERIAL_NK_PATHS)}, got {name!r}")
    return Material.from_nk_file(name, str(MATERIAL_NK_PATHS[name]), NK_WAVELENGTH_UNIT)


def build_geometry(num_slices=None, period=None, tcd=None, bcd=None, depth=None):
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

    ridge = _material(RIDGE_MATERIAL)
    groove = _material(GROOVE_MATERIAL)
    substrate = _material(SUBSTRATE_MATERIAL)
    air = Material("air", 1.0)
    lattice = Lattice1D(period)
    layers = staircase_slab_layers(
        center_x=0.0,
        top_halfwidth=0.5 * tcd,
        bottom_halfwidth=0.5 * bcd,
        thickness=depth,
        num_slices=num_slices,
        shape_material=ridge,
        background_material=groove,
    )
    return layers, lattice, air, substrate


def main() -> None:
    print(f"Period (pitch)     = {PERIOD * 1e6:.4f} um  (= TCD + SPACING)")
    print(f"TCD / BCD          = {TCD * 1e6:.4f} / {BCD * 1e6:.4f} um")
    print(f"Depth              = {DEPTH * 1e6:.4f} um")
    print(f"Ridge / Substrate / Groove = {RIDGE_MATERIAL} / {SUBSTRATE_MATERIAL} / {GROOVE_MATERIAL}")
    print()

    num_orders = 2 * NUM_ORD + 1

    reflectance = np.zeros((len(SLICE_COUNTS), len(WAVELENGTHS)))
    transmittance = np.zeros((len(SLICE_COUNTS), len(WAVELENGTHS)))
    absorptance = np.zeros((len(SLICE_COUNTS), len(WAVELENGTHS)))

    for si_idx, num_slices in enumerate(SLICE_COUNTS):
        layers, lattice, air, substrate = build_geometry(num_slices=num_slices)
        sim = Simulation(lattice, layers, num_orders=num_orders, incidence=air, transmission=substrate)

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
        ridge_material=RIDGE_MATERIAL,
        substrate_material=SUBSTRATE_MATERIAL,
        groove_material=GROOVE_MATERIAL,
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
