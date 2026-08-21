"""Phase 4a example: circular air via array in a Si substrate on a square lattice.

Stack (top to bottom):
    air              (incidence, semi-infinite)
    Si substrate with air-filled circular via (finite thickness, 2D-periodic)
    Si              (exit, semi-infinite)

The geometry mirrors a through-silicon-via (TSV) scatterometry target:
Si bulk with air holes etched through it, as in the vendored
`EMTutorial/Scatterometry/ThroughSiliconVia` JCMsuite reference case.

Si is real dispersive `n,k` data (`NK_FILE/si_KLA.txt`, via
`Material.from_nk_file`) rather than a flat constant, and the wavelength
range is 0.4-0.8um/401 points -- matching `pillar_array.py`'s own updated
convention (decisions.md ADR-039) and `structures/trench/tapered_trench.py`'s
range, per the project owner's direct request to stop using a simplified
constant index. The corresponding Lumerical build needs the *same*
imported `si_KLA.txt` data as a custom material (not the built-in
"Si (Silicon) - Palik" dataset, a materials-model mismatch already found
and fixed once, ADR-036).

Run with:  python structures/via/via_array.py
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Circle, Lattice, Pattern
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.output_paths import run_output_dir, write_run_metadata
from sougata_solver.simulation import Simulation
from sougata_solver.sweep import avoid_rayleigh_wood_anomalies

# Si_KLA.txt has columns: wavelength [nm], n, k -- same file/loader
# structures/trench/tapered_trench.py and pillar_array.py already use.
NK_DIR = Path(__file__).resolve().parents[3] / "NK_FILE"
SI_NK_PATH = NK_DIR / "si_KLA.txt"

# ============================================================================
# EDIT: via geometry and materials
# ============================================================================
PERIOD = 0.7e-6          # lattice period (meters)
VIA_RADIUS = 0.18e-6     # via (hole) radius (meters)
THICKNESS = 0.46e-6      # via depth (meters)
N_VIA = 1.0              # via fill index (air)

# ============================================================================
# EDIT: incident light -- angle, polarization, order truncation
# ============================================================================
INCIDENT_ANGLE_DEG = 0.0
# Re-measured after switching to real dispersive Si and the new 0.4-0.8um
# range (see decisions.md ADR-039's addendum for the exact convergence scan
# -- do not assume the old constant-index/0.5-1.5um value still applies).
NUM_ORDERS = 81            # 2D Fourier order truncation parameter -- PLACEHOLDER, being re-measured
S_AMPLITUDE = 1.0         # 1.0/0.0 = s-pol; 0.0/1.0 = p-pol
P_AMPLITUDE = 0.0

# ============================================================================
# EDIT: wavelength sweep (meters)
# ============================================================================
# 0.4-0.8um/401 points -- matches pillar_array.py's updated range and
# structures/trench/tapered_trench.py's own convention. Any grid point
# landing exactly on a Rayleigh/Wood's-anomaly wavelength for this
# PERIOD/NUM_ORDERS/angle (troubleshooting.md's documented q==0 divide-by-
# zero) is nudged automatically -- no manual per-range recomputation needed.
WAVELENGTHS = avoid_rayleigh_wood_anomalies(
    np.linspace(0.40e-6, 0.80e-6, 401), period=PERIOD, num_orders=NUM_ORDERS, theta=math.radians(INCIDENT_ANGLE_DEG)
)

OUTPUT_CSV = "output_via_RT.csv"


def build_geometry(period=None, via_radius=None, thickness=None):
    """Returns (layers, lattice, incidence, transmission)."""
    period = period if period is not None else PERIOD
    via_radius = via_radius if via_radius is not None else VIA_RADIUS
    thickness = thickness if thickness is not None else THICKNESS

    air = Material("air", N_VIA**2)
    # Real dispersive Si, not a flat n=3.48 constant -- same NK_FILE/loader
    # tapered_trench.py and pillar_array.py already use.
    substrate = Material.from_nk_file("Si", str(SI_NK_PATH), "nm")

    # Via (air hole) centered in the unit cell at (period/2, period/2)
    pattern = Pattern(
        background=substrate,
        shapes=[Circle(center=(period / 2, period / 2), radius=via_radius, material=air)],
    )
    lattice = Lattice(a=(period, 0.0), b=(0.0, period))
    layers = [Layer("via_layer", thickness, pattern=pattern)]
    return layers, lattice, air, substrate


def main() -> None:
    layers, lattice, air, substrate = build_geometry()
    sim = Simulation(lattice, layers, num_orders=NUM_ORDERS, incidence=air, transmission=substrate)

    reflectance = np.zeros(len(WAVELENGTHS))
    transmittance = np.zeros(len(WAVELENGTHS))

    print(f"{'wavelength (nm)':>16}  {'R':>8}  {'T':>8}  {'R+T':>8}")
    for i, wavelength in enumerate(WAVELENGTHS):
        excitation = PlaneWaveExcitation(
            wavelength=wavelength,
            theta=math.radians(INCIDENT_ANGLE_DEG),
            phi=0.0,
            s_amplitude=S_AMPLITUDE,
            p_amplitude=P_AMPLITUDE,
        )
        result = sim.solve(excitation)
        reflectance[i] = result.reflectance()
        transmittance[i] = result.transmittance()
        rt = reflectance[i] + transmittance[i]
        print(f"{wavelength * 1e9:16.1f}  {reflectance[i]:8.4f}  {transmittance[i]:8.4f}  {rt:8.4f}")

    out = run_output_dir("via_array")
    csv_path = out / OUTPUT_CSV
    table = np.column_stack([WAVELENGTHS * 1e9, reflectance, transmittance])
    np.savetxt(csv_path, table, delimiter=",", header="wavelength_nm,R,T", comments="")
    write_run_metadata(
        out,
        __file__,
        period_m=PERIOD,
        via_radius_m=VIA_RADIUS,
        thickness_m=THICKNESS,
        substrate_material="Si (dispersive, NK_FILE/si_KLA.txt)",
        n_via=N_VIA,
        incident_angle_deg=INCIDENT_ANGLE_DEG,
        num_orders=NUM_ORDERS,
        s_amplitude=S_AMPLITUDE,
        p_amplitude=P_AMPLITUDE,
        wavelength_range_m=(WAVELENGTHS[0], WAVELENGTHS[-1], len(WAVELENGTHS)),
    )
    print(f"\nSaved {len(WAVELENGTHS)} rows to {csv_path}")
    print(f"Run metadata: {out / 'run_metadata.txt'}")


if __name__ == "__main__":
    main()
