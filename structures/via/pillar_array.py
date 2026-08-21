"""Phase 4a example: circular Si pillar array on a square lattice (air
background), grown on a Si substrate.

Stack (top to bottom):
    air              (incidence, semi-infinite)
    Si circular pillar / air background (finite thickness, 2D-periodic)
    Si substrate     (exit, semi-infinite)

Substrate reuses the pillar's own material (same `Material` instance, not
an independently chosen one) -- a pillar grown from/etched into the same
wafer material it stands on, the common case absent a reason to pick a
different substrate. Requested directly by the project owner after
reviewing a free-standing (air/pillar/air) render and pointing out that a
real pillar device is essentially never fully free-standing -- see
decisions.md's ADR-039 addendum for the full account, including the
re-run convergence check this change required (the substrate changes the
sharp resonance behavior the free-standing case showed).

Si is real dispersive `n,k` data (`NK_FILE/si_KLA.txt`, via
`Material.from_nk_file`) rather than a flat constant, and the wavelength
range is 0.4-0.8um/401 points -- matching `structures/trench/tapered_trench.py`'s
own convention exactly, per the project owner's direct request to stop
using a simplified constant index. The corresponding Lumerical build
needs the *same* imported `si_KLA.txt` data as a custom material (not the
built-in "Si (Silicon) - Palik" dataset, a materials-model mismatch
already found and fixed once, ADR-036), the same way the project owner
already did for the trench cross-validation.

Run with:  python structures/via/pillar_array.py
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
# structures/trench/tapered_trench.py already uses.
NK_DIR = Path(__file__).resolve().parents[3] / "NK_FILE"
SI_NK_PATH = NK_DIR / "si_KLA.txt"

# ============================================================================
# EDIT: pillar geometry and materials
# ============================================================================
PERIOD = 0.7e-6           # lattice period (meters)
PILLAR_RADIUS = 0.18e-6  # pillar radius (meters)
THICKNESS = 0.46e-6      # pillar height (meters)
N_BG = 1.0                # background index (air), within the patterned layer only

# ============================================================================
# EDIT: incident light -- angle, polarization, order truncation
# ============================================================================
INCIDENT_ANGLE_DEG = 0.0
# Re-measured after switching to real dispersive Si (Si's index varies
# substantially over 0.4-0.8um -- higher and lossier at the short-wavelength
# end than the flat n=3.48 this file used before) and the new 0.4-0.8um
# range: see decisions.md's ADR-039 addendum for the exact convergence scan.
NUM_ORDERS = 121          # 2D Fourier order truncation parameter
S_AMPLITUDE = 1.0         # 1.0/0.0 = s-pol; 0.0/1.0 = p-pol
P_AMPLITUDE = 0.0

# ============================================================================
# EDIT: wavelength sweep (meters)
# ============================================================================
# 0.4-0.8um/401 points -- matches structures/trench/tapered_trench.py's own
# range exactly, per the project owner's direct request (real dispersive Si
# needs a wavelength range where that dispersion is actually being probed
# meaningfully, not an arbitrary choice). Any grid point landing exactly on
# a Rayleigh/Wood's-anomaly wavelength for this PERIOD/NUM_ORDERS/angle
# (troubleshooting.md's documented q==0 divide-by-zero) is nudged
# automatically -- no manual per-range recomputation needed.
WAVELENGTHS = avoid_rayleigh_wood_anomalies(
    np.linspace(0.40e-6, 0.80e-6, 401), period=PERIOD, num_orders=NUM_ORDERS, theta=math.radians(INCIDENT_ANGLE_DEG)
)

OUTPUT_CSV = "output_pillar_RT.csv"


def build_geometry(period=None, pillar_radius=None, thickness=None):
    """Returns (layers, lattice, incidence, transmission)."""
    period = period if period is not None else PERIOD
    pillar_radius = pillar_radius if pillar_radius is not None else PILLAR_RADIUS
    thickness = thickness if thickness is not None else THICKNESS

    air = Material("air", N_BG**2)
    # Real dispersive Si, not a flat n=3.48 constant -- same NK_FILE/loader
    # tapered_trench.py already uses. Pillar and substrate share the exact
    # same Material instance (grown from/etched into the same wafer).
    pillar = Material.from_nk_file("Si", str(SI_NK_PATH), "nm")
    substrate = pillar

    # Pillar centered in the unit cell at (period/2, period/2)
    pattern = Pattern(
        background=air,
        shapes=[Circle(center=(period / 2, period / 2), radius=pillar_radius, material=pillar)],
    )
    lattice = Lattice(a=(period, 0.0), b=(0.0, period))
    layers = [Layer("pillar_layer", thickness, pattern=pattern)]
    return layers, lattice, air, substrate


def main() -> None:
    layers, lattice, air, transmission = build_geometry()
    sim = Simulation(lattice, layers, num_orders=NUM_ORDERS, incidence=air, transmission=transmission)

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

    out = run_output_dir("pillar_array")
    csv_path = out / OUTPUT_CSV
    table = np.column_stack([WAVELENGTHS * 1e9, reflectance, transmittance])
    np.savetxt(csv_path, table, delimiter=",", header="wavelength_nm,R,T", comments="")
    write_run_metadata(
        out,
        __file__,
        period_m=PERIOD,
        pillar_radius_m=PILLAR_RADIUS,
        thickness_m=THICKNESS,
        pillar_substrate_material="Si (dispersive, NK_FILE/si_KLA.txt)",
        n_bg=N_BG,
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
