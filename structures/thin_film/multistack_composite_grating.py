"""Laterally-alternating composite of two multilayer stacks (1D-periodic).

EXCEPTION to this folder's usual "no in-plane pattern" rule (see
`structures/README.md`'s `thin_film/` table) -- this script lives here at
the project owner's request even though it uses real `Lattice1D`/`Slab`
patterning (the Phase 3 capability every other file under `trench/` uses).
It is not a groove/ridge grating: every z-slice below alternates between
two *complete* materials across the period, not "material vs. air".

Reproduces a structure built in a commercial RCWA tool (Lumerical RCWA
solver): period 2 um in x, invariant in y. The two sides (left/right) each
carry their own independent stack of (material, thickness) layers -- they
need not match in layer count or in any individual thickness, including
the substrate. Both sides are aligned at one shared reference plane, z=0
(substrate top / bottom of the first film), matching the original
Lumerical model where both substrates' top interface was flush at z=0:

    LEFT_SUBSTRATE_LAYERS / LEFT_FILM_LAYERS    (Si substrate, SiO2 film)
    RIGHT_SUBSTRATE_LAYERS / RIGHT_FILM_LAYERS  (Ni substrate, SiO film)

`build_geometry()` slices the z-axis at every interface height contributed
by *either* side (not just the shared ones), and for a side that has
nothing at a given height, fills in the ambient half-space material on
that side (`TRANSMISSION_MATERIAL` below either side's own stack,
`INCIDENCE_MATERIAL` above it) -- this is how the original example's
SiO2(0.2um)/SiO(0.3um) height mismatch produces a real step, and the same
mechanism now also covers a substrate-thickness mismatch (e.g. Ni deeper
than Si) with no separate code path.

n,k data for Si/Ni/SiO2/SiO comes from the same `NK_FILE/*_KLA.txt` files
`custom_multistack.py` already uses.

Run with:  python structures/thin_film/multistack_composite_grating.py
"""

import math
from pathlib import Path

import numpy as np

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice1D, Pattern, Slab
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.output_paths import run_output_dir, write_run_metadata
from sougata_solver.simulation import Simulation
from sougata_solver.sweep import avoid_rayleigh_wood_anomalies

# KLA material files have columns: wavelength [nm], n, k.
NK_DIR = Path(__file__).resolve().parents[3] / "NK_FILE"
NK_WAVELENGTH_UNIT = "nm"

# ============================================================================
# EDIT (1): materials
# ============================================================================
air = Material("air", 1.0)
si = Material.from_nk_file("Si", str(NK_DIR / "si_KLA.txt"), NK_WAVELENGTH_UNIT)
ni = Material.from_nk_file("Ni", str(NK_DIR / "ni_KLA.txt"), NK_WAVELENGTH_UNIT)
sio2 = Material.from_nk_file("SiO2", str(NK_DIR / "sio2_KLA.txt"), NK_WAVELENGTH_UNIT)
sio = Material.from_nk_file("SiO", str(NK_DIR / "sio_KLA.txt"), NK_WAVELENGTH_UNIT)

INCIDENCE_MATERIAL = air
TRANSMISSION_MATERIAL = air

# ============================================================================
# EDIT (2): grating period and the width of the "left" region within it;
# the remainder of the period is the "right" region.
# ============================================================================
PERIOD = 2.0e-6        # grating period along x (meters)
LEFT_WIDTH = 1.0e-6    # width of the left region within one period (meters)

# ============================================================================
# EDIT (3): each side's own independent stack. Lists are (name, material,
# thickness) tuples. *_SUBSTRATE_LAYERS go deepest-first (bottom to top,
# ending at the shared z=0 reference plane); *_FILM_LAYERS go
# nearest-substrate-first (bottom to top, starting at z=0). Add/remove
# tuples or change any thickness freely, independently per side -- nothing
# else in this file needs to change.
# ============================================================================
LEFT_SUBSTRATE_LAYERS = [("Si", si, 0.5e-6)]
LEFT_FILM_LAYERS = [("SiO2", sio2, 0.2e-6)]

RIGHT_SUBSTRATE_LAYERS = [("Ni", ni, 0.5e-6)]
RIGHT_FILM_LAYERS = [("SiO", sio, 0.3e-6)]

# ============================================================================
# EDIT (4): incident light -- angle, polarization, order truncation
# ============================================================================
INCIDENT_ANGLE_DEG = 0.0
NUM_ORD = 15              # orders per side; total Fourier orders = 2*NUM_ORD+1
S_AMPLITUDE = 1.0         # 1.0/0.0 = pure TE (s); 0.0/1.0 = pure TM (p)
P_AMPLITUDE = 0.0

# ============================================================================
# EDIT (5): wavelength sweep (meters). At this PERIOD/NUM_ORD/angle a few
# grid points land exactly on a Rayleigh/Wood's-anomaly threshold
# (troubleshooting.md's documented q==0 divide-by-zero) -- nudged
# automatically, no manual per-range recomputation needed.
# ============================================================================
WAVELENGTHS = avoid_rayleigh_wood_anomalies(
    np.linspace(0.4e-6, 0.8e-6, 401), period=PERIOD, num_orders=2 * NUM_ORD + 1, theta=math.radians(INCIDENT_ANGLE_DEG)
)

# ============================================================================
# EDIT (6): where to save results (set OUTPUT_CSV_PATH to None to skip)
# ============================================================================
RUN_NAME = "multistack_composite_grating"
OUTPUT_CSV_PATH = "output_multistack_composite_grating_RT.csv"


def _left_region_slab(material: Material) -> Slab:
    """The left region: x in [-PERIOD/2, -PERIOD/2 + LEFT_WIDTH]."""
    center_x = -PERIOD / 2 + LEFT_WIDTH / 2
    return Slab(center_x=center_x, halfwidth=LEFT_WIDTH / 2, material=material)


def _side_intervals(substrate_layers, film_layers):
    """One side's full stack as (z_bottom, z_top, material) tuples, bottom
    to top, z=0 at the substrate/film reference plane both sides share."""
    intervals = []
    z = -sum(thickness for _, _, thickness in substrate_layers)
    for _name, material, thickness in substrate_layers:
        intervals.append((z, z + thickness, material))
        z += thickness
    for _name, material, thickness in film_layers:  # z is 0.0 here
        intervals.append((z, z + thickness, material))
        z += thickness
    return intervals


def _material_at(intervals, z_bottom, z_top, ambient_below, ambient_above):
    """Which material one side has across [z_bottom, z_top]. That interval
    is guaranteed (by construction, see `build_geometry`) to lie entirely
    within one of `intervals` or entirely outside all of them -- in the
    latter case the side is exposed to whichever ambient half-space is on
    that side of its own stack."""
    z_mid = 0.5 * (z_bottom + z_top)
    for lo, hi, material in intervals:
        if lo <= z_mid <= hi:
            return material
    return ambient_below if z_mid < intervals[0][0] else ambient_above


def build_geometry():
    """Returns (layers, lattice, incidence, transmission). Merges the two
    sides' independent stacks into one ordered list of patterned Layers by
    slicing z at every interface height from either side."""
    left_intervals = _side_intervals(LEFT_SUBSTRATE_LAYERS, LEFT_FILM_LAYERS)
    right_intervals = _side_intervals(RIGHT_SUBSTRATE_LAYERS, RIGHT_FILM_LAYERS)

    raw_breakpoints = sorted({z for lo, hi, _ in left_intervals + right_intervals for z in (lo, hi)})
    breakpoints = [raw_breakpoints[0]]
    for z in raw_breakpoints[1:]:
        if z - breakpoints[-1] > 1e-15:  # merge near-duplicates from independent cumulative sums
            breakpoints.append(z)

    layers = []
    for z_bottom, z_top in zip(breakpoints[:-1], breakpoints[1:]):
        left_material = _material_at(left_intervals, z_bottom, z_top, TRANSMISSION_MATERIAL, INCIDENCE_MATERIAL)
        right_material = _material_at(right_intervals, z_bottom, z_top, TRANSMISSION_MATERIAL, INCIDENCE_MATERIAL)
        pattern = Pattern(background=right_material)
        pattern.add(_left_region_slab(left_material))
        name = f"{z_bottom * 1e9:.0f}to{z_top * 1e9:.0f}nm"
        layers.append(Layer(name, z_top - z_bottom, pattern=pattern))

    layers.reverse()  # top to bottom (highest z first), matching this file's Layer-list convention
    lattice = Lattice1D(PERIOD)
    return layers, lattice, INCIDENCE_MATERIAL, TRANSMISSION_MATERIAL


def main():
    layers, lattice, incidence, transmission = build_geometry()
    num_orders = 2 * NUM_ORD + 1
    sim = Simulation(lattice, layers, num_orders=num_orders, incidence=incidence, transmission=transmission)

    print("Layers (top to bottom):")
    for layer in layers:
        print(f"  {layer.name}: {layer.thickness * 1e9:.1f} nm")
    print()

    reflectance = np.zeros(len(WAVELENGTHS))
    transmittance = np.zeros(len(WAVELENGTHS))

    print(f"{'wavelength (nm)':>16}  {'R':>8}  {'T':>8}  {'A':>8}")
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
        absorptance = 1 - reflectance[i] - transmittance[i]
        print(f"{wavelength * 1e9:16.1f}  {reflectance[i]:8.4f}  {transmittance[i]:8.4f}  {absorptance:8.4f}")

    if OUTPUT_CSV_PATH:
        output_dir = run_output_dir(RUN_NAME)
        write_run_metadata(
            output_dir,
            __file__,
            period_m=PERIOD,
            left_width_m=LEFT_WIDTH,
            left_substrate_layers=[(name, t) for name, _mat, t in LEFT_SUBSTRATE_LAYERS],
            left_film_layers=[(name, t) for name, _mat, t in LEFT_FILM_LAYERS],
            right_substrate_layers=[(name, t) for name, _mat, t in RIGHT_SUBSTRATE_LAYERS],
            right_film_layers=[(name, t) for name, _mat, t in RIGHT_FILM_LAYERS],
            merged_layers=[(layer.name, layer.thickness) for layer in layers],
            incidence_material=INCIDENCE_MATERIAL.name,
            transmission_material=TRANSMISSION_MATERIAL.name,
            incident_angle_deg=INCIDENT_ANGLE_DEG,
            num_orders=num_orders,
            s_amplitude=S_AMPLITUDE,
            p_amplitude=P_AMPLITUDE,
            wavelength_range_m=(WAVELENGTHS[0], WAVELENGTHS[-1], len(WAVELENGTHS)),
        )
        absorptance = 1.0 - reflectance - transmittance
        table = np.column_stack([WAVELENGTHS * 1e9, reflectance, transmittance, absorptance])
        output_path = output_dir / OUTPUT_CSV_PATH
        np.savetxt(output_path, table, delimiter=",", header="wavelength_nm,R,T,A", comments="")
        print(f"\nSaved {len(WAVELENGTHS)} rows to {output_path}")
        print(f"Run metadata: {output_dir / 'run_metadata.txt'}")

    return reflectance, transmittance


if __name__ == "__main__":
    main()
