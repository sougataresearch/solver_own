"""Template: an arbitrary-length multilayer thin-film stack, any materials.

Copy this file (within structures/thin_film/) for any new multistack (thin
film / DBR / anti-reflection coating / etc.) instead of editing
sio2_on_si_thin_film.py in place. Everything you're likely to change is in
the numbered EDIT blocks below.

LIMITATION (see phases.md / troubleshooting.md): only *uniform* (unpatterned)
layers work today -- each layer has a thickness (z-direction) only, and is
treated as extending infinitely in x/y. Trench/via/pillar patterning (real
x/y dimensions: line width, radius, pitch) is not implemented yet -- once
Phase 3/4 land, those get their own structures/trench/ and structures/via/
folders and templates.

Run with:  python structures/thin_film/custom_multistack.py
"""

import cmath
import math
from pathlib import Path

import numpy as np

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.output_paths import run_output_dir, write_run_metadata
from sougata_solver.simulation import Simulation
from sougata_solver.sweep import sweep_polarization

# KLA material files have columns: wavelength [nm], n, k.
NK_DIR = Path(__file__).resolve().parents[3] / "NK_FILE"
NK_WAVELENGTH_UNIT = "nm"


# ============================================================================
# EDIT (1): define every material you need, one of three ways
# ============================================================================
air = Material("air", 1.0)                                                          # constant, lossless
sio2 = Material.from_nk_file("SiO2", str(NK_DIR / "sio2_KLA.txt"), NK_WAVELENGTH_UNIT)
sio = Material.from_nk_file("SiO", str(NK_DIR / "sio_KLA.txt"), NK_WAVELENGTH_UNIT)
ni = Material.from_nk_file("Ni", str(NK_DIR / "ni_KLA.txt"), NK_WAVELENGTH_UNIT)
si = Material.from_nk_file("Si", str(NK_DIR / "si_KLA.txt"), NK_WAVELENGTH_UNIT)     # dispersive, from file

# ============================================================================
# EDIT (2): the stack itself -- as many Layer(name, thickness, material=...)
# entries as you want, top to bottom. incidence/transmission below are the
# semi-infinite half-spaces above/below this list (not part of it).
# ============================================================================
INCIDENCE_MATERIAL = air     # what light travels through before hitting the stack
TRANSMISSION_MATERIAL = si   # semi-infinite substrate below the stack
# ^ if you want the substrate to have a *finite* thickness instead of being
#   semi-infinite, add it as a normal Layer(...) in the list below and set
#   TRANSMISSION_MATERIAL back to whatever is truly underneath it (e.g. air).

# Stack: air / SiO2 (200 nm) / SiO (300 nm) / Ni (10 nm) / SiO2 (500 nm) / semi-infinite Si
layers = [
    Layer("SiO2", 200e-9, material=sio2),
    Layer("SiO", 300e-9, material=sio),
    Layer("Ni", 10e-9, material=ni),
    Layer("SiO2", 500e-9, material=sio2),
]

# ============================================================================
# EDIT (3): incident light -- angle (degrees), polarization.
#
# Every polarization state -- TE/TM, any linear angle, RCP/LCP, any ellipse
# -- is one formula, per CONVENTIONS.md's "Worked polarization examples"
# table (Category 6 target 6.1): s=sin(alpha), p=cos(alpha)*exp(i*delta).
# alpha=0 -> pure p (TM), alpha=90 -> pure s (TE) -- matched to the commercial
# RCWA tool's own polarization-angle convention (0=P, 90=S; confirmed against
# a Lumerical FDTD "grating_power" Rs_power/Rp_power export,
# `R_linear = sin^2(alpha)*Rs_power + cos^2(alpha)*Rp_power`) so a solver
# `linear_Xdeg` state means the same physical input as the commercial tool's
# "Linear X deg" with no angle conversion needed (see decisions.md ADR-033).
# Linear/circular are just special-case (alpha, delta) values of that same
# formula, not separate physics -- so POLARIZATION_STATES_DEG below is a
# single flat list of (name, alpha_deg, delta_deg) triples, all run through
# the one `_jones_state` function. Add/remove/edit rows freely -- as many
# linear angles (delta=0) or general ellipses as you want, in one place.
# ============================================================================
INCIDENT_ANGLE_DEG = 0.0
AZIMUTHAL_ANGLE_DEG = 0.0

POLARIZATION_STATES_DEG = [
    ("TE", 90.0, 0.0),             # alpha=90 -> pure s; delta irrelevant (cos(90)=0 kills p)
    ("TM", 0.0, 0.0),              # alpha=0 -> pure p; delta irrelevant (sin(0)=0 kills s)
    ("RCP", 45.0, 90.0),           # equal split, +90 deg phase -> right-hand circular
    ("LCP", 45.0, -90.0),          # equal split, -90 deg phase -> left-hand circular
    ("linear_15deg", 15.0, 0.0),   # delta=0 -> linear, at whatever angle alpha is
    ("linear_30deg", 30.0, 0.0),
    ("elliptical_a30_d45", 30.0, 45.0),    # general (alpha, delta) -> a genuine ellipse
    ("elliptical_a45_d90", 45.0, 90.0),
    ("elliptical_a60_d135", 60.0, 135.0),
]


def _jones_state(alpha_deg: float, delta_deg: float) -> tuple[complex, complex]:
    """The one formula behind every polarization state in this file --
    s=sin(alpha), p=cos(alpha)*exp(i*delta) -- per CONVENTIONS.md's table
    (ADR-033: alpha=0=P/alpha=90=S, matched to the commercial RCWA tool's
    convention)."""
    alpha, delta = math.radians(alpha_deg), math.radians(delta_deg)
    return math.sin(alpha), math.cos(alpha) * cmath.exp(1j * delta)


POLARIZATION_STATES = {name: _jones_state(a, d) for name, a, d in POLARIZATION_STATES_DEG}

# Pick exactly one state (a POLARIZATION_STATES key, printed by main() if
# unsure of the exact generated name) for the main R/T-vs-wavelength run.
POLARIZATION_STATE = "RCP"
S_AMPLITUDE, P_AMPLITUDE = POLARIZATION_STATES[POLARIZATION_STATE]

# Optionally compare EVERY generated state side-by-side at one wavelength
# (set to None to skip). Per CONVENTIONS.md, at exactly INCIDENT_ANGLE_DEG=0
# every state gives identical R/T (no preferred in-plane direction at normal
# incidence) -- this comparison is only informative away from theta=0.
COMPARE_POLARIZATIONS_AT_M = 550e-9

# ============================================================================
# EDIT (4): wavelength sweep (meters)
# ============================================================================
WAVELENGTHS = np.linspace(0.4e-6, 0.8e-6, 401)

# ============================================================================
# EDIT (5): where to save results (set OUTPUT_CSV_PATH to None to skip
# saving); RUN_NAME tags the output subfolder -- rename it if you copy this
# file for a new stack, so its outputs don't get labeled "custom_multistack".
# ============================================================================
RUN_NAME = "sio2_sio_ni_sio2_on_semi_infinite_si"
OUTPUT_CSV_PATH = "output_multistack_RT.csv"


def build_geometry():
    """Returns (layers, lattice, incidence, transmission)."""
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))  # unused -- only matters for patterned layers
    return layers, lattice, INCIDENCE_MATERIAL, TRANSMISSION_MATERIAL


def main():
    layers, lattice, incidence, transmission = build_geometry()
    sim = Simulation(
        lattice, layers, num_orders=1,
        incidence=incidence, transmission=transmission,
    )

    print(f"Polarization: {POLARIZATION_STATE}  (s_amplitude={S_AMPLITUDE}, p_amplitude={P_AMPLITUDE})")
    print(f"Incidence angle: {INCIDENT_ANGLE_DEG:g} deg\n")

    if COMPARE_POLARIZATIONS_AT_M is not None:
        state_names = list(POLARIZATION_STATES)
        jones_states = [POLARIZATION_STATES[name] for name in state_names]
        comparison = sweep_polarization(
            sim,
            wavelength=COMPARE_POLARIZATIONS_AT_M,
            theta=math.radians(INCIDENT_ANGLE_DEG),
            phi=math.radians(AZIMUTHAL_ANGLE_DEG),
            jones_states=jones_states,
        )
        r_cmp, t_cmp = comparison.reflectance(), comparison.transmittance()
        print(f"Polarization comparison at {COMPARE_POLARIZATIONS_AT_M * 1e9:.1f} nm:")
        print(f"{'state':>12}  {'R':>8}  {'T':>8}  {'A':>8}")
        for name, r_i, t_i in zip(state_names, r_cmp, t_cmp):
            print(f"{name:>12}  {r_i:8.4f}  {t_i:8.4f}  {1 - r_i - t_i:8.4f}")
        print()

    reflectance = np.zeros(len(WAVELENGTHS))
    transmittance = np.zeros(len(WAVELENGTHS))

    print(f"{'wavelength (nm)':>16}  {'R':>8}  {'T':>8}  {'A':>8}")
    for i, wavelength in enumerate(WAVELENGTHS):
        excitation = PlaneWaveExcitation(
            wavelength=wavelength,
            theta=math.radians(INCIDENT_ANGLE_DEG),
            phi=math.radians(AZIMUTHAL_ANGLE_DEG),
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
            layers=[(layer.name, layer.thickness) for layer in layers],
            incidence_material=INCIDENCE_MATERIAL.name,
            transmission_material=TRANSMISSION_MATERIAL.name,
            incident_angle_deg=INCIDENT_ANGLE_DEG,
            azimuthal_angle_deg=AZIMUTHAL_ANGLE_DEG,
            polarization_state=POLARIZATION_STATE,
            s_amplitude=S_AMPLITUDE,
            p_amplitude=P_AMPLITUDE,
            polarization_states_deg=POLARIZATION_STATES_DEG,
            wavelength_range_m=(WAVELENGTHS[0], WAVELENGTHS[-1], len(WAVELENGTHS)),
        )
        absorptance = 1.0 - reflectance - transmittance
        table = np.column_stack([WAVELENGTHS * 1e9, reflectance, transmittance, absorptance])
        output_path = output_dir / OUTPUT_CSV_PATH
        np.savetxt(output_path, table, delimiter=",", header="wavelength_nm,R,T,A", comments="")
        print(f"\nSaved {len(WAVELENGTHS)} rows to {output_path}")

    return reflectance, transmittance


if __name__ == "__main__":
    main()
