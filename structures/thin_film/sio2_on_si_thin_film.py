"""One SiO2 or SiO film on a semi-infinite Si or Ni substrate.

Stack (top to bottom):

    air (semi-infinite) / selected finite film / selected substrate (semi-infinite)

This is the appropriate stack for the KLA/Filmetrics Reflectance Calculator
when its substrate is set to Si or Ni.  The substrate is deliberately passed
as ``transmission=substrate`` rather than added to ``layers``: transmission
media in the solver are semi-infinite half-spaces.
"""

import math
from pathlib import Path

import numpy as np

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.output_paths import run_output_dir, write_run_metadata
from sougata_solver.simulation import Simulation


# KLA material files have columns: wavelength [nm], n, k.
NK_DIR = Path(__file__).resolve().parents[3] / "NK_FILE"
SUBSTRATE_NK_PATHS = {
    "Si": NK_DIR / "si_KLA.txt",
    "Ni": NK_DIR / "ni_KLA.txt",
}
FILM_NK_PATHS = {
    "SiO": NK_DIR / "sio_KLA.txt",
    "SiO2": NK_DIR / "sio2_KLA.txt",
}
NK_WAVELENGTH_UNIT = "nm"

# Select exactly one semi-infinite substrate.  Use "Si" or "Ni".
SUBSTRATE_MATERIAL = "Si"

# Select exactly one finite film.  Use "SiO" or "SiO2".
FILM_MATERIAL = "SiO2"
FILM_THICKNESS_M = 500e-9

# Match these three settings to the KLA calculator before comparing curves.
INCIDENT_ANGLE_DEG = 40.0
AZIMUTHAL_ANGLE_DEG = 0.0
# "s", "p", "mixed" (equal-power average of s and p), "rcp"/"lcp" (circular),
# or "elliptical" (uses ELLIPTICAL_ALPHA_DEG/ELLIPTICAL_DELTA_DEG below).
POLARIZATION = "elliptical"

# Only used when POLARIZATION == "elliptical": s_amplitude = cos(alpha),
# p_amplitude = sin(alpha) * exp(1j * delta), per CONVENTIONS.md's "Worked
# polarization examples" table (delta != 0, pi, or it degenerates to linear).
ELLIPTICAL_ALPHA_DEG = 20.0
ELLIPTICAL_DELTA_DEG = 50.0

# KLA export supports this same 400--800 nm, 1 nm grid.  Adjust if needed.
WAVELENGTHS = np.linspace(400e-9, 800e-9, 401)
OUTPUT_CSV_PATH = "output_R.csv"


def _polarization_amplitudes(polarization: str) -> tuple[complex, complex]:
    """Return one pure polarization state's `(s_amplitude, p_amplitude)`;
    "mixed" light is calculated separately (equal-power average of two
    solves, not a single Jones vector).

    RCP/LCP/elliptical values are exactly `CONVENTIONS.md`'s "Worked
    polarization examples" table (also exercised by
    `tests/test_polarization_states.py`, Category 6 targets 6.2/6.3):
    `PlaneWaveExcitation.s_amplitude`/`p_amplitude` are already complex, so
    circular/elliptical states need no solver changes -- only the phase
    relationship between the two amplitudes.
    """
    if polarization == "s":
        return 1.0, 0.0
    if polarization == "p":
        return 0.0, 1.0
    if polarization == "rcp":
        return 1 / math.sqrt(2), 1j / math.sqrt(2)
    if polarization == "lcp":
        return 1 / math.sqrt(2), -1j / math.sqrt(2)
    if polarization == "elliptical":
        alpha = math.radians(ELLIPTICAL_ALPHA_DEG)
        delta = math.radians(ELLIPTICAL_DELTA_DEG)
        return math.cos(alpha), math.sin(alpha) * complex(np.exp(1j * delta))
    raise ValueError("polarization must be 's', 'p', 'rcp', 'lcp', or 'elliptical'")


def _solve_spectrum(sim: Simulation, polarization: str) -> np.ndarray:
    s_amplitude, p_amplitude = _polarization_amplitudes(polarization)
    reflectance = np.zeros(len(WAVELENGTHS))
    for i, wavelength in enumerate(WAVELENGTHS):
        excitation = PlaneWaveExcitation(
            wavelength=wavelength,
            theta=math.radians(INCIDENT_ANGLE_DEG),
            phi=math.radians(AZIMUTHAL_ANGLE_DEG),
            s_amplitude=s_amplitude,
            p_amplitude=p_amplitude,
        )
        result = sim.solve(excitation)
        reflectance[i] = result.reflectance()
    return reflectance


def build_geometry():
    """Returns (layers, lattice, incidence, transmission)."""
    if SUBSTRATE_MATERIAL not in SUBSTRATE_NK_PATHS:
        raise ValueError(
            f"SUBSTRATE_MATERIAL must be one of {tuple(SUBSTRATE_NK_PATHS)}, got {SUBSTRATE_MATERIAL!r}"
        )
    if FILM_MATERIAL not in FILM_NK_PATHS:
        raise ValueError(f"FILM_MATERIAL must be one of {tuple(FILM_NK_PATHS)}, got {FILM_MATERIAL!r}")

    substrate_nk_path = SUBSTRATE_NK_PATHS[SUBSTRATE_MATERIAL]
    substrate = Material.from_nk_file(SUBSTRATE_MATERIAL, str(substrate_nk_path), NK_WAVELENGTH_UNIT)
    film_nk_path = FILM_NK_PATHS[FILM_MATERIAL]
    film = Material.from_nk_file(FILM_MATERIAL, str(film_nk_path), NK_WAVELENGTH_UNIT)
    air = Material("air", 1.0)

    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))
    layers = [Layer(FILM_MATERIAL, FILM_THICKNESS_M, material=film)]
    # The substrate is semi-infinite, not a finite layer in ``layers``.
    return layers, lattice, air, substrate


def main():
    if POLARIZATION not in ("s", "p", "mixed", "rcp", "lcp", "elliptical"):
        raise ValueError("POLARIZATION must be 's', 'p', 'mixed', 'rcp', 'lcp', or 'elliptical'")

    layers, lattice, air, substrate = build_geometry()
    substrate_nk_path = SUBSTRATE_NK_PATHS[SUBSTRATE_MATERIAL]
    film_nk_path = FILM_NK_PATHS[FILM_MATERIAL]
    sim = Simulation(lattice, layers, num_orders=1, incidence=air, transmission=substrate)

    if POLARIZATION == "mixed":
        reflectance = (_solve_spectrum(sim, "s") + _solve_spectrum(sim, "p")) / 2.0
    else:
        reflectance = _solve_spectrum(sim, POLARIZATION)

    print(f"Stack: air / {FILM_MATERIAL} ({FILM_THICKNESS_M * 1e9:g} nm) / semi-infinite {SUBSTRATE_MATERIAL}")
    print(f"Angle: {INCIDENT_ANGLE_DEG:g} deg; polarization: {POLARIZATION}")
    print(f"{'wavelength (nm)':>16}  {'R':>8}")
    for wavelength, r in zip(WAVELENGTHS, reflectance):
        print(f"{wavelength * 1e9:16.1f}  {r:8.4f}")

    if OUTPUT_CSV_PATH:
        output_dir = run_output_dir(f"{FILM_MATERIAL.lower()}_on_semi_infinite_{SUBSTRATE_MATERIAL.lower()}")
        write_run_metadata(
            output_dir,
            __file__,
            stack=f"air / {FILM_MATERIAL} / semi-infinite {SUBSTRATE_MATERIAL}",
            substrate_material=SUBSTRATE_MATERIAL,
            substrate_is_semi_infinite=True,
            substrate_nk_path=str(substrate_nk_path),
            film_material=FILM_MATERIAL,
            film_nk_path=str(film_nk_path),
            film_thickness_m=FILM_THICKNESS_M,
            incident_angle_deg=INCIDENT_ANGLE_DEG,
            azimuthal_angle_deg=AZIMUTHAL_ANGLE_DEG,
            polarization=POLARIZATION,
            wavelength_range_m=(WAVELENGTHS[0], WAVELENGTHS[-1], len(WAVELENGTHS)),
        )
        table = np.column_stack([WAVELENGTHS * 1e9, reflectance])
        output_path = output_dir / OUTPUT_CSV_PATH
        np.savetxt(output_path, table, delimiter=",", header="wavelength_nm,R", comments="")
        print(f"\nSaved {len(WAVELENGTHS)} rows to {output_path}")
        print(f"Run metadata: {output_dir / 'run_metadata.txt'}")

    return reflectance


if __name__ == "__main__":
    main()
