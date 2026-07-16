"""Post-processing: load the raw reflected-field data written by
structures/thin_film/sio2_on_si_ellipsometry_run.py, assemble the Jones
reflection matrix, the Mueller matrix, and the standard ellipsometric angles
(Psi, Delta) -- no re-solving of the physics, just extracting these
derived quantities from already-computed raw field data.

Run structures/thin_film/sio2_on_si_ellipsometry_run.py first to produce the
input CSV, then:

Run with:  python postprocessing/jones_mueller_ellipsometry.py
"""

import csv
import math
from collections import defaultdict

from sougata_solver.output_paths import find_latest_output
from sougata_solver.polarimetry import decompose_sp, jones_to_mueller

# ============================================================================
# EDIT (1): filename of the raw CSV produced by the matching "structures"
# script -- looked up automatically under outputs/YYYY-MM-DD/, most recent
# date first, so this doesn't need editing if you ran that script today.
# ============================================================================
INPUT_CSV_FILENAME = "sio2_on_si_ellipsometry_raw.csv"


def _load_raw_fields(csv_path: str):
    """Group raw (Ex, Ey) rows by (wavelength, theta, phi), each group
    holding the 's' and 'p' polarization runs needed to build one Jones
    matrix."""
    groups: dict[tuple[float, float, float], dict[str, complex]] = defaultdict(dict)
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            key = (float(row["wavelength_m"]), float(row["theta_deg"]), float(row["phi_deg"]))
            ex = complex(float(row["Ex_re"]), float(row["Ex_im"]))
            ey = complex(float(row["Ey_re"]), float(row["Ey_im"]))
            groups[key][row["polarization"]] = (ex, ey)
    return groups


def jones_matrix_from_raw(ex_ey_by_pol: dict[str, tuple[complex, complex]], theta: float, phi: float):
    """Assemble the 2x2 Jones reflection matrix `[[rss, rsp], [rps, rpp]]`
    from the raw (Ex, Ey) reflected field of an s-incidence run and a
    p-incidence run, using the same `decompose_sp` convention the solver's
    own `polarimetry.jones_reflection_matrix` uses internally (see
    `src/sougata_solver/polarimetry.py`) -- reused here, not re-derived, so this
    post-processing step can't silently drift from the solver's convention.
    """
    cos_theta = math.cos(theta)
    jones = [[0j, 0j], [0j, 0j]]
    for column, polarization in enumerate(("s", "p")):
        ex, ey = ex_ey_by_pol[polarization]
        e_s, e_p = decompose_sp(ex, ey, phi, cos_theta)  # reflected: +cos(theta)
        jones[0][column] = e_s
        jones[1][column] = e_p
    return jones


def main():
    input_path = find_latest_output(INPUT_CSV_FILENAME)
    print(f"Reading {input_path}")
    groups = _load_raw_fields(input_path)

    for (wavelength, theta_deg, phi_deg), ex_ey_by_pol in sorted(groups.items()):
        if "s" not in ex_ey_by_pol or "p" not in ex_ey_by_pol:
            print(f"Skipping wavelength={wavelength * 1e9:.1f} nm theta={theta_deg} deg: missing s or p run")
            continue

        theta = math.radians(theta_deg)
        phi = math.radians(phi_deg)
        jones = jones_matrix_from_raw(ex_ey_by_pol, theta, phi)
        rss, rsp = jones[0][0], jones[0][1]
        rps, rpp = jones[1][0], jones[1][1]

        print(f"\nwavelength = {wavelength * 1e9:.1f} nm, theta = {theta_deg} deg")
        print("Jones reflection matrix [[rss, rsp], [rps, rpp]]:")
        print(f"  rss = {rss:.6f}   |rss| = {abs(rss):.6f}")
        print(f"  rsp = {rsp:.6f}")
        print(f"  rps = {rps:.6f}")
        print(f"  rpp = {rpp:.6f}   |rpp| = {abs(rpp):.6f}")

        # Standard ellipsometric angles: rho = rpp / rss = tan(Psi) * exp(i*Delta)
        rho = rpp / rss
        psi = math.degrees(math.atan(abs(rho)))
        delta = math.degrees(math.atan2(rho.imag, rho.real))
        print(f"Ellipsometric angles: Psi = {psi:.3f} deg, Delta = {delta:.3f} deg")

        mueller = jones_to_mueller(jones)
        print("Mueller reflection matrix:")
        for row in mueller:
            print("  " + "  ".join(f"{v:9.5f}" for v in row))


if __name__ == "__main__":
    main()
