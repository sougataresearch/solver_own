"""Post-processing: load the raw per-order reflected-field data written by
structures/trench/trench_grating_ellipsometry_run.py and assemble a Jones
reflection matrix, a Mueller matrix, and the ellipsometric angles (Psi,
Delta) for *every* diffraction order -- no re-solving of the physics, just
extracting these derived quantities from already-computed raw field data.

Reuses `decompose_sp`/`jones_to_mueller` from `sougata_solver.polarimetry`,
the same functions `postprocessing/jones_mueller_ellipsometry.py` uses for
the zeroth-order-only (uniform-stack) case, so this can't silently drift
from the solver's own convention.

Run structures/trench/trench_grating_ellipsometry_run.py first to produce
the input CSV, then:

Run with:  python postprocessing/jones_mueller_per_order.py
"""

import csv
import math
from collections import defaultdict

from sougata_solver.output_paths import find_latest_output
from sougata_solver.polarimetry import decompose_sp, jones_to_mueller

# ============================================================================
# EDIT (1): filename of the raw CSV produced by the matching "structures"
# script -- looked up automatically under outputs/YYYY_MM_DD/, most recent
# date first, so this doesn't need editing if you ran that script today.
# ============================================================================
INPUT_CSV_FILENAME = "trench_grating_ellipsometry_raw.csv"


def _load_raw_fields(csv_path: str):
    """Group raw (Ex, Ey) rows by (wavelength, theta, phi, order), each
    group holding the 's' and 'p' polarization runs needed to build one
    per-order Jones matrix."""
    groups: dict[tuple[float, float, float, int, int], dict[str, tuple[complex, complex]]] = defaultdict(dict)
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            key = (
                float(row["wavelength_m"]), float(row["theta_deg"]), float(row["phi_deg"]),
                int(row["order_g1"]), int(row["order_g2"]),
            )
            ex = complex(float(row["Ex_re"]), float(row["Ex_im"]))
            ey = complex(float(row["Ey_re"]), float(row["Ey_im"]))
            groups[key][row["polarization"]] = (ex, ey)
    return groups


def jones_matrix_from_raw(ex_ey_by_pol: dict[str, tuple[complex, complex]], theta: float, phi: float):
    """Assemble the 2x2 Jones reflection matrix `[[rss, rsp], [rps, rpp]]`
    for one diffraction order from its raw (Ex, Ey) s-incidence and
    p-incidence runs, using the same `decompose_sp` convention the
    solver's own `polarimetry.jones_reflection_matrix_by_order` uses
    internally."""
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

    for (wavelength, theta_deg, phi_deg, g1, g2), ex_ey_by_pol in sorted(groups.items()):
        if "s" not in ex_ey_by_pol or "p" not in ex_ey_by_pol:
            print(f"Skipping wavelength={wavelength * 1e9:.1f} nm theta={theta_deg} deg order=({g1},{g2}): "
                  "missing s or p run")
            continue

        theta = math.radians(theta_deg)
        phi = math.radians(phi_deg)
        jones = jones_matrix_from_raw(ex_ey_by_pol, theta, phi)
        rss, rsp = jones[0][0], jones[0][1]
        rps, rpp = jones[1][0], jones[1][1]

        print(f"\nwavelength = {wavelength * 1e9:.1f} nm, theta = {theta_deg} deg, order = ({g1}, {g2})")
        print("Jones reflection matrix [[rss, rsp], [rps, rpp]]:")
        print(f"  rss = {rss:.6f}   |rss| = {abs(rss):.6f}")
        print(f"  rsp = {rsp:.6f}")
        print(f"  rps = {rps:.6f}")
        print(f"  rpp = {rpp:.6f}   |rpp| = {abs(rpp):.6f}")

        # Standard ellipsometric angles: rho = rpp / rss = tan(Psi) * exp(i*Delta)
        # Only meaningful when rss is non-negligible -- off-specular orders
        # can have rss close to zero, where rho blows up; that's a physical
        # feature of this order, not a bug, so it's printed as-is.
        rho = rpp / rss if abs(rss) > 1e-12 else complex("nan")
        psi = math.degrees(math.atan(abs(rho)))
        delta = math.degrees(math.atan2(rho.imag, rho.real))
        print(f"Ellipsometric angles: Psi = {psi:.3f} deg, Delta = {delta:.3f} deg")

        mueller = jones_to_mueller(jones)
        print("Mueller reflection matrix:")
        for row in mueller:
            print("  " + "  ".join(f"{v:9.5f}" for v in row))


if __name__ == "__main__":
    main()
