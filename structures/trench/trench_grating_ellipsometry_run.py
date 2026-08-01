"""Build the 1D lamellar grating (see trench_grating.py), run the solver at
one or more (wavelength, angle) points for both s- and p-polarized
incidence, and save the *raw* reflected transverse field data -- for every
diffraction order, not just the zeroth -- to a CSV.

This is the "build & run" half of a per-order ellipsometry measurement: it
produces raw field data only. Turning that raw data into a per-order Jones
matrix, Mueller matrix, and ellipsometric angles (Psi, Delta) is
post-processing -- see postprocessing/jones_mueller_per_order.py, which
reads the CSV this script writes.

Per-order isolation reuses the identical masking technique already
validated in `SimulationResult.diffraction_efficiencies` (`simulation.py`)
and `polarimetry.jones_reflection_matrix_by_order` -- zero every other
order's mode amplitude before reading the tangential field, rather than a
new per-order formula.

Run with:  python structures/trench/trench_grating_ellipsometry_run.py
"""

import math

import numpy as np

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.fields import tangential_e_field
from sougata_solver.geometry import Lattice1D, Pattern, Slab
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.output_paths import run_output_dir, write_run_metadata
from sougata_solver.simulation import Simulation

# ============================================================================
# EDIT (1): grating geometry and materials (same defaults as trench_grating.py)
# ============================================================================
PERIOD = 0.7e-6          # grating period (meters)
FILL_FACTOR = 0.3        # fraction of the period occupied by the ridge
THICKNESS = 0.46e-6      # grating (groove) depth (meters)
N_RIDGE = 3.48           # ridge index (e.g. Si)

# ============================================================================
# EDIT (2): measurement points and order truncation
# ============================================================================
WAVELENGTHS = [0.8e-6]           # meters; add more points for a spectral sweep
INCIDENT_ANGLES_DEG = [0.0]      # degrees from surface normal
AZIMUTHAL_ANGLE_DEG = 0.0
NUM_ORD = 5                      # orders per side; total Fourier orders = 2*NUM_ORD+1

# ============================================================================
# EDIT (3): where to save the raw field data
# ============================================================================
OUTPUT_CSV_PATH = "trench_grating_ellipsometry_raw.csv"  # saved under outputs/YYYY_MM_DD/


def main():
    air = Material("air", 1.0)
    ridge = Material("ridge", N_RIDGE**2)

    pattern = Pattern(background=air)
    pattern.add(Slab(center_x=-PERIOD * (1 - FILL_FACTOR) / 2, halfwidth=0.5 * FILL_FACTOR * PERIOD, material=ridge))

    lattice = Lattice1D(PERIOD)
    layers = [Layer("grating", THICKNESS, pattern=pattern)]
    num_orders = 2 * NUM_ORD + 1
    sim = Simulation(lattice, layers, num_orders=num_orders, incidence=air, transmission=air)

    rows = []
    for wavelength in WAVELENGTHS:
        for theta_deg in INCIDENT_ANGLES_DEG:
            theta = math.radians(theta_deg)
            phi = math.radians(AZIMUTHAL_ANGLE_DEG)
            for polarization, (s_amp, p_amp) in [("s", (1.0, 0.0)), ("p", (0.0, 1.0))]:
                excitation = PlaneWaveExcitation(wavelength, theta, phi, s_amplitude=s_amp, p_amplitude=p_amp)
                result = sim.solve(excitation)
                modes_inc = result.all_modes[0]
                omega = excitation.omega()
                zeros = np.zeros_like(result.a0)
                block = result.b_reflected.shape[0] // 2

                for i in range(result.num_orders):
                    g1, g2 = int(result.g[i, 0]), int(result.g[i, 1])
                    b_masked = np.zeros_like(result.b_reflected)
                    idx = [i, i + block]
                    b_masked[idx] = result.b_reflected[idx]
                    ex, ey = tangential_e_field(omega, modes_inc.q, modes_inc.kp, modes_inc.phi, zeros, b_masked)
                    rows.append(
                        [wavelength, theta_deg, AZIMUTHAL_ANGLE_DEG, g1, g2, polarization,
                         ex[i].real, ex[i].imag, ey[i].real, ey[i].imag]
                    )
                print(f"wavelength={wavelength * 1e9:.1f} nm  theta={theta_deg} deg  pol={polarization}  "
                      f"solved {result.num_orders} orders")

    output_dir = run_output_dir("trench_grating_ellipsometry_run")
    write_run_metadata(
        output_dir,
        __file__,
        wavelengths_m=WAVELENGTHS,
        incident_angles_deg=INCIDENT_ANGLES_DEG,
        azimuthal_angle_deg=AZIMUTHAL_ANGLE_DEG,
        period_m=PERIOD,
        fill_factor=FILL_FACTOR,
        thickness_m=THICKNESS,
        n_ridge=N_RIDGE,
        num_orders=num_orders,
    )
    output_path = output_dir / OUTPUT_CSV_PATH
    with open(output_path, "w") as f:
        f.write("wavelength_m,theta_deg,phi_deg,order_g1,order_g2,polarization,Ex_re,Ex_im,Ey_re,Ey_im\n")
        for row in rows:
            f.write(",".join(str(v) for v in row) + "\n")
    print(f"\nSaved {len(rows)} raw field rows ({num_orders} orders x 2 pol x "
          f"{len(WAVELENGTHS)} wavelengths x {len(INCIDENT_ANGLES_DEG)} angles) to {output_path}")


if __name__ == "__main__":
    main()
