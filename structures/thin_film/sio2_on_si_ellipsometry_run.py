"""Build the SiO2/Si stack (see sio2_on_si_thin_film.py), run the solver at
one or more (wavelength, angle) points for both s- and p-polarized
incidence, and save the *raw* reflected transverse field data to a CSV.

This is the "build & run" half of an ellipsometry measurement: it produces
raw field data only. Turning that raw data into a Jones matrix, a Mueller
matrix, or ellipsometric angles (Psi, Delta) is post-processing -- see
postprocessing/jones_mueller_ellipsometry.py, which reads the CSV this
script writes.

Run with:  python structures/thin_film/sio2_on_si_ellipsometry_run.py
"""

import math

import numpy as np

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.fields import tangential_e_field
from sougata_solver.geometry import Lattice
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.output_paths import run_output_dir, write_run_metadata
from sougata_solver.simulation import Simulation

# ============================================================================
# EDIT (1): materials -- swap in material_from_csv(...) (see
# custom_material_from_nk_data.py) for your real Si/SiO2 n,k data.
# ============================================================================
air = Material("air", 1.0)
si = Material("Si", (3.9 + 0.02j) ** 2)   # placeholder -- use your real Si n,k
sio2 = Material("SiO2", 1.46**2)          # placeholder -- use your real SiO2 n,k

# ============================================================================
# EDIT (2): the stack (thickness in meters, top to bottom)
# ============================================================================
layers = [
    Layer("SiO2", 50e-9, material=sio2),
    Layer("Si", 12e-6, material=si),
]

# ============================================================================
# EDIT (3): measurement points -- ellipsometry is most sensitive away from
# normal incidence, hence the default 65 degrees.
# ============================================================================
WAVELENGTHS = [0.55e-6]          # meters; add more points for a spectral sweep
INCIDENT_ANGLES_DEG = [65.0]     # degrees from surface normal
AZIMUTHAL_ANGLE_DEG = 0.0

# ============================================================================
# EDIT (4): where to save the raw field data
# ============================================================================
OUTPUT_CSV_PATH = "sio2_on_si_ellipsometry_raw.csv"  # filename only; saved under outputs/YYYY_MM_DD/


def main():
    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))  # unused for uniform (unpatterned) layers
    sim = Simulation(lattice, layers, num_orders=1, incidence=air, transmission=air)

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
                ex, ey = tangential_e_field(
                    omega, modes_inc.q, modes_inc.kp, modes_inc.phi, zeros, result.b_reflected
                )
                i = result.zeroth_order_index
                rows.append(
                    [wavelength, theta_deg, AZIMUTHAL_ANGLE_DEG, polarization, ex[i].real, ex[i].imag, ey[i].real, ey[i].imag]
                )
                print(
                    f"wavelength={wavelength * 1e9:.1f} nm  theta={theta_deg} deg  "
                    f"pol={polarization}  Ex={ex[i]:.6f}  Ey={ey[i]:.6f}"
                )

    output_dir = run_output_dir("sio2_on_si_ellipsometry_run")
    write_run_metadata(
        output_dir,
        __file__,
        wavelengths_m=WAVELENGTHS,
        incident_angles_deg=INCIDENT_ANGLES_DEG,
        azimuthal_angle_deg=AZIMUTHAL_ANGLE_DEG,
        layers=[(layer.name, layer.thickness) for layer in layers],
    )
    output_path = output_dir / OUTPUT_CSV_PATH
    with open(output_path, "w") as f:
        f.write("wavelength_m,theta_deg,phi_deg,polarization,Ex_re,Ex_im,Ey_re,Ey_im\n")
        for row in rows:
            f.write(f"{row[0]},{row[1]},{row[2]},{row[3]},{row[4]},{row[5]},{row[6]},{row[7]}\n")
    print(f"\nSaved {len(rows)} raw field rows to {output_path}")


if __name__ == "__main__":
    main()
