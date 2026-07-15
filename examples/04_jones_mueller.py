"""Extract the reflection Jones matrix (rss, rsp, rps, rpp) and the
corresponding Mueller matrix for the SiO2/Si stack from 03_sio2_on_si.py.

Run with:  python examples/04_jones_mueller.py
"""

import math

from pyrcwa.geometry import Lattice
from pyrcwa.layer import Layer
from pyrcwa.materials import Material
from pyrcwa.polarimetry import jones_reflection_matrix, jones_to_mueller
from pyrcwa.simulation import Simulation

WAVELENGTH = 0.55e-6
INCIDENT_ANGLE_DEG = 65.0  # ellipsometry is most sensitive away from normal incidence

# Same stack as 03_sio2_on_si.py -- swap in material_from_csv(...) for your
# real Si/SiO2 n,k data as shown there.
air = Material("air", 1.0)
si = Material("Si", (3.9 + 0.02j) ** 2)     # placeholder -- use your real Si n,k
sio2 = Material("SiO2", 1.46**2)             # placeholder -- use your real SiO2 n,k

lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))
layers = [
    Layer("SiO2", 50e-9, material=sio2),
    Layer("Si", 12e-6, material=si),
]
sim = Simulation(lattice, layers, num_orders=1, incidence=air, transmission=air)

theta = math.radians(INCIDENT_ANGLE_DEG)
jones = jones_reflection_matrix(sim, WAVELENGTH, theta, phi=0.0)
rss, rsp = jones[0, 0], jones[0, 1]
rps, rpp = jones[1, 0], jones[1, 1]

print(f"theta = {INCIDENT_ANGLE_DEG} deg, wavelength = {WAVELENGTH * 1e9:.0f} nm\n")
print("Jones reflection matrix [[rss, rsp], [rps, rpp]]:")
print(f"  rss = {rss:.6f}   |rss| = {abs(rss):.6f}   phase = {math.degrees(math.atan2(rss.imag, rss.real)):.2f} deg")
print(f"  rsp = {rsp:.6f}")
print(f"  rps = {rps:.6f}")
print(f"  rpp = {rpp:.6f}   |rpp| = {abs(rpp):.6f}   phase = {math.degrees(math.atan2(rpp.imag, rpp.real)):.2f} deg")

# Standard ellipsometric angles, derived directly from the Jones matrix:
#   rho = rpp / rss = tan(Psi) * exp(i*Delta)
rho = rpp / rss
psi = math.degrees(math.atan(abs(rho)))
delta = math.degrees(math.atan2(rho.imag, rho.real))
print(f"\nEllipsometric angles: Psi = {psi:.3f} deg, Delta = {delta:.3f} deg")

mueller = jones_to_mueller(jones)
print("\nMueller reflection matrix:")
for row in mueller:
    print("  " + "  ".join(f"{v:9.5f}" for v in row))
