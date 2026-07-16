"""Minimal usage example: reflectance/transmittance of a uniform multilayer
stack (an anti-reflection coating), the Phase 1 capability of sougata_solver.

Run with:  python structures/thin_film/anti_reflection_coating.py
"""

import math

from sougata_solver.excitation import PlaneWaveExcitation
from sougata_solver.geometry import Lattice
from sougata_solver.layer import Layer
from sougata_solver.materials import Material
from sougata_solver.simulation import Simulation

WAVELENGTH = 0.55e-6  # 550 nm, green light

# The in-plane lattice only matters once layers are patterned (Phase 2+);
# for a uniform stack any nonzero lattice works, and num_orders=1 (just the
# zeroth order) is all that's needed.
lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))

air = Material("air", 1.0)
glass = Material("glass", 1.5**2)

# Single quarter-wave MgF2 anti-reflection coating on glass.
mgf2_thickness = WAVELENGTH / (4 * 1.38)
layers = [Layer("MgF2", mgf2_thickness, material=Material("MgF2", 1.38**2))]

sim = Simulation(lattice, layers, num_orders=1, incidence=air, transmission=glass)

for theta_deg in (0.0, 30.0, 60.0):
    excitation = PlaneWaveExcitation(
        wavelength=WAVELENGTH,
        theta=math.radians(theta_deg),
        phi=0.0,
        s_amplitude=1.0,  # s-polarized (TE) incident light
        p_amplitude=0.0,
    )
    result = sim.solve(excitation)
    r, t = result.reflectance(), result.transmittance()
    print(f"theta={theta_deg:5.1f} deg   R={r:.4f}   T={t:.4f}   R+T={r + t:.4f}")
