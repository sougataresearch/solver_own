"""Template: build a Material from your own measured (wavelength, n, k)
dataset instead of a constant refractive index.

Edit the three places marked EDIT below: (1) how you load your data,
(2) the substrate/incidence materials, (3) the layer stack and thicknesses.

Run with:  python examples/02_custom_nk_dataset.py
"""

import math

import numpy as np
from scipy.interpolate import interp1d

from pyrcwa.excitation import PlaneWaveExcitation
from pyrcwa.geometry import Lattice
from pyrcwa.layer import Layer
from pyrcwa.materials import Material
from pyrcwa.simulation import Simulation


def _parse_refractiveindex_csv(csv_path: str):
    """Parse the refractiveindex.info CSV export format: an 'n' block
    (header 'wl,n') and, optionally, a separate 'k' block (header 'wl,k')
    after a blank line -- the two blocks can even be on different
    wavelength grids. If there's no 'k' block, k is treated as zero
    (lossless) at every wavelength.
    """
    with open(csv_path) as f:
        lines = [line.strip() for line in f if line.strip()]

    blocks: list[list[str]] = []
    current: list[str] = []
    for line in lines:
        try:
            float(line.split(",")[0])
        except ValueError:
            if current:
                blocks.append(current)
            current = []
        else:
            current.append(line)
    if current:
        blocks.append(current)
    if not blocks:
        raise ValueError(f"No numeric data found in {csv_path!r}")

    def parse_block(block):
        arr = np.array([[float(x) for x in line.split(",")] for line in block])
        return arr[:, 0], arr[:, 1]

    wl_n, n_vals = parse_block(blocks[0])
    if len(blocks) > 1:
        wl_k, k_vals = parse_block(blocks[1])
    else:
        wl_k, k_vals = wl_n, np.zeros_like(wl_n)
    return wl_n, n_vals, wl_k, k_vals


def material_from_csv(name: str, csv_path: str, wavelength_unit: str = "um") -> Material:
    """Build a dispersive Material from a refractiveindex.info-style CSV
    (an 'n' block and an optional 'k' block).

    `wavelength_unit` is the unit used in the CSV's first column ("um" or
    "nm"); pyrcwa itself always works in meters internally.
    """
    wl_n, n_vals, wl_k, k_vals = _parse_refractiveindex_csv(csv_path)
    scale = {"um": 1e-6, "nm": 1e-9, "m": 1.0}[wavelength_unit]

    n_interp = interp1d(wl_n * scale, n_vals, bounds_error=False, fill_value="extrapolate")
    k_interp = interp1d(wl_k * scale, k_vals, bounds_error=False, fill_value="extrapolate")

    def n_of_wavelength(wavelength: float) -> float:
        return float(n_interp(wavelength))

    def k_of_wavelength(wavelength: float) -> float:
        return float(k_interp(wavelength))

    return Material.from_nk(name, n_of_wavelength, k_of_wavelength)


# You can also build a dispersive Material directly from in-memory arrays,
# without a CSV file, the same way:
def material_from_arrays(name: str, wavelengths_m: np.ndarray, n_vals: np.ndarray, k_vals: np.ndarray) -> Material:
    n_interp = interp1d(wavelengths_m, n_vals, bounds_error=False, fill_value="extrapolate")
    k_interp = interp1d(wavelengths_m, k_vals, bounds_error=False, fill_value="extrapolate")
    return Material.from_nk(
        name,
        lambda wl: float(n_interp(wl)),
        lambda wl: float(k_interp(wl)),
    )


if __name__ == "__main__":
    # ---- EDIT (1): load your sample's n,k data --------------------------
    # Option A, from a CSV file (wavelength_um, n, k):
    #   sample = material_from_csv("my_sample", r"C:\path\to\my_nk_data.csv", wavelength_unit="um")
    # Option B, from arrays you already have in Python:
    wavelengths_um = np.array([0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70])
    n_vals = np.array([2.10, 2.05, 2.00, 1.98, 1.95, 1.93, 1.90])
    k_vals = np.array([0.30, 0.20, 0.12, 0.08, 0.05, 0.03, 0.02])
    sample = material_from_arrays("my_sample", wavelengths_um * 1e-6, n_vals, k_vals)

    # ---- EDIT (2): incidence medium and substrate ------------------------
    air = Material("air", 1.0)
    substrate = Material("substrate", 1.5**2)  # replace with your real substrate index

    # ---- EDIT (3): layer stack (thickness in meters, incidence -> substrate order)
    thickness = 100e-9  # 100 nm film of your sample
    layers = [Layer("sample_film", thickness, material=sample)]

    lattice = Lattice((1e-6, 0.0), (0.0, 1e-6))  # unused for uniform (unpatterned) layers
    sim = Simulation(lattice, layers, num_orders=1, incidence=air, transmission=substrate)

    for wavelength_um in [0.45, 0.55, 0.65]:
        excitation = PlaneWaveExcitation(
            wavelength=wavelength_um * 1e-6,
            theta=math.radians(0.0),
            phi=0.0,
            s_amplitude=1.0,
            p_amplitude=0.0,
        )
        result = sim.solve(excitation)
        r, t = result.reflectance(), result.transmittance()
        print(f"wavelength={wavelength_um:.3f} um   R={r:.4f}   T={t:.4f}   A={1 - r - t:.4f}")
