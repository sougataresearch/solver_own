"""Reciprocal lattice (G-vector) truncation for Fourier order selection."""

from __future__ import annotations

import numpy as np

from sougata_solver.geometry import Lattice


def truncate_fourier_orders(lattice: Lattice, num_orders: int, method: str = "circular") -> np.ndarray:
    """Return an `(n, 2)` integer array of `(g1, g2)` reciprocal lattice
    indices, sorted by increasing `|k|` and truncated to `num_orders`.

    `method="circular"` selects indices whose reciprocal-space magnitude is
    smallest, matching S4's default G-vector selection (`gsel.c`).
    """
    if method != "circular":
        raise NotImplementedError(f"truncation method {method!r} not implemented")

    Lk = lattice.reciprocal_vectors()
    # A (2*radius+1)^2 integer box always contains more than num_orders
    # points once radius >= sqrt(num_orders); pad generously for safety.
    radius = max(1, int(np.ceil(np.sqrt(num_orders))) + 2)
    g1, g2 = np.meshgrid(np.arange(-radius, radius + 1), np.arange(-radius, radius + 1), indexing="ij")
    g1 = g1.ravel()
    g2 = g2.ravel()
    k = g1[:, None] * Lk[0] + g2[:, None] * Lk[1]
    kmag = np.linalg.norm(k, axis=1)
    order = np.argsort(kmag, kind="stable")

    selected = order[:num_orders]
    return np.stack([g1[selected], g2[selected]], axis=1)
