"""Reciprocal lattice (G-vector) truncation for Fourier order selection."""

from __future__ import annotations

import numpy as np

from sougata_solver.geometry import Lattice, Lattice1D


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


def truncate_fourier_orders_1d(lattice: Lattice1D, num_orders: int) -> np.ndarray:
    """1D analogue of `truncate_fourier_orders`: same "sort candidate
    G-indices by `|k|`, take the smallest `num_orders`" selection, restricted
    to a single reciprocal direction (`g2` always `0`, matching S4's 1D
    G-vector convention, `S4/S4/S4.cpp:236-251`/`1031-1044`). Returned shape
    is the same `(n, 2)` int array as `truncate_fourier_orders`, so it plugs
    into `fourier_factorization.py`'s `pattern_epsilon_hat`/`toeplitz_matrix`
    unmodified.
    """
    Lk = lattice.reciprocal_vectors()
    radius = num_orders + 2
    g1 = np.arange(-radius, radius + 1)
    kmag = np.abs(g1 * Lk[0, 0])
    order = np.argsort(kmag, kind="stable")

    selected = order[:num_orders]
    g1_sel = g1[selected]
    return np.stack([g1_sel, np.zeros_like(g1_sel)], axis=1)
