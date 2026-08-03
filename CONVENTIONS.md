# Solver conventions

## Scope and status

This document records the conventions that are already implemented by the
isotropic uniform, 1D, and 2D patterned RCWA paths. It is descriptive only:
it does not add anisotropic physics or change any numerical result. Future
tensor work must retain these conventions or explicitly document a validated
conversion at its interface.

## Coordinates, units, and excitation

- Layers are laterally periodic in `x`/`y` and are stacked along `z`; a
  `LayerStack` is ordered from the incidence half-space at the left/top of
  the stack to the transmission half-space at the right/bottom.
- Length inputs, including wavelength, layer thickness, lattice vectors, and
  in-plane wavevectors, use SI metres and radians per metre.
- `PlaneWaveExcitation.omega()` returns `omega = 2*pi/wavelength`. This is a
  vacuum angular **wavenumber**, using the solver's natural-unit convention
  `c = 1`, rather than an angular frequency in s^-1.
- `theta` is the polar angle from `+z`; `phi` is the azimuth from `+x` toward
  `+y`. For incident refractive index `n`, the zeroth-order in-plane vector is
  `kx0 = omega*n*sin(theta)*cos(phi)` and
  `ky0 = omega*n*sin(theta)*sin(phi)`.
- The reciprocal-lattice Fourier orders shift those incident components to
  form each order's `(kx, ky)`.

## Phasor and propagation convention

- The implemented frequency-domain convention uses `d/dt -> -i*omega` and a
  forward layer propagation factor `exp(+i*q*z)`. This is the convention
  used by the S-matrix and eigenmode paths; it is not interchangeable with a
  textbook convention that uses the opposite time sign without changing the
  associated spatial signs.
- `q` is the z-directed modal propagation constant. `_select_q_branch` uses
  the outgoing/decaying branch: real propagating roots are non-negative,
  purely evanescent roots have positive imaginary part, and complex roots are
  sign-flipped when needed so `Im(q) >= 0`.
- Consequently, propagation through a finite positive thickness uses
  `exp(+i*q*thickness)`. A mode with positive imaginary `q` decays in the
  forward `+z` direction.

The branch rule and propagation factor are transcribed from
`S4/S4/rcwa.cpp::SolveLayerEigensystem_uniform` (lines 422–502),
`SolveLayerEigensystem` (lines 684–827), and the S4 propagation convention
used by `smatrix.py::propagation_smatrix`.

## Modal vectors and normalization

- Every layer has `2*n` modes for `n` retained Fourier orders. `q` has shape
  `(2*n,)`; `phi` and `kp` have shape `(2*n, 2*n)`.
- `phi` is an internal eigenvector basis. It is **not** individually
  power-normalized, and users must not interpret `phi @ a` directly as the
  tangential electric field.
- For forward amplitude `a` and backward amplitude `b`, the implemented
  tangential electric-field reconstruction is

  ```text
  u  = kp @ phi @ ((a - b) / (omega*q))
  Ex = u[n:2*n]
  Ey = -u[0:n]
  ```

  Thus the internal transverse ordering is `u = [-Ey; Ex]` per Fourier-order
  block. The corresponding tangential field-continuity matrices are
  `A = phi` and `B = kp @ phi @ diag(1/q)`.
- Incident amplitudes are normalized operationally: `incident_mode_amplitude`
  solves for `a0` that produces the caller-requested zeroth-order `(Ex, Ey)`
  with all other orders zero. Reported R/T is then normalized by the incident
  z-Poynting flux, not by a unit-amplitude or unit-power eigenvector rule.

These field/amplitude relations are transcribed from
`S4/S4/rcwa.cpp::GetInPlaneFieldVector` (lines 1959–1995),
`GetZPoyntingFlux` (lines 1846–1897), and `GetSMatrix` (lines 936–1096).

## S-matrix direction convention

With `a` denoting forward-going and `b` backward-going amplitudes, each
scattering matrix is stored as

```text
[a_right; b_left] = S @ [a_left; b_right]
```

The inputs are incoming amplitudes and the outputs are outgoing amplitudes.
The full stack is assembled using the Redheffer star product. This convention
is fixed across all layers and Fourier orders; a new solver path must produce
`q`, `phi`, and `kp` compatible with it rather than introducing a separate
amplitude normalization.

The interface and star-product formulas are verified against
`S4/S4/rcwa.cpp::GetSMatrix` (lines 936–1096) and
`S4/S4r/StarProduct.hpp` (`T2Sblocks`, lines 51–65; `StarProduct`,
lines 83–110).

## Polarization convention

- The incident electric field is built from the internal basis
  `s_hat = (-sin(phi), cos(phi), 0)` and transverse
  `p_hat_xy = -cos(theta)*(cos(phi), sin(phi))`.
- Complex `s_amplitude` and `p_amplitude` therefore represent TE/TM, linear,
  circular, and elliptical input states.
- This sign/phase convention is internally self-consistent and validated for
  power quantities. It has **not** yet been externally matched to an S4 or
  EMpy polarization-phase convention, so phase-sensitive comparisons must
  include an explicit convention check.

## Permittivity tensor ordering

`Material.epsilon_tensor(wavelength)` returns a complex `(3, 3)` array in
Cartesian order:

```text
[[epsilon_xx, epsilon_xy, epsilon_xz],
 [epsilon_yx, epsilon_yy, epsilon_yz],
 [epsilon_zx, epsilon_zy, epsilon_zz]]
```

Today, only scalar-isotropic material values enter `Simulation.solve`;
non-isotropic layers intentionally raise `NotImplementedError`. The tensor
ordering is documented now so the later anisotropy milestones can be tested
without ambiguity. Storing a tensor must not be mistaken for solving it.

## Known limitation

This document freezes the existing isotropic RCWA conventions. It does not
claim that a full anisotropic formulation, exact degenerate-mode treatment,
or public propagating/evanescent classification has been implemented; those
are separate Category-1 targets.
