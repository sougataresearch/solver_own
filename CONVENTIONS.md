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

### Worked polarization examples (Category 6 target 6.1)

`PlaneWaveExcitation(s_amplitude, p_amplitude)` constructs any fully-
polarized input state via these amplitude pairs (`|s_amplitude|^2 +
|p_amplitude|^2` sets total incident power, per `incident_mode_amplitude`'s
docstring):

| State | `s_amplitude` | `p_amplitude` |
|---|---|---|
| TE (s-pol) | `1.0` | `0.0` |
| TM (p-pol) | `0.0` | `1.0` |
| Linear at 45 deg | `1/sqrt(2)` | `1/sqrt(2)` |
| Linear at an arbitrary angle `alpha` | `cos(alpha)` | `sin(alpha)` |
| Right-hand circular (RCP) | `1/sqrt(2)` | `1j/sqrt(2)` |
| Left-hand circular (LCP) | `1/sqrt(2)` | `-1j/sqrt(2)` |
| Elliptical | `cos(alpha)` | `sin(alpha) * exp(1j*delta)`, `delta != 0, pi` |

These are the states `tests/test_polarization_states.py` (Category 6
targets 6.2/6.3) exercises: normal-incidence regression checks that, for
an isotropic (non-birefringent) stack, `R`/`T` are identical across every
row above at fixed total power (a direct consequence of the stack having
no preferred in-plane direction at `theta=0` — see that file for the
proof-by-symmetry argument); oblique-incidence regression additionally
varies `phi` and checks energy conservation for mixed/elliptical states.

## Real-space field reconstruction (Category 9 target 9.1, Phase 7)

`fields.py`'s `modal_field_components`/`reconstruct_field_at_points` extend
the already-established modal-vector conventions above (`u = [-Ey; Ex]`,
`Ex = u[n:2n]`, `Ey = -u[0:n]`) to a full 6-component `(Ex,Ey,Ez,Hx,Hy,Hz)`
per-order Fourier representation, then sum that Fourier series onto a
real-space `(x, y)` point or grid. Transcribed from
`S4/S4/rcwa.cpp::GetInPlaneFieldVector` (lines 1959-1995, the transverse
components — this project's existing `tangential_e_field` already
implements the `E` half of this) and `GetFieldAtPoint` (lines 1997-2074,
the longitudinal components and the real-space phase sum).

- **Transverse components** (per Fourier order `i`, forward amplitude `a`,
  backward amplitude `b`, both referenced at the same `z`-plane):
  ```text
  Hx, Hy = split(phi @ (a + b))          # NOT E, despite the naive-looking
                                          # resemblance -- see troubleshooting.md
  Ex, Ey = tangential_e_field(...)        # already implemented, unchanged
  ```
- **Longitudinal components** (from the source-free Maxwell curl equations,
  i.e. `del x H = -i*omega*eps*E` restricted to its z-component and
  `del x E = i*omega*mu0*H` similarly, evaluated per Fourier order so
  `del -> i*(kx, ky, d/dz)` in-plane):
  ```text
  Ez = epsilon_inv @ (ky*Hx - kx*Hy) / omega   # epsilon_inv: scalar (uniform
                                                 # isotropic) or (n,n) matrix
                                                 # (patterned/anisotropic),
                                                 # same dispatch as
                                                 # eigenmodes.build_kp_matrix
  Hz = (kx*Ey - ky*Ex) / omega
  ```
- **Real-space phase sum** (the actual "reconstruction" step; `kx`, `ky` are
  this project's already-angular per-order wavevectors, the same arrays
  `eigenmodes.py`/`simulation.py` already use — **not** `geometry.py`'s
  separate "cycles per unit length" convention for shape Fourier
  transforms, a different, already-documented convention used only inside
  `Shape.fourier_transform`/`fourier_factorization.py`):
  ```text
  F(x, y) = sum_i F_i * exp(i*(kx_i*x + ky_i*y))
  ```
  for each field component `F in {Ex,Ey,Ez,Hx,Hy,Hz}`.
- **Depth dependence within a layer** (independently derived — see
  `decisions.md` ADR-015 for why this wasn't transcribed and how it was
  validated instead — algebraically consistent with `smatrix.propagation_smatrix`'s
  already-established `exp(+i*q*thickness)` convention, confirmed directly,
  not assumed): given forward/backward amplitudes `(a_top, b_top)` at a
  layer's top reference plane (`z=0` local to that layer, exactly what
  `smatrix.interior_amplitudes` recovers), the amplitudes at depth `z`
  (`0 <= z <= thickness`) are
  ```text
  a(z) = a_top * exp(+i*q*z)
  b(z) = b_top * exp(-i*q*z)
  ```

**Real-space Poynting flux, and the missing factor of `0.5` (Category 9
target 9.6, found this session, confirmed directly not assumed)**: the
textbook time-averaged Poynting flux is `Sz = 0.5*Re(Ex*conj(Hy) -
Ey*conj(Hx))`. `fields.z_poynting_flux`'s own already-oracle-validated
modal quadratic form does **not** include that `0.5` -- confirmed by
direct comparison for a single-order uniform-layer case:
`z_poynting_flux`'s output is exactly `2x` the textbook formula evaluated
on `modal_field_components`' own output for the same mode. This never
affected any existing result because `reflectance()`/`transmittance()`
only ever use *ratios* of `z_poynting_flux` outputs (the factor of 2
cancels), but it matters for `tests/test_field_reconstruction.py`'s target
9.6 real-space flux integral, which must therefore use
`Sz = Re(Ex*conj(Hy) - Ey*conj(Hx))` (no `0.5`) to match this project's
established convention -- verified to reproduce the solver's own `R`/`T`
to full double precision once the missing factor is accounted for, not
merely "close."

## Bottom (reverse-side) illumination (Category 6 target 6.6)

There is no separate "illuminate from below" mode. Reverse the layer list
and swap which material plays `incidence`/`transmission`:
`Simulation(lattice, list(reversed(layers)), num_orders, incidence=<old
transmission material>, transmission=<old incidence material>)` — a
`Layer`'s `thickness`/`pattern` carry no inherent z-direction, so this is a
complete, correct description of the mirrored problem. See `decisions.md`
ADR-014 and `tests/test_bottom_incidence.py` for the reciprocity check that
validates this recipe.

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
