# API Reference — `sougata_solver`

Target 18.4 of `COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 18. This
expands `src/sougata_solver/README.md`'s existing Module Map into a full
per-symbol reference (signatures, units, exceptions) — see that file for
the one-line responsibility summary this builds on, and
[`theory.md`](theory.md) for *why* each stage exists. Only the **public**
surface is listed (no leading-underscore helpers); `design.md`'s Public
API Inventory (Category 15 target 15.1) is the authoritative "what's
stable to depend on" table if you need that distinction.

**Units, once, not repeated per function**: lengths (wavelength,
thickness, lattice vectors) are SI metres; angles are radians;
`omega = 2*pi/wavelength` is a vacuum wavenumber in the solver's natural
units (`c=1`), not an angular frequency in `s^-1` — see `CONVENTIONS.md`.

**Exceptions, once, not repeated per function**: see `design.md`'s
**Failure Contract** table for the complete, tested inventory of which
condition raises `ValueError` vs. `NotImplementedError` vs.
`numpy.linalg.LinAlgError` vs. only a `WARNING` log. Only unusual,
function-specific exception behavior is called out below.

## `materials.py` — dispersive and tensor materials

| Symbol | Signature | Notes |
|---|---|---|
| `Material` | `Material(name, eps, source=None)` | `eps` is a scalar, callable(wavelength), `(3,3)` array, or callable returning one. Validated at construction (probe wavelength `1.0`) and on every `epsilon_tensor()` call (Category 5 target 5.1). |
| `Material.is_isotropic` / `.is_diagonal` | properties | Dispatch flags `simulation.py` uses to route to the correct eigensolver. |
| `Material.epsilon_tensor(wavelength)` | method | Always returns a `(3,3)` complex array (see `CONVENTIONS.md`'s tensor ordering), even for isotropic materials. |
| `Material.from_nk(name, wavelengths, n, k, source=None)` | classmethod | Tabulated dispersive data from in-memory arrays. |
| `Material.from_nk_file(name, path, wavelength_unit="um", source=None)` | classmethod | Same, from a refractiveindex.info-style CSV file. |
| `Material.from_refractiveindex_formula_file(name, path, ...)` | classmethod | Parses a refractiveindex.info analytic-formula file; only `"formula 4"` supported today (else `NotImplementedError`). |
| `Material.from_sellmeier(name, B, C, source=None)` | classmethod | Analytic Sellmeier dispersion. |
| `Material.from_cauchy(name, A, B, source=None)` | classmethod | Analytic Cauchy dispersion. |
| `Material.from_lorentz(name, ...)` | classmethod | Lorentz oscillator model (absorbing). |
| `Material.from_drude(name, eps_inf, omega_p_ev, gamma_ev, source=None)` | classmethod | Free-electron Drude model (metals); includes `RAKIC_GOLD`/`RAKIC_SILVER`/`RAKIC_ALUMINUM`/`RAKIC_TITANIUM` published presets (target 5.6). |
| `Material.from_drude_lorentz(name, ...)` | classmethod | Combined Drude + Lorentz oscillators. |
| `Material.from_permittivity_tensor(name, eps3x3, source=None)` | classmethod | Direct anisotropic construction from a `(3,3)` tensor (or callable). |

## `geometry.py` — lattices, shapes, patterns

| Symbol | Signature | Notes |
|---|---|---|
| `Lattice(a, b)` | 2D lattice, basis vectors `a`/`b` as `(x,y)` tuples | Raises `ValueError` for a non-finite or degenerate (zero-area) lattice. |
| `Lattice1D(period)` | 1D (trench/grating) lattice | `period` must be positive and finite. |
| `Shape` (ABC) | — | `Circle`, `Rectangle`, `Ellipse`, `Polygon`, `Slab` all implement `.area`, `.bounding_radius`, and a closed-form `.fourier_transform(kx, ky)` (no raster/FFT — see `theory.md` Stage 1). |
| `Circle(center, radius, material)` | | |
| `Rectangle(center, halfwidth, material)` | `halfwidth` is `(hx, hy)` | |
| `Ellipse(center, halfwidth, material)` | | Anisotropically rescaled `jinc` transform (`design.md` §3b). |
| `Polygon(vertices, material)` | | Closed-form boundary-integral Fourier transform, not raster/FFT. |
| `Slab(center, halfwidth, material)` | 1D-lattice shape | `center` is `(x, 0.0)` by the `Shape` ABC's tuple contract. |
| `Pattern(shapes, background)` | | `.containment_tree()` resolves nested-shape double-counting (S4 subtraction-rule convention). |
| `validate_pattern_fits_lattice(pattern, lattice)` | function | Called automatically from `Simulation.__init__`; raises `ValueError` if a shape could overlap its own periodic image. |

## `geometry_io.py` — JSON pattern import (parser only)

| Symbol | Notes |
|---|---|
| `pattern_from_dict(data)` / `pattern_from_json_string(text)` / `pattern_from_json_file(path)` | Minimal, safe schema — stdlib `json` only, never `eval`/`exec` (Category 4 target 4.6). Not wired into `Simulation`/`Layer` automatically; construct a `Pattern` from the result yourself. |

## `fourier_basis.py` — which Fourier orders are kept

| Symbol | Notes |
|---|---|
| `truncate_fourier_orders(lattice, num_orders, method="circular")` | 2D circular G-vector truncation. `NotImplementedError` for an unrecognized `method`. |
| `truncate_fourier_orders_1d(lattice, num_orders)` | 1D analogue. |

## `fourier_factorization.py` — Toeplitz permittivity construction

See `theory.md` Stage 1 and `design.md` Algorithm 3a before using these directly — which function you want depends on the direct-rule/inverse-rule/numerical-inverse distinction, not just isotropic-vs-tensor.

| Symbol | Notes |
|---|---|
| `pattern_epsilon_hat(pattern, g_vectors, wavelength, inverse=False)` | Scalar (isotropic) Fourier coefficients. |
| `toeplitz_matrix(pattern, g_vectors, wavelength, inverse=False)` | Scalar Toeplitz matrix from the above. |
| `pattern_epsilon_hat_component(pattern, g_vectors, wavelength, i, j)` | Per-tensor-component `(i,j)` Fourier coefficients, for anisotropic patterns. |
| `toeplitz_matrix_component(pattern, g_vectors, wavelength, i, j)` | Per-component Toeplitz matrix. |

## `layer.py` — the stack data model

| Symbol | Notes |
|---|---|
| `Layer(name, thickness, material=None, pattern=None)` | Exactly one of `material`/`pattern` must be given (`ValueError` otherwise). `thickness=math.inf` marks a semi-infinite half-space. |
| `LayerStack(layers, incidence, transmission)` | Prepends/appends the semi-infinite incidence/transmission half-spaces automatically. |
| `LayerEigenmodes` | Container: `q`, `phi`, `kp` (see `theory.md` Stage 2) plus `.diagnostics`. |
| `EigenmodeDiagnostics` | Per-solve conditioning info (`cond(epsilon_hat)`, `cond(phi)`, min pairwise eigenvalue gap). |

## `staircase.py` — tapered/sloped sidewalls

| Symbol | Notes |
|---|---|
| `slice_profile(...)` | General geometry-to-layer-slices interface (Category 4 target 4.7). |
| `staircase_circle_layers(...)` / `staircase_rectangle_layers(...)` / `staircase_slab_layers(...)` | Tapered-sidewall generators built on `slice_profile`. `num_slices < 1` raises `ValueError`. |

## `eigenmodes.py` — per-layer eigenmode solve

See `theory.md` Stage 2 for which function applies to which layer type; this is the highest-risk module in the project (`design.md`), so don't guess — check the theory doc's dispatch table first.

| Symbol | Notes |
|---|---|
| `solve_layer_eigenmodes_uniform(omega, kx, ky, eps)` | Closed-form, uniform isotropic. |
| `solve_layer_eigenmodes_uniform_diagonal(...)` / `_uniform_inplane(...)` | Uniform anisotropic (diagonal / in-plane-coupled tensor). |
| `solve_layer_eigenmodes_1d(...)` | 1D-patterned (trench/grating); raises `ValueError` for any nonzero `ky` (conical mounting, out of scope). |
| `solve_layer_eigenmodes_patterned(...)` | 2D-patterned, isotropic. |
| `solve_layer_eigenmodes_patterned_inplane(...)` | 2D-patterned, anisotropic (diagonal/in-plane only — target 1.5 longitudinal coupling still `NotImplementedError`). |
| `build_kp_matrix(omega, kx, ky, epsilon_inv)` | Shared k-parallel operator builder; accepts scalar or full `(n,n)` `epsilon_inv`. |
| `classify_propagating(q, tol=1e-9)` | Propagating vs. evanescent classification, reusing `_select_q_branch`'s convention. |
| `svd_diagnostics(phi, relative_threshold=1e-6)` | Opt-in singular-value spectrum diagnostic (Category 12 target 12.4). |
| `ILL_CONDITIONED_THRESHOLD` (`1e4`), `DEGENERATE_GAP_THRESHOLD` (`1e-6`) | Module-level thresholds controlling the `WARNING`s in `design.md`'s Failure Contract. |

## `excitation.py` — incident plane wave

| Symbol | Notes |
|---|---|
| `PlaneWaveExcitation(s_amplitude, p_amplitude, theta=0.0, phi=0.0)` | See `CONVENTIONS.md`'s worked polarization-state table (TE/TM/linear/circular/elliptical) for how to set `s_amplitude`/`p_amplitude`. |
| `PlaneWaveExcitation.omega(wavelength)` | Returns `2*pi/wavelength` (natural units, `c=1` — not `s^-1`). |
| `PlaneWaveExcitation.incident_mode_amplitude(...)` | Solves for the mode-amplitude vector `a0` producing the requested zeroth-order `(Ex,Ey)`; `LinAlgError` if `kp @ phi` is exactly singular. |

## `smatrix.py` — S-matrix cascading

Full derivation in [`s_matrix_method.md`](s_matrix_method.md); see `theory.md` Stage 3 for the summary.

| Symbol | Notes |
|---|---|
| `interface_smatrix(modes_l, modes_lp1)` | Single-interface S-matrix. `LinAlgError` (via `scipy.linalg.lu_factor`) at exact grazing incidence (`theta=90 deg`) — see `design.md`'s Failure Contract for why this is a genuine, documented edge case, not a bug. |
| `propagation_smatrix(q, thickness)` | Single-layer propagation S-matrix, `exp(+i*q*thickness)` convention. |
| `star_product(n2, a, b)` | Redheffer star product of two S-matrices. |
| `SMatrixStack(thicknesses, all_modes)` | Cascades a full `LayerStack`'s S-matrices; `ValueError` if `len(thicknesses) != len(all_modes)`. |
| `SMatrixStack.partial_smatrix_up_to(...)` | S-matrix up to an intermediate interface — the basis for `interior_amplitudes`. |
| `interior_amplitudes(...)` | Forward/backward mode amplitudes at an interior interface (Category 9 target 9.3; independently derived, `decisions.md` ADR-015). |

## `fields.py` — power and real-space fields

See `theory.md` Stage 4.

| Symbol | Notes |
|---|---|
| `z_poynting_flux(...)` | Reflected/transmitted power from mode amplitudes. **No `0.5` time-average factor** (harmless in ratios; matters for absolute real-space flux — see `CONVENTIONS.md`). |
| `tangential_e_field(...)` | Tangential `E` from mode amplitudes — `E = phi @ (a-b)` with an index swap/sign flip, **not** the naive-looking `phi @ (a+b)` (that combination is `H`). |
| `modal_field_components(...)` | Full per-order `(Ex,Ey,Ez,Hx,Hy,Hz)`. |
| `propagate_amplitudes(q, z, a_top, b_top)` | Amplitudes at depth `z` within a layer. |
| `reconstruct_field_at_points(...)` | Sums the Fourier series to real-space `(x,y)` points/grids. |
| `save_field_grid_npz(...)` | NumPy `.npz` field export (Phase 7). |

## `polarimetry.py` — Jones/Mueller polarimetry

| Symbol | Notes |
|---|---|
| `decompose_sp(ex, ey, phi, cos_theta_signed)` | Cartesian → s/p decomposition; also reused by `postprocessing/`. |
| `jones_reflection_matrix(sim, wavelength, theta, phi)` | Total (order-summed) Jones reflection matrix. |
| `jones_reflection_matrix_by_order(...)` | Per-diffraction-order Jones matrices. |
| `jones_to_mueller(jones)` | Jones → Mueller matrix conversion. |

## `sweep.py` — parameter sweeps and convergence

| Symbol | Notes |
|---|---|
| `SweepResult` | Typed container for a one-parameter sweep. |
| `sweep_wavelength` / `sweep_theta` / `sweep_phi` / `sweep_polarization` / `sweep_thickness` | Thin wrappers repeating `Simulation.solve()` across one parameter. |
| `harmonic_study(...)` | Reflectance/transmittance vs. `num_orders`, for convergence studies. |
| `find_convergence_index(values, tolerance)` | First index where successive values fall within `tolerance`. |
| `auto_select_num_orders(...)` | Conservative auto-stopping criterion (`decisions.md` ADR-018). |
| `rayleigh_wood_wavelengths(...)` / `avoid_rayleigh_wood_anomalies(...)` | Identify/avoid the diffraction-order Rayleigh threshold where `q=0` (see `design.md`'s Failure Contract). |

## `simulation.py` — top-level orchestration

| Symbol | Notes |
|---|---|
| `Simulation(lattice, layers, num_orders, incidence, transmission)` | Validates every patterned layer via `geometry.validate_pattern_fits_lattice` at construction; caches per-pattern Toeplitz matrices and per-layer eigenmode solves (`decisions.md` ADR-016/ADR-022). |
| `Simulation.solve(excitation, wavelength)` | Returns a `SimulationResult`. Raises `NotImplementedError` for any layer/material with nonzero longitudinal anisotropic coupling (target 1.5). |
| `SimulationResult.reflectance()` / `.transmittance()` | Scalar R/T. |
| `SimulationResult.diffraction_efficiencies()` | Per-order efficiencies. |
| `SimulationResult.complex_amplitudes()` | Raw Cartesian `(Ex,Ey)` per order — deliberately not an s/p basis (target 10.1; see `CONVENTIONS.md`'s note on the p-polarization sign-convention ambiguity). |
| `SimulationResult.diffraction_angles()` | `theta=None` for evanescent orders; `phi` always geometric. |
| `SimulationResult.energy_balance()` | Incident/reflected/transmitted/loss + residual, composed from the above (no new formula). |
| `SimulationResult.layer_absorption()` | Per-layer absorbed power. |
| `SimulationResult.order_classification()` | Propagating/evanescent per order (target 1.8). |

## `vectorized.py` — batched wavelength sweep (opt-in)

| Symbol | Notes |
|---|---|
| `sweep_wavelength_vectorized(...)` | Narrowly scoped: uniform-isotropic-only stacks, `num_orders=1`. Every batched function is a formula-identical re-expression of an already-cited scalar function (`decisions.md` ADR-023) — confirmed to `1e-12` agreement with `sweep.sweep_wavelength`, ~31x faster. Raises if the stack isn't uniform-isotropic (`_require_uniform_isotropic_stack`). |

## `ocd.py` — semiconductor OCD parameterization

| Symbol | Notes |
|---|---|
| `OCDTrapezoidParams` | Validated CD-first (critical-dimension) trapezoid parameters; `.sidewall_angle_deg` property. |
| `trapezoid_trench_layers(...)` | Thin wrapper around `staircase.staircase_slab_layers`. |
| `rounded_rectangle_polygon(...)` | Arc-sampled corner rounding built on `geometry.Polygon`. |

## `config.py` — JSON simulation configuration

| Symbol | Notes |
|---|---|
| `simulation_from_dict(data)` | Builds `(Simulation, PlaneWaveExcitation)` from a dict; only ever constructs objects, never calls `.solve()`. |
| `simulation_from_json_string(text)` / `simulation_from_json_file(path)` | Same, from JSON text/file. |

## `cli.py` — command-line entry point

| Symbol | Notes |
|---|---|
| `main(argv=None)` | `sougata-solver run <config.json>` entry point (`pyproject.toml`'s `[project.scripts]`). Exit codes: `0` success, `1` solver failure, `2` invalid config. |

## `export.py` — sweep serialization

| Symbol | Notes |
|---|---|
| `export_sweep_npz(sweep, path)` / `load_sweep_npz(path)` | Serializes a `SweepResult` to/from a NumPy `.npz` archive; metadata JSON-encoded as a string array (no `allow_pickle=True` needed on load). |

## `plotting.py` — matplotlib visualization (lazy import)

Every function takes plain arrays/already-computed result objects — never a `Simulation` — and returns `(fig, ax)`; none calls `.solve()`.

| Symbol | Notes |
|---|---|
| `plot_unit_cell(pattern, lattice, resolution=200, ax=None)` | |
| `plot_layer_stack(thicknesses, labels, ax=None)` | |
| `plot_structure_3d(layer_stack, lattice, resolution=40, extrusion_length=None, ax=None)` | Category 16 3D preview. |
| `plot_rt_spectrum(wavelengths, reflectance, transmittance=None, metadata=None, ax=None)` | |
| `plot_harmonic_convergence(num_orders_values, values, convergence_index=None, value_label="Reflectance", ax=None)` | |
| `plot_diffraction_orders(diffraction_efficiencies, kind="reflected", ax=None)` | |
| `plot_field_intensity(...)` / `plot_field_phase(...)` / `plot_poynting_vector(...)` | Real-space field/Poynting visualizations. |

## `output_paths.py` — run output folders

| Symbol | Notes |
|---|---|
| `run_output_dir(run_name)` | `outputs/YYYY_MM_DD/HH_MM_SS_<run_name>/`. |
| `write_run_metadata(output_dir, script_path, **params)` | Writes `run_metadata.txt` (explicit UTF-8). |
| `find_latest_output(filename)` | Used by `postprocessing/` scripts to locate a `structures/` run's output without an explicit path. |
