# Tutorials — `sougata_solver`

Targets 18.5 ("thin film"), 18.6 ("grating"), and 18.7 ("via/taper") of
`COMMERCIAL_RCWA_ATOMIC_TARGETS.md` Category 18. Each tutorial below walks
through an existing, already-oracle-validated example script — no new
example code was written for this; the scripts themselves already state
what they demonstrate and what they're checked against (see each script's
own module docstring). New here is the surrounding pedagogy: why the
example is built the way it is, and how to read its output.

Before starting, make sure the environment is set up per
[`GETTING_STARTED.md`](GETTING_STARTED.md). All commands below assume an
activated virtual environment and a working directory of the repo root
(`sougata_solver/`).

**Sample outputs below were captured directly by re-running each script on
2026-08-20** (not remembered or assumed) — you should get the same numbers
(these are deterministic solves, no randomness anywhere in the pipeline).

## 18.5 — Thin film: `structures/thin_film/sio2_on_si_thin_film.py`

**What it builds**: `air (semi-infinite) / SiO2 (500 nm) / Si (semi-infinite)`
— read as `layers=[Layer("SiO2", 500e-9, material=film)]`,
`incidence=air`, `transmission=substrate` in
[`api_reference.md`](api_reference.md)'s `simulation.py` section's
`Simulation` signature. The substrate is passed as `transmission=`, not
appended to `layers`, because transmission media are semi-infinite
half-spaces (`theory.md` Stage 2 doesn't apply a finite thickness to
them at all).

**Why this is the right first example**: it's the one geometry where
`sougata_solver`'s answer can be checked against a completely independent,
closed-form Fresnel/TMM calculation — no Fourier factorization or general
eigendecomposition involved (`num_orders=1`, `theory.md` Stage 1 is a
no-op for a uniform layer). If this doesn't match, nothing downstream can
be trusted either.

**Run it**:

```bash
python structures/thin_film/sio2_on_si_thin_film.py
```

**What you'll see** (elliptical polarization, 40° incidence, 400–800 nm
sweep — script's own `main()` prints every row; first few shown here):

```
Stack: air / SiO2 (500 nm) / semi-infinite Si
Angle: 40 deg; polarization: elliptical
 wavelength (nm)         R
           400.0    0.2901
           401.0    0.2935
           402.0    0.2968
           ...
```

**What to look at**: `R` varies smoothly and shows thin-film interference
fringes as wavelength increases (constructive/destructive interference
between the air/SiO2 and SiO2/Si reflections) — a real, physically
expected pattern, not noise. Try changing `POLARIZATION` (line 50) to
`"s"` or `"p"` and re-running to see the polarization-dependent Fresnel
reflectance at oblique incidence instead.

**How this is validated** (see [`validation_guide.md`](validation_guide.md)
for the full account): `tests/oracles/fresnel.py`, an independent
Fresnel/TMM implementation from Born & Wolf/Macleod, and
`tests/oracles/empy_tmm.py`, transcribed from the vendored `EMpy` — both
cross-checked against this exact `sio2_on_si_thin_film.py` structure in
`tests/test_thin_film_empy_cross_check.py`. This script's own docstring
also names the KLA/Filmetrics Reflectance Calculator as the real-world
target it's built to match.

**Try it yourself**: edit `SUBSTRATE_MATERIAL`/`FILM_MATERIAL` (lines
39/42) to swap in `"Ni"`/`"SiO"`, or `FILM_THICKNESS_M` (line 43) to see
the interference fringe spacing change — thinner films show fewer, wider
fringes across the same wavelength range.

## 18.6 — Grating: `structures/trench/trench_grating.py`

**What it builds**: a 1D-periodic lamellar grating (line/space trench) —
`air / [Si ridge + air groove, 1D-periodic] / air`, using `Lattice1D`,
`Pattern`, and `Slab` (`api_reference.md`'s `geometry.py` section).

**Why this is the right second example**: it's the simplest case that
actually exercises `theory.md`'s full pipeline — Fourier factorization
(Stage 1, since the layer is now patterned) feeding a real eigenmode
solve (Stage 2), not the uniform closed-form shortcut Stage 1's a no-op
for.

**Run it**:

```bash
python structures/trench/trench_grating.py
```

**What you'll see** (TE polarization, normal incidence, 500–1500 nm
sweep; 15 orders per side = 31 total Fourier orders — first several rows
shown):

```
 wavelength (nm)         R         T       R+T
           500.0    0.1646    0.8354    1.0000
           510.0    0.0037    0.9963    1.0000
           520.0    0.1774    0.8226    1.0000
           ...
           580.0    0.8203    0.1797    1.0000
           ...
```

**What to look at**: `R+T` should equal `1.0000` at every wavelength —
this grating has no absorption (both materials are lossless at these
wavelengths), so this is a live energy-conservation check, not just a
printed column. Notice how sharply `R` swings between adjacent rows near
580 nm — that's a genuine diffraction-order threshold (a new order
starting to propagate or a Rayleigh–Wood anomaly, `theory.md`'s Stage 2
note on mode classification), not numerical noise.

**How this is validated**: `tests/oracles/rcwa_1d_gaylord.py`, hand-
transcribed from the vendored `Rigorous-Coupled-Wave-Analysis` project
citing Moharam, Grann, Pommet & Gaylord (1995) — this exact geometry
(period, fill factor, ridge index, thickness) is the Phase 3 system-test
benchmark, cross-checked in `tests/test_1d_grating.py`. A second,
structurally-different oracle (`tests/oracles/rcwa_1d_pyrcwa.py`, from the
vendored `PyRCWA` project) checks the same case independently.

**Try it yourself**: increase `NUM_ORD` (line 43, e.g. to 25) and re-run —
`R`/`T` should barely change if 15 orders was already enough (a quick,
informal convergence check; `sweep.harmonic_study` in
[`api_reference.md`](api_reference.md)'s `sweep.py` section
automates this properly).

## 18.7 — Via/taper: `structures/via/tapered_pillar.py`

**What it builds**: a linearly-tapered square Si pillar array — top side
length `TCD=0.42 µm` narrowing to bottom side length `BCD=0.18 µm` over a
`0.46 µm` height, approximated by a stack of `num_slices` uniform-cross-
section "staircase" layers (`staircase.staircase_rectangle_layers`,
`theory.md` doesn't have its own section for this — see `design.md`
Phase 5 and `phases.md`'s Phase 5 entry).

**Why this is the right third example**: unlike the first two, **no
external oracle exists for a genuinely tapered structure** (per
`staircase.py`'s own docstring) — a discretization scheme's correctness
evidence has to come from somewhere else. This script *is* that evidence:
it sweeps `num_slices` at fixed wavelength/angle and shows the answer
converging as the staircase gets finer, which is the only correctness
argument available for this geometry (also exactly what
`tests/test_staircase.py`'s tapered-rectangle convergence checks pin as a
regression).

**Run it**:

```bash
python structures/via/tapered_pillar.py
```

**What you'll see** (actual output, captured 2026-08-20):

```
num_slices           R           T         R+T
         1    0.263198    0.736802    1.000000
         2    0.701689    0.298311    1.000000
         4    0.707004    0.292996    1.000000
         8    0.701246    0.298754    1.000000
        16    0.700712    0.299288    1.000000
        32    0.700621    0.299379    1.000000
```

**What to look at — this is the whole point of the tutorial**:
`num_slices=1` (a single untapered "average" pillar, not tapered at all)
gives `R=0.263`, a **qualitatively wrong** answer. By `num_slices=8` the
result has already settled to `R≈0.701` and changes by less than `0.001`
per doubling thereafter — that convergence trend, not agreement with any
external reference, is what tells you the staircase approximation is
trustworthy for this geometry at this thickness. If you ever see `R`
still changing by more than ~1% between `num_slices=16` and `32` for a
different geometry, that's a signal to increase `SLICE_COUNTS` (line 64)
before trusting the result, not a bug.

**How this is validated**: convergence-vs-`num_slices` (shown above), plus
a zero-taper regression to the already-oracle-validated single-uniform-
layer result (`tests/test_staircase.py`) — i.e. when `TCD == BCD`
(no taper), this must reduce to Phase 4a's already-validated uniform
pillar answer, which it does.

**Try it yourself**: set `TCD = BCD` (both `0.3e-6`, say) and re-run — you
should see `R`/`T` become essentially independent of `num_slices`
immediately (a true prism has nothing to discretize), confirming the
staircase machinery correctly reduces to the untapered case.
