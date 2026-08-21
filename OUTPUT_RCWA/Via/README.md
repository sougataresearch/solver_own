# `OUTPUT_RCWA/Via/` — Lumerical RCWA cross-validation for pillar/via

Unlike `Trench/` and `Thin_Film/Complex_Multi/`, no real Lumerical `.fsp`
file existed for a pillar/via structure to transcribe (per the project
owner's direct choice — see `decisions.md`'s Phase 4a Lumerical
cross-validation ADR, ADR-039). These two structures were instead built
**from scratch in Lumerical to match `sougata_solver`'s own already-shipped
`structures/via/pillar_array.py`/`via_array.py` constants exactly**, the
same approach ADR-034 used for `multistack_composite_grating.py`.

## Workflow — build-only script, then a separate export-only script

Per the project owner's direct preference (they run Lumerical themselves
and don't want a combined build+solve+export script): build and export
are two separate `.lsf` files per structure.

1. Run `build_pillar_structure.lsf` (or `build_via_structure.lsf`) — builds
   the geometry and RCWA solver region, no `run;`, no export.
2. Click **Run** in Lumerical yourself.
3. Run the matching `export_*_grating_power.lsf` script (edit its `fname`
   save path first; delete any existing output file, since `write()`
   appends rather than overwrites).

## Pillar (`build_pillar_structure.lsf` + `export_pillar_grating_power.lsf` ↔ `structures/via/pillar_array.py`)

| Setting | Value | Note |
|---|---|---|
| Lattice | square, period 0.7 µm (x and y) | `PERIOD` |
| Shape | circle, radius 0.18 µm, centered in unit cell | `PILLAR_RADIUS` |
| Shape material | **real dispersive Si** (`NK_FILE/si_KLA.txt`) | not a flat constant — see "Materials" below |
| Background material | `"etch"` (confirmed constant, index=1 / air) | `N_BG` |
| Layer thickness (z-extent) | 0.46 µm | `THICKNESS` |
| Incidence medium | air, semi-infinite | |
| Transmission medium | same dispersive Si as the pillar, **genuinely semi-infinite** | growth substrate — see "Semi-infinite substrate" below; not free-standing (changed from an earlier free-standing air/pillar/air version after reviewing the render) |
| Incidence angle | 0° (normal) | `INCIDENT_ANGLE_DEG` |
| Wavelength range | 0.40–0.80 µm, 401 points | `WAVELENGTHS` |

## Via (`build_via_structure.lsf` + `export_via_grating_power.lsf` ↔ `structures/via/via_array.py`)

| Setting | Value | Note |
|---|---|---|
| Lattice | square, period 0.7 µm (x and y) | `PERIOD` |
| Shape | circle (air-filled hole), radius 0.18 µm, centered in unit cell | `VIA_RADIUS` |
| Background material | same real dispersive Si as the pillar | not a flat constant — see "Materials" below |
| Hole fill material | `"etch"` (confirmed constant, index=1 / air) | `N_VIA` |
| Layer thickness (via depth) | 0.46 µm | `THICKNESS` |
| Incidence medium | air, semi-infinite | |
| Transmission medium | same dispersive Si, **genuinely semi-infinite** | see "Semi-infinite substrate" below |
| Incidence angle | 0° (normal) | |
| Wavelength range | 0.40–0.80 µm, 401 points | |

Polarization: no explicit polarization property is set on the RCWA
object — `grating_power`'s export already returns both `Rs`/`Ts` (s-pol)
and `Rp`/`Tp` (p-pol) regardless, so the solve isn't tied to one
polarization. `sougata_solver`'s own `S_AMPLITUDE=1.0`/`P_AMPLITUDE=0.0`
just means the overlay scripts compare against the `Rs`/`Ts` columns.

## Object structure convention — one Rectangle + one Structure Group

The project owner's standard way of building an RCWA structure in
Lumerical is one plain Rectangle (the background/bulk material) plus one
Structure Group (holding the patterned/etched feature via a construction
script) — confirmed directly by extracting the actual embedded
construction script from `my_trench_0.3.fsp`'s binary (no Lumerical
install needed for this — same technique as `decisions.md` ADR-035's
`.fsp` spot-check). That file's object tree is `::model::Si_slab` (a
plain Rectangle) + `::model::Etch` (a Structure Group whose `setupscript`
property runs `deleteall; for(...) { addrect; ... }`, drawing 32 named
`etch_i` rectangles — the staircase taper).

`build_pillar_structure.lsf`/`build_via_structure.lsf` mirror this
exactly: one plain background Rectangle + one Structure Group whose
construction script draws one `addcircle` child (no loop needed — neither
structure is tapered, so there's only one slice). The geometry `set(...)`
property names (`"x span"`/`"z span"`/etc., with a space, and
`"material"`) are copied verbatim from the real trench file's own script
text, not guessed. The Structure Group's script-*enabling* properties
(`"construction group"`/`"script"`) were a best-effort reconstruction
from the internal storage keys found (`constructionflag`/`setupscript`) —
**confirmed working**: the built object tree matched exactly
(`Air_bg`/`Si_bg`, `Si_substrate`, `Pillar`/`Via` > `pillar_0`/`via_0`,
`RCWA`), no errors.

## Materials — real dispersive Si (`Si_KLA`), named database material

`pillar_array.py`/`via_array.py` both load real dispersive `n,k` data via
`Material.from_nk_file` (`NK_FILE/si_KLA.txt`) — not a flat constant —
matching `structures/trench/tapered_trench.py`'s own convention. Changed
from an earlier flat `n=3.48` constant per the project owner's direct
request.

The Lumerical scripts reference a named database material, `"Si_KLA"`
(the project owner's own confirmed naming), the same way
`my_trench_0.3.fsp`'s own construction script references `"etch"` without
ever defining it inline. **You need `Si_KLA` already imported in your
Material Database** (`si_KLA.txt`'s `n` *and* `k` columns both loaded, as
a Sampled-data material, wavelength unit nm) before running either build
script — it must **not** be the built-in dispersive "Si (Silicon) -
Palik" material, a different dataset that reproduces ADR-036's
already-found-and-fixed mismatch, and must **not** be a plain constant
"Dielectric" type (no `k` at all).

**Real mistake found and fixed this way once already**: an earlier
version of the pillar build had its `Si_substrate`/`pillar_0` objects
still referencing an old constant `Si_n3p48` material rather than the
newly-imported `Si_KLA`, so the solve never touched real absorption data
at all — confirmed by comparing energy conservation on both sides
(`sougata_solver`'s `R+T` ranged 0.46–0.96 as expected from real Si
absorption; Lumerical's `Rs+Ts` was exactly `1.000000` everywhere, the
signature of a lossless material). If a future overlay ever shows
`Rs+Ts=1` across the whole spectrum again, check the object's actual
`material` property first — that was the real cause last time, not a bad
import.

## Semi-infinite substrate (pillar and via)

Per ADR-034's finding: Lumerical infers whether a medium is semi-infinite
from whether the drawn object's z-extent reaches past the RCWA
computation region's own z-boundary — not from an "Interfaces" list
entry. The substrate Rectangle extends to `z min=-1e-6`, comfortably past
the RCWA region's own `z min=-0.05e-6`, the same way the trench/
composite-grating cross-validations confirmed a substrate reads as
genuinely semi-infinite.

## Harmonic order count — a genuine truncation-scheme difference, not directly translatable

`sougata_solver` truncates its 2D Fourier order set **circularly** by
`|k|`-magnitude (`fourier_basis.truncate_fourier_orders`, matching S4's
`gsel.c` convention) — not a square per-axis box. Lumerical's RCWA solver
tab instead exposes a per-axis `"max number ku"`/`"max number kv"` (a
square `(2*ku+1) x (2*kv+1)` truncation for a 2D lattice) — for the 1D
trench case these two schemes coincide exactly (`ADR-036`'s ku=15 matched
`NUM_ORD=31` exactly), but for a genuinely 2D lattice they do not, and no
exact ku/kv-to-`NUM_ORDERS` conversion formula has been independently
confirmed here (an honest gap, per `rules.md` AI Coding Rule 1/5 — not
fabricated). Both build scripts set `ku=kv=17` (1225 total), a generous
over-provision.

**`sougata_solver`'s own `NUM_ORDERS` is still an open question as of
this writing** — see `decisions.md` ADR-040's addendum: after switching
to real dispersive Si, a re-measured convergence check found the
short-wavelength end (~420nm, where Si absorbs strongly) has not settled
even at `num_orders=289`. Don't treat any specific `NUM_ORDERS` value in
`pillar_array.py`/`via_array.py` as final until that's resolved.

## Known feature to expect, not a bug

An earlier, free-standing (air/pillar/air) version of the pillar
structure showed sharp, narrow guided-mode/Fano-type resonances across
the swept spectrum. Adding the Si growth substrate measurably damps that
behavior. Real dispersive Si's strong short-wavelength absorption is a
separate, also-expected effect — expect `R+T` well below 1 (down to
roughly 0.46 in the worst case measured so far) near 400nm, rising toward
~1 near 700-800nm where Si is nearly lossless. This is real physics
(`si_KLA.txt`'s own `k` column), not a bug.
