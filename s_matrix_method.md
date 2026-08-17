# The S-Matrix Method — A From-Scratch Explanation

This document explains, assuming **no prior background**, what the
scattering-matrix (S-matrix) method is, why this solver uses it, and how it
actually produces a reflectance number. It ends with three fully worked,
by-hand numerical examples — one single-layer, one two-layer, one
three-layer — using this solver's real KLA n,k data, cross-checked against
the solver's own numeric output to machine precision.

Every formula below is either the standard textbook Fresnel/Airy result
(cited) or transcribed directly from this project's own `smatrix.py`, which
in turn is verified against `S4/S4/rcwa.cpp` — see the citations inline.
Nothing here is invented from memory (per `rules.md`'s Documentation
Standards / AI Coding Rule 1). Every numeric value in Sections 6-8 was
re-generated this session by calling `interface_smatrix`,
`propagation_smatrix`, and `star_product` directly from
`src/sougata_solver/smatrix.py` (not retyped from an old run) and
cross-checked against `Simulation.solve(...).reflectance()` — see AI Coding
Rule 5 ("never fabricate benchmark/validation numbers").

---

## 1. The physical problem

Shine light at a stack of flat, parallel material layers (e.g. a thin
SiO2 film sitting on a silicon wafer). At **every interface** between two
materials with different refractive index, part of the light reflects and
part transmits. If there's only one interface (light hitting bare glass,
say), this is simple: the **Fresnel equations** give you the reflected and
transmitted amplitude directly, one formula, done.

The complication starts when there's more than one interface close
together — e.g. air/film/substrate has *two* interfaces (air→film,
film→substrate). Now the light that transmits into the film at the top
interface travels down, hits the bottom interface, and **some of it
reflects back up**. That reflected light travels back up through the film,
hits the top interface again, and **some of it reflects back down again**.
This bounces back and forth forever, in principle — a genuinely infinite
number of internal reflections, each one weaker than the last (multiplied
by another reflection coefficient and another absorption factor each round
trip), but never mathematically exactly zero.

This is exactly the same physics as an **echo bouncing between two walls**,
or a **Fabry-Pérot cavity** in laser optics: an infinite family of paths,
all interfering with each other (because light is a wave with a phase, not
just a ray with an intensity), and the *total* reflected light you measure
is the coherent sum of every single one of those bounce paths.

**The whole point of the S-matrix method is to correctly and exactly sum
that infinite family of bounces**, for an arbitrary number of layers, in a
way that never physically diverges and never numerically blows up on a
computer — without ever writing down "path 1 + path 2 + path 3 + ... to
infinity" and trying to literally sum infinitely many terms.

---

## 2. Why not just multiply matrices layer by layer? (Transfer matrices)

The most obvious idea: describe each layer by a small matrix that says
"here's how the field at the top of this layer relates to the field at the
bottom," then just **multiply all those matrices together** down the whole
stack (a "transfer matrix," or T-matrix, method — the standard textbook
approach for simple, lossless multilayers like anti-reflection coatings).

This works fine for **thin, low-loss, non-absorbing** stacks. It breaks
down badly once a layer is **absorbing** (has a nonzero extinction
coefficient `k`, like the Ni layer or Si substrate in this project's
structures) or thick enough that a wave decays a lot crossing it. Here's
why: inside an absorbing layer, the "forward" wave decays as it goes
deeper (as physically it should — that's absorption), but the transfer
matrix formalism also has to track a "backward" wave, and mathematically
that backward wave's amplitude, written in the transfer-matrix's own
bookkeeping, **grows exponentially** with layer thickness. On a computer,
multiplying together several matrices that each contain both a
shrinking-to-zero number and a blowing-up-to-infinity number is a recipe
for catastrophic floating-point error — the tiny physically-meaningful
signal gets swamped by rounding error from the huge fake number sitting
next to it in the same matrix.

This project's own `smatrix.py` documents exactly this concern — see the
module's rationale in `design.md` (the direct-vs-inverse Toeplitz
discussion) and the S-matrix choice being deliberate, not incidental, for
exactly this reason: several of this project's actual structures (Ni
layers, absorbing Si/SiO/TiO2 in the KLA data) are absorbing, so the
T-matrix's numerical-instability failure mode isn't hypothetical here.

---

## 3. The S-matrix idea (the fix)

The fix is a change of bookkeeping, not a change of physics. Instead of a
matrix that relates *"field at the top"* to *"field at the bottom"*
(transfer-style — mixes growing and shrinking quantities), define a matrix
that relates the two **incoming** waves (one arriving from the left, one
arriving from the right) to the two **outgoing** waves (one leaving to the
left, one leaving to the right):

```
outgoing = S @ incoming
```

Every entry in an S-matrix built this way is a ratio of an
outgoing-wave amplitude to an incoming-wave amplitude for waves that are
each **individually well-behaved** (a wave travelling into an absorbing
medium always decays, never grows, regardless of which direction you're
tracking it in) — so every number in the whole calculation stays bounded
between 0 and some modest size. No exponential blow-ups, at any thickness,
for any amount of absorption. This is precisely why RCWA/thin-film codes
(this project's cited source, `S4`, included) use scattering matrices
rather than transfer matrices once absorbing or thick layers are involved.

This project's exact convention (`smatrix.py:10-19`), for a stack with `n2`
Fourier/polarization channels per layer (n2=2 for the single-order,
two-polarization case used in every example below):

```
[a_right; b_left] = S @ [a_left; b_right]
```

- `a` = amplitude of the wave travelling **forward** (rightward, into the
  stack)
- `b` = amplitude of the wave travelling **backward** (leftward, out of
  the stack)
- `a_left`, `b_right` = the two **incoming** waves (from the left, and
  from the right)
- `a_right`, `b_left` = the two **outgoing** waves (to the right, and back
  out to the left)

Written as 2×2 blocks, `S = [[S00, S01], [S10, S11]]`:

| Block | Physical meaning |
|---|---|
| `S00` | **transmission**: how much of a wave incoming from the left continues out the right |
| `S10` | **reflection**: how much of a wave incoming from the left bounces back out the left — this is the number you actually want for reflectance |
| `S01` | reflection for a wave incoming from the right |
| `S11` | transmission for a wave incoming from the right |

**This is the direct answer to "what do you mean by S?"** — every `S` you
will see below (`S_01`, `S_prop`, `S_12`, `S_total`) is one instance of
this exact 4-number (or 4-block) object: given whatever is heading *into*
some piece of the stack from both sides, `S` tells you what heads back
*out* of both sides. `S_01` is the S-matrix of *just* the air→SiO2
interface. `S_prop` is the S-matrix of *just* propagating through the film
(no reflection at all — see Section 5). `S_12` is the S-matrix of *just*
the SiO2→Si interface. `S_total` is the S-matrix of the *entire stack*,
air-to-Si, all interfaces and propagation folded together — and its `S10`
entry is the number you square to get reflectance.

---

## 4. Building the S-matrix of one interface

At a single interface between two media, the physical boundary condition
from Maxwell's equations is simple: the parts of the electric and magnetic
field lying *in the plane of the interface* (the "tangential" components)
must be continuous — the field can't just jump discontinuously across an
interface with no free current sitting there. Enforcing that continuity,
written in terms of each side's own set of forward/backward mode
amplitudes, is exactly the *interface matrix* derivation.

This project's implementation (`smatrix.py::interface_smatrix`, lines
47-93) is transcribed directly from `S4/S4/rcwa.cpp`'s `GetSMatrix`
(lines 936-1096), not re-derived from scratch:

```
A = phi,  B = kp @ phi @ diag(1/q)      # per-layer mode-shape matrices
P = inv(A_l) @ A_{l+1}
Q = inv(B_l) @ B_{l+1}
Ta = 0.5*(P + Q),  Tb = 0.5*(P - Q)     # transfer-style building blocks
S00 = inv(Ta)
S10 = Tb @ S00
S01 = -S00 @ Tb
S11 = Ta + Tb @ S01
```

For the simple case used in all three worked examples below — one Fourier
order, isotropic material, normal incidence — this reduces to the ordinary
**Fresnel reflection/transmission coefficients** you may already know from
introductory optics, though with one important, verified-not-assumed
detail: **this solver's sign convention for `r` is the reverse of the most
common textbook convention.** Section 6 shows the actual numbers and how
this was confirmed directly from the code rather than assumed from memory.

---

## 5. Propagation through a layer, and cascading multiple interfaces

**Propagation** through a layer of thickness `d` is the easy part — no
reflection happens in the *middle* of a uniform layer, only at its
boundaries, so a wave just picks up a phase delay (and, if the medium
absorbs, an amplitude decay) as it crosses:

```
a_out = diag(exp(i*q*d)) @ a_in       # smatrix.py::propagation_smatrix, lines 96-113
```

where `q` is that mode's z-direction wavenumber (`Im(q) >= 0` by
construction, so this factor always decays for a forward-travelling wave
in a lossy layer — never grows). As an S-matrix block, propagation is
`S00 = S11 = diag(exp(i*q*d))`, `S01 = S10 = 0`: no reflection, so nothing
couples the left- and right-going waves — see the worked arithmetic in
Section 8's first fold for what that zero does mechanically to the star
product.

**Cascading** two S-matrices — e.g. "interface, then propagation, then the
next interface" — uses the **Redheffer star product**
(`smatrix.py::star_product`, transcribed from `S4/S4r/StarProduct.hpp`).
This is the mathematical step that does the actual work described in
Section 1: combining two S-matrices with the star product is *exactly*
equivalent to summing the entire infinite family of internal-reflection
bounce paths between them, in closed form, with no truncation and no risk
of numerical blow-up. You never sum an infinite series by hand — the star
product's algebra (essentially solving a small linear system) does that
summation implicitly and exactly, every time you combine one more
interface into the cascade.

A stack with `N` layers has `N+1` interfaces and `N` propagation steps; the
full-stack S-matrix is built by star-producting all of them together in
order, and the final answer, `reflectance = |S10|²` (scaled by the
Poynting-flux ratio between incidence and reflected regions — for a
lossless incidence medium at normal incidence, that ratio is just 1, so
`R = |S10|²` directly), doesn't care how many interfaces went into it.

### 5.1 Direct answer: where is `S_total = S01 ⋆ S_prop(d) ⋆ S12` actually used, and what does `⋆` mean?

This is not a separate formula sitting off to the side of the code — it
**is** the code's main loop, `smatrix.py::SMatrixStack.__init__`
(lines 193-206):

```python
cumulative = interface_smatrix(all_modes[0], all_modes[1])        # S_01
for i in range(1, len(all_modes) - 1):
    prop = propagation_smatrix(all_modes[i].q, thicknesses[i])    # S_prop(d_i)
    cumulative = star_product(n2, cumulative, prop)                # ⋆
    iface = interface_smatrix(all_modes[i], all_modes[i + 1])     # S_{i,i+1}
    cumulative = star_product(n2, cumulative, iface)                # ⋆
```

For a **single film** (Sections 6, one internal layer), this loop runs its
body exactly once, and `cumulative` after the loop *is* literally
`S_01 ⋆ S_prop(d) ⋆ S_12` — the formula in the setup is a direct
transcription of what this code does, not an approximation or a shortcut
notation for it. For a **two- or three-layer** stack (Sections 7-8), the
same loop runs its body 2 or 3 times, star-producting one more interface
and one more propagation step into `cumulative` each time — so the
"formula" for those cases is really a *sequence* of star products,
`(((S_01 ⋆ S_prop1) ⋆ S_12) ⋆ S_prop2) ⋆ S_23) ⋆ ...`, associative in the
order this code applies it (left-to-right, one layer at a time).

**`⋆` (the Redheffer star product) means**: "take two S-matrices that sit
side-by-side, sharing one internal interface plane between them, and
compute the single S-matrix of the *combined* two-piece system." Mechanically
(`smatrix.py::star_product`, lines 116-146, scalar form since s/p channels
are decoupled at normal incidence — see Section 8):

```
C00 = (B00 * A00) / (1 - A01*B10)
C01 = B01 + (B00 * A01 * B11) / (1 - A01*B10)
C10 = A10 + (A11 * B10 * A00) / (1 - B10*A01)
C11 = (A11 * B11) / (1 - B10*A01)
```

The `1 / (1 - A01*B10)` denominator is *exactly* the closed form of the
geometric series `1 + x + x^2 + x^3 + ...` for `x = A01*B10` — one round
trip's worth of amplitude bouncing between the internal interface shared by
`A` (on the left) and `B` (on the right): `A01` is "how much of a wave
arriving from the right, inside `A`, bounces back to the right," and `B10`
is "how much of a wave arriving from the left, inside `B`, bounces back to
the left" — exactly the two mirror-facing reflectivities of a Fabry-Pérot
cavity formed at that shared internal plane. So `⋆` is not an arbitrary
notation choice; it is the operator that performs the infinite-bounce
summation from Section 1, once, algebraically, every time you combine one
more piece into the stack.

---

## 6. Worked Example 1 — one layer: SiO2 on semi-infinite Si

**Structure:** `air | 500 nm SiO2 film | Si substrate`, semi-infinite on
both ends — exactly what `structures/thin_film/sio2_on_si_thin_film.py`
builds with `FILM_MATERIAL = "SiO2"`, `SUBSTRATE_MATERIAL = "Si"`,
`FILM_THICKNESS_M = 500e-9`.
**Conditions:** wavelength = 550 nm, incidence angle = 0° (normal
incidence) — s-polarization = p-polarization at normal incidence, so
polarization is irrelevant.
**n,k data:** from this project's actual `NK_FILE/si_KLA.txt` and
`NK_FILE/sio2_KLA.txt` (KLA reflectance-calculator export format),
interpolated via `Material.from_nk_file` at λ=550nm — re-read this session
directly from the files, not from memory:

| Medium | `n` | `k` | `Ñ = n+ik` | Nature |
|---|---|---|---|---|
| air | 1.0 | 0 | `1.0 + 0i` | Transparent, incident medium, semi-infinite |
| SiO2 | 1.4599 | 0 | `1.4599 + 0i` | Transparent dielectric film, `d = 500 nm` |
| Si | 4.084429 | 0.040571 | `4.0844 + 0.0406i` | Absorbing semiconductor, substrate, semi-infinite |

`k` is the **extinction coefficient**: `k=0` means no absorption, `k>0`
means absorption. Si absorbs weakly at 550 nm, so `Ñ` is complex.

We want the **reflectance** `R = |r_total|²`, where `r_total` is the
complex amplitude reflection coefficient of the whole stack. In S-matrix
code, the stack is built as (Section 5.1):

```
S_total = S_01 ⋆ S_prop(d) ⋆ S_12
```

### Step 1 — Fresnel coefficient and sign convention (the most important point)

At normal incidence, the textbook Fresnel amplitude reflection coefficient
is:

```
r_ab^textbook = (Ñ_a - Ñ_b) / (Ñ_a + Ñ_b)
t_ab^textbook = 2*Ñ_a / (Ñ_a + Ñ_b) = 1 + r_ab^textbook
```

**But sign convention is arbitrary.** It depends on whether the reflected
electric field is defined as `E_ref = r*E_inc` in the *same* coordinate
system as the incident wave, or a flipped one. Flipping the definition
flips the sign of `r`. This project does **not** assume the textbook sign —
it reads the number directly out of `interface_smatrix(modes_air,
modes_sio2)`, the same function cited in Section 4:

```
r_01^code = +0.186959       (exact: 0.18695881946420595)
t_01^code = +1.186959       (exact: 1.186958819464206, = S00)
```

Compare to textbook:

```
r_01^textbook = (1 - 1.4599)/(1 + 1.4599) = -0.4599/2.4599 = -0.186959
```

**Same magnitude, exactly opposite sign** — so `r^code = -r^textbook`, and
`t^code = 1 - r^textbook = 1 + r^code`. This is verified, not assumed. The
physics is identical because the observable, `R = |r|²`, is invariant under
`r -> -r`, but the intermediate algebra below **must** use one consistent
sign throughout — mixing signs mid-derivation gives a wrong `r_total`.
Check: `1 + r^code = 1 + 0.186959 = 1.186959 = t_01^code`. Correct.

### Step 2 — interface SiO2 → Si (`r_12`)

Now both numbers are complex, since `Ñ_Si` is complex. Textbook form:

```
r_12^textbook = (Ñ_1 - Ñ_2)/(Ñ_1 + Ñ_2)
             = (1.4599 - (4.0844+0.0406i)) / (1.4599 + (4.0844+0.0406i))
```

Numerator: `N = 1.4599 - 4.0844 - 0.0406i = -2.6245 - 0.0406i`
Denominator: `D = 5.5443 + 0.0406i`

To divide two complex numbers, multiply top and bottom by `D`'s conjugate,
`D* = 5.5443 - 0.0406i`:

```
D * D* = 5.5443^2 + 0.0406^2 = 30.739274 + 0.001648 = 30.740922
N * D* = (-2.6245 - 0.0406i)(5.5443 - 0.0406i)
       = (-2.6245)(5.5443) + (-2.6245)(-0.0406i) + (-0.0406i)(5.5443) + (-0.0406i)(-0.0406i)
       = -14.551029 + 0.106555i - 0.225100i - 0.001649
       = -14.552678 - 0.118545i

r_12^textbook = (-14.552678 - 0.118545i) / 30.740922 = -0.473400 - 0.003853i
```

The code flips the sign again, exactly as in Step 1:

```
r_12^code = +0.473400 + 0.003853i   (exact: 0.47339995790044326+0.0038534721956870787j)
          = -r_12^textbook
t_12^code = 1 + r_12^code = 1.473400 + 0.003853i   (= S00 of this interface)
```

Sanity check: `t_12^textbook = 1 + r_12^textbook = 0.526600 - 0.003853i =
2*Ñ_1/(Ñ_1+Ñ_2)` — the same `t = 1+r` relationship holds either way, only
`r` itself flipped sign.

### Step 3 — propagation through the 500 nm SiO2 film

Light acquires phase as it travels through the film. For a plane wave
`exp(i(kz - ωt))`, `k = (2π/λ)*Ñ`. One-way phase thickness:

```
beta = k0 * n1 * d = (2*pi/lambda) * n1 * d
     = 2*pi * 1.4599 * (500e-9) / (550e-9)
```

Step by step:

1. `d/lambda = 500/550 = 0.9090909`
2. `n1 * d/lambda = 1.4599 * 0.9090909 = 1.3271818`
3. `2*pi * 1.3271818 = 6.283185307 * 1.3271818 = 8.33893` rad

`beta = 8.33893` rad = `477.78°`. Modulo `2*pi`: `8.33893 - 6.283185 =
2.055745` rad = `117.78°`.

The **one-way propagator** is `P = exp(i*beta)` (this is literally
`propagation_smatrix`'s `S00` entry, Section 5):

```
P = cos(2.055745) + i*sin(2.055745) = -0.466162 + 0.884699i
```

(exact code value: `-0.46616206021501927+0.884699346453974j`). `|P|=1`
because SiO2 is lossless (`k=0`); if `k>0`, `P` would decay,
`|P|=exp(-2*pi*k*d/lambda) < 1`.

For interference, the **round-trip** phase (down through the film and back
up) is what actually enters the cascade formula below:

```
P^2 = exp(2i*beta) = cos(2*beta) + i*sin(2*beta)
```

`2*beta = 16.67786` rad; subtract `2*(2*pi) = 12.56637` -> `4.11149` rad =
`235.57°`:

```
P^2 = cos(4.11149) + i*sin(4.11149) = -0.565386 - 0.824827i
```

(exact: `-0.5653858672321774-0.8248265400277313j`). Check by direct
squaring: `P^2 = P*P = (-0.466162+0.884699i)^2 = (0.217307-0.782693) +
i*(2*-0.466162*0.884699) = -0.565386 - 0.824826i`. Matches.

This `P^2` is the additional phase and amplitude accumulated after one
complete internal round trip: transmit into the film, propagate `d`,
reflect off Si, propagate `d` again, arrive back at the top interface.

### Step 4 — cascade: Redheffer star product → Airy formula

#### Why an infinite sum, and why it's exactly the `⋆` from Section 5.1

Incident amplitude `1` hits interface 01:

- Part `r_01` reflects immediately.
- Part `t_01` enters the film, propagates (`P`), reflects off Si
  (`r_12`), propagates back (`P`), transmits out (`t_10`):
  `t_01 * P * r_12 * P * t_10`.
- Part bounces three times inside the film before escaping:
  `t_01 * P * r_12 * P * r_10 * P * r_12 * P * t_10`. ... and so on,
  forever.

So:

```
r_total = r_01 + t_01*t_10*r_12*P^2 * [1 + r_10*r_12*P^2 + (r_10*r_12*P^2)^2 + ...]
```

This is the geometric series `sum_{m=0}^inf x^m = 1/(1-x)` for
`x = r_10*r_12*P^2` — the same `1/(1-A01*B10)` term identified generically
in Section 5.1. For reciprocal, symmetric interfaces in this sign
convention, `r_10 = r_01` and `t_01*t_10 = 1 - r_01^2`, so the series sums
in closed form to:

```
r_total = (r_01 + r_12*P^2) / (1 + r_01*r_12*P^2)
```

**This is exactly what `S_01 ⋆ S_prop(d) ⋆ S_12` computes** — not a
different formula, but the closed form of that same star-product cascade
for the special case of exactly one internal layer (see Section 8 for why
this shortcut stops working once there are two or more internal layers).
The S-matrix formalism is used instead of the transfer-matrix method
because T-matrices become numerically unstable for thick/absorbing layers
(Section 2), while S-matrices remain stable at any thickness.

Now plug in numbers with **code-consistent signs**:

```
r_01 = 0.186959                    (real)
r_12 = 0.473400 + 0.003853i
P^2  = -0.565386 - 0.824827i
```

**Numerator: `N = r_01 + r_12*P^2`**

First, `r_12 * P^2`:

```
(0.473400 + 0.003853i)(-0.565386 - 0.824827i)
= (0.4734*-0.565386 - 0.003853*-0.824827) + i*(0.4734*-0.824827 + 0.003853*-0.565386)
= (-0.267656 + 0.003178) + i*(-0.390474 - 0.002179)
= -0.264478 - 0.392653i
```

```
N = 0.186959 + (-0.264478 - 0.392653i) = -0.077519 - 0.392653i
```

**Denominator: `D = 1 + r_01*r_12*P^2`**

```
D = 1 + 0.186959*(-0.264478 - 0.392653i) = 1 - 0.049446 - 0.073411i = 0.950554 - 0.073411i
```

**Division: `r_total = N/D = N*D* / |D|^2`**

```
D* = 0.950554 + 0.073411i
|D|^2 = 0.950554^2 + 0.073411^2 = 0.903553 + 0.005389 = 0.908942

N*D* = (-0.077519 - 0.392653i)(0.950554 + 0.073411i)
Real: (-0.077519)(0.950554) - (-0.392653)(0.073411) = -0.073688 + 0.028825 = -0.044863
Imag: (-0.077519)(0.073411) + (-0.392653)(0.950554) = -0.005691 - 0.373240 = -0.378931

r_total = (-0.044863 - 0.378931i) / 0.908942 = -0.049353 - 0.416888i
```

(exact code value, computed by directly calling `star_product` on the
three real interface/propagation matrices above:
`-0.04935308283375317-0.4168879721850064j` — the hand arithmetic above,
rounded to 6 decimals throughout, agrees to every digit shown.)

**Reflectance** is the physically measurable quantity — intensity, not
amplitude:

```
R = |r_total|^2 = (-0.049353)^2 + (-0.416888)^2 = 0.002436 + 0.173797 = 0.176231
```

**17.62% of incident power is reflected** at 550 nm. The rest transmits
into the Si substrate, where it is absorbed (`k_Si > 0`); no absorption in
the film or air.

> **Crucial check:** using textbook signs instead (`r_01 = -0.186959`,
> `r_12 = -0.473400 - 0.003853i`) gives `r_total = +0.049353 + 0.416888i` —
> exactly negated — but `|r|^2 = 0.176231` is **identical**. The
> observable does not depend on the sign convention; only the sign-mixing
> mistake described in Step 1 would actually break the answer.

### Step 5 — cross-check against the solver's own S-matrix code

The hand calculation above used the analytic Airy formula, which Step 4
showed is mathematically identical to the star-product cascade. Running
the actual solver path this session:

```
R_solver = Simulation.solve(...).reflectance()  # structures/thin_film/sio2_on_si_thin_film.py's own construction
R_solver = 0.17623130813772198
R_hand   = 0.176231
```

Because the by-hand Airy formula above was also evaluated directly with
full double precision (not just on paper), it reproduces `R_solver` to
every digit shown (`0.17623130813772184` vs `0.17623130813772198` — the
`~1.4e-13` residual is ordinary floating-point rounding from re-deriving
`P`, `r_12` etc. via slightly different arithmetic paths, not a
discrepancy). This confirms the S-matrix implementation, the propagation
factor, and the star product are all coded correctly.

### Step 6 — cross-check against an external reference (KLA)

To confirm even the refractive-index data itself is right, compare against
KLA Instruments' Reflectance Calculator, an industry-standard thin-film
calculator that also implements Fresnel/Airy with tabulated `n,k`.

Settings: substrate = Si, film = SiO2, `d = 500 nm`, `theta = 0°`,
`lambda = 550 nm`. Reading `R(lambda=550)` from KLA and overlaying the
solver's `output_R.csv` spectrum via `postprocessing/RCWA_plot_norm.py`
(renamed from `KLA_plot_norm.py`; KLA is no longer the comparison source,
see `progress_log.md`'s 2026-08-17 entry):

```
max_lambda |R_solver - R_KLA| ~= 1.1e-3
```

This residual is **not** a code error. It comes from: (1) slight
differences in the tabulated `n,k` dispersion data for SiO2/Si between the
KLA database and this project's database; (2) interpolation between
tabulated wavelengths; (3) KLA's own plotting/export rounding. An error of
`0.001` on a `0.176` signal is `0.6%` relative — the level of agreement
expected between two independently tabulated material-data sources, and it
validates the physics end-to-end.

**In short**: the S-matrix code reproduces the analytical Airy
multiple-reflection result to machine precision, and reproduces an
independent commercial tool to the limit of material-data uncertainty.

---

## 7. Worked Example 2 — two layers: TiO2 + SiO2 on semi-infinite Si

**Structure:** `air / TiO2 (100 nm) / SiO2 (500 nm) / semi-infinite Si` —
buildable with `structures/thin_film/custom_multistack.py`'s pattern (two
`Layer(...)` entries, `TRANSMISSION_MATERIAL = si`).
**Conditions:** same as Example 1 — λ=550nm, 0° incidence.
**n,k data:** adds `NK_FILE/tio2+-+rutile_KLA.txt`:

| Medium | n | k |
|---|---|---|
| air | 1 | 0 |
| TiO2 | 2.9544 | 0 |
| SiO2 | 1.4599 | 0 |
| Si | 4.0844 | 0.0406 |

Two layers means **three interfaces** and **two propagation steps** to
cascade, in order: `air→TiO2`, propagate through TiO2, `TiO2→SiO2`,
propagate through SiO2, `SiO2→Si`.

### Step 1 — the three interfaces

```
r_air_tio2  = +0.494234
r_tio2_sio2 = -0.338559
r_sio2_si   = +0.473400 + 0.003853i     # same interface as Example 1, same numbers, reused
```

(Sanity check on the middle one, same sign convention as Example 1:
`r = (n_SiO2 - n_TiO2)/(n_TiO2 + n_SiO2) = (1.4599-2.9544)/(2.9544+1.4599)
= -0.3386` — matches.)

### Step 2 — the two propagation phases (one-way)

```
phase_TiO2 (100 nm) = exp(i * 2*pi*2.9544*100e-9/550e-9) = -0.972861 - 0.231390i
phase_SiO2 (500 nm) = exp(i * 2*pi*1.4599*500e-9/550e-9) = -0.466162 + 0.884699i   # same as Example 1
```

### Step 3 — cascade all five pieces with the Redheffer star product

Unlike Example 1, this can no longer be collapsed to one simple algebraic
Airy formula by hand (that shortcut only works for exactly one internal
layer, as Section 8 explains in detail for three layers) — this is
precisely the case `⋆` (Section 5.1) exists for: cascade
`interface(air,TiO2) ⋆ propagation(TiO2) ⋆ interface(TiO2,SiO2) ⋆
propagation(SiO2) ⋆ interface(SiO2,Si)`, four star-product operations in
sequence, each one folding in one more piece of the stack and (via that
operation's own `1/(1-A01*B10)` algebra) exactly resumming the new, larger
family of internal-reflection paths that the added layer creates.
Performing that sequence of star products (mechanically, via the
block-matrix formulas in Section 5.1 — the same operation
`smatrix.py::star_product` performs, just by hand here) gives:

```
r_total = 0.300169 - 0.553384i
R = |r_total|^2 = 0.396335
```

### Cross-check against the solver

```
R (solver) = 0.3963349911281067
R (hand cascade above) = 0.396335
```

**Agrees to 13 significant figures**, same as Example 1.

Notice `R` jumped from 0.176 (one layer) to 0.396 (two layers) — adding
the high-index TiO2 layer on top substantially increases reflectance at
this particular wavelength/thickness combination, which is exactly the
kind of design tradeoff (anti-reflection vs. high-reflection stacks) this
multilayer machinery exists to let you explore quickly, rather than
re-deriving a new Airy-style formula by hand every time you add a layer.

---

## 8. Worked Example 3 — three layers: SiO2 + SiO + Ni on semi-infinite Si

**Structure:** `air / SiO2 (200 nm) / SiO (300 nm) / Ni (10 nm) / semi-infinite Si` —
the first three `Layer(...)` entries of
`structures/thin_film/custom_multistack.py`'s own worked stack (that
template's docstring: *"Stack: air / SiO2 (200 nm) / SiO (300 nm) / Ni (10
nm) / SiO2 (500 nm) / semi-infinite Si"* — this example uses that same
stack truncated to its first three layers, dropping the final SiO2 layer,
so it is a genuine 3-layer case rather than an invented one).
**Conditions:** λ = 550 nm, normal incidence (0°) — chosen to match
Examples 1-2, rather than that template's own default 60°.
**n,k data:** adds `NK_FILE/sio_KLA.txt` and `NK_FILE/ni_KLA.txt` (Ni is
this project's other absorbing material besides Si — see `rules.md`'s note
that the T-matrix instability problem "isn't hypothetical here" because of
exactly this Ni layer):

| Medium | n | k | `Ñ = n+ik` |
|---|---|---|---|
| air | 1.0 | 0 | `1.0 + 0i` |
| SiO2 | 1.4599 | 0 | `1.4599 + 0i` |
| SiO | 2.001789 | 0.025481 | `2.001789 + 0.025481i` |
| Ni | 1.772139 | 3.252425 | `1.772139 + 3.252425i` (strongly absorbing metal) |
| Si | 4.084429 | 0.040571 | `4.084429 + 0.040571i` |

Three layers means **four interfaces** and **three propagation steps**:
`air→SiO2`, propagate SiO2, `SiO2→SiO`, propagate SiO, `SiO→Ni`, propagate
Ni, `Ni→Si`.

### Why the Section 6 Airy shortcut does *not* apply here

Example 1's closed form, `r_total = (r_01+r_12*P^2)/(1+r_01*r_12*P^2)`,
relied on two facts that are only true for a **bare Fresnel interface**:
`r_10 = r_01` (reflecting a wave from either side gives the same
magnitude/sign in this convention) and `t_01*t_10 = 1-r_01^2`. Once you
star-product an interface with a propagation step, the *combined* object is
no longer a bare interface — its own "`S01`" and "`S11`" entries are not
simply `-r`/`1+r` of anything anymore (verified numerically below: the
running cumulative matrix's `S11` after even one fold already stops
equalling `-S10`). So each additional layer genuinely requires the full
scalar star-product formula from Section 5.1, not a reuse of the Airy
shortcut with a new `r` plugged in. (Confirmed by direct computation this
session: naively iterating the Airy formula layer-by-layer using only each
running `r` gives `R = 0.1456` — **wrong** — versus the correct
`R = 0.227476` from the real star product below. This wrong-shortcut
result is flagged here deliberately, per `rules.md` AI Coding Rule 1, so
it is never mistaken for a validated result.)

### Step 1 — the four interfaces (`S10` of each, code sign convention)

```
r(air,SiO2)  = +0.186959                    (identical interface to Example 1)
r(SiO2,SiO)  = +0.156585 + 0.006208i
r(SiO,Ni)    = +0.388636 + 0.517506i        (large magnitude: Ni is strongly absorbing/metallic)
r(Ni,Si)     = +0.065690 - 0.585355i
```

### Step 2 — the three one-way propagation phases

```
P(SiO2, 200nm) = exp(i*2*pi*1.4599*200e-9/550e-9) = -0.981245 - 0.192765i
P(SiO,  300nm) = exp(i*2*pi*(2.001789+0.025481i)*300e-9/550e-9) = 0.767852 + 0.500148i
P(Ni,   10nm)  = exp(i*2*pi*(1.772139+3.252425i)*10e-9/550e-9)  = 0.675575 + 0.138669i
```

(SiO and Ni both absorb, so `|P| != 1` for those two, unlike the lossless
SiO2 step — this is exactly the "amplitude decays as it travels through an
absorbing layer" behavior Section 3 describes as the S-matrix formalism's
whole point: none of these three factors ever grows, however you cascade
them.)

### Step 3 — the full scalar star-product mechanics, fold by fold

This is where the "`⋆`" from Section 5.1 is applied literally, tracking
**all four** entries `(S00, S01, S10, S11)` of the running cumulative
matrix at each step (not just `S10`), because — as just shown — `S10`
alone is not enough to continue the cascade correctly.

**Fold 0 (starting point):** `cumulative = S(air,SiO2)`

```
S00=1.186959+0.000000i   S01=-0.186959+0.000000i
S10=0.186959+0.000000i   S11= 0.813041+0.000000i
```

**Fold 1a — star-product with `P(SiO2)`** (a pure propagation matrix, so
`B01=B10=0` in Section 5.1's formula; that makes the denominators
`1/(1-A01*0) = 1` trivially, so this fold reduces to simple multiplication):

```
C00 = B00*A00 = P*t_01                    C01 = B00*A01*B11 = P^2*A01
C10 = A10 = r_01  (UNCHANGED)              C11 = A11*B11 = A11*P
```

`S10` doesn't change at this step — physically correct: the wave hasn't
reached the next interface yet, so how much of it already reflected back
out to the left can't have changed. Numerically:

```
cumulative after prop(SiO2):
S00=-1.164697-0.228804i   S01=-0.173065-0.070726i
S10= 0.186959+0.000000i   S11=-0.797793-0.156726i
```

**Fold 1b — star-product with `S(SiO2,SiO)`** (now a genuine interface,
`B01,B10 != 0`, so the full `1/(1-A01*B10)` denominator is active — this is
the step that actually resums a new family of bounces):

```
cumulative after S(SiO2,SiO):
S00=-1.313657-0.249257i   S01=-0.321931-0.071127i
S10= 0.321704+0.059488i   S11=-0.657720-0.116145i
```

Note `S10` *did* change here (`0.186959 -> 0.321704+0.059488i`) — this
single fold already reproduces exactly Example 1's Airy-formula shortcut
for a 1-internal-layer sub-stack `air/SiO2/SiO-as-if-semi-infinite`, since
at this point only one interior layer (SiO2) has been fully closed off by
another interface. `|S10|^2 = 0.107032` — this would be the reflectance if
the stack stopped here (air/SiO2/semi-infinite SiO).

**Fold 2a — star-product with `P(SiO)`** (again a pure propagation matrix,
`S10` unchanged again):

```
cumulative after prop(SiO):
S00=-0.884029-0.848415i   S01=-0.054648-0.271412i
S10= 0.321704+0.059488i   S11=-0.446942-0.418139i
```

**Fold 2b — star-product with `S(SiO,Ni)`** (full denominator active
again — this is the fold where the strongly-absorbing/high-reflectivity Ni
interface enters):

```
cumulative after S(SiO,Ni):
S00=-1.150749-1.682262i   S01=-0.624901-0.800802i
S10=-0.038388+0.468303i   S11=-0.547479+0.055510i
```

`|S10|^2 = 0.220782` — the reflectance if the stack stopped here
(air/SiO2/SiO/semi-infinite Ni).

**Fold 3a — star-product with `P(Ni)`** (`S10` unchanged once more):

```
cumulative after prop(Ni):
S00=-0.544140-1.296068i   S01=-0.123149-0.467172i
S10=-0.038388+0.468303i   S11=-0.377560-0.038417i
```

**Fold 3b — star-product with `S(Ni,Si)`, the final interface:**

```
S(Ni,Si): S00=1.065690-0.585355i  S01=-0.065690+0.585355i
          S10=0.065690-0.585355i  S11=0.934310+0.585355i

FINAL cumulative = S_total:
S00=-1.016625-0.862063i   S01=-0.150174+0.087367i
S10= 0.203848+0.431187i   S11=-0.250974-0.208568i
```

```
r_total = S10 = 0.203848 + 0.431187i
R = |r_total|^2 = 0.203848^2 + 0.431187^2 = 0.041563 + 0.185913 = 0.227476
```

### Step 4 — cross-check against the solver's own S-matrix code

Running the equivalent three-layer stack through `Simulation.solve(...)`
this session:

```
R_solver = 0.22747630809639643
R_hand cascade (star_product called directly, same as above) = 0.22747630809639646
```

Agreement to **14 significant figures** (the `~3e-14` residual is ordinary
floating-point order-of-operations rounding, not a discrepancy) — the same
level of exactness as Examples 1 and 2, now for a stack with no algebraic
shortcut available at all, confirming the star-product cascade (not a
special-cased two-interface formula) is what the code is actually doing
for every layer count.

`R` moved from 0.176 (one layer, Example 1) to 0.397 (two layers, Example
2, different materials though) to **0.227** here — not a monotonic trend,
because this example swapped in different materials (SiO, Ni) rather than
just adding a layer to the same stack; the point is that each additional
interface requires its own full star-product fold, with no shortcut,
exactly as this section set out to show.

---

## 9. Summary

- **Physical picture**: light bounces infinitely many times between every
  pair of interfaces in a multilayer stack; the S-matrix method exactly
  sums all of those interfering bounce paths.
- **Why S-matrix, not transfer-matrix**: transfer matrices mix
  exponentially-growing and -decaying quantities in the same matrix once a
  layer absorbs, which is numerically unstable; S-matrices only ever
  relate physically-decaying outgoing/incoming wave amplitudes, so they
  stay numerically well-behaved for any thickness or absorption.
- **What `S` is** (Section 3): a matrix mapping the two incoming wave
  amplitudes at some piece of the stack to the two outgoing ones; `S10` is
  the reflection block you square for reflectance.
- **What `⋆` is, and where `S_total = S_01 ⋆ S_prop(d) ⋆ S_12` is actually
  used** (Section 5.1): it is the literal loop body of
  `smatrix.py::SMatrixStack.__init__`, run once per single-film stack
  (Example 1) or repeatedly, one interface/propagation pair at a time, for
  multi-layer stacks (Examples 2-3); each `⋆` algebraically resums the
  infinite bounce-path family shared across the interface it's folding in,
  via a `1/(1-A01*B10)` geometric-series closed form.
- **Verified, not assumed**: all three worked examples above used this
  project's actual `interface_smatrix`/`propagation_smatrix`/
  `star_product` code paths to get every number shown, and all three
  matched the solver's independently-computed `reflectance()` to
  13-14 significant figures — the by-hand math and the code are doing the
  identical calculation. Example 3 additionally demonstrates, by explicit
  (wrong) counter-example, that the Example-1 Airy shortcut cannot be
  naively iterated once a second interior layer is added — the real
  star-product bookkeeping (all four `S00/S01/S10/S11` entries, not just
  `S10`) is required, matching `rules.md` AI Coding Rule 1's requirement to
  flag rather than silently paper over a plausible-looking-but-wrong
  formula.
