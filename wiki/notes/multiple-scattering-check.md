---
type: note
projects: [godt]
created: 2026-09-22
updated: 2026-10-20
status: current
---

# Multiple scattering: what holds of the design solution

Re-tests [[multishear-fusion-design-solution]] against a forward model it did not
generate. That note's assumption A1 was first-order [[rytov-approximation|Rytov]];
worse, exp03–05 *synthesised* the measurement from the very model they then
inverted, so the only error present was photon noise. Here the data comes from
[[brenner1993]]'s wave-propagation method (WPM) — GODT's own forward engine, and
the multi-slice method that survives the 60° illumination inclination beam
propagation cannot — and the **estimator is unchanged**.

Everything below is **computed**, not sourced: `lab/exp06_multiple_scattering.py`
on shared math in `lab/shearlib.py` (extended with WPM, an exact first-Born
reference, and the phantom, with five new selftest invariants). Nothing from a
source enters except [[s42005-025-02416-3]]'s optical parameters.

## Headline

1. **The transfer nulls are exactly preserved.** This follows from linearity: shearing
   acts on the measured Rytov field *after* all propagation, so
   $\mathcal F[u_R(\mathbf r+\boldsymbol\Delta)-u_R(\mathbf r)]=G(\mathbf K_\perp)\,\mathcal F[u_R]$
   is an identity whatever $u_R$'s relation to $\tilde O$. The **4.2×10⁻¹⁶** check on WPM
   data at full contrast confirms the implementation (on a periodic grid it is the DFT
   shift theorem). Multiple scattering cannot put signal where $G=0$, for monochromatic
   light, a uniform shear and an unwrapped phase difference.
2. **So the *geometry* of the design study is untouched and the *accuracy* is not.**
   Which frequencies are estimable, where the comb sits, why incommensurate
   magnitudes fill it, the rule $\Delta\approx2\pi/K_{\max}$ — all survive exactly.
   What multiple scattering corrupts is the *value* at each $\mathbf K$: a bias no
   design can fuse away.
3. **The design margins collapse, and at one contrast the ordering inverts.**
   Fusion-beats-averaging is 9.6× under the ideal linear model, 1.8× once gridding
   is realistic, and **1.1× once multiple scattering is in**, all measured at
   [[chen1998]]'s bound $a\,\delta n=0.23\lambda$. At $0.46\lambda$ the margin is
   gone entirely and naive averaging is in fact marginally *better* than coherent
   fusion (0.6025 vs 0.6109) — which is what a common bias does, since averaging
   partly cancels it and coherent fusion adds it in phase. Past $\approx0.23\lambda$
   the error budget stops being about shear design.
4. **Nearest-neighbour $K_z$ gridding is the largest *correctable* error in the
   pipeline.** Assumption A3 was assumed away in exp03–05; it is worth a factor
   **6.2** on the best design (0.0171 → 0.1058) before any physics changes, against
   the further ×2.7 multiple scattering adds at the lowest contrast tested
   (0.1058 → 0.2882). Multiple scattering is the bigger absolute error; gridding is
   the bigger multiplier, and unlike the physics it is free to fix.
5. **WPM vs. first Born (provisional).** For a single slice WPM carries the obliquity
   factor of the *incident* ray, $1/k_{iz}$, where Born carries the *scattered* one,
   $1/k_z$ (a known property of multi-slice models). On a thick volume this
   implementation shows a larger residual gap (0.47–0.69) that does not vanish as
   $\delta z\to0$. That residual may come from the simulation window: steep rays walk
   ~11 μm sideways at 60°, comparable to the 19.2 μm window and its absorber. It needs
   re-running on a wider window before being read as physics.
6. **Shear data need no unwrapping step.** Against an ODT baseline that is *not*
   unwrapped, in a purely single-scattering pipeline at $a\,\delta n=0.92\lambda$, the ODT
   benchmark collapses to RMSE 2.57 while the shear designs sit unchanged at 0.11 —
   entirely because the absolute Rytov phase reaches 4.96 rad and wraps at 5.7 % of
   pixels, and the shear difference does not. This is [[guizar-sicairos2011]]'s
   argument demonstrated end-to-end. A comparison against *unwrapped* ODT has not been
   run, so the 24× is an upper bound on the benefit. The other caveat: the shear
   difference carries its own $|\Delta\varphi|<\pi$ limit across the shear distance,
   which grows with $\Delta$ (2.2 % of pixels at 0.8 μm, 4.2 % at 1.6 μm at that
   contrast) — it cost nothing measurable here, but it is what would bite first at
   larger shear.

---

## 1. Five forward models, one estimator

Every column below feeds the **same** Rytov/Ewald estimator — same caps, same
$h=G\,k_m/k_z$, same Wiener/BLUE/naive fusion, same noise realisation at equal
photon budget. Only the data changes, so every difference is forward-model error.

| | data | what it adds |
|---|---|---|
| **A** ideal-linear | $h\,\tilde O$ at the binned voxel | exp05's model, exactly |
| **B** Born-exact | first Born/Rytov, continuum $K_z$, real pupil, real field | realistic gridding |
| **B′** WPM-linearised | WPM at contrast $\varepsilon$, scaled up | multi-slice propagation |
| **C** WPM | the real thing | multiple scattering |
| **D** Born + MS | $d_B+(d_C-d_{B'})$ | multiple scattering *only* |

**D is what makes the question answerable.** B′ turns out to sit at an RMSE floor of
0.635 — the multi-slice-vs-Born gap of §3 — which ruins the reconstruction before
multiple scattering is even asked about, so C−B′ cannot be read at the design level.
Since depositing is linear in the data, subtracting WPM's own linearisation and
adding the difference to the Born baseline cancels that gap exactly while keeping
the genuine nonlinear perturbation. **B→D is then multiple scattering, measured on
a reconstruction that works.**

Setup: [[optical-diffraction-tomography|GODT_SIM]] ($\lambda$ 0.5 μm, $n_0$ 1.33,
NA 1.33, 60° cone), $128^3$ at 0.1 μm isotropic, 180 illuminations snapped to the
propagation lattice (realised inclination 60.16 ± 0.57°) so the incident plane wave
is exactly periodic under WPM's FFTs. WPM runs on a 19.2 μm padded grid with an
absorbing edge and is cropped back to 12.8 μm. Cell phantom, cytoplasm radius
2.98 μm, $\delta n\in[-0.0086,+0.0769]$ at contrast 1. $\rho=100$, aggregates on the
covered support above $|K|=2$ rad/μm, as in exp05.

*(The 180 projections cover 92.1 % of the support the geometry reaches — exp03's
angular-starvation finding, here paid by every model equally.)*

## 2. The nulls are exactly preserved (T1′)

This result is structural. The shear difference is a *transverse linear operation on
the measured Rytov field* $u_R=\ln(u/u_i)$; in the field $u$ itself the interferogram is
bilinear. Whatever nonlinear path the light took to reach the detector, the measurement
is $u_R$, and

$$\mathcal F\big[u_R(\mathbf r+\boldsymbol\Delta)-u_R(\mathbf r)\big](\mathbf K_\perp)
=\big(e^{i\mathbf K_\perp\cdot\boldsymbol\Delta}-1\big)\,\mathcal F[u_R](\mathbf K_\perp).$$

Checked on the WPM field at contrast 1, against the digital shear. With whole-pixel shears on a periodic grid this is the DFT shift theorem, so it confirms the implementation rather than the physics:

| shear | in-band nulls | $\max\lvert\mathcal F[\text{sheared}]-G\,\mathcal F[u_R]\rvert/\max\lvert\mathcal F[u_R]\rvert$ | $\lvert$data$\rvert$ on the null columns |
|---|---|---|---|
| 0.40 μm | 1 | 4.16e-16 | 2.1e-15 |
| 0.80 μm | 3 | 4.16e-16 | 1.1e-15 |
| 1.60 μm | 7 | 4.21e-16 | 5.6e-16 |

**Consequence.** Every geometric conclusion of [[multishear-fusion-design-solution]]
is immune: the null comb at $K=2\pi m/\Delta$, the blind plane
$\mathbf K_\perp\!\cdot\!\boldsymbol\Delta=0$, the 83.6 / 16.4 split of the lost
volume, the incommensurate-ratio rule, and the separation
$J=\rho\,S_{\mathbf t}\,M$ — because $S_{\mathbf t}$ is built from $G$ alone and $M$
from the cap geometry alone, and multiple scattering changes neither. It changes
$\tilde O\mapsto\tilde O+\text{bias}$. For a linear, first-order reconstruction,
estimability is a property of the design and accuracy is not. The scope is monochromatic
light, a uniform shear, no wrapping of the phase difference and a periodic or unbounded
field of view; a nonlinear inverse could in principle draw weak information about the
null planes from data off them.

## 3. WPM vs. first Born (provisional)

> **Caveat.** The thick-volume residual below may be a window/absorber effect rather
> than physics: at 60° the scattered light walks ~11 μm sideways through the object,
> against a 19.2 μm WPM window with an absorbing edge. The single-slice obliquity
> factor is solid (and known for multi-slice models); the rest needs a wider-window
> re-run. Models B′, C and D below inherit this.

Measured on a **band-limited random probe** ($\delta n_{\rm rms}=4.1\times10^{-4}$,
correlation 0.15 μm, smoothly windowed) rather than on the cell phantom, whose
sharp edges alias at the band edge and would be read as a forward-model gap that is
really a sampling one. Against an exact first-Born reference with the same pupil
and the same continuum $K_z$:

| $k_z/k_{iz}$ | rays | median $\lvert$WPM/Born$\rvert$ | ratio to $k_z/k_{iz}$ |
|---|---|---|---|
| 0.244 | 29 | 0.110 | 0.47 |
| 0.408 | 106 | 0.210 | 0.52 |
| 0.614 | 152 | 0.321 | 0.52 |
| 0.808 | 195 | 0.489 | 0.61 |
| 1.003 | 256 | 0.767 | 0.76 |
| 1.637 | 2203 | 1.719 | 1.02 |

The leading term is exact and derivable. For a single slice WPM reduces
analytically to $t(\mathbf r)=\exp[i\,k_z(n(\mathbf r),\mathbf k_{i\perp})\,\delta z]$
— the axial wavenumber evaluated at the **incident** transverse frequency — whose
scattered amplitude carries $1/k_{iz}$, while first Born carries $1/k_z(\mathbf k_s)$.
`shearlib.selftest` checks WPM against that transmission to **1e-12** and the
resulting ratio against $k_z/k_{iz}$ to **1e-3**.

On a thick volume the law only holds where the bulk of the rays are — the
$k_z/k_{iz}=1.64$ bin, 2203 of 2941 rays, lands at 1.02 — and WPM additionally
*under*-propagates everything steeper, by a further ~25 % already at
$k_z/k_{iz}\approx1$ and ~2× by 0.24. So the obliquity factor is the identified
leading term, not the whole story. Overall relative gap 0.69, still 0.47 after
dividing $k_z/k_{iz}$ out — and it is not a step-size artefact:

| $\delta z$ (μm) | slices | relative difference to Born |
|---|---|---|
| 0.1000 | 128 | 0.6922 |
| 0.0500 | 256 | 0.6736 |
| 0.0250 | 512 | 0.6671 |

It converges to a floor near 0.667, not toward zero. (The 0.1 μm step used
throughout §5 is within 4 % of the converged value, so nothing below rests on the
discretisation.) On the **cell phantom** the same comparison gives a raw gap of
0.475, and 0.467 after the obliquity factor — smaller than the probe's only because
the phantom's spectrum is concentrated at low $K$, where the two models agree.

**Two things follow.** First, a caveat that belongs on [[brenner1993]]: multi-slice
forward models and Born/Rytov inverses disagree most at *large scattering angle*,
which is exactly where the high-$K_\perp$ end of a multi-shear design lives.
Second, a practical point for validating any Rytov/Ewald pipeline, including this
project's: when synthetic data come from WPM, the comparison tests the pipeline and
the forward-model gap together, so a forward model whose linear limit matches the
inverse gives a cleaner test. This does not affect GODT's own Shepp–Logan check, which
compares GODT with differentiated ODT on the same WPM data.

## 4. Where the linearisation actually breaks

WPM against **its own linearisation**, so §3's gap cancels and what is left is
multiple scattering plus Rytov breakdown. One illumination.

| contrast | $\delta n_{\max}$ | $a\,\delta n/\lambda$ | nonlin | $\lvert K\rvert<10$ | $\lvert K\rvert>20$ |
|---|---|---|---|---|---|
| 0.25 | 0.0192 | 0.114 | 0.070 | 0.069 | 0.228 |
| 0.50 | 0.0384 | 0.229 | 0.149 | 0.146 | 0.481 |
| 1.00 | 0.0769 | 0.458 | 0.330 | 0.323 | 1.077 |
| 2.00 | 0.1537 | 0.915 | 1.045 | 1.010 | 4.724 |
| 4.00 | 0.3075 | 1.830 | 1.234 | 1.224 | 2.699 |

Two points. While the linearisation still means anything, the nonlinearity grows
slightly faster than linearly in contrast — ratios 2.12, 2.22 for the first two
doublings, as a second-order term should — then 3.17 through the breakdown and
1.18 after it, where the measure saturates because the "error" is already the size
of the signal. Read the last two rows as *broken*, not as a trend.

And **it is concentrated at the band edge**: the $|K|>20$ column runs 3.3× the
$|K|<10$ one across the three rows where the comparison is meaningful (4.7× at
contrast 2, then 2.2× once both have saturated). Multiple scattering damages
precisely the high transverse frequencies that multi-shear design exists to
recover.

[[chen1998]]'s first-Rytov criterion is $a\,\delta n\le0.24\lambda$. The phantom at
contrast 0.5 sits at 0.229 — and §5 finds a reconstruction error of 0.288 there
against chen1998's ~25 % at the bound. Different geometry (3D vs. 2D cylinders) and
metric, so this is consistency in magnitude rather than a precise cross-check.

## 5. The design comparison under each model

RMSE against the band-limited truth, equal photon budget, identical noise.
A, B and B′ are linear models with noise scaled to the signal, so they should be
contrast-independent, and they reproduce to three digits across the sweep — a free
consistency check on the machinery. **One row breaks it and is meant to**: B's ODT
benchmark jumps to 2.5705 at the highest contrast, which is phase wrapping, not
scattering. See below and §6.

**Contrast 0.5** ($a\,\delta n/\lambda=0.229$, at chen1998's bound):

| design | fusion | A ideal | B Born | B′ WPM-lin | C WPM | **D Born+MS** |
|---|---|---|---|---|---|---|
| ODT benchmark | wiener | 0.0221 | 0.1067 | 0.6353 | 0.6825 | **0.2883** |
| 0.1 x | wiener | 0.2172 | 0.2395 | 0.6531 | 0.6967 | 0.3541 |
| 0.8 x | wiener | 0.2453 | 0.2651 | 0.6596 | 0.7018 | 0.3704 |
| 0.8 x | blue | diverges | diverges | diverges | diverges | diverges |
| 0.8 x + 0.8 y | wiener | 0.0408 | 0.1120 | 0.6354 | 0.6824 | 0.2898 |
| 0.8 x + 0.8 y | naive | 0.1804 | 0.2090 | 0.6442 | 0.6857 | 0.3287 |
| {0.8,0.55} x + y | wiener | **0.0171** | **0.1058** | 0.6354 | 0.6827 | **0.2882** |
| {0.8,0.55} x + y | naive | 0.1642 | 0.1952 | 0.6398 | 0.6818 | 0.3197 |

**Across contrast**, for the key rows (the script prints the full tables):

| $a\,\delta n/\lambda$ | | A ideal | B Born | B′ WPM-lin | C WPM | **D Born+MS** |
|---|---|---|---|---|---|---|
| 0.229 | best design | 0.0171 | 0.1058 | 0.6354 | 0.6827 | **0.2882** |
| 0.458 | best design | 0.0172 | 0.1057 | 0.6349 | 0.8526 | **0.6109** |
| 0.915 | best design | 0.0172 | 0.1057 | 0.6340 | 1.7008 | **1.6349** |
| 0.915 | ODT benchmark | 0.0222 | *2.5705* | 0.6339 | 2.5271 | *4.2380* |

Two readings. The multiple-scattering column **D** degrades smoothly and steeply:
0.29 → 0.61 → 1.63, crossing RMSE 1 (the estimate is no better than zero) between
0.46λ and 0.92λ. And the italicised ODT entries at 0.915 are *not* a
multiple-scattering effect at all — they are phase wrapping, which is §6's subject
and which hits the absolute-phase channel only.

**The margins, across models and contrasts:**

| | fusion beats averaging | best pair beats single 0.8 |
|---|---|---|
| A ideal-linear | **9.6×** | 14.3× |
| B Born-exact | 1.8× | 2.5× |
| D Born+MS, contrast 0.5 | 1.1× | 1.3× |
| D Born+MS, contrast 1.0 | 1.0× | 1.0× |
| D Born+MS, contrast 2.0 | 1.1× | 1.3× |

At contrast 1.0 the ordering actually inverts: naive averaging is marginally
*better* than coherent fusion (0.6025 vs 0.6109). That is not a defect in the
fusion rule — it is what happens when a common bias dominates, since averaging
per-channel reconstructions partly cancels the bias while coherent fusion adds it
in phase. (At contrast 2.0 fusion is ahead again by 1.1×; the crossover is a
coincidence of two large numbers, not a trend.) It is the clearest single sign
that the error budget has stopped being about design.

**How to read this.** The exp05 conclusions are not refuted; they are *bounded*. Two
of the three ways the error grows — gridding and the multi-slice/Born gap — are
implementation choices with known fixes ([[pan1983]]'s bilinear interpolation for
the first; a forward model matched to the inverse for the second). The third,
multiple scattering, has a known fix too: [[ma2017]]'s nonlinear WPM forward
operator inside the reconstruction loop. **The design study's numbers are what you
get once those are dealt with, not instead of dealing with them.**

## 6. Wrapping: what the shearing geometry buys, and what it costs

Two different things get measured, on the same multiply-scattered fields.

**`unwrap` / `wrap Δ`** — the two routes to the same shear difference,
$\log[u(\mathbf r+\boldsymbol\Delta)/u(\mathbf r)]$ (what the interferogram
delivers) versus $\log u(\mathbf r+\boldsymbol\Delta)-\log u(\mathbf r)$ (what you
hold if the absolute phase is formed first), disagree exactly where a $2\pi$ jump
separates the two sampled points. At one pixel that is the residue density standard
Rytov ODT must unwrap; at 0.8 and 1.6 μm it is wrapping the shearing construction
never has to see.

**`self Δ`** — whether the shear difference *itself* exceeds $\pi$, found by
chaining one-pixel wrap-safe differences (each far below $\pi$) into the
unambiguous difference over $\Delta$ and comparing with the direct principal value.
This is the one wrapping the interferogram cannot escape, because $\arg[\cdot]$ is a
principal value too.

| $a\,\delta n/\lambda$ | $\varphi_{\rm lin}$ | unwrap@1px | wrap 0.8 | wrap 1.6 | **self 0.8** | **self 1.6** |
|---|---|---|---|---|---|---|
| 0.114 | 0.32 | 0 | 0 | 0 | 0 | 0 |
| 0.229 | 0.64 | 0 | 0 | 0 | 0 | 0 |
| 0.458 | 1.27 | 0 | 0 | 0 | 0 | 0 |
| 0.915 | 2.55 | 0.0045 | 0.0101 | 0.0179 | **0.0218** | **0.0422** |
| 1.830 | 5.10 | 0.0067 | 0.0087 | 0.0084 | **0.0358** | **0.0706** |

### What it is worth, in the reconstruction

The design table of §5 contains the demonstration, in the **B Born-exact** column —
a strictly single-scattering model, so nothing here is about multiple scattering:

| $a\,\delta n/\lambda$ | ODT benchmark | {0.8,0.55} x + y |
|---|---|---|
| 0.229 | 0.1067 | 0.1058 |
| 0.458 | 0.1066 | 0.1057 |
| **0.915** | **2.5705** | **0.1057** |

The ODT row, which is *not* unwrapped here, does not degrade — it *falls over*, by a factor 24, between two
contrasts at which the shear rows do not move at all. The cause is confirmed
directly: the true Born Rytov phase reaches 1.21, 2.44 and **4.96** rad at the
three contrasts, so at the third it exceeds $\pi$ and wraps at 5.7 % of pixels, in
100 % of illuminations. The shear channels never form that quantity.

This is [[guizar-sicairos2011]]'s argument carried through to a reconstruction. It
shows what a gradient measurement avoids, but against unwrapped ODT the margin would be
smaller; that comparison is a natural next step.

### The caveat that remains

The surviving condition is $|\Delta\varphi|<\pi$ **across the shear distance**, not
across a pixel, and unlike the absolute-phase problem it gets *worse* with
$\Delta$ — 2.2 % of pixels at 0.8 μm against 4.2 % at 1.6 μm, at
$a\,\delta n=0.92\lambda$. It cost nothing measurable in §5, and below
$a\,\delta n\approx0.46\lambda$ — which covers the whole regime where the design
margins survive anyway — nothing wraps at all. But it is real, it scales with
$\Delta$, and it is a third independent argument for small shear alongside
$\Delta\approx2\pi/K_{\max}$ and the coherence limit of
[[multishear-fusion-design-solution]] §4.

In short: shear data need no global unwrap, so a violated pixel stays a local error
instead of propagating through an unwrap across the field.

## 7. What this does *not* establish

- **Model D is a construction, not a simulation.** $d_B+(d_C-d_{B'})$ cancels WPM's
  linear gap exactly, but the nonlinear excess it transplants is still computed
  inside WPM and inherits WPM's angular accuracy. Its *magnitude* is corroborated
  independently by chen1998 (§4) and its *contrast scaling* is clean, but a
  Lippmann–Schwinger or convergent-Born-series forward model — whose linear limit is
  Born by construction — would settle it properly, and is a natural follow-up.
- **One phantom, one geometry.** A single cell-like object at one size; $a\,\delta n$
  is swept through contrast only, so size-dependence (which is the *other* half of
  chen1998's criterion) is untested here.
- **WPM at $\delta z=0.1$ μm.** The step-size floor of §3 is real but B′ and C share
  it exactly, so §5's D column is unaffected; §3's absolute gap is not a converged
  number for the steepest rays.
- **Still no spatial-coherence treatment.** This is a coherent monochromatic
  multiple-scattering check. [[ledwig-robles2023]]'s Ewald-*plane* question — whether
  $G$ survives as a multiplicative factor under a genuinely partially coherent
  source — is untouched, and §2's theorem does *not* cover it: that identity is about
  a single coherent realisation of $u_R$.
- **Amplitude is kept.** GODT discards it in practice ([[krauze2016]] prices that at
  0.44 quality-index points); here $u_R$ is complex throughout.
- **E6, spatially varying shear, is still open** and is now the cheapest remaining
  item, as it was before.

## Reproduction

```bash
cd lab
python shearlib.py                     # selftest, includes the WPM invariants
python exp06_multiple_scattering.py    # ~28 min
python exp06_multiple_scattering.py --quick   # ~10 min, 60 illuminations
```

Figures land in `lab/out/` (gitignored): `ms_linearity.png`, `ms_slices.png`.
