---
type: note
projects: [godt]
created: 2026-11-25
updated: 2026-11-25
status: current
---

# The design margins without exp06's implementation errors

[[multiple-scattering-check]] (exp06) measured how much of the multi-shear design
advantage survives realistic data, and found it collapsing: the best design's lead over
the published single 0.8 μm shear fell from 14.3× (ideal linear model) to 2.5× (realistic
gridding) to 1.3× (multiple scattering at $a\,\delta n=0.23\lambda$). Two of those three
steps were implementation choices rather than physics, and exp06 said so. This note
removes both and measures what is left.

- **Gridding** (exp10, `lab/exp10_exact_gridding.py`): exp06 deposits every measured ray
  at the *nearest* $K_z$ voxel (assumption A3 of [[multishear-fusion-design-problem]]).
  Here the same data go through an estimator with no gridding at all.
- **The forward model** (exp11, `lab/exp11_exact_design.py`): exp06's model D spliced
  WPM's nonlinear excess onto first Born. [[full-wave-referee]] (exp09) checked that
  excess against an exact Helmholtz solver and found it dominated by the simulation
  window and by WPM's own errors. Here a forward model that the referee validates
  replaces it.

Everything is **computed**; nothing is sourced except [[s42005-025-02416-3]]'s optics and,
for the gridding idea, [[pan1983]]'s framing of direct inversion. Both scripts reproduce
exp06's A and B columns to its four printed digits before changing anything, and assert it.

## Headline

1. **Most of exp06's "realistic gridding" loss was the gridding.** With an exact estimator
   the best design's RMSE on continuum first-Born data falls from 0.1058 to 0.0280. The
   A→B factor drops from 6.18 to 1.64, so 73 % of it (in log terms) was $K_z$ rounding.
2. **The margins come back.** On the same data, the best design beats the single 0.8 μm
   shear by **8.8×** instead of 2.5×, and per-frequency fusion beats equal-weight
   averaging by **5.9×** instead of 1.8×. The ideal-model values are 14.3× and 9.6×.
3. **Multiple scattering then takes them away again, gridding or no gridding.** With a
   forward model checked against the exact referee, at the edge of the weak-scattering
   regime ($a\,\delta n=0.23\lambda$) the best design's RMSE is 0.27–0.28, its lead over the
   single shear is 1.3×, and fusion beats averaging by 1.1×. At $0.46\lambda$ both are
   1.0–1.1×. At that point multiple scattering, not the reconstruction or the shear
   design, sets the error.
4. **exp06's multiple-scattering conclusion stands, and its numbers were pessimistic.**
   Its spliced model D gave 1.3× and 1.1× at $0.23\lambda$, the same margins, although
   [[full-wave-referee]] shows D's correction has the wrong structure. Its absolute errors
   were too high by 7 % ($0.23\lambda$) and 24 % ($0.46\lambda$).

## 1. Why nearest-voxel gridding costs so much here

A ray of illumination $i$ measures $\tilde O$ at a continuous $K_z=k_z-k_{iz}$. exp06
assigns it to the nearest voxel of a grid with pitch $\delta k=2\pi/12.8\ \mu$m$^{-1}$
= 0.491 rad/μm. Over exp06's 649 080 rays the rounding is 0.142 rad/μm rms (0.245
maximum, half a voxel). Because $\tilde O(K_z)=\int O(z)e^{-iK_zz}dz$, a shift $\delta K_z$
multiplies the contribution of a slice at height $z$ by $e^{-i\,\delta K_z z}$: at the
cell's edge ($|z|=3$ μm) that is a phase error of 0.43 rad rms, enough to scramble the
high-$|K|$ detail the multi-shear design exists to recover.

## 2. An estimator with no gridding

Along one transverse frequency column $\mathbf K_\perp$ the continuum datum is an exact
discrete-time Fourier transform of the $n_z$ slices $x_k$ of the 2D-transformed object:

$$y_r=h_r\,c_r\sum_k x_k\,e^{-iK_{z,r}(k-n_z/2)\,\delta z},$$

with $h_r$ the channel factor and $c_r$ the centring phase exp06 already applies. The MAP
estimate under exp06's own prior ($\tilde O_j$ independent with variance $S_j$) solves, per
column, the $n_z\times n_z$ system

$$\Big[\sum_r w_r\,a_ra_r^{H}+F^{H}\operatorname{diag}(1/S)\,F\Big]x=\sum_r\frac{h_r^*y_r}{\sigma^2}\,a_r,$$

whose data term is Toeplitz and whose prior term is circulant
(`shearlib.ColumnWiener`). Hermitian partners of a real object fold in as extra rays.
Nothing is approximated: no kernel, no oversampling, no support constraint. With every ray
placed at its voxel's $K_z$ it **is** exp06's per-voxel rule $\text{num}/(J+1/S)$; the
selftest checks this to $10^{-12}$ relative. Cost on $128^3$: about one minute per
design, against seconds for the per-voxel rule.

## 3. Results, contrast 0.5 ($a\,\delta n=0.229\lambda$)

RMSE against the band-limited truth, exp06's data, prior, noise realisation and equal
photon budget:

| design | fusion | A ideal | B nearest voxel | **B exact** |
|---|---|---|---|---|
| ODT benchmark | Wiener | 0.0221 | 0.1067 | **0.0328** |
| 0.1 μm $x$ | Wiener | 0.2172 | 0.2395 | 0.2189 |
| 0.8 μm $x$ | Wiener | 0.2453 | 0.2651 | 0.2464 |
| 0.8 $x$ + 0.8 $y$ | Wiener | 0.0408 | 0.1120 | **0.0468** |
| 0.8 $x$ + 0.8 $y$ | naive average | 0.1804 | 0.2090 | 0.1820 |
| {0.8, 0.55} $x$ + $y$ | Wiener | 0.0171 | 0.1058 | **0.0280** |
| {0.8, 0.55} $x$ + $y$ | naive average | 0.1642 | 0.1952 | 0.1661 |

| | best design | fusion vs averaging | best vs single 0.8 μm |
|---|---|---|---|
| A ideal linear | 0.0171 | 9.6× | 14.3× |
| B, nearest voxel (exp06) | 0.1058 | 1.8× | 2.5× |
| **B, exact estimator** | **0.0280** | **5.9×** | **8.8×** |

The designs that gain most are the ones that rely on high-$|K|$ detail (the fused pairs
and the ODT benchmark); the single-channel rows barely move, because their error is
dominated by their nulls and blind plane, which no estimator can fill.

**What the remaining factor 1.64 is.** It is not shear-specific: the ODT benchmark keeps a
similar factor (0.0221 → 0.0328, 1.48×), and its data are exactly the continuum DTFT the
estimator assumes. The remaining cost is therefore that of sampling at off-grid $K_z$
itself under this prior. One off-grid sample constrains a whole column through the
Dirichlet kernel, so voxels next to the missing cone have to be inferred rather than
measured *(inference; not isolated by a separate experiment)*.

## 4. Multiple scattering with a validated forward model (exp11)

exp06's stage 4 again, everything identical (and asserted so for A and B), with two new
forward models in place of its WPM-based C and D:

- **C**: multi-slice first Born ([[multi-slice-forward-models]]) on a 38.4 μm window,
  measured on the whole window at the reconstruction lattice. Its linear limit is first
  Born by construction, and its measurement matches periodic Born's to 1.7 %. On the rays
  where the exact referee is converged it reproduces the exact excess to 4 % of the signal
  at $0.23\lambda$ and 10 % at $0.46\lambda$ ([[full-wave-referee]] §3). No splice.
- **D**: $B+e_\text{MSB}$, MSB's excess on exactly Born's linear part, the analogue of
  exp06's D.

The 38.4 μm window contains both exp06's 19.2 μm illumination lattice and the 12.8 μm
reconstruction lattice, so nothing is re-snapped. Every model is reconstructed with both
estimators: nearest voxel and exact (x).

**Contrast 0.5** ($a\,\delta n/\lambda=0.229$):

| design | fusion | A | B | **C** | D | exp06 D | Bx | **Cx** | Dx |
|---|---|---|---|---|---|---|---|---|---|
| ODT benchmark | Wiener | 0.0221 | 0.1067 | **0.2690** | 0.2495 | 0.2883 | 0.0328 | **0.2714** | 0.2481 |
| 0.1 $x$ | Wiener | 0.2172 | 0.2395 | 0.3406 | 0.3259 | 0.3541 | 0.2189 | 0.3337 | 0.3175 |
| 0.8 $x$ | Wiener | 0.2453 | 0.2651 | 0.3572 | 0.3435 | 0.3704 | 0.2464 | 0.3663 | 0.3458 |
| 0.8 $x$+0.8 $y$ | Wiener | 0.0408 | 0.1120 | 0.2707 | 0.2515 | 0.2898 | 0.0468 | 0.2838 | 0.2570 |
| 0.8 $x$+0.8 $y$ | naive | 0.1804 | 0.2090 | 0.3161 | 0.2996 | 0.3287 | 0.1820 | 0.3164 | 0.2939 |
| {0.8,0.55} $x$+$y$ | Wiener | 0.0171 | 0.1058 | **0.2686** | 0.2491 | 0.2882 | 0.0280 | **0.2802** | 0.2538 |
| {0.8,0.55} $x$+$y$ | naive | 0.1642 | 0.1952 | 0.3080 | 0.2909 | 0.3197 | 0.1661 | 0.3024 | 0.2816 |

**The margins across contrast** (best design = {0.8, 0.55} $x$+$y$, Wiener):

| $a\,\delta n/\lambda$ | model | best design | vs single 0.8 μm | fusion vs averaging |
|---|---|---|---|---|
| any | A ideal linear | 0.0171 | 14.3× | 9.6× |
| any | B Born, nearest voxel | 0.1058 | 2.5× | 1.8× |
| any | B Born, exact | 0.0280 | 8.8× | 5.9× |
| 0.229 | C MSB, nearest voxel | 0.2686 | 1.33× | 1.15× |
| 0.229 | C MSB, exact | 0.2802 | 1.31× | 1.08× |
| 0.229 | exp06's D | 0.2882 | 1.29× | 1.11× |
| 0.458 | C MSB, nearest voxel | 0.4917 | 1.09× | 1.02× |
| 0.458 | C MSB, exact | 0.5353 | 1.07× | 0.96× |
| 0.458 | exp06's D | 0.6109 | — | 0.99× |
| 0.915 | C MSB, nearest / exact | 0.886 / 1.003 | 1.01× / 1.00× | 0.99× / 0.91× |

Four readings.

1. **Multiple scattering, not gridding, is the floor.** On single-scattering data the
   exact estimator recovers most of what realistic data cost (B → Bx). With multiple
   scattering in, both estimators land near 0.27 at $0.23\lambda$ and near 0.5 at
   $0.46\lambda$, and every design lands with them, so the design margins go.
2. **The exact estimator is slightly worse under multiple scattering** (0.2802 against
   0.2686, and 0.5353 against 0.4917). The data are now biased, and an estimator that fits
   them more faithfully fits the bias more faithfully. Nearest-voxel rounding scrambles the
   phase of high-$|K|$ terms and so partly averages the bias away *(inference; not
   isolated)*.
3. **exp06's verdict survives a better forward model.** exp06's D and the validated C give
   the same margins at $0.23\lambda$ (1.3×, 1.1×) and the same tie at $0.46\lambda$,
   although D's excess points the wrong way on the rays the referee can judge. The
   design-level verdict depends on how big the bias is, which D estimated about right.
   D's absolute RMSEs were too high: 0.2882 against 0.2686 at $0.23\lambda$, and 0.6109
   against 0.4917 at $0.46\lambda$.
4. **At $0.92\lambda$ nothing is validated.** MSB fails against the referee there
   ([[full-wave-referee]] §3), so the last row is MSB's answer, not a checked one. The
   ODT benchmark's absolute-phase wrapping (exp06 §6) is visible again at that contrast
   (B: 2.57).

**For the project's design question** this narrows where multi-shear GODT pays off: the
14× and 9.6× of the linear model are mostly recoverable on single-scattering data if the
reconstruction does not grid, and mostly lost to multiple scattering by $a\,\delta n\approx0.25\lambda$
unless the inverse models it. For thin or weakly scattering samples the design study
applies almost undiminished. For a 6 μm cell at $\delta n\approx0.04$ it does not, without a
nonlinear inverse such as [[ma2017]]'s *(inference about real samples; their
$a\,\delta n/\lambda$ is not measured here)*.

## 5. What this does *not* establish

- **One phantom, one geometry, one noise realisation**, exactly as exp06: the comparison is
  paired, so differences between columns are not noise, but the absolute numbers belong to
  this cell phantom.
- **The prior is exp06's** shell-averaged spectrum of the true object, which a real
  reconstruction does not have. Both estimators use the same prior, so the comparison is
  fair; the absolute RMSEs are optimistic.
- **The exact estimator is the MAP estimate for this linear model**, not a regulariser of
  the kind practical ODT uses (TV, non-negativity; [[lim2015]], [[krauze2016]]). It shows
  what the data support, not what a given practical pipeline will reach.
- **Cheaper gridding** ([[pan1983]]'s bilinear interpolation after zero-padding, or a
  NUFFT) was not tested. The exact solve bounds what any of them can recover.
- **The multiple-scattering model is validated narrowly**: against the exact referee on
  one illumination, on near-axis rays, up to $a\,\delta n=0.46\lambda$
  ([[full-wave-referee]] §3). Its steep rays are its own answer, and its $0.92\lambda$ row is
  not validated at all. A design-level run of the exact referee (540 solves of several
  minutes each) would close this.

## Reproduction

```bash
cd lab
python exp10_exact_gridding.py    # ~10 min, ~4 GB
python exp11_exact_design.py      # ~45 min, ~6 GB; --quick: contrast 0.5 only
```

Both print every number above, and exp11 writes `lab/out/realistic_margins.png`, the
figure in the README.
