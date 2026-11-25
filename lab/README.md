# lab/: analysis code

Numerical experiments, literature-graph queries and ingest helpers. Anything
that needs third-party dependencies lives here; `tools/` stays stdlib-only.

The results are written up in `wiki/`; this folder holds the code. Each
experiment script names the wiki page it supports, and each page that rests on
a computation names its script.

## Setup

From the repository root:

```bash
pip install -r lab/requirements.txt
cd lab && python shearlib.py      # selftest, should print "shearlib selftest OK"
```

Tested with Python 3.10–3.13 (CI runs 3.10 and 3.13). `scipy` is required from
`exp03` onward: the design solvers use SLSQP and HiGHS, and the visibility model uses
`scipy.special.j1`. On a Windows console, set `PYTHONIOENCODING=utf-8` (OpenAlex paper
titles are not ASCII). Run the tests from the repository root with `pytest -q`
(`pytest -q -m "not slow"` skips the end-to-end script checks).

## Layout

| file | what it does | supports |
|---|---|---|
| `shearlib.py` | shared maths: shear transfer factor, ODT band geometry, Wiener error, cell phantom; 3D Ewald caps, cap density (binned + closed form), convex design solvers, fusion rules; WPM multiple-scattering propagator and an exact first-Born reference; scalar Mie series, convergent-Born-series Helmholtz solver (`CBS`), multi-slice first Born (`MSB`), full-window lattice sampling, gridding-free column estimator (`ColumnWiener`). Has a selftest. | — |
| `exp01_shear_nulls.py` | transfer-null structure, conditioning, directions vs magnitudes | `wiki/notes/finite-shear-null-experiment.md` |
| `exp02_shear_recovery.py` | Monte Carlo validation, per-frequency error, image demo | same |
| `exp03_ewald_fisher.py` | 3D Fisher maps, closed-form cap density, P1, true lost volume, angular coverage | `wiki/notes/multishear-fusion-design-solution.md` |
| `exp04_shear_design.py` | convex optimal designs for three objectives (P-A, P-D, P-W), conjectures C1–C3, calibration error (E1), chromatic shear (E2) | same |
| `exp05_fusion_validation.py` | phantom validation of Wiener/BLUE/naive fusion against the predicted MSE spectrum | same |
| `exp05b_weighted_average.py` | audit of exp05: precision-weighted averaging equals fusion; the low-K shell | same, and `multishear-fusion-derivation.md` |
| `exp06_multiple_scattering.py` | WPM forward model vs the linear one: null preservation, WPM-vs-Born obliquity, linearisation breakdown, design ranking under five forward models, phase wrapping. Needs more than 6 GB of RAM | `wiki/notes/multiple-scattering-check.md` |
| `exp07_spatially_varying_shear.py` | spatially varying shear Delta(x): null space (chirps), calibrated vs naive estimator, windowed null map, Monte Carlo, 2D CG demo, sample translation as diversity. ~90 s | `wiki/notes/spatially-varying-shear.md` |
| `exp08_colour_shear.py` | wavelength as a shear channel: per-colour Ewald geometry, colour shares, dispersion, mosaic crosstalk. numpy + matplotlib only, ~30 s | `wiki/notes/colour-multiplexed-shear.md` |
| `exp09_full_wave_referee.py` | exact Helmholtz referee (convergent Born series) vs scalar Mie; WPM's linear law on a matched window; WPM / multi-slice Born / exp06's model D against the exact nonlinear excess. ~45 min, ~6 GB | `wiki/notes/full-wave-referee.md` |
| `exp10_exact_gridding.py` | exp06's data reconstructed with no gridding (per-column exact solve) vs nearest-voxel deposition. ~10 min, ~4 GB | `wiki/notes/gridding-and-forward-model.md` |
| `exp11_exact_design.py` | exp06's design comparison with both implementation errors removed: multi-slice Born on a 38.4 um window, both estimators, three contrasts; writes the README figure. ~45 min, ~6 GB | same |
| `realsample_table.py` | arithmetic for the literature table: a*dn/lambda per sample from sourced inputs; `--phantom` gives exp06's max- vs mean-contrast readings. Instant | `wiki/notes/real-samples-on-the-scattering-axis.md` |
| `citation_sweep.py` | OpenAlex forward citation sweep over the project's anchors | `wiki/notes/forward-citation-sweep-2026-07.md` |
| `report.py` | renders `out/citations.json` | same |
| `extract_pdfs.py` | PDF → text for ingests (PyMuPDF text layer; `pip install -r lab/requirements-ingest.txt`) | ingest workflow |
| `out/` | generated figures and data — **gitignored**, regenerate by re-running | — |

Import shared maths from `shearlib`. It is the single source of truth for the
transfer factor and band geometry, so a wrong formula is wrong in exactly one
place. New invariants go into its `selftest()`.

## Established so far

`exp01` + `exp02`, both reproducing exactly on re-run:

- `|H_d(K)| = 2|sin(K d/2)| / d`, verified to machine precision, `→ |iK|` as `d → 0`.
- Nulls at `K = 2πm/d`. At the published experimental shear (0.8 μm) three fall
  inside the measurable band (7.85, 15.71, 23.56 rad/μm), costing 3.56 % of it at
  20 dB (the null frequencies themselves at any SNR, for monochromatic light). At
  the 0.1 μm simulation shear there are none.
- Wiener error `[ρG² + 1]^(-1/2)` validated against Monte Carlo to 3–4 decimals.
- Commensurate shear ratios fail (0.8 + 0.4 keeps a shared null at 15.71);
  incommensurate ones work (0.8 + 0.55 clears the band). exp01/02 normalise by
  1/d, which favours small partner shears; the physical ranking is exp04's.
- In the 2D image test a second magnitude along x gains 1.6× in RMSE over the
  single 0.8 shear, and a second direction (x + y at 0.8) gains 3.4×, because the
  y data also cover most of the x nulls as well as the `K·d = 0` blind line. Both
  together: 8.5×.

`exp03` + `exp04` + `exp05`, all reproducing on re-run (exp04 ~4 min, exp05 ~2.5
min, exp03 ~20 s):

- `J(K) = rho * S_t(K_perp) * M(K)` separates exactly. Shear design is 2D and the
  tomography is a per-voxel weight, so the two can be optimised separately. `M(K)` has a closed form
  (`cap_density_closed`), verified against binned digital caps to sub-percent on
  totals. The missing cone falls out of `|c(K)| = 1` with half-angle `90 - theta`.
- The true lost volume at 0.8 um is 3.02 % of the sampled support (the 1D proxy
  said 3.56 %), split 83.6 % null planes / 16.4 % blind plane. A second direction
  at 0.8 um cuts it to 0.09 %, a second magnitude along x to 0.50 %, both to 0.00 %.
  The 0.1 um shear has no in-band nulls and still loses 2.54 %, all near the blind plane.
- Design rule `d ~ 2*pi/K_max` (0.201 um here); the optimiser picks 0.20 um
  as the best single channel and 0.15 um as the best orthogonal pair.
- Per-frequency (precision-weighted) fusion beats equal-weight averaging of
  reconstructions 21x at equal photons; BLUE diverges for a single shear and is
  fine for an incommensurate pair.
- Optimal designs are not sparse (conjecture C1 is refuted); what makes them
  sparse is falling coherence (C2, confirmed). Avoiding small-integer shear ratios
  (C3) matters for two- and three-channel designs only.
- 180 projections reach only 76.6 % of the reachable support on a 25.6 um field.

`exp06`, reproducing on re-run (~20–30 min, more than 6 GB of RAM):

- The null planes survive multiple scattering exactly. Shearing is linear in the
  measured Rytov field `u_R`, so `FT[shear] = G(K_perp) FT[u_R]` holds whatever the
  scattering; on WPM data at full contrast it holds to 4e-16, which checks the
  implementation (on a periodic grid it is the DFT shift theorem). The null geometry
  behind exp03/exp04 therefore does not depend on the forward model; the size of the
  gains does (below).
- WPM vs. first Born (provisional). For one slice WPM takes the obliquity
  factor at the incident ray (`1/k_iz`) where Born takes it at the scattered one
  (`1/k_z`), checked in the selftest. The larger thick-object gap (~0.67) may be
  partly a window/absorber effect and needs a wider-window re-run.
  [Resolved by exp09: it was the window mismatch.]
- Binned -> continuum K_z with nearest-neighbour deposition (assumption A3) costs a
  factor 6.2 on the best design, against the further 2.7x multiple scattering adds at
  the lowest contrast tested. A bilinear ablation is the next step.
  [Done more strongly in exp10: no gridding at all; 73 % of the factor was rounding.]
- Margins compress, and at one contrast the ordering inverts: fusion beats
  averaging 9.6x (ideal linear), 1.8x (realistic gridding), 1.1x (multiple
  scattering in) at chen1998's `0.23 lambda` bound; by `0.46 lambda` naive
  averaging is marginally ahead, which is what a common bias does.
- Wrapping: shear data skip global unwrapping, but the shear difference has its
  own `< pi` limit across the shear distance, which affects ~5x more pixels at 0.8 um
  than the absolute phase's one-pixel residues do.

`exp08`, reproducing on re-run (~30 s):

- Colour acts mainly as a shear magnitude. B+G+R in x and y scores 0.00103
  against 0.00104 for the mechanical {0.8, 0.55} pair; three colours at the *same*
  shear score the same as one. With all colours forced onto the 0.5 um geometry the
  script reproduces exp03/exp04's rows exactly.
- Blue + green is optimal (0.00089 at shares 0.55/0.45; 0.00090 at equal shares) and
  red gets zero share. Supports of different
  colours are not nested, so score every design on one reference support.
- `J = rho S_t M` does not factorise across colours: each has its own `M_c(K)`.
- At the design-rule shear the idea buys nothing (530 nm alone: 0.00094).

`exp09` + `exp10` + `exp11`, run once each (~45 min, ~10 min, ~45 min); exp11's
`--quick` re-run reproduces its contrast-0.5 table:

- `CBS` (convergent Born series) matches the scalar Mie series for a 1.5 um sphere
  to 0.4 % (normal) / 0.9 % (60 deg); first Born is off by 16 % / 80 %.
- On ONE window, WPM = Born x k_z/k_iz x (slice phase) to 0.3 % at every dz, and
  |WPM/Born| / (k_z/k_iz) = 1.000 +- 0.001 per angle bin. exp06's extra steep-ray
  deficit was its window mismatch.
- Measure on the whole window at the reconstruction lattice (`lattice_samples`):
  a 12.8 um crop moves even CBS's linear part by 10 % (54 % near-axis).
- On near-axis rays (the only ones where the referee is window-converged, 2.8 %),
  MSB's excess is within 4 % / 10 % of the signal at a dn = 0.23 / 0.46 lambda;
  exp06's D is off by 46 % / 85 %, more than the excess itself (36 % / 68 %).
  At 0.92 lambda every cheap model fails.
- `ColumnWiener` (no gridding) on exp06's B data: best design 0.1058 -> 0.0280,
  margins 2.5x -> 8.8x (vs single 0.8) and 1.8x -> 5.9x (fusion vs averaging).
- With MSB as the forward model: 1.33x / 1.15x at 0.23 lambda, 1.09x / 1.02x at
  0.46 lambda, either estimator. exp06's verdict holds; its RMSEs were 7-24 % high.

## Next: an inverse that models multiple scattering

~~exp07, spatially varying shear~~: done (`wiki/notes/spatially-varying-shear.md`).
A calibrated estimator neither gains nor loses at 20 dB (the null is smeared over
`K_m eps` but its area is conserved); ignoring the variation costs 0.21 -> 0.28 at
a 3 % ramp; translating the sample across the shear gradient between two exposures
is a real source of shear diversity (0.208 -> 0.124 at a 10 % ramp).

~~A second candidate, raised by `exp06`: a **Born-series / Lippmann-Schwinger**
forward model.~~ Done in exp09 (convergent Born series as an exact referee) and
exp11 (multi-slice Born as the validated design-level model).

What exp11 changes about priorities: at `a dn = 0.23 lambda` the reconstruction
error is now set by multiple scattering, not by gridding or by the design, so
the item with the largest expected payoff is an **inverse that models multiple
scattering**: `MSB` inside the reconstruction loop, in the spirit of `ma2017`'s
WPM-in-the-loop. `MSB` is cheap (two FFTs per slice) and its linear limit is the
current estimator, so the existing fusion is its first iterate. Two smaller items:
a design-level run of the exact referee on a cluster (540 solves of ~5 min each
at 512^2), and a wider window or a near-to-far-field transform so that steep rays
(tan theta_s > 2) can be adjudicated at all.

## Conventions

- Seed anything stochastic (`SEED` at module top) and say so in the docstring.
- Write artifacts via `shearlib.out_path()` so they land in `out/`.
- Print the numbers that end up in the wiki, so a re-run visibly confirms or
  contradicts the page.
- Keep the "what this does **not** establish" honesty in the wiki page, not
  buried in code comments.

## Gotchas hit already

- **The centring phase must use the DEPOSITED VOXEL's `K_z`, not the ray's exact
  one.** `fftn` of a centred array carries `e^{-i K.r_c}`, and with `z_c = 6.4 um`
  that convention flips by ~pi between neighbouring `K_z` voxels. Using the exact
  `K_z` hands each ray a phase error up to pi/2 and destroys the reconstruction
  (0.519 vs 0.104) while looking exactly like a physics result. Cost `exp06` an
  hour. `cap_geometry` now returns `Kz_vox` for this.
- **Scale absorbing-boundary strength per unit LENGTH, not per step.** Otherwise
  halving `dz` doubles the total absorption and a step-size convergence test
  measures the absorber rather than the propagator — which inverted the trend and
  briefly suggested WPM diverges as `dz -> 0`.
- **A "chain" of shear differences that telescopes is not a measurement.**
  `sum_j [ws((j+1)d) - ws(j d)]` collapses to `ws(md)` identically, so it can only
  ever return zero. The real chain sums one-pixel differences shifted along the
  shear axis: `sum_j roll(ws(dx), -j)`. The selftest now pins both branches.

- **Don't test null recovery on a low-frequency-dominated phantom without a
  bandpass.** Recovering `O` from gradient data is ill-posed as `K → 0` for
  *every* shear (that is gradient integration, which GODT deliberately avoids).
  Un-bandpassed, that effect buries the nulls and inverts the apparent ranking
  — the first attempt at `exp02` concluded the opposite of the truth.
- **Compare forward models on ONE window.** exp06 ran WPM on a 19.2 um window with
  an absorbing edge, cropped it to 12.8 um, and compared it with first Born on a
  12.8 um periodic window. Steep rays walk 10-25 um sideways through a 6 um cell
  at 60 deg, so the two windows disagreed about them, and that looked exactly like
  a forward-model gap (the "~0.67, provisional" of exp06 section 3). On a matched
  window WPM = Born x k_z/k_iz to 0.3 %. Even the exact solver's answer moves by
  tens of percent under a 12.8 um crop; measure on the whole window at the
  reconstruction lattice (`shearlib.lattice_samples`) instead. exp06 also printed
  "absorbing edge outside |x| = 7.9 um", but `WPM` sizes the taper on the FULL
  width, so it starts at 6.14 um, inside the 6.4 um crop.
- **A periodic box stalls the Helmholtz solver.** Without absorption at the sides,
  near-grazing waves and guided modes of the periodic array never reach the z
  layers, and CBS converges only algebraically. Give `CBS` a `frame`. Where frame
  and z layers meet, ADD them: taking the max lowers eps but left a trapped mode in
  small boxes (1500 iterations unconverged vs ~200).
- **Sub-slicing WPM changes the object.** Repeating each slice s times makes a
  z-staircase, whose first Born is (1/s) sum_j e^{-i K_z j dz/s} times the
  rectangle rule, not the rectangle rule. Comparing against the wrong one briefly
  made a converging method look divergent.
- **Interpolate at the positions the FFT assumes.** A trigonometric-interpolation
  operator `sum_k P[k] e^{i k (x + Delta)}` must use x = i*dx, numpy's sample
  positions, not the centred coordinates the shear field is defined on. With
  centred x every frequency picks up (-1)^q and the operator also shifts the image
  by half the field. exp07's first draft did this in 2D; its calibrated estimator
  inverted the same wrong operator and looked fine, and only the naive filter
  (error ~2.25 at every amplitude) gave it away. Check any such operator against
  the Fourier shift for a uniform Delta.
- **Python strings eat LaTeX.** "\times" in a normal string is a TAB plus "imes".
  Write wiki pages with an editor or a raw string.
- Windows consoles default to cp1252. Use the `safe()` helper in `report.py`, or set
  `PYTHONIOENCODING=utf-8`.
- OpenAlex wants a contact address: `OPENALEX_MAILTO=you@example.com`.
- `chen1998` does not resolve in OpenAlex (DOI lookup fails, title fallback
  matches an unrelated paper). It is commented out in `citation_sweep.py`
  rather than left to pollute results; its forward citations are unswept.
