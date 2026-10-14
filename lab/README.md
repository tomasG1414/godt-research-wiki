# lab/ — analysis code

Numerical experiments, literature-graph queries, and ingest helpers. Anything
that needs third-party dependencies lives here; `tools/` stays stdlib-only.

Findings belong in `wiki/`, not here. Every script names the wiki page it
supports, and every page that rests on a computation names the script. Code
without a wiki page is unfinished work.

## Setup

From the repository root:

```bash
pip install -r lab/requirements.txt
cd lab && python shearlib.py      # selftest, should print "shearlib selftest OK"
```

Tested with Python 3.10–3.13. `scipy` is required from `exp03` onward — the design
solvers use SLSQP and HiGHS, and the visibility model uses `scipy.special.j1`. On a
Windows console, set `PYTHONIOENCODING=utf-8` (OpenAlex paper titles are not ASCII).

## Layout

| file | what it does | supports |
|---|---|---|
| `shearlib.py` | **shared math** — shear transfer factor, ODT band geometry, Wiener error, cell phantom; 3D Ewald caps, cap density (binned + closed form), convex design solvers, fusion rules; WPM multiple-scattering propagator and an exact first-Born reference. Has a selftest. | — |
| `exp01_shear_nulls.py` | transfer-null structure, conditioning, directions vs magnitudes | `wiki/notes/finite-shear-null-experiment.md` |
| `exp02_shear_recovery.py` | Monte Carlo validation, per-frequency error, image demo | same |
| `exp03_ewald_fisher.py` | 3D Fisher maps, closed-form cap density, P1, true lost volume, angular coverage | `wiki/notes/multishear-fusion-design-solution.md` |
| `exp04_shear_design.py` | convex (P-A)/(P-D)/(P-W) optimal designs, C1–C3, calibration (E1), chromatic shear (E2) | same |
| `exp05_fusion_validation.py` | phantom validation of Wiener/BLUE/naive fusion against the predicted MSE spectrum | same |
| `exp05b_weighted_average.py` | audit of exp05: precision-weighted averaging equals fusion; the low-K shell | same, and `multishear-fusion-derivation.md` |
| `exp06_multiple_scattering.py` | WPM forward model vs the linear one: null preservation, WPM-vs-Born obliquity, linearisation breakdown, design ranking under five forward models, phase wrapping | `wiki/notes/multiple-scattering-check.md` |
| `exp08_colour_shear.py` | wavelength as a shear channel: per-colour Ewald geometry, colour shares, dispersion, mosaic crosstalk. numpy + matplotlib only, ~30 s | `wiki/notes/colour-multiplexed-shear.md` |
| `citation_sweep.py` | OpenAlex forward citation sweep over the project's anchors | `wiki/notes/forward-citation-sweep-2026-07.md` |
| `report.py` | renders `out/citations.json` | same |
| `extract_pdfs.py` | PDF → text for ingests (PyMuPDF text layer) | ingest workflow |
| `out/` | generated figures and data — **gitignored**, regenerate by re-running | — |

**Import shared math from `shearlib`.** It is the single source of truth for
the transfer factor and band geometry; if a formula is wrong it should be
wrong in exactly one place. Add new invariants to its `selftest()`.

## Established so far

`exp01` + `exp02`, both reproducing exactly on re-run:

- `|H_d(K)| = 2|sin(K d/2)| / d`, verified to machine precision, `→ |iK|` as `d → 0`.
- Nulls at `K = 2πm/d`. **At the published experimental shear (0.8 μm) three
  fall inside the measurable band** (7.85, 15.71, 23.56 rad/μm) costing 3.56 %
  of it at 20 dB (the null frequencies themselves at any SNR); **at the 0.1 μm simulation shear, none.**
- Wiener error `[ρG² + 1]^(-1/2)` validated against Monte Carlo to 4 decimals.
- Commensurate shear ratios fail (0.8 + 0.4 keeps a shared null at 15.71);
  incommensurate ones work (0.8 + 0.55 clears the band). exp01/02 normalise by
  1/d, which favours small partner shears; the physical ranking is exp04's.
- Magnitudes fill the `m≠0` comb; directions address the `K·d = 0` blind line.
  Both together: 8.5× RMSE improvement over single-shear 0.8.

`exp03` + `exp04` + `exp05`, all reproducing on re-run (exp04 ~4 min, exp05 ~2.5
min, exp03 ~20 s):

- **`J(K) = rho * S_t(K_perp) * M(K)` separates exactly.** Shear design is 2D, the
  tomography is a per-voxel weight, and they never mix. `M(K)` has a closed form
  (`cap_density_closed`), verified against binned digital caps to sub-percent on
  totals. The missing cone falls out of `|c(K)| = 1` with half-angle `90 - theta`.
- **True lost volume at 0.8 um is 3.02 % of the sampled support**, not the 1D
  proxy's 3.56 %, and split 83.6 % null planes / 16.4 % blind plane. The 0.1 um
  shear has no in-band nulls and still loses 2.54 %, all blind plane.
- **Design rule `d ~ 2*pi/K_max`** (0.201 um here); the optimiser picks 0.20 um
  as the best single channel and 0.15 um as the best orthogonal pair.
- **Per-frequency (precision-weighted) fusion beats equal-weight averaging of
  reconstructions 21x** at equal photons; BLUE diverges for a single shear and is
  fine for an incommensurate pair.
- **Sparsity (C1) is refuted**; coherence is what sparsifies a design. C2 and C3
  confirmed, with C3 scoped to small channel sets.
- 180 projections reach only 76.6 % of the reachable support on a 25.6 um field.

`exp06`, reproducing on re-run (~20–30 min):

- **The nulls survive multiple scattering exactly.** Shearing is linear in the
  MEASURED field, so `FT[shear] = G(K_perp) FT[u_R]` is an identity however the
  light got there: verified to 4e-16 on WPM data at full contrast. Every
  geometric conclusion of exp03/exp04 is therefore unconditional.
- **WPM vs. first Born (provisional).** For one slice WPM takes the obliquity
  factor at the incident ray (`1/k_iz`) where Born takes it at the scattered one
  (`1/k_z`), checked in the selftest. The larger thick-object gap (~0.67) may be
  partly a window/absorber effect and needs a wider-window re-run.
- **Binned -> continuum K_z with nearest-neighbour deposition (A3) costs a factor
  6.2** on the best design, against the further 2.7x multiple scattering adds at
  the lowest contrast tested. A bilinear ablation is the next step.
- **Margins compress and at one contrast the ordering inverts**: fusion beats
  averaging 9.6x (ideal linear), 1.8x (realistic gridding), 1.1x (multiple
  scattering in) at chen1998's `0.23 lambda` bound; by `0.46 lambda` naive
  averaging is marginally ahead, which is what a common bias does.
- **Wrapping**: GODT escapes global unwrapping but the shear difference has its
  own `< pi` limit across the shear distance, which affects ~5x more pixels at 0.8 um
  than the absolute phase's one-pixel residues do.

`exp08`, reproducing on re-run (~1 min):

- **Colour acts mainly as a shear magnitude.** B+G+R in x and y scores 0.00103
  against 0.00104 for the mechanical {0.8, 0.55} pair; three colours at the *same*
  shear score the same as one. With all colours forced onto the 0.5 um geometry the
  script reproduces exp03/exp04's rows exactly.
- **Blue + green is optimal (0.00090); red gets zero share.** Supports of different
  colours are not nested, so score every design on one reference support.
- **`J = rho S_t M` does not factorise across colours**: each has its own `M_c(K)`.
- At the design-rule shear the idea buys nothing (530 nm alone: 0.00094).

## Next: exp07, spatially varying shear

The cheapest remaining item. Everything above assumes a
single scalar `Delta` per channel (assumption A5), and the whole design rests on
null planes sitting at fixed, known frequencies. A real two-direction instrument
has `Delta` varying across the field from tube-lens distortion. `exp04` §6 covers
the *random* case (10 nm rms already costs 20 dB at the band edge); a smooth
spatial variation is a different object — deterministic, potentially calibratable,
and one reading says it averages the nulls away for free.

Build it by making `Delta` a function of position and computing a windowed
local-Fourier information map, then comparing against the uniform-shear design.
Reuses `shearlib`'s `G2_vec` and cap machinery unchanged.

A second candidate, raised by `exp06`: a **Born-series / Lippmann-Schwinger**
forward model, whose linear limit *is* first Born by construction. `exp06` had
to splice one (model D) because WPM's linear limit is not Born, and a proper
solver would remove that splice.

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
- Windows consoles default to cp1252. Use the `safe()` helper in `report.py`, or set
  `PYTHONIOENCODING=utf-8`.
- OpenAlex wants a contact address: `OPENALEX_MAILTO=you@example.com`.
- `chen1998` does not resolve in OpenAlex (DOI lookup fails, title fallback
  matches an unrelated paper). It is commented out in `citation_sweep.py`
  rather than left to pollute results; its forward citations are unswept.
