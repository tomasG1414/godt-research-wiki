---
type: note
projects: [godt]
created: 2026-11-25
updated: 2026-11-25
status: current
---

# An exact referee for the forward models

[[multiple-scattering-check]] (exp06) left two things provisional. Its forward model, the
[[brenner1993]] wave-propagation method (WPM), disagreed with first Born by about 50 % even
in the weak-scattering limit, with an unexplained extra deficit at steep angles. And its
multiple-scattering column, model D, transplanted WPM's nonlinear excess onto first Born
with no way to check that excess. This note adds an exact referee and settles both.

The referee is a [[convergent-born-series]] (CBS) solver of the full scalar Helmholtz
equation, with no weak-scattering, forward-only, paraxial or thin-slice approximation.
Everything here is **computed** by `lab/exp09_full_wave_referee.py` on `lab/shearlib.py`
(`CBS`, `MSB`, `mie_scalar`, `lattice_samples`, all covered by its selftest). Nothing is
sourced except [[s42005-025-02416-3]]'s optics and the CBS iteration itself
(Osnabrugge et al. 2016, cited on the concept page).

## Headline

1. **The referee is exact where it can be checked.** Against the exact scalar Mie series
   for a 1.5 μm sphere it is within 0.4 % at normal incidence and 0.9 % at 60°, for
   spheres where first Born is off by 16 % and 80 % (§1).
2. **exp06's provisional WPM–Born gap was the simulation window plus a known factor.**
   With WPM and Born on the same window, WPM's linear limit is Born $\times\,k_z/k_{iz}$
   to 0.1 % in every angle bin, at every step size. The extra steep-angle deficit came
   from comparing a cropped, absorbing-edge window with a periodic one (§2).
3. **Model D got the size of the multiple-scattering correction roughly right, not its
   structure.** On the rays where the referee itself is window-converged, D's error is
   *larger* than the excess it stands in for: 0.46 against 0.36 of the measured signal at
   $a\,\delta n=0.23\lambda$. Modelling no multiple scattering at all would have been closer
   there. Its design-level margins nonetheless match a validated model's
   ([[gridding-and-forward-model]] §4), because that verdict depends on the bias's size (§3).
4. **Multi-slice first Born is the better cheap model, up to $0.46\lambda$.** On the same
   rays it reproduces the exact excess to 4 % of the signal at $0.23\lambda$ and 10 % at
   $0.46\lambda$, and its full-window measurement matches first Born's to 1–2 %. At
   $0.92\lambda$ every cheap model fails, MSB included (§3).
5. **A 12.8 μm crop is not a neutral way to measure.** It moves even the exact solution's
   linear part by 10 % overall and 54 % on near-axis rays. Measuring on the whole window at
   the reconstruction lattice converges on periodic first Born instead (§3a).

## 1. The referee against the exact Mie series

A sphere of radius 1.5 μm, index 1.34 or 1.38 in $n_0=1.33$, at $\lambda=0.5$ μm, on a
25.6 μm window with an absorbing frame outside $|x|=8.2$ μm. The far-field amplitude comes
from the angular spectrum of the scattered field, $f=k_zA(\mathbf k_\perp)/(2\pi i)$, and is
compared with `mie_scalar`, the scalar partial-wave series. That series is checked
separately against the optical theorem ($10^{-4}$) and against first Born in the weak limit.
Only rays whose sideways walk to the measurement plane stays clear of the frame are scored.
Relative error, all clear rays:

| sphere | phase delay | incidence | **CBS** | MSB | WPM | first Born |
|---|---|---|---|---|---|---|
| 1.34 | 0.38 rad | normal | **0.0038** | 0.0076 | 0.0315 | 0.1553 |
| 1.34 | | 60° | **0.0093** | 0.0191 | 0.1458 | 0.1536 |
| 1.38 | 1.88 rad | normal | **0.0039** | 0.0420 | 0.0351 | 0.7993 |
| 1.38 | | 60° | **0.0094** | 0.1288 | 0.1554 | 0.7899 |

The forward amplitude agrees to 0.1–0.2 % at normal incidence (1.38: Mie
$16.782+14.211i$ μm, CBS $16.806+14.170i$, first Born $24.07$) and to 0.5–0.6 % at 60°. CBS's error grows toward the
edge of the scored range (about 6 % at 40–55°), where MSB's grows in step. That suggests the
measurement geometry rather than the solver *(inference)*. WPM's error concentrates at
20–55°, up to 55 % at normal incidence and 95 % at 60°, as the obliquity law predicts (§2).

## 2. WPM's linear limit, on a window matched to Born's

exp06 compared WPM, run on a 19.2 μm window with an absorbing edge and cropped to 12.8 μm,
against first Born evaluated on a 12.8 μm periodic window. Steep scattered light walks
10–25 μm sideways through a 6 μm cell at 60° illumination, so the two windows disagree
about exactly the rays in question. Here the same probe as exp06 (band-limited,
smoothly windowed, $\delta n$ amplitude $2\times10^{-3}$) is run both ways.
$|\text{WPM}/\text{Born}|$ per obliquity bin, and divided by $k_z/k_{iz}$:

| $k_z/k_{iz}$ | rays | exp06's windows | ÷ obliquity | one window | ÷ obliquity |
|---|---|---|---|---|---|
| 0.262 | 55 | 0.0735 | 0.3265 | 0.2128 | **1.0012** |
| 0.414 | 146 | 0.1764 | 0.4182 | 0.4092 | **1.0001** |
| 0.609 | 222 | 0.3117 | 0.5203 | 0.6128 | **1.0005** |
| 0.805 | 291 | 0.4620 | 0.5707 | 0.8053 | **1.0003** |
| 1.001 | 376 | 0.7374 | 0.7370 | 0.9970 | **0.9994** |
| 1.594 | 2473 | 1.6305 | 1.0022 | 1.5586 | **1.0004** |

(The matched column uses amplitude $2\times10^{-4}$; the next paragraph says why.) The
complex law, including where WPM places each thin screen ($s$ sub-slices per voxel):

$$\text{WPM}=\text{Born}\times\frac{k_z}{k_{iz}}\times\frac1s\sum_{j=1}^{s}e^{-iK_z\,j\,\delta z/s}$$

| $\delta z$ (μm) | amplitude $2\times10^{-4}$ | amplitude $2\times10^{-3}$ |
|---|---|---|
| 0.1 | 0.0034 | 0.0400 |
| 0.05 | 0.0033 | 0.0508 |
| 0.025 | 0.0033 | 0.0744 |
| 0.0125 | 0.0033 | 0.1157 |

At the weaker amplitude the residual is second order and does not move with $\delta z$: the
law is exact. At exp06's probe amplitude the steep rays are already nonlinear, and that
nonlinear part grows as the slices get thinner. **This retires exp06 §3's provisional
reading.** The single-slice obliquity factor was the whole linear story; the further
"25 %–2×" steep-angle deficit and the gap's floor near 0.67 were the windows and the
probe's own nonlinearity. See [[multi-slice-forward-models]].

**An audit finding on exp06's printout.** exp06 printed "absorbing edge outside
$|x|$ = 7.9 μm", but `WPM` sizes its taper on the full width: 18 % of 19.2 μm, so the
absorber starts at 6.14 μm, inside the 6.4 μm crop. The label is now computed
(`lab/exp06_multiple_scattering.py`); no number changes.

## 3. What model D transplants, against the exact excess

exp06's cell phantom (cytoplasm radius 2.98 μm), one illumination at 61.78°, on a lattice
shared by every window below; 3593 collected rays. For each model the measured Rytov
spectrum $d(c)$ on the pupil is taken two ways:

- **full**: the FFT of the whole simulation window at the 12.8 μm reconstruction
  lattice's frequencies (`shearlib.lattice_samples`). This is the infinite-field-of-view
  value, the measurement periodic first Born (exp06's model B) makes;
- **crop**: the FFT of the central 12.8 μm, which is how exp06 measured WPM.

$e(c)=d(c)-(c/\varepsilon)\,d(\varepsilon)$ is the nonlinear excess, so exp06's
$D=B+e(\text{WPM, 19.2 μm, crop})$. The referee is CBS on a 51.2 μm window (absorbing frame
outside $|x|=20.5$ μm), measured full. Rays are split by their angle to the axis: near-axis
($\tan\theta_s<0.8$, 1433 rays, whose sideways walk stays clear of any frame) and steep
($\tan\theta_s>2$, 676 rays).

### (a) Linear parts against periodic first Born, contrast 0.5

| model | all rays | near-axis | steep |
|---|---|---|---|
| CBS, 51.2 μm, full | 0.041 | 0.015 | 0.139 |
| CBS, 38.4 μm, full | 0.061 | 0.015 | 0.207 |
| MSB, 51.2 μm, full | 0.011 | 0.016 | 0.038 |
| MSB, 38.4 μm, full | 0.017 | 0.016 | 0.057 |
| WPM, 51.2 μm, full | 0.095 | 1.313 | 0.242 |
| CBS, 51.2 μm, **crop** | 0.104 | 0.544 | 0.346 |
| WPM, 19.2 μm, **crop** (exp06) | 0.231 | 1.362 | 0.535 |

Measured on the whole window, both MSB and the exact solver converge on Born as the window
grows. A crop does not converge, even for the exact solver, and its error is worst on the
near-axis rays. Those rays carry the weak high-$|K|$ signal (at 60° illumination they are
the back-toward-the-axis scattering), so clipping the strong forward lobe, which walks
about 11 μm sideways, leaks into them *(inference about the mechanism)*. WPM's near-axis
1.31 is its obliquity factor ($k_z/k_{iz}\approx2$ there).

### (b) The nonlinear excess against the referee

$\|e-e_\text{ref}\|/\|e_\text{ref}\|$; the first row of each block is the exact excess's size,
$\|e_\text{ref}\|/\|\text{lin}_\text{ref}\|$.

| contrast ($a\,\delta n/\lambda$) | model | all rays | near-axis | steep |
|---|---|---|---|---|
| 0.5 (0.229) | *exact excess size* | *0.042* | *0.367* | *0.099* |
| | CBS, 38.4 μm (referee vs itself) | 0.264 | **0.028** | 0.384 |
| | MSB, 38.4 μm | 0.694 | **0.110** | 0.575 |
| | WPM, 38.4 μm | 0.510 | 1.250 | 0.520 |
| | WPM, 19.2 μm crop = **D** | 1.141 | **1.272** | 1.217 |
| 1.0 (0.458) | *exact excess size* | *0.085* | *0.724* | *0.198* |
| | MSB, 38.4 μm | 0.733 | **0.152** | 0.581 |
| | WPM, 38.4 μm | 0.500 | 1.245 | 0.517 |
| | WPM, 19.2 μm crop = **D** | 1.117 | **1.243** | 1.283 |
| 2.0 (0.915) | *exact excess size* | *1.046* | *5.141* | *2.059* |
| | MSB, 38.4 μm | 0.502 | 1.436 | 0.462 |
| | WPM, 38.4 μm | 0.183 | 1.018 | 0.157 |
| | WPM, 19.2 μm crop = **D** | 0.258 | 1.046 | 0.312 |

**Read the near-axis column only.** There the referee agrees with itself across windows to
2.8 %, so the excess is a property of the cell. On steep rays it moves by 38 % between
windows, so nothing can be scored there. (MSB at 51.2 μm gives the same numbers as at
38.4 μm to two digits.)

### (c) At the level of the measurement, near-axis rays

| contrast | exact excess, share of signal | D's error | MSB's error |
|---|---|---|---|
| 0.5 | 0.357 | **0.455** | 0.039 |
| 1.0 | 0.683 | **0.849** | 0.104 |
| 2.0 | 1.051 | 1.099 | 1.509 |

D's error exceeds the excess it models at both contrasts where any model works. That is
the sense in which D gets the size right and the structure wrong: its excess has about the
right norm but points the wrong way on the rays where the referee can judge it. MSB, whose
only approximations are forward-only propagation, no evanescent coupling and an 84° cut on
radiated angles, is within 4–10 %. Neither survives $a\,\delta n=0.92\lambda$, where the
exact excess is larger than the linear signal itself; [[chen1998]] puts first Rytov's
limit at $0.24\lambda$.

## 4. What this does *not* establish

- **Scalar waves throughout.** The referee is exact for the scalar Helmholtz equation, which
  is the model GODT and this project use; it says nothing about polarization effects.
- **One cell, one illumination** for the excess comparison in §3; the design-level
  consequence is measured on all 180 illuminations in [[gridding-and-forward-model]], but
  with multi-slice Born, not with this referee.
- **Steep rays stay unadjudicated.** On rays with $\tan\theta_s>2$ even the exact solver's
  excess moves with the window (§3), so neither model can be scored there.
- **Nothing is validated at $a\,\delta n=0.92\lambda$.** exp11's top contrast therefore
  reports what multi-slice Born says there, not a checked answer.
- **The referee shares the voxelised object.** CBS solves the discretised problem exactly;
  the $0.1$ μm voxels and the phantom's sub-voxel texture are common to every model, so
  discretisation error relative to a continuous cell is not part of any comparison here.
- **The Mie check is in a framed window** with a finite measurement plane, which is why
  errors grow toward 55°. A near-to-far-field transformation over a closed surface would
  extend it; it was not needed for the referee's role here.

## Reproduction

```bash
cd lab
python exp09_full_wave_referee.py    # ~45 min on 4 cores, ~6 GB
```

Prints every number above and writes `lab/out/referee_mie.png` (CBS, MSB, WPM and first
Born against the Mie series for the $n=1.38$ sphere).
