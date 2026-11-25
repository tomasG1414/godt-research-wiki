---
type: concept
projects: [godt]
created: 2026-11-25
updated: 2026-11-25
status: current
---

# Multi-slice forward models

Forward models that march a field through a 3D object slice by slice, propagating
through homogeneous medium between slices by angular spectrum. They are cheap (a few
2D FFTs per slice) and include forward multiple scattering, which first Born does not.
They drop backward scattering and, usually, evanescent coupling. Three variants appear
in this wiki:

| | each slice does | linear limit | used in |
|---|---|---|---|
| **BPM** | phase screen of the *averaged* index, applied to all angles | paraxial | not used: fails beyond ~10° ([[brenner1993]]) |
| **WPM** | phase screen $e^{ik_z(n(\mathbf r),\mathbf k_\perp)\delta z}$ per plane-wave component ([[brenner1993]]) | first Born × $k_z/k_{iz}$ | GODT's simulations; exp06 |
| **MSB** (multi-slice first Born) | radiates through the exact Green's function $-i\,\delta z/(2k_z)$ of the *scattered* wave | first Born, exactly | exp09, exp11 |

MSB is related to the multi-layer Born model of Chen, Ren, Liu, Chowdhury & Waller,
*Optica* **7**, 394 (2020) *(not ingested; cited from its abstract page)*. The version
here is forward-only and lives in `lab/shearlib.py` (`MSB`).

## The obliquity law (WPM)

For one slice WPM reduces analytically to the transmission
$t(\mathbf r)=\exp[i\,k_z(n(\mathbf r),\mathbf k_{i\perp})\,\delta z]$: the axial wavenumber
is taken at the **incident** transverse frequency, so the scattered amplitude carries
$1/k_{iz}$ where first Born carries $1/k_z$ of the scattered wave. On a thick object, with
WPM and Born sharing one simulation window, the law holds exactly ([[full-wave-referee]] §2):

$$\text{WPM}=\text{Born}\times\frac{k_z}{k_{iz}}\times\frac1s\sum_{j=1}^{s}e^{-iK_z\,j\,\delta z/s},$$

the sum being where WPM places each of the $s$ thin screens in a voxel. Measured residual:
0.0034 at every step size from 0.1 to 0.0125 μm, which is second order. With GODT's 60°
illumination $k_z/k_{iz}$ runs from about 0.2 to 2, so WPM data fed to a Rytov/Ewald
inverse carry an angle-dependent error of up to a factor 5 before any multiple scattering.
An earlier, larger steep-angle gap in exp06 was the simulation window, not this law
([[full-wave-referee]] §2, [[multiple-scattering-check]] §3).

## How accurate they are beyond Born

Against the exact scalar Mie series for a sphere ([[full-wave-referee]] §1; $a=1.5$ μm,
normal incidence, angles up to 55°):

| sphere index (medium 1.33) | phase delay | first Born | WPM | MSB |
|---|---|---|---|---|
| 1.34 | 0.38 rad | 15.5 % | 3.2 % | 0.8 % |
| 1.38 | 1.88 rad | 79.9 % | 3.5 % | 4.2 % |

WPM's error concentrates at 20–55° (20–55 %), where the obliquity law bites; MSB's is
spread evenly. At 60° illumination MSB's error rises to 13 % on the strong sphere and WPM's
to 15 %.

On exp06's cell phantom, what matters for a multi-shear design is the **nonlinear
excess** $e(c)=d(c)-(c/\varepsilon)\,d(\varepsilon)$, the part of the data first Born
misses. On the rays where an exact referee can pin it down, MSB reproduces it to 11–15 %
of its size (4–10 % of the measured signal) up to $a\,\delta n=0.46\lambda$, while WPM's
excess is off by more than its own size. Neither holds at $0.92\lambda$
([[full-wave-referee]] §3). At the design level the two give similar margins anyway,
because what decides them is how large the bias is ([[gridding-and-forward-model]] §4).

## Two practical traps

- **The window.** A multi-slice model and the Born formula it is compared with must see
  the same transverse window. Steep scattered light walks $\delta z\tan\theta$ sideways,
  10–25 μm through a 6 μm cell at 60°, so an absorbing edge or a crop removes different
  light from each. Measure on the whole window at the reconstruction lattice
  (`shearlib.lattice_samples`), which is what periodic first Born does.
- **Step size at finite contrast.** WPM's linear limit converges in $\delta z$, but its
  nonlinear response at steep angles does not: at $\delta n$ amplitude $2\times10^{-3}$ the
  deviation from the law grows from 4 % to 12 % as $\delta z$ goes from 0.1 to
  0.0125 μm ([[full-wave-referee]] §2).
