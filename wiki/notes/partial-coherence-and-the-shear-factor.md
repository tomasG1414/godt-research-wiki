---
type: note
projects: [godt]
created: 2026-11-25
updated: 2026-11-25
status: current
---

# Partial coherence and the shear factor

Every result in this project treats shear data as the coherent Rytov spectrum multiplied by
$G(\mathbf K_\perp)=e^{i\mathbf K_\perp\cdot\boldsymbol\Delta}-1$, with spatial coherence
entering only as a scalar fringe visibility $|\mu(\Delta)|^2$ per channel
([[multishear-fusion-design-solution]] §4, where the visibility model is labelled
unsourced). [[multiple-scattering-check]] §2 proves the multiplicative identity for one
coherent realisation only. GODT's experiment used an LED (530 nm, 35 nm bandwidth;
[[s42005-025-02416-3]]), and partially coherent imaging theory says the transfer
function then becomes bilinear ([[partially-coherent-transfer-function]],
[[mehta-sheppard2008]]). This note asks whether the separable model survives, and where it
fails.

Everything is **computed** by `lab/exp13_partial_coherence.py`: a 2D (x–z) simulation at
GODT's simulation optics (0.5 μm, $n_0=1.33$, 60° illumination), an extended source of
41 mutually incoherent tilted plane waves (uniform in transverse wavevector over a cone of
half-angle $\alpha$), multi-slice Born per source point, five-step phase shifting of the
sheared interferogram with intensities summed over the source, and the empty-field
fringe divided out. The script was drafted by a subagent; the review here fixed a
units bug in its forward-model check and added the pupil-interior analysis that the
conclusions rest on.

## Headline

1. **The derived model is right.** The measured phase difference is a source average of
   the per-source coherent shear data, each on its own Ewald cap, with complex weights
   $e^{i d_s\Delta}/\mu$ (§1). It matches the simulation to 0.1–1.4 % for every cone,
   shear and object tested.
2. **The separable model holds for small sources, and fails first for out-of-focus
   structure.** Away from the pupil edge, $G\times$ (coherent phase) is within 1 % of the
   simulation up to a cone of about 0.5–1° for a weakly scattering cell-like object, and up
   to about 4° for an in-focus thin slice. A slice 2 μm from focus is off by 3.8 % at 1°,
   13 % at 2° and 36 % at 4°. The scale is $|z|\,|v|\,h\ll1$, with $h$ the source's spread in
   transverse wavevector and $v$ the ray-slope difference (§2).
3. **Partial coherence does not fill the nulls usefully.** Data appear on $G$'s null
   frequencies only through out-of-focus structure: 21 % of the band's rms for a slice 2 μm
   from focus at a 4° cone, 0.2 % in focus. That is depth-mixed signal, not information about
   the null frequencies of one slice (§3).
4. **The visibility model needs the plane of incidence.** For a cone about the 60°
   direction, $|\mu|$ for a shear perpendicular to the plane of incidence is the jinc
   `shearlib.visibility` uses. For a shear in that plane the argument carries an extra
   $\cos60°$: at an 8° cone and 0.8 μm, $|\mu|$ is 0.895, not 0.625 (§4).
5. **GODT's bandwidth is harmless to the separable model; a hard pupil edge is not.** With
   a polarization-grating shear $\propto\lambda$, a point source and 35/530 fractional
   bandwidth, the empty-field visibility stays 1.000 (0.951 for a wavelength-independent
   0.8 μm shear), and away from the pupil edge the monochromatic model holds to 0.4–0.8 %
   in focus. Over the full band, the blurring of the hard pupil edge by the source spread or
   by the spectrum costs 6–12 %, even at a 0.25° cone, for objects with content at the band
   edge (§5).

## 1. The derived model

A source point $s$ is a plane wave of transverse wavevector $q_s=q_0+d_s$, weight $w_s$.
Five-step phase shifting returns the fringe phasor, which is the mutual intensity of the
image field across the shear:

$$\Gamma(x)=\sum_s w_s\,u_s(x+\Delta)\,u_s^*(x),\qquad
\Gamma_0=e^{iq_0\Delta}\mu(\Delta),\qquad \mu(\Delta)=\sum_s w_s e^{i d_s\Delta}$$

(the second is van Cittert–Zernike for the empty field). For a weak object
$u_s=e^{iq_sx}(1+\psi_s)$, and to first order

$$\delta\varphi(x)=\arg\frac{\Gamma}{\Gamma_0}
=\operatorname{Im}\Big\{\frac1\mu\sum_s w_s\,e^{id_s\Delta}\big[\psi_s(x+\Delta)+\psi_s^*(x)\big]\Big\},$$

where each $\psi_s$ lives on its own Ewald cap, $K_z=k_z(\kappa+q_s)-k_z(q_s)$ *(derivation,
checked numerically below)*. For a point source this is $G$ times the coherent phase. Two
things change for an extended one. Each source point samples a slightly different cap, so
the coherent phases decorrelate with depth. And the weights $e^{id_s\Delta}/\mu$ are
complex, which leaks log-amplitude into the phase (§2). Linearising $K_z$ in $d_s$ gives
a closed form per slice at depth $z$,

$$G_\text{pc}(\kappa;z)=\frac{e^{i\kappa\Delta}\mu(\Delta-zv)-\mu(\Delta+zv)}{\mu(\Delta)},
\qquad v=\tan\theta_i-\tan\theta_s(\kappa),$$

which is $G$ at $z=0$ or for a point source, and is not zero on $G$'s nulls. It is used
here for orientation (§3); its own accuracy was not established separately.

**Checks.** Multi-slice Born against first Born at contrast $10^{-4}$: $5\times10^{-5}$.
A point source against $G\times$ the coherent phase: 0.29–0.31 % at the weak cell's contrast,
ten times smaller at ten times lower contrast, so it is the weak-object linearisation.

## 2. When the separable model holds

Relative error of $G\times$ coherent phase (M0), of an unweighted source average (M1) and of
the derived model (M2) against the simulation, shear 0.8 μm, with a band of width
$h+2$ bins at each pupil edge removed from both images:

| object | cone 0.25° | 0.5° | 1° | 2° | 4° | 8° | M2, worst |
|---|---|---|---|---|---|---|---|
| thin slice in focus | 0.67 % | 0.68 % | 0.83 % | 1.87 % | 4.09 % | 11.8 % | 0.87 % |
| thin slice 2 μm from focus | 0.76 % | 1.18 % | 3.83 % | 13.2 % | 35.8 % | 114 % | 1.40 % |
| weak cell, 10 × 4.8 μm | 0.36 % | 0.71 % | 1.63 % | 4.89 % | 9.20 % | 26.3 % | 0.21 % |

The cone half-angle is in the medium; 1° corresponds to a condenser NA of 0.023. The
validity length $1/(h\,|v|)$ for the median ray slope is 15 μm at 0.25°, 3.8 μm at 1° and
0.96 μm at 4°, which tracks where M0 fails for structure 2 μm from focus. Shear magnitude
matters less than the cone: at a fixed cone the three shears (0.2, 0.55, 0.8 μm) differ by
less than a factor 2 in the full-band table the script prints.

**The amplitude leak.** The complex weights also mix log-amplitude into the measured
phase: 0.5–2 % of the signal at a 1° cone, 4–10 % at 4° (15–22 % for the slice 2 μm from
focus), 12–28 % at 8° (55 % out of focus).
GODT discards amplitude, so this is an error source for a phase-only pipeline
*(inference)*.

## 3. The nulls

$|\mathcal F[\delta\varphi]|$ on the three null bins of the 0.8 μm shear, relative to the rms
over the band (zero by construction for the coherent model):

| object | cone 1° (m = 1, 2, 3) | 4° | 8° |
|---|---|---|---|
| thin slice in focus | 0.0001, 0.0000, 0.0000 | 0.0019, 0.0001, 0.0002 | 0.0088, 0.0004, 0.0008 |
| thin slice 2 μm from focus | 0.0128, 0.0059, 0.0014 | 0.2134, 0.0774, 0.0121 | 0.5267, 0.0445, 0.0147 |
| weak cell | 0.0013, 0.0009, 0.0003 | 0.0243, 0.0126, 0.0030 | 0.1633, 0.0381, 0.0050 |

The derived model reproduces every entry to within a few percent of its value. What
appears on the nulls is out-of-focus structure seen through a depth-dependent factor
$G_\text{pc}$, so it does not supply the missing in-focus information that a second shear
does ([[finite-shear-null-experiment]]).

## 4. The visibility model

Sum over a circular cone of half-angle $\alpha$ about the 60° direction, against the
small-cone formula $\mu=2J_1(r)/r$ with
$r=k_m\sin\alpha\,\sqrt{(\Delta_\parallel\cos60°)^2+\Delta_\perp^2}$ and against
`shearlib.visibility` (the same jinc with $r=\pi\Delta/l_c$, $l_c=\lambda/(2n_0\sin\alpha)$,
no $\cos60°$):

| cone | shear | $|\mu|$, shear in the plane of incidence | perpendicular | shearlib |
|---|---|---|---|---|
| 2° | 0.8 μm | 0.9932 | 0.9730 | 0.9730 |
| 4° | 0.8 μm | 0.9730 | 0.8951 | 0.8951 |
| 8° | 0.55 μm | 0.9494 | 0.8086 | 0.8089 |
| 8° | 0.8 μm | 0.8950 | 0.6247 | 0.6253 |

The formula matches the direct sum to $10^{-3}$. Under GODT's conical scan the plane of
incidence rotates with the illumination azimuth, so a fixed shear direction sees both
values over one scan *(inference)*. The coherence argument of
[[multishear-fusion-design-solution]] §4 (falling coherence favours smaller shears) keeps
its direction but uses the perpendicular, worst-case value throughout.

## 5. Bandwidth, and the pupil edge

Gaussian spectrum with GODT's 35/530 fractional bandwidth, 15 wavelengths, geometry scaled
per wavelength (`shearlib.at_wavelength`), grating shear $\Delta\propto\lambda$, achromatic
phase steps. Point source, relative error of the monochromatic model:

| object | full band, 0.2 / 0.55 / 0.8 μm | pupil interior |
|---|---|---|
| thin slice in focus | 7.3 / 10.7 / 12.4 % | 0.80 / 0.71 / 0.83 % |
| thin slice 2 μm from focus | 9.7 / 13.9 / 16.7 % | 3.81 / 1.86 / 1.67 % |
| weak cell | 9.0 / 9.1 / 8.8 % | 0.73 / 0.40 / 0.41 % |

The spectral average of $G$ (`shearlib.chromatic_G`) changes these by at most 0.25
percentage points: the per-wavelength Ewald caps and pupil edges matter, not the spread of
the shear. Adding a 2° cone takes the weak cell to 3.9–4.4 % in the interior.

**The pupil edge.** The large full-band numbers come from a few frequencies at the
collection pupil's edge, which the source spread $h$ (or the spectrum) moves between
source points: compare the full-band and interior columns, or §2's interior table with the
6–10 % full-band errors the script prints at a 0.25° cone. This is the classic apodisation of a
hard pupil by partial coherence ([[partially-coherent-transfer-function]]). It is not a
property of the shear factor, but a phase-only pipeline that assumes a hard pupil will see
it.

At full contrast (projected phase 0.87 rad) the weak-object model itself fails by 11–27 %,
the coherent measurement already differs by 6 % at a 0.25° cone over the full band, and the
question becomes the multiple-scattering one of [[full-wave-referee]].

## 6. What this does *not* establish

- **2D, scalar, one illumination direction.** A 3D cone was used only for the visibility
  table. The conical scan, in which the plane of incidence rotates, was not simulated.
- **Uniform slit or cone sources.** A real LED behind a condenser has a non-uniform angular
  spectrum; only the shape of $\mu$ changes in the derivation, but the numbers do.
- **The closed form $G_\text{pc}$ is orientation.** The derived source-average model (M2) is
  what was checked against the simulation.
- **No tomographic consequence was computed.** Whether depth-dependent decorrelation at
  GODT's actual condenser aperture costs reconstruction accuracy needs this model inside
  the 3D pipeline (`lab/exp11_exact_design.py`).
- **GODT's condenser aperture was not found** in the sources ingested here, so where its
  illumination sits on the cone axis is open.

## Reproduction

```bash
cd lab
python exp13_partial_coherence.py    # ~35 s
```

Prints every number above and writes `lab/out/partial_coherence.png`.
