---
type: note
projects: [godt]
created: 2026-10-05
updated: 2026-10-05
status: current
---

# Colour as a shear channel

Tests one idea: the polarization-grating shear scales with wavelength
($\Delta\propto\lambda$, confirmed on [[s42005-025-02416-3]]'s own grating numbers in
[[multishear-fusion-design-solution]] §6 E2), so resolving colours instead of averaging
them turns each colour into its own shear magnitude with no moving parts. Everything
numerical here is **computed** by `lab/exp08_colour_shear.py` on `lab/shearlib.py`.
The prior art is **web-sourced at abstract level**.

## Headline

1. **It works as a substitute for a mechanical shear change, and no better.** Three LED colours with $x$ and $y$ shear give the same error as the mechanical {0.8, 0.55} μm pair in $x$ and $y$ (0.00103 against 0.00104). The gain comes from the shear diversity, not the colour: three colours forced to the *same* shear gain nothing (§2).
2. **Blue + green is the right pair, and red is dead weight.** The optimal time share for red is zero. A longer wavelength has a smaller band, less information per photon, and reaches only 58.5 % of the 530 nm support (§4).
3. **It only matters for a large-shear instrument.** A single colour at the design-rule shear $\Delta\approx2\pi/K_{\max}$ already scores 0.00094. Colour multiplexing helps at a large shear such as 0.8 μm; at the design-rule shear it adds nothing (§7).
4. **It does not fill the missing cone.** The cone's angular opening is wavelength-free (§3).
5. **The hardware is already demonstrated in 2D.** The GODT authors' group posted a three-wavelength, single-colour-frame polarization-grating shearing preprint on 2026-09-07, for 2D replica removal. This note studies the tomographic side (§8).
6. **The separation theorem breaks across colours.** $J=\rho\,S_{\mathbf t}M$ no longer factorises, because each colour has its own $M_c(\mathbf K)$ (§1).

## 1. Model

Each colour $c$ is the same instrument at wavelength $\lambda_c$. Angles and apertures are
geometric and stay fixed (60° cone, NA 1.33, $n_0=1.33$, GODT simulation optics), so:

- every wavenumber scales as $1/\lambda_c$: Ewald radius, band edge, cap density $M_c(\mathbf K)$, Rytov weight;
- the shear scales the other way, $\Delta_c=0.8\,\mu\text{m}\times\lambda_c/530$ nm;
- the estimand is the dispersion-free contrast $f=n_m^2-n^2$, with $O=k_0^2f$. A photon at a shorter wavelength therefore carries more information about $f$: $\rho_c=\rho\,(\lambda_\text{ref}/\lambda_c)^2$ *(inference; follows from the $k_0^4/k_m^2$ in [[multishear-fusion-design-solution]] §1)*.

$$J(\mathbf K)=\sum_c\rho_c\,M_c(\mathbf K)\sum_j t_{cj}\,\big|e^{i\boldsymbol\Delta_{cj}\cdot\mathbf K_\perp}-1\big|^2$$

This is a sum of separable terms, not a separable product. The problem stays convex,
but the 2D shortcut of [[optimal-experimental-design]] is gone, and designs must be
scored on the 3D volume.

| colour | shear | $K_{\max}$ (rad/μm) | $\rho_c$ | support vs 530's | in-band nulls (rad/μm) |
|---|---|---|---|---|---|
| 450 nm | 0.679 μm | 34.65 | 123.5 | 162.2 % | 9.25, 18.50, 27.75 |
| 530 nm | 0.800 μm | 29.42 | 89.0 | 100 % | 7.85, 15.71, 23.56 |
| 635 nm | 0.958 μm | 24.56 | 62.0 | 58.5 % | 6.56, 13.11, 19.67 |

All scores below are the mean relative Wiener MSE **on the 530 nm support**, flat prior,
equal total photon count, $\rho=100$ at 0.5 μm as in exp03 and exp04. "Lost" is the share
of that support with relative MSE above 0.25. The grid is $320^2\times96$ at 0.08 μm,
which has the same $K$ lattice as exp03's and room for the blue band.

**Consistency check.** With every colour forced onto the 0.5 μm geometry, the script
reproduces [[multishear-fusion-design-solution]]'s rows exactly (0.03420 and 3.02 % for
the single 0.8 μm shear, 0.00700, 0.00275, 0.00097). In that fixed-geometry model two
colours in $x$ score 0.00654, three score 0.00628, and three in $x$ and $y$ score
0.00087. Those are an upper bound on what per-wavelength geometry gives.

## 2. What colour buys, with per-wavelength geometry

| design (equal shares) | mean rel. MSE | lost | unsampled |
|---|---|---|---|
| 530 nm single, $x$ | 0.03534 | 3.10 % | 0 |
| B+G+R, $x$, all at the **same** 0.8 μm shear (control) | 0.03567 | 3.10 % | 0 |
| mechanical 0.8 + 0.55 μm at 530 nm, $x$ | 0.00741 | 0.52 % | 0 |
| blue + red, $x$ | 0.01798 | 1.82 % | 0.25 % |
| blue + green, $x$ | 0.00691 | 0.53 % | 0 |
| blue + green + red, $x$ | 0.00702 | 0.53 % | 0 |
| 530 nm, $x$ + $y$ | 0.00291 | 0.10 % | 0 |
| mechanical {0.8, 0.55} μm at 530 nm, $x$ + $y$ | 0.00104 | 0.00 % | 0 |
| blue + red, $x$ + $y$ | 0.00467 | 0.30 % | 0.25 % |
| **blue + green, $x$ + $y$** | **0.00090** | 0.00 % | 0 |
| blue + green + red, $x$ + $y$ | 0.00103 | 0.00 % | 0 |

- **The control row is the mechanism.** An achromatic shear with three colours scores the same as one colour. Colour helps only because the grating makes the shear chromatic.
- **Colour replaces a second magnitude one for one.** In one direction it is a 5× gain over the single shear, the same as the mechanical pair.
- **It does not replace the second direction.** The residual 0.53 % is the blind plane $K_x=0$, which no magnitude change touches ([[multishear-fusion-design-solution]] §2).
- **Blue + red fails, although it has the largest shear ratio (1.41).** Red covers 58.5 % of the 530 nm support, so over the rest blue is alone and keeps its own nulls. A quarter of a percent of the support is sampled by neither.
- **Blue + green beats the mechanical pair by 14 %.** That margin is blue's higher $\rho_c$, not a better null fill.

Blue also adds voxels outside the 530 nm support, worth 62 % of it. That is the
resolution gain of a shorter wavelength and any blue instrument gets it. Out there blue
is alone and keeps its null comb.

Sample of where the error sits, three colours against the mechanical pair, both in $x$ and $y$:

| $K_\perp$ band (rad/μm) | share of support | three colours | mechanical pair |
|---|---|---|---|
| 0 – 10 | 11.9 % | 0.00062 | 0.00055 |
| 10 – 20 | 41.2 % | 0.00107 | 0.00107 |
| 20 – 24.6 | 22.6 % | 0.00101 | 0.00136 |
| 24.6 – 29.4 (beyond red's band) | 24.3 % | 0.00118 | 0.00094 |

## 3. The missing cone is untouched

Every Ewald cap passes through the origin at the same angle whatever its radius, so the
cone boundary near the origin is $K_\perp=K_z\cot\theta$ for every colour
([[missing-cone-problem]]). The inner boundary at fixed $K_z$:

| $K_z$ (rad/μm) | 450 nm | 530 nm | 635 nm |
|---|---|---|---|
| 0.98 | 0.49 | 0.49 | 0.49 |
| 1.96 | 0.98 | 0.98 | 0.88 |
| 3.93 | 1.72 | 1.65 | 1.47 |
| 5.89 | 2.20 | 1.98 | out of range |

At second order the boundary bends inward as $K_z^2/(2k_m\sin^3\theta)$, more for a
longer wavelength *(inference from the cap geometry; the table agrees)*. So supports of
different colours are **not nested**: 0.26 % of the 530 nm support is unsampled at
450 nm. The three colours together sample 52.7 % of the nominal cone
($K_\perp<|K_z|\cot\theta$, $|K_z|\le7.88$) against 30.7 % for 530 nm alone. That
difference is 0.2 % of the support and comes from blue's longer axial reach
($K_z$ up to 9.29 against 7.88), not from any angular fill. Only 0.4 % of the voxels
that blue and red add lie in the cone.

## 4. Which colours, and how much of each

Simplex scan of the colour shares for the $x$ + $y$ design, step 0.05:

| shares (B, G, R) | mean rel. MSE |
|---|---|
| **0.55, 0.45, 0 (optimum)** | **0.00089** |
| ½, ½, 0 | 0.00090 |
| ¼, ½, ¼ (Bayer pixel areas) | 0.00101 |
| ⅓, ⅓, ⅓ | 0.00103 |
| 0, ½, ½ | 0.00265 |
| 0, 1, 0 | 0.00291 |
| ½, 0, ½ | 0.00467 |
| 1, 0, 0 | 0.00470 |

The rule is to use the **two shortest wavelengths available**. Blue alone scores worse
than green alone on this support mainly because of the 0.26 % sliver it does not sample,
which alone adds 0.0026.
The comparison is at equal photon *counts*. LED power, camera quantum efficiency and
phototoxicity in the blue are not modelled, and all three push the other way.

## 5. Sample dispersion

The model assumes the object is the same at every colour. If the contrast $f$ differs by
an rms fraction $\varepsilon$ between a colour and 530 nm, and that is treated as
coherent model error (the construction of E1 in [[multishear-fusion-design-solution]]),
each off-reference channel is capped at an SNR of $1/\varepsilon^2$ per voxel. The
photon-limited level here is 33.9 dB.

| $\varepsilon$ | ceiling | B+G, $x$+$y$ | B+G+R, $x$+$y$ |
|---|---|---|---|
| 0 | ∞ | 0.00090 | 0.00103 |
| 1 % | 40 dB | 0.00092 | 0.00105 |
| 2 % | 34 dB | 0.00097 | 0.00110 |
| 5 % | 26 dB | 0.00126 | 0.00140 |
| 10 % | 20 dB | 0.00175 | 0.00201 |
| 20 % | 14 dB | 0.00241 | 0.00295 |

Below 2 % dispersion is invisible. At 20 % the colour design is still ahead of a single
colour in $x$ and $y$ (0.00291). The real $\varepsilon$ of cells between 450 and 530 nm
is **not in this wiki**, and neither is a source for it. This is the conservative
reading: a one-parameter dispersion fit per voxel would recover most of the loss
*(inference)*. It also means geometry and composition compete for the same data, as the
hyperspectral ODT literature in §8 treats them.

## 6. Single-shot with a colour camera

A mosaic sensor costs nothing in total photons by itself: each colour gets its pixel-area
share of one exposure. Crosstalk does cost, because unmixing amplifies shot noise
(`shearlib.unmix_shares`, a **model, not sourced**; adjacent-channel leakage $x$,
next-nearest $x^2$):

| leakage $x$ | effective shares B, G, R | sum | three colours, $x$+$y$ |
|---|---|---|---|
| 0 | 0.250, 0.500, 0.250 | 1.000 | 0.00101 |
| 0.1 | 0.220, 0.386, 0.220 | 0.825 | 0.00122 |
| 0.2 | 0.182, 0.269, 0.182 | 0.632 | 0.00160 |
| 0.3 | 0.142, 0.172, 0.142 | 0.456 | 0.00223 |

So one colour exposure per shear direction replaces two or three monochrome ones, at a
20 to 60 % error penalty for a leakage of 0.1 to 0.2. Two physical configurations
($x$ and $y$) are still needed. **Not modelled:** the mosaic samples red and blue at
twice the pixel pitch, so the image needs 2× more magnification or the colour planes
alias.

Each LED's own linewidth does not matter once colours are fused: with 20, 35 and 17 nm
FWHM (the 35 nm is GODT's; the others are typical values, general knowledge) the
three-colour scores move by under 4 %.

## 7. At the design-rule shear the idea is not needed

Grating separation set for 0.2 μm at 530 nm, so no colour has an in-band null:

| design, $x$ + $y$ | mean rel. MSE |
|---|---|
| 530 nm only | 0.00094 |
| three colours | 0.00097 |
| three colours, fixed geometry | 0.00092 |
| blue + green | 0.00079 |

With no comb to fill, shear diversity is worth nothing. The blue + green gain is again
blue's information per photon. And 530 nm alone at 0.2 μm (0.00094) already beats three
colours at 0.8 μm (0.00103). Colour multiplexing is therefore a **repair for an
instrument that must run at large shear**, for instance because fringe visibility or
phase sensitivity demands it. Whether GODT is such an instrument depends on
$|\mu(\Delta)|$ and on assumption A4, both still open.

## 8. Prior art (web search, abstracts only)

- **Closest: from the GODT authors' group.** Zdańkowski, Piekarska, Rogalski, Stefaniuk & Trusiak, "Single-Shot Replica-Free Common-Path Quantitative Phase Imaging with Doubled Field of View," Research Square preprint, doi:10.21203/rs.3.rs-10652545/v1, posted 2026-09-07. The abstract describes wavelength-dependent diffraction in polarization gratings moving "shear control entirely into the spectral domain", three discrete wavelengths at once, and one multichannel colour hologram, for replica separation in 2D total-shear imaging. It is the single-shot successor of R2D-QPI ([[forward-citation-sweep-2026-07]]).
- Toto-Arellano, Reyes & Escamilla Angeles, *Appl. Opt.* **65**, 1437 (2026), doi:10.1364/AO.586119. White light, RGB channel separation, phase derivative in one exposure with a cyclic-path interferometer. Its blue channel was unusable for lack of fringe visibility and camera response. That is a practical warning, since §4 makes blue the most valuable colour.
- Wavelength-resolved ODT for composition: Jung, Kim, Yoon & Park, *Opt. Express* **24**, 2006 (2016); Sung, *Phys. Rev. Applied* **19**, 014064 (2023); Hörmann et al., arXiv:2401.07391.
- Wavelength in place of angle: Wu et al., *Light Sci. Appl.* (2024), doi:10.1038/s41377-024-01568-1, lens-free and on-chip.
- Reported by the earlier search in this thread and **not re-checked**: *Opt. Lett.* **49**, 3066 (2024) and Matlis et al., *Nat. Commun.* (2012), on wavelength-encoded angles in one shot.

**What this note adds.** The hardware idea, a polarization-grating module giving a
different shear per colour read by a colour camera, is published. This note contributes
the quantitative side for gradient-regime tomography: which colours, how many, what
dispersion and crosstalk cost, and that the correct shear makes the trick unnecessary.

## 9. What this does *not* establish

- **Fringe visibility per colour.** Unit visibility is assumed. GODT's $|\mu(\Delta)|$ is unmeasured, and here it would differ per colour.
- **Polarization-grating efficiency off its design wavelength.** A grating built for 530 nm leaks zero order at 450 and 635 nm (general knowledge, not checked). That lowers visibility and adds an unsheared background.
- **Chromatic aberration and immersion dispersion.** $n_0$ is held at 1.33 for every colour and the objective is taken as perfectly achromatic.
- **The real dispersion of cells**, as §5 says.
- **Multiple scattering.** [[multiple-scattering-check]] showed the null positions survive it exactly and the margins shrink with contrast. Per colour the contrast variable $a\,\delta n/\lambda$ differs, so blue reaches the Rytov bound first *(inference)*.
- **Nearest-neighbour gridding (A3)** and **sequential-channel noise (A6)** are inherited unchanged.
- **Spatial coherence.** Like E2, this treats each colour as coherent and monochromatic.

## Reproduction

```bash
cd lab
python exp08_colour_shear.py     # ~1 min, numpy + matplotlib only
```

Writes `lab/out/colour_shear.png` and prints every number above. The two helpers it
adds to `shearlib` (`at_wavelength`, `unmix_shares`) are covered by the selftest.
