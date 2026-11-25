---
type: note
projects: [godt]
created: 2026-11-25
updated: 2026-11-25
status: current
---

# Spatially varying shear

Relaxes assumption A5 of [[multishear-fusion-design-problem]]: one uniform shear per
channel. In real two-direction instruments the shear varies across the field, from
tube-lens distortion in [[2605-04297-double-sagnac-ls-qpm]]. A smooth variation is
deterministic and can be calibrated, unlike the random shear error priced in
[[multishear-fusion-design-solution]] §6 (E1). Two questions: does the variation smear
the null comb of [[finite-shear-null-experiment]] so that the data carry information
there, a "free" multi-shear effect? And what does it cost an estimator that assumes a
uniform shear?

Everything is **computed** by `lab/exp07_spatially_varying_shear.py` at projection level
(1D and 2D, as exp01–02; no Ewald caps, no Rytov factor). The data model is
$\psi(x)=\varphi(x+\Delta(x))-\varphi(x)$ + noise with $\Delta(x)=0.8\,\mu\text{m}\,
(1+\varepsilon f(x))$, $f$ a linear ramp or a quadratic bowl of unit range across the
field, $\varepsilon$ = 1 %, 3 %, 10 %; band 2–31.2 rad/μm (GODT simulation optics),
flat prior, 20 dB per frequency unless stated, dimensionless shear power as in exp03+.
Two estimators: **calibrated** (knows $\Delta(x)$, exact MAP for the space-variant
operator) and **naive** (exp02's per-frequency Wiener filter for a uniform shear).
The script was drafted by a subagent and reviewed here. The review found and fixed a
sign-convention bug in its 2D operator; the fix and its check are recorded in the script
and in `lab/README.md`.

## Headline

1. **A calibrated estimator gets neither a gain nor a loss at 20 dB.** The null is
   smeared over a band of width $\approx K_m\varepsilon$ (measured 2.24 rad/μm against a
   predicted 2.36 for the third null and a 10 % ramp), so the error per Fourier mode at the
   null drops from 0.70 to 0.40. But the area is conserved: each null still costs 2.5–3.1
   modes, against 3.06 for a uniform shear, and the overall RMSE stays at 0.21–0.22
   (uniform: 0.208). A smooth shear field is not a free multi-shear design. The
   {0.8, 0.55} μm pair at the same exposure reaches 0.099.
2. **Ignoring the variation is expensive.** With the naive filter a 3 % ramp (6.9 nm rms
   shear deviation) raises the RMSE from 0.21 to 0.28 at 20 dB, and 10 % to 0.66; the
   model error passes the noise above $K\approx14$ rad/μm at 3 %. At 40 dB even a 1 %
   ramp (2.3 nm rms) costs 0.11 → 0.14. The shear field has to be calibrated to a few
   nanometres, the same scale as the 10 nm random error of E1.
3. **Real diversity comes from moving the sample, not from the variation itself.** Two
   exposures at half the photons each, with the sample translated across the shear
   gradient between them, see every point at two local shears. With a 10 % ramp and the
   two looks 40 nm apart, the RMSE falls from 0.208 to 0.124 and the third null's cost from
   3.06 to 0.53 modes: about three quarters of the way to the deliberate {0.8, 0.55} μm
   pair. With a
   3 % ramp (12 nm apart) it is 0.179.
4. **The local picture is exact.** Windowed over 6 μm, the worst-recovered component
   sits at $2\pi m/\Delta(x_0)$ to within 0.04 rad/μm, and a phase-space (Husimi)
   smoothing of the local Wiener error predicts the windowed error to within 6 %.
5. **The 2D images agree.** On a band-passed Shepp–Logan phantom, the calibrated
   estimator scores 0.223–0.226 against the uniform shear's 0.222, and the naive one
   0.227–0.229 at 3 % and 0.27–0.29 at 10 %, the damage concentrated in the third null band.

## 1. The null moves, it does not vanish

For a uniform shear the lost functions are the plane waves at $K_m=2\pi m/\Delta_0$. For a
smooth $\Delta(x)$ the natural candidates are chirps $e^{2\pi i m\alpha(x)}$ with
$\alpha'=(1+\Delta'/2)/\Delta$ (second-order WKB, exact for the ramp), whose local
frequency is $2\pi m/\Delta(x)$. Rayleigh quotients $\|B c\|^2/\|c\|^2$ of Hann-windowed
test functions:

| profile | $\varepsilon$ | m | chirp | plane wave at $K_m$ |
|---|---|---|---|---|
| uniform | 0 | 1, 3 | $3.21\times10^{-3}$ | $3.21\times10^{-3}$ |
| ramp | 3 % | 3 | $3.21\times10^{-3}$ | $9.60\times10^{-3}$ |
| ramp | 10 % | 3 | $3.21\times10^{-3}$ | $7.31\times10^{-2}$ |
| bowl | 10 % | 3 | $3.40\times10^{-3}$ | $6.11\times10^{-2}$ |

The chirp is as nearly null as the uniform shear's plane wave; the floor of
$3.2\times10^{-3}$ is the window. The plane wave at the old null frequency is now measured.
The lost information has moved to a different function.

## 2. Calibrated estimator

51.2 μm field, exact posterior. "Lost" is the sum of relative MSE over $\pm1.5$ rad/μm of
the null, in modes, per sign of $K$:

| case | m = 1 lost | m = 2 lost | m = 3 lost | overall RMSE |
|---|---|---|---|---|
| uniform 0.8 μm | 3.06 | 3.06 | 3.06 | 0.2077 |
| ramp 3 % | 3.06 | 3.06 | 3.05 | 0.2103 |
| ramp 10 % | 3.05 | 3.01 | 2.94 | 0.2116 |
| bowl 10 % | 2.97 | 2.91 | 2.46 | 0.2182 |
| {0.8, 0.55} pair, half exposure each | 0.17 | 0.14 | 1.30 | 0.0987 |

The equivalent width (area / peak) of the error bump over the whole field approaches the
prediction $K_m\varepsilon$ for a ramp once that exceeds the intrinsic 0.394 rad/μm: m = 3 at
10 %, 2.24 against 2.36; m = 2, 1.54 against 1.57; at 1 % it stays at the intrinsic width.
The per-mode error is the field average of the local Wiener error, to a maximum deviation
of 0.05–0.13 (ramp) and 0.05–0.22 (bowl) in rms relative error.

**At higher SNR the picture changes.** At 40 dB the uniform shear loses 3.10 modes over
the three nulls per 51.2 μm, a 10 % ramp 1.34 and a 10 % bowl 2.67 (RMSE 0.114, 0.076,
0.110). So the variation does recover information once noise is low. The ramp's larger
gain may partly come from its jump at the periodic boundary of the simulated field, which
the bowl does not have *(inference; not isolated)*.

## 3. The naive estimator

Relative RMSE at 20 dB (40 dB in brackets), calibrated vs naive, and the rms shear
deviation the naive filter ignores:

| case | rms $\delta\Delta$ | calibrated | naive | bias alone |
|---|---|---|---|---|
| ramp 1 % | 2.3 nm | 0.208 (0.113) | 0.217 (0.144) | 0.158 (0.142) |
| ramp 3 % | 6.9 nm | 0.210 (0.109) | 0.281 (0.286) | 0.239 (0.286) |
| ramp 10 % | 23.1 nm | 0.212 (0.076) | 0.658 (0.871) | 0.641 (0.871) |
| bowl 3 %, nominal shear | 10.7 nm | 0.217 (0.058) | 0.365 (0.424) | 0.333 (0.423) |
| bowl 10 %, nominal shear | 35.8 nm | 0.218 (0.110) | 0.939 (1.274) | 0.928 (1.274) |

The bias grows with $K$: per mode it is $\approx(Ks)^2G_0^2/(G_0^2+1/\rho)^2$ with $s$ the rms
shear deviation, so it passes the 20 dB noise floor above $K\approx14$ rad/μm at 3 %
and above $K\approx4$ rad/μm at 10 %. Using the field-mean shear instead of the nominal one
helps at 20 dB and can hurt at 40 dB (bowl 3 %: 0.75 instead of 0.42). The Monte Carlo
check (200 objects, real noise) reproduces the exact posterior to within 2.5 % in every
row tested.

## 4. Where diversity does come from: translating the sample

Two exposures at half the photons each; between them the sample moves by $s$ across a
ramp field, so every point is seen at two local shears differing by $\varepsilon D_0 s/L$:

| case | shear difference | m = 1 lost | m = 2 lost | m = 3 lost | RMSE |
|---|---|---|---|---|---|
| uniform, single exposure | — | 3.06 | 3.06 | 3.06 | 0.2077 |
| ramp 3 %, $s=L/2$ | 12 nm | 2.75 | 2.18 | 1.71 | 0.1787 |
| ramp 10 %, $s=L/4$ | 20 nm | 1.98 | 1.26 | 0.87 | 0.1446 |
| ramp 10 %, $s=L/2$ | 40 nm | 1.58 | 0.84 | 0.53 | 0.1244 |
| {0.8, 0.55} pair, half exposure each | 250 nm | 0.17 | 0.14 | 1.30 | 0.0987 |

The higher nulls benefit most, because the same shear difference moves them further. This
turns a calibrated field distortion into shear diversity with no change to the optics. It
needs a sample stage and a calibrated $\Delta(x)$, and the simulated field is periodic, so
every point has a partner position *(inference about practice)*.

## 5. What this does *not* establish

- **Projection level only.** No Ewald caps, Rytov factor or missing cone; the 3D lost-volume
  metric of [[multishear-fusion-design-solution]] was not recomputed. The separation
  $J=\rho S_tM$ fails for a space-variant shear, so a 3D version needs the space-variant
  operator inside the tomographic model *(inference)*.
- **Shear magnitude only.** $\Delta(x)$ varies in length along a fixed direction. A rotating
  shear direction, which distortion also produces, was not modelled.
- **Smooth, known profiles.** Linear and quadratic; real distortion maps, and the precision
  to which they can be calibrated, were not measured.
- **A periodic field.** The ramp jumps at the boundary; §2's 40 dB result may depend on it,
  and §4 assumes every point has a partner position.
- **Flat prior and one SNR pair** (20 and 40 dB).

## Reproduction

```bash
cd lab
python exp07_spatially_varying_shear.py    # ~90 s, < 1.5 GB
```

Prints every number above and writes `lab/out/varying_shear.png` (per-mode error of the
calibrated and naive estimators, the null's local track) and
`lab/out/varying_shear_images.png` (2D reconstructions).
