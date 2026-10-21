---
type: source
projects: [godt]
created: 2026-07-19
updated: 2026-10-20
status: current
source: "https://doi.org/10.1038/s41598-025-86370-8"
---

# A Novel Multiple Directional Shearing Interferometry System with Synchronous Polarization Phase Shifting

**Citation:** Yahui Zhu, Ailing Tian, Hongjun Wang, Bingcai Liu, Kai Wang, *Scientific Reports* **15**, 2775 (2025). DOI `10.1038/s41598-025-86370-8`. Received 29 Sep 2024, published 22 Jan 2025.

**Gist:** The companion instrument paper to [[zhu2024-sensors]] from the same group — this one builds the hardware and validates it against a ZYGO interferometer, where the other applies it to aspheres. Its value here is a single number: going from **4 shear directions to 16 changes the reconstructed RMS by 0.0000 λ**. That is the sharpest available statement that direction diversity saturates, and it matters for how the [[multi-shear-acquisition|multi-shear]] extension is argued.

## Key results

**Direction count saturates, hard** (Table 1, spherical mirror, $\lambda$ = 632.8 nm):

| | ZYGO | 4 directional derivatives | 16 directional derivatives |
|---|---|---|---|
| PV (λ) | 0.2480 | 0.2957 | 0.2922 |
| RMS (λ) | 0.0390 | 0.0390 | 0.0390 |

Quadrupling the number of shear directions moves PV by 0.0035 λ and RMS not at all. The authors state it plainly: PV and RMS "remain almost unchanged after sixteen or more directional derivatives." Agreement with ZYGO is PV error λ/20, RMS error λ/100 (ΔPV 0.0448 λ, ΔRMS 0.002 λ).

**Single-frame four-step phase shifting, without a pixelated polarization camera.** A checkerboard 2D phase grating (20 μm period, 692 nm etch, >60 % into the four first orders) makes four copies of the sheared field; a 2×2 array of ordinary polarizers at 0°/45°/90°/135° imposes phase shifts of 0, π/2, π, 3π/2 (Eq. 15), read out by one CMOS frame, and $\Delta W = \arctan[(I_{135}-I_{45})/(I_0-I_{90})]$. Multi-direction comes from rotating the whole displacer module — polarizer + crystal + quarter-wave plate as one unit, deliberately, so the polarization state does not change with rotation. Eight directions demonstrated (0° to 315° in 45° steps).

**Polarizer azimuth error is amplified threefold**: $(\Delta\varphi)_{\max} \approx 3\Delta\theta$ (Eq. 28), "comparable in magnitude to the measurement error." Calibrated by fitting Lissajous figures to adjacent-polarizer intensity curves — a perfect circle means exactly 90° apart. Measured azimuth differences 89.1°, 90.8°, 91.2°, 89.3°, giving phase-shift error below 0.02 λ.

**Shear from a single uniaxial crystal**: $s = d\tan\theta_t$, maximised when the optic axis sits at 45° to the surface. With $n_o$ = 1.65578, $n_e$ = 1.48520, $d$ = 9.28 mm the shear is 1.005 mm — the same crystal as [[zhu2024-sensors]]. Measured shear in the interferogram: 36 pixels, recovered by binarising and circle-fitting the two sheared apertures. Sub-frame registration of the four quadrants uses phase correlation, verified against known injected offsets of (30,10), (−10,−20), (5,8) pixels.

## Method

Jones-calculus derivation of the phase-shifting scheme (Eqs. 1–16), crystal-optics design of the shear (Eqs. 17–19), then experiment: He-Ne 632.992 nm through spatial filter and collimator, 50 mm spherical mirror of 400 mm curvature radius under test, differential Zernike polynomial fitting to 36 coefficients, coefficients averaged across direction pairs. The paper evaluates each direction against a reference direction and selects the best, rather than fusing all directions into one estimate.

## Relevance

Per [[godt]] and [[multi-shear-acquisition]]:

- **Together with [[zhu2024-sensors]], this closes the question of what direction count alone buys.** Two independent experiments from the same group, on different optics, both find saturation. The multi-shear argument should therefore rest on **magnitude diversity** — which fills the transfer nulls quantified in [[finite-shear-null-experiment]] — and treat direction diversity as the narrower fix for the $K_x=0$ blind line, not as a second independent win. Both papers hold shear magnitude fixed and vary only direction, so neither is evidence against magnitude diversity; they simply do not test it. **[audit 2026-10-20]** Superseded in 3D: a second direction also covers most of the first direction's null planes (0.8 x + 0.8 y loses 0.09 % of the support vs 3.02 % for 0.8 x); see [[multi-shear-acquisition]].
- **Same caveat as its companion:** a smooth spherical wavefront fitted in 36 Zernikes cannot exhibit transfer-null damage, because a global modal basis is insensitive to a zero at one spatial frequency. The saturation result is about direction redundancy, not about nulls. *(inference)*
- **The rotating-module design is directly transferable.** Rotating polarizer + crystal + QWP as one rigid unit, rather than rotating components individually, addresses the same problem any multi-direction polarization-grating setup faces.
- **Alignment tolerance is a real cost.** The 3× azimuth-error amplification, and the fact that a 0.2° polarizer error is already comparable to the measurement error, sets expectations for any polarization-multiplexed multi-shear scheme — consistent with the 0.1-pixel shear-error sensitivity reported by [[2509-15382-osi-flex]]. Multi-shear is cheap in photons and expensive in calibration.
- As with [[zhu2024-sensors]], the He-Ne source and carrier-free four-step demodulation are not directly portable to GODT's LED, though the polarizer-array phase-shifting scheme itself needs no spatial coherence and could be.

## Raises

Three references here bear on the **reconstruction** question:

- Velghe, Primot, Guérineau et al., *Opt. Lett.* **30**(3), 245 (2005) — "Wave-front reconstruction from multidirectional phase derivatives generated by multilateral shearing interferometers." The canonical multi-direction fusion reference.
- Legarda-Saenz & Espinosa-Romero, *Opt. Eng.* **50**(4), 040501 (2011) — "Wavefront reconstruction using multiple directional derivatives and Fourier transform." A *Fourier-domain* multi-direction fusion, i.e. the same domain GODT works in, unlike the modal fitting used here.
- Legarda-Saenz, Flores & Lakatos, *Opt. Laser Technol.* **182**(A), 112053 (2025) — the same line extended to discontinuous phase objects.

Also re-flagged from [[zhu2024-sensors]]: Shanmugam et al., *Light: Adv. Manuf.* **3**(2), 193 (2022), variable shearing holography — the only source in either paper's bibliography that varies shear *magnitude*.
