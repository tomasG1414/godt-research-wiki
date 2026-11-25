---
type: project
projects: [godt]
created: 2026-07-18
updated: 2026-11-25
status: current
---


# GODT: Gradient Optical Diffraction Tomography

> Reading on GitHub? Names in double square brackets, such as [[finite-shear-null-experiment]],
> are links between wiki pages. They work in Obsidian, and here each one is the file of that
> name under `wiki/` (in this case [`notes/finite-shear-null-experiment.md`](../notes/finite-shear-null-experiment.md)).
> The notes behind the results are [blind spots](../notes/finite-shear-null-experiment.md),
> [choosing and combining shears](../notes/multishear-fusion-design-solution.md),
> [realistic physics](../notes/multiple-scattering-check.md), and its follow-ups with
> [an exact referee](../notes/full-wave-referee.md) and
> [no gridding](../notes/gridding-and-forward-model.md).

## Goal

Understand gradient ODT ([[s42005-025-02416-3]], *Communications Physics* 8, 510 (2025);
preprint [[2411-08423-godt]]) in depth, then study one extension numerically: acquiring
and fusing several shear channels (different shear magnitudes and directions) in a
single tomographic reconstruction. Everything here runs on synthetic data, with no lab needed.

## Background

GODT takes common-path [[common-path-shearing-interferometry|shearing-interferometry]] phase-gradient data and feeds it straight into a first-order [[rytov-approximation|Rytov]] tomographic reconstruction (the Ewald-sphere / direct-inversion pipeline of [[optical-diffraction-tomography]]). It does not integrate the gradient into a phase map first. Its central result is $\tilde{O}_\Delta = iK_x\tilde{O}$: reconstructing from shear data gives the exact derivative of the ODT scattering-potential spectrum, so ODT's diffraction compensation carries over unchanged. The published version extends this to finite shear with the exact factor $\frac{e^{i\Delta K_x}-1}{\Delta}$ (the paper writes $j$ for $i$), which tends to $iK_x$ as $\Delta\to0$. The authors validated it numerically on a [[shepp-logan-phantom]] and experimentally on a [[nano-printed-cell-phantom]] and hippocampal neurons, using a polarization-grating common-path system.

Total-shear common-path ODT needs a sample-free region as its reference, so it fails on dense, confluent samples. At small shear GODT needs no such region and works with low-coherence LED illumination, which is how it images confluent samples (neurons) that total-shear mode cannot.

The gradient *measurement* is shared with the whole [[gradient-phase-imaging]] family ([[nguyen2017]], [[kandel2019]], [[kwon2019]], [[optica-8-1-6]]). Reconstructing a 3D index derivative straight from gradient data also has precedent with straight-ray tomography ([[vishnyakov2016]]) and in X-ray phase tomography ([[guizar-sicairos2011]]). GODT does it inside a diffraction framework, which is the distinction its published version draws itself.

## Why several shear channels

See [[multi-shear-acquisition]] for the full picture. There are two reasons.

1. **Directional contrast.** A single shear sees index variation along one direction only ([[kandel2019]]; [[mehta-sheppard2008]] shows the transfer perpendicular to the shear is brightfield). Measured gains from extra directions saturate quickly, at about one orthogonal pair ([[zhu2024-appl-phys-b]], [[zhu2024-sensors]], [[zhu2025-sci-rep]]). Combining shear orientations is also the future work the GODT paper names.
2. **Transfer nulls.** The finite-shear factor vanishes at $K=2\pi m/\Delta$. The null comb, and filling it with a second shear in a coprime ratio, are known in 2D lateral-shearing interferometry ([[servin2007]], [[mehta-sheppard2008]]; Elster & Weingärtner 1999, summarised in Tomioka et al. 2023, arXiv:2309.16434). This project quantifies the comb in the 3D Ewald setting ([[finite-shear-null-experiment]]). In 3D a second shear direction already covers most of the first direction's null planes, and a second, incommensurate magnitude removes what both directions miss ([[multishear-fusion-design-solution]] §3). The null frequencies carry no data at any SNR, so this is recovery of missing information rather than noise averaging.

Acquisition hardware for two shear directions already exists in several forms ([[nguyen2025-deep-learning-lsi]], [[sciadv-adl0501]], [[2605-04297-double-sagnac-ls-qpm]]), and colour can act as a shear magnitude ([[colour-multiplexed-shear]]). The question studied here is the reconstruction side: how to fuse shear channels in $\tilde O$-space, and how to choose them.

## Results of the simulation thread

Code lives in [`lab/`](../../lab/README.md). The results below come from deterministic or seeded scripts there. All of them use the GODT paper's simulation optics ($\lambda=0.5$ μm, NA 1.33, 60° illumination) and monochromatic light unless stated.

### Step 1: transfer nulls ([[finite-shear-null-experiment]])

At 0.8 μm shear, three nulls fall inside the measurable band and cost 3.56 % of it along the shear axis. At 0.1 μm there are none. Two design rules follow: avoid commensurate magnitude ratios, and use both kinds of diversity, since a second direction covers the blind plane and most of the null planes and a second magnitude removes what both directions miss.

### Step 2: fusion and optimal design in 3D

Notes: [[multishear-fusion-design-problem]], [[multishear-fusion-design-solution]], derivation in [[multishear-fusion-derivation]].

1. The information factorises exactly as $J(\mathbf K)=\rho\,S_{\mathbf t}(\mathbf K_\perp)\,M(\mathbf K)$, a consequence of the paper's Eq. 20 (the shear factor multiplies the spectrum point by point), and the cap density $M$ has a closed form. The shear design enters only through $S_{\mathbf t}$ and the tomographic geometry only through $M$, so the whole study takes minutes of compute.
2. Transmission ODT never samples the $K_z$ axis, so the fact that in-plane shears cannot see it costs nothing extra. What a single shear does lose is its blind *plane*, 0.50 % of the support, and a second shear direction recovers it.
3. At 0.8 μm the lost volume is 3.02 %: 83.6 % on the null planes and 16.4 % on the blind plane. A second direction at the same magnitude cuts it to 0.09 %, a second magnitude along the same direction to 0.50 %, and both together to 0.00 %. A larger shear narrows the blind slab but lets nulls into the band, so single-shear quality is not monotone in $\Delta$.
4. The best single shear is $\Delta\approx2\pi/K_{\max}=\lambda/(\mathrm{NA_{det}}+\mathrm{NA_{ill}})$ (0.20 μm for the simulation optics), the largest shear that keeps the null comb out of the band and roughly the lateral resolution limit. This assumes coherent, monochromatic light; partial coherence pulls the optimum to smaller shears.
5. In the linear model and at equal photon budget, weighting each frequency by its information gives 21× lower error (textured cell phantom) than equal-weight averaging of the reconstructions. Fusion is the same as precision-weighted averaging, so the estimator is standard; the 21× measures what the per-frequency weights are worth against plain averaging.
6. Calibration may limit the high frequencies. With an assumed 10 nm rms random shear error, the band-edge SNR ceiling is about 20 dB below the photon limit. A broadband LED already partly fills the higher-order nulls.

### Colour as a shear channel ([[colour-multiplexed-shear]])

Three LED colours act like a mechanical pair of shear magnitudes, blue + green is the best pair, and at the design-rule shear the trick is unnecessary. The separation theorem breaks across colours, and the missing cone is untouched.

### Step 3: multiple scattering ([[multiple-scattering-check]])

The same estimator, run on data from a [[brenner1993]] WPM forward model:

1. The null planes are set by the shear, whatever the scattering. The shear difference acts linearly on the measured Rytov field after all propagation, so the data vanish on the same planes for any scattering. The 4e-16 check confirms the implementation (on a periodic grid it is the DFT shift theorem). Wrapping of the phase difference, a finite field of view, broadband (dispersive) shear and spatially varying shear all blur the nulls, and a nonlinear inverse could in principle recover weak information there.
2. The margins shrink with contrast. Fusion beats averaging by 9.6× under the ideal linear model (on that note's cell phantom), 1.8× with realistic gridding and 1.1× with multiple scattering at [[chen1998]]'s validity bound. Over the same steps, the best design's advantage over a single 0.8 μm shear goes from 14.3× to 2.5× and 1.3×. (Revisited in step 4.)
3. Shear data never form the absolute phase, so they avoid the global unwrapping step whose failures end first-Rytov ODT's range in [[chen1998]]. They do not avoid the Rytov linearisation, and the shear difference has its own wrap limit (2.2 % of pixels at 0.8 μm vs 0.45 % one-pixel residues). A comparison against properly unwrapped ODT has not been run.

### Step 4: an exact referee, and the two implementation errors removed

Notes: [[full-wave-referee]] (exp09), [[gridding-and-forward-model]] (exp10, exp11). Step 3 left the forward-model gap provisional and the gridding error unfixed. Both are now settled.

1. **An exact referee.** A [[convergent-born-series]] solver of the full scalar Helmholtz equation matches the exact Mie series for a 1.5 μm sphere to 0.4 % (0.9 % at 60°), where first Born is off by 16–80 %.
2. **The WPM–Born gap was the window.** On one simulation window WPM's linear limit is first Born × $k_z/k_{iz}$ to 0.1 % in every angle bin ([[multi-slice-forward-models]]). Step 3's extra steep-angle deficit came from comparing a cropped, absorbing-edge window with a periodic one. A 12.8 μm crop moves even the exact solution's linear part by 10 % (54 % on near-axis rays).
3. **Most of the realistic-gridding loss was the gridding.** An estimator with no Fourier-space gridding (an exact solve per $K_z$ column) takes the best design from 0.1058 to 0.0280 on continuum Born data. The margins come back to 8.8× over the single shear and 5.9× for fusion over averaging (ideal model: 14.3×, 9.6×).
4. **Multiple scattering then sets the error.** Step 3's model D had the size of the multiple-scattering correction about right but not its structure. Multi-slice first Born reproduces the exact excess to 4–10 % of the signal up to $a\,\delta n=0.46\lambda$. With it, the margins are 1.3× and 1.1× at $0.23\lambda$ and about 1× at $0.46\lambda$, with either estimator. Step 3's verdict stands; its absolute errors were 7–24 % too high.

## Limitations

- Everything is simulated. The experiment in the paper used a 530 nm LED with 35 nm bandwidth and an NA 1.45 objective, not the simulation optics used here, and the LED's bandwidth partly fills the higher-order nulls ([[multishear-fusion-design-solution]] §3, §6).
- The fringe-visibility model $\mu(\Delta)$ and the 10 nm shear-calibration error are assumptions, not measurements.
- The multiple-scattering numbers rest on multi-slice first Born, checked against an exact solver on one illumination and only on near-axis rays; steep rays cannot be adjudicated with the windows used, and nothing is validated at $a\,\delta n=0.92\lambda$ ([[full-wave-referee]]).
- The exact estimator is a MAP estimate under an oracle prior (the true object's shell spectrum), not a practical regularised reconstruction.
- Validation reports RMSE on cell-like phantoms only; SSIM, axial elongation, the Shepp–Logan phantom and the nano-printed phantom are not reported.
- Each experiment note lists the rest under "What this does not establish".

## Next technical steps

- **An inverse that models multiple scattering.** Step 4 makes this the item with the largest payoff: at $a\,\delta n\approx0.25\lambda$ multiple scattering, not gridding or design, sets the error. Multi-slice Born inside the reconstruction loop ([[ma2017]]'s route with a Born-consistent model) has the current fusion estimator as its first iterate.
- An unwrapped-ODT baseline for the wrapping comparison.
- ~~Spatially varying shear~~: done in [[spatially-varying-shear]]. A smooth, calibrated $\Delta(x)$ moves the null rather than removing it (no gain or loss at 20 dB); an estimator that ignores it degrades quickly (3 % ramp: RMSE 0.21 → 0.28), so the field must be calibrated to nanometres; translating the sample across the shear gradient gives genuine shear diversity.
- ~~Bilinear (or NUFFT) gridding~~ and ~~a Born-series forward model~~: done, step 4.
- Partial coherence: whether the shear factor stays a simple multiplier under [[ledwig-robles2023]]'s broadband Ewald-plane treatment.
- ~~Real samples on the $a\,\delta n/\lambda$ axis~~: done from the literature in [[real-samples-on-the-scattering-axis]]. Bacteria and isolated organelles sit below the thresholds, single cells straddle them (and the axis has a max- and a mean-contrast reading a factor 3.6 apart), GODT's printed phantom sits at 0.15–0.30, and its neurons cannot be placed from the paper. A uniform-sphere run would pin the threshold for nearly uniform objects.
