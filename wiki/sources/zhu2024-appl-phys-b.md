---
type: source
projects: [godt]
created: 2026-07-19
updated: 2026-07-22
status: current
source: "https://doi.org/10.1007/s00340-024-08353-3"
---

# Wavefront Reconstruction Based on Multi-Directional Orthogonal Lateral Shearing Interferometry

**Citation:** Yahui Zhu, Ailing Tian, Hongjun Wang, Bingcai Liu, *Applied Physics B* **130**, 218 (2024). DOI `10.1007/s00340-024-08353-3`. Received 25 Jul 2024, published 13 Nov 2024.

**Gist:** The reconstruction-theory paper of the Zhu group trio, alongside its two application papers, [[zhu2024-sensors]] and [[zhu2025-sci-rep]]. Its contribution is a differential-Zernike least-squares formulation that accepts shear directions at **arbitrary** angles rather than requiring alignment to $x$ and $y$, which is the practically useful part; and it supplies the middle data point showing how multi-direction gains scale.

## Key results

**Direction-independence is the actual claim.** Simulation on a shallow parabolic wavefront (PV 3.1222 λ, RMS 1.3211 λ) reconstructed from any two orthogonal directions — 0°/90°, 15°/105°, 30°/120°, 45°/135°, 60°/150° — gives residual ΔRMS between $1.676\times10^{-5}$ and $1.681\times10^{-5}$ λ (Table 1). Identical to three digits regardless of orientation. The point is that misaligned shear no longer costs accuracy, not that more directions help.

**Where the gain saturates.** Four directions give ΔRMS $9.894\times10^{-6}$ λ against $1.676\times10^{-5}$ λ for two — the paper's claimed **40.73 %** RMS improvement (Table 2). Set beside its companions this completes a consistent picture across all three papers:

| directions | source | RMS gain |
|---|---|---|
| 2 → 4 | this paper (simulation) | 40.7 % |
| 2 → 3 | [[zhu2024-sensors]] (experiment) | ~0.7 % |
| 4 → 16 | [[zhu2025-sci-rep]] (experiment) | 0 % |

Gains are real up to about two orthogonal pairs and vanish after that. The mechanism is error averaging — the paper says so directly, "effectively mitigates the impact of random system errors" — not recovery of missing information.

**How you fuse matters more than how many you take.** Averaging *Zernike coefficients* across direction sets beats averaging *differential wavefronts*: relative ΔRMS against a ZYGO reference is 0.0001 λ versus 0.0006 λ (Table 5), a sixfold difference, from identical raw data. Experimental reconstruction of a 50 mm parabola (0.5732 μm asphericity, 400 mm radius) gives PV 0.2590 λ / RMS 0.0398 λ by coefficient averaging against ZYGO's 0.2486 λ / 0.0397 λ.

**The authors note the high-frequency limitation.** The method "primarily addresses the high-precision reconstruction of aspherical wavefronts with **low asphericity**"; for wavefronts "which involve extensive high-frequency information, it is crucial to consider the retrieval of both sheared-phase and un-sheared phases," left to future work on "the frequency response of the shearing and un-shearing operators." That question is closely related to the transfer nulls of [[finite-shear-null-experiment]].

## Method

Differential Zernike polynomial fitting (DZF, after Liu 2003), which avoids the shear-matrix computations that the elliptic and numerical orthogonal transforms require. Each shear direction $\theta_i$ contributes rows $\Delta Z^{x,\theta_i}_j = Z_j(x - s\cos\theta_i, y - s\sin\theta_i) - Z_j(x,y)$ (Eq. 16); direction sets are stacked into one overdetermined system $\Delta W = \Delta Z a$ and solved by least squares $\hat a = (\Delta Z^T \Delta Z)^{-1}\Delta Z^T \Delta W$; coefficients are then averaged over sets (Eq. 35). Phase per direction comes from FFT carrier demodulation. Shear magnitude and direction are calibrated from the interferogram itself by shear-wavefront feature extraction — edge detection on the two sheared apertures, averaged over many sampling points (Eqs. 39–41). Spectral leakage from the finite aperture is handled by extending the fringe pattern 5 % beyond the interference region before filtering. Hardware is the same DBCs-BD two-calcite displacer described in [[zhu2024-sensors]].

## Relevance

Per [[godt]] and [[multi-shear-acquisition]]:

- **It settles the direction-count question quantitatively.** Together with its two companions: two orthogonal pairs are worth having, more is not. A GODT multi-shear argument resting on direction count would be arguing for a ~40 % noise-averaging gain that saturates immediately; the argument resting on **magnitude diversity filling transfer nulls** is arguing for recovery of information that is otherwise structurally absent. These are not comparable claims and should not be presented as two halves of one motivation.
- **The coefficient-averaging-beats-wavefront-averaging result is a transferable design hint.** Fusing in the representation you are solving for, rather than in the measurement domain, gave 6× here. The analogue for GODT would be fusing shear channels in $\tilde O$-space rather than averaging reconstructed gradient volumes. *(inference — but it is the same structural argument, and it costs nothing to respect.)*
- **Its stated future work is close to this project's question.** "The frequency response of the shearing and un-shearing operators" for high-frequency objects is what [[finite-shear-null-experiment]] computed for GODT's filter. It is encouraging that a metrology group considers the same question important from the 2D modal side.
- **Confirms the modal caveat.** All three Zhu papers reconstruct smooth, low-order, Zernike-representable wavefronts. A global modal basis is insensitive to a transfer zero at a single spatial frequency, so null damage would not show up in these experiments, as the authors' own note on high-frequency objects suggests. Their saturation results therefore bound direction diversity but say nothing about magnitude diversity.
- Same LED incompatibility as its companions: FFT carrier demodulation needs spatial coherence.

## Raises

Two references bear most directly on the reconstruction question, both from Servín's group:

- **Servín, Cywiak & Dávila, *Opt. Express* **15**(26), 17805 (2007) — "Lateral shearing interferometry: theoretical limits with practical consequences."** Cited here as the reference for the frequency response of shearing/un-shearing operators; contains the 2D transfer-null analysis ([[servin2007]]).
- Servín, Malacara & Marroquín, *Appl. Opt.* **35**(22), 4343 (1996) — "Wave-front recovery from two orthogonal sheared interferograms." The original two-orthogonal-shear inversion.
- Rimmer & Wyant, *Appl. Opt.* **14**(1), 142 (1975) — "Evaluation of large aberrations using a lateral-shear interferometer having **variable shear**." The oldest variable-magnitude shearing reference found so far; worth checking whether magnitude diversity was already exploited in 1975 and for what.
- Jeong, Park, Ghim et al., *Opt. Laser Eng.* **154**, 107020 (2022) — flexible lateral shearing interferometry based on **polarization gratings** for surface metrology. Same component family as GODT, in the metrology lineage.
