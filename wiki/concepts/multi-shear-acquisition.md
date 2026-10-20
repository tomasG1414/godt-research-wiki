---
type: concept
projects: [godt]
created: 2026-07-18
updated: 2026-10-20
status: current
---

# Multi-Shear Acquisition

Capturing shearing-interferometry data along **more than one shear direction** (and/or at more than one shear magnitude). Central to this project's simulation thread: [[s42005-025-02416-3|GODT]] operates at a single shear and names combining multiple shear orientations as future work, since single-shear reconstruction gives orientation-dependent contrast.

## Why multi-shear

Two distinct motivations. The GODT paper states the first; the second comes from this project's simulations:

1. **Directional contrast (source-backed).** A single shear direction makes GODT sensitive to RI variation along that direction only, so features oriented parallel to the shear show weakly. [[kandel2019]] independently confirms the underlying cause — "the shear in GLIM renders the transfer function non-radially symmetric."

   A **second, distinct argument for direction diversity** arrives with [[2605-04297-double-sagnac-ls-qpm]]: integrating a single-direction differential phase map "introduces line artifacts in the reconstructed images," which their two-direction system removes. This is a *reconstruction*-side motivation — 1D integration accumulates error along its integration path, and a second orthogonal gradient constrains it — rather than a transfer-function one, and it is orthogonal to the saturation result below, which measured only accuracy in a global modal basis. It is also the weaker kind of claim for this project, since GODT does not integrate at all; but it is the reason the group itself gives for building two directions, and their single-Sagnac predecessor (arXiv:2604.20739) is where the artifact is documented.
2. **Filling transfer-function nulls — tested, confirmed** ([[finite-shear-null-experiment]]). $|H_\Delta(K)| = 2|\sin(K\Delta/2)|/\Delta$, zeros at $K = 2\pi m/\Delta$. At GODT's published experimental shear of 0.8 μm, **three nulls fall inside the measurable band** (7.85, 15.71, 23.56 rad/μm); with their surroundings that is 3.56 % of the band at 20 dB, and the null frequencies themselves are lost at any SNR. At the 0.1 μm shear used in the paper's *simulations* there are none.

   Two design rules fall out. **Magnitude ratios must not be ratios of small integers**: 0.8 + 0.4 barely helps because the nulls coincide, while 0.8 + 0.55 clears the band (the 3D study puts the best partner for 0.8 μm near 0.15–0.2 μm). And **magnitudes and directions fix different things** — magnitudes fill the $m\neq0$ comb, directions address the $K_x = 0$ blind line that no magnitude choice touches. Combining both gave an 8.5× RMSE improvement over single-shear 0.8 in a band-passed image test.

   **Redone in 3D against the real Ewald sampling, 2026-09-21** ([[multishear-fusion-design-solution]]). The true volume lost at the published 0.8 μm is **3.02 % of the sampled support**, not the 1D proxy's 3.56 %, and it is a different 3.02 %: **83.6 % of it sits on the null planes and 16.4 % on the blind plane**, which no magnitude choice touches. Two things the 1D calculation could not see. The 0.1 μm simulation shear has *no in-band nulls and still loses 2.54 %*, all blind plane, because a small shear has a wider low-$|G|$ slab around $K_x=0$; and 0.33 μm beats 0.1 μm despite having a null, for the same reason. **Single-shear quality is a trade between a null comb and a blind slab, not a monotone function of $\Delta$.**

   The incommensurate-ratio rule survives and was *re-derived* rather than imposed: scanning a second shear pair against a fixed 0.8 μm pair, the worst choices land exactly on ratios $\{1, \tfrac23, \tfrac34, \tfrac32, 2, 3, \tfrac35\}$. But it is **a rule for small channel sets only** — an optimised 12-magnitude design contains ratios 2, 3 and 4 and is unharmed, because a null shared by two channels is covered by the other ten.

### The shear magnitude should be about $2\pi/K_{\max}$

The cleanest design rule the 3D study produced, and the one a reader can act on. $2\pi/K_{\max}$ is the largest shear whose first null $2\pi/\Delta$ stays outside the measurable band; for GODT's simulation geometry that is **0.201 μm**, and an optimiser over a 0.05 μm candidate grid picks 0.20 μm as the best single channel. Larger shears admit nulls and smaller ones widen the blind slab, so the optimum sits *at* the boundary, $\Delta\approx2\pi/K_{\max}$. How this maps onto the published 0.8 μm experiment depends on the experimental illumination NA. The best orthogonal pair sits at 0.15 μm and stays there under every coherence model tested.

   Also from that experiment: the $K=0$ zero is a *separate* problem from the comb — it is the DC loss inherent to every gradient method, independent of $\Delta$, and is what [[s42005-025-02416-3]] presents as background suppression. Only the comb is fixable by magnitude diversity.

## Direction diversity saturates — the two motivations are not equal

Three papers from one metrology group measured what extra shear *directions* actually buy, holding magnitude fixed:

| directions | source | RMS gain |
|---|---|---|
| 2 → 4 | [[zhu2024-appl-phys-b]] (simulation) | 40.7 % |
| 2 → 3 | [[zhu2024-sensors]] (experiment) | ~0.7 % |
| 4 → 16 | [[zhu2025-sci-rep]] (experiment) | 0.0 % |

Gains are modest, saturate by about two orthogonal pairs, and are **error averaging** rather than recovery of missing information — [[zhu2024-appl-phys-b]] says so directly ("effectively mitigates the impact of random system errors"). That is theoretically unsurprising: one orthogonal pair already spans the full transverse gradient, so further directions are largely redundant.

The consequence for this project's argument: motivations 1 and 2 above are **not comparable claims**, and framing them as two halves of one case is a mistake. Directional contrast is a real but bounded effect fixable by one orthogonal pair. Null-filling recovers information that is otherwise structurally absent at any SNR. **[audit 2026-10-20]** An earlier version credited that to magnitude diversity alone and called direction diversity the narrower fix for the $K_x=0$ blind line. In 3D that is wrong: a second direction's data also cover most of the first direction's null planes (0.8 $x$ + 0.8 $y$ loses 0.09 % of the support against 3.02 % for 0.8 $x$; [[multishear-fusion-design-solution]] §3), and a second magnitude removes the residue where both directions are blind. The saturation discussed here is about adding directions beyond one orthogonal pair.

**Two reasons not to over-read the saturation result.** All three papers reconstruct smooth, low-order wavefronts in a **global Zernike basis**, which is largely immune to a transfer zero at any single spatial frequency, so null damage would not show up there. And none of them varies shear *magnitude*. They bound direction diversity; they say nothing about magnitude diversity. *(inference; as the authors note, their method "primarily addresses... low asphericity", and high-frequency objects are left to future work on "the frequency response of the shearing and un-shearing operators".)*

**The saturation is now reproduced from theory, in 3D, for a different reason than theirs** ([[multishear-fusion-design-solution]]). Opening an optimal design from two directions to four (adding 45° and 135°) improves the mean-error objective by **4 %**, even though the diagonals take half the time budget. The mechanism is not error averaging: one orthogonal pair already covers every blind plane, so further directions are close to redundant in the information sense too. Both routes lead to the same conclusion: beyond one orthogonal pair, extra directions add little, and the remaining gain comes from shear magnitude.

## The transfer nulls have 2D prior art — now confirmed

**Updated 2026-09-11 after reading [[servin2007]].** The null comb is established knowledge from multiple independent lines:

- **Servín, Cywiak & Dávila (2007)** — [[servin2007]] — the explicit Fourier-domain analysis: forward filter $2i\sin(dx\,\omega_x)\,\Phi$ (Eq. 6), inverse poles at $\sin(dx\,\omega_x)=0$ (Eq. 7), and a regularised LS inverse (Eq. 11) that is stable but cannot recover the lost information. Their $n$-direction generalisation (§7) allows variable magnitudes and directions, but addresses **spatial support** (boundary erosion on bounded pupils), not frequency-null filling. The authors call the frequency response "well known."
- **Mehta & Sheppard (2008)** — [[mehta-sheppard2008]] — derives the DIC effective object spectrum as $T_N = 2iT\sin(2\pi m\Delta - \phi)$, the same zero comb, tracing it to Cogswell & Sheppard (1992) and naming contrast reversal at zero crossings as an artifact source.

What [[finite-shear-null-experiment]] adds is the quantification at GODT's specific shear values (3 nulls in band at 0.8 μm vs. none at the 0.1 μm simulation shear) and the design rules in the tomographic setting. Combining coprime shear magnitudes is itself known in 2D (Elster & Weingärtner 1999 and follow-ups, reviewed in Tomioka et al., *Appl. Opt.* 2023, arXiv:2309.16434). In short: known 2D observations, quantified here in 3D.

Two things follow from Mehta-Sheppard that Servín does not cover. The DIC bias $\phi$ translates the whole comb ($2\pi m\Delta - \phi = k\pi$) — GODT discards this degree of freedom because its phase-shifting recovers the complex gradient directly, and whether bias diversity could substitute for or cheapen magnitude diversity is unexamined *(inference)*. And perpendicular to the shear, the DIC transmission cross-coefficient is *even*, i.e. brightfield — meaning a single shear transfers **zero** in-band phase information along that axis, not merely reduced contrast. That is firmer quantitative backing for motivation 1 than anything previously on this page.

## The acquisition-cost problem

GODT already needs 5 phase-shifted frames per illumination angle across ~180 angles. Naively adding shear directions multiplies this. Four strategies appear across the ingested sources:

| Strategy | Source | Cost / limitation |
|---|---|---|
| Sequential acquisition | baseline | Multiplies frame count by number of shears |
| Spatial multiplexing | [[sciadv-adl0501]] | Costs field of view or sensor area |
| Carrier-frequency multiplexing (QWLSI-style) | discussed in [[nguyen2025-deep-learning-lsi]]; **reopened** by [[2605-04297-double-sagnac-ls-qpm]] | Incompatible with GODT's **LED**, not with low spatial coherence as such — see below |
| Learned demultiplexing | [[nguyen2025-deep-learning-lsi]] | Needs training data; introduces a non-analytic step |
| Double birefringent crystal displacer (DBCs-BD) | [[zhu2025-sci-rep]], [[zhu2024-sensors]] | Passive and cheap, both orthogonal directions per exposure; but demodulation as published needs spatial coherence |
| Wavelength multiplexing (colour = shear magnitude) | [[colour-multiplexed-shear]] (computed); hardware in a Trusiak-group preprint, 2026-09 (abstract only) | Magnitudes only, not directions; needs resolved colours, so it forfeits the single broadband LED; assumes low sample dispersion; crosstalk costs 20–60 % |

**Wavelength multiplexing, added 2026-10-05.** The polarization-grating shear scales with $\lambda$, so two or three LED colours read by a colour camera are two or three shear magnitudes in one exposure. [[colour-multiplexed-shear]] finds it matches a mechanical magnitude pair, wants the two shortest wavelengths, and is unnecessary if the shear can simply be set to $2\pi/K_{\max}$. It supplies the *magnitude* diversity that the direction-only hardware in this table lacks, so the two combine.

**A possible constraint on mixing magnitudes.** [[mehta-sheppard2008]] shows that in DIC the coherence ratio $S = \mathrm{NA}_\text{cond}/\mathrm{NA}_\text{obj}$ must fall as shear grows — "for bigger shear, the condenser will need to be stopped down further." If that applies to GODT, a scheme mixing shear magnitudes might need different illumination apertures per channel. This bites only for Köhler-like geometries where shear is applied at the objective pupil; Nomarski-DIC escapes it via a matched condenser prism.

**Calibration is the real cost, not photons.** Converging evidence: shear error of 0.1 pixel already produces noticeable 2D phase error ([[2509-15382-osi-flex]]); polarizer azimuth error is amplified threefold, $(\Delta\varphi)_{\max}\approx3\Delta\theta$, with 0.2° already comparable to the measurement error ([[zhu2025-sci-rep]]). Both groups calibrate shear from the interferogram itself rather than trusting nominal values — [[zhu2024-appl-phys-b]] by edge-extraction on the sheared apertures, OSI-flex by joint estimation in the optimisation.

**Quantified from the information side, 2026-09-21** ([[multishear-fusion-design-solution]]). An rms shear error $s$ caps the per-voxel SNR at $|G|^2/(s K_\perp)^2$ no matter how many photons are collected. At exactly the 0.1 px error OSI-flex calls noticeable (10 nm at a 0.1 μm pixel — *inference*: OSI-flex's 0.1 px is a systematic offset on a 5.3 px shear with no stated physical pixel size), that ceiling is **16.5 dB at the band edge against a photon-limited 37.1 dB** — calibration is 20 dB below photons and is what sets high-frequency quality. Two consequences for a multi-shear design. Mean error roughly doubles between perfect calibration and 0.1 px. And the *optimal response* is to spread the time budget across every available channel rather than concentrate it, which is the opposite of the photon-limited optimum — so a calibration-limited instrument wants more channels, not fewer.

**And shear is not one number per channel.** [[2605-04297-double-sagnac-ls-qpm]] is the corpus's clearest measurement of how badly this misbehaves in a real two-direction instrument, and all three failure modes are worth expecting:

- **Anisotropy between channels.** The zero-shear camera position differed by 0.7 mm between the $x$ and $y$ axes (5.7 vs 6.4 mm), so at any single camera position the two channels have genuinely unequal shear.
- **It propagates into reconstruction.** 2D integration with equal $dx$ and $dy$ elongated the reconstructed beads; a $dx\!:\!dy$ ratio of **1.5:1** was needed for correct sphericity, and every result in the paper is processed that way. The integration step must be told the true shear ratio — assuming isotropy is a reconstruction error, not just a calibration one.
- **It varies across the field of view.** Tube-lens distortion made shear non-uniform spatially, visibly degrading one side of the image, and this — not nominal resolution — drove their choice of operating point. A single scalar $\Delta$ per channel is therefore an approximation even after calibration.

The last point is the one that bites a *tomographic* multi-shear scheme hardest: [[finite-shear-null-experiment]]'s null positions $K = 2\pi m/\Delta$ move with $\Delta$, so a spatially varying shear smears the null comb across the field rather than placing it at known frequencies. Whether that smearing helps (nulls average out) or hurts (no clean inverse filter anywhere) is unexamined and is a cheap thing to simulate in `lab/`.

**The carrier-frequency row was too strongly stated (corrected after reading [[2605-04297-double-sagnac-ls-qpm]]).** That paper runs single-shot Fourier-transform demodulation on a *two-carrier mesh* interferogram encoding both orthogonal shear directions, under low spatial coherence, and still reaches 11.7/13.3 mrad spatial sensitivity. The binding requirement for dense carrier fringes is **temporal** coherence, not spatial: the two interfering beams are sheared copies of the same field, so their fringe visibility is the source's degree of coherence at the shear distance, $|\mu(\Delta)|$, which stays high while the shear is inside the transverse coherence length. Their [[pseudo-thermal-illumination|pseudo-thermal source]] (laser + rotating diffuser) lowers spatial coherence only, keeping speckle suppression *and* fringe density.

So the route is closed for GODT's broadband LED specifically, and available with a source swap. Two things stop this from being a free win, both on [[pseudo-thermal-illumination]]: a quasi-monochromatic source forfeits the transfer-function flattening [[ledwig-robles2023]] attributes to broadband illumination, and the separate obstacle that FT demodulation breaks down under high-NA *conical tomographic* illumination ([[common-path-shearing-interferometry]]) is untouched by the source change — [[2605-04297-double-sagnac-ls-qpm]] never meets it because it illuminates on-axis. Single-shot multi-shear *tomography* therefore is not unlocked simply by adopting this source.

Two further single-shot routes exist that reduce per-shear cost without addressing multiplexing directly: pixelated polarization cameras (noted as viable in [[s42005-025-02416-3]] itself, and exploited by [[sciadv-adl0501]]) and Hilbert spiral transform demodulation ([[2024-06-14-598984-polarization-hilbert-holotomography]], demonstrated in total-shear mode only).

## Hardware feasibility

[[nguyen2025-deep-learning-lsi]] establishes that two shear directions can be built from **polarization gratings** — the same component family GODT uses — by adding a second PG + flat-mirror arm. Its modules work in double-pass reflection, whereas GODT's PG1/PG2 pair is transmissive; whether two transmissive PG pairs at different orientations can cascade without the first corrupting the second is unexamined.

**Cascading two shearing stages does work, as [[2605-04297-double-sagnac-ls-qpm]] demonstrates.** Not with polarization gratings, but with two Sagnac interferometers (polarizing beam splitter + two mirrors each) in series at orthogonal orientations, retrofitted to a bright-field microscope's output port. The enabling trick is a **45° polarizer between the stages**: stage 1 emits two orthogonally polarized $x$-sheared copies, the polarizer projects them onto a common state so they can interfere *and* so stage 2's PBS sees a 45° input it can split afresh. That resets the polarization basis and is what stops the first stage from consuming the degree of freedom the second one needs. Mirror tilt in each stage sets that stage's shear independently.

Two limits on how far this transfers. The architecture is **polarization-multiplexed and therefore two-stage-limited** — each stage consumes the polarization basis and hands it on exactly once, so a third orthogonal shear direction would need a different multiplexing axis. And a serial cascade may also produce diagonal cross-terms in the spectrum, at $(\pm d_1/2, \pm d_2/2)$ (*inference*, see that page's Raises).

[[sciadv-adl0501]] achieves both orthogonal directions plus all four phase shifts in a **single exposure** using one spatial- and polarization-multiplexed metalens with a polarization camera — the most acquisition-efficient demonstration in the corpus, though 2D and non-tomographic.

## Existing tooling

[[2509-15382-osi-flex]] — from the GODT authors — is an open-source (github.com/QCI-LAB/OSI-flex) optimization framework that already handles **arbitrary numbers, magnitudes, and orientations of shear vectors**, jointly estimating phase and shear via automatic differentiation + ADAM. It integrates gradients to 2D phase rather than reconstructing in 3D; its architecture (data-fidelity summed over shear channels, differentiable forward model, joint shear estimation) is a plausible template for a multi-shear tomographic reconstruction with the forward model swapped for the Ewald-sphere relation.

Practical warning from that work: a shear error of **0.1 pixel** already causes noticeable phase error in 2D, and cross-correlation shear estimation is only good to ~0.5 pixels. A multi-shear tomographic scheme should probably estimate shear jointly rather than trust nominal values.

## Precedent for combining shear data

- **2D**: DIC-family work (Cogswell, Shribak, cited in [[nguyen2017]]) long established that two orthogonal shears recover the full transverse gradient. [[sciadv-adl0501]] uses exactly this to make 2D phase reconstruction *deterministic*, where single-direction systems need iterative algorithms plus object priors.
- **3D**: no ingested source combines multiple shears tomographically. Whether the 2D orthogonal-shear math carries into the Ewald-space setting is open — in 3D each shear applies a distinct $\frac{e^{j\Delta K}-1}{\Delta}$ filter to the *same* scattering-potential spectrum, rather than sampling independent gradient components.

## How to fuse the channels — three options now on the table

The reconstruction question is what this project studies, and the corpus contains three distinct answers to "how do you combine shear channels", worth naming explicitly:

1. **Modal** — stack each direction's differential-Zernike design matrix as extra rows in one global least-squares problem ([[zhu2024-appl-phys-b]]). Null-immune, because a global basis compensates a zero at any single frequency, but capped at $N$ coefficients and therefore incapable of high-frequency content. Probably wrong for textured biological samples; right for smooth optics.
2. **Optimisation** — sum a data-fidelity term over shear channels in a differentiable forward model, estimating shear jointly with phase ([[2509-15382-osi-flex]]). Currently 2D, but architecturally the closest template.
3. **Spectral** — combine in $\tilde O$-space, which is where GODT natively works and where the $(e^{j\Delta K}-1)/\Delta$ filters live. None of the ingested sources does this in 3D. Formulated in [[multishear-fusion-design-problem]] and **solved in [[multishear-fusion-design-solution]]**: the per-voxel Fisher information separates as $J = \rho\,S_{\mathbf t}(\mathbf K_\perp)\,M(\mathbf K)$, so the optimal weights are inverse-variance in $\tilde O$-space and the shear design decouples from the tomography entirely.

One transferable empirical hint from [[zhu2024-appl-phys-b]]: averaging *Zernike coefficients* across direction sets beat averaging *differential wavefronts* by 6× (ΔRMS 0.0001 λ vs 0.0006 λ) on identical raw data. Fuse in the representation you are solving for, not in the measurement domain — the analogue here is fusing in $\tilde O$-space rather than averaging reconstructed gradient volumes.

**Demonstrated in 3D: 21× over equal-weight averaging.** On a textured cell phantom at equal photon budget, spectral (Wiener) fusion of a {0.8, 0.55} μm $x$+$y$ design reaches 0.0078 relative RMSE against **0.1607** for an equal-weight average of the four per-channel reconstructions ([[multishear-fusion-design-solution]]). A per-voxel precision-weighted average of the same reconstructions gives 0.0079: fusion *is* inverse-variance averaging in $\tilde O$-space ([[multishear-fusion-derivation]] §6.2). So the 21× measures the value of per-frequency weights, not a new estimator, and on [[multiple-scattering-check]]'s phantom the margin is 9.6× (linear) → 1.8× → 1.1× as realism is added.

The estimator-side statement of the same fact: the prior-free BLUE **diverges** for a single shear, because it is undefined where $J=0$, and is perfectly well-conditioned for an incommensurate pair, where it agrees with Wiener to four digits. Multi-shear is the difference between an estimable and a non-estimable parameter, not a noise-averaging improvement.

## Relevance

Per [[godt]]: this page collects what the corpus says about the project's core extension. The reconstruction question (how to fuse multi-shear data in $\tilde O$-space) is the question this project studies; the acquisition question is largely solved in principle by existing hardware, with the caveats above.
