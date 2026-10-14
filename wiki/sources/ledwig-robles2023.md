---
type: source
projects: [godt]
created: 2026-07-19
updated: 2026-07-19
status: current
source: "https://doi.org/10.1063/5.0123206"
---

# Partially Coherent Broadband 3D Optical Transfer Functions

**Citation:** Patrick Ledwig, Francisco E. Robles, *APL Photonics* **8**, 041301 (2023). DOI `10.1063/5.0123206`.

**Gist:** The theoretical follow-on to [[optica-8-1-6]] from the same Georgia Tech group. It gives a closed method for computing the non-paraxial 3D optical transfer function of a microscope under an *arbitrary* source — any angular distribution, any spectrum, no symmetry assumed — and makes it cheap: $O(n^3)$ instead of $O(mn^4)$, 20 seconds instead of 1 h 20 min at $256\times256\times128$. For [[godt]] this is the tool that would answer what GODT's transfer function actually is under its real broadband LED, rather than under the monochromatic coherent idealisation its [[rytov-approximation|Rytov]] reconstruction assumes.

## Key results

**Spatial and temporal coherence unified into one 3D angular spectrum.** The paper's conceptual move is to map temporal frequency onto the *radial* coordinate of a 3D reciprocal space $\vec u$, so direction and wavelength become one object $S(\vec u)$. That yields a 3D broadband analogue of the van Cittert–Zernike theorem,

$$W(\vec x, \nu) = \int S(\vec u)\, e^{i2\pi \vec x\cdot\vec u}\,\delta(|\vec u| - \nu)\, d^3\vec u \qquad \text{(Eq. 1)}$$

in which the cross-spectral density depends on three spatial coordinates rather than six, because of shift invariance in the local frame.

**The Ewald sphere becomes an Ewald plane.** The 3D phase/amplitude transfer function (Eq. 3) sifts the source with delta functions whose argument, after the coordinate change, describes a *plane* rather than a sphere:

$$\tilde T_{\alpha/\phi}(\vec q) = \int S(\vec u)\left[|P(\vec u - \vec q)|^2\delta(|\vec u| - |\vec u - \vec q|) \pm |P(\vec u+\vec q)|^2\delta^*(|\vec u+\vec q| - |\vec u|)\right] d^3\vec u$$

The plane collects every source field pair that can contribute incoherently to image spatial frequency $\vec q$, *regardless of wavelength*. Amplitude signal is the sum of the two sampling planes, phase signal the difference. Compare the monochromatic coherent case in [[optical-diffraction-tomography]], where the sampling surface is the Ewald sphere cap.

**The speed comes from analytic geometry, not approximation.** Source and pupil acceptance are cones in $\vec u$; spectral limits are spheres. Ewald planes cut these in hyperbolas and circles, so the integration bounds can be written down in closed form (Eqs. 4–5) and only non-zero regions are ever summed (Eq. 6). No 3D FFTs, no convolutions. Cost decouples from image resolution and re-couples to the smoothness of the source in k-space — and most real sources are smooth there.

**Broadband illumination flattens the transfer function.** In the demonstration — 6000 K white light Monte-Carlo-scattered through brain tissue, arriving with angle-dependent colour "like a dawn-lit sky" — the broadband OTF's main lobes are noticeably flatter than the single-wavelength OTF at the mean wavelength (Fig. 4c vs 4d). The paper argues this is advantageous for deconvolution: more uniform frequency sensitivity, less damping of the object's high-frequency content.

**Assumes slowly varying, not weak, objects.** Stated as nearly equivalent in practice to the weak-object transfer function, but additionally admitting dark-field terms where source NA exceeds pupil NA, and justified on the [[rytov-approximation|Rytov]] approximation as more valid for biological samples — the same reasoning GODT uses.

## Method

Analytic derivation (details in supplementary material), then GPU-parallelised numerical integration on the sampling planes. Benchmarked against the conventional convolutional approach of discretising the source and accumulating spherical-cap convolutions per incident direction, on the same hardware (i7-7800X + GTX 1080 Ti). Validated against the quasi-monochromatic case in the supplement. Example code promised on the lab website — **not** a named public repository, unlike [[2509-15382-osi-flex]].

## Relevance

Per [[godt]]:

- **GODT's illumination is broadband and partially coherent; its reconstruction theory is not.** GODT uses an LED specifically because low coherence suppresses speckle and multiple-scattering background, then reconstructs with a coherent monochromatic Rytov model. That mismatch is normal practice and probably benign, but it is unquantified in this wiki, and this paper supplies the machinery to quantify it — compute the true broadband 3D transfer function for GODT's actual source and pupil and compare against the idealisation. That is a well-defined, tractable numerical experiment for `lab/`, and a natural sequel to [[finite-shear-null-experiment]].
- **It may partly answer the shear-versus-coherence tension** raised by [[mehta-sheppard2008]]. If the illumination cone must narrow to keep coherence across a large shear, this framework computes the resulting 3D transfer function directly, including the [[missing-cone-problem|missing cone]] consequences, rather than requiring the two effects to be reasoned about separately.
- **The flattening result cuts in GODT's favour.** If broadband illumination genuinely evens out the transfer function, GODT's LED is buying more than speckle suppression. Not demonstrated for a shearing geometry — the demonstration here is epi-mode qOBM-like — so this is a hypothesis, not a result. *(inference)*
- Confirms [[optica-8-1-6]] was not a one-off: this group is systematically building non-paraxial partially coherent 3D transfer-function theory, and is the corpus's other main programme in this theory.

## Raises

- Streibl, *Opt. Acta* **31**, 1233 (1984) — named here as the landmark this result mirrors; the origin of the 3D partially coherent phase/amplitude transfer function, and a foundational citation alongside [[wolf1969]].
- Soto, Rodrigo & Alieva, *Opt. Lett.* **43**, 4699 (2018), "Partially coherent illumination engineering for enhanced refractive index tomography" — cited (ref 31) for the flattening claim above. If illumination engineering can shape the effective transfer function, it is a *third* lever alongside shear magnitude and shear direction.
- Does the Ewald-plane formulation compose with a shear filter at all? GODT's $(e^{j\Delta K_x}-1)/\Delta$ is derived in a coherent monochromatic setting; whether it survives as a simple multiplicative factor under this broadband treatment is unexamined and is the crux of using this machinery for GODT.
