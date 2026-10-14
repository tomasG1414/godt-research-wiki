---
type: source
projects: [godt]
created: 2026-07-18
updated: 2026-07-18
status: current
source: "https://doi.org/10.1126/sciadv.adl0501"
---

# Single-Shot Deterministic Complex Amplitude Imaging with a Single-Layer Metalens

**Citation:** L. Li, S. Wang, F. Zhao, Y. Zhang, S. Wen, H. Chai, Y. Gao, W. Wang, L. Cao, Y. Yang, "Single-shot deterministic complex amplitude imaging with a single-layer metalens," Science Advances 10, eadl0501 (2024). Tsinghua University.

**Gist:** A single-layer metalens with spatial- and polarization-multiplexed point spread functions that records **four shearing interference patterns along both orthogonal in-plane directions in a single exposure**, using a commercial polarization camera. This is the strongest demonstration in the ingested corpus of simultaneous multi-direction shear acquisition — precisely the acquisition problem a multi-shear [[s42005-025-02416-3|GODT]] extension would face — though the reconstruction here is 2D and still integrates the gradient.

## Key results

- **Acquisition scheme (the transferable part).** A single metalens simultaneously produces two spatially-separated shearing interference images — one sheared along $x$, one along $y$ — each consisting of overlapped LCP and RCP copies. A polarization camera then splits each into four polarization channels (0°, 45°, 90°, 135°), which act as the four phase-shift steps. **One exposure therefore yields both shear directions and all four phase shifts**: $\nabla_x\varphi$ and $\nabla_y\varphi$ are recovered deterministically by four-step phase shifting, with no sequential acquisition, no moving parts, and no iterative algorithm.
- Deterministic 2D phase reconstruction from the full 2D gradient via higher-order finite-difference least-squares integration — determinism follows from having *both* gradient components, in contrast to single-direction systems ([[kwon2019]] and others) that need iterative algorithms and object priors to get 2D phase from a 1D gradient, "with limited reconstruction fidelity and speed."
- Performance: smallest measurable phase gradient 42 mrad/μm (compare [[kwon2019]]'s 92.3 mrad/μm), relative phase accuracy 0.0021λ. Metalens D = 2 mm, f = 1.5 cm, NA = 0.066; shearing distance Δs = 25 μm at magnification M = 1.5.
- **LED-compatible**, avoiding laser speckle — with the explicit condition that the *spatial coherence length must exceed the shearing distance* Δs. This is the same low-coherence/small-shear coupling that motivates GODT's design. Illumination is linearly polarized to avoid a random polarization phase difference between LCP and RCP.
- Design constraint worth noting: Δs must exceed twice the polarization camera's macro-pixel size (Nyquist), while smaller Δs is preferable for reconstruction sampling — a concrete tension between shear magnitude and sampling that any polarization-camera-based multi-shear scheme inherits.
- Demonstrated on surface metrology and dynamic complex-amplitude imaging of living cells.

## Method

Single-step-lithography single-layer metalens encoding spatial multiplexing (two shear directions at separated image locations) plus polarization multiplexing (LCP/RCP shearing pairs), read out by a commercial polarization camera. Requires prior calibration with a collimated beam to remove phase gradients from non-planar incident wavefront and alignment error.

## Relevance

Per [[godt]]: this is the most directly useful source in the corpus for the **acquisition side** of the multi-shear extension. The core difficulty with multi-shear GODT is cost: GODT already needs 5-frame temporal phase shifting per illumination angle, over 180 angles; naively adding shear directions multiplies that. This paper shows that both orthogonal shear directions *and* the full phase-shift set can be captured in a single exposure through combined spatial + polarization multiplexing on a polarization camera.

That is compatible with GODT's existing hardware in principle: [[s42005-025-02416-3]] already notes that "phase shifted interferograms can be also obtained in a single shot employing a pixelated polarization camera." The open piece is whether polarization-grating optics can additionally spatially multiplex *two shear directions*, as the metalens does. *(inference)*

The reconstruction side does **not** transfer: this paper integrates the gradient (least-squares integration) to recover 2D phase, which is exactly the step GODT is designed to avoid, and there is no tomography here at all. The contribution to borrow is acquisition, not inversion.

## Raises

- Can polarization gratings spatially multiplex two shear directions the way this metalens does, or would multi-shear GODT need a fundamentally different shearing element (metasurface, or two cascaded PG pairs at different orientations)? Unresolved by any ingested source; a hardware-feasibility question for the extension.
- The Nyquist constraint (Δs > 2× macro-pixel) interacts with GODT's shear-dependent transfer nulls — if multi-shear GODT wants specific shear *values* to place nulls favorably, the polarization camera's pixel geometry bounds what's achievable.
