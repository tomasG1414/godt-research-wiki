---
type: note
projects: [godt]
created: 2026-09-21
updated: 2026-09-21
status: current
---

# Optimal fusion and design of multi-shear GODT channels — problem formulation

> **Solved 2026-09-21 — see [[multishear-fusion-design-solution]].** Every conjecture
> and proposition below has been computed. Headline: P1 holds but is minor (the
> direction-*dependent* blind plane is the real loss), C1 (sparsity) is **false**,
> C2 and C3 are confirmed, and fusing in $\tilde O$-space beats equal-weight averaging
> of reconstructions by 21× (precision-weighted averaging is equivalent to fusion; see
> [[multishear-fusion-derivation]]). This page is retained as the formulation; the answers,
> the numbers and the revised claim list live on the solution page.

Goal: derive the per-voxel information that a set of
shear channels delivers under [[optical-diffraction-tomography|Ewald]] sampling, use it
to build the optimal fusion estimator, and turn "which shears, for how long" into an
optimal-design problem. Validated on phantoms, no lab needed. Extends
[[finite-shear-null-experiment]] and the "spectral" fusion option in
[[multi-shear-acquisition]].

Status of claims: the forward model restates [[s42005-025-02416-3]]; estimation and
design theory are standard (general knowledge, not source-backed in this wiki);
their combination for GODT and every conjecture below are *(inference)* until
computed.

---

## 1. Concepts

### A. Tomographic forward model
- **Scattering potential** $O(\mathbf r)=k_0^2[n_m^2-n^2(\mathbf r)]$ (sign/normalisation as in [[optical-diffraction-tomography]]); the unknown. RI from $n=n_m\sqrt{1-O/k_m^2}$-type relation.
- **Rytov field** $u_R=\ln(u/u_i)$; first-order [[rytov-approximation]] linearises $u_R$ in $O$.
- **Ewald cap** $C_i$: the spherical shell of object frequencies $\mathbf K=\mathbf k-\mathbf k_i$ reached by illumination $\mathbf k_i$, truncated by the detection NA.
- **Rytov mapping factor** $k_0^2/(2ik_z)$: links detected spectrum to $\tilde O$ on the cap; noise in $\tilde O$ scales with $k_z^2$.
- **Cap multiplicity** $m(\mathbf K)$: number of illumination caps passing through voxel $\mathbf K$.
- **Measurable band / missing cone**: union of caps; the unfilled cone around $K_z$ is the [[missing-cone-problem]].
- **Direct inversion (DI) / gridding**: placing cap samples on a Cartesian grid ([[pan1983]]).

### B. Shearing
- **Phase difference** $\psi=\varphi(\mathbf x+\boldsymbol\Delta)-\varphi(\mathbf x)$: what the interferometer measures ([[common-path-shearing-interferometry]]).
- **Shear vector** $\boldsymbol\Delta_j$ (transverse): magnitude and direction of channel $j$.
- **Shear transfer function** $G_j(\mathbf K)=e^{i\boldsymbol\Delta_j\cdot\mathbf K_\perp}-1$; $|G_j|^2=4\sin^2(\boldsymbol\Delta_j\cdot\mathbf K_\perp/2)$. Its $\Delta\to0$ limit gives GODT's $i\Delta K_x$.
- **Null comb**: zeros at $\boldsymbol\Delta_j\cdot\mathbf K_\perp=2\pi m$, $m\neq0$ — planes in 3D. Prior art in 2D: [[servin2007]], [[mehta-sheppard2008]].
- **Blind plane**: $m=0$ member, $\boldsymbol\Delta_j\cdot\mathbf K_\perp=0$; contains the DC loss common to all gradient methods.
- **Magnitude vs direction diversity**: magnitudes move the $m\neq0$ comb; directions rotate the blind plane ([[multi-shear-acquisition]]).
- **Incommensurate-ratio rule**: shear-magnitude ratios must avoid small-integer ratios or nulls coincide ([[finite-shear-null-experiment]]).
- **Chromatic shear**: polarization-grating diffraction angle scales with $\lambda$, so under a broadband LED $\Delta$ is a distribution, not a number *(inference)*.
- **Spatially varying shear**: $\Delta(\mathbf x)$ from optics distortion ([[2605-04297-double-sagnac-ls-qpm]]); breaks shift invariance.

### C. Noise and coherence
- **Interferometric phase noise** $\sigma_j^2\propto 1/(V_j^2 N_\text{ph})$: depends on visibility and photons, not on $|\Delta|$ (general knowledge).
- **Fringe visibility** $V_j=|\mu(\boldsymbol\Delta_j)|$, with $\mu$ the **complex degree of spatial coherence** of the illumination at separation $\Delta_j$. Couples shear size to noise under partial coherence ([[partially-coherent-transfer-function]]).
- **Photon / time budget**: total acquisition $T$; **time share** $t_j$ per channel, $\sum t_j=1$.
- **Channel covariance** $\Sigma$: nonzero off-diagonals for single-shot multiplexed acquisition (cross-talk).

### D. Estimation
- **Channel**: one (shear, illumination) measurement set; generalises to (shear, illumination, rotation, wavelength).
- **Fisher information per voxel** $J(\mathbf K)$: inverse of the best achievable variance at $\mathbf K$.
- **BLUE / inverse-variance weighting**: optimal linear unbiased fusion; weights $\propto |h|^2/\sigma^2$.
- **Wiener / MAP fusion** with spectral prior $S_O(\mathbf K)$: MSE $(J+S_O^{-1})^{-1}$.
- **Unmeasured set** $\mathcal U=\{J=0\}$: only a prior can fill it.
- **Analogues to borrow** (general knowledge, not yet in wiki): SIM generalised Wiener filter (Gustafsson 2000, 2008); multiview deconvolution (Preibisch 2014); Fourier gradient integration (Frankot–Chellappa 1988; Velghe 2005; Legarda-Saenz 2011).

### E. Design
- **Optimal experimental design**: choose measurements to optimise a scalar function of $J$.
- **A-optimal** (mean error), **D-optimal** (log-information / capacity), **E-/worst-case** (max–min information).
- **Convexity**: $J$ linear in $\mathbf t$ ⇒ A-, D-, worst-case problems are convex.
- **Design sparsity**: optimal designs often use few support points (Carathéodory-type results) — here a conjecture.

### F. Robustness and model mismatch
- **Shear calibration error** $\delta\Delta$: $\delta G_j = i(\mathbf K_\perp\!\cdot\!\hat{\boldsymbol\Delta}_j)e^{i\boldsymbol\Delta_j\cdot\mathbf K_\perp}\delta\Delta$; relative error diverges at nulls. 0.1 px already matters ([[2509-15382-osi-flex]]).
- **Uncertainty-aware weights**: signal-dependent effective noise.
- **Rytov breakdown**, **phase-only data** (amplitude neglected; cost quantified in [[krauze2016]]), **broadband vs monochromatic model** ([[ledwig-robles2023]]).

### G. Validation
- **Phantoms**: textured cell-like (high-$K$ content, where nulls live); [[shepp-logan-phantom]] (baseline, smooth); [[nano-printed-cell-phantom]] (GODT's experimental target).
- **Forward models**: linear Rytov (theory check); WPM ([[brenner1993]]) for realism.
- **Key diagnostic**: empirical error spectrum vs predicted $(J+S_O^{-1})^{-1}$.

---

## 2. Problem formulation

### 2.1 Notation

| symbol | meaning |
|---|---|
| $k_m=n_m k_0$ | wavenumber in medium |
| $\mathcal I=\{\mathbf k_i\}$ | illumination wavevectors, $|\mathbf k_i|=k_m$, $|\mathbf k_{i\perp}|\le k_m\mathrm{NA_{ill}}/n_m$ |
| $C_i$ | cap: $\{\mathbf K=\mathbf k-\mathbf k_i:\ |\mathbf k|=k_m,\ k_z>0,\ |\mathbf k_\perp|\le k_0\mathrm{NA}\}$ |
| $\mathcal I(\mathbf K)$ | illuminations whose cap contains $\mathbf K$ |
| $k_z(\mathbf K;i)$ | $\sqrt{k_m^2-|\mathbf K_\perp+\mathbf k_{i\perp}|^2}$ |
| $\mathcal D=\{\boldsymbol\Delta_1,\dots,\boldsymbol\Delta_N\}$ | candidate shear channels |
| $G_j(\mathbf K)$ | $e^{i\boldsymbol\Delta_j\cdot\mathbf K_\perp}-1$ |
| $t_j$ | time share of channel $j$ |
| $\mu(\boldsymbol\Delta)$ | complex degree of spatial coherence |
| $S_O(\mathbf K)$ | prior power spectrum of $O$ |

### 2.2 Assumptions (baseline; each relaxed later)

- **A1** First-order Rytov, weak scattering.
- **A2** Phase-only data with log-amplitude neglected (as GODT). Relaxing it couples $\mathbf K$ and $-\mathbf K$ (Hermitian mixing) → $2\times2$ blocks per voxel pair.
- **A3** Exact gridding: cap samples land on the voxel grid (interpolation error ignored).
- **A4** Additive circular Gaussian phase-difference noise, white in the detector plane; no wrapping.
- **A5** Shears known exactly and spatially uniform.
- **A6** Channels acquired sequentially ⇒ independent noise.
- **A7** Monochromatic illumination.

### 2.3 Measurement model

For illumination $i$ and shear channel $j$, at every $\mathbf K\in C_i$:

$$y_{ij}(\mathbf K)=h_{ij}(\mathbf K)\,\tilde O(\mathbf K)+n_{ij}(\mathbf K),\qquad h_{ij}(\mathbf K)=G_j(\mathbf K)\,\frac{k_0^2}{2i\,k_z(\mathbf K;i)},$$

$$n_{ij}\sim\mathcal{CN}(0,\sigma_j^2),\qquad \sigma_j^2=\frac{\sigma_0^2}{t_j\,T\,N_\text{ph}\,|\mu(\boldsymbol\Delta_j)|^2}.$$

Single-shear GODT is $N=1$ with $\Delta\to0$.

### 2.4 Estimation subproblem (design fixed)

Per-voxel Fisher information

$$J(\mathbf K)=\sum_{j=1}^{N}\ \sum_{i\in\mathcal I(\mathbf K)}\frac{|h_{ij}(\mathbf K)|^2}{\sigma_j^2}
=\sum_j t_j\,a_j(\mathbf K),\qquad
a_j(\mathbf K)=\frac{T N_\text{ph}|\mu(\boldsymbol\Delta_j)|^2}{\sigma_0^2}\,|G_j(\mathbf K)|^2\sum_{i\in\mathcal I(\mathbf K)}\frac{k_0^4}{4k_z(\mathbf K;i)^2}.$$

BLUE, defined where $J>0$:

$$\hat{\tilde O}(\mathbf K)=\frac{1}{J(\mathbf K)}\sum_{i,j}\frac{h_{ij}^*(\mathbf K)\,y_{ij}(\mathbf K)}{\sigma_j^2},\qquad \operatorname{Var}=\frac1{J(\mathbf K)}.$$

Wiener (MAP under Gaussian prior), defined everywhere:

$$\hat{\tilde O}_W=\frac{\sum_{i,j}h_{ij}^*y_{ij}/\sigma_j^2}{J+S_O^{-1}},\qquad \mathrm{MSE}(\mathbf K)=\frac{1}{J(\mathbf K)+S_O(\mathbf K)^{-1}}.$$

Unmeasured set

$$\mathcal U=\Big(\mathbb R^3\setminus\textstyle\bigcup_i C_i\Big)\ \cup\ \Big\{\mathbf K:\ G_j(\mathbf K)=0\ \ \forall j\Big\}.$$

**Proposition to verify (P1, inference).** The universal blind set $\{\mathbf K_\perp=0\}$ (the $K_z$ axis) lies inside the missing cone for illumination scanning, since a cap reaches $\mathbf K_\perp=0$ only at $\mathbf k_\perp=\mathbf k_{i\perp}$, i.e. $\mathbf K=0$. Hence with two non-parallel shears the *extra* loss relative to plain ODT is confined to null-plane intersections $\{\boldsymbol\Delta_j\cdot\mathbf K_\perp\in2\pi\mathbb Z\ \forall j\}\setminus\{\mathbf K_\perp=0\}$. (Addresses Step 2 of [[godt]]'s simulation thread.)

### 2.5 Design subproblem

Choose $\mathbf t\in\Delta^{N-1}=\{t_j\ge0,\ \sum t_j=1\}$ over a dense candidate grid $\mathcal D$ (magnitudes × directions). Let $B$ be the measurable band and $w(\mathbf K)\ge0$ a task weight.

- **(P-A) Bayesian A-optimal**
$$\min_{\mathbf t\in\Delta^{N-1}}\ \int_B\frac{w(\mathbf K)}{\sum_j t_j a_j(\mathbf K)+S_O(\mathbf K)^{-1}}\,d\mathbf K$$
- **(P-D) D-optimal / capacity**
$$\max_{\mathbf t}\ \int_B w\,\log\!\Big(1+S_O\textstyle\sum_j t_j a_j\Big)\,d\mathbf K$$
- **(P-W) worst-case**, on $B'=B\setminus\mathcal N_\varepsilon(\{\mathbf K_\perp=0\})$
$$\max_{\mathbf t}\ \min_{\mathbf K\in B'}\ \sum_j t_j a_j(\mathbf K)\quad(\text{an LP after discretisation})$$

All three are convex in $\mathbf t$. Optional constraints: cardinality $\|\mathbf t\|_0\le k$ (non-convex; use $\ell_1$/reweighting or enumeration for small $k$); hardware constraint of at most two orthogonal directions per polarization cascade ([[multi-shear-acquisition]]).

**Conjectures to test.**
- **C1 (sparsity).** Optimal $\mathbf t$ has few non-zeros, giving a principled answer to "how many shears", to compare with the Zhu group's empirical saturation.
- **C2 (coherence-limited optimum).** With $|\mu(\Delta)|$ decreasing, the optimal magnitude balances $|G|^2$ growth (largest near $\Delta K_\text{max}\approx\pi$) against visibility loss.
- **C3 (incommensurability emerges).** (P-W) optima reproduce the incommensurate-ratio rule without imposing it.

### 2.6 Extensions (relaxing assumptions)

- **E1 Calibration uncertainty (relax A5).** With $\delta\Delta_j\sim\mathcal N(0,s_j^2)$ and the prior replacing $|\tilde O|^2$:
$$J_\text{rob}(\mathbf K)=\sum_j\sum_{i}\frac{|h_{ij}|^2}{\sigma_j^2(t_j)+s_j^2(\mathbf K_\perp\!\cdot\hat{\boldsymbol\Delta}_j)^2\,\big|\tfrac{k_0^2}{2k_z}\big|^2 S_O(\mathbf K)}.$$
Each term has the form $t/(c+dt)$, which is concave, so $J_\text{rob}$ is concave and (P-A) stays convex.
- **E2 Chromatic shear (relax A7).** Spectrum $s(\lambda)$; $\boldsymbol\Delta_j(\lambda)\approx\boldsymbol\Delta_{j0}\lambda/\lambda_0$ for small grating angles *(inference)*; $k_m(\lambda)$ changes cap radius. Treat $(j,\lambda)$ as sub-channels with fixed ratios: $J=\sum_j t_j\int s(\lambda)a_j(\mathbf K;\lambda)\,d\lambda$. Prediction: spectral averaging lifts the nulls. Rigorous version needs the partially coherent 3D transfer function ([[ledwig-robles2023]]); whether $G_j$ stays multiplicative there is open.
- **E3 Correlated channels (relax A6).** Stack $\mathbf y(\mathbf K)$; $J=\mathbf h^H\Sigma^{-1}\mathbf h$, estimator $(\mathbf h^H\Sigma^{-1}\mathbf h)^{-1}\mathbf h^H\Sigma^{-1}\mathbf y$. Compares sequential vs single-shot schemes on equal footing.
- **E4 Sample rotation.** Rotations $R_r$: $G_{jr}(\mathbf K)=e^{i(R_r^\top\boldsymbol\Delta_j)\cdot\mathbf K}-1$ on rotated caps $R_r^\top C_i$. Rotation supplies direction diversity and fills the missing cone at once *(inference)*; gridding becomes non-Cartesian (NUFFT, iterative LS).
- **E5 Phase-only coupling (relax A2).** $2\times2$ information matrix per $\{\mathbf K,-\mathbf K\}$ pair.
- **E6 Spatially varying shear (relax A5 differently).** Breaks voxel-wise diagonality; local-Fourier (windowed) version of $J$. Open.

### 2.7 Validation protocol

- Phantoms: textured cell-like (primary), Shepp–Logan, nano-printed geometry.
- Designs compared at fixed $T N_\text{ph}$: single shear 0.1 μm; single 0.8 μm; orthogonal pair; equal-weight averaging of per-channel reconstructions; BLUE/Wiener fusion; optimal design (P-A, P-W).
- Metrics: MSE spectrum vs theory; RMSE and SSIM of $O$ (or $n$ after integration) on the measurable band; axial elongation on spherical features.
- Code: extend `lab/shearlib.py`; new script to be named here once written.
