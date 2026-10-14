---
type: note
projects: [godt]
created: 2026-09-25
updated: 2026-09-25
status: current
---

# Multi-shear fusion and design: derivation from first principles

The theory behind [[multishear-fusion-design-solution]], derived step by step from
the scalar Helmholtz equation. Its companion [[multishear-fusion-design-problem]]
states the model; this page shows where every factor comes from, which assumption
each step spends, and which "results" are theorems rather than numerics. Everything
here is derivation *(inference, checked against `lab/shearlib.py`)* unless a source is
cited; the only source-backed inputs are [[s42005-025-02416-3]]'s optics
($\lambda=0.5$ μm, $n_m=1.33$, NA 1.33, 60° conical scan, 180 frames; LED 530 nm,
$\Delta\lambda=35$ nm, $\Gamma=6.29$ μm).

Prerequisites, each with its own page: [[optical-diffraction-tomography]],
[[rytov-approximation]], [[common-path-shearing-interferometry]],
[[gradient-phase-imaging]], [[missing-cone-problem]], [[optimal-experimental-design]],
[[partially-coherent-transfer-function]].

**Conventions.** Fourier transforms $\hat f(\mathbf K)=\int f(\mathbf r)e^{-i\mathbf K\cdot\mathbf r}d^3r$
(3D) and $\hat\psi(\boldsymbol\kappa)=\int\psi(\boldsymbol\rho)e^{-i\boldsymbol\kappa\cdot\boldsymbol\rho}d^2\rho$
(2D, detector plane). $\mathbf r=(\boldsymbol\rho,z)$, $z$ the optical axis. $k_0=2\pi/\lambda$,
$k_m=n_mk_0$. Transverse parts carry $\perp$. Assumptions A1–A7 are those of
[[multishear-fusion-design-problem]] §2.2.

---

## 1. Scattering: from Helmholtz to the Fourier diffraction theorem

### 1.1 The wave equation with a scattering potential

A monochromatic scalar field in a medium of index $n(\mathbf r)$ obeys
$\nabla^2U+k_0^2n^2U=0$. Moving the inhomogeneity to the right-hand side,

$$(\nabla^2+k_m^2)\,U(\mathbf r)=O(\mathbf r)\,U(\mathbf r),\qquad O(\mathbf r)=k_0^2\big(n_m^2-n(\mathbf r)^2\big)\approx-2k_0^2n_m\,\delta n(\mathbf r),$$

with $O$ the **scattering potential** (units μm⁻², sign as in `shearlib.cell_phantom`).
It vanishes outside the sample. This is the unknown, and a tomogram is an estimate of
$\hat O(\mathbf K)$ on as much of $\mathbf K$-space as the instrument reaches.

The outgoing Green's function of $(\nabla^2+k_m^2)$ is
$g(\mathbf r)=-e^{ik_m|\mathbf r|}/(4\pi|\mathbf r|)$, so formally
$U=U_i+\int g(\mathbf r-\mathbf r')O(\mathbf r')U(\mathbf r')\,d^3r'$ with incident plane
wave $U_i=e^{i\mathbf k_i\cdot\mathbf r}$, $|\mathbf k_i|=k_m$.

### 1.2 First-order Rytov (A1)

Write $U=U_ie^{\psi}$. Substituting into the wave equation gives the exact Riccati
equation $2\nabla\psi_s\!\cdot\!\nabla\ln U_i+\nabla^2\psi_s+(\nabla\psi_s)^2=O$ for the
scattered complex phase $\psi_s$. Dropping $(\nabla\psi_s)^2$, the only nonlinear
term, is the first-order [[rytov-approximation]]:

$$\psi_s(\mathbf r)\approx\frac{1}{U_i(\mathbf r)}\int g(\mathbf r-\mathbf r')\,O(\mathbf r')\,U_i(\mathbf r')\,d^3r'.$$

It is valid when the *phase gradient* is small, not the phase itself. That is why it
outlives Born for large, gentle cells ([[chen1998]] gives the bounds). The price of
dropping $(\nabla\psi_s)^2$ is priced in [[multiple-scattering-check]].

### 1.3 Weyl expansion and the cap

Downstream of the sample ($z>z'$ for all $\mathbf r'$ in the support), expand the
Green's function in plane waves (Weyl):

$$g(\mathbf r)=\frac{-i}{8\pi^2}\int\frac{e^{i(\mathbf q\cdot\boldsymbol\rho+q_z z)}}{q_z}\,d^2q,\qquad q_z=\sqrt{k_m^2-|\mathbf q|^2}.$$

Evanescent components ($|\mathbf q|>k_m$) die before the detector and are dropped.
Inserting this into §1.2 and doing the $\mathbf r'$ integral turns it into a 3D Fourier
transform of $O$, evaluated at $\mathbf k-\mathbf k_i$ with $\mathbf k=(\mathbf q,q_z)$. The
$\boldsymbol\rho$ transform then collapses the $\mathbf q$ integral onto
$\mathbf q=\boldsymbol\kappa+\mathbf k_{i\perp}$ (the $(2\pi)^2$ cancels the $8\pi^2$
down to 2):

$$\boxed{\ \hat\psi_s(\boldsymbol\kappa;z)=\frac{\hat O(\mathbf K)}{2i\,k_z}\,e^{i(k_z-k_{iz})z},\qquad \mathbf K=\mathbf k-\mathbf k_i,\ \ \mathbf k_\perp=\boldsymbol\kappa+\mathbf k_{i\perp},\ \ k_z=\sqrt{k_m^2-|\mathbf k_\perp|^2}\ }$$

This is the **Fourier diffraction theorem** (Wolf 1969). Three consequences carry the
rest of the page:

1. **One 2D measurement fills one spherical cap.** As $\boldsymbol\kappa$ ranges over
   the detection pupil $|\mathbf k_\perp|\le k_0\mathrm{NA}$, $\mathbf K$ sweeps the cap
   $C_i=\{\mathbf k-\mathbf k_i:\ |\mathbf k|=k_m,\ k_z>0\}$: an Ewald sphere through
   the origin, shifted by $-\mathbf k_i$.
2. **The cap's transverse coordinate is the detector frequency.** Since
   $\mathbf K_\perp=\mathbf k_\perp-\mathbf k_{i\perp}=\boldsymbol\kappa$, a filter applied
   to the detected phase at frequency $\boldsymbol\kappa$ acts on $\hat O$ at
   $\mathbf K_\perp=\boldsymbol\kappa$, *whichever illumination produced the cap*. This
   is the root of the separation theorem (§3).
3. **The Rytov gain is $1/(2k_z)$.** A sample at scattered direction $k_z$ is weighted by
   $1/(2k_z)$. The propagation factor $e^{i(k_z-k_{iz})z}$ is unimodular and carries no
   information; refocusing removes it.

### 1.4 Phase-only data (A2)

The instrument measures the phase $\phi_s=\operatorname{Im}\psi_s$, not the complex
$\psi_s$. Because $\phi_s$ is real,
$\hat\phi_s(\boldsymbol\kappa)=[\hat\psi_s(\boldsymbol\kappa)-\hat\psi_s^*(-\boldsymbol\kappa)]/2i$,
which couples $\hat O(\mathbf K)$ on cap $i$ to $\hat O^*$ at the mirrored point
$\mathbf K'$ of the same cap ($\mathbf K'_\perp=-\mathbf K_\perp$, same $k_z$). The whole
design study takes A2: model $\hat\phi_s$ as if it were $\hat\psi_s$, one complex
observation of $\hat O(\mathbf K)$ per cap point. Relaxing A2 turns every scalar below
into a $2\times2$ block on $\{\mathbf K,\mathbf K'\}$ (E5 of the problem note). For the
design question this matters little, because $|G|^2$ is even in $\mathbf K_\perp$ and
both members of the pair see the same shear factor *(inference)*.

---

## 2. Shearing: the measurement is a Fourier multiplier

### 2.1 Phase difference

A lateral shearing interferometer ([[common-path-shearing-interferometry]]) interferes
the image with a copy displaced by $\boldsymbol\Delta_j$. Its fringe phase is

$$\delta_j\phi(\boldsymbol\rho)=\phi(\boldsymbol\rho+\boldsymbol\Delta_j)-\phi(\boldsymbol\rho),\qquad\phi=\mathbf k_{i\perp}\!\cdot\boldsymbol\rho+\phi_s.$$

The illumination tilt adds only the constant $\mathbf k_{i\perp}\cdot\boldsymbol\Delta_j$,
a piston that is removed. By the shift theorem,

$$\widehat{\delta_j\phi_s}(\boldsymbol\kappa)=G_j(\boldsymbol\kappa)\,\hat\phi_s(\boldsymbol\kappa),\qquad G_j(\boldsymbol\kappa)=e^{i\boldsymbol\Delta_j\cdot\boldsymbol\kappa}-1,$$

$$|G_j|^2=4\sin^2\!\big(\tfrac12\boldsymbol\Delta_j\!\cdot\boldsymbol\kappa\big)\in[0,4].$$

With §1.3, the $(i,j)$ datum at cap point $\mathbf K\in C_i$ is

$$y_{ij}(\mathbf K)=h_{ij}(\mathbf K)\,\hat O(\mathbf K)+n_{ij}(\mathbf K),\qquad h_{ij}(\mathbf K)=\frac{G_j(\mathbf K_\perp)}{2i\,k_z(\mathbf K;i)}.$$

*Units.* [[multishear-fusion-design-problem]] writes $h=G\,k_0^2/(2ik_z)$ with $\hat O$
there being the dimensionless $\widehat{n^2-n_m^2}$. That is the same model with $k_0^2$
moved between $h$ and $\hat O$. Only $\rho$ in §3 changes.

### 2.2 Why shear commutes with everything else

Shear, free-space propagation (refocusing) and the pupil are all shift-invariant, so all
three are Fourier multipliers in $\boldsymbol\kappa$, and multipliers commute. Hence:

- It does not matter optically where the shear is applied, as long as that plane is
  conjugate to the detector. This is exact for the linearised model *(inference)*.
- **It survives multiple scattering in geometry.** Whatever the true scattered phase is,
  the shear multiplies *its* spectrum by $G_j$. The nulls $G_j=0$ are therefore a
  property of the instrument, not of the scattering model. [[multiple-scattering-check]]
  confirmed this numerically to 4e-16.

### 2.3 The dimensionless normalisation

The derivative-style transfer function $H=G/|\boldsymbol\Delta|$ of
[[finite-shear-null-experiment]] tends to $i\,\hat{\boldsymbol\Delta}\cdot\boldsymbol\kappa$
as $\Delta\to0$, which is the gradient. Its $1/\Delta$ is a convention. What matters
for information is the noise on the quantity actually recorded, which is the phase
difference $\delta_j\phi$, not $\delta_j\phi/\Delta$ (§3.1). So $|G|^2$, peak 4, is the
honest weight.

---

## 3. Noise, Fisher information and the separation theorem

### 3.1 Noise model (A4, A6)

In a two-beam interferometer with $N$ detected photons per pixel and fringe
visibility $V$, the Cramér–Rao bound on the fringe phase is
$\operatorname{Var}\phi\ge 2/(NV^2)$ for four-step phase shifting, and of the same order
for off-axis demodulation. A shearing interferometer's two beams are the *same* field
at two points, so its visibility is the complex degree of spatial coherence across the
shear: $V=|\mu(\boldsymbol\Delta_j)|$. With a time share $t_j$ of a total exposure
$T N_\text{ph}$,

$$\operatorname{Var}\big[\delta_j\phi(\boldsymbol\rho)\big]=\frac{\sigma_0^2}{t_j\,T N_\text{ph}\,|\mu(\boldsymbol\Delta_j)|^2}\equiv\sigma_j^2,$$

white across pixels. It is *independent of $|\boldsymbol\Delta_j|$ at fixed visibility*;
that is assumption A4, and it is what lets $|G|^2\le4$ exceed the direct-phase weight
1. White Gaussian noise stays white Gaussian under a unitary DFT, so each Fourier
sample carries independent noise of variance $\sigma_j^2$, and any constant from the
DFT normalisation is absorbed in $\sigma_0$. A6 (sequential channels) makes different
$j$ independent. Hermitian symmetry of real noise pairs $\boldsymbol\kappa$ with
$-\boldsymbol\kappa$; a factor 2 is absorbed in $\sigma_0$ *(inference)*. Phase wrapping
is excluded (A4). It is an independent argument for small $\Delta$, analysed in
[[multiple-scattering-check]] §6.

### 3.2 Per-voxel Fisher information

A3 (exact gridding) puts every cap sample on a voxel of the $\mathbf K$ grid. The
observations of $\hat O(\mathbf K)$ are then $\{y_{ij}(\mathbf K)\}$ for
$i\in\mathcal I(\mathbf K)$, the illuminations whose caps pass through that voxel, and
every channel $j$. They are independent complex Gaussians with means
$h_{ij}\hat O(\mathbf K)$. Different voxels are decoupled: the forward operator is
diagonal in $\mathbf K$. The Fisher information for $\hat O(\mathbf K)$ is therefore a
scalar:

$$J(\mathbf K)=\sum_j\sum_{i\in\mathcal I(\mathbf K)}\frac{|h_{ij}(\mathbf K)|^2}{\sigma_j^2}
=\frac{TN_\text{ph}}{\sigma_0^2}\sum_j t_j|\mu(\boldsymbol\Delta_j)|^2|G_j(\mathbf K_\perp)|^2\sum_{i\in\mathcal I(\mathbf K)}\frac{1}{4k_z(\mathbf K;i)^2}.$$

### 3.3 The separation theorem (T1)

The double sum factorises because the $j$-dependence ($G_j$, $\mu$, $t_j$) and the
$i$-dependence ($k_z$) are carried by different factors:

$$\boxed{\ J(\mathbf K)=\rho\;S_{\mathbf t}(\mathbf K_\perp)\;M(\mathbf K),\quad
S_{\mathbf t}=\sum_jt_j|\mu(\boldsymbol\Delta_j)|^2|G_j(\mathbf K_\perp)|^2,\quad
M=\sum_{i\in\mathcal I(\mathbf K)}\Big(\frac{k_m}{k_z(\mathbf K;i)}\Big)^2,\quad
\rho=\frac{TN_\text{ph}}{4\sigma_0^2k_m^2}\ }$$

**What the factorisation actually needs.** Only two things. (a) The shear factor must
not depend on the illumination index. It does not, by §1.3 point 2: $G_j$ is a
function of $\mathbf K_\perp=\boldsymbol\kappa$ whichever cap carried the sample. (b)
Every channel must be acquired at every illumination with the same share $t_j$, a full
factorial design. Conical illumination is **not** needed for the factorisation. It is
needed for §4's *closed form* of $M$, and for $k_z$ to depend on $K_z$ alone.
[[multishear-fusion-design-solution]] §1 lists conical illumination as one of the "two
facts" behind exactness; strictly it is only a fact behind the closed form. If (b) is
violated, for example by different shear sets at different illuminations, $J$ becomes
$\sum_{ij}t_{ij}(\cdot)$ and the design problem is genuinely 3D.

**Why it makes the design cheap.** The design $\mathbf t$ enters only through the 2D
function $S_{\mathbf t}(\mathbf K_\perp)$, which is constant along every $K_z$ column. So
any objective $\int_B w\,F(J,\mathbf K)\,d^3K$ can be precomputed per column and
optimised over a 2D field.

---

## 4. The geometry factor $M(\mathbf K)$ in closed form

### 4.1 Which caps pass through a voxel

Conical illumination: $\mathbf k_i=k_m(\sin\theta\cos\phi,\ \sin\theta\sin\phi,\ \cos\theta)$,
with $N_\text{fr}$ azimuths $\phi$ uniform on $[0,2\pi)$. The point $\mathbf K$ lies on
cap $\phi$ iff $|\mathbf K+\mathbf k_i|=k_m$, i.e. $|\mathbf K|^2+2\mathbf K\cdot\mathbf k_i=0$.
Write $\mathbf K_\perp=K_\perp(\cos\psi,\sin\psi)$:

$$|\mathbf K|^2+2k_m\big(K_\perp\sin\theta\cos(\phi-\psi)+K_z\cos\theta\big)=0
\iff\cos(\phi-\psi)=c(\mathbf K)\equiv-\frac{|\mathbf K|^2+2k_mK_z\cos\theta}{2k_mK_\perp\sin\theta}.$$

For $|c|<1$ there are two roots $\phi_\pm=\psi\pm\arccos c$, and for $|c|>1$ none. Every
cap has the same $k_{iz}=k_m\cos\theta$, so the scattered
$k_z=K_z+k_m\cos\theta$ is the same for both roots and depends on $K_z$ alone.

### 4.2 Counting digital samples

Digital caps assign one $K_z$ per detector pixel. A voxel of edge $\delta k$ at
$(\mathbf K_\perp,K_z)$ collects the azimuths whose cap height at that $\mathbf K_\perp$ falls
within $\pm\delta k/2$ of $K_z$. Along the $K_\perp$ column, cap $\phi$ has height
$K_z(\phi)=\sqrt{k_m^2-|\mathbf K_\perp+\mathbf k_{i\perp}(\phi)|^2}-k_m\cos\theta$, with

$$\Big|\frac{dK_z}{d\phi}\Big|=\frac{|(\mathbf K_\perp+\mathbf k_{i\perp})\cdot\partial_\phi\mathbf k_{i\perp}|}{k_z}
=\frac{k_m\sin\theta\,K_\perp\,|\sin(\phi-\psi)|}{k_z}=\frac{k_m\sin\theta\,K_\perp\sqrt{1-c^2}}{k_z}.$$

With azimuth density $N_\text{fr}/2\pi$ and two roots, the count is
$2\cdot\frac{N_\text{fr}}{2\pi}\cdot\delta k/|dK_z/d\phi|$. Each sample carries weight
$(k_m/k_z)^2$, which gives

$$\boxed{\ M(\mathbf K)=\frac{N_\text{fr}\,\delta k}{\pi}\cdot\frac{k_m}{k_z\,\sin\theta\,K_\perp\sqrt{1-c(\mathbf K)^2}}\quad\text{on }\{|c|<1,\ k_z\in[k_{z,\text{lo}},k_m]\},\ \text{else }0\ }$$

where $k_{z,\text{lo}}=\sqrt{k_m^2-(k_0\mathrm{NA})^2}$ enforces the detection pupil. This
matches `shearlib.cap_density_closed`, which the solution note checks against direct
binning to 0.5 % in total. The $1/\sqrt{1-c^2}$ divergence at $|c|=1$ is the **caustic**:
the envelope of the cap family, where two roots merge and $K_z(\phi)$ is stationary. It
is integrable, $\int d\phi$ being finite.

### 4.3 The missing cone, derived

Near the origin $|\mathbf K|^2\ll k_m|\mathbf K|$, so $c\approx-K_z\cot\theta/K_\perp$.
Sampling requires $|c|<1$, i.e. $|K_z|<K_\perp\tan\theta$. The unsampled region is a
double cone around the $K_z$ axis with half-angle $\arctan(\cot\theta)=90^\circ-\theta$,
which is **30° for θ = 60°**. This is the [[missing-cone-problem]] as a consequence of
§4.1 rather than a picture. At finite $|\mathbf K|$ the $|\mathbf K|^2$ term curves the
boundary, and the cone narrows away from the origin.

**P1 ($K_z$ axis unsampled).** $\mathbf K_\perp=0$ on cap $i$ means
$\mathbf k_\perp=\mathbf k_{i\perp}$. With $|\mathbf k|=|\mathbf k_i|$ and both $k_z>0$, that
forces $\mathbf k=\mathbf k_i$, so $\mathbf K=0$. Hence the $K_z$ axis carries no samples
apart from the origin, for *any* illumination pattern, not only conical ones.

### 4.4 Band limit

The largest transverse reach is $K_{\max}=k_0\mathrm{NA}+k_m\sin\theta$ (detection plus
illumination). For GODT's simulation optics that is $16.71+14.47=31.19$ rad/μm, and
$2\pi/K_{\max}=0.2014$ μm. This takes the 60° as measured in the medium, which is
*(inference)*: [[s42005-025-02416-3]] does not say.

---

## 5. Nulls, blind sets and incommensurability

From §2.1, $G_j=0\iff\boldsymbol\Delta_j\cdot\mathbf K_\perp\in2\pi\mathbb Z$. $G_j$ does
not depend on $K_z$, so in 3D each null is a **plane** containing the $K_z$ direction,
perpendicular to $\boldsymbol\Delta_j$, at spacing $2\pi/|\boldsymbol\Delta_j|$:

- $m=0$, the **blind plane** $\mathbf K_\perp\cdot\boldsymbol\Delta_j=0$. It exists for every
  magnitude and meets the sampled support in a set of positive area. No choice of
  $|\boldsymbol\Delta_j|$ removes it.
- $m\neq0$, the **null comb**, known in 2D as the Servín null structure ([[servin2007]]).
  It enters the band iff $2\pi/|\boldsymbol\Delta_j|<K_{\max}$.
- Near the blind plane $|G_j|^2\approx(\boldsymbol\Delta_j\cdot\mathbf K_\perp)^2$. A small
  shear has a *wide* low-information slab $|K_\parallel|\lesssim1/|\boldsymbol\Delta_j|$.

**Common nulls of a set.** A voxel is lost for a design iff $G_j=0$ for all $j$ with
$t_j|\mu_j|>0$.

- *Parallel pair* $d_1,d_2$. Nulls coincide at $K_\parallel=2\pi m/d_1=2\pi n/d_2$, which
  needs $d_1/d_2=m/n$. For $d_1/d_2=p/q$ in lowest terms, the first common null is at
  $2\pi p/d_1$. It lies in band iff $p<d_1K_{\max}/2\pi$, so **small-integer ratios are
  worst**. Irrational ratios have no common null except the blind plane, but they have
  near-coincidences at the continued-fraction convergents. That is the incommensurate
  rule (C3), derived.
- *Non-parallel pair*. The blind planes meet in the $K_z$ axis, which P1 says is already
  empty. The combs meet in vertical **lines** through the reciprocal-lattice points
  $\mathbf K_\perp=2\pi B^{-\top}\mathbf m$ with $B=[\boldsymbol\Delta_1\ \boldsymbol\Delta_2]$.
  For an orthogonal pair of magnitude $d$ these sit at $(2\pi/d)(m_1,m_2)$, i.e. every
  7.85 rad/μm for 0.8 μm, which puts many inside a 31 rad/μm band. These are the
  "doubly-null voxels" of the solution note's §5.
- *Parallel and orthogonal pairs combined*. $\{d_1,d_2\}\,x+y$ with $d_1/d_2$
  incommensurate has no common null in the band apart from the empty axis, so $J>0$ on
  the whole support. That is the condition for BLUE to exist (§6).

**The single-shear design rule, stated correctly.** For one channel, information across
the band grows with $|\boldsymbol\Delta|$: the blind slab narrows and $|G|^2$ rises
towards its peak at $\boldsymbol\Delta\cdot\mathbf K_\perp=\pi$. It keeps growing until the
first comb plane $2\pi/|\boldsymbol\Delta|$ enters the band. The optimum therefore sits
*at* the boundary, $|\boldsymbol\Delta|\approx2\pi/K_{\max}$, not anywhere below it. The
solution note's headline originally wrote this as "$\Delta\lesssim2\pi/K_{\max}$". That wording is
too loose: its own §3 shows 0.1 μm, which satisfies the inequality, doing *worse* than
0.33 μm, which violates it. The rule is "as large as possible without admitting the
comb", i.e. $\Delta\approx2\pi/K_{\max}$.

---

## 6. Estimation and fusion

### 6.1 BLUE and Wiener per voxel

Voxels decouple (§3.2), so estimation is a scalar problem per $\mathbf K$. Define the
back-projected data $b(\mathbf K)=\sum_{ij}h_{ij}^*y_{ij}/\sigma_j^2$. Then:

- **BLUE / maximum likelihood** (where $J>0$): $\hat O_\text{BLUE}=b/J$, unbiased, with
  variance $1/J$. It is undefined on $J=0$ and has unbounded variance as $J\to0$, which
  is why it diverges for a single shear.
- **Wiener / MAP**, under the Gaussian prior
  $\hat O(\mathbf K)\sim\mathcal{CN}(0,S_O(\mathbf K))$: the posterior is Gaussian with
  precision $J+S_O^{-1}$, so

$$\hat O_W=\frac{b}{J+S_O^{-1}},\qquad\mathrm{MSE}(\mathbf K)=\frac{1}{J(\mathbf K)+S_O(\mathbf K)^{-1}},\qquad\frac{\mathrm{MSE}}{S_O}=\frac{1}{1+S_O J}.$$

Both expressions are exact under A1–A7, not asymptotic, because the model is linear
Gaussian. That is why the solution note's §5 finds empirical/predicted ratios of 1.00
on a Gaussian-random-field phantom.

### 6.2 Fusion *is* precision-weighted averaging

Per channel, $b=\sum_jb_j$ and $J=\sum_jJ_j$, with $J_j=\rho\,t_j|\mu_j|^2|G_j|^2M$. So

$$\hat O_W^{\text{fused}}=\frac{\sum_jb_j}{\sum_jJ_j+S_O^{-1}}=\frac{\sum_jJ_j\,\hat O_{j,\text{BLUE}}}{\sum_jJ_j+S_O^{-1}}.$$

The fused Wiener estimate is exactly the **per-frequency inverse-variance average** of
the per-channel estimates, shrunk by the prior. Two consequences:

- An *equal-weight* average of per-channel Wiener reconstructions,
  $\frac1N\sum_jb_j/(J_j+S_O^{-1})$, is biased toward zero by the factor
  $\frac1N\sum_jS_OJ_j/(1+S_OJ_j)$. On a null of channel $j$ that channel contributes 0,
  and the average is diluted by $1/N$ at every such frequency. Averaging in real space
  equals averaging in $\mathbf K$-space, since the FT is linear, so the "naive average" is
  this estimator. Numerically, a precision-weighted average (weights $J_j+S_O^{-1}$)
  gives 0.0079 against fused 0.0078 and naive 0.1607
  (`lab/exp05b_weighted_average.py`). The solution note's 21× (and [[multiple-scattering-check]]'s 9.6×
  model-A ladder) measure **the value of per-frequency weights**, not of a new estimator.
- Anyone who keeps the per-channel reconstructions in $\mathbf K$-space and weights each
  voxel by $J_j$ reproduces fusion exactly. What is specific to this
  study is the *design* and the *3D quantification*, and showing the weighting matters
  this much. The estimator is the standard information-form combination used in
  multi-frame deconvolution and SIM Wiener reconstruction. Legarda-Saenz 2011's
  Fourier-domain multi-direction fusion is flagged on [[zhu2025-sci-rep]] as possible
  direct prior art for the 2D case *(inference; not yet checked against the paper)*.

### 6.3 What no design can recover

The MSE formula shows that a voxel with $J=0$ returns the prior mean (zero) with error
$S_O$. That holds for the missing cone ($M=0$) and for common nulls ($S_{\mathbf t}=0$).
Low $|\mathbf K|$ is the dangerous region for gradient data. There, $|G_j|^2\approx(\boldsymbol\Delta_j\cdot\mathbf K_\perp)^2$
is small for *small* shears, while $S_O$ of a cell is largest: most of a cell's spectral
energy sits at the lowest frequencies. For large shears (0.8 μm at $K=2$ rad/μm,
$|G|^2\approx2$) the low-$K$ shell is *not* intrinsically lost. The one voxel that no
shear design can reach is **DC**, $\mathbf K=0$, where $G_j(0)=0$ for every $j$: the mean
index offset, which gradient data cannot see and an absolute reference must supply. On
exp05's cell phantom DC holds 19.4 % of the power. Every non-degenerate pair then
recovers the rest of the $|K|<2$ shell (shell error 0.274, exactly DC's share), while
a single shear loses about half of it (0.55). Source: `lab/exp05b_weighted_average.py`,
discussed in [[multishear-fusion-design-solution]] §5.

---

## 7. Design as optimal experimental design

### 7.1 Convexity

Let $a_j(\mathbf K)=\rho|\mu_j|^2|G_j|^2M$ so that $J=\sum_jt_ja_j$ is linear in $\mathbf t$ on
the simplex $\{t_j\ge0,\sum t_j=1\}$. Then:

- (P-A) $f_A=\int_Bw/(J+S_O^{-1})$ is convex, since $x\mapsto1/(x+c)$ is convex and
  decreasing on $x\ge0$ and composition with an affine map preserves convexity.
- (P-D) $f_D=\int_Bw\log(1+S_OJ)$ is concave (log of affine).
- (P-W) $\max_{\mathbf t}\min_{\mathbf K\in B'}J$ is the maximum of a concave piecewise
  linear function, which is an LP.

See [[optimal-experimental-design]].

### 7.2 Optimality certificate (equivalence theorem)

For a differentiable convex $f$ on the simplex with gradient $g_j=\partial f/\partial t_j$,
convexity gives, for the optimum $\mathbf t^*$,
$f(\mathbf t^*)\ge f(\mathbf t)+g\cdot(\mathbf t^*-\mathbf t)\ge f(\mathbf t)+\min_jg_j-\sum_jt_jg_j$.
Hence

$$0\le f(\mathbf t)-f^*\le\sum_jt_jg_j-\min_jg_j,$$

which is zero iff no vertex is a descent direction. For (P-A),
$g_j=-\int_Bw\,a_j/(J+S_O^{-1})^2$. The condition becomes the Kiefer–Wolfowitz form
$\max_j\varphi_j\le\sum_jt_j\varphi_j$ with $\varphi_j=\int w a_j/(J+S_O^{-1})^2$. The
solution note's "gap $\le3\times10^{-4}$" is this bound divided by $|f|$, and it matches
`DesignProblem.certificate`.

### 7.3 Why the optimum need not be sparse (C1)

Carathéodory-type sparsity (an optimal design with at most $p(p+1)/2+1$ support points)
holds when the information is a $p\times p$ matrix. Here the "information" is the whole
function $J(\mathbf K)$ on $B$, whose dimension is the number of voxels. The theorem
therefore gives no useful bound, and with strongly overlapping $|G_j|^2$ the optimum can
spread over many channels. Refuting C1 is consistent with the theory, not a surprise.

### 7.4 Why coherence sparsifies (C2)

Visibility multiplies each channel's efficiency: $a_j\propto|\mu(\boldsymbol\Delta_j)|^2$.
For a uniform circular incoherent source, van Cittert–Zernike gives
$\mu(\Delta)=2J_1(x)/x$ with $x=\pi\Delta/l_c$ (see
[[partially-coherent-transfer-function]]; the $l_c$ is a **model parameter**). Large
shears are penalised, and past the first zero ($\Delta=1.22\,l_c$) they are useless. The
effective candidate set shrinks, and the optimum concentrates on the few magnitudes that
balance the growth of $|G|^2$ against the loss of $|\mu|^2$. The balance point moves
toward small $\Delta$ as $l_c$ falls.

---

## 8. Robustness

### 8.1 Calibration error (E1)

Let the true shear be $\boldsymbol\Delta_j+\boldsymbol\delta_j$ with a small unknown
$\boldsymbol\delta_j$. Then

$$G_j^{\text{true}}=G_j+i(\boldsymbol\delta_j\cdot\mathbf K_\perp)e^{i\boldsymbol\Delta_j\cdot\mathbf K_\perp}+O(\delta^2).$$

The reconstruction assumes $G_j$, so the datum carries an extra term
$i(\boldsymbol\delta_j\cdot\mathbf K_\perp)e^{i\boldsymbol\Delta_j\cdot\mathbf K_\perp}\hat O/(2ik_z)$.
Averaging over $\boldsymbol\delta_j$ with component rms $s$ along
$\hat{\boldsymbol\Delta}_j$, and replacing $|\hat O|^2$ by its prior $S_O$, gives an
error term of variance $s^2(\mathbf K_\perp\cdot\hat{\boldsymbol\Delta}_j)^2S_O/(4k_z^2)$
added to $\sigma_j^2/t_j$:

$$J_{\text{rob}}=\sum_j\frac{t_ja_j}{1+t_j\,a_j\,s^2(\mathbf K_\perp\!\cdot\hat{\boldsymbol\Delta}_j)^2S_O/|G_j|^2}\ \xrightarrow{t_j\to\infty}\ \frac{|G_j|^2}{s^2(\mathbf K_\perp\!\cdot\hat{\boldsymbol\Delta}_j)^2}\,S_O^{-1}.$$

Each term has the form $t/(c+dt)$, which is concave, so (P-A) stays convex. The
**ceiling** $|G_j|^2/(sK_\parallel)^2$ on the per-voxel SNR does not depend on photons:
more exposure cannot beat a wrong model. The solution note's $\lvert G\rvert^2/(sK_\perp)^2$
is this ceiling evaluated along $\hat{\boldsymbol\Delta}$.

Two approximations are hidden here.

1. $\boldsymbol\delta_j$ is **one** number per channel, shared by every voxel and every
   illumination. Treating it as independent noise per voxel overstates how much fusion
   across channels averages it down. With $N$ channels carrying independent
   $\boldsymbol\delta_j$, the fused error averages over $N$ draws, not over $N\times$
   (voxels) *(inference)*.
2. Linearisation requires $sK_\perp\ll1$. At $K_\perp=30$ rad/μm and $s=50$ nm,
   $sK_\perp=1.5$, so the last row of the solution note's E1 table is outside the
   derivation's validity.

### 8.2 Chromatic shear (E2)

A polarization-grating pair of period $\Gamma$ and separation $L$ shears by
$\Delta=L\tan\theta_\lambda$ with $\sin\theta_\lambda=\lambda/\Gamma$. For small angles
$\Delta(\lambda)\approx\Delta_0\lambda/\lambda_0$, so the fractional shear spread equals the
fractional bandwidth $\varepsilon=(\lambda-\lambda_0)/\lambda_0$, of rms $\sigma_\varepsilon$. For a
Gaussian with 35 nm FWHM at 530 nm, $\sigma_\varepsilon=0.0280$. [[s42005-025-02416-3]] gives
only "Δλ = 35 nm"; the Gaussian FWHM reading is a modelling choice.

Wavelengths add in intensity. For a weak object the demodulated fringe term is linear
in the object, so the effective multiplier is the spectral average
$G_\text{eff}(\mathbf K_\perp)=\langle e^{i\boldsymbol\Delta(\lambda)\cdot\mathbf K_\perp}\rangle_\lambda-1$.
At the $m$-th null $\Delta_0K=2\pi m$ this is

$$G_\text{eff}=\langle e^{2\pi im\varepsilon}\rangle-1=e^{-2\pi^2m^2\sigma_\varepsilon^2}-1\approx-2\pi^2m^2\sigma_\varepsilon^2,\qquad\frac{|G_\text{eff}|^2}{4}\approx\pi^4m^4\sigma_\varepsilon^4.$$

That gives $6.0\times10^{-5}$, $9.6\times10^{-4}$ and $4.9\times10^{-3}$ for $m=1,2,3$
(the exact Gaussian values are $5.9\times10^{-5}$, $9.1\times10^{-4}$ and $4.3\times10^{-3}$).
*Why higher nulls fill:* the $m$-th null's position spreads by $2\pi m\sigma_\varepsilon/\Delta_0$,
growing with $m$, while its width in $K$ does not. The residue therefore grows as
$m^4$.

This model holds three things fixed that the real instrument does not, all
*(inference)*:

- the object phase $\phi_s\propto k_0$, which scales as $1/\lambda$;
- the Ewald radius $k_m(\lambda)$, so each wavelength samples a slightly different cap;
- the Rytov weight $1/k_z$.

For $\sigma_\varepsilon\approx3\,\%$ these are few-percent corrections to $|h|$ but move the caps by
$\sim\!\sigma_\varepsilon K$. The spectral average is then over *different* $\mathbf K$, not
only different $G$, which is the [[ledwig-robles2023]] question left open.

---

## 9. Status of each result

| result | status |
|---|---|
| Fourier diffraction theorem, $h=G/(2ik_z)$ | derived (A1, A2) |
| Separation $J=\rho S_{\mathbf t}M$ | **theorem** under A1–A6 plus a full-factorial acquisition; conical not required |
| Closed form $M(\mathbf K)$, caustic | **theorem** under conical illumination and A3; checked numerically |
| Missing-cone half-angle $90^\circ-\theta$ | theorem (small-$\lvert\mathbf K\rvert$ limit) |
| P1: $K_z$ axis empty | theorem, any illumination pattern |
| Null planes, blind plane, lattice lines, rational-ratio coincidences | theorem |
| Single-shear optimum $\Delta\approx2\pi/K_{\max}$ | argued here (monotone gain until the comb enters); located numerically |
| BLUE/Wiener MSE $1/(J+S_O^{-1})$ | exact under A1–A7 |
| Fusion = per-frequency precision-weighted average | **identity** |
| Convexity of P-A/P-D/P-W, equivalence gap | theorem |
| C1 false, C2, C3 | numerical; theory explains *why* (§7.3–7.4, §5) |
| E1 ceiling $\lvert G\rvert^2/(sK_\parallel)^2$ | first order in $sK$; treats a per-channel constant as per-voxel noise |
| E2 residue $\pi^4m^4\sigma_\varepsilon^4$ | second order; fixed-cap, fixed-object-phase model |

Numerics: `lab/shearlib.py` (functions `cap_density_closed`, `DesignProblem`,
`chromatic_G`, `visibility`) and the scripts listed in
[[multishear-fusion-design-solution]] § Reproduction.
