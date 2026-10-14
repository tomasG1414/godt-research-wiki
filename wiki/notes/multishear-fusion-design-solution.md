---
type: note
projects: [godt]
created: 2026-09-21
updated: 2026-10-05
status: current
---

# Multi-shear fusion and design: solved

Answers the problem posed in [[multishear-fusion-design-problem]]. Everything here
is **computed**, not sourced: `lab/exp03_ewald_fisher.py`, `lab/exp04_shear_design.py`,
`lab/exp05_fusion_validation.py`, all on shared math in `lab/shearlib.py` (extended
with the 3D Ewald and design machinery, selftest updated). The only inputs taken from
a source are [[s42005-025-02416-3]]'s Eq. 20 shear factor and its optical parameters;
every consequence below is derived here. Assumptions A1–A7 of the problem note hold
unless a section says otherwise.

> **Follow-up, 2026-09-22.** [[multiple-scattering-check]] re-ran this study's
> design comparison against a [[brenner1993]] WPM forward model. The geometry
> below survives exactly; the margins do not. Read §7 and §8 alongside it.
>
> **Audit, 2026-09-25.** The theory is now derived from the Helmholtz equation in
> [[multishear-fusion-derivation]]. Every number below was reproduced by rerunning
> exp03–05. The corrections are to framing, not arithmetic:
> - the 21× is a margin over *equal-weight* averaging, and precision-weighted averaging equals fusion (§5);
> - the design rule is $\Delta\approx2\pi/K_{\max}$, not $\lesssim$;
> - the $|K|<2$ shell is a mask, not a physical loss;
> - E1 runs on a different grid from §4;
> - C3's ratio is $0.8/d_2$;
> - "0.1 px = 10 nm" is an inference.
>
> The affected passages are marked **[audit]**.

## Headline

1. **The problem separates exactly.** $J(\mathbf K)=\rho\,S_{\mathbf t}(\mathbf K_\perp)\,M(\mathbf K)$ — the shear design enters only through a 2D function, the tomography only through a cap density with a closed form. Design and geometry never mix. This is what makes the whole study cheap.
2. **P1 is true, and it was the wrong thing to worry about.** The direction-*independent* blind set $\{\mathbf K_\perp=0\}$ is unsampled by transmission ODT anyway. But the direction-*dependent* blind plane $\mathbf K_\perp\!\cdot\!\boldsymbol\Delta=0$ of a single shear **is** measured, and is 0.50 % of the support — a real, avoidable loss.
3. **The single-shear magnitude should be $\Delta\approx2\pi/K_{\max}$: as large as possible without admitting the null comb.** The optimal single channel is 0.20 μm, and $2\pi/K_{\max}=0.201$ μm. Smaller is *not* safe: the 0.1 μm simulation shear loses more than 0.33 μm (§3), because a small shear widens the blind slab. These numbers use the *simulation* optics; the experiment's illumination is 58° in air per [[2024-06-14-598984-polarization-hilbert-holotomography]], and the rule has not been re-run on that geometry.
4. **Fusing in $\tilde O$-space beats *equal-weight* averaging of reconstructions by 21×** on identical data at equal photon budget. **[audit]** A per-frequency precision-weighted average of the same per-channel reconstructions reproduces fusion (0.0079 against 0.0078). The fused Wiener estimator *is* that weighted average ([[multishear-fusion-derivation]] §6.2). The 21× therefore measures the value of per-frequency weighting, not a new estimator. On [[multiple-scattering-check]]'s phantom the margin is 9.6× (linear) → 1.8× (realistic gridding) → 1.1× (multiple scattering).
5. **Optimal designs are not sparse** — C1 is false. What makes them sparse is *coherence*: as visibility falls, the optimum collapses from 24 channels to 6.
6. **Calibration can cap the high frequencies.** Assuming 10 nm rms random shear error, the calibration-limited SNR at the band edge is 16.5 dB against a photon-limited 37.1 dB on E1's grid. The 10 nm is what "0.1 px" becomes on a 0.1 μm simulation pixel *(inference)*.

---

## 1. The separation theorem (T1)

For illumination azimuth $\phi$ and shear channel $j$, assembling §2.4 of the problem
note over the conical scan gives, at every voxel of the sampled support,

$$J(\mathbf K)=\underbrace{\frac{T N_\text{ph}k_0^4}{4\sigma_0^2 k_m^2}}_{\rho}\ \underbrace{\sum_j t_j\,|\mu(\boldsymbol\Delta_j)|^2\,\big|e^{i\boldsymbol\Delta_j\cdot\mathbf K_\perp}-1\big|^2}_{S_{\mathbf t}(\mathbf K_\perp)\ \text{— the design}}\ \underbrace{\sum_{i\in\mathcal I(\mathbf K)}\Big(\frac{k_m}{k_z(\mathbf K;i)}\Big)^{2}}_{M(\mathbf K)\ \text{— the geometry}}$$

The factorisation is exact because the shear factor does not depend on which
illumination delivered the sample. It depends on $\mathbf K_\perp$ only, the detector
frequency, and is constant along the $K_z$ column. Every channel is also acquired at
every illumination. **[audit]** Conical illumination is *not* needed for the separation
([[multishear-fusion-derivation]] §3.3). It is needed for the closed form of $M$ below,
where it makes the scattered $k_z=K_z+k_m\cos\theta$ depend on $K_z$ alone. The design
problem is therefore 2D and the tomography is a fixed weight. That is why a dense
shear-candidate grid over a $256^2\times96$ volume takes minutes, not days.

Note the natural normalisation: the **dimensionless** shear power
$|G_j|^2=4\sin^2(\boldsymbol\Delta_j\!\cdot\!\mathbf K_\perp/2)$, peak 4, not
$|H|^2=|G|^2/\Delta^2$. Under A4 the measured phase difference is one measurement
whose noise does not grow with $|\Delta|$, so $|G|^2$ is the right weight. The $1/\Delta$
in [[finite-shear-null-experiment]]'s $H$ leaves the null positions unchanged but does
change how designs *rank*: it favours small shears, which is why that page's 1D ranking
of second shears is superseded by §4 here.

### $M(\mathbf K)$ in closed form

Writing $\psi$ for the azimuth of $\mathbf K_\perp$, the caps through $\mathbf K$ are
the roots of $\cos(\phi-\psi)=c(\mathbf K)$ with

$$c(\mathbf K)=-\frac{|\mathbf K|^2+2k_mK_z\cos\theta}{2k_mK_\perp\sin\theta},\qquad
M(\mathbf K)=\frac{N_\text{fr}\,\delta k}{\pi}\cdot\frac{k_m}{k_z\sin\theta\,K_\perp\sqrt{1-c^2}}$$

on $\{|c|<1\}$ and zero elsewhere. Two samples per voxel for $|c|<1$, none outside;
the integrable $1/\sqrt{\cdot}$ divergence on $|c|=1$ is the caustic where the cap
family folds. Checked against direct binning of digital caps: total within 0.5 %,
voxel-wise median deviation 0.6 % and 90th percentile 3.2 % on a $128^2\times48$ grid
with 11 520 azimuths (the residual is the $K_z$ rounding of the digital caps, and it
shrinks with azimuth count). **[audit]** On the $256^2\times96$ grid that §3–4 actually
use, exp03 prints a total ratio of 0.979, a median deviation of 7.6 % and a 90th
percentile of 26 %. The voxel-wise agreement is much looser there, and it is the
closed form, not the binning, that feeds the design.

**The missing cone falls out of $|c|=1$.** Near the origin the boundary is
$K_z=\pm K_\perp\tan\theta$: half-angle $90^\circ-\theta=30^\circ$ for GODT's 60° cone,
narrowing to 19° at $K_z=5.9$ rad/μm as the caps curve. That is the
[[missing-cone-problem]] derived rather than asserted.

### Angular sampling is a separate budget

180 projections reach only **76.6 %** of the support the geometry allows on a 25.6 μm
field (88.4 % on 12.8 μm), under nearest-neighbour direct inversion. [[multiple-scattering-check]]
§1 reports 92.1 % for a 12.8 μm field. The difference is unreconciled; its axial
extent and pupil floor may differ *(inference)*. This is
orthogonal to the shear design — every channel shares the same caps — but it means
a 180-projection acquisition is angle-limited before shear enters, and it is the argument
for [[pan1983]]'s bilinear interpolation over nearest-neighbour.

## 2. P1, and the correction it needs (T2)

**The $K_z$ axis carries exactly zero samples except at $\mathbf K=0$.** A cap reaches
$\mathbf K_\perp=0$ only where $\mathbf k_\perp=\mathbf k_{i\perp}$, i.e. at
$\mathbf K=0$. Confirmed exactly on the digital caps. So the blind set common to
*every* in-plane shear set is unmeasured by transmission ODT regardless — the shear
blind spot and the missing cone coincide there, and that part of the shear penalty
is already paid. P1 holds.

**But P1 was the wrong worry.** The blind set of a *single* shear is not the axis, it
is the whole plane $\mathbf K_\perp\!\cdot\!\boldsymbol\Delta=0$, and that plane is
measured:

| set | share of the sampled support |
|---|---|
| plane $K_x=0$ exactly | 0.50 % |
| slab $\lvert K_x\rvert<1$ rad/μm | 4.52 % |
| within 1 rad/μm of the $K_z$ axis | 0.019 % |

The second shear direction is what buys the plane back; its residual blind set (the
axis) is free. So direction diversity has a *structural* justification in 3D, not
only the error-averaging one the Zhu group measured — but it is worth 0.5 % of
the support, which is the right order for the modest gains they report.

## 3. What the nulls actually cost in 3D

Fraction of the sampled support with relative Wiener MSE > 0.25, at $\rho=100$ per
frame-sample and a flat prior, every design at **equal photon budget**:

| design | lost volume | mean rel. MSE |
|---|---|---|
| ODT benchmark ($\lvert G\rvert^2\!=\!1$) | 0.00 % | 0.0016 |
| 0.1 μm, $x$ (GODT simulation; paper shears in $y$, equivalent by symmetry) | 2.54 % | 0.0230 |
| 0.33 μm, $x$ (preprint experiment) | 2.13 % | 0.0192 |
| **0.8 μm, $x$ (published experiment)** | **3.02 %** | 0.0342 |
| 0.8 + 0.4, $x$ (ratio 2) | 1.40 % | 0.0173 |
| 0.8 + 0.55, $x$ | 0.50 % | 0.0070 |
| 0.8 $x$ + 0.8 $y$ | 0.09 % | 0.0027 |
| {0.8, 0.55} $x$ + $y$ | 0.00 % | 0.0010 |

Of the 3.02 % lost by the published single shear, **83.6 % sits on the $m\neq0$ null
planes and 16.4 % on the blind plane**; nothing is lost elsewhere. The 1D proxy in
[[finite-shear-null-experiment]] said 3.56 %. **The 3D answer is 3.02 %, and it is a
different 3.02 %** — one sixth of it is the blind plane, which no magnitude choice
fixes and which the 1D calculation could not see.

All rows use GODT's *simulation* optics ($\lambda=0.5$ μm, NA 1.33, 60°). The 0.33 and
0.8 μm rows put experimental shears onto that geometry. The experiment itself ran at
530 nm with an NA 1.45 objective and 58° illumination in air, so their lost
volumes are indicative.

Two things the 1D proxy also missed. The 0.1 μm simulation shear has **no in-band
nulls yet still loses 2.54 %**. All of it sits in the low-$|G|$ slab around the blind
plane, and only 0.50 % on the plane itself. The reason is that a small shear has a *wider*
low-$|G|$ slab around $K_x=0$ because $|G|^2\approx(K_x\Delta)^2$ there. And 0.33 μm
beats 0.1 μm despite having a null, for the same reason. Single-shear quality is a
trade between a null comb and a blind slab, not a monotone function of $\Delta$.

## 4. Optimal designs (T3)

Candidate grid: 32 magnitudes (0.05–1.6 μm) × 2 directions = 64 channels, GODT
simulation geometry, $\rho=100$. Convexity holds as the problem note argued;
optimality is certified by the equivalence-theorem gap
$(\sum_jt_jg_j-\min_jg_j)/|f|$ (derived in [[multishear-fusion-derivation]] §7.2).
It is $\le3\times10^{-4}$ for the (P-A)/(P-D) rows of this table. (P-W) is solved as an
LP. The E1 solves in §6 use a looser tolerance of $10^{-3}$ and do not print their gap.

| design | $f_A$ (mean rel. MSE) | worst-case efficiency |
|---|---|---|
| ODT benchmark | 0.00157 | — |
| best single channel, $\Delta=0.20$ μm | 0.01224 | 0 |
| 0.8 μm single (published) | 0.03420 | 0 |
| best $x\!+\!y$ pair, $\Delta=0.15$ μm | 0.00089 | — |
| {0.8, 0.55} $x$ + $y$ | 0.00097 | 0.007 |
| {0.8, 0.55, 0.35} $x$ + $y$ | 0.00086 | 0.070 |
| uniform over all 64 | 0.00083 | 0.175 |
| **(P-A) optimum, 24 channels** | **0.00076** | 0.096 |
| (P-D) optimum, 20 channels | 0.00076 | — |
| (P-W) optimum, 28 channels | 0.00081 | **0.233** |

Read the benchmark carefully: $|G|^2$ peaks at 4, so on this model a shear design can
score *better* than direct ODT. That is assumption A4 plus unit visibility talking,
not a claim that gradient beats direct ODT. §6 prices the visibility.

**The design rule worth remembering.** $2\pi/K_{\max}=0.201$ μm is the largest shear
whose first null $2\pi/\Delta$ stays outside the measurable band. On a 0.05 μm
candidate grid the largest such value is 0.20 μm — and that is exactly what the
optimiser picks as the best single channel. Choose $\Delta\lesssim2\pi/K_{\max}$ and
the comb never enters. **[audit]** "$\lesssim$" should read "$\approx$". Below the line,
information falls again as the blind slab widens: the 0.1 μm simulation shear, "comfortably
inside", scores 0.0230 against 0.0122 at the line. The experimental 0.8 μm sits ~4× above it
on these optics; the paper gives practical reasons for a larger shear (noise, hardware)
that assumption A4 does not model. *(The grid resolves to 0.05 μm, so this identifies the
mechanism — largest null-free shear — not three significant figures. The code's
reachable edge, clipped by the pupil floor, is 31.10 rad/μm, giving 0.202 μm.)*

**The gains are lopsided.** Going from the published
single shear to a well-chosen pair is 35×. Going from a well-chosen hand-picked set to
the certified optimum is 12 %, and a *uniform* design over all 64 candidates is
within 9 % of optimal. The objective is flat near its minimum. In short,
multi-shear matters and the optimisation is a refinement — except in the
worst-case metric, where (P-W) buys 2.4× over (P-A) for 7 % of the mean.

### C1 (sparsity) is false

The optima use 24, 20 and 28 of 64 candidates and stay spread at $\rho=10$ and
$\rho=1000$. No Carathéodory-style collapse. None was guaranteed: the information here is a
function over every voxel, not a small matrix, so the theorem gives no useful bound
([[multishear-fusion-derivation]] §7.3). The mechanism is visible in the flat
objective: with $|G_j|^2$ overlapping heavily between neighbouring magnitudes, many
designs are near-optimal and the optimiser has no incentive to concentrate.

### Extra directions are nearly free and nearly useless

Adding 45° and 135° (128 candidates) improves $f_A$ from 0.00076 to 0.00073 — **4 %**
— while the diagonals take half the weight. A theory-side confirmation of the
Zhu group's experimental saturation ([[zhu2024-sensors]] ~0.7 % for a third
direction, [[zhu2025-sci-rep]] no change to 4-digit precision from 4 to 16). It is
reached independently and for a different reason: one orthogonal pair already covers
every blind plane. The analogy is loose *(inference)*. The Zhu group's Zernike fits are
null-immune and their gains are error averaging. [[zhu2024-appl-phys-b]] also reports a
40.7 % simulated gain from 2 to 4 directions, which this model does not reproduce.

### C2 (coherence-limited optimum) confirmed — and it is what creates sparsity

With $|\mu(\Delta)|^2$ a van Cittert–Zernike jinc of coherence length $l_c$ (a
**model, not measured** — see §8 and [[partially-coherent-transfer-function]]):

| $l_c$ | channels in the (P-A) optimum | largest magnitude used | $f_A$ |
|---|---|---|---|
| ∞ | 24 | 1.60 μm | 0.00076 |
| 4 μm | 14 | 0.80 μm | 0.00079 |
| 2 μm | 8 | 0.55 μm | 0.00082 |
| 1 μm | 6 | 0.35 μm | 0.00090 |

Coherence, not information structure, is what makes an optimal design sparse and
small-shear. The best *pair* magnitude stays 0.15 μm throughout, so the headline
design rule is coherence-robust even though the full design is not.

### C3 (incommensurability emerges) confirmed, with a scope limit

Scanning a second $x\!+\!y$ pair at $d_2\in[0.05,1.6]$ μm (step 0.01) against a fixed
0.8 μm pair, the local *maxima* of $f_A$ — the worst choices — land on small-integer
ratios $0.8/d_2\in\{3,2,\tfrac32,1,\tfrac34,\tfrac23,\tfrac35\}$ (**[audit]** the ratio is
$0.8/d_2$, not $d_2/0.8$). That is $d_2=0.27,0.40,0.53,0.80,1.07,1.20,1.33$ μm, each
within one grid step of the exact ratio. $\min S_{\mathbf t}=0$ at ratios 1, 2 and 2/3.
Nothing imposed the rule; it fell out, and the lattice argument in
[[multishear-fusion-derivation]] §5 predicts it. Best second pair:
0.18 μm (ratio 4.44), $f_A$ 0.00080 against 0.00275 for the 0.8 pair alone.

**But the rule is a rule for small sets.** The (P-A) optimum contains ratios 2, 3 and
4 among its twelve magnitudes and still has $\min S_{\mathbf t}=0.44$, because a null
shared by two channels is covered by the other ten. State the incommensurate-ratio
rule for two- and three-channel designs, which is where it was derived and where it
bites.

## 5. Validation on a phantom (T4)

Textured cell phantom (cytoplasm, nucleus, nucleolus, 25 organelles, 0.06 μm
correlated texture; $\delta n\in[-0.013,+0.085]$), $128^3$ at 0.1 μm isotropic,
720 azimuths, digital caps, complex Gaussian noise per sample, seeded.

Two losses are quoted once. The **missing cone** costs 0.431 of the truth's norm. The
**$|K|<2$ rad/μm shell** holds 0.842 of it, and every estimate is masked there,
including the ODT benchmark (`K_LO` in exp05). All aggregates below therefore exclude
it. This is the trap that inverted the first attempt at `exp02`, and it is recorded in
`lab/README.md`.

**[audit]** An earlier version said the shell is lost because "$|G|^2\to0$ for every
shear". That is wrong. The mask is a choice, and the shell is only partly lost:

- The shell has 241 voxels. The **DC voxel** holds 19.4 % of the total power and is blind
  for every shear.
- The rest have $K_\perp\in[0.49,1.96]$. There, 0.8 μm has median $|G|^2=0.15$ (max 2.0),
  and only 0.1 μm is near zero.
- Unmasked Wiener error in the shell: ODT 0.000, either single shear 0.55, both pairs
  0.274. The 0.274 is exactly DC's share, so a pair recovers the whole shell except
  DC.

So gradient data lose DC (the mean index, recoverable from any absolute reference) and
a single shear loses about half the low shell. The mask hides a real ODT advantage and
a real pair advantage alike.

| design | fusion | RMSE vs band-limited truth |
|---|---|---|
| ODT benchmark | Wiener | 0.0100 |
| 0.1 μm $x$ | Wiener | 0.2194 |
| 0.8 μm $x$ | Wiener | 0.2564 |
| 0.8 μm $x$ | BLUE | **diverges** |
| 0.8 $x$ + 0.8 $y$ | Wiener | 0.0424 |
| 0.8 $x$ + 0.8 $y$ | naive average | 0.1824 |
| {0.8, 0.55} $x$ + $y$ | Wiener | **0.0078** |
| {0.8, 0.55} $x$ + $y$ | BLUE | 0.0078 |
| {0.8, 0.55} $x$ + $y$ | naive average | 0.1607 |
| {0.8, 0.55} $x$ + $y$ | precision-weighted average **[audit]** | 0.0079 |
| 0.8 $x$ + 0.8 $y$ | precision-weighted average **[audit]** | 0.0424 |
| (P-A) optimum (exp05's own: 22 channels, shell prior; not §4's) | Wiener | 0.0070 |

The precision-weighted rows use weights $J_j+S_O^{-1}$ per voxel on each channel's Wiener
estimate (`lab/exp05b_weighted_average.py`, same seed and data as exp05).

**Three results, in order of how much they matter.**

*Fusion beats equal-weight averaging by 21×* (0.0078 against 0.1607) on identical data at
equal photon budget. **[audit] Rewritten.** The earlier text claimed that averaging
per-channel reconstructions *cannot* recover what each channel lacks. That is false.
Each per-channel Wiener estimate is shrunk to zero at its own nulls, and an
*equal-weight* average dilutes it there by $1/N$. A per-voxel precision-weighted
average of the same reconstructions matches fusion to 1 %. Algebraically, the fused
Wiener estimate *is* the inverse-variance average of per-channel BLUEs
([[multishear-fusion-derivation]] §6.2). What the 21× measures is how much a practitioner
loses by averaging reconstructions or phase maps without per-frequency weights. That is
a useful, quotable number, but it is not a new estimator: the estimator is standard
information-form fusion. The 3D tomographic extension of the [[zhu2024-appl-phys-b]]
hint (6× for coefficients over wavefronts) stands, and [[multiple-scattering-check]]
shows how it shrinks with contrast.

*BLUE diverges for a single shear and is fine for an incommensurate pair.* The
prior-free unbiased estimator is undefined where $J=0$ and blows up beside it; with a
pair, $J>0$ across the band and BLUE and Wiener agree to four digits. This is the
estimator-side statement of why multi-shear is not merely noise averaging — it is the
difference between an estimable and a non-estimable parameter.

*The theory is right.* For a Gaussian random field drawn from the prior — where
$(J+S_O^{-1})^{-1}$ should be exact, not approximate — empirical and predicted error
agree to a shell-median ratio of 1.003–1.011, with 96–98 % of shells within 5 % for
designs that leave no near-null voxels. Where a design does leave them (0.8 $x$ + 0.8
$y$ has doubly-null voxels), one shell scatters to 2.25 because 99 % of its predicted
error sits on the top 1 % of its voxels — a small-sample estimate of an exponential
variable, not a model failure.

## 6. Robustness (T5)

### E1 — shear calibration error

With rms shear error $s$ treated as coherent model error, the per-voxel SNR is capped
at $|G|^2/(s\,K_\parallel)^2$ regardless of photons, where $K_\parallel$ is the
projection of $\mathbf K_\perp$ on the shear direction. It is quoted below along
$\hat{\boldsymbol\Delta}$, where $K_\parallel=K_\perp$ (derivation in
[[multishear-fusion-derivation]] §8.1).

**[audit] Different setup from §4.** The grid is $128^2\times48$ (a 12.8 μm field), with
a binned density from 1440 azimuths. The candidates are 16 magnitudes (0.1–1.6 μm) ×
$x,y$ = 32, which is what "32 (all)" means. $\rho$ is per frame-sample and $M\propto\delta k$, so the coarser $K$ pitch doubles
$\rho M$ and halves $f_A$. The $s=0$ value 0.00038 therefore corresponds to §4's
0.00076 and is not a better design. Compare rows within this table only.

| $s$ | ceiling at $K_\perp=30$ | $f_A$ | channels in the optimum |
|---|---|---|---|
| 0 | ∞ | 0.00038 | 20 |
| 5 nm | 22.5 dB | 0.00048 | 32 (all) |
| **10 nm = 0.1 px** | **16.5 dB** | 0.00069 | 32 |
| 20 nm | 10.5 dB | 0.00135 | 32 |
| 50 nm | 2.5 dB | 0.00536 | 32 |

The 50 nm row has $sK_\perp=1.5$, outside the first-order linearisation, and is
indicative only.

Photon-limited SNR on the same grid is 37.1 dB (about 34 dB on §4's grid). **At 10 nm
rms error, calibration is 20 dB below the photon limit at the band edge and is what sets
high-frequency quality.** **[audit]** [[2509-15382-osi-flex]]'s "0.1 px" (§3.2) is a
*systematic* offset on a 5.3-pixel shear, about 1.9 %, and its pixels have no stated
physical size. "10 nm" converts it on GODT's 0.1 μm simulation pixel *(inference)*.
The source reports a systematic error, and this table models an rms random one.

The model also treats each channel's single calibration error as independent noise at
every voxel. That overstates how much spreading over channels averages it down,
since with $N$ channels there are only $N$ independent draws *(inference)*. This is
[[multi-shear-acquisition]]'s "calibration is the real cost, not photons" arriving
independently from the information side. The optimal response is to *spread* the
budget over every available channel rather than concentrate it — the opposite of what
the photon-limited optimum does.

### E2 — chromatic shear under GODT's real LED

The polarization-grating shear scales with wavelength. This is now **confirmed from
the source's own numbers** rather than assumed: [[s42005-025-02416-3]] reports
$\Gamma=6.29$ μm and $\theta\approx4.83^\circ$ at $\lambda=530$ nm, and
$\arcsin(0.53/6.29)=4.834^\circ$ — the grating equation $\sin\theta=\lambda/\Gamma$
holds exactly, so $\boldsymbol\Delta(\lambda)\approx\boldsymbol\Delta_0\lambda/\lambda_0$.

The source gives only "$\Delta\lambda=35$ nm". Reading it as a Gaussian of 35 nm FWHM is a
modelling choice, and the design geometry stays at the simulation's 0.5 μm. Under that
spectrum the wavelengths add in intensity, so the
effective filter is the spectral *average* of $G$. At the null
$K=2\pi m/\Delta_0$ the residue is second order in the fractional bandwidth,
$|G_\text{eff}|\approx2\pi^2m^2\sigma_\varepsilon^2$:

| null | $\lvert G_\text{eff}\rvert^2/4$ at 35 nm | 2nd-order prediction | rel. MSE at 28 dB |
|---|---|---|---|
| $m=1$, 7.85 rad/μm | $5.7\times10^{-5}$ | $6.0\times10^{-5}$ | 0.872 |
| $m=2$, 15.71 | $8.7\times10^{-4}$ | $9.6\times10^{-4}$ | 0.307 |
| $m=3$, 23.56 | $4.1\times10^{-3}$ | $4.9\times10^{-3}$ | 0.086 |

**Partial answer to the open question on broadband illumination.** The nulls survive —
the first one is still 87 % lost — but the higher-order ones are substantially filled.
The $m$-th null's position spreads by $2\pi m\sigma_\varepsilon/\Delta_0$, which grows
with $m$, while its width does not, so the residue grows as $m^4$
([[multishear-fusion-derivation]] §8.2). Net effect on the published configuration:
single-shear 0.8 μm improves from $f_A$ 0.03420 to 0.02115. That is a **1.6× gain that comes for free
from GODT's broadband LED**. The optimal design barely moves (0.00078 against 0.00076).

Caveat, and it is the same one the wiki already carries: this treats the broadband
source as a spectral average of the coherent monochromatic model. Whether $G_j$
remains multiplicative under [[ledwig-robles2023]]'s Ewald-*plane* treatment of a
genuinely partially coherent source is still open — this computes the chromatic half
of the question, not the spatial-coherence half. The model also holds three things
fixed across the spectrum: the Ewald radius $k_m(\lambda)$, the object phase
($\propto1/\lambda$) and the Rytov weight. Each wavelength really samples a slightly
different cap, displaced by $\sim\sigma_\varepsilon K$ *(inference)*.

**Follow-up, 2026-10-05.** [[colour-multiplexed-shear]] takes the "incoherent bound" of
`chromatic_G` seriously: resolve the colours instead of averaging them, with
per-wavelength geometry. Three LED colours then act as three shear magnitudes and match
the mechanical {0.8, 0.55} μm pair of §3.

## 7. What this does *not* establish

- **$\mu(\Delta)$ is a model.** The jinc in §4 is van Cittert–Zernike for a circular incoherent source; GODT's real LED-plus-condenser visibility against shear is not measured anywhere in this wiki. The *direction* of the C2 effect is robust; the numbers are not.
- **A4 is doing a lot of work.** That phase-difference noise is independent of $|\Delta|$ is what lets $|G|^2\le4$ beat the ODT benchmark. It is right for a common-path shearing interferometer at fixed visibility and wrong the moment visibility enters, which is exactly what C2 shows.
- **A6 (sequential channels).** E3 of the problem note — correlated channels under single-shot multiplexing — is not computed. It is the acquisition-side question and [[2605-04297-double-sagnac-ls-qpm]]'s cascade is the case that matters.
- **Spatially varying shear (E6)** is untouched and remains the cheapest open item; the null planes move with $\Delta(\mathbf x)$ and the design is built on them being fixed.
- ~~**Linear Rytov only.**~~ **Done 2026-09-22 — see [[multiple-scattering-check]].** The WPM comparison the problem note's T4 asked for now exists. Short version: the *geometry* of this note is exact under multiple scattering (the shear factor multiplies the measured field, so the nulls survive to 4e-16), the *accuracy* is not (on its phantom the fusion-beats-averaging margin falls from 9.6x (linear) to 1.1x at [[chen1998]]'s validity bound, and to nothing beyond it). And nearest-neighbour gridding (A3, below) turns out to be the single largest *correctable* error in the pipeline: a factor 6.2, against the further 2.7x multiple scattering adds at the lowest contrast tested.
- **Nearest-neighbour gridding (A3).** The digital caps *are* the model, so interpolation error is assumed away in theory and simulation alike. The 76.6 % angular coverage number is the visible edge of this. **Now priced** in [[multiple-scattering-check]]: against data carrying the true continuum $K_z$, nearest-neighbour deposition costs a factor **6.2** on the best design (0.0171 -> 0.1058 RMSE), against the further 2.7x multiple scattering adds at the lowest contrast tested — exactly what [[pan1983]] warned nearest-neighbour would do.

- **The 60° illumination angle is taken as measured in the medium** ($k_{\rm ill}=k_m\sin\theta$). [[s42005-025-02416-3]] says only "with respect to the optical axis" *(inference)*. $K_{\max}$ and the design rule inherit this.
- **Rotation (E4) and phase-only coupling (E5/A2)** are not computed. A2 is argued harmless for the design ($|G|^2$ is even in $\mathbf K_\perp$; [[multishear-fusion-derivation]] §1.4), but the $2\times2$ version has not been run ([[krauze2016]]).
- **Wrapping** is excluded by A4. [[multiple-scattering-check]] §6 gives a third, independent argument for small $\Delta$ from the wrapping limit across the shear. It reinforces the $\Delta\approx2\pi/K_{\max}$ rule and should be cited with it.
- **Validation protocol only partly followed.** The problem note's §2.8 asked for SSIM, axial elongation, a Shepp–Logan phantom and a nano-printed phantom. Only RMSE on a textured cell is reported.
- **Real samples.** Whether any of these margins apply to GODT's actual samples depends on their $a\,\delta n/\lambda$.

## Reproduction

```bash
pip install -r lab/requirements.txt          # numpy, scipy, matplotlib
cd lab
python shearlib.py                # selftest, includes the 3D invariants
python exp03_ewald_fisher.py      # §1-3, ~20 s, caches M(K) in out/
python exp04_shear_design.py      # §4, 6, ~4 min
python exp05_fusion_validation.py # §5, ~2.5 min
python exp05b_weighted_average.py # §5 [audit] rows and low-K shell, ~1 min
```

Figures land in `lab/out/` (gitignored): `ewald_maps.png`, `lost_volume.png`,
`design_shares.png`, `design_scan.png`, `design_maps.png`,
`validation_spectra.png`, `validation_slices.png`.
