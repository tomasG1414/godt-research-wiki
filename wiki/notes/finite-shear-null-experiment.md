---
type: note
projects: [godt]
created: 2026-07-18
updated: 2026-10-20
status: current
---

# Finite-shear transfer nulls: experiment

Tests the inference that GODT's finite-shear transfer factor nulls a comb of spatial frequencies, and that combining shears could fill those nulls. **Result: confirmed, with two corrections to the original framing.**

Everything below is computed, not sourced. The one thing taken from a source is [[s42005-025-02416-3]]'s Eq. 20; all consequences are derived here. **Note:** the forward shearing filter and its nulls are established 2D prior art — [[servin2007]] derives $2i\sin(dx\,\omega_x)\,\Phi$ (their Eq. 6), plots the inverse's poles (Fig. 2), and states the information loss explicitly — and combining shear magnitudes with a coprime ratio is also known in 2D lateral shearing (Elster & Weingärtner 1999 and follow-ups, summarised in Tomioka et al., *Appl. Opt.* 2023, arXiv:2309.16434). What this page adds is the quantification at GODT's shear values in the tomographic setting.

## Why this is testable without the tomographic pipeline

GODT's relation is $\tilde{O}_\Delta(K) = H_\Delta(K)\,\tilde{O}(K)$. The shear factor multiplies the object spectrum **pointwise**. The Ewald/DI machinery determines *which* $K$ get sampled; it does not change $H$. So the conditioning of the shear inversion is a property of $H$ alone and separates cleanly from the tomography. That is what makes this a cheap experiment rather than a full reconstruction build.

## 1. The transfer factor, in closed form

From Eq. 20, $H_\Delta(K) = \frac{e^{jK\Delta}-1}{\Delta}$. Using $|e^{i\theta}-1| = 2|\sin(\theta/2)|$:

$$\boxed{\;|H_\Delta(K)| = \frac{2\left|\sin\!\left(K\Delta/2\right)\right|}{\Delta}\;}$$

Verified against direct numerical evaluation to about $10^{-14}$, and confirmed to converge to $|iK|$ as $\Delta\to0$ (relative deviation $10^{-3}$ at $\Delta=10^{-4}\,\mu$m). Zeros at $K = 2\pi m/\Delta$.

### Correction 1: there are two different zero structures

The original inference conflated them.

- **$m=0$, i.e. $K=0$.** The DC/low-frequency loss inherent to *any* gradient method. **Independent of $\Delta$** — no choice of shear magnitude touches it. This is the thing gradient methods trade away, and which [[s42005-025-02416-3]] presents as an advantage ("inherent suppression of slowly varying backgrounds"). Recovering it is gradient integration, which GODT is designed not to need.
- **$m\neq0$, i.e. $K = 2\pi m/\Delta$.** The finite-shear nulls. These **move with $\Delta$** and are the ones multi-shear can fill.

Only the second was the real content of the original claim.

## 2. The finding: the nulls depend strongly on the shear value

Object-spectrum band for the GODT simulation configuration (λ=0.5 μm, n₀=1.33, NA=1.33, 60° illumination): $|K_x| \le k_{det} + k_{ill} = 16.71 + 14.47 = 31.19$ rad/μm, which matches the sampling Nyquist of 31.42 rad/μm at dx=0.1 μm — critically sampled.

| shear | source of value | first null | nulls inside band | positions (rad/μm) |
|---|---|---|---|---|
| 0.1 μm | GODT **simulation** (1-pixel shift) | 62.83 | **0** | none |
| 0.33 μm | GODT **preprint** experiment | 19.04 | 1 | 19.04 |
| 0.8 μm | GODT **published** experiment | 7.85 | **3** | 7.85, 15.71, 23.56 |

**At the 0.1 μm simulation shear there are no in-band nulls; at the 0.8 μm experimental shear there are three.** The published paper notes that Δ=0.8 μm "slightly smooths high-frequency details". That smoothing is the roll-off of the finite-difference filter: $|H_\Delta|/|K|=|\mathrm{sinc}(K\Delta/2)|$ falls to $2/\pi$ at $K=\pi/\Delta$, with sign-reversed side lobes beyond. The nulls are where that roll-off crosses zero, and they matter for any reconstruction that tries to *invert* the filter (see the caveat in "What this does not establish").

## 3. What the nulls cost, quantitatively

For the optimal linear (Wiener) estimator with white noise and per-frequency SNR $\rho$, the relative recovery error is

$$\varepsilon(K) = \left[\rho\,G^2(K) + 1\right]^{-1/2}, \qquad G^2(K) = \sum_n |H_n(K)|^2$$

Validated against Monte Carlo (40 000 trials/frequency) to 3–4 decimals. At a null $G=0$, so $\varepsilon = 1$ and the frequency is lost at any SNR (with monochromatic light). Better hardware cannot fix this: the transfer factor is exactly zero there, so the data carry no information about that frequency.

At ρ = 100 (20 dB), errors quoted for K > 2 rad/μm to exclude the shared low-frequency rise:

| shear set | worst error | at K | % of band with ε>0.5 | mean ε |
|---|---|---|---|---|
| 0.1 (sim) | 0.050 | — | 0.00 % | 0.010 |
| 0.33 (preprint exp) | 0.9998 | 19.04 | 1.18 % | 0.045 |
| **0.8 (published exp)** | **0.9999** | 23.56 | **3.56 %** | 0.108 |
| 0.8 + 0.4 (ratio 2) | 0.9994 | 15.71 | 0.85 % | 0.044 |
| 0.8 + 0.55 | 0.195 | 23.20 | **0.00 %** | 0.038 |
| 0.8 + 0.08 | 0.037 | — | 0.00 % | 0.009 |
| 0.8 + 0.55 + 0.35 | 0.031 | — | 0.00 % | 0.020 |

Rows with small shears are flattered by the $H=G/\Delta$ normalisation used on this page; see the note on ranking below, and [[multishear-fusion-design-solution]] §3 for the physical ranking.

### Correction 2: commensurate shear ratios do not help

**0.8 + 0.4 barely improves on 0.8 alone** (0.85 % vs 3.56 % of band lost) because a ratio of exactly 2 makes the nulls coincide at K = 15.71. The original claim — "combine different shear magnitudes" — is too loose. The design rule is:

> **Choose shear magnitudes whose ratio is not a ratio of small integers.** 0.8 + 0.55 (ratio ≈ 1.45) removes every in-band null; 0.8 + 0.4 does not.

A sweep of second shears against a fixed 0.8 μm found the worst choices are 0.4 and 0.8 (coincident nulls, $G \to 0.0002$). **Note on ranking:** this page uses $H=G/\Delta$ with noise independent of $\Delta$, which weights small shears up, so its preference for very small partners (0.05–0.09 μm) is an artifact of that normalisation. Under the physical noise model of the 3D study the best partner for 0.8 μm is about 0.15–0.2 μm ([[multishear-fusion-design-solution]] §4). The null positions and the commensurability rule are unaffected.

## 4. Image-domain check

512², dx=0.1 μm, band-passed 2–31 rad/μm so the shared DC loss is excluded and the nulls are what shows. Shepp-Logan object, 20 dB, Wiener recovery. Relative RMSE:

| set | RMSE | vs single 0.8 |
|---|---|---|
| 0.8 only | 0.1533 | — |
| 0.8 + 0.4 | 0.1052 | 1.5× better |
| 0.8 + 0.55 | 0.0964 | 1.6× |
| x + y at 0.8 | 0.0448 | 3.4× |
| x + y at 0.8 and 0.55 | **0.0180** | **8.5×** |

Visible artifact in the single-0.8 reconstruction is directional banding along the shear axis, consistent with losing discrete $K_x$ frequencies.

## 5. Directions vs magnitudes — and the missing-cone overlap

$H$ depends on $K\!\cdot\!\mathbf{d}$, so a shear along $x$ is blind wherever $K_x=0$ regardless of magnitude. Fraction of the 2D band with $G<0.5$:

| shear set | min G (\|K\|>0.5) | % band with G<0.5 |
|---|---|---|
| x only, 0.8 | 0.000 | 12.13 % |
| x only, 0.8+0.55 | 0.000 | 1.15 % |
| x+y, 0.8 | 0.010 | 1.14 % |
| **x+y, 0.8+0.55 each** | **0.505** | **0.00 %** |
| 3 dirs (0/60/120°), 0.8 | 0.013 | 0.64 % |

So **directions and magnitudes complement each other**: a second direction covers the $K_x=0$ blind line and most of the first direction's nulls (12.13 % → 1.14 %), a second magnitude fills the $m\neq0$ comb along one direction (→ 1.15 %), and only both together clear the band. **[audit 2026-10-20]** An earlier reading credited directions with the blind line only; see [[multishear-fusion-design-solution]] §2–3.

**The finding worth carrying forward:** no set of *in-plane* shears recovers the $K_x=K_y=0$ axis — and that axis is exactly where the [[missing-cone-problem]] already sits. The shear blind spot and the missing cone overlap, which partly defuses the shear blind spot (those frequencies were already unmeasured). Checked against the full 3D Ewald sampling in [[multishear-fusion-design-solution]] §2.

## What this does and does not establish

**Does:** the null structure is real, exactly where the closed form predicts, falls inside the measurable band at GODT's published experimental shear, costs frequencies irrecoverably at any SNR (for monochromatic light), and is removed by adding an incommensurate second shear. Multi-shear now has a second quantitative motivation independent of directional contrast.

**Does not:** this models the shear filter in isolation. It does not include the tomographic sampling (which $K$ are actually visited by the Ewald caps), the [[rytov-approximation]]'s validity limits, amplitude effects, or realistic experimental noise structure. The percentages are along the 1D shear axis; in 3D the null set is a family of planes perpendicular to the shear direction, so the volume fraction differs. And recovering $O$ is a different task from what GODT publishes — GODT displays the gradient field itself, so for single-shear GODT-as-published the nulls mean the displayed image simply lacks those frequencies rather than that an inversion is ill-conditioned.

**Follow-up:** done in 3D with the actual Ewald sampling; see [[multishear-fusion-design-solution]].

## Reproduction

Minimal self-contained reproducer for the load-bearing numbers (needs only numpy):

```python
import numpy as np
LAM, N0, NA, TH = 0.5, 1.33, 1.33, 60.0                     # GODT sim config
KMAX = 2*np.pi*NA/LAM + 2*np.pi*N0*np.sin(np.deg2rad(TH))/LAM   # 31.19 rad/um

aH = lambda K, d: 2*np.abs(np.sin(K*d/2))/d                 # |H| closed form
G2 = lambda K, ds: sum(aH(K, d)**2 for d in ds)             # combined power
err = lambda K, ds, rho: np.sqrt(1/(rho*G2(K, ds)+1))       # Wiener rel. error

for d in (0.1, 0.33, 0.8):                                  # in-band nulls
    n = [2*np.pi*m/d for m in range(1, 99) if 2*np.pi*m/d <= KMAX]
    print(f"d={d:4.2f}  {len(n)} null(s): {[round(v,2) for v in n]}")

K = np.linspace(0.01, KMAX, 6000); m = K > 2.0              # skip shared DC rise
for lbl, ds in [("0.8", [0.8]), ("0.8+0.4", [0.8, 0.4]),
                ("0.8+0.55", [0.8, 0.55]), ("0.8+0.08", [0.8, 0.08])]:
    e = err(K, ds, 100.0)                                   # rho=100 => 20 dB
    print(f"{lbl:9s} worst={e[m].max():.4f}  %band>0.5={100*np.mean(e[m]>0.5):.2f}")
```

Expected output: `d=0.80  3 null(s): [7.85, 15.71, 23.56]`, and `0.8` worst≈1.0 / 3.56 % versus `0.8+0.55` worst≈0.19 / 0.00 %.

The full experiments now live in `lab/` (added as schema v2):

```bash
pip install -r lab/requirements.txt
cd lab
python shearlib.py                # selftest on the shared math
python exp01_shear_nulls.py       # sections 1-5 above, writes out/transfer.png
python exp02_shear_recovery.py    # Monte Carlo, error curves, image demo
```

Both reproduce every number on this page exactly. `lab/README.md` records what
is established and the gotcha that sank the first attempt (see "Gotchas hit already" there).
