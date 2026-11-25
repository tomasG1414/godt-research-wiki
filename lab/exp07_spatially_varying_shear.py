"""exp07 -- Spatially varying shear: does a smooth, known Delta(r) fill the nulls?

Relaxes assumption A5 (one uniform shear per channel). The data model is

    psi(r) = phi(r + Delta(r)) - phi(r) + noise,     Delta(r) = Delta0 (1 + eps f(r)) x^

with f a smooth profile of unit range across the field: a linear ramp, or a
quadratic (distortion-like) bowl. eps = 1 %, 3 %, 10 % of 0.8 um. Projection
level (1D/2D) like exp01/exp02: no Ewald caps, no Rytov factor.

Two estimators:
  calibrated -- knows Delta(r); MAP/Wiener for the true space-variant operator
                (dense linear algebra in 1D, conjugate gradient in 2D);
  naive      -- assumes a uniform shear and applies the per-frequency Wiener
                filter of exp02.

Establishes (all for a flat band-limited prior, per-frequency SNR rho):
  0. For uniform shear the dense machinery reproduces wiener_rel_err exactly.
  1. The null moves rather than vanishes: chirps exp(2 pi i m alpha(x)) with
     local frequency 2 pi m / Delta(x) (alpha = Abel function of x -> x +
     Delta(x)) are as nearly null as the uniform shear's plane waves, to the
     window-limited precision of the test. (At 40 dB the periodic ramp does
     lose fewer modes than the uniform shear -- section 2b -- which the jump
     of the ramp at the periodic wrap may explain; not isolated.)
  2. Calibrated: the error per FOURIER frequency near a null drops, the number
     of lost modes does not. Per-frequency error = field average of the local
     Wiener error (checked), i.e. the null is smeared over K_m * eps.
  3. A windowed (Gabor) information map: the null follows 2 pi m / Delta(x0);
     its width over a window W grows as K_m |grad Delta| W / Delta0.
  4. Naive: model error, exact bias/noise split, versus amplitude and SNR.
  5. Monte Carlo check of 2 and 4 with real objects and real noise.
  6. 2D image demo (band-passed Shepp-Logan), CG for the calibrated estimator.
  7. Where diversity does come from: the same object seen at two field
     positions (sample translated across the shear gradient).

Convention: the DIMENSIONLESS shear power |e^{iK Delta} - 1|^2 with noise
independent of Delta (assumption A4, as exp03+), not exp02's H = G/Delta.
Uniform-shear numbers therefore differ slightly from exp02's table.
All "pair" rows state whether the two channels share the exposure (rho/2
each, as exp03/exp04) or get a full exposure each (as exp02).

Stochastic parts are seeded (SEED). Sections 0-4 and 7 are deterministic
linear algebra with no random numbers.

Writes: out/varying_shear.png, out/varying_shear_images.png
Wiki:   wiki/notes/spatially-varying-shear.md

Run: python exp07_spatially_varying_shear.py        (~1-2 min, < 1.5 GB)
"""
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from shearlib import GODT_SIM, G2_vec, nulls_in_band, shepp_logan_2d, out_path

SEED = 7
RHO = 100.0            # per-frequency SNR, 20 dB (as exp02/exp03)
D0 = 0.8               # nominal shear, um (published experiment)
D2 = 0.55              # incommensurate partner used throughout the project
DX = 0.1               # pixel, um (GODT simulation sampling)
KMIN = 2.0             # band-pass lower edge, rad/um (README gotcha)
KMAX = GODT_SIM.kmax
EPS = (0.01, 0.03, 0.10)
CORE, WIDE = 0.25, 1.5     # half-widths of the null bands, rad/um
RULE = "=" * 78


# ------------------------------------------------------------------ 1D model
class Field1D:
    """Periodic 1D field of N pixels; unknowns are the in-band Fourier modes."""

    def __init__(self, n, dx=DX, kmin=KMIN, kmax=KMAX):
        self.n, self.dx, self.L = n, dx, n * dx
        self.x = (np.arange(n) - n // 2) * dx
        k = 2 * np.pi * np.fft.fftfreq(n, dx)
        self.k = np.sort(k[(np.abs(k) >= kmin) & (np.abs(k) <= kmax)])
        self.dk = 2 * np.pi / self.L
        self.U = np.exp(1j * np.outer(self.x, self.k)) / np.sqrt(n)   # unitary columns

    def profile(self, name):
        """Unit-range profile f(x): Delta = D0 (1 + eps f)."""
        u = self.x / self.L                       # [-1/2, 1/2)
        if name == "ramp":                        # tilt: -eps/2 ... +eps/2
            return u
        if name == "quad":                        # bowl: 0 at centre, 1 at the edges
            return (2 * u) ** 2
        raise ValueError(name)

    def B(self, delta):
        """Operator band coefficients -> data, psi = phi(x + Delta(x)) - phi(x)."""
        return self.U * (np.exp(1j * delta[:, None] * self.k[None, :]) - 1.0)

    def band(self, K, hw):
        return np.abs(np.abs(self.k) - K) < hw


def posterior_relmse(Bs, rhos):
    """diag of (sum_j rho_j B_j^H B_j + I)^-1: relative MSE per Fourier mode of
    the MAP estimate, flat unit prior. Also returns the Fisher matrix."""
    F = sum(r * (B.conj().T @ B) for B, r in zip(Bs, rhos))
    C = np.linalg.inv(F + np.eye(F.shape[0]))
    return np.real(np.diag(C)), F, C


def naive_relmse(fld, B, d_assumed, rho):
    """Exact relative MSE per mode of the uniform-shear Wiener filter applied to
    data from the true operator B. Returns (bias part, noise part)."""
    G0 = np.exp(1j * fld.k * d_assumed) - 1.0
    W = np.conj(G0) / (np.abs(G0) ** 2 + 1.0 / rho)
    M = W[:, None] * (fld.U.conj().T @ B) - np.eye(len(fld.k))
    return np.sum(np.abs(M) ** 2, axis=1), np.abs(W) ** 2 / rho


def local_average(fld, delta, rho):
    """Prediction: field average of the local Wiener relative MSE."""
    G2 = 4.0 * np.sin(0.5 * np.outer(delta, fld.k)) ** 2
    return np.mean(1.0 / (rho * G2 + 1.0), axis=0)


def summarise(fld, rel):
    """(per-null core rms, wide rms, peak, lost modes), edge sum, overall, %>0.5."""
    out = []
    for Km in nulls_in_band(D0, KMAX):
        c, w = fld.band(Km, CORE), fld.band(Km, WIDE)
        out.append((np.sqrt(rel[c].mean()), np.sqrt(rel[w].mean()),
                    np.sqrt(rel[w].max()), rel[w].sum() / 2))
    edge = rel[np.abs(fld.k) > 29.5].sum() / 2
    return out, edge, np.sqrt(rel.mean()), 100 * np.mean(np.sqrt(rel) > 0.5)


def row(lbl, fld, rel):
    nb, edge, tot, frac = summarise(fld, rel)
    s = f"  {lbl:30s}"
    for c, w, p, n in nb:
        s += f" {c:5.3f} {w:5.3f} {n:5.2f} |"
    print(s + f" {edge:5.2f} {tot:7.4f} {frac:6.2f}%")
    return tot


def header():
    print(f"  {'':30s} {'m=1':^17s} | {'m=2':^17s} | {'m=3':^17s} | "
          f"{'edge':>4s} {'overall':>7s} {'>0.5':>7s}")
    print(f"  {'case':30s}" + " core  wide  lost |" * 3 + "  lost    RMSE    band")


# ------------------------------------------------------------- sections 0, 1
def sec0_consistency(fld):
    print(RULE)
    print("0. CONSISTENCY: uniform shear, dense machinery vs closed form")
    print(RULE)
    print(f"  field {fld.n} px x {fld.dx} um = {fld.L:.1f} um, dk = {fld.dk:.4f} rad/um, "
          f"band {KMIN}-{KMAX:.2f} rad/um, {len(fld.k)} modes, rho = {RHO:g}")
    print(f"  in-band nulls of {D0} um: "
          + ", ".join(f"{v:.2f}" for v in nulls_in_band(D0, KMAX)) + " rad/um;"
          f" m=4 at {4 * 2 * np.pi / D0:.2f}, {4 * 2 * np.pi / D0 - KMAX:.2f} above the band edge")
    rel, _, _ = posterior_relmse([fld.B(np.full(fld.n, D0))], [RHO])
    ana = 1.0 / (RHO * G2_vec(fld.k, 0.0, (D0, 0.0)) + 1.0)
    print(f"  max |dense - 1/(rho G^2 + 1)| over the band = {np.max(np.abs(rel - ana)):.2e}")
    b, nz = naive_relmse(fld, fld.B(np.full(fld.n, D0)), D0, RHO)
    print(f"  max |naive(bias+noise) - same|              = {np.max(np.abs(b + nz - ana)):.2e}")
    print(f"  intrinsic null width (rho G^2 < 1): FWHM = 2/(D0 sqrt(rho)) = "
          f"{2 / (D0 * np.sqrt(RHO)):.3f} rad/um;")
    print(f"  Lorentzian area pi/(D0 sqrt(rho)) = {np.pi / (D0 * np.sqrt(RHO)):.3f} rad/um "
          f"= {np.pi / (D0 * np.sqrt(RHO)) / fld.dk:.2f} modes per null and sign of K")


def sec1_nullspace(fld):
    print()
    print(RULE)
    print("1. NULL SPACE OF A SMOOTH SHEAR FIELD: chirps exp(2 pi i m alpha(x))")
    print(RULE)
    print("  alpha' = (1 + Delta'/2)/Delta (second-order WKB; exact for the ramp).")
    print("  Rayleigh quotient |B c|^2/|c|^2 of a Hann-windowed test function,")
    print("  projected on the band. Lost at this SNR when it is below 1/rho = "
          f"{1 / RHO:g}.")
    print(f"  {'profile':8s} {'eps':>5s} {'m':>2s} {'chirp':>10s} {'plane wave':>11s} "
          f"{'in band':>8s}   (plane wave = exp(i K_m x), same window)")
    hann = 0.5 * (1 + np.cos(2 * np.pi * fld.x / fld.L))
    for name in ("ramp", "quad"):
        for eps in (0.0, 0.03, 0.10):
            if eps == 0.0 and name == "quad":
                continue
            delta = D0 * (1 + eps * fld.profile(name))
            dprime = np.gradient(delta, fld.dx)
            alpha = np.cumsum((1 + dprime / 2) / delta) * fld.dx
            alpha -= alpha[fld.n // 2]
            Bm = fld.B(delta)
            for m in (1, 3):
                q = []
                for ph in (2 * np.pi * m * alpha, 2 * np.pi * m / D0 * fld.x):
                    full = np.fft.fft(hann * np.exp(1j * ph)) / np.sqrt(fld.n)
                    c = fld.U.conj().T @ (hann * np.exp(1j * ph))
                    inband = np.sum(np.abs(c) ** 2) / np.sum(np.abs(full) ** 2)
                    q.append(np.sum(np.abs(Bm @ c) ** 2) / np.sum(np.abs(c) ** 2))
                print(f"  {name if eps else 'uniform':8s} {100 * eps:4.0f}% {m:2d} {q[0]:10.2e} "
                      f"{q[1]:11.2e} {100 * inband:7.2f}%")
    print("  => the chirp stays (numerically) null at every amplitude; the plane")
    print("     wave at K_m does not. The null has moved to another function, not gone.")


# ---------------------------------------------------------------- section 2
def sec2_calibrated(fld, store):
    print()
    print(RULE)
    print(f"2. CALIBRATED ESTIMATOR, EXACT  (L = {fld.L:.1f} um, rho = {RHO:g})")
    print(RULE)
    print("  core/wide = rms relative error over |K| within "
          f"+-{CORE}/+-{WIDE} rad/um of the null;")
    print("  lost = sum of relative MSE over the wide band, in modes (per sign of K);")
    print("  edge lost = same for |K| > 29.5 (skirt of the m=4 null, outside the band).")
    header()
    uni = fld.B(np.full(fld.n, D0))
    rel, _, _ = posterior_relmse([uni], [RHO])
    store["uniform"] = rel
    base = row(f"uniform {D0}", fld, rel)
    for name in ("ramp", "quad"):
        for eps in EPS:
            delta = D0 * (1 + eps * fld.profile(name))
            rel, F, C = posterior_relmse([fld.B(delta)], [RHO])
            store[(name, eps)] = rel
            store[(name, eps, "C")] = C
            row(f"{name} {100 * eps:.0f}%  [{delta.min():.3f},{delta.max():.3f}]", fld, rel)
    p2 = fld.B(np.full(fld.n, D2))
    rel, _, _ = posterior_relmse([uni, p2], [RHO / 2, RHO / 2])
    store["pair"] = rel
    pair = row(f"pair {{{D0},{D2}}}, rho/2 each", fld, rel)
    rel, _, _ = posterior_relmse([uni, p2], [RHO, RHO])
    row(f"pair {{{D0},{D2}}}, rho each", fld, rel)
    print(f"  overall RMSE: uniform {base:.4f}; pair at equal exposure {pair:.4f} "
          f"({base / pair:.2f}x better)")

    print()
    print("  Peak error per Fourier mode near each null: exact, the local-average")
    print("  prediction <1/(rho G^2(K, Delta(x)) + 1)>_x, and the diagonal-Fisher")
    print("  lower bound 1/(rho <G^2>_x + 1).  (rms relative error)")
    print(f"  {'case':12s}" + "".join(f" | m={m}: exact  local  bound" for m in (1, 2, 3))
          + " | max |exact-local|")
    for name in ("ramp", "quad"):
        for eps in EPS:
            delta = D0 * (1 + eps * fld.profile(name))
            rel = store[(name, eps)]
            loc = local_average(fld, delta, RHO)
            G2 = 4.0 * np.sin(0.5 * np.outer(delta, fld.k)) ** 2
            bound = 1.0 / (RHO * G2.mean(axis=0) + 1.0)
            store[(name, eps, "loc")] = loc
            s = f"  {name} {100 * eps:3.0f}%   "
            for Km in nulls_in_band(D0, KMAX):
                w = fld.band(Km, WIDE)
                i = np.flatnonzero(w)[np.argmax(rel[w])]
                s += f" |      {np.sqrt(rel[i]):5.3f}  {np.sqrt(loc[i]):5.3f}  {np.sqrt(bound[i]):5.3f}"
            print(s + f" |  {np.max(np.abs(np.sqrt(rel) - np.sqrt(loc))):.3f}")
    print()
    print("  Smear width of the m-th null on the whole field. Equivalent width =")
    print("  (sum of relative MSE) * dk / peak, to be compared with the intrinsic width")
    print("  of the uniform shear (below) and, once larger, with the prediction for a")
    print("  ramp, K_m * eps (the range of 2 pi m/Delta(x)).")
    ew0 = {}
    for Km in nulls_in_band(D0, KMAX):
        w = fld.band(Km, 4.0) & (fld.k > 0)
        ew0[Km] = store["uniform"][w].sum() * fld.dk / store["uniform"][w].max()
    print(f"  intrinsic equivalent width (uniform): {ew0[nulls_in_band(D0, KMAX)[0]]:.3f} rad/um")
    print(f"  {'case':12s}" + "".join(f" | m={m}: width  K_m*eps" for m in (1, 2, 3)))
    for eps in EPS:
        s = f"  ramp {100 * eps:3.0f}%   "
        for Km in nulls_in_band(D0, KMAX):
            w = fld.band(Km, 4.0) & (fld.k > 0)
            rel = store[("ramp", eps)]
            ew = rel[w].sum() * fld.dk / rel[w].max()
            s += f" |     {ew:6.3f}  {Km * eps:6.3f}"
        print(s)


def sec2b_field_and_snr():
    print()
    print("  Dependence on field size and SNR (calibrated, overall RMSE and lost")
    print("  modes summed over the three nulls, per sign of K and per 51.2 um of field):")
    print(f"  {'L (um)':>7s} {'rho':>7s} {'case':10s} {'RMSE':>8s} {'lost/51.2um':>12s} "
          f"{'peak m=3':>9s}")
    for n in (512, 1024):
        f = Field1D(n)
        for rho in (100.0, 1e4):
            for name, eps in (("uniform", 0.0), ("ramp", 0.10), ("quad", 0.10)):
                delta = D0 * (1 + eps * f.profile("ramp" if name == "uniform" else name))
                rel, _, _ = posterior_relmse([f.B(delta)], [rho])
                lost = sum(rel[f.band(Km, WIDE)].sum() / 2 for Km in nulls_in_band(D0, KMAX))
                pk = np.sqrt(rel[f.band(nulls_in_band(D0, KMAX)[2], WIDE)].max())
                print(f"  {f.L:7.1f} {rho:7.0f} {name + (f' {100*eps:.0f}%' if eps else ''):10s} "
                      f"{np.sqrt(rel.mean()):8.4f} {lost * 51.2 / f.L:12.2f} {pk:9.3f}")


# ---------------------------------------------------------------- section 3
def gabor_profile(fld, C, x0, sig_x, Ks):
    """Posterior/prior variance ratio of Gabor atoms centred (x0, K)."""
    out = np.empty(len(Ks))
    for i, K in enumerate(Ks):
        g = np.exp(-((fld.k - K) * sig_x) ** 2) * np.exp(-1j * (fld.k - K) * x0)
        out[i] = np.real(g.conj() @ C @ g) / np.real(g.conj() @ g)
    return out


def husimi_profile(fld, delta, x0, sig_x, Ks, rho):
    """Same, predicted by smoothing the local Wiener error 1/(rho G^2 + 1) with
    the atom's phase-space footprint (std sig_x in x, 1/(2 sig_x) in K)."""
    wx = np.exp(-0.5 * ((fld.x - x0) / sig_x) ** 2)
    wx /= wx.sum()
    sym = wx @ (1.0 / (rho * 4.0 * np.sin(0.5 * np.outer(delta, fld.k)) ** 2 + 1.0))
    sk = 1.0 / (2 * sig_x)
    out = np.empty(len(Ks))
    for i, K in enumerate(Ks):
        v = np.exp(-0.5 * ((fld.k - K) / sk) ** 2)
        out[i] = v @ sym / v.sum()
    return out


def sec3_local_fourier(fld, store):
    print()
    print(RULE)
    print("3. WINDOWED LOCAL-FOURIER MAP  (Gabor atoms, exact posterior)")
    print(RULE)
    m = 3
    Km = 2 * np.pi * m / D0
    Ks = np.arange(Km - 4.0, Km + 4.0, 0.02)
    f2s = 2 * np.sqrt(2 * np.log(2))
    print(f"  (a) where the m={m} null sits, window FWHM 6.0 um: location of the")
    print("      worst-recovered atom vs 2 pi m / Delta(x0)")
    print(f"  {'case':10s} {'x0 (um)':>8s} {'Delta(x0)':>10s} {'predicted':>10s} {'found':>8s}")
    for name in ("ramp", "quad"):
        delta = D0 * (1 + 0.10 * fld.profile(name))
        C = store[(name, 0.10, "C")]
        for x0 in (-0.25 * fld.L, 0.0, 0.25 * fld.L):
            d0 = np.interp(x0, fld.x, delta)
            pr = gabor_profile(fld, C, x0, 6.0 / f2s, Ks)
            print(f"  {name + ' 10%':10s} {x0:8.1f} {d0:10.4f} {2 * np.pi * m / d0:10.2f} "
                  f"{Ks[np.argmax(pr)]:8.2f}")
    print()
    print(f"  (b) width of the m={m} null seen through a window, ramp, x0 = 0.")
    print("      W = FWHM of the window intensity. EW = equivalent width of the")
    print("      variance-ratio bump (area/peak). Scales: intrinsic FWHM "
          f"{2 / (D0 * np.sqrt(RHO)):.2f},")
    print("      window resolution = FWHM of the atom in K, smear = K_m eps W / L")
    print("      = K_m |grad Delta| W / Delta0.")
    print(f"  {'eps':>5s} {'W (um)':>7s} {'resol.':>7s} {'smear':>7s} {'EW exact':>9s} "
          f"{'EW Husimi':>10s} {'peak exact':>11s} {'peak Husimi':>12s}")
    for eps in (0.0, 0.03, 0.10):
        delta = D0 * (1 + eps * fld.profile("ramp"))
        C = store[("ramp", eps, "C")] if eps else posterior_relmse([fld.B(delta)], [RHO])[2]
        for W in (1.5, 3.0, 6.0, 12.0, 24.0):
            sx = W / f2s
            pe = gabor_profile(fld, C, 0.0, sx, Ks)
            ph = husimi_profile(fld, delta, 0.0, sx, Ks, RHO)
            dK = Ks[1] - Ks[0]
            print(f"  {100 * eps:4.0f}% {W:7.1f} {f2s / (2 * sx):7.2f} {Km * eps * W / fld.L:7.2f} "
                  f"{pe.sum() * dK / pe.max():9.2f} {ph.sum() * dK / ph.max():10.2f} "
                  f"{pe.max():11.3f} {ph.max():12.3f}")
    # map for the figure
    delta = D0 * (1 + 0.10 * fld.profile("ramp"))
    x0s = np.linspace(-0.4 * fld.L, 0.4 * fld.L, 41)
    Kg = np.arange(Km - 3.0, Km + 3.0, 0.05)
    gmap = np.array([gabor_profile(fld, store[("ramp", 0.10, "C")], x0, 6.0 / f2s, Kg)
                     for x0 in x0s])
    return x0s, Kg, gmap, 2 * np.pi * m / np.interp(x0s, fld.x, delta)


# ---------------------------------------------------------------- section 4
def sec4_naive(fld, store):
    print()
    print(RULE)
    print(f"4. NAIVE ESTIMATOR (uniform-shear Wiener filter), EXACT  (L = {fld.L:.1f} um)")
    print(RULE)
    print("  'nominal' assumes Delta0 = 0.8 um; 'mean' assumes the field-mean shear")
    print("  (identical for the ramp). bias = model-error part alone (noise-free data,")
    print("  same filter). hi-K = rms relative error over |K| > 25 rad/um.")
    print(f"  {'case':22s} {'rho':>6s} {'calibrated':>11s} {'naive':>8s} {'bias only':>10s} "
          f"{'hi-K cal.':>10s} {'hi-K naive':>11s} {'rms dDelta':>11s}")
    hi = np.abs(fld.k) > 25.0
    for rho in (RHO, 1e4):
        for name in ("ramp", "quad"):
            for eps in EPS:
                delta = D0 * (1 + eps * fld.profile(name))
                Bm = fld.B(delta)
                cal = store[(name, eps)] if rho == RHO else posterior_relmse([Bm], [rho])[0]
                for tag, da in (("nominal", D0), ("mean", delta.mean())):
                    if name == "ramp" and tag == "mean":
                        continue
                    b, nz = naive_relmse(fld, Bm, da, rho)
                    store[(name, eps, tag, rho)] = b + nz
                    print(f"  {name + f' {100*eps:.0f}% ' + tag:22s} {rho:6.0f} "
                          f"{np.sqrt(cal.mean()):11.4f} {np.sqrt((b + nz).mean()):8.4f} "
                          f"{np.sqrt(b.mean()):10.4f} {np.sqrt(cal[hi].mean()):10.4f} "
                          f"{np.sqrt((b + nz)[hi].mean()):11.4f} "
                          f"{1e3 * np.sqrt(np.mean((delta - da) ** 2)):8.1f} nm")
    print()
    print("  Null bands, naive (nominal) vs calibrated, rho = 100:")
    header()
    row(f"uniform {D0}", fld, store["uniform"])
    for name in ("ramp", "quad"):
        for eps in EPS:
            row(f"{name} {100 * eps:.0f}% calibrated", fld, store[(name, eps)])
            row(f"{name} {100 * eps:.0f}% naive", fld, store[(name, eps, "nominal", RHO)])
    print()
    print("  Small-error estimate of the naive model error: per mode, relative MSE")
    print("  ~ (K s)^2 / (G0^2 + 1/rho)^2 * G0^2, s = rms shear error; it passes the")
    print("  noise term where K s ~ 1/sqrt(rho):")
    for eps in EPS:
        s = D0 * eps / np.sqrt(12)
        print(f"    ramp {100 * eps:3.0f}%: s = {1e3 * s:5.1f} nm, K_max s = {KMAX * s:.3f} "
              f"(1/sqrt(rho) = {1 / np.sqrt(RHO):.3f}); model error exceeds 20 dB noise "
              f"above K = {min(1 / (np.sqrt(RHO) * s), 999):.1f} rad/um")


# ---------------------------------------------------------------- section 5
def sec5_montecarlo(fld, store, rng, trials=200):
    print()
    print(RULE)
    print(f"5. MONTE CARLO CHECK  ({trials} real objects, real noise, seed {SEED})")
    print(RULE)
    print(f"  {'case':14s} {'estimator':11s} {'exact RMSE':>11s} {'MC RMSE':>9s} "
          f"{'exact m=3 wide':>15s} {'MC m=3 wide':>12s}")
    nk = len(fld.k)
    pos = fld.k > 0
    w3 = fld.band(nulls_in_band(D0, KMAX)[2], WIDE)
    for name, eps in (("ramp", 0.03), ("ramp", 0.10), ("quad", 0.10)):
        delta = D0 * (1 + eps * fld.profile(name))
        Bm = fld.B(delta)
        F = Bm.conj().T @ Bm
        solve = np.linalg.inv(RHO * F + np.eye(nk)) @ (RHO * Bm.conj().T)
        G0 = np.exp(1j * fld.k * D0) - 1.0
        W = np.conj(G0) / (np.abs(G0) ** 2 + 1.0 / RHO)
        e_cal, e_nai = np.zeros(nk), np.zeros(nk)
        for _ in range(trials):
            c = np.zeros(nk, complex)
            half = (rng.standard_normal(pos.sum()) + 1j * rng.standard_normal(pos.sum())) / np.sqrt(2)
            c[pos] = half
            c[~pos] = np.conj(half[::-1])          # Hermitian: phi real
            psi = np.real(Bm @ c) + rng.standard_normal(fld.n) / np.sqrt(RHO)
            e_cal += np.abs(solve @ psi - c) ** 2
            e_nai += np.abs(W * (fld.U.conj().T @ psi) - c) ** 2
        for tag, mc, ex in (("calibrated", e_cal / trials, store[(name, eps)]),
                            ("naive", e_nai / trials, store[(name, eps, "nominal", RHO)])):
            print(f"  {name + f' {100*eps:.0f}%':14s} {tag:11s} {np.sqrt(ex.mean()):11.4f} "
                  f"{np.sqrt(mc.mean()):9.4f} {np.sqrt(ex[w3].mean()):15.4f} "
                  f"{np.sqrt(mc[w3].mean()):12.4f}")


# ---------------------------------------------------------------- section 6
class Shear2D:
    """psi(x, y) = phi(x + Delta(x, y), y) - phi(x, y) on an n x n periodic grid.
    Exact for band-limited phi: row-wise trigonometric interpolation."""

    def __init__(self, n, delta, dx=DX):
        self.n = n
        kx = 2 * np.pi * np.fft.fftfreq(n, dx)
        # The phase must use the sample positions numpy's FFT assumes, i*dx,
        # not the centred coordinates Delta(x, y) is defined on: with centred
        # x every frequency picks up (-1)^q and the operator shifts the image
        # by half the field as well as by Delta. (Caught in review: the
        # calibrated estimator, inverting the same operator, hid it.)
        xs = np.arange(n) * dx
        # T[y, x, kx] = exp(i kx (x + Delta(x, y))) / n ; Nyquist column dropped
        ph = (xs[:, None] + delta)[..., None] * kx[None, None, :]     # (x, y, kx)
        self.T = (np.exp(1j * ph) / n).transpose(1, 0, 2).copy()
        self.T[:, :, n // 2] = 0.0

    def fwd(self, phi):
        P = np.fft.fft(phi, axis=0)                                    # (kx, y)
        return np.real(np.matmul(self.T, P.T[:, :, None])[:, :, 0].T) - phi

    def adj(self, r):
        Y = np.matmul(self.T.conj().transpose(0, 2, 1), r.T[:, :, None])[:, :, 0].T   # (kx, y)
        return np.real(np.fft.ifft(Y, axis=0)) * self.n - r


def cg(apply, b, tol=1e-7, maxiter=600):
    x = np.zeros_like(b)
    r = b.copy()
    p = r.copy()
    rs = np.vdot(r, r).real
    b2 = rs
    for it in range(1, maxiter + 1):
        Ap = apply(p)
        a = rs / np.vdot(p, Ap).real
        x += a * p
        r -= a * Ap
        rn = np.vdot(r, r).real
        if rn < tol ** 2 * b2:
            break
        p = r + (rn / rs) * p
        rs = rn
    return x, it


def sec6_images(rng, n=256):
    print()
    print(RULE)
    print(f"6. 2D IMAGE DEMO  ({n}^2 px, dx = {DX} um, field {n * DX:.1f} um, band-pass "
          f"{KMIN}-{KMAX:.1f} rad/um, 20 dB)")
    print(RULE)
    x = (np.arange(n) - n // 2) * DX
    X, Y = np.meshgrid(x, x, indexing="ij")
    k = 2 * np.pi * np.fft.fftfreq(n, DX)
    KX, KY = np.meshgrid(k, k, indexing="ij")
    KR = np.hypot(KX, KY)
    band = (KR <= KMAX) & (KR >= KMIN)
    proj = lambda a: np.real(np.fft.ifft2(np.fft.fft2(a) * band))
    O = proj(shepp_logan_2d(X, Y, scale=n * DX / 2.56))
    Ok = np.fft.fft2(O)
    Ssig = np.mean(np.abs(Ok[band]) ** 2) / n ** 2        # per-mode power, unitary DFT
    sig = np.sqrt(Ssig / RHO)                             # per-pixel noise std
    lam = 1.0 / RHO
    noise = [sig * rng.standard_normal((n, n)) for _ in range(2)]
    masks = [(np.abs(np.abs(KX) - Km) < WIDE) & band for Km in nulls_in_band(D0, KMAX)]

    def score(rec):
        Ek = np.fft.fft2(rec - O)
        return (np.linalg.norm(rec - O) / np.linalg.norm(O),
                [np.linalg.norm(Ek[m]) / np.linalg.norm(Ok[m]) for m in masks])

    def wiener_uniform(psis, ds, lam_):
        num, den = 0.0, 0.0
        for psi, d in zip(psis, ds):
            G = np.exp(1j * KX * d) - 1.0
            num = num + np.conj(G) * np.fft.fft2(psi)
            den = den + np.abs(G) ** 2
        return proj(np.real(np.fft.ifft2(num / (den + lam_))))

    def shear_uniform(d):
        return np.real(np.fft.ifft2(Ok * (np.exp(1j * KX * d) - 1.0)))

    print("  null bands: relative error of the part of the image spectrum with")
    print(f"  |K_x| within +-{WIDE} rad/um of each null. Same noise realisation in every row.")
    print(f"  {'case':28s} {'estimator':11s} {'rel RMSE':>9s} {'m=1':>6s} {'m=2':>6s} {'m=3':>6s} "
          f"{'CG its':>7s}")
    recs = {}

    def show(lbl, est, rec, its=""):
        e, nb = score(rec)
        print(f"  {lbl:28s} {est:11s} {e:9.4f} " + " ".join(f"{v:6.3f}" for v in nb)
              + f" {its:>7}")
        return e

    op0 = Shear2D(n, np.full((n, n), D0))
    print(f"  operator check: uniform Delta through Shear2D vs Fourier shift, max |diff| / max "
          f"= {np.abs(op0.fwd(O) - shear_uniform(D0)).max() / np.abs(shear_uniform(D0)).max():.1e}")
    del op0
    rec = wiener_uniform([shear_uniform(D0) + noise[0]], [D0], lam)
    base = show(f"uniform {D0}", "Wiener", rec)
    recs["uniform 0.8"] = rec
    rec = wiener_uniform([shear_uniform(D0) + np.sqrt(2) * noise[0],
                          shear_uniform(D2) + np.sqrt(2) * noise[1]], [D0, D2], 2 * lam)
    pair = show(f"pair {{{D0},{D2}}}, rho/2 each", "Wiener", rec)
    recs["pair {0.8, 0.55}"] = rec
    rec = wiener_uniform([shear_uniform(D0) + noise[0], shear_uniform(D2) + noise[1]],
                         [D0, D2], lam)
    show(f"pair {{{D0},{D2}}}, rho each", "Wiener", rec)

    u, v = X / (n * DX), Y / (n * DX)
    profiles = {"ramp along x": u, "ramp along y": v,
                "radial bowl": (u ** 2 + v ** 2) / 0.5}
    t_cg, worst_adj = 0.0, 0.0
    for name, f in profiles.items():
        for eps in (0.03, 0.10):
            delta = D0 * (1 + eps * f)
            op = Shear2D(n, delta)
            a, b = proj(rng.standard_normal((n, n))), rng.standard_normal((n, n))
            worst_adj = max(worst_adj, abs(np.vdot(op.fwd(a), b) - np.vdot(a, op.adj(b)))
                            / abs(np.vdot(op.fwd(a), b)))
            psi = op.fwd(O) + noise[0]
            t0 = time.time()
            rec, its = cg(lambda p: proj(op.adj(op.fwd(p))) + lam * p, proj(op.adj(psi)))
            t_cg += time.time() - t0
            lbl = f"{name} {100 * eps:.0f}%"
            show(lbl, "calibrated", rec, str(its))
            show(lbl, "naive", wiener_uniform([psi], [D0], lam))
            if name == "radial bowl":
                show(lbl, "naive mean", wiener_uniform([psi], [delta.mean()], lam))
            if eps == 0.10 and name != "ramp along y":
                recs[f"{lbl}, calibrated"] = rec
                recs[f"{lbl}, naive"] = wiener_uniform([psi], [D0], lam)
            del op
    print(f"  adjoint test, worst relative mismatch <A a, b> vs <a, A^T b>: {worst_adj:.1e}")
    print(f"  uniform / pair(rho/2 each) RMSE ratio: {base / pair:.2f}x   (CG total {t_cg:.0f} s)")
    return O, recs


# ---------------------------------------------------------------- section 7
def sec7_translation(fld, store):
    print()
    print(RULE)
    print("7. WHERE DIVERSITY DOES COME FROM: the object seen at two field positions")
    print(RULE)
    print("  Two exposures at rho/2 each; between them the sample is translated by s")
    print("  across the ramp field (periodic object, so every pixel is seen twice,")
    print("  with shears differing by eps*D0*s/L away from the wrap). Exact posterior.")
    header()
    row(f"uniform {D0}, single", fld, store["uniform"])
    for eps in EPS:
        f = fld.profile("ramp")
        d1 = D0 * (1 + eps * f)
        for frac in (0.25, 0.5):
            sh = int(round(frac * fld.n))
            d2 = np.roll(d1, sh)
            rel, _, _ = posterior_relmse([fld.B(d1), fld.B(d2)], [RHO / 2, RHO / 2])
            row(f"ramp {100 * eps:.0f}%, s = {frac:.2f} L", fld, rel)
            store[("trans", eps, frac)] = rel
    row(f"pair {{{D0},{D2}}}, rho/2 each", fld, store["pair"])
    print("  (shear difference between the two looks at s = L/2: "
          + ", ".join(f"{100*e:.0f}% -> {1e3 * D0 * e / 2:.0f} nm" for e in EPS) + ")")


# -------------------------------------------------------------------- figure
C_ = ["#257de0", "#f47230", "#00b894", "#f2a600", "#6a4ce2"]
INK, MUTED = "#191919", "#666666"


def style(ax, title, xl, yl):
    ax.set_title(title, fontsize=10, color=INK, loc="left")
    ax.set_xlabel(xl, color=MUTED)
    ax.set_ylabel(yl, color=MUTED)
    ax.grid(alpha=0.25, lw=0.6)
    ax.tick_params(colors=MUTED, labelsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def figures(fld, store, gab, O, recs):
    pos = fld.k > 0
    K = fld.k[pos]
    fig, ax = plt.subplots(2, 2, figsize=(13.5, 8.6))
    a = ax[0, 0]
    a.plot(K, np.sqrt(store["uniform"][pos]), color=C_[0], lw=1.6, label="uniform 0.8")
    for eps, c in zip((0.03, 0.10), C_[1:3]):
        a.plot(K, np.sqrt(store[("ramp", eps)][pos]), color=c, lw=1.6,
               label=f"ramp {100 * eps:.0f} %, calibrated")
    a.plot(K, np.sqrt(store["pair"][pos]), color=C_[4], lw=1.6, label="pair {0.8, 0.55}, equal exposure")
    style(a, "Calibrated: the null spreads, its area stays", "K (rad/um)", "relative error per Fourier mode")
    a.set_ylim(0, 1.05)
    a.legend(fontsize=8, frameon=False)

    a = ax[0, 1]
    Km = 3 * 2 * np.pi / D0
    z = pos & (np.abs(fld.k - Km) < 4)
    for (name, eps), c in zip((("ramp", 0.03), ("ramp", 0.10), ("quad", 0.10)), C_[1:4]):
        a.plot(fld.k[z], np.sqrt(store[(name, eps)][z]), color=c, lw=1.8,
               label=f"{name} {100 * eps:.0f} %: exact")
        a.plot(fld.k[z], np.sqrt(store[(name, eps, "loc")][z]), color=c, lw=1.0, ls="--")
    a.plot([], [], color=MUTED, lw=1.0, ls="--", label="field average of the local error")
    style(a, "Null m = 3: exact vs local-average prediction", "K (rad/um)", "relative error per Fourier mode")
    a.set_ylim(0, 1.05)
    a.legend(fontsize=8, frameon=False)

    a = ax[1, 0]
    for eps, c in zip((0.01, 0.03, 0.10), (C_[0], C_[1], C_[2])):
        a.plot(K, np.sqrt(store[("ramp", eps, "nominal", RHO)][pos]), color=c, lw=1.4,
               label=f"ramp {100 * eps:.0f} %, naive")
    a.plot(K, np.sqrt(store[("ramp", 0.10)][pos]), color=INK, lw=1.0, ls="--",
           label="ramp 10 %, calibrated")
    style(a, "Naive filter: model error grows with K", "K (rad/um)", "relative error per Fourier mode")
    a.set_ylim(0, 1.6)
    a.legend(fontsize=8, frameon=False)

    a = ax[1, 1]
    x0s, Kg, gmap, Kpred = gab
    im = a.imshow(gmap.T, origin="lower", aspect="auto", cmap="Blues", vmin=0, vmax=1,
                  extent=[x0s[0], x0s[-1], Kg[0], Kg[-1]])
    a.plot(x0s, Kpred, color=C_[1], lw=1.5, label="2 pi m / Delta(x0)")
    style(a, "Ramp 10 %: lost Gabor atoms (6 um window) follow the local null",
          "window centre x0 (um)", "K (rad/um)")
    a.grid(False)
    a.legend(fontsize=8, frameon=False, loc="upper right")
    cb = fig.colorbar(im, ax=a)
    cb.set_label("posterior / prior variance", color=MUTED, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path("varying_shear.png"), dpi=115)

    fig2, axes = plt.subplots(2, 4, figsize=(15, 7.6))
    items = [("ground truth (band-passed)", O)] + list(recs.items())
    for a, (lbl, r) in zip(axes.ravel(), items):
        a.imshow(r.T, cmap="gray", vmin=O.min(), vmax=O.max(), origin="lower")
        a.set_title(lbl, fontsize=9, color=INK)
        a.axis("off")
    for a in axes.ravel()[len(items):]:
        a.axis("off")
    fig2.tight_layout()
    fig2.savefig(out_path("varying_shear_images.png"), dpi=115)
    print()
    print(f"  wrote {out_path('varying_shear.png')}")
    print(f"  wrote {out_path('varying_shear_images.png')}")


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    fld = Field1D(512)
    store = {}
    sec0_consistency(fld)
    sec1_nullspace(fld)
    sec2_calibrated(fld, store)
    sec2b_field_and_snr()
    gab = sec3_local_fourier(fld, store)
    sec4_naive(fld, store)
    sec5_montecarlo(fld, store, rng)
    O, recs = sec6_images(rng)
    sec7_translation(fld, store)
    figures(fld, store, gab, O, recs)
    print(f"  total run time {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
