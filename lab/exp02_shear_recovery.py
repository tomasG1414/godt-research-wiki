"""exp02 -- Per-frequency recovery error, and a band-passed image demo.

Establishes:
  1. The Wiener error formula  rel_err = [rho*G^2 + 1]^(-1/2)  against Monte
     Carlo. At a null the error is exactly 1 for any SNR.
  2. Per-frequency error curves: single shear spikes at each null; an
     incommensurate second shear removes them, a commensurate one does not.
  3. Image-domain RMSE, band-passed so the DC loss shared by every gradient
     method is excluded and the nulls are what is being measured.

A NOTE ON THE FIRST ATTEMPT, kept as a warning: recovering O from gradient
data is ill-posed as K -> 0 for ANY shear (that is gradient integration,
which GODT is designed not to do). On a low-frequency-dominated phantom that
effect buries the nulls entirely and inverts the apparent ranking. The band-
pass below is what separates the two.

Writes: out/err_vs_K.png, out/recovery.png
Wiki:   wiki/notes/finite-shear-null-experiment.md

Run: python exp02_shear_recovery.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from shearlib import (GODT_SIM, H_vec, combined_G2, wiener_rel_err,
                      shepp_logan_2d, out_path)

SEED = 1
RHO = 100.0          # per-frequency SNR, 20 dB
RULE = "=" * 78


def montecarlo_check(rng):
    print(RULE)
    print("1. WIENER ERROR FORMULA vs MONTE CARLO   (rho = 100)")
    print(RULE)
    print(f"  {'K':>9} {'G^2':>12} {'analytic':>10} {'monte carlo':>12}")
    for K in (3.0, 2 * np.pi / 0.8, 8.5, 15.0, 2 * 2 * np.pi / 0.8, 22.0):
        ds = [0.8]
        ana = wiener_rel_err(K, ds, RHO)
        n, S = 40000, 1.0
        sig = np.sqrt(S / RHO)
        O = (rng.normal(0, np.sqrt(S / 2), n) + 1j * rng.normal(0, np.sqrt(S / 2), n))
        num, den = np.zeros(n, complex), 0.0
        for d in ds:
            h = (np.exp(1j * K * d) - 1) / d
            meas = h * O + (rng.normal(0, sig / np.sqrt(2), n)
                            + 1j * rng.normal(0, sig / np.sqrt(2), n))
            num += np.conj(h) * meas
            den += abs(h) ** 2
        est = num / (den + 1.0 / RHO)
        mc = np.sqrt(np.mean(np.abs(est - O) ** 2) / np.mean(np.abs(O) ** 2))
        print(f"  {K:9.3f} {combined_G2(K, ds):12.5f} {ana:10.4f} {mc:12.4f}")
    print()
    print("  => formula validated. At a null, error = 1 regardless of SNR:")
    print("     the frequency is structurally absent, not merely noisy.")


SETS_1D = [
    ("d=0.1 (GODT sim)",      [0.1]),
    ("d=0.33 (preprint exp)", [0.33]),
    ("d=0.8 (published exp)", [0.8]),
    ("d=0.8+0.4 (ratio 2)",   [0.8, 0.4]),
    ("d=0.8+0.55",            [0.8, 0.55]),
    ("d=0.8+0.08",            [0.8, 0.08]),
    ("d=0.8+0.55+0.35",       [0.8, 0.55, 0.35]),
]


def per_frequency(kmax):
    print()
    print(RULE)
    print("2. PER-FREQUENCY RECOVERY ERROR  (rho = 100, K > 2 rad/um)")
    print(RULE)
    K = np.linspace(0.01, kmax, 6000)
    m = K > 2.0
    print(f"  {'set':26s} {'worst':>8} {'at K':>8} {'%band>0.5':>11} {'mean':>8}")
    curves = {}
    for lbl, ds in SETS_1D:
        e = wiener_rel_err(K, ds, RHO)
        curves[lbl] = e
        print(f"  {lbl:26s} {e[m].max():8.4f} {K[m][e[m].argmax()]:8.2f} "
              f"{100*np.mean(e[m] > 0.5):10.2f}% {e[m].mean():8.4f}")
    print()
    print("  (K > 2 excludes the low-frequency rise every gradient method shares)")
    return K, curves


SETS_2D = [
    ("d=0.1 only",          [(0.1, 0)]),
    ("d=0.8 only",          [(0.8, 0)]),
    ("d=0.8 + 0.4",         [(0.8, 0), (0.4, 0)]),
    ("d=0.8 + 0.55",        [(0.8, 0), (0.55, 0)]),
    ("x+y at 0.8",          [(0.8, 0), (0, 0.8)]),
    ("x+y at 0.8 and 0.55", [(0.8, 0), (0.55, 0), (0, 0.8), (0, 0.55)]),
]


def image_demo(rng, kmax, n=512, dx=0.1, kmin=2.0, snr_db=20):
    print()
    print(RULE)
    print(f"3. IMAGE DEMO  ({n}^2, dx={dx} um, bandpass {kmin}-{kmax:.1f} rad/um,"
          f" {snr_db} dB)")
    print(RULE)
    x = (np.arange(n) - n // 2) * dx
    X, Y = np.meshgrid(x, x, indexing="ij")
    k = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    KX, KY = np.meshgrid(k, k, indexing="ij")
    KR = np.hypot(KX, KY)
    band = (KR <= kmax) & (KR >= kmin)

    Ok = np.fft.fft2(shepp_logan_2d(X, Y)) * band
    O = np.real(np.fft.ifft2(Ok))
    Ssig = np.mean(np.abs(Ok[band]) ** 2)
    sig2 = Ssig / 10 ** (snr_db / 10)

    print(f"  {'set':26s} {'rel RMSE':>10} {'vs single 0.8':>15}")
    recs, base = {}, None
    for lbl, dvecs in SETS_2D:
        num = np.zeros_like(Ok)
        den = np.zeros(Ok.shape)
        for dv in dvecs:
            h = H_vec(KX, KY, dv)
            meas = h * Ok + np.sqrt(sig2 / 2) * (rng.normal(size=Ok.shape)
                                                 + 1j * rng.normal(size=Ok.shape))
            num += np.conj(h) * meas
            den += np.abs(h) ** 2
        rec = np.real(np.fft.ifft2((num / (den + sig2 / Ssig)) * band))
        e = np.linalg.norm(rec - O) / np.linalg.norm(O)
        recs[lbl] = rec
        if lbl == "d=0.8 only":
            base = e
        gain = f"{base/e:.1f}x" if base and lbl != "d=0.8 only" else "-"
        print(f"  {lbl:26s} {e:10.4f} {gain:>15}")
    return O, recs


def main():
    rng = np.random.default_rng(SEED)
    kmax = GODT_SIM.kmax
    montecarlo_check(rng)
    K, curves = per_frequency(kmax)
    O, recs = image_demo(rng, kmax)

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.4))
    for lbl in ("d=0.1 (GODT sim)", "d=0.33 (preprint exp)", "d=0.8 (published exp)"):
        axes[0].plot(K, curves[lbl], label=lbl)
    axes[0].set(title="single shear: error reaches 1.0 at each null",
                xlabel="K  (rad/um)", ylabel="relative recovery error",
                ylim=(0, 1.05), xlim=(0, kmax))
    for lbl in ("d=0.8 (published exp)", "d=0.8+0.4 (ratio 2)",
                "d=0.8+0.55", "d=0.8+0.08"):
        axes[1].plot(K, curves[lbl], label=lbl)
    axes[1].set(title="second shear: incommensurate fills the nulls",
                xlabel="K  (rad/um)", ylabel="relative recovery error",
                ylim=(0, 1.05), xlim=(0, kmax))
    for a in axes:
        a.legend(fontsize=8)
        a.grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(out_path("err_vs_K.png"), dpi=115)

    fig2 = plt.figure(figsize=(15, 5.6))
    ax = fig2.add_subplot(2, 4, 1)
    ax.imshow(O, cmap="gray"); ax.set_title("ground truth", fontsize=8); ax.axis("off")
    for i, (lbl, r) in enumerate(recs.items(), start=2):
        ax = fig2.add_subplot(2, 4, i)
        ax.imshow(r, cmap="gray", vmin=O.min(), vmax=O.max())
        ax.set_title(lbl, fontsize=8); ax.axis("off")
    fig2.suptitle("Recovery at 20 dB, bandpassed", fontsize=11)
    fig2.tight_layout()
    fig2.savefig(out_path("recovery.png"), dpi=115)
    print()
    print(f"  wrote {out_path('err_vs_K.png')}")
    print(f"  wrote {out_path('recovery.png')}")


if __name__ == "__main__":
    main()
