"""exp01 -- Finite-shear transfer nulls: structure and conditioning.

Establishes:
  1. |H_d(K)| = 2|sin(K d/2)|/d exactly; H -> iK as d -> 0.
  2. At the published experimental shear (0.8 um) three nulls fall inside the
     measurable band; at the 0.1 um simulation shear, none.
  3. Conditioning G = sqrt(sum|H|^2): commensurate shear ratios fail.
     (H carries a 1/d factor, which favours small partner shears; for the
     physical ranking of designs use exp04, which weights by |e^{iKd} - 1|^2.)
  4. Magnitudes fill the m!=0 comb; directions address the K.d = 0 blind line.
     No in-plane shear set recovers K_x = K_y = 0.

Writes: out/transfer.png
Wiki:   wiki/notes/finite-shear-null-experiment.md

Run: python exp01_shear_nulls.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from shearlib import (GODT_SIM, GODT_EXP, H, absH, combined_G, nulls_in_band, out_path)

SHEARS = [0.1, 0.33, 0.8]
RULE = "=" * 78


def main():
    kmax = GODT_SIM.kmax

    # ---------------------------------------------------- 1. transfer factor
    print(RULE)
    print("1. TRANSFER FACTOR")
    print(RULE)
    K = np.linspace(-40, 40, 20001)
    for d in SHEARS:
        err = np.max(np.abs(np.abs(H(K, d)) - absH(K, d)))
        print(f"  d={d:5.2f} um   max | |H| - 2|sin(Kd/2)|/d | = {err:.3e}")
    for d in (1e-2, 1e-3, 1e-4):
        Kt = np.array([0.5, 5.0, 20.0])
        rel = np.max(np.abs(H(Kt, d) - 1j * Kt) / np.abs(1j * Kt))
        print(f"  d={d:.0e} um    max rel. deviation of H from iK = {rel:.3e}")
    print()
    print("  => |H_d(K)| = 2|sin(K d/2)|/d.  Zeros at K = 2*pi*m/d.")
    print("     m == 0 : DC loss, inherent to any gradient method, d-independent.")
    print("     m != 0 : finite-shear comb, moves with d, fillable by multi-shear.")

    # -------------------------------------------------- 2. nulls vs the band
    print()
    print(RULE)
    print("2. DO THE NULLS FALL INSIDE THE MEASURABLE BAND?")
    print(RULE)
    print(f"  {GODT_SIM.name}: k_det={GODT_SIM.k_det:.2f}  k_ill={GODT_SIM.k_ill:.2f}"
          f"  -> band |K_x| <= {kmax:.2f} rad/um")
    print(f"    sampling Nyquist at dx={GODT_SIM.dx} um: "
          f"{GODT_SIM.k_nyquist:.2f} rad/um (critically sampled)")
    print(f"  {GODT_EXP.name}: band |K_x| <= {GODT_EXP.kmax:.2f} rad/um")
    print()
    print(f"  {'shear':>7} {'first null':>12} {'in sim band':>13} {'positions':>30}")
    for d in SHEARS:
        ns = nulls_in_band(d, kmax)
        pos = ", ".join(f"{v:.2f}" for v in ns) if ns else "none"
        print(f"  {d:6.2f}  {2*np.pi/d:9.2f}    {len(ns):9d}    {pos:>30}")
    print()
    print("  => The 0.1 um simulation shear is null-free in band;")
    print("     the 0.8 um experimental shear is not.")

    # ------------------------------------------------------ 3. conditioning
    print()
    print(RULE)
    print("3. CONDITIONING  G(K) = sqrt(sum_n |H_n(K)|^2)")
    print(RULE)
    print("  Noise amplification of the least-squares inversion is 1/G.")
    print("  G = 0 means irrecoverable at any SNR.  K > 0.5 only, so the")
    print("  shared low-frequency rise does not dominate.")
    print()
    Kg = np.linspace(-kmax, kmax, 40001)
    m = np.abs(Kg) > 0.5

    def report(label, shears):
        g = combined_G(Kg, shears)[m]
        i = g.argmin()
        print(f"  {label:36s} min G = {g[i]:8.4f} at K = {Kg[m][i]:7.2f}"
              f"   median = {np.median(g):7.3f}")
        return g[i]

    print("  single:")
    for d in SHEARS:
        report(f"d = {d}", [d])
    print()
    print("  pairs:")
    for a, b in [(0.8, 0.4), (0.8, 0.55), (0.8, 0.08), (0.1, 0.33)]:
        ratio = a / b
        tag = " (ratio 2 - nulls coincide)" if abs(ratio - 2) < 1e-9 else ""
        report(f"d = {a} + {b}{tag}", [a, b])
    print()
    print("  triple:")
    report("d = 0.8 + 0.55 + 0.35", [0.8, 0.55, 0.35])

    # ------------------------------- 4. best partner for a fixed large shear
    print()
    print(RULE)
    print("4. IF HARDWARE FIXES d1 = 0.8 um, WHAT SECOND SHEAR IS BEST?")
    print("   (1/d-normalised H favours small shears; exp04 gives the physical ranking)")
    print(RULE)
    cands = np.arange(0.05, 0.81, 0.005)
    scored = sorted(((combined_G(Kg[m], [0.8, d2]).min(), d2) for d2 in cands),
                    reverse=True)
    print("  best 6:")
    for v, d2 in scored[:6]:
        print(f"    d2 = {d2:5.3f} -> min G = {v:8.4f}")
    print("  worst 4 (avoid, nulls coincide):")
    for v, d2 in scored[-4:]:
        print(f"    d2 = {d2:5.3f} -> min G = {v:8.4f}")

    # ------------------------------------------- 5. directions vs magnitudes
    print()
    print(RULE)
    print("5. DIRECTIONS vs MAGNITUDES")
    print(RULE)
    kx = np.linspace(-kmax, kmax, 601)
    KX, KY = np.meshgrid(kx, kx, indexing="ij")

    def G2_vec(sets):
        tot = np.zeros_like(KX)
        for dx_, dy_ in sets:
            mag = np.hypot(dx_, dy_)
            tot += (2 * np.abs(np.sin((KX * dx_ + KY * dy_) / 2.0)) / mag) ** 2
        return np.sqrt(tot)

    R = np.hypot(KX, KY) > 0.5
    for label, sets in [
        ("x only, 0.8", [(0.8, 0)]),
        ("x only, 0.8 + 0.55", [(0.8, 0), (0.55, 0)]),
        ("x + y, 0.8", [(0.8, 0), (0, 0.8)]),
        ("x + y, 0.8 and 0.55 each", [(0.8, 0), (0.55, 0), (0, 0.8), (0, 0.55)]),
        ("3 dirs 0/60/120 deg, 0.8", [(0.8 * np.cos(a), 0.8 * np.sin(a))
                                      for a in np.deg2rad([0, 60, 120])]),
    ]:
        g = G2_vec(sets)
        print(f"  {label:28s} min G = {g[R].min():7.4f}   "
              f"band with G<0.5: {100*np.mean(g[R] < 0.5):5.2f}%")
    print()
    print("  => Magnitudes fill the comb; directions address the K.d=0 blind")
    print("     line. No in-plane set recovers K_x=K_y=0 -- where the missing")
    print("     cone already sits. Confirmed in 3D by exp03.")

    # -------------------------------------------------------------- figure
    Kp = np.linspace(0, kmax, 4000)
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.2))
    ax = axes[0]
    for d, c in zip(SHEARS, ["C0", "C1", "C2"]):
        ax.plot(Kp, absH(Kp, d), c, label=f"d = {d} um")
        for nu in nulls_in_band(d, kmax):
            ax.axvline(nu, color=c, ls=":", alpha=.5)
    ax.plot(Kp, Kp, "k--", lw=1, alpha=.6, label="|iK|  (d -> 0)")
    ax.set(xlabel="K  (rad/um)", ylabel="|H|", xlim=(0, kmax),
           title="single shear; dotted = in-band nulls")
    ax.legend(fontsize=8)

    ax = axes[1]
    for lbl, ds in [("0.8 only", [0.8]), ("0.8 + 0.4", [0.8, 0.4]),
                    ("0.8 + 0.55", [0.8, 0.55]), ("0.8 + 0.08", [0.8, 0.08])]:
        ax.semilogy(Kp, np.maximum(combined_G(Kp, ds), 1e-6), label=lbl)
    ax.set(xlabel="K  (rad/um)", ylabel="G = sqrt(sum |H|^2)", xlim=(0, kmax),
           title="combined conditioning (higher is better)")
    ax.legend(fontsize=8)
    ax.grid(alpha=.3)
    fig.tight_layout()
    p = out_path("transfer.png")
    fig.savefig(p, dpi=115)
    print()
    print(f"  wrote {p}")


if __name__ == "__main__":
    main()
