"""exp09 -- An exact Helmholtz referee for exp06's forward models.

exp06 left two things provisional. Its WPM forward model disagreed with first
Born by ~50 % even in the weak-scattering limit, with an unexplained excess at
steep angles; and its multiple-scattering column D transplanted WPM's
nonlinear excess onto Born without any way to check that excess. This script
adds an exact referee and settles both.

  1. THE REFEREE. shearlib.CBS (convergent Born series, Osnabrugge et al. 2016)
     solves the full scalar Helmholtz equation. Checked against the exact
     scalar partial-wave (Mie) series for a sphere far outside first Born,
     alongside first Born, WPM and multi-slice Born (shearlib.MSB).
  2. WPM'S LINEAR LIMIT, IN A WINDOW MATCHED TO BORN'S. exp06 compared WPM on a
     19.2 um window with an absorbing edge, cropped to 12.8 um, against Born on
     a 12.8 um periodic window. Steep rays walk 10-25 um sideways, so the two
     windows disagree about them. On one window the law is exact:
     WPM / Born = (k_z/k_iz) x (slice-placement phase).
  3. WHAT D TRANSPLANTS. exp06's cell phantom, one illumination, contrasts 0.5
     to 2: the nonlinear excess e(c) = d(c) - (c/eps) d(eps) of WPM and MSB
     against the exact one. Measured on the whole simulation window at the
     reconstruction lattice's frequencies (the infinite-field-of-view value
     periodic Born gives), because a 12.8 um crop changes even the exact
     solver's answer by tens of percent.

Deterministic. ~45 min on 4 cores, ~6 GB.
Wiki: wiki/notes/full-wave-referee.md
Run:  python exp09_full_wave_referee.py
"""
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

from shearlib import (GODT_SIM, KGrid, WPM, MSB, CBS, mie_scalar, born_sphere,
                      snap_illumination, cap_geometry, born_rytov_ratio, meas_scale,
                      cell_phantom, lattice_samples, out_path)
from exp06_multiple_scattering import SEED, PHANTOM_SCALE, EPS

RULE = "=" * 86
CBS_KW = dict(nb=30, atten=8.0, frame=0.18, frame_atten=8.0)
Z0, NZC = 28, 72          # the z slices of exp06's 128 that hold the whole cell


def embed(v, pad, fill=0.0):
    out = np.full((pad, pad, v.shape[2]), fill, dtype=v.dtype)
    o = (pad - v.shape[0]) // 2
    out[o:o + v.shape[0], o:o + v.shape[0]] = v
    return out


def crop(f, n=128):
    o = (f.shape[0] - n) // 2
    return f[o:o + n, o:o + n]


def rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


# ======================================================================= 1
def stage_mie(cfg):
    print(RULE)
    print("1. THE REFEREE AGAINST THE EXACT SCALAR MIE SERIES")
    print(RULE)
    lam, n0 = cfg.lam, cfg.n0
    k0, km = 2 * np.pi / lam, cfg.k_m
    n, dx, nz, a = 256, 0.1, 40, 1.5
    x = (np.arange(n) - n // 2) * dx
    z = (np.arange(nz) - nz // 2) * dx
    X, Y, Z = np.meshgrid(x, x, z, indexing="ij")
    frac = np.zeros(X.shape)
    ss = (np.arange(5) - 2) / 5 * dx               # partial-volume sphere, 5^3 subsamples
    for sx in ss:
        for sy in ss:
            for sz in ss:
                frac += (X + sx) ** 2 + (Y + sy) ** 2 + (Z + sz) ** 2 <= a * a
    frac /= 125
    k = 2 * np.pi * np.fft.fftfreq(n, dx)
    KX, KY = np.meshgrid(k, k, indexing="ij")
    kzs = np.sqrt(np.maximum(km ** 2 - KX ** 2 - KY ** 2, 0.0))
    shift = np.exp(-1j * (KX + KY) * x[0])
    KI = snap_illumination(cfg, KGrid(n, dx, n), 8)[0]
    print(f"  sphere a = {a} um on a {n*dx:.1f} um window, absorbing frame outside "
          f"|x| = {0.5*n*dx*(1-2*0.18):.1f} um; far-field amplitude f(theta) from the angular")
    print("  spectrum of the scattered field, f = k_z A(k_perp) / (2 pi i). Relative error vs Mie:")
    fig, ax0 = plt.subplots(figsize=(6.4, 4.2))
    rows = []
    for n1 in (1.34, 1.38):
        m = n1 / n0
        O = frac * k0 ** 2 * (n0 ** 2 - n1 ** 2)
        nvol = np.sqrt(n0 ** 2 - O / k0 ** 2)
        models = {"CBS": None, "MSB": MSB(O, dx, dx, km, absorb_frac=0.18),
                  "WPM": WPM(nvol, dx, dx, lam, n0, n_levels=12, absorb_frac=0.18,
                             dtype=np.complex128)}
        cb = CBS(O, dx, dx, km, **CBS_KW)
        print(f"  n_sphere = {n1}: phase delay 2 k0 a dn = "
              f"{2*k0*a*(n1-n0):.2f} rad")
        print(f"  {'incidence':>9} {'band':>15} {'rays':>6} {'CBS':>8} {'MSB':>8} {'WPM':>8} {'Born':>8}")
        for inc, ki in (("normal", (0.0, 0.0)), ("60 deg", tuple(KI[0]))):
            kiz = np.sqrt(km ** 2 - ki[0] ** 2 - ki[1] ** 2)
            ramp = np.exp(1j * (ki[0] * X[:, :, 0] + ki[1] * Y[:, :, 0]))
            ps, info = cb.solve(ki, tol=1e-7)
            ratios = {"CBS": cb.ratio(ps, ki), "MSB": models["MSB"].ratio(ki)}
            ratios["WPM"] = models["WPM"].to_centre(models["WPM"].run(ki), ki)
            f = {mm: (kzs * dx * dx * np.fft.fft2(r * ramp - ramp) * shift / (2j * np.pi))
                 for mm, r in ratios.items()}
            mu = (KX * ki[0] + KY * ki[1] + kzs * kiz) / km ** 2
            polar = np.degrees(np.arccos(np.clip(kzs / km, -1, 1)))
            fm = mie_scalar(km, m, a, mu)
            fb = born_sphere(km, m, a, mu)
            th = np.degrees(np.arccos(np.clip(mu, -1, 1)))
            # rays whose sideways walk to the measurement plane stays clear of the frame
            clear = (kzs > 0.1 * km) & (polar < (55 if inc == "normal" else 62))
            for band, (lo, hi) in (("all clear", (0, 180)), ("theta < 20", (0, 20)),
                                   ("20-40 deg", (20, 40)), ("40-55 deg", (40, 55))):
                s = clear & (th >= lo) & (th < hi)
                if s.sum() < 20:
                    continue
                e = [rel(f[mm][s], fm[s]) for mm in ("CBS", "MSB", "WPM")] + [rel(fb[s], fm[s])]
                rows.append((n1, inc, band, *e))
                print(f"  {inc:>9} {band:>15} {s.sum():6d} " + " ".join(f"{v:8.4f}" for v in e))
            i0 = np.unravel_index(np.argmax(np.where(clear, mu, -2)), mu.shape)
            print(f"  {inc:>9}   forward f: Mie {fm[i0]:.3f}  CBS {f['CBS'][i0]:.3f}  "
                  f"Born {fb[i0]:.3f}   (CBS {info['iterations']} iterations, "
                  f"Helmholtz residual {info['residual']:.1e})")
            if inc == "normal" and n1 == 1.38:
                ax = ax0
                s = clear & (KY == 0) & (KX >= 0)
                o = np.argsort(th[s])
                tt = np.linspace(0, 55, 400)
                ax.semilogy(tt, np.abs(mie_scalar(km, m, a, np.cos(np.radians(tt)))), "k-",
                            lw=2, label="exact (scalar Mie)")
                ax.semilogy(tt, np.abs(born_sphere(km, m, a, np.cos(np.radians(tt)))), "k:",
                            label="first Born")
                for mm, mk in (("CBS", "o"), ("MSB", "s"), ("WPM", "^")):
                    ax.semilogy(th[s][o], np.abs(f[mm][s][o]), mk, ms=3.5, label=mm)
                ax.set(xlabel="scattering angle (deg)", ylabel="|f| (um)",
                       title=f"sphere a = {a} um, n = {n1} in {n0}, normal incidence")
                ax.legend(fontsize=8)
                ax.grid(alpha=.3, which="both")
    print("  => CBS reproduces the exact series where first Born is off by tens of percent;")
    print("     that is the referee every comparison below is scored against.")
    return fig


# ======================================================================= 2
def stage_window(cfg):
    print()
    print(RULE)
    print("2. WPM'S LINEAR LIMIT, ON A WINDOW MATCHED TO BORN'S")
    print(RULE)
    grid = KGrid(128, 0.1, 128)
    k0, nmed = 2 * np.pi / cfg.lam, cfg.n0
    sc = meas_scale(cfg, grid)
    a_um = 4.8 * PHANTOM_SCALE
    xs = (np.arange(grid.n) - grid.n // 2) * grid.dx
    Xp, Yp, Zp = np.meshgrid(xs, xs, xs, indexing="ij")
    Rp = np.sqrt(Xp ** 2 + Yp ** 2 + Zp ** 2)
    pr = gaussian_filter(np.random.default_rng(5).standard_normal(Rp.shape), 0.15 / grid.dx)
    shape = (pr / pr.std()) * 0.5 * (1 - np.tanh((Rp - a_um) / 0.25))   # exp06's probe, unit amp

    def gap(amp, pad, absorb, sub=1):
        probe = amp * shape
        KI = snap_illumination(cfg, grid, 8, dk=2 * np.pi / (pad * grid.dx))[0]
        ki = KI[0]
        cap = cap_geometry(cfg, grid, ki)
        ix, iy = np.nonzero(cap["ok"])
        kz, kiz = cap["kz"][ix, iy], cap["kiz"]
        w = WPM(np.repeat(embed(nmed + probe, pad, nmed), sub, axis=2), grid.dx,
                grid.dz / sub, cfg.lam, nmed, n_levels=24, absorb_frac=absorb,
                dtype=np.complex128)
        ap = (sc * np.fft.fft2(np.log(crop(w.to_centre(w.run(ki), ki), grid.n))))[ix, iy]
        nq = nmed + w.step * np.rint(probe / w.step)        # Born of what WPM actually sees
        bp = (sc * np.fft.fft2(np.log(born_rytov_ratio(k0 ** 2 * (nmed ** 2 - nq ** 2),
                                                       cfg, grid, ki, cap))))[ix, iy]
        Kz = kz - kiz
        law = (kz / kiz) * sum(np.exp(-1j * Kz * j * grid.dz / sub) for j in range(1, sub + 1)) / sub
        return ap, bp, kz / kiz, law

    print("  The probe is exp06's (band-limited, smoothly windowed, dn amplitude 2e-3 inside a")
    print(f"  {a_um:.2f} um sphere). |WPM / Born| per k_z/k_iz bin, and the same divided by k_z/k_iz:")
    print(f"  {'':12} {'exp06: 19.2 um + absorber, cropped':>36} {'matched: 12.8 um periodic':>28}")
    print(f"  {'k_z/k_iz':>10} {'rays':>6} {'|W/B|':>8} {'/ obliquity':>12} {'|W/B|':>9} "
          f"{'/ obliquity':>12}")
    A6, B6, o6, _ = gap(2e-3, 192, 0.18)
    Am, Bm, om, Lm = gap(2e-4, 128, 0.0)
    for lo, hi in ((0.1, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 1.1), (1.1, 2.0)):
        s6 = (o6 >= lo) & (o6 < hi) & (np.abs(B6) > 0.05 * np.median(np.abs(B6)))
        sm = (om >= lo) & (om < hi) & (np.abs(Bm) > 0.05 * np.median(np.abs(Bm)))
        if s6.sum() > 20 and sm.sum() > 20:
            r6, rm = np.abs(A6 / B6), np.abs(Am / Bm)
            print(f"  {np.median(o6[s6]):10.3f} {s6.sum():6d} {np.median(r6[s6]):8.4f} "
                  f"{np.median(r6[s6] / o6[s6]):12.4f} {np.median(rm[sm]):9.4f} "
                  f"{np.median(rm[sm] / om[sm]):12.4f}")
    print(f"  (the matched column uses amplitude 2e-4; the law below says why)")
    print()
    print("  Complex law  WPM = Born x (k_z/k_iz) x (1/s) sum_j e^{-i K_z j dz/s}  on the matched")
    print("  window (s sub-slices per voxel; the sum is where WPM places each thin screen):")
    print(f"  {'dz (um)':>8} {'amp 2e-4':>10} {'amp 2e-3':>10}   relative difference to the law")
    for sub in (1, 2, 4, 8):
        e = []
        for amp in (2e-4, 2e-3):
            ap, bp, _, law = gap(amp, 128, 0.0, sub)
            e.append(rel(ap, bp * law))
        print(f"  {0.1/sub:8.4f} {e[0]:10.4f} {e[1]:10.4f}")
    print("  => On one window, WPM's linear limit is first Born times k_z/k_iz, exactly: the")
    print("     residual at 2e-4 is second order and does not move with dz. exp06's extra deficit")
    print("     of steep rays was its window (absorber plus crop), not propagation physics.")
    print("     exp06's probe amplitude (2e-3) is already nonlinear at steep angles, and that")
    print("     nonlinear part is NOT converged in dz -- it grows as the slices get thinner.")


# ======================================================================= 3
def stage_cell(cfg):
    print()
    print(RULE)
    print("3. WHAT exp06's MODEL D TRANSPLANTS, AGAINST THE EXACT EXCESS")
    print(RULE)
    grid = KGrid(128, 0.1, 128)
    k0, nmed, km = 2 * np.pi / cfg.lam, cfg.n0, cfg.k_m
    sc = meas_scale(cfg, grid)
    O1, dn1 = cell_phantom(grid, np.random.default_rng(SEED), scale=PHANTOM_SCALE)
    assert not np.any(dn1[:, :, :Z0]) and not np.any(dn1[:, :, Z0 + NZC:])
    # one illumination on a lattice that the 128, 192, 384 and 512 windows share
    KI, th = snap_illumination(cfg, grid, 1, dk=2 * np.pi / 6.4)
    ki = KI[0]
    cap = cap_geometry(cfg, grid, ki)
    ix, iy = np.nonzero(cap["ok"])
    kz = cap["kz"][ix, iy]
    tan_s = np.sqrt(np.maximum(km ** 2 - kz ** 2, 0)) / kz
    classes = (("all rays", tan_s >= 0), ("near-axis, tan < 0.8", tan_s < 0.8),
               ("steep, tan > 2", tan_s > 2))
    print(f"  exp06's cell, one illumination at {th[0]:.2f} deg; {ix.size} collected rays.")
    print("  Each model's measured Rytov spectrum d(c) on the pupil is taken two ways:")
    print("    full -- FFT of the whole simulation window at the 12.8 um lattice's frequencies")
    print("            (shearlib.lattice_samples): the infinite-field-of-view measurement that")
    print("            exp06's model B, periodic first Born, makes;")
    print("    crop -- FFT of the central 12.8 um, as exp06 measured WPM.")
    print("  e(c) = d(c) - (c/eps) d(eps) is the nonlinear excess; exp06's D = B + e(WPM, 19.2 um, crop).")
    print("  Referee: CBS on a 51.2 um window (frame outside |x| = 20.5 um), full.")

    def full(r):
        return sc * lattice_samples(np.log(r), grid.n, ix, iy)

    def crop_m(r):
        return (sc * np.fft.fft2(np.log(crop(r, grid.n))))[ix, iy]

    def O_of(c, pad):
        return embed(k0 ** 2 * (nmed ** 2 - (nmed + c * dn1[:, :, Z0:Z0 + NZC]) ** 2), pad)

    def n_of(c, pad):
        return embed(nmed + c * dn1[:, :, Z0:Z0 + NZC], pad, nmed)

    D = {}
    t0 = time.time()
    contrasts = (0.5, 1.0, 2.0)
    for c in (EPS,) + contrasts:
        O128 = k0 ** 2 * (nmed ** 2 - (nmed + c * dn1) ** 2)
        D[("Born", c)] = (sc * np.fft.fft2(np.log(born_rytov_ratio(O128, cfg, grid, ki, cap))))[ix, iy]
        for pad in ((384, 512) if c in (EPS, 0.5) else (512,)):
            cb = CBS(O_of(c, pad), grid.dx, grid.dz, km, nb=30, atten=8.0, frame=0.1,
                     frame_atten=8.0, dtype=np.complex64)
            ps, info = cb.solve(ki, tol=1e-6)
            r = cb.ratio(ps, ki)
            D[(f"CBS {pad}", c)], D[(f"CBS {pad} crop", c)] = full(r), crop_m(r)
            print(f"    CBS {pad} contrast {c}: {info['iterations']} iterations  "
                  f"[{time.time()-t0:.0f} s]", flush=True)
        for pad in (384, 512):
            r = MSB(O_of(c, pad), grid.dx, grid.dz, km, absorb_frac=0.1).ratio(ki)
            D[(f"MSB {pad}", c)] = full(r)
            w = WPM(n_of(c, pad), grid.dx, grid.dz, cfg.lam, nmed, n_levels=10, absorb_frac=0.1)
            D[(f"WPM {pad}", c)] = full(w.to_centre(w.run(ki), ki))
        w = WPM(n_of(c, 192), grid.dx, grid.dz, cfg.lam, nmed, n_levels=10)   # exp06's setup
        D[("WPM 192 crop = D", c)] = crop_m(w.to_centre(w.run(ki), ki))

    def lin(m, c):
        return (c / EPS) * D[(m, EPS)]

    def exc(m, c):
        return D[(m, c)] - lin(m, c)

    ref = "CBS 512"
    print()
    print("  (a) Linear parts against periodic first Born, ||lin - Born|| / ||Born||, contrast 0.5:")
    print(f"  {'model':>18}" + "".join(f"{lab:>22}" for lab, _ in classes))
    for m in ("CBS 512", "CBS 384", "MSB 512", "MSB 384", "WPM 512", "CBS 512 crop",
              "WPM 192 crop = D"):
        print(f"  {m:>18}" + "".join(f"{rel(lin(m, 0.5)[s], lin('Born', 0.5)[s]):22.3f}"
                                     for _, s in classes))
    print("  The full-window measurement converges on Born's (MSB and CBS to a few percent);")
    print("  a 12.8 um crop does not, even for the exact solver. WPM's own gap is k_z/k_iz.")
    print()
    print("  (b) Nonlinear excess against the referee, ||e - e_ref|| / ||e_ref||, and its size")
    print("      ||e_ref|| / ||lin_ref|| (top row):")
    for c in contrasts:
        eR = exc(ref, c)
        print(f"  contrast {c} (a dn/lambda = {4.8*PHANTOM_SCALE*c*np.abs(dn1).max()/cfg.lam:.3f})")
        print(f"  {'exact excess size':>18}" + "".join(
            f"{np.linalg.norm(eR[s])/np.linalg.norm(lin(ref, c)[s]):22.3f}" for _, s in classes))
        for m in ("CBS 384", "MSB 512", "MSB 384", "WPM 512", "WPM 384", "WPM 192 crop = D"):
            if (m, c) in D:
                print(f"  {m:>18}" + "".join(f"{rel(exc(m, c)[s], eR[s]):22.3f}" for _, s in classes))
    print("  'CBS 384' is the referee against itself on a smaller window: where it is small")
    print("  the excess is a property of the cell, where it is not (steep rays) nothing below")
    print("  can be adjudicated by this referee.")
    print()
    print("  (c) Model D's error at the level of the whole measurement, ||e_D - e_ref|| / ||d_ref||,")
    print("      next to the exact excess it stands in for, ||e_ref|| / ||d_ref||, near-axis rays:")
    s = classes[1][1]
    for c in contrasts:
        eR, dR = exc(ref, c)[s], D[(ref, c)][s]
        print(f"    c = {c:4.2f}: exact excess {np.linalg.norm(eR)/np.linalg.norm(dR):.3f}, "
              f"D's error {np.linalg.norm(exc('WPM 192 crop = D', c)[s]-eR)/np.linalg.norm(dR):.3f}, "
              f"MSB's error {np.linalg.norm(exc('MSB 384', c)[s]-eR)/np.linalg.norm(dR):.3f}")
    return D


def main():
    t_all = time.time()
    cfg = GODT_SIM
    fig = stage_mie(cfg)
    stage_window(cfg)
    stage_cell(cfg)
    fig.tight_layout()
    fig.savefig(out_path("referee_mie.png"), dpi=130)
    print()
    print(f"  wrote {out_path('referee_mie.png')}")
    print(f"  total {time.time()-t_all:.0f} s")


if __name__ == "__main__":
    main()
