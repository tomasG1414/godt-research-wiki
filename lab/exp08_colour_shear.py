"""exp08 -- Wavelength as a shear channel: colour-multiplexed multi-shear GODT.

The polarization-grating shear scales with wavelength, d(lam) = d0 lam/lam0
(solution note E2). Under one broadband LED the colours are averaged and the
nulls survive; resolved, each colour is its own shear magnitude with no moving
parts. This script prices that, with the per-wavelength geometry the first
back-of-envelope run held fixed:

  1. FIXED GEOMETRY (the upper-bound model): every colour on the lam = 0.5 um
     caps at equal rho. Ties this script to exp03/exp04's published rows.
  2. PER-WAVELENGTH GEOMETRY: each colour gets its own Ewald radius, band edge,
     cap density M_c(K) and Rytov weight, and its own per-photon information
     rho_c = rho (lam_ref/lam_c)^2  (O = k0^2 f, so a photon at shorter
     wavelength says more about the dispersion-free contrast f = n_m^2 - n^2).
     The separation J = rho S_t(K_perp) M(K) no longer holds across colours:
         J(K) = sum_c rho_c M_c(K) sum_j t_cj |G_cj(K_perp)|^2.
     Scored on the 530 nm support so every design answers the same question;
     what a colour adds outside it is reported separately.
  3. Does colour diversity touch the missing cone?  (No.)
  4. Colour time shares: simplex scan for the three-colour x + y design.
  5. Sample dispersion as coherent model error: SNR ceiling 1/eps^2 per
     off-reference channel (same construction as E1's calibration error).
  6. Single-shot mosaic sensor: pixel-area shares and crosstalk unmixing
     (shearlib.unmix_shares, a MODEL), and each LED's own bandwidth.
  7. At the design-rule shear (no in-band nulls) colour buys nothing as shear
     diversity; what is left is the geometry.

Stated models, none sourced: LED FWHMs, the crosstalk matrix, the dispersion
levels. Assumptions A1-A7 of the design problem note otherwise.

Writes: out/colour_shear.png
Wiki:   wiki/notes/colour-multiplexed-shear.md

Run: python exp08_colour_shear.py     (~1 min; numpy + matplotlib only, no cache)
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from shearlib import (GODT_SIM, KGrid, G2_vec, hermitian_sym, cap_density, chromatic_G,
                      at_wavelength, unmix_shares, nulls_in_band, out_path)

RHO = 100.0                       # per frame-sample SNR at lam = LAM_RHO, as exp03/exp04
LAM_RHO = 0.50                    # um; the wavelength exp03/exp04's rho = 100 refers to
N_AZ, N_FRAMES = 1440, 180
LAM0, D0 = 0.53, 0.8              # um; published shear at the LED centre wavelength
B, G, R = 0.45, 0.53, 0.635       # um; three common LED lines
FWHM = {B: 0.020, G: 0.035, R: 0.017}     # um; G is GODT's, B and R are typical (model)
X, Y = (1.0, 0.0), (0.0, 1.0)
LOST = 0.25                       # relative MSE above which a voxel counts as lost
RULE = "=" * 78


def shear(lam, d0=D0):
    """Grating shear at wavelength lam for a module that gives d0 at LAM0."""
    return d0 * lam / LAM0


def ch(lam, u, t=1.0, d=None, geom=None, fwhm=0.0):
    """One channel: colour, direction, share, shear (default: the grating's),
    the wavelength whose geometry it is scored on (default: its own)."""
    return {"lam": lam, "u": u, "t": t, "d": shear(lam) if d is None else d,
            "geom": lam if geom is None else geom, "fwhm": fwhm}


class Stack:
    """Per-wavelength cap densities on one K grid, flattened over the union
    of their supports so a design evaluates in a few vector operations."""

    def __init__(self, grid, lams):
        self.grid = grid
        self.cfg = {l: at_wavelength(GODT_SIM, l) for l in lams}
        M = {}
        for l, c in self.cfg.items():
            grid.check(c)
            M[l] = hermitian_sym(cap_density(c, grid, N_AZ, N_FRAMES))
        self.supp3 = {l: m > 0 for l, m in M.items()}
        self.idx = np.flatnonzero(np.any(list(self.supp3.values()), axis=0))
        self.ixy = self.idx // grid.nz
        # A_c = rho_c M_c with rho_c = RHO (LAM_RHO / lam_c)^2
        self.A = {l: RHO * (LAM_RHO / l) ** 2 * m.ravel()[self.idx] for l, m in M.items()}
        self.supp = {l: a > 0 for l, a in self.A.items()}
        self._E = {}

    def E(self, c):
        """Dimensionless shear power of a channel on the flattened voxels."""
        key = (round(c["d"], 9), c["u"], c["lam"], c["fwhm"])
        if key not in self._E:
            g = self.grid
            if c["fwhm"] > 0:        # the LED's own spectrum, averaged coherently (E2)
                Kpar = (g.KX * c["u"][0] + g.KY * c["u"][1]).ravel()
                Ge, _ = chromatic_G(Kpar, c["d"], c["lam"], c["fwhm"], n_lam=31)
                e = np.abs(Ge) ** 2
            else:
                e = G2_vec(g.KX, g.KY, (c["d"] * c["u"][0], c["d"] * c["u"][1])).ravel()
            self._E[key] = e[self.ixy]
        return self._E[key]

    def relmse(self, chans, eps=0.0, ref=G):
        """Relative Wiener MSE 1/(J + 1) per voxel, flat prior, shares
        normalised to 1. eps: rms relative dispersion of every channel whose
        colour is not `ref`, as coherent model error (ceiling 1/eps^2)."""
        tot = sum(c["t"] for c in chans)
        J = np.zeros(self.idx.size)
        for c in chans:
            term = (c["t"] / tot) * self.A[c["geom"]] * self.E(c)
            if eps > 0 and c["lam"] != ref:
                term = term / (1.0 + term * eps ** 2)
            J += term
        return 1.0 / (J + 1.0)

    def to3d(self, flat, fill=np.nan):
        out = np.full(self.grid.size, fill)
        out[self.idx] = flat
        return out.reshape(self.grid.shape)


def score(st, chans, ref, **kw):
    """(mean rel. MSE, % lost) on the support of `ref`, the voxels recovered
    outside it as a % of it, and the % of it no channel of the design samples
    at all (supports of different colours are not nested, see section 3)."""
    r = st.relmse(chans, **kw)
    m = st.supp[ref]
    seen = np.any([st.supp[c["geom"]] for c in chans], axis=0)
    return (r[m].mean(), 100 * np.mean(r[m] > LOST), 100 * np.sum(r[~m] <= LOST) / m.sum(),
            100 * np.mean(~seen[m]))


def main():
    grid = KGrid(320, 0.08, 96)       # same K lattice as exp03's 256 x 0.1, wider for blue
    st = Stack(grid, (LAM_RHO, B, G, R))
    dB, dG, dR = shear(B), shear(G), shear(R)

    # ------------------------------------------------ 1. fixed geometry
    print(RULE)
    print("1. FIXED GEOMETRY (lam = 0.5 um caps and rho for every colour): the upper bound")
    print(RULE)
    print(f"  grating shear d(lam) = {D0} um x lam/{1e3*LAM0:.0f} nm: "
          f"{1e3*B:.0f} nm -> {dB:.3f}, {1e3*G:.0f} -> {dG:.3f}, {1e3*R:.0f} -> {dR:.3f} um "
          f"(ratio R/B {dR/dB:.3f})")
    f = LAM_RHO
    fixed = [
        ("0.8 x single (published)", [ch(G, X, geom=f)]),
        ("2 colours B+R, x", [ch(B, X, geom=f), ch(R, X, geom=f)]),
        ("3 colours B+G+R, x", [ch(l, X, geom=f) for l in (B, G, R)]),
        ("mechanical 0.8 + 0.55, x", [ch(G, X, geom=f), ch(G, X, d=0.55, geom=f)]),
        ("0.8 x + 0.8 y", [ch(G, X, geom=f), ch(G, Y, geom=f)]),
        ("2 colours, x + y", [ch(l, u, geom=f) for l in (B, R) for u in (X, Y)]),
        ("3 colours, x + y", [ch(l, u, geom=f) for l in (B, G, R) for u in (X, Y)]),
        ("mechanical {0.8, 0.55}, x + y", [ch(G, u, d=d, geom=f) for d in (0.8, 0.55) for u in (X, Y)]),
    ]
    print(f"  {'design':32s} {'mean relMSE':>12} {'% lost':>8}")
    for lbl, chans in fixed:
        m, lost = score(st, chans, f)[:2]
        print(f"  {lbl:32s} {m:12.5f} {lost:8.2f}")
    print("  (rows 1, 4, 5, 8 are exp03/exp04's 0.0342 / 0.0070 / 0.0027 / 0.00097)")

    # ------------------------------------------- 2. per-wavelength geometry
    print()
    print(RULE)
    print("2. PER-WAVELENGTH GEOMETRY, scored on the 530 nm support")
    print(RULE)
    for l in (B, G, R):
        c = st.cfg[l]
        lo, hi = c.Kz_range
        print(f"  {1e3*l:.0f} nm: shear {shear(l):.3f} um, kmax {c.kmax:5.2f} rad/um, K_z in "
              f"[{lo:5.2f}, {hi:5.2f}], rho {RHO*(LAM_RHO/l)**2:5.1f}, support "
              f"{100*st.supp[l].sum()/st.supp[G].sum():5.1f} % of 530's, in-band nulls at "
              + ", ".join(f"{k:.2f}" for k in nulls_in_band(shear(l), c.kmax)))
    designs = [
        ("530 single, x", [ch(G, X)]),
        ("B+G+R x, SAME shear 0.8 (control)", [ch(l, X, d=0.8) for l in (B, G, R)]),
        ("mechanical 0.8 + 0.55 at 530, x", [ch(G, X), ch(G, X, d=0.55)]),
        ("2 colours B+R, x", [ch(B, X), ch(R, X)]),
        ("2 colours B+G, x", [ch(B, X), ch(G, X)]),
        ("3 colours B+G+R, x", [ch(l, X) for l in (B, G, R)]),
        ("530, x + y", [ch(G, X), ch(G, Y)]),
        ("mechanical {0.8, 0.55} at 530, x + y", [ch(G, u, d=d) for d in (0.8, 0.55) for u in (X, Y)]),
        ("2 colours B+R, x + y", [ch(l, u) for l in (B, R) for u in (X, Y)]),
        ("2 colours B+G, x + y", [ch(l, u) for l in (B, G) for u in (X, Y)]),
        ("3 colours B+G+R, x + y", [ch(l, u) for l in (B, G, R) for u in (X, Y)]),
    ]
    D = dict(designs)
    print()
    print(f"  {'design (equal shares)':38s} {'mean relMSE':>12} {'% lost':>8} {'% unsampled':>12} {'% gained outside':>17}")
    res = {}
    for lbl, chans in designs:
        res[lbl] = score(st, chans, G)
        print(f"  {lbl:38s} {res[lbl][0]:12.5f} {res[lbl][1]:8.2f} {res[lbl][3]:12.2f} {res[lbl][2]:17.2f}")
    b = res["mechanical {0.8, 0.55} at 530, x + y"][0]
    for lbl in ("3 colours B+G+R, x + y", "2 colours B+G, x + y"):
        print(f"  {lbl} vs the mechanical pair x + y: {res[lbl][0]/b:.2f}x the error on the same support")
    r3 = st.relmse(D["3 colours B+G+R, x + y"])
    rm = st.relmse(D["mechanical {0.8, 0.55} at 530, x + y"])
    KP = grid.KP.ravel()[st.ixy]
    m = st.supp[G]
    print("  mean relMSE by |K_perp| band (rad/um):   3 colours x+y   mechanical x+y")
    edges = [0, 10, 20, st.cfg[R].kmax, st.cfg[G].kmax + 1e-9]
    for lo_, hi_ in zip(edges[:-1], edges[1:]):
        s = m & (KP >= lo_) & (KP < hi_)
        print(f"    {lo_:5.1f} - {hi_:5.1f} ({100*s.sum()/m.sum():4.1f} % of support): "
              f"{r3[s].mean():12.5f} {rm[s].mean():15.5f}")

    # ------------------------------------------------ 3. the missing cone
    print()
    print(RULE)
    print("3. DOES COLOUR DIVERSITY REACH INTO THE MISSING CONE?")
    print(RULE)
    KPg = grid.KP
    print("  smallest sampled |K_perp| at fixed K_z (the cone's inner boundary), rad/um:")
    print(f"  {'K_z':>6} " + " ".join(f"{1e3*l:>7.0f}nm" for l in (B, G, R)) + f" {'union':>9}")
    for kzv in (1.0, 2.0, 4.0, 6.0, 8.0):
        iz = int(np.rint(kzv / grid.dk))
        row = []
        for l in (B, G, R):
            sl = st.supp3[l][:, :, iz]
            row.append(KPg[sl].min() if sl.any() else np.nan)
        print(f"  {grid.kz[iz]:6.2f} " + " ".join(f"{v:9.2f}" for v in row) + f" {np.nanmin(row):9.2f}")
    KX, KY, KZ = grid.K3()
    KP3 = np.hypot(KX, KY)
    cG = st.cfg[G]
    kzmax = max(abs(cG.Kz_range[0]), cG.Kz_range[1])
    th = np.deg2rad(cG.theta_ill)
    cone = (KP3 < np.abs(KZ) / np.tan(th)) & (np.abs(KZ) <= kzmax) & (np.abs(KZ) > 0)
    union = st.supp3[B] | st.supp3[G] | st.supp3[R]
    print(f"  nominal cone K_perp < |K_z| cot(theta), |K_z| <= {kzmax:.2f} ({cone.sum()} voxels): sampled by")
    print("    " + ", ".join(f"{1e3*l:.0f} nm {100*np.mean(st.supp3[l][cone]):.1f} %" for l in (B, G, R))
          + f", all three {100*np.mean(union[cone]):.1f} %")
    nG = st.supp3[G].sum()
    new = union & ~st.supp3[G]
    print(f"  voxels B and R add to the 530 support: {100*new.sum()/nG:.1f} % of it, of which "
          f"{100*(new & cone).sum()/new.sum():.1f} % lie in that cone")
    for l in (B, R):
        sl = st.supp3[G] & ~st.supp3[l]
        print(f"  530 support NOT sampled at {1e3*l:.0f} nm: {100*sl.sum()/nG:.2f} % of it "
              f"({100*(sl & cone).sum()/max(sl.sum(), 1):.0f} % of that in the cone, "
              f"{100*(sl & (KP3 > st.cfg[l].kmax)).sum()/max(sl.sum(), 1):.0f} % beyond that colour's band edge)")
    print("  => to first order in K_z the boundary K_perp = K_z cot(theta) is wavelength-free:")
    print("     every cap passes through the origin at the same angle. At second order the")
    print("     boundary bends inward as K_z^2 / (2 k_m sin^3 theta), MORE for longer wavelengths,")
    print("     so supports are not nested: a longer wavelength reaches slightly further into")
    print("     the cone at a given K_z, a shorter one reaches further out and further up.")

    # ------------------------------------------------ 4. colour shares
    print()
    print(RULE)
    print("4. COLOUR TIME SHARES for B+G+R x + y (simplex, step 0.05; x and y equal)")
    print(RULE)
    Bc = {l: st.A[l] * 0.5 * (st.E(ch(l, X)) + st.E(ch(l, Y))) for l in (B, G, R)}
    m = st.supp[G]
    best = (np.inf, None)
    grid_t = np.round(np.arange(0, 1.0001, 0.05), 2)
    named = {}
    for tb in grid_t:
        for tg in grid_t:
            tr = round(1 - tb - tg, 2)
            if tr < -1e-9:
                continue
            val = (1.0 / (tb * Bc[B] + tg * Bc[G] + tr * Bc[R] + 1.0))[m].mean()
            named[(tb, tg, tr)] = val
            if val < best[0]:
                best = (val, (tb, tg, tr))
    print(f"  best shares (B, G, R) = ({best[1][0]:.2f}, {best[1][1]:.2f}, {abs(best[1][2]):.2f}): mean relMSE {best[0]:.5f}")
    for lbl, t in (("equal thirds", None), ("Bayer areas (1/4, 1/2, 1/4)", (0.25, 0.5, 0.25)),
                   ("blue only", (1.0, 0.0, 0.0)), ("green only", (0.0, 1.0, 0.0)),
                   ("blue + red halves", (0.5, 0.0, 0.5)), ("blue + green halves", (0.5, 0.5, 0.0)),
                   ("green + red halves", (0.0, 0.5, 0.5))):
        v = (1.0 / ((Bc[B] + Bc[G] + Bc[R]) / 3 + 1.0))[m].mean() if t is None else named[t]
        print(f"    {lbl:30s} {v:.5f}")

    # ------------------------------------------------ 5. dispersion
    print()
    print(RULE)
    print("5. SAMPLE DISPERSION as coherent model error on every non-530 channel")
    print(RULE)
    print("  eps = rms relative change of the contrast f between that colour and 530 nm;")
    print("  per-voxel SNR of the channel is capped at 1/eps^2 whatever the photon count.")
    print(f"  photon-limited rho M |G|^2 at 530, median M, |G|^2 = 4: "
          f"{10*np.log10(4*np.median(st.A[G][st.supp[G]])):.1f} dB")
    print(f"  {'eps':>6} {'ceiling':>8} {'B+G x':>10} {'B+G+R x':>10} {'B+G x+y':>10} {'B+G+R x+y':>10}")
    for eps in (0.0, 0.01, 0.02, 0.05, 0.10, 0.20):
        v = [score(st, D[k], G, eps=eps)[0] for k in ("2 colours B+G, x", "3 colours B+G+R, x",
                                                     "2 colours B+G, x + y", "3 colours B+G+R, x + y")]
        ceil = "inf" if eps == 0 else f"{-20*np.log10(eps):.0f} dB"
        print(f"  {eps:6.2f} {ceil:>8} " + " ".join(f"{x:10.5f}" for x in v))
    print(f"  references on the same support: 530 single x {res['530 single, x'][0]:.5f}, "
          f"530 x + y {res['530, x + y'][0]:.5f}")

    # ------------------------------------------------ 6. single shot
    print()
    print(RULE)
    print("6. SINGLE-SHOT RGB: mosaic area shares, crosstalk unmixing, LED bandwidths")
    print(RULE)
    area = np.array([0.25, 0.5, 0.25])            # B, G, R pixel fractions (Bayer)
    print(f"  {'adjacent crosstalk':>19} {'eff. shares B, G, R':>24} {'sum':>6} {'3 col x':>10} {'3 col x+y':>10}")
    for x in (0.0, 0.1, 0.2, 0.3, 0.4):
        Xm = np.array([[1, x, x * x], [x, 1, x], [x * x, x, 1]], dtype=float)
        t = unmix_shares(Xm, area)
        # shares are fractions of ONE exposure; do not renormalise (relmse would)
        Jx = sum(t[i] * st.A[l] * st.E(ch(l, X)) for i, l in enumerate((B, G, R)))
        Jxy = sum(t[i] * Bc[l] for i, l in enumerate((B, G, R)))
        print(f"  {x:19.1f} {'  '.join(f'{v:.3f}' for v in t):>24} {t.sum():6.3f} "
              f"{(1/(Jx+1))[m].mean():10.5f} {(1/(Jxy+1))[m].mean():10.5f}")
    print("  (x + y here is two exposures of half the time each, i.e. two physical configurations)")
    print("  NOT modelled: the mosaic samples each colour on a coarser lattice (R, B at twice the")
    print("  pixel pitch), so the image must be oversampled 2x or the colour planes alias.")
    print()
    print("  each LED's own bandwidth (coherent average inside a colour, as E2):")
    for lbl, lams, dirs in (("3 colours, x", (B, G, R), (X,)), ("3 colours, x + y", (B, G, R), (X, Y)),
                            ("530 single, x", (G,), (X,))):
        mono = score(st, [ch(l, u) for l in lams for u in dirs], G)
        led = score(st, [ch(l, u, fwhm=FWHM[l]) for l in lams for u in dirs], G)
        print(f"    {lbl:18s} monochromatic {mono[0]:.5f} / {mono[1]:.2f} % lost   "
              f"LED lines {led[0]:.5f} / {led[1]:.2f} % lost")

    # ------------------------------------------------ 7. design-rule shear
    print()
    print(RULE)
    print("7. AT THE DESIGN-RULE SHEAR (0.2 um at 530 nm: no in-band nulls)")
    print(RULE)
    for lbl, chans in (("530, x + y", [ch(G, u, d=0.2) for u in (X, Y)]),
                       ("3 colours, x + y", [ch(l, u, d=shear(l, 0.2)) for l in (B, G, R) for u in (X, Y)]),
                       ("3 colours x + y, fixed geometry", [ch(l, u, d=shear(l, 0.2), geom=G) for l in (B, G, R) for u in (X, Y)]),
                       ("B+G, x + y", [ch(l, u, d=shear(l, 0.2)) for l in (B, G) for u in (X, Y)]),
                       ("blue only, x + y", [ch(B, u, d=shear(B, 0.2)) for u in (X, Y)])):
        s = score(st, chans, G)
        print(f"  {lbl:34s} {s[0]:.5f}   {s[1]:.2f} % lost ({s[3]:.2f} % unsampled)   {s[2]:.1f} % gained outside")
    print("  => with no comb to fill, the shear diversity is worth nothing; what moves the")
    print("     number is which colour's band and per-photon information you bought.")

    # ----------------------------------------------------------- figure
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
    K = np.linspace(0, st.cfg[B].kmax, 2000)
    for l, col in ((B, "tab:blue"), (G, "tab:green"), (R, "tab:red")):
        ax[0].plot(K, G2_vec(K, 0 * K, (shear(l), 0)), color=col, lw=1,
                   label=f"{1e3*l:.0f} nm, d = {shear(l):.2f} um")
        ax[0].axvline(st.cfg[l].kmax, color=col, ls=":", lw=1)
    ax[0].plot(K, sum(G2_vec(K, 0 * K, (shear(l), 0)) for l in (B, G, R)) / 3, "k", lw=1.6,
               label="equal-share sum")
    ax[0].set(xlabel="K_x (rad/um)", ylabel="|G|^2", title="shear power per colour (dotted: band edges)")
    ax[0].legend(fontsize=7)
    ax[0].grid(alpha=.3)
    ext = [grid.kx.min(), -grid.kx.min(), grid.kz.min(), -grid.kz.min()]
    for a, (lbl, chans) in zip(ax[1:], (designs[0], ("3 colours B+G+R, x", D["3 colours B+G+R, x"]))):
        img = np.fft.fftshift(st.to3d(st.relmse(chans))[:, 0, :]).T
        im = a.imshow(img, origin="lower", extent=ext, aspect="auto", cmap="magma", vmin=0, vmax=1)
        a.set(title=f"rel. MSE, K_y = 0: {lbl}", xlabel="K_x (rad/um)", ylabel="K_z (rad/um)",
              ylim=(-10, 10))
        plt.colorbar(im, ax=a, fraction=0.046)
    fig.tight_layout()
    fig.savefig(out_path("colour_shear.png"), dpi=110)
    print()
    print(f"  wrote {out_path('colour_shear.png')}")


if __name__ == "__main__":
    main()
