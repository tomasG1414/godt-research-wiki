"""exp06 -- Do the exp03-05 design conclusions survive multiple scattering?

exp03-05 assumed A1: first-order Rytov, weak scattering. The measurement was
SYNTHESISED from the model it was then inverted with, so the only error was
photon noise. This script replaces the forward model with Brenner & Singer's
wave-propagation method (brenner1993) -- GODT's own forward engine, and the
one multi-slice method valid at GODT's 60 deg illumination -- and re-runs the
same estimator on data it did not generate.

Five forward models feed ONE unchanged Rytov/Ewald estimator, so every
difference between columns is forward-model error and nothing else:

  A  ideal-linear   data = h O~ at the binned voxel  (exp05's model, exactly)
  B  Born-exact     first Born/Rytov, continuum K_z, real pupil, real field
  B' WPM-linearised WPM run at contrast eps and scaled up: WPM's own physics
                    with every nonlinearity removed
  C  WPM            the real thing
  D  Born + MS      d_B + (d_C - d_B'): WPM's nonlinear EXCESS transplanted onto
                    the Born baseline

A->B is gridding (assumption A3). B->B' is the forward-model gap between
multi-slice propagation and Born, which turns out to dominate everything and is
a property of WPM rather than of the sample. B'->C is multiple scattering and
Rytov breakdown, isolated from the other two -- but on top of a reconstruction
B' has already ruined, so it cannot be read at the design level.

D is what makes the design-level question answerable. Deposition is linear in
the data, so subtracting WPM's own linearisation and adding the difference to
Born cancels the multi-slice gap exactly while keeping the genuine
multiple-scattering perturbation. B->D is then multiple scattering measured on
a reconstruction that works. (Caveat: the residual is still computed inside
WPM, so it inherits WPM's angular accuracy -- but it is the nonlinear part
only, which is what is wanted.)

Establishes (argument and numbers in the wiki note):
  1. The shear factor G(K_perp) multiplies the measured Rytov spectrum EXACTLY,
     whatever the scattering physics, so the null comb and the whole design
     geometry are structurally immune to multiple scattering. Verified to
     machine precision on WPM data at full contrast, not assumed.
  2. For a single slice, WPM's linear limit is first Born with the obliquity factor taken at the
     INCIDENT angle: |WPM/Born| = k_z/k_iz, exactly. Reconstructing WPM data
     with a Rytov/Ewald inverse therefore carries an angle-dependent error
     before any multiple scattering.
  3. Where the linearisation breaks, in a*dn/lambda, against chen1998.
  4. What multiple scattering costs the exp05 design comparison: a bias floor
     no design can remove, so absolute errors rise and the gaps compress.
  5. The shearing measurement is wrap-immune where standard Rytov ODT is not
     (guizar-sicairos2011 vs chen1998), quantified.

Equal photon budget throughout, as in exp05: a design's channels share one
exposure. Noise is accumulated once per channel at unit share and rescaled
(J(t) = t J, num(t) = t*signal + sqrt(t)*noise), so every design and every
forward model sees the SAME noise realisation and the comparison is paired.

Seeded (SEED). Writes: out/ms_linearity.png, out/ms_designs.png, out/ms_slices.png
Wiki:   wiki/notes/multiple-scattering-check.md

Run: python exp06_multiple_scattering.py          (~28 min)
     python exp06_multiple_scattering.py --quick  (~10 min, 60 illuminations)
"""
import sys
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from shearlib import (GODT_SIM, KGrid, WPM, hermitian_sym, cap_density,
                      snap_illumination, cap_geometry, born_rytov_ratio,
                      meas_scale, shear_rytov, cell_phantom, shell_spectrum,
                      hsym_c, fuse, nulls_in_band, out_path)

SEED = 11
RHO = 100.0             # per-frame-sample SNR, calibrated exactly as in exp05
N_ILL = 180             # GODT's own projection count; also the exposure budget
K_LO = 2.0              # every aggregate is on supp & |K| >= K_LO (exp05's trap)
PHANTOM_SCALE = 0.62    # cytoplasm radius 3.0 um, so the field stays inside the pad
PAD = 192               # WPM transverse grid (19.2 um); measurement cropped to 128
N_LEVELS = 10           # WPM index quantisation
EPS = 0.02              # contrast of the linearised-WPM control (model B')
MODELS = ("A ideal-linear", "B Born-exact", "B' WPM-linear", "C WPM", "D Born+MS")
SIMULATED = MODELS[:4]          # D is a linear combination of B, B' and C
RULE = "=" * 86

QUICK = "--quick" in sys.argv


def embed(vol, pad):
    """Centre a (n, n, ns) index volume in a (pad, pad, ns) background."""
    n = vol.shape[0]
    out = np.full((pad, pad, vol.shape[2]), vol.flat[0], dtype=vol.dtype)
    o = (pad - n) // 2
    out[o:o + n, o:o + n] = vol
    return out


def crop(field, n):
    o = (field.shape[0] - n) // 2
    return field[o:o + n, o:o + n]


def shear_factor(grid, dvec):
    if dvec is None:
        return np.ones(grid.KX.shape, complex)
    return np.exp(1j * (grid.KX * dvec[0] + grid.KY * dvec[1])) - 1.0


def channel_meas(logratio, ratio, grid, dvec):
    """2D spectrum of what a channel measures, in `meas_scale` units.

    dvec None -> the ODT benchmark: the Rytov field itself on its principal
    branch, which is what standard Rytov ODT holds before unwrapping.
    Otherwise the shear difference, which the interferogram delivers wrap-free.
    """
    return np.fft.fft2(logratio if dvec is None else
                       shear_rytov(ratio, grid.dx, dvec, wrap_safe=True))


def bin_add(dst, lin, vals):
    n = dst.size
    dst += (np.bincount(lin, weights=vals.real, minlength=n)
            + 1j * np.bincount(lin, weights=vals.imag, minlength=n)).astype(dst.dtype)


def main():
    t_all = time.time()
    rng = np.random.default_rng(SEED)
    cfg = GODT_SIM
    grid = KGrid(128, 0.1, 128)
    grid.check(cfg)
    Kabs = grid.Kabs()
    k0, nmed = 2 * np.pi / cfg.lam, cfg.n0
    n_ill = 60 if QUICK else N_ILL

    KI, th_eff = snap_illumination(cfg, grid, n_ill, dk=2 * np.pi / (PAD * grid.dx))
    caps = [cap_geometry(cfg, grid, ki) for ki in KI]
    sc = meas_scale(cfg, grid)

    print(RULE)
    print("SETUP")
    print(RULE)
    print(f"  {cfg!r}")
    print(f"  {grid!r}")
    print(f"  {n_ill} illuminations snapped to the PROPAGATION lattice, so the incident plane")
    print(f"  wave is exactly periodic under WPM's FFTs: realised inclination "
          f"{th_eff.mean():.2f} +- {th_eff.std():.2f} deg")
    print(f"  (nominal {cfg.theta_ill:.1f}); the realised value is used everywhere below.")
    print(f"  WPM grid {PAD}x{PAD}x{grid.nz} = {PAD*grid.dx:.1f} um, {N_LEVELS} index levels,")
    print(f"  absorbing edge outside |x| = 7.9 um; measurement cropped to the central "
          f"{grid.n*grid.dx:.1f} um.")

    # ------------------------------------------------------ phantom + prior
    O1, dn1 = cell_phantom(grid, rng, scale=PHANTOM_SCALE)
    a_um = 4.8 * PHANTOM_SCALE
    dn_max = float(np.abs(dn1).max())
    Ot1 = np.fft.fftn(O1)
    cover = np.zeros(grid.size, np.int32)
    for cap in caps:
        ixc, iyc = np.nonzero(cap["ok"])
        cover += np.bincount((ixc * grid.n + iyc) * grid.nz + cap["iz"][ixc, iyc],
                             minlength=grid.size).astype(np.int32)
    supp = hermitian_sym(cover.reshape(grid.shape)) > 0
    reach = hermitian_sym(cap_density(cfg, grid, 720, n_ill)) > 0
    ref = supp & (Kabs >= 8) & (Kabs <= 24)
    meas = supp & (Kabs >= K_LO)
    print(f"  phantom: cell, cytoplasm radius a = {a_um:.2f} um, dn in "
          f"[{dn1.min():+.4f}, {dn1.max():+.4f}] at contrast 1.0")
    print(f"  support actually covered by these {n_ill} projections: {supp.sum()} voxels "
          f"= {100*supp.sum()/reach.sum():.1f} % of what the geometry reaches")
    print("  (angular starvation, exp03's finding; every model below is scored on the covered")
    print(f"  set only, so it costs every model equally). Aggregates on supp & |K| >= {K_LO}: "
          f"{meas.sum()} voxels; rho = {RHO:.0f}")

    # ================================================================ STAGE 1
    print()
    print(RULE)
    print("1. WHAT WPM'S LINEAR LIMIT IS -- the forward-model gap, before any nonlinearity")
    print(RULE)
    ki0, cap0 = KI[0], caps[0]
    ix0, iy0 = np.nonzero(cap0["ok"])
    kz0, kiz0 = cap0["kz"][ix0, iy0], cap0["kiz"]
    Ka0 = np.sqrt(grid.KX[ix0, iy0] ** 2 + grid.KY[ix0, iy0] ** 2 + cap0["Kz"][ix0, iy0] ** 2)
    # A band-limited probe: the cell phantom's sharp edges alias at the band edge
    # and would be read as a forward-model gap that is really a sampling one.
    from scipy.ndimage import gaussian_filter
    xs = (np.arange(grid.n) - grid.n // 2) * grid.dx
    zs = (np.arange(grid.nz) - grid.nz // 2) * grid.dz
    Xp, Yp, Zp = np.meshgrid(xs, xs, zs, indexing="ij")
    Rp = np.sqrt(Xp ** 2 + Yp ** 2 + Zp ** 2)
    pr = gaussian_filter(np.random.default_rng(5).standard_normal(Rp.shape), 0.15 / grid.dx)
    probe = 2e-3 * (pr / pr.std()) * 0.5 * (1 - np.tanh((Rp - a_um) / 0.25))
    w_p = WPM(embed(nmed + probe, PAD), grid.dx, grid.dz, cfg.lam, nmed,
              n_levels=24, dtype=np.complex128)
    ap = (sc * np.fft.fft2(np.log(crop(w_p.to_centre(w_p.run(ki0), ki0), grid.n))))[ix0, iy0]
    bp = (sc * np.fft.fft2(np.log(born_rytov_ratio(
        k0 ** 2 * (nmed ** 2 - (nmed + probe) ** 2), cfg, grid, ki0, cap0))))[ix0, iy0]
    obl = kz0 / kiz0
    keep = np.abs(bp) > 0.05 * np.median(np.abs(bp))
    print(f"  Band-limited probe (dn_rms {probe.std():.1e}, deep in the linear regime) vs exact")
    print("  first Born -- same pupil, same continuum K_z, same illumination:")
    print(f"    relative difference                       "
          f"{np.linalg.norm(ap-bp)/np.linalg.norm(bp):.4f}")
    print(f"    after dividing out k_z/k_iz               "
          f"{np.linalg.norm(ap-bp*obl)/np.linalg.norm(bp*obl):.4f}")
    print(f"    median |WPM/Born| / (k_z/k_iz)            "
          f"{np.median(np.abs((ap/bp)[keep]) / obl[keep]):.4f}   (1.0 = exact)")
    print(f"  {'k_z/k_iz':>10} {'rays':>7} {'median |WPM/Born|':>19} {'ratio to k_z/k_iz':>19}")
    for lo, hi in ((0.1, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 1.1), (1.1, 2.0)):
        m = keep & (obl >= lo) & (obl < hi)
        if m.sum() > 20:
            r = np.abs((ap / bp)[m])
            print(f"  {np.median(obl[m]):10.3f} {m.sum():7d} {np.median(r):19.4f} "
                  f"{np.median(r/obl[m]):19.4f}")
    print("  Is that gap a step-size artefact? Same probe, finer WPM steps:")
    print(f"  {'dz (um)':>9} {'slices':>7} {'rel. difference to Born':>25}")
    for sub in (1, 2, 4):
        wv = WPM(np.repeat(embed(nmed + probe, PAD), sub, axis=2), grid.dx,
                 grid.dz / sub, cfg.lam, nmed, n_levels=24, dtype=np.complex128)
        av = (sc * np.fft.fft2(np.log(crop(wv.to_centre(wv.run(ki0), ki0), grid.n))))[ix0, iy0]
        print(f"  {grid.dz/sub:9.4f} {grid.nz*sub:7d} "
              f"{np.linalg.norm(av-bp)/np.linalg.norm(bp):25.4f}")
    print("  (A step-size artefact would fall toward zero here. Stage 4's B' and C share one")
    print("   step size, so their difference is unaffected either way.)")

    # the same gap, measured on the phantom the rest of the script uses
    w_e = WPM(embed(nmed + EPS * dn1, PAD), grid.dx, grid.dz, cfg.lam, nmed,
              n_levels=N_LEVELS, dtype=np.complex128)
    logr_eps = np.log(crop(w_e.to_centre(w_e.run(ki0), ki0), grid.n))
    a = (sc * np.fft.fft2(logr_eps))[ix0, iy0]
    b = (sc * np.fft.fft2(np.log(born_rytov_ratio(
        k0 ** 2 * (nmed ** 2 - (nmed + EPS * dn1) ** 2), cfg, grid, ki0, cap0))))[ix0, iy0]
    print(f"  On the cell phantom at contrast {EPS}: raw gap "
          f"{np.linalg.norm(a-b)/np.linalg.norm(b):.4f}, "
          f"after k_z/k_iz {np.linalg.norm(a-b*obl)/np.linalg.norm(b*obl):.4f}")
    print("  => WPM carries the obliquity factor of the INCIDENT ray, 1/k_iz, where first Born")
    print("     carries the SCATTERED one, 1/k_z. They agree for forward scattering and part")
    print("     company as k_z/k_iz leaves 1 -- i.e. exactly at the band edge the multi-shear")
    print("     design is about. This is a property of multi-slice propagation, not of this")
    print("     implementation: for a single slice WPM reduces analytically to")
    print("     exp(i k_z(n(x), k_i) dz), which the shearlib selftest checks to 1e-12.")
    print("  => so any validation that inverts WPM data with a Rytov/Ewald model")
    print("     includes this factor before any multiple scattering enters.")

    # ================================================================ STAGE 2
    print()
    print(RULE)
    print("2. THE NULLS ARE EXACTLY PRESERVED -- shearing is linear in the MEASURED field")
    print(RULE)
    print("  Shearing acts on u_R after all propagation, so")
    print("      FT[u_R(r+D) - u_R(r)](K_perp) = (e^{i K_perp.D} - 1) FT[u_R](K_perp)")
    print("  is an identity in the measured field, whatever u_R's relation to O. Multiple")
    print("  scattering cannot put signal where G = 0. Checked on the WPM field at contrast 1:")
    w_f = WPM(embed(nmed + dn1, PAD), grid.dx, grid.dz, cfg.lam, nmed,
              n_levels=N_LEVELS, dtype=np.complex128)
    r_ms = crop(w_f.to_centre(w_f.run(ki0), ki0), grid.n)
    uRk = np.fft.fft2(np.log(r_ms))
    print()
    print(f"  {'shear':>9} {'in-band nulls':>14} {'max|FT[shear]-G FT[u_R]|/max|FT[u_R]|':>39} "
          f"{'|data| on null cols':>20}")
    for d in (0.4, 0.8, 1.6):                 # nulls of these land exactly on the grid
        shk = np.fft.fft2(shear_rytov(r_ms, grid.dx, (d, 0.0)))
        G = shear_factor(grid, (d, 0.0))
        resid = np.max(np.abs(shk - G * uRk)) / np.max(np.abs(uRk))
        cols = np.zeros(grid.KX.shape, bool)
        for nu in nulls_in_band(d, cfg.kmax):
            cols |= np.abs(np.abs(grid.KX) - nu) < 1e-9
        on = np.max(np.abs(shk[cols])) / np.max(np.abs(shk))
        print(f"  {d:6.2f} um {len(nulls_in_band(d, cfg.kmax)):14d} {resid:39.2e} {on:20.2e}")
    print("  => machine precision at full scattering contrast. Which K are estimable, where the")
    print("     comb sits, and why incommensurate magnitudes fill it are not approximations that")
    print("     multiple scattering can invalidate. What it corrupts is the VALUE at each K.")

    # ================================================================ STAGE 3
    print()
    print(RULE)
    print("3. WHERE THE LINEARISATION BREAKS, AND WHETHER THE SHEAR PHASE WRAPS")
    print(RULE)
    contrasts = (0.5, 1.0, 2.0) if QUICK else (0.25, 0.5, 1.0, 2.0, 4.0)
    print("  nonlin = || d_WPM(c) - (c/eps) d_WPM(eps) || / || (c/eps) d_WPM(eps) ||: WPM against")
    print("  its OWN linearisation, so stage 1's forward-model gap cancels and what is left is")
    print("  multiple scattering plus Rytov breakdown. chen1998: first Rytov needs")
    print("  a*dn <= 0.24 lambda. The wrap columns are defined under the table.")
    print()
    print(f"  {'contrast':>8} {'dn_max':>7} {'a dn/lam':>9} {'nonlin':>7} {'|K|<10':>7} "
          f"{'|K|>20':>7} {'phi_lin':>8} {'unwrap@1px':>11} {'wrap0.8':>8} {'wrap1.6':>8}"
          f" {'self0.8':>9} {'self1.6':>9}")
    nonlin_tab = []
    for c in contrasts:
        wc = WPM(embed(nmed + c * dn1, PAD), grid.dx, grid.dz, cfg.lam, nmed,
                 n_levels=N_LEVELS, dtype=np.complex128)
        rc = crop(wc.to_centre(wc.run(ki0), ki0), grid.n)
        ac = (sc * np.fft.fft2(np.log(rc)))[ix0, iy0]
        bc = (sc * np.fft.fft2((c / EPS) * logr_eps))[ix0, iy0]
        lo, hi = Ka0 < 10, Ka0 > 20
        e = np.linalg.norm(ac - bc) / np.linalg.norm(bc)
        e_lo = np.linalg.norm(ac[lo] - bc[lo]) / np.linalg.norm(bc[lo])
        e_hi = np.linalg.norm(ac[hi] - bc[hi]) / np.linalg.norm(bc[hi])
        phi_lin = np.abs(np.imag((c / EPS) * logr_eps)).max()
        # Every wrap column is the SAME exact test on the SAME actual field: the two
        # ways of forming a shear difference disagree exactly where a 2 pi jump sits
        # between the two sampled points. At one pixel that is the residue test, i.e.
        # the unwrapping load standard Rytov ODT faces; at 0.8 and 1.6 um it is how
        # much wrapping the shearing construction sidesteps by never forming phi.
        wraps = []
        for d in (grid.dx, 0.8, 1.6):
            sw = shear_rytov(rc, grid.dx, (d, 0.0), wrap_safe=True)
            su = shear_rytov(rc, grid.dx, (d, 0.0), wrap_safe=False)
            wraps.append(np.mean(np.abs(sw - su) > 1.0))
        # Does the SHEAR DIFFERENCE itself exceed pi? That is the sampling condition
        # guizar-sicairos2011 leaves in place, and the only wrapping the
        # interferogram cannot escape. Chain one-pixel wrap-safe differences (each
        # far below pi) to get the unambiguous difference over D, and compare with
        # the direct principal value.
        # one-pixel wrap-safe steps, telescoped along x: sum_j s(r + j dx) with
        # s = log[u(r+dx)/u(r)], which sums to the UNAMBIGUOUS difference over D
        # because every term is a principal value of a one-pixel difference.
        s1 = shear_rytov(rc, grid.dx, (grid.dx, 0.0), wrap_safe=True)
        selfwrap = []
        for d in (0.8, 1.6):
            m = int(round(d / grid.dx))
            chained = sum(np.roll(s1, -j, axis=0) for j in range(m))
            direct = shear_rytov(rc, grid.dx, (m * grid.dx, 0.0), wrap_safe=True)
            selfwrap.append(np.mean(np.abs(np.imag(chained - direct)) > 1.0))
        nonlin_tab.append((c, e, e_lo, e_hi))
        print(f"  {c:8.2f} {c*dn_max:7.4f} {a_um*c*dn_max/cfg.lam:9.4f} {e:7.4f} {e_lo:7.4f} "
              f"{e_hi:7.4f} {phi_lin:8.2f} {wraps[0]:11.4f} {wraps[1]:8.4f} {wraps[2]:8.4f} "
              f"{selfwrap[0]:9.4f} {selfwrap[1]:9.4f}")
    print()
    print("  phi_lin is the linear-theory absolute Rytov phase excursion (radians), a scale.")
    print("  'unwrap@1px' is the fraction of neighbouring pixel pairs separated by a 2 pi jump")
    print("  in the principal-branch phase: exactly the residues standard Rytov ODT must unwrap")
    print("  before it can use the data, and where chen1998 says first Rytov collapses.")
    print("  'wrap 0.8' / 'wrap 1.6' apply the identical test at the working shears. They are")
    print("  NOT failures of those channels -- the interferogram hands over arg[u(r+D)u*(r)]")
    print("  directly, so that fraction is what the shearing geometry gets for free.")
    print("  'self D' IS a failure of the channel: the fraction where the shear difference")
    print("  itself exceeds pi, found by chaining one-pixel wrap-safe differences and comparing")
    print("  with the direct principal value. This is the one wrapping the interferogram cannot")
    print("  escape, and it bounds D from above independently of the design rule.")

    # ================================================================ STAGE 4
    print()
    print(RULE)
    print("4. THE exp05 DESIGN COMPARISON, RE-RUN ON EACH FORWARD MODEL (equal photon budget)")
    print(RULE)
    X, Y = (1.0, 0.0), (0.0, 1.0)
    CH = {"odt": None, "0.1x": (0.1, 0.0), "0.8x": (0.8, 0.0),
          "0.8y": (0.0, 0.8), "0.55x": (0.55, 0.0), "0.55y": (0.0, 0.55)}
    designs = [
        ("ODT benchmark", ["odt"], ["wiener"]),
        ("0.1 x", ["0.1x"], ["wiener"]),
        ("0.8 x", ["0.8x"], ["wiener", "blue"]),
        ("0.8 x + 0.8 y", ["0.8x", "0.8y"], ["wiener", "naive"]),
        ("{0.8,0.55} x + y", ["0.8x", "0.55x", "0.8y", "0.55y"], ["wiener", "blue", "naive"]),
    ]
    ms_contrasts = (1.0, 2.0) if QUICK else (0.5, 1.0, 2.0)
    results, recs = {}, {}
    logr_cache = {}          # model B' is contrast-independent: cache log u_R(eps)

    for c in ms_contrasts:
        dnc = c * dn1
        Oc = k0 ** 2 * (nmed ** 2 - (nmed + dnc) ** 2)
        Otc = np.fft.fftn(Oc)
        Sc, _, _ = shell_spectrum(Otc, Kabs, cfg.kmax, grid.dk)
        sig0 = np.sqrt(np.mean(np.abs(Otc[ref]) ** 2) / RHO)
        Ob = np.real(np.fft.ifftn(Otc * meas))
        wpm_c = WPM(embed(nmed + dnc, PAD), grid.dx, grid.dz, cfg.lam, nmed, n_levels=N_LEVELS)
        wpm_e = None if logr_cache else WPM(embed(nmed + EPS * dn1, PAD), grid.dx, grid.dz,
                                            cfg.lam, nmed, n_levels=N_LEVELS)
        # signal accumulators per (model, channel); noise and J per channel only
        numS = {(m, k): np.zeros(grid.size, np.complex64) for m in SIMULATED for k in CH}
        numN = {k: np.zeros(grid.size, np.complex64) for k in CH}
        J1 = {k: np.zeros(grid.size, np.float32) for k in CH}
        t0 = time.time()
        for j, (ki, cap) in enumerate(zip(KI, caps)):
            ixj, iyj = np.nonzero(cap["ok"])
            lin = (ixj * grid.n + iyj) * grid.nz + cap["iz"][ixj, iyj]
            phase = np.exp(-1j * cap["Kz_vox"][ixj, iyj] * cap["z_c"])
            hbase = cfg.k_m / cap["kz"][ixj, iyj]
            if j not in logr_cache:
                logr_cache[j] = np.log(crop(wpm_e.to_centre(wpm_e.run(ki), ki),
                                            grid.n)).astype(np.complex64)
            fields = {}
            r = crop(wpm_c.to_centre(wpm_c.run(ki), ki), grid.n)
            fields["C WPM"] = (np.log(r), r)
            lg = (c / EPS) * logr_cache[j].astype(np.complex128)
            fields["B' WPM-linear"] = (lg, np.exp(lg))
            rb = born_rytov_ratio(Oc, cfg, grid, ki, cap)
            fields["B Born-exact"] = (np.log(rb), rb)
            gen = np.random.default_rng((SEED, j))
            for k, dvec in CH.items():
                h = shear_factor(grid, dvec)[ixj, iyj] * hbase
                J1[k] += np.bincount(lin, weights=np.abs(h) ** 2 / sig0 ** 2,
                                     minlength=grid.size).astype(np.float32)
                nz_ = gen.standard_normal((2, lin.size))
                noise = np.sqrt(sig0 ** 2 / 2) * (nz_[0] + 1j * nz_[1])
                bin_add(numN[k], lin, np.conj(h) * noise / sig0 ** 2)
                for m in SIMULATED:
                    if m.startswith("A"):
                        y = h * Otc.reshape(-1)[lin]
                    else:
                        lg_, r_ = fields[m]
                        y = (sc * channel_meas(lg_, r_, grid, dvec))[ixj, iyj] * phase
                    bin_add(numS[(m, k)], lin, np.conj(h) * y / sig0 ** 2)
            if j == 0 or (j + 1) % 30 == 0:
                print(f"    contrast {c}: illumination {j+1}/{n_ill}  [{time.time()-t0:.0f} s]")

        # D: Born baseline + WPM's nonlinear excess. Free, because depositing is
        # linear in the data, so the combination can be formed on the accumulators.
        for k in CH:
            numS[("D Born+MS", k)] = (numS[("B Born-exact", k)]
                                      + numS[("C WPM", k)] - numS[("B' WPM-linear", k)])

        def chan(m, k, t):
            """(num, J) for channel k of a design giving it time share t."""
            num = (t * numS[(m, k)] + np.sqrt(t) * numN[k]).reshape(grid.shape)
            return hsym_c(num.astype(np.complex128)), hermitian_sym(
                (t * J1[k]).astype(np.float64).reshape(grid.shape))

        print()
        print(f"  --- contrast {c:.2f}   dn_max {c*dn_max:.4f}, a dn/lam "
              f"{a_um*c*dn_max/cfg.lam:.3f}   |   RMSE vs band-limited truth")
        print(f"  {'design':20s} {'fusion':7s} " + " ".join(f"{m:>14s}" for m in MODELS))
        for name, ks, hows in designs:
            t = 1.0 / len(ks)
            built = {m: [chan(m, k, t) for k in ks] for m in MODELS}
            for how in hows:
                row = []
                for m in MODELS:
                    est, _ = fuse(built[m], Sc, how)
                    Oe = np.real(np.fft.ifftn(est * meas))
                    r_ = np.linalg.norm(Oe - Ob) / np.linalg.norm(Ob)
                    results[(c, name, how, m)] = r_
                    if c == 1.0:
                        recs[(name, how, m)] = Oe
                    row.append(r_)
                print(f"  {name:20s} {how:7s} " +
                      " ".join(f"{v:14.4f}" if v < 100 else f"{'diverges':>14s}" for v in row))
            del built
        del numS, numN, J1, wpm_c, wpm_e

    # ------------------------------------------------------------- summary
    print()
    print(RULE)
    print("5. WHAT SURVIVES")
    print(RULE)
    key = "{0.8,0.55} x + y"
    for c in ms_contrasts:
        print(f"  contrast {c:.2f}  (a dn/lam = {a_um*c*dn_max/cfg.lam:.3f}):")
        for m in MODELS:
            f = results[(c, key, "wiener", m)]
            n_ = results[(c, key, "naive", m)]
            s_ = results[(c, "0.8 x", "wiener", m)]
            print(f"    {m:16s} fused {f:8.4f} | averaged {n_:8.4f} -> fusion wins "
                  f"{n_/f:6.1f}x | single 0.8 {s_:8.4f} -> pair wins {s_/f:6.1f}x")

    # ------------------------------------------------------------- figures
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.4))
    xs = [a_um * t[0] * dn_max / cfg.lam for t in nonlin_tab]
    for i, lbl in ((1, "all K"), (2, "|K| < 10"), (3, "|K| > 20")):
        ax[0].loglog(xs, [t[i] for t in nonlin_tab], "o-", label=lbl)
    ax[0].axvline(0.24, color="k", ls=":", label="chen1998 first-Rytov bound")
    ax[0].set(xlabel="a dn / lambda", ylabel="relative nonlinearity",
              title="WPM against its own linearisation")
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=.3, which="both")
    names = [d[0] for d in designs]
    w = 0.16
    for i, m in enumerate(MODELS):
        v = [results.get((1.0, n, "wiener", m), np.nan) for n in names]
        ax[1].bar(np.arange(len(names)) + (i - 2) * w, v, w, label=m)
    ax[1].set_yscale("log")
    ax[1].set_xticks(np.arange(len(names)))
    ax[1].set_xticklabels(names, rotation=20, ha="right", fontsize=7)
    ax[1].set(ylabel="RMSE", title="Wiener fusion at contrast 1.0")
    ax[1].legend(fontsize=7)
    ax[1].grid(alpha=.3, axis="y")
    fig.tight_layout()
    fig.savefig(out_path("ms_linearity.png"), dpi=110)

    O_band = np.real(np.fft.ifftn(np.fft.fftn(
        k0 ** 2 * (nmed ** 2 - (nmed + dn1) ** 2)) * meas))
    keys = [("truth", None, None), ("0.8 x", "wiener", "D Born+MS"),
            (key, "naive", "D Born+MS"), (key, "wiener", "A ideal-linear"),
            (key, "wiener", "B Born-exact"), (key, "wiener", "D Born+MS")]
    fig2, ax2 = plt.subplots(2, len(keys), figsize=(3.0 * len(keys), 6.2))
    cz, cy = grid.nz // 2, grid.n // 2
    vmin, vmax = O_band.min(), O_band.max()
    for j, (name, how, m) in enumerate(keys):
        vol = O_band if name == "truth" else recs.get((name, how, m))
        if vol is None:
            continue
        ax2[0, j].imshow(vol[:, :, cz].T, cmap="gray", vmin=vmin, vmax=vmax, origin="lower")
        ax2[1, j].imshow(vol[:, cy, :].T, cmap="gray", vmin=vmin, vmax=vmax, origin="lower",
                         aspect=grid.dz / grid.dx)
        ax2[0, j].set_title("band-limited truth" if name == "truth" else f"{name}\n{how} / {m}",
                            fontsize=7)
        for aa in ax2[:, j]:
            aa.axis("off")
    fig2.suptitle("Reconstructions at contrast 1.0, equal photon budget", fontsize=10)
    fig2.tight_layout()
    fig2.savefig(out_path("ms_slices.png"), dpi=110)
    print()
    print(f"  wrote {out_path('ms_linearity.png')}")
    print(f"  wrote {out_path('ms_slices.png')}")
    print(f"  total {time.time()-t_all:.0f} s")


if __name__ == "__main__":
    main()
