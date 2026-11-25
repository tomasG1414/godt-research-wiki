"""exp13 -- Partial coherence and the shear factor.

Question: with a spatially partially coherent (extended, incoherent) source,
is the phase-difference datum of a common-path lateral shearing interferometer
still  G(K_perp) x (coherent Rytov phase spectrum),  G = e^{i K_perp.D} - 1,
with the coherence entering only as a scalar visibility |mu(D)|?

Derivation (wiki note, section 2).  Source point s = plane wave with transverse
wavevector q_s = q0 + d_s, weight w_s.  Five-step phase shifting returns the
fringe phasor, which is the mutual intensity of the image field at the shear:

    Gamma(x) = sum_s w_s u_s(x + D) u_s*(x),     Gamma_0 = e^{i q0 D} mu(D),
    mu(D)    = sum_s w_s e^{i d_s D}             (van Cittert-Zernike).

For a weak object, u_s = e^{i q_s x} (1 + psi_s), and with c_s = e^{i d_s D}:

    dphi(x) = arg(Gamma / Gamma_0)
            = Im { (1/mu) sum_s w_s c_s [ psi_s(x + D) + psi_s*(x) ] }      (M2)

psi_s lives on source point s's OWN Ewald cap, K_z = k_z(kappa + q_s) - k_z(q_s).
Linearising K_z in d_s gives, per object slice at depth z, a generalised shear
factor in closed form (M3):

    G_pc(kappa; z) = [ e^{i kappa D} mu(D - z v) - mu(D + z v) ] / mu(D),
    v(kappa) = tan(theta_i) - tan(theta_s(kappa))      (ray-slope difference)

which is G(kappa) for z = 0 or a point source, and is NOT zero on G's nulls.

Numerical check, 2D (x, z), GODT simulation optics (0.5 um, n0 = 1.33, 60 deg):
N_s mutually incoherent tilted plane waves, uniform in q over a cone of
half-angle alpha about 60 deg; each is propagated through the object by a
multi-slice Born scheme (forward multiple scattering, exact first Born in its
linear limit -- checked below); the sheared copy is added with five phase
steps, INTENSITIES are summed over the source, and the five-step formula
recovers the phase difference.  Compared against

    M0  G x coherent phase of the nominal source point (the pipeline's model)
    M1  G x unweighted source average of the per-source phases (real weights)
    M2  the derived linear model above (complex weights c_s / mu)
    M3  the closed form with mu(D -+ z v)

Part 6 repeats this for a broadband source (grating shear D ~ lambda).
Parts 2 and 6 also report errors with the pupil-edge band removed: an extended
source (or a spectrum) blurs a hard pupil edge over ~2h, which dominates the
full-band numbers and is not a property of the shear factor. (The interior
analysis and a units fix in the part-0 check were added in review.)

Stated models, none sourced: the source is a uniform cone about the nominal
direction; the phase shift is achromatic; the empty-field fringe Gamma_0 is
measured and divided out; the phantoms.

Writes: out/partial_coherence.png
Wiki:   wiki/notes/partial-coherence-and-the-shear-factor.md

Run: python exp13_partial_coherence.py     (~1 min; numpy + scipy + matplotlib)
"""
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import j1

from shearlib import GODT_SIM, at_wavelength, chromatic_G, visibility, out_path

RULE = "=" * 78
N, DX = 2048, 0.05                # transverse grid; L = 102.4 um (0.8 um nulls on-grid)
NZ, DZ = 65, 0.1                  # slices; z in [-3.2, 3.2] um, slice 32 at z = 0
NS = 41                           # source points (odd: the nominal one is included)
SHEARS = (0.2, 0.55, 0.8)         # um
ALPHAS = (0.25, 0.5, 1.0, 2.0, 4.0, 8.0)   # cone half-angle in the medium, deg
WEAK = 0.01                       # contrast factor of the "weak" thick phantom

X = (np.arange(N) - N // 2) * DX
Z = (np.arange(NZ) - NZ // 2) * DZ
KAP = 2 * np.pi * np.fft.fftfreq(N, DX)          # kappa = k_x - q_s (envelope frequency)


# ------------------------------------------------------------------ phantoms
def _smooth(a, sig_x, sig_z=None):
    """Periodic Gaussian low-pass (FFT), a has shape (..., N) or (NZ, N)."""
    F = np.fft.fft(a, axis=-1) * np.exp(-0.5 * (KAP * sig_x) ** 2)
    a = np.fft.ifft(F, axis=-1).real
    if sig_z is not None:
        kz = 2 * np.pi * np.fft.fftfreq(a.shape[0], DZ)
        a = np.fft.ifft(np.fft.fft(a, axis=0) * np.exp(-0.5 * (kz * sig_z) ** 2)[:, None],
                        axis=0).real
    return a


def thin_profile(rng):
    """1D index profile dn(x) of a thin (one-slice, 0.1 um) phase object."""
    win = 0.5 * (1 - np.tanh((np.abs(X) - 10.0) / 0.5))
    tex = _smooth(rng.standard_normal(N), 0.12)
    tex /= tex.std()
    f = win * (0.6 + 0.4 * tex)
    for x0, s, a in ((-4.0, 0.2, 1.0), (1.5, 0.12, -0.8), (5.5, 0.4, 0.7)):
        f += a * np.exp(-0.5 * ((X - x0) / s) ** 2)
    return 0.01 * f / np.abs(f).max()


def thin_object(prof, z0):
    dn = np.zeros((NZ, N))
    dn[int(np.argmin(np.abs(Z - z0)))] = prof
    return dn


def thick_object(rng):
    """Cell-like 2D index map dn(x, z): 10 x 4.8 um body, nucleus, granules, texture."""
    XX, ZZ = np.meshgrid(X, Z)
    s = np.hypot(XX / 5.0, ZZ / 2.4)
    body = 0.5 * (1 - np.tanh((s - 1) / 0.06))
    dn = 0.010 * body
    sn = np.hypot((XX + 1.0) / 1.6, (ZZ - 0.3) / 1.2)
    dn += 0.008 * 0.5 * (1 - np.tanh((sn - 1) / 0.08))
    for _ in range(8):
        x0, z0 = rng.uniform(-3.5, 3.5), rng.uniform(-1.6, 1.6)
        dn += rng.choice([-1, 1]) * 0.02 * np.exp(-0.5 * ((XX - x0) ** 2 + (ZZ - z0) ** 2) / 0.15 ** 2)
    tex = _smooth(rng.standard_normal((NZ, N)), 0.12, 0.12)
    dn += 0.003 * tex / tex.std() * body
    return dn


# ------------------------------------------------------- fields per source point
def geom(cfg, q):
    """Collected-pupil mask, scattered k_z and K_z on the envelope grid for tilt q."""
    kx = KAP + q
    ok = np.abs(kx) <= cfg.k_perp_max
    kz = np.sqrt(np.maximum(cfg.k_m ** 2 - kx ** 2, 0.0))
    kzi = np.sqrt(cfg.k_m ** 2 - q ** 2)
    return ok, np.where(ok, kz, 1.0), kzi


def potential(cfg, dn):
    """O = k0^2 (n^2 - n0^2), Helmholtz (lap + k_m^2) u = -O u."""
    k0 = 2 * np.pi / cfg.lam
    return k0 ** 2 * ((cfg.n0 + dn) ** 2 - cfg.n0 ** 2)


def msb(cfg, O, q):
    """Multi-slice Born: envelope a(x) = u e^{-i q x} refocused to the plane z = 0.

    Per slice: half step, add the slab's first-Born scattering of the CURRENT
    field with the outgoing obliquity i dz / (2 k_z), half step. Components
    outside the collection pupil are dropped at every step. Forward multiple
    scattering only; its first-order term is `born` exactly.
    """
    ok, kz, kzi = geom(cfg, q)
    half = np.where(ok, np.exp(0.5j * kz * DZ), 0.0)
    scat = np.where(ok, 1j * DZ / (2 * kz), 0.0)
    z_top, z_bot = Z[0] - DZ / 2, Z[-1] + DZ / 2
    A = np.zeros(N, complex)
    A[0] = N * np.exp(1j * kzi * z_top)
    for j in range(NZ):
        A = A * half
        if O[j].any():
            A = A + scat * np.fft.fft(O[j] * np.fft.ifft(A))
        A = A * half
    return np.where(ok, A * np.exp(-1j * kz * z_bot), 0.0)     # spectrum of a(x) at z = 0


def born(cfg, Ot, q):
    """First-Born psi~_s(kappa) = i/(2 k_z) sum_z dz O~(kappa; z) e^{-i K_z z}."""
    ok, kz, kzi = geom(cfg, q)
    ph = np.exp(-1j * np.outer(Z, kz - kzi))                   # (NZ, N)
    return np.where(ok, 1j * DZ / (2 * kz) * (Ot * ph).sum(axis=0), 0.0)


def source(cfg, alpha_deg, ns=NS):
    """Uniform slit source in q about the nominal tilt: offsets d_s, weights."""
    h = cfg.k_m * np.cos(np.deg2rad(cfg.theta_ill)) * np.sin(np.deg2rad(alpha_deg))
    if h == 0:
        return np.zeros(1), np.ones(1), 0.0
    d = h * (2 * (np.arange(ns) + 0.5) / ns - 1)               # midpoints; d = 0 included
    return d, np.full(ns, 1.0 / ns), h


def mu_slit(r, h):
    """Continuum degree of coherence of the uniform slit source: sin(h r)/(h r)."""
    return np.sinc(h * np.asarray(r) / np.pi)


# ----------------------------------------------------- the instrument (5 frames)
def five_step(frames):
    """Fringe phasor from five frames at phase steps (n - 3) pi/2 (Hariharan)."""
    I1, I2, I3, I4, I5 = frames
    return ((2 * I3 - I1 - I5) + 2j * (I2 - I4)) / 8.0


def interferogram_frames(specs, piston, shear_of):
    """Sum over mutually incoherent components of the five phase-stepped
    sheared interferograms. `specs` = list of (weight, envelope spectrum A,
    key) and shear_of(key) the shear that component sees."""
    frames = [np.zeros(N) for _ in range(5)]
    for w, A, key in specs:
        D = shear_of(key)
        a = np.fft.ifft(A)
        ash = np.fft.ifft(A * np.exp(1j * KAP * D)) * piston(key, D)
        for n in range(5):
            frames[n] += w * np.abs(ash * np.exp(1j * (n - 2) * np.pi / 2) + a) ** 2
    return frames


def measure(specs_obj, specs_empty, piston, shear_of):
    """Phase difference as the instrument delivers it: five-step phasor of the
    object frames, referenced to the empty-field phasor. Also the visibility."""
    g = five_step(interferogram_frames(specs_obj, piston, shear_of))
    f0 = interferogram_frames(specs_empty, piston, shear_of)
    g0 = five_step(f0)
    vis = 2 * np.abs(g0).mean() / np.mean(f0[2] + f0[0]) * 2   # 2|Gamma0| / mean intensity
    return np.angle(g * np.conj(g0)), vis, g0.mean()


def shift(F, D):
    return np.fft.ifft(F * np.exp(1j * KAP * D))


def rel(a, b):
    return 100.0 * np.linalg.norm(a - b) / np.linalg.norm(a)


# ------------------------------------------------------------------- models
def model_set(cfg, Ot, d, w, h, D):
    """M0..M3 phase differences (real arrays) for a monochromatic slit source."""
    q0 = cfg.k_ill
    c = np.exp(1j * d * D)
    mu = (w * c).sum()
    psis = [born(cfg, Ot, q0 + ds) for ds in d]
    pair = [shift(P, D) + np.conj(np.fft.ifft(P)) for P in psis]      # psi(x+D) + psi*(x)
    i0 = int(np.argmin(np.abs(d)))
    M0 = pair[i0].imag
    M1 = sum(wi * p.imag for wi, p in zip(w, pair))
    M2 = (sum(wi * ci * p for wi, ci, p in zip(w, c, pair)) / mu).imag
    M1c = sum(wi * ci.real * p.imag for wi, ci, p in zip(w, c, pair)) / mu.real
    # ---- M3: closed form, K_z linearised in the source offset
    ok, kz, kzi = geom(cfg, q0)
    v = np.where(ok, q0 / kzi - (KAP + q0) / kz, 0.0)                  # dK_z/dq
    Psi = np.where(ok, 1j * DZ / (2 * kz), 0.0) * Ot * np.exp(-1j * np.outer(Z, kz - kzi))
    zv = np.outer(Z, v)
    m0 = mu_slit(D, h)
    t1 = np.exp(1j * KAP * D) * (mu_slit(D - zv, h) * Psi).sum(axis=0) / m0
    s2 = (mu_slit(D + zv, h) * Psi).sum(axis=0) / m0
    M3 = (np.fft.ifft(t1) + np.conj(np.fft.ifft(s2))).imag
    return {"M0": M0, "M1": M1, "M2": M2, "M3": M3, "M1c": M1c, "mu": mu}


def run_case(cfg, dn, alpha, shears=SHEARS):
    """Simulate + model one object under one cone. Returns per-shear dicts."""
    O = potential(cfg, dn)
    Ot = np.fft.fft(O, axis=1)
    d, w, h = source(cfg, alpha)
    q0 = cfg.k_ill
    A = [msb(cfg, O, q0 + ds) for ds in d]
    E = [msb(cfg, np.zeros_like(O), q0 + ds) for ds in d]
    specs = [(wi, Ai, ds) for wi, Ai, ds in zip(w, A, d)]
    empty = [(wi, Ei, ds) for wi, Ei, ds in zip(w, E, d)]
    out = {}
    for D in shears:
        piston = lambda ds, DD: np.exp(1j * (q0 + ds) * DD)
        dphi, vis, g0 = measure(specs, empty, piston, lambda key: D)
        m = model_set(cfg, Ot, d, w, h, D)
        m.update(sim=dphi, vis=vis, h=h, mu_cont=float(mu_slit(D, h)))
        out[D] = m
    return out


# ============================================================================
def main():
    t_start = time.time()
    cfg = GODT_SIM
    rng = np.random.default_rng(13)
    prof = thin_profile(rng)
    thick = thick_object(rng)
    objects = {
        "thin, z = 0": thin_object(prof, 0.0),
        "thin, z = +2 um": thin_object(prof, 2.0),
        "thick (weak)": WEAK * thick,
    }
    q0, km = cfg.k_ill, cfg.k_m
    print(RULE)
    print("exp13 -- partial coherence and the shear factor (2D x-z, GODT simulation optics)")
    print(RULE)
    print(f"lambda = {cfg.lam} um, n0 = {cfg.n0}, k_m = {km:.3f}, q0 = k_m sin 60 = {q0:.3f}, "
          f"k_perp_max = {cfg.k_perp_max:.3f} rad/um")
    print(f"grid: N = {N}, dx = {DX} um (L = {N * DX:.1f} um), {NZ} slices of {DZ} um; "
          f"{NS} source points")
    print(f"envelope band kappa in [{-cfg.k_perp_max - q0:.2f}, {cfg.k_perp_max - q0:.2f}] rad/um")
    k0 = 2 * np.pi / cfg.lam
    print(f"thin object: peak dn = {np.abs(prof).max():.3f} over {DZ} um -> peak phase "
          f"{k0 * np.abs(prof).max() * DZ:.4f} rad (normal incidence)")
    print(f"thick object: peak dn = {np.abs(thick).max():.4f}, body 10 x 4.8 um; projected phase "
          f"(normal incidence) peak {k0 * np.abs(thick.sum(axis=0)).max() * DZ:.3f} rad at contrast 1, "
          f"{WEAK * k0 * np.abs(thick.sum(axis=0)).max() * DZ:.4f} rad at contrast {WEAK}")

    # ---------------------------------------------------------------- part 0
    print("\n" + RULE)
    print("0. Forward-model checks")
    print(RULE)
    O = potential(cfg, 1e-4 * thick)
    A = msb(cfg, O, q0)
    psi_msb = A.copy()               # both spectra in numpy FFT units: remove the
    psi_msb[0] -= N                  # incident wave (N at kappa = 0), do not divide by N
    psi_b = born(cfg, np.fft.fft(O, axis=1), q0)
    e = np.linalg.norm(psi_msb - psi_b) / np.linalg.norm(psi_b)
    print(f"MSB -> first Born as contrast -> 0 (contrast 1e-4, nominal tilt): "
          f"|psi_MSB - psi_Born| / |psi_Born| = {e:.2e}")
    # coherent sanity: one source point, weak object -> data = G x phase (M0)
    r = run_case(cfg, objects["thick (weak)"], 0.0)
    for D in SHEARS:
        print(f"coherent source, thick (weak), shear {D:4.2f}: |sim - M0| / |sim| = "
              f"{rel(r[D]['sim'], r[D]['M0']):.3f} %   (linearisation floor)")
    r2 = run_case(cfg, 0.1 * objects["thick (weak)"], 0.0, shears=(0.8,))
    print(f"   same at 10x lower contrast, shear 0.80: {rel(r2[0.8]['sim'], r2[0.8]['M0']):.4f} %  "
          f"(floor scales with contrast -> it is the weak-object linearisation)")

    # ---------------------------------------------------------------- part 1
    print("\n" + RULE)
    print("1. Fringe visibility |mu(D)|: simulated (five-step, empty field) vs van Cittert-Zernike")
    print(RULE)
    print("2D slit source, uniform in q over q0 +- h, h = k_m cos(60) sin(alpha):  mu = sin(hD)/(hD)")
    print(f"{'alpha':>6} {'NA_c':>6} {'h':>7} {'1st zero':>9} | " +
          " | ".join(f"D={D:4.2f}: sim   sinc " for D in SHEARS))
    vis_tab = {}
    empty_dn = np.zeros((NZ, N))
    for al in ALPHAS:
        d, w, h = source(cfg, al)
        E = [(wi, msb(cfg, empty_dn, q0 + ds), ds) for wi, ds in zip(w, d)]
        row = []
        for D in SHEARS:
            f0 = interferogram_frames(E, lambda ds, DD: np.exp(1j * (q0 + ds) * DD), lambda k: D)
            g0 = five_step(f0)
            vis = 2 * np.abs(g0).mean() / (0.5 * np.mean(f0[0] + f0[2]))
            vis_tab[(al, D)] = (vis, float(mu_slit(D, h)))
            row.append(f"       {vis:.4f} {mu_slit(D, h):.4f}")
        print(f"{al:6.2f} {cfg.n0 * np.sin(np.deg2rad(al)):6.3f} {h:7.4f} {np.pi / h:9.2f} | " + " | ".join(row))
    print("(1st zero in um. Fringe MODULATION is |mu|; the Fisher weight is |mu|^2.)")

    print("\n3D: circular cone of half-angle alpha about the 60 deg direction (direct sum over")
    print("directions, uniform in solid angle), shear along x (in the plane of incidence) and")
    print("along y (perpendicular). Small-cone formula: mu = 2 J1(r)/r,")
    print("r = k_m sin(alpha) sqrt((D_par cos 60)^2 + D_perp^2); shearlib.visibility is that")
    print("jinc with r = pi D / lc, i.e. lc = lambda / (2 n0 sin alpha), and NO cos(60).")
    th0 = np.deg2rad(cfg.theta_ill)
    print(f"{'alpha':>6} {'lc':>7} {'D':>5} | {'|mu| x: sum':>11} {'formula':>8} | {'|mu| y: sum':>11} "
          f"{'formula':>8} | {'shearlib |mu|':>13} | {'|mu|^2: x':>9} {'y':>6} {'shearlib':>8}")
    for al in ALPHAS:
        a = np.deg2rad(al)
        nb, ng = 200, 256
        beta = (np.arange(nb) + 0.5) / nb * a
        gam = (np.arange(ng) + 0.5) / ng * 2 * np.pi
        B, Gm = np.meshgrid(beta, gam, indexing="ij")
        wgt = np.sin(B)
        wgt /= wgt.sum()
        sx, sy, sz = np.sin(B) * np.cos(Gm), np.sin(B) * np.sin(Gm), np.cos(B)
        qx = km * (sx * np.cos(th0) + sz * np.sin(th0))        # rotate the cone axis to 60 deg
        qy = km * sy
        lc = cfg.lam / (2 * cfg.n0 * np.sin(a))
        for D in SHEARS:
            mx = abs((wgt * np.exp(1j * qx * D)).sum())
            my = abs((wgt * np.exp(1j * qy * D)).sum())
            rx, ry = km * np.sin(a) * D * np.cos(th0), km * np.sin(a) * D
            fx, fy = abs(2 * j1(rx) / rx), abs(2 * j1(ry) / ry)
            sl = float(np.sqrt(visibility(D, lc)))
            print(f"{al:6.2f} {lc:7.2f} {D:5.2f} | {mx:11.4f} {fx:8.4f} | {my:11.4f} {fy:8.4f} | "
                  f"{sl:13.4f} | {mx ** 2:9.4f} {my ** 2:6.4f} {sl ** 2:8.4f}")

    # ---------------------------------------------------------------- part 2
    print("\n" + RULE)
    print("2. Phase-difference data vs models: relative L2 error, % of |sim|")
    print(RULE)
    print("M0 = G x coherent phase (nominal source point)   M1 = G x unweighted source average")
    print("M2 = derived linear model (weights c_s / mu)      M3 = closed form mu(D -+ z v)")
    res = {}
    for name, dn in objects.items():
        print(f"\n--- {name}")
        print(f"{'alpha':>6} | " + " | ".join(f"D = {D:4.2f}:   M0     M1     M2     M3" for D in SHEARS))
        for al in ALPHAS:
            r = run_case(cfg, dn, al)
            res[(name, al)] = r
            print(f"{al:6.2f} | " + " | ".join(
                "          " + " ".join(f"{rel(r[D]['sim'], r[D][m]):6.2f}" for m in ("M0", "M1", "M2", "M3"))
                for D in SHEARS))

    # the same comparison away from the pupil edge
    print("\nThe same at D = 0.8 um, with a band of width h + 2 bins at each pupil edge removed")
    print("from both images (an extended source blurs a hard pupil edge over ~2h; this")
    print("separates that from the shear factor proper):")
    dkk = 2 * np.pi / (N * DX)
    print(f"{'object':>18} {'alpha':>6} |    M0     M1     M2   | full-band M0")
    for name in objects:
        for al in ALPHAS:
            r = res[(name, al)][0.8]
            m = r["h"] + 2 * dkk
            edge = (np.abs(np.abs(KAP + q0) - cfg.k_perp_max) < m) | \
                   (np.abs(np.abs(-KAP + q0) - cfg.k_perp_max) < m)
            Fs = np.fft.fft(r["sim"]) * ~edge
            errs = [100 * np.linalg.norm(Fs - np.fft.fft(r[k]) * ~edge) / np.linalg.norm(Fs)
                    for k in ("M0", "M1", "M2")]
            print(f"{name:>18} {al:6.2f} | " + " ".join(f"{v:6.2f}" for v in errs)
                  + f"   | {rel(r['sim'], r['M0']):6.2f}")

    # amplitude -> phase leak
    print("\nAmplitude-to-phase leak: |M2 - M1c| / |M2|, M1c = cos(d_s D)/mu-weighted PHASE differences only")
    print(f"{'object':>18} {'alpha':>6} | " + "  ".join(f"D={D:4.2f}" for D in SHEARS))
    for name in objects:
        for al in (1.0, 4.0, 8.0):
            r = res[(name, al)]
            print(f"{name:>18} {al:6.2f} | " + "  ".join(
                f"{100 * np.linalg.norm(r[D]['M2'] - r[D]['M1c']) / np.linalg.norm(r[D]['M2']):6.2f}%"
                for D in SHEARS))

    # ---------------------------------------------------------------- part 3
    print("\n" + RULE)
    print("3. The nulls of G at 0.8 um shear (kappa = 7.854, 15.708, 23.562 rad/um; on-grid)")
    print(RULE)
    print("|F[dphi]| on the three null bins, relative to the rms |F[dphi]| over the band")
    print("0.5 < |kappa| < 31 rad/um.  M0 is zero there by construction.")
    dk = 2 * np.pi / (N * DX)
    nb_idx = [int(round(m * 2 * np.pi / 0.8 / dk)) for m in (1, 2, 3)]
    band = (np.abs(KAP) > 0.5) & (np.abs(KAP) < 31.0)
    print(f"{'object':>18} {'alpha':>6} | {'sim: m=1':>9} {'m=2':>7} {'m=3':>7} | {'M2: m=1':>9} {'m=2':>7} "
          f"{'m=3':>7} | {'M0 max':>8}")
    null_tab = {}
    for name in objects:
        for al in ALPHAS:
            r = res[(name, al)][0.8]
            Fs, F2, F0 = (np.fft.fft(r[k]) for k in ("sim", "M2", "M0"))
            rms = np.sqrt(np.mean(np.abs(Fs[band]) ** 2))
            vs = [abs(Fs[i]) / rms for i in nb_idx]
            v2 = [abs(F2[i]) / rms for i in nb_idx]
            null_tab[(name, al)] = vs
            print(f"{name:>18} {al:6.2f} | {vs[0]:9.4f} {vs[1]:7.4f} {vs[2]:7.4f} | {v2[0]:9.4f} {v2[1]:7.4f} "
                  f"{v2[2]:7.4f} | {max(abs(F0[i]) for i in nb_idx) / rms:8.1e}")

    # closed-form generalised shear factor at the nulls, for orientation
    ok, kz, kzi = geom(cfg, q0)
    v = np.where(ok, q0 / kzi - (KAP + q0) / kz, np.nan)
    print("\nClosed form |G_pc(kappa; z)| = |e^{i kappa D} mu(D - z v) - mu(D + z v)| / |mu(D)| at D = 0.8 um,")
    print("on the negative-kappa side of each null (the side inside the envelope band); |G| <= 2.")
    print(f"{'kappa':>8} {'v':>6} | " + " ".join(f"z={z:3.1f},a={al:3.1f}" for z in (1.0, 2.0) for al in (1.0, 4.0)))
    for m in (1, 2, 3):
        i = (-nb_idx[m - 1]) % N
        row = []
        for z in (1.0, 2.0):
            for al in (1.0, 4.0):
                h = source(cfg, al)[2]
                g = (np.exp(1j * KAP[i] * 0.8) * mu_slit(0.8 - z * v[i], h) - mu_slit(0.8 + z * v[i], h)) / mu_slit(0.8, h)
                row.append(f"{abs(g):11.4f}")
        print(f"{KAP[i]:8.3f} {v[i]:6.3f} | " + " ".join(row))

    # depth of coherence: where the walk-off z v reaches the coherence scale
    print("\nValidity of the separable model: |z| v h << 1, i.e. |z| << 1/(h |v|).  rms |v| over the")
    vv = v[ok & (np.abs(KAP) > 0.5)]
    vr = float(np.sqrt(np.nanmean(vv ** 2)))
    print(f"collected band = {vr:.2f} (median |v| = {np.nanmedian(np.abs(vv)):.2f}, max {np.nanmax(np.abs(vv)):.1f}).")
    print(f"{'alpha':>6} {'h':>7} {'1/(h v_med)':>12} {'1/(h v_rms)':>12}   (um)")
    for al in ALPHAS:
        h = source(cfg, al)[2]
        print(f"{al:6.2f} {h:7.4f} {1 / (h * np.nanmedian(np.abs(vv))):12.2f} {1 / (h * vr):12.2f}")

    # ---------------------------------------------------------------- part 4
    print("\n" + RULE)
    print("4. Full contrast (thick phantom, contrast 1): does the weak-object conclusion survive?")
    print(RULE)
    print("coh = the SIMULATED coherent measurement (nominal source point, same MSB, same")
    print("five-step); M2 = derived linear model. Errors in % of |sim|.")
    print(f"{'alpha':>6} | " + " | ".join(f"D = {D:4.2f}:  coh     M0     M2" for D in SHEARS)
          + " |  max |dphi| (rad) at 0.8")
    rc = run_case(cfg, thick, 0.0)
    full = {}
    for al in ALPHAS:
        r = run_case(cfg, thick, al)
        full[al] = r
        print(f"{al:6.2f} | " + " | ".join(
            "         " + " ".join(f"{rel(r[D]['sim'], x):6.2f}" for x in (rc[D]['sim'], r[D]['M0'], r[D]['M2']))
            for D in SHEARS) + f" |  {np.abs(r[0.8]['sim']).max():.3f}")
    print("same, weak thick phantom (contrast %.2g), coh column only:" % WEAK)
    rcw = run_case(cfg, objects["thick (weak)"], 0.0)
    for al in ALPHAS:
        r = res[("thick (weak)", al)]
        print(f"{al:6.2f} | " + " | ".join(f"D = {D:4.2f}: {rel(r[D]['sim'], rcw[D]['sim']):6.2f}" for D in SHEARS))

    # ---------------------------------------------------------------- part 5
    print("\n" + RULE)
    print("5. Where in the band the coherent model fails (thick weak, D = 0.8 um)")
    print(RULE)
    edges = [0.5, 4, 8, 12, 16, 20, 24, 28, 31]
    print("ratio rms|F[sim]| / rms|F[M0]| and error |F[sim] - F[M0]| / |F[sim]| per |kappa| bin")
    print(f"{'|kappa| bin':>12} | " + " | ".join(f"a={al:4.2f}: ratio  err%" for al in (1.0, 2.0, 4.0, 8.0)))
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (np.abs(KAP) >= lo) & (np.abs(KAP) < hi)
        row = []
        for al in (1.0, 2.0, 4.0, 8.0):
            r = res[("thick (weak)", al)][0.8]
            Fs, F0 = np.fft.fft(r["sim"]), np.fft.fft(r["M0"])
            row.append(f"        {np.linalg.norm(Fs[sel]) / np.linalg.norm(F0[sel]):6.3f} "
                       f"{100 * np.linalg.norm(Fs[sel] - F0[sel]) / np.linalg.norm(Fs[sel]):5.1f}")
        print(f"{lo:5.1f}-{hi:5.1f}  | " + " | ".join(row))

    # ---------------------------------------------------------------- part 6
    print("\n" + RULE)
    print("6. Broadband source (secondary): Gaussian spectrum, grating shear D(lam) = D0 lam/lam0")
    print(RULE)
    lam0 = cfg.lam
    fwhm = lam0 * 35.0 / 530.0
    sig = fwhm / (2 * np.sqrt(2 * np.log(2)))
    nl = 15
    lams = lam0 + sig * np.linspace(-2.5, 2.5, nl)
    sw = np.exp(-0.5 * ((lams - lam0) / sig) ** 2)
    sw /= sw.sum()
    print(f"lam0 = {lam0} um, FWHM = {fwhm * 1e3:.1f} nm (GODT's 35/530 fractional bandwidth), "
          f"{nl} wavelengths over +-2.5 sigma")
    print("Geometry per wavelength: angles and NA fixed (shearlib.at_wavelength), no dispersion,")
    print("O = k0(lam)^2 (n^2 - n0^2). Phase steps taken as achromatic (geometric phase).")
    print("\nEmpty-field fringe visibility under the spectrum alone (point source):")
    for D in SHEARS:
        ph_gr = np.array([at_wavelength(cfg, l).k_ill * D * l / lam0 for l in lams])
        ph_fx = np.array([at_wavelength(cfg, l).k_ill * D for l in lams])
        print(f"  D0 = {D:4.2f}: grating shear (D ~ lam) |sum s e^(i q D)| = {abs((sw * np.exp(1j * ph_gr)).sum()):.6f}"
              f"   wavelength-independent shear: {abs((sw * np.exp(1j * ph_fx)).sum()):.4f}"
              f"   (piston q0 D0 = {cfg.k_ill * D:.2f} rad)")

    def broadband(dn, alpha, D0):
        specs, empty, pairs = [], [], []
        norm_mu = 0.0
        for l, s in zip(lams, sw):
            c = at_wavelength(cfg, l)
            O = potential(c, dn)
            Ot = np.fft.fft(O, axis=1)
            d, w, h = source(c, alpha, ns=11 if alpha > 0 else 1)
            Dl = D0 * l / lam0
            for ds, wi in zip(d, w):
                key = (c.k_ill + ds, Dl)
                specs.append((s * wi, msb(c, O, c.k_ill + ds), key))
                empty.append((s * wi, msb(c, np.zeros_like(O), c.k_ill + ds), key))
                P = born(c, Ot, c.k_ill + ds)
                cs = np.exp(1j * ds * Dl)
                pairs.append(s * wi * cs * (shift(P, Dl) + np.conj(np.fft.ifft(P))))
                norm_mu += s * wi * cs
        # common piston e^{i q0(lam) D(lam)} is wavelength-free for grating shear; keep it exact anyway
        dphi, vis, _ = measure(specs, empty, lambda key, DD: np.exp(1j * key[0] * key[1]), lambda key: key[1])
        M2 = (sum(pairs) / norm_mu).imag
        # monochromatic nominal (M0) and shearlib's chromatic_G on the nominal cap (Mc)
        Ot0 = np.fft.fft(potential(cfg, dn), axis=1)
        P0 = born(cfg, Ot0, cfg.k_ill)
        M0 = (shift(P0, D0) + np.conj(np.fft.ifft(P0))).imag
        phi0 = np.fft.fft(np.fft.ifft(P0).imag)
        Mc = np.fft.ifft(chromatic_G(KAP, D0, lam0, fwhm)[0] * phi0).real
        return dphi, M0, Mc, M2, vis

    print("\nErrors in % of |sim|.  M0 = monochromatic G x coherent phase at lam0;")
    print("Mc = shearlib.chromatic_G (spectral average of G) x the lam0 phase on the lam0 cap;")
    print("M2 = derived linear model, each wavelength on its own cap with its own shear.")
    print("'interior': the band swept by the pupil edges over the spectrum (+ source spread")
    print("+ 2 bins) removed from both images, as in part 2.")
    dkk = 2 * np.pi / (N * DX)
    hi = [at_wavelength(cfg, l).k_perp_max - at_wavelength(cfg, l).k_ill for l in lams]
    lo = [-at_wavelength(cfg, l).k_perp_max - at_wavelength(cfg, l).k_ill for l in lams]

    def interior(al):
        pad = source(cfg, al)[2] + 2 * dkk
        keep = np.ones(N, bool)
        for ed in (hi, lo):
            for sgn in (1, -1):
                keep &= ~((sgn * KAP >= min(ed) - pad) & (sgn * KAP <= max(ed) + pad))
        return keep

    def rel_in(a, b, keep):
        A = np.fft.fft(a) * keep
        return 100 * np.linalg.norm(A - np.fft.fft(b) * keep) / np.linalg.norm(A)

    print(f"{'object':>18} {'alpha':>6} {'D0':>5} | {'M0':>6} {'Mc':>6} {'M2':>6} | interior: "
          f"{'M0':>6} {'Mc':>6} {'M2':>6} | visibility")
    bb = {}
    for name in ("thin, z = 0", "thin, z = +2 um", "thick (weak)"):
        for al in (0.0, 2.0):
            keep = interior(al)
            for D in SHEARS:
                dphi, M0, Mc, M2, vis = broadband(objects[name], al, D)
                bb[(name, al, D)] = (rel(dphi, M0), rel(dphi, Mc), rel(dphi, M2))
                print(f"{name:>18} {al:6.2f} {D:5.2f} | {rel(dphi, M0):6.2f} {rel(dphi, Mc):6.2f} "
                      f"{rel(dphi, M2):6.2f} | interior: {rel_in(dphi, M0, keep):6.2f} "
                      f"{rel_in(dphi, Mc, keep):6.2f} {rel_in(dphi, M2, keep):6.2f} | {vis:.4f}")

    # ------------------------------------------------------------------ figure
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    dd = np.linspace(1e-3, 3.0, 300)
    for al in ALPHAS:
        h = source(cfg, al)[2]
        ln, = ax[0].plot(dd, np.abs(mu_slit(dd, h)), lw=1.2, label=f"{al}°")
        ax[0].plot(SHEARS, [vis_tab[(al, D)][0] for D in SHEARS], "o", color=ln.get_color(), ms=4)
    ax[0].set(xlabel="shear Δ (µm)", ylabel="|µ(Δ)|", title="Fringe visibility: slit source about 60°\n"
              "lines sin(hΔ)/(hΔ), dots five-step simulation")
    ax[0].legend(title="cone half-angle", fontsize=8)
    r = res[("thick (weak)", 4.0)][0.8]
    ks = np.fft.fftshift(KAP)
    for key, lab, st in (("sim", "simulated (5-step, intensities summed)", dict(color="k", lw=1.0)),
                         ("M0", "M0: G × coherent phase", dict(color="C3", lw=0.8, alpha=0.8)),
                         ("M2", "M2: derived (per-source caps)", dict(color="C0", lw=0.8, ls="--"))):
        ax[1].semilogy(ks, np.abs(np.fft.fftshift(np.fft.fft(r[key]))) + 1e-12, label=lab, **st)
    for m in (1, 2, 3):
        ax[1].axvline(-2 * np.pi * m / 0.8, color="0.7", lw=0.6)
    ax[1].set(xlim=(-32, 0), xlabel="κ (rad/µm)", ylabel="|F[δφ]|",
              title="Thick weak object, Δ = 0.8 µm, cone 4°\n(grey: nulls of G)")
    ax[1].legend(fontsize=8)
    for D, col in zip(SHEARS, ("C0", "C1", "C2")):
        ax[2].loglog(ALPHAS, [rel(res[("thick (weak)", al)][D]["sim"], res[("thick (weak)", al)][D]["M0"]) for al in ALPHAS],
                     "o-", color=col, label=f"M0, Δ = {D} µm")
        ax[2].loglog(ALPHAS, [rel(res[("thick (weak)", al)][D]["sim"], res[("thick (weak)", al)][D]["M2"]) for al in ALPHAS],
                     "s--", color=col, lw=0.8, ms=3, label=f"M2, Δ = {D} µm")
    ax[2].set(xlabel="source cone half-angle (deg, in medium)", ylabel="relative L2 error (%)",
              title="Error of the separable model (M0) and of the\nderived model (M2), thick weak object")
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    p = out_path("partial_coherence.png")
    fig.savefig(p, dpi=130)
    print(f"\nwrote {p}")
    print(f"run time {time.time() - t_start:.1f} s")


if __name__ == "__main__":
    main()
