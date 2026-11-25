"""Shared math for the GODT shearing / tomography experiments.

Single source of truth for the shear transfer factor and the ODT band
geometry. Experiments import from here rather than re-deriving, so a fix
lands everywhere at once.

References
----------
Winnik et al., "Gradient optical diffraction tomography", Communications
Physics 8:510 (2025), doi:10.1038/s42005-025-02416-3.
  Eq. 16  -> infinitesimal-shear factor  j*K
  Eq. 20  -> finite-shear factor  (exp(j*K*d) - 1)/d
Wolf, Opt. Commun. 1, 153 (1969) -> Ewald/band geometry.
Brenner & Singer, Appl. Opt. 32, 4984 (1993) -> wave-propagation method (WPM),
  the angle-unrestricted multiple-scattering forward model used by exp06.
Osnabrugge, Leedumrongwatthanakun & Vellekoop, J. Comput. Phys. 322, 113
  (2016), arXiv:1601.05997 -> convergent Born series (CBS), the exact
  Helmholtz solver exp09 uses as a referee. Not ingested as a wiki source.

Wiki: wiki/notes/finite-shear-null-experiment.md          (H, absH, ... : exp01/02)
      wiki/notes/multishear-fusion-design-problem.md      (the model this implements)
      wiki/notes/multishear-fusion-design-solution.md     (KGrid, cap_*, Design*: exp03/04/05)
      wiki/notes/multiple-scattering-check.md              (cell_phantom, wpm_*: exp06)
      wiki/notes/colour-multiplexed-shear.md               (at_wavelength, unmix_shares: exp08)
      wiki/notes/full-wave-referee.md                      (mie_scalar, CBS, MSB: exp09)
"""
import numpy as np

__all__ = [
    "GODT_SIM", "GODT_EXP",
    "H", "absH", "combined_G2", "combined_G",
    "band_limit", "nulls_in_band", "wiener_rel_err",
    "shepp_logan_2d",
    # 3D Ewald sampling / Fisher information / design (exp03+)
    "G2_vec", "visibility", "KGrid", "flipK", "hermitian_sym", "azimuths",
    "cap_samples", "cap_density", "cap_density_closed", "cached_density",
    "DesignProblem", "solve_design", "solve_maximin", "chromatic_G",
    # multiple-scattering forward model (exp06+)
    "cell_phantom", "snap_illumination", "WPM", "angular_propagate",
    "shear_rytov", "meas_scale", "cap_geometry", "born_rytov_ratio", "deposit",
    "shell_spectrum", "hsym_c", "fuse",
    # wavelength as a shear channel (exp08+)
    "at_wavelength", "unmix_shares",
    # exact referee and a Born-consistent multiple-scattering model (exp09+)
    "mie_scalar", "born_sphere", "CBS", "MSB", "lattice_samples",
    # gridding-free reconstruction (exp10+)
    "ColumnWiener",
]


# ---------------------------------------------------------------- configs
class Config:
    """Optical configuration of an illumination-scanning ODT setup."""

    def __init__(self, name, lam, n0, na_det, theta_ill_deg, dx=None, shear=None):
        self.name = name
        self.lam = lam                      # vacuum wavelength, um
        self.n0 = n0                        # immersion index
        self.na_det = na_det                # detection NA
        self.theta_ill = theta_ill_deg      # illumination inclination, deg
        self.dx = dx                        # sampling pitch, um (optional)
        self.shear = shear                  # effective shear, um (optional)

    @property
    def k_det(self):
        return 2 * np.pi * self.na_det / self.lam

    @property
    def k_ill(self):
        return 2 * np.pi * self.n0 * np.sin(np.deg2rad(self.theta_ill)) / self.lam

    @property
    def kmax(self):
        """Max |K_x| reachable in the object spectrum.

        K = k_scattered - k_incident, so transverse reach is the sum of the
        detection and illumination transverse wavenumbers.
        """
        return self.k_det + self.k_ill

    @property
    def k_nyquist(self):
        return None if self.dx is None else np.pi / self.dx

    # ---- 3D geometry (exp03+). Convention: K = k_scattered - k_incident,
    # transmission (k_z > 0), illumination k_i = k_m (sin th cos phi, sin th sin phi, cos th).
    kz_floor_frac = 0.1
    """Samples with scattered k_z < kz_floor_frac * k_m are discarded. The
    Rytov mapping factor k0^2/(2 i k_z) diverges at grazing scattered angles,
    which only bites when NA = n0 (GODT_SIM); for GODT_EXP it never acts."""

    @property
    def k_m(self):
        """Wavenumber in the immersion medium."""
        return 2 * np.pi * self.n0 / self.lam

    @property
    def k_iz(self):
        """Axial component of every illumination wavevector (conical scan)."""
        return self.k_m * np.cos(np.deg2rad(self.theta_ill))

    @property
    def kz_floor(self):
        return self.kz_floor_frac * self.k_m

    @property
    def k_perp_max(self):
        """Largest scattered transverse wavenumber kept: the detection NA,
        or the k_z floor if that bites first."""
        return min(self.k_det, np.sqrt(self.k_m ** 2 - self.kz_floor ** 2))

    @property
    def Kz_range(self):
        """(min, max) of K_z = k_z - k_iz over all cap samples."""
        kz_lo = max(self.kz_floor, np.sqrt(max(self.k_m ** 2 - self.k_perp_max ** 2, 0.0)))
        return kz_lo - self.k_iz, self.k_m - self.k_iz

    def __repr__(self):
        return (f"<Config {self.name}: lam={self.lam} n0={self.n0} "
                f"NA={self.na_det} ill={self.theta_ill:.4g}deg "
                f"kmax={self.kmax:.2f} rad/um>")


# Numerical-simulation configuration used in the GODT paper.
GODT_SIM = Config("GODT simulation", lam=0.5, n0=1.33, na_det=1.33,
                  theta_ill_deg=60.0, dx=0.1, shear=0.1)

# Experimental configuration (published version). NOTE: the arXiv preprint
# reported an effective shear of 0.33 um for the same setup; the published
# Communications Physics version reports 0.8 um. Use the published value.
# `theta_ill` is the inclination INSIDE the immersion medium. The system paper
# (2024-06-14-598984-polarization-hilbert-holotomography, ref. 84 of the
# GODT paper) gives "illumination incidence in air at 58 deg", i.e. an
# illumination NA of sin 58 = 0.848, which is 34.04 deg in n0 = 1.515.
# (Corrected 2026-09-30: this used to read 58.0, putting 58 deg inside the oil
# and overstating k_ill 15.2 -> 10.05 rad/um. Only exp01's band printout used it.)
GODT_EXP = Config("GODT experiment", lam=0.53, n0=1.515, na_det=1.45,
                  theta_ill_deg=float(np.degrees(np.arcsin(np.sin(np.radians(58.0)) / 1.515))),
                  shear=0.8)


# ------------------------------------------------------- shear transfer
def H(K, d):
    """Complex finite-difference shear transfer factor (Eq. 20).

    K may be an array of scalar frequencies (1D, shear assumed along K) or,
    for vector shear, pass the projection K.d_hat * |d| via `K` and the
    magnitude via `d` -- see `H_vec` usage in the 2D experiments.
    """
    return (np.exp(1j * np.asarray(K) * d) - 1.0) / d


def absH(K, d):
    """|H| in closed form.  |exp(ix)-1| = 2|sin(x/2)|.

    Zeros at K = 2*pi*m/d.  Two distinct structures:
      m == 0 : the DC loss inherent to any gradient method; independent of d.
      m != 0 : the finite-shear comb; moves with d, fillable by multi-shear.
    """
    return 2.0 * np.abs(np.sin(np.asarray(K) * d / 2.0)) / d


def H_vec(KX, KY, dvec):
    """Complex transfer factor for a 2D shear vector dvec = (dx, dy)."""
    mag = np.hypot(*dvec)
    proj = KX * dvec[0] + KY * dvec[1]
    return (np.exp(1j * proj) - 1.0) / mag


def combined_G2(K, shears):
    """sum_n |H_n(K)|^2 for a set of scalar shear magnitudes."""
    return sum(absH(K, d) ** 2 for d in shears)


def combined_G(K, shears):
    return np.sqrt(combined_G2(K, shears))


def band_limit(cfg):
    """(kmax, k_det, k_ill) for a Config."""
    return cfg.kmax, cfg.k_det, cfg.k_ill


def nulls_in_band(d, kmax):
    """Positions of the m != 0 transfer nulls inside |K| <= kmax."""
    if d <= 0:
        return []
    first = 2 * np.pi / d
    n = int(np.floor(kmax / first))
    return [m * first for m in range(1, n + 1)]


def wiener_rel_err(K, shears, rho):
    """Relative recovery error of the optimal linear estimator.

        rel_err(K) = [ rho * G^2(K) + 1 ]^(-1/2)

    rho is the per-frequency SNR (signal power / noise power). At a null
    G = 0 so the error is 1 regardless of rho: the frequency is lost
    outright, not merely noisy. Validated against Monte Carlo in
    exp02_shear_recovery.py.
    """
    return np.sqrt(1.0 / (rho * combined_G2(K, shears) + 1.0))


# ------------------------------------------------------------- phantom
def shepp_logan_2d(X, Y, scale=20.0):
    """Modified Shepp-Logan phantom evaluated on a coordinate grid in um.

    `scale` sets the phantom radius in um (default 20 um, matching the
    40 um extent of the phantom used in the GODT paper's simulations).
    """
    E = [  # (intensity, a, b, x0, y0, phi_deg)
        (1.00, .69, .92, 0, 0, 0), (-.80, .6624, .874, 0, -.0184, 0),
        (-.20, .11, .31, .22, 0, -18), (-.20, .16, .41, -.22, 0, 18),
        (.10, .21, .25, 0, .35, 0), (.10, .046, .046, 0, .1, 0),
        (.10, .046, .046, 0, -.1, 0), (.10, .046, .023, -.08, -.605, 0),
        (.10, .023, .023, 0, -.606, 0), (.10, .023, .046, .06, -.605, 0),
    ]
    out = np.zeros_like(X, dtype=float)
    xs, ys = X / scale, Y / scale
    for I, a, b, x0, y0, phi in E:
        t = np.deg2rad(phi)
        xr = (xs - x0) * np.cos(t) + (ys - y0) * np.sin(t)
        yr = -(xs - x0) * np.sin(t) + (ys - y0) * np.cos(t)
        out[(xr / a) ** 2 + (yr / b) ** 2 <= 1] += I
    return out


# ---------------------------------------- dimensionless shear power (exp03+)
def G2_vec(KX, KY, dvec):
    """Dimensionless shear power |e^{i d.K} - 1|^2 = 4 sin^2(d.K / 2).

    This, not |H|^2, enters the Fisher information: the measured phase
    difference psi = phi(x+d) - phi(x) carries noise that does not scale with
    |d| (assumption A4 of the design note), so the 1/d of `H` is a
    convention that cancels. Peak 4 at d.K = pi (mod 2 pi); zeros at
    d.K = 2 pi m, m = 0 (blind plane) and m != 0 (finite-shear comb).
    """
    proj = np.asarray(KX) * dvec[0] + np.asarray(KY) * dvec[1]
    return 4.0 * np.sin(proj / 2.0) ** 2


def visibility(d, lc=None):
    """|mu(d)|^2: fringe-visibility factor of a shear of magnitude d.

    MODEL, not sourced: circular incoherent source, van Cittert-Zernike
    jinc  mu = 2 J1(pi d/lc) / (pi d/lc), first zero at d = 1.22 lc.
    lc=None or inf -> fully coherent (1). GODT's real mu(d) is unknown.
    """
    d = np.asarray(d, dtype=float)
    if lc is None or not np.isfinite(lc):
        return np.ones_like(d)
    from scipy.special import j1
    x = np.pi * d / lc
    with np.errstate(divide="ignore", invalid="ignore"):
        mu = np.where(x == 0, 1.0, 2.0 * j1(x) / x)
    return mu ** 2


def chromatic_G(K, d0, lam0, fwhm, n_lam=61, nsig=3.0):
    """Broadband effective shear factor under a Gaussian source spectrum.

    Polarization-grating shear scales with wavelength (diffraction angle
    sin th = lam / Gamma, GODT: Gamma = 6.29 um), so d(lam) = d0 lam/lam0.
    For a weak object the interference term is linear in the object and
    the wavelengths add in intensity, so the effective filter is the
    spectral AVERAGE of G (coherent-average model, inference):
        G_eff(K) = sum_i s_i (e^{i d(lam_i) K} - 1).
    Returns (G_eff, G2_incoherent) where G2_incoherent = sum_i s_i |G_i|^2
    is the (unreachable) bound if each wavelength were detected separately.
    """
    sig = fwhm / (2 * np.sqrt(2 * np.log(2)))
    lam = np.linspace(lam0 - nsig * sig, lam0 + nsig * sig, n_lam)
    s = np.exp(-0.5 * ((lam - lam0) / sig) ** 2)
    s /= s.sum()
    K = np.asarray(K, dtype=float)
    ph = np.exp(1j * np.outer(lam / lam0 * d0, K))            # (n_lam, K)
    G = ph - 1.0
    return (s[:, None] * G).sum(axis=0), (s[:, None] * np.abs(G) ** 2).sum(axis=0)


# ------------------------------------- wavelength as a shear channel (exp08+)
def at_wavelength(cfg, lam):
    """The same instrument run at vacuum wavelength `lam` (um).

    Angles and apertures are geometric and stay put (theta_ill, NA, n0 --
    immersion dispersion ignored), so every wavenumber scales as 1/lam: the
    Ewald radius, the band edge and the cap family are a scaled copy. The
    polarization-grating shear scales the other way, d = d0 lam/lam0
    (`chromatic_G`), which is what makes a colour a distinct shear channel.
    """
    shear = None if cfg.shear is None else cfg.shear * lam / cfg.lam
    return Config(f"{cfg.name} {1e3 * lam:.0f}nm", lam=lam, n0=cfg.n0,
                  na_det=cfg.na_det, theta_ill_deg=cfg.theta_ill, dx=cfg.dx, shear=shear)


def unmix_shares(X, area):
    """Effective time shares of colour channels read through a mosaic sensor.

    MODEL, not sourced. Pixel class p covers a fraction area[p] of the sensor
    and responds to colour c with relative sensitivity X[p, c] (X[p, p] = 1,
    off-diagonals = crosstalk), all colours at equal irradiance. Shot noise in
    class p then scales as sum_c X[p, c] / area[p]; unmixing with X^-1 gives
    colour c the variance [X^-1 diag(v) X^-T]_cc, and its reciprocal is the
    share of a monochrome full-sensor exposure that colour is worth. X = I
    returns `area`: a mosaic by itself costs nothing in total photons.
    """
    X = np.asarray(X, dtype=float)
    v = X.sum(axis=1) / np.asarray(area, dtype=float)
    Xi = np.linalg.inv(X)
    return 1.0 / np.einsum("cp,p,cp->c", Xi, v, Xi)


# ---------------------------------------------- 3D Ewald sampling (exp03+)
class KGrid:
    """FFT-ordered Cartesian grid of the object spectrum.

    n transverse voxels at real-space pitch dx (um) -> K pitch dk = 2 pi/(n dx);
    nz axial voxels at the SAME K pitch, so real-space dz = n dx / nz. All
    axes are in numpy FFT order so that `flipK` maps K -> -K exactly and a
    real-space volume of shape (n, n, nz) FFTs straight onto the grid.
    """

    def __init__(self, n, dx, nz):
        self.n, self.dx, self.nz = int(n), float(dx), int(nz)
        self.dk = 2 * np.pi / (self.n * self.dx)
        self.dz = self.n * self.dx / self.nz
        self.kx = self.dk * np.fft.fftfreq(self.n, 1.0 / self.n)
        self.ky = self.kx.copy()
        self.kz = self.dk * np.fft.fftfreq(self.nz, 1.0 / self.nz)
        self.KX, self.KY = np.meshgrid(self.kx, self.ky, indexing="ij")
        self.KP = np.hypot(self.KX, self.KY)

    @property
    def shape(self):
        return (self.n, self.n, self.nz)

    @property
    def size(self):
        return self.n * self.n * self.nz

    @property
    def kz_max(self):
        """Largest |K_z| representable without wrap-around."""
        return self.dk * (self.nz // 2)

    def check(self, cfg):
        lo, hi = cfg.Kz_range
        assert -self.kz_max < lo - self.dk and hi + self.dk < self.kz_max, \
            "K_z axis too short for the caps"
        assert cfg.kmax <= self.dk * (self.n // 2) + 1e-9, "transverse axis too short"

    def K3(self):
        return np.meshgrid(self.kx, self.ky, self.kz, indexing="ij")

    def Kabs(self):
        KX, KY, KZ = self.K3()
        return np.sqrt(KX ** 2 + KY ** 2 + KZ ** 2)

    def __repr__(self):
        return (f"<KGrid {self.n}x{self.n}x{self.nz}, dx={self.dx} um, "
                f"dz={self.dz:.3f} um, dk={self.dk:.4f} rad/um>")


def flipK(A):
    """Map voxel K -> -K on an FFT-ordered 3D grid (Hermitian partner)."""
    return np.roll(A[::-1, ::-1, ::-1], 1, axis=(0, 1, 2))


def hermitian_sym(A):
    """A(K) + A(-K), without double-counting the voxels that are their own
    partner (the origin and, on even axes, the Nyquist planes). For a real
    object the samples at K and -K inform the same parameter."""
    S = A + flipK(A)
    selfs = []
    for n in A.shape:
        m = np.zeros(n, dtype=bool)
        m[0] = True
        if n % 2 == 0:
            m[n // 2] = True
        selfs.append(m)
    self3 = selfs[0][:, None, None] & selfs[1][None, :, None] & selfs[2][None, None, :]
    S[self3] = A[self3]
    return S


def azimuths(n_az):
    return 2 * np.pi * np.arange(n_az) / n_az


def cap_samples(cfg, grid, phi):
    """Digital Ewald cap of one illumination azimuth `phi` (rad).

    One sample per transverse grid frequency K_perp inside the shifted
    detection disc |K_perp + k_i,perp| <= k_perp_max, deposited at the
    nearest K_z voxel of K_z = k_z - k_iz (assumption A3: the cap IS this
    digital surface, so theory and simulation share one discretisation).
    Returns (ix, iy, iz, kz) with kz the scattered axial wavenumber.
    """
    kix, kiy = cfg.k_ill * np.cos(phi), cfg.k_ill * np.sin(phi)
    kxp, kyp = grid.KX + kix, grid.KY + kiy
    kp2 = kxp * kxp + kyp * kyp
    ok = kp2 <= cfg.k_perp_max ** 2
    ix, iy = np.nonzero(ok)
    kz = np.sqrt(cfg.k_m ** 2 - kp2[ok])
    iz = np.rint((kz - cfg.k_iz) / grid.dk).astype(np.int64) % grid.nz
    return ix, iy, iz, kz


def cap_density(cfg, grid, n_az, n_frames=180, chunk=32):
    """M(K): Rytov-weighted cap multiplicity per voxel, by direct binning.

        M(K) = (n_frames / n_az) * sum_az sum_{samples in voxel} (k_m / k_z)^2

    = expected weighted number of samples a voxel receives from an
    n_frames-projection acquisition, estimated with n_az dense azimuths.
    (k_m/k_z)^2 is the squared Rytov mapping factor k0^2/(2 i k_z),
    normalised to 1 for paraxial scattering: a sample at scattered k_z
    carries (k_m/k_z)^2 times the information about O~ of a paraxial one.
    Not Hermitian-symmetrised; callers add flipK(M) for a real object.
    """
    n, nz = grid.n, grid.nz
    M = np.zeros(grid.size)
    phis = azimuths(n_az)
    for start in range(0, n_az, chunk):
        lins, wts = [], []
        for phi in phis[start:start + chunk]:
            ix, iy, iz, kz = cap_samples(cfg, grid, phi)
            lins.append((ix * n + iy) * nz + iz)
            wts.append((cfg.k_m / kz) ** 2)
        M += np.bincount(np.concatenate(lins), weights=np.concatenate(wts),
                         minlength=grid.size)
    return M.reshape(grid.shape) * (n_frames / n_az)


def cap_density_closed(cfg, grid, n_frames=180):
    """Closed-form continuum limit of `cap_density` (design note, T1).

    Under conical illumination every cap is a translate of one spherical
    cap, so the scattered k_z = K_z + k_m cos(th) depends on K_z alone, and
    the azimuths whose cap passes through K solve

        cos(phi - psi) = c(K) = -(|K|^2 + 2 k_m K_z cos th) / (2 k_m K_perp sin th)

    (psi = azimuth of K_perp): two roots for |c| < 1, none otherwise. With
    |dK_z/dphi| = (k_m sin th / k_z) K_perp sqrt(1 - c^2) the weighted sample
    count in a voxel of edge dk is

        M(K) = (n_frames dk / pi) * k_m / (k_z sin th K_perp sqrt(1 - c^2))

    on {|c| < 1, k_z in [kz_lo, k_m]}, zero elsewhere. The 1/sqrt divergence
    on the caustic |c| = 1 (envelope of the cap family) is integrable.
    Near the origin |c| = 1 reduces to K_z = +-K_perp tan(th): the missing
    cone of half-angle 90 deg - th. Verified against binning in selftest.
    """
    KX, KY, KZ = grid.K3()
    KP = np.hypot(KX, KY)
    th = np.deg2rad(cfg.theta_ill)
    s, c = np.sin(th), np.cos(th)
    km = cfg.k_m
    kz = KZ + km * c
    kz_lo = max(cfg.kz_floor, np.sqrt(max(km ** 2 - cfg.k_perp_max ** 2, 0.0)))
    with np.errstate(divide="ignore", invalid="ignore"):
        cc = -(KX ** 2 + KY ** 2 + KZ ** 2 + 2 * km * KZ * c) / (2 * km * KP * s)
        inside = (kz >= kz_lo) & (kz <= km) & (KP > 0) & (np.abs(cc) < 1)
        M = np.where(inside,
                     n_frames * grid.dk / np.pi * km / (kz * s * KP * np.sqrt(1 - cc ** 2)),
                     0.0)
    return M


def cached_density(cfg, grid, n_az, n_frames=180):
    """Hermitian-symmetrised M(K) + flipK(M(K)), cached under out/."""
    import os
    tag = cfg.name.replace(" ", "_")
    p = out_path(f"M_{tag}_{grid.n}_{grid.nz}_{grid.dx}_{n_az}_{n_frames}.npy")
    if os.path.exists(p):
        return np.load(p)
    M = hermitian_sym(cap_density(cfg, grid, n_az, n_frames))
    np.save(p, M)
    return M


# ----------------------------------------------- optimal design (exp04+)
class DesignProblem:
    """Convex time-share design over N candidate shear channels.

    A   : 3D array  s(K) rho M(K)  (prior shape x per-frame SNR x symmetrised
          density), zero off the support.  S(K) J(K) = A(K) S_t(K_perp).
    G2  : (N, n, n) dimensionless shear powers |G_j(K_perp)|^2
    v   : (N,) visibility factors |mu(d_j)|^2
    P   : optional (N, n, n) squared projections (K_perp . d_hat_j)^2, only
          needed for the calibration-robust objective (E1).
    w   : optional 3D weight (default 1). Objectives are per-support-voxel means.
    """

    def __init__(self, A, G2, v, P=None, w=None):
        self.A = np.asarray(A, dtype=np.float64)
        self.supp = self.A > 0
        self.nsupp = int(self.supp.sum())
        self.E = np.asarray(G2, dtype=np.float64) * np.asarray(v, dtype=np.float64)[:, None, None]
        self.N = self.E.shape[0]
        self.P = None if P is None else np.asarray(P, dtype=np.float64)
        self.w = None if w is None else np.asarray(w, dtype=np.float64)

    def S_t(self, t):
        return np.tensordot(np.asarray(t, dtype=np.float64), self.E, axes=(0, 0))

    def relmse(self, t):
        """3D relative Wiener MSE 1/(A S_t + 1); equals 1 off the support."""
        return 1.0 / (self.A * self.S_t(t)[:, :, None] + 1.0)

    def _wsum(self, X):
        return (X if self.w is None else self.w * X)

    def f_A(self, t):
        """(P-A): mean relative MSE over the support, and its gradient."""
        den = self.A * self.S_t(t)[:, :, None] + 1.0
        val = self._wsum(1.0 / den)[self.supp].sum() / self.nsupp
        D = self._wsum(self.A / den ** 2).sum(axis=2)
        grad = -np.tensordot(self.E, D, axes=([1, 2], [0, 1])) / self.nsupp
        return val, grad

    def f_D(self, t):
        """(P-D): minus the mean log(1 + A S_t) over the support (to minimise)."""
        x = self.A * self.S_t(t)[:, :, None]
        val = -self._wsum(np.log1p(x))[self.supp].sum() / self.nsupp
        D = self._wsum(self.A / (1.0 + x)).sum(axis=2)
        grad = -np.tensordot(self.E, D, axes=([1, 2], [0, 1])) / self.nsupp
        return val, grad

    def f_A_robust(self, t, s):
        """(P-A) with an independent shear-calibration error of rms s (um)
        per channel, treated as coherent model error (E1 of the design note):

            S J = sum_j t_j A E_j / (1 + t_j A s^2 P_j).

        Each term is t/(c + d t): concave, so the problem stays convex.
        Ceiling as t_j -> inf:  E_j / (s^2 P_j) -- calibration-limited SNR.
        """
        t = np.asarray(t, dtype=np.float64)
        J = np.zeros(self.A.shape)
        grad = np.zeros(self.N)
        dens = []
        for j in range(self.N):
            if t[j] <= 0 and s > 0:
                dens.append(None)
                continue
            d = 1.0 + t[j] * self.A * (s ** 2 * self.P[j])[:, :, None]
            J += t[j] * self.A * self.E[j][:, :, None] / d
            dens.append(d)
        tot = J + 1.0
        val = self._wsum(1.0 / tot)[self.supp].sum() / self.nsupp
        w2 = self._wsum(1.0 / tot ** 2)
        for j in range(self.N):
            d = dens[j]
            if d is None:
                d = 1.0
            dJ = self.A * self.E[j][:, :, None] / d ** 2
            grad[j] = -(w2 * dJ)[self.supp].sum() / self.nsupp
        return val, grad

    def certificate(self, t, f):
        """Equivalence-theorem gap. For convex f on the simplex, t is optimal
        iff min_j g_j(t) >= sum_j t_j g_j(t) (no vertex is a descent
        direction). Returns (sum_j t_j g_j - min_j g_j)/|f(t)| >= 0, zero at
        the optimum. For (P-A) this is the Kiefer-Wolfowitz-type statement
        max_j phi_j(t) <= sum_j t_j phi_j(t), phi_j = sum_K w a_j/(J + 1/S)^2.
        """
        val, g = f(t)
        return (np.dot(t, g) - g.min()) / abs(val)


def solve_design(prob, f, t0=None, maxiter=400, ftol=1e-12, prune=1e-4,
                 gap_tol=1e-4, rounds=12):
    """Minimise the convex objective f (f_A, f_D, or a lambda) over the simplex.

    SLSQP with analytic gradient on the full candidate set, then
    Fedorov-Wynn-style rounds: prune shares below `prune`, add the vertex
    with the most negative equivalence-theorem gap, re-solve on that support
    until the gap is below `gap_tol`. Returns (t, value, gap).
    """
    from scipy.optimize import minimize
    N = prob.N
    t = np.full(N, 1.0 / N) if t0 is None else np.asarray(t0, dtype=float)

    def run(idx, tstart):
        def fm(x):
            tt = np.zeros(N)
            tt[idx] = x
            v, g = f(tt)
            return v, g[idx]
        x0 = np.clip(tstart[idx], 0, None)
        x0 = x0 / x0.sum() if x0.sum() > 0 else np.full(len(idx), 1.0 / len(idx))
        res = minimize(fm, x0, jac=True, method="SLSQP",
                       bounds=[(0.0, 1.0)] * len(idx),
                       constraints=[{"type": "eq", "fun": lambda x: x.sum() - 1.0,
                                     "jac": lambda x: np.ones_like(x)}],
                       options={"maxiter": maxiter, "ftol": ftol})
        tt = np.zeros(N)
        tt[idx] = np.clip(res.x, 0, None)
        return tt / tt.sum()

    t = run(np.arange(N), t)
    best_val = f(t)[0]
    for _ in range(rounds):
        val, g = f(t)
        gap = (np.dot(t, g) - g.min()) / abs(val)
        if gap < gap_tol:
            break
        keep = t > prune
        keep[np.argmin(g)] = True
        idx = np.nonzero(keep)[0]
        tstart = np.where(keep, np.maximum(t, prune), 0.0)
        t_new = run(idx, tstart)
        v_new = f(t_new)[0]
        if v_new <= best_val + 1e-14:
            t, best_val = t_new, v_new
        else:
            break
    t[t < prune] = 0.0
    t /= t.sum()
    val, g = f(t)
    gap = (np.dot(t, g) - g.min()) / abs(val)
    return t, val, gap


def solve_maximin(G2, v, mask2d, subsample=2):
    """(P-W) maximin-efficiency design:

        max_t min_{K_perp in mask}  S_t(K_perp) / max_j v_j |G_j(K_perp)|^2

    i.e. the worst-case fraction of the information the best single channel
    could deliver at that transverse frequency. Scale-free, so no low-K
    cutoff is needed (the ratio stays finite as K_perp -> 0). Solved as an LP
    (scipy HiGHS) on every `subsample`-th grid point, then evaluated on the
    full mask. Returns (t, z_lp, z_full).
    """
    from scipy.optimize import linprog
    E = np.asarray(G2, dtype=float) * np.asarray(v, dtype=float)[:, None, None]
    Emax = E.max(axis=0)
    pts = mask2d & (Emax > 0)
    eff = E / np.where(Emax > 0, Emax, 1.0)
    sub = np.zeros_like(pts)
    sub[::subsample, ::subsample] = True
    Pm = eff[:, pts & sub].T
    N = E.shape[0]
    c = np.zeros(N + 1)
    c[-1] = -1.0
    A_ub = np.hstack([-Pm, np.ones((Pm.shape[0], 1))])
    res = linprog(c, A_ub=A_ub, b_ub=np.zeros(Pm.shape[0]),
                  A_eq=np.r_[np.ones(N), 0.0][None, :], b_eq=[1.0],
                  bounds=[(0, None)] * N + [(None, None)], method="highs")
    t = np.clip(res.x[:N], 0, None)
    t /= t.sum()
    z_full = np.tensordot(t, eff, axes=(0, 0))[pts].min()
    return t, -res.fun, z_full



# ------------------------------- multiple-scattering forward model (exp06+)
#
# Everything above linearises: the measurement is assumed to be
# h(K) O~(K) + noise. This section builds the nonlinear alternative --
# Brenner & Singer's wave-propagation method (WPM) -- so that the design
# conclusions can be re-tested on data the linear model did not generate.
# WPM is GODT's own forward engine (its ref 78) and, unlike BPM, is valid at
# the 60 deg illumination inclination GODT uses.
#
# Conventions fixed here and used by exp06 (derivation in the wiki note):
#   O = k0^2 (n_m^2 - n^2)          so that  (lap + k_m^2) u = O u
#   first Born, plane z = 0:        U~_s(k_perp) = -i O~(K) / (2 k_z)
#   first Rytov:                    U~_R(K_perp) = -i O~(K) / (2 k_z)
#   arrays are centred (index n//2 <-> x = 0), exactly as cell_phantom builds
#   them, so the DFT of a centred array carries e^{-i K.r_c}; the transverse
#   ramp cancels between the two sides and the axial one does not, leaving
#
#       meas_scale * fft2(u_R) * e^{-i K_z z_c}  ==  (k_m/k_z) * fftn(O)
#
#   with meas_scale = 2i k_m / dz and z_c = (nz//2) dz. The right-hand side is
#   exactly the "h * Ot" that exp05's `acquire` synthesises, so WPM data and
#   linear data are interchangeable inputs to the same estimator.


def meas_scale(cfg, grid):
    """2i k_m / dz -- puts a Rytov field into the units exp05's `acquire` uses."""
    return 2j * cfg.k_m / grid.dz


def cell_phantom(grid, rng, scale=1.0, contrast=1.0, cfg=None):
    """Textured cell-like scattering potential O = k0^2 (n_m^2 - n^2), centred.

    Cytoplasm sphere, nucleus, nucleolus, 25 organelles of either sign, and a
    fine random texture (0.06 um correlation, NOT rescaled) so the spectrum
    keeps power out to the band edge. Returns (O, dn).

    `scale` multiplies every length (radii and centres), `contrast` every index
    step. scale=contrast=1 is the exp05 phantom, bit-for-bit: the rng is drawn
    in the same order regardless of the parameters.
    """
    from scipy.ndimage import gaussian_filter
    cfg = cfg or GODT_SIM
    n, nz = grid.n, grid.nz
    x = (np.arange(n) - n // 2) * grid.dx
    z = (np.arange(nz) - nz // 2) * grid.dz
    X, Y, Z = np.meshgrid(x, x, z, indexing="ij")
    dn = np.zeros(X.shape)
    s = float(scale)

    def sphere(c, r):
        c = (c[0] * s, c[1] * s, c[2] * s)
        return (X - c[0]) ** 2 + (Y - c[1]) ** 2 + (Z - c[2]) ** 2 <= (r * s) ** 2
    cyto = sphere((0, 0, 0), 4.8)
    dn[cyto] += 0.020
    nc = (1.0, 0.5, -0.3)
    dn[sphere(nc, 2.2)] += 0.015
    dn[sphere((1.6, 0.5, 0.0), 0.6)] += 0.030
    k = 0
    while k < 25:
        c = rng.uniform(-4.2, 4.2, 3)
        if np.linalg.norm(c) < 4.2 and np.linalg.norm(c - nc) > 2.7:
            dn[sphere(c, rng.uniform(0.25, 0.5))] += rng.choice([-1, 1]) * rng.uniform(0.008, 0.015)
            k += 1
    tex = gaussian_filter(rng.standard_normal(X.shape),
                          sigma=[0.06 / grid.dx, 0.06 / grid.dx, 0.06 / grid.dz])
    tex *= 0.006 / tex[cyto].std()
    dn += tex * cyto
    dn = dn * contrast
    nm, k0 = cfg.n0, 2 * np.pi / cfg.lam
    return k0 ** 2 * (nm ** 2 - (nm + dn) ** 2), dn


def snap_illumination(cfg, grid, n_az, tol=0.75, dk=None):
    """Incident transverse wavevectors snapped to the grid lattice.

    WPM propagates by FFT, so the incident plane wave must be exactly periodic
    on the transverse grid: k_i,perp has to be an integer multiple of dk. This
    picks, for each of `n_az` target azimuths, the nearest unused lattice point
    within `tol` grid units of the circle |k_perp| = k_m sin(theta_ill).

    `dk` defaults to the grid's own pitch; pass the PROPAGATION grid's pitch
    when WPM runs on a padded grid, since that is the one periodicity is
    required on. The reconstruction never needs k_i on its own lattice.

    Returns (KI, theta_eff) with KI shape (n_az, 2) in rad/um and theta_eff the
    per-azimuth illumination inclination actually realised (deg) -- snapping
    perturbs it by a fraction of a degree, and exp06 uses the realised value,
    not the nominal one, everywhere.
    """
    dk = grid.dk if dk is None else float(dk)
    R = cfg.k_ill / dk
    m = int(np.ceil(R + tol + 1))
    p, q = np.meshgrid(np.arange(-m, m + 1), np.arange(-m, m + 1), indexing="ij")
    r = np.hypot(p, q)
    ok = np.abs(r - R) <= tol
    P, Q, RR = p[ok], q[ok], r[ok]
    if P.size < n_az:
        raise ValueError(f"only {P.size} lattice points within {tol} dk of the "
                         f"illumination circle; raise tol or the grid size")
    ang = np.arctan2(Q, P)
    used = np.zeros(P.size, bool)
    KI = np.empty((n_az, 2))
    for j, phi in enumerate(azimuths(n_az)):
        d = np.abs(np.angle(np.exp(1j * (ang - phi))))
        d = d + 0.05 * np.abs(RR - R)      # angle first, radius as tiebreak
        d[used] = np.inf
        i = int(np.argmin(d))
        used[i] = True
        KI[j] = dk * P[i], dk * Q[i]
    kp = np.hypot(KI[:, 0], KI[:, 1])
    return KI, np.rad2deg(np.arcsin(np.clip(kp / cfg.k_m, -1, 1)))


def _fft2(a):
    import scipy.fft as sfft
    return sfft.fft2(a, workers=-1)


def _ifft2(a):
    import scipy.fft as sfft
    return sfft.ifft2(a, workers=-1)


def angular_propagate(U, dx, k_m, dz):
    """Exact angular-spectrum propagation through the homogeneous medium.

    Evanescent components are dropped (WPM drops them too, and for dz < 0 --
    back-propagation to the object centre -- they would blow up).
    """
    n = U.shape[0]
    k = 2 * np.pi * np.fft.fftfreq(n, dx)
    KX, KY = np.meshgrid(k, k, indexing="ij")
    kt2 = KX ** 2 + KY ** 2
    prop = np.where(kt2 < k_m ** 2,
                    np.exp(1j * np.sqrt(np.maximum(k_m ** 2 - kt2, 0.0)) * dz), 0.0)
    return _ifft2(_fft2(U) * prop)


class WPM:
    """Wave-propagation method (brenner1993) on a fixed index volume.

    The volume is quantised into index levels once; each slice then costs one
    forward FFT plus one inverse FFT per level *present in that slice*, so the
    free-space slices above and below the object cost a single FFT pair. The
    quantisation is built on dn = n - n0 with a uniform step through zero, so
    the background is represented exactly and free space propagates exactly.

    Per Brenner & Singer, the local phase is applied per plane-wave component
    with its own angle (kz = sqrt(n_l^2 k0^2 - k_perp^2)) rather than to an
    averaged index, which is what removes BPM's ~10 deg validity ceiling --
    mandatory here, since GODT illuminates at 60 deg.

    Absorbing boundary: a taper applied to the SCATTERED part only (u_i is known
    analytically and is left exact), so that light leaving the field of view
    sideways does not wrap around the periodic grid. `absorb_rate` is per
    micron of propagation, so the total absorption is independent of `dz`.
    """

    def __init__(self, nvol, dx, dz, lam, n0, n_levels=8,
                 absorb_frac=0.18, absorb_rate=0.6, dtype=np.complex64):
        self.dx, self.dz, self.lam, self.n0 = float(dx), float(dz), float(lam), float(n0)
        self.dtype = dtype
        self.nx, self.ny, self.ns = nvol.shape
        self.k0 = 2 * np.pi / self.lam
        self.k_m = self.k0 * self.n0
        k = 2 * np.pi * np.fft.fftfreq(self.nx, self.dx)
        self.KX, self.KY = np.meshgrid(k, k, indexing="ij")
        kt2 = self.KX ** 2 + self.KY ** 2
        x = (np.arange(self.nx) - self.nx // 2) * self.dx
        self.X, self.Y = np.meshgrid(x, x, indexing="ij")

        # ---- quantise dn on a uniform step through zero (background exact)
        d = nvol - self.n0
        span = float(np.abs(d).max())
        self.step = span / max(n_levels, 1) if span > 0 else 1.0
        qi = np.rint(d / self.step).astype(np.int32)
        self.levels = np.unique(qi)
        self.quant_err = float(np.abs(d - qi * self.step).max())
        # ---- one propagator per level (evanescent dropped, as in BPM/WPM)
        self.P = []
        for lv in self.levels:
            kz2 = ((self.n0 + lv * self.step) * self.k0) ** 2 - kt2
            self.P.append(np.where(kz2 > 0,
                                   np.exp(1j * np.sqrt(np.maximum(kz2, 0)) * self.dz),
                                   0.0).astype(dtype))
        # ---- per-slice level presence (computed once, reused every azimuth)
        self.present = [np.unique(qi[:, :, s]) for s in range(self.ns)]
        self.masks = {}
        for s in range(self.ns):
            if len(self.present[s]) > 1:
                self.masks[s] = {int(lv): (qi[:, :, s] == lv).astype(dtype)
                                 for lv in self.present[s]}
        self.lvindex = {int(lv): i for i, lv in enumerate(self.levels)}
        # ---- absorbing boundary (separable edge taper), applied to u - u_i
        if absorb_frac > 0:
            w = absorb_frac * self.nx * self.dx
            half = 0.5 * self.nx * self.dx
            t = np.clip((np.abs(x) - (half - w)) / w, 0.0, 1.0)
            # per unit LENGTH, not per step: otherwise halving dz doubles the
            # total absorption and a step-size convergence test measures the
            # absorber instead of the propagator.
            f = np.exp(-absorb_rate * self.dz * t ** 2)
            self.absorb = (f[:, None] * f[None, :]).astype(dtype)
        else:
            self.absorb = None

    @property
    def z0(self):
        """Entrance-face position. The volume is centred on z = 0, exactly as
        `cell_phantom` centres it, so that the incident plane wave carries its
        true absolute phase exp(i k_i.r) at every plane and u/u_i needs no
        bookkeeping constant."""
        return -0.5 * self.ns * self.dz

    def run(self, ki):
        """Propagate a unit plane wave with transverse wavevector ki = (kix, kiy)
        from the entrance face to the exit face. Returns the exit field."""
        kix, kiy = float(ki[0]), float(ki[1])
        kiz = np.sqrt(max(self.k_m ** 2 - kix ** 2 - kiy ** 2, 0.0))
        ramp = np.exp(1j * (kix * self.X + kiy * self.Y)).astype(self.dtype)
        U = (ramp * np.exp(1j * kiz * self.z0)).astype(self.dtype)
        for s in range(self.ns):
            Uk = _fft2(U)
            pres = self.present[s]
            if len(pres) == 1:
                U = _ifft2(Uk * self.P[self.lvindex[int(pres[0])]])
            else:
                V = np.zeros_like(U)
                msk = self.masks[s]
                for lv in pres:
                    lv = int(lv)
                    V += msk[lv] * _ifft2(Uk * self.P[self.lvindex[lv]])
                U = V
            if self.absorb is not None:
                Ui = ramp * np.exp(1j * kiz * (self.z0 + (s + 1) * self.dz))
                U = Ui + (U - Ui) * self.absorb
        return U.astype(np.complex128)

    def to_centre(self, U, ki):
        """Back-propagate an exit field to the plane through the object centre
        and divide out the incident wave there, giving u/u_i (whose log is u_R).

        Back-propagating the *exact* field to the object centre before taking
        the logarithm is chen1998's hybrid prescription, and is what makes
        first Rytov's validity range ~3x first Born's; taking the log at the
        exit face instead would handicap the linear model for the wrong reason.
        """
        kix, kiy = float(ki[0]), float(ki[1])
        Uc = angular_propagate(U, self.dx, self.k_m, self.z0)   # exit -> centre
        Ui = np.exp(1j * (kix * self.X + kiy * self.Y))         # u_i at z = 0
        return Uc / Ui


def shell_spectrum(Ot, Kabs, kmax, dk):
    """Shell-averaged |O~|^2 -> isotropic 3D prior S(K) by |K| lookup.
    Returns (S, edges, profile)."""
    edges = np.arange(0, kmax + 2 * dk, dk)
    idx = np.clip(np.digitize(Kabs.ravel(), edges) - 1, 0, len(edges) - 2)
    pw = np.abs(Ot.ravel()) ** 2
    num = np.bincount(idx, weights=pw, minlength=len(edges) - 1)
    cnt = np.bincount(idx, minlength=len(edges) - 1)
    prof = num / np.maximum(cnt, 1)
    prof[cnt == 0] = prof[cnt > 0].min()
    return prof[idx].reshape(Ot.shape), edges, prof


def hsym_c(num):
    """num(K) + conj(num(-K)); self-partner voxels keep Re(num).

    The complex counterpart of `hermitian_sym`: for a real object the samples
    at K and -K inform the same parameter."""
    S = num + np.conj(flipK(num))
    m = [np.zeros(n, bool) for n in num.shape]
    for k, n in enumerate(num.shape):
        m[k][0] = True
        if n % 2 == 0:
            m[k][n // 2] = True
    self3 = m[0][:, None, None] & m[1][None, :, None] & m[2][None, None, :]
    S[self3] = num[self3].real
    return S


def fuse(chans, S, how):
    """Combine per-channel (num, J) into an estimate of O~.

    wiener -- MAP with the prior S; blue -- unbiased and prior-free, undefined
    where J = 0; naive -- mean of the per-channel Wiener reconstructions, the
    "average the reconstructions" baseline. Returns (estimate, total J)."""
    if how == "wiener":
        num = sum(c[0] for c in chans)
        J = sum(c[1] for c in chans)
        return num / (J + 1.0 / S), J
    if how == "blue":
        num = sum(c[0] for c in chans)
        J = sum(c[1] for c in chans)
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(J > 0, num / np.where(J > 0, J, 1), 0), J
    if how == "naive":
        ests = [c[0] / (c[1] + 1.0 / S) for c in chans]
        return sum(ests) / len(ests), sum(c[1] for c in chans)
    raise ValueError(how)


def cap_geometry(cfg, grid, ki):
    """Cap geometry for one illumination with transverse wavevector `ki`.

    Unlike `cap_samples`, the illumination is given explicitly (WPM needs it
    snapped to the grid, so the realised inclination differs slightly from the
    nominal one and must be used consistently). Returns a dict with the
    scattered k_z, the exact continuum K_z, the nearest-voxel K_z index, and
    the collected-pupil mask -- everything the deposition needs.
    """
    kix, kiy = float(ki[0]), float(ki[1])
    kiz = np.sqrt(max(cfg.k_m ** 2 - kix ** 2 - kiy ** 2, 0.0))
    kxp, kyp = grid.KX + kix, grid.KY + kiy
    kp2 = kxp ** 2 + kyp ** 2
    ok = kp2 <= cfg.k_perp_max ** 2
    kz = np.sqrt(np.maximum(cfg.k_m ** 2 - kp2, 0.0))
    Kz = np.where(ok, kz - kiz, 0.0)
    iz = np.rint(Kz / grid.dk).astype(np.int64) % grid.nz
    # The centring phase must use the VOXEL's K_z, not the ray's exact one.
    # fftn of a centred array carries e^{-i K.r_c}, and with z_c = 6.4 um the
    # convention flips by ~pi between neighbouring K_z voxels; using the exact
    # K_z here would hand each deposited ray a phase error of up to pi/2 and
    # destroy the reconstruction. Rounding is then applied to the smooth,
    # un-ramped transform, which varies on the scale 2 pi / (object depth).
    return {"ok": ok, "kz": kz, "kiz": kiz, "Kz": Kz, "iz": iz,
            "Kz_vox": grid.dk * np.where(iz < grid.nz // 2 + 1, iz, iz - grid.nz),
            "z_c": (grid.nz // 2) * grid.dz}


def born_rytov_ratio(O, cfg, grid, ki, cap=None):
    """Exact first-Rytov field u/u_i at the object-centre plane, for the SAME
    illumination and pupil as `WPM.run` + `WPM.to_centre`.

    This is the linear control. It evaluates O~ at the *continuum* K_z of each
    collected ray -- no nearest-voxel rounding -- so that when it is pushed
    through the same deposition as the WPM data, the only thing that differs
    between the two is the physics, not the gridding.

        U~_R(K_perp) = -i O~(K) / (2 k_z),   u/u_i = exp(u_R)

    Components outside the collection pupil are dropped, as they are for WPM.
    """
    cap = cap or cap_geometry(cfg, grid, ki)
    Oxy = np.fft.fft2(O, axes=(0, 1))                       # (n, n, nz)
    # O~ at the exact K_z, with the transverse centring ramp cancelling. The
    # z-transform is a polynomial in w = e^{-i K_z dz}, so Horner costs one
    # complex exponential of size n^2 instead of one of size n^2 nz.
    w = np.exp(-1j * cap["Kz"] * grid.dz)
    Ohat = np.zeros(w.shape, complex)
    for k in range(grid.nz - 1, -1, -1):
        Ohat *= w
        Ohat += Oxy[:, :, k]
    Ohat *= np.exp(1j * cap["Kz"] * (grid.nz // 2) * grid.dz)
    with np.errstate(divide="ignore", invalid="ignore"):
        URk = np.where(cap["ok"], -1j * grid.dz * Ohat / (2 * np.maximum(cap["kz"], 1e-12)), 0.0)
    return np.exp(np.fft.ifft2(URk))


def deposit(d2d, cap, grid, use_exact_Kz=False):
    """Collected 2D measurement -> (linear voxel index, value) on the K grid.

    `d2d` is in the units `meas_scale` produces. The e^{-i K_z z_c} factor
    converts the centred-array phase convention of the field to that of
    `fftn(O)`; see the derivation at the head of this section.
    """
    ok = cap["ok"]
    ix, iy = np.nonzero(ok)
    val = d2d[ix, iy] * np.exp(-1j * cap["Kz_vox"][ix, iy] * cap["z_c"])
    return (ix * grid.n + iy) * grid.nz + cap["iz"][ix, iy], val, ix, iy


def shear_rytov(ratio, dx, dvec, wrap_safe=True):
    """Shear difference of the Rytov field, as a shearing interferometer forms it.

    `ratio` is u/u_i at the object-centre plane (so its log is u_R). The shift
    is done in Fourier space: exact for the band-limited periodic field, and
    for non-integer shears alike. u_i cancels from the ratio up to its own
    known, unambiguous phase ramp, which the caller never needs.

    wrap_safe=True  -- returns log[ u(r+D) / u(r) ], whose phase is
        arg[u(r+D) u*(r)]: the principal value of a SMALL quantity. This is
        what the interferogram delivers, and per guizar-sicairos2011 it is
        structurally immune to 2 pi wrapping of the absolute phase.
    wrap_safe=False -- returns log u(r+D) - log u(r) with both logs on the
        principal branch: what you get if the absolute Rytov phase is formed
        first and left unwrapped, i.e. standard Rytov ODT's failure mode
        (chen1998). The difference between the two settings IS the wrapping.
    """
    n = ratio.shape[0]
    k = 2 * np.pi * np.fft.fftfreq(n, dx)
    KX, KY = np.meshgrid(k, k, indexing="ij")
    shifted = np.fft.ifft2(np.fft.fft2(ratio) * np.exp(1j * (KX * dvec[0] + KY * dvec[1])))
    if wrap_safe:
        return np.log(shifted / ratio)
    return np.log(shifted) - np.log(ratio)


# ------------------- exact referee + Born-consistent multiple scattering (exp09+)
#
# exp06 had to splice its multiple-scattering model (D = B + C - B') because
# WPM's linear limit is not first Born. This section adds
#   mie_scalar / born_sphere -- the exact scalar solution for a homogeneous
#       sphere (partial waves), and its first-Born counterpart;
#   CBS -- the convergent Born series of Osnabrugge et al. (2016): an iterative
#       solver of the full scalar Helmholtz equation, forward and backward
#       scattering and evanescent coupling included. The exact referee.
#   MSB -- multi-slice first Born, forward only: as cheap as WPM, but each
#       slice scatters through the exact Green's function -i/(2 k_z) of the
#       SCATTERED wave, so its linear limit is first Born by construction
#       (related to the multi-layer Born model of Chen et al., Optica 7, 394
#       (2020), not ingested here).
# Conventions are those of the exp06 section: O = k0^2 (n_m^2 - n^2), so
# (lap + k_m^2) u = O u, centred arrays, slice k at z_k = (k - nz//2) dz.


def mie_scalar(k, m, a, cos_theta, lmax=None):
    """Far-field amplitude f(theta) of a homogeneous sphere, SCALAR waves.

    Exact partial-wave solution of (lap + k(r)^2) psi = 0 for a sphere of
    radius a and relative index m = n_sphere / n_medium in a medium of
    wavenumber k, with psi and d psi / dr continuous at r = a (the scalar
    analogue of Mie theory; GODT's whole model is scalar). Outside,
    psi = e^{ikz} + f(theta) e^{ikr} / r as r -> inf, with

        f = (-i/k) sum_l (2l+1) b_l P_l(cos theta),
        b_l = [k1 j_l'(k1 a) j_l(k a) - k j_l(k1 a) j_l'(k a)]
              / [k j_l(k1 a) h_l'(k a) - k1 j_l'(k1 a) h_l(k a)],   k1 = m k.

    Checked in selftest against first Born in the weak limit and against the
    optical theorem (sigma_ext = sigma_sca for a lossless sphere).
    """
    from scipy.special import spherical_jn, spherical_yn, eval_legendre
    k1 = m * k
    if lmax is None:
        x = max(k, k1) * a
        lmax = int(x + 4 * x ** (1 / 3) + 10)
    ct = np.asarray(cos_theta, dtype=float)
    f = np.zeros(ct.shape, complex)
    for l in range(lmax + 1):
        j, jp = spherical_jn(l, k * a), spherical_jn(l, k * a, derivative=True)
        j1, j1p = spherical_jn(l, k1 * a), spherical_jn(l, k1 * a, derivative=True)
        h = j + 1j * spherical_yn(l, k * a)
        hp = jp + 1j * spherical_yn(l, k * a, derivative=True)
        b = (k1 * j1p * j - k * j1 * jp) / (k * j1 * hp - k1 * j1p * h)
        f += (2 * l + 1) * b * eval_legendre(l, ct)
    return -1j / k * f


def born_sphere(k, m, a, cos_theta):
    """First-Born far-field amplitude of the same sphere:
    f_B = (k1^2 - k^2) a^3 j1(q a) / (q a),  q = 2 k sin(theta/2)."""
    from scipy.special import spherical_jn
    q = k * np.sqrt(np.maximum(2.0 - 2.0 * np.asarray(cos_theta, float), 0.0))
    qa = np.maximum(q * a, 1e-12)
    return ((m * k) ** 2 - k ** 2) * a ** 3 * spherical_jn(1, qa) / qa


def _kz_prop(KX, KY, k_m):
    kt2 = KX ** 2 + KY ** 2
    return np.sqrt(np.maximum(k_m ** 2 - kt2, 0.0)), kt2 < k_m ** 2


class CBS:
    """Convergent Born series for the scattered field (osnabrugge2016).

    Solves the full scalar Helmholtz equation for the scattered field
    psi_s of a plane wave psi_i = exp(i k_i . r) hitting the object O:

        (lap + k_m^2 - O + i alpha) psi_s = O psi_i,

    on an n x n x nz box that is periodic transversely (k_i must lie on the
    transverse lattice, as for WPM) with 2*nb absorbing slices appended in z.
    alpha (> 0 only in those slices) acts on the scattered field alone, the
    same trick WPM's edge taper uses, so the incident wave is never absorbed.
    Pseudo-spectral Laplacian; no paraxial, forward-only or weak-scattering
    approximation. The periodic transverse box is the SAME one first Born is
    evaluated on (`born_rytov_ratio`), so the two share a window exactly.

    Iteration (the paper's Eqs.): with k(r)^2 = k0^2 + i eps + V(r),
    k0^2 = midpoint of Re k^2, eps >= max |k^2 - k0^2|, gamma = (i/eps) V,
    g0(p) = 1/(|p|^2 - k0^2 - i eps):
        psi <- psi + gamma [ g0 * (V psi + S) - psi ],   S = -O psi_i.
    Each iteration advances the field by about 2 k0 / eps, so cost grows with
    the contrast and with the absorber strength.

    Absorber: Im k^2 = alpha_max * b(u)^2 over the 2 nb appended slices,
    b rising linearly from the object box to the middle of the layer and back
    (the z axis is periodic, so the two halves are one contiguous layer).
    `alpha_max` defaults to the value that attenuates a normally travelling
    wave by exp(-atten) across the whole layer.
    """

    def __init__(self, O, dx, dz, k_m, nb=40, atten=12.0, alpha_max=None,
                 frame=0.0, frame_atten=8.0, dtype=np.complex128):
        import scipy.fft as sfft
        self.sfft = sfft
        O = np.asarray(O)
        self.n, _, self.nz = O.shape
        self.nb, self.dx, self.dz, self.k_m = int(nb), float(dx), float(dz), float(k_m)
        self.dtype = dtype
        nzt = self.nz + 2 * self.nb
        self.nzt = nzt
        L = 2 * self.nb * self.dz
        if alpha_max is None:
            # normal-incidence amplitude decay: int alpha/(2 k_m) dz = alpha_max L / (6 k_m)
            alpha_max = 6.0 * self.k_m * atten / L if self.nb > 0 else 0.0
        self.alpha_max = float(alpha_max)
        u = (np.arange(2 * self.nb) + 0.5) / max(2 * self.nb, 1)
        prof = self.alpha_max * (1.0 - np.abs(2 * u - 1.0)) ** 2
        delta = np.zeros((self.n, self.n, nzt), complex)      # k^2 - k_m^2
        delta[:, :, :self.nz] = -O
        delta[:, :, self.nz:] = 1j * prof[None, None, :]
        # Optional transverse frame: the box stays periodic, but scattered
        # light walking sideways is absorbed in the outer `frame` fraction of
        # each edge instead of circulating. Without it a periodic array of
        # scatterers is simulated, whose near-grazing and guided waves never
        # reach the z layers (slow, algebraic convergence) and are not the
        # physics of an isolated object. Same profile shape as WPM's taper.
        self.frame = float(frame)
        if self.frame > 0:
            xs = (np.arange(self.n) - self.n // 2) * self.dx
            w = self.frame * self.n * self.dx
            half = 0.5 * self.n * self.dx
            t = np.clip((np.abs(xs) - (half - w)) / w, 0.0, 1.0)
            # lateral amplitude decay across one edge for a grazing wave:
            # int alpha/(2 k_m) dx = a_max w / (6 k_m)  (quadratic profile)
            a_xy = 6.0 * self.k_m * frame_atten / w
            Fxy = a_xy * np.maximum(t[:, None] ** 2, t[None, :] ** 2)
            if np.any(O[Fxy > 0, :] != 0):
                raise ValueError("object overlaps the absorbing frame")
            # Where frame and z layers meet the two add. (Taking the max
            # instead lowers eps but left a weakly damped mode in small boxes:
            # 1500 iterations without converging, against ~200 with the sum.)
            delta = delta + 1j * Fxy[:, :, None]
            self.frame_alpha = a_xy
            self.clear_half = half - w
        re = delta.real
        self.shift = 0.5 * (min(re.min(), 0.0) + max(re.max(), 0.0))
        self.k0sq = self.k_m ** 2 + self.shift
        self.eps = max(float(np.abs(delta - self.shift).max()), 1e-6 * self.k_m ** 2)
        self.V = (delta - self.shift - 1j * self.eps).astype(dtype)
        self.gamma = ((1j / self.eps) * self.V).astype(dtype)
        kx = 2 * np.pi * np.fft.fftfreq(self.n, self.dx)
        kz = 2 * np.pi * np.fft.fftfreq(nzt, self.dz)
        p2 = kx[:, None, None] ** 2 + kx[None, :, None] ** 2 + kz[None, None, :] ** 2
        self.g0 = (1.0 / (p2 - self.k0sq - 1j * self.eps)).astype(dtype)
        self.p2 = p2
        self.O = O
        x = (np.arange(self.n) - self.n // 2) * self.dx
        self.x = x
        self.z = (np.arange(self.nz) - self.nz // 2) * self.dz
        self.X, self.Y = np.meshgrid(x, x, indexing="ij")

    def incident(self, ki, z):
        kix, kiy = float(ki[0]), float(ki[1])
        kiz = np.sqrt(self.k_m ** 2 - kix ** 2 - kiy ** 2)
        return np.exp(1j * (kix * self.X + kiy * self.Y))[:, :, None] * \
            np.exp(1j * kiz * np.asarray(z))[None, None, :]

    def solve(self, ki, tol=1e-7, maxiter=5000, verbose=False):
        """Scattered field on the object box (n, n, nz) and an info dict
        (iterations, final update size, Helmholtz residual)."""
        f, fi = self.sfft.fftn, self.sfft.ifftn
        S = np.zeros((self.n, self.n, self.nzt), self.dtype)
        S[:, :, :self.nz] = -self.O * self.incident(ki, self.z)
        g0S = fi(self.g0 * f(S, workers=-1), workers=-1)
        psi = (self.gamma * g0S).astype(self.dtype)          # first iterate
        nS = np.linalg.norm(S)
        it, upd = 1, np.inf
        while it < maxiter:
            d = self.gamma * (fi(self.g0 * f(self.V * psi, workers=-1), workers=-1) + g0S - psi)
            psi += d
            it += 1
            upd = np.linalg.norm(d) / max(np.linalg.norm(psi), 1e-300)
            if verbose and it % 50 == 0:
                print(f"      CBS iteration {it}: update {upd:.2e}")
            if upd < tol:
                break
        # Helmholtz residual of the converged field, inside the object box
        lap = fi(-self.p2 * f(psi, workers=-1), workers=-1)
        k2 = self.k_m ** 2 + (self.V + self.shift + 1j * self.eps)
        res = (lap + k2 * psi + S)[:, :, :self.nz]
        info = {"iterations": it, "update": upd, "eps": self.eps,
                "residual": float(np.linalg.norm(res) / max(nS, 1e-300))}
        return psi[:, :, :self.nz].astype(np.complex128), info

    def ratio(self, psi_s, ki):
        """u / u_i at the object-centre plane, from the field on the top
        slice of the object box: exactly what `WPM.to_centre` returns."""
        zt = self.z[-1]
        U = self.incident(ki, [zt])[:, :, 0] + psi_s[:, :, -1]
        Uc = angular_propagate(U, self.dx, self.k_m, -zt)
        return Uc / self.incident(ki, [0.0])[:, :, 0]


class MSB:
    """Multi-slice first Born, forward only.

    March the total field slice by slice. Slice k (centred at z_k) is lit by
    the field propagated to its centre, and radiates through the exact
    angular-spectrum Green's function of the SCATTERED wave,

        dU~(k_perp) = -i dz / (2 k_z(k_perp)) * FT[O_k U(z_k)],

    so its linear limit is first Born, slice placement included, by
    construction (WPM's is Born times k_z/k_iz). Approximations, all measured
    against CBS in exp09: no backward scattering, no evanescent coupling, a
    slice does not scatter its own scattered field, and components with
    k_z < kz_floor (grazing beyond ~84 deg in GODT_SIM) are not radiated --
    on a periodic window they would circulate sideways for ever.
    Two FFTs per non-empty slice, one complex multiply per empty one.
    """

    def __init__(self, O, dx, dz, k_m, kz_floor_frac=0.1, absorb_frac=0.0,
                 absorb_rate=0.6):
        import scipy.fft as sfft
        self.sfft = sfft
        self.O = np.asarray(O)
        self.n, _, self.nz = self.O.shape
        self.dx, self.dz, self.k_m = float(dx), float(dz), float(k_m)
        k = 2 * np.pi * np.fft.fftfreq(self.n, self.dx)
        KX, KY = np.meshgrid(k, k, indexing="ij")
        self.kz, prop = _kz_prop(KX, KY, self.k_m)
        rad = prop & (self.kz >= kz_floor_frac * self.k_m)
        self.Ph = np.where(prop, np.exp(1j * self.kz * self.dz / 2), 0.0)
        self.P = self.Ph ** 2
        with np.errstate(divide="ignore"):
            self.kern = np.where(rad, -1j * self.dz / (2 * np.maximum(self.kz, 1e-12)), 0.0) * self.Ph
        self.active = [bool(np.any(self.O[:, :, s] != 0)) for s in range(self.nz)]
        x = (np.arange(self.n) - self.n // 2) * self.dx
        self.X, self.Y = np.meshgrid(x, x, indexing="ij")
        self.z = (np.arange(self.nz) - self.nz // 2) * self.dz
        # optional edge taper on the scattered part, identical to WPM's
        self.absorb = None
        if absorb_frac > 0:
            w = absorb_frac * self.n * self.dx
            t = np.clip((np.abs(x) - (0.5 * self.n * self.dx - w)) / w, 0.0, 1.0)
            fe = np.exp(-absorb_rate * self.dz * t ** 2)
            self.absorb = fe[:, None] * fe[None, :]

    def ratio(self, ki):
        """u / u_i at the object-centre plane, as `WPM.to_centre` returns."""
        f, fi = self.sfft.fft2, self.sfft.ifft2
        kix, kiy = float(ki[0]), float(ki[1])
        kiz = np.sqrt(self.k_m ** 2 - kix ** 2 - kiy ** 2)
        ramp = np.exp(1j * (kix * self.X + kiy * self.Y))
        zlo = self.z[0] - self.dz / 2
        Uk = f(ramp * np.exp(1j * kiz * zlo), workers=-1)
        for s in range(self.nz):
            if self.active[s]:
                Uc = fi(Uk * self.Ph, workers=-1)
                Uk = Uk * self.P + self.kern * f(self.O[:, :, s] * Uc, workers=-1)
            else:
                Uk = Uk * self.P
            if self.absorb is not None:
                Ui = ramp * np.exp(1j * kiz * (zlo + (s + 1) * self.dz))
                Uk = f(Ui + (fi(Uk, workers=-1) - Ui) * self.absorb, workers=-1)
        zex = self.z[-1] + self.dz / 2
        back = np.where(self.kz > 0, np.exp(-1j * self.kz * zex), 0.0)
        return fi(Uk * back, workers=-1) / ramp


def lattice_samples(F, n, ix, iy):
    """2D FFT of a centred N x N field (N a multiple of n, same pitch) at the
    frequencies of the n-point lattice, indexed by (ix, iy) on that lattice,
    in the phase convention of the centred n x n crop.

    For a field contained in the crop this equals fft2(crop)[ix, iy] exactly.
    For one that is not -- steep scattered light, refocused to the centre
    plane, spreads tens of microns sideways -- it is the infinite-field-of-
    view value, which is what first Born on a periodic n grid computes
    (Poisson summation). Cropping instead truncates the field before the
    transform, so every model, Born included, loses its steep rays differently;
    exp09 shows that this, not physics, drove exp06's forward-model gaps.
    """
    N = F.shape[0]
    m, r = divmod(N, n)
    assert r == 0, "the big box must be a whole number of small boxes"
    q = np.fft.fftfreq(n, 1.0 / n).astype(np.int64)
    o = (N - n) // 2
    U = np.fft.fft2(F)
    return U[(q[ix] * m) % N, (q[iy] * m) % N] * np.exp(2j * np.pi * (q[ix] + q[iy]) * o / n)


# ------------------------------------- gridding-free reconstruction (exp10+)
#
# exp06's estimator deposits every ray at the NEAREST K_z voxel (assumption
# A3) and that alone costs a factor ~6 in RMSE. The fix needs no interpolation
# kernel: along one K_perp column the continuum datum is an exact DTFT of the
# nz slices of the 2D-FFT'd object, x_k = FFT_xy[O](K_perp, z_k), so
#
#     y_r = h_r c_r sum_k x_k exp(-i K_r (k - nz//2) dz),   c_r = e^{-i K_vox z_c}
#
# holds exactly (c_r is the centring phase exp06's `deposit` applies). The MAP
# estimate under exp06's prior (O~_j independent, variance S_j) solves, per
# column, the nz x nz system
#
#     [ sum_r w_r a_r a_r^H + F^H diag(1/S) F ] x = sum_r (h_r^* y_r / s^2) a_r,
#
# a Toeplitz data term plus a circulant prior term. With K_r = K_vox it
# collapses to exp06's per-voxel  num / (J + 1/S)  exactly (selftest).


def _partner_cols(cols, n):
    ix, iy = np.divmod(cols, n)
    return ((-ix) % n) * n + ((-iy) % n)


class ColumnWiener:
    """Exact (gridding-free) Wiener / fusion estimator over K_perp columns.

    Built once per ray geometry: `cols` (ix * n + iy), continuum `Kz` and the
    centring phase `c` of every ray of every illumination. Hermitian partners
    of a real object are folded in, so each ray also measures column -K_perp
    at -K_z with conjugated coefficients, exactly as `hsym_c` does per voxel.
    Channel data enter through `data_term` (weights) and `rhs` (data), and
    `solve` combines any set of channels with any time shares.
    """

    def __init__(self, grid, cols, Kz, c):
        self.grid, self.n, self.nz, self.dz = grid, grid.n, grid.nz, grid.dz
        cols = np.asarray(cols, np.int64)
        pc = _partner_cols(cols, self.n)
        Kz = np.asarray(Kz, float)
        c = np.asarray(c, complex)
        # a ray whose partner is itself (self-partner column AND K_z = 0)
        # must not be counted twice -- the hsym_c rule for self-partner voxels
        twin = (pc != cols) | (np.abs(Kz) > 1e-12)
        self.cols = np.concatenate([cols, pc[twin]])
        self.K = np.concatenate([Kz, -Kz[twin]])
        self.c = np.concatenate([c, np.conj(c[twin])])
        self.twin = twin
        self.nray = cols.size
        self.ucols, self.inv = np.unique(self.cols, return_inverse=True)
        self.e1 = np.exp(1j * self.K * self.dz)          # one lag step per ray

    def _lagsum(self, wts, lags):
        """sum over each column's rays of wts * e1^lag, for lag in `lags`,
        as a (n_used_cols, len(lags)) array."""
        out = np.empty((self.ucols.size, len(lags)), complex)
        p = wts * self.e1 ** lags[0]
        for j in range(len(lags)):
            out[:, j] = (np.bincount(self.inv, weights=p.real, minlength=self.ucols.size)
                         + 1j * np.bincount(self.inv, weights=p.imag, minlength=self.ucols.size))
            p = p * self.e1
        return out

    def _fold(self, v):
        """Per-ray values of the original rays -> with the partner copies."""
        return np.concatenate([v, v[: self.nray][self.twin]])

    def data_term(self, w):
        """Toeplitz lags t(D) = sum_r w_r e^{i K_r D dz}, D = 0..nz-1, for
        ray weights w_r = |h_r|^2 / s^2 (one channel at unit share)."""
        return self._lagsum(self._fold(np.asarray(w, float)), np.arange(self.nz))

    def rhs(self, g):
        """b(k) = sum_r g_r conj(c_r) e^{i K_r (k - nz//2) dz}, with
        g_r = conj(h_r) y_r / s^2 per original ray (partners get conj(g))."""
        g = np.asarray(g, complex)
        # a self-partner ray (K = 0 in a self-partner column) measures a real
        # number: keep Re(g), as hsym_c keeps Re(num) on self-partner voxels
        g = np.where(self.twin, g, g.real)
        gf = np.concatenate([g, np.conj(g[self.twin])])
        return self._lagsum(gf * np.conj(self.c), np.arange(self.nz) - self.nz // 2)

    def solve(self, T_list, b_list, S, chunk=256):
        """Wiener estimate O~ on the full grid from channels combined as
        sum_j T_j (data term, already scaled by its share) with right-hand
        sides sum_j b_j. Columns no ray reaches are left at zero (no data).

        `b_list` may instead be a list of such lists -- several data sets
        (forward models) seen through the same channels -- which share one
        factorisation; a list of estimates is then returned."""
        n, nz = self.n, self.nz
        multi = not isinstance(b_list[0], np.ndarray)
        sets = b_list if multi else [b_list]
        T = sum(T_list)
        B = np.stack([sum(bl) for bl in sets], axis=-1)          # (cols, nz, nrhs)
        d = np.arange(nz)[:, None] - np.arange(nz)[None, :]
        est = np.zeros((n * n, nz, B.shape[-1]), complex)
        Sc = np.asarray(S).reshape(n * n, nz)[self.ucols]
        pcirc = nz * np.fft.ifft(1.0 / Sc, axis=1)       # circulant prior lags
        for s0 in range(0, self.ucols.size, chunk):
            sl = slice(s0, s0 + chunk)
            Tl = T[sl]
            N = np.where(d[None] >= 0, Tl[:, np.abs(d)], np.conj(Tl[:, np.abs(d)]))
            N = N + pcirc[sl][:, d % nz]
            x = np.linalg.solve(N, B[sl])
            est[self.ucols[sl]] = np.fft.fft(x, axis=1)
        out = [est[:, :, r].reshape(n, n, nz) for r in range(B.shape[-1])]
        return out if multi else out[0]


# ------------------------------------------------------------- utilities
def out_path(fname):
    """Path inside lab/out/ (created on demand). Keeps artifacts together."""
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    d = os.path.join(here, "out")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, fname)


def selftest():
    """Cheap invariants. Run: python shearlib.py"""
    K = np.linspace(-40, 40, 20001)
    for d in (0.1, 0.33, 0.8):
        assert np.max(np.abs(np.abs(H(K, d)) - absH(K, d))) < 1e-12, "closed form"
    # H -> iK as d -> 0
    Kt = np.array([0.5, 5.0, 20.0])
    assert np.max(np.abs(H(Kt, 1e-4) - 1j * Kt) / np.abs(1j * Kt)) < 2e-3, "limit"
    # nulls where predicted
    for d in (0.33, 0.8):
        for nu in nulls_in_band(d, GODT_SIM.kmax):
            assert absH(nu, d) < 1e-9, "null position"
    # published-config expectations, guarding against silent regressions
    assert abs(GODT_SIM.kmax - 31.19) < 0.01, "sim band"
    assert len(nulls_in_band(0.8, GODT_SIM.kmax)) == 3, "0.8 um -> 3 in-band nulls"
    assert len(nulls_in_band(0.1, GODT_SIM.kmax)) == 0, "0.1 um -> none"
    # error formula sanity
    assert abs(wiener_rel_err(2 * np.pi / 0.8, [0.8], 1e6) - 1.0) < 1e-9, "null err=1"

    # ---- 3D Ewald sampling / design (exp03+) ----
    # dimensionless shear power == (d |H|)^2
    Kt = np.linspace(-30, 30, 1001)
    assert np.max(np.abs(G2_vec(Kt, 0 * Kt, (0.8, 0)) - (0.8 * absH(Kt, 0.8)) ** 2)) < 1e-10, "G2"
    # flipK maps K -> -K exactly on an FFT-ordered grid
    g = KGrid(64, 0.1, 24)
    g.check(GODT_SIM)
    KX, KY, KZ = g.K3()
    # (the Nyquist element of an even axis is its own partner, so skip it)
    nyq = (np.abs(KX) >= g.dk * (g.n // 2)) | (np.abs(KY) >= g.dk * (g.n // 2)) | \
          (np.abs(KZ) >= g.kz_max)
    for A_, name in ((KX, "x"), (KY, "y"), (KZ, "z")):
        assert np.allclose(flipK(A_)[~nyq], -A_[~nyq]), f"flipK {name}"
    # P1: the K_z axis is never sampled except at K = 0 (missing cone)
    M = hermitian_sym(cap_density(GODT_SIM, g, n_az=720))
    assert M[0, 0, 1:].max() == 0.0, "P1: axis must be unsampled"
    assert abs(M[0, 0, 0] - 720 * 4 * 180 / 720) < 1e-9, "origin: every cap, weight (1/cos60)^2, once"
    # closed form vs binning: totals and pointwise away from the caustic
    Mc = hermitian_sym(cap_density_closed(GODT_SIM, g))
    assert abs(M.sum() / Mc.sum() - 1) < 0.05, "closed-form total"
    th = np.deg2rad(GODT_SIM.theta_ill)
    KP = np.hypot(KX, KY)
    with np.errstate(divide="ignore", invalid="ignore"):
        cc = -(KX ** 2 + KY ** 2 + KZ ** 2 + 2 * GODT_SIM.k_m * KZ * np.cos(th)) / (
            2 * GODT_SIM.k_m * KP * np.sin(th))
    safe = (np.abs(cc) < 0.8) & (Mc > 0) & (M > 0) & (KP > 3)
    dev = np.median(np.abs(M[safe] / Mc[safe] - 1))
    assert dev < 0.1, f"closed-form pointwise (median dev {dev:.3f})"
    # design solver: 2 channels, brute force vs SLSQP + certificate
    G2 = np.stack([G2_vec(g.KX, g.KY, (0.8, 0)), G2_vec(g.KX, g.KY, (0.55, 0))])
    prob = DesignProblem(100.0 * M, G2, np.ones(2))
    ts = np.linspace(0, 1, 201)
    vals = [prob.f_A(np.array([a, 1 - a]))[0] for a in ts]
    t_bf = ts[int(np.argmin(vals))]
    t, val, gap = solve_design(prob, prob.f_A)
    assert abs(t[0] - t_bf) < 0.01 and gap < 1e-3, f"design solver ({t[0]:.3f} vs {t_bf:.3f}, gap {gap:.1e})"
    # chromatic factor reduces to monochromatic as the bandwidth -> 0
    Ge, G2i = chromatic_G(np.array([5.0, 7.85]), 0.8, 0.53, 1e-6)
    assert np.allclose(np.abs(Ge) ** 2, G2_vec(np.array([5.0, 7.85]), 0.0, (0.8, 0)), atol=1e-6), "chromatic"

    # ---- multiple-scattering forward model (exp06+) ----
    gm = KGrid(64, 0.1, 64)
    cf = GODT_SIM
    KI, th = snap_illumination(cf, gm, 8)
    # snapped k_i must sit exactly on the lattice (that is the whole point)
    assert np.allclose(np.round(KI / gm.dk) * gm.dk, KI, atol=1e-12), "snap on lattice"
    assert np.abs(th - cf.theta_ill).max() < 6.0, "snapped inclination"  # coarse 64 grid
    # free space: WPM must reproduce exact angular-spectrum propagation
    w0 = WPM(np.full((64, 64, 10), cf.n0), gm.dx, gm.dz, cf.lam, cf.n0,
             absorb_frac=0.0, dtype=np.complex128)
    ki = KI[0]
    kiz = np.sqrt(cf.k_m ** 2 - ki[0] ** 2 - ki[1] ** 2)
    U = w0.run(ki)
    exact = np.exp(1j * (ki[0] * w0.X + ki[1] * w0.Y + kiz * (w0.z0 + 10 * gm.dz)))
    assert np.max(np.abs(U - exact)) < 1e-11, "WPM free space"
    # a single slice must equal the analytic thin-screen transmission exactly
    rs = np.random.default_rng(0).standard_normal((64, 64)) * 1e-5
    w1 = WPM((cf.n0 + rs)[:, :, None], gm.dx, gm.dz, cf.lam, cf.n0, n_levels=4096,
             absorb_frac=0.0, dtype=np.complex128)
    nq = cf.n0 + w1.step * np.rint(rs / w1.step)
    k0_ = 2 * np.pi / cf.lam
    t_an = np.exp(1j * np.sqrt(nq ** 2 * k0_ ** 2 - ki[0] ** 2 - ki[1] ** 2) * gm.dz)
    Ui0 = np.exp(1j * (ki[0] * w1.X + ki[1] * w1.Y + kiz * w1.z0))
    assert np.max(np.abs(w1.run(ki) - Ui0 * t_an)) < 1e-12, "WPM thin screen"
    # and that transmission's scattered spectrum is first Born times k_z/k_iz,
    # which is WPM's obliquity factor: taken at the incident, not scattered, ray
    cap = cap_geometry(cf, gm, ki)
    Us = np.fft.fft2(w1.run(ki) / Ui0 - 1.0)
    Ot2 = np.fft.fft2(k0_ ** 2 * (cf.n0 ** 2 - nq ** 2)) * gm.dz
    Born = np.where(cap["ok"], -1j * Ot2 / (2 * np.maximum(cap["kz"], 1e-12)), 0.0)
    sel = cap["ok"] & (np.abs(Born) > 1e-6 * np.abs(Born[cap["ok"]]).max())
    assert abs(np.median(np.abs(Us[sel] / Born[sel]) /
                         (cap["kz"][sel] / cap["kiz"])) - 1) < 1e-3, "WPM obliquity"
    # the shear factor multiplies the measured Rytov field exactly, however it arose
    rr = np.exp(0.3 * (np.random.default_rng(1).standard_normal((64, 64))
                       + 1j * np.random.default_rng(2).standard_normal((64, 64))))
    for d in (0.4, 0.8):
        sh = np.fft.fft2(shear_rytov(rr, gm.dx, (d, 0.0)))
        Gd = np.exp(1j * gm.KX * d) - 1.0
        ref = np.max(np.abs(np.fft.fft2(np.log(rr))))
        assert np.max(np.abs(sh - Gd * np.fft.fft2(np.log(rr)))) < 1e-9 * ref, \
            "shear acts as G on the measured Rytov field"
    # the wrap-safe branch must be the unambiguous one: chaining one-pixel
    # differences has to reproduce a large shear exactly on a smooth field, and
    # to disagree with it on a field whose absolute phase wraps
    rg = np.random.default_rng(7)
    for amp, want_wrap in ((0.01, False), (0.4, True)):
        uu = np.exp(1j * amp * np.cumsum(rg.standard_normal((64, 64)), axis=0))
        s1 = shear_rytov(uu, gm.dx, (gm.dx, 0.0), wrap_safe=True)
        ch = sum(np.roll(s1, -j, axis=0) for j in range(8))
        di = shear_rytov(uu, gm.dx, (8 * gm.dx, 0.0), wrap_safe=True)
        wrapped = np.mean(np.abs(np.imag(ch - di)) > 1.0)
        assert (wrapped > 0.001) == want_wrap, f"wrap detector (amp {amp}: {wrapped})"
        if not want_wrap:
            assert np.max(np.abs(ch - di)) < 1e-10, "chained == direct when unwrapped"

    # ---- wavelength as a shear channel (exp08+) ----
    c2 = at_wavelength(GODT_SIM, 0.25)
    assert abs(c2.kmax - 2 * GODT_SIM.kmax) < 1e-9 and abs(c2.shear - 0.05) < 1e-12, "1/lam scaling"
    # the cap family is a scaled copy: M_lam(K) on a grid of pitch dk equals
    # M_lam0 on the grid scaled by lam0/lam (same indices, M is per voxel)
    g1, g2 = KGrid(64, 0.1, 24), KGrid(64, 0.05, 24)
    assert np.allclose(cap_density(c2, g2, 90), cap_density(GODT_SIM, g1, 90)), "cap scaling"
    assert np.allclose(unmix_shares(np.eye(3), [.25, .5, .25]), [.25, .5, .25]), "mosaic, no crosstalk"
    assert unmix_shares([[1, .2], [.2, 1]], [.5, .5]).max() < 0.5, "crosstalk costs"

    # ---- exact referee + Born-consistent multiple scattering (exp09+) ----
    km_ = GODT_SIM.k_m
    mu = np.linspace(-1.0, 1.0, 8001)
    fm = mie_scalar(km_, 1.38 / 1.33, 1.0, mu)
    sca = 2 * np.pi * np.sum(0.5 * (np.abs(fm[1:]) ** 2 + np.abs(fm[:-1]) ** 2) * np.diff(mu))
    ext = 4 * np.pi / km_ * mie_scalar(km_, 1.38 / 1.33, 1.0, np.array([1.0]))[0].imag
    assert abs(sca / ext - 1) < 1e-4, f"Mie optical theorem ({sca / ext:.6f})"
    fw, fb = mie_scalar(km_, 1.0001, 1.0, mu), born_sphere(km_, 1.0001, 1.0, mu)
    assert np.linalg.norm(fw - fb) / np.linalg.norm(fb) < 5e-3, "Mie -> first Born"
    # MSB: free space is exact, and the linear limit IS first Born on the pupil
    gb = KGrid(64, 0.1, 64)
    kib = snap_illumination(GODT_SIM, gb, 4)[0][1]
    assert np.max(np.abs(MSB(np.zeros(gb.shape), gb.dx, gb.dz, km_).ratio(kib) - 1)) < 1e-12, "MSB free space"
    xb = (np.arange(64) - 32) * 0.1
    Xb, Yb, Zb = np.meshgrid(xb, xb, xb, indexing="ij")
    Ob = 1e-4 * np.exp(-(Xb ** 2 + Yb ** 2 + Zb ** 2) / 2.0) * np.cos(9 * Xb + 4 * Zb)
    capb = cap_geometry(GODT_SIM, gb, kib)
    okb = capb["ok"]
    ub = np.fft.fft2(np.log(MSB(Ob, gb.dx, gb.dz, km_).ratio(kib)))[okb]
    bb = np.fft.fft2(np.log(born_rytov_ratio(Ob, GODT_SIM, gb, kib, capb)))[okb]
    assert np.linalg.norm(ub - bb) / np.linalg.norm(bb) < 1e-5, "MSB linear limit = Born"
    # CBS: a framed weak sphere converges geometrically and solves Helmholtz
    gc = (np.arange(48) - 24) * 0.1
    Xc, Yc, Zc = np.meshgrid(gc, gc, gc[12:36], indexing="ij")
    Oc = np.where(Xc ** 2 + Yc ** 2 + Zc ** 2 <= 0.25, (2 * np.pi / 0.5) ** 2 * (1.33 ** 2 - 1.36 ** 2), 0.0)
    cb = CBS(Oc, 0.1, 0.1, km_, nb=16, atten=6, frame=0.25, frame_atten=6)
    _, info = cb.solve((0.0, 0.0), tol=1e-6, maxiter=1500)
    assert info["iterations"] < 1500 and info["residual"] < 1e-5, f"CBS convergence {info}"
    # full-box lattice sampling == crop FFT for a field the crop contains
    Fq = np.zeros((96, 96), complex)
    Fq[40:56, 38:60] = np.random.default_rng(4).standard_normal((16, 22))
    ixq, iyq = np.nonzero(np.ones((32, 32), bool))
    assert np.allclose(lattice_samples(Fq, 32, ixq, iyq),
                       np.fft.fft2(Fq[32:64, 32:64])[ixq, iyq], atol=1e-10), "lattice_samples"

    # ---- gridding-free reconstruction (exp10+) ----
    # with every ray on its voxel's K_z the column solver IS exp06's per-voxel
    # Wiener rule, Hermitian folding and self-partner voxels included
    g = KGrid(32, 0.1, 32)
    rq = np.random.default_rng(3)
    cols, Kv, cc, lin, hh = [], [], [], [], []
    for ki in snap_illumination(GODT_SIM, g, 6)[0]:
        cp = cap_geometry(GODT_SIM, g, ki)
        ix, iy = np.nonzero(cp["ok"])
        cols.append(ix * g.n + iy)
        Kv.append(cp["Kz_vox"][ix, iy])
        cc.append(np.exp(-1j * cp["Kz_vox"][ix, iy] * cp["z_c"]))
        lin.append((ix * g.n + iy) * g.nz + cp["iz"][ix, iy])
        hh.append(rq.standard_normal(ix.size) + 1j * rq.standard_normal(ix.size))
    cols, Kv, cc, lin, hh = map(np.concatenate, (cols, Kv, cc, lin, hh))
    gy = np.conj(hh) * (rq.standard_normal(hh.size) + 1j * rq.standard_normal(hh.size))
    Sq = rq.uniform(0.5, 2.0, g.shape)
    Sq = 0.5 * (Sq + flipK(Sq))
    numq = (np.bincount(lin, gy.real, g.size) + 1j * np.bincount(lin, gy.imag, g.size)).reshape(g.shape)
    Jq = np.bincount(lin, np.abs(hh) ** 2, g.size).reshape(g.shape)
    ref_pv = fuse([(hsym_c(numq), hermitian_sym(Jq))], Sq, "wiener")[0]
    cw = ColumnWiener(g, cols, Kv, cc)
    est_cw = cw.solve([cw.data_term(np.abs(hh) ** 2)], [cw.rhs(gy)], Sq)
    assert np.abs(est_cw - ref_pv).max() < 1e-12 * np.abs(ref_pv).max(), "ColumnWiener == per-voxel"
    print("shearlib selftest OK")
    print(f"  {GODT_SIM!r}")
    print(f"  {GODT_EXP!r}")


if __name__ == "__main__":
    selftest()
