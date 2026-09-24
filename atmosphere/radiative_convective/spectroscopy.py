"""Molecular absorption and scattering cross-sections on a uniform wavenumber grid.

- Lines: HITRAN line parameters with the standard temperature and pressure
  scaling, Voigt profiles cut 25 cm-1 from line centre. Water lines have the
  profile value at 25 cm-1 subtracted (the pedestal), which is the convention the
  MT_CKD continuum assumes. Profiles are summed exactly (Voigt) within about
  1 cm-1 of each centre, on a ten-times coarser mesh in the near wings, where the
  Lorentzian equals the Voigt profile, and by FFT convolution beyond a few cm-1;
  the tests bound the error of this multi-mesh line-by-line method.
- Water-vapour continuum: MT_CKD 4.3 self and foreign coefficients with AER's
  density scaling, temperature dependence, radiation term and four-point
  interpolation.
- Collision-induced absorption: HITRAN CIA tables, linear in temperature between
  tabulated sets, clamped outside them.
- Rayleigh scattering of air: Bodhaine et al. (1999).

Cross-sections are cm^2 per molecule; CIA coefficients cm^5 per molecule^2.
"""
from __future__ import annotations
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from scipy.io import netcdf_file
from scipy.signal import fftconvolve
from scipy.special import voigt_profile

from shared.constants import AVOGADRO, BOLTZMANN, PLANCK, SPEED_OF_LIGHT
from atmosphere.radiative_convective import fetch_inputs

C2 = PLANCK * SPEED_OF_LIGHT / BOLTZMANN * 100.0   # second radiation constant, cm K
T_REF = 296.0
VOIGT_CORE = 50.0   # Gaussian standard deviations within which the full Voigt profile is used
KERNEL_VERSION = 3  # bump when line_absorption changes, so cached cross-sections are rebuilt
CACHE = Path(__file__).resolve().parent / 'cache'

# (HITRAN molecule id, local isotopologue code) -> (global id, molar mass g/mol),
# from https://hitran.org/docs/iso-meta/ (retrieved 2026-09-24).
ISOTOPOLOGUES = {
    (1, '1'): (1, 18.010565), (1, '2'): (2, 20.014811), (1, '3'): (3, 19.014780),
    (1, '4'): (4, 19.016740), (1, '5'): (5, 21.020985), (1, '6'): (6, 20.020956),
    (1, '7'): (129, 20.022915),
    (2, '1'): (7, 43.989830), (2, '2'): (8, 44.993185), (2, '3'): (9, 45.994076),
    (2, '4'): (10, 44.994045), (2, '5'): (11, 46.997431), (2, '6'): (12, 45.997400),
    (2, '7'): (13, 47.998320), (2, '8'): (14, 46.998291), (2, '9'): (121, 45.998262),
    (2, '0'): (15, 49.001675), (2, 'A'): (120, 48.001646), (2, 'B'): (122, 47.001618),
    (3, '1'): (16, 47.984745), (3, '2'): (17, 49.988991), (3, '3'): (18, 49.988991),
    (3, '4'): (19, 48.988960), (3, '5'): (20, 48.988960),
    (7, '1'): (36, 31.989830), (7, '2'): (37, 33.994076), (7, '3'): (38, 32.994045),
}

LINE_FILES = {
    'H2O': ('hitran_h2o_0-3500.par', 'hitran_h2o_3500-25000.par'),
    'CO2': ('hitran_co2_0-3500.par', 'hitran_co2_3500-10000.par'),
    'O2': ('hitran_o2_0-20000.par',),
    'O3': ('hitran_o3_0-3500.par',),
}


@dataclass(frozen=True)
class Grid:
    """Uniform wavenumber grid, cm^-1."""
    start: float
    step: float
    size: int

    @classmethod
    def span(cls, start, stop, step):
        return cls(float(start), float(step), int(round((stop - start) / step)) + 1)

    @property
    def nu(self):
        return self.start + self.step * np.arange(self.size)

    @property
    def stop(self):
        return self.start + self.step * (self.size - 1)


# ---------------------------------------------------------------- line data

@dataclass(frozen=True)
class LineList:
    molecule: str
    nu: np.ndarray
    strength: np.ndarray      # cm/molecule at 296 K, abundance-weighted (HITRAN)
    gamma_air: np.ndarray     # cm^-1/atm HWHM at 296 K
    gamma_self: np.ndarray
    elower: np.ndarray        # cm^-1
    n_air: np.ndarray
    delta_air: np.ndarray     # cm^-1/atm
    gid: np.ndarray           # HITRAN global isotopologue id
    mass: np.ndarray          # g/mol

    def select(self, lo, hi):
        m = (self.nu >= lo) & (self.nu <= hi)
        return LineList(self.molecule, *(getattr(self, f)[m] for f in
                        ('nu', 'strength', 'gamma_air', 'gamma_self', 'elower', 'n_air', 'delta_air', 'gid', 'mass')))

    def __len__(self):
        return len(self.nu)


def _parse_par(path: Path, s_min: float):
    fields = {k: [] for k in ('nu', 'strength', 'gamma_air', 'gamma_self', 'elower', 'n_air', 'delta_air', 'gid', 'mass')}
    with open(path, 'r') as handle:
        for row in handle:
            if len(row) < 67:
                continue
            s = float(row[15:25])
            if s < s_min:
                continue
            gid, mass = ISOTOPOLOGUES[(int(row[0:2]), row[2])]
            fields['nu'].append(float(row[3:15])); fields['strength'].append(s)
            fields['gamma_air'].append(float(row[35:40])); fields['gamma_self'].append(float(row[40:45]))
            fields['elower'].append(float(row[45:55])); fields['n_air'].append(float(row[55:59]))
            fields['delta_air'].append(float(row[59:67])); fields['gid'].append(gid); fields['mass'].append(mass)
    return {k: np.asarray(v, dtype=np.int64 if k == 'gid' else float) for k, v in fields.items()}


def load_lines(molecule: str, s_min: float) -> LineList:
    """Lines of a molecule with 296-K intensity >= s_min, cached per input hash."""
    manifest = {f['name']: f['sha256'] for f in fetch_inputs.manifest()['files']}
    names = LINE_FILES[molecule]
    key = hashlib.sha256(json.dumps([[manifest[n] for n in names], s_min]).encode()).hexdigest()[:16]
    cache = CACHE / f'lines_{molecule}_{key}.npz'
    if cache.is_file():
        data = dict(np.load(cache))
    else:
        parts = [_parse_par(fetch_inputs.path(n), s_min) for n in names]
        data = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
        order = np.argsort(data['nu'], kind='stable')
        data = {k: v[order] for k, v in data.items()}
        CACHE.mkdir(exist_ok=True)
        np.savez(cache, **data)
    # HITRAN gives no air broadening (0.0000) for 271 weak minor-isotopologue O3
    # lines near 1020 cm-1 (4e-5 of the O3 band strength); they take the
    # molecule's median value so every line has a finite Lorentz width.
    missing = data['gamma_air'] <= 0
    if missing.any():
        data['gamma_air'] = np.where(missing, np.median(data['gamma_air'][~missing]), data['gamma_air'])
    return LineList(molecule, *(data[f] for f in ('nu', 'strength', 'gamma_air', 'gamma_self', 'elower',
                                                   'n_air', 'delta_air', 'gid', 'mass')))


class PartitionSums:
    """HITRAN total internal partition sums, tabulated at 1-K steps."""
    def __init__(self, gids):
        self.tables = {}
        for gid in sorted(set(int(g) for g in gids)):
            data = np.loadtxt(fetch_inputs.path(f'hitran_q{gid}.txt'))
            self.tables[gid] = (data[:, 0], data[:, 1])

    def q(self, gid, t):
        tt, qq = self.tables[int(gid)]
        return float(np.interp(t, tt, qq))

    def ratio(self, gids, t):
        """Q(296 K)/Q(T) for each line."""
        out = np.empty(len(gids))
        for gid in np.unique(gids):
            out[gids == gid] = self.q(gid, T_REF) / self.q(gid, t)
        return out


def _ramp(d, inner, outer):
    """Smooth step from 0 at |d| <= inner to 1 at |d| >= outer."""
    x = np.clip((np.abs(d) - inner) / (outer - inner), 0.0, 1.0)
    return x * x * (3 - 2 * x)


def line_absorption(grid: Grid, lines: LineList, sums: PartitionSums, t: float, p_pa: float,
                    p_self_pa: float = 0.0, pedestal: bool = False, cutoff: float = 25.0,
                    fine_halfwidth: float = 1.0, coarse_factor: int = 10, far_start: float = 4.0,
                    far_full: float = 6.0, batch: int = 4000):
    """Absorption cross-section (cm^2/molecule) of one molecule's lines at (T, p).

    Three parts per line, which sum to its cut-off Voigt profile:
    - far wing (smoothly from far_start to the cutoff): the Lorentzian expanded
      as (gamma/pi)(1/d^2 - gamma^2/d^4), summed for all lines at once by FFT
      convolution on the coarse mesh; the neglected term is (gamma/d)^4 < 3e-4;
    - near wing (within far_full): the Lorentzian on the coarse mesh;
    - core (within fine_halfwidth plus one coarse step): the exact profile on
      the fine mesh minus the coarse interpolant of the same line; the Voigt
      function is evaluated within VOIGT_CORE Gaussian widths of centre and
      the Lorentzian beyond, where the two agree to 3/VOIGT_CORE^2.
    """
    lines = lines.select(grid.start - cutoff, grid.stop + cutoff)
    out_fine = np.zeros(grid.size)
    h, cf = grid.step, int(coarse_factor)
    big = h * cf
    ncoarse = (grid.size - 1) // cf + 2
    coarse = np.zeros(ncoarse)
    if len(lines) == 0:
        return out_fine
    if not fine_halfwidth + big < far_start < far_full < cutoff:
        raise ValueError('Need core window < far_start < far_full < cutoff')
    p_atm, ps_atm = p_pa / 101325.0, p_self_pa / 101325.0
    q_ratio = sums.ratio(lines.gid, t)
    strength = (lines.strength * q_ratio * np.exp(-C2 * lines.elower * (1 / t - 1 / T_REF))
                * (-np.expm1(-C2 * lines.nu / t)) / (-np.expm1(-C2 * lines.nu / T_REF)))
    centre = lines.nu + lines.delta_air * (p_atm - ps_atm)
    lorentz = (T_REF / t) ** lines.n_air * (lines.gamma_air * (p_atm - ps_atm) + lines.gamma_self * ps_atm)
    gauss = lines.nu * math.sqrt(BOLTZMANN * t * AVOGADRO * 1e3) / SPEED_OF_LIGHT / np.sqrt(lines.mass)
    # Beyond ~1 cm^-1 from centre the Voigt profile equals the Lorentzian to within
    # (Doppler width / distance)^2 < 1e-5, so the wings use the Lorentzian.
    base = lorentz / (math.pi * (cutoff**2 + lorentz**2)) if pedestal else np.zeros(len(lines))

    # Far wings: sticks split linearly between the two nearest coarse nodes.
    kfar = int(math.ceil(cutoff / big))
    pad = kfar + 2
    pos = (centre - grid.start) / big
    k0 = np.floor(pos).astype(np.int64)
    frac = pos - k0
    nstick = ncoarse + 2 * pad
    sticks = []
    for weight in (strength * lorentz / math.pi, strength * lorentz**3 / math.pi, strength * base):
        a = np.bincount(k0 + pad, weights=(1 - frac) * weight, minlength=nstick)
        a += np.bincount(k0 + 1 + pad, weights=frac * weight, minlength=nstick)
        sticks.append(a[:nstick])
    d = np.arange(-kfar, kfar + 1) * big
    ramp = np.where(np.abs(d) <= cutoff, _ramp(d, far_start, far_full), 0.0)
    safe = np.where(d == 0, 1.0, d)
    kernels = (ramp / safe**2, ramp / safe**4, ramp)
    far = (fftconvolve(sticks[0], kernels[0]) - fftconvolve(sticks[1], kernels[1])
           - fftconvolve(sticks[2], kernels[2]))
    coarse += far[kfar + pad: kfar + pad + ncoarse]

    window = fine_halfwidth + big
    kf = int(math.ceil(window / h))
    kc = int(math.ceil(far_full / big))
    kw = int(math.ceil(window / big)) + 2
    off_f = np.arange(-kf, kf + 1)
    off_c = np.arange(-kc, kc + 1)
    off_w = np.arange(0, 2 * kw + 1)
    for a in range(0, len(lines), batch):
        sl = slice(a, a + batch)
        s_, c_, l_, g_, b_ = strength[sl, None], centre[sl, None], lorentz[sl, None], gauss[sl, None], base[sl, None]

        def wing(dnu):
            return s_ * (l_ / (math.pi * (dnu * dnu + l_ * l_)) - b_) * (1 - _ramp(dnu, far_start, far_full))

        def profile(dnu):
            shape = l_ / (math.pi * (dnu * dnu + l_ * l_))
            core = np.abs(dnu) < VOIGT_CORE * g_
            if core.any():
                gb, lb = np.broadcast_to(g_, dnu.shape)[core], np.broadcast_to(l_, dnu.shape)[core]
                shape[core] = voigt_profile(dnu[core], gb, lb)
            return s_ * (shape - b_)

        # Near wings on the coarse mesh.
        jc = np.rint((c_[:, 0] - grid.start) / big).astype(np.int64)
        jj = jc[:, None] + off_c[None, :]
        vals = wing(grid.start + jj * big - c_)
        ok = (jj >= 0) & (jj < ncoarse)
        coarse += np.bincount(np.where(ok, jj, 0).ravel(), weights=np.where(ok, vals, 0.0).ravel(), minlength=ncoarse)
        # The same line on the coarse points around its centre, for the correction.
        jw0 = jc - kw
        wvals = wing(grid.start + (jw0[:, None] + off_w[None, :]) * big - c_)
        # Core on the fine mesh: exact profile minus its coarse interpolant.
        ic = np.rint((c_[:, 0] - grid.start) / h).astype(np.int64)
        ii = ic[:, None] + off_f[None, :]
        nu_f = grid.start + ii * h
        exact = profile(nu_f - c_)
        jf = np.floor_divide(ii, cf)
        fr = (ii - jf * cf) / cf
        q0 = jf - jw0[:, None]
        q0c = np.clip(q0, 0, 2 * kw - 1)
        interp = ((1 - fr) * np.take_along_axis(wvals, q0c, axis=1)
                  + fr * np.take_along_axis(wvals, q0c + 1, axis=1))
        interp = np.where((q0 >= 0) & (q0 < 2 * kw), interp, 0.0)
        near = np.abs(nu_f - c_) <= window
        okf = near & (ii >= 0) & (ii < grid.size)
        out_fine += np.bincount(np.where(okf, ii, 0).ravel(), weights=np.where(okf, exact - interp, 0.0).ravel(),
                                minlength=grid.size)
    idx = np.arange(grid.size)
    j = idx // cf
    fr = (idx - j * cf) / cf
    out_fine += (1 - fr) * coarse[j] + fr * coarse[np.minimum(j + 1, ncoarse - 1)]
    return np.maximum(out_fine, 0.0)


# ---------------------------------------------------------- water continuum

class WaterContinuum:
    """MT_CKD 4.3 water-vapour continuum (Mlawer et al. 2012, 2023), AER coefficients."""
    def __init__(self):
        with netcdf_file(fetch_inputs.path('mt_ckd_4.3_absco-ref_wv-mt-ckd.nc'), 'r', mmap=False) as f:
            self.wn = f.variables['wavenumbers'].data.astype(float).copy()
            self.self_ref = f.variables['self_absco_ref'].data.astype(float).copy()
            self.for_ref = f.variables['for_absco_ref'].data.astype(float).copy()
            self.self_texp = f.variables['self_texp'].data.astype(float).copy()
            self.ref_press_mb = float(f.variables['ref_press'].data)
            self.ref_temp = float(f.variables['ref_temp'].data)
        self.dv = self.wn[1] - self.wn[0]

    @staticmethod
    def radiation_term(nu, t):
        """AER radiation field term nu tanh(c2 nu / 2T), cm^-1."""
        x = C2 * nu / t
        return np.where(x <= 0.01, 0.5 * x * nu, nu * np.tanh(0.5 * x))

    def _xint(self, coeff, nu):
        """AER four-point interpolation (myxint) from the 10 cm^-1 coefficient grid."""
        j = np.floor((nu - self.wn[0]) / self.dv).astype(np.int64)
        j = np.clip(j, 1, len(self.wn) - 3)
        p = (nu - self.wn[j]) / self.dv
        c = (3 - 2 * p) * p * p
        b = 0.5 * p * (1 - p)
        b1, b2 = b * (1 - p), b * p
        return -coeff[j - 1] * b1 + coeff[j] * (1 - c + b2) + coeff[j + 1] * (c + b1) - coeff[j + 2] * b2

    def cross_sections(self, nu, t, p_pa, x_h2o):
        """Self and foreign continuum cross-sections per water molecule (cm^2)."""
        rho = (p_pa / 100.0 / self.ref_press_mb) * (self.ref_temp / t)
        rad = self.radiation_term(self.wn, t)
        self_c = self.self_ref * (self.ref_temp / t) ** self.self_texp * x_h2o * rho * rad
        for_c = self.for_ref * (1 - x_h2o) * rho * rad
        nu = np.asarray(nu, dtype=float)
        inside = (nu >= self.wn[1]) & (nu <= self.wn[-3])
        return (np.where(inside, self._xint(self_c, nu), 0.0), np.where(inside, self._xint(for_c, nu), 0.0))


# ---------------------------------------------------- collision-induced absorption

class CollisionInduced:
    """One HITRAN CIA file: bands of tabulated k(nu) at several temperatures.

    Sets from one reference, or complementary in temperature, form a band and
    are interpolated linearly in temperature. Where a later band from a
    different reference overlaps an earlier one at the same temperatures, the
    earlier band is used and the later one only fills wavenumbers it does not
    cover, so no absorption is counted twice.
    """
    def __init__(self, name):
        rows = fetch_inputs.path(name).read_text().splitlines()
        sets, i = [], 0
        while i < len(rows):
            head = rows[i]
            n = int(head[40:47])
            data = np.array([[float(v) for v in r.split()[:2]] for r in rows[i + 1:i + 1 + n]])
            sets.append(dict(lo=float(head[20:30]), hi=float(head[30:40]), t=float(head[47:54]),
                             ref=head[-3:].strip(), nu=data[:, 0], k=np.maximum(data[:, 1], 0.0)))
            i += n + 1
        bands = []
        for s in sets:
            for b in bands:
                overlap = min(b['hi'], s['hi']) - max(b['lo'], s['lo'])
                if overlap <= 0.5 * min(b['hi'] - b['lo'], s['hi'] - s['lo']):
                    continue
                t_overlap = min(b['tmax'], s['t']) >= max(b['tmin'], s['t']) and s['ref'] not in b['refs']
                if s['ref'] in b['refs'] or not t_overlap:
                    b['sets'].append(s); b['refs'].add(s['ref'])
                    b['lo'], b['hi'] = min(b['lo'], s['lo']), max(b['hi'], s['hi'])
                    b['tmin'], b['tmax'] = min(b['tmin'], s['t']), max(b['tmax'], s['t'])
                    break
            else:
                bands.append(dict(sets=[s], refs={s['ref']}, lo=s['lo'], hi=s['hi'], tmin=s['t'], tmax=s['t']))
        for b in bands:
            seen, uniq = set(), []
            for s in sorted(b['sets'], key=lambda s: s['t']):
                if s['t'] not in seen:
                    seen.add(s['t']); uniq.append(s)
            b['sets'] = uniq
        # Earlier bands take precedence where later bands overlap them.
        for n_, b in enumerate(bands):
            b['exclude'] = [(e['lo'], e['hi']) for e in bands[:n_]
                            if min(e['hi'], b['hi']) > max(e['lo'], b['lo'])]
        self.name, self.bands = name, bands

    def coefficient(self, nu, t):
        nu = np.asarray(nu, dtype=float)
        out = np.zeros_like(nu)
        for b in self.bands:
            sets = b['sets']
            ts = np.array([s['t'] for s in sets])
            tc = min(max(t, ts[0]), ts[-1])
            k = int(np.searchsorted(ts, tc, side='right')) - 1
            k = min(max(k, 0), len(ts) - 1)
            def on(s):
                return np.interp(nu, s['nu'], s['k'], left=0.0, right=0.0)
            if k == len(ts) - 1 or ts[k] == tc:
                val = on(sets[k])
            else:
                w = (tc - ts[k]) / (ts[k + 1] - ts[k])
                val = (1 - w) * on(sets[k]) + w * on(sets[k + 1])
            for lo, hi in b['exclude']:
                val = np.where((nu >= lo) & (nu <= hi), 0.0, val)
            out += val
        return out


CIA_PAIRS = (('N2', 'N2', 'hitran_cia_N2-N2_2021.cia'), ('O2', 'N2', 'hitran_cia_O2-N2_2024.cia'),
             ('O2', 'O2', 'hitran_cia_O2-O2_2024.cia'), ('N2', 'H2O', 'hitran_cia_N2-H2O_2018.cia'))


# ------------------------------------------------------------------- Rayleigh

def rayleigh_cross_section(nu, co2_fraction=4e-4):
    """Rayleigh scattering cross-section of dry air (cm^2/molecule), Bodhaine et al. (1999).

    Peck & Reeder refractivity of standard air scaled for CO2 (their eq. 19) and
    the depolarisation (King) factor of the N2/O2/Ar/CO2 mixture (eq. 23).
    Composition differences between Earth and the lunar case (a few per cent
    less O2) change this by under one per cent and are ignored.
    """
    lam_um = 1e4 / np.asarray(nu, dtype=float)
    inv2 = lam_um ** -2
    n300 = 1e-8 * (8060.51 + 2480990.0 / (132.274 - inv2) + 17455.7 / (39.32957 - inv2))
    n = 1 + n300 * (1 + 0.54 * (co2_fraction - 0.0003))
    c_pct = co2_fraction * 100
    f_n2 = 1.034 + 3.17e-4 * inv2
    f_o2 = 1.096 + 1.385e-3 * inv2 + 1.448e-4 * inv2 ** 2
    king = (78.084 * f_n2 + 20.946 * f_o2 + 0.934 * 1.0 + c_pct * 1.15) / (78.084 + 20.946 + 0.934 + c_pct)
    ns = 2.546899e19
    lam_cm = lam_um * 1e-4
    return 24 * math.pi ** 3 * (n * n - 1) ** 2 / (lam_cm ** 4 * ns ** 2 * (n * n + 2) ** 2) * king
