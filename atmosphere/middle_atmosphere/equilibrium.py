"""Radiative-convective equilibrium of a column with its ozone chemistry, at a fixed surface temperature.

The troposphere follows the saturated pseudoadiabat from the surface with the
Manabe-Wetherald humidity profile; above the tropopause every layer is in
radiative equilibrium (absorbed sunlight equals net thermal emission) and water
vapour keeps its tropopause (cold-trap) mixing ratio. The tropopause is the
lowest level at which radiative equilibrium is stable against the moist
adiabat. Thermal radiation comes from the correlated-k scheme
(atmosphere/radiative_convective/ck.py, with ozone); sunlight below 500 nm from
photolysis.py and from 500 nm to 5 um from the line-by-line solar model with the
ozone Chappuis band added. Ozone and the other trace gases come from the
steady-state chemistry (chemistry.py), and the temperature, chemistry and
sunlight are iterated to mutual consistency.

Holding the surface temperature fixed makes the top-of-atmosphere imbalance a
result: it is the extra forcing (from clouds, surface albedo or greenhouse gases)
the column would need to stay at that temperature. The calculation is a global
and diurnal mean; it has no clouds, circulation, waves or tides, and it treats
the thermal emission of CO2 as in local thermodynamic equilibrium, which fails
above roughly 1-10 Pa.
"""
from __future__ import annotations
from dataclasses import dataclass, field, replace
import math
import numpy as np

from shared.constants import SYNODIC_MONTH_DAYS
from atmosphere.radiative_convective import thermodynamics as th, optics, shortwave as sw, climate as cl, ck
from atmosphere.middle_atmosphere import photolysis as ph, chemistry as ch

SECONDS_PER_DAY = 86400.0
SW_CONFIG = optics.SHORTWAVE          # line-by-line solar configuration (run.py may set its worker count)


@dataclass(frozen=True)
class Case:
    name: str
    planet: th.Planet
    dry_pressure_pa: float
    shield: str = 'none'
    co2_ppm: float = 400.0
    surface_temperature_k: float = 288.0
    surface_albedo: float = 0.13
    n2o_ppb: float = 330.0
    h2_ppm: float = 0.53
    kzz_scale: float = 1.0
    precipitation_m_per_year: float = 1.0
    stratospheric_h2o: float | None = None     # None: the cold-trap value at the tropopause
    top_pressure_pa: float = 0.1
    layers: int = 100
    sw_stride: int = 2                         # line-by-line sunlight on every second level
    day_days: float = SYNODIC_MONTH_DAYS       # length of the solar day (1 for the Earth control)
    nonlte: bool = False                       # CO2 15-um emission out of LTE above 50 Pa (ck.NonLTE)
    nir_thermalisation: str = 'full'           # near-infrared sunlight absorbed above 50 Pa: 'full' heat,
                                               # or 'collisional': the fraction eps/(1+eps) of ck.NonLTE

    def air(self):
        return th.earthlike_air(self.dry_pressure_pa, self.co2_ppm)

    def boundary(self):
        fixed = (('N2O', self.n2o_ppb * 1e-9), ('H2', self.h2_ppm * 1e-6))
        return ch.Boundary(fixed=fixed, precipitation_m_per_year=self.precipitation_m_per_year)

    def mixing(self):
        return ch.Mixing(scale=self.kzz_scale)


class Profile:
    """Levels (surface first) with an adiabatic troposphere below level `trop` and free temperatures above."""

    def __init__(self, case: Case):
        self.case = case
        self.air = case.air()
        self.rh = th.manabe_wetherald()
        ts = case.surface_temperature_k
        self.ps = case.dry_pressure_pa + float(self.rh(np.array([1.0]), 1.0)[0]) * float(th.saturation_pressure(ts))
        self.p = np.exp(np.linspace(math.log(self.ps), math.log(case.top_pressure_pa), case.layers + 1))
        adiabat = th.moist_adiabat(ts, self.ps, self.air, t_stop=80.0, x_max=math.log(self.ps / case.top_pressure_pa) + 0.1)
        self.t_ad = th.adiabat_temperature(adiabat, self.ps, self.p)
        self.t_ad[0] = ts

    def initial(self, t_strat=200.0):
        t = np.maximum(self.t_ad, t_strat)
        trop = int(np.argmax(self.t_ad <= t_strat))
        t[:trop + 1] = self.t_ad[:trop + 1]
        return t, trop

    def water(self, t, trop):
        x = self.rh(self.p, self.ps) * th.saturation_pressure(np.maximum(t, 100.0)) / self.p
        x_trop = self.case.stratospheric_h2o
        if x_trop is None:
            x_trop = float(x[trop])
        x = np.where(np.arange(self.p.size) > trop, x_trop, x)
        return np.clip(x, 0.0, 1.0)

    def column(self, t, trop):
        col = th.column_from_levels(self.case.planet, self.air, self.p, t, self.water(t, trop))
        col['tropopause_pa'] = float(self.p[trop])
        return col


def near_infrared(col, air, o3_layer_columns, shield, surface_albedo, stride=2, n_mu0=6):
    """Line-by-line solar absorption from 500 nm to 5 um (W/m^2): per fine layer, at the surface, and reflected.

    Computed on a column coarsened to every `stride`-th level and interpolated back to the
    fine layers as heating per unit mass (column total conserved).
    """
    idx = np.arange(0, col['p_pa'].size, stride)
    if idx[-1] != col['p_pa'].size - 1:
        idx = np.append(idx, col['p_pa'].size - 1)
    coarse = th.column_from_levels(col['planet'], air, col['p_pa'][idx], col['t_k'][idx], col['x_h2o'][idx])
    nu, tau_abs = optics.optical_depth(coarse, SW_CONFIG)
    # Ozone Chappuis and Wulf absorption on the line-by-line grid (1-nm bins, room temperature).
    sigma_bins = ph.CrossSections().o3_room
    lam = 1e7 / nu
    sigma = np.interp(lam, ph.CENTRES_NM, sigma_bins, left=0.0, right=0.0)
    o3_coarse = np.add.reduceat(o3_layer_columns, idx[:-1])
    tau_abs = tau_abs + o3_coarse[:, None] * sigma[None, :]
    tau_r = optics.rayleigh_optical_depth(coarse, nu)
    spectrum, longward = sw.solar_spectrum(nu, shield=shield)
    radius = (coarse['planet'].radius_m + coarse['z_m'])[::-1]
    # As shortwave.planetary_mean, keeping each zenith-cosine node's local net-flux profile.
    mu, w = sw.disk_nodes(n_mu0)
    weights = np.full(nu.size, nu[1] - nu[0]); weights[0] *= 0.5; weights[-1] *= 0.5
    f = spectrum * weights
    nets = np.zeros((n_mu0, len(idx)))
    reflected = 0.0
    chunk = 40000
    for a in range(0, nu.size, chunk):
        sl = slice(a, a + chunk)
        for i, m in enumerate(mu):
            direct, diffuse, up = sw.column_fluxes(tau_abs[::-1, sl], tau_r[::-1, sl], np.zeros_like(tau_r[::-1, sl]),
                                                   surface_albedo, m, radius)
            nets[i] += ((direct + diffuse - up) @ f[sl])[::-1]          # surface first
            reflected += 0.5 * w[i] * float(up[0] @ f[sl])
    mass = col['layer_mass_kg_m2']
    lp_coarse = np.log(coarse['layer_p_pa'])[::-1]
    lp_fine = np.log(col['layer_p_pa'])[::-1]

    def to_fine(net):
        # Heating per unit mass interpolated smoothly in log pressure, then scaled so the column
        # total is kept. Sharing each coarse layer's absorption evenly within it would give a
        # stepwise heating profile, and in the thin upper layers, where each layer's temperature
        # follows its own heating, those steps force an alternating level temperature pattern.
        coarse_abs = net[1:] - net[:-1]
        specific = coarse_abs / coarse['layer_mass_kg_m2']
        fine = np.interp(lp_fine, lp_coarse, specific[::-1])[::-1] * mass
        total = fine.sum()
        return fine * (coarse_abs.sum() / total) if abs(total) > 0 else fine

    net = 0.5 * (w[:, None] * nets).sum(axis=0)
    fine = to_fine(net)
    nodes = np.array([to_fine(n_) for n_ in nets]).T               # local heating (layers, nodes)
    # Sunlight beyond 5 um (about 1 W/m^2) is absorbed in the lowest layer.
    tail = 0.25 * (longward or 0.0)
    fine[0] += tail
    incident = 0.25 * float(spectrum @ weights)
    return dict(layers=fine, nodes=nodes, surface=float(net[0]), reflected=reflected, incident=incident + tail)


@dataclass
class Result:
    case: Case
    t: np.ndarray
    trop: int
    col: dict
    state: ch.State
    sw_layers: np.ndarray
    lw: dict
    budget: dict
    history: list = field(default_factory=list)
    day_night: dict = field(default_factory=dict)


class Solver:
    def __init__(self, tables=None, verbose=False):
        self.lw = ck.CKLongwave(tables)
        self.verbose = verbose

    def say(self, *a):
        if self.verbose:
            print(*a, flush=True)

    def longwave_absorption(self, col, o3_layers):
        if getattr(self, 'nonlte', False):
            f = self.lw.fluxes_nonlte(col, extra={'O3': o3_layers}, o_mixing=self.o_mixing)
        else:
            f = self.lw.fluxes(col, extra={'O3': o3_layers})
        net_up = f['up'] - f['down']
        return net_up[:-1] - net_up[1:], f

    def heating_rate(self, col, air, sw_layers, o3_layers):
        lw_abs, f = self.longwave_absorption(col, o3_layers)
        cp = col_cp(col, air)
        return (sw_layers + lw_abs) / (col['layer_mass_kg_m2'] * cp) * SECONDS_PER_DAY, f

    def radiative_equilibrium(self, prof, t, trop, sw_layers, o3_mixing, tol=1e-2, max_iter=40, smoothing=1e-3):
        """Solve for the level temperatures above `trop` at fixed absorbed sunlight and ozone mixing ratio.

        Energy balance is posed on layers and temperatures live on levels, so an alternating
        (odd-even) pattern of level temperatures leaves every layer's heating almost unchanged:
        the Jacobian is nearly singular in that one direction and plain Newton steps let the
        pattern grow. Each step therefore solves the least-squares problem |J d + r|^2 +
        s^2 |D2 (t + d)|^2, where D2 is the second difference among the free levels (not
        across the tropopause, whose kink is physical) and s = 1e-3 K/day per K. The pattern's
        singular value is about 1.5e-3 just above the tropopause, where radiative relaxation
        takes 100-400 days, against 4e-3 to 0.8 for the physical modes. Layer means, and so
        every flux, heating rate and chemical rate, do not see the pattern; results are
        reported as layer means. The exact discrete equilibrium carries the pattern (about
        +-5 K in levels just above the tropopause); the smooth profile found here leaves a
        heating residual of up to ~0.01 K/day that alternates in sign between layers and
        cancels over neighbouring pairs, hence the 1e-2 K/day tolerance.
        """
        free = np.arange(trop + 1, t.size)
        layers = np.arange(trop, t.size - 1)
        n = free.size
        d2 = np.zeros((n - 2, n + 1))                          # on [t_trop, free levels]
        for i in range(n - 2):
            d2[i, i + 1:i + 4] = (1.0, -2.0, 1.0)

        def resid(tt):
            col = prof.column(tt, trop)
            q, f = self.heating_rate(col, prof.air, sw_layers, o3_mixing * col['layer_column_cm2'])
            return q[layers], q, col, f

        def jacobian(tt, r):
            jac = np.empty((layers.size, n))
            for j, lev in enumerate(free):
                tp = tt.copy(); tp[lev] += 0.5
                jac[:, j] = (resid(tp)[0] - r) / 0.5
            return jac

        r, q, col, f = resid(t)
        jac = None
        self.jacobian = (trop, None)
        best, stalled = np.inf, 0
        for it in range(max_iter):
            size = np.max(np.abs(r))
            if size < tol:
                break
            stalled = 0 if size < 0.8 * best else stalled + 1
            best = min(best, size)
            if jac is None or stalled >= 3:
                stalled = 0
                jac = jacobian(t, r)
            a = np.vstack([jac, smoothing * d2[:, 1:]])
            b = np.concatenate([-r, -smoothing * (d2 @ t[trop:])])
            step = np.linalg.lstsq(a, b, rcond=None)[0]
            scale = min(1.0, 15.0 / max(np.max(np.abs(step)), 1e-12))
            for _ in range(6):
                tn = t.copy(); tn[free] = np.clip(t[free] + scale * step, 60.0, 400.0)
                rn, qn, coln, fn = resid(tn)
                if np.max(np.abs(rn)) < np.max(np.abs(r)) or scale < 0.05:
                    break
                scale *= 0.5
            t, r, q, col, f = tn, rn, qn, coln, fn
            self.say(f'    RE iteration {it}: max |heating| {np.max(np.abs(r)):.2e} K/day')
        if jac is None:
            jac = jacobian(t, r)
        self.jacobian = (trop, jac)
        return t, q, col, f

    def day_night(self, prof, t, trop, col, local_nodes, n_mu0=6, samples=720):
        """Linear, radiative-only response of the levels above the tropopause to the day at the equator.

        The local solar heating at the Gauss zenith-cosine nodes (layers x nodes, W/m^2) is
        interpolated along the equatorial day (mu0 = cos of the hour angle, zero at night);
        its first harmonic at the day length drives dT/dt = J T' + Q'(t), with J the
        radiative-equilibrium Jacobian (radiative exchange between levels included). Returns
        the amplitude of the first harmonic and the local radiative relaxation time.
        """
        t_jac, jac = self.jacobian
        if jac is None or t_jac != trop:
            raise RuntimeError('day_night needs the Jacobian of the final radiative-equilibrium solve')
        mu, _ = sw.disk_nodes(n_mu0)
        phase = 2 * np.pi * (np.arange(samples) + 0.5) / samples
        mu_t = np.maximum(np.cos(phase), 0.0)
        xs = np.concatenate([[0.0], mu])
        layers = np.arange(trop, t.size - 1)
        heat = np.array([np.interp(mu_t, xs, np.concatenate([[0.0], local_nodes[k]])) for k in layers])
        cp = col_cp(col, prof.air)[layers]
        rate = heat / (col['layer_mass_kg_m2'][layers] * cp)[:, None] * SECONDS_PER_DAY      # K/day
        q1 = 2.0 * (rate * np.exp(-1j * phase)[None, :]).mean(axis=1)
        omega = 2 * np.pi / prof.case.day_days
        response = np.linalg.solve(1j * omega * np.eye(layers.size) - jac, q1)
        return dict(p_pa=prof.p[trop + 1:], amplitude_k=np.abs(response),
                    relaxation_days=-1.0 / np.minimum(np.diag(jac), -1e-12),
                    mean_heating_k_day=rate.mean(axis=1), peak_heating_k_day=rate.max(axis=1))

    def solve(self, case: Case, max_outer=10, state=None):
        prof = Profile(case)
        t, trop = prof.initial()
        self.nonlte, self.o_mixing = case.nonlte, None
        shield = cl.load_shield(case.shield)
        sun = ph.Sunlight(shield)
        history = []
        o3_mixing = None
        sw_layers = None
        for outer in range(max_outer):
            col = prof.column(t, trop)
            state = ch.solve(col, sun, case.surface_albedo, col['tropopause_pa'], case.mixing(), case.boundary(),
                             state=state)
            m = ch.fixed_densities(col)['M']
            new_o3 = state.n[ch.INDEX['O3']] / m
            du = ch.column_du(col, state.n)
            o3_layers = new_o3 * col['layer_column_cm2']
            self.o_mixing = state.n[ch.INDEX['O']] / m
            nir = near_infrared(col, prof.air, o3_layers, shield, case.surface_albedo, stride=case.sw_stride)
            if case.nir_thermalisation == 'collisional':
                # Only the collisionally deactivated fraction of absorbed near-infrared sunlight heats the air
                # above 50 Pa; the rest is re-emitted. A lower bound, since some energy reaches heat by other paths.
                eps = ck.NonLTE().epsilon(col, self.o_mixing)
                keep = np.where(col['layer_p_pa'] < ck.NonLTE().start_pa, eps / (1 + eps), 1.0)
                nir['layers'] = nir['layers'] * keep
                nir['nodes'] = nir['nodes'] * keep[:, None]
            sw_layers = state.heating + nir['layers']
            local_nodes = state.info['light']['heating_nodes'] + nir['nodes']
            t_old = t.copy()
            # Radiative equilibrium above the tropopause, then move the tropopause until it is consistent.
            seen = set()
            while True:
                t, q, col, f = self.radiative_equilibrium(prof, t, trop, sw_layers, new_o3)
                seen.add(trop)
                # Stability of the lowest radiative layer against the adiabat, judged on layer means.
                if (trop + 1 < t.size - 1 and t[trop] + t[trop + 1] < prof.t_ad[trop] + prof.t_ad[trop + 1] - 0.1
                        and trop + 1 not in seen):
                    trop += 1
                    t[:trop + 1] = prof.t_ad[:trop + 1]
                elif trop > 1 and q[trop - 1] > 0 and trop - 1 not in seen:
                    trop -= 1
                    t[:trop + 1] = prof.t_ad[:trop + 1]
                else:
                    break
            d_t = float(np.max(np.abs(t - t_old)))
            # Relative to the peak mixing ratio, with a 10-ppb floor so an ozone-free column converges.
            d_o3 = 1.0 if o3_mixing is None else float(np.max(np.abs(new_o3 - o3_mixing)) / max(np.max(new_o3), 1e-8))
            o3_mixing = new_o3
            history.append(dict(outer=outer, max_dt_k=d_t, o3_change=d_o3, o3_du=du, trop=trop,
                                tropopause_pa=float(prof.p[trop]), t_top=float(t[-1])))
            self.say(f'  outer {outer}: dT {d_t:.2f} K, dO3 {d_o3:.3f}, O3 {du:.1f} DU, tropopause '
                     f'{prof.p[trop]:.0f} Pa ({t[trop]:.1f} K), T(top) {t[-1]:.1f} K')
            if d_t < 0.3 and d_o3 < 0.01:
                break
        col = prof.column(t, trop)
        state = ch.solve(col, sun, case.surface_albedo, col['tropopause_pa'], case.mixing(), case.boundary(),
                         state=state)
        lw_abs, f = self.longwave_absorption(col, state.n[ch.INDEX['O3']] / ch.fixed_densities(col)['M']
                                             * col['layer_column_cm2'])
        light = state.info['light']
        uv_below = ph.CENTRES_NM < ph.HEATING_LIMIT_NM
        budget = dict(olr=f['olr'], surface_down_lw=f['surface_down'],
                      asr_uv=float(state.info['light']['toa_net'][uv_below].sum()),
                      asr_nir=nir['incident'] - nir['reflected'],
                      surface_sw=float(light['surface_net'][uv_below].sum()) + nir['surface'])
        budget['asr'] = budget['asr_uv'] + budget['asr_nir']
        budget['net_toa'] = budget['asr'] - budget['olr']
        cycle = self.day_night(prof, t, trop, col, local_nodes)
        return Result(case=case, t=t, trop=trop, col=col, state=state, sw_layers=sw_layers,
                      lw=dict(absorbed=lw_abs, fluxes=f), budget=budget, history=history, day_night=cycle)


def col_cp(col, air):
    """Specific heat of moist air in each layer (J/kg/K)."""
    x = col['layer_x_h2o']
    m_dry, m_w = air.molar_mass, th.MOLAR_MASS['H2O']
    q = x * m_w / (x * m_w + (1 - x) * m_dry)
    return (1 - q) * air.cp + q * th.CP_VAPOUR
