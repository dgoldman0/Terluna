"""Single-column model of the month-long lunar day: ground or sea, boundary layer and air column.

    python -m climate.lunar_day                 # every case; writes climate/results/lunar_day/
    python -m climate.lunar_day --only moon_titania_equator_land earth_equator_land

A stopped run resumes where it left off: finished cases are skipped, and an
unfinished case restarts from the end of its last simulated day (checkpoints in
results/lunar_day/checkpoints/, removed when the case finishes).

A column of air over land or sea is followed through repeated solar days (29.5
days on the Moon, one day for the Earth control) until it repeats, and the
last day is recorded. It starts from a radiative-convective equilibrium of
atmosphere/middle_atmosphere (composition, ozone and temperatures at a 288 K
surface), regridded so the lowest kilometres are resolved: layers 10 m thick at
the ground, growing upward.

Each 30-minute step (10 minutes for the Earth control):

- Thermal radiation from the correlated-k scheme (atmosphere/radiative_convective/ck.py)
  with the ground's own temperature for its emission, recomputed every 2 hours
  on the Moon (every 30 minutes for the Earth control). Sunlight absorbed in each
  layer and at the ground is interpolated in the solar zenith cosine from a
  table computed once, line by line from 500 nm to 5 um and with the
  ultraviolet-visible model below 500 nm, for the case's shield and
  composition.
- The ground exchanges heat with the lowest layer by bulk formulas with a
  stability correction (Louis 1979) and a prescribed wind, in 60-s substeps.
  Land conducts heat from a massless skin into 24 soil layers 10 m deep (the
  top one 2.6 mm thick); the sea is a mixed layer.
  Evaporation meets an aerodynamic and a surface resistance in series (100 s/m
  for well-watered vegetation, 10,000 s/m for arid ground, none for the sea);
  dew meets only the first.
- Water vapour is followed in the lowest 3 km; vapour beyond saturation
  condenses, releasing its heat, and is removed (dew, fog or drizzle, whose
  optical effects are not modelled). Higher up, relative humidity is capped at
  the Manabe-Wetherald profile: vapour that convection carries beyond it
  condenses there, releasing its heat, and falls out as rain (clouds are not
  modelled).
- Unstable layers mix (dry convective adjustment in potential temperature, with
  their vapour); a weak background diffusivity mixes the lowest 3 km.
- Exchange with the rest of the Moon is a relaxation toward the global-mean
  equilibrium. It stands in for the circulation; it is the model's main
  assumption, and the cases vary it. Above 3 km, temperature relaxes over
  12 hours. The Moon turns so slowly that its Rossby radius far exceeds its
  size, so gravity waves (about 50 m/s, as on Earth) erase temperature
  differences in the free atmosphere in about half a day, as they do in Earth's
  tropics: the weak-temperature-gradient approximation (Sobel and Bretherton
  2000; for slowly rotating planets, Mills and Abbot 2013). Heat that
  convection and condensation deliver above 3 km is thereby carried away.
  Vapour relaxes over 20 days, below 3 km toward the global-mean humidity and
  above it toward the reference relative humidity, so drier air regains vapour
  no faster than the circulation brings it.

A few simulated days cannot bring 10 m of soil or 50 m of sea into their
periodic state, so between days (not before the recorded one) the model moves
them there. Every soil layer is shifted so its mean over the day equals the
skin's: the day-mean temperature of a periodically heated conducting ground is
the same at every depth. The sea is held to the reference surface temperature
on average: the heat that currents would have to bring or carry away is set to
the day's mean uptake, and is reported. Otherwise an equatorial sea would keep
warming for many lunar days, because the single column has no way to export the
heat that ocean and air currents carry poleward. The air of the lowest 3 km is
not relaxed, so the column still takes several days to settle, longest near
the pole and over dry ground: each Moon case runs eight lunar days and the
Earth control 60 days. ground_drift_k reports the change left between the
last two days.

The model has no clouds (except that condensed water is removed), no horizontal
advection beyond the relaxation, no topography, a single prescribed wind and
no vegetation or soil maps. It bounds how far the long day and night can swing
temperatures at a place and shows which processes set the swing. It does not
predict a lunar climate.
"""
from __future__ import annotations
import argparse
import csv
from dataclasses import dataclass, field, replace
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np

from shared.constants import SYNODIC_MONTH_DAYS, STEFAN_BOLTZMANN
from atmosphere.radiative_convective import thermodynamics as th, optics, shortwave as sw, ck, climate as cl
from atmosphere.middle_atmosphere import equilibrium as eq, photolysis as ph, run as ma_run, chemistry as ch

HERE = Path(__file__).resolve().parent
RESULTS = HERE / 'results' / 'lunar_day'
MA_RESULTS = HERE.parent / 'atmosphere' / 'middle_atmosphere' / 'results'
DAY_S = 86400.0
MU_NODES = np.array([0.02, 0.05, 0.1, 0.17, 0.25, 0.35, 0.5, 0.65, 0.8, 0.92, 1.0])
BL_HEIGHTS = np.array([0, 10, 25, 50, 100, 175, 275, 400, 575, 800, 1100, 1500, 2000, 2700, 3500, 4500, 6000,
                       8000, 10500, 13500, 17000, 21000], dtype=float)
BL_TOP_M = 3000.0
KARMAN = 0.4


@dataclass(frozen=True)
class Surface:
    kind: str = 'land'                   # 'land' or 'sea'
    roughness_m: float = 0.1             # land 0.1 m (short vegetation); sea 1e-4 m
    surface_resistance_s_m: float = 100.0   # to evaporation: well-watered vegetation ~100, dry soil ~500, sea 0
    soil_diffusivity_m2_s: float = 5e-7
    soil_heat_capacity_j_m3_k: float = 2.0e6
    mixed_layer_m: float = 50.0


@dataclass(frozen=True)
class Case:
    name: str
    equilibrium: str                     # middle-atmosphere case supplying the column
    latitude_deg: float = 0.0
    surface: Surface = field(default_factory=Surface)
    wind_m_s: float = 5.0
    stable_tail: str = 'short'           # 'short': 1/(1+10 Ri)^2; 'long': 1/(1+10 Ri)
    background_k_m2_s: float = 1.0
    relaxation_days: float = 0.5           # temperature above 3 km (weak temperature gradient)
    vapour_relaxation_days: float = 20.0   # vapour below 3 km
    days: int = 8                        # solar days simulated; the last is recorded
    step_s: float = 1800.0
    radiation_every: int = 4             # thermal radiation recomputed every this many steps


SEA = Surface(kind='sea', roughness_m=1e-4, surface_resistance_s_m=0.0)


def cases():
    base = Case('moon_titania_equator_land', 'moon_1.2atm_titania_stack')
    out = [base,
           replace(base, name='moon_titania_equator_sea', surface=SEA),
           replace(base, name='moon_titania_lat60_land', latitude_deg=60.0),
           replace(base, name='moon_titania_lat60_sea', latitude_deg=60.0, surface=SEA),
           replace(base, name='moon_titania_lat80_land', latitude_deg=80.0),
           replace(base, name='moon_edge200_equator_land', equilibrium='moon_1.2atm_edge_200nm'),
           replace(base, name='moon_titania_equator_land_calm', wind_m_s=2.0),
           replace(base, name='moon_titania_equator_land_windy', wind_m_s=8.0),
           replace(base, name='moon_titania_equator_land_long_tail', stable_tail='long'),
           replace(base, name='moon_titania_equator_land_weak_mixing', background_k_m2_s=0.1),
           replace(base, name='moon_titania_equator_land_dry_soil', surface=replace(Surface(), surface_resistance_s_m=500.0)),
           replace(base, name='moon_titania_equator_land_arid', surface=replace(Surface(), surface_resistance_s_m=1e4)),
           replace(base, name='moon_titania_equator_land_exchange_3h', relaxation_days=0.125),
           replace(base, name='moon_titania_equator_land_exchange_5d', relaxation_days=5.0),
           # The 20-day vapour relaxation needs about 50 Earth days to settle the column.
           Case('earth_equator_land', 'earth_control', days=60, step_s=600.0, radiation_every=3),
           Case('earth_equator_sea', 'earth_control', surface=SEA, days=60, step_s=600.0, radiation_every=3),
           Case('earth_equator_land_arid', 'earth_control', surface=replace(Surface(), surface_resistance_s_m=1e4),
                days=60, step_s=600.0, radiation_every=3)]
    return {c.name: c for c in out}


def _profile(name):
    with open(MA_RESULTS / 'profiles' / f'{name}.csv', newline='') as handle:
        rows = list(csv.DictReader(handle))
    return {k: np.array([float(r[k]) if r[k] != '' else np.nan for r in rows]) for k in rows[0]}


class Column:
    """The equilibrium column regridded with a resolved boundary layer, and its sunlight table."""

    def __init__(self, case: Case):
        self.case = case
        self.eq_case = ma_run.cases()[case.equilibrium]
        self.day_s = self.eq_case.day_days * DAY_S
        prof = eq.Profile(self.eq_case)
        self.air = prof.air
        self.planet = self.eq_case.planet
        src = _profile(case.equilibrium)
        # Heights of the equilibrium layers; new levels: resolved boundary layer, then the equilibrium levels.
        lp_src = np.log(src['p_pa'])
        z_src = src['z_km'] * 1e3
        ps = prof.ps
        p_bl = np.exp(np.interp(BL_HEIGHTS[1:], np.concatenate([[0.0], z_src]),
                                np.concatenate([[math.log(ps)], lp_src])))
        upper = prof.p[prof.p < p_bl[-1] * 0.999]
        self.p = np.concatenate([[ps], p_bl, upper])
        lp_mid = np.log(np.sqrt(self.p[1:] * self.p[:-1]))
        # The surface closes the table, so layers below the lowest equilibrium layer (3.6 km up on the Moon)
        # follow the profile down to the reference surface temperature and humidity.
        self.ts_ref = self.eq_case.surface_temperature_k
        x_surface = float(prof.rh(np.array([ps]), ps)[0] * th.saturation_pressure(self.ts_ref) / ps)
        lp_tab = np.concatenate([[math.log(ps)], lp_src])
        interp = lambda values, surface: np.interp(-lp_mid, -lp_tab, np.concatenate([[surface], values]))
        self.t_ref = interp(src['t_k'], self.ts_ref)
        self.x_ref = interp(src['h2o'], x_surface)
        self.mix = {s: interp(src[s], src[s][0]) for s in ('O3', 'NO2', 'NO3', 'N2O5', 'N2O', 'HNO3', 'H2O2')}
        self.trop_p = json.loads((MA_RESULTS / 'cases' / f'{case.equilibrium}.json').read_text())['tropopause_pa']
        self.albedo = self.eq_case.surface_albedo
        col = self.column(self.t_ref, self.x_ref)
        self.mass = col['layer_mass_kg_m2']
        self.heights = col['layer_z_m']
        self.thickness = np.diff(col['z_m'])
        self.bl = col['layer_z_m'] < BL_TOP_M
        self._solar_table(col)

    def o3_columns(self, col):
        return self.mix['O3'] * col['layer_column_cm2']

    def levels(self, t_layers):
        lp = np.log(self.p)
        lp_mid = np.log(np.sqrt(self.p[1:] * self.p[:-1]))
        inner = np.interp(-lp[1:-1], -lp_mid, t_layers)
        return np.concatenate([[t_layers[0]], inner, [t_layers[-1]]])

    def column(self, t_layers, x_layers):
        t_lev = self.levels(t_layers)
        x_lev = np.concatenate([[x_layers[0]], 0.5 * (x_layers[1:] + x_layers[:-1]), [x_layers[-1]]])
        return th.column_from_levels(self.planet, self.air, self.p, t_lev, np.clip(x_lev, 0.0, 0.5))

    def _solar_table(self, col):
        """Layer and ground absorption (W/m^2) of sunlight at each zenith-cosine node."""
        m = ch.fixed_densities(col)['M']
        dens = {s: v * m for s, v in self.mix.items()}
        shield = cl.load_shield(self.eq_case.shield)
        sun = ph.Sunlight(shield)
        tau_uv, sca_uv, _ = sun.optics(col, dens)
        below = ph.CENTRES_NM < ph.HEATING_LIMIT_NM
        nu, tau = optics.optical_depth(col, eq.SW_CONFIG)
        tau = tau + self.o3_columns(col)[:, None] * np.interp(1e7 / nu, ph.CENTRES_NM, ph.CrossSections().o3_room,
                                                              left=0.0, right=0.0)[None, :]
        ray = optics.rayleigh_optical_depth(col, nu)
        spectrum, longward = sw.solar_spectrum(nu, shield=shield)
        weights = np.full(nu.size, nu[1] - nu[0]); weights[0] *= 0.5; weights[-1] *= 0.5
        radius = (self.planet.radius_m + col['z_m'])[::-1]
        nl = col['layer_p_pa'].size
        self.table_layers = np.zeros((MU_NODES.size, nl))
        self.table_ground = np.zeros(MU_NODES.size)
        for i, mu in enumerate(MU_NODES):
            d, f, u = sw.column_fluxes(tau_uv[::-1][:, below], sca_uv[::-1][:, below],
                                       np.zeros_like(sca_uv[::-1][:, below]), self.albedo, mu, radius)
            net = ((d + f - u) @ sun.energy[below])[::-1]
            for a in range(0, nu.size, 40000):
                sl = slice(a, a + 40000)
                d, f, u = sw.column_fluxes(tau[::-1, sl], ray[::-1, sl], np.zeros_like(ray[::-1, sl]), self.albedo,
                                           mu, radius)
                net = net + ((d + f - u) @ (spectrum[sl] * weights[sl]))[::-1]
            self.table_layers[i] = net[1:] - net[:-1]
            self.table_ground[i] = net[0]
            self.table_layers[i, 0] += (longward or 0.0) * mu        # beyond 5 um: absorbed near the ground
        self.incident_top = float(spectrum @ weights) + (longward or 0.0)

    def sunlight(self, mu0):
        """Layer and ground absorption for a zenith cosine, linear between nodes and to zero at the horizon."""
        if mu0 <= 0:
            return np.zeros(self.table_layers.shape[1]), 0.0
        xs = np.concatenate([[0.0], MU_NODES])
        layers = np.vstack([np.zeros(self.table_layers.shape[1]), self.table_layers])
        ground = np.concatenate([[0.0], self.table_ground])
        i = min(int(np.searchsorted(xs, mu0, side='right')) - 1, xs.size - 2)
        w = (mu0 - xs[i]) / (xs[i + 1] - xs[i])
        return (1 - w) * layers[i] + w * layers[i + 1], float((1 - w) * ground[i] + w * ground[i + 1])


def saturation_mole_fraction(t, p):
    return np.minimum(th.saturation_pressure(np.clip(t, 150.0, 400.0)) / p, 1.0)


def stability(ri, tail):
    ri = np.asarray(ri, dtype=float)
    stable = 1.0 / (1 + 10 * np.maximum(ri, 0.0)) ** (2 if tail == 'short' else 1)
    unstable = 1 + 10 * np.sqrt(np.maximum(-ri, 0.0))
    return np.where(ri > 0, stable, np.minimum(unstable, 4.0))


def _tridiagonal(a, b, c, d):
    n = b.size
    cp, dp = np.empty(n), np.empty(n)
    cp[0], dp[0] = c[0] / b[0], d[0] / b[0]
    for i in range(1, n):
        den = b[i] - a[i] * cp[i - 1]
        cp[i] = c[i] / den if i < n - 1 else 0.0
        dp[i] = (d[i] - a[i] * dp[i - 1]) / den
    x = np.empty(n)
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


class Model:
    def __init__(self, case: Case, column: Column):
        self.case, self.c = case, column
        self.lw = ck.CKLongwave()
        self.t = column.t_ref.copy()
        self.x = column.x_ref.copy()
        s = case.surface
        if s.kind == 'land':
            edges = np.concatenate([[0.0], 0.01 * (1.35 ** np.arange(25) - 1) / 0.35])
            self.dz_soil = np.diff(edges * 10.0 / edges[-1])
            self.soil = np.full(self.dz_soil.size, column.ts_ref)
        else:
            self.soil = np.array([column.ts_ref])
        self.ocean_heat = 0.0                    # W/m^2 that currents bring to the sea (negative: carry away)
        self.condensed = 0.0                     # kg/m^2 of vapour condensed in the air since the last record
        self.cp_air = column.air.cp
        self.r_air = column.air.gas_constant
        self.g = column.planet.surface_gravity
        self.kappa = self.r_air / self.cp_air

    @property
    def ts(self):
        return float(self.soil[0])

    def mu0(self, time_s):
        lat = math.radians(self.case.latitude_deg)
        hour = 2 * math.pi * time_s / self.c.day_s
        return max(0.0, math.cos(lat) * math.cos(hour))

    def ground_heat_capacity(self):
        s = self.case.surface
        return s.mixed_layer_m * 1025.0 * 3990.0 if s.kind == 'sea' else None

    def exchange(self, dt, lw_down, sw_ground):
        """Ground-air exchange in 60-s substeps; returns mean sensible and latent fluxes (W/m^2) and dew (kg/m^2).

        Each substep solves the ground implicitly with its fluxes linearised about the current skin
        temperature, then the lowest layer implicitly toward the new ground state, then mixes any
        instability upward, so heat entering the thin lowest layer reaches the mixed layer at once.
        """
        c, s = self.c, self.case.surface
        p1 = math.sqrt(c.p[0] * c.p[1])
        z1 = max(c.heights[0], 1.0)                  # mid-height of the lowest layer
        m1 = c.mass[0]
        cn = (KARMAN / math.log(max(z1, 2 * s.roughness_m) / s.roughness_m)) ** 2
        sub = 60.0
        n = int(round(dt / sub))
        h_sum = le_sum = dew = 0.0
        for _ in range(n):
            ts, ta = self.ts, self.t[0]
            exner = (c.p[0] / p1) ** self.kappa
            theta_a = ta * exner
            rho = p1 / (self.r_air * ta)
            u = max(self.case.wind_m_s, 0.5)
            ri = self.g * z1 * (theta_a - ts) / (0.5 * (theta_a + ts) * u * u)
            conduct = rho * cn * float(stability(ri, self.case.stable_tail)) * u    # kg m^-2 s^-1
            lv = float(th.latent_heat(max(ts, 150.0)))
            q_sat = 0.622 * float(saturation_mole_fraction(ts, c.p[0]))
            q_air = 0.622 * self.x[0]
            # Evaporation meets the surface resistance in series with the air's; dew only the air's.
            r_air = 1.0 / max(conduct / rho, 1e-9)
            wet = conduct if q_sat <= q_air else rho / (r_air + s.surface_resistance_s_m)
            e = wet * (q_sat - q_air)
            h = conduct * self.cp_air * (ts - theta_a)
            le = lv * e
            emit = STEFAN_BOLTZMANN * ts ** 4
            net = sw_ground + lw_down - emit - h - le
            dnet = 4 * STEFAN_BOLTZMANN * ts ** 3 + conduct * self.cp_air + wet * lv * q_sat * lv / (th.R_VAPOUR * ts * ts)
            if s.kind == 'sea':
                cap = self.ground_heat_capacity()
                self.soil[0] += (net + self.ocean_heat) * sub / (cap + dnet * sub)
            else:
                self._conduct(net, sub, dnet)
            # Lowest layer, implicitly toward the new ground state (potential temperature and vapour).
            a = conduct * sub / m1
            theta_new = (theta_a + a * self.ts) / (1 + a)
            self.t[0] = theta_new / exner
            q_sat_new = 0.622 * float(saturation_mole_fraction(self.ts, c.p[0]))
            aw = wet * sub / m1
            q_new = (q_air + aw * q_sat_new) / (1 + aw)
            e = (q_new - q_air) * m1 / sub
            self.x[0] = max(q_new / 0.622, 0.0)
            if e < 0:
                dew -= e * sub
            h_sum += conduct * self.cp_air * (self.ts - theta_new)
            le_sum += lv * e
            self.adjust(boundary_layer_only=True)
        return h_sum / n, le_sum / n, dew

    def _conduct(self, flux, dt, dflux_dt):
        """Implicit soil conduction; the surface flux linearised about the current skin temperature."""
        s = self.case.surface
        k = s.soil_diffusivity_m2_s * s.soil_heat_capacity_j_m3_k
        dz = self.dz_soil
        cap = s.soil_heat_capacity_j_m3_k * dz
        gap = 0.5 * (dz[1:] + dz[:-1])
        cond = k / gap
        n = dz.size
        a = np.zeros(n); b = cap / dt; c = np.zeros(n); d = cap / dt * self.soil
        a[1:] = -cond; c[:-1] = -cond
        b = b.copy(); b[1:] += cond; b[:-1] += cond
        b[0] += dflux_dt
        d = d.copy(); d[0] += flux + dflux_dt * self.soil[0]
        self.soil = _tridiagonal(a, b, c, d)

    def adjust(self, boundary_layer_only=False):
        """Dry convective adjustment in potential temperature (mixing vapour too), then condensation."""
        c = self.c
        p_mid = np.sqrt(c.p[1:] * c.p[:-1])
        theta = self.t * (c.p[0] / p_mid) ** self.kappa
        m = c.mass
        top = int(np.argmax(p_mid < c.trop_p))
        if boundary_layer_only:
            top = min(top, int(np.argmax(~c.bl)) + 1)
        # Mix each unstable pair until none is left; mixing conserves the pair's heat (and vapour).
        for _ in range(4 * top * top):
            unstable = np.where(theta[1:top] < theta[:top - 1] - 1e-6)[0]
            if unstable.size == 0:
                break
            for k in unstable:
                sl = slice(k, k + 2)
                theta[sl] = (theta[sl] * m[sl]).sum() / m[sl].sum()
                self.x[sl] = (self.x[sl] * m[sl]).sum() / m[sl].sum()
        self.t = theta * (p_mid / c.p[0]) ** self.kappa
        sat = saturation_mole_fraction(self.t, p_mid)
        over = np.where(c.bl, np.maximum(self.x - sat, 0.0), 0.0)
        if np.any(over > 0):
            lv = th.latent_heat(np.clip(self.t, 150.0, 400.0))
            dq = 0.622 * over
            self.t = self.t + lv * dq / self.cp_air
            self.x = self.x - over
            self.condensed += float((dq * m).sum())

    def diffuse(self, dt):
        """Background mixing of potential temperature and vapour through the lowest 3 km (implicit)."""
        c = self.c
        k_bg = self.case.background_k_m2_s
        idx = np.where(c.bl)[0]
        if k_bg <= 0 or idx.size < 2:
            return
        p_mid = np.sqrt(c.p[1:] * c.p[:-1])
        z = c.heights[idx]
        m = c.mass[idx]
        rho = m / c.thickness[idx]
        coef = 0.5 * (rho[1:] + rho[:-1]) * k_bg / np.diff(z)      # kg m^-2 s^-1 between layer centres
        n = idx.size
        for arr, is_t in ((self.t, True), (self.x, False)):
            v = arr[idx] * ((c.p[0] / p_mid[idx]) ** self.kappa if is_t else 1.0)
            a = np.zeros(n); cc = np.zeros(n); b = m / dt
            a[1:] = -coef; cc[:-1] = -coef
            b = b.copy(); b[1:] += coef; b[:-1] += coef
            new = _tridiagonal(a, b, cc, m / dt * v)
            arr[idx] = new * ((p_mid[idx] / c.p[0]) ** self.kappa if is_t else 1.0)

    def relax(self, dt):
        c = self.c
        above = ~c.bl
        # Exact exponential decay, stable for time scales shorter than the step.
        keep_t = math.exp(-dt / (self.case.relaxation_days * DAY_S))
        keep_x = math.exp(-dt / (self.case.vapour_relaxation_days * DAY_S))
        self.t[above] = c.t_ref[above] + (self.t[above] - c.t_ref[above]) * keep_t
        self.x[c.bl] = c.x_ref[c.bl] + (self.x[c.bl] - c.x_ref[c.bl]) * keep_x
        p_mid = np.sqrt(c.p[1:] * c.p[:-1])
        # Above the boundary layer, relative humidity is capped at the reference profile's: vapour that
        # convection carried beyond it condenses there and falls out, releasing its heat. Drier air
        # regains vapour only through the exchange, never faster than the circulation would bring it.
        rh_ref = c.x_ref / saturation_mole_fraction(c.t_ref, p_mid)
        target = np.minimum(rh_ref[above] * saturation_mole_fraction(self.t[above], p_mid[above]), 0.5)
        dq = 0.622 * np.maximum(self.x[above] - target, 0.0)
        self.t[above] += th.latent_heat(np.clip(self.t[above], 150.0, 400.0)) * dq / self.cp_air
        self.condensed += float((dq * c.mass[above]).sum())
        capped = np.minimum(self.x[above], target)
        self.x[above] = target + (capped - target) * keep_x

    def spin_up(self, start, mean):
        """Between days, move the ground toward its periodic state (see the module notes)."""
        if self.case.surface.kind == 'land':
            self.soil += mean[0] - mean
        else:
            self.ocean_heat -= self.ground_heat_capacity() * (self.soil[0] - start[0]) / self.c.day_s
            self.soil[0] += self.c.ts_ref - mean[0]

    def run(self, record_last=True, checkpoint=None):
        """Simulate the case's days and return the last one's series. With checkpoint (a path) the
        state is saved after every day, and a run that was stopped resumes from its last saved day."""
        c, case = self.c, self.case
        dt = case.step_s
        steps_per_day = int(round(c.day_s / dt))
        series = []
        self.day_means = []
        first = self._resume(checkpoint) if checkpoint else 0
        for day in range(first, case.days):
            start = self.soil.copy()
            soil_sum = np.zeros_like(self.soil)
            for i in range(steps_per_day):
                time_s = (i + 0.5) * dt
                mu0 = self.mu0(time_s)
                if i % case.radiation_every == 0:
                    col = c.column(self.t, self.x)
                    f = self.lw.fluxes(col, t_surface=self.ts, extra={'O3': c.mix['O3'] * col['layer_column_cm2']})
                    net_up = f['up'] - f['down']
                    lw_layers = net_up[:-1] - net_up[1:]
                sw_layers, sw_ground = c.sunlight(mu0)
                self.t += dt * (sw_layers + lw_layers) / (self.cp_air * c.mass)
                h, le, dew = self.exchange(dt, f['surface_down'], sw_ground)
                self.diffuse(dt)
                self.adjust()
                self.relax(dt)
                cond, self.condensed = self.condensed, 0.0
                if record_last and day == case.days - 1:
                    p_mid = np.sqrt(c.p[1:] * c.p[:-1])
                    theta = self.t * (c.p[0] / p_mid) ** self.kappa
                    mixed = c.heights[np.argmax(theta > theta[0] + 1.0)]
                    rh0 = self.x[0] / float(saturation_mole_fraction(self.t[0], p_mid[0]))
                    series.append(dict(hour=time_s / 3600.0, mu0=mu0, ground_k=self.ts, air_k=float(self.t[0]),
                                       air_300m_k=float(np.interp(300.0, c.heights, self.t)),
                                       column_mean_k=float((self.t * c.mass).sum() / c.mass.sum()),
                                       sensible_w_m2=h, latent_w_m2=le, sw_ground_w_m2=sw_ground,
                                       lw_down_w_m2=f['surface_down'], olr_w_m2=f['olr'],
                                       mixed_layer_m=float(mixed), rh_ground=rh0, dew_kg_m2=dew,
                                       condensed_kg_m2=cond))
                soil_sum += self.soil
            mean = soil_sum / steps_per_day
            self.day_means.append(float(mean[0]))
            if day < case.days - 1:
                self.spin_up(start, mean)
                if checkpoint:
                    self._save(checkpoint, day + 1)
        return series

    def _state_key(self):
        """Identifies the case, this code and the column, so a checkpoint from anything else is ignored."""
        h = hashlib.sha256(repr(self.case).encode() + Path(__file__).read_bytes())
        for a in (self.c.t_ref, self.c.x_ref, self.c.table_layers, self.c.table_ground):
            h.update(np.ascontiguousarray(a, dtype=float).tobytes())
        return h.hexdigest()

    def _save(self, path, day):
        path.parent.mkdir(parents=True, exist_ok=True)
        part = path.with_name(path.stem + '.part.npz')
        np.savez(str(part), key=self._state_key(), day=day, t=self.t, x=self.x, soil=self.soil,
                 ocean_heat=self.ocean_heat, day_means=np.array(self.day_means))
        part.replace(path)

    def _resume(self, path):
        """Restore the state saved after the last finished day; returns the day to start from."""
        if not path.is_file():
            return 0
        with np.load(path) as saved:
            if str(saved['key']) != self._state_key():
                return 0
            self.t, self.x, self.soil = saved['t'].copy(), saved['x'].copy(), saved['soil'].copy()
            self.ocean_heat = float(saved['ocean_heat'])
            self.day_means = saved['day_means'].tolist()
            return int(saved['day'])


def summarise(case, series, column, model):
    arr = {k: np.array([r[k] for r in series]) for k in series[0]}
    night = arr['mu0'] <= 0
    out = dict(case=case.name, equilibrium=case.equilibrium, latitude_deg=case.latitude_deg, surface=case.surface.kind,
               wind_m_s=case.wind_m_s, stable_tail=case.stable_tail, background_k_m2_s=case.background_k_m2_s,
               relaxation_days=case.relaxation_days, vapour_relaxation_days=case.vapour_relaxation_days,
               surface_resistance_s_m=case.surface.surface_resistance_s_m,
               day_length_days=column.day_s / DAY_S)
    for k in ('ground_k', 'air_k', 'air_300m_k', 'column_mean_k'):
        out[k.replace('_k', '') + '_max_k'] = float(arr[k].max())
        out[k.replace('_k', '') + '_min_k'] = float(arr[k].min())
        out[k.replace('_k', '') + '_mean_k'] = float(arr[k].mean())
    out['ground_range_k'] = out['ground_max_k'] - out['ground_min_k']
    out['air_range_k'] = out['air_max_k'] - out['air_min_k']
    out['max_mixed_layer_m'] = float(arr['mixed_layer_m'].max())
    out['night_hours_saturated'] = float((arr['rh_ground'][night] > 0.999).sum() * case.step_s / 3600.0)
    out['night_dew_kg_m2'] = float(arr['dew_kg_m2'][night].sum())
    out['mean_sensible_w_m2'] = float(arr['sensible_w_m2'].mean())
    out['mean_latent_w_m2'] = float(arr['latent_w_m2'].mean())
    out['mean_olr_w_m2'] = float(arr['olr_w_m2'].mean())
    out['ocean_heat_w_m2'] = model.ocean_heat
    out['ground_drift_k'] = model.day_means[-1] - model.day_means[-2]
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--only', nargs='*')
    parser.add_argument('--out', type=Path, default=RESULTS)
    parser.add_argument('--processes', type=int, default=0)
    args = parser.parse_args(argv)
    if args.processes:
        eq.SW_CONFIG = replace(optics.SHORTWAVE, processes=args.processes)
    plan = cases()
    (args.out / 'series').mkdir(parents=True, exist_ok=True)
    (args.out / 'cases').mkdir(parents=True, exist_ok=True)
    columns = {}
    for name, case in plan.items():
        if (args.only and name not in args.only) or (args.out / 'cases' / f'{name}.json').is_file():
            continue
        key = case.equilibrium
        if key not in columns:
            columns[key] = Column(case)
        column = columns[key]
        column.case = case
        model = Model(case, column)
        checkpoint = args.out / 'checkpoints' / f'{name}.npz'
        series = model.run(checkpoint=checkpoint)
        summary = summarise(case, series, column, model)
        with open(args.out / 'series' / f'{name}.csv', 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(series[0]))
            writer.writeheader()
            writer.writerows({k: (f'{v:.6g}' if isinstance(v, float) else v) for k, v in r.items()} for r in series)
        (args.out / 'cases' / f'{name}.json').write_text(json.dumps(summary, indent=1) + '\n')
        checkpoint.unlink(missing_ok=True)
        print(f"{name:44s} ground {summary['ground_min_k']:6.1f}-{summary['ground_max_k']:6.1f} K  air "
              f"{summary['air_min_k']:6.1f}-{summary['air_max_k']:6.1f} K  column {summary['column_mean_min_k']:6.1f}-"
              f"{summary['column_mean_max_k']:6.1f} K  mixed layer {summary['max_mixed_layer_m']:6.0f} m", flush=True)
    rows = [json.loads((args.out / 'cases' / f'{n}.json').read_text()) for n in plan
            if (args.out / 'cases' / f'{n}.json').is_file()]
    if rows:
        with open(args.out / 'lunar_day.csv', 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        code = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()[:16] for p in (Path(__file__),)}
        meta = dict(schema='terluna.climate.lunar-day/1', producer=dict(domain='climate', files=code,
                    column_source='atmosphere/middle_atmosphere/results'),
                    evidence=' '.join(part.replace('\n', ' ') for part in __doc__.split('\n\n')[1:]),
                    reading_rule=('Temperatures of the last simulated solar day; series hours count from local '
                                  'noon. ground is the land skin or sea mixed layer; air is the lowest layer (about '
                                  '5 m); column is the mass-weighted mean of the whole air column. The sea is held '
                                  'to the reference surface temperature on average: ocean_heat_w_m2 is the heat '
                                  'currents must bring (negative: carry away) to hold it (zero for land). '
                                  'ground_drift_k is the change in the day-mean ground temperature between the last '
                                  'two days. Compare cases, not absolute values, where the relaxation, wind and '
                                  'mixing assumptions differ.'),
                    table='lunar_day.csv', series='series/<case>.csv')
        (args.out / 'lunar_day.json').write_text(json.dumps(meta, indent=1) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
