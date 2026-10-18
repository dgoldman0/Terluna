"""The waters' oxygen in the three target climates: what the seas and lakes hold at saturation, and a lake model
through the lunar month for stratification, mixing and dissolved oxygen in deep crater lakes and productive shallows.

    python -m research.studies.lake_oxygen.lakes                                  # the lakes (a few seconds)
    climate/gcm/.venv/bin/python -m research.studies.lake_oxygen.extract           # the GCM's forcing (a minute)
    python -m research.studies.lake_oxygen.run                                    # the runs (about an hour)
    python -m research.studies.lake_oxygen.run --quick                            # two years, to research/runs

It reads results/lakes.json, research/runs/lake_oxygen/forcing.npz with forcing.json, and the committed products of
climate_targets (the targets, their seas' temperatures and the share of the way between the two runs), sea_chemistry
(the held fresh seas' salinity) and geography's groundwater uptake (the heat flow through the floor), and writes
results/lake_oxygen.json. Each lake run is cached in research/runs/lake_oxygen/cases under a hash of the code, the
inputs and the case, so a run stopped part way resumes where it stopped.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
import time
from pathlib import Path

import numpy as np

from atmosphere.radiative_convective.thermodynamics import earthlike_air
from biosphere.ecology import lake as L
from shared.constants import GAS_CONSTANT, STEFAN_BOLTZMANN, SYNODIC_MONTH_DAYS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUNS = ROOT / 'research' / 'runs' / 'lake_oxygen'
FORCING = RUNS / 'forcing.npz'
FORCING_META = RUNS / 'forcing.json'
LAKES = HERE / 'results' / 'lakes.json'
TARGETS = ROOT / 'research' / 'studies' / 'climate_targets' / 'results' / 'climate_targets.json'
SEA_CHEMISTRY = ROOT / 'research' / 'studies' / 'sea_chemistry' / 'results' / 'sea_chemistry.json'
UPTAKE = ROOT / 'geography' / 'results' / 'groundwater_uptake.json'
OUT = HERE / 'results' / 'lake_oxygen.json'
SCHEMA = 'terluna.research.lake-oxygen/2'

DESIGN_PRESSURE_PA = 121590.0          # decisions register: 1.2 atm at the water level
DESIGN_CO2_PPM = 400.0
HELD_FRESH = dict(corner='central', delivery='fresh', sinks='managed', sediment='kept_open')
TARGET_NAMES = dict(hot='no further dimming', resort='300 K', lock='295 K')
DEEP_M = 352.0                         # the rain-fed lakes' median local depth (results/lakes.json)
AREA_MEDIAN_DEEPEST_M = 1925.0         # half the lake area lies in lakes whose deepest point passes this
SHALLOW_M, SHALLOW_LAYER_M = 3.0, 0.25
STEPS_PER_MEAN = 72                    # model steps in each 3-day mean, of 3575 s
LOWEST_FLOOR_MS = 1.0                  # PlaSim's minimum wind for the surface exchange (fluxmod.f90 zumin)
SOLAR_HOUR_BINS = 10
CENTRAL = '/renewal/supply'           # the columns the variant runs keep: surface renewal, demand fixed by supply
SITE_BANDS = ('0-15', '15-30', '30-45')  # the chosen lakes' bands, in results/lakes.json order
HYPOXIC, ANOXIC = 2.0, 0.5             # g/m3
G_C_TO_O2 = 32.0 / 12.0                # photosynthetic quotient of one
COARSE = 4                             # the spin-up's steps join this many model steps
WARM_START_K = 1.0                     # the column starts this much above the air's mean and cools from above
SETTLED_K, SETTLED_O2 = 0.005, 0.02    # a spin-up pass settles when the deep water's mean moves less than these
MAX_PASSES = dict(deep=20, deepest=0, deepest_long=30, shallow=1, perturbed=20)
# After the coarse passes settle, the deep lakes pass over the record again at the analysis step until the settle test
# passes there too: the 8-hour steps settle a slightly different deep water, and the review of 18 October found the
# analysed decade drifting from it (its bottom water 4-53% older than the coarse spin-up's last pass)
FINE_PASSES = dict(deep=10, deepest=0, deepest_long=10, shallow=0, perturbed=10)
AFTER_OVERTURN_DAYS = 30.0             # the window in which an overturn's oxygen at the surface is read

# Oxygen demand from oligotrophic to eutrophic in g O2 per m2 per day, fixed by each state's supply of organic
# matter (sources.json):
# gross production from Wetzel's (2001) ranges (oligotrophic 50-300, mesotrophic 250-1,000, eutrophic above
# 1,000 mg C per m2 a day) at 200, 600 and 2,000 mg C; the light's attenuation from Carlson's (1977) Secchi depths
# (oligotrophic beyond 4 m, mesotrophic 2-4 m, eutrophic 0.5-2 m) at 6, 3 and 1 m, with Kd = 1.7/SD (Poole & Atkins
# 1929); the sediment and the water below the euphotic zone from deep European lakes (Steinsberger et al. 2020:
# sediment uptake 7.9-10.6 mmol, water-column mineralisation 2.1-27.4 mmol per m2 a day) and eutrophic ones
# (Muller et al. 2012: areal hypolimnetic mineralisation 0.47-1.31, mean 0.90, of which reduced substances from the
# sediment 0.37 +- 0.12 g O2 per m2 a day). New reservoirs run eutrophic for six to fourteen years (Paterson et al.
# 2019), here with Muller's highest mineralisation. The euphotic water respires four fifths of production where it is
# made, so that eutrophic production exports about what Muller's eutrophic lakes consume below.
TROPHIC = {
    'oligotrophic': dict(gpp=0.200 * G_C_TO_O2, secchi_m=6.0, sediment=0.25, deep_water=0.07),
    'mesotrophic': dict(gpp=0.600 * G_C_TO_O2, secchi_m=3.0, sediment=0.30, deep_water=0.30),
    'eutrophic': dict(gpp=2.000 * G_C_TO_O2, secchi_m=1.0, sediment=0.37, deep_water=0.53),
    'new_reservoir': dict(gpp=2.000 * G_C_TO_O2, secchi_m=1.0, sediment=0.49, deep_water=0.82),
}
TRANSFER = ('renewal', 'cole_caraco')
KINETICS = dict(supply=0.0, warm=0.62)  # Arrhenius energy (eV): ecosystem respiration's 0.62 (Yvon-Durocher 2012)
# The demand is fixed by the trophic state's supply of organic matter, whatever the water's warmth. The 'warm' columns
# let the deep water's and floor's rates follow respiration's warm kinetics from the temperature of the hypolimnia
# they were measured in, 4-8 C (Muller et al. 2012; Steinsberger et al. 2020): about 3-8 times the deep rates and
# 5-13 times the shallow floors' in these lakes, which no settling supply sustains. They show the weeks to months
# after a warming while a stock of settled matter lasts, a transient kept out of the bounds.
SOURCE_C = 6.0
# The euphotic water respires its share of the day's production from a labile pool that decays with an e-folding
# time: a plankton community's dark respiration falls through the night at 0.04 an hour towards a steady level
# (Markager et al. 1992), and over the 354-hour lunar night the pool's 2-7 days bracket that run-down (design values,
# 4 days central).
LABILE_DAYS = dict(central=4.0, fast=2.0, slow=7.0)
# Anaerobic mineralisation's products reach the density through the debt: its methane and CO2 at Lake Kivu's
# coefficients, with half its oxygen equivalents ending as methane, lake_anoxia's central share for sulphate-poor water
# (Crowe et al. 2011). The iron bracket adds iron reduction's dissolved iron and bicarbonate for the other half, up to
# the 120 uM of Fe(II) ferruginous Lake Matano holds (Crowe et al. 2008), at the model's 0.8 kg/m3 per g/kg: each gram
# of oxygen equivalents reduces 0.125 mol of iron to Fe(II) and leaves 0.25 mol of bicarbonate, 22.2 g of salts
# (CH2O + 4 Fe(OH)3 + 7 CO2 -> 4 Fe2+ + 8 HCO3- + 3 H2O), consuming 7/32 mol of CO2 and making none of the 1/32 mol
# the debt's other products count, so 0.25 mol of CO2 is netted off at Lake Kivu's coefficient (2 mol for each mol of
# iron at the cap). The cap acts on the debt, so the bracket overstates where capped and uncapped water mix.
METHANE_SHARE = 0.5
IRON = dict(per_debt=(1.0 - METHANE_SHARE) * ((0.125 * 55.845 + 0.25 * 61.017) * 0.8e-3 - 0.25 * 44.01 * 0.284e-3),
            cap=(120e-3 * 55.845 + 240e-3 * 61.017) * 0.8e-3 - 240e-3 * 44.01 * 0.284e-3)
# Below the mixed layer the wind's energy sets the dissipation, the treatment the author chose on 18 October: the
# share of the wind's energy flux at 10 m (P10 = rho_air C10 U10^3, with Wuest & Lorke's 2003 drag) that reaches the
# stratified water as turbulent kinetic energy, 0.3 +/- 0.1% across Earth's lakes (Wuest, Piepke & Van Senden 2000;
# 0.42% in Lake Alpnach), spread over the water below the mixed layer, with Osborn's rule for the diffusivity. The cap
# is raised to 0.1 m2/s, above what a step mixes across the thickest layers, so that it holds no interface back. The
# fixed dissipation of the first round (1e-10 and 1e-9 W/kg, capped at 1e-4 m2/s) stays as sensitivity rows, and
# Yeates & Imberger's (2003) diffusivity from the Lake number, calibrated in DYRESM, as a bracket on the 0-15 lakes.
WIND_SHARE = dict(low=0.002, central=0.003, high=0.004)
ETA_BRACKET = ('low', 'high')
DEEP_MIXING = 'wind'                   # the lake model's deep-mixing treatment (lake.DEEP_MIXING), central
CAP_WIND = 1e-1                        # m2/s: the cap under the wind and the Lake number, raised to hold none back
CAP_WIND_RAISED = 1.0                  # m2/s: ten times that, the check that the cap holds nothing back
DISSIPATION = dict(low=1e-11, central=1e-10, high=1e-9)   # W/kg: the fixed dissipation's sensitivity rows
CAP = dict(low=1e-5, central=1e-4, high=1e-3)              # m2/s: the fixed dissipation's cap and its brackets
UNSTABLE_DIFFUSIVITY = 1e-4            # m2/s across statically unstable interfaces, raised to the cap where higher
LAKE_NUMBER_BELOW_M = dict(day=2.0, month=10.0)   # z_T: the strongest stratification below this depth, the day's
# warm layer (2 m) or the lunar month's thermocline (10 m); the review of 18 October found the 2 m rule taking the
# day's warm layer on 64-65% of days
SETTLED_AGE, SETTLED_SOLUTE = 0.02, 0.02   # the deep water's age (share) and solute (g/m3) must also settle a pass
LAKE_NUMBER_EVERY_DAYS = 3                 # the Lake number from the column every three days
# What aerobic respiration's products add to the water's density, per g/m3 of O2 respired: the CO2 alone, 1.375 g a
# gram of O2 at 0.284e-3 per g/kg (Lake Kivu's coefficient), or the CO2 dissolving calcite into calcium bicarbonate,
# 5.07 g a gram of O2 at the model's salinity term of 0.8 kg/m3 per g/kg; 'iron' is the CO2 with the iron bracket.
SOLUTE = dict(none=0.0, co2=44.01 / 32.0 * 0.284e-3, carbonate=162.11 / 32.0 * 0.8e-3, iron=44.01 / 32.0 * 0.284e-3)
COUPLED = dict(deep='mesotrophic/renewal/supply', shallow='eutrophic/renewal/supply')   # whose products the
# physics carries: the deep lakes' mesotrophic water, the productive shallows' eutrophic water
STEPS_PER_MEAN_DEEP = 36               # the deep lakes' steps of 7150 s


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def kd(secchi_m):
    return 1.7 / secchi_m


def scenarios():
    """Each trophic state's oxygen columns: two transfer velocities by two kinetics with the central labile pool, and
    the pool's fast and slow brackets at the central transfer and kinetics."""
    out = []
    for trophic, p in TROPHIC.items():
        common = dict(gpp=p['gpp'], attenuation=kd(p['secchi_m']), sediment=p['sediment'], deep_water=p['deep_water'])
        for transfer in TRANSFER:
            for kinetics, ev in KINETICS.items():
                out.append(L.Scenario(f'{trophic}/{transfer}/{kinetics}', transfer=transfer, activation_ev=ev,
                                      reference_c=SOURCE_C if ev else 20.0, labile_days=LABILE_DAYS['central'],
                                      **common))
        for label in ('fast', 'slow'):
            out.append(L.Scenario(f'{trophic}/renewal/supply/labile_{LABILE_DAYS[label]:g}d', transfer='renewal',
                                  labile_days=LABILE_DAYS[label], **common))
    return out


def dumps(obj) -> str:
    """JSON with one key a line and short lists of numbers on one line."""
    text = json.dumps(obj, indent=1)
    return re.sub(r'\[\s+([-0-9.eE,\s]+?|null|true|false)\s+\]', lambda m: '[' + ' '.join(m.group(1).split()) + ']',
                  text)


def weighted_quantile(x, w, q):
    i = np.argsort(x, kind='stable')
    c = np.cumsum(w[i])
    return float(np.interp(q * c[-1], c - 0.5 * w[i], x[i]))


# --- the setting ----------------------------------------------------------------------------------------------------

class Setting:
    def __init__(self, quick=False):
        self.quick = quick
        self.lakes = json.loads(LAKES.read_text())
        self.meta = json.loads(FORCING_META.read_text())
        if self.meta['npz_sha256'] != sha256(FORCING):
            raise SystemExit('forcing.npz does not match forcing.json; run the extract again')
        if self.meta['producer']['read']['research/studies/lake_oxygen/results/lakes.json'] != sha256(LAKES):
            raise SystemExit('the forcing was extracted for another lakes.json; run the extract again')
        self.f = dict(np.load(FORCING))
        self.targets = json.loads(TARGETS.read_text())
        if self.targets['schema'] != 'terluna.research.climate-targets/2':
            raise SystemExit(f"climate_targets is {self.targets['schema']}")
        chem = json.loads(SEA_CHEMISTRY.read_text())
        held = next(r for r in chem['runs'] if all(r[k] == v for k, v in HELD_FRESH.items()))
        self.salinity = float(held['rows'][-1]['salinity_g_kg'])
        self.held_ph = float(held['rows'][-1]['ph'])
        uptake = json.loads(UPTAKE.read_text())
        self.heat_flow = float(uptake['cases']['central']['params']['heat_flow'])
        dry = earthlike_air(DESIGN_PRESSURE_PA, DESIGN_CO2_PPM)
        self.o2_fraction = float(dry.fractions['O2'])
        self.air = L.Air(gas_constant=GAS_CONSTANT / dry.molar_mass, cp=float(dry.cp))
        self.gravity = float(self.meta['runs']['lock']['gravity'])
        self.window_s = self.meta['runs']['lock']['steps_per_mean'] * self.meta['runs']['lock']['timestep_min'] * 60.0
        self.dt = self.window_s / STEPS_PER_MEAN
        self.share = {k: next(t for t in self.targets['targets'] if t['name'] == v)['interpolated_share_of_the_way']
                      for k, v in TARGET_NAMES.items()}
        self.lapse = {k: self._lapse(k) for k in TARGET_NAMES}
        self.phase = {run: self._phase(run) for run in ('lock', 'hot', 'eddies')}
        self.eddy_factor = self._eddy_factor()
        self.inputs = {str(p.relative_to(ROOT)): sha256(p) for p in (LAKES, FORCING, FORCING_META, TARGETS,
                                                                       SEA_CHEMISTRY, UPTAKE)}

    def _lapse(self, name):
        """K per m by which the air warms going down, from the target's own column over land (0.9 to 5 km)."""
        t = next(x for x in self.targets['targets'] if x['name'] == TARGET_NAMES[name])['climate']['ladder']
        by = {r['height_km']: r['air_c'] for r in t}
        return (by[0.9] - by[5]) / 4100.0

    def _phase(self, run):
        """The subsolar longitude at time zero (degrees east), on the synodic clock, from each mean's Sun."""
        r = self.meta['runs'][run]
        w = r['steps_per_mean'] * r['timestep_min'] * 60.0
        sub = self.f[f'{run}__subsolar_lon']
        t = (np.arange(sub.size) + 0.5) * w
        omega = 360.0 / (SYNODIC_MONTH_DAYS * 86400.0)
        est = np.radians(sub + omega * t)
        return float(np.degrees(np.arctan2(np.sin(est).mean(), np.cos(est).mean())) % 360.0)

    def _eddy_factor(self):
        """By how much the mean speed of the 6.4-hour means exceeds that of their 3-day mean vector, in its cube
        (the wind's stirring), at each site, from the two years of the design case written every 6.4 hours."""
        u, v = self.f['eddies__ua'], self.f['eddies__va']
        r = self.meta['runs']['eddies']
        per = r['writes_per_lunar_day'] // self.meta['runs']['lock']['writes_per_lunar_day']
        n = (u.shape[0] // per) * per
        uu, vv = u[:n].reshape(-1, per, u.shape[1]), v[:n].reshape(-1, per, v.shape[1])
        speed = np.hypot(uu, vv)
        vector = np.hypot(uu.mean(1), vv.mean(1))
        return ((speed ** 3).mean((0, 1)) / (vector ** 3).mean(0)) ** (1.0 / 3.0)

    def mu(self, lat, lon, t, run):
        """Cosine of the Sun's zenith angle (clipped at the horizon) at a place and time on the run's clock; the
        Moon's 1.5 degree tilt left out."""
        omega = 360.0 / (SYNODIC_MONTH_DAYS * 86400.0)
        hour = np.radians(lon - (self.phase[run] - omega * t))
        return np.maximum(0.0, math.cos(math.radians(lat)) * np.cos(hour))


# --- forcing --------------------------------------------------------------------------------------------------------

def site_info(st: Setting, i, water=False):
    s = st.meta['sites'][i]
    lake = st.lakes['chosen'][i]
    col = i + 3 if water else i
    return dict(index=col, band=s['band'], lat=s['lat'], lon=s['water_lon'] if water else s['lon'],
                elevation=s['water_elevation_m'] if water else s['elevation_m'], lake=lake,
                eddy_factor=float(st.eddy_factor[col]))


def window_humidity(hus):
    """Each window's humidity from PlaSim's records. The output's specific humidity is the value at the step that
    opens the window, held through it (plasim.f90 computes sqout only when the step is a multiple of the output
    interval, and outmod.f90 accumulates that held value), where the other fields are means over the window; so each
    window takes the mean of its opening record and the next, and the last window its own."""
    hus = np.asarray(hus, dtype=float)
    out = hus.copy()
    out[:-1] = 0.5 * (hus[:-1] + hus[1:])
    return out


def windows(st: Setting, run, site, at_lake=True, target=None, eddy=True):
    """The GCM's means at a cell as forcing: sunlight and sky infrared at the surface, the lowest level's potential
    temperature at the surface pressure (PlaSim's tas, its dtsa: the level, about 0.86 km up, brought down its dry
    adiabat), that level's humidity and wind and its height, which the bulk exchange pairs as PlaSim does, the wind at
    10 m for Cole & Caraco's transfer, and the pressure, carried to the lake's level when at_lake.

    The carry-down moves the potential temperature by the target's free-air lapse rate over the lake's depth below the
    cell. The air the lake's exchange takes is the free air at the lowest level's height above the lake; that air is
    warmer than the cell's lowest level by the lapse rate times the depth, and its potential temperature at the lake's
    own surface pressure is its temperature times sigma^-kappa, as tas is the cell's. Carrying tas by the lapse rate is
    that to within the factor sigma^-kappa = 1.005 on the carry, a few thousandths of a kelvin."""
    f, c = st.f, site['index']
    g = st.gravity
    rd = st.air.gas_constant
    tas, ta = f[f'{run}__tas'][:, c], f[f'{run}__ta'][:, c]
    q = window_humidity(f[f'{run}__hus'][:, c])
    ps = f[f'{run}__ps'][:, c] * 100.0
    alb, rss, rls, ts, lsm = (f[f'{run}__{k}'][:, c] for k in ('alb', 'rss', 'rls', 'ts', 'lsm'))
    sw = np.where(alb < 0.999, rss / np.maximum(1.0 - alb, 1e-3), 0.0)
    eps = np.where(lsm > 0.5, 1.0, L.WATER_EMISSIVITY)
    lw = rls / eps + STEFAN_BOLTZMANN * ts ** 4
    sigma = st.meta['runs'][run]['lowest_sigma']
    z_low = rd * 0.5 * (ta + tas) / g * math.log(1.0 / sigma)
    speed = np.hypot(f[f'{run}__ua'][:, c], f[f'{run}__va'][:, c]) * (site['eddy_factor'] if eddy else 1.0)
    u_low = np.maximum(speed, LOWEST_FLOOR_MS)                       # PlaSim's zumin on its lowest level's wind
    u10 = np.array([L.reduce_wind(u, z, 10.0, g) for u, z in zip(u_low, z_low)])
    air = tas.copy()
    p = ps.copy()
    if at_lake:
        dz = site['elevation'] - site['lake']['level_above_sea_m']
        air = tas + st.lapse[target] * dz
        lw = lw * (air / tas) ** 4
        p = ps * np.exp(g * dz / (rd * 0.5 * (air + tas)))
    qs = np.array([L.saturation_humidity(a, pp) for a, pp in zip(air, p)])
    r = st.meta['runs'][run]
    w = r['steps_per_mean'] * r['timestep_min'] * 60.0
    t = (np.arange(tas.size) + 0.5) * w
    omega = 360.0 / (SYNODIC_MONTH_DAYS * 86400.0)
    hour = ((site['lon'] - (st.phase[run] - omega * t) + 180.0) % 360.0) - 180.0
    return dict(run=run, window_s=w, sw=sw, lw=lw, air=air, q=np.minimum(q, qs), u10=u10, u_low=u_low, z_low=z_low,
                p=p, hour=hour, gcm_ts=ts, gcm_alb=alb,
                gcm_turbulent=f[f'{run}__hfss'][:, c] + f[f'{run}__hfls'][:, c])


def resort_windows(st: Setting, lock_w, hot_w, from_side='lock'):
    """The resort between the runs by the ground's temperature: one run's means moved by the share of the way of the
    difference between the two runs' mean cycles, binned by the solar hour. It mixes the two runs' physics."""
    share = st.share['resort']
    base, factor = (lock_w, share) if from_side == 'lock' else (hot_w, share - 1.0)
    bins = lambda h: np.clip(((h + 180.0) / 360.0 * SOLAR_HOUR_BINS).astype(int), 0, SOLAR_HOUR_BINS - 1)

    def composite(w, k):
        b = bins(w['hour'])
        return np.array([w[k][b == i].mean() for i in range(SOLAR_HOUR_BINS)])

    out = dict(base)
    b = bins(base['hour'])
    for k in ('sw', 'lw', 'air', 'q', 'u10', 'u_low', 'z_low', 'p'):
        out[k] = base[k] + factor * (composite(hot_w, k) - composite(lock_w, k))[b]
    out['sw'] = np.maximum(out['sw'], 0.0)
    out['u10'] = np.maximum(out['u10'], 0.0)
    out['u_low'] = np.maximum(out['u_low'], LOWEST_FLOOR_MS)
    qs = np.array([L.saturation_humidity(a, pp) for a, pp in zip(out['air'], out['p'])])
    out['q'] = np.clip(out['q'], 0.0, qs)
    return out


def hourly(st: Setting, w, site, first=0, count=None, steps_per_mean=STEPS_PER_MEAN, paired=True):
    """Step forcing from means: linear in time between the means' centres, the sunlight shaped within each mean by
    the Sun's height there and scaled to the mean. paired gives the bulk exchange the lowest level's wind and height,
    as PlaSim pairs them; without it the exchange takes the 10 m wind and coefficient, the pairing before 18 October."""
    n_w = w['air'].size
    count = n_w - first if count is None else count
    sl = slice(first, first + count)
    window = w['window_s']
    dt = window / steps_per_mean
    centres = (np.arange(n_w) + 0.5) * window
    t = first * window + (np.arange(count * steps_per_mean) + 0.5) * dt
    keys = ('lw', 'air', 'q', 'u10', 'p') + (('u_low', 'z_low') if paired and 'u_low' in w else ())
    interp = {k: np.interp(t, centres, w[k]) for k in keys}
    mu = st.mu(site['lat'], site['lon'], t, w['run']).reshape(count, steps_per_mean)
    mean = mu.mean(1, keepdims=True)
    shape = np.where(mean > 1e-6, mu / np.maximum(mean, 1e-12), 1.0)
    sw = (w['sw'][sl][:, None] * shape).ravel()
    return L.Forcing(dt=dt, sw_down=sw, lw_down=interp['lw'], air_k=interp['air'], q_air=interp['q'],
                     wind10=interp['u10'], pressure=interp['p'], o2_fraction=st.o2_fraction,
                     wind_level=interp.get('u_low'), level_height=interp.get('z_low')), t, mu.ravel()


# --- the lakes ------------------------------------------------------------------------------------------------------

def basin(kind, site, depth=None):
    lake = site['lake']
    if kind == 'shallow':
        return L.Basin.pan(SHALLOW_M, SHALLOW_LAYER_M, 1e6)
    h = lake['hypsography']
    depth = depth or DEEP_M
    return L.Basin.from_hypsography(depth, lake['area_km2'] * 1e6, h['depth_fraction'], h['area_fraction'],
                                    max_m=8.0 if depth <= 500.0 else 25.0)


def physics(st: Setting, kind, mixing=DEEP_MIXING, eta='central', dissipation='central', full_mixing=False,
            albedo=None, solute='co2', coupled=None, cap=None):
    """The column's physics: the coupled trophic state's water clarity and its respiration's products in the density,
    aerobic CO2 and the debt's anaerobic CO2 and methane ('none' carries neither; 'carbonate' takes the aerobic CO2 as
    calcite's bicarbonate; 'iron' adds iron reduction's capped dissolved load); the deep mixing, the wind's share eta
    of P10, the fixed dissipation or the Lake number's, and the cap (CAP_WIND under the wind and the Lake number,
    CAP['central'] under the fixed dissipation, unless given by CAP's name or in m2/s)."""
    coupled = coupled or COUPLED[kind]
    trophic = coupled.split('/')[0]
    names = [s.name for s in scenarios()]
    if cap is None:
        cap_m2_s = CAP['central'] if mixing == 'dissipation' else CAP_WIND
    else:
        cap_m2_s = CAP[cap] if isinstance(cap, str) else float(cap)
    return L.Physics(gravity=st.gravity, salinity=st.salinity, attenuation=kd(TROPHIC[trophic]['secchi_m']),
                     deep_mixing=mixing, wind_share=WIND_SHARE[eta], dissipation=DISSIPATION[dissipation],
                     max_diffusivity=cap_m2_s, unstable_diffusivity=UNSTABLE_DIFFUSIVITY, heat_flow=st.heat_flow,
                     full_mixing=full_mixing, albedo=0.069 if albedo is None else albedo,
                     solute_density=SOLUTE[solute], solute_column=names.index(coupled),
                     debt_density=0.0 if solute == 'none' else L.anaerobic_density(METHANE_SHARE),
                     iron_density_per_debt=IRON['per_debt'] if solute == 'iron' else 0.0,
                     iron_density_cap=IRON['cap'] if solute == 'iron' else 0.0)


def nights(mu):
    """Index ranges of each night: from a sunset to the next sunrise."""
    up = mu > 0.0
    sunset = np.flatnonzero(up[:-1] & ~up[1:]) + 1
    sunrise = np.flatnonzero(~up[:-1] & up[1:]) + 1
    out = []
    for s in sunset:
        r = sunrise[sunrise > s]
        if r.size:
            out.append((int(s), int(r[0])))
    return out


def days(mu):
    """Index ranges of each day: from a sunrise to the next sunset."""
    up = mu > 0.0
    sunrise = np.flatnonzero(~up[:-1] & up[1:]) + 1
    sunset = np.flatnonzero(up[:-1] & ~up[1:]) + 1
    out = []
    for r in sunrise:
        s = sunset[sunset > r]
        if s.size:
            out.append((int(r), int(s[0])))
    return out


def cycles(mu):
    """Index ranges of each lunar cycle at the site: from one sunrise to the next."""
    up = mu > 0.0
    sunrise = np.flatnonzero(~up[:-1] & up[1:]) + 1
    return [(int(a), int(b)) for a, b in zip(sunrise[:-1], sunrise[1:])]


def quantiles(x, qs=(0.1, 0.5, 0.9)):
    """The 10th, 50th and 90th percentiles as a list, or None for no values."""
    x = np.asarray(x, dtype=float)
    return [round(float(np.quantile(x, q)), 3) for q in qs] if x.size else None


def summarise(kind, out, mu, b: L.Basin, scen, dt, keep=None):
    hours = dt / 3600.0
    res = dict(surface_t_c=dict(mean=round(float(out['surface_t'].mean()), 2),
                                p1=round(float(np.quantile(out['surface_t'], 0.01)), 2),
                                p99=round(float(np.quantile(out['surface_t'], 0.99)), 2)),
               bottom_t_c=dict(mean=round(float(out['bottom_t'].mean()), 3),
                               min=round(float(out['bottom_t'].min()), 3), max=round(float(out['bottom_t'].max()), 3)),
               ice_max_m=round(float(out['ice'].max()), 3),
               saturation_mg_l=dict(mean=round(float(out['saturation'].mean()), 3),
                                    min=round(float(out['saturation'].min()), 3),
                                    max=round(float(out['saturation'].max()), 3)))
    cyc = cycles(mu)
    reach = np.array([out['mixed_depth'][a:b].max() for a, b in cyc])
    bottom_top = float(b.z[-2])
    full = reach >= bottom_top - 1e-6
    res['cycles'] = len(cyc)
    res['mixing'] = dict(complete_overturn_share=round(float(full.mean()), 4) if len(cyc) else None,
                         deepest_share_of_depth=quantiles(reach / b.depth),
                         deepest_m=quantiles(reach),
                         bottom_age_days=dict(mean=round(float(out['bottom_age'].mean()), 2),
                                              max=round(float(out['bottom_age'].max()), 2)))
    # complete overturns: the first step of each cycle at which the mixed layer reaches the bottom layer
    events = []
    for (a, b_), hit in zip(cyc, full):
        if hit:
            events.append(a + int(np.argmax(out['mixed_depth'][a:b_] >= bottom_top - 1e-6)))
    steps = out['mixed_depth'].size
    res['mixing']['complete_overturns'] = len(events)
    res['mixing']['years_between'] = (round(float(np.diff(events).mean() * dt / (365.25 * 86400.0)), 2)
                                      if len(events) > 1 else None)
    after = int(AFTER_OVERTURN_DAYS * 86400.0 / dt)
    ns = nights(mu)
    per = {}
    for m, s in enumerate(scen):
        if keep and not s.name.endswith(keep):
            continue
        so, bo, do = out['surface_o2'][:, m], out['bottom_o2'][:, m], out['deep_o2'][:, m]
        entry = dict(surface=dict(min=round(float(so.min()), 3), mean=round(float(so.mean()), 3),
                                  max=round(float(so.max()), 3)),
                     bottom=dict(min=round(float(bo.min()), 3), mean=round(float(bo.mean()), 3),
                                 hypoxic_share=round(float((bo < HYPOXIC).mean()), 4),
                                 anoxic_share=round(float((bo < ANOXIC).mean()), 4)),
                     deepest_tenth=dict(min=round(float(do.min()), 3), mean=round(float(do.mean()), 3),
                                        p5=round(float(np.quantile(do, 0.05)), 3),
                                        hypoxic_share=round(float((do < HYPOXIC).mean()), 4),
                                        anoxic_share=round(float((do < ANOXIC).mean()), 4)),
                     whole_lake_mean=round(float(out['mean_o2'][:, m].mean()), 3),
                     bottom_debt_max=round(float(out['bottom_debt'][:, m].max()), 3))
        if events:
            ends = events[1:] + [steps]
            hyp, anox = [], []
            for e, nxt in zip(events, ends):
                seg = do[e:nxt]
                h_i, a_i = np.flatnonzero(seg < HYPOXIC), np.flatnonzero(seg < ANOXIC)
                hyp.append(float(h_i[0] * dt / 86400.0) if h_i.size else None)
                anox.append(float(a_i[0] * dt / 86400.0) if a_i.size else None)
            reached = lambda xs: dict(share=round(sum(x is not None for x in xs) / len(xs), 3),
                                      days=quantiles([x for x in xs if x is not None]))
            recover = []
            for e, nxt in zip(events, ends):
                low = e + int(np.argmin(so[e:min(e + after, nxt)])) if nxt > e else e
                back = np.flatnonzero(so[low:nxt] >= HYPOXIC)
                recover.append(float(back[0] * dt / 86400.0) if back.size else None)
            entry['overturns'] = dict(
                lake_mean_after=quantiles([out['mean_o2'][min(e + 1, steps - 1), m] for e in events]),
                surface_minimum_after=quantiles([so[e:e + after].min() for e in events]),
                surface_back_above_hypoxic=reached(recover),
                to_hypoxia=reached(hyp), to_anoxia=reached(anox))
        h_i, a_i = np.flatnonzero(do < HYPOXIC), np.flatnonzero(do < ANOXIC)
        entry['deepest_tenth']['days_to_hypoxia_from_start'] = round(float(h_i[0] * dt / 86400.0), 1) \
            if h_i.size else None
        entry['deepest_tenth']['days_to_anoxia_from_start'] = round(float(a_i[0] * dt / 86400.0), 1) \
            if a_i.size else None
        if ns:
            low_s = np.array([so[a:b].min() for a, b in ns])
            low_b = np.array([bo[a:b].min() for a, b in ns])
            sat_at = np.array([out['saturation'][a + int(np.argmin(so[a:b]))] for a, b in ns])
            dawn_s = np.array([so[b - 1] for a, b in ns])
            when = np.array([np.argmin(so[a:b]) / max(b - a - 1, 1) for a, b in ns])
            hyp_h = np.array([(bo[a:b] < HYPOXIC).sum() * hours for a, b in ns])
            anox_h = np.array([(bo[a:b] < ANOXIC).sum() * hours for a, b in ns])
            entry['nights'] = dict(count=len(ns), surface_minimum=quantiles(low_s), bottom_minimum=quantiles(low_b),
                                   surface_minimum_share_of_saturation=quantiles(low_s / sat_at),
                                   surface_before_dawn=quantiles(dawn_s),
                                   surface_minimum_share_of_the_night=quantiles(when),
                                   hypoxic_share=round(float((low_s < HYPOXIC).mean()), 3),
                                   anoxic_share=round(float((low_s < ANOXIC).mean()), 3),
                                   hours_bottom_hypoxic=quantiles(hyp_h), hours_bottom_anoxic=quantiles(anox_h))
            ds = days(mu)
            if ds:
                entry['days'] = dict(count=len(ds),
                                     bottom_minimum=quantiles([bo[a:b].min() for a, b in ds]),
                                     hours_bottom_anoxic=quantiles([(bo[a:b] < ANOXIC).sum() * hours for a, b in ds]),
                                     surface_maximum_share_of_saturation=quantiles(
                                         [so[a:b].max() / out['saturation'][a + int(np.argmax(so[a:b]))]
                                          for a, b in ds]))
        per[s.name] = entry
    res['oxygen'] = per
    return res


def coarse(forcing: L.Forcing, k=COARSE):
    """The forcing in means of k steps, for the spin-up."""
    n = len(forcing) // k * k
    agg = lambda x: np.asarray(x, dtype=float)[:n].reshape(-1, k).mean(1)
    return L.Forcing(dt=forcing.dt * k, sw_down=agg(forcing.sw_down), lw_down=agg(forcing.lw_down),
                     air_k=agg(forcing.air_k), q_air=agg(forcing.q_air), wind10=agg(forcing.wind10),
                     pressure=agg(forcing.pressure), o2_fraction=forcing.o2_fraction,
                     wind_level=None if forcing.wind_level is None else agg(forcing.wind_level),
                     level_height=None if forcing.level_height is None else agg(forcing.level_height))


def state_from(saved, bsn: L.Basin, scen):
    """A saved column carried to another basin of the same lake: each new layer takes the saved layer at its centre's
    depth, and below the saved column's floor the floor's water."""
    z = np.asarray(saved['centre'])
    idx = np.clip(np.searchsorted(z, bsn.centre), 0, z.size - 1)
    pool = np.asarray(saved['pool'])[idx].copy() if 'pool' in saved else None
    state = L.State(np.asarray(saved['t'])[idx].copy(), np.asarray(saved['oxygen'])[idx].copy(),
                    np.asarray(saved['debt'])[idx].copy(), np.asarray(saved['age'])[idx].copy(),
                    solute=np.asarray(saved['solute'])[idx].copy(), pool=pool)
    state.mixed_depth = float(saved['mixed_depth'])
    return state


def forcing_for(st: Setting, target, site, resort_from='lock', years=None, steps_per_mean=STEPS_PER_MEAN):
    if target == 'resort':
        w = resort_windows(st, windows(st, 'lock', site, target='resort'), windows(st, 'hot', site, target='resort'),
                           from_side=resort_from)
    else:
        w = windows(st, target, site, target=target)
    n_w = w['air'].size if years is None else years * 120
    forcing, t, mu = hourly(st, w, site, 0, n_w, steps_per_mean)
    return w, forcing, mu


def lake_number(bsn: L.Basin, temps, u2, gravity, below_m=LAKE_NUMBER_BELOW_M['month']):
    """The Lake number (Imberger & Patterson 1990, in the form the review's check takes it from rLakeAnalyzer) of each
    saved column: L_N = g S z_T / (rho u*^2 A0^(1/2) z_cv), with S the column's Schmidt work per unit g, z_T the depth
    of the strongest stratification below below_m (2 m takes the day's warm layer, 10 m the lunar month's
    thermocline) and z_cv the centre of volume. Under 1 the wind's stress tilts the thermocline to the surface and
    stirs the deep water, basin-scale motion a one-dimensional column does not hold."""
    v = bsn.volume
    zcv = float((bsn.centre * v).sum() / v.sum())
    out = []
    below = bsn.centre[1:] > below_m
    for t, uu in zip(temps, u2):
        rho = L.density(t)
        rbar = (rho * v).sum() / v.sum()
        work = ((bsn.centre - zcv) * (rho - rbar) * v).sum() / bsn.area[0]
        n2 = np.diff(rho) / np.diff(bsn.centre)
        zt = float(bsn.centre[1:][below][np.argmax(n2[below])]) if below.any() else float(bsn.centre[-1])
        out.append(gravity * work * zt / (L.RHO0 * max(uu, 1e-12) * math.sqrt(bsn.area[0]) * zcv))
    return np.array(out)


def settled(history):
    """A pass settles when the deep water's mean temperature, the coupled column's deepest-tenth oxygen, the bottom
    water's age and its solute each move less than their tolerance from the pass before; a column that has stopped
    exchanging keeps ageing and so does not pass."""
    if len(history) < 2:
        return False
    (t1, o1, a1, s1), (t2, o2, a2, s2) = history[-2], history[-1]
    return (abs(t2 - t1) < SETTLED_K and abs(o2 - o1) < SETTLED_O2
            and abs(a2 - a1) <= SETTLED_AGE * max(a1, a2, 30.0) and abs(s2 - s1) < SETTLED_SOLUTE)


def deep_mixing_record(ph: L.Physics, out):
    """What set the mixing below the mixed layer through the analysed record: the treatment, the wind's share and the
    dissipation it gave (W/kg), and the share of the interfaces below the mixed layer at the cap, held back by it
    (the cap slower than a step's own mixing across the interface) and statically unstable."""
    d = np.asarray(out['dissipation'], dtype=float)
    d = d[np.isfinite(d)]
    below = max(int(out['interfaces_below']), 1)
    return dict(treatment=ph.deep_mixing, wind_share=ph.wind_share if ph.deep_mixing == 'wind' else None,
                dissipation_w_kg=dict(mean=float(d.mean()), p10=float(np.quantile(d, 0.1)),
                                      median=float(np.median(d)), p90=float(np.quantile(d, 0.9))) if d.size else None,
                cap_m2_s=ph.max_diffusivity, unstable_m2_s=max(ph.unstable_diffusivity, ph.max_diffusivity),
                share_at_the_cap=round(out['interfaces_capped'] / below, 5),
                share_held_by_the_cap=round(out['interfaces_held_by_cap'] / below, 5),
                share_unstable=round(out['interfaces_unstable'] / below, 5))


def budget_errors(bsn: L.Basin, state, h0, o0, out, seconds):
    """How closely a run keeps its heat and oxygen: what the column gained less what crossed its surface and floor
    (for oxygen, less production and plus respiration), against the column's own heat (above freezing, at least a
    kelvin's worth) and oxygen, and as a mean flux over the lake's surface through the run. A repeating record takes
    in almost no net heat, so the net is no scale for the error; what is left is the implicit diffusion's round-off,
    which the wind's large diffusivities across near-neutral water raise."""
    heat = L.heat_content(bsn, state) - h0 - out['heat_in']
    ox = L.oxygen_content(bsn, state) - o0 - out['oxygen_in'] - out['produced'] + out['respired']
    scale = max(abs(h0), L.RHO_CP * float(bsn.volume.sum()))
    area_s = float(bsn.area[0]) * seconds
    return dict(heat_relative_error=float(heat / scale), heat_error_w_m2=float(heat / area_s),
                oxygen_relative_error=float(np.max(np.abs(ox) / np.maximum(np.abs(o0), 1.0))),
                oxygen_error_g_m2_day=float(np.max(np.abs(ox)) / area_s * 86400.0))


def settle_passes(lake, forcing, state, col, bsn, limit):
    """Passes over the record with the given forcing until the settle test passes or the limit is reached; the
    deep water's mean temperature, the coupled column's deepest-tenth oxygen, the bottom water's age and solute by
    pass, and whether each pass overturned. The settle test needs two passes at the same step."""
    history, overturned = [], []
    for _ in range(limit):
        o = lake.run(forcing, state)
        history.append((float(o['bottom_t'].mean()), float(o['deep_o2'][:, col].mean()),
                        float(o['bottom_age'].mean()), float(o['bottom_solute'].mean())))
        overturned.append(bool((o['mixed_depth'] >= bsn.z[-2] - 1e-6).any()))
        if settled(history):
            break
    return dict(passes=len(history), settled=settled(history) if limit > 1 else None, limit=limit,
                steps_per_pass=len(forcing), step_s=forcing.dt,
                bottom_t_by_pass=[round(h[0], 4) for h in history], deep_o2_by_pass=[round(h[1], 3) for h in history],
                bottom_age_days_by_pass=[round(h[2], 1) for h in history],
                bottom_solute_by_pass=[round(h[3], 3) for h in history], overturned_by_pass=overturned)


def run_case(st: Setting, kind, target, site_i, dissipation='central', resort_from='lock', depth=None, years=None,
             start_from=None, solute='co2', coupled=None, keep=None, cap=None, long_spin=False, mixing=DEEP_MIXING,
             eta='central'):
    """One lake through a target's forcing: a spin-up over the whole record (in 8-hour steps for the deep lakes, a
    year for the shallows) from a warm, saturated column or a saved one, until it settles or reaches the passes'
    limit; for the deep lakes, passes over the record again at the analysis step until the settle test passes there
    as well (or their limit); then the record analysed step by step. A deeper lake taken from a saved column starts
    oxygenated, as if just overturned: with long_spin it passes over the record until it settles, else it runs the
    record from there."""
    site = site_info(st, site_i)
    spm = STEPS_PER_MEAN_DEEP if kind == 'deep' else STEPS_PER_MEAN
    w, forcing, mu = forcing_for(st, target, site, resort_from, years, spm)
    bsn = basin(kind, site, depth)
    scen = scenarios()
    lake = L.Lake(bsn, physics(st, kind, mixing, eta, dissipation, solute=solute, coupled=coupled, cap=cap), scen,
                  st.air)
    lake.calibrate_production(forcing)
    start = float(forcing.air_k.mean()) - 273.15 + WARM_START_K
    if start_from is not None:
        state = state_from(start_from, bsn, scen)
        if depth:                                       # the deeper lake starts oxygenated, as if just overturned
            state.oxygen[:] = float(L.oxygen_saturation(state.t.mean(), DESIGN_PRESSURE_PA / 101325.0, st.salinity,
                                                        st.o2_fraction))
            state.debt[:], state.age[:], state.solute[:], state.pool[:] = 0.0, 0.0, 0.0, 0.0
    else:
        state = L.initial_state(bsn, scen, start, float(L.oxygen_saturation(start, DESIGN_PRESSURE_PA / 101325.0,
                                                                             st.salinity, st.o2_fraction)))
    spin = coarse(forcing) if kind == 'deep' else hourly(st, w, site, 0, min(120, w['air'].size), spm)[0]
    tic = time.time()
    col = lake.physics.solute_column
    which = 'deepest_long' if long_spin else 'deepest' if depth else \
        ('perturbed' if start_from is not None and kind == 'deep' else kind)
    passes = MAX_PASSES[which]
    spun = settle_passes(lake, spin, state, col, bsn, passes)
    spun['start_c'] = None if start_from is not None else round(start, 2)
    # the deep lakes settle again at the analysis step, so the analysed record starts from that step's own column
    spun['at_the_analysis_step'] = settle_passes(lake, forcing, state, col, bsn, FINE_PASSES[which]) \
        if kind == 'deep' and FINE_PASSES[which] else None
    h0, o0 = L.heat_content(bsn, state), L.oxygen_content(bsn, state)
    every = int(round(LAKE_NUMBER_EVERY_DAYS * 86400.0 / forcing.dt)) if kind == 'deep' else 0
    out = lake.run(forcing, state, profiles_at=range(0, len(forcing), every) if every else ())
    errors = budget_errors(bsn, state, h0, o0, out, len(forcing) * forcing.dt)
    res = summarise(kind, out, mu, bsn, scen, forcing.dt, keep)
    res['bottom_solute_g_m3'] = dict(mean=round(float(out['bottom_solute'].mean()), 3),
                                     max=round(float(out['bottom_solute'].max()), 3))
    res['deep_mixing'] = deep_mixing_record(lake.physics, out)
    if every:
        day = int(round(86400.0 / forcing.dt))
        steps_at = sorted(out['profiles'])
        u2 = [float((out['u_star'][k:k + day] ** 2).mean()) for k in steps_at]
        temps = [out['profiles'][k]['t'] for k in steps_at]
        res['lake_number'] = {}
        for rule, below_m in LAKE_NUMBER_BELOW_M.items():          # the month's thermocline, and the day's warm layer
            ln = lake_number(bsn, temps, u2, st.gravity, below_m)
            res['lake_number'][rule] = dict(z_t_below_m=below_m, median=round(float(np.median(ln)), 3),
                                            p10=round(float(np.quantile(ln, 0.1)), 3),
                                            p90=round(float(np.quantile(ln, 0.9)), 3),
                                            share_under_1=round(float((ln < 1.0).mean()), 3), samples=int(ln.size))
    out.pop('profiles', None)
    year = 120 * spm
    first, last = slice(0, min(2 * year, len(forcing))), slice(max(0, len(forcing) - 2 * year), len(forcing))
    res['settling'] = dict(spin_up=spun,
                           bottom_t_first_two_years=round(float(out['bottom_t'][first].mean()), 3),
                           bottom_t_last_two_years=round(float(out['bottom_t'][last].mean()), 3),
                           deep_o2_first_two_years=[round(float(x), 3) for x in out['deep_o2'][first].mean(0)],
                           deep_o2_last_two_years=[round(float(x), 3) for x in out['deep_o2'][last].mean(0)])
    days = len(forcing) * forcing.dt / 86400.0
    res['budget'] = dict(errors,
                         production_g_m2_day=[round(float(x) / bsn.area[0] / days, 3) for x in out['produced']],
                         respiration_g_m2_day=[round(float(x) / bsn.area[0] / days, 3) for x in out['respired']],
                         from_the_air_g_m2_day=[round(float(x) / bsn.area[0] / days, 3) for x in out['oxygen_in']])
    res['case'] = dict(kind=kind, target=target, site=site['band'], lat=site['lat'], lon=site['lon'],
                       lake_level_above_sea_m=site['lake']['level_above_sea_m'], cell_elevation_m=site['elevation'],
                       depth_m=bsn.depth, layers=bsn.n, area_km2=round(bsn.area[0] / 1e6, 1),
                       deep_mixing=mixing, wind_share=WIND_SHARE[eta] if mixing == 'wind' else None,
                       dissipation=DISSIPATION[dissipation] if mixing == 'dissipation' else None,
                       cap=lake.physics.max_diffusivity,
                       resort_from=resort_from if target == 'resort' else None,
                       solute=solute, solute_density=SOLUTE[solute], coupled=coupled or COUPLED[kind],
                       debt_density=lake.physics.debt_density, pairing='lowest level' if forcing.wind_level is not None
                       else '10 m', labile_days=LABILE_DAYS, warm_reference_c=SOURCE_C,
                       step_s=forcing.dt,
                       years=round(days / 365.25, 2), seconds=round(time.time() - tic, 1),
                       forcing=dict(air_c=round(float(forcing.air_k.mean()) - 273.15, 2),
                                    sunlight_w_m2=round(float(forcing.sw_down.mean()), 1),
                                    sky_infrared_w_m2=round(float(forcing.lw_down.mean()), 1),
                                    wind10_m_s=round(float(forcing.wind10.mean()), 2),
                                    pressure_atm=round(float(forcing.pressure.mean()) / 101325.0, 4)))
    res['final_state'] = dict(centre=bsn.centre.tolist(), t=state.t.tolist(), oxygen=state.oxygen.tolist(),
                              debt=state.debt.tolist(), age=state.age.tolist(), solute=state.solute.tolist(),
                              pool=state.pool.tolist(), mixed_depth=state.mixed_depth)
    return res


def slab_check(st: Setting, run, site_i, years=None):
    """A 50-m mixed slab at the nearest water cell, forced as the lakes are, against the GCM's own slab ocean: with
    the exchange paired as PlaSim pairs it (the lowest level's wind and height) and as it was paired before 18 October
    (the 10 m wind and coefficient with the same air), each with the free convection over water carried to lunar
    gravity, as the lakes take it, and with Earth's constant, as PlaSim takes it."""
    import dataclasses
    site = site_info(st, site_i, water=True)
    w = windows(st, run, site, at_lake=False)
    if years:
        w = {k: (v[:years * 120] if isinstance(v, np.ndarray) else v) for k, v in w.items()}
    bsn = L.Basin.pan(50.0, 50.0, 1.0)
    scen = [L.Scenario('none', 0.0, 0.5, 0.0, 0.0)]
    albedo = float(np.median(w['gcm_alb'][w['sw'] > 50.0])) if np.any(w['sw'] > 50.0) else 0.069
    gcm = w['gcm_ts'] - 273.15
    res = dict(site=site['band'], run=run, lon=site['lon'], water_elevation_m=site['elevation'], albedo=albedo,
               gcm_mean_c=round(float(gcm.mean()), 2), gcm_range_k=round(float(np.ptp(gcm)), 2))
    for pairing, paired in (('lowest_level', True), ('ten_metres', False)):
        forcing, t, mu = hourly(st, w, site, paired=paired)
        spin, _, _ = hourly(st, w, site, 0, 120, paired=paired)
        res[pairing] = {}
        for name, air in (('lunar_free_convection', st.air),
                          ('earth_free_convection', dataclasses.replace(st.air, free_convection_gravity=9.80665))):
            ph = dataclasses.replace(physics(st, 'deep', full_mixing=True, albedo=albedo, solute='none'),
                                     solute_column=0)              # the slab's one column
            lake = L.Lake(bsn, ph, scen, air)
            state = L.initial_state(bsn, scen, float(gcm[0]), 8.0)
            lake.run(spin, state)
            out = lake.run(forcing, state)
            model = out['surface_t'].reshape(-1, STEPS_PER_MEAN).mean(1)
            res[pairing][name] = dict(model_mean_c=round(float(model.mean()), 2),
                                      bias_k=round(float((model - gcm).mean()), 2),
                                      rms_k=round(float(np.sqrt(((model - gcm) ** 2).mean())), 2),
                                      correlation=round(float(np.corrcoef(model, gcm)[0, 1]), 3),
                                      model_range_k=round(float(np.ptp(model)), 2))
    return res


def resolution_check(st: Setting, site_i, starts=None):
    """The design case's two years written every 6.4 hours, its years 20-21, drive the lakes from those means at a
    half-hour step, and from their own 3-day means with the wind factor, as the main runs see the GCM, at the same
    half-hour step and at the main runs' own step (two hours for the deep lake, one for the shallows)."""
    site = site_info(st, site_i)
    fine = windows(st, 'eddies', site, target='lock', eddy=False)
    per = st.meta['runs']['eddies']['writes_per_lunar_day'] // st.meta['runs']['lock']['writes_per_lunar_day']
    n = (fine['air'].size // per) * per
    coarse_w = {k: (v[:n].reshape(-1, per).mean(1) if isinstance(v, np.ndarray) and v.size == fine['air'].size
                    else v) for k, v in fine.items()}
    # the coarse means' wind as the main runs take it: the speed of the mean vector, raised by the factor
    f, c = st.f, site['index']
    u = f['eddies__ua'][:n, c].reshape(-1, per).mean(1)
    v = f['eddies__va'][:n, c].reshape(-1, per).mean(1)
    ta = f['eddies__ta'][:n, c].reshape(-1, per).mean(1)
    tas = f['eddies__tas'][:n, c].reshape(-1, per).mean(1)
    sigma = st.meta['runs']['eddies']['lowest_sigma']
    z_low = st.air.gas_constant * 0.5 * (ta + tas) / st.gravity * math.log(1.0 / sigma)
    coarse_w['u_low'] = np.maximum(np.hypot(u, v) * site['eddy_factor'], LOWEST_FLOOR_MS)
    coarse_w['z_low'] = z_low
    coarse_w['u10'] = np.array([L.reduce_wind(s, z, 10.0, st.gravity) for s, z in zip(coarse_w['u_low'], z_low)])
    coarse_w['window_s'] = fine['window_s'] * per
    coarse_w['hour'] = fine['hour'][per // 2::per][:coarse_w['air'].size]
    results = {}
    for name, w, spm in (('six_hourly', fine, 13), ('three_day', coarse_w, 13 * per), ('three_day_main_step',
                                                                                         coarse_w, None)):
        for kind in ('shallow', 'deep'):
            k_spm = spm if spm else (STEPS_PER_MEAN_DEEP if kind == 'deep' else STEPS_PER_MEAN)
            forcing, t, mu = hourly(st, w, site, steps_per_mean=k_spm)
            bsn = basin(kind, site)
            scen = scenarios()
            lake = L.Lake(bsn, physics(st, kind), scen, st.air)
            lake.calibrate_production(forcing)
            state = state_from(starts[kind], bsn, scen)
            lake.run(forcing, state)                    # the first pass settles the column on this forcing
            out = lake.run(forcing, state)
            results.setdefault(name, {})[kind] = summarise(kind, out, mu, bsn, scen, forcing.dt, CENTRAL)
    return results


# --- saturation -----------------------------------------------------------------------------------------------------

def saturation(st: Setting):
    p_design = DESIGN_PRESSURE_PA / 101325.0
    lunar = lambda t, p=p_design, s=None: float(L.oxygen_saturation(t, p, st.salinity if s is None else s,
                                                                      st.o2_fraction))
    earth = lambda t, s=0.0: float(L.oxygen_saturation(t, 1.0, s))
    out = dict(checks=dict(earth_fresh_20c=round(earth(20.0), 3), earth_fresh_35c=round(earth(35.0), 3),
                           expected=dict(earth_fresh_20c=9.09, earth_fresh_35c=6.95),
                           at_earlier_medians={str(t): round(lunar(t), 3) for t in (34.0, 27.7, 23.0)},
                           # the integration history review of 17 October (its finding F04, on the research
                           # drive), from the same fit at the sea medians climate_targets gave then
                           integration_history_review_f04={'34.0': 7.14, '27.7': 7.93, '23.0': 8.63}),
               pressure_atm=round(p_design, 4), o2_fraction=st.o2_fraction, salinity_g_kg=st.salinity, targets={})
    runs = st.targets['runs']
    weights = {}
    for name in ('lock', 'hot'):
        ts = st.f[f'{name}__waters_ts'] - 273.15
        ps = st.f[f'{name}__waters_ps'] * 100.0 / 101325.0
        c = L.oxygen_saturation(ts, ps, st.salinity, st.o2_fraction)
        w = np.broadcast_to(st.f[f'{name}__waters_weight'], c.shape)
        lat = np.broadcast_to(np.abs(st.f[f'{name}__waters_lat']), c.shape)
        sea = np.broadcast_to(st.f[f'{name}__waters_surface_m'] < 50.0, c.shape)
        stats = dict(all=dict(median=weighted_quantile(c.ravel(), w.ravel(), 0.5),
                              p1=weighted_quantile(c.ravel(), w.ravel(), 0.01),
                              p99=weighted_quantile(c.ravel(), w.ravel(), 0.99)),
                     seas=weighted_quantile(c[sea], w[sea], 0.5), raised_lakes=weighted_quantile(c[~sea], w[~sea], 0.5),
                     pressure_seas_atm=float(np.average(ps[sea], weights=w[sea])),
                     pressure_raised_lakes_atm=float(np.average(ps[~sea], weights=w[~sea])),
                     temperature_median_c=weighted_quantile(ts.ravel(), w.ravel(), 0.5), bands={})
        for lo, hi in ((0, 15), (15, 30), (30, 45), (45, 60), (60, 90)):
            m = (lat >= lo) & (lat < hi)
            stats['bands'][f'{lo}-{hi}'] = weighted_quantile(c[m], w[m], 0.5)
        weights[name] = stats
    for name, ct_name in TARGET_NAMES.items():
        tgt = next(t for t in st.targets['targets'] if t['name'] == ct_name)
        share = tgt['interpolated_share_of_the_way']
        mix = lambda a, b: a + share * (b - a)
        design, test = runs['design']['seas'], runs['test']['seas']
        median_c = tgt['climate']['sea_median_c']
        entry = dict(sea_median_c=round(median_c, 2), mg_l=round(lunar(median_c), 3),
                     seawater_mg_l=round(lunar(median_c, s=35.0), 3), earth_air_fresh_mg_l=round(earth(median_c), 3),
                     warm_99th_c=round(mix(design['all']['p99'], test['all']['p99']), 2),
                     warm_99th_mg_l=round(lunar(mix(design['all']['p99'], test['all']['p99'])), 3),
                     bands={})
        for band in design['bands']:
            t_band = mix(design['bands'][band]['median'], test['bands'][band]['median'])
            entry['bands'][band] = dict(median_c=round(t_band, 2), mg_l=round(lunar(t_band), 3))
        lw, hw = weights['lock'], weights['hot']
        entry['three_day_means'] = dict(
            median=round(mix(lw['all']['median'], hw['all']['median']), 3),
            p1=round(mix(lw['all']['p1'], hw['all']['p1']), 3), p99=round(mix(lw['all']['p99'], hw['all']['p99']), 3),
            seas_median=round(mix(lw['seas'], hw['seas']), 3),
            raised_lakes_median=round(mix(lw['raised_lakes'], hw['raised_lakes']), 3),
            bands={b: round(mix(lw['bands'][b], hw['bands'][b]), 3) for b in lw['bands']})
        out['targets'][name] = entry
    out['waters_pressure_atm'] = dict(seas=round(weights['lock']['pressure_seas_atm'], 4),
                                      raised_lakes=round(weights['lock']['pressure_raised_lakes_atm'], 4))
    out['waters_temperature_check_c'] = dict(lock=round(weights['lock']['temperature_median_c'], 2),
                                             hot=round(weights['hot']['temperature_median_c'], 2))
    return out


# --- the whole ------------------------------------------------------------------------------------------------------

def cache_key(st: Setting, name, **case):
    h = hashlib.sha256()
    for p in (Path(__file__), Path(L.__file__)):
        h.update(sha256(p).encode())
    h.update(json.dumps(st.inputs, sort_keys=True).encode())
    h.update(json.dumps(dict(name=name, quick=st.quick, **case), sort_keys=True).encode())
    return h.hexdigest()[:24]


def cached(st: Setting, name, fn, **case):
    folder = RUNS / 'cases'
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f'{name}_{cache_key(st, name, **case)}.json'
    if path.exists():
        return json.loads(path.read_text())
    tic = time.time()
    res = fn(st, **case)
    path.write_text(json.dumps(res) + '\n')
    shown = {k: v for k, v in case.items() if k not in ('start_from', 'starts')}
    print(f'  {name} {shown} done in {time.time() - tic:.0f} s', flush=True)
    return res


def strip(runs):
    """The runs without their saved columns."""
    return {group: {k: {kk: vv for kk, vv in v.items() if kk != 'final_state'} for k, v in cases.items()}
            for group, cases in runs.items()}


def answers(runs, group='deep'):
    """The three questions in each target, from one share's deep runs: oxygen with the surface-renewal transfer and
    demand fixed by the trophic state; the bracket spans the other transfer and the labile pool's e-folding, and the
    warm kinetics' transient stands beside it."""
    out = {}
    trophic = list(TROPHIC)
    for target in TARGET_NAMES:
        deep = {v['case']['site']: v for v in runs[group].values() if v['case']['target'] == target}
        shallow = {v['case']['site']: v for v in runs['shallow'].values() if v['case']['target'] == target}
        a = dict(overturn={}, deep_oxygen={}, overturn_oxygen={}, shallows_late_night={})
        for site, v in deep.items():
            m = v['mixing']
            a['overturn'][site] = dict(complete_overturns=m['complete_overturns'], cycles=v['cycles'],
                                       years=v['case']['years'], typical_mixing_m=m['deepest_m'][1],
                                       mixing_90th_m=m['deepest_m'][2], bottom_age_days=m['bottom_age_days'],
                                       lake_number=v.get('lake_number'))
            a['deep_oxygen'][site] = {t: dict(mean=v['oxygen'][f'{t}/renewal/supply']['deepest_tenth']['mean'],
                                              p5=v['oxygen'][f'{t}/renewal/supply']['deepest_tenth']['p5'],
                                              anoxic_share=v['oxygen'][f'{t}/renewal/supply']['deepest_tenth'][
                                                  'anoxic_share']) for t in trophic}
            a['overturn_oxygen'][site] = {t: v['oxygen'][f'{t}/renewal/supply'].get('overturns') for t in trophic}
        for site, v in shallow.items():
            entry = {}
            for t in trophic:
                n = v['oxygen'][f'{t}/renewal/supply']['nights']
                labile = {name.split('/')[-1]: o['nights']['surface_minimum'][1]
                          for name, o in v['oxygen'].items() if name.startswith(f'{t}/renewal/supply/labile')}
                lows = [v['oxygen'][f'{t}/{tr}/supply']['nights']['surface_minimum'][1] for tr in TRANSFER] \
                    + list(labile.values())
                warm = [v['oxygen'][f'{t}/{tr}/warm']['nights']['surface_minimum'][1] for tr in TRANSFER]
                entry[t] = dict(median=n['surface_minimum'][1], p10=n['surface_minimum'][0],
                                share_of_saturation=n['surface_minimum_share_of_saturation'][1],
                                share_of_the_night=n['surface_minimum_share_of_the_night'][1],
                                before_dawn=n['surface_before_dawn'][1],
                                hypoxic_nights=n['hypoxic_share'], anoxic_nights=n['anoxic_share'],
                                bottom_hypoxic_hours=n['hours_bottom_hypoxic'][1],
                                bracket_of_medians=[min(lows), max(lows)], labile_brackets=labile,
                                warm_transient=[min(warm), max(warm)])
            a['shallows_late_night'][site] = entry
        out[target] = a
    return out


_SETTINGS = {}


def _setting(quick):
    if quick not in _SETTINGS:
        _SETTINGS[quick] = Setting(quick=quick)
    return _SETTINGS[quick]


FUNCTIONS = dict(run_case=run_case, slab_check=slab_check, resolution_check=resolution_check)


def _work(item):
    """One cached case in a worker process: (group, key, cache name, function name, case, quick)."""
    group, key, name, fn, case, quick = item
    return group, key, cached(_setting(quick), name, FUNCTIONS[fn], **case)


def main(argv=None) -> int:
    from concurrent.futures import ProcessPoolExecutor
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true', help='two years of forcing; writes to research/runs')
    ap.add_argument('--jobs', type=int, default=2, help='worker processes, at most 2 on the shared machine')
    args = ap.parse_args(argv)
    st = _setting(args.quick)
    years = 2 if args.quick else None
    jobs = max(1, min(args.jobs, 2))
    tic = time.time()
    sat = saturation(st)
    groups = ('deep', 'deep_eta_low', 'deep_eta_high', 'deep_lake_number', 'deep_fixed', 'shallow', 'deepest',
              'deepest_long', 'deepest_long_eta_low', 'deepest_long_eta_high', 'wind_old_cap', 'wind_cap_raised',
              'sensitivity', 'resort_from_hot')
    runs = {g: {} for g in groups}
    eta_group = dict(central='deep', low='deep_eta_low', high='deep_eta_high')
    first = []
    for eta, group in eta_group.items():
        for target in ('lock', 'resort', 'hot'):
            for site_i in range(3):
                case = dict(kind='deep', target=target, site_i=site_i, years=years)
                if eta != 'central':
                    case['eta'] = eta
                first.append((group, f'{target}/{SITE_BANDS[site_i]}', 'deep', 'run_case', case, args.quick))
    for target in ('lock', 'resort', 'hot'):
        # Yeates & Imberger's diffusivity from the Lake number, and the fixed dissipation of the first round with its
        # cap, each from the warm start, at the 0-15 lake
        first.append(('deep_lake_number', f'{target}/0-15', 'deep', 'run_case',
                      dict(kind='deep', target=target, site_i=0, years=years, mixing='lake_number'), args.quick))
        for diss in ('central', 'high'):
            first.append(('deep_fixed', f'{target}/0-15/{diss}', 'deep', 'run_case',
                          dict(kind='deep', target=target, site_i=0, years=years, mixing='dissipation',
                               dissipation=diss), args.quick))
        for site_i in range(3):
            first.append(('shallow', f'{target}/{SITE_BANDS[site_i]}', 'shallow', 'run_case',
                          dict(kind='shallow', target=target, site_i=site_i, years=years), args.quick))
    slabs = [('slab', f'{run}/{SITE_BANDS[i]}', 'slab', 'slab_check', dict(run=run, site_i=i, years=years),
              args.quick) for run in ('lock', 'hot') for i in range(3)]
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        slab = {}
        for group, key, res in pool.map(_work, first + slabs):
            (slab if group == 'slab' else runs[group])[key] = res
        second = []
        for eta, group in eta_group.items():
            suffix = '' if eta == 'central' else f'_eta_{eta}'
            extra = {} if eta == 'central' else dict(eta=eta)
            for target in ('lock', 'resort', 'hot'):
                base = runs[group][f'{target}/0-15']['final_state']
                if eta == 'central':
                    second.append(('deepest', target, 'deepest', 'run_case',
                                   dict(kind='deep', target=target, site_i=0, depth=AREA_MEDIAN_DEEPEST_M,
                                        years=years, start_from=base, keep=CENTRAL), args.quick))
                second.append((f'deepest_long{suffix}', target, 'deepest', 'run_case',
                               dict(kind='deep', target=target, site_i=0, depth=AREA_MEDIAN_DEEPEST_M, years=years,
                                    start_from=base, keep=CENTRAL, long_spin=True, **extra), args.quick))
        for target in ('lock', 'resort', 'hot'):
            base = runs['deep'][f'{target}/0-15']['final_state']
            second.append(('wind_old_cap', target, 'deep', 'run_case',
                           dict(kind='deep', target=target, site_i=0, years=years, start_from=base, keep=CENTRAL,
                                cap='central'), args.quick))
            second.append(('wind_cap_raised', target, 'deep', 'run_case',
                           dict(kind='deep', target=target, site_i=0, years=years, start_from=base, keep=CENTRAL,
                                cap=CAP_WIND_RAISED), args.quick))
        for target in ('lock', 'hot'):
            base = runs['deep'][f'{target}/0-15']['final_state']
            common = dict(kind='deep', target=target, site_i=0, years=years, start_from=base, keep=CENTRAL)
            for solute in ('none', 'carbonate', 'iron'):
                second.append(('sensitivity', f'{target}/solute_{solute}', 'deep', 'run_case',
                               dict(common, solute=solute), args.quick))
            second.append(('sensitivity', f'{target}/eutrophic', 'deep', 'run_case',
                           dict(common, coupled='eutrophic/renewal/supply'), args.quick))
        for kind in ('shallow', 'deep'):
            base = runs[kind]['resort/0-15']['final_state']
            second.append(('resort_from_hot', kind, kind, 'run_case',
                           dict(kind=kind, target='resort', site_i=0, resort_from='hot', years=years,
                                start_from=base, keep=CENTRAL), args.quick))
        starts = {k: runs[k]['lock/0-15']['final_state'] for k in ('shallow', 'deep')}
        second.append(('resolution', 'resolution', 'resolution', 'resolution_check',
                       dict(site_i=0, starts=starts), args.quick))
        resolution = None
        for group, key, res in pool.map(_work, second):
            if group == 'resolution':
                resolution = res
            else:
                runs[group][key] = res
    checks = dict(slab=[slab[f'{run}/{SITE_BANDS[i]}'] for run in ('lock', 'hot') for i in range(3)],
                  resolution=resolution)
    out = dict(
        schema=SCHEMA,
        producer=dict(files={str(p.resolve().relative_to(ROOT)): sha256(p)
                             for p in (Path(__file__), Path(L.__file__), HERE / 'lakes.py', HERE / 'extract.py')},
                      inputs=st.inputs,
                      constants=dict(STEFAN_BOLTZMANN=STEFAN_BOLTZMANN, SYNODIC_MONTH_DAYS=SYNODIC_MONTH_DAYS,
                                     GAS_CONSTANT=GAS_CONSTANT, MOON_SURFACE_GRAVITY=L.MOON_SURFACE_GRAVITY,
                                     EXOPLASIM_CHARNOCK=L.EXOPLASIM_CHARNOCK,
                                     EXOPLASIM_SEA_ROUGHNESS_FLOOR=L.EXOPLASIM_SEA_ROUGHNESS_FLOOR,
                                     EXOPLASIM_VON_KARMAN=L.EXOPLASIM_VON_KARMAN, BOLTZMANN=L.BOLTZMANN,
                                     ELEMENTARY_CHARGE=L.ELEMENTARY_CHARGE)),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        setting=dict(quick=args.quick, design_pressure_pa=DESIGN_PRESSURE_PA, o2_fraction=st.o2_fraction,
                     salinity_g_kg=st.salinity, held_ph=st.held_ph, heat_flow_w_m2=st.heat_flow,
                     air_gas_constant=st.air.gas_constant, air_cp=st.air.cp, gravity=st.gravity,
                     step_s=st.dt, share_of_the_way=st.share, lapse_k_per_km={k: round(v * 1000.0, 3)
                                                                             for k, v in st.lapse.items()},
                     subsolar_longitude_at_zero=st.phase,
                     eddy_wind_factor=[round(float(x), 3) for x in st.eddy_factor],
                     deep_m=DEEP_M, area_median_deepest_m=AREA_MEDIAN_DEEPEST_M, shallow_m=SHALLOW_M,
                     spin_up=dict(coarse_steps=COARSE, warm_start_k=WARM_START_K, max_passes=MAX_PASSES,
                                  passes_at_the_analysis_step=FINE_PASSES,
                                  settled=dict(k=SETTLED_K, o2_g_m3=SETTLED_O2, age_share=SETTLED_AGE,
                                               solute_g_m3=SETTLED_SOLUTE)),
                     trophic=TROPHIC, transfer=list(TRANSFER), kinetics_ev=KINETICS, warm_reference_c=SOURCE_C,
                     labile_days=LABILE_DAYS, methane_share=METHANE_SHARE, iron=IRON, solute_density=SOLUTE,
                     debt_density=L.anaerobic_density(METHANE_SHARE), coupled=COUPLED,
                     steps_per_mean=dict(deep=STEPS_PER_MEAN_DEEP, shallow=STEPS_PER_MEAN),
                     deep_mixing=DEEP_MIXING, wind_share=WIND_SHARE, cap_wind_m2_s=CAP_WIND,
                     cap_wind_raised_m2_s=CAP_WIND_RAISED,
                     fixed_dissipation_w_kg=DISSIPATION, fixed_cap_m2_s=CAP,
                     unstable_diffusivity_m2_s=UNSTABLE_DIFFUSIVITY, unstable_rule='at least the cap',
                     lake_number_constant=L.Physics().lake_number_constant,
                     lake_number_z_t_below_m=LAKE_NUMBER_BELOW_M,
                     pairing='the lowest level\'s wind and height (PlaSim\'s), with its potential temperature and '
                             'humidity carried to the lake', hypoxic_g_m3=HYPOXIC, anoxic_g_m3=ANOXIC),
        saturation=sat, runs=strip(runs), checks=checks,
        answers=dict(central=answers(runs, 'deep'), eta_low=answers(runs, 'deep_eta_low'),
                     eta_high=answers(runs, 'deep_eta_high')),
        seconds=round(time.time() - tic, 1))
    path = OUT if not args.quick else RUNS / 'lake_oxygen_quick.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(dumps(out) + '\n')
    print('wrote', path.relative_to(ROOT) if path.is_relative_to(ROOT) else path, f'in {time.time() - tic:.0f} s')
    return 0


EVIDENCE = (
    'Saturation is the Benson and Krause fit (Garcia & Gordon 1992) with Benson & Krause\'s pressure factor, extended '
    'from its stated 1.1 atm to the Moon\'s 1.2, for the air\'s 17.46% O2 and the held fresh seas\' salinity. The '
    'lakes are a one-dimensional model (biosphere/ecology/lake.py) checked against analytic diffusion, convective '
    'and wind deepening, reaeration and conduction-limited ice, with heat and oxygen conserved to round-off, forced '
    'from the GCM\'s 3-day means at the T21 cells holding three crater lakes: the design case (the lock) and the '
    'test with the longwave refitted by height and ice optics (the hot world), years 10-19, the resort interpolated '
    'between them by the ground\'s temperature, which mixes two physics. The surface exchange is paired as PlaSim '
    'pairs it, at the lowest level\'s height with its wind, potential temperature and humidity, as corrected after '
    'the independent review of 18 October and its re-check. Below the mixed layer the dissipation is the share of '
    'the wind\'s energy flux at 10 m that Earth\'s lakes pass to their stratified water (Wuest, Piepke & Van Senden '
    '2000), spread over the water\'s mass there and bracketed by the share\'s measured range, with the flux on Wuest '
    '& Lorke\'s lake drag as Simstrat pairs it (the drag the measured share was taken against is open; the lake '
    'model\'s own stress at lunar gravity gives about half to four fifths of that flux at these light winds); the '
    'fixed dissipation and Yeates & Imberger\'s Lake-number diffusivity are sensitivity rows, and the cap on the '
    'diffusivity holds no interface back. Each deep lake settles at the analysis step after its coarse spin-up. Oxygen '
    'demand is bounded scenarios from published lakes on Earth (design values where labelled), fixed by the supply '
    'of organic matter; the warm kinetics are a transient beside the bounds. The lunar gravity enters every buoyancy '
    'term. The runs are calculations on one GCM\'s climate at T21; no lunar lake has been measured.')
READING_RULE = (
    'Oxygen in g/m3 (mg/L); quantiles are lists of the 10th, 50th and 90th percentiles. saturation.targets gives '
    'each target\'s seas at their median and band temperatures at 1.2 atm, and three_day_means the distribution over '
    'the GCM\'s water cells at each cell\'s own pressure, the resort interpolated by the share of the way as '
    'climate_targets does. runs.deep and runs.shallow are keyed target/band, under the wind\'s deep mixing at its '
    'central share of P10 (0.3%); runs.deep_eta_low and deep_eta_high hold the deep lakes at its bracket (0.2% and '
    '0.4%); runs.deep_lake_number holds the 0-15 lakes under Yeates & Imberger\'s diffusivity from the Lake number, '
    'and runs.deep_fixed under the first round\'s fixed dissipation (central 1e-10 and high 1e-9 W/kg, capped at '
    '1e-4 m2/s), keyed target/band/dissipation; runs.wind_old_cap the wind\'s mixing under the old cap of 1e-4 m2/s, '
    'and runs.wind_cap_raised under a cap ten times the wind\'s, 1 m2/s, the check that its cap holds nothing back. '
    'Each case\'s deep_mixing gives the treatment, the dissipation the wind gave (W/kg) and the shares of the '
    'interfaces below the mixed layer at the cap, held back by it (the cap slower than a step\'s own mixing) and '
    'unstable. Each case\'s oxygen is keyed trophic/transfer/kinetics: transfer \'renewal\' (Tedford et al. 2014, '
    'central) or \'cole_caraco\'; kinetics \'supply\' (demand fixed by the trophic state\'s supply, the bounds) or '
    '\'warm\' (the deep water\'s and floor\'s demand with 0.62 eV from the source hypolimnia\'s 6 C, production and '
    'the euphotic water\'s labile pool held: a transient no settling supply sustains, kept out of the bounds), and '
    'trophic/renewal/supply/labile_2d and labile_7d carry the labile pool\'s brackets about its central 4 days. The '
    'variant runs keep the central renewal/supply columns. mixing.complete_overturns counts the lunar cycles, sunrise '
    'to sunrise at the site, in which the mixed layer reached the bottom layer, and deepest_m each cycle\'s deepest '
    'mixing; deepest_tenth is the deepest tenth of the lake\'s volume; overturns reads the oxygen after each complete '
    'overturn; lake_number is the Lake number from the column every three days with z_T at the strongest '
    'stratification below 10 m (month, the lunar month\'s thermocline) and below 2 m (day, the day\'s warm layer). '
    'nights run from sunset to sunrise and days from sunrise to sunset; surface_minimum is each night\'s lowest '
    'surface oxygen and surface_minimum_share_of_the_night when it falls. settling.spin_up records each spin-up pass '
    'in coarse steps and, under at_the_analysis_step, each pass at the analysis step after them (the deep water\'s '
    'temperature, oxygen, age and solute, and whether the pass overturned). runs.sensitivity varies the solutes and '
    'the trophic state the physics follows at the 0-15 lake; runs.deepest takes that lake to the 1,925 m deepest '
    'point that half the lake area exceeds, oxygenated at the start as if just overturned, and runs the record, and '
    'runs.deepest_long (with its eta bracket) first passes over the record until it settles (up to 30 passes, then '
    'up to 10 at the analysis step); runs.resort_from_hot builds the resort from the hot world\'s weather. '
    'checks.slab compares a 50-m mixed slab with the GCM\'s own slab ocean, paired at the lowest level and at 10 m, '
    'each with the free convection at lunar gravity and at Earth\'s; checks.resolution drives the lakes from the '
    'design case\'s 6.4-hour means of its years 20-21 at a half-hour step (six_hourly), from their 3-day means at '
    'that step (three_day) and at the main runs\' steps, two hours for the deep lake and one for the shallows '
    '(three_day_main_step). answers holds the central runs and the eta bracket\'s.')


if __name__ == '__main__':
    sys.exit(main())
