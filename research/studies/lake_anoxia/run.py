"""What the deep lakes' anoxic water does, and what design can do about it, in the three target climates.

    python -m research.studies.lake_anoxia.profiles      # the deep lakes' columns through time (about nine minutes)
    python -m research.studies.lake_anoxia.run           # the analysis (about three minutes)
    python -m research.studies.lake_anoxia.run --cache-only   # the same, stopping unless the columns are current

It reads the columns that profiles.py saved from the lake model with the wind's deep mixing at its central share and
its bracket (research/runs/lake_anoxia), the waters'
oxygen study's product and lakes, geography's rain-fed lakes on LOLA, the domain synthesis's methane ledger, the deep
column's reduced-gas forcing, the climate targets as the lake studies ran on them (lake_oxygen's
climate_targets_as_run.json), the sea chemistry's lakes and the land's phosphorus, and writes
results/lake_anoxia.json:
- methane: what the anoxic depths make, by trophic state and target, the bubbles that leave the sediment at the
  Moon's hydrostatic pressure and how much of them reaches the air, the lakes' flux Moon-wide, and what it does to
  the air's methane and to the dimming that answers the gases;
- the gas build-up against the hydrostatic bubble point, the criterion behind Lake Nyos's eruption and Lake Kivu's
  hazard assessments: in the 352 m lakes between overturns, in the 1,925 m lake left unmixed, and for any isolated
  hypolimnion by depth;
- phosphorus released from anoxic floors at Earth's rates, with the Moon's low sulphate, and its feedback on the
  lakes' trophic state through Vollenweider's loading relation, its settling bracketed for lunar gravity and the
  internal load bounded by the steady state;
- the remedies, sized per km2 and for the Moon's lakes.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import sys
import time
from pathlib import Path

import numpy as np

from biosphere.ecology import lake as L
from research.studies.lake_anoxia import gases as G
from research.studies.lake_oxygen import run as R
from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUNS = ROOT / 'research' / 'runs' / 'lake_anoxia'
PROFILES = RUNS / 'profiles.json'
OXYGEN = ROOT / 'research' / 'studies' / 'lake_oxygen' / 'results' / 'lake_oxygen.json'
LAKES = ROOT / 'research' / 'studies' / 'lake_oxygen' / 'results' / 'lakes.json'
LEDGERS = ROOT / 'research' / 'studies' / 'domain_synthesis' / 'results' / 'ledgers.json'
GAS_FORCING = ROOT / 'atmosphere' / 'radiative_convective' / 'results' / 'reduced_gases.json'
TARGETS = R.TARGETS                    # the climate targets as the lake studies ran on them (lake_oxygen's run.py)
CHEMISTRY = ROOT / 'research' / 'studies' / 'sea_chemistry' / 'results' / 'sea_chemistry.json'
DRAINAGE = ROOT / 'geography' / 'results' / 'drainage.json'
OUT = HERE / 'results' / 'lake_anoxia.json'
SCHEMA = 'terluna.research.lake-anoxia/2'

TRAITS = ('oligotrophic', 'mesotrophic', 'eutrophic', 'new_reservoir')
TARGET_ORDER = ('lock', 'resort', 'hot')
CENTRAL = 'renewal/supply'
# The share of anaerobic mineralisation (in oxygen equivalents) that ends as methane. Low-sulphate fresh water leaves
# methanogenesis the main terminal path once iron and manganese oxides are spent: over half of the organic matter in
# ferruginous Lake Matano (Crowe et al. 2011); the Moon's lakes hold about 1.4 uM of sulphate (sea chemistry's export).
METHANE_SHARE = dict(low=0.3, central=0.5, high=0.8)
HALF_SATURATION = 0.5          # g/m3: the lake model's respiration slows toward it, the unmet share K/(K + O2)
ANOXIC = 0.5                   # g/m3: anaerobic mineralisation is counted in water below it; between 0.5 and 2 the
                               # model's Monod limit owes a fifth to half of the demand and repays it in place
ETAS = ('central', 'low', 'high')                 # the wind's share of P10 below the mixed layer and its bracket
GIANT_M = 1000.0               # lakes deeper than this take the oxycline's flux alone in the giants' bracket
FILL_TO_M = 12.0               # the depth the nights ventilate, to which the fill brings a lake's floor
# Share of the methane reaching oxygenated water at an overturn that methanotrophs oxidise: at Earth's autumn overturns
# 54-98%, most estimates 75-94% (Kankaala et al. 2007; Schubert et al. 2012).
OXIDISED = dict(low=0.5, central=0.85, high=0.98)
# Share of the methane crossing a lasting oxycline between overturns that escapes: Kivu's 0.2 of the 35 g C m-2 a-1
# reaching its oxycline (Pasche et al. 2011) central, the bracket from all of it oxidised to the 15% Earth's overturns
# let through; the methane study refines it.
OXYCLINE_ESCAPE = dict(low=0.0, central=0.2 / 35.0, high=0.15)
OVERTURN_UNMET_ESCAPE = 0.5    # stored methane the overturned lake's oxygen cannot meet: half escapes (a guess)
RELEASE_RULES = ('net_rise', 'stock')   # what an overturn releases: the lower and upper bounds (deep_lake())
# The deep water's age at overturn and the water that last left the surface saturated: dissolved N2 and Ar at the
# surface's saturation, CO2 at equilibrium with the design air.
DESIGN_CO2_PPM = 400.0
SURFACE_ATM = R.DESIGN_PRESSURE_PA / 101325.0
LAKE_IONIC = 1.0e-3            # mol/kg, a soft lake (sea chemistry: 0.03-0.04 g/kg)
REDFIELD_P_PER_O2 = 31.0 / (138.0 * 32.0) * 1000.0      # mg P per g O2 respired
# Bubbles: Earth's lake bubbles are 2-8 mm across (DelSontro et al. 2015); lunar gravity's weaker buoyancy lets them
# grow larger before they break free of the mud (fracture scaling, a guess), so 5-20 mm brackets them.
BUBBLE_MM = dict(small=5.0, central=10.0, large=20.0)
BUBBLE_CH4_FRACTION = 0.8      # Langenegger et al. (2019): sediment bubbles 75-90% methane
SEDIMENT_DECAY_PER_M = 30.0    # Langenegger et al. (2019): production falls off at 10-43 per m, mean 30
# Phosphorus from anoxic sediments (Nurnberg 1997, as cited in the Canadian review and later work), mg P m-2 d-1: the
# means of her trophic classes, taken as inputs, so release switches on at the class boundary by construction. New
# reservoirs take the eutrophic class with their production, and the hypereutrophic 25.6 as the bracket for their
# flooded soils. Low sulphate leaves iron to bind phosphate as vivianite (Caraco et al. 1989), a factor of 0.1-1.
P_RELEASE = dict(oligotrophic=0.0, mesotrophic=6.1, eutrophic=10.8, new_reservoir=10.8)
P_RELEASE_FLOODED_SOILS = 25.6
# Vollenweider's retention, TP = L / (q_s (1 + tau^0.5)), stands for settling at Earth's gravity: an apparent settling
# speed of q_s tau^0.5. Carried to the Moon it slows with gravity, as Stokes settling does (g/6) or as aggregates do
# (about its square root). The internal load is added outside the retention, as Nurnberg's (1984) mass balance has it,
# for lakes whose floors hold phosphorus beyond their steady state; at steady state the floors return at most what
# settled, so the lake keeps at most what flows in, L / q_s: anoxic floors can at most cancel the retention, and
# burial, the vivianite among it, keeps part of it.
SETTLING = dict(earth=1.0, aggregates=(MOON_SURFACE_GRAVITY / STANDARD_GRAVITY) ** 0.5,
                stokes=MOON_SURFACE_GRAVITY / STANDARD_GRAVITY)
EARTH_AIR = dict(o2=0.20946, n2=0.78084, ar=0.00934, co2_ppm=420.0)   # dry air, mole fractions
LOW_SULPHATE_FACTOR = (0.1, 1.0)
ANOXIC_FLOOR = 1.0             # g/m3 above the floor counted as anoxic for the release
HYPOLIMNION_M = 20.0           # the bubble point is read below this depth and the euphotic zone, where the
                               # day's photosynthesis supersaturates the water with oxygen
OECD_TP = dict(ultraoligotrophic=4.0, oligotrophic=10.0, mesotrophic=35.0, eutrophic=100.0)   # mg/m3, upper bounds
# Remedies (sources.json): oxygen by vacuum swing adsorption at 0.29-0.34 kWh per Nm3 (1.429 kg); compression to the
# diffusers isothermal at 60% efficiency; Lorenzen & Fast's (1977) 9.2 m3 of air a minute per km2 for circulation;
# oxygenation's demand raised 1-3 times by the stirring (Moore et al. 2015 used three).
O2_KWH_PER_KG = (0.29 / 1.429, 0.34 / 1.429)
COMPRESSOR_EFFICIENCY = 0.6
CIRCULATION_AIR_M3_MIN_KM2 = 9.2
INDUCED_DEMAND = (1.0, 3.0)
EARTH_INSTALLATIONS_T_DAY_KM2 = dict(Camanche=9.0 / 31.16, Richard_B_Russell=200.0 / 107.85,
                                     Upper_San_Leandro=9.0 / 2.509, Spring_Hollow=0.3 / 0.639, Carvins_Cove=2.0 / 3.237)
CH4_ENERGY_KJ_MOL = 890.0
TG = 1e12


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def dumps(obj) -> str:
    text = json.dumps(obj, indent=1)
    return re.sub(r'\[\s+([-0-9.eE,\s]+?|null|true|false)\s+\]', lambda m: '[' + ' '.join(m.group(1).split()) + ']',
                  text)


def rnd(x, sig=4):
    if isinstance(x, dict):
        return {k: rnd(v, sig) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [rnd(v, sig) for v in x]
    if isinstance(x, (np.floating, np.integer)):
        x = x.item()
    if isinstance(x, bool) or not isinstance(x, float) or x == 0.0 or not math.isfinite(x):
        return x
    return round(x, sig - 1 - int(math.floor(math.log10(abs(x)))))


# --- the setting ------------------------------------------------------------------------------------------------

class Setting:
    def __init__(self):
        self.oxygen = json.loads(OXYGEN.read_text())
        self.meta = json.loads(PROFILES.read_text())
        for name, case in self.meta['cases'].items():
            if sha256(RUNS / f'{name}.npz') != case['npz_sha256']:
                raise SystemExit(f'{name}.npz does not match profiles.json; run profiles.py again')
        s = self.oxygen['setting']
        self.gravity, self.salinity, self.o2_fraction = s['gravity'], s['salinity_g_kg'], s['o2_fraction']
        self.scenarios = R.scenarios()
        self.names = [x.name for x in self.scenarios]
        self.ledgers = json.loads(LEDGERS.read_text())
        self.forcing = json.loads(GAS_FORCING.read_text())['moon_at_design']
        self.targets = json.loads(TARGETS.read_text())
        chem = json.loads(CHEMISTRY.read_text())
        held = next(r for r in chem['runs'] if all(r[k] == v for k, v in R.HELD_FRESH.items()))
        self.lake_ph = float(held['rows'][-1]['lake_ph'])
        self.p_loads = chem['setting']['p_kg_yr']
        self.pco2 = DESIGN_CO2_PPM * 1e-6 * SURFACE_ATM
        self.ta = G.alkalinity_at(self.lake_ph, self.pco2, 21.15, LAKE_IONIC)[0]
        drainage = json.loads(DRAINAGE.read_text())
        self.lake_area_km2 = json.loads(LAKES.read_text())['lakes']['lake_area_km2']
        self.lake_volume_km3 = drainage['lakes']['volume_km3']
        self.residence_yr = 445.0                        # sea chemistry: the lakes pass the runoff through 445 years
        self.n2_frac = 0.8153104493254115                # earthlike_air at 1.2 atm and 400 ppm (lake_oxygen's air)
        self.ar_frac = 0.009711049153081852
        self.inputs = {str(p.relative_to(ROOT)): sha256(p) for p in (OXYGEN, LAKES, LEDGERS, GAS_FORCING, TARGETS,
                                                                       CHEMISTRY, DRAINAGE, PROFILES)}
        self.absorbed = {}
        share = {t['name']: t['interpolated_share_of_the_way'] for t in self.targets['targets']}
        for k, name in R.TARGET_NAMES.items():
            # absorbed sunlight between the design case's 211.2 W/m2 and the warm Moon's 225.4 (climate_targets)
            self.absorbed[k] = 211.2 + share[name] * (225.4 - 211.2)

    def load(self, name):
        return dict(np.load(RUNS / f'{name}.npz'))


def air_gases(st, t_c):
    """N2 and Ar (mol/m3) in water that left the surface saturated at t_c, and the vapour's pressure."""
    dry = SURFACE_ATM - L.vapour_pressure_atm(t_c)
    return G.solubility('N2', t_c) * st.n2_frac * dry, G.solubility('Ar', t_c) * st.ar_frac * dry


# --- methane and gas in the deep lakes --------------------------------------------------------------------------

def column_gases(st, rec, bsn, m, share):
    """Methane, CO2 and the gas saturation ratio through the record for oxygen column m: methane the debt's methane
    share, CO2 what the water respired since it left the surface (its oxygen deficit and its debt), below the
    euphotic zone; N2 and Ar at the surface's saturation."""
    t = rec['t'].astype(float)
    ox = rec['oxygen'][:, :, m].astype(float)
    debt = rec['debt'][:, :, m].astype(float)
    sat = L.oxygen_saturation(t, SURFACE_ATM, st.salinity, st.o2_fraction)
    z_eu = 4.6 / st.scenarios[m].attenuation
    below = bsn.centre >= z_eu
    ch4 = share * debt / 64.0
    respired = np.where(below, np.maximum(sat - ox, 0.0), 0.0) / 32.0
    dic_added = respired + debt * (share / 64.0 + (1.0 - share) / 32.0)
    base = equilibrium_dic(st, t)
    _, co2 = G.speciate(base + dic_added / 1000.0, st.ta, t, LAKE_IONIC)
    n2, ar = air_gases(st, t)
    p_gas = G.gas_pressure_atm(t, ch4=ch4, co2_aq=co2, n2=n2, ar=ar, o2=ox / 32.0)
    p_local = G.local_pressure_atm(bsn.centre, SURFACE_ATM, st.gravity)
    return dict(ch4=ch4, dic_added=dic_added, ratio=p_gas / p_local[None, :], p_gas=p_gas, p_local=p_local)


def equilibrium_dic(st, t_c, pco2=None):
    """DIC (mol/kg) of lake water at the lake alkalinity in equilibrium with the design air at t_c (or with pco2)."""
    t = np.asarray(t_c, dtype=float)
    k1, k2, kw, _ = G.carbonate_constants(t)
    g1, g2 = G.gammas(21.15, LAKE_IONIC)
    k1, k2, kw = k1 / (g1 * g1), k2 / g2, kw / (g1 * g1)
    co2 = G.co2_kh(t) * (st.pco2 if pco2 is None else pco2)
    lo, hi = np.full(t.shape, -11.0), np.full(t.shape, -4.0)
    for _ in range(70):
        mid = 0.5 * (lo + hi)
        h = 10.0 ** mid
        alk = co2 * (k1 / h + 2.0 * k1 * k2 / (h * h)) + kw / h - h
        high = alk > st.ta
        lo, hi = np.where(high, mid, lo), np.where(high, hi, mid)
    h = 10.0 ** (0.5 * (lo + hi))
    return co2 * (1.0 + k1 / h + k1 * k2 / (h * h))


def anoxic_methane(st, rec, bsn, m, share, target_t=None):
    """Methane made in anoxic water and on anoxic floors (mol per m2 of lake surface a year): the model's demand in
    each layer times the share respiration cannot meet, K/(K + O2), averaged over the record's daily columns."""
    s = st.scenarios[m]
    demand = L.demand_profile(bsn, s)                       # g O2 m-3 d-1 as the scenario gives it
    t = rec['t'].astype(float)
    if s.activation_ev:                                     # the warm bracket, from the scenario's reference
        demand = demand[None, :] * np.exp(s.activation_ev / L.BOLTZMANN_EV
                                          * (1.0 / (s.reference_c + 273.15) - 1.0 / (t + 273.15)))
    else:
        demand = np.broadcast_to(demand, t.shape)
    ox = np.maximum(rec['oxygen'][:, :, m].astype(float), 0.0)
    # the model's Monod form leaves part of the demand unmet in oxygenated and hypoxic water and repays it in place
    # within a day, which is aerobic respiration in a lake and, on a floor, methane oxidised at the sediment's
    # surface: only anoxic water counts
    unmet = np.where(ox < ANOXIC, demand * HALF_SATURATION / (HALF_SATURATION + ox), 0.0)
    per_day = (unmet * bsn.volume[None, :]).sum(1).mean() / bsn.area[0]       # g O2 m-2 d-1
    return share * per_day * 365.25 / 64.0, per_day


def overturns(rec, bsn, dt, every):
    """Steps at which the mixed layer reaches the bottom layer, grouped into events a lunar day apart."""
    md = rec['mixed_depth']
    hit = np.flatnonzero(md >= bsn.z[-2] - 1e-6)
    if not hit.size:
        return []
    events, start = [], hit[0]
    for a, b in zip(hit[:-1], hit[1:]):
        if (b - a) * dt > 15 * 86400.0:
            events.append((start, a))
            start = b
    events.append((start, hit[-1]))
    return events


def deep_lake(st, name, case):
    rec = st.load(name)
    z, area = rec['z'], rec['area']
    bsn = L.Basin(z, area)
    dt, every = float(rec['dt']), case['every']
    days = len(rec['mixed_depth']) * dt / 86400.0
    years = days / 365.25
    events = overturns(rec, bsn, dt, every)
    prof_step = rec['step']
    out = dict(target=case['target'], band=case['band'], years=round(years, 2), overturn_events=len(events),
               traits={})
    for trophic in TRAITS:
        m = st.names.index(f'{trophic}/{CENTRAL}')
        mw = st.names.index(f'{trophic}/renewal/warm')
        gas = column_gases(st, rec, bsn, m, METHANE_SHARE['central'])
        made = {k: anoxic_methane(st, rec, bsn, m, v)[0] for k, v in METHANE_SHARE.items()}
        made_warm = anoxic_methane(st, rec, bsn, mw, METHANE_SHARE['central'])[0]
        inventory = (gas['ch4'] * bsn.volume[None, :]).sum(1) / bsn.area[0]          # mol m-2 through time
        before, after = [], []
        for a, b in events:
            before.append(float(inventory[max(int(np.searchsorted(prof_step, a)) - 1, 0)]))
            after.append(float(inventory[min(int(np.searchsorted(prof_step, b)) + 1, len(prof_step) - 1)]))
        # What an overturn releases, by two rules until the methane tracer settles it (the record repeating, the
        # first event's gap runs from the record's last). The net rise: the inventory's rise since the previous event,
        # which in a repeating record is what the model's own oxygen removed during the overturns, and so a lower
        # bound. The stock: the standing stock the overturn brings up, no more than the methane made since the
        # previous event, the study's own assumption taken whole and so an upper bound. The share of the stock the
        # restratifying mixed layer degasses lies between them.
        gaps, net_rise, stock = [], [], []
        for i, (a, _) in enumerate(events):
            gap = ((a - events[i - 1][1]) if i else (len(rec['mixed_depth']) - events[-1][1] + a)) * dt / 86400.0
            gaps.append(gap)
            net_rise.append(min(max(before[i] - after[i - 1], 0.0), made['central'] * gap / 365.25))
            stock.append(min(before[i], made['central'] * gap / 365.25))
        released = dict(net_rise=net_rise, stock=stock)
        deep_rows = np.flatnonzero(bsn.centre >= max(HYPOLIMNION_M, 4.6 / st.scenarios[m].attenuation))
        sub = gas['ratio'][:, deep_rows]
        worst = np.unravel_index(np.argmax(sub), sub.shape)
        worst = (worst[0], deep_rows[worst[1]])
        ox_after = []
        for a, _ in events:
            k = min(int(np.searchsorted(prof_step, a)) + 1, len(prof_step) - 1)
            ox_after.append(float((rec['oxygen'][k, :, m] * bsn.volume).sum() / bsn.volume.sum()))
        need = {rule: [64.0 * b for b in x] for rule, x in released.items()}   # g O2 m-2 to oxidise it
        # g O2 m-2 in the overturned lake, read from the saved column after the overturn, when the model has spent a
        # little of it on the debt already: a slight underestimate of what meets the methane
        have = [o * bsn.volume.sum() / bsn.area[0] for o in ox_after]
        out['traits'][trophic] = dict(
            methane_made_mol_m2_yr=made, methane_made_warm_transient_mol_m2_yr=made_warm,
            inventory_before_overturns_mol_m2=before, inventory_after_overturns_mol_m2=after,
            days_since_the_previous_overturn=gaps, released_at_overturns_mol_m2=released,
            inventory_mean_mol_m2=float(inventory.mean()), inventory_max_mol_m2=float(inventory.max()),
            bottom_ch4_max_umol_l=float(gas['ch4'][:, -1].max() * 1000.0),
            bottom_co2_added_max_mmol_l=float(gas['dic_added'][:, -1].max()),
            saturation_ratio_max=float(gas['ratio'][worst]), saturation_ratio_depth_m=float(bsn.centre[worst[1]]),
            saturation_ratio_bottom_max=float(gas['ratio'][:, -1].max()),
            gas_pressure_bottom_max_atm=float(gas['p_gas'][:, -1].max()),
            bottom_pressure_atm=float(gas['p_local'][-1]),
            oxygen_after_overturns_g_m2=have, oxygen_to_oxidise_stored_methane_g_m2=need,
            anoxic_floor=anoxic_floor(rec, bsn, m, years))
    out['bottom_t_c'] = float(rec['t'][:, -1].mean())
    out['surface_t_c'] = float(rec['surface_t'].mean())
    out['mixed_depth_m'] = [float(np.quantile(rec['mixed_depth'], q)) for q in (0.5, 0.9, 0.99)]
    out['area_km2'] = float(area[0] / 1e6)
    out['mean_depth_m'] = float(bsn.volume.sum() / area[0])
    out['schmidt_stability_j_m2'] = schmidt_stability(rec, bsn, st.gravity)
    # the fill that would bring the floor up to 12 m everywhere, per m2 of the lake: the water below that plane
    below = np.clip((bsn.z[1:] - FILL_TO_M) / bsn.dz, 0.0, 1.0)
    out['fill_to_12m_m'] = float((bsn.volume * below).sum() / area[0])
    return out


def anoxic_floor(rec, bsn, m, years):
    """Nurnberg's anoxic factor: days a year that a floor area equal to the lake's surface lies under water below
    ANOXIC_FLOOR, from the daily columns."""
    ox = rec['oxygen'][:, :, m].astype(float)
    anoxic = (ox < ANOXIC_FLOOR) @ bsn.sediment / bsn.area[0]          # share of the lake's area, each day
    return dict(days_per_year=float(anoxic.mean() * 365.25), floor_share=float(anoxic.mean()))


def schmidt_stability(rec, bsn, gravity):
    """The column's Schmidt stability (J/m2), its mean through the record: the work to mix it fully."""
    rho = L.density(rec['t'].astype(float))
    v = bsn.volume
    rbar = (rho * v).sum(1, keepdims=True) / v.sum()
    zg = (bsn.centre * v).sum() / v.sum()
    work = gravity * ((bsn.centre - zg)[None, :] * (rho - rbar) * v[None, :]).sum(1) / bsn.area[0]
    return float(work.mean())


def deepest(st, name, case):
    rec = st.load(name)
    bsn = L.Basin(rec['z'], rec['area'])
    dt = float(rec['dt'])
    years = rec['step'] * dt / (365.25 * 86400.0)
    deep_mask = np.cumsum(bsn.volume) - 0.5 * bsn.volume >= 0.9 * bsn.volume.sum()
    w = np.where(deep_mask, bsn.volume, 0.0) / bsn.volume[deep_mask].sum()
    out = dict(target=case['target'], years=float(years[-1]), traits={})
    for trophic in TRAITS:
        m = st.names.index(f'{trophic}/{CENTRAL}')
        ox = rec['oxygen'][:, :, m].astype(float) @ w
        debt = rec['debt'][:, :, m].astype(float) @ w
        gas = column_gases(st, rec, bsn, m, METHANE_SHARE['central'])
        # the oxygen lost and the debt gained over the run, net of what the lake's overturns bring down; the lake
        # turns over, so its deep water gathers no gas toward a bubble point, and the years it would take isolated
        # follow the hazard map's rule: the trophic state's supply-fixed floor and deep-water demand over a 50 m
        # hypolimnion at the lake's floor
        span = float(years[-1] - years[0])
        net = ((ox[0] - ox[-1]) + (debt[-1] - debt[0])) / span
        t_deep = float(rec['t'][-1][deep_mask].mean())
        p = R.TROPHIC[trophic]
        rate = (p['sediment'] + p['deep_water']) / 50.0 * 365.25                     # g O2 m-3 a-1
        out['traits'][trophic] = dict(deepest_tenth_o2_start=float(ox[0]), deepest_tenth_o2_end=float(ox[-1]),
                                      net_loss_over_the_overturns_g_m3_yr=float(net), debt_end_g_m3=float(debt[-1]),
                                      saturation_ratio_max=float(gas['ratio'][:, bsn.centre >= HYPOLIMNION_M].max()),
                                      bottom_pressure_atm=float(gas['p_local'][-1]),
                                      isolated_demand_g_m3_yr=rate,
                                      years_to_bubble_point_if_isolated=isolated_years(st, bsn.z[-1], t_deep, rate))
    out['mixing_reached_bottom'] = bool((rec['mixed_depth'] >= bsn.z[-2] - 1e-6).any())
    out['overturn_events'] = len(overturns(rec, bsn, dt, None))
    return out


def isolated_years(st, depth_m, t_c, demand_g_m3_yr, gravity=None, surface_atm=SURFACE_ATM, share=None,
                   n2_atm=None, earth=False):
    """Years for water isolated at depth_m to reach the bubble point when it loses oxygen at demand_g_m3_yr and its
    respiration then goes on anaerobically, methane taking share of it and the rest CO2 (with the lake's alkalinity).
    earth takes Earth's dry air at surface_atm (its oxygen, nitrogen, argon and CO2) for the comparison."""
    gravity = gravity or st.gravity
    share = METHANE_SHARE['central'] if share is None else share
    dry = surface_atm - float(L.vapour_pressure_atm(t_c))
    if earth:
        sat = float(L.oxygen_saturation(t_c, surface_atm, st.salinity, EARTH_AIR['o2']))
        n2 = float(G.solubility('N2', t_c)) * EARTH_AIR['n2'] * dry
        ar = float(G.solubility('Ar', t_c)) * EARTH_AIR['ar'] * dry
        base = float(equilibrium_dic(st, np.array(t_c), pco2=EARTH_AIR['co2_ppm'] * 1e-6 * surface_atm))
    else:
        sat = float(L.oxygen_saturation(t_c, surface_atm, st.salinity, st.o2_fraction))
        n2 = float(G.solubility('N2', t_c)) * (st.n2_frac if n2_atm is None else n2_atm / dry) * dry
        ar = float(G.solubility('Ar', t_c)) * st.ar_frac * dry
        base = float(equilibrium_dic(st, np.array(t_c)))
    p_local = float(G.local_pressure_atm(depth_m, surface_atm, gravity))

    def ratio(debt):
        ch4 = share * debt / 64.0
        dic = sat / 32.0 + debt * (share / 64.0 + (1.0 - share) / 32.0)
        _, co2 = G.speciate(base + dic / 1000.0, st.ta, t_c, LAKE_IONIC)
        return float(G.gas_pressure_atm(t_c, ch4=ch4, co2_aq=co2, n2=n2, ar=ar)) / p_local

    if ratio(0.0) >= 1.0:
        return 0.0
    lo, hi = 0.0, 1.0
    while ratio(hi) < 1.0:
        hi *= 2.0
        if hi > 1e7:
            return None
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if ratio(mid) < 1.0 else (lo, mid)
    return (sat + hi) / demand_g_m3_yr


def hazard_map(st, temps):
    """Years for an isolated hypolimnion to reach the bubble point at its floor, on the Moon and on Earth, by floor
    depth and trophic state: its demand the trophic state's floor and deep-water demand over a hypolimnion reaching
    from the floor to 15 m below the surface, or 50 m at most."""
    out = {}
    for depth in (20.0, 30.0, 50.0, 100.0, 200.0, 352.0, 1000.0, 1925.0):
        thick = min(max(depth - 15.0, 5.0), 50.0)
        row = {}
        for trophic in TRAITS:
            p = R.TROPHIC[trophic]
            rate = (p['sediment'] + p['deep_water']) / thick * 365.25            # g O2 m-3 a-1
            row[trophic] = dict(
                moon={k: isolated_years(st, depth, t, rate) for k, t in temps.items()},
                earth_at_25c=isolated_years(st, depth, 25.0, rate, gravity=STANDARD_GRAVITY, surface_atm=1.0,
                                            earth=True))
        out[f'{depth:g}'] = dict(hypolimnion_m=thick, traits=row)
    return out


# --- ebullition and the Moon's lakes ----------------------------------------------------------------------------

def survival_table(st, temps):
    """Share of a sediment bubble's methane reaching the air by release depth, bubble size and interface, at each
    target's deep temperature; the water around it holds N2 at saturation, no oxygen and no methane."""
    depths = (1.0, 3.0, 10.0, 30.0, 60.0, 100.0, 200.0, 352.0)
    out = {}
    for target, t_c in temps.items():
        n2, _ = air_gases(st, t_c)
        rows = {}
        for size, mm in BUBBLE_MM.items():
            for iface in ('rigid', 'clean'):
                rows[f'{size}/{iface}'] = [G.bubble_survival(z, mm / 1000.0, t_c, SURFACE_ATM, st.gravity,
                                                             dict(N2=float(n2)), BUBBLE_CH4_FRACTION, iface)[0]
                                           for z in depths]
        out[target] = rows
    n2e = float(G.solubility('N2', 20.0)) * 0.78 * (1.0 - float(L.vapour_pressure_atm(20.0)))
    earth = {f'5mm/{iface}': [G.bubble_survival(z, 0.005, 20.0, 1.0, STANDARD_GRAVITY, dict(N2=n2e), 0.8, iface)[0]
                              for z in depths] for iface in ('rigid', 'clean')}
    return dict(depths_m=list(depths), moon=out, earth_20c=earth)


def moon_lakes():
    """The rain-fed lakes' cells (local depth and area) and each body's deepest point, from geography's drainage
    product on LOLA, as lake_oxygen's lakes.py reads them."""
    from research.studies.lake_oxygen import lakes as K
    d = K.lake_depths()
    b = K.bodies(d)
    deepest = b['deepest'][b['body']]
    return dict(depth=np.maximum(d['depth'], 0.1), area=d['area'], deepest=deepest)


DEPTH_EDGES = np.array([0.0, 2.0, 5.0, 10.0, 20.0, 40.0, 70.0, 100.0, 150.0, 250.0, 400.0, 700.0, 1e5])
CLASS_EDGES = (('pans, under 10 m', 0.0, 10.0), ('10-100 m', 10.0, 100.0), ('100-1,000 m', 100.0, 1000.0),
               ('over 1,000 m', 1000.0, 1e6))


def lake_classes(cells):
    a = cells['area']
    total = a.sum()
    by_depth = np.histogram(cells['depth'], DEPTH_EDGES, weights=a)[0] / total
    classes = {}
    for name, lo, hi in CLASS_EDGES:
        m = (cells['deepest'] >= lo) & (cells['deepest'] < hi)
        classes[name] = dict(area_share=float(a[m].sum() / total),
                             by_local_depth=(np.histogram(cells['depth'][m], DEPTH_EDGES, weights=a[m])[0]
                                             / max(a[m].sum(), 1.0)).tolist())
    giant = cells['deepest'] > GIANT_M
    return dict(area_km2=float(total / 1e6), depth_edges_m=DEPTH_EDGES.tolist(), by_local_depth=by_depth.tolist(),
                classes=classes, giant_area_share=float(a[giant].sum() / total))


def sediment_ebullition(st, trophic, depth, t_c, share, size, iface, warm=False):
    """Methane (mol m-2 of floor a year) that a floor at depth makes, sends up as bubbles, and gets to the air."""
    p = R.TROPHIC[trophic]
    # a floor within 25 m of the euphotic zone's base takes the deep water's demand too, as the lake model has it
    rate = p['sediment'] + (p['deep_water'] if depth < 4.6 / R.kd(p['secchi_m']) + 25.0 else 0.0)
    if warm:                                    # lake_oxygen's warm bracket, from the source hypolimnia's 6 C
        rate *= math.exp(R.KINETICS['warm'] / L.BOLTZMANN_EV * (1.0 / (R.SOURCE_C + 273.15) - 1.0 / (t_c + 273.15)))
    production = share * rate / 64.0 / 86400.0                       # mol m-2 s-1
    n2, _ = air_gases(st, t_c)
    dry_n2 = float(n2 / G.solubility('N2', t_c))
    f_e = G.ebullition_share(production, depth, t_c, SURFACE_ATM, st.gravity, dry_n2,
                             decay_per_m=SEDIMENT_DECAY_PER_M)
    surv = G.bubble_survival(depth, BUBBLE_MM[size] / 1000.0, t_c, SURFACE_ATM, st.gravity, dict(N2=float(n2)),
                             BUBBLE_CH4_FRACTION, iface)[0] if f_e > 0 else 0.0
    year = 365.25 * 86400.0
    return dict(made=production * year, bubbles=production * f_e * year, to_air=production * f_e * surv * year,
                share=f_e, survival=surv)


def budget(st, deep, pans, classes, temps):
    """The lakes' methane Moon-wide by target and trophic state: bubbles from every floor by its local depth, and the
    pelagic flux the anoxic water sends up: what each overturn releases, by each of the two release rules (the net
    rise and the stock), oxidised at the overturn's share OXIDISED gives as far as the overturned lake's oxygen
    reaches, and the rest crossing a lasting oxycline at the share OXYCLINE_ESCAPE gives. The giants' bracket gives
    the lakes deeper than GIANT_M the oxycline's flux alone, as if they never overturned."""
    mids = 0.5 * (DEPTH_EDGES[:-1] + DEPTH_EDGES[1:])
    mids[-1] = 1000.0
    area_m2 = classes['area_km2'] * 1e6
    out = {}
    for target in TARGET_ORDER:
        runs = [v for v in deep.values() if v['target'] == target]
        surf_t = float(np.mean([v['surface_t_c'] for v in runs]))
        deep_t = float(np.mean([v['bottom_t_c'] for v in runs]))
        entry = {}
        for trophic in TRAITS:
            ebul = {}
            # low: small bubbles with clean surfaces, which lose the most on the way up; high: large contaminated
            # ones, which keep the most
            for case, (size, iface) in dict(central=('central', 'rigid'), low=('small', 'clean'),
                                            high=('large', 'rigid')).items():
                rows = [sediment_ebullition(st, trophic, z, surf_t if z < 15.0 else deep_t,
                                            METHANE_SHARE['central'], size, iface) for z in mids]
                ebul[case] = float(sum(r['to_air'] * s for r, s in zip(rows, classes['by_local_depth'])))
            made = float(np.mean([v['traits'][trophic]['methane_made_mol_m2_yr']['central'] for v in runs]))
            made_hi = float(np.mean([v['traits'][trophic]['methane_made_mol_m2_yr']['high'] for v in runs]))
            made_lo = float(np.mean([v['traits'][trophic]['methane_made_mol_m2_yr']['low'] for v in runs]))
            oxycline_only = {k: (made_lo if k == 'low' else made_hi if k == 'high' else made) * OXYCLINE_ESCAPE[k]
                             for k in ('low', 'central', 'high')}
            giant_share = classes['giant_area_share']
            rules = {}
            for rule in RELEASE_RULES:
                # what each overturn releases meets only the oxygen of the overturned lake: what that cannot oxidise
                # leaves or is oxidised as the air brings oxygen in, taken as half each; the rest of the methane made
                # crosses a lasting oxycline
                released, met = [], []
                for v in runs:
                    tr = v['traits'][trophic]
                    released.append(sum(tr['released_at_overturns_mol_m2'][rule]) / v['years'])
                    need = sum(tr['oxygen_to_oxidise_stored_methane_g_m2'][rule])
                    have = sum(tr['oxygen_after_overturns_g_m2'])
                    met.append(min(1.0, have / need) if need > 0 else 1.0)
                released_m, met_m = float(np.mean(released)), float(np.mean(met))
                through = max(made - released_m, 0.0)

                def pelagic_flux(scale, oxidised, oxycline):
                    over = met_m * (1.0 - oxidised) + (1.0 - met_m) * OVERTURN_UNMET_ESCAPE
                    return scale * (released_m * over + through * oxycline)

                pelagic = dict(central=pelagic_flux(1.0, OXIDISED['central'], OXYCLINE_ESCAPE['central']),
                               low=pelagic_flux(made_lo / made if made else 0.0, OXIDISED['high'],
                                                OXYCLINE_ESCAPE['low']),
                               high=pelagic_flux(made_hi / made if made else 0.0, OXIDISED['low'],
                                                 OXYCLINE_ESCAPE['high']))
                per_m2 = {k: ebul[k] + pelagic[k] for k in ('low', 'central', 'high')}
                per_m2_giants = {k: ebul[k] + (1.0 - giant_share) * pelagic[k] + giant_share * oxycline_only[k]
                                 for k in ('low', 'central', 'high')}
                rules[rule] = dict(
                    released_at_overturns_mol_m2_yr=released_m, overturn_oxygen_met_share=met_m,
                    pelagic_mol_m2_yr=pelagic, bubble_share_of_the_flux=ebul['central'] / per_m2['central']
                    if per_m2['central'] > 0 else None,
                    flux_g_m2_yr={k: x * 16.04 for k, x in per_m2.items()},
                    moon_tg_yr={k: x * 16.04 * area_m2 / TG for k, x in per_m2.items()},
                    moon_tg_yr_giants_on_the_oxycline={k: x * 16.04 * area_m2 / TG for k, x in per_m2_giants.items()})
            entry[trophic] = dict(
                ebullition_mol_m2_yr=ebul, oxycline_only_mol_m2_yr=oxycline_only,
                anoxic_methane_made_mol_m2_yr=dict(low=made_lo, central=made, high=made_hi), by_release_rule=rules)
        out[target] = dict(surface_t_c=surf_t, deep_t_c=deep_t, traits=entry)
    return out


def forcing_at(st, ch4_ppm):
    """All-sky forcing (W/m2) of methane with N2O at the held 0.1 ppm: the deep column's paired grid in ln(ppm), and
    beyond it the methane-alone grid's increments on the grid's top."""
    both = st.forcing['both']
    xs = sorted(float(k.split('+')[0]) for k in both if k.endswith('+0.1'))
    ys = [both[f'{x:g}+0.1']['allsky_w_m2'] for x in xs]
    alone = st.forcing['CH4']
    ax = sorted(float(k) for k in alone)
    ay = [alone[f'{x:g}']['allsky_w_m2'] for x in ax]
    if ch4_ppm <= xs[-1]:
        return float(np.interp(math.log(max(ch4_ppm, xs[0])), np.log(xs), ys))
    return float(ys[-1] + np.interp(math.log(ch4_ppm), np.log(ax), ay) - np.interp(math.log(xs[-1]), np.log(ax), ay))


def air(st, budget_rows):
    """What the lakes' flux does to the air's methane: the domain synthesis's ledger with its lakes replaced, the
    sources cut by biology as the held pair assumes, held at 1.87 ppm by removal or left to the soils."""
    m = st.ledgers['gases']['methane']
    held = st.ledgers['gases']['held']
    tg_per_ppm = m['tg_per_ppm']
    lakes_ledger = 10.0 * 0.1034 * 4 * math.pi * 1737.4e3 ** 2 / TG       # the screen's central 10 g m-2 a-1
    sink = m['central'][0]['soil_sink_tg_yr_per_ppm']
    base_low, base_high = (c['sources_tg_yr']['total'] for c in m['central'])
    out = dict(ledger_lakes_tg_yr=lakes_ledger, ledger_range_tg_yr=m['wild_sources_tg_yr']['lakes'],
               soil_sink_tg_yr_per_ppm=sink, tg_per_ppm=tg_per_ppm, held_ppm=held['methane_ppm'],
               held_forcing_w_m2=held['forcing_w_m2']['allsky'], targets={})
    f_held = forcing_at(st, held['methane_ppm'])
    for target in TARGET_ORDER:
        rows = {}
        for trophic in TRAITS:
            lake = budget_rows[target]['traits'][trophic]['moon_tg_yr']
            row = {}
            for case in ('low', 'central', 'high'):
                sources = (base_low - lakes_ledger + lake[case], base_high - lakes_ledger + lake[case])
                after_cut = tuple(0.1 * s for s in sources)
                steady = tuple(s / sink for s in after_cut)
                removal = tuple(max(s - sink * held['methane_ppm'], 0.0) for s in after_cut)
                df = tuple(forcing_at(st, x) - f_held for x in steady)
                uncut = tuple(0.1 * (s - lake[case]) + lake[case] for s in sources)
                uncut_ppm = tuple(x / sink for x in uncut)
                uncut_df = tuple(forcing_at(st, x) - f_held for x in uncut_ppm)
                row[case] = dict(lakes_tg_yr=lake[case], sources_tg_yr=sources, left_after_90pct_cut_tg_yr=after_cut,
                                 lakes_uncut=dict(sources_tg_yr=uncut, steady_ppm_by_soils=uncut_ppm,
                                                  removal_to_hold_tg_yr=tuple(max(x - sink * held['methane_ppm'], 0.0)
                                                                              for x in uncut),
                                                  air_passed_per_year_share=tuple(
                                                      max(x - sink * held['methane_ppm'], 0.0)
                                                      / (held['methane_ppm'] * tg_per_ppm) for x in uncut),
                                                  forcing_beyond_held_w_m2=uncut_df,
                                                  added_dimming_percent=tuple(100.0 * x / st.absorbed[target]
                                                                              for x in uncut_df)),
                                 steady_ppm_by_soils=steady, removal_to_hold_tg_yr=removal,
                                 air_passed_per_year_share=tuple(r / (held['methane_ppm'] * tg_per_ppm)
                                                                 for r in removal),
                                 forcing_beyond_held_w_m2=df,
                                 added_dimming_percent=tuple(100.0 * x / st.absorbed[target] for x in df),
                                 lakes_alone_uncut_ppm=lake[case] / sink)
            rows[trophic] = row
        out['targets'][target] = rows
    return out


# --- phosphorus -------------------------------------------------------------------------------------------------

def phosphorus(st, deep):
    """Internal loading from anoxic floors (Nurnberg's release rate times the anoxic factor) and its feedback on the
    lakes' phosphorus: Vollenweider's relation for the external load, TP = L / (q_s (1 + f tau^0.5)), with its settling
    at Earth's gravity (f = 1) and carried to the Moon's (f = g/6 or its square root), and the internal load added
    outside the retention, TP_int = L_int / q_s, as Nurnberg's (1984) mass balance has it for lakes drawing on a
    floor's stock. At steady state the phosphorus a floor releases is what settled on it, so whatever the floors do
    the lake keeps at most what flows in: TP <= L / q_s, (1 + f tau^0.5) times the oxic lake's, the retention at most
    cancelled; burial keeps part of it."""
    mean_depth = st.lake_volume_km3 * 1e9 / (st.lake_area_km2 * 1e6)
    q_s = mean_depth / st.residence_yr
    area = st.lake_area_km2 * 1e6
    dissolved = [x / area * 1000.0 for x in st.p_loads['without_lakes']]               # g P m-2 a-1 (kg -> g)
    particulate = [x / area * 1000.0 for x in st.p_loads['particulate_without_lakes']]
    tp = lambda load_g, qs, tau, f: 1000.0 * load_g / (qs * (1.0 + f * math.sqrt(tau)))   # mg/m3, external load
    out = dict(mean_depth_m=mean_depth, hydraulic_load_m_yr=q_s, residence_yr=st.residence_yr, settling=SETTLING,
               loads_g_m2_yr=dict(dissolved=dissolved, particulate=particulate), targets={})
    deep_classes = dict(mean_depth_m=float(np.mean([v['mean_depth_m'] for v in deep.values()])))
    tau_352 = deep_classes['mean_depth_m'] / q_s
    out['apparent_settling_m_yr'] = {k: f * q_s * math.sqrt(tau_352) for k, f in SETTLING.items()}
    out['tp_from_dissolved_load_mg_m3'] = {k: dict(all_lakes=[tp(x, q_s, st.residence_yr, f) for x in dissolved],
                                                   lakes_352m=[tp(x, q_s, tau_352, f) for x in dissolved])
                                           for k, f in SETTLING.items()}
    out['tp_from_all_load_mg_m3'] = {k: dict(all_lakes=[tp(x + y, q_s, st.residence_yr, f)
                                                        for x, y in zip(dissolved, particulate)],
                                             lakes_352m=[tp(x + y, q_s, tau_352, f)
                                                         for x, y in zip(dissolved, particulate)])
                                     for k, f in SETTLING.items()}
    out['permissible_load_g_m2_yr'] = {k: {c: v * q_s * (1.0 + f * math.sqrt(tau_352)) / 1000.0
                                           for c, v in OECD_TP.items()} for k, f in SETTLING.items()}
    # at steady state anoxic floors can at most cancel the retention: the lake keeps no more than what flows in
    out['without_retention'] = dict(
        tp_from_dissolved_load_mg_m3=[1000.0 * x / q_s for x in dissolved],
        tp_from_all_load_mg_m3=[1000.0 * (x + y) / q_s for x, y in zip(dissolved, particulate)],
        permissible_load_g_m2_yr={c: v * q_s / 1000.0 for c, v in OECD_TP.items()},
        over_the_oxic_lake={k: 1.0 + f * math.sqrt(tau_352) for k, f in SETTLING.items()})
    for target in TARGET_ORDER:
        runs = [v for v in deep.values() if v['target'] == target]
        rows = {}
        for trophic in TRAITS:
            af = float(np.mean([v['traits'][trophic]['anoxic_floor']['days_per_year'] for v in runs]))
            internal = [P_RELEASE[trophic] * af * f / 1000.0 for f in LOW_SULPHATE_FACTOR]     # g P m-2 a-1
            with_internal = {k: [tp(dissolved[0], q_s, tau_352, f) + 1000.0 * i / q_s for i in internal]
                             + [tp(dissolved[1], q_s, tau_352, f) + 1000.0 * internal[1] / q_s]
                             for k, f in SETTLING.items()}
            ceiling = [1000.0 * dissolved[0] / q_s] * len(internal) + [1000.0 * dissolved[1] / q_s]
            row = dict(anoxic_factor_days_yr=af, internal_load_g_m2_yr=internal,
                       internal_over_dissolved_load=[i / dissolved[1] for i in internal] + [internal[1] / dissolved[0]],
                       tp_with_internal_mg_m3=with_internal,
                       tp_with_internal_at_steady_state_mg_m3={k: [min(a, c) for a, c in zip(v, ceiling)]
                                                               for k, v in with_internal.items()})
            if trophic == 'new_reservoir':
                row['flooded_soils_internal_load_g_m2_yr'] = [P_RELEASE_FLOODED_SOILS * af * f / 1000.0
                                                              for f in LOW_SULPHATE_FACTOR]
            rows[trophic] = row
        out['targets'][target] = rows
    return out


# --- remedies ---------------------------------------------------------------------------------------------------

def remedies(st, deep, temps):
    area = st.lake_area_km2
    out = dict(lake_area_km2=area, targets={})
    q_s = st.lake_volume_km3 * 1e9 / (st.lake_area_km2 * 1e6) / st.residence_yr     # m/yr of throughflow
    for target in TARGET_ORDER:
        runs = [v for v in deep.values() if v['target'] == target]
        t_deep = float(np.mean([v['bottom_t_c'] for v in runs]))
        sat = float(L.oxygen_saturation(t_deep, SURFACE_ATM, st.salinity, st.o2_fraction))
        # respiration's warm kinetics from the source hypolimnia's temperature (lake_oxygen's warm columns): a
        # transient no settling supply sustains, given beside the bounds
        warm = math.exp(R.KINETICS['warm'] / L.BOLTZMANN_EV * (1.0 / (R.SOURCE_C + 273.15) - 1.0 / (t_deep + 273.15)))
        rows = {}
        for trophic in TRAITS:
            p = R.TROPHIC[trophic]
            demand = p['sediment'] + p['deep_water']                       # g O2 m-2 d-1, fixed by supply
            o2 = (demand * INDUCED_DEMAND[0], demand * INDUCED_DEMAND[1])   # t O2 km-2 d-1 (g/m2 = t/km2): the demand,
            delivery = []                                                   # and as the devices' stirring raises it
            for depth in (352.0, 1925.0):
                ratio = float(G.local_pressure_atm(depth, SURFACE_ATM, st.gravity)) / SURFACE_ATM
                kwh = 1000.0 / 32.0 * G.GAS_R * (t_deep + 273.15) * math.log(ratio) / COMPRESSOR_EFFICIENCY / 3.6e6
                delivery.append(kwh)
            kw = [o2[0] * 1000.0 * (O2_KWH_PER_KG[0] + delivery[0]) / 24.0,
                  o2[1] * 1000.0 * (O2_KWH_PER_KG[1] + delivery[1]) / 24.0]
            transient = demand * warm * INDUCED_DEMAND[1]
            rows[trophic] = dict(
                oxygenation=dict(o2_t_km2_day=o2, kwh_per_kg=[O2_KWH_PER_KG[0] + delivery[0],
                                                              O2_KWH_PER_KG[1] + delivery[1]],
                                 kw_per_km2=kw, gw_all_lakes=[k * area / 1e6 for k in kw]),
                warm_transient=dict(o2_t_km2_day=transient,
                                    kw_per_km2=transient * 1000.0 * (O2_KWH_PER_KG[1] + delivery[1]) / 24.0,
                                    note='the stirred demand at warm kinetics from 6 C, which no settling supply '
                                         'sustains: the weeks to months after a warming, beside the bounds'),
                inflow_routing=dict(o2_g_m2_day=q_s * sat / 365.25, share_of_demand=q_s * sat / 365.25 / demand))
        circulation = {}
        for depth in (352.0, 1925.0):
            ratio = float(G.local_pressure_atm(depth, SURFACE_ATM, st.gravity)) / SURFACE_ATM
            q = CIRCULATION_AIR_M3_MIN_KM2 / 60.0
            kw_moon = SURFACE_ATM * G.ATM * q * math.log(ratio) / COMPRESSOR_EFFICIENCY / 1000.0
            ratio_e = (1.0 + 1000.0 * STANDARD_GRAVITY * depth / G.ATM)
            kw_earth = G.ATM * q * math.log(ratio_e) / COMPRESSOR_EFFICIENCY / 1000.0
            circulation[f'{depth:g}'] = dict(kw_per_km2=kw_moon, earth_kw_per_km2=kw_earth,
                                             gw_all_lakes=kw_moon * area / 1e6,
                                             basis='Lorenzen & Fast\'s per-area air flow, Earth\'s guideline for lakes '
                                                   'tens of metres deep, compressed to the depth; no plume is computed')
        stab = float(np.mean([v['schmidt_stability_j_m2'] for v in runs]))
        out['targets'][target] = dict(deep_t_c=t_deep, saturation_mg_l=sat, warm_factor=warm, traits=rows,
                                      circulation=circulation, schmidt_stability_j_m2=stab,
                                      mixed_depth_m=[float(np.mean([v['mixed_depth_m'][i] for v in runs]))
                                                     for i in range(3)])
    out['earth_installations_t_o2_km2_day'] = EARTH_INSTALLATIONS_T_DAY_KM2
    return out


def degassing(st, deep):
    """Pipes as at Nyos: the water each km2 must lift a year to carry away its methane, whether the lifted water
    would bubble on its own, and the methane's energy against the pumping."""
    out = {}
    for target in TARGET_ORDER:
        runs = [v for v in deep.values() if v['target'] == target]
        rows = {}
        for trophic in TRAITS:
            made = float(np.mean([v['traits'][trophic]['methane_made_mol_m2_yr']['central'] for v in runs]))
            ch4 = float(np.mean([v['traits'][trophic]['bottom_ch4_max_umol_l'] for v in runs])) / 1000.0  # mol/m3
            p_gas = float(np.mean([v['traits'][trophic]['gas_pressure_bottom_max_atm'] for v in runs]))
            flow = made * 1e6 / max(0.9 * ch4, 1e-9) / (365.25 * 86400.0)           # m3/s per km2
            pump_kw = flow * 5e4 / 0.7 / 1000.0                                    # 50 kPa through pipe and degasser
            rows[trophic] = dict(methane_mol_m2_yr=made, deep_ch4_mol_m3=ch4, water_m3_s_per_km2=flow,
                                 self_lifting=bool(p_gas > SURFACE_ATM), deep_gas_pressure_atm=p_gas,
                                 pumping_kw_per_km2=pump_kw,
                                 methane_energy_kw_per_km2=made * 1e6 * CH4_ENERGY_KJ_MOL / (365.25 * 86400.0))
        out[target] = rows
    return out


def analyse(st, cases, classes):
    """Everything the study reads from one wind share's saved columns."""
    deep = {name: deep_lake(st, name, c) for name, c in cases.items() if c['kind'] == 'deep'}
    deepest_runs = {name: deepest(st, name, c) for name, c in cases.items() if c['kind'] == 'deepest'}
    temps = {t: float(np.mean([v['bottom_t_c'] for v in deep.values() if v['target'] == t])) for t in TARGET_ORDER}
    b = budget(st, deep, None, classes, temps)
    by_rule = {rule: {t: dict(traits={tr: v['by_release_rule'][rule] for tr, v in b[t]['traits'].items()})
                      for t in b} for rule in RELEASE_RULES}
    return dict(deep_lakes=deep, deepest=deepest_runs, hazard_map=hazard_map(st, temps),
                survival=survival_table(st, temps), methane=b,
                air={rule: air(st, rows) for rule, rows in by_rule.items()}, phosphorus=phosphorus(st, deep),
                remedies=remedies(st, deep, temps), degassing=degassing(st, deep))


def stale_columns(st):
    """The saved columns whose starting run the waters' oxygen study's cache does not hold under the key its current
    code and setting give it: profiles.py started each 352 m lake from that study's cached deep run, whose file it
    records, and each 1,925 m lake from the 0-15 lake's; empty when the cache holds every one with the same contents."""
    lo = R.Setting()
    out = []
    for name, case in st.meta['cases'].items():
        extra = {} if case['eta'] == 'central' else dict(eta=case['eta'])
        run = dict(kind='deep', target=case['target'], site_i=R.SITE_BANDS.index(case['band']), years=None, **extra)
        now = R.RUNS / 'cases' / f"deep_{R.cache_key(lo, 'deep', **run)}.json"
        saved = R.RUNS / 'cases' / case['from_case'] if 'from_case' in case else None
        if not now.exists() or (saved is not None and sha256(saved) != sha256(now)):
            out.append(name)
    return out


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--cache-only', action='store_true',
                    help="stop unless the saved columns start from lake_oxygen's cached runs (it computes none)")
    args = ap.parse_args(argv)
    tic = time.time()
    st = Setting()
    if args.cache_only or R.cache_only():
        stale = stale_columns(st)
        if stale:
            raise SystemExit(f"cache only: lake_oxygen's cache does not hold the runs {len(stale)} saved columns "
                             f'start from ({", ".join(stale[:4])})')
    cells = moon_lakes()
    classes = lake_classes(cells)
    by = {e: analyse(st, {n: c for n, c in st.meta['cases'].items() if c.get('eta', 'central') == e}, classes)
          for e in ETAS}
    out = dict(
        schema=SCHEMA,
        producer=dict(files={str(p.resolve().relative_to(ROOT)): sha256(p)
                             for p in (Path(__file__), HERE / 'gases.py', HERE / 'profiles.py', Path(L.__file__))},
                      inputs=st.inputs,
                      constants=dict(MOON_SURFACE_GRAVITY=MOON_SURFACE_GRAVITY, STANDARD_GRAVITY=STANDARD_GRAVITY)),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        setting=dict(methane_share=METHANE_SHARE, oxidised=OXIDISED, oxycline_escape=OXYCLINE_ESCAPE,
                     overturn_unmet_escape=OVERTURN_UNMET_ESCAPE, methane_made_below_g_m3=ANOXIC,
                     lake_alkalinity_meq_kg=st.ta * 1000.0, lake_ph=st.lake_ph, bubble_mm=BUBBLE_MM,
                     bubble_ch4_fraction=BUBBLE_CH4_FRACTION, sediment_decay_per_m=SEDIMENT_DECAY_PER_M,
                     p_release_mg_m2_d=P_RELEASE, p_release_flooded_soils_mg_m2_d=P_RELEASE_FLOODED_SOILS,
                     settling=SETTLING, low_sulphate_factor=LOW_SULPHATE_FACTOR, anoxic_floor_g_m3=ANOXIC_FLOOR,
                     o2_kwh_per_kg=O2_KWH_PER_KG, induced_demand=INDUCED_DEMAND, warm_reference_c=R.SOURCE_C,
                     circulation_air_m3_min_km2=CIRCULATION_AIR_M3_MIN_KM2, absorbed_sunlight_w_m2=st.absorbed,
                     surface_atm=SURFACE_ATM, gravity=st.gravity, giant_m=GIANT_M, earth_air=EARTH_AIR,
                     deep_mixing=R.DEEP_MIXING, wind_share={e: R.WIND_SHARE[e] for e in ETAS}),
        moon_lakes=classes, **by['central'], eta_low=by['low'], eta_high=by['high'],
        seconds=round(time.time() - tic, 1))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(dumps(rnd(out)) + '\n')
    print('wrote', OUT.relative_to(ROOT), f'in {time.time() - tic:.0f} s')
    return 0


EVIDENCE = (
    'Arithmetic and small reference calculations on one lake model\'s columns: the waters\' oxygen study\'s deep '
    'lakes (biosphere/ecology/lake.py as corrected after the review of 18 October and its re-check) rerun from their '
    'own settled columns with the column saved daily, with the wind\'s deep mixing at its central share of P10 and at '
    'its bracket. Methane is the model\'s debt of reduced substances in anoxic water times a methane share from '
    'Earth\'s low-sulphate lakes; what each overturn releases is bounded by two rules, the inventory\'s net rise '
    'since the last and the stock it brings up within the methane made since the last, until a methane tracer '
    'settles it; CO2 '
    'is what the water respired since it left the surface, in the lake\'s carbonate system; the gas saturation ratio '
    'sets them, with N2 and Ar at the surface\'s saturation, against the local pressure at lunar gravity. Sediment '
    'ebullition is the steady porewater balance of Bazhin as Langenegger et al. (2019) apply it, checked against Lake '
    'Soppen; bubbles rise with a contaminated-bubble drag and Froessling or Higbie transfer at lunar gravity. '
    'Oxidation and escape shares, phosphorus release rates (Nurnberg\'s class means, inputs), loads, Vollenweider\'s '
    'retention with its settling bracketed for lunar gravity, and the remedies\' energy are Earth\'s published '
    'values; the oxygenation\'s bounds take the demand fixed by supply, the warm kinetics a transient beside them; '
    'the circulation\'s air flow is Lorenzen & Fast\'s guideline from Earth, with no plume computed. No lunar lake '
    'has been measured, and the model carries the methane and CO2 of its debt only in the density, as their net '
    'effect, and no phosphorus of its own.')
READING_RULE = (
    'Concentrations in mol/m3 (mmol/L) or umol/L as named; fluxes in mol or g per m2 of lake surface a year unless '
    'named per m2 of floor; Tg of CH4 a year Moon-wide over geography\'s 3.92 million km2 of rain-fed lakes. The '
    'top level holds the runs with the wind\'s deep mixing at its central share of P10 (0.3%), and eta_low and '
    'eta_high the same at its bracket (0.2% and 0.4%). deep_lakes are keyed deep_<target>_<band> (with _eta_low or '
    '_eta_high), each trophic state at the central '
    'renewal/supply column; methane_made is the anoxic mineralisation\'s methane at the low, central and high '
    'methane share (methane_made_warm_transient at the warm kinetics\' transient); released_at_overturns gives what '
    'each overturn releases by each rule: net_rise, the inventory\'s rise since the previous overturn (a lower '
    'bound), and stock, the stock the overturn brings up within the methane made since the previous overturn (an '
    'upper bound). saturation_ratio is the dissolved gases\' pressure over the local pressure: 1 is the bubble '
    'point. hazard_map gives years for an isolated hypolimnion to reach it at its floor, earth_at_25c under Earth\'s '
    'own air; deepest holds the 1,925 m lake\'s fifty years from an oxygenated start, which turns over '
    '(overturn_events) and so gathers no gas, with the years it would take isolated by the hazard map\'s rule. '
    'survival lists by depths_m the share of a bubble\'s methane reaching the air. '
    'methane.<target>.traits gives each trophic state\'s bubbles and, under by_release_rule, its flux by each rule '
    'with low/central/high brackets and moon_tg_yr_giants_on_the_oxycline with the lakes deeper than 1,000 m sending '
    'up only what crosses a lasting oxycline; air.<rule> gives the steady state and dimming with the lakes in the '
    'ledger, each pair the low- and high-cost ends of the domain synthesis. phosphorus is keyed by settling (earth, '
    'aggregates at the square root of g/6, stokes at g/6), loads in g P m-2 a-1; tp_with_internal_mg_m3 adds the '
    'internal load outside the retention (Nurnberg 1984), for floors drawing on a stock beyond their steady state, '
    'and tp_with_internal_at_steady_state_mg_m3 holds it to phosphorus.without_retention, what flows in over the '
    'hydraulic load, the most a lake keeps at steady state whatever its floors release; burial keeps it below that. '
    'remedies per km2 of lake: oxygenation\'s bounds at the demand and at three times it, as the devices\' stirring '
    'raises it, and warm_transient beside them.')


if __name__ == '__main__':
    sys.exit(main())
