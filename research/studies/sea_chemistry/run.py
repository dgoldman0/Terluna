"""The seas' make-up through time: the water's delivered salt and phosphorus, weathering, and the states it leads to.

    python -m research.studies.sea_chemistry.run      # results/sea_chemistry.json, a few seconds

Resources' question 9 asks how the seas' composition changes through time, with the delivered water's salts a design
lever. On 13 October 2026 the author asked for the model with fresh, salted and phosphorus-added water, and asked
that it look past the first thousands of years to the stable or dynamic state the seas might hold over the world's
billion-year life. This runner follows the common-level waters and the rain-fed lakes that feed them from the
delivery to a billion years:

- the delivered water: fresh, brackish (Earth's seawater diluted tenfold), salted (Earth's seawater) or fresh with
  Earth's mean phosphorus;
- weathering on the carbon cycle's terms (research/studies/carbon_cycle): Dessert's basalt law over the land's
  discharge, the developed soils ageing toward the catchment rate, which erosion keeps renewing, and the highland
  soil's calcium, magnesium and sodium with potassium, with the silica and sulphate of its export (lunar cycle
  ecology), of which the share the springs study finds reaching the waters arrives, with the springs' make-up; the
  lakes take 96.5% of it and pass their water to the seas;
- the floors: while the crust fills (stage 1) the lakes' and seas' floors draw off water with its dissolved load, and
  delivered fresh water replaces it;
- the carbonate system of each box in equilibrium with the air's 400 ppm, through the carbon cycle's constants and
  solver, and calcite where a water stands above its threshold (three times saturation in fresh water, the
  delivered saturation in a salted sea);
- magnesium's fate: kept in solution, taken into magnesium clays once a sea reaches pH 9, or steered into magnesite
  where the rivers enter the lakes, as the carbon cycle's managed route does;
- silica: held at amorphous silica's solubility without silica-building life, or at 0.1 mmol/kg with it;
- sediment: the seas kept open, the land's sediment returned to it and the waters' own deposits taken out, or the
  sediment and deposits left to settle, with the salts the pore water buries, until the lakes and then the seas are
  full;
- phosphorus: the sea phosphorus study's box, with the lakes' trap ending when they fill.

The air's CO2 is held at 400 ppm, as the planned world holds it; the CO2 the waters take and give back is reported.
"""
from __future__ import annotations
import hashlib
import json
import math
import sys
from pathlib import Path
from scipy.optimize import brentq

from research.studies.carbon_cycle import run as cc

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'sea_chemistry.json'
SCHEMA = 'terluna.research.sea-chemistry/1'
FILES = ('research/studies/sea_chemistry/run.py', 'research/studies/carbon_cycle/run.py')
SCREEN = 'resources/results/first_screen.json'
SEA_P = 'resources/results/sea_phosphorus.json'
CO2_STUDY = 'research/studies/atmospheric_co2/results/atmospheric_co2.json'
LUNAR_CYCLE = 'research/studies/lunar_cycle_ecology/results/lunar_cycle_ecology.json'
GATES = 'research/studies/conservation/results/gates.json'
ATLAS = 'geography/results/atlas.json'
DRAINAGE = 'geography/results/drainage.json'
CARBON = 'research/studies/carbon_cycle/results/carbon_cycle.json'
INPUTS = {SCREEN: 'terluna.resources.first-screen/1', SEA_P: 'terluna.resources.sea-phosphorus/1',
          CO2_STUDY: 'terluna.research.atmospheric-co2/1', LUNAR_CYCLE: 'terluna.research.lunar-cycle-ecology/1',
          GATES: None, ATLAS: 'terluna.geography.atlas/1', DRAINAGE: 'terluna.geography.drainage/1',
          CARBON: 'terluna.research.carbon-cycle/1', cc.SPRINGS: 'terluna.research.springs/1',
          cc.STAGE1: 'terluna.geography.groundwater-uptake/1'}
FLOOR_SPANS = ('10', '100', '1000', '10000', '100000', '1e+06', '1e+07')   # stage 1's spans, to ten million years

YEAR_S = cc.YEAR_S
# Earth's reference seawater at a salinity of 35.165 (Millero et al. 2008), mol/kg; alkalinity eq/kg.
SEAWATER = dict(na=0.46906, k=0.01021, mg=0.05282, ca=0.01028, cl=0.54586, so4=0.02824, alkalinity=0.0023)
SEAWATER_G_KG = 35.165
BRACKISH_SHARE = 0.1                     # the brackish delivery: Earth's seawater diluted tenfold
EARTH_P_MG_M3 = 70.0                     # Earth's ocean mean phosphorus (the sea phosphorus study)
M_P, M_CL, M_SO4, M_SIO2, M_HCO3 = 30.974, 35.453, 96.06, 60.084, 61.017
# Amorphous silica's solubility as H4SiO4 (Gunnarsson & Arnorsson 2000) and silicic acid's first dissociation.
SILICA_LOG_K_25C = -2.71
SILICA_PKA = 9.82
SILICA_LIFE_MOL_KG = 1e-4                # design guess: silica-building life holds dissolved silica near Earth's ocean's
PH_CLAY = 9.0                            # design guess: magnesium clays form from alkaline, silica-rich water above this
SI_PER_MG_CLAY = 1.5                     # sepiolite, Mg4Si6O15(OH)2.6H2O: one and a half silicons for each magnesium
# Sediment: dry bulk density of what settles (the sea phosphorus study's design guess) and the grains' own density.
DRY_DENSITY = (500.0, 1500.0)            # kg/m3
GRAIN_DENSITY = 2700.0
# What the waters deposit: the minerals' own densities (kg/m3) and the porosity of the deposits they build, a design
# guess between compacted carbonate and the oozes of silica-building life.
MINERAL_DENSITY = dict(calcite=2710.0, magnesite=3010.0, magnesium_clays=2100.0, opal=2050.0, soda=2200.0,
                       sodium_sulphate=2660.0)
DEPOSIT_POROSITY = (0.5, 0.8)
MOLAR_MASS = dict(calcite=100.087, magnesium_clays=0.25 * 647.8, opal=60.084)   # g per mole of Ca, Mg or Si
P_RESIDENCE_YR = (1000.0, 20000.0)       # the sea phosphorus study's range
P_REACTIVE_SHARE = (0.1, 0.3)            # Earth's reactive share of river particulate phosphorus, freed by salt water
FRESH_LIMIT_G_KG = 1.0
VALID_IONIC_STRENGTH = 0.5               # Davies activities hold to about this; beyond it only the masses are given
VALID_SALINITY_G_KG = (19.0, 43.0)       # the seawater constants' range (Lueker et al. 2000)
TRIAL_MARGIN = 2.0                       # a step's own inputs may carry a box this far past the range before its minerals
SEGMENTS = cc.SEGMENTS                   # the carbon cycle's steps to 10,000 years
LOG_STEPS_PER_DECADE = 40
END_YEARS = 1e9
REPORT_YEARS = (0, 10, 100, 1000, 3000, 10000, 30000, 100000, 300000, 1e6, 3e6, 1e7, 3e7, 1e8, 3e8, 1e9)
REPORT_SET = {round(float(y), 6) for y in REPORT_YEARS}
DELIVERIES = ('fresh', 'brackish', 'salted', 'fresh_phosphorus')
BRINE_CEILING_G_KG = 300.0              # design guess: soda and salt minerals deposit from brines of this order
SINKS = {
    'natural': dict(magnesium='clays', silica='abiotic', soda='kept'),
    'natural_without_clays': dict(magnesium='kept', silica='abiotic', soda='kept'),
    'managed': dict(magnesium='steered', silica='life', soda='works'),
}
SEDIMENT = ('kept_open', 'settling')
EVIDENCE = ('A box model on committed products: the rain-fed lakes and the common-level waters as well-mixed boxes in '
            'CO2 equilibrium with the air at 400 ppm, at the design case\'s temperatures, fed by the carbon cycle\'s '
            'weathering of the highland soil, with the carbon cycle\'s carbonate constants and solver (fresh water '
            'with Davies activities; Earth\'s seawater with its stoichiometric constants). Calcite, magnesium clays, '
            'magnesite and opal leave the water at stated thresholds (design guesses and bounds); the springs study '
            'and stage 1 give the share of the weathering that reaches the waters and the floors\' uptake; sodium, '
            'potassium, chloride and sulphate leave only with pore water and the water the floors take. Ion pairs are '
            'not followed. Beyond an ionic strength '
            'of 0.5 mol/kg the activities fail, and the masses alone are reported.')
READING_RULE = ('Concentrations in mmol/kg of water unless a key names g/kg or mg/m3; alkalinity in meq/kg. Years from '
                'the water\'s delivery. corners are the carbon cycle\'s low, central and high weathering. Each run is a '
                'delivery, a set of sinks and a fate for the sediment; pH is null where the ionic strength exceeds 0.5. '
                'co2 keys are Gt C: taken from the air into the waters and their carbonate, and given back by '
                'calcite and clays.')


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path, schema=None):
    data = json.loads((ROOT / path).read_text())
    if schema and data.get('schema') != schema:
        raise ValueError(f'{path}: unexpected schema {data.get("schema")}')
    return data


def rnd(x, sig=6):
    if isinstance(x, dict):
        return {k: rnd(v, sig) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [rnd(v, sig) for v in x]
    if isinstance(x, bool) or not isinstance(x, float) or x == 0.0 or not math.isfinite(x):
        return x
    return round(x, sig - 1 - int(math.floor(math.log10(abs(x)))))


# --- The setting -------------------------------------------------------------------------------------------------

def setting(data):
    atlas, drainage, gates = data[ATLAS], data[DRAINAGE], data[GATES]
    discharge = drainage['rivers']['discharge_to_sea_m3s'] * YEAR_S
    law = cc.DESSERT_A * math.exp(cc.DESSERT_B * (cc.LAND_AIR_K - 273.15))
    limb = data[LUNAR_CYCLE]['return_limb']
    rock = limb['rock_stock_top_metre']['highland_soil']
    mol = dict(ca=rock['Ca_kg_m2'] / cc.M_CA, mg=rock['Mg_kg_m2'] / cc.M_MG,
               na=rock['Na_kg_m2'] / cc.M_NA, k=rock['K_kg_m2'] / cc.M_K)
    eq = 2.0 * mol['ca'] + 2.0 * mol['mg'] + mol['na'] + mol['k']
    export = limb['soils']['highland_soil']['export_g_m2_yr']
    export_eq = (2.0 * export['Ca'][0] / cc.M_CA + 2.0 * export['Mg'][0] / cc.M_MG + export['Na'][0] / cc.M_NA
                 + export['K'][0] / cc.M_K)
    sea_p, loss = data[SEA_P]['partition'], data[SEA_P]['land_loss']
    # The land's erosion: the conservation study's 10-100 m/Myr of rock at the low end, the soil loss that the
    # phosphorus product's particles imply at the high (the sea phosphorus study).
    erosion = [min(loss['conservation_erosion_kg_yr']) / 1e12,
               max(max(loss['conservation_erosion_kg_yr']), max(loss['implied_erosion_kg_yr'])) / 1e12]
    dry_pa = data[CO2_STUDY]['design_air']['pressure_pa']
    return dict(
        sea_kg=atlas['water_mass_kg'], lake_kg=drainage['lakes']['volume_km3'] * 1e12,
        discharge_kg_yr=discharge * 1e3, weathering_eq_yr=law * discharge / 1e3,
        shares=dict(ca=2.0 * mol['ca'] / eq, mg=2.0 * mol['mg'] / eq, na=mol['na'] / eq, k=mol['k'] / eq),
        si_per_eq=export['Si'][0] / 28.086 / export_eq, so4_per_eq=export['S'][0] / 32.06 / export_eq,
        lake_share=sea_p['lake_passage_share_of_land_runoff']['tau_ge_0_yr'],
        particle_pass=sea_p['particulate_survival'],
        erosion_gt_yr=erosion,
        pco2_atm=cc.DESIGN_CO2_PPM * 1e-6 * dry_pa / cc.ATM_PA,
        lake_const=cc.fresh_constants(cc.LAND_AIR_K), fresh_sea_const=cc.fresh_constants(cc.SEA_SURFACE_K),
        p_kg_yr=dict(after_lakes=sea_p['to_seas_dissolved_kg_yr'], without_lakes=sea_p['if_lakes_passed_all_dissolved_kg_yr'],
                     particulate_after_lakes=sea_p['to_seas_particulate_kg_yr'],
                     particulate_without_lakes=loss['particulate_kg_yr']),
        lake_volume_m3=sea_p['lake_volume_m3'], groundwater=cc.ground(data), floors=floors(data))


def floors(data):
    """The lakes' and seas' floors' uptake (kg a year) over each of stage 1's spans, for each crust."""
    stage1 = data[cc.STAGE1]
    areas = {b['band']: b for b in stage1['setting']['bands']}
    out = {}
    for case in set(cc.CRUST_CASES.values()):
        rows, start = [], 0.0
        for end in FLOOR_SPANS:
            row = dict(start=start, end=float(end), seas_kg_yr=0.0, lakes_kg_yr=0.0)
            for name, band in stage1['cases'][case]['bands'].items():
                cols = band['columns']
                for key, col, area in (('seas_kg_yr', 'sea_floor', 'sea_m2'), ('lakes_kg_yr', 'lake_floor', 'lake_m2')):
                    if col in cols:
                        row[key] += cols[col]['by_year'][end]['uptake_m_yr_since_last'] * areas[name][area] * 1e3
            rows.append(row)
            start = float(end)
        out[case] = rows
    return out


def floor_at(rows, t, key):
    for row in rows:
        if t < row['end']:
            return row[key]
    return 0.0


def arriving(st, case, t):
    """The share of the catchment law's weathering that reaches the waters at year t, and its shares by ion: the
    springs study's calcium, magnesium and alkali shares, the alkalis split as the highland soil's."""
    g = st['groundwater'][case]
    ratio = cc.at_year(g['years'], g['ratio'], t)
    sh = {k: cc.at_year(g['years'], [x[k] for x in g['shares']], t) for k in ('ca', 'mg', 'mono')}
    rock = st['shares']
    alk = rock['na'] + rock['k']
    return ratio, dict(ca=sh['ca'], mg=sh['mg'], na=sh['mono'] * rock['na'] / alk, k=sh['mono'] * rock['k'] / alk)


def sediment_times(st, works):
    """When the lakes, then the seas, are full: the lakes keep all but the particles that pass them, the seas take
    those and their own deposits (calcite, magnesium minerals and opal) and, once the lakes are full, everything the
    land loses. The earliest pairs the high erosion, the loose dry density and the porous deposits of the natural
    seas; the latest the low erosion, the dense sediment and the managed seas' compact deposits. The same without the
    deposits gives the clastic sediment's share of the filling."""
    sea_m3 = st['sea_kg'] / 1e3

    def fill(e, rho, passing, dep):
        clastic = e / rho
        lakes = st['lake_volume_m3'] / ((1.0 - passing) * clastic)
        before = (passing * clastic + dep) * lakes
        seas = sea_m3 / (passing * clastic + dep) if before >= sea_m3 else lakes + (sea_m3 - before) / (clastic + dep)
        return lakes, seas
    out = {}
    for name, e, rho, passing, dep in (
            ('earliest', st['erosion_gt_yr'][1] * 1e12, DRY_DENSITY[0], st['particle_pass'][1],
             works['deposits_in_the_seas']['natural']['deposit_km3_yr'][1] * 1e9),
            ('latest', st['erosion_gt_yr'][0] * 1e12, DRY_DENSITY[1], st['particle_pass'][0],
             works['deposits_in_the_seas']['managed']['deposit_km3_yr'][0] * 1e9)):
        lakes, seas = fill(e, rho, passing, dep)
        out[name] = dict(erosion_gt_yr=e / 1e12, dry_density_kg_m3=rho, deposits_km3_yr=dep / 1e9,
                         lakes_full_yr=lakes, seas_full_yr=seas, seas_full_without_deposits_yr=fill(e, rho, passing, 0.0)[1])
    return dict(lakes_full_yr=[out['earliest']['lakes_full_yr'], out['latest']['lakes_full_yr']],
                seas_full_yr=[out['earliest']['seas_full_yr'], out['latest']['seas_full_yr']],
                seas_full_without_deposits_yr=[out['earliest']['seas_full_without_deposits_yr'],
                                               out['latest']['seas_full_without_deposits_yr']],
                ends=out)


# --- The waters ----------------------------------------------------------------------------------------------------

def box(mass, delivery, kind):
    b = dict(kind=kind, mass=mass, ta=0.0, ca=0.0, mg=0.0, na=0.0, k=0.0, cl=0.0, so4=0.0, si=0.0, dic=0.0)
    share = dict(brackish=BRACKISH_SHARE, salted=1.0).get(delivery, 0.0)
    for key in ('na', 'k', 'mg', 'ca', 'cl', 'so4'):
        b[key] = share * SEAWATER[key] * mass
    b['ta'] = share * SEAWATER['alkalinity'] * mass
    return b


def salinity(b):
    """Dissolved salts, g/kg, with the alkalinity counted as bicarbonate and silica as SiO2."""
    g = (b['ca'] * cc.M_CA + b['mg'] * cc.M_MG + b['na'] * cc.M_NA + b['k'] * cc.M_K + b['cl'] * M_CL
         + b['so4'] * M_SO4 + max(b['ta'], 0.0) * M_HCO3 + b['si'] * M_SIO2)
    return g / b['mass']


def ion_salinity(b):
    """The salinity the seawater constants read: the ions alone, without the uncharged silica."""
    return salinity(b) - b['si'] * M_SIO2 / b['mass']


def ionic(b):
    return 0.5 * (4.0 * b['ca'] + 4.0 * b['mg'] + b['na'] + b['k'] + b['cl'] + 4.0 * b['so4'] + max(b['ta'], 0.0)) / b['mass']


def valid(b, salted, margin=1.0):
    """Whether the box lies where its constants hold: Davies activities in fresh water, the seawater constants'
    salinity range in a salted sea; margin widens the range for a step's inputs before its minerals form."""
    if salted:
        return VALID_SALINITY_G_KG[0] <= ion_salinity(b) <= VALID_SALINITY_G_KG[1] * margin
    return ionic(b) <= VALID_IONIC_STRENGTH * margin


def proxy(b):
    """The box as the carbon cycle's solver reads it, whose ionic strength takes the monovalent ions as 'mono': here
    sodium, potassium, chloride and four times the sulphate, so that its Davies activities see every ion."""
    return dict(mass=b['mass'], ta=b['ta'], ca=b['ca'], mg=b['mg'], dic=b['dic'],
                mono=b['na'] + b['k'] + b['cl'] + 4.0 * b['so4'])


def constants(st, b, salted):
    if salted:
        return cc.sea_constants(cc.SEA_SURFACE_K, max(ion_salinity(b), 1.0))
    return st['lake_const'] if b['kind'] == 'lake' else st['fresh_sea_const']


def speciate(st, b, salted):
    return cc.speciate(proxy(b), st['pco2_atm'], constants(st, b, salted))


def calcite(st, b, salted, omega):
    """Calcite down to the threshold. The carbon cycle's solver holds the activities at the water's starting ionic
    strength, which a long step's inputs raise, so it is repeated until the saturation settles at the threshold."""
    total = 0.0
    for _ in range(12):
        p = proxy(b)
        x = cc.precipitate(p, st['pco2_atm'], constants(st, b, salted), omega)
        b['ca'], b['ta'] = p['ca'], p['ta']
        total += x
        if x <= 0.0 or speciate(st, b, salted)['omega'] <= omega * (1.0 + 1e-6):
            break
    return total


def clays(st, b, salted):
    """Magnesium clays that bring a sea down to PH_CLAY: each magnesium takes two equivalents of alkalinity and one
    and a half silicons, and its carbon returns to the air."""
    if b['mg'] <= 0.0 or speciate(st, b, salted)['ph'] <= PH_CLAY:
        return 0.0

    def above(x):
        t = dict(b, mg=b['mg'] - x, ta=b['ta'] - 2.0 * x, si=max(b['si'] - SI_PER_MG_CLAY * x, 0.0))
        return speciate(st, t, salted)['ph'] - PH_CLAY
    top = min(b['mg'], 0.5 * b['ta']) * (1.0 - 1e-12)
    x = top if above(top) > 0.0 else brentq(above, 0.0, top, xtol=1e-3, rtol=1e-12)
    b['mg'] -= x
    b['ta'] -= 2.0 * x
    b['si'] = max(b['si'] - SI_PER_MG_CLAY * x, 0.0)
    return x


def silica_cap(st, b, salted, rule):
    if rule == 'life':
        limit = SILICA_LIFE_MOL_KG
    else:
        ph = speciate(st, b, salted)['ph']
        limit = 10 ** SILICA_LOG_K_25C * (1.0 + 10 ** (ph - SILICA_PKA))
    x = max(b['si'] - limit * b['mass'], 0.0)
    b['si'] -= x
    return x


def decay(b, inputs, k, dt):
    """Add a year's inputs over dt while pore water takes a share k a year of everything dissolved: x relaxes to
    inputs/k with time constant 1/k (k = 0 is plain accumulation)."""
    f = math.exp(-k * dt) if k > 0 else 1.0
    g = (1.0 - f) / k if k > 0 else dt
    lost = {}
    for key, rate in inputs.items():
        before = b[key]
        b[key] = before * f + rate * g
        lost[key] = before + rate * dt - b[key]
    return lost


def times():
    """The carbon cycle's steps to 10,000 years, then LOG_STEPS_PER_DECADE a decade, with every report year."""
    out = [0.0]
    for start, end, dt in SEGMENTS:
        n = round((end - start) / dt)
        out += [start + dt * (i + 1) for i in range(n)]
    t = out[-1]
    while t < END_YEARS * (1 - 1e-9):
        t = min(t * 10 ** (1.0 / LOG_STEPS_PER_DECADE), END_YEARS)
        out.append(t)
    return sorted(set(round(x, 6) for x in out + [float(y) for y in REPORT_YEARS]))


def molar_grams(d):
    g = dict(ca=cc.M_CA, mg=cc.M_MG, na=cc.M_NA, k=cc.M_K, cl=M_CL, so4=M_SO4, ta=M_HCO3, si=M_SIO2)
    return sum(v * g[k] for k, v in d.items() if k in g)


def run(st, corner, delivery, sinks, sediment):
    """One run from the delivery to END_YEARS, or until the seas are full of sediment. Settling sediment takes the
    middle of the land's erosion range at the denser dry density, so its lakes and seas fill at central times."""
    i = cc.CORNERS.index(corner)
    p = dict(fresh_factor=cc.FRESH_FACTOR[i], soil_age_yr=cc.SOIL_AGE_AT_START_YR[i])
    salted = delivery == 'salted'
    lake, sea = box(st['lake_kg'], 'fresh', 'lake'), box(st['sea_kg'], delivery, 'sea')
    for b in (lake, sea):
        b['dic'] = b['mass'] * speciate(st, b, salted and b is sea)['dic']
    held0 = lake['dic'] + sea['dic']
    omega_sea = cc.OMEGA_FRESH if not salted else speciate(st, sea, True)['omega']
    sh, mg_rule, si_rule, soda_rule = (st['shares'], SINKS[sinks]['magnesium'], SINKS[sinks]['silica'],
                                       SINKS[sinks]['soda'])
    tau_lake = st['lake_kg'] / (st['discharge_kg_yr'] * st['lake_share'])
    erosion = 0.5 * sum(st['erosion_gt_yr']) * 1e12
    passing = 0.5 * sum(st['particle_pass'])
    lake_full = st['lake_volume_m3'] * DRY_DENSITY[1] / ((1.0 - passing) * erosion)
    porosity = 1.0 - DRY_DENSITY[1] / GRAIN_DENSITY
    cum = {k: 0.0 for k in ('weathering_eq', 'calcite', 'clays', 'magnesite', 'opal', 'pore_salts_g', 'sediment_m3',
                            'soda_works', 'sulphate_works', 'soda_minerals', 'deposit_m3', 'floor_alkalinity_eq',
                            'floor_carbon', 'floor_salts_g')}
    case = cc.CRUST_CASES[corner]
    floor_rows = st['floors'][case]
    # Each mole the sea deposits builds this much deposit (m3), at the middle porosity.
    build = {k: MOLAR_MASS[k] * 1e-3 / MINERAL_DENSITY[k] / (1.0 - 0.5 * sum(DEPOSIT_POROSITY))
             for k in ('calcite', 'magnesium_clays', 'opal')}
    rows, milestones, filled = [state(st, 0.0, sea, lake, salted, cum, held0)], {}, None
    grid = times()
    last = dict(sea)
    active_calcite = active_clays = False
    for t0, t1 in zip(grid[:-1], grid[1:]):
        dt = t1 - t0
        tm = 0.5 * (t0 + t1)
        ratio, sh = arriving(st, case, tm)
        w = st['weathering_eq_yr'] * cc.weathering_factor(tm, p) * ratio
        cum['weathering_eq'] += w * dt
        k_lake_floor = floor_at(floor_rows, tm, 'lakes_kg_yr') / lake['mass']
        k_sea_floor = floor_at(floor_rows, tm, 'seas_kg_yr') / sea['mass']
        cum['floor_carbon'] += (k_lake_floor * lake['dic'] + k_sea_floor * sea['dic']) * dt
        so4 = w * st['so4_per_eq']
        inflow = dict(ca=0.5 * w * sh['ca'], mg=0.5 * w * sh['mg'], na=w * sh['na'], k=w * sh['k'], so4=so4,
                      ta=w - 2.0 * so4, si=w * st['si_per_eq'])
        if mg_rule == 'steered':                    # the magnesite route: the rivers' magnesium made into MgCO3
            cum['magnesite'] += inflow['mg'] * dt
            inflow['ta'] -= 2.0 * inflow['mg']
            inflow['mg'] = 0.0
        if soda_rule == 'works':                    # soda works take the rivers' sodium, potassium and sulphate
            cum['soda_works'] += (inflow['na'] + inflow['k']) * dt
            cum['sulphate_works'] += inflow['so4'] * dt
            inflow['ta'] -= inflow['na'] + inflow['k'] - 2.0 * inflow['so4']
            inflow['na'] = inflow['k'] = inflow['so4'] = 0.0
        # The lakes: their share of the weathering, relaxing toward steady flow-through, and their outflow.
        to_lake = {k: v * st['lake_share'] for k, v in inflow.items()}
        k_lake = 1.0 / tau_lake + k_lake_floor
        lost = decay(lake, to_lake, k_lake, dt)
        outflow = {k: v / (tau_lake * k_lake) for k, v in lost.items()}
        down = {k: v * k_lake_floor / k_lake for k, v in lost.items()}
        cum['floor_alkalinity_eq'] += down['ta']
        cum['floor_salts_g'] += molar_grams(down)
        cum['calcite'] += calcite(st, lake, False, cc.OMEGA_FRESH)
        # The seas: the rest of the weathering, the lakes' outflow and, if sediment settles, pore water's burial.
        k_pw = 0.0
        if sediment == 'settling':
            v_rate = erosion * (1.0 if t1 > lake_full else passing) / DRY_DENSITY[1]
            k_pw = porosity * v_rate * 1e3 / sea['mass']
            cum['sediment_m3'] += v_rate * dt
        to_sea = {k: inflow[k] * (1.0 - st['lake_share']) + outflow[k] / dt for k in inflow}
        to_sea['cl'] = 0.0
        # The water delivered to replace what the floors take carries the delivery's own salts.
        share = dict(brackish=BRACKISH_SHARE, salted=1.0).get(delivery, 0.0)
        if share > 0.0 and k_sea_floor > 0.0:
            top = k_sea_floor * sea['mass'] * share
            for key in ('na', 'k', 'mg', 'ca', 'cl', 'so4'):
                to_sea[key] = to_sea.get(key, 0.0) + top * SEAWATER[key]
            to_sea['ta'] += top * SEAWATER['alkalinity']
        k_sea = k_pw + k_sea_floor
        lost = decay(sea, to_sea, k_sea, dt)
        pore = {k: v * k_pw / k_sea for k, v in lost.items()} if k_sea > 0 else lost
        down = {k: v * k_sea_floor / k_sea for k, v in lost.items()} if k_sea > 0 else {}
        cum['pore_salts_g'] += molar_grams(pore)
        cum['floor_alkalinity_eq'] += down.get('ta', 0.0)
        cum['floor_salts_g'] += molar_grams(down)
        # The minerals by the constants where the step ends inside their range (a step's inputs may carry the box
        # past it until calcite and clays take them), by the frozen caps below where it does not.
        trial, x_ca, x_mg, x_si = dict(sea), 0.0, 0.0, 0.0
        if not valid(trial, salted, TRIAL_MARGIN):
            # A long step's calcium (and, with clays forming, magnesium) leaves as it did in the last step before the
            # constants see the rest: the sea stood at its calcite threshold then.
            for key, active in (('ca', active_calcite), ('mg', active_clays), ('si', True)):
                gained = trial[key] - last[key]
                if active and gained > 0.0:
                    trial[key] = last[key]
                    if key == 'si':
                        cum['opal'] += gained
                        cum['deposit_m3'] += gained * build['opal']
                        continue
                    trial['ta'] -= 2.0 * gained
                    if key == 'ca':
                        x_ca += gained
                    else:
                        x_mg += gained
        if valid(trial, salted, TRIAL_MARGIN):
            x_si = silica_cap(st, trial, salted, si_rule)          # first, since the salinity sets the constants
            x_ca += calcite(st, trial, salted, omega_sea)
            if mg_rule == 'clays':
                x_mg += clays(st, trial, salted)
        ok = valid(trial, salted)
        if ok:
            sea.update(trial)
            if x_ca > 0.0:
                milestones.setdefault('calcite_in_the_seas', t1)
            cum['calcite'] += x_ca
            cum['clays'] += x_mg
            cum['opal'] += x_si
            cum['deposit_m3'] += x_ca * build['calcite'] + x_mg * build['magnesium_clays'] + x_si * build['opal']
            sea['dic'] = sea['mass'] * speciate(st, sea, salted)['dic']
            sp = speciate(st, sea, salted)
            active_calcite = x_ca > 0.0 or sp['omega'] >= omega_sea * (1.0 - 1e-3)
            active_clays = x_mg > 0.0 or (mg_rule == 'clays' and sp['ph'] >= PH_CLAY - 1e-3)
        else:
            # Beyond the constants' range the minerals keep calcium, magnesium (with clays) and silica where they
            # stood, taking what arrives, each calcium and magnesium with two equivalents of alkalinity.
            for key, keep in (('ca', True), ('mg', mg_rule == 'clays'), ('si', True)):
                gained = sea[key] - last[key]
                if keep and gained > 0.0:
                    sea[key] = last[key]
                    if key != 'si':
                        sea['ta'] -= 2.0 * gained
                    cum[dict(ca='calcite', mg='clays', si='opal')[key]] += gained
                    cum['deposit_m3'] += gained * build[dict(ca='calcite', mg='magnesium_clays', si='opal')[key]]
        if salinity(sea) > BRINE_CEILING_G_KG:
            # The brine ceiling: soda, sodium sulphate and magnesium carbonates deposit what arrives, each sodium and
            # potassium with its equivalent of alkalinity, each sulphate with two sodiums and each magnesium with two.
            gained = {key: max(sea[key] - last[key], 0.0) for key in ('na', 'k', 'so4', 'mg')}
            for key, x in gained.items():
                sea[key] -= x
            sea['ta'] -= gained['na'] + gained['k'] - 2.0 * gained['so4'] + 2.0 * gained['mg']
            cum['soda_minerals'] += gained['na'] + gained['k']
            cum['clays'] += gained['mg']                    # as magnesium carbonates, counted with the clays
            milestones.setdefault('brine_ceiling', t1)
        cum['opal'] += silica_cap(st, lake, False, si_rule)
        lake['dic'] = lake['mass'] * speciate(st, lake, False)['dic']
        last = dict(sea)
        s = salinity(sea)
        for name, hit in (('fresh_limit', s >= FRESH_LIMIT_G_KG and delivery == 'fresh'),
                          ('earth_salinity', s >= SEAWATER_G_KG and not salted),
                          ('salinity_up_a_tenth', salted and s >= 1.1 * SEAWATER_G_KG),
                          ('beyond_the_constants', not ok)):
            if hit and name not in milestones:
                milestones[name] = t1
        if ok and 'ph_9' not in milestones and speciate(st, sea, salted)['ph'] >= PH_CLAY - 1e-6:
            milestones['ph_9'] = t1
        if sediment == 'settling' and filled is None and cum['sediment_m3'] + cum['deposit_m3'] >= sea['mass'] / 1e3:
            filled = t1
        if t1 in REPORT_SET or t1 == grid[-1] or filled == t1:
            rows.append(state(st, t1, sea, lake, salted, cum, held0))
        if filled is not None:
            break
    return dict(corner=corner, delivery=delivery, sinks=sinks, sediment=sediment, rows=rows,
                milestones_yr=milestones, seas_full_yr=filled, lakes_full_yr=lake_full if sediment == 'settling' else None,
                calcite_threshold=omega_sea)


def state(st, t, sea, lake, salted, cum, held0):
    sp = speciate(st, sea, salted) if valid(sea, salted) else None
    m = sea['mass']
    held = sea['dic'] + lake['dic']
    took = cum['weathering_eq']
    buried = cum['calcite'] + cum['magnesite']
    return dict(
        year=t, salinity_g_kg=salinity(sea), ionic_strength=ionic(sea), ph=sp['ph'] if sp else None,
        calcite_saturation=sp['omega'] if sp else None, alkalinity_meq_kg=1e3 * sea['ta'] / m,
        ca=1e3 * sea['ca'] / m, mg=1e3 * sea['mg'] / m, na_k=1e3 * (sea['na'] + sea['k']) / m, cl=1e3 * sea['cl'] / m,
        so4=1e3 * sea['so4'] / m, si=1e3 * sea['si'] / m, lake_ph=speciate(st, lake, False)['ph'],
        lake_salinity_g_kg=salinity(lake),
        buried=dict(calcite_gt_c=cum['calcite'] * cc.M_C / 1e15, magnesite_gt_c=cum['magnesite'] * cc.M_C / 1e15,
                    magnesium_clays_gt=cum['clays'] * 0.25 * 647.8 / 1e15, opal_gt=cum['opal'] * M_SIO2 / 1e15,
                    soda_minerals_mol=cum['soda_minerals'], soda_works_mol=cum['soda_works'],
                    sulphate_works_mol=cum['sulphate_works'],
                    pore_water_salts_gt=cum['pore_salts_g'] / 1e15, sediment_km3=cum['sediment_m3'] / 1e9,
                    deposits_km3=cum['deposit_m3'] / 1e9),
        floors=dict(alkalinity_eq=cum['floor_alkalinity_eq'], salts_gt=cum['floor_salts_g'] / 1e15,
                    carbon_gt_c=cum['floor_carbon'] * cc.M_C / 1e15),
        co2_gt_c=dict(weathering_took=took * cc.M_C / 1e15, waters_gained=(held - held0) * cc.M_C / 1e15,
                      carbonate=buried * cc.M_C / 1e15, floors=cum['floor_carbon'] * cc.M_C / 1e15,
                      returned=(took - (held - held0) - buried - cum['floor_carbon']) * cc.M_C / 1e15))


# --- Phosphorus --------------------------------------------------------------------------------------------------

def phosphorus(st, sed, residence, start_mg_m3, extra_kg_yr=0.0, salted=False, end=END_YEARS):
    """The sea phosphorus study's mixed box, dP/dt = I - P/tau, with the dissolved input after the lakes until the
    lakes are full and the land's whole dissolved loss after; a salted sea also frees the reactive share of the
    particles. Returns mg/m3 at REPORT_YEARS for the low and high inputs."""
    vol = st['sea_kg'] / 1e3
    out = {}
    for e in (0, 1):
        before = st['p_kg_yr']['after_lakes'][e] + extra_kg_yr
        after = st['p_kg_yr']['without_lakes'][e] + extra_kg_yr
        if salted:
            before += P_REACTIVE_SHARE[e] * st['p_kg_yr']['particulate_after_lakes'][e]
            if st['p_kg_yr']['particulate_without_lakes']:
                after += P_REACTIVE_SHARE[e] * st['p_kg_yr']['particulate_without_lakes'][e]
        t_switch = sed['lakes_full_yr'][1 - e]               # low input with late filling, high with early
        series = []
        for y in REPORT_YEARS:
            p0 = start_mg_m3 * vol * 1e-6                     # kg
            if y <= t_switch:
                p = before * residence + (p0 - before * residence) * math.exp(-y / residence)
            else:
                ps = before * residence + (p0 - before * residence) * math.exp(-t_switch / residence)
                p = after * residence + (ps - after * residence) * math.exp(-(y - t_switch) / residence)
            series.append(p / vol * 1e6)
        out['low' if e == 0 else 'high'] = series
    return out


def works(st):
    """What each sink takes a year once the soils have aged to the catchment rate, with the mature Moon's share of it
    reaching the waters at the springs' make-up (the central crust), Mt a year: calcite from the
    calcium, magnesite or magnesium clays from the magnesium, soda from the sodium and potassium (as bicarbonates) and
    sodium sulphate from the sulphate."""
    ratio, sh = arriving(st, 'central', 1.0e6)
    w = st['weathering_eq_yr'] * ratio
    so4 = w * st['so4_per_eq']
    mt = dict(calcite=0.5 * w * sh['ca'] * 100.087 / 1e12, magnesite=0.5 * w * sh['mg'] * 84.314 / 1e12,
              magnesium_clays=0.5 * w * sh['mg'] * 0.25 * 647.8 / 1e12,
              soda=(w * sh['na'] * 84.007 + w * sh['k'] * 100.115) / 1e12, sodium_sulphate=so4 * 142.04 / 1e12,
              opal=w * st['si_per_eq'] * M_SIO2 / 1e12)
    km3 = {k: v * 1e9 / MINERAL_DENSITY[k] / 1e9 for k, v in mt.items()}
    sea_km3 = st['sea_kg'] / 1e12
    # The natural seas deposit calcite, magnesium clays and opal; the managed ones calcite and opal, with the magnesite
    # and the soda made at the rivers.
    deposits = {}
    for name, minerals in (('natural', ('calcite', 'magnesium_clays', 'opal')), ('managed', ('calcite', 'opal'))):
        solid = sum(km3[m] for m in minerals)
        built = [solid / (1.0 - f) for f in DEPOSIT_POROSITY]
        deposits[name] = dict(solids_km3_yr=solid, deposit_km3_yr=built,
                              seas_full_yr=[sea_km3 / built[1], sea_km3 / built[0]])
    return dict(weathering_eq_yr=w, share_reaching_the_waters=ratio, shares=sh, co2_gt_c_yr=w * cc.M_C / 1e15,
                minerals_mt_yr=mt, minerals_km3_yr=km3,
                carbonate_burial_gt_c_yr=0.5 * w * sh['ca'] * cc.M_C / 1e15,
                magnesite_burial_gt_c_yr=0.5 * w * sh['mg'] * cc.M_C / 1e15, deposits_in_the_seas=deposits)


def results():
    data = {p: load(p, s) for p, s in INPUTS.items()}
    st = setting(data)
    wk = works(st)
    sed = sediment_times(st, wk)
    runs = []
    for delivery in DELIVERIES:
        for sinks in SINKS:
            for corner in cc.CORNERS:
                for sediment in SEDIMENT:
                    if delivery == 'fresh_phosphorus':
                        continue                              # its major ions are the fresh delivery's
                    runs.append(run(st, corner, delivery, sinks, sediment))
    p = {}
    for residence in P_RESIDENCE_YR:
        key = f'{residence:.0f}_yr'
        p[key] = dict(fresh=phosphorus(st, sed, residence, 0.0), fresh_phosphorus=phosphorus(st, sed, residence, EARTH_P_MG_M3),
                      salted=phosphorus(st, sed, residence, 0.0, salted=True))
    carbon = data[CARBON]
    return dict(
        schema=SCHEMA,
        producer=dict(files={f: sha256(ROOT / f) for f in FILES}, inputs={k: sha256(ROOT / k) for k in INPUTS}),
        evidence=EVIDENCE, reading_rule=READING_RULE, report_years=list(REPORT_YEARS),
        setting=rnd(dict({k: v for k, v in st.items() if not k.endswith('_const') and k != 'groundwater'}), 5),
        sediment=rnd(sed), works_at_the_catchment_rate=rnd(wk), runs=rnd(runs), phosphorus_mg_m3=rnd(p))


def dumps(product) -> str:
    return json.dumps(product, indent=1, allow_nan=False) + '\n'


def main(argv=None) -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(dumps(results()))
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
