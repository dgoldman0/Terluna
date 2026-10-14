"""The Open Moon's carbon cycle as a box model, under the parameters the project now holds.

    python -m research.studies.carbon_cycle.run
    # -> research/studies/carbon_cycle/results/carbon_cycle.json

Resources' question 7 and the CO2 study's next work. From the seeding of the biosphere to year 10,000 the model tracks:
- the air's CO2 and methane, as carbon;
- the land's plants and soils over C''s land uses with the food blend, from Earth's biome inventories placed on the
  products' rain regions, with peat in the wet margins and carbonate in the dry soils;
- the sky's aerophytes and sky reefs at the sky ecology's covers;
- the rain-fed lakes and the seas: alkalinity, calcium, magnesium and dissolved inorganic carbon, in CO2 equilibrium
  with the air at the design temperatures, and the calcite they precipitate;
- organic burial in the lakes and seas, and the lake deltas dredged for phosphorus;
- people's stocks: the food chain, bodies and built stocks;
- the mineral strand: the calcium, magnesium, sodium and potassium weathering frees, in the waters, the carbonates and
  the kilns' lime;
- the groundwater (the springs study and stage 1): the share of the catchment law's weathering that reaches the waters
  and its make-up, the calcite the rock keeps, and the carbon the filling crust takes down with the water it draws from
  the land, the lakes and the seas.

Two modes run. Held: the managed return (kilns on the carbonate as it forms, the magnesite route, dredging) and an
import keep CO2 at 400 ppm, and the model reports what each must supply. Free: once the biosphere is grown, the kilns,
the magnesite route and the import stop, and CO2 follows the sinks; weathering slows with the cooling and the plants'
weaker help, and the land's and the sky's stocks follow the plants' production as CO2 falls.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

from scipy.optimize import brentq

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_RADIUS
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'carbon_cycle.json'
SCHEMA = 'terluna.research.carbon-cycle/1'
FILES = ('research/studies/carbon_cycle/run.py',)

SCREEN = 'resources/results/first_screen.json'
SEA_P = 'resources/results/sea_phosphorus.json'
INDUSTRY = 'resources/results/industry.json'
CO2_STUDY = 'research/studies/atmospheric_co2/results/atmospheric_co2.json'
LEDGERS = 'research/studies/domain_synthesis/results/ledgers.json'
FOOD = 'research/studies/food_systems/results/food_systems.json'
SKY = 'research/studies/sky_ecology/results/sky_ecology.json'
LUNAR_CYCLE = 'research/studies/lunar_cycle_ecology/results/lunar_cycle_ecology.json'
ECO_SCREEN = 'biosphere/ecology/results/first_screen.json'
ECO_PEOPLE = 'biosphere/ecology/results/people.json'
ATLAS = 'geography/results/atlas.json'
DRAINAGE = 'geography/results/drainage.json'
FORCING = 'atmosphere/radiative_convective/results/inverse_climate.csv'
SPRINGS = 'research/studies/springs/results/springs.json'
STAGE1 = 'geography/results/groundwater_uptake.json'
INPUTS = {
    SCREEN: 'terluna.resources.first-screen/1',
    SEA_P: 'terluna.resources.sea-phosphorus/1',
    INDUSTRY: 'terluna.resources.industry/1',
    CO2_STUDY: 'terluna.research.atmospheric-co2/1',
    LEDGERS: 'terluna.research.domain-synthesis-ledgers/1',
    FOOD: 'terluna.research.food-systems/1',
    SKY: 'terluna.research.sky-ecology/1',
    LUNAR_CYCLE: 'terluna.research.lunar-cycle-ecology/1',
    ECO_SCREEN: 'terluna.biosphere.ecology-first-screen/1',
    ECO_PEOPLE: 'terluna.biosphere.ecology-people/1',
    ATLAS: 'terluna.geography.atlas/1',
    DRAINAGE: 'terluna.geography.drainage/1',
    FORCING: None,
    SPRINGS: 'terluna.research.springs/1',
    STAGE1: 'terluna.geography.groundwater-uptake/1',
}

YEAR_S = JULIAN_DAY * JULIAN_YEAR_DAYS
MOON_AREA = 4.0 * math.pi * MOON_RADIUS ** 2
M_C, M_CO2, M_O2, M_CH4 = 12.011, 44.009, 31.998, 16.043             # g/mol
M_CA, M_MG, M_NA, M_K = 40.078, 24.305, 22.990, 39.098
M_CAO, M_MGO, M_FEO = 56.077, 40.304, 71.844
MOL_PER_GT_C = 1e15 / M_C
ATM_PA = 101325.0

# --- The planned world -------------------------------------------------------------------------------------------
DESIGN_CO2_PPM = 400.0          # CO2 held at the level the planned climate was run with
# The corrected design case, A28_dim5_moon years 20-29 (climate/gcm/README.md, "The corrected design case"). No product
# carries these means, so they are named here: the global surface mean, the air over land (weathering and the lakes) and
# the sea surface (the seas).
DESIGN_SURFACE_K = 294.9
LAND_AIR_K = 294.3
SEA_SURFACE_K = 297.0
# The held methane level: the author's choice of 13 October, Earth's present methane with a third of its N2O, the
# balanced middle of three pairs weighed (domain/ecology, research/studies/domain_synthesis, gases.held). Its working
# case: biology at the sources cuts 90% of methane's sources, and the soils and contactors in open air oxidize the rest
# to CO2.
HELD_METHANE_PPM = 1.8663
HELD_SOURCE_CUT = 0.9

# --- Weathering --------------------------------------------------------------------------------------------------
# Dessert et al. (2003), basaltic catchments: CO2 consumed (mol km-2 yr-1) = runoff (mm/yr) x 323.44 exp(0.0642 T), T in
# C, as the CO2 study uses it; the total is the law times the land's discharge, whatever area it is spread over.
DESSERT_A, DESSERT_B = 323.44, 0.0642
AGE_EXPONENT = 0.61             # field rates fall as t^-0.61 (White & Brantley 2003)
# Design guess: the developed soils start at F0 times the catchment rate (the CO2 study's 1-10 for fine, glassy ground)
# and age along White & Brantley's law from their age at the start, falling toward the catchment rate, which erosion
# keeps renewing: F(t) = 1 + (F0 - 1)(1 + t/age)^-0.61.
FRESH_FACTOR = (1.0, math.sqrt(10.0), 10.0)
SOIL_AGE_AT_START_YR = (100.0, 100.0, 250.0)
HIGHLAND_CAO_WT = 0.157         # Apollo 16 soil (Lunar Sourcebook Table 7.15), to recover the rock stock's mass basis
HIGHLAND_FEO_WT = 0.051         # the same soil's FeO, oxidized as the rock weathers (2 FeO + 1/2 O2 -> Fe2O3)

# --- Groundwater -------------------------------------------------------------------------------------------------
# The springs study and stage 1. Only a share of the catchment law's weathering reaches the waters: nearly half of the
# rivers' flow is rain on the lakes, which touches no rock, groundwater lays part of its calcium in the rock, and while
# the crust fills it takes much of the recharge with its load. What arrives has the springs' calcium-poor make-up. The
# rock's calcite keeps one carbon for each calcium it holds. While the crust fills, the land's water carries the soil's
# CO2 down with it, and the lakes' and seas' floors draw off their dissolved load, which delivered fresh water
# replaces. Each corner takes stage 1's crust of the same end; before 1,000 years the 1,000-year state holds.
CRUST_CASES = dict(low='slow', central='central', high='fast')
SPRING_STATES = ((1000.0, 'young_1000_yr'), (1.0e4, 'young_10000_yr'), (1.0e6, 'mature'))
CRUST_SPANS = ((0.0, '10'), (10.0, '100'), (100.0, '1000'), (1000.0, '10000'), (1.0e4, '100000'))
SPAN_STATE = {'10': 'young_1000_yr', '100': 'young_1000_yr', '1000': 'young_1000_yr', '10000': 'young_10000_yr',
              '100000': 'mature'}
CENTRAL_SOILS = dict(fresh_factor=FRESH_FACTOR[1], soil_age_yr=SOIL_AGE_AT_START_YR[1])

# --- Land --------------------------------------------------------------------------------------------------------
# Earth's carbon in plants and the top metre of soil (IPCC TAR WG1 Table 3.2 and Watson et al. 2000, the WBGU column,
# as 't Mannetje 2007 prints it): area 10^6 km2, plants and soils Gt C.
BIOMES = {'tropical_forests': (17.6, 212.0, 216.0), 'tropical_savannas': (22.5, 66.0, 264.0),
          'deserts_semideserts': (45.5, 8.0, 191.0), 'wetlands': (3.5, 15.0, 225.0), 'croplands': (16.0, 3.0, 128.0)}
FOREST_RAIN_MM_DAY = 2.0        # design mapping: wet land in bands of at least this mean rain holds forest
WILD_GROWTH_YEARS = 50.0        # plants reach 90% of their stock in about 50 years (the CO2 study; Poorter et al. 2016)
MANAGED_GROWTH_YEARS = 5.0      # design guess for fields, settlements and pasture
SOIL_BUILD_KG_C_M2_YR = 0.033   # soil carbon builds at 33 g C m-2 a year on establishment (Post & Kwon 2000)
SLOW_SOIL_FACTOR = 5.0          # the slow-soil variant: soils build five times slower (design guess)
SEED_SHARE = 0.01               # plants and the sky seeded at 1% of their stock (design guess)
SOIL_START_SHARE = 0.05         # developed soils start with 5% of their organic carbon (design guess)
PEAT_G_C_M2_YR = (18.0, 23.0, 28.0)   # Holocene peat accumulation, 23 +/- 2, 18-19 at its least and 25-28 when warmest
                                      # (Loisel et al. 2014), applied to the wet margins
SOIL_CARBONATE_G_C_M2_YR = (0.12, 0.27, 0.42)   # pedogenic carbonate in Mojave Aridisols (Schlesinger 1985)

# --- Waters ------------------------------------------------------------------------------------------------------
OMEGA_FRESH = 3.0               # design guess: calcifiers hold fresh waters at three times calcite saturation
SEA_SALINITY = 35.0             # the salted case: Earth's seawater delivered with the water (textbook composition)
SEAWATER = dict(alkalinity=2.30e-3, ca=10.28e-3, mg=52.8e-3)   # eq/kg, mol/kg; the salted seas start here and their
                                # calcifiers hold the saturation they arrive with
PISTON_M_DAY = 1.0              # an air-sea piston velocity at the low end of Earth's waters, for the exchange check

# --- Burial ------------------------------------------------------------------------------------------------------
LAKE_BURIAL_G_C_M2_YR = (19.0, 33.5, 48.0)       # lakes and reservoirs (Mendonca et al. 2017)
# Earth's seas by setting (Duarte, Middelburg & Caraco 2005, Table 1): unvegetated shelves bury 17 g C m-2 a year on the
# tenth of the shelf that is depositional; the deep sea buries 6 Tg C a year (Berner 1982) over Earth's 3.6e14 m2 of
# ocean less 2.66e13 m2 of shelf and 2.7e12 m2 of estuaries and vegetated coast.
SHELF_G_C_M2_YR = 17.0 * 0.10
DEEP_G_C_M2_YR = 6.0e12 / (3.6e14 - 2.66e13 - 2.7e12)
# Organic carbon of deltaic and shelf sediment, 0.7-1.5% (Hedges & Keil 1995) with Berner's 0.75% in the middle, as
# Duarte et al. (2005) give them, applied to the dredged lake deltas.
DELTA_OC = (0.007, 0.0075, 0.015)

# --- People ------------------------------------------------------------------------------------------------------
FOOD_RESIDENCE_YR = 0.5         # design guess: half a year's harvest in stores, kitchens and compost
BODY_KG, BODY_CARBON = 62.0, 0.18   # an adult's mean mass (Walpole et al. 2012) and the body's carbon share (textbook)
BUILT_T_C_PER_PERSON = (0.2, 1.0, 5.0)   # design guess: wood and polymers in buildings and goods; the domain products
                                         # give the built stocks as lunar rock

# --- The managed return ------------------------------------------------------------------------------------------
CALCITE_GJ_T_CO2 = 4.08         # calcination at its minimum, 3.2 GJ a tonne of lime (Schorcht et al. 2013)
BEST_KILN_GJ_T_LIME = (3.2, 4.2)   # the best lime kilns (Schorcht et al. 2013)
MAGNESITE_GJ_T_CO2 = 2.57       # magnesite at its minimum (Schorcht et al. 2013)
MAGNESIUM_STEERED = 1.0         # the magnesite route: the share of the lakes' magnesium made into MgCO3 (design lever)
LIME_RECARBONATING = (0.25, 0.5)   # shares of lime that might take CO2 back, for the kiln load 1/(1 - f)

# --- The free run ------------------------------------------------------------------------------------------------
FEEDBACK_W_M2_K = (1.0, 1.3, 1.6)   # the GCM's 1.0-1.6 W m-2 K-1 (the CO2 study, from climate/gcm/compare.py --settle)
FERT = 0.4                      # plant-assisted weathering against CO2, [2R/(1 + R)]^FERT (Berner & Kothavala 2001)
# As CO2 falls the land's and the sky's stocks follow the CO2 study's light-saturated C3 leaf rate at 25 C, relative to
# the design's 48.6 Pa, down to zero at its compensation point: stocks taken to scale with production (design guess).
GROWN_SHARE = 0.95              # the biosphere counts as grown when its plants and soils reach 95% of their stock
THRESHOLDS_PA = {'comfortable_43_pa': 43.0, 'productive_c3_farming_35_pa': 35.0,
                 'diverse_ecosystems_29_pa': 29.0, 'hard_floor_20_pa': 20.0}   # the CO2 study's plant table

CORNERS = ('low', 'central', 'high')
SEGMENTS = ((0.0, 200.0, 0.5), (200.0, 2000.0, 2.0), (2000.0, 10000.0, 5.0))
OUT_YEARS = (0, 10, 20, 50, 100, 200, 300, 500, 700, 1000, 1500, 2000, 3000, 5000, 7000, 10000)
REPORT_YEARS = (100, 1000, 2000, 5000, 10000)

EVIDENCE = ('A box model on committed products: the air and its 400 ppm of CO2, the design case\'s temperatures, '
            'the drainage\'s runoff and lakes, the atlas\'s seas, C\'\'s land uses with the food blend, the gas '
            'ledger\'s methane, the sky ecology\'s aerophytes, resources\' dredged deltas and its kilns. Its rates are '
            'Earth measurements carried to the Moon (the basalt weathering law, field ageing, biome inventories, lake '
            'and marine burial, peat and desert carbonate) beside design guesses, each labelled in the assumptions; '
            'the groundwater\'s share of the weathering, the rock\'s calcite and the filling crust\'s uptake come '
            'from the springs study and stage 1; '
            'its carbonate chemistry uses published constants at the design temperatures. A runnable calculation '
            'with checked conservation and step convergence: ecosystems, sediments and the seas enter as well-mixed '
            'boxes with Earth rates, and no lunar measurement stands behind any rate.')
READING_RULE = ('Carbon in Gt C (1e12 kg of carbon) and Gt C a year unless a key names Gt CO2, mol or kg; gases in ppm '
                'of the dry air; waters in meq or mmol per kg; power in GW; years from the seeding of the biosphere. '
                'Corners pair the low ends of every input (the smallest stocks and sinks), their middles and their '
                'high ends; the middle of the fresh-soil factor is geometric. Held runs keep CO2 at 400 ppm: their '
                'import is what the managed return leaves, and their net sink is what the Moon\'s reservoirs take from '
                'the air before any return. Free runs start from the held run at the year the biosphere is grown. '
                'Rates at a year are the mean over the step that ends there; the first year carries the first step\'s.')


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path, schema=None):
    data = json.loads((ROOT / path).read_text())
    if schema and data.get('schema') != schema:
        raise ValueError(f'{path}: unexpected schema {data.get("schema")}')
    return data


def rnd(x, sig=4):
    if isinstance(x, dict):
        return {k: rnd(v, sig) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [rnd(v, sig) for v in x]
    if isinstance(x, bool) or not isinstance(x, float) or x == 0.0 or not math.isfinite(x):
        return x
    return round(x, sig - 1 - int(math.floor(math.log10(abs(x)))))


def gt(mol):
    return mol / MOL_PER_GT_C


# --- Carbonate chemistry -----------------------------------------------------------------------------------------

def fresh_constants(t_k):
    """Fresh water: Plummer & Busenberg (1982) for CO2's solubility (mol kg-1 atm-1), carbonic acid's K1 and K2 and
    calcite's Ksp, as activities; Harned & Owen's water product; the Davies A from its linear fit in temperature."""
    lt = math.log10(t_k)
    return dict(
        kh=10 ** (108.3865 + 0.01985076 * t_k - 6919.53 / t_k - 40.45154 * lt + 669365.0 / t_k ** 2),
        k1=10 ** (-356.3094 - 0.06091964 * t_k + 21834.37 / t_k + 126.8339 * lt - 1684915.0 / t_k ** 2),
        k2=10 ** (-107.8871 - 0.03252849 * t_k + 5151.79 / t_k + 38.92561 * lt - 563713.9 / t_k ** 2),
        ksp=10 ** (-171.9065 - 0.077993 * t_k + 2839.319 / t_k + 71.595 * lt),
        kw=10 ** -(4470.99 / t_k - 6.0875 + 0.01706 * t_k),
        kb=1.0, tb=0.0, davies=0.4883 + 8.074e-4 * (t_k - 273.15), fresh=True)


def sea_constants(t_k, s):
    """Seawater, stoichiometric on the total pH scale: Weiss (1974) K0, Lueker et al. (2000) K1 and K2, Dickson (1990)
    KB with Uppstrom's (1974) boron, Millero's (1995) Kw as the DOE handbook gives it, Mucci (1983) calcite."""
    x, ln, rs = t_k / 100.0, math.log, math.sqrt(s)
    k0 = math.exp(-60.2409 + 93.4517 / x + 23.3585 * ln(x) + s * (0.023517 - 0.023656 * x + 0.0047036 * x * x))
    pk1 = 3633.86 / t_k - 61.2172 + 9.67770 * ln(t_k) - 0.011555 * s + 0.0001152 * s * s
    pk2 = 471.78 / t_k + 25.9290 - 3.16967 * ln(t_k) - 0.01781 * s + 0.0001122 * s * s
    lnkb = ((-8966.90 - 2890.53 * rs - 77.942 * s + 1.728 * s * rs - 0.0996 * s * s) / t_k + 148.0248 + 137.1942 * rs
            + 1.62142 * s - (24.4344 + 25.085 * rs + 0.2474 * s) * ln(t_k) + 0.053105 * rs * t_k)
    lnkw = 148.96502 - 13847.26 / t_k - 23.6521 * ln(t_k) + (118.67 / t_k - 5.977 + 1.0495 * ln(t_k)) * rs - 0.01615 * s
    lksp = (-171.9065 - 0.077993 * t_k + 2839.319 / t_k + 71.595 * math.log10(t_k)
            + (-0.77712 + 0.0028426 * t_k + 178.34 / t_k) * rs - 0.07711 * s + 0.0041249 * s * rs)
    return dict(kh=k0, k1=10 ** -pk1, k2=10 ** -pk2, ksp=10 ** lksp, kw=math.exp(lnkw), kb=math.exp(lnkb),
                tb=0.0004157 * s / 35.0, davies=0.0, fresh=False)


def gammas(c, box):
    """Activity coefficients (monovalent, divalent): Davies in fresh water at an ionic strength from the box's cations
    and alkalinity; one in seawater, whose constants are stoichiometric."""
    if not c['fresh']:
        return 1.0, 1.0
    i = 0.5 * (4.0 * box['ca'] + 4.0 * box['mg'] + box['mono'] + max(box['ta'], 0.0)) / box['mass']
    r = math.sqrt(i)
    f = r / (1.0 + r) - 0.3 * i
    return 10 ** (-c['davies'] * f), 10 ** (-4.0 * c['davies'] * f)


def carbonate(ta, co2, c, g1, g2):
    """[H+], HCO3- and CO3-- (mol/kg) at alkalinity ta and dissolved CO2 co2: the charge balance
    ta = HCO3 + 2 CO3 + B(OH)4 + OH - H solved for H by Newton's method in ln H, kept inside its bracket."""
    k1, k2, kw = c['k1'] / (g1 * g1), c['k2'] / g2, c['kw'] / (g1 * g1)
    a, b, tb, kb = k1 * co2, k1 * k2 * co2, c['tb'], c['kb']
    lo, hi = -46.0, 0.0
    if ta > 1e-7:                                    # the carbonate-only root, a close start in alkaline water
        h = (a + math.sqrt(a * a + 8.0 * b * ta)) / (2.0 * ta)
    else:
        h = math.sqrt(a + kw)
    u = math.log(min(max(h, 1e-19), 0.5))
    for _ in range(200):
        h = math.exp(u)
        borate = tb * kb / (kb + h)
        f = a / h + 2.0 * b / (h * h) + kw / h + borate - h - ta
        if f > 0.0:
            lo = u
        else:
            hi = u
        move = f / (a / h + 4.0 * b / (h * h) + kw / h + borate * h / (kb + h) + h)
        if abs(move) < 1e-12:
            u += move
            break
        v = u + move
        u = v if lo < v < hi else 0.5 * (lo + hi)
    h = math.exp(u)
    return h, a / h, b / (h * h)


def speciate(box, pco2_atm, c):
    """The box's carbonate system in equilibrium with the air: pH, DIC (mol/kg) and calcite saturation."""
    g1, g2 = gammas(c, box)
    co2 = c['kh'] * pco2_atm
    h, hco3, co3 = carbonate(box['ta'] / box['mass'], co2, c, g1, g2)
    omega = g2 * g2 * (box['ca'] / box['mass']) * co3 / c['ksp']
    return dict(ph=-math.log10(h * g1), dic=co2 + hco3 + co3, co3=co3, omega=omega, co2=co2, hco3=hco3)


def precipitate(box, pco2_atm, c, omega_c):
    """Calcite (mol) that brings the box down to omega_c if it stands above: each mole takes a calcium, two
    equivalents of alkalinity and a carbon; the second carbon returns to the air as the water re-equilibrates."""
    g1, g2 = gammas(c, box)
    m, co2 = box['mass'], c['kh'] * pco2_atm

    def excess(x):
        _, _, co3 = carbonate((box['ta'] - 2.0 * x) / m, co2, c, g1, g2)
        return g2 * g2 * (box['ca'] - x) / m * co3 / c['ksp'] - omega_c
    if excess(0.0) <= 0.0:
        return 0.0
    top = min(box['ca'], 0.5 * box['ta']) * (1.0 - 1e-12)
    x = brentq(excess, 0.0, top, xtol=1e-6, rtol=1e-13)
    box['ca'] -= x
    box['ta'] -= 2.0 * x
    box['dic'] -= x
    return x


def dic_mol(box, pco2_atm, c):
    return box['mass'] * speciate(box, pco2_atm, c)['dic']


# --- The setting from the products -------------------------------------------------------------------------------

def forcing_sweep(path):
    """Net flux at the top of the lunar 1.2-atm column against CO2 (the radiative-convective carbon_dioxide sweep, with
    its 400 ppm base), at the fixed surface temperature the sweep uses."""
    rows = list(csv.DictReader((ROOT / path).read_text().splitlines()))
    pick = {}
    for r in rows:
        if (r['scenario'] == 'moon_1.2atm' and r['ts_k'] == '290.0' and r['shield'] == 'none'
                and r['humidity'] == 'manabe_wetherald' and r['stratosphere_k'] == '200.0'
                and r['surface_albedo'] == '0.13' and r['sweep'] in ('carbon_dioxide', 'earthlike_humidity')):
            pick[float(r['co2_ppm'])] = float(r['net_toa'])
    if DESIGN_CO2_PPM not in pick or len(pick) < 3:
        raise ValueError('The CO2 sweep lacks its points')
    return dict(sorted(pick.items()))


def forcing(ppm, sweep):
    """W/m2 against 400 ppm, linear in ln(ppm) between the sweep's points and along its end segments beyond them."""
    xs, ys = [math.log(k) for k in sweep], list(sweep.values())
    x = math.log(ppm)
    j = 1
    while j < len(xs) - 1 and x > xs[j]:
        j += 1
    y = ys[j - 1] + (ys[j] - ys[j - 1]) * (x - xs[j - 1]) / (xs[j] - xs[j - 1])
    return y - sweep[DESIGN_CO2_PPM]


def setting(data):
    """Everything the corners share, from the products."""
    screen, atlas, drainage = data[SCREEN], data[ATLAS], data[DRAINAGE]
    inv = screen['inventories']
    air_mol = inv['air_mol']
    co2 = data[CO2_STUDY]['design_air']
    dry_pa = co2['pressure_pa']
    discharge = drainage['rivers']['discharge_to_sea_m3s'] * YEAR_S                     # m3 a year to the seas
    land_m2 = (atlas['main_land_share'] + atlas['island_share']) * MOON_AREA
    # Dessert's law over the land's discharge: mol/km2/yr per mm of runoff, times mm km2 of water a year.
    law = DESSERT_A * math.exp(DESSERT_B * (LAND_AIR_K - 273.15))
    weathering_eq_yr = law * discharge / 1e3
    rock = data[LUNAR_CYCLE]['return_limb']['rock_stock_top_metre']['highland_soil']
    mol = dict(ca=rock['Ca_kg_m2'] / M_CA * 1e3, mg=rock['Mg_kg_m2'] / M_MG * 1e3,
               mono=rock['Na_kg_m2'] / M_NA * 1e3 + rock['K_kg_m2'] / M_K * 1e3)
    eq = 2.0 * mol['ca'] + 2.0 * mol['mg'] + mol['mono']
    rock_kg_m2 = rock['Ca_kg_m2'] / (HIGHLAND_CAO_WT * M_CA / M_CAO)
    o2_per_eq = rock_kg_m2 * HIGHLAND_FEO_WT / (M_FEO * 1e-3) / 4.0 / eq                # 4 FeO take one O2
    shares = dict(ca=2.0 * mol['ca'] / eq, mg=2.0 * mol['mg'] / eq, mono=mol['mono'] / eq)
    sea_p = data[SEA_P]
    lake_share = sea_p['partition']['lake_passage_share_of_land_runoff']['tau_ge_0_yr']
    shelves = sea_p['gathering']['shelves']['all_common_level_water']
    sea_m2 = atlas['water_share'] * MOON_AREA
    lake_m2 = drainage['lakes']['area_share'] * MOON_AREA
    shelf = shelves['under_200_m']
    sea_burial = (SHELF_G_C_M2_YR * shelf + DEEP_G_C_M2_YR * (1.0 - shelf)) * sea_m2
    industry = data[INDUSTRY]['places']['moon_surface']
    stocks = data[INDUSTRY]['people']
    n = data[LEDGERS]['choices']['moon_billion'] * 1e9
    renewal = [n * stocks['stock_t_per_person'][0] / stocks['lifetime_years'][1] / 1e9,
               n * stocks['stock_t_per_person'][1] / stocks['lifetime_years'][0] / 1e9]    # Gt a year
    leaf = data[CO2_STUDY]['photosynthesis']['25C']
    curve = sorted([(leaf['compensation_gamma_star_pa'], 0.0)]
                   + [(float(k), v) for k, v in leaf['relative_to_earth_today'].items()])
    return dict(
        air_mol=air_mol, dry_pa=dry_pa, o2_mol=inv['o2_kg'] / (M_O2 * 1e-3), o2_pa=co2['o2_pa'],
        held_co2_mol=air_mol * DESIGN_CO2_PPM * 1e-6, pco2_atm=DESIGN_CO2_PPM * 1e-6 * dry_pa / ATM_PA,
        tg_ch4_per_ppm=air_mol * 1e-6 * M_CH4 / 1e12,
        discharge_m3_yr=discharge, land_m2=land_m2, runoff_mm_yr=discharge / land_m2 * 1e3,
        dessert_mol_km2_mm=law, weathering_eq_yr=weathering_eq_yr, cation_shares=shares, o2_per_eq=o2_per_eq,
        rock_kg_per_eq=rock_kg_m2 / eq, lake_share=lake_share,
        sea_kg=atlas['water_mass_kg'], sea_m2=sea_m2, lake_kg=drainage['lakes']['volume_km3'] * 1e12, lake_m2=lake_m2,
        shelf_share=shelves['under_200_m'], sea_burial_mol_yr=sea_burial / M_C,
        lake_const=fresh_constants(LAND_AIR_K), fresh_sea_const=fresh_constants(SEA_SURFACE_K),
        salted_sea_const=sea_constants(SEA_SURFACE_K, SEA_SALINITY),
        heat_w_m2_per_tw=industry['heat_w_m2_per_tw_used'], sunlight_w_m2=industry['sunlight_w_m2'],
        sweep=forcing_sweep(FORCING), leaf_curve=curve, built_renewal_gt_yr=renewal,
        small_waters_m2=atlas['small_lakes_share'] * MOON_AREA, ground=ground(data))


def ground(data):
    """The groundwater's terms for each of stage 1's crusts: by year, the share of the catchment law's weathering that
    reaches the waters, its cations' shares and the calcite the rock keeps (the springs study's central chemistry); by
    span, the filling crust's uptake from the land, the lakes and the seas, with the carbon the land's water carries
    at the water table."""
    springs, stage1 = data[SPRINGS], data[STAGE1]
    areas = {b['band']: b for b in stage1['setting']['bands']}
    out = {}
    for case in set(CRUST_CASES.values()):
        states = springs['cases'][case]['states']
        years, ratio, shares, rock = [], [], [], []
        for year, name in SPRING_STATES:
            moon = states[name]['variants']['central']['moon']
            w = moon['weathering']
            c = w['cations_mol_yr']
            eq = 2.0 * c['ca'] + 2.0 * c['mg'] + c['mono']
            years.append(year)
            ratio.append(w['over_carbon_cycle_catchment'])
            shares.append(dict(ca=2.0 * c['ca'] / eq, mg=2.0 * c['mg'] / eq, mono=c['mono'] / eq))
            laid = moon['calcite_laid_gt_c_yr']
            rock.append((laid['saturated_rock'] + laid['unsaturated_rock']) * MOL_PER_GT_C)
        spans = []
        for start, end in CRUST_SPANS:
            entry = states[SPAN_STATE[end]]['variants']['central']['by_band']
            row = dict(start=start, end=float(end), seas_kg_yr=0.0, lakes_kg_yr=0.0, land_kg_yr=0.0, land_mol_c_yr=0.0)
            for name, band in stage1['cases'][case]['bands'].items():
                cols = band['columns']
                up = lambda key: cols[key]['by_year'][end]['uptake_m_yr_since_last'] if key in cols else 0.0
                row['seas_kg_yr'] += up('sea_floor') * areas[name]['sea_m2'] * 1e3
                row['lakes_kg_yr'] += up('lake_floor') * areas[name]['lake_m2'] * 1e3
                kg = up('land') * areas[name]['land_m2'] * 1e3
                row['land_kg_yr'] += kg
                if name in entry:
                    row['land_mol_c_yr'] += kg * entry[name]['water_table_entry']['dic_mmol_kg'] * 1e-3
            spans.append(row)
        out[case] = dict(years=years, ratio=ratio, shares=shares, rock_mol_yr=rock, spans=spans)
    return out


def at_year(years, values, t):
    """Linear in the logarithm of the year between the springs study's states, held before the first."""
    x = math.log10(max(t, years[0]))
    xs = [math.log10(y) for y in years]
    for (x0, v0), (x1, v1) in zip(zip(xs, values), zip(xs[1:], values[1:])):
        if x <= x1:
            return v0 + (v1 - v0) * (x - x0) / (x1 - x0)
    return values[-1]


def span_of(spans, t):
    for row in spans:
        if t < row['end']:
            return row
    return spans[-1]


def leaf_rate(pa, curve):
    """The CO2 study's light-saturated C3 rate against CO2 partial pressure, linear between its points, zero at its
    compensation point and flat beyond its last point."""
    if pa <= curve[0][0]:
        return 0.0
    for (x0, y0), (x1, y1) in zip(curve, curve[1:]):
        if pa <= x1:
            return y0 + (y1 - y0) * (pa - x0) / (x1 - x0)
    return curve[-1][1]


def production(ppm, st):
    """Production relative to the design's, from the leaf rate (1 at 400 ppm)."""
    rate = lambda q: leaf_rate(q * 1e-6 * st['dry_pa'], st['leaf_curve'])
    return rate(ppm) / rate(DESIGN_CO2_PPM)


def plant_weathering(ppm):
    """GEOCARB III's vascular-plant weathering against CO2, relative to 400 ppm."""
    r = ppm / DESIGN_CO2_PPM
    return (2.0 * r / (1.0 + r)) ** FERT


def land_classes(data, i):
    """The dry land's classes at corner i: Earth's biome stocks on the products' rain regions, C''s land uses with the
    blend, the wet margins and the parks with their region's wild stock."""
    ledger, people = data[LEDGERS]['land'], data[ECO_PEOPLE]['land']
    dry = ledger['dry_land_km2'] * 1e6
    regions = people['regions']
    rows = people['by_latitude']
    wet = sum(r['dry_land_share'] * r['wet_share'] for r in rows)
    if not math.isclose(wet, regions['wet_land'], rel_tol=2e-3):
        raise ValueError('The bands do not add up to the wet land')
    forest = sum(r['dry_land_share'] * r['wet_share'] for r in rows if r['rain_mm_day'] >= FOREST_RAIN_MM_DAY)
    total = sum(regions.values())
    area = dict(forest=dry * forest / total, savanna=dry * (regions['wet_land'] - forest) / total,
                desert=dry * (regions['fog_desert'] + regions['polar_dry_land']) / total)
    # Pair the corners: the low corner takes the most land used (the blend's high-cost end, 2 kg of heritage meat).
    def ends(pair):
        return (pair[1], 0.5 * (pair[0] + pair[1]), pair[0])[i]
    fields, settled = ends(ledger['blend_fields_km2']) * 1e6, ends(ledger['settled_km2']) * 1e6
    pasture = ledger['with_heritage'][('2_kg', '1_kg', '0.5_kg')[i]]['heritage_pasture_km2'] * 1e6
    wet_share = data[ECO_SCREEN]['trace_gases']['methane']
    wetland = (wet_share['wetland_share_scenario'][0], wet_share['central']['inputs']['wetland_share'],
               wet_share['wetland_share_scenario'][1])[i] * MOON_AREA
    dens = {k: (v[1] / v[0], v[2] / v[0]) for k, v in BIOMES.items()}                  # kg C/m2: plants, soils
    classes = {
        'forest': (area['forest'] - fields - settled - wetland, 'tropical_forests', WILD_GROWTH_YEARS),
        'savanna': (area['savanna'] - pasture, 'tropical_savannas', WILD_GROWTH_YEARS),
        'desert': (area['desert'], 'deserts_semideserts', WILD_GROWTH_YEARS),
        'wetland': (wetland, 'wetlands', WILD_GROWTH_YEARS),
        'fields': (fields, 'croplands', MANAGED_GROWTH_YEARS),
        'settlements': (settled, 'croplands', MANAGED_GROWTH_YEARS),
        'pasture': (pasture, 'tropical_savannas', MANAGED_GROWTH_YEARS),
    }
    out = {}
    for name, (a, biome, years) in classes.items():
        if a <= 0.0:
            raise ValueError(f'{name} has no land left')
        plants, soils = dens[biome]
        out[name] = dict(area_m2=a, biome=biome, plants_mol=a * plants * 1e3 / M_C,
                         soils_mol=0.0 if name == 'wetland' else a * soils * 1e3 / M_C, growth_years=years)
    return out, area, dict(fields=fields, settlements=settled, pasture=pasture, wetland=wetland)


def sky_stock(data, i):
    """The sky's standing carbon and growth rate at corner i: the sky ecology's cover, class shares and bodies at the
    reference light."""
    sky = data[SKY]
    cover = (min(sky['cover_basis']['range']), sky['cover_basis']['central'], max(sky['cover_basis']['range']))[i]
    shares = sky['phosphorus']['design_guesses']['class_shares_of_cover']
    rows = [r for r in sky['food_web']['production'] if r['gpp_c_kg_m2_year'] == 2.0]
    standing, rate = 0.0, 0.0
    for cls, share in shares.items():
        vals = sorted({r['standing_host_c_kg_m2'] for r in rows if r['cls'] == cls})
        standing += share * (vals[0], sum(vals) / len(vals), vals[-1])[i]
        g = sky['headlines']['growth_rate_per_year'][cls]
        rate += share * (g[0], 0.5 * (g[0] + g[1]), g[1])[i]
    return dict(cover=cover, standing_kg_c_m2=standing, growth_rate_per_year=rate,
                stock_mol=cover * MOON_AREA * standing * 1e3 / M_C)


def parameters(data, st, corner, methane='held', magnesium='natural', seas='fresh', soils='post_kwon'):
    i = CORNERS.index(corner)
    classes, regions, uses = land_classes(data, i)
    ledger = data[LEDGERS]
    m = ledger['gases']['methane']
    wild = m['wild_sources_tg_yr']['central']
    people = [c['sources_tg_yr']['farms'] + c['sources_tg_yr']['heritage_cattle'] for c in m['central']]
    people_tg = (people[0], 0.5 * (people[0] + people[1]), people[1])[i]
    food = data[FOOD]
    diet = food['basis']['diet']
    c_per_kcal = 6 * M_C / 180.156 / food['derived']['fermentation']['glucose_kcal_per_g']   # g C a kcal, as glucose
    n = ledger['choices']['moon_billion'] * 1e9
    loss = (diet['loss_factor'][0], 0.5 * sum(diet['loss_factor']), diet['loss_factor'][1])[i]
    grown_g_yr = n * diet['intake_kcal_day'] * JULIAN_YEAR_DAYS * loss * c_per_kcal
    lost = ledger['phosphorus']['lost']
    p_need = (lost[0], 0.5 * (lost[0] + lost[1]), lost[1])[i] * 1e9                     # kg P a year dredged back
    sed = data[SEA_P]['routes']['dredging']['sediment_kg_per_kg_p']
    sediment = (sed[0], 0.5 * (sed[0] + sed[1]), sed[1])[i] * p_need                    # kg of delta sediment a year
    fuel = data[SEA_P]['routes']['dredging']['fuel_mj_per_kg_p']                        # resources' dredging fuel
    ends = (fuel[0] / sed[0], fuel[1] / sed[1])
    dredge_mj = (ends[0], 0.5 * (ends[0] + ends[1]), ends[1])                           # MJ a kg of sediment
    soil_tau = BIOMES['tropical_forests'][2] / BIOMES['tropical_forests'][0] / SOIL_BUILD_KG_C_M2_YR
    if soils == 'slow':
        soil_tau *= SLOW_SOIL_FACTOR
    sky = sky_stock(data, i)
    return dict(
        corner=corner, methane=methane, magnesium=magnesium, seas=seas, soils=soils,
        classes=classes, regions_m2=regions, uses_m2=uses, soil_tau_yr=soil_tau,
        fresh_factor=FRESH_FACTOR[i], soil_age_yr=SOIL_AGE_AT_START_YR[i],
        lake_burial_mol_yr=LAKE_BURIAL_G_C_M2_YR[i] * st['lake_m2'] / M_C,
        peat_mol_yr=PEAT_G_C_M2_YR[i] * uses['wetland'] / M_C,
        soil_carbonate_mol_yr=SOIL_CARBONATE_G_C_M2_YR[i] * regions['desert'] / M_C,
        dredged_mol_yr=sediment * DELTA_OC[i] * 1e3 / M_C, dredged_sediment_kg_yr=sediment, p_need_kg_yr=p_need,
        delta_oc=DELTA_OC[i], dredge_mj_per_kg=dredge_mj[i],
        wild_ch4_mol_yr=wild * 1e12 / M_CH4, people_ch4_mol_yr=people_tg * 1e12 / M_CH4,
        soil_ch4_per_yr=m['central'][0]['soil_sink_tg_yr_per_ppm'] / st['tg_ch4_per_ppm'],
        sky=sky, people=n, food_mol=grown_g_yr * FOOD_RESIDENCE_YR / M_C + n * BODY_KG * BODY_CARBON * 1e3 / M_C,
        food_grown_mol_yr=grown_g_yr / M_C, built_mol=n * BUILT_T_C_PER_PERSON[i] * 1e6 / M_C,
        power_made_share=ledger['choices']['blend_weights']['power_made'],
        feedback_w_m2_k=FEEDBACK_W_M2_K[i], crust_case=CRUST_CASES[corner])


# --- The model ---------------------------------------------------------------------------------------------------

def weathering_factor(t, p):
    return 1.0 + (p['fresh_factor'] - 1.0) * (1.0 + t / p['soil_age_yr']) ** -AGE_EXPONENT


def sea_const(st, p):
    return st['salted_sea_const'] if p['seas'] == 'salted' else st['fresh_sea_const']


def initial_state(st, p):
    lake = dict(mass=st['lake_kg'], ta=0.0, ca=0.0, mg=0.0, mono=0.0, dic=0.0)
    sea = dict(mass=st['sea_kg'], ta=0.0, ca=0.0, mg=0.0, mono=0.0, dic=0.0)
    if p['seas'] == 'salted':
        sea.update(ta=SEAWATER['alkalinity'] * st['sea_kg'], ca=SEAWATER['ca'] * st['sea_kg'],
                   mg=SEAWATER['mg'] * st['sea_kg'])
    lake['dic'] = dic_mol(lake, st['pco2_atm'], st['lake_const'])
    sea['dic'] = dic_mol(sea, st['pco2_atm'], sea_const(st, p))
    s = dict(t=0.0, air_co2=st['held_co2_mol'], air_ch4=0.0,
             plants={k: SEED_SHARE * v['plants_mol'] for k, v in p['classes'].items()},
             soils={k: SOIL_START_SHARE * v['soils_mol'] for k, v in p['classes'].items()},
             peat=0.0, soil_carbonate=0.0, sky=SEED_SHARE * p['sky']['stock_mol'], food=p['food_mol'],
             built=p['built_mol'], lake=lake, sea=sea, lake_org=0.0, sea_org=0.0, calcite_sed=0.0, mgco3_sed=0.0,
             rock_calcite=0.0, crust=0.0,
             cao=0.0, mgo=0.0, delivered=dict(ca=lake['ca'] + sea['ca'], mg=lake['mg'] + sea['mg'],
                                              mono=lake['mono'] + sea['mono'],
                                              offset=sea['ta'] - 2 * sea['ca'] - 2 * sea['mg'] - sea['mono']),
             cum={k: 0.0 for k in ('import', 'kiln_calcite', 'kiln_mg', 'dredged', 'calcite', 'calcite_lake',
                                   'mgco3', 'rock_ca', 'rock_mg', 'rock_mono', 'rock_ca_dry', 'ch4_made', 'ch4_cut',
                                   'ch4_removed', 'ch4_soil', 'o2_burial', 'o2_rock', 'lake_burial', 'sea_burial',
                                   'weathering', 'floor_ta', 'floor_ca', 'floor_mg', 'floor_mono', 'floor_dic',
                                   'rock_ca_crust', 'crust_land')})
    s['initial_organic'] = organic(s)
    s['initial_total'] = carbon_total(s)
    s['omega_sea'] = OMEGA_FRESH if p['seas'] == 'fresh' else speciate(sea, st['pco2_atm'], sea_const(st, p))['omega']
    return s


def organic(s):
    """Organic carbon in the living and built stocks the biosphere's growth makes from CO2 (peat counted as burial)."""
    return sum(s['plants'].values()) + sum(s['soils'].values()) + s['sky'] + s['food'] + s['built']


def carbon_total(s):
    return (s['air_co2'] + s['air_ch4'] + sum(s['plants'].values()) + sum(s['soils'].values()) + s['peat']
            + s['soil_carbonate'] + s['sky'] + s['food'] + s['built'] + s['lake']['dic'] + s['sea']['dic']
            + s['lake_org'] + s['sea_org'] + s['calcite_sed'] + s['mgco3_sed'] + s['rock_calcite'] + s['crust'])


def land_organic(s):
    return sum(s['plants'].values()) + sum(s['soils'].values())


def mature_land(p):
    return sum(v['plants_mol'] + v['soils_mol'] for v in p['classes'].values())


def steps(t_start, scale=1.0):
    out = []
    for a, b, dt in SEGMENTS:
        n = int(round((b - a) / (dt * scale)))
        for k in range(n):
            t0 = a + (b - a) * k / n
            if t0 >= t_start - 1e-9:
                out.append((t0, (b - a) / n))
    return out


def logistic(x, x_max, rate, dt):
    if x_max <= 0.0 or x <= 0.0:
        return 0.0
    return x_max / (1.0 + (x_max / x - 1.0) * math.exp(-rate * dt))


def growth_rate(years):
    """Logistic rate that takes a stock from the seed share to 90% in the given years."""
    return math.log(9.0 * (1.0 - SEED_SHARE) / SEED_SHARE) / years


def step(s, st, p, dt, held, managed):
    """One step of length dt (years). Returns its fluxes (mol)."""
    t0 = s['t']
    f = {}
    mature_plants = sum(v['plants_mol'] for v in p['classes'].values())
    dev = sum(s['plants'].values()) / mature_plants                    # the land's development, burial and wild methane
    wet_dev = s['plants']['wetland'] / p['classes']['wetland']['plants_mol']
    ppm = s['air_co2'] / st['air_mol'] * 1e6
    prod = 1.0 if held else production(ppm, st)                       # stocks follow production as CO2 falls
    air = 0.0                                                          # carbon into the air this step
    # Land plants and soils, the wet margins' peat, the dry soils' carbonate, the sky.
    land0 = land_organic(s)
    for k, v in p['classes'].items():
        s['plants'][k] = logistic(s['plants'][k], prod * v['plants_mol'], growth_rate(v['growth_years']), dt)
        target = prod * v['soils_mol']
        s['soils'][k] = target + (s['soils'][k] - target) * math.exp(-dt / p['soil_tau_yr'])
    f['land'] = land_organic(s) - land0
    f['peat'] = p['peat_mol_yr'] * wet_dev * dt
    s['peat'] += f['peat']
    f['soil_carbonate'] = p['soil_carbonate_mol_yr'] * dt
    s['soil_carbonate'] += f['soil_carbonate']
    s['cum']['rock_ca_dry'] += f['soil_carbonate']
    sky0 = s['sky']
    s['sky'] = logistic(s['sky'], prod * p['sky']['stock_mol'], p['sky']['growth_rate_per_year'], dt)
    f['sky'] = s['sky'] - sky0
    air -= f['land'] + f['peat'] + f['soil_carbonate'] + f['sky']
    # Organic burial in the lakes and seas, and the deltas dredged back to the farms.
    f['lake_burial'] = p['lake_burial_mol_yr'] * dev * dt
    f['sea_burial'] = st['sea_burial_mol_yr'] * dev * dt
    f['dredged'] = min(p['dredged_mol_yr'] * dt, f['lake_burial'])
    s['lake_org'] += f['lake_burial'] - f['dredged']
    s['sea_org'] += f['sea_burial']
    air -= f['lake_burial'] + f['sea_burial'] - f['dredged']
    # Methane: made from the biosphere's carbon, a share oxidized at the sources when held, the soils' uptake, and the
    # removal that holds the level. Every oxidized molecule returns to CO2.
    made = (p['wild_ch4_mol_yr'] * dev + p['people_ch4_mol_yr']) * dt
    cut = HELD_SOURCE_CUT if p['methane'] == 'held' else 0.0
    emitted = made * (1.0 - cut)
    k = p['soil_ch4_per_yr']
    m0 = s['air_ch4']
    eq = emitted / dt / k
    m1 = eq + (m0 - eq) * math.exp(-k * dt)
    hold = HELD_METHANE_PPM * 1e-6 * st['air_mol'] if p['methane'] == 'held' else float('inf')
    removed = max(0.0, m1 - hold)
    m1 -= removed
    s['air_ch4'] = m1
    soil_ox = emitted - (m1 - m0) - removed
    f['methane'] = m1 - m0
    s['cum']['ch4_made'] += made
    s['cum']['ch4_cut'] += made * cut
    s['cum']['ch4_removed'] += removed
    s['cum']['ch4_soil'] += soil_ox
    air -= f['methane']
    # Weathering: the developed soils' factor, slowed in a free run by the cooling and the plants' weaker help as CO2
    # falls.
    cool = 0.0 if held else forcing(ppm, st['sweep']) / p['feedback_w_m2_k']
    plants_help = 1.0 if held else plant_weathering(ppm)
    tm = t0 + 0.5 * dt
    g = st['ground'][p['crust_case']]
    intensity = math.exp(DESSERT_B * cool) * plants_help
    w = st['weathering_eq_yr'] * weathering_factor(tm, p) * intensity * dt * at_year(g['years'], g['ratio'], tm)
    f['weathering'] = w
    sh = {k: at_year(g['years'], [x[k] for x in g['shares']], tm) for k in ('ca', 'mg', 'mono')}
    # The rock's calcite, carried from the springs study's central soils to the corner's, and the land water's carbon
    # into the filling crust, which follows the plants' production.
    x_rock = (at_year(g['years'], g['rock_mol_yr'], tm) * weathering_factor(tm, p)
              / weathering_factor(tm, CENTRAL_SOILS) * intensity * dt)
    span = span_of(g['spans'], tm)
    c_land = span['land_mol_c_yr'] * prod * dt
    s['rock_calcite'] += x_rock
    s['crust'] += c_land
    s['cum']['rock_ca_crust'] += x_rock
    s['cum']['crust_land'] += c_land
    f['rock_calcite'] = x_rock
    air -= x_rock + c_land
    pco2 = s['air_co2'] / st['air_mol'] * st['dry_pa'] / ATM_PA
    lake, sea = s['lake'], s['sea']
    lc, sc = st['lake_const'], sea_const(st, p)
    dic0 = lake['dic'] + sea['dic']
    for box, share in ((lake, st['lake_share']), (sea, 1.0 - st['lake_share'])):
        box['ta'] += w * share
        box['ca'] += 0.5 * w * share * sh['ca']
        box['mg'] += 0.5 * w * share * sh['mg']
        box['mono'] += w * share * sh['mono']
        box['dic'] += w * share
    s['cum']['rock_ca'] += 0.5 * w * sh['ca']
    s['cum']['rock_mg'] += 0.5 * w * sh['mg']
    s['cum']['rock_mono'] += w * sh['mono']
    air -= w
    # The magnesite route: the lakes' new magnesium made into MgCO3 where the rivers enter.
    x_mg = 0.0
    if managed and p['magnesium'] == 'magnesite':
        x_mg = MAGNESIUM_STEERED * 0.5 * w * st['lake_share'] * sh['mg']
        lake['mg'] -= x_mg
        lake['ta'] -= 2.0 * x_mg
        lake['dic'] -= x_mg
    # The lakes' outflow to the seas.
    share = min(1.0, st['discharge_m3_yr'] * 1e3 * st['lake_share'] * dt / lake['mass'])
    for key in ('ta', 'ca', 'mg', 'mono', 'dic'):
        moved = lake[key] * share
        lake[key] -= moved
        sea[key] += moved
    # The floors draw off the waters' dissolved load while the crust fills; delivered fresh water replaces the water.
    floors = 0.0
    for box, key in ((lake, 'lakes_kg_yr'), (sea, 'seas_kg_yr')):
        frac = min(1.0, span[key] * dt / box['mass'])
        for k in ('ta', 'ca', 'mg', 'mono', 'dic'):
            moved = box[k] * frac
            box[k] -= moved
            s['cum'][f'floor_{k}'] += moved
            if k == 'dic':
                floors += moved
        if box is sea:
            d = s['delivered']
            d['offset'] *= 1.0 - frac
            if p['seas'] == 'salted':                    # the delivered top-up is seawater, as the seas were
                top = frac * sea['mass']
                add = dict(ta=top * SEAWATER['alkalinity'], ca=top * SEAWATER['ca'], mg=top * SEAWATER['mg'])
                for k, v in add.items():
                    sea[k] += v
                d['ca'] += add['ca']
                d['mg'] += add['mg']
                d['offset'] += add['ta'] - 2.0 * add['ca'] - 2.0 * add['mg']
    s['crust'] += floors
    f['crust'] = c_land + floors
    # Calcite where the waters stand above their saturation, then CO2 equilibrium with the air.
    x_lake = precipitate(lake, pco2, lc, OMEGA_FRESH)
    x_sea = precipitate(sea, pco2, sc, s['omega_sea'])
    s['calcite_sed'] += x_lake + x_sea
    s['mgco3_sed'] += x_mg
    s['cum']['calcite'] += x_lake + x_sea
    s['cum']['calcite_lake'] += x_lake
    s['cum']['mgco3'] += x_mg
    f['calcite_lake'], f['calcite_sea'], f['mgco3'] = x_lake, x_sea, x_mg
    # Kilns (held runs) take the carbonate as it forms; the lime and magnesia are kept from the air.
    f['kiln_calcite'] = f['kiln_mg'] = 0.0
    if managed:
        f['kiln_calcite'], f['kiln_mg'] = x_lake + x_sea, x_mg
        s['calcite_sed'] -= x_lake + x_sea
        s['mgco3_sed'] -= x_mg
        s['cao'] += x_lake + x_sea
        s['mgo'] += x_mg
        s['cum']['kiln_calcite'] += x_lake + x_sea
        s['cum']['kiln_mg'] += x_mg
        air += x_lake + x_sea + x_mg
    s['air_co2'] += air
    before = lake['dic'] + sea['dic']
    if held:
        lake['dic'] = dic_mol(lake, pco2, lc)
        sea['dic'] = dic_mol(sea, pco2, sc)
        s['air_co2'] -= lake['dic'] + sea['dic'] - before            # the air-water exchange
        f['import'] = st['held_co2_mol'] - s['air_co2']
        s['air_co2'] = st['held_co2_mol']
        s['cum']['import'] += f['import']
    else:
        total = s['air_co2'] + before

        def mismatch(x):
            q = x * st['dry_pa'] / ATM_PA
            return x * st['air_mol'] + dic_mol(lake, q, lc) + dic_mol(sea, q, sc) - total
        x = brentq(mismatch, 1e-14, 0.05, xtol=1e-20, rtol=1e-12)
        q = x * st['dry_pa'] / ATM_PA
        lake['dic'] = dic_mol(lake, q, lc)
        sea['dic'] = dic_mol(sea, q, sc)
        s['air_co2'] = total - lake['dic'] - sea['dic']
        f['import'] = 0.0
    # The waters' dissolved growth: what the air gave them beyond the carbonate that formed.
    f['dissolved'] = lake['dic'] + sea['dic'] - dic0
    # Oxygen: one O2 for each organic carbon buried and kept, four FeO for each O2 the rock takes.
    s['cum']['o2_burial'] += f['lake_burial'] + f['sea_burial'] + f['peat'] - f['dredged']
    s['cum']['o2_rock'] += (w + 2.0 * x_rock) * st['o2_per_eq']
    s['cum']['dredged'] += f['dredged']
    s['cum']['lake_burial'] += f['lake_burial']
    s['cum']['sea_burial'] += f['sea_burial']
    s['cum']['weathering'] += w
    s['t'] = t0 + dt
    return f


def snapshot(s, st, p, f, dt):
    """The state at a year and the mean rates over the step that ended there."""
    pco2 = s['air_co2'] / st['air_mol'] * st['dry_pa'] / ATM_PA
    lake = speciate(s['lake'], pco2, st['lake_const'])
    sea = speciate(s['sea'], pco2, sea_const(st, p))
    water = lambda b, sp: dict(alkalinity_meq_kg=b['ta'] / b['mass'] * 1e3, ca_mmol_kg=b['ca'] / b['mass'] * 1e3,
                               mg_mmol_kg=b['mg'] / b['mass'] * 1e3, dic_mmol_kg=b['dic'] / b['mass'] * 1e3,
                               ph=sp['ph'], omega_calcite=sp['omega'])
    out = dict(
        year=s['t'], co2_ppm=s['air_co2'] / st['air_mol'] * 1e6, ch4_ppm=s['air_ch4'] / st['air_mol'] * 1e6,
        weathering_factor=weathering_factor(s['t'], p),
        stocks_gt_c=dict(air_co2=gt(s['air_co2']), air_ch4=gt(s['air_ch4']), land_plants=gt(sum(s['plants'].values())),
                         land_soils=gt(sum(s['soils'].values())), peat=gt(s['peat']),
                         soil_carbonate=gt(s['soil_carbonate']), sky=gt(s['sky']), people=gt(s['food'] + s['built']),
                         lake_dic=gt(s['lake']['dic']), sea_dic=gt(s['sea']['dic']),
                         organic_sediments=gt(s['lake_org'] + s['sea_org']),
                         carbonate_sediments=gt(s['calcite_sed'] + s['mgco3_sed']),
                         rock_calcite=gt(s['rock_calcite']), filling_crust=gt(s['crust'])),
        waters=dict(lakes=water(s['lake'], lake), seas=water(s['sea'], sea)),
        cumulative_gt_c=dict(import_=gt(s['cum']['import']),
                             kilns=gt(s['cum']['kiln_calcite'] + s['cum']['kiln_mg']), dredged=gt(s['cum']['dredged']),
                             carbonate_formed=gt(s['cum']['calcite'] + s['cum']['mgco3'])),
        lime_gt=dict(cao=s['cao'] * M_CAO / 1e15, mgo=s['mgo'] * M_MGO / 1e15),
        oxygen=oxygen(s, st))
    if f is not None:
        r = {k: gt(v) / dt for k, v in f.items()}
        r['net_sink'] = (r['land'] + r['peat'] + r['soil_carbonate'] + r['sky'] + r['lake_burial'] + r['sea_burial']
                         + r['methane'] + r['dissolved'] + r['calcite_lake'] + r['calcite_sea'] + r['mgco3']
                         + r['rock_calcite'] + r['crust'])
        r['kilns'] = r['kiln_calcite'] + r['kiln_mg']
        out['rates_gt_c_yr'] = r
        out['kilns'] = kiln_power(f['kiln_calcite'] / dt, f['kiln_mg'] / dt, st)
    return out


def oxygen(s, st):
    """O2 the carbon cycle leaves in the air (mol, then shares): burial net of the dredged return, the biosphere's
    stocks if made from imported CO2, methane's two O2 a carbon, and the weathering rock's FeO."""
    burial, rock = s['cum']['o2_burial'], -s['cum']['o2_rock']
    stocks = organic(s) - s['initial_organic']
    methane = 2.0 * s['air_ch4']
    net = burial + rock
    return dict(burial_kg=burial * M_O2 * 1e-3, rock_kg=rock * M_O2 * 1e-3,
                burial_less_rock_share=net / st['o2_mol'], burial_less_rock_kpa=net / st['o2_mol'] * st['o2_pa'] / 1e3,
                biosphere_from_co2_kg=stocks * M_O2 * 1e-3, methane_kg=methane * M_O2 * 1e-3,
                all_share=(net + stocks + methane) / st['o2_mol'])


def kiln_power(calcite_mol_yr, mg_mol_yr, st):
    """GW to calcine the carbonate as it forms: at the minimum and in the best lime kilns (the magnesite's best
    scaled by the lime kilns' ratio, a design guess); the heat on the Moon if the energy is brought in."""
    t_co2 = lambda mol: mol * M_CO2 / 1e6                                     # tonnes CO2 a year
    best = [g / BEST_KILN_GJ_T_LIME[0] for g in BEST_KILN_GJ_T_LIME]          # the best kilns over the minimum
    gw = lambda tco2, gj: tco2 * gj * 1e9 / YEAR_S / 1e9
    theory = gw(t_co2(calcite_mol_yr), CALCITE_GJ_T_CO2) + gw(t_co2(mg_mol_yr), MAGNESITE_GJ_T_CO2)
    kilns = [gw(t_co2(calcite_mol_yr), CALCITE_GJ_T_CO2 * b) + gw(t_co2(mg_mol_yr), MAGNESITE_GJ_T_CO2 * b)
             for b in best]
    as_calcite = gw(t_co2(calcite_mol_yr + mg_mol_yr), CALCITE_GJ_T_CO2)
    return dict(gt_co2_yr=dict(calcite=t_co2(calcite_mol_yr) / 1e9, magnesite=t_co2(mg_mol_yr) / 1e9),
                gw_minimum=theory, gw_best_kilns=kilns, gw_if_all_calcite=as_calcite,
                heat_w_m2_if_brought_in=[k / 1e3 * st['heat_w_m2_per_tw'] for k in kilns],
                share_of_sunlight=[k * 1e9 / (st['sunlight_w_m2'] * MOON_AREA) for k in kilns],
                lime_gt_yr=dict(cao=calcite_mol_yr * M_CAO / 1e15, mgo=mg_mol_yr * M_MGO / 1e15))


def checks_of(s, st):
    """Carbon and the cations, against their sources, as shares of the stock they balance."""
    c = carbon_total(s) - s['initial_total'] - s['cum']['import']
    d = s['delivered']
    ca = (d['ca'] + s['cum']['rock_ca'] + s['cum']['rock_ca_dry'] + s['cum']['rock_ca_crust']
          - (s['lake']['ca'] + s['sea']['ca'] + s['cum']['calcite'] + s['soil_carbonate'] + s['rock_calcite']
             + s['cum']['floor_ca']))
    mg = d['mg'] + s['cum']['rock_mg'] - (s['lake']['mg'] + s['sea']['mg'] + s['cum']['mgco3'] + s['cum']['floor_mg'])
    mono = d['mono'] + s['cum']['rock_mono'] - (s['lake']['mono'] + s['sea']['mono'] + s['cum']['floor_mono'])
    ta = (s['lake']['ta'] + s['sea']['ta'] - 2 * (s['lake']['ca'] + s['sea']['ca'] + s['lake']['mg'] + s['sea']['mg'])
          - s['lake']['mono'] - s['sea']['mono'] - d['offset'])
    lime = s['cao'] - s['cum']['kiln_calcite'] + s['mgo'] - s['cum']['kiln_mg']
    scale_c = carbon_total(s)
    scale_ca = d['ca'] + s['cum']['rock_ca'] + s['cum']['rock_ca_dry'] + s['cum']['rock_ca_crust'] or 1.0
    scale_mg = d['mg'] + s['cum']['rock_mg'] or 1.0
    scale_mono = d['mono'] + s['cum']['rock_mono'] or 1.0
    scale_ta = abs(s['lake']['ta']) + abs(s['sea']['ta']) or 1.0
    return dict(carbon=c / scale_c, calcium=ca / scale_ca, magnesium=mg / scale_mg, sodium_potassium=mono / scale_mono,
                alkalinity=ta / scale_ta, lime=lime / (s['cao'] + s['mgo'] or 1.0))


def simulate(st, p, mode='held', start=None, step_scale=1.0):
    """Integrate a run. Held runs keep CO2 at 400 ppm with the managed return and an import; free runs start from a held
    state and let CO2 follow the sinks with the managed return and the import stopped."""
    held = mode == 'held'
    s = copy.deepcopy(start) if start is not None else initial_state(st, p)
    t_start = s['t']
    out, rates_for_start = {}, None
    snaps_at = set(OUT_YEARS) | {t_start}
    mature = mature_land(p)
    grown = None
    crossings, calcite_from, peak = {}, {}, (0.0, 0.0)
    worst = dict(carbon=0.0, calcium=0.0, magnesium=0.0, sodium_potassium=0.0, alkalinity=0.0, lime=0.0)
    prev_ppm = s['air_co2'] / st['air_mol'] * 1e6
    for t0, dt in steps(t_start, step_scale):
        f = step(s, st, p, dt, held, managed=held)
        if rates_for_start is None:
            rates_for_start = (copy.deepcopy(f), dt)
        for k, v in checks_of(s, st).items():
            worst[k] = max(worst[k], abs(v))
        t1 = round(s['t'], 6)
        for box in ('lake', 'sea'):
            if box not in calcite_from and f[f'calcite_{box}'] > 0.0:
                calcite_from[box] = t0
        if f['import'] / dt > peak[0]:
            peak = (f['import'] / dt, t0 + 0.5 * dt)
        if not held:
            ppm = s['air_co2'] / st['air_mol'] * 1e6
            for name, pa in THRESHOLDS_PA.items():
                limit = pa / st['dry_pa'] * 1e6
                if name not in crossings and ppm <= limit < prev_ppm:
                    crossings[name] = t0 + dt * (prev_ppm - limit) / (prev_ppm - ppm) - t_start
            prev_ppm = ppm
        if any(abs(t1 - y) < 1e-6 for y in snaps_at):
            out[int(round(t1))] = snapshot(s, st, p, f, dt)
        if held and grown is None and land_organic(s) >= GROWN_SHARE * mature:
            grown = (t1, copy.deepcopy(s))
    first = int(round(t_start))                     # the first year carries the first step's rates
    out[first] = snapshot(copy.deepcopy(start) if start is not None else initial_state(st, p), st, p,
                          rates_for_start[0], rates_for_start[1])
    return dict(series=dict(sorted(out.items())), grown=grown, crossings=crossings, checks=worst, final=s,
                calcite_from=calcite_from, peak_import=dict(gt_c_yr=gt(peak[0]), year=peak[1]))


def columns(series):
    """A run's snapshots as columns: every leaf a list over the years."""
    snaps = list(series.values())

    def merge(items):
        if isinstance(items[0], dict):
            return {k: merge([x[k] for x in items]) for k in items[0]}
        return list(items)
    return merge(snaps)


# --- The runs and what they answer -------------------------------------------------------------------------------

HELD_RUNS = [(c, m, g) for c in CORNERS for m in ('held', 'soils') for g in ('natural', 'magnesite')]
VARIANTS = {'central/held/natural/salted': dict(seas='salted'),
            'central/held/magnesite/salted': dict(seas='salted'),
            'central/held/natural/slow_soils': dict(soils='slow')}
FREE_RUNS = [f'{c}/held/natural' for c in CORNERS] + ['central/soils/natural']


def name_of(corner, methane, magnesium):
    return f'{corner}/{methane}/{magnesium}'


def parse(name):
    corner, methane, magnesium = name.split('/')[:3]
    return dict(corner=corner, methane=methane, magnesium=magnesium)


def first_year(series, test):
    for y, sn in series.items():
        if test(sn):
            return y
    return None


def held_summary(run, st, p):
    """A held run's answers: the build's carbon when the biosphere is grown, the managed return and the import."""
    t_g, s_g = run['grown']
    grown = snapshot(s_g, st, p, None, 1.0)
    final = run['final']
    at = {}
    for y in REPORT_YEARS:
        r, k = run['series'][y]['rates_gt_c_yr'], run['series'][y]['kilns']
        at[str(y)] = dict(
            net_sink=r['net_sink'], kilns=r['kilns'], dredged=r['dredged'], import_=r['import'],
            import_kg_c_s=r['import'] * 1e12 / YEAR_S,
            managed_share_of_net_sink=(r['kilns'] + r['dredged']) / r['net_sink'],
            kiln_gt_co2_yr=k['gt_co2_yr'], kiln_gw_minimum=k['gw_minimum'], kiln_gw_best_kilns=k['gw_best_kilns'],
            lime_gt_yr=k['lime_gt_yr'],
            sinks=dict(land=r['land'], peat=r['peat'], soil_carbonate=r['soil_carbonate'], sky=r['sky'],
                       lake_burial=r['lake_burial'], sea_burial=r['sea_burial'], methane=r['methane'],
                       waters_dissolved=r['dissolved'], carbonate=r['calcite_lake'] + r['calcite_sea'] + r['mgco3'],
                       rock_calcite=r['rock_calcite'], filling_crust=r['crust']))
    total_g = gt(carbon_total(s_g))
    w = grown['stocks_gt_c']
    filled = sum(w[k] for k in ('air_co2', 'air_ch4', 'land_plants', 'land_soils', 'sky', 'people'))
    burial_10k = run['series'][10000]['rates_gt_c_yr']['lake_burial']
    sediment = burial_10k * 1e12 / p['delta_oc']                               # kg of delta sediment a year
    cum, cum_g = final['cum'], s_g['cum']
    return dict(
        grown_year=t_g,
        build=dict(carbon_on_the_moon=total_g, at_seeding=gt(s_g['initial_total']),
                   imported_after_seeding=gt(cum_g['import']), where=w, beyond_the_air=total_g - w['air_co2'],
                   stocks_filled=filled, taken_by_the_sinks=total_g - filled),
        peak_import=run['peak_import'],
        after_grown=dict(imported=gt(cum['import'] - cum_g['import']),
                         kilns=gt(cum['kiln_calcite'] + cum['kiln_mg'] - cum_g['kiln_calcite'] - cum_g['kiln_mg']),
                         dredged=gt(cum['dredged'] - cum_g['dredged'])),
        methane_after_grown=methane_rates(final, s_g, st, p),
        at=at, oxygen={str(y): run['series'][y]['oxygen'] for y in (1000, 5000, 10000)},
        calcite_from_year=run['calcite_from'],
        carbonate_formed_gt_c=dict(lakes=gt(cum['calcite_lake']), seas=gt(cum['calcite'] - cum['calcite_lake']),
                                   magnesite_route=gt(cum['mgco3'])),
        lake_burial_returned_by_dredging_at_10000=dict(
            sediment_gt_yr=sediment / 1e12, fuel_gw=sediment * p['dredge_mj_per_kg'] * 1e6 / YEAR_S / 1e9,
            basis='delta sediment at its organic content, enough to return all of the lakes\' burial'),
        seas_saturated_year=first_year(run['series'], lambda sn: sn['waters']['seas']['omega_calcite'] >= 1.0),
        lakes_saturated_year=first_year(run['series'], lambda sn: sn['waters']['lakes']['omega_calcite'] >= 1.0))


def methane_rates(final, grown, st, p):
    """Methane's mean flows once the biosphere is grown (Tg CH4 a year), its stock's growth, and the share of the air
    the contactors pass a year at full capture when the level is held."""
    years = final['t'] - grown['t']
    tg = lambda key: (final['cum'][key] - grown['cum'][key]) * M_CH4 / 1e12 / years
    out = dict(made=tg('ch4_made'), cut_at_the_sources=tg('ch4_cut'), soils=tg('ch4_soil'), removed=tg('ch4_removed'),
               stock_growth=(final['air_ch4'] - grown['air_ch4']) * M_CH4 / 1e12 / years)
    if p['methane'] == 'held':
        out['air_share_passed_a_year'] = out['removed'] / (st['tg_ch4_per_ppm'] * HELD_METHANE_PPM)
    return out


def free_summary(run, start_year):
    series = run['series']
    first, final = series[min(series)], series[max(series)]
    land = lambda sn: sn['stocks_gt_c']['land_plants'] + sn['stocks_gt_c']['land_soils']
    return dict(start_year=start_year, years_to=run['crossings'],
                co2_ppm={str(y): sn['co2_ppm'] for y, sn in series.items()},
                land_released_gt_c=land(first) - land(final),
                net_sink_gt_c_yr=dict(start=first['rates_gt_c_yr']['net_sink'],
                                      end=final['rates_gt_c_yr']['net_sink']))


def exchange_check(st, s):
    """Years for a box to come within 1/e of CO2 equilibrium with the air at PISTON_M_DAY: its mean depth over the
    piston velocity, times the DIC it takes up for each unit of dissolved CO2 at fixed alkalinity."""
    out = {}
    for name, box, c, area in (('lakes', s['lake'], st['lake_const'], st['lake_m2']),
                               ('seas', s['sea'], st['fresh_sea_const'], st['sea_m2'])):
        q = st['pco2_atm']
        d0, d1 = speciate(box, q, c), speciate(box, q * 1.01, c)
        buffer = (d1['dic'] - d0['dic']) / (d1['co2'] - d0['co2'])
        depth = box['mass'] / 1000.0 / area
        out[name] = dict(mean_depth_m=depth, dic_per_co2=buffer,
                         years=depth * buffer / (PISTON_M_DAY * JULIAN_YEAR_DAYS))
    return out


def law_check(data):
    """The weathering law against the CO2 study's own cases."""
    cases = data[CO2_STUDY]['weathering']['cases']
    worst = 0.0
    for c in cases:
        mine = c['runoff_mm_yr'] * DESSERT_A * math.exp(DESSERT_B * c['temperature_c']) * c['fresh_glass_factor']
        worst = max(worst, abs(mine - c['gross_mol_km2_yr']) / c['gross_mol_km2_yr'])
    return dict(cases=len(cases), worst_relative_difference=worst)


def convergence(st, params, runs):
    """The central held run and its free run again at half the step: the largest relative change in what the answers
    read, and the change in the free run's threshold years."""
    name = name_of('central', 'held', 'natural')
    fine = simulate(st, params[name], 'held', step_scale=0.5)
    fine_free = simulate(st, params[name], 'free', start=fine['grown'][1], step_scale=0.5)
    worst = 0.0
    for y in REPORT_YEARS:
        a, b = runs[name]['series'][y], fine['series'][y]
        for key in ('import', 'net_sink', 'kilns', 'dissolved', 'lake_burial'):
            x, z = a['rates_gt_c_yr'][key], b['rates_gt_c_yr'][key]
            if max(abs(x), abs(z)) > 1e-6:
                worst = max(worst, abs(x - z) / max(abs(x), abs(z)))
        for key in ('sea_dic', 'organic_sediments', 'land_soils'):
            x, z = a['stocks_gt_c'][key], b['stocks_gt_c'][key]
            worst = max(worst, abs(x - z) / max(abs(x), abs(z)))
    coarse = runs[name + '/free']['crossings']
    return dict(run=name, held_worst_relative_change=worst,
                free_threshold_years_change={k: abs(coarse[k] - fine_free['crossings'][k]) for k in coarse},
                grown_year=dict(step=runs[name]['grown'][0], half_step=fine['grown'][0]))


def comparisons(data, st, held):
    """This model's numbers beside the records they revise."""
    industry = data[INDUSTRY]['after_build']['carbon']
    co2 = data[CO2_STUDY]
    case = [c for c in co2['weathering']['cases'] if c['runoff_mm_yr'] == 349 and c['temperature_c'] == 26.0
            and c['fresh_glass_factor'] == 1.0][0]
    gross = gt(st['weathering_eq_yr']) * M_CO2 / M_C
    lakes = [LAKE_BURIAL_G_C_M2_YR[0] * st['lake_m2'] / 1e15, LAKE_BURIAL_G_C_M2_YR[2] * st['lake_m2'] / 1e15]
    where = {c: held[name_of(c, 'held', 'natural')]['build']['where'] for c in CORNERS}
    steady = {c: held[name_of(c, 'held', 'natural')]['at']['10000']['import_kg_c_s'] for c in CORNERS}
    beyond = data[SCREEN]['carbon_beyond_the_biosphere']
    return dict(
        co2_study=dict(gross_weathering_gt_co2_yr=dict(its_case=case['gross_gt_yr'], design_catchment_rate=gross,
                                                       ratio=gross / case['gross_gt_yr']),
                       earth_like_biosphere_kg_co2=co2['biosphere']['earth_like_land_carbon_kg_co2'],
                       air_kg_co2_constant_gravity=co2['design_air']['inventory_kg_co2'],
                       this_land_and_peat_at_grown_gt_c={c: w['land_plants'] + w['land_soils'] + w['peat']
                                                         for c, w in where.items()}),
        resources_industry=dict(
            burial_gt_co2_yr=industry['burial_gt_co2_yr'], its_basis=industry['burial_basis'],
            this_sea_burial_gt_co2_yr=gt(st['sea_burial_mol_yr']) * M_CO2 / M_C,
            this_lake_burial_gt_co2_yr=[x * M_CO2 / M_C for x in lakes],
            kilns_gw=industry['kilns_gw'], burial_import_kg_c_s=industry['burial_import_as_carbon_kg_s'],
            this_steady_import_kg_c_s_at_10000=steady),
        resources_screen=dict(seas_and_lakes_dic_kg_co2_at_earths_dic=beyond['seas_and_lakes_dic_co2_kg'],
                              air_co2_kg=data[SCREEN]['inventories']['co2_kg']))


def assumptions(data, st):
    m = data[LEDGERS]['gases']['methane']
    sea_p = data[SEA_P]['routes']['dredging']
    rows = [
        ('co2_ppm', DESIGN_CO2_PPM, 'design choice',
         'held at the level the planned climate was run with: 48.6 Pa in the design\'s 121.6 kPa of dry air'),
        ('methane_held_ppm', HELD_METHANE_PPM, 'design choice',
         'the author\'s choice of 13 October, Earth\'s present methane, the balanced middle of three pairs weighed '
         '(domain/ecology, domain synthesis, gases.held); 90% of the sources cut at the sources, the rest oxidized '
         'by the soils and by contactors in open air'),
        ('methane_soils_case_ppm', [c['steady_ppm'] for c in m['central']], 'project product',
         'no cut and no removal: the gas ledger\'s central sources at the blend\'s two ends over the soils\' central '
         'sink; the central corner takes their mean'),
        ('temperatures_k', dict(global_surface=DESIGN_SURFACE_K, land_air=LAND_AIR_K, sea_surface=SEA_SURFACE_K),
         'project result', 'the corrected design case A28_dim5_moon, years 20-29 (climate/gcm/README.md); '
         'weathering and the lakes at the air over land, the seas at the sea surface'),
        ('runoff_mm_yr', st['runoff_mm_yr'], 'project product',
         'the drainage product\'s discharge to the seas over the atlas\'s land, lakes included, as '
         'resources/industry.py runoff_mm() reads it'),
        ('weathering_law', dict(a=DESSERT_A, b=DESSERT_B), 'Earth measurement',
         'Dessert et al. (2003) for basaltic catchments, CO2 consumed = runoff x a exp(b T); the highlands weather '
         'like basalt through their anorthite (resources finding 5)'),
        ('groundwater', dict(springs=SPRINGS, stage_1=STAGE1, crusts=CRUST_CASES), 'project product',
         'the springs study\'s share of the catchment law\'s weathering that reaches the waters, its cations and the '
         'calcite the rock keeps, at 1,000 and 10,000 years and a million (the 1,000-year values held before), linear '
         'in the logarithm of the year; stage 1\'s uptake by the floors and the land by span, the land\'s water with '
         'the springs study\'s carbon at the water table; the waters\' load drawn off with the floors\' water, which '
         'delivered fresh water replaces'),
        ('fresh_soil_factor', list(FRESH_FACTOR), 'design guess',
         'the developed soils start at this multiple of the catchment rate, within the CO2 study\'s 1-10 for fine '
         'glassy ground'),
        ('soil_age_at_start_yr', list(SOIL_AGE_AT_START_YR), 'design guess',
         'the soils\' age when the biosphere is seeded, made in the build\'s last centuries; they age as '
         '(1 + t/age)^-0.61 toward the catchment rate'),
        ('age_exponent', AGE_EXPONENT, 'Earth measurement', 'field rates fall as t^-0.61 (White & Brantley 2003)'),
        ('cation_shares', st['cation_shares'], 'lunar measurement',
         'the highland soil\'s Ca, Mg, Na and K (the lunar-cycle register, from the Lunar Sourcebook), released in '
         'proportion; its 5.1 wt% FeO oxidized as the rock weathers'),
        ('biomes_gt_c', BIOMES, 'Earth inventory',
         'plants and the top metre of soil by biome (IPCC TAR WG1 Table 3.2, WBGU column): area 10^6 km2, plants, '
         'soils'),
        ('rain_regions', dict(forest_rain_mm_day=FOREST_RAIN_MM_DAY), 'design mapping',
         'wet land in bands of at least 2 mm a day holds tropical forest, the rest of the wet land savanna, the fog '
         'desert and polar dry land desert; the kept wild and the parks hold their region\'s wild stock'),
        ('land_uses', 'C\' with the blend', 'project product with a design mapping',
         'the synthesis ledgers: fields and settlements on the forest region at cropland stocks, heritage pasture on '
         'the savanna region (the middle-wet homestead belt) at savanna stocks; the low corner takes the blend\'s '
         'high-cost end with 2 kg of meat a person, the high corner the low-cost end with 0.5 kg'),
        ('plant_growth_years', dict(wild=WILD_GROWTH_YEARS, managed=MANAGED_GROWTH_YEARS),
         'Earth measurement and design guess',
         'logistic growth from a 1% seed to 90% of the stock: forest regrowth in about 50 years (the CO2 study, '
         'Poorter et al. 2016), five years for fields, settlements and pasture'),
        ('soil_build_kg_c_m2_yr', SOIL_BUILD_KG_C_M2_YR, 'Earth measurement',
         'Post & Kwon (2000) on establishment: the soils approach their stock on the forest\'s time, 372 years, in '
         'every class; the slow-soil variant takes five times as long'),
        ('starting_stocks', dict(plants=SEED_SHARE, soils=SOIL_START_SHARE, sky=SEED_SHARE), 'design guess',
         'shares of the mature stocks at the seeding'),
        ('wet_margins_share_of_moon', data[ECO_SCREEN]['trace_gases']['methane']['wetland_share_scenario'],
         'project product', 'the ecology screen\'s wetland scenario, central 2.6%, on the forest region'),
        ('peat_g_c_m2_yr', list(PEAT_G_C_M2_YR), 'Earth measurement',
         'Holocene peat accumulation (Loisel et al. 2014), following the wet margins\' plants'),
        ('soil_carbonate_g_c_m2_yr', list(SOIL_CARBONATE_G_C_M2_YR), 'Earth measurement',
         'pedogenic carbonate in Mojave Aridisols (Schlesinger 1985) over the fog desert and polar dry land; its '
         'calcium comes from silicate, so each carbon is a sink'),
        ('lake_burial_g_c_m2_yr', list(LAKE_BURIAL_G_C_M2_YR), 'Earth measurement',
         'lakes and reservoirs (Mendonca et al. 2017) over the rain-fed lakes, following the land\'s plants'),
        ('sea_burial_g_c_m2_yr', dict(shelf=SHELF_G_C_M2_YR, deep=DEEP_G_C_M2_YR), 'Earth measurement',
         'Earth\'s unvegetated shelves and deep sea (Duarte et al. 2005) on the seas\' 14% under 200 m and the rest; '
         'the rivers bring the seas little phosphorus, so these are upper bounds'),
        ('dredging', dict(p_need_tg_yr=data[LEDGERS]['phosphorus']['lost'],
                          sediment_kg_per_kg_p=sea_p['sediment_kg_per_kg_p'], organic_carbon=list(DELTA_OC)),
         'project product and Earth measurement',
         'the people\'s phosphorus losses dredged back from the lake deltas (the synthesis ledger, resources\' '
         'routes), the sediment at Earth\'s deltaic and shelf organic content (Hedges & Keil 1995, Berner 1982) and '
         'respired on the farms'),
        ('waters', 'well mixed, in CO2 equilibrium with the air', 'design simplification',
         'checks.exchange gives the time to equilibrium'),
        ('lake_passage_share', st['lake_share'], 'project product',
         'the share of the land\'s runoff that crosses a rain-fed lake (the routing); the rest reaches the seas '
         'directly'),
        ('calcite_saturation_held', OMEGA_FRESH, 'design guess',
         'calcifiers precipitate calcite wherever fresh water stands above three times saturation; salted seas '
         'keep the saturation they arrive with'),
        ('salted_seas', dict(salinity=SEA_SALINITY, **SEAWATER), 'textbook composition',
         'Earth\'s seawater delivered with the water'),
        ('carbonate_constants', 'at the design temperatures', 'textbook constant',
         'Plummer & Busenberg (1982), Harned & Owen and Davies in fresh water; Weiss (1974), Lueker et al. (2000), '
         'Dickson (1990), Uppstrom (1974), Millero (1995) and Mucci (1983) in seawater; no ion pairs'),
        ('sky', 'the sky ecology', 'project product',
         'cover 0.1%, 1% and 2%, its class shares and its bodies\' standing host carbon at the reference light; '
         'logistic growth at the class-weighted growth rates'),
        ('people', dict(billion=data[LEDGERS]['choices']['moon_billion'], food_residence_yr=FOOD_RESIDENCE_YR,
                        body_kg=BODY_KG, body_carbon=BODY_CARBON, built_t_c=list(BUILT_T_C_PER_PERSON)),
         'design guess',
         'C\'\'s 8 billion from the seeding; the diet\'s food grown with its losses (the food study) held for half '
         'a year; the built stocks\' wood and polymers'),
        ('kilns', 'all carbonate calcined as it forms', 'requirement',
         'the managed return that keeps the carbonate sediments from growing; the lime and magnesia kept from the '
         'air, sealed or built in; a share f that takes CO2 back raises the load by 1/(1 - f)'),
        ('kiln_energy_gj_t_co2', dict(calcite=CALCITE_GJ_T_CO2, magnesite=MAGNESITE_GJ_T_CO2,
                                      best_kilns_gj_t_lime=list(BEST_KILN_GJ_T_LIME)),
         'Earth measurement and design guess',
         'calcination at its minimum and the best lime kilns (Schorcht et al. 2013); magnesite\'s kilns scaled by '
         'the lime kilns\' ratio'),
        ('magnesite_route', MAGNESIUM_STEERED, 'design lever',
         'the lakes\' new magnesium made into MgCO3 where the rivers enter, then calcined; the route\'s own energy '
         'is outside the model'),
        ('free_run', 'from the grown biosphere', 'design simplification',
         'kilns, the magnesite route and the import stop; dredging and the methane hold go on; weathering cools with '
         'the forcing over the GCM\'s feedback and loses the plants\' help as CO2 falls; the land\'s and the sky\'s '
         'stocks follow the leaf rate'),
        ('feedback_w_m2_k', list(FEEDBACK_W_M2_K), 'project result',
         'the GCM\'s 1.0-1.6 W m-2 K-1 as the CO2 study gives it; the forcing from the radiative-convective CO2 '
         'sweep'),
        ('fert', FERT, 'literature model parameter', 'GEOCARB III, eq 5 (Berner & Kothavala 2001)'),
    ]
    return [dict(name=n, value=v, evidence=e, basis=b) for n, v, e, b in rows]


def groundwater_setting(st):
    """The groundwater's terms as the runs take them, for each of stage 1's crusts."""
    out = {}
    for case, g in sorted(st['ground'].items()):
        out[case] = dict(
            corners=[c for c, v in CRUST_CASES.items() if v == case],
            reaching_the_waters={f'{y:g}': r for y, r in zip(g['years'], g['ratio'])},
            cation_shares={f'{y:g}': sh for y, sh in zip(g['years'], g['shares'])},
            rock_calcite_gt_c_yr={f'{y:g}': gt(m) for y, m in zip(g['years'], g['rock_mol_yr'])},
            filling_crust={f'{r["start"]:g}-{r["end"]:g}': dict(
                seas_km3_yr=r['seas_kg_yr'] / 1e12, lakes_km3_yr=r['lakes_kg_yr'] / 1e12,
                land_km3_yr=r['land_kg_yr'] / 1e12, land_water_gt_c_yr=gt(r['land_mol_c_yr'])) for r in g['spans']})
    return out


def corner_summary(p):
    c = p['classes']
    return dict(
        land=dict(classes_km2={k: v['area_m2'] / 1e6 for k, v in c.items()},
                  plants_gt_c={k: gt(v['plants_mol']) for k, v in c.items()},
                  soils_gt_c={k: gt(v['soils_mol']) for k, v in c.items()},
                  regions_km2={k: v / 1e6 for k, v in p['regions_m2'].items()}),
        fresh_factor=p['fresh_factor'], soil_age_yr=p['soil_age_yr'], soil_tau_yr=p['soil_tau_yr'],
        lake_burial_gt_c_yr=gt(p['lake_burial_mol_yr']), peat_gt_c_yr=gt(p['peat_mol_yr']),
        soil_carbonate_gt_c_yr=gt(p['soil_carbonate_mol_yr']), dredged_gt_c_yr=gt(p['dredged_mol_yr']),
        dredged_sediment_gt_yr=p['dredged_sediment_kg_yr'] / 1e12, p_need_tg_yr=p['p_need_kg_yr'] / 1e9,
        methane_sources_tg_yr=dict(wild=p['wild_ch4_mol_yr'] * M_CH4 / 1e12,
                                   people=p['people_ch4_mol_yr'] * M_CH4 / 1e12),
        sky=dict(cover=p['sky']['cover'], standing_kg_c_m2=p['sky']['standing_kg_c_m2'],
                 growth_rate_per_year=p['sky']['growth_rate_per_year'], stock_gt_c=gt(p['sky']['stock_mol'])),
        people_gt_c=dict(food_and_bodies=gt(p['food_mol']), built=gt(p['built_mol']),
                         food_grown_gt_c_yr=gt(p['food_grown_mol_yr']),
                         power_made_food_gt_c_yr=gt(p['food_grown_mol_yr']) * p['power_made_share']),
        feedback_w_m2_k=p['feedback_w_m2_k'])


def questions(held, free):
    """The four answers, gathered from the runs' summaries (methane held at 1.9 ppm unless named)."""
    holds = {c: dict(start_year=free[f'{c}/held/natural']['start_year'],
                     years_to=free[f'{c}/held/natural']['years_to'],
                     co2_ppm_at_10000=free[f'{c}/held/natural']['co2_ppm']['10000'],
                     land_released_gt_c=free[f'{c}/held/natural']['land_released_gt_c']) for c in CORNERS}
    holds['central_with_methane_at_the_soils_steady_state'] = dict(years_to=free['central/soils/natural']['years_to'])
    ret = {}
    for c in CORNERS:
        for g in ('natural', 'magnesite'):
            ret[f'{c}/{g}'] = {y: dict(net_sink=v['net_sink'], kilns=v['kilns'], dredged=v['dredged'],
                                       import_=v['import_'], kiln_gt_co2_yr=v['kiln_gt_co2_yr'],
                                       kiln_gw_minimum=v['kiln_gw_minimum'],
                                       kiln_gw_best_kilns=v['kiln_gw_best_kilns'], lime_gt_yr=v['lime_gt_yr'])
                               for y, v in held[name_of(c, 'held', g)]['at'].items()}
    imports = {name: dict(grown_year=h['grown_year'], carbon_on_the_moon=h['build']['carbon_on_the_moon'],
                          beyond_the_air=h['build']['beyond_the_air'], peak_import=h['peak_import'],
                          import_gt_c_yr={y: h['at'][y]['import_'] for y in ('2000', '5000', '10000')},
                          imported_after_grown=h['after_grown']['imported'])
               for name, h in held.items()}
    return dict(holds=holds,
                managed_return=dict(by_corner=ret,
                                    kiln_load_if_lime_recarbonates={f'{f:g}': 1.0 / (1.0 - f)
                                                                    for f in LIME_RECARBONATING}),
                imports=imports, oxygen={c: held[name_of(c, 'held', 'natural')]['oxygen'] for c in CORNERS})


def results(step_scale=1.0):
    data = {p: (load(p, s) if not p.endswith('.csv') else None) for p, s in INPUTS.items()}
    st = setting(data)
    runs, params, held = {}, {}, {}
    for corner, methane, magnesium in HELD_RUNS:
        name = name_of(corner, methane, magnesium)
        params[name] = parameters(data, st, corner, methane=methane, magnesium=magnesium)
        runs[name] = simulate(st, params[name], 'held', step_scale=step_scale)
    for name, change in VARIANTS.items():
        base = parse(name)
        params[name] = parameters(data, st, base['corner'], methane=base['methane'], magnesium=base['magnesium'],
                                  **change)
        runs[name] = simulate(st, params[name], 'held', step_scale=step_scale)
    for name in list(runs):
        held[name] = held_summary(runs[name], st, params[name])
    free = {}
    for name in FREE_RUNS:
        run = simulate(st, params[name], 'free', start=runs[name]['grown'][1], step_scale=step_scale)
        runs[name + '/free'] = run
        free[name] = free_summary(run, runs[name]['grown'][0])
    central = name_of('central', 'held', 'natural')
    worst = {k: max(r['checks'][k] for r in runs.values()) for k in runs[central]['checks']}
    return dict(
        schema=SCHEMA,
        producer=dict(files={f: sha256(ROOT / f) for f in FILES}, inputs={p: sha256(ROOT / p) for p in INPUTS},
                      constants=constants_used(FILES)),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        assumptions=assumptions(data, st),
        setting=dict(air_mol=st['air_mol'], co2_gt_c_at_400_ppm=gt(st['held_co2_mol']),
                     o2_kg=st['o2_mol'] * M_O2 * 1e-3, runoff_mm_yr=st['runoff_mm_yr'],
                     weathering_gt_c_yr_at_catchment_rate=gt(st['weathering_eq_yr']),
                     cation_shares=st['cation_shares'], groundwater=groundwater_setting(st),
                     lakes=dict(area_m2=st['lake_m2'], water_kg=st['lake_kg']),
                     seas=dict(area_m2=st['sea_m2'], water_kg=st['sea_kg'], shelf_share=st['shelf_share']),
                     sea_burial_gt_c_yr=gt(st['sea_burial_mol_yr']), tg_ch4_per_ppm=st['tg_ch4_per_ppm'],
                     forcing_w_m2=st['sweep'], leaf_curve=st['leaf_curve'],
                     built_stocks_renewal_gt_yr=st['built_renewal_gt_yr'],
                     small_common_level_waters=dict(
                         area_m2=st['small_waters_m2'],
                         at_lake_rates_gt_c_yr=[LAKE_BURIAL_G_C_M2_YR[0] * st['small_waters_m2'] / 1e15,
                                                LAKE_BURIAL_G_C_M2_YR[2] * st['small_waters_m2'] / 1e15])),
        corners={c: corner_summary(params[name_of(c, 'held', 'natural')]) for c in CORNERS},
        questions=questions(held, free),
        held=held, free=free,
        runs={name: dict(series=columns(r['series']), crossings=r['crossings']) for name, r in runs.items()},
        checks=dict(conservation_worst_share=worst, law=law_check(data),
                    exchange=exchange_check(st, runs[central]['final']),
                    convergence=convergence(st, params, runs) if step_scale == 1.0 else None),
        comparisons=comparisons(data, st, held))


def to_json(x, depth=0):
    """JSON with one key a line and lists of numbers kept on one line."""
    pad = ' ' * depth
    if isinstance(x, dict):
        if not x:
            return '{}'
        body = ',\n'.join(f'{pad} {json.dumps(str(k))}: {to_json(v, depth + 1)}' for k, v in x.items())
        return '{\n' + body + '\n' + pad + '}'
    if isinstance(x, (list, tuple)):
        if all(not isinstance(v, (dict, list, tuple)) for v in x):
            return json.dumps(list(x), allow_nan=False)
        return '[\n' + ',\n'.join(f'{pad} {to_json(v, depth + 1)}' for v in x) + '\n' + pad + ']'
    return json.dumps(x, allow_nan=False)


def dumps(product) -> str:
    """The product with its numbers to four figures, the producer's hashes and constants as they are."""
    return to_json({k: (v if k == 'producer' else rnd(v)) for k, v in product.items()}) + '\n'


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out', type=Path, default=OUT)
    args = parser.parse_args(argv)
    product = results()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(dumps(product))
    for c in CORNERS:
        h = product['questions']['holds'][c]
        print(f"{c}: grown in year {h['start_year']:g}; left alone, CO2 falls to 29 Pa in "
              f"{h['years_to']['diverse_ecosystems_29_pa']:.0f} years")
    print(f'wrote {args.out}')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
