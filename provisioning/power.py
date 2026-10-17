"""The power system for the Open Moon's beamed power: a screen of the design chosen on 12 October 2026.

    python -m provisioning.power
    # -> provisioning/results/power.json

On 12 October 2026 the author chose beamed power from the ring fleet for food and the sky towns, on the design set out
in the domain synthesis (research/studies/domain_synthesis/README.md, "Power from the array: the chosen
distribution"): 2.45 GHz beams to the sky towns' crowns, closed receivers at the ground's loads joined by a grid,
several beams for every receiver, retrodirective beams, and exposure above 10 W/m2 kept to closed ground, decks and
columns. This screen answers provisioning question 3 for that design:

- absorption: ITU-R P.676-13 Annex 1 (equations 1-9: 44 oxygen lines, the Debye spectrum and pressure-induced
  nitrogen, 35 water-vapour lines; equations 13-19b for slant paths through exponentially thickening layers, with
  ITU-R P.453's refractivity), checked against ITU-R Study Group 3's own validation examples on Earth and then
  integrated through the atmosphere domain's moon_1.2atm column at 295 K, with oxygen at Earth's partial pressure;
- rain: ITU-R P.838-3 at both bands, checked against its own table, with the Moon's slower-falling rain and the
  storms of the cloud-resolving rings;
- the ring: the held matched stack of 1,859 rings, its tiles in view by latitude through the lunar cycle and the
  libration of the ring plane, the passes and handovers, the power a tile's collector can beam, the aperture each
  beam needs by the 1978 reference's tau, and the tiles the Moon's intake takes;
- the grid: receivers at the loads as the distributions study places C' with the food blend, the high latitudes,
  lines and their losses (US EIA 2018), and what backs up a receiver;
- reserves: storms, sky-reef crossings, single tiles lost, and a fault across the whole ring, with the critical loads
  and their storage;
- the summit metropolis's power plan.

Closed-form arithmetic and one-dimensional radiative transfer; no beam, receiver or grid is designed or simulated.
"""
from __future__ import annotations
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

from atmosphere.radiative_convective import thermodynamics as th
from shared.constants import (EARTH_RADIUS, JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_EQUATOR_TO_ECLIPTIC_DEG, MOON_GM,
                              MOON_ORBIT_INCLINATION_DEG, MOON_RADIUS, MOON_SURFACE_GRAVITY, SOLAR_CONSTANT,
                              SPEED_OF_LIGHT, SYNODIC_MONTH_DAYS)
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SCHEMA = 'terluna.provisioning.power/1'
OUT = HERE / 'results' / 'power.json'
SOURCES = HERE / 'power_sources.json'
CODE = ('provisioning/power.py', 'atmosphere/radiative_convective/thermodynamics.py')

STUDIES = 'research/studies'
LEDGERS = f'{STUDIES}/domain_synthesis/results/ledgers.json'
POPULATION = 'provisioning/results/population.json'
LINKS = f'{STUDIES}/solar_shield_array/results/electromagnetic.json'
LAYOUT = f'{STUDIES}/solar_shield_array/results/ring_layout.json'
LEDGER = f'{STUDIES}/solar_shield_array/results/integrated_ledger.json'
ARRAY_HEAT = f'{STUDIES}/solar_shield_array/results/array_heat.json'
DISTRIBUTIONS = f'{STUDIES}/distributions/results/distributions.json'
BLEND = f'{STUDIES}/food_systems/results/food_blend.json'
WAYS = f'{STUDIES}/ways_of_living/results/ways_of_living.json'
SHIPS = f'{STUDIES}/sky_ships/results/sky_ships.json'
VOLUME = f'{STUDIES}/inhabited_volume/results/inhabited_volume.json'
SKY = f'{STUDIES}/sky_ecology/results/sky_ecology.json'
TOWER = f'{STUDIES}/summit_tower/results/summit_tower.json'
PEOPLE = 'biosphere/ecology/results/people.json'
SCREEN = 'provisioning/results/first_screen.json'
COLUMN_TABLE = 'atmosphere/radiative_convective/results/inverse_climate.csv'
RINGS = {'equator': 'climate/results/crm/ring_ring_equator.json',
         '45n': 'climate/results/crm/ring_ring_45n_lsw.json', '45s': 'climate/results/crm/ring_ring_45s_lsw.json',
         '80n': 'climate/results/crm/ring_ring_80n_lsw.json', '80s': 'climate/results/crm/ring_ring_80s_lsw.json'}
INPUTS = (LEDGERS, POPULATION, LINKS, LAYOUT, LEDGER, ARRAY_HEAT, DISTRIBUTIONS, BLEND, WAYS, SHIPS, VOLUME, SKY, TOWER,
          PEOPLE, SCREEN, COLUMN_TABLE, *RINGS.values(), 'provisioning/power_sources.json')

BANDS_GHZ = (2.45, 5.8)
KM_PER_DEG = math.pi * MOON_RADIUS / 180.0 / 1e3
HOURS_PER_YEAR = JULIAN_YEAR_DAYS * 24.0

# ---------------------------------------------------------------------------------------------------------------
# ITU-R P.676-13 (08/2022), Annex 1, Table 1: oxygen lines f0 (GHz), a1..a6; Table 2: water-vapour lines f0, b1..b6.
OXYGEN = (
    (50.474214, 0.975, 9.651, 6.690, 0.0, 2.566, 6.850), (50.987745, 2.529, 8.653, 7.170, 0.0, 2.246, 6.800),
    (51.503360, 6.193, 7.709, 7.640, 0.0, 1.947, 6.729), (52.021429, 14.320, 6.819, 8.110, 0.0, 1.667, 6.640),
    (52.542418, 31.240, 5.983, 8.580, 0.0, 1.388, 6.526), (53.066934, 64.290, 5.201, 9.060, 0.0, 1.349, 6.206),
    (53.595775, 124.600, 4.474, 9.550, 0.0, 2.227, 5.085), (54.130025, 227.300, 3.800, 9.960, 0.0, 3.170, 3.750),
    (54.671180, 389.700, 3.182, 10.370, 0.0, 3.558, 2.654), (55.221384, 627.100, 2.618, 10.890, 0.0, 2.560, 2.952),
    (55.783815, 945.300, 2.109, 11.340, 0.0, -1.172, 6.135), (56.264774, 543.400, 0.014, 17.030, 0.0, 3.525, -0.978),
    (56.363399, 1331.800, 1.654, 11.890, 0.0, -2.378, 6.547), (56.968211, 1746.600, 1.255, 12.230, 0.0, -3.545, 6.451),
    (57.612486, 2120.100, 0.910, 12.620, 0.0, -5.416, 6.056), (58.323877, 2363.700, 0.621, 12.950, 0.0, -1.932, 0.436),
    (58.446588, 1442.100, 0.083, 14.910, 0.0, 6.768, -1.273), (59.164204, 2379.900, 0.387, 13.530, 0.0, -6.561, 2.309),
    (59.590983, 2090.700, 0.207, 14.080, 0.0, 6.957, -0.776), (60.306056, 2103.400, 0.207, 14.150, 0.0, -6.395, 0.699),
    (60.434778, 2438.000, 0.386, 13.390, 0.0, 6.342, -2.825), (61.150562, 2479.500, 0.621, 12.920, 0.0, 1.014, -0.584),
    (61.800158, 2275.900, 0.910, 12.630, 0.0, 5.014, -6.619), (62.411220, 1915.400, 1.255, 12.170, 0.0, 3.029, -6.759),
    (62.486253, 1503.000, 0.083, 15.130, 0.0, -4.499, 0.844), (62.997984, 1490.200, 1.654, 11.740, 0.0, 1.856, -6.675),
    (63.568526, 1078.000, 2.108, 11.340, 0.0, 0.658, -6.139), (64.127775, 728.700, 2.617, 10.880, 0.0, -3.036, -2.895),
    (64.678910, 461.300, 3.181, 10.380, 0.0, -3.968, -2.590), (65.224078, 274.000, 3.800, 9.960, 0.0, -3.528, -3.680),
    (65.764779, 153.000, 4.473, 9.550, 0.0, -2.548, -5.002), (66.302096, 80.400, 5.200, 9.060, 0.0, -1.660, -6.091),
    (66.836834, 39.800, 5.982, 8.580, 0.0, -1.680, -6.393), (67.369601, 18.560, 6.818, 8.110, 0.0, -1.956, -6.475),
    (67.900868, 8.172, 7.708, 7.640, 0.0, -2.216, -6.545), (68.431006, 3.397, 8.652, 7.170, 0.0, -2.492, -6.600),
    (68.960312, 1.334, 9.650, 6.690, 0.0, -2.773, -6.650), (118.750334, 940.300, 0.010, 16.640, 0.0, -0.439, 0.079),
    (368.498246, 67.400, 0.048, 16.400, 0.0, 0.000, 0.000), (424.763020, 637.700, 0.044, 16.400, 0.0, 0.000, 0.000),
    (487.249273, 237.400, 0.049, 16.000, 0.0, 0.000, 0.000), (715.392902, 98.100, 0.145, 16.000, 0.0, 0.000, 0.000),
    (773.839490, 572.300, 0.141, 16.200, 0.0, 0.000, 0.000), (834.145546, 183.100, 0.145, 14.700, 0.0, 0.000, 0.000))
WATER = (
    (22.235080, .1079, 2.144, 26.38, .76, 5.087, 1.00), (67.803960, .0011, 8.732, 28.58, .69, 4.930, .82),
    (119.995940, .0007, 8.353, 29.48, .70, 4.780, .79), (183.310087, 2.273, .668, 29.06, .77, 5.022, .85),
    (321.225630, .0470, 6.179, 24.04, .67, 4.398, .54), (325.152888, 1.514, 1.541, 28.23, .64, 4.893, .74),
    (336.227764, .0010, 9.825, 26.93, .69, 4.740, .61), (380.197353, 11.67, 1.048, 28.11, .54, 5.063, .89),
    (390.134508, .0045, 7.347, 21.52, .63, 4.810, .55), (437.346667, .0632, 5.048, 18.45, .60, 4.230, .48),
    (439.150807, .9098, 3.595, 20.07, .63, 4.483, .52), (443.018343, .1920, 5.048, 15.55, .60, 5.083, .50),
    (448.001085, 10.41, 1.405, 25.64, .66, 5.028, .67), (470.888999, .3254, 3.597, 21.34, .66, 4.506, .65),
    (474.689092, 1.260, 2.379, 23.20, .65, 4.804, .64), (488.490108, .2529, 2.852, 25.86, .69, 5.201, .72),
    (503.568532, .0372, 6.731, 16.12, .61, 3.980, .43), (504.482692, .0124, 6.731, 16.12, .61, 4.010, .45),
    (547.676440, .9785, .158, 26.00, .70, 4.500, 1.00), (552.020960, .1840, .158, 26.00, .70, 4.500, 1.00),
    (556.935985, 497.0, .159, 30.86, .69, 4.552, 1.00), (620.700807, 5.015, 2.391, 24.38, .71, 4.856, .68),
    (645.766085, .0067, 8.633, 18.00, .60, 4.000, .50), (658.005280, .2732, 7.816, 32.10, .69, 4.140, 1.00),
    (752.033113, 243.4, .396, 30.86, .68, 4.352, .84), (841.051732, .0134, 8.177, 15.90, .33, 5.760, .45),
    (859.965698, .1325, 8.055, 30.60, .68, 4.090, .84), (899.303175, .0547, 7.914, 29.85, .68, 4.530, .90),
    (902.611085, .0386, 8.429, 28.65, .70, 5.100, .95), (906.205957, .1836, 5.110, 24.08, .70, 4.700, .53),
    (916.171582, 8.400, 1.441, 26.73, .70, 5.150, .78), (923.112692, .0079, 10.293, 29.00, .70, 5.000, .80),
    (970.315022, 9.009, 1.919, 25.50, .64, 4.940, .67), (987.926764, 134.6, .257, 29.85, .68, 4.550, .90),
    (1780.000000, 17506., .952, 196.3, 2.00, 24.15, 5.00))
# ITU-R Study Group 3, validation examples rev. 8.3.0 (CG-3M3J-13-ValEx-Rev8.3.0.xlsx), sheet "P.676-13 SpAtt":
# gamma_o and gamma_w (dB/km) at a dry pressure of 1,013.25 hPa, 288.15 K and 7.5 g/m3 of water vapour.
ITU_STANDARD = dict(p_hpa=1013.25, t_k=288.15, rho_g_m3=7.5)
ITU_SPECIFIC = {1.0: (5.3886581679065499E-3, 5.0904617324964401E-5), 2.0: (6.7160384744085E-3, 2.0438137865037601E-4),
                3.0: (7.0759601327262398E-3, 4.62792856493146E-4), 5.0: (7.4004262562906101E-3, 1.31345892117269E-3),
                6.0: (7.5372120829487396E-3, 1.92157180441783E-3), 10.0: (8.2244167027098799E-3, 5.9741252454767199E-3),
                12.0: (8.6982640687735709E-3, 9.5353882202459297E-3),
                20.0: (1.18835504778076E-2, 9.7047304815111701E-2),
                22.0: (1.3130222965391699E-2, 0.17420703333692), 60.0: (14.6234747964861, 0.15484184063624701),
                90.0: (3.8869711072423499E-2, 0.341973394422181), 118.0: (1.13486620187051, 0.605921952594212),
                130.0: (4.1509083599522802E-2, 0.75184470364612899), 183.0: (1.2733908835870901E-2, 27.665008314166499),
                350.0: (3.6805627050457097E-2, 10.1226731395126)}
# The same file, sheet "P.676-13 A_Gas_A1_2.2.1a", example 1: 28 GHz from the ground at an apparent elevation of 30
# degrees through ITU-R P.835-6's mean annual global reference atmosphere, 922 layers.
ITU_SLANT = dict(f_ghz=28.0, elevation_deg=30.0, layers=922, attenuation_db=0.47081173472870474)
# P.676's line strengths take Earth's dry air: its oxygen and nitrogen mole fractions (the column model's Earth air).
EARTH_O2, EARTH_N2 = 0.2095, 0.7808

# ITU-R P.838-3 (03/2005), Tables 1-4: (a_j, b_j, c_j), then m and c.
P838 = dict(
    kh=(((-5.33980, -0.10008, 1.13098), (-0.35351, 1.26970, 0.45400), (-0.23789, 0.86036, 0.15354),
         (-0.94158, 0.64552, 0.16817)), -0.18961, 0.71147),
    kv=(((-3.80595, 0.56934, 0.81061), (-3.44965, -0.22911, 0.51059), (-0.39902, 0.73042, 0.11899),
         (0.50167, 1.07319, 0.27195)), -0.16398, 0.63297),
    ah=(((-0.14318, 1.82442, -0.55187), (0.29591, 0.77564, 0.19822), (0.32177, 0.63773, 0.13164),
         (-5.37610, -0.96230, 1.47828), (16.1721, -3.29980, 3.43990)), 0.67849, -1.95537),
    av=(((-0.07771, 2.33840, -0.76284), (0.56727, 0.95545, 0.54039), (-0.20238, 1.14520, 0.26809),
         (-48.2991, 0.791669, 0.116226), (48.5833, 0.791459, 0.116479)), -0.053739, 0.83433))
# Its Table 5 (kH, alphaH, kV, alphaV), printed to four or five figures, and the validation file's sheet "P.838-3
# Sp.Att" at full precision.
P838_TABLE5 = {2.0: (0.0000847, 1.0664, 0.0000998, 0.9490), 2.5: (0.0001321, 1.1209, 0.0001464, 1.0085),
               3.0: (0.0001390, 1.2322, 0.0001942, 1.0688), 5.0: (0.0002162, 1.6969, 0.0002428, 1.5317),
               5.5: (0.0003909, 1.6499, 0.0003115, 1.5882), 6.0: (0.0007056, 1.5900, 0.0004878, 1.5728)}
P838_VALIDATION = {14.25: (3.9186742018373702E-2, 1.1352800825355001, 4.3451215355383301E-2, 1.05905291969976),
                   29.0: (0.22240103375669601, 0.95802573203141295, 0.212395484704367, 0.92032488857384398)}

# ---------------------------------------------------------------------------------------------------------------
# The design column: the atmosphere domain's moon_1.2atm scenario (atmosphere/radiative_convective/run.py) at the design
# case's surface temperature, with its Manabe-Wetherald humidity and 200 K stratosphere, on the product's levels.
COLUMN = dict(dry_pressure_pa=121590.0, co2_ppm=400.0, stratosphere_k=200.0, levels=(70, 20))
COLUMN_TS_K = 295.0
COLUMN_SENSITIVITY_K = (290.0, 300.0)
ELEVATIONS_DEG = tuple(float(e) for e in range(1, 91))

# Earth reference facts the screen uses.
G_MORRISON = 9.81                 # the gravity Morrison's fall-speed laws were fitted for (climate/crm)
RAIN_B = 0.8                      # Morrison's raindrop fall-speed exponent, V = a D^b, so V ~ g^((b+1)/3) (climate/crm)
HVDC_PER_1000_KM = (0.02, 0.035)  # line losses per 1,000 km (US EIA 2018)
HVDC_CONVERTER = 0.01             # each converter station (US EIA 2018)
PUBLIC_LIMIT_W_M2 = 10.0          # ICNIRP 2020, general public, whole body, 2-300 GHz
EMERGENCY_DAYS = 7.0              # on-site fuel for safety-related emergency power (US NRC Regulatory Guide 1.137)
H2_LHV_KWH_KG = 33.3              # hydrogen's lower heating value, 120 MJ/kg
H2_ROUND_TRIP = 0.40              # hydrogen storage, electricity back over electricity in (Headley & Schoenung 2020)

# Design guesses, each labelled in the product.
USABLE_ELEVATION_DEG = 30.0       # the lowest elevation a beam is planned from
BEAMS_PER_RECEIVER = 3            # each rated for half the load, so any two carry it
BEAM_RATING_SHARE = 0.5
RECEIVER_SIDE_KM = 3.0            # a ground receiver's flat-filled patch, a few kilometres across
TILE_UTILISATION = 0.5            # share of a transmitting tile's rating in use, its loads lying in part of its view
SMALL_COLLECTOR_SHARE = 0.1       # the variant: collectors over a tenth of each transmitting tile, ten times the tiles
TILE_FAULTS_PER_YEAR = 0.1        # whole-tile beam outages a transmitting tile suffers each year, any cause
REASSIGN_S = 600.0                # time to put a fresh beam on a receiver that lost one
FUEL_CELL = 0.5                   # electricity per unit of hydrogen's lower heating value
BATTERY_KWH_KG = (0.15, 0.3)      # pack specific energy
BRIDGE_HOURS = 1.0                # batteries carry the critical loads while the hydrogen plant starts
CRITICAL_SHARE = dict(residents=0.10, sky_control=0.20, water=1.0, cultures=0.10, lamps=0.05, close_computing=0.25,
                      hydrogen_remade=0.0, trips=0.0)
DRY_HIGHLAND_KM = 4.0             # height of the dry highland cities' receivers (C' places them on ground above 4 km)
REEF_LAYER_KM = (10.0, 20.0)      # sky reefs drift at 10-20 km (sky ecology)
REEF_ABOVE_TOWN_KM = (5.0, 10.0)  # the synthesis's reefs 5-10 km above a town, for the offset a switch needs

EPS_DEG = MOON_EQUATOR_TO_ECLIPTIC_DEG + MOON_ORBIT_INCLINATION_DEG   # the ring plane's pole from the lunar pole
OMEGA_SUN_FRAME = 2 * math.pi / (SYNODIC_MONTH_DAYS * JULIAN_DAY)       # the Moon's turn under the Sun-fixed rings

EVIDENCE = ("Screen-level arithmetic and one-dimensional radiative transfer on committed products. The gaseous "
            "absorption is ITU-R P.676-13 Annex 1 line by line, reproduced to 1e-12 against ITU-R Study Group 3's "
            "validation examples on Earth, then integrated through the atmosphere domain's prescribed moon_1.2atm "
            "column at 295 K with oxygen at Earth's partial pressure: a computed result for a prescribed clear-sky "
            "column, on line parameters fitted to Earth's air. Rain uses ITU-R P.838-3 with "
            "Earth's drop spectra, its rate scaled by the Moon's slower fall (an estimate). The ring's view is "
            "geometry on the held matched stack, circular rings at their layout radii and tilts. Receivers, beams, "
            "grid, faults and storage are sized from the products and labelled design guesses; nothing is "
            "simulated through time beyond sampling the lunar cycle and the libration of the ring plane.")
READING_RULE = ("Pairs are [low-cost end, high-cost end] of C' with the food blend as the domain synthesis gives them, "
                "unless labelled. Intake is the beam power arriving at receivers and crowns, which the synthesis "
                "counts as heat on the Moon; delivered power is what reaches the loads. Absorption in dB and as the "
                "share of the beam lost; elevations in degrees; heights in km above the column's surface (sea level) "
                "unless named; latitudes are lunar, with the north and south alike. Tile counts are of the fleet's "
                "28 million tiles; times in hours unless named; storage in TWh of electricity delivered.")


def sha256(path) -> str:
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def load(path):
    return json.loads((ROOT / path).read_text())


def sig(x, n=4):
    """x rounded to n significant figures through dicts and lists."""
    if isinstance(x, dict):
        return {k: sig(v, n) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [sig(v, n) for v in x]
    if isinstance(x, (bool, str)) or x is None:
        return x
    if isinstance(x, (int, np.integer)):
        return int(x)
    x = float(x)
    if x == 0.0 or not math.isfinite(x):
        return x
    return float(f'{x:.{n}g}')


def pair(f):
    return [f(0), f(1)]


def loss(db):
    """Share of the power lost to an attenuation in dB."""
    return 1.0 - 10.0 ** (-np.asarray(db, float) / 10.0)


# ---------------------------------------------------------------------------------------------------------------
# Gaseous absorption.

def _shape(f, fi, df, delta):
    """P.676-13 equation 5, the line shape factor."""
    return f / fi * ((df - delta * (fi - f)) / ((fi - f) ** 2 + df ** 2)
                     + (df - delta * (fi + f)) / ((fi + f) ** 2 + df ** 2))


def specific_attenuation(f, p_dry_hpa, e_hpa, t_k, o2_scale=1.0, n2_scale=1.0):
    """P.676-13 equations 1-9: gamma_o and gamma_w (dB/km) at f GHz, dry pressure and vapour pressure (hPa) and T (K).

    o2_scale multiplies the oxygen line strengths and the Debye term, which go with oxygen's number density; the line
    widths and the Debye width go with the total pressure as P.676 gives them. n2_scale multiplies the pressure-induced
    nitrogen term, which goes as the square of nitrogen's pressure. Both are 1 for Earth's air."""
    p, e, t = (np.asarray(v, float) for v in (p_dry_hpa, e_hpa, t_k))
    th_ = 300.0 / t
    ox = np.asarray(OXYGEN)
    fi, a1, a2, a3, a4, a5, a6 = (c[None, :] for c in ox.T)
    P, E, TH = p[..., None], e[..., None], th_[..., None]
    s = a1 * 1e-7 * P * TH ** 3 * np.exp(a2 * (1 - TH)) * o2_scale
    df = a3 * 1e-4 * (P * TH ** (0.8 - a4) + 1.1 * E * TH)
    df = np.sqrt(df ** 2 + 2.25e-6)
    delta = (a5 + a6 * TH) * 1e-4 * (P + E) * TH ** 0.8
    n_ox = (s * _shape(f, fi, df, delta)).sum(-1)
    d = 5.6e-4 * (p + e) * th_ ** 0.8
    n_dry = f * p * th_ ** 2 * (6.14e-5 / (d * (1 + (f / d) ** 2)) * o2_scale
                                + 1.4e-12 * p * th_ ** 1.5 / (1 + 1.9e-5 * f ** 1.5) * n2_scale)
    wv = np.asarray(WATER)
    fi, b1, b2, b3, b4, b5, b6 = (c[None, :] for c in wv.T)
    s = b1 * 1e-1 * E * TH ** 3.5 * np.exp(b2 * (1 - TH))
    df = b3 * 1e-4 * (P * TH ** b4 + b5 * E * TH ** b6)
    df = 0.535 * df + np.sqrt(0.217 * df ** 2 + 2.1316e-12 * fi ** 2 / TH)
    n_w = (s * _shape(f, fi, df, 0.0)).sum(-1)
    return (0.1820 * f * (n_ox + n_dry)).reshape(p.shape), (0.1820 * f * n_w).reshape(p.shape)


def refractivity(p_dry_hpa, e_hpa, t_k):
    """ITU-R P.453 equations 1-2, N-units."""
    return 77.6 * p_dry_hpa / t_k + 72.0 * e_hpa / t_k + 3.75e5 * e_hpa / t_k ** 2


def layer_edges(h_low, h_top, imax=None):
    """P.676-13 equations 14-15: layers 10 cm thick at the ground, thickening by e^(1/100) each, cut at h_low and the
    column's top (or the first imax layers)."""
    i = np.arange(1, 4001)
    edges = 1e-4 * (np.exp((i - 1) / 100.0) - 1) / (math.exp(0.01) - 1)
    if imax is not None:
        return edges[:imax + 1]
    edges = edges[(edges > h_low) & (edges < h_top)]
    return np.concatenate([[h_low], edges, [h_top]])


def slant_path(f, elevation_deg, profile, radius_km, edges, o2_scale=1.0, n2_scale=1.0, refraction=True):
    """P.676-13 equations 13, 17 and 19b: attenuation (dB) from the lowest edge to the top at an apparent elevation,
    split into oxygen and water vapour; profile(h_km) gives T (K), dry pressure and vapour pressure (hPa)."""
    hb, delta = edges[:-1], np.diff(edges)
    t, p, e = profile(hb + delta / 2)
    n = 1 + refractivity(p, e, t) * 1e-6 if refraction else np.ones_like(t)
    r = radius_km + hb
    beta = np.arcsin(np.clip(n[0] * r[0] * math.sin(math.radians(90.0 - elevation_deg)) / (n * r), -1.0, 1.0))
    a = -r * np.cos(beta) + np.sqrt(r ** 2 * np.cos(beta) ** 2 + 2 * r * delta + delta ** 2)
    go, gw = specific_attenuation(f, p, e, t, o2_scale, n2_scale)
    return float(a @ go), float(a @ gw)


def p835_global(h_km):
    """ITU-R P.835-6 mean annual global reference atmosphere: T (K), total pressure (hPa), vapour pressure (hPa) at
    geometric heights (km); water vapour 7.5 exp(-h/2) g/m3 held to a mixing ratio of at least 2e-6."""
    h = np.asarray(h_km, float)
    hp = 6356.766 * h / (6356.766 + h)
    t, p = np.empty_like(h), np.empty_like(h)
    segments = (
        (0.0, 11.0, lambda x: 288.15 - 6.5 * x, lambda x: 1013.25 * (288.15 / (288.15 - 6.5 * x)) ** (-34.1632 / 6.5)),
        (11.0, 20.0, lambda x: 216.65 + 0 * x, lambda x: 226.3226 * np.exp(-34.1632 * (x - 11) / 216.65)),
        (20.0, 32.0, lambda x: 216.65 + (x - 20), lambda x: 54.74980 * (216.65 / (216.65 + (x - 20))) ** 34.1632),
        (32.0, 47.0, lambda x: 228.65 + 2.8 * (x - 32),
         lambda x: 8.680422 * (228.65 / (228.65 + 2.8 * (x - 32))) ** (34.1632 / 2.8)),
        (47.0, 51.0, lambda x: 270.65 + 0 * x, lambda x: 1.109106 * np.exp(-34.1632 * (x - 47) / 270.65)),
        (51.0, 71.0, lambda x: 270.65 - 2.8 * (x - 51),
         lambda x: 0.6694167 * (270.65 / (270.65 - 2.8 * (x - 51))) ** (-34.1632 / 2.8)),
        (71.0, 84.852, lambda x: 214.65 - 2.0 * (x - 71),
         lambda x: 0.03956649 * (214.65 / (214.65 - 2.0 * (x - 71))) ** (-34.1632 / 2.0)))
    done = np.zeros(h.shape, bool)
    for lo, hi, ft, fp in segments:
        m = (hp >= lo) & (hp < hi) & (h < 86.0)
        t[m], p[m] = ft(hp[m]), fp(hp[m])
        done |= m
    m, hh = ~done, h[~done]
    t[m] = np.where(hh <= 91.0, 186.8673,
                    263.1905 - 76.3232 * np.sqrt(np.clip(1 - ((hh - 91.0) / 19.9429) ** 2, 0, None)))
    p[m] = np.exp(95.571899 - 4.011801 * hh + 6.424731e-2 * hh ** 2 - 4.789660e-4 * hh ** 3 + 1.340543e-6 * hh ** 4)
    e = 7.5 * np.exp(-h / 2.0) * t / 216.7
    e = np.where(e / p < 2e-6, 2e-6 * p, e)
    return t, p, e


def earth_profile(h_km):
    t, p, e = p835_global(h_km)
    return t, p - e, e


def earth_checks():
    """The model against ITU-R Study Group 3's own Earth values, and Earth's zenith at the two bands."""
    e_std = ITU_STANDARD['rho_g_m3'] * ITU_STANDARD['t_k'] / 216.7          # P.676-13 equation 4
    rows, worst = {}, 0.0
    for f, (go, gw) in ITU_SPECIFIC.items():
        mo, mw = specific_attenuation(f, ITU_STANDARD['p_hpa'], e_std, ITU_STANDARD['t_k'])
        rel = max(abs(float(mo) / go - 1), abs(float(mw) / gw - 1))
        worst = max(worst, rel)
        rows[f'{f:g}'] = dict(itu_gamma_o=go, model_gamma_o=float(mo), itu_gamma_w=gw, model_gamma_w=float(mw))
    edges = layer_edges(0.0, 100.0, imax=ITU_SLANT['layers'])
    ao, aw = slant_path(ITU_SLANT['f_ghz'], ITU_SLANT['elevation_deg'], earth_profile, EARTH_RADIUS / 1e3, edges)
    zenith = {}
    for f in BANDS_GHZ:
        zo, zw = slant_path(f, 90.0, earth_profile, EARTH_RADIUS / 1e3, edges)
        zenith[f'{f:g}_ghz'] = dict(dry_db=zo, wet_db=zw, total_db=zo + zw, share_lost=float(loss(zo + zw)))
    return dict(specific_attenuation=rows, worst_relative_error=worst,
                slant=dict(**ITU_SLANT, model_db=ao + aw, relative_error=(ao + aw) / ITU_SLANT['attenuation_db'] - 1),
                earth_zenith_p835_global=zenith)


def moon_column(ts_k):
    air = th.earthlike_air(COLUMN['dry_pressure_pa'], COLUMN['co2_ppm'])
    spec = th.ColumnSpec(th.MOON, air, ts_k, COLUMN['stratosphere_k'], th.manabe_wetherald(),
                         troposphere_levels=COLUMN['levels'][0], upper_levels=COLUMN['levels'][1])
    col = th.build_column(spec)
    z = col['z_m'] / 1e3
    lnp, lnx = np.log(col['p_pa']), np.log(np.maximum(col['x_h2o'], 1e-30))

    def profile(h_km):
        p = np.exp(np.interp(h_km, z, lnp)) / 100.0
        e = np.exp(np.interp(h_km, z, lnx)) * p
        return np.interp(h_km, z, col['t_k']), p - e, e
    water = float(np.sum(col['layer_mass_kg_m2'] * col['layer_x_h2o'] * th.MOLAR_MASS['H2O']
                         / (th.MOLAR_MASS['H2O'] * col['layer_x_h2o'] + air.molar_mass * (1 - col['layer_x_h2o']))))
    trop = float(np.interp(math.log(col['tropopause_pa']), np.log(col['p_pa'][::-1]), col['z_m'][::-1]) / 1e3)
    freezing = float(np.interp(-273.15, -col['t_k'], z))                 # the 0 C level, T falling with height
    return dict(profile=profile, air=air, top_km=float(z[-1]), surface_pressure_pa=float(col['surface_pressure_pa']),
                water_column_kg_m2=water, tropopause_km=trop, freezing_level_km=freezing,
                air_column_t_m2=float(col['layer_mass_kg_m2'].sum() / 1e3))


def column_check(table_rows):
    """The column built here against the atmosphere product's rows at 290 and 300 K."""
    out = {}
    for ts in COLUMN_SENSITIVITY_K:
        row = next(r for r in table_rows if r['scenario'] == 'moon_1.2atm' and r['humidity'] == 'manabe_wetherald'
                   and float(r['ts_k']) == ts and float(r['co2_ppm']) == COLUMN['co2_ppm']
                   and float(r['stratosphere_k']) == COLUMN['stratosphere_k'] and r['shield'] == 'none')
        col = moon_column(ts)
        out[f'{ts:g}K'] = dict(product=dict(water_column_kg_m2=float(row['water_column_kg_m2']),
                                            tropopause_height_km=float(row['tropopause_height_km']),
                                            surface_pressure_pa=float(row['surface_pressure_pa'])),
                               here=dict(water_column_kg_m2=col['water_column_kg_m2'],
                                         tropopause_height_km=col['tropopause_km'],
                                         surface_pressure_pa=col['surface_pressure_pa']))
    return out


def receiver_heights(volume, tower):
    """Receiver heights (km above sea level): lowland homes at sea level, the dry highland cities, the sky towns at the
    inhabited-volume study's reference height and the summit metropolis on its ground."""
    return dict(lowland=0.0, highland=DRY_HIGHLAND_KM, towns=float(volume['assumptions']['reference_height_km']),
                summit=tower['site']['ground_above_sea_m'] / 1e3)


def absorption_tables(col, heights):
    """Attenuation (dB) against apparent elevation for each band and receiver height, through the Moon's column."""
    o2_scale = col['air'].fractions['O2'] / EARTH_O2
    n2_scale = (col['air'].fractions['N2'] / EARTH_N2) ** 2
    tables = {}
    for f in BANDS_GHZ:
        for name, h0 in heights.items():
            edges = layer_edges(h0, col['top_km'])
            rows = [slant_path(f, e, col['profile'], MOON_RADIUS / 1e3, edges, o2_scale, n2_scale)
                    for e in ELEVATIONS_DEG]
            tables[(f, name)] = np.array([r[0] + r[1] for r in rows]), np.array([r[0] for r in rows])
    return tables, o2_scale, n2_scale


def absorption(col, tables, o2_scale, n2_scale, checks, heights):
    out = {}
    for f in BANDS_GHZ:
        rows = {}
        for name in heights:
            total, dry = tables[(f, name)]
            rows[name] = {f'{e:g}_deg': dict(db=float(total[int(e) - 1]), dry_db=float(dry[int(e) - 1]),
                                             share_lost=float(loss(total[int(e) - 1])))
                          for e in (90.0, 60.0, 45.0, 30.0, 20.0, 10.0, 5.0)}
        out[f'{f:g}_ghz'] = rows
    sens = {}
    for ts in COLUMN_SENSITIVITY_K:
        c = moon_column(ts)
        for f in BANDS_GHZ:
            edges = layer_edges(0.0, c['top_km'])
            zo, zw = slant_path(f, 90.0, c['profile'], MOON_RADIUS / 1e3, edges, o2_scale, n2_scale)
            sens[f'{ts:g}K_{f:g}_ghz_zenith_db'] = zo + zw
    earth = checks['earth_zenith_p835_global']
    edges0 = layer_edges(0.0, col['top_km'])
    zo_e, zw_e = slant_path(BANDS_GHZ[0], 90.0, col['profile'], MOON_RADIUS / 1e3, edges0, 1.0, n2_scale)
    return dict(
        model='ITU-R P.676-13 Annex 1, equations 1-9 (line by line) and 13-19b (slant paths), refractivity by '
              'ITU-R P.453; layers 10 cm thick at the receiver, thickening by e^(1/100), to the column top at 1 Pa',
        composition=dict(o2_mole_fraction=col['air'].fractions['O2'], o2_strength_scale=o2_scale,
                         n2_pressure_induced_scale=n2_scale,
                         rule='oxygen line strengths and the Debye term scaled by the oxygen mole fraction over '
                              "Earth's 0.2095; widths kept as P.676 gives them for air at the total pressure"),
        column=dict(scenario='moon_1.2atm', surface_temperature_k=COLUMN_TS_K, **{k: col[k] for k in (
            'surface_pressure_pa', 'water_column_kg_m2', 'tropopause_km', 'freezing_level_km', 'top_km',
            'air_column_t_m2')}),
        receiver_heights_km=heights, by_elevation=out,
        zenith_sensitivity_to_surface_temperature_db=sens,
        zenith_2_45_with_earths_oxygen_fraction_db=zo_e + zw_e,
        moon_over_earth_zenith={f'{f:g}_ghz': float(tables[(f, 'lowland')][0][-1]) / earth[f'{f:g}_ghz']['total_db']
                                for f in BANDS_GHZ})


# ---------------------------------------------------------------------------------------------------------------
# The ring and what a receiver sees.

def tile_samples(layout, fleet_tiles, n_rings=93, du_deg=1.0):
    """The held matched stack as n_rings circular rings, heights spread evenly between the flown extremes, radius
    interpolated in height from the flown sample, tilt from height over radius; nodes at right angles to the Sun."""
    rings = sorted(layout['flown']['matched']['rings'], key=lambda r: r['height_km'])
    h = np.array([r['height_km'] for r in rings])
    rad = np.array([r['radius_km'] for r in rings])
    hs = np.linspace(h[0], h[-1], n_rings)
    rs = np.interp(hs, h, rad)
    tilt = np.arcsin(hs / rs)
    u = np.radians(np.arange(0.0, 360.0, du_deg))
    pos = np.stack([rs[:, None] * np.cos(u)[None, :] * np.cos(tilt)[:, None],
                    rs[:, None] * np.sin(u)[None, :] * np.ones((n_rings, 1)),
                    rs[:, None] * np.cos(u)[None, :] * np.sin(tilt)[:, None]], -1)
    motion = np.sqrt(MOON_GM / (rs * 1e3) ** 3)
    return dict(pos=pos, weight=fleet_tiles / (n_rings * u.size), radius_km=rs, tilt_deg=np.degrees(tilt),
                height_km=hs, mean_motion=motion, du_deg=du_deg, rings=layout['stack']['rings'])


def receivers(lat_deg, hours_deg, alpha_deg, height_km):
    """Receiver positions (km) in the Sun-fixed ring frame: x to the Sun, z the ring plane's pole. The lunar pole
    stands EPS_DEG from it toward azimuth alpha; the hour angle runs from local noon."""
    eps, a = math.radians(EPS_DEG), math.radians(alpha_deg)
    s = np.array([math.sin(eps) * math.cos(a), math.sin(eps) * math.sin(a), math.cos(eps)])
    e1 = np.array([1.0, 0.0, 0.0]) - s * s[0]
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(s, e1)
    phi, hh = math.radians(lat_deg), np.radians(np.asarray(hours_deg, float))
    r = MOON_RADIUS / 1e3 + height_km
    return r * (math.cos(phi) * np.cos(hh)[:, None] * e1 + math.cos(phi) * np.sin(hh)[:, None] * e2
                + math.sin(phi) * s)


def runs(mask):
    """Circular runs of True along the last axis: for each row, the sum of run lengths squared, their sum and the
    longest."""
    rows, n = mask.shape
    full = mask.all(1)
    start = np.argmin(mask, 1)
    idx = (start[:, None] + np.arange(n)[None, :]) % n
    rolled = np.take_along_axis(mask, idx, 1)
    flat = np.concatenate([rolled, np.zeros((rows, 1), bool)], 1).ravel()
    d = np.diff(np.concatenate([[0], flat.astype(np.int8), [0]]))
    s, e = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
    length = (e - s).astype(float)
    row = s // (n + 1)
    sq = np.bincount(row, length ** 2, rows)
    tot = np.bincount(row, length, rows)
    longest = np.zeros(rows)
    np.maximum.at(longest, row, length)
    sq[full], tot[full], longest[full] = n * n, n, n
    return sq, tot, longest


def view(tiles, lat_deg, tables, heights=('lowland', 'highland', 'towns'), hours=np.arange(0.0, 360.0, 15.0),
         alphas=np.arange(0.0, 360.0, 45.0)):
    """What a receiver at a latitude sees through the lunar cycle and the ring plane's libration: tiles above 0, 20,
    30 and 45 degrees, the highest tile, its passes above the usable elevation, and the absorption on its beams for
    each receiver height. Elevations are taken at the surface; a receiver 12 km up sees the tiles within 0.04 degrees
    of that."""
    pos, w = tiles['pos'], tiles['weight']
    n_h = len(hours)
    rel_rate = tiles['mean_motion'] + OMEGA_SUN_FRAME                  # retrograde rings under the turning Moon
    hours_per_sample = math.radians(tiles['du_deg']) / np.tile(rel_rate, n_h) / 3600.0
    acc = {k: [] for k in ('n0', 'n20', 'n30', 'n45', 'emax', 'longest_h', 'mean_pass_h')}
    absorb = {(f, hname, kind): [] for f in BANDS_GHZ for hname in heights for kind in ('best', 'usable')}
    for alpha in alphas:
        p = receivers(lat_deg, hours, alpha, 0.0)
        d = pos[None] - p[:, None, None, :]
        up = p / np.linalg.norm(p, axis=1)[:, None]
        el = np.degrees(np.arcsin(np.einsum('hrud,hd->hru', d, up) / np.linalg.norm(d, axis=-1)))
        flat = el.reshape(n_h, -1)
        usable = el >= USABLE_ELEVATION_DEG
        use_flat = usable.reshape(n_h, -1)
        emax = flat.max(1)
        for key, thr in (('n0', 0.0), ('n20', 20.0), ('n30', USABLE_ELEVATION_DEG), ('n45', 45.0)):
            acc[key].append((flat >= thr).sum(1) * w)
        acc['emax'].append(emax)
        sq, tot, longest = runs(usable.reshape(-1, usable.shape[-1]))
        acc['longest_h'].append((longest * hours_per_sample).reshape(n_h, -1).max(1))
        total = tot.reshape(n_h, -1).sum(1)
        acc['mean_pass_h'].append(np.where(total > 0, (sq * hours_per_sample).reshape(n_h, -1).sum(1)
                                           / np.maximum(total, 1), 0.0))
        idx = np.clip(np.rint(flat).astype(int), 1, len(ELEVATIONS_DEG)) - 1
        n_use = use_flat.sum(1)
        for f in BANDS_GHZ:
            for hname in heights:
                table = tables[(f, hname)][0]
                best = np.interp(np.maximum(emax, ELEVATIONS_DEG[0]), ELEVATIONS_DEG, table)
                mean_use = np.where(n_use > 0, (table[idx] * use_flat).sum(1) / np.maximum(n_use, 1), best)
                absorb[(f, hname, 'best')].append(best)
                absorb[(f, hname, 'usable')].append(mean_use)
    a = {k: np.concatenate(v) for k, v in acc.items()}
    seen = a['n30'] > 0
    out = dict(tiles_above_horizon=dict(mean=a['n0'].mean(), least=a['n0'].min()),
               tiles_above_20_deg=dict(mean=a['n20'].mean(), least=a['n20'].min()),
               tiles_above_usable=dict(mean=a['n30'].mean(), least=a['n30'].min()),
               tiles_above_45_deg=dict(mean=a['n45'].mean(), least=a['n45'].min()),
               highest_tile_deg=dict(mean=a['emax'].mean(), least=a['emax'].min()),
               share_of_time_with_no_usable_tile=float((~seen).mean()),
               share_of_time_with_no_tile_above_20_deg=float((a['n20'] == 0).mean()),
               longest_usable_pass_h=dict(mean=a['longest_h'][seen].mean() if seen.any() else 0.0,
                                          least=a['longest_h'][seen].min() if seen.any() else 0.0),
               typical_usable_pass_h=float(a['mean_pass_h'][seen].mean()) if seen.any() else 0.0)
    ab = {}
    for f in BANDS_GHZ:
        for hname in heights:
            best = np.concatenate(absorb[(f, hname, 'best')])
            use = np.concatenate(absorb[(f, hname, 'usable')])
            ab[f'{f:g}_ghz_{hname}'] = dict(best_tile_db=best.mean(), beams_db=use.mean(),
                                            beams_share_lost=float(loss(use).mean()))
    out['absorption'] = ab
    return out


# ---------------------------------------------------------------------------------------------------------------
# Loads by latitude.

def area_share(lo, hi):
    return math.sin(math.radians(hi)) - math.sin(math.radians(lo))


def ground_groups(dist):
    """The ground's people of C' with the blend (1 kg of heritage meat, the compromise rule), central composition, by
    latitude band and receiver height: the metropolis, the dry highland cities by band, and the scattered communities,
    47.5% within 30 degrees and the rest on the high-latitude lowland homes poleward of 48.6 degrees."""
    full = dist['moon']['cprime_with_blend']['1_kg']['wild_30_hybrid_20']['full']
    mid = lambda v: (v[0] + v[1]) / 2.0
    dry = {k: mid(v) for k, v in full['dry_people_by_band'].items()}
    metropolis = mid(full['people']['metropolis'])
    ground = 4.0
    scattered = ground - metropolis - sum(dry.values())
    low = scattered * mid(full['scattered_share_0_30'])
    high = scattered - low
    bands = dist['basis']['land_bands']['dry_land_km2']
    lowland_high = {(48.6, 60.0): bands['45-60'] * area_share(48.6, 60.0) / area_share(45.0, 60.0),
                    (60.0, 75.0): bands['60-75'], (75.0, 90.0): bands['75-90']}
    total_high = sum(lowland_high.values())
    groups = [dict(name='metropolis', lat=(5.0, 10.0), height='summit', billion=metropolis),
              dict(name='scattered_0_30', lat=(0.0, 30.0), height='lowland', billion=low),
              dict(name='dry_30_45', lat=(30.0, 45.0), height='highland', billion=dry['30-45']),
              dict(name='dry_45_60', lat=(45.0, 60.0), height='highland', billion=dry['45-60']),
              dict(name='dry_60_90', lat=(60.0, 90.0), height='highland', billion=dry['60-90'])]
    for (a, b), km2 in lowland_high.items():
        groups.append(dict(name=f'high_lowland_{a:g}_{b:g}', lat=(a, b), height='lowland',
                           billion=high * km2 / total_high))
    return groups, dict(metropolis=metropolis, dry_highland=dry, scattered=scattered, scattered_0_30=low,
                        scattered_high_latitude=high)


def town_shares(ways):
    """The sky towns' shares by absolute latitude at 10 km, from the design winds' gathering (both hemispheres)."""
    shares = {}
    for row in ways['winds']['gathering_by_lat']:
        a, b = row['lat_deg']
        key = (min(abs(a), abs(b)), max(abs(a), abs(b)))
        shares[key] = shares.get(key, 0.0) + row['share']
    total = sum(shares.values())
    return {k: v / total for k, v in sorted(shares.items())}


def bin_weights(lo, hi, edges):
    """Area weights of a latitude band over 5-degree bins."""
    w = np.zeros(len(edges) - 1)
    for i in range(len(edges) - 1):
        a, b = max(lo, edges[i]), min(hi, edges[i + 1])
        if b > a:
            w[i] = area_share(a, b)
    return w / w.sum()


# ---------------------------------------------------------------------------------------------------------------
# Rain.

def p838(f):
    """ITU-R P.838-3 equations 2-3: kH, alphaH, kV, alphaV, and the circular-polarization k and alpha (equations
    4-5 with a tilt of 45 degrees, independent of elevation)."""
    x = math.log10(f)

    def fit(table):
        terms, m, c = table
        return sum(a * math.exp(-((x - b) / cc) ** 2) for a, b, cc in terms) + m * x + c
    kh, kv = 10 ** fit(P838['kh']), 10 ** fit(P838['kv'])
    ah, av = fit(P838['ah']), fit(P838['av'])
    k = (kh + kv) / 2.0
    return dict(kh=kh, ah=ah, kv=kv, av=av, k=k, alpha=(kh * ah + kv * av) / (2.0 * k))


def rain_db_per_km(f, rate_mm_h, fall):
    c = p838(f)
    return c['k'] * (fall * rate_mm_h) ** c['alpha']


def ring_storms(ring):
    lat = abs(ring['latitude_deg']) if ring.get('tilt_deg') in (None, 0.0) else 0.0
    length_km = 2 * math.pi * MOON_RADIUS / 1e3 * math.cos(math.radians(lat))
    st = ring['storms']
    rr = ring['rain_rate_mm_h']
    wet = [rr['land']['wet_share'], rr['water']['wet_share']]
    out = dict(latitude_deg=lat, ring_km=length_km, wet_share_above_0_1_mm_h=wet,
               p99_mm_h=[rr['land']['p99'], rr['water']['p99']], max_mm_h=[rr['land']['max'], rr['water']['max']],
               storms_per_snapshot=st.get('per_snapshot', 0.0) or 0.0)
    if st.get('width_km'):
        out.update(width_km=st['width_km'], rain_max_mm_h=st['rain_max_mm_h'], cloud_top_km=st['cloud_top_km'],
                   coverage_above_1_mm_h=out['storms_per_snapshot'] * st['width_km']['median'] / length_km)
    else:
        out.update(coverage_above_1_mm_h=0.0)
    return out


def edge_buffer_km(spot_m, level_w_m2, limit_w_m2):
    """How far beyond a flat-filled patch's edge the beam stays above limit_w_m2, when the edge is the patch's level
    blurred by a Gaussian spot that holds 95% of its power within spot_m: the level times erfc(sqrt(2) x / w) / 2."""
    w = spot_m / 2.0 / math.sqrt(-math.log(0.05) / 2.0)
    lo, hi = 0.0, 10.0 * w
    for _ in range(80):
        mid = (lo + hi) / 2.0
        if level_w_m2 * 0.5 * math.erfc(math.sqrt(2.0) * mid / w) > limit_w_m2:
            lo = mid
        else:
            hi = mid
    return lo / 1e3


def path_in_rain(h0_km, top_km, width_km, elevation_deg):
    """Path length (km) inside a rain cell standing from the ground to top_km, width_km across, from a receiver at
    h0_km looking up at an elevation through its centre."""
    depth = max(top_km - h0_km, 0.0)
    e = math.radians(elevation_deg)
    up = depth / math.sin(e)
    across = width_km / max(math.cos(e), 1e-9)
    return min(up, across)


def crossing_factor(depth_km, width_km, elevation_deg):
    """How much more often a slant path through rain of the given depth meets cells of the given width than a
    vertical one: cells as discs, the path's ground track depth/tan(e) long."""
    if depth_km <= 0:
        return 1.0
    return 1.0 + 4.0 * depth_km / (math.pi * width_km * math.tan(math.radians(elevation_deg)))


# ---------------------------------------------------------------------------------------------------------------
# The screen.

def results():
    data = {p: load(p) for p in INPUTS if p.endswith('.json')}
    with open(ROOT / COLUMN_TABLE, newline='') as handle:
        table_rows = list(csv.DictReader(handle))
    ledgers, population, links = data[LEDGERS], data[POPULATION], data[LINKS]
    layout, ledger, aheat = data[LAYOUT], data[LEDGER], data[ARRAY_HEAT]
    dist, blend, ways, ships = data[DISTRIBUTIONS], data[BLEND], data[WAYS], data[SHIPS]
    volume, sky, tower, people, screen = data[VOLUME], data[SKY], data[TOWER], data[PEOPLE], data[SCREEN]
    rings = {k: data[v] for k, v in RINGS.items()}
    # The sky ships' storms are the equatorial ring's, and the synthesis's reef share is sky ecology's.
    ring_storms_eq, ship_storms = rings['equator']['storms'], ships['storms']
    if not (math.isclose(ship_storms['per_snapshot'], ring_storms_eq['per_snapshot'], rel_tol=1e-2)
            and math.isclose(ship_storms['width_km']['median'], ring_storms_eq['width_km']['median'], rel_tol=1e-2)):
        raise ValueError("The sky ships' storms are not the equatorial ring's")
    reef_share = sky['phosphorus']['design_guesses']['class_shares_of_cover']['sky_reef']
    if not math.isclose(reef_share, ledgers['heights']['sky_towns']['reef_share_of_cover'], rel_tol=1e-3):
        raise ValueError("The synthesis's reef share is not sky ecology's")

    # 1. The deep air's absorption.
    checks = earth_checks()
    col = moon_column(COLUMN_TS_K)
    heights = receiver_heights(volume, tower)
    tables, o2_scale, n2_scale = absorption_tables(col, heights)
    absorb = absorption(col, tables, o2_scale, n2_scale, checks, heights)
    air_check = [dict(height_km=r['height_km'], volume_pressure_atm=r['pressure_atm'],
                      column_pressure_atm=float(np.exp(np.interp(r['height_km'], *_column_lnp())) / 101325.0),
                      volume_temperature_c=r['temperature_c'],
                      column_temperature_c=float(col['profile'](np.array([r['height_km']]))[0][0] - 273.15))
                 for r in volume['air']]

    # 3. The ring: tiles in view by latitude.
    fleet_tiles = float(np.mean(aheat['fleet']['tiles']))
    tiles = tile_samples(layout, fleet_tiles)
    edges = np.arange(0.0, 91.0, 5.0)
    centres = (edges[:-1] + edges[1:]) / 2.0
    by_lat = {f'{c:g}': view(tiles, c, tables) for c in centres}
    summit_lat = tower['site']['lat_deg']
    summit_view = view(tiles, summit_lat, tables, heights=('summit',))

    # Loads by latitude.
    groups, composition = ground_groups(dist)
    towns_by_lat = town_shares(ways)
    ap = ledgers['array_power']
    sky_tw, ground_tw = ap['sky_towns']['intake_tw'], ap['ground']['intake_tw']
    towns = ledgers['heights']['sky_towns']['towns']
    town_w = np.zeros(len(centres))
    for (a, b), s in towns_by_lat.items():
        town_w += s * bin_weights(a, b, edges)
    ground_w = {}
    for g in groups:
        ground_w[g['name']] = (g, g['billion'] / 4.0 * bin_weights(*g['lat'], edges))

    def weighted(f, hname, key, weights):
        vals = np.array([by_lat[f'{c:g}']['absorption'][f'{f:g}_ghz_{hname}'][key] for c in centres])
        return float(vals @ weights)

    def ground_mean(f, key, beam=True):
        tot = 0.0
        for name, (g, w) in ground_w.items():
            if g['height'] == 'summit':
                tot += g['billion'] / 4.0 * summit_view['absorption'][f'{f:g}_ghz_summit'][key]
                continue
            vals = np.array([by_lat[f'{c:g}']['absorption'][f'{f:g}_ghz_{g["height"]}'][key] for c in centres])
            tot += float(vals @ w)
        return tot

    link = next(c for c in links['link_cases'] if c['kind'] == 'microwave')
    chain = dict(tx=link['tx_efficiency'], capture=link['capture'], rx=link['rx_efficiency'],
                 bus=link['bus_efficiency'])
    heat_per_w_no_air = population['earth_service']['heat_on_earth_per_w_delivered']
    load_mean = {}
    for f in BANDS_GHZ:
        towns_db = weighted(f, 'towns', 'beams_db', town_w)
        ground_db = ground_mean(f, 'beams_db')
        towns_loss, ground_loss = float(loss(towns_db)), float(loss(ground_db))
        load_mean[f'{f:g}_ghz'] = dict(
            towns_db=towns_db, towns_share_lost=towns_loss, ground_db=ground_db, ground_share_lost=ground_loss,
            towns_best_tile_db=weighted(f, 'towns', 'best_tile_db', town_w),
            ground_best_tile_db=ground_mean(f, 'best_tile_db'))
    a_t = load_mean['2.45_ghz']['towns_share_lost']
    a_g = load_mean['2.45_ghz']['ground_share_lost']
    a_t_best = float(loss(load_mean['2.45_ghz']['towns_best_tile_db']))
    a_g_best = float(loss(load_mean['2.45_ghz']['ground_best_tile_db']))
    intake = pair(lambda e: sky_tw[e] + ground_tw[e])
    absorbed_tw = pair(lambda e: sky_tw[e] * a_t / (1 - a_t) + ground_tw[e] * a_g / (1 - a_g))
    absorbed_best_tw = pair(lambda e: sky_tw[e] * a_t_best / (1 - a_t_best) + ground_tw[e] * a_g_best / (1 - a_g_best))
    a_mean = pair(lambda e: absorbed_tw[e] / (intake[e] + absorbed_tw[e]))
    a_best = pair(lambda e: absorbed_best_tw[e] / (intake[e] + absorbed_best_tw[e]))
    sheat = ledgers['heat']
    dim = sheat['factors']['dimming_percent_per_tw']
    with_computing = sheat['array']['with_close_computing']['heat_tw']
    link_effect = dict(
        chain=chain, dc_efficiency_without_air=link['dc_efficiency'],
        dc_efficiency_with_air=pair(lambda e: link['dc_efficiency'] * (1 - a_mean[e])),
        heat_per_w_delivered_without_air=heat_per_w_no_air,
        heat_per_w_delivered_with_air=pair(lambda e: heat_per_w_no_air / (1 - a_mean[e])),
        towns_heat_per_w=heat_per_w_no_air / (1 - a_t), ground_heat_per_w=heat_per_w_no_air / (1 - a_g),
        absorbed_in_the_air_tw=absorbed_tw, intake_tw=intake, share_of_beams_lost=a_mean,
        array_heat_with_close_computing_tw=dict(synthesis=with_computing,
                                                with_air=pair(lambda e: with_computing[e] + absorbed_tw[e])),
        added_dimming_percent=pair(lambda e: absorbed_tw[e] * dim[e]),
        rise_in_the_arrays_heat=pair(lambda e: absorbed_tw[e] / with_computing[e]),
        with_the_best_tiles=dict(towns_share_lost=a_t_best, ground_share_lost=a_g_best, share_of_beams_lost=a_best,
                                 absorbed_in_the_air_tw=absorbed_best_tw,
                                 heat_per_w_delivered=pair(lambda e: heat_per_w_no_air / (1 - a_best[e])),
                                 rise_in_the_arrays_heat=pair(lambda e: absorbed_best_tw[e] / with_computing[e])),
        rule=('the main figures take each beam from a tile drawn at random from the arc above the usable elevation; '
              'the best tiles take every beam from the highest tile in view, the floor of what the air takes'))

    # 2. Rain.
    fall = (G_MORRISON / MOON_SURFACE_GRAVITY) ** ((RAIN_B + 1) / 3)
    rain_coeffs = {f'{f:g}_ghz': p838(f) for f in BANDS_GHZ}
    table5 = {f'{f:g}': dict(model=list(p838(f)[k] for k in ('kh', 'ah', 'kv', 'av')), table=list(v))
              for f, v in P838_TABLE5.items()}
    validation = {f'{f:g}': dict(model=list(p838(f)[k] for k in ('kh', 'ah', 'kv', 'av')), itu=list(v))
                  for f, v in P838_VALIDATION.items()}
    rates = (1.0, 5.0, 10.0, 25.0, 50.0, 100.0)
    specific = {f'{f:g}_ghz': {f'{r:g}_mm_h': dict(earth_fall_db_km=rain_db_per_km(f, r, 1.0),
                                                   lunar_fall_db_km=rain_db_per_km(f, r, fall)) for r in rates}
                for f in BANDS_GHZ}
    storms = {k: ring_storms(v) for k, v in rings.items()}
    eq = storms['equator']
    freeze = col['freezing_level_km']
    e_ground = float(np.array([by_lat[f'{c:g}']['highest_tile_deg']['mean'] for c in centres]) @
                     sum(w for _, w in ground_w.values()))
    e_towns = float(np.array([by_lat[f'{c:g}']['highest_tile_deg']['mean'] for c in centres]) @ town_w)
    cases = {
        'p99_rate': dict(rate=eq['p99_mm_h'][0], top=eq['cloud_top_km']['median'], width=eq['width_km']['median']),
        'median_storm': dict(rate=eq['rain_max_mm_h']['median'], top=eq['cloud_top_km']['median'],
                             width=eq['width_km']['median']),
        'p90_storm': dict(rate=eq['rain_max_mm_h']['p90'], top=eq['cloud_top_km']['p90'], width=eq['width_km']['p90']),
        'largest_storm': dict(rate=eq['rain_max_mm_h']['max'], top=eq['cloud_top_km']['max'],
                              width=eq['width_km']['max'])}
    looks = lambda e: (('operating', e), ('usable_floor', USABLE_ELEVATION_DEG))
    fades = {}
    for name, c in cases.items():
        top = min(c['top'], freeze)
        row = dict(rate_mm_h=c['rate'], rain_top_km=top, width_km=c['width'])
        for who, h0, e in (('ground', 0.0, e_ground), ('towns', heights['towns'], e_towns)):
            row[who] = {}
            for label, el in looks(e):
                length = path_in_rain(h0, top, c['width'], el)
                bands = {}
                for f in BANDS_GHZ:
                    db = rain_db_per_km(f, c['rate'], fall) * length
                    bands[f'{f:g}_ghz'] = dict(db=db, share_lost=float(loss(db)))
                row[who][label] = dict(elevation_deg=el, path_km=length, **bands)
        fades[name] = row
    crossing = {}
    coverage_45 = float(np.mean([storms['45n']['coverage_above_1_mm_h'], storms['45s']['coverage_above_1_mm_h']]))
    for who, h0, e in (('ground', 0.0, e_ground), ('towns', heights['towns'], e_towns)):
        depth = min(eq['cloud_top_km']['median'], freeze) - h0
        crossing[who] = {}
        for label, el in looks(e):
            k = crossing_factor(depth, eq['width_km']['median'], el)
            crossing[who][label] = dict(elevation_deg=el, factor=k,
                                        storm_share_equator=eq['coverage_above_1_mm_h'] * k,
                                        above_p99_share_equator=0.01 * k, storm_share_45_deg=coverage_45 * k)
    cloud = [(h, d, n) for band in people['sky']['bands'].values()
             for h, d, n in zip(band['cloud_fraction_ring']['heights_km'], band['cloud_fraction_ring']['day'],
                                band['cloud_fraction_ring']['night'])]

    def cloud_on_path(h0):
        rows = [(h, d, n) for h, d, n in cloud if h >= h0]
        return dict(day=max(r[1] for r in rows), night=max(r[2] for r in rows),
                    day_height_km=max(rows, key=lambda r: r[1])[0])
    rain = dict(model='ITU-R P.838-3 equations 1-5, circular polarization (tilt 45 degrees)',
                cloud_on_the_path=dict(ground=cloud_on_path(0.0), towns=cloud_on_path(heights['towns']),
                                       rule='the largest share of time the equatorial ring holds cloud at any height '
                                            'the path crosses, from the ecology work\'s sky bands'),
                coefficients=rain_coeffs, check_table5=table5, check_validation_file=validation,
                lunar_fall_factor=fall,
                lunar_fall_rule=('at a sixth of the gravity raindrops fall about 2.9 times slower (CM1, '
                                 'V ~ g^((b+1)/3) '
                                 'with b = 0.8), so a lunar rain rate holds 2.9 times the water of the same rate on '
                                 'Earth; attenuation below 10 GHz goes with the water held, so P.838-3 is read at '
                                 'the Earth rate that holds the same water'),
                specific_attenuation=specific, storms_by_ring=storms, freezing_level_km=freeze,
                operating_elevation_deg=dict(ground=e_ground, towns=e_towns), fades=fades, crossing=crossing)

    # 3, continued: power per tile, apertures, beams and tiles.
    tile_side_km = ledger['requirements']['tile_side_m'] / 1e3
    clear_km = ledger['requirements']['clear_side_m'] / 1e3
    eta = aheat['design']['collector_conversion']
    heat_per_electric = sorted({c['conversion']: c['heat_per_electric'] for c in aheat['collectors']}.items())
    collector_heat = [heat_per_electric[1][1], heat_per_electric[0][1]]           # 30% then 20% conversion
    per_km2_electric = [SOLAR_CONSTANT * 1e6 * x / 1e9 for x in sorted(eta, reverse=True)]   # GW per km2
    per_km2_rf = [x * chain['tx'] for x in per_km2_electric]
    tile_rf = [x * clear_km ** 2 for x in per_km2_rf]
    tau = population['earth_service']['tau_1978']
    r_min, r_max = float(tiles['radius_km'].min()), float(tiles['radius_km'].max())
    rr = MOON_RADIUS / 1e3
    e_u = math.radians(USABLE_ELEVATION_DEG)
    d_min = r_min - rr
    d_max = math.sqrt(r_max ** 2 - (rr * math.cos(e_u)) ** 2) - rr * math.sin(e_u)
    crown_m2 = math.pi * (dist['basis']['sky']['envelope_radius_km'] * 1e3) ** 2
    patch_m2 = (RECEIVER_SIDE_KM * 1e3) ** 2
    apertures = {}
    for f in BANDS_GHZ:
        lam = SPEED_OF_LIGHT / (f * 1e9)
        side = lambda area, d: tau * lam * d * 1e3 / math.sqrt(area) / 1e3
        spot = lambda d: 2 * math.sqrt((tau * lam * d * 1e3) ** 2 / (clear_km * 1e3) ** 2 / math.pi)
        apertures[f'{f:g}_ghz'] = dict(town_crown_side_km=[side(crown_m2, d_min), side(crown_m2, d_max)],
                                       ground_patch_side_km=[side(patch_m2, d_min), side(patch_m2, d_max)],
                                       whole_tile_spot_m=[spot(d_min), spot(d_max)])
    peak = ap['intensities_w_m2']['rectenna_1978_peak']
    patch_power_gw = peak * patch_m2 / 1e9
    receivers_n = pair(lambda e: ground_tw[e] * 1e3 / patch_power_gw)
    rf_in_use = pair(lambda e: sky_tw[e] / (1 - a_t) + ground_tw[e] / (1 - a_g))
    rated = BEAM_RATING_SHARE * BEAMS_PER_RECEIVER                 # each receiver's beams together carry 1.5 loads
    rf_installed = pair(lambda e: rf_in_use[e] * rated)
    electric = pair(lambda e: rf_in_use[e] / chain['tx'])
    collector_km2 = [electric[0] * 1e3 / per_km2_electric[0], electric[1] * 1e3 / per_km2_electric[1]]
    transmitting = [rf_installed[0] * 1e3 / tile_rf[0] / TILE_UTILISATION,
                    rf_installed[1] * 1e3 / tile_rf[1] / TILE_UTILISATION]
    share = [t / fleet_tiles for t in transmitting]
    beams = pair(lambda e: (towns + receivers_n[e]) * BEAMS_PER_RECEIVER)
    in_orbit_tw = [electric[0] * (1 - chain['tx']) + electric[0] * collector_heat[0],
                   electric[1] * (1 - chain['tx']) + electric[1] * collector_heat[1]]
    orbit_w_m2 = aheat['placement']['moon_W_m2_per_TW_released_in_orbit']
    glow = ledgers['heat']['fleet_glow']['w_m2']
    dimmer = ledger['electricity']['dimmer_electric_TW']
    town_beam_mw = pair(lambda e: sky_tw[e] * 1e6 / towns / BEAMS_PER_RECEIVER)
    installed_km2 = [t * clear_km ** 2 for t in transmitting]
    lam_245 = SPEED_OF_LIGHT / (BANDS_GHZ[0] * 1e9)
    tau_tile = [math.sqrt(clear_km ** 2 * 1e6 * crown_m2) / (lam_245 * d * 1e3) for d in (d_max, d_min)]
    r_mid = 20000.0
    n_mid = math.sqrt(MOON_GM / (r_mid * 1e3) ** 3) + OMEGA_SUN_FRAME
    psi = lambda e: math.acos(rr / r_mid * math.cos(math.radians(e))) - math.radians(e)
    overhead_pass = dict(horizon_h=2 * psi(0.0) / n_mid / 3600.0,
                         usable_h=2 * psi(USABLE_ELEVATION_DEG) / n_mid / 3600.0)
    per_ring = fleet_tiles / tiles['rings']
    tile_gap_s = [2 * math.pi / float(m) / per_ring for m in (tiles['mean_motion'].max(), tiles['mean_motion'].min())]
    typical = [by_lat[f'{c:g}']['typical_usable_pass_h'] for c in centres]
    usable_least = {f'{c:g}': by_lat[f'{c:g}']['tiles_above_usable']['least'] for c in centres}
    ring_out = dict(
        stack=dict(rings=tiles['rings'], sampled_rings=len(tiles['radius_km']), radius_km=[r_min, r_max],
                   tilt_deg=[float(tiles['tilt_deg'].min()), float(tiles['tilt_deg'].max())], fleet_tiles=fleet_tiles,
                   orbit_days=[2 * math.pi / float(m) / JULIAN_DAY for m in (tiles['mean_motion'].max(),
                                                                             tiles['mean_motion'].min())],
                   ring_plane_from_lunar_pole_deg=EPS_DEG),
        by_latitude={k: {kk: vv for kk, vv in v.items() if kk != 'absorption'} for k, v in by_lat.items()},
        middle_ring_alone_hidden_within_deg_of_its_poles=math.degrees(math.asin(rr / 20000.0)),
        tile=dict(side_km=tile_side_km, clear_side_km=clear_km, sunlight_gw=SOLAR_CONSTANT * clear_km ** 2 * 1e6 / 1e9,
                  electric_gw_per_km2=per_km2_electric, rf_gw_per_km2=per_km2_rf, whole_tile_rf_gw=tile_rf,
                  dimmer_share_rf_gw=[SOLAR_CONSTANT * clear_km ** 2 * 1e6 * ledger['requirements']['dimmer_share']
                                      * x * chain['tx'] / 1e9 for x in sorted(eta, reverse=True)]),
        tau=tau, distance_km=[d_min, d_max], apertures=apertures,
        beams=dict(per_receiver=BEAMS_PER_RECEIVER, rating_share=BEAM_RATING_SHARE, towns=towns,
                   ground_receivers=receivers_n, patch_side_km=RECEIVER_SIDE_KM, patch_gw=patch_power_gw,
                   total=beams, town_beam_mw_in_use=town_beam_mw, ground_beam_gw_in_use=patch_power_gw /
                   BEAMS_PER_RECEIVER),
        rf_in_use_tw=rf_in_use, rf_installed_tw=rf_installed, electricity_at_the_tiles_tw=electric,
        dimmer_electric_tw=dimmer, collector_km2_in_use=collector_km2,
        transmitting_tiles=transmitting, share_of_fleet=share, collector_and_aperture_km2_installed=installed_km2,
        whole_tile_tau_on_a_crown=tau_tile, overhead_pass_middle_ring=overhead_pass, tiles_per_ring=per_ring,
        seconds_between_tiles_on_a_ring=tile_gap_s,
        handover_every_h_by_latitude={f'{c:g}': t / BEAMS_PER_RECEIVER for c, t in zip(centres, typical)},
        beams_per_transmitting_tile=[beams[0] / transmitting[0], beams[1] / transmitting[1]],
        transmitting_tiles_in_usable_view_least={k: [v * s for s in share] for k, v in usable_least.items()},
        small_collector_variant=dict(collector_share_of_tile=SMALL_COLLECTOR_SHARE,
                                     rf_gw_per_tile=[x * SMALL_COLLECTOR_SHARE for x in tile_rf],
                                     transmitting_tiles=[t / SMALL_COLLECTOR_SHARE for t in transmitting],
                                     share_of_fleet=[x / SMALL_COLLECTOR_SHARE for x in share],
                                     beams_per_transmitting_tile=[beams[0] * SMALL_COLLECTOR_SHARE / transmitting[0],
                                                                  beams[1] * SMALL_COLLECTOR_SHARE / transmitting[1]],
                                     transmitting_tiles_in_usable_view_least={
                                         k: [v * x / SMALL_COLLECTOR_SHARE for x in share]
                                         for k, v in usable_least.items()}),
        heat_released_in_orbit_tw=in_orbit_tw,
        moon_w_m2_from_that_heat=[x * orbit_w_m2 for x in in_orbit_tw],
        share_of_the_fleets_glow=[in_orbit_tw[0] * orbit_w_m2 / glow[0], in_orbit_tw[1] * orbit_w_m2 / glow[1]])

    # 4. The grid.
    settled_land = dist['moon']['cprime_with_blend']['1_kg']['wild_30_hybrid_20']['full']['land']['settled_km2']
    beam_loss_by_bin = {}
    for c in centres:
        v = by_lat[f'{c:g}']['absorption']
        beam_loss_by_bin[f'{c:g}'] = dict(lowland=v['2.45_ghz_lowland']['beams_share_lost'],
                                          highland=v['2.45_ghz_highland']['beams_share_lost'],
                                          no_usable_tile=by_lat[f'{c:g}']['share_of_time_with_no_usable_tile'])
    always = [c for c in centres if by_lat[f'{c:g}']['share_of_time_with_no_usable_tile'] == 0.0]
    edge_lat = max(always)
    lines = {}
    for c in centres:
        if c <= edge_lat:
            continue
        km = (c - edge_lat) * KM_PER_DEG
        line = [min(1.0, x * km / 1000.0 + 2 * HVDC_CONVERTER) for x in HVDC_PER_1000_KM]
        ref = beam_loss_by_bin[f'{edge_lat:g}']['highland']
        via_line = [1 - (1 - ref) * (1 - x) for x in line]
        lines[f'{c:g}'] = dict(km=km, line_loss=line, via_line_loss=via_line,
                               beam_loss_highland=beam_loss_by_bin[f'{c:g}']['highland'],
                               beam_loss_lowland=beam_loss_by_bin[f'{c:g}']['lowland'],
                               line_better=bool(via_line[0] < beam_loss_by_bin[f'{c:g}']['highland']))
    spots = apertures['2.45_ghz']['whole_tile_spot_m']
    buffer_km = [edge_buffer_km(x, peak, PUBLIC_LIMIT_W_M2) for x in spots]
    closed_patch = [(RECEIVER_SIDE_KM + 2 * b) ** 2 for b in buffer_km]
    closed_land = [receivers_n[0] * closed_patch[0], receivers_n[1] * closed_patch[1]]
    exposure = dict(public_limit_w_m2=PUBLIC_LIMIT_W_M2, patch_level_w_m2=peak, edge_buffer_km=buffer_km,
                    closed_km2_per_patch=closed_patch, closed_land_km2=closed_land,
                    closed_over_settled_land=pair(lambda e: closed_land[e] / settled_land[e]),
                    rule='a patch flat-filled at the 1978 peak with edges blurred by a whole tile\'s spot; the '
                         'closed ground runs out to where the beam falls below the public level')
    load_by_bin = sum(w for _, w in ground_w.values())
    poleward_share = float(sum(load_by_bin[i] for i, c in enumerate(centres) if c > edge_lat))
    poleward_line_share = float(sum(load_by_bin[i] for i, c in enumerate(centres)
                                    if c > edge_lat and lines[f'{c:g}']['line_better']))
    settled = settled_land
    spacing = pair(lambda e: math.sqrt(settled[e] / receivers_n[e]))
    receiver_land = ap['ground']['receiver_km2_at_1978_peak']
    band_load = {}
    for name, (g, w) in ground_w.items():
        band_load[name] = dict(lat_deg=list(g['lat']), height=g['height'], billion=g['billion'],
                               intake_tw=pair(lambda e: ground_tw[e] * g['billion'] / 4.0),
                               receiver_km2=pair(lambda e: ground_tw[e] * g['billion'] / 4.0 * 1e12 / peak / 1e6))
    grid = dict(rule=('receivers at the loads, the ground\'s intake spread as its people; poleward of the last band '
                      'that always sees a usable tile, the loads take power by line when that loses less than '
                      'their own low beams'),
                composition_billion=composition, groups=band_load, beam_loss_by_latitude_2_45=beam_loss_by_bin,
                always_usable_to_deg=edge_lat + 2.5, poleward_lines=lines,
                load_share_poleward=poleward_share, load_share_taking_lines=poleward_line_share,
                receiver_spacing_km=spacing, receiver_km2=receiver_land, settled_km2=settled,
                receivers_over_settled_land=pair(lambda e: receiver_land[e] / settled[e]), exposure=exposure,
                backup=dict(rule='each receiver tied to two neighbours by lines rated for its whole load',
                            line_gw=patch_power_gw,
                            loss_through_a_neighbour=[2 * HVDC_CONVERTER + x * s / 1000.0
                                                      for x, s in zip(HVDC_PER_1000_KM, spacing)]),
                hvdc=dict(per_1000_km=list(HVDC_PER_1000_KM), per_converter=HVDC_CONVERTER))

    # 5. Reserves.
    skytowns = ledgers['heights']['sky_towns']
    share_reef = skytowns['reef_share_of_cover']
    widths = skytowns['reef_width_km']
    drift = skytowns['relative_drift_m_s']
    sky_light = people['sky']['bands']['storm_layer']
    top_of_air = ledger['requirements']['climate_top_of_air_W_m2']
    sun = [sky_light['light_at_bottom']['direct_sun_overhead'] * top_of_air,
           sky_light['light_at_top']['direct_sun_overhead'] * top_of_air]
    reef = {}
    layer = REEF_LAYER_KM[1] - REEF_LAYER_KM[0]
    for who, width_km, e in (('towns', 2 * dist['basis']['sky']['envelope_radius_km'], e_towns),
                             ('ground', RECEIVER_SIDE_KM, e_ground)):
        rows = {}
        for label, el in (('operating', e), ('usable_floor', USABLE_ELEVATION_DEG)):
            along = width_km + layer / math.tan(math.radians(el))
            mean_width = 2 * (width_km + along) / math.pi
            per_cover = {}
            for cover in skytowns['reef_passes_per_town_day']:
                rates, blocked, hours = [], [], []
                for d in widths:
                    n = float(cover) * share_reef / (math.pi * (d / 2) ** 2)
                    r = d / 2 + 0.05
                    area = width_km * along + 2 * (width_km + along) * r + math.pi * r * r   # the strip grown by a reef
                    blocked.append(1.0 - math.exp(-n * area))
                    for v in drift:
                        rates.append(n * (mean_width + d + 0.1) * v * JULIAN_DAY / 1e3)
                        hours.append((mean_width + d) * 1e3 / v / 3600.0)
                per_cover[cover] = dict(crossings_per_beam_day=[min(rates), max(rates)],
                                        share_of_time_blocked_without_switching=[min(blocked), max(blocked)],
                                        crossing_hours=[min(hours), max(hours)],
                                        switches_per_receiver_day=[BEAMS_PER_RECEIVER * min(rates),
                                                                   BEAMS_PER_RECEIVER * max(rates)])
            rows[label] = dict(elevation_deg=el, strip_km=[width_km, along], mean_width_km=mean_width,
                               by_cover=per_cover)
        reef[who] = rows
    s_e = math.sin(math.radians(e_towns)) ** 2
    offsets = [math.degrees((min(widths) + 1.0) * s_e / REEF_ABOVE_TOWN_KM[1]),
               math.degrees((max(widths) + 1.0) * s_e / REEF_ABOVE_TOWN_KM[0])]
    town_intensity = ap['sky_towns']['w_m2']
    reef_dose = dict(town_beam_w_m2=[x / BEAMS_PER_RECEIVER for x in town_intensity],
                     town_beam_rated_w_m2=[x * BEAM_RATING_SHARE for x in town_intensity],
                     ground_beam_w_m2=peak / BEAMS_PER_RECEIVER, ground_beam_rated_w_m2=peak * BEAM_RATING_SHARE,
                     overhead_sun_in_the_reef_layer_w_m2=sun)
    lam = TILE_FAULTS_PER_YEAR
    tau_y = REASSIGN_S / (HOURS_PER_YEAR * 3600.0)
    per_receiver = BEAMS_PER_RECEIVER * lam * (BEAMS_PER_RECEIVER - 1) * lam * tau_y
    faults = dict(tile_faults_per_year=lam, reassign_s=REASSIGN_S,
                  beam_losses_per_receiver_year=BEAMS_PER_RECEIVER * lam,
                  double_losses_per_receiver_year=per_receiver,
                  double_losses_per_year_moonwide=pair(lambda e: per_receiver * (towns + receivers_n[e])),
                  after_one_loss='the other two beams carry the whole load at their rating',
                  after_two_losses='the last beam carries half the load until a fresh beam is placed',
                  transmitting_tiles_failing_per_day=[t * lam / JULIAN_YEAR_DAYS for t in transmitting])
    parts = sheat['array']['heat_by_component_tw']
    h_w = heat_per_w_no_air
    aloft = ledgers['choices']['aloft_billion'] / ledgers['choices']['moon_billion']
    bp = blend['blend'][str(ledgers['choices']['moon_billion'])]['wild_30_hybrid_20']['power_tw']
    cultures = pair(lambda e: bp['power_made_food'][e] + bp['microbial'][e])
    lamps = pair(lambda e: bp['lamps'][e])
    computing = [x / h_w for x in sheat['array']['close_computing_heat_tw']]
    crit = dict(
        residents_aloft=pair(lambda e: parts['residents'][e] / h_w * aloft * CRITICAL_SHARE['residents']),
        residents_ground=pair(lambda e: parts['residents'][e] / h_w * (1 - aloft) * CRITICAL_SHARE['residents']),
        sky_control=pair(lambda e: parts['sky_control'][e] / h_w * CRITICAL_SHARE['sky_control']),
        water=pair(lambda e: parts['water'][e] / h_w * CRITICAL_SHARE['water']),
        food_cultures=pair(lambda e: cultures[e] * CRITICAL_SHARE['cultures']),
        lamp_farms=pair(lambda e: lamps[e] * CRITICAL_SHARE['lamps']),
        close_computing=pair(lambda e: computing[e] * CRITICAL_SHARE['close_computing']))
    total_crit = pair(lambda e: sum(v[e] for v in crit.values()))
    delivered = pair(lambda e: (sky_tw[e] + ground_tw[e]) / h_w)
    towns_crit = pair(lambda e: crit['residents_aloft'][e] + crit['sky_control'][e])
    ground_crit = pair(lambda e: total_crit[e] - towns_crit[e])
    hours = EMERGENCY_DAYS * 24.0
    h2_lift = dist['moon']['cprime_with_blend']['1_kg']['wild_30_hybrid_20']['full']['sky']['hydrogen_gt']
    reserve_lift = volume['assumptions']['lift_allocation_shares']['reserve']
    kwh_m3 = tower['night_storage']['delivered_kwh_per_m3']
    town_h2_mt = pair(lambda e: towns_crit[e] * 1e9 * hours / (FUEL_CELL * H2_LHV_KWH_KG) / 1e9)
    ground_twh = pair(lambda e: ground_crit[e] * hours)
    ring_fault = dict(
        days=EMERGENCY_DAYS, critical_share=CRITICAL_SHARE, critical_tw=crit, critical_total_tw=total_crit,
        delivered_total_tw=delivered, critical_share_of_delivered=pair(lambda e: total_crit[e] / delivered[e]),
        energy_twh=pair(lambda e: total_crit[e] * hours), per_day_twh=pair(lambda e: total_crit[e] * 24.0),
        towns=dict(critical_tw=towns_crit, hydrogen_mt=town_h2_mt,
                   share_of_lifting_hydrogen=pair(lambda e: town_h2_mt[e] / (h2_lift[e] * 1e3)),
                   share_of_reserve_lift=pair(lambda e: town_h2_mt[e] / (h2_lift[e] * 1e3) / reserve_lift)),
        ground=dict(critical_tw=ground_crit, energy_twh=ground_twh,
                    hydrogen_mt=pair(lambda e: ground_twh[e] * 1e9 / (FUEL_CELL * H2_LHV_KWH_KG) / 1e9),
                    electricity_to_make_it_twh=pair(lambda e: ground_twh[e] / H2_ROUND_TRIP),
                    pumped_hydro_km3_at_korolevs_head=pair(lambda e: ground_twh[e] / kwh_m3),
                    korolev_level_m=pair(lambda e: ground_twh[e] / kwh_m3 * 1e9
                                         / (tower['site']['lakes'][1]['area_km2'] * 1e6)),
                    battery_bridge_twh=pair(lambda e: total_crit[e] * BRIDGE_HOURS),
                    battery_bridge_mt=[total_crit[0] * BRIDGE_HOURS * 1e9 / BATTERY_KWH_KG[1] / 1e9,
                                       total_crit[1] * BRIDGE_HOURS * 1e9 / BATTERY_KWH_KG[0] / 1e9]),
        food=dict(rule='the cities hold about a month of food (first screen), so production stops while the '
                       'living cultures are kept alive'),
        towns_hydrogen_unremade_share=[x * hours / HOURS_PER_YEAR
                                       for x in dist['basis']['sky']['hydrogen_loss_per_year']])
    reserves = dict(storms=dict(fades=fades, rating_headroom_db=float(10 * math.log10(BEAMS_PER_RECEIVER
                                                                                       * BEAM_RATING_SHARE))),
                    reefs=dict(beam_columns=reef, switch_offset_deg=offsets, what_a_reef_takes=reef_dose),
                    single_tiles=faults, ring_fault=ring_fault)

    # 6. The summit metropolis.
    met = composition['metropolis']
    moon_b = ledgers['choices']['moon_billion']
    res_d = pair(lambda e: parts['residents'][e] / h_w * met / moon_b)
    comp_d = pair(lambda e: computing[e] * met / moon_b)
    food_d = pair(lambda e: bp['total'][e] * met / moon_b)
    total_d = pair(lambda e: res_d[e] + comp_d[e] + food_d[e])
    a_s = summit_view['absorption']['2.45_ghz_summit']['beams_share_lost']
    arriving = pair(lambda e: total_d[e] * h_w)
    rec_km2 = pair(lambda e: arriving[e] * 1e12 / peak / 1e6)
    city_r = 32.0
    annulus = pair(lambda e: math.sqrt(city_r ** 2 + rec_km2[e] / math.pi) - city_r)
    crit_s = pair(lambda e: res_d[e] * CRITICAL_SHARE['residents'] + comp_d[e] * CRITICAL_SHARE['close_computing']
                  + food_d[e] * (cultures[e] * CRITICAL_SHARE['cultures'] + lamps[e] * CRITICAL_SHARE['lamps'])
                  / bp['total'][e])
    store = pair(lambda e: crit_s[e] * hours)
    summit = dict(
        people_billion=met, latitude_deg=summit_lat, ground_km=heights['summit'],
        brief_gw=screen['heat']['human_use']['metropolis_gw'],
        delivered_tw=dict(residents=res_d, close_computing=comp_d, food=food_d, total=total_d),
        intake_tw=arriving, receiver_km2=rec_km2, patches=pair(lambda e: arriving[e] * 1e3 / patch_power_gw),
        receiver_ring_width_km=annulus, city_radius_km=city_r,
        view=summit_view, beams_share_lost=a_s, heat_per_w_delivered=h_w / (1 - a_s),
        critical_tw=crit_s, store_twh=store,
        korolev_km3=pair(lambda e: store[e] / kwh_m3),
        hydrogen_mt=pair(lambda e: store[e] * 1e9 / (FUEL_CELL * H2_LHV_KWH_KG) / 1e9),
        night_store_without_beams_twh=screen['night']['metropolis_night_twh'])

    synthesis_check = dict(
        absorption_estimate_share=[0.05, 0.08],
        absorption_computed_share=dict(zenith_sea_level=float(loss(tables[(2.45, 'lowland')][0][-1])),
                                       zenith_10_km=float(loss(tables[(2.45, 'towns')][0][-1])),
                                       beams_load_weighted=a_mean, best_tiles_load_weighted=a_best),
        heat_per_w=dict(synthesis=heat_per_w_no_air, with_air=link_effect['heat_per_w_delivered_with_air']),
        dimmer_against_tile_electricity=dict(dimmer_tw=dimmer, needed_tw=electric),
        polar=dict(synthesis='the ring below the horizon within about 5 degrees of its plane\'s poles',
                   computed='the stack\'s tilted rings keep tiles above every horizon',
                   tiles_above_the_horizon_at_the_pole=by_lat['87.5']['tiles_above_horizon'],
                   highest_tile_at_the_pole_deg=by_lat['87.5']['highest_tile_deg']),
        rain_10_km_50_100_mm_h_2_45_db=dict(earth_fall=[10 * rain_db_per_km(2.45, r, 1.0) for r in (50.0, 100.0)],
                                            lunar_fall=[10 * rain_db_per_km(2.45, r, fall) for r in (50.0, 100.0)]),
        rain_10_km_50_mm_h_5_8_share_lost=dict(earth_fall=float(loss(10 * rain_db_per_km(5.8, 50.0, 1.0))),
                                               lunar_fall=float(loss(10 * rain_db_per_km(5.8, 50.0, fall)))))

    assumptions = dict(
        design_guesses=dict(
            usable_elevation_deg=dict(value=USABLE_ELEVATION_DEG, basis='beams planned from 30 degrees up, where the '
                                      'air takes about 11% at sea level against 6% overhead'),
            beams_per_receiver=dict(value=BEAMS_PER_RECEIVER, rating_share=BEAM_RATING_SHARE,
                                    basis='three beams from different parts of the ring, each rated for half the load, '
                                          'switched make-before-break'),
            receiver_side_km=dict(value=RECEIVER_SIDE_KM, basis='a flat-filled patch a few kilometres across, as the '
                                  'domain synthesis sizes the ground receivers'),
            storage_days=dict(value=EMERGENCY_DAYS, basis='the critical loads carried through a fault across the '
                              'whole ring for a week, after the 7-day fuel that Earth keeps for safety-related '
                              'emergency power (US NRC Regulatory Guide 1.137)'),
            tile_utilisation=dict(value=TILE_UTILISATION, basis='half a transmitting tile\'s rating in use on average'),
            small_collector_share=dict(value=SMALL_COLLECTOR_SHARE, basis='the variant with collectors over a tenth '
                                       'of each transmitting tile'),
            transmitter=dict(value='a sandwich of collector and phased array over a tile\'s clear face, sun-pointing '
                                   'collectors', basis='each beam may use the whole aperture'),
            tile_faults_per_year=dict(value=TILE_FAULTS_PER_YEAR, basis='one whole-tile outage a decade per tile'),
            reassign_s=dict(value=REASSIGN_S),
            fuel_cell=dict(value=FUEL_CELL, basis='electricity per unit of hydrogen\'s lower heating value'),
            battery_kwh_kg=dict(value=list(BATTERY_KWH_KG), basis='pack specific energy'),
            bridge_hours=dict(value=BRIDGE_HOURS,
                              basis='batteries carry the critical loads while hydrogen plant starts'),
            critical_shares=dict(value=CRITICAL_SHARE, basis='residents: ventilation, water and sewage pumps, '
                                 'refrigeration of food and medicine, safety lighting, communications and health care; '
                                 'sky control: towns drift with the wind and steer only to keep clear; water: all; '
                                 'cultures: the power-made food\'s living cultures kept alive at a tenth of their '
                                 'feed; lamps: the lamp-lit farms dark, with air and water kept moving; close '
                                 'computing: the millisecond control kept on the Moon; the towns\' hydrogen remade '
                                 'and trips shed'),
            dry_highland_height_km=dict(value=DRY_HIGHLAND_KM, basis='the dry highland cities on ground 4 km or more '
                                        'above the sea, at the low end of it'),
            ground_loads_follow_people=dict(value=True,
                                            basis='food plants, lamp farms and close computing at the cities'),
            composition=dict(value='the middle of each C\' pair, with the ground\'s 4 billion'),
            reef_above_town_km=dict(value=list(REEF_ABOVE_TOWN_KM),
                                    basis='the synthesis\'s reefs 5-10 km above a town')),
        estimates=dict(
            lunar_fall=dict(value=fall, basis='attenuation with the water held, Earth\'s drop spectra'),
            rain_cell=dict(value='uniform at the storm\'s peak surface rate from the ground to its top or the 0 C '
                                 'level, as wide as the storm'),
            circular_rings=dict(value='the layout\'s forced eccentricities of 0.11-0.16 are left out; distances '
                                      'swing by that share about each orbit')),
        published=dict(hvdc=dict(per_1000_km=list(HVDC_PER_1000_KM), per_converter=HVDC_CONVERTER,
                                 source='US EIA 2018'),
                       public_limit_w_m2=dict(value=PUBLIC_LIMIT_W_M2, source='ICNIRP 2020'),
                       emergency_fuel_days=dict(value=EMERGENCY_DAYS, source='US NRC Regulatory Guide 1.137'),
                       hydrogen_round_trip=dict(value=H2_ROUND_TRIP, source='Headley & Schoenung 2020'),
                       hydrogen_lhv_kwh_kg=dict(value=H2_LHV_KWH_KG)))

    return dict(
        schema=SCHEMA,
        producer=dict(domain='provisioning', files={f: sha256(f) for f in CODE},
                      inputs={p: sha256(p) for p in INPUTS}, constants=constants_used(list(CODE))),
        evidence=EVIDENCE, reading_rule=READING_RULE, assumptions=sig(assumptions),
        absorption=sig(dict(earth_checks=checks, column_check=column_check(table_rows), air_check=air_check,
                            **absorb, by_latitude={k: v['absorption'] for k, v in by_lat.items()},
                            load_weighted=load_mean, link=link_effect)),
        rain=sig(rain), ring=sig(ring_out), grid=sig(grid), reserves=sig(reserves), summit=sig(summit),
        synthesis_check=sig(synthesis_check))


def _column_lnp():
    """Heights (km) and log pressure (Pa) of the design column, for comparing with other profiles."""
    air = th.earthlike_air(COLUMN['dry_pressure_pa'], COLUMN['co2_ppm'])
    spec = th.ColumnSpec(th.MOON, air, COLUMN_TS_K, COLUMN['stratosphere_k'], th.manabe_wetherald(),
                         troposphere_levels=COLUMN['levels'][0], upper_levels=COLUMN['levels'][1])
    col = th.build_column(spec)
    return col['z_m'] / 1e3, np.log(col['p_pa'])


def dumps(product) -> str:
    return json.dumps(product, indent=1, allow_nan=False) + '\n'


def main(argv=None) -> int:
    product = results()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(dumps(product))
    a = product['absorption']
    print(json.dumps(dict(zenith=a['by_elevation']['2.45_ghz']['lowland']['90_deg'], load=a['load_weighted'],
                          link=a['link']), indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
