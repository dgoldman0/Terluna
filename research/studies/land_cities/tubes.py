"""Lava tubes on the Open Moon: the candidates from the literature, how wide a roof spans, where each candidate lies
against the 28% atlas, the water table it meets, what keeping it dry takes, what it holds, and its two lives, as the
first phase's shelter and as a lasting city.

The inventory is a literature compilation (sources.json gives each source and how far it was read). Everything
else is closed-form physics, grid lookups and one small solve: Darcy flow with lunar gravity, the water table as a
Dupuit-Forchheimer mound in two dimensions on the atlas (groundwater.py), seepage through a sealing ring, the levels
a vaulted section holds, conduction into the roof and view factors. No structural or radiation-transport model runs
here.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.special import erfcinv

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

from research.studies.land_cities import grid, groundwater as gw

WATER_DENSITY = 1000.0          # kg/m3, fresh water
SECONDS_PER_YEAR = JULIAN_YEAR_DAYS * JULIAN_DAY
GRAVITY_RATIO = MOON_SURFACE_GRAVITY / STANDARD_GRAVITY

# ---------------------------------------------------------------------------------------------------------------
# Rock and water (values with their sources in sources.json; ranges are [low, central, high]).
BASALT_BULK_DENSITY = (3010.0, 3150.0, 3270.0)      # kg/m3, mare basalt samples (Kiefer et al. 2012, LPSC table)
# Hydraulic conductivity of basalt aquifers on Earth (m/s): Columbia River basalt bulk low end (Burns et al. 2016,
# after Hansen et al. 1994); Kahle et al. 2011's regional median for the basalt units (70 ft/d); Hawaii's dike-free
# lava flows at the top of their usual range (5,000 ft/d; Hunt 1996).
K_EARTH = (3.0e-7, 2.5e-4, 1.8e-2)
RECHARGE_SHARE = (0.1, 0.36, 0.5)                   # of rainfall: Hawaii's 10-50%, Oahu's island-wide 36%
SPECIFIC_YIELD = (0.04, 0.07, 0.10)                 # Hawaiian basalt, "most of the values ... about 5 to 10 percent"
BASALT_THICKNESS_M = (100.0, 740.0, 1620.0)         # nearside mare basalt, 1-sigma low, mean, high (Gong et al. 2016)
ROCK_DIFFUSIVITY = (0.5e-6, 1.0e-6)                 # m2/s, a design range for basalt's thermal diffusivity (guess)

# Sealing. A grouted ring of the stated thickness carries the head; Earth practice reaches K of about 1e-8 m/s
# (Katibeh & Aalianvari 2012: "the practical minimum for tunnels is around 10^-6 cm/s").
GROUTED_K_EARTH = (1e-9, 1e-8, 1e-7)                # m/s on Earth; the same rock passes a sixth at lunar gravity
GROUTED_RING_M = 5.0
PUMP_EFFICIENCY = 0.7                               # the inhabited-volume product's pump efficiency
# Allowable leakage into tunnels and underground public spaces (NASEM 2024): 0.002 and 0.001 gal/ft2/day.
LEAKAGE_ALLOWED_L_M2_DAY = dict(tunnel=0.0815, public_space=0.0407)
# ITA (1988): waterproofing sheets where the water stands more than 15 m above the crown (an Earth head).
ITA_MEMBRANE_HEAD_EARTH_M = 15.0
# A tunnel-waterproofing guide's exposure classes for mined tunnels by water head on Earth (Sika): below 5 m low;
# 5-30 m moderate, where membrane systems serve; 30-60 m high, membranes with watertight concrete; above 60 m very
# high, double-layer actively controlled membranes with watertight concrete.
MEMBRANE_CLASS_HEADS_EARTH_M = (5.0, 30.0, 60.0)

# Capacity.
STOREY_M = (4.0, 5.0)                               # floor to floor inside a tube
USABLE_SHARE = (0.5, 0.7)                           # floor kept for homes, work and services; the rest open (guess)
FLOOR_M2_PER_PERSON = (35.0, 70.0, 140.0)           # half, the metropolis brief's 70, and twice

# Shielding: Earth's atmosphere as a shield.
EARTH_AIR_COLUMN_KG_M2 = 101325.0 / STANDARD_GRAVITY


def r(x, n=3):
    if x is None or isinstance(x, (bool, int, str)):
        return x
    x = float(x)
    if x == 0.0 or not math.isfinite(x):
        return x
    return float(f'{x:.{n}g}')


def rr(seq, n=3):
    return [r(x, n) for x in seq]


def metres(x):
    """A height to the decimetre (heights from the solve carry no meaning below it)."""
    return round(float(x), 1) + 0.0


# ---------------------------------------------------------------------------------------------------------------
# The inventory. Longitudes east, 0-360. Sources by key in sources.json.

# The mare pits of the LROC Lunar Pit Atlas (Wagner & Robinson 2022's 21 pits outside impact melt, the 16 the atlas
# calls mare), with the two near-mare pits Wagner & Robinson 2014 counted among their eight mare pits. inner and
# funnel are mouth diameters (m), depth the atlas's inner depth (m), floor_m the floor's depth below the
# surrounding surface where the sources give it, void_m the void seen beyond the wall at floor level and its height.
PITS = [
    dict(key='marius_hills_pit', name='Marius Hills pit', lat=14.0917, lon=303.2299, inner=(49, 58), funnel=(79, 92),
         depth=40, floor_m=(45, 60), void=dict(reach_m=12, ceiling_m=17),
         setting='in a sinuous rille on the Marius Hills volcanic plateau, Oceanus Procellarum',
         sources=('lroc_pit_atlas', 'ashley_2011', 'haruyama_2009', 'wagner_robinson_2014')),
    dict(key='west_marius_hills_pit', name='West Marius Hills pit', lat=13.5507, lon=301.8267, inner=(47, 47),
         funnel=(70, 95), depth=16, floor_m=(16, 16), void=None,
         setting='Marius Hills; a linear collapse 2 km NW may mark a partly intact tube', sources=('lroc_pit_atlas',)),
    dict(key='mare_tranquillitatis_pit', name='Mare Tranquillitatis pit', lat=8.3355, lon=33.222, inner=(88, 100),
         funnel=(140, 146), depth=105, floor_m=(125, 175), void=dict(reach_m=20, ceiling_m=60),
         setting='near the centre of the Cauchy shield volcano, Mare Tranquillitatis; radar shows a conduit '
                 'running west',
         sources=('lroc_pit_atlas', 'ashley_2011', 'wagner_robinson_2014', 'carrer_2024', 'horvath_2022')),
    dict(key='sw_mare_tranquillitatis_pit', name='Southwest Mare Tranquillitatis pit (with its inner pit B)',
         lat=4.1438, lon=24.6871, inner=(26, 32), funnel=(80, 100), depth=25, floor_m=(25, 38), void=None,
         setting='Mare Tranquillitatis, 2.3 km from a wrinkle ridge', sources=('lroc_pit_atlas',)),
    dict(key='mare_ingenii_pit', name='Mare Ingenii pit', lat=-35.9494, lon=166.0559, inner=(71, 104),
         funnel=(157, 180), depth=55, floor_m=(75, 100), void=dict(reach_m=20, ceiling_m=None),
         setting='Mare Ingenii in the South Pole-Aitken basin, far side, near a swirl',
         sources=('lroc_pit_atlas', 'haruyama_2010', 'wagner_robinson_2014')),
    dict(key='runge_pit', name='Runge pit (Mare Smythii)', lat=-2.7024, lon=86.78, inner=(12, 14), funnel=(28, 31),
         depth=5, floor_m=(5, 6), void=None, setting='Mare Smythii, in Runge crater', sources=('lroc_pit_atlas',)),
    dict(key='lacus_mortis_pit', name='Lacus Mortis pit', lat=44.9608, lon=25.6119, inner=(110, 140),
         funnel=(210, 280), depth=60, floor_m=(100, 115), void=dict(reach_m=None, ceiling_m=None),
         setting='Lacus Mortis, beside the Rima Burg graben; a ramp reaches the floor',
         sources=('lroc_pit_atlas', 'wagner_robinson_2014')),
    dict(key='compton_pit', name='Compton pit', lat=56.2247, lon=106.1955, inner=(92, 110), funnel=(140, 230),
         depth=35, floor_m=(35, 40), void=None, setting='mare on the floor of Compton crater, far side',
         sources=('lroc_pit_atlas',)),
    dict(key='north_procellarum_1_pit', name='North Procellarum 1 pit', lat=35.4097, lon=314.3602, inner=(108, 157),
         funnel=(175, 195), depth=54, floor_m=(54, 58), void=None, setting='northern Oceanus Procellarum',
         sources=('lroc_pit_atlas',)),
    dict(key='north_procellarum_2_pit', name='North Procellarum 2 pit', lat=35.3433, lon=314.3447, inner=(155, 160),
         funnel=(170, 175), depth=43, floor_m=(43, 45), void=None, setting='2 km from North Procellarum 1',
         sources=('lroc_pit_atlas',)),
    dict(key='sinus_iridum_pit', name='Sinus Iridum pit', lat=45.623, lon=331.189, inner=(32, 66), funnel=(125, 145),
         depth=16, floor_m=(16, 21), void=None, setting='Sinus Iridum', sources=('lroc_pit_atlas',)),
    dict(key='mare_insularum_pit', name='Mare Insularum pit', lat=8.6861, lon=330.9665, inner=(159, 159),
         funnel=(255, 315), depth=63, floor_m=(63, 63), void=None, setting='Mare Insularum, mostly filled',
         sources=('lroc_pit_atlas',)),
    dict(key='mare_moscoviense_pit', name='Mare Moscoviense pit', lat=29.5025, lon=150.4647, inner=(160, 175),
         funnel=(205, 205), depth=40, floor_m=(40, 40), void=None, setting='Mare Moscoviense, far side',
         sources=('lroc_pit_atlas',)),
    dict(key='mare_serenitatis_pit', name='Mare Serenitatis pit', lat=35.1041, lon=17.4017, inner=(17, 21),
         funnel=(55, 60), depth=25, floor_m=(25, 30), void=None,
         setting='Mare Serenitatis, 500 m from a possible vent', sources=('lroc_pit_atlas',)),
    dict(key='central_mare_fecunditatis_pit', name='Central Mare Fecunditatis pit', lat=-0.9182, lon=48.6595,
         inner=(100, 130), funnel=(162, 198), depth=42, floor_m=(42, 57), void=None,
         setting='central Mare Fecunditatis; no overhang', sources=('lroc_pit_atlas', 'wagner_robinson_2014')),
    dict(key='sw_mare_fecunditatis_pit', name='Southwest Mare Fecunditatis pit', lat=-6.7521, lon=42.7595,
         inner=(15, 19), funnel=(79, 90), depth=51, floor_m=(51, 81), void=dict(reach_m=7, ceiling_m=None),
         setting='a highland kipuka 1 km from the mare edge', sources=('lroc_pit_atlas', 'wagner_robinson_2014')),
    dict(key='schluter_pit', name='Schlüter pit', lat=-5.8395, lon=276.95, inner=(23, 40), funnel=(50, 65),
         depth=57, floor_m=(57, 75), void=dict(reach_m=None, ceiling_m=None),
         setting='highland 4 km from a mare patch, Schlüter crater',
         sources=('lroc_pit_atlas', 'wagner_robinson_2014')),
]

# Roofed segments along sinuous rilles (Coombs & Hawke 1992, Table 1): the rille's position to the half degree,
# each measured segment's length and width (km) and minimum roof (m, from the largest superposed crater; None where
# no crater was measured), and the rank for a base (A prime, B good, C possible).
RILLES = [
    dict(key='rille_A', name='Rille A, Oceanus Procellarum', lat=22.0, lon=329.5,
         segments=[(2.20, 0.55, 110, 'B'), (1.10, 0.55, 110, 'A'), (0.55, 0.33, 110, 'C')]),
    dict(key='rille_B', name='Rille B, Oceanus Procellarum', lat=22.0, lon=331.0,
         segments=[(0.88, 0.44, 110, 'B'), (1.10, 0.55, 276, 'C'), (4.40, 0.03, 110, 'B'), (0.55, 0.33, 220, 'C'),
                   (1.10, 0.33, 56, 'C'), (5.50, 0.44, 220, 'B')]),
    dict(key='rille_C', name='Rille C, south of Gruithuisen K (the Oberbeck tube)', lat=35.0, lon=317.0,
         segments=[(0.90, 0.72, 30, 'A'), (1.20, 0.97, 46, 'A'), (0.66, 0.60, 136, 'B'), (0.54, 0.60, 60, 'A'),
                   (0.60, 0.45, 46, 'A'), (0.60, 0.82, 46, 'A'), (0.51, 0.82, 46, 'A')], rille_km=60),
    dict(key='rille_D', name='Rille D, Oceanus Procellarum', lat=36.0, lon=320.0, segments=[(0.90, 0.15, 30, 'C')]),
    dict(key='rille_E1', name='Rille E (E1), Oceanus Procellarum', lat=36.0, lon=320.0,
         segments=[(4.92, 1.20, 60, 'B')]),
    dict(key='rille_E2', name='Rille E (E2-E4), northern Oceanus Procellarum', lat=47.5, lon=303.0,
         segments=[(1.20, 0.36, None, 'B'), (2.40, 0.48, 120, 'B'), (1.32, 0.60, 120, 'B')]),
    dict(key='rille_F', name='Rille F, Oceanus Procellarum', lat=30.0, lon=311.5, segments=[(2.40, 1.20, 300, 'C')]),
    dict(key='rille_G', name='Rille G, Oceanus Procellarum', lat=27.0, lon=318.25,
         segments=[(2.37, 0.32, 190, 'B'), (1.87, 0.23, 375, 'C'), (0.88, 1.25, 100, 'C')]),
    dict(key='rille_H', name='Rille H, Oceanus Procellarum', lat=15.5, lon=313.0, segments=[(1.65, 0.88, 110, 'C')]),
    dict(key='rille_I', name='Rille I, Oceanus Procellarum', lat=10.5, lon=311.0,
         segments=[(1.32, 0.88, None, 'B'), (1.32, 0.77, 110, 'C'), (6.60, 0.66, None, 'C'), (4.60, 0.55, 165, 'C')]),
    dict(key='rille_J', name='Rille J (Ja and Jb), Marius Hills', lat=15.0, lon=303.0,
         segments=[(0.25, 0.33, 40, 'B'), (0.25, 0.34, 125, 'B'), (0.25, 0.25, 15, 'B'), (0.13, 0.33, 25, 'C'),
                   (0.35, 0.37, 40, 'B'), (0.25, 0.32, 25, 'B'), (1.63, 0.33, 75, 'B'), (1.63, 0.38, 40, 'C'),
                   (7.37, 0.27, 125, 'C'), (0.37, 0.28, 40, 'C'), (0.75, 0.48, 40, 'C')], rille_km=57.5),
    dict(key='rille_K', name='Rille K, Marius Hills region', lat=12.0, lon=307.0,
         segments=[(6.60, 0.55, 165, 'C'), (10.45, 1.03, 525, 'C'), (8.80, 0.88, 330, 'C'), (2.75, 0.66, 55, 'B'),
                   (1.10, 0.33, None, 'B'), (5.50, 0.55, 110, 'C')]),
    dict(key='rille_L', name='Rille L, Mare Tranquillitatis', lat=2.0, lon=44.0,
         segments=[(3.30, 1.32, None, 'C'), (6.60, 1.32, None, 'C')]),
    dict(key='rille_M', name='Rille M, Mare Tranquillitatis', lat=4.0, lon=28.0,
         segments=[(7.70, 0.88, None, 'C'), (4.40, 0.55, None, 'C')]),
    dict(key='rille_N', name='Rille N, Mare Tranquillitatis', lat=11.0, lon=20.0,
         segments=[(17.60, 0.77, 495, 'C'), (2.53, 1.10, 165, 'C'), (4.62, 0.99, None, 'C')]),
    dict(key='rille_O', name='Rille O, southern Mare Serenitatis', lat=18.0, lon=26.0,
         segments=[(0.77, 0.51, 35, 'B'), (0.32, 0.51, 25, 'B'), (0.32, 0.48, 25, 'B'), (0.88, 0.48, 15, 'B'),
                   (0.77, 0.46, 40, 'B')], rille_km=30),
    dict(key='rille_P', name='Rille P, northern Mare Imbrium by the Alpine Valley', lat=51.0, lon=352.0,
         segments=[(0.84, 0.21, 65, 'B'), (1.26, 0.21, 85, 'C'), (1.68, 0.29, 65, 'C')]),
    dict(key='rille_Q', name='Rille Q, southern Oceanus Procellarum north of Mare Humorum', lat=-16.0, lon=323.0,
         segments=[(1.65, 0.44, None, 'C')]),
    dict(key='rille_R', name='Rille R, eastern Mare Serenitatis', lat=29.0, lon=29.0,
         segments=[(1.32, 1.10, None, 'A'), (7.15, 0.55, 330, 'A')], rille_km=15),
]
RILLE_SOURCES = ('coombs_hawke_1992', 'horz_1985')
# What the sources give for the named rilles beyond their IAU positions and lengths.
NAMED_RILLE_NOTES = {
    'Rima Hadley': dict(length_km=135.0, width_m=900.0, depth_m=370.0,
                        basis='a lava channel partly roofed into a tube whose roof has since collapsed',
                        source='greeley_1971'),
    'Vallis Schroteri': dict(basis='GRAIL suggests a 1-4 km wide empty tube continuing the rille', source='blair_2015'),
    'Rima Sharp': dict(basis='GRAIL suggests a 1-4 km wide empty tube continuing the rille; a buried mass deficit '
                             'lies south of it', source='blair_2015'),
}

# Collapse chains measured on Kaguya-LOLA terrain (Sauro et al. 2020): the chain's start, its length and intact
# share, the collapses' mean width and depth (m), and the tube those give with the terrestrial collapse factor f = 1.
# The authors put the conduits a half to a third smaller if inflation dominates.
CHAINS = [
    dict(key='marius_hills_south_chain', name='Marius Hills South collapse chain', lat=13.096, lon=302.944,
         chain_km=16.9, intact_km=6.4, collapse_w_m=574, collapse_d_m=87),
    dict(key='marius_hills_north_chain', name='Marius Hills North collapse chain', lat=13.603, lon=301.953,
         chain_km=8.5, intact_km=4.1, collapse_w_m=527, collapse_d_m=96),
    dict(key='gruithuisen_chain', name='Gruithuisen collapse chain', lat=34.618, lon=316.533, chain_km=41.0,
         intact_km=27.2, collapse_w_m=858, collapse_d_m=190),
]
CHAIN_SOURCES = ('sauro_2020', 'kaku_2017')
TUBE_SCALE = (1.0 / 3.0, 1.0)       # tube size against the collapses: a third (inflation) to equal (f = 1)


# ---------------------------------------------------------------------------------------------------------------
# Water at lunar gravity.

def pressure_pa(head_m):
    """Water pressure at a depth of head_m at lunar gravity."""
    return WATER_DENSITY * MOON_SURFACE_GRAVITY * np.asarray(head_m, float)


def earth_head_m(head_m):
    """The head on Earth that gives the same pressure."""
    return np.asarray(head_m, float) * GRAVITY_RATIO


def k_moon(k_earth):
    """The same rock and water at lunar gravity: K = k rho g / mu, so K scales with g."""
    return np.asarray(k_earth, float) * GRAVITY_RATIO


def mound_m(x_m, divide_m, recharge_m_s, k_m_s, base_m):
    """Dupuit-Forchheimer water table above the shore's level at distance x from the shore, with the groundwater
    divide at divide_m, uniform recharge, and a saturated base base_m below the shore's level:
    (h + b)^2 = b^2 + (R/K) x (2L - x)."""
    x = np.minimum(np.asarray(x_m, float), divide_m)
    return np.sqrt(base_m ** 2 + recharge_m_s / k_m_s * x * (2.0 * divide_m - x)) - base_m


GROUND_TOLERANCE_M = 1.0        # a water table within this of the ground counts as reaching it


def mound_ratio_limits():
    """How much higher the Moon's mound stands than Earth's for the same rock and rain. K falls with g, so the ratio
    runs from 1 / (g ratio) where the mound is small against the saturated thickness beneath it to the square root
    of that where the mound dwarfs it."""
    return dict(mound_small_against_base=r(1.0 / GRAVITY_RATIO),
                mound_large_against_base=r(math.sqrt(1.0 / GRAVITY_RATIO)))


def groundwater_site(at: grid.Atlas, lat, lon, rain_mm_day, runoff_mm_day):
    """The steady water table about a place in three cases, solved in two dimensions on the atlas window about it
    with the seas and lakes at their levels and the ground capping the mound (groundwater.py), and the central case's
    rise. Low: the least recharge (a tenth of the rain, never more than the GCM's rain less evaporation) into
    permeable, thick basalt; central: Oahu's 36% of the rain into basalt of the Columbia Plateau's median
    conductivity at the nearside's mean basalt thickness; high: half the rain into tight, thin basalt. Recharge never
    exceeds the rock's conductivity. The rain and the rain less evaporation are the GCM fields carried onto the atlas.
    The rise starts from the water table at the sea's level under the land, at the central specific yield. Raw values
    at the place's atlas cell."""
    win = gw.Window.about(at, lat, lon)
    ci, cj = win.centre
    sea = at.sea_level_m
    rain = np.asarray(rain_mm_day)[win.ix] / 1e3 / JULIAN_DAY
    net = np.maximum(np.asarray(runoff_mm_day)[win.ix], 0.0) / 1e3 / JULIAN_DAY
    land = ~win.wet
    cases = {}
    for i, case in enumerate(('low', 'central', 'high')):
        k = float(k_moon(K_EARTH[2 - i]))
        b = BASALT_THICKNESS_M[2 - i]
        recharge = RECHARGE_SHARE[i] * rain
        if case == 'low':
            recharge = np.minimum(recharge, net)
        recharge = np.minimum(recharge, k)
        free = gw.steady(win, recharge, k, b, cap=False)
        held = gw.steady(win, recharge, k, b)
        c = dict(k=k, b=b, recharge=float(recharge[ci, cj]), uncapped=float(free['z'][ci, cj]),
                 z=float(held['z'][ci, cj]), cell_pinned=bool(held['pinned'][ci, cj]),
                 waterlogged_share=float(held['pinned'][land].mean()), base=held['base'])
        if case == 'central':
            c['earth_uncapped'] = float(gw.steady(win, recharge, float(K_EARTH[1]), b, cap=False)['z'][ci, cj])
            c['without_rain'] = float(gw.steady(win, np.zeros(win.shape), k, b, cap=False)['z'][ci, cj])
            if c['recharge'] > 0.0:
                c['series'] = gw.rise(win, recharge, k, SPECIFIC_YIELD[1], sea, held)
        cases[case] = c
    return dict(sea=sea, cases=cases)


def water_table(site, ground_m, floor_below_m):
    """A land candidate's own ground (m above the geoid, LOLA 16 px/deg) and floor against the water table of its
    atlas cell (7.6 km across) in the three cases, and how fast the water comes. Rise times scale with the specific
    yield; the central case's come from the solve, the high case's at the first rate, the recharge over the
    specific yield, which its tight rock barely drains sideways. From dry basalt the water first needs at most
    S (sea level - base) / R more to reach the sea's level. Returns the cases and each case's unrounded depth of the
    water table below the candidate's ground (negative above it)."""
    sea = site['sea']
    floor_z = ground_m - floor_below_m[0] if floor_below_m else None
    s_low, s_mid, s_high = SPECIFIC_YIELD
    out, depths = {}, {}
    for case, c in site['cases'].items():
        depth = ground_m - c['z']
        depths[case] = depth
        per_year = c['recharge'] * SECONDS_PER_YEAR
        row = dict(recharge_mm_year=r(per_year * 1e3), k_moon_m_s=r(c['k']), base_m=c['b'],
                   uncapped_mound_above_sea_m=metres(c['uncapped'] - sea), water_table_above_sea_m=metres(c['z'] - sea),
                   depth_below_ground_m=r(max(depth, 0.0)), reaches_ground=bool(depth <= GROUND_TOLERANCE_M),
                   cell_waterlogged=c['cell_pinned'], waterlogged_share_of_window_land=r(c['waterlogged_share']))
        if case == 'central':
            rain_part = c['uncapped'] - c['without_rain']
            row.update(uncapped_without_rain_above_sea_m=metres(c['without_rain'] - sea),
                       uncapped_mound_at_earth_gravity_above_sea_m=metres(c['earth_uncapped'] - sea),
                       rain_mound_moon_over_earth=r(rain_part / (c['earth_uncapped'] - c['without_rain']))
                       if rain_part > 1.0 else None)
        if case != 'low' and per_year > 0.0:
            rise = dict(specific_yield=[s_low, s_high], first_rise_m_per_year=rr([per_year / s_high, per_year / s_low]),
                        years_to_floor_at_first_rate=None if floor_z is None else rr(
                            [max(floor_z - sea, 0.0) * s / per_year for s in (s_low, s_high)]),
                        dry_basalt_adds_at_most_years=rr([(sea - c['base']) * s / per_year for s in (s_low, s_high)]))
            series = c.get('series')
            if series is not None:
                start = float(series['levels'][0])
                ninety = gw.first_reached(series, start + 0.9 * (c['z'] - start)) if c['z'] > start + 1e-6 else 0.0
                floor = gw.first_reached(series, floor_z) if floor_z is not None else None
                rise['years_to_floor'] = None if floor is None else rr([floor * s / s_mid for s in (s_low, s_high)])
                rise['years_to_90_percent'] = None if ninety is None else rr(
                    [ninety * s / s_mid for s in (s_low, s_high)])
            row['rise'] = rise
        out[case] = row
    return out, depths


# ---------------------------------------------------------------------------------------------------------------
# Keeping a tube dry.

def seepage_m3_s(head_m, area_m2, k_earth=GROUTED_K_EARTH[1], ring_m=GROUTED_RING_M):
    """Water through a grouted ring of thickness ring_m carrying a mean head over the tube's envelope area, at lunar
    gravity: flux = K_moon head / ring."""
    return float(k_moon(k_earth) * head_m / ring_m * area_m2)


def pumping_mw(q_m3_s, lift_m, efficiency=PUMP_EFFICIENCY):
    return float(WATER_DENSITY * MOON_SURFACE_GRAVITY * q_m3_s * lift_m / efficiency / 1e6)


def envelope_m2_per_m(width_m, height_m):
    """Floor, walls and vault of a semi-elliptical section over a flat floor, per metre of tube (Ramanujan's
    perimeter for the half ellipse)."""
    a, b = width_m / 2.0, height_m
    h = ((a - b) / (a + b)) ** 2
    half = math.pi * (a + b) * (1.0 + 3.0 * h / (10.0 + math.sqrt(4.0 - 3.0 * h))) / 2.0
    return width_m + half


def sealing(head_m, width_m, height_m, length_km, lift_m):
    """Seepage and pumping for a tube carrying head_m of water on its envelope (mean over the section), with the
    grouted ring's conductivity at the low, central and high ends of Earth practice."""
    area = envelope_m2_per_m(width_m, height_m) * length_km * 1e3
    rows = {}
    for label, k in zip(('tight', 'central', 'loose'), GROUTED_K_EARTH):
        q = seepage_m3_s(head_m, area, k)
        rows[label] = dict(grouted_k_earth_m_s=k, seepage_m3_s=r(q), pumping_mw=r(pumping_mw(q, lift_m)))
    return dict(head_m=r(head_m), pressure_kpa=r(float(pressure_pa(head_m)) / 1e3),
                earth_head_for_same_pressure_m=r(float(earth_head_m(head_m))), envelope_km2=r(area / 1e6),
                lift_m=r(lift_m), by_ring=rows)


# ---------------------------------------------------------------------------------------------------------------
# What a tube holds.

def level_widths(width_m, height_m, storey_m):
    """Usable widths of stacked levels in a semi-elliptical vault over a flat floor: a level counts where the vault
    leaves a full storey of headroom over it."""
    n = int(height_m // storey_m)
    z = (np.arange(n) + 1.0) * storey_m
    widths = width_m * np.sqrt(np.clip(1.0 - (z / height_m) ** 2, 0.0, 1.0))
    return widths[widths > 0.0]


def capacity(width_m, height_m, length_km):
    """Floor and people for a tube section, across storey heights, usable shares and floor per person."""
    out = {}
    for s in STOREY_M:
        widths = level_widths(width_m, height_m, s)
        floor_km2 = float(widths.sum()) * length_km * 1e3 / 1e6
        out[f'{s:g}_m_storeys'] = dict(
            levels=int(widths.size), gross_floor_km2=r(floor_km2),
            full_width_floor_km2=r(width_m * widths.size * length_km * 1e3 / 1e6),
            people_million={f'{a:g}_m2': [r(floor_km2 * 1e6 * u / a / 1e6) for u in USABLE_SHARE]
                            for a in FLOOR_M2_PER_PERSON})
    return out


def check_estimate(width_m=1000.0, length_km=50.0, floors=10, floor_m2=70.0, height_m=50.0):
    """The rough estimate the author was given, its arithmetic, and the same tube as a vault."""
    naive_km2 = width_m * length_km * 1e3 * floors / 1e6
    vault = capacity(width_m, height_m, length_km)
    # The vault's height whose levels add up to the estimate's ten full-width floors, for each storey height.
    heights = np.arange(height_m, 4.0 * height_m + 0.5, 0.5)
    matching = {f'{s:g}_m_storeys': r(float(next(h for h in heights
                                                  if level_widths(width_m, h, s).sum() >= floors * width_m)))
                for s in STOREY_M}
    return dict(width_m=width_m, length_km=length_km, floors=floors, floor_m2_per_person=floor_m2,
                floor_km2=r(naive_km2), people_million=r(naive_km2 * 1e6 / floor_m2 / 1e6),
                as_vault_of_height_m=height_m, vault=vault, vault_height_matching_ten_floors_m=matching)


# ---------------------------------------------------------------------------------------------------------------
# The first phase: shelter before the air.

def roof_shielding(roof_m):
    """A roof's areal mass against Earth's atmosphere, which shields Earth's surface."""
    lo, hi = BASALT_BULK_DENSITY[0], BASALT_BULK_DENSITY[2]
    return [dict(roof_m=t, areal_mass_t_m2=rr([t * d / 1e3 for d in (lo, hi)]),
                 earth_atmospheres=rr([t * d / EARTH_AIR_COLUMN_KG_M2 for d in (lo, hi)]))
            for t in roof_m]


def air_matching_earth_column():
    """Surface pressure at which the Open Moon's air column equals Earth's atmosphere's mass per area."""
    return dict(earth_column_t_m2=r(EARTH_AIR_COLUMN_KG_M2 / 1e3),
                lunar_surface_pressure_kpa=r(EARTH_AIR_COLUMN_KG_M2 * MOON_SURFACE_GRAVITY / 1e3))


def roof_to_hold_air(pressure_kpa):
    """Roof thickness whose weight balances an inside pressure, at lunar gravity, for the basalt's density range."""
    densest, lightest = BASALT_BULK_DENSITY[2], BASALT_BULK_DENSITY[0]
    return dict(pressure_kpa=pressure_kpa,
                roof_m=rr([pressure_kpa * 1e3 / (d * MOON_SURFACE_GRAVITY) for d in (densest, lightest)]))


def warming_lag_years(depth_m, fraction=0.5):
    """Years after a step in mean surface temperature until the rock at depth_m has moved by `fraction` of it
    (conduction, erfc solution), for the diffusivity range."""
    arg = float(erfcinv(fraction))
    return rr([(depth_m / (2.0 * arg)) ** 2 / k / SECONDS_PER_YEAR for k in (ROCK_DIFFUSIVITY[1], ROCK_DIFFUSIVITY[0])])


# ---------------------------------------------------------------------------------------------------------------
# A lasting city: daylight, air and floods.

def beam_transmission(zenith_deg, haze_daylight):
    """Direct transmission of sunlight through the clear air at a zenith angle, interpolated in air mass (1/cos)
    between the haze product's three zenith angles."""
    zs = np.array([float(k) for k in haze_daylight])
    ts = np.array([haze_daylight[k]['clear']['direct'] for k in haze_daylight])
    m = 1.0 / np.cos(np.radians(zs))
    order = np.argsort(m)
    mz = 1.0 / max(math.cos(math.radians(zenith_deg)), 1e-3)
    return float(np.exp(np.interp(mz, m[order], np.log(ts[order]), right=np.nan)))


def heliostat_share(lat_deg, above_air_lux, haze_daylight, target_lux, delivered):
    """Mirror area per m2 of lit floor for heliostats sending the beam down a vertical shaft into light guides, as a
    mean over the daylight hours at a latitude (the Sun on the equator's path; zero declination). Each heliostat
    sends the beam straight up to a fixed mirror over the shaft, which turns it down; the heliostat's mirror then
    stands at half the Sun's zenith angle to the beam, cos((90 - e) / 2). The fixed mirror's and the guides' losses
    sit inside the delivered share."""
    hours = np.linspace(-89.5, 89.5, 180)                  # hour angle, degrees
    sin_e = math.cos(math.radians(lat_deg)) * np.cos(np.radians(hours))
    e = np.degrees(np.arcsin(np.clip(sin_e, 0.0, 1.0)))
    up = e > 10.0                                         # below 10 degrees the beam is thin; mirrors stand idle
    beam = np.array([beam_transmission(90.0 - x, haze_daylight) if u else 0.0 for x, u in zip(e, up)])
    beam = np.nan_to_num(beam) * above_air_lux
    on_mirror = beam * np.cos(np.radians((90.0 - e) / 2.0))
    mean = float(on_mirror.mean())
    return dict(lat_deg=r(lat_deg), mean_beam_on_mirror_lux=r(mean), target_lux=target_lux,
                delivered_share=list(delivered),
                mirror_m2_per_floor_m2=rr([target_lux / (mean * d) for d in delivered]))


def ventilation(people, litres_per_person_s, pressure_pa_drop, fan_efficiency):
    q = people * litres_per_person_s / 1e3
    return dict(people=people, outdoor_air_m3_s=r(q),
                fan_mw=rr([q * dp / fan_efficiency / 1e6 for dp in pressure_pa_drop]))


def stack_pressure_pa(height_m, delta_t_k, air_density, air_t_k, gravity=MOON_SURFACE_GRAVITY):
    """Stack-effect pressure across a shaft: (rho_out - rho_in) g h, with the density difference from the
    temperature difference at constant pressure."""
    return air_density * delta_t_k / air_t_k * gravity * height_m


def fill_years(volume_m3, inflow_m3_s):
    return volume_m3 / inflow_m3_s / SECONDS_PER_YEAR if inflow_m3_s > 0 else None


def breach_inflow_m3_s(head_m, area_m2, discharge_coefficient=0.6):
    """Water through an opening under a head at lunar gravity (Torricelli, with an orifice coefficient)."""
    return discharge_coefficient * area_m2 * math.sqrt(2.0 * MOON_SURFACE_GRAVITY * head_m)


# ---------------------------------------------------------------------------------------------------------------
# Where each candidate lies on the Open Moon.

SEA_NAMES = {'Mare Smythii': 'the Smythii-Marginis Sea', 'Mare Marginis': 'the Smythii-Marginis Sea',
             'Mare Humboldtianum': 'the Humboldtianum Sea', 'Mare Moscoviense': 'the Moscoviense lake',
             'Mare Orientale': 'the Orientale lake'}


def water_name(body):
    """The atlas body's working name, as the geography README names the seas and lakes."""
    share = body.get('share', 0.0)
    if share > 0.1:
        return 'the Near-side Sea'
    if share > 0.05:
        return 'the South Pole-Aitken Sea'
    for name in body.get('water_names') or []:
        if name in SEA_NAMES:
            return SEA_NAMES[name]
    return 'a crater lake'


def locate(cands, at: grid.Atlas, lola16, geoid_at, lake16, level16, bodies, rain, runoff, water_km, nearest):
    """Ground above the geoid at each candidate from LOLA's 16 px/deg grid less the degree-200 geoid, against the
    atlas's sea level and the rain-fed lakes; the class; the distance to standing water and its level; the rain.
    A candidate whose ground lies below the common level is under water by the atlas's own rule (every closed
    depression below the level holds water), read at 16 px/deg: under the sea where its atlas cell or one of the
    eight around it is sea, else under a lake at the common level. It is under a lake where a rain-fed lake covers
    it, coastal where its atlas cell or one of the eight around it holds standing water, and dry beyond. Each row
    carries its unrounded ground (m above the geoid) and water depth under 'raw' for the water table and the heads."""
    out = []
    ppd = grid.LAKE_PPD
    sea = at.sea_level_m
    names = {b['label']: b for b in bodies}
    wet = at.wet()
    level_flat = at.water_level_m().ravel()
    for c, geo in zip(cands, geoid_at):
        i = int(np.clip(np.floor((90.0 - c['lat']) * ppd), 0, lola16.shape[0] - 1))
        j = int(np.floor((c['lon'] % 360.0) * ppd)) % lola16.shape[1]
        ground = float(lola16[i, j]) - geo
        cols = [(j + d) % lola16.shape[1] for d in range(-3, 4)]
        win = lola16[max(i - 3, 0):i + 4][:, cols] - geo
        ia, ja = at.cell(c['lat'], c['lon'])
        around = np.ix_([x for x in (ia - 1, ia, ia + 1) if 0 <= x < at.lat.size],
                        [(ja + d) % at.lon.size for d in (-1, 0, 1)])
        row = dict(key=c['key'], name=c['name'], lat=c['lat'], lon=c['lon'],
                   ground_above_sea_m=r(ground - sea, 4),
                   window_13km_above_sea_m=rr([float(win.min()) - sea, float(win.max()) - sea], 4),
                   atlas_cell_above_sea_m=r(float(at.height_m[ia, ja]) - sea, 4),
                   rain_mm_day=r(float(rain[ia, ja])), runoff_mm_day=r(float(runoff[ia, ja])))
        raw = dict(ground_m=ground, water_depth_m=0.0)
        if ground < sea:
            labels = at.water_label[around]
            label = int(at.water_label[ia, ja]) or (int(np.bincount(labels[labels > 0]).argmax())
                                                    if (labels > 0).any() else 0)
            raw['water_depth_m'] = sea - ground
            if label:
                row.update({'class': 'under the sea', 'water_depth_m': r(sea - ground, 4),
                            'water': water_name(names.get(label, {}))})
                if not at.water[ia, ja]:
                    row['note'] = 'below the common level behind the atlas coast, on LOLA\'s 16 px/deg ground'
            else:
                row.update({'class': 'under a lake', 'water_depth_m': r(sea - ground, 4),
                            'water': 'a hollow below the common level, too small for the atlas grid'})
        elif lake16[i, j] and ground < float(level16[i, j]):
            raw['water_depth_m'] = float(level16[i, j]) - ground
            row.update({'class': 'under a lake', 'water_depth_m': r(float(level16[i, j]) - ground, 4),
                        'water': 'a rain-fed lake'})
        else:
            level = float(level_flat[nearest[ia, ja]])
            row.update({'class': 'coastal' if wet[around].any() else 'dry', 'water_depth_m': 0.0,
                        'km_to_standing_water': r(float(water_km[ia, ja])),
                        'nearest_water_above_sea_m': r(level - sea, 4)})
        box = c.get('box')
        if box is None and c.get('type') == 'roofed rille segments':
            box = ((c['lat'] - RILLE_BOX_DEG, c['lat'] + RILLE_BOX_DEG),
                   (c['lon'] - RILLE_BOX_DEG, c['lon'] + RILLE_BOX_DEG))
        if box:
            row['box'] = dict(lat_deg=list(box[0]), lon_deg=list(box[1]), **box_shares(at, box))
        row['raw'] = raw
        out.append(row)
    return out


RILLE_BOX_DEG = 0.5     # Coombs & Hawke give each rille's position to the half degree


def box_shares(at: grid.Atlas, box):
    """Shares of a latitude-longitude box under the seas and atlas lakes, under rain-fed lakes and dry, and the dry
    land's mean height above the sea."""
    (la0, la1), (lo0, lo1) = box
    lon, w, e = at.LO % 360.0, lo0 % 360.0, lo1 % 360.0
    in_lon = ((lon >= w) & (lon <= e)) if w <= e else ((lon >= w) | (lon <= e))     # a box across 0 E wraps
    m = (at.LA >= la0) & (at.LA <= la1) & in_lon
    area = float((at.area * m).sum())
    if area <= 0.0:
        return dict(sea_share=None, lake_share=None, dry_share=None, dry_mean_above_sea_m=None)
    dry = float((at.dry * m).sum())
    return dict(sea_share=r(float((at.area * at.water * m).sum()) / area),
                lake_share=r(float((at.land * at.f_lake * m).sum()) / area), dry_share=r(dry / area),
                dry_mean_above_sea_m=r(float((at.hsl_km * at.dry * m).sum()) / dry * 1e3, 4) if dry > 0 else None)


# ---------------------------------------------------------------------------------------------------------------
# Tube sizes for the candidates that have them, and what each holds.

HEIGHT_TO_WIDTH = (0.15, 0.22)      # the collapses' depth over width (Sauro et al. 2020): Marius Hills to Gruithuisen


def tube_cases(width_m, length_km, height_m=None, roof_m=None):
    """A small and a large tube for a candidate: the small a third of the width (and height) the collapses or the
    rille give, at the lower depth-to-width ratio; the large at their full size and the higher ratio."""
    w_lo, w_hi = TUBE_SCALE[0] * width_m, TUBE_SCALE[1] * width_m
    if height_m is None:
        h_lo, h_hi = HEIGHT_TO_WIDTH[0] * w_lo, HEIGHT_TO_WIDTH[1] * w_hi
    else:
        h_lo, h_hi = TUBE_SCALE[0] * height_m, TUBE_SCALE[1] * height_m
    return dict(small=dict(width_m=w_lo, height_m=h_lo, length_km=length_km, roof_m=roof_m),
                large=dict(width_m=w_hi, height_m=h_hi, length_km=length_km, roof_m=roof_m))


def holds(case):
    """People a tube section holds at the metropolis brief's 70 m2, the levels in it and its floor; low end at
    5 m storeys and half the floor usable, high end at 4 m storeys and 70% usable."""
    w, h, L = case['width_m'], case['height_m'], case['length_km']
    tall, short = level_widths(w, h, STOREY_M[1]), level_widths(w, h, STOREY_M[0])
    lo = tall.sum() * L * 1e3 * USABLE_SHARE[0]
    hi = short.sum() * L * 1e3 * USABLE_SHARE[1]
    return dict(width_m=r(w), height_m=r(h), length_km=r(L), levels=[int(tall.size), int(short.size)],
                usable_floor_km2=rr([lo / 1e6, hi / 1e6]),
                people_million_at_70_m2=rr([lo / FLOOR_M2_PER_PERSON[1] / 1e6, hi / FLOOR_M2_PER_PERSON[1] / 1e6]))


def stable_width_m(roof_m, state='lithostatic'):
    """The widest stable span for a roof thickness in Blair et al.'s maps (interpolated in log roof thickness,
    held at the table's ends)."""
    t = np.log([x['roof_m'] for x in STABILITY])
    w = np.array([x[f'{state}_m'] for x in STABILITY], float)
    return float(np.interp(math.log(max(roof_m, 1e-3)), t, w))


def theinat_min_roof_m(width_m):
    """The thinnest stable roof for a width in Theinat et al.'s map (stepwise: the next tabulated width up)."""
    for row in STABILITY_THEINAT:
        if width_m <= row['width_m']:
            return float(row['min_roof_m'])
    return None


# ---------------------------------------------------------------------------------------------------------------
# Gravity and radar candidates. GRAIL sees only tubes with "minimum widths of several kilometers, heights of
# hundreds of meters, and lengths of tens of kilometers" (Chappaz et al. 2017), so a candidate it finds is a tube
# of that size if it is a tube at all; no candidate's own size has been published. Their capacity is given for the
# author's ten levels, the lowest in a semi-elliptical vault as high as the collapses' depth-to-width ratio makes it,
# for widths of 1-4 km (Blair et al. 2015 on the Vallis Schroteri and Rima Sharp candidates) and lengths of 10-50 km.
# box gives the search region where the sources give one.
GRAIL_SIZE = dict(width_m=(1000.0, 4000.0), length_km=(10.0, 50.0), levels=10)
GRAVITY = [
    dict(key='grail_marius_hills', name='GRAIL low along the Marius Hills skylight rille', lat=14.0917, lon=303.2299,
         box=((11.0, 16.0), (301.0, 307.0)), region='marius_hills',
         basis='a Bouguer gravity low surrounding the rille along which the Marius Hills skylight lies',
         sources=('chappaz_2017', 'sood_2016')),
    dict(key='grail_rima_sharp_mairan', name='GRAIL low south of Rima Sharp and west of Rima Mairan', lat=35.892,
         lon=313.235, box=((35.0, 40.0), (311.0, 316.0)), region='rima_sharp_mairan',
         basis='a buried mass deficit in both methods and in free-air and Bouguer maps; the point is the radar '
               'candidate P0 inside it (Kaku et al. 2018)', sources=('chappaz_2014', 'chappaz_2017', 'kaku_2018')),
    dict(key='grail_sinus_iridum', name='GRAIL low by the Sinus Iridum pit', lat=45.623, lon=331.189,
         box=None, region='sinus_iridum', basis='an underground mass deficit in the vicinity of the pit',
         sources=('sood_2016',)),
    dict(key='grail_vallis_schroteri', name='GRAIL continuation of Vallis Schroteri', lat=26.16, lon=308.42,
         box=((24.16, 26.28), (306.39, 310.70)), region='vallis_schroteri',
         basis='GRAIL suggests 1-4 km wide empty tubes continuing the sinuous rilles of Vallis Schroteri and Rima '
               'Sharp; the point is the valley\'s IAU centre and the box its IAU extent',
         sources=('blair_2015', 'chappaz_2014', 'chappaz_2017')),
]
RADAR = [
    dict(key='lrs_marius_hills', name='Radar-sounded tube along Marius Hills rille A', lat=14.0917, lon=303.2299,
         box=((13.0, 15.0), (301.85, 304.01)), region='marius_hills', length_km=(30.0, 30.0),
         floor_depth_m=(75.0, 150.0),
         basis='echo patterns along rille A and its southwest underground extension (Kaku et al. 2017); a length '
               'of 30 km (Naito et al. 2020, citing Kaku) and a depth of 75-150 m (Sauro et al. 2020, citing Kaku)',
         sources=('kaku_2017', 'naito_2020', 'sauro_2020')),
    dict(key='lrs_rima_mairan', name='Radar candidate west of Rima Mairan (P0)', lat=35.892, lon=313.235,
         box=((35.0, 37.0), (312.5, 314.5)), region='rima_sharp_mairan',
         basis='a two-peak echo suggesting "a wide cave" at depths of tens to a few hundred metres (Kaku et al. 2018)',
         sources=('kaku_2018',)),
    dict(key='minirf_tranquillitatis', name='Radar conduit below the Mare Tranquillitatis pit', lat=8.3355, lon=33.222,
         box=None, region='mare_tranquillitatis', floor_depth_m=(130.0, 175.0),
         basis='a cave conduit "tens of metres long" running west from the pit floor (Carrer et al. 2024); reported '
               'at 30-80 m long, at least 45-60 m wide and reaching 130-175 m below the surface',
         sources=('carrer_2024',)),
]
# Candidates that describe the same tube: the Marius Hills rille, its collapse chains, radar tube and gravity low;
# the Gruithuisen rille and chain; the Rima Sharp-Mairan gravity and radar candidates. The total counts each once.
REGIONS = {'rille_J': 'marius_hills', 'marius_hills_south_chain': 'marius_hills',
           'marius_hills_north_chain': 'marius_hills', 'marius_hills_pit': 'marius_hills',
           'west_marius_hills_pit': 'marius_hills', 'rille_C': 'gruithuisen', 'gruithuisen_chain': 'gruithuisen',
           'sinus_iridum_pit': 'sinus_iridum', 'mare_tranquillitatis_pit': 'mare_tranquillitatis'}
# Blair et al. 2017's finite-element stability maps (their Fig. 3, read from its reproduction in Williams 2021,
# to one 250 m column): the widest stable span for a roof thickness, for rock emplaced in one thick layer
# (lithostatic stress) and in many thin layers (Poisson stress). Half-elliptical tubes three times as wide as high.
STABILITY = [dict(roof_m=t, lithostatic_m=a, poisson_m=b) for t, a, b in (
    (1, 750, 750), (2, 1000, 1250), (5, 1500, 2250), (10, 2250, 3000), (20, 3000, 3250), (50, 3500, 3500),
    (100, 3750, 3500), (200, 4250, 2000), (500, 5000, 2000))]
# Theinat et al. 2020 (their Fig. 13 via Williams 2021): the thinnest stable roof for a width, for a fractured rock
# mass; tubes come out less stable than in Blair et al.
STABILITY_THEINAT = [dict(width_m=w, min_roof_m=t) for w, t in (
    (300, 2), (500, 5), (750, 50), (1000, 100), (1750, 100), (2000, 200), (3000, 200), (3250, 1000), (4000, 1000))]
STABILITY_SOURCES = ('blair_2017', 'blair_2015', 'williams_2021', 'theinat_2020')
QUAKES = dict(
    shallow_per_year=dict(value=5.0, basis='28 shallow (high-frequency teleseismic) moonquakes in 1971-1976 on the '
                                           'Apollo network, "averaging about five per year"', source='nakamura_1979'),
    largest_mb=dict(above=5.0, energy_j=6e10, basis='the largest in eight years of observation',
                    source='oberst_nakamura_1992'),
    site_mb_4_5_within_100_km=dict(once_in_years=400.0, basis='a lunar base at a random site',
                                   source='oberst_nakamura_1992'),
    energy_j_per_year=dict(moon=1e14, earth=1e18, source='oberst_nakamura_1992'),
    largest_observed_mb=dict(value=4.8, basis='the largest shallow moonquake observed, on a b-value of 0.55; '
                                              'extrapolated, magnitude 7-7.5 somewhere on the Moon once in 50-100 '
                                              'years', source='nakamura_1980'),
    near_scarps=dict(count=8, of=28, within_km=30.0, basis='epicentres within the distance of strong shaking from '
                                                           'young thrust-fault scarps', source='watters_2019'),
    analogue='a rate like the central United States\' intraplate seismicity, where rare large events (New Madrid, '
             'mb 7.1-7.4) cannot be ruled out (Oberst & Nakamura 1992)')
NORWEGIAN_INFLOW_L_MIN_100M = dict(urban=(2.0, 4.0), subsea=30.0, grouted_fjellinjen=2.0,
                                   unsealed_30_m_deep=(10.0, 40.0))
TEMPERATURES = dict(
    airless_deep=dict(value_c=-20.0, basis='a constant -20 C inside tubes (Hörz 1985, after Mendell); 253 K',
                      sources=('horz_1985', 'wagner_robinson_2014')),
    near_skylight=dict(value_k=290.0, basis='beyond the opening in permanent shadow, about 290 K and nearly '
                                             'constant, as a blackbody cavity at 1 AU (Horvath et al. 2022)',
                       sources=('horvath_2022',)),
    interiors_from_albedo=dict(values_c={'Marius Hills': 18.0, 'Mare Tranquillitatis': 19.0, 'Mare Ingenii': 4.0},
                               basis='thermal equilibrium from the holes\' albedo (Haruyama et al. 2010)',
                               sources=('haruyama_2010',)))
# Radiation before the air: the airless Moon's cosmic-ray dose and its fall with depth.
DOSE = dict(
    surface_gcr_msv_year=dict(solar_minimum=416.0, solar_maximum=160.5, source='naito_2020'),
    largest_sep_event_msv=dict(value=2190.0, source='naito_2020'),
    effective_dose_by_depth_msv_year=dict(values={'0': 275.0, '1': 270.0, '2': 200.0, '3': 120.0, '4': 60.0,
                                                  '5': 20.0, '6': 0.0},
                                          basis='GCR at solar minimum under 5 m of the airless surface\'s '
                                                'loose layer over 3.3 g/cm3 basalt; read from De Angelis et al. 2002, '
                                                'Fig. 3; none observable past 6 m',
                                          source='de_angelis_2002'),
    pit_floor_43_m_msv_year=dict(below=30.0, source='naito_2020'),
    tube_beyond_75_m_msv_year=dict(below=1.0, source='naito_2020'),
    earth_cosmic_msv_year=dict(value=0.39, source='naito_2020'),
    open_moon_limit_msv_day=0.027)
STORM_RAIN_MM_H = (17.0, 87.0)      # the peak beneath flashing storm cores (infrastructure review, from the CM1 box)


def floor_depth_m(c):
    """The floor's depth below the surrounding ground, [low, high] in m, where the sources give or imply it."""
    if 'floor_m' in c:
        return list(c['floor_m'])
    if c.get('type') == 'roofed rille segments':
        rille = next(x for x in RILLES if x['key'] == c['key'])
        roofs = [s[2] for s in rille['segments'] if s[2] is not None]
        widths = [s[1] for s in rille['segments']]
        if not roofs:
            return None
        h = [HEIGHT_TO_WIDTH[0] * TUBE_SCALE[0] * min(widths) * 1e3, HEIGHT_TO_WIDTH[1] * max(widths) * 1e3]
        return [min(roofs) + h[0], max(roofs) + h[1]]
    if c.get('type') == 'collapse chain':
        return [TUBE_SCALE[0] * c['collapse_d_m'], c['collapse_d_m']]
    if 'floor_depth_m' in c:
        return list(c['floor_depth_m'])
    return None


def sizes(c):
    """Tube cases (small and large) with what each holds, for the candidates that have a size."""
    out = []
    if c.get('type') == 'roofed rille segments':
        rille = next(x for x in RILLES if x['key'] == c['key'])
        for n, (length, width, roof, rank) in enumerate(rille['segments'], start=1):
            cases = tube_cases(width * 1e3, length, roof_m=roof)
            out.append(dict(segment=n, rank=rank, small=holds(cases['small']), large=holds(cases['large'])))
    elif c.get('type') == 'collapse chain':
        cases = tube_cases(c['collapse_w_m'], c['intact_km'], height_m=c['collapse_d_m'])
        out.append(dict(segment='intact length', rank=None, small=holds(cases['small']), large=holds(cases['large'])))
    elif c.get('type') == 'gravity candidate':
        out.append(dict(segment='at the detection threshold', rank=None,
                        small=holds_levels(GRAIL_SIZE['width_m'][0], GRAIL_SIZE['length_km'][0], GRAIL_SIZE['levels'],
                                           HEIGHT_TO_WIDTH[0]),
                        large=holds_levels(GRAIL_SIZE['width_m'][1], GRAIL_SIZE['length_km'][1], GRAIL_SIZE['levels'],
                                           HEIGHT_TO_WIDTH[1])))
    elif c.get('type') == 'radar void' and c.get('length_km'):
        chain = next(x for x in CHAINS if x['key'] == 'marius_hills_south_chain')
        cases = tube_cases(chain['collapse_w_m'], c['length_km'][0], height_m=chain['collapse_d_m'])
        out.append(dict(segment='30 km at the Marius Hills collapses\' size', rank=None,
                        small=holds(cases['small']), large=holds(cases['large'])))
    return out


def holds_levels(width_m, length_km, levels, height_to_width):
    """People a tube holds in its lowest `levels` levels, as in the author's estimate of ten, in the semi-elliptical
    vault the other candidates take, as high as height_to_width of its width: low end at 5 m storeys and half the
    floor usable, high end at 4 m storeys and 70%, at 70 m2 a person."""
    h = height_to_width * width_m
    tall = level_widths(width_m, h, STOREY_M[1])[:levels]
    short = level_widths(width_m, h, STOREY_M[0])[:levels]
    lo = tall.sum() * length_km * 1e3 * USABLE_SHARE[0]
    hi = short.sum() * length_km * 1e3 * USABLE_SHARE[1]
    return dict(width_m=r(width_m), height_m=r(h), length_km=r(length_km), levels=[int(tall.size), int(short.size)],
                usable_floor_km2=rr([lo / 1e6, hi / 1e6]),
                people_million_at_70_m2=rr([lo / FLOOR_M2_PER_PERSON[1] / 1e6, hi / FLOOR_M2_PER_PERSON[1] / 1e6]))


def rille_summary():
    out = []
    for x in RILLES:
        seg = x['segments']
        roofs = [s[2] for s in seg if s[2] is not None]
        out.append(dict(key=x['key'], name=x['name'], lat=x['lat'], lon=x['lon'], segments=len(seg),
                        measured_km=r(sum(s[0] for s in seg)),
                        widths_km=[min(s[1] for s in seg), max(s[1] for s in seg)],
                        min_roof_m=[min(roofs), max(roofs)] if roofs else None,
                        prime_segments=sum(1 for s in seg if s[3] == 'A'), rille_km=x.get('rille_km'),
                        sources=RILLE_SOURCES))
    return out


def chain_summary():
    out = []
    for c in CHAINS:
        cases = tube_cases(c['collapse_w_m'], c['intact_km'], height_m=c['collapse_d_m'])
        out.append(dict(c, tube_width_m=[r(cases['small']['width_m']), r(cases['large']['width_m'])],
                        tube_height_m=[r(cases['small']['height_m']), r(cases['large']['height_m'])],
                        sources=CHAIN_SOURCES))
    return out


def stability_checks():
    """Each roofed rille segment with a measured roof: its width against Blair et al.'s widest stable spans and
    its roof against Theinat et al.'s thinnest stable roof."""
    out = []
    for x in RILLES:
        for n, (length, width, roof, rank) in enumerate(x['segments'], start=1):
            if roof is None:
                continue
            w = width * 1e3
            need = theinat_min_roof_m(w)
            out.append(dict(key=x['key'], segment=n, width_m=w, roof_m=roof,
                            blair_lithostatic_m=r(stable_width_m(roof, 'lithostatic')),
                            blair_poisson_m=r(stable_width_m(roof, 'poisson')),
                            theinat_min_roof_m=need,
                            stable_in=dict(blair_lithostatic=bool(w <= stable_width_m(roof, 'lithostatic')),
                                           blair_poisson=bool(w <= stable_width_m(roof, 'poisson')),
                                           theinat=bool(need is not None and roof >= need))))
    return out
