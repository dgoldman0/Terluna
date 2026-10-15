"""Distributions of people across Earth, the Moon, the array and the colonies: step six of the population work.

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 nice -n 10 \
        python3 -m research.studies.distributions.run          # -> results/distributions.json

A coupled screen on committed products. It builds the Moon's capacity at 4, 6 and 8 billion residents with 2-4 billion
aloft and the ground's people mostly in the drier highlands, checks every band of the land with the food fields and the
settlement together, sets out the synergies among the Moon's regions, prices Earth, the array and the Solar System's
colonies, and assembles four distributions of 25-30 billion people for the author to choose among. Its inputs are the
food-systems, land-cities, ways-of-living, daily-life, inhabited-volume, conservation, population, computing
and zonal-wind products and the shared population cases, and seven files from the former domain branches, read in-tree
with their hashes checked (BOUND). Two primary sources on free-space and Mars settlements add their
published figures (LIT). Arithmetic and closed-form screens only: no crop, climate, hydrology, ecosystem, structural
or traffic model runs here.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
import re
import subprocess
from pathlib import Path

from shared.constants import AU, JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_RADIUS, MOON_SURFACE_GRAVITY, SOLAR_CONSTANT
from shared.constants import SPEED_OF_LIGHT
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.distributions/1'
OUT = HERE / 'results' / 'distributions.json'
CODE = ['run.py']
SECONDS_PER_YEAR = JULIAN_YEAR_DAYS * JULIAN_DAY


class MissingInput(RuntimeError):
    """An input this run needs is not readable here."""


# Committed products on this branch: path -> schema.
INPUTS = {
    'research/studies/food_systems/results/food_systems.json': 'terluna.research.food-systems/1',
    'research/studies/land_cities/results/land_cities.json': 'terluna.research.land-cities/1',
    'research/studies/ways_of_living/results/ways_of_living.json': 'terluna.research.ways-of-living/1',
    'research/studies/daily_life/results/daily_life.json': 'terluna.research.daily-life/1',
    'research/studies/inhabited_volume/results/inhabited_volume.json': 'terluna.research.inhabited-volume/1',
    'research/studies/conservation/results/heritage_register.json': 'terluna.research.conservation-register/1',
    'provisioning/results/population.json': 'terluna.provisioning.population/1',
    'provisioning/results/computing.json': 'terluna.provisioning.computing/1',
    'climate/results/gcm/zonal_winds_A28_dim5_moon.json': 'terluna.climate.gcm-zonal-winds/1',
    'shared/scenarios/population.json': 'terluna.scenario.population/1',
}
KEYS = dict(food='research/studies/food_systems/results/food_systems.json',
            land='research/studies/land_cities/results/land_cities.json',
            wol='research/studies/ways_of_living/results/ways_of_living.json',
            daily='research/studies/daily_life/results/daily_life.json',
            iv='research/studies/inhabited_volume/results/inhabited_volume.json',
            heritage='research/studies/conservation/results/heritage_register.json',
            pop='provisioning/results/population.json',
            comp='provisioning/results/computing.json',
            winds='climate/results/gcm/zonal_winds_A28_dim5_moon.json',
            cases='shared/scenarios/population.json')

# Files from the former domain branches, read in-tree since the joint integration, their hashes checked.
BOUND = {
    'sea_phosphorus': dict(
        path='resources/results/sea_phosphorus.json', schema='terluna.resources.sea-phosphorus/1',
        branch='main', worktree='terluna', commit=None,
        blob='649cd9a8cc4bdc943f70d79156bc4019c70ecc98',
        sha256='c5d12abbf5c18dadf016b9d664c5a71fa160beda75856933ef07db3c71141005',
        note="The land's phosphorus through the lakes to the seas, people's needs, dredging the lake deltas, sewage "
             'recovery and the natural cycle kept beside them.'),
    'industry_md': dict(
        path='resources/industry.md', schema=None,
        branch='main', worktree='terluna', commit=None,
        blob='763c9ef5d48e78ac0679944887e6867c28311de7',
        sha256='394f707ab1478370c12af288e742fd440057aa6e226cd30da44ca788b3e85200',
        note="The array's habitats as the largest continuing material stream, the places' sunlight and the Moon's "
             "people's stocks; the values used are transcribed in TRANSCRIBED and checked against these bytes."),
    'phosphorus': dict(
        path='resources/results/phosphorus.json', schema='terluna.resources.phosphorus/1',
        branch='main', worktree='terluna', commit=None,
        blob='bab3177a4c9119a3f57fb4a6fc01e62307a2f444',
        sha256='4089899e2b70e8b7888bb83e705e1773f316a1b3c220a43313a48b2578b615a5',
        note='The energy to lift material to the sky towns at 10 km.'),
    'ecology_people': dict(
        path='biosphere/ecology/results/people.json', schema='terluna.biosphere.ecology-people/1',
        branch='main', worktree='terluna', commit=None,
        blob='e76556670bf282fa37635f5a4f3976a1f9adf0f1',
        sha256='758a8dcf68c6eb58082e5749e22c35e2da7ffd530f7483622e6c1b3011779ee8',
        note="The kept shares and the lunar design's land, the basis of the food and lighting studies."),
    'living_light': dict(
        path='research/studies/living_light/results/living_light.json', schema='terluna.research.living-light/1',
        branch='main', worktree='terluna', commit=None,
        blob='7479768679bda3bb5ef2178f141f48271571485a',
        sha256='24941b1f4e33eb849aeaaaffdd80df7428b37727daf661b4de5a97195b454521',
        note='Land under a light-polluted sky a person by home and lighting design, the metropolis and residents '
             'aloft; lamp hours by home; the dark belt.'),
    'sky_ecology': dict(
        path='research/studies/sky_ecology/results/sky_ecology.json', schema='terluna.research.sky-ecology/1',
        branch='main', worktree='terluna', commit=None,
        blob='0723a9375a78fba45bf86fbf7545ae691c5fd915',
        sha256='2f8156f73806aab57002a07a6a6ec79230d728f02c891ccbab15c764ca6bebf3',
        note="Aerophyte cover and the encounter rule for bodies on the same wind tracks."),
    'places': dict(
        path='research/studies/joint_synthesis/results/places.json',
        schema='terluna.research.joint-synthesis-places/1',
        branch='main', worktree='terluna', commit=None,
        blob='3aa6d0939faf502433c03105b831a9994ab3315d',
        sha256='9707e4d1c7332ba352b90dbf77c66fde080ec5e0fe24e57b54067c61ecdcfc8d',
        note="Main's corrected per-place night light and lightning at the summit, the highland box and the polar "
             'plain.'),
}

# Values read from the bound industry note; each phrase must occur in its bytes (test_run).
TRANSCRIBED = {
    'array_habitat_gt_per_billion': dict(value=1005.0, phrase='Each billion residents is 1,005 Gt',
                                         unit='Gt of habitat a billion residents, the Stanford torus'),
    'array_structure_t_per_person': dict(value=15.0, phrase='with 15 t of structure', unit='t a person'),
    'array_growth_kg_s_per_billion': dict(value=[6400.0, 32000.0], phrase='6,400–32,000 kg/s',
                                          unit='kg/s a billion residents grown over 1,000-5,000 years'),
    'array_growth_years': dict(value=[1000.0, 5000.0], phrase='grown over 1,000–5,000 years', unit='years'),
    'water_shield_share_of_build': dict(value=[0.003, 0.1], phrase='0.003–0.1% of what the build delivers',
                                        unit='percent of the build’s delivered mass, the whole habitat as water'),
    'moon_stocks_t_per_person': dict(value=[43.0, 335.0], phrase='each billion people hold 43–335 Gt',
                                     unit='t a person of stocks in use, decent living to industrial countries'),
    'main_belt': dict(value=[2.7, 187.0], phrase='| Main belt, 2.7 AU |', unit='AU; W/m2 (the row\'s 187 W/m²)'),
    'main_belt_sunlight': dict(value=187.0, phrase='| 187 W/m² |', unit='W/m2 at 2.7 AU'),
    'trojans': dict(value=[5.2, 50.0], phrase="| Jupiter's trojans, 5.2 AU |", unit='AU; W/m2 (the row\'s 50 W/m²)'),
    'trojans_sunlight': dict(value=50.0, phrase='| 50 W/m² |', unit='W/m2 at 5.2 AU'),
    'kuiper_belt': dict(value=[40.0, 0.85], phrase='| Kuiper belt, 40 AU |', unit='AU; W/m2 (the row\'s 0.85 W/m²)'),
    'kuiper_belt_sunlight': dict(value=0.85, phrase='| 0.85 W/m² |', unit='W/m2 at 40 AU'),
}
PLACES_AU = ('main_belt', 'trojans', 'kuiper_belt')

# Published figures from the two primary sources read for this study (sources.json).
LIT = {
    'kalpana_one': dict(
        value=dict(residents=3000.0, radius_m=250.0, length_m=325.0, shield_t_m2=10.0, total_mass_t=7.0e6,
                   living_m2_per_resident=170.0, kw_per_resident=60.0, agriculture_kw_per_resident=50.0,
                   other_kw_per_resident=10.0, farm_m2_per_person=50.0),
        unit='a revised orbital settlement design: people, m, t/m2, t, m2 and kW a resident as its authors count it',
        source=['globus2007_kalpana']),
    'city_npc': dict(
        value=dict(colony_kwe_per_member=2.0, city_kwe_per_capita=5.46, city_gwe=5.0, people=1.0e6,
                   greenhouse_ha=15000.0, diet_kcal_day=1600.0, tunnel_people=25000.0, tunnel_km=25.0,
                   tunnel_diameter_m=100.0),
        unit='a million-person Mars city plan: kWe a person, GWe, ha of greenhouse, kcal a day, a lava-tube district',
        source=['donofrio2021_city_npc']),
}

# Design guesses, each a range: every Moon case runs over all their combinations.
GUESSES = {
    'highland_share': dict(value=[0.6, 0.75],
                           unit="share of the ground's people in the highlands, the summit metropolis included",
                           note="the author's 'mostly', read as three-fifths to three-quarters"),
    'tower_share': dict(value=[1.0 / 3.0, 0.5],
                        unit="share of the dry highlands' people in tower districts of 30 storeys",
                        note='the rest in highland towns at 10,000 a km2'),
    'relative_drift_m_s': dict(value=[0.1, 1.0], unit='m/s between neighbours drifting at one height',
                               note="the sky-ecology study's design guess for aerophytes on the same tracks, taken "
                                    'here for sky towns'),
    'homestead_people_per_km2': dict(value=[2.0, 20.0],
                                     unit="homestead residents a km2 of field among the middle latitudes' farms (C')",
                                     note="a household of four on each 20-200 ha of field, the span of Earth's "
                                          'family farms from mixed-farming country to the mechanised grain belts'),
}
HOMESTEAD_HOUSEHOLD = 4.0             # people in a homestead household, with the guess above

MOON_BILLION = (4, 6, 8)
ALOFT_BILLION = (2, 3, 4)


class LandRule(float):
    """A land rule, as the food study defines it. Its value is the share of each zone, and of the sky, closed to
    farms, grazing, the homesteads' homes, cities, towns and sky towns; `hybrid` is the part of each zone held as parks
    and gardens, the sky's hybrid part holding no towns, and `key` names the rule in the products. Two rules are equal
    only when both parts are, so a kept half and the compromise stay apart, and a rule with a hybrid part needs its own
    key."""

    def __new__(cls, closed, hybrid=0.0, key=None):
        if hybrid and not key:
            raise ValueError('a land rule with a hybrid part needs a key')
        rule = super().__new__(cls, closed)
        rule.hybrid = hybrid
        rule.key = key or f'kept_{round(closed * 100)}'
        return rule

    def __eq__(self, other):
        if isinstance(other, LandRule):
            return float(self) == float(other) and self.hybrid == other.hybrid
        return float(self) == other and not self.hybrid

    def __ne__(self, other):
        return not self == other

    __hash__ = float.__hash__


# The author's rule of 11 October 2026 (research/decisions.md, People on the Moon and beyond): at least half of the
# land and of the dark sky wild is preferred, applied here to every zone and to the sky (kept_50); it can come down to
# 30% wild with another 20% hybrid, a compromise between keeping 30% and keeping half. The hybrid land is heavily
# manicured megaparks and great gardens for leisure, with the small shops, restrooms and places to rest they need, so
# farms, grazing, the homesteads' homes and towns keep to the rest. At the author's request the earlier reading, half
# wild with another 20% hybrid, which leaves farms and towns 30%, runs for comparison. The food study's
# basis.land_rules carry the same rules, and checks.land_rules compares them.
KEPT = (LandRule(0.3), LandRule(0.5), LandRule(0.5, 0.2, 'wild_30_hybrid_20'),
        LandRule(0.7, 0.2, 'wild_50_hybrid_20'))
CHOSEN_MIX = 'power_heavy'
SUPPLIES = ('fusion', 'array')
DESIGNS = ('shielded', 'twilight', 'living')
DISTRICT_T = ('20.0', '50.0', '100.0')
# Placement combinations: (highland share, tower share, scattered land a person), each index 0 = low, 1 = high.
COMBOS = tuple((h, t, s) for h in (0, 1) for t in (0, 1) for s in (0, 1))
TWILIT_HOMES = ('nearside_towns', 'limb_towns', 'far_side_places', 'high_latitude_towns')
SCATTERED_HOMES = ('coastal_and_lake_towns', 'wetland_edge_towns', 'canopy_towns', 'farming_country',
                   'nearside_towns', 'limb_towns', 'far_side_places', 'high_latitude_towns')
HOME_REGION = dict(coastal_and_lake_towns='big_water_shore', wetland_edge_towns='wet_margin',
                   canopy_towns='rain_belt', farming_country='rain_belt', nearside_towns='middle_nearside',
                   limb_towns='middle_limb', far_side_places='middle_far_side',
                   high_latitude_towns='high_latitude')
FOOD_BANDS = ('0-15', '15-30', '30-45', '45-60', '60-75', '75-90')
DRY_BANDS = ('30-45', '45-60', '60-90')
WATER_BANDS = {'30-45': (30.0, 45.0), '45-60': (45.0, 60.0), '60-90': (60.0, 90.0)}
LAND_BANDS = ('all_land_0_30', 'highlands_30_45', 'highlands_45_90', 'middle_wet_zone', 'fog_desert_and_polar')
WET_RULES = ('area', 'dry_side')       # the middle band's scattered land: spread by area, or off full middle-wet fields
TIE_KM = 0.5                          # tube candidates this close in distance to a site count as tied
JOIN_KM = 250.0                       # an undersea community joins a tube within this distance
EARTH_SENSITIVITY = (8.0, 11.0)       # beside the distributions' 9 billion

DISTRIBUTIONS = {
    'A': dict(name='The Moon at 4 billion, the array and the colonies carrying more',
              earth=9.0, moon=4.0, aloft=2.0, array=7.0, colonies=5.0),
    'B': dict(name='The Moon at 6 billion', earth=9.0, moon=6.0, aloft=3.0, array=6.0, colonies=4.0),
    'C': dict(name='The Moon at 8 billion', earth=9.0, moon=8.0, aloft=4.0, array=5.0, colonies=3.0),
    'D': dict(name='30 billion, the Moon at 8 billion mostly on the land',
              earth=9.0, moon=8.0, aloft=2.0, array=8.0, colonies=5.0),
    "C'": dict(name="C', the author's distribution: C with the middle latitudes' lowland people moved",
               earth=9.0, moon=8.0, aloft=4.0, array=5.0, colonies=3.0, case='cprime'),
}

EVIDENCE = ('A screen on committed and bound products: the food study\'s mixes and their feasible ranges, the '
            'land-cities highlands, tubes and water lift, the ways-of-living settings and sky towns, the '
            'inhabited-volume envelopes and roaming case, the zonal-wind product\'s spreads and gathering, the '
            'heritage register, the population and computing screens, ecology\'s living light and sky ecology, '
            'resources\' phosphorus routes and main\'s places, with two primary sources on free-space and Mars '
            'settlements. People are placed by rule; every cost is arithmetic on those products. No settlement, '
            'farm, tube, sky town, traffic system, habitat or colony has been designed or tested, and no capacity is '
            'credited to the colonies.')
READING_RULE = ('People in billions unless named; areas km2; water km3 a year; power, energy and heat TW (TW a '
                'billion is kW a person); phosphorus kt a year unless a key names Mt; mass Gt; shares are '
                'fractions and _percent values percentage points. moon.cases is keyed by lunar residents and those '
                'aloft ("6_3") and then by land rule. kept_30 and kept_50 keep 30% and half of every zone and of the '
                'sky wild; wild_30_hybrid_20 is the author\'s compromise, 30% wild and another 20% hybrid, and '
                'wild_50_hybrid_20 half wild and another 20% hybrid, kept for comparison at the author\'s request '
                '(assumptions.land_rules: closed is the share closed to farms, grazing, the homesteads\' homes, '
                'cities, towns and sky towns, hybrid the part of each zone held as the author\'s megaparks and great '
                'gardens, the sky\'s hybrid part holding no towns; checks.land_rules compares them with the food '
                'study\'s). Capacities, traffic and the wild shares use the land and sky outside the closed share, '
                'which differs by rule, so used_share_of_dry_land compares the rules directly; '
                'sky.towns_grid_share_of_their_sky gives the sky towns\' grid against that sky, and land.parks_km2 and '
                'land.parks_share_of_dry_land give the parks. Each case runs '
                'over 16 corners, every combination of the placement guesses with every input at its favourable and '
                'unfavourable end; full gives [low, high] over them all, feasible over the corners whose land bands '
                'all fit (null when all do), and remedy the same case with the dry highland towns poleward of 45 '
                'degrees first, the scattered communities on the twilit and Earthlit homes and the middle band\'s '
                'scattered communities on the fog-desert side of their homes, where some corner overfills '
                '(buildable_corners_by_parts counts the corners each part fits). bands gives each corner range\'s '
                'use, capacity and overfill by land band; the middle-wet zone is a band of its own, which the '
                'scattered communities\' middle-band land (30-48.6 degrees) covers by area, and its fields are held '
                'to the food fill\'s room against the food product\'s rounding. moon.cprime gives the results of '
                'C\' (C prime) by rule, and each distribution\'s moon.land_rules its own. Energy: '
                'primary_tw counts the residents\' primary energy and the electricity at fusion\'s input, 2.0-3.0 W '
                'a watt (one basis); electricity_tw is the electricity alone; heat_by_component_tw splits the heat. '
                'fusion and array name the supply whose heat the Moon takes. The chosen food mix is the food study\'s '
                'power-heavy mix; food_menu gives every mix. distributions holds the four for the author to choose '
                'among; earth_sensitivity moves Earth to 8 and 11 billion beside them.')
MARKS = dict(
    screen='Numbers computed here: the placement of people and its land, the land bands, energy, heat, dimming, '
           'water, phosphorus, lit sky, wild land and sky, sky traffic, the tube pairs, the synergies, and every '
           'distribution.',
    derived='Arithmetic on cited values: the published settlement figures carried per person and per billion (LIT), '
            'and the values transcribed from the bound industry note.')


# ---------------------------------------------------------------------------------------------------------------
# Reading inputs.

def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(path) -> str:
    return sha256(path)[:16]


def load_json(rel, schema):
    path = ROOT / rel
    if not path.exists():
        raise MissingInput(f'{rel} is absent')
    data = json.loads(path.read_text())
    if data.get('schema') != schema:
        raise ValueError(f'Unexpected schema in {rel}: {data.get("schema")}')
    return data


def bound_bytes(key) -> bytes:
    """A bound input's bytes: the in-tree file, which since the joint integration of 12 October 2026 holds every
    former branch's inputs; its hash is checked against the pin."""
    b = BOUND[key]
    local = ROOT / b['path']
    if not local.exists():
        raise MissingInput(f"{b['path']} is not in this tree")
    raw = local.read_bytes()
    if hashlib.sha256(raw).hexdigest() != b['sha256']:
        raise ValueError(f"{b['path']} does not match its pinned hash")
    return raw


def load_bound(key):
    data = json.loads(bound_bytes(key))
    if data.get('schema') != BOUND[key]['schema']:
        raise ValueError(f'Unexpected schema in bound input {key}')
    return data


def transcribed_check(text):
    """Every transcribed phrase occurs in the bound industry note."""
    return {k: v['phrase'] in text for k, v in TRANSCRIBED.items()}


# ---------------------------------------------------------------------------------------------------------------
# Small helpers.

def r(x, n=3):
    """Round a float to n significant figures."""
    if x is None or isinstance(x, (bool, int, str)):
        return x
    x = float(x)
    if x == 0.0 or not math.isfinite(x):
        return x
    return float(f'{x:.{n}g}')


def rd(o, n=3):
    if isinstance(o, dict):
        return {k: rd(v, n) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [rd(v, n) for v in o]
    return r(o, n)


def merge(rows):
    """[min, max] of every numeric leaf over a list of like-shaped dicts; other leaves taken from the first."""
    first = rows[0]
    if isinstance(first, dict):
        return {k: merge([row[k] for row in rows]) for k in first}
    if isinstance(first, bool) or first is None or isinstance(first, str):
        return first
    if isinstance(first, (int, float)):
        values = [float(v) for v in rows]
        return [min(values), max(values)]
    if isinstance(first, list):
        return first
    raise TypeError(type(first))


def gc_km(lat1, lon1, lat2, lon2):
    """Great-circle distance on the Moon's sphere, km (the haversine form)."""
    p1, l1, p2, l2 = map(math.radians, (lat1, lon1, lat2, lon2))
    h = math.sin((p2 - p1) / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin((l2 - l1) / 2) ** 2
    return MOON_RADIUS / 1e3 * 2.0 * math.asin(min(1.0, math.sqrt(h)))


def dumps(product):
    """JSON with one key per line and lists of plain values kept on one line."""
    text = json.dumps(product, indent=1, ensure_ascii=False)
    return re.sub(r'\[([^\[\]{}]*)\]', lambda m: '[' + re.sub(r'\s+', ' ', m.group(1)).strip() + ']', text)


# ---------------------------------------------------------------------------------------------------------------
# The basis: every product value this study uses, gathered once with where it comes from.

def wind_spread(winds, height_km):
    """Area- and coverage-weighted spread of the 3-day-mean wind at one height, from the 10th and 90th percentiles
    (2.563 standard deviations apart), and the relative speed of two towns whose winds are uncorrelated."""
    g, h, lats = winds['ground_frame'], winds['axes']['height_km'], winds['axes']['lat_deg']
    i = h.index(height_km)
    su = sv = w_sum = 0.0
    for j, lat in enumerate(lats):
        vals = [g[k][i][j] for k in ('u_p10_m_s', 'u_p90_m_s', 'v_p10_m_s', 'v_p90_m_s')]
        if None in vals:
            continue
        w = math.cos(math.radians(lat)) * g['coverage'][i][j]
        su += w * ((vals[1] - vals[0]) / 2.563) ** 2
        sv += w * ((vals[3] - vals[2]) / 2.563) ** 2
        w_sum += w
    sigma_u, sigma_v = math.sqrt(su / w_sum), math.sqrt(sv / w_sum)
    return dict(height_km=height_km, sigma_u_m_s=sigma_u, sigma_v_m_s=sigma_v,
                uncorrelated_relative_m_s=math.sqrt(2.0 * (sigma_u ** 2 + sigma_v ** 2)))


# The homes the middle latitudes' people can move to under C', and the study's measures of what sustains people
# there. Dry bands take the land-cities rows for their latitudes; lowland homes take their own rows.
DESTINATIONS = dict(
    high_latitude_lowland=dict(name='high-latitude lowland homes poleward of 48.6 degrees', kind='lowland',
                               band='fog_desert_and_polar', homes=('high_latitude_towns',), lat=(60, 90),
                               rain_bands=('high',), places=('polar plain',)),
    lowland_0_30=dict(name='lowland homes equatorward of 30 degrees', kind='lowland', band='all_land_0_30',
                      homes=('coastal_and_lake_towns', 'wetland_edge_towns', 'canopy_towns', 'farming_country'),
                      lat=(0, 30), rain_bands=('equatorial', 'tropical'),
                      places=('southern Mare Nubium', 'eastern Smythii headland')),
    dry_30_45=dict(name='dry highland cities at 30-45 degrees', kind='dry', band='30-45', lat=(30, 60),
                   light_home='far_side_places', rain_bands=('middle',), places=('highland box',)),
    dry_45_60=dict(name='dry highland cities at 45-60 degrees', kind='dry', band='45-60', lat=(30, 60),
                   light_home='high_latitude_towns', rain_bands=('middle', 'high'),
                   places=('highland box', 'polar plain')),
    dry_60_90=dict(name='dry highland cities at 60-90 degrees', kind='dry', band='60-90', lat=(60, 90),
                   light_home='high_latitude_towns', rain_bands=('high',), places=('polar plain',)),
)
# Where C' moves people from, measured the same way for comparison.
ORIGIN = dict(name='lowland homes at 30-48.6 degrees', kind='lowland', band='middle', lat=(30, 60),
              homes=('nearside_towns', 'limb_towns', 'far_side_places', 'coastal_and_lake_towns', 'wetland_edge_towns'),
              rain_bands=('middle',), places=('Oceanus Procellarum by Russell', 'South Pole-Aitken coast by Mare Ingenii'))


def destination_measures(hl, wol, ll, places, water_rows, high_homes_km2):
    """Each destination's water energy with nothing recovered, comfortable hours under CM1 and the GCM, rain,
    lightning, walking-lamp hours and lit sky a person, read from the products for its latitudes and height."""
    rings = hl['comfort']['rings_corrected_for_height']

    def ring(lat, height_km):
        return next(r for r in rings if r['abs_lat_deg'][0] == lat[0] and r['abs_lat_deg'][1] == lat[1]
                    and r['height_m'][0] <= height_km * 1e3 < r['height_m'][1])

    lowland_gw = wol['settings']['high_latitude_towns']['per_billion']['domestic_gross_lift_from_seas_gw']
    bands = wol['moon']['bands']
    pp = ll['headlines']['per_person_m2']['shielded']['homes']
    hours = {h: v['needs']['walking']['hours_below_mean'] for h, v in ll['needs']['homes'].items()}
    out = {}
    for key, spec in list(DESTINATIONS.items()) + [('origin_middle_lowland', ORIGIN)]:
        if spec['kind'] == 'dry':
            band_lat = WATER_BANDS[spec['band']]
            water = water_rows[band_lat]['gw_per_billion_from_nearest_p50']
            rows = [ring(spec['lat'], hl['land']['poleward_of_30']['mean_height_above_sea_km'])]
            light_homes = (spec['light_home'],)
        else:
            water = lowland_gw
            rows = [ring(spec['lat'], wol['settings'][h]['physical']['height_above_sea_km_mean'])
                    for h in spec['homes']]
            light_homes = spec['homes']
        rain = [bands[b]['rain_mm_day_all_surfaces'] for b in spec['rain_bands']]
        flashes = [places['places'][x]['lightning']['central'] for x in spec['places']]
        out[key] = dict(
            name=spec['name'], kind=spec['kind'], band=spec['band'],
            water_gw_per_billion_nothing_recovered=list(water),
            comfort_cm1_strict_hours=[min(r['cm1_strict_hours'] for r in rows),
                                      max(r['cm1_strict_hours'] for r in rows)],
            comfort_gcm_strict_hours=[min(r['gcm_strict_hours'] for r in rows),
                                      max(r['gcm_strict_hours'] for r in rows)],
            rain_mm_day=[min(rain), max(rain)], flashes_per_km2_year=[min(flashes), max(flashes)],
            walking_lamp_hours=[min(hours[h] for h in light_homes), max(hours[h] for h in light_homes)],
            lit_m2_per_person=[min(pp[h][0] for h in light_homes), max(pp[h][1] for h in light_homes)],
            homes=list(light_homes) if spec['kind'] == 'lowland' else None,
            land_km2=high_homes_km2 if key == 'high_latitude_lowland' else None)
    return out


def build_basis(d, b):
    food, land, wol, pop, comp, winds, daily, iv, cases = (d['food'], d['land'], d['wol'], d['pop'], d['comp'],
                                                           d['winds'], d['daily'], d['iv'], d['cases'])
    ll, se, sp, ph, eco, places = (b['living_light'], b['sky_ecology'], b['sea_phosphorus'], b['phosphorus'],
                                   b['ecology_people'], b['places'])
    hl = land['highlands']
    metro = wol['settings']['summit_metropolis']
    roam = wol['settings']['roaming_sky_towns']['per_billion']
    aer = wol['aerial_reference']
    pb = wol['per_billion_common']
    water_rows = {tuple(x['abs_lat_deg']): x for x in hl['water']['by_latitude']}
    fb = {f"{x['abs_lat_deg'][0]:g}-{x['abs_lat_deg'][1]:g}": x for x in food['basis']['zones']['bands']}
    high_bands = {k: fb[k]['highland_dry_land_km2'] for k in FOOD_BANDS}
    dry_bands = {'30-45': high_bands['30-45'], '45-60': high_bands['45-60'],
                 '60-90': high_bands['60-75'] + high_bands['75-90']}
    zones = food['basis']['zones']['zones']
    lift_row = next(t for t in ph['transport'] if t['altitude_m'] == 10000 and t['p_mass_fraction'] == 0.1)
    fl = food['checks']['food_lift_to_10_km']
    designs = ll['skyglow']['designs']['lunar_design']
    near_side = sp['partition']['by_receiving_water']['Oceanus Procellarum']
    near_side_mg = sp['seas']['bodies']['near_side_sea']['table_mg_m3']['dissolved_after_lakes']['1000_yr']
    bands = wol['moon']['bands']
    sky_clear = {k: v['sky_clear_of_terrain_km2'] for k, v in bands.items()}
    sky_total = sum(sky_clear.values())
    occupancy = wol['winds']['occupancy_by_band']
    gather = {k: occupancy[k] / (sky_clear[k] / sky_total) for k in occupancy}
    fixed10 = winds['trajectories']['fixed_points']['10']
    gathering10 = se['tenants']['contact']['gathering']['10.0']
    systems = food['systems']
    fog_glass = systems['fog_desert_glass_vented']
    irr = systems['irrigated_fields']
    heritage_sites = wol['settings']['undersea_communities']['design']
    regions = wol['moon']['regions']
    home_0_30 = {h: (regions[HOME_REGION[h]]['bands_dry_land_km2']['equatorial']
                     + regions[HOME_REGION[h]]['bands_dry_land_km2']['tropical'])
                 / regions[HOME_REGION[h]]['dry_land_km2'] for h in SCATTERED_HOMES}
    # the ways-of-living middle band, 30-48.6 degrees, holds the food study's middle-wet zone within the scattered
    # homes; with no overlay of homes on zones, the homes' middle-band land is taken to cover the zone evenly
    home_middle = {h: regions[HOME_REGION[h]]['bands_dry_land_km2']['middle'] / regions[HOME_REGION[h]]['dry_land_km2']
                   for h in SCATTERED_HOMES}
    home_high = {h: regions[HOME_REGION[h]]['bands_dry_land_km2']['high'] / regions[HOME_REGION[h]]['dry_land_km2']
                 for h in SCATTERED_HOMES}
    middle_homes_km2 = sum(hl['scattered']['homes_km2'][h] * home_middle[h] for h in SCATTERED_HOMES)
    high_homes_km2 = sum(hl['scattered']['homes_km2'][h] * home_high[h] for h in SCATTERED_HOMES)
    measured = destination_measures(hl, wol, ll, places, water_rows, high_homes_km2)
    origin = measured.pop('origin_middle_lowland')
    destinations = measured
    access = hl['access']
    zone_lat = zones['middle_wet']['abs_lat_deg_mean']
    roaming = iv['roaming_control_sensitivity']
    sky_town = next(x for x in daily['volume']['heights'] if x['label'] == 'a roaming sky town')
    case_earth = [c['earth'] for c in cases['cases'].values()]
    case_elsewhere = [c['elsewhere'] for c in cases['cases'].values()]
    return dict(
        moon=dict(area_km2=wol['moon']['area_km2'], dry_land_km2=wol['moon']['dry_land_km2'],
                  sky_clear_of_terrain_km2=sky_clear, sky_clear_total_km2=sky_total),
        kept_shares=eco['capacity']['kept_shares'],
        land_bands=dict(dry_land_km2={k: fb[k]['dry_land_km2'] for k in FOOD_BANDS}, highland_km2=high_bands,
                        zones_km2={k: zones[k]['dry_land_km2'] for k in zones},
                        middle_band_homes_km2=middle_homes_km2,
                        middle_wet_share_of_middle_homes=zones['middle_wet']['dry_land_km2'] / middle_homes_km2),
        metropolis=dict(billion=metro['limit']['holds_billion'],
                        dry_land_km2=metro['home_dry_land_km2'],
                        people_per_km2=metro['limit']['numbers']['density_per_km2'],
                        heat_w_m2=metro['per_billion']['heat_over_town_w_m2'],
                        home_runoff_km3_year=metro['per_billion']['home_runoff_km3_year'],
                        rain_mm_day=hl['summit_rise'][0]['rain_mm_day'],
                        lake_share=hl['summit_rise'][0]['lake_share'],
                        within_100_km=hl['summit_rise'][2]),
        towers=dict(people_per_km2=hl['tower_district']['people_per_km2'],
                    building_t_per_person=hl['tower_district']['building_t_per_person'],
                    heat_w_m2=wol['settings']['elevated_districts']['per_billion']['heat_over_town_w_m2']),
        towns=dict(people_per_km2=land['assumptions']['people']['densities_per_km2']['towns']),
        scattered=dict(m2_per_person=[wol['assumptions']['settlement_m2_per_person']['town'][0],
                                      wol['assumptions']['settlement_m2_per_person']['spread'][1]],
                       homes_km2={h: hl['scattered']['homes_km2'][h] for h in SCATTERED_HOMES},
                       home_share_0_30=home_0_30, home_share_middle=home_middle, home_share_high=home_high,
                       high_homes_km2=high_homes_km2,
                       mix_weight=wol['assumptions']['mix_weight']),
        highlands=dict(all_km2=hl['land']['all_highland']['dry_land_km2'],
                       dry_km2=sum(dry_bands.values()),
                       dry_bands_km2=dry_bands,
                       dry_runoff_km3_year=hl['land']['poleward_of_30']['runoff_km3_year'],
                       dry_rain_mm_day=hl['land']['poleward_of_30']['rain_mm_day'],
                       dry_far_side_share=hl['land']['poleward_of_30']['far_side_share'],
                       dry_flashes=hl['land']['poleward_of_30']['flashes_per_km2_year'],
                       dry_mean_height_km=hl['land']['poleward_of_30']['mean_height_above_sea_km'],
                       wet_km2=hl['land']['equatorward_of_30']['dry_land_km2'],
                       wet_rain_mm_day=hl['land']['equatorward_of_30']['rain_mm_day']),
        undersea=dict(billion=heritage_sites['residents_all_candidates'] / 1e9,
                      per_dome=heritage_sites['residents_per_dome'], sites=heritage_sites['candidate_sites'],
                      depth_m=heritage_sites['depth_m']),
        energy=dict(band_kw=[pb['energy_tw']['efficient'], pb['energy_tw']['high_income']],
                    european_kw=pb['energy_tw']['european'],
                    heat_tw_per_billion_fusion=[pb['heat_on_moon_tw']['efficient'],
                                                pb['heat_on_moon_tw']['high_income']],
                    placement_power_from_orbit=pop['heat']['placement_factor']['industry_in_orbit_power_from_orbit'],
                    heat_per_w=food['basis']['heat']['heat_per_w'],
                    w_m2_per_tw=food['basis']['heat']['moon_w_m2_per_tw'],
                    dimming_percent_per_tw=food['basis']['heat']['dimming_percent_per_tw'],
                    close_computing_tw_per_billion=comp['moon_side']['cases']['per_billion_people_x100']['tw'],
                    fleet_glow_dimming_percent=[x / pop['heat']['design_sunlight_w_m2'] * 100.0
                                                for x in pop['heat']['fleet_glow_w_m2']]),
        water=dict(domestic_km3_per_billion=hl['water']['domestic_km3_per_billion'],
                   metropolis_gw_per_billion=water_rows[(0.0, 15.0)]['gw_per_billion_from_nearest_p50'],
                   dry_gw_per_billion_by_band={b: water_rows[WATER_BANDS[b]]['gw_per_billion_from_nearest_p50']
                                               for b in DRY_BANDS},
                   dry_lift_km_by_band={b: water_rows[WATER_BANDS[b]]['lift_km_p50_p90'][0] for b in DRY_BANDS},
                   scattered_makeup_gw_per_billion=pb['domestic_makeup_lift_from_seas_gw'],
                   recovery_cut=hl['water']['per_billion_by_recovery']['0.9']['from_nearest_gw'][0]
                   / hl['water']['per_billion_by_recovery']['0']['from_nearest_gw'][0],
                   makeup_l_per_resident_day=[x * 1e12 / 1e9 / JULIAN_YEAR_DAYS
                                              for x in pb['domestic_water_makeup_km3_year']],
                   sky_catchment_l_per_resident_per_mm=roam['rain_catchment_l_per_resident_per_mm_day']),
        sky=dict(envelopes_per_billion={k: 1e9 / v for k, v in aer['residents_per_envelope_by_district'].items()},
                 projected_share_per_billion_50t=roam['projected_share_of_moon'],
                 envelope_radius_km=aer['envelope_radius_m'] / 1e3,
                 spacing_diameters=wol['assumptions']['spacing_diameters'],
                 hydrogen_gt_per_billion=roam['hydrogen_gt'], gas_gt_per_billion=roam['lifting_gas_mixture_gt'],
                 control_tw_per_billion=[roam['control_tw_by_duty']['0.01'], roam['control_tw_by_duty']['0.1']],
                 duty=[min(x['powered_duty_fraction'] for x in roaming),
                       max(x['powered_duty_fraction'] for x in roaming)],
                 powered_m_s=min(x['relative_air_speed_m_s'] for x in roaming),
                 trips_tw_per_billion=roam['trips_down_tw'],
                 hydrogen_remake_tw_per_billion_per_percent=roam['hydrogen_upkeep_per_percent_a_year']['remake_tw'],
                 hydrogen_ppm_per_billion_per_percent=roam['hydrogen_upkeep_per_percent_a_year'][
                     'air_hydrogen_ppm_held'],
                 hydrogen_loss_per_year=aer['hydrogen_upkeep']['loss_bracket_per_year'],
                 food_lift_tw_per_billion=fl['gw_at_efficiency'] / fl['residents_aloft_billion'] / 1e3,
                 lift_j_kg_10_km=lift_row['ideal_j_kg'], lift_efficiency=lift_row['efficiency'],
                 storm_width_km=aer['storm_width_km'],
                 gather_by_band=gather, occupancy_by_band=occupancy,
                 wind_10_km=wind_spread(winds, 10.0),
                 sink_10_km_e_folding_days=fixed10[0]['e_folding_days'] if fixed10 else None,
                 routes_sharing_a_cell_10_km=[min(gathering10['within_cell_share'][1:]),
                                              max(gathering10['within_cell_share'][1:])],
                 aerophyte_cover=se['cover_basis']['central'], aerophyte_cover_range=[
                     se['cover_basis']['range'][0], se['cover_basis']['range'][-1]],
                 aerophyte_passes_check=next(e for e in se['tenants']['contact']['encounters']
                                             if e['diameter_m'] == 100.0 and e['cover'] == 0.01),
                 breathes_like_earth_m=sky_town['earth_equivalent_m']),
        food=dict(settlement_m2=food['basis']['diet']['settlement_m2'],
                  power_made_kw=systems['power_made_food']['total_kw_per_resident'],
                  microbial_animal_part_kw=systems['microbial_protein_animal_part']['power_tw'],
                  fog_glass=dict(glass_km2=fog_glass['glass_km2'], share_of_zone=fog_glass['share_of_zone'],
                                 from_seas_km3=fog_glass['from_seas_km3_year'],
                                 power_tw=fog_glass['power_tw']['total'],
                                 pipe_tw=fog_glass['power_tw']['water'],
                                 n2o_kt_n=fog_glass['nitrogen']['n2o']['kt_n_year']),
                  watered_fog_fields=dict(land_km2=irr['fog_desert']['land_km2'],
                                          withdrawal_km3=irr['fog_desert']['withdrawal_km3_year']),
                  eaten_kg_day=fl['food_kg_per_person_day'],
                  fog_desert_mean_height_km=zones['fog_desert']['height_above_sea_km_mean'],
                  fog_desert_km2=zones['fog_desert']['dry_land_km2'],
                  fog_desert_sea_km=zones['fog_desert']['km_to_sea_or_great_lake_mean']),
        light=dict(per_person=ll['headlines']['per_person_m2'],
                   metropolis_town_lit={des: designs[des]['metropolis_m2_per_person_lit_like_a_town_of_1e4']
                                        for des in DESIGNS},
                   dark_belt_km2=ll['night_life']['dark_night_requirement']['dark_belt_km2'],
                   walking_hours={h: v['needs']['walking']['hours_below_mean']
                                  for h, v in ll['needs']['homes'].items()},
                   never_below_walking={h: v['needs']['walking']['share_never_below']
                                        for h, v in ll['needs']['homes'].items()},
                   check_spread=ll['headlines']['shares_all_on_the_land']['shielded']['spread'],
                   summit_darkest_lux=places['places']['summit port']['night']['natural_darkest_lux'],
                   highland_box_darkest_lux=places['places']['highland box']['night']['natural_darkest_lux'],
                   polar_plain_darkest_lux=places['places']['polar plain']['night']['natural_darkest_lux'],
                   summit_flashes=places['places']['summit port']['lightning']['central'],
                   highland_box_flashes=places['places']['highland box']['lightning']['central']),
        phosphorus=dict(design_kg_per_person=sp['people']['ecology_design_target']['lost_kg_per_person'],
                        sediment_kg_per_kg=sp['routes']['dredging']['sediment_kg_per_kg_p'],
                        fuel_mj_per_kg=sp['routes']['dredging']['fuel_mj_per_kg_p'],
                        delta_particles_kg_yr=sp['routes']['dredging']['supply_kg_yr']['lake_delta_particles'],
                        sea_mouth_particles_kg_yr=sp['routes']['dredging']['supply_kg_yr']['sea_mouth_particles'],
                        struvite_kg_per_kg=sp['distribution']['carrier_kg_per_kg_p']['struvite'][0],
                        lift_mj_per_kg_per_km=sp['distribution']['lift_mj_per_kg_p_per_km'],
                        sewage_recoverable_kg_per_billion=sp['routes']['wastewater']['recoverable_per_billion_kg'],
                        town_return_needed=sp['people']['ecology_design_target']['town_return_needed_at_low_end'],
                        natural=sp['natural_cycle']['decent'],
                        waters_life_ceiling_kg_yr=sp['routes']['waters_life_ceiling_kg_yr']['base'],
                        near_side_sea=dict(particulate_kg_yr=near_side['particulate_kg_yr'],
                                           dissolved_kg_yr=near_side['dissolved_kg_yr'],
                                           mg_m3_1000_years=near_side_mg),
                        land_reach=sp['distribution']['land_reach']),
        diet_kcal_day=food['basis']['diet']['intake_kcal_day'],
        destinations=destinations, origin_middle_lowland=origin,
        crews=dict(sky_boat_seats=access['sky_boat']['seats'],
                   sky_boat_range_with_climb_km=access['sky_boat']['range_with_climb_km'],
                   sky_boat_range_km=access['sky_boat']['range_km'],
                   sky_boat_kwh_50_km_with_climb=next(x['sky_boat']['kwh_with_climb'] for x in access['trips']
                                                      if x['km'] == 50.0),
                   sky_ferry_route_km=access['sky_ferry']['route_km'],
                   regional_ship_route_km=access['regional_ship']['route_km'],
                   regional_ship_wh_per_passenger_km=access['regional_ship']['wh_per_passenger_km'],
                   middle_wet_mean_lat_deg=zone_lat,
                   km_to_45_degrees=(45.0 - zone_lat) * math.pi / 180.0 * MOON_RADIUS / 1e3),
        tubes=dict(totals=land['tubes']['capacity']['totals'], sealing=land['tubes']['sealing'],
                   reference_people_million=land['tubes']['sealing']['residents_energy']['people_million'],
                   by_candidate=land['tubes']['capacity']['by_candidate'],
                   rows={x['key']: dict(lat=x['lat'], lon=x['lon'], cls=x['class'], type=x['type'], name=x['name'])
                         for x in land['tubes']['on_the_atlas']['rows']}),
        earth=dict(within_boundaries_billion=pop['food']['within_boundaries_billion'],
                   farmland_km2_per_billion=wol['settings']['earth']['per_billion']['plant_rich_farmland_km2'],
                   farmland_share_per_billion=wol['settings']['earth']['per_billion']['share_of_todays_farmland'],
                   k_per_billion_thermal=[pop['heat']['earth']['efficient']['k_per_billion_thermal'],
                                          pop['heat']['earth']['high_income']['k_per_billion_thermal']],
                   rectenna_km2_per_tw=pop['earth_service']['rectenna_km2_per_tw'],
                   collector_km2_per_tw_delivered=pop['earth_service']['collector_km2_per_tw_delivered'],
                   heat_per_w_delivered=pop['earth_service']['heat_on_earth_per_w_delivered'],
                   thermal_plant_heat_per_w=pop['earth_service']['thermal_plant_heat_per_w'],
                   earth_w_m2_per_tw=pop['heat']['earth_w_m2_per_tw'],
                   earth_k_per_w_m2=pop['heat']['earth_k_per_w_m2'],
                   transfer_days=daily['earth']['transfer_days'],
                   shared_cases_billion=[min(case_earth), max(case_earth)],
                   shared_cases_elsewhere_billion=[min(case_elsewhere), max(case_elsewhere)]),
        array=dict(shield_gt_per_billion=wol['settings']['array_habitats']['per_billion']['shield_gt'],
                   farm_km2_per_billion=wol['settings']['array_habitats']['per_billion']['farm_km2'],
                   power_tw_per_billion=wol['settings']['array_habitats']['per_billion']['power_tw'],
                   commute_kw=wol['settings']['array_habitats']['per_billion']['commute_kw_per_commuter'],
                   moon_w_m2_per_tw=pop['array']['moon_w_m2_per_tw_in_orbit'],
                   fleet_mass_gt=comp['supply']['fleet_mass_gt'],
                   matter_stream_gt_per_year=pop['array']['matter_stream_gt_per_year'],
                   night_half_pw=comp['supply']['night_half_intercepts_pw'],
                   collector_km2_per_tw=pop['array']['collector_km2_per_tw'],
                   shield_t_per_person=pop['array']['shield_t_per_person']),
        transcribed={k: v['value'] for k, v in TRANSCRIBED.items()},
        literature={k: v['value'] for k, v in LIT.items()},
        guesses={k: v['value'] for k, v in GUESSES.items()},
    )


# ---------------------------------------------------------------------------------------------------------------
# The Moon: placement, land bands, food, energy, water, phosphorus, light, wild land and sky, traffic.

def kept_key(kept):
    """A land rule's key in the products; a plain kept share is kept_30 or kept_50."""
    return getattr(kept, 'key', None) or f'kept_{round(kept * 100)}'


def hybrid_share(kept):
    """The share of each zone a land rule holds as parks and gardens."""
    return getattr(kept, 'hybrid', 0.0)


def food_block(food, mix, moon, kept):
    """A mix's case for `moon` billion residents at a kept share, at its full range when both ends fit and at its
    feasible range otherwise; None when even its low-cost end overflows."""
    case = food['mixes'][mix]['cases'][str(int(moon))][kept_key(kept)]
    feasible = case['feasible']
    if not feasible[0]:
        return None
    if feasible[1]:
        block, reach = case, 1.0
    else:
        block, reach = case['feasible_range'], case['feasible_to']
    return dict(block=block, reach=reach, both_ends=bool(feasible[1]))


def food_town_km2(B, surface_billion, end, reach):
    """The food study's town allowance on the highlands, 25-100 m2 a person, moved with its other inputs."""
    lo, hi = B['food']['settlement_m2']
    m2 = lo + (reach if end == 1 else 0.0) * (hi - lo)
    return surface_billion * 1e9 * m2 / 1e6


def placement(B, ground, combo, rule='spread', dry_rule='land', kept=0.3):
    """The ground's people by setting for one combination of the guesses (each index 0 low, 1 high). The dry
    highland people spread over the dry high ground's bands by land, or fill the bands poleward of 45 degrees first
    (dry_rule 'poleward'); the scattered communities spread over the other homes by land, or four times as densely on
    the twilit and Earthlit homes (rule 'twilit')."""
    gs = B['guesses']
    hs = gs['highland_share'][combo[0]]
    ts = gs['tower_share'][combo[1]]
    sm2 = B['scattered']['m2_per_person'][combo[2]]
    met = B['metropolis']['billion']
    domes = B['undersea']['billion']
    high = hs * ground
    dry = high - met
    tower = ts * dry
    towns = dry - tower
    scattered = ground - high - domes
    km2_per_billion_dry = 1e9 * (ts / B['towers']['people_per_km2'] + (1.0 - ts) / B['towns']['people_per_km2'])
    land = B['highlands']['dry_bands_km2']
    north_land = land['45-60'] + land['60-90']
    if dry_rule == 'land':
        share = {k: v / sum(land.values()) for k, v in land.items()}
    else:
        need = dry * km2_per_billion_dry
        north = min(need, north_land * (1.0 - kept))
        frac_north = north / need if need > 0 else 0.0
        share = {'30-45': 1.0 - frac_north, '45-60': frac_north * land['45-60'] / north_land,
                 '60-90': frac_north * land['60-90'] / north_land}
    dry_people = {k: dry * share[k] for k in DRY_BANDS}
    dry_area = {k: dry_people[k] * km2_per_billion_dry for k in DRY_BANDS}
    homes = B['scattered']['homes_km2']
    weight = {h: (B['scattered']['mix_weight'] if (rule == 'twilit' and h in TWILIT_HOMES) else 1.0) * homes[h]
              for h in homes}
    w_sum = sum(weight.values())
    by_home = {h: scattered * weight[h] / w_sum for h in homes}
    share_0_30 = sum(weight[h] * B['scattered']['home_share_0_30'][h] for h in homes) / w_sum
    share_middle = sum(weight[h] * B['scattered']['home_share_middle'][h] for h in homes) / w_sum
    areas = dict(metropolis=B['metropolis']['dry_land_km2'],
                 tower_districts=tower * 1e9 / B['towers']['people_per_km2'],
                 highland_towns=towns * 1e9 / B['towns']['people_per_km2'],
                 scattered=scattered * 1e9 * sm2 / 1e6, undersea_and_tubes=0.0)
    people = dict(metropolis=met, tower_districts=tower, highland_towns=towns, scattered=scattered,
                  undersea_domes=domes)
    land = dict(low_0_30=areas['scattered'] * share_0_30, middle=areas['scattered'] * share_middle)
    land['high'] = areas['scattered'] - land['low_0_30'] - land['middle']
    return dict(highland_share=hs, tower_share=ts, scattered_m2=sm2, people=people, dry_highland=dry,
                dry_people_by_band=dry_people, dry_area_by_band=dry_area, scattered_by_home=by_home,
                scattered_share_0_30=share_0_30, scattered_share_middle=share_middle, areas_km2=areas,
                scattered_land_km2=land, homestead_land_km2=0.0, homestead_rule='zone',
                settled_km2=sum(areas.values()))


def land_bands(B, blk, reach, place, e, kept, wet_rule='area'):
    """Every band of the dry land outside the kept share, with the food fields and this study's settlement together.
    The fields' band split follows the food study's fill: wettest band first, and in the 30-45 degree band the
    highland share of what is available there, its own towns spread over the highland bands by land. Rain-fed fields
    in the middle-wet zone cannot move to the fog desert, so the zone is a band of its own: the scattered communities'
    middle-band land (30-48.6 degrees) covers it by area (wet_rule 'area'), or keeps off its fields where they fill it
    and settles the fog-desert side of the homes (wet_rule 'dry_side')."""
    lb = B['land_bands']
    hb, dry_b, zones = lb['highland_km2'], lb['dry_land_km2'], lb['zones_km2']
    out_share = 1.0 - kept
    f = blk.get('fields')
    rain_belt = f['field_km2_by_zone']['rain_belt'][e] if f else 0.0
    middle_wet = middle_wet_fields(B, blk, e, kept)
    highland_fields = f['field_km2_by_zone']['highlands'][e] if f else 0.0
    fog_watered = f['fog_desert_watered_km2'][e] if f else 0.0
    if blk['land_km2']['glass'][e] > 0.0:
        raise ValueError('the band check is built for mixes without glass')
    towns_fs = food_town_km2(B, blk['people_billion']['surface'][e], e, reach)
    reserved = {k: towns_fs * hb[k] / sum(hb.values()) for k in FOOD_BANDS}
    high_30_45 = max(hb['30-45'] * out_share - reserved['30-45'], 0.0)
    fields_30_45 = middle_wet * high_30_45 / (zones['middle_wet'] * out_share) if middle_wet > 0 else 0.0
    # the fields stay within what the food study's fill had in the 0-30 degree bands, which its product's three
    # significant figures can overstate by a few thousand km2
    fill_0_30 = sum(max(dry_b[k] * out_share - reserved[k], 0.0) for k in ('0-15', '15-30'))
    fields_0_30 = min(rain_belt + highland_fields - fields_30_45, fill_0_30)
    a = place['areas_km2']
    sl = place['scattered_land_km2']
    middle_land = sl['middle']
    mw_cap = zones['middle_wet'] * out_share
    on_mw = middle_land * lb['middle_wet_share_of_middle_homes']
    if wet_rule == 'dry_side':
        on_mw = min(on_mw, max(mw_cap - middle_wet, 0.0))
    elif wet_rule != 'area':
        raise ValueError(f'unknown wet rule {wet_rule}')
    homes = place['homestead_land_km2']          # C': homesteads live among the middle-wet fields
    rule = place['homestead_rule']
    if rule == 'zone':                           # beside the fields, on the zone
        homes_on_mw, homes_elsewhere = homes, 0.0
    elif rule == 'within_fields':                # within the fields they work, taking a share of the fields' land
        homes_on_mw, homes_elsewhere = 0.0, 0.0
    else:
        raise ValueError(f'unknown homestead rule {rule}')
    grazing = blk['grazing_km2'][e]
    use = dict(all_land_0_30=fields_0_30 + a['metropolis'] + sl['low_0_30'],
               highlands_30_45=fields_30_45 + place['dry_area_by_band']['30-45'],
               highlands_45_90=place['dry_area_by_band']['45-60'] + place['dry_area_by_band']['60-90'],
               middle_wet_zone=middle_wet + on_mw + homes_on_mw,
               fog_desert_and_polar=(fog_watered + grazing + middle_land - on_mw + sl['high'] + homes_elsewhere))
    cap = dict(all_land_0_30=(dry_b['0-15'] + dry_b['15-30']) * out_share,
               highlands_30_45=hb['30-45'] * out_share,
               highlands_45_90=(hb['45-60'] + hb['60-75'] + hb['75-90']) * out_share,
               middle_wet_zone=mw_cap,
               fog_desert_and_polar=(zones['fog_desert'] + zones['polar']) * out_share)
    overfill = {k: max(use[k] - cap[k], 0.0) for k in LAND_BANDS}
    dry_used = (place['dry_area_by_band']['30-45'] + use['highlands_45_90'] + fields_30_45)
    dry_side_cap = (lb['middle_band_homes_km2'] - zones['middle_wet']) * out_share
    return dict(use_km2=use, capacity_km2=cap, use_share={k: use[k] / cap[k] for k in LAND_BANDS},
                overfill_km2=overfill, fields_on_dry_high_ground_km2=fields_30_45,
                fields_on_middle_wet_km2=middle_wet, scattered_middle_band_km2=middle_land,
                scattered_on_middle_wet_km2=on_mw, homestead_land_on_middle_wet_km2=homes_on_mw,
                homestead_share_of_middle_wet_fields=homes / middle_wet if middle_wet > 0.0 else 0.0,
                middle_homes_dry_side_used_share=(middle_land - on_mw) / dry_side_cap,
                dry_high_ground_used_share=dry_used / (B['highlands']['dry_km2'] * out_share),
                dry_high_ground_settled_share=(dry_used - fields_30_45) / (B['highlands']['dry_km2'] * out_share))


def lit_sky(B, place, aloft, end, design):
    """Land under a light-polluted sky at the night's darkest moment (living light's per-person areas for one
    lighting design), each home capped at its land. The dry highland people at 30-45 degrees take the far side's
    twilit-night places' figure, the home at the same latitudes and side; poleward of 45 degrees towns light none."""
    pp = B['light']['per_person'][design]
    homes = B['scattered']['homes_km2']
    dry = place['dry_people_by_band']['30-45'] * 1e9 * pp['homes']['far_side_places'][end] / 1e6
    scattered = sum(min(homes[h], place['scattered_by_home'][h] * 1e9 * pp['homes'][h][end] / 1e6) for h in homes)
    up = aloft * 1e9 * pp['aloft'][end] / 1e6
    met_town = B['metropolis']['billion'] * 1e9 * B['light']['metropolis_town_lit'][design][end] / 1e6
    met_earth = B['metropolis']['billion'] * 1e9 * pp['metropolis'][end] / 1e6
    dry_land = B['moon']['dry_land_km2']
    total_town, total_earth = met_town + dry + scattered + up, met_earth + dry + scattered + up
    return dict(metropolis_lit_as_a_town=met_town, metropolis_at_earth_area=met_earth, dry_highlands=dry,
                scattered=scattered, aloft=up, total_metropolis_as_a_town=total_town,
                total_metropolis_at_earth_area=total_earth,
                share_of_dry_land_metropolis_as_a_town=total_town / dry_land,
                share_of_dry_land_metropolis_at_earth_area=total_earth / dry_land,
                metropolis_share_of_dark_belt_as_a_town=met_town / B['light']['dark_belt_km2'],
                metropolis_share_of_dark_belt_at_earth_area=met_earth / B['light']['dark_belt_km2'])


def moon_corner(B, food, moon, aloft, kept, supply, combo, end, rule='spread', dry_rule='land', wet_rule='area',
                mix=CHOSEN_MIX, middle=None):
    """Every account of one Moon case at one corner: one combination of the guesses, every input at `end`. With
    `middle` (C'), the middle latitudes' lowland people keep homesteads among the farms there and the rest move."""
    ground = moon - aloft
    place = placement(B, ground, combo, rule, dry_rule, kept)
    fb = food_block(food, mix, moon, kept)
    blk, reach = fb['block'], fb['reach']
    e = end
    if middle is not None:
        place = move_middle(B, place, blk, reach, e, kept, wet_rule, middle)
    dry_land = B['moon']['dry_land_km2']
    towns_fs = food_town_km2(B, blk['people_billion']['surface'][e], e, reach)
    food_land = blk['land_km2']['total'][e] - towns_fs
    food_high = blk['land_km2']['highlands'][e] - towns_fs
    grazing = blk['grazing_km2'][e]
    settled = place['settled_km2']
    used = food_land + grazing + settled
    a = place['areas_km2']
    bands = land_bands(B, blk, reach, place, e, kept, wet_rule)
    land = dict(
        settled_km2=settled, settled_by_setting_km2=a, food_km2=food_land, grazing_km2=grazing,
        used_km2=used, used_share_of_dry_land=used / dry_land,
        wild_share_outside_kept=1.0 - used / (dry_land * (1.0 - kept)),
        parks_km2=dry_land * hybrid_share(kept), parks_share_of_dry_land=hybrid_share(kept),
        highlands_used_share_outside_kept=(food_high + a['metropolis'] + a['tower_districts'] + a['highland_towns'])
        / (B['highlands']['all_km2'] * (1.0 - kept)),
        dry_high_ground_settled_share=bands['dry_high_ground_settled_share'],
        dry_high_ground_used_share=bands['dry_high_ground_used_share'],
        fields_on_dry_high_ground_km2=bands['fields_on_dry_high_ground_km2'])
    # energy, on one basis and as electricity, and heat by component
    en = B['energy']
    band = en['band_kw'][e]
    resident_primary = moon * band
    fusion_input = en['heat_per_w']['fusion'][e]
    if supply == 'fusion':
        resident_heat = moon * en['heat_tw_per_billion_fusion'][e]
        elec_factor = fusion_input
    else:
        resident_heat = moon * band * en['placement_power_from_orbit']
        elec_factor = en['heat_per_w']['array']
    sky = B['sky']
    loss = sky['hydrogen_loss_per_year'][e] * 100.0
    sky_parts = dict(control=aloft * sky['control_tw_per_billion'][e], trips_down=aloft * sky['trips_tw_per_billion'],
                     hydrogen_remade=aloft * sky['hydrogen_remake_tw_per_billion_per_percent'][e] * loss,
                     food_lift=aloft * sky['food_lift_tw_per_billion'])
    sky_elec = sum(sky_parts.values())
    w = B['water']
    rec = w['recovery_cut'] if e == 0 else 1.0
    dry_gw = sum(place['dry_people_by_band'][k] * w['dry_gw_per_billion_by_band'][k][e] for k in DRY_BANDS)
    water_gw = (B['metropolis']['billion'] * w['metropolis_gw_per_billion'][e] * rec + dry_gw * rec
                + place['people']['scattered'] * w['scattered_makeup_gw_per_billion'][e])
    water_elec = water_gw / 1e3
    food_elec = blk['power_tw']['total'][e]
    electricity = food_elec + sky_elec + water_elec
    heat_parts = dict(residents=resident_heat, food=blk['heat'][supply]['heat_tw'][e],
                      sky_control=sky_parts['control'] * elec_factor,
                      sky_hydrogen_remade=sky_parts['hydrogen_remade'] * elec_factor,
                      sky_trips_and_food_lift=(sky_parts['trips_down'] + sky_parts['food_lift']) * elec_factor,
                      water=water_elec * elec_factor)
    heat = sum(heat_parts.values())
    dim = en['dimming_percent_per_tw'][e]
    close = moon * en['close_computing_tw_per_billion']
    energy = dict(resident_primary_tw=resident_primary, food_tw=food_elec, sky_tw=sky_elec, sky_parts_tw=sky_parts,
                  water_lift_tw=water_elec, electricity_tw=electricity,
                  primary_tw=resident_primary + electricity * fusion_input,
                  heat_tw=heat, heat_by_component_tw=heat_parts,
                  heat_share={k: v / heat for k, v in heat_parts.items()},
                  heat_w_m2=heat * en['w_m2_per_tw'], added_dimming_percent=heat * dim,
                  close_computing_tw=close, close_computing_dimming_percent=close * elec_factor * dim,
                  hydrogen_ppm_held=aloft * sky['hydrogen_ppm_per_billion_per_percent'] * loss)
    # water
    dom = w['domestic_km3_per_billion'][e]
    dry_use = place['dry_highland'] * dom
    water = dict(domestic_gross_km3=moon * dom, lift_gw=water_gw, dry_highlands_lift_gw=dry_gw * rec,
                 dry_highlands_domestic_km3=dry_use, dry_highlands_runoff_km3=B['highlands']['dry_runoff_km3_year'],
                 dry_highlands_runoff_over_use=B['highlands']['dry_runoff_km3_year'] / dry_use,
                 food_rain_km3=blk['water_km3_year']['rain_used'][e],
                 food_from_seas_km3=blk['water_km3_year']['withdrawn_from_seas'][e])
    # phosphorus
    p = B['phosphorus']
    need_kt = blk['phosphorus_mt_year']['lost'][e] * 1e3
    catch_kt = p['delta_particles_kg_yr'][1 - e] / 1e6
    struvite_people = place['dry_highland'] + aloft
    phosphorus = dict(
        need_kt=need_kt, eaten_kt=blk['phosphorus_mt_year']['eaten'][e] * 1e3,
        dredged_sediment_gt=need_kt * 1e6 * p['sediment_kg_per_kg'][e] / 1e12,
        dredging_fuel_gw=need_kt * 1e6 * p['fuel_mj_per_kg'][e] * 1e6 / SECONDS_PER_YEAR / 1e9,
        share_of_delta_catch=need_kt / catch_kt,
        struvite_to_highlands_and_sky_mt=need_kt * struvite_people / moon * p['struvite_kg_per_kg'] / 1e3,
        struvite_lift_mw=(need_kt * 1e6 * place['dry_highland'] / moon * p['lift_mj_per_kg_per_km']['struvite'][e]
                          * B['highlands']['dry_mean_height_km'] * 1e6 / SECONDS_PER_YEAR / 1e6),
        sewage_recoverable_kt=moon * p['sewage_recoverable_kg_per_billion'][e] / 1e6,
        microbial_biomass_kt=blk['phosphorus_mt_year']['in_microbial_biomass'][e] * 1e3)
    light = {des: lit_sky(B, place, aloft, e, des) for des in DESIGNS}
    # sky
    moon_km2 = B['moon']['area_km2']
    shade = aloft * sky['projected_share_per_billion_50t']
    n_env = aloft * sky['envelopes_per_billion']['50.0']
    cell = (sky['spacing_diameters'] * 2.0 * sky['envelope_radius_km']) ** 2
    grid = n_env * cell / moon_km2
    sky_out = dict(towns_50t=n_env, shade_share_of_moon=shade, grid_share_of_moon=grid,
                   wild_sky_outside_kept_by_shade=1.0 - shade / (1.0 - kept),
                   wild_sky_outside_kept_by_grid=1.0 - grid / (1.0 - kept),
                   towns_grid_share_of_their_sky=grid / (1.0 - kept),
                   hydrogen_gt=aloft * sky['hydrogen_gt_per_billion'],
                   lifting_gas_gt=aloft * sky['gas_gt_per_billion'])
    buildable = all(v <= 1.0 for v in bands['overfill_km2'].values())
    return dict(buildable=buildable, people=place['people'],
                people_unplaced_billion=ground - sum(place['people'].values()),
                n2o_kt_n_per_billion=blk['n2o']['kt_n_year'][e] / moon,
                dry_highland_billion=place['dry_highland'], dry_people_by_band=place['dry_people_by_band'],
                scattered_share_0_30=place['scattered_share_0_30'],
                highland_share=place['highland_share'], tower_share=place['tower_share'],
                scattered_m2_per_person=place['scattered_m2'], bands=bands, land=land, energy=energy, water=water,
                phosphorus=phosphorus, light=light, sky=sky_out, moved=place.get('moved'))


def destination_spare_km2(B, blk, reach, place, e, kept, wet_rule='area'):
    """Spare land outside the kept share, after the fields and the settlement, in each destination's own land."""
    bands = land_bands(B, blk, reach, place, e, kept, wet_rule)
    hb, out = B['land_bands']['highland_km2'], 1.0 - kept
    use, cap = bands['use_km2'], bands['capacity_km2']
    return dict(high_latitude_lowland=B['scattered']['high_homes_km2'] * out - place['scattered_land_km2']['high'],
                lowland_0_30=cap['all_land_0_30'] - use['all_land_0_30'],
                dry_30_45=cap['highlands_30_45'] - use['highlands_30_45'],
                dry_45_60=hb['45-60'] * out - place['dry_area_by_band']['45-60'],
                dry_60_90=(hb['60-75'] + hb['75-90']) * out - place['dry_area_by_band']['60-90'])


def middle_wet_fields(B, blk, e, kept):
    """The fields on the middle-wet zone, held to the room the food study's fill gives them in the 30-45 degree band:
    the zone's land outside the kept share. The fill reaches that room exactly when it fills the zone, and the food
    product's three significant figures can put its field area up to 0.5% above the room computed here."""
    f = blk.get('fields')
    if not f:
        return 0.0
    lb = B['land_bands']
    hb, dry_b, zones = lb['highland_km2'], lb['dry_land_km2'], lb['zones_km2']
    room = (min(dry_b['30-45'], zones['middle_wet'] + hb['30-45']) - hb['30-45']) * (1.0 - kept)
    return min(max(f['field_km2_by_zone']['middle_wet'][e], 0.0), room)


def move_middle(B, place, blk, reach, e, kept, wet_rule, spec):
    """C': the scattered communities of the 30-48.6 degree lowland homes leave, except homesteads among the farms
    there, at the guessed density on the middle-wet zone's fields; the rest go to the destinations in `spec['order']`,
    each taking what its spare land outside the kept share holds at its own density."""
    fields_mw = middle_wet_fields(B, blk, e, kept)
    density = GUESSES['homestead_people_per_km2']['value'][spec['h_end']]
    middle_people = place['people']['scattered'] * place['scattered_share_middle']
    homesteads = min(density * fields_mw / 1e9, middle_people)
    moved = middle_people - homesteads
    m2 = place['scattered_m2']
    new = copy.deepcopy(place)
    new['scattered_land_km2']['middle'] = 0.0
    new['homestead_land_km2'] = homesteads * 1e9 * m2 / 1e6
    new['homestead_rule'] = spec.get('homes_on', 'zone')
    share_mid = B['scattered']['home_share_middle']
    by_home = new['scattered_by_home']
    for h in by_home:
        by_home[h] *= 1.0 - share_mid[h]
    by_home['far_side_places'] += homesteads   # lit like the far side's twilit-night places, as at 30-45 degrees
    ts = new['tower_share']
    dry_km2_per_b = 1e9 * (ts / B['towers']['people_per_km2'] + (1.0 - ts) / B['towns']['people_per_km2'])
    lowland_km2_per_b = 1e9 * m2 / 1e6
    homes_0_30 = DESTINATIONS['lowland_0_30']['homes']
    land_0_30 = {h: B['scattered']['homes_km2'][h] * B['scattered']['home_share_0_30'][h] for h in homes_0_30}
    placed = {k: 0.0 for k in DESTINATIONS}
    left = moved
    for dest in spec['order']:
        if left <= 0.0:
            break
        kind = DESTINATIONS[dest]['kind']
        per_b = dry_km2_per_b if kind == 'dry' else lowland_km2_per_b
        spare = destination_spare_km2(B, blk, reach, new, e, kept, wet_rule)[dest]
        take = min(left, max(spare, 0.0) / per_b)
        if take <= 0.0:
            continue
        if dest == 'high_latitude_lowland':
            new['scattered_land_km2']['high'] += take * per_b
            by_home['high_latitude_towns'] += take
        elif dest == 'lowland_0_30':
            new['scattered_land_km2']['low_0_30'] += take * per_b
            for h in homes_0_30:
                by_home[h] += take * land_0_30[h] / sum(land_0_30.values())
        else:
            b = DESTINATIONS[dest]['band']
            new['dry_people_by_band'][b] += take
            new['dry_area_by_band'][b] += take * per_b
            new['dry_highland'] += take
            new['people']['tower_districts'] += ts * take
            new['people']['highland_towns'] += (1.0 - ts) * take
            new['people']['scattered'] -= take
            a = new['areas_km2']
            a['tower_districts'] += ts * take * 1e9 / B['towers']['people_per_km2']
            a['highland_towns'] += (1.0 - ts) * take * 1e9 / B['towns']['people_per_km2']
            a['scattered'] -= take * lowland_km2_per_b
        placed[dest] = take
        left -= take
    new['settled_km2'] = sum(new['areas_km2'].values())
    scattered_land = new['areas_km2']['scattered']
    if scattered_land > 0.0:
        new['scattered_share_0_30'] = new['scattered_land_km2']['low_0_30'] / scattered_land
        new['scattered_share_middle'] = new['homestead_land_km2'] / scattered_land
    new['moved'] = dict(middle_billion=middle_people, homesteads_billion=homesteads, moved_billion=moved,
                        placed_billion=placed, unplaced_billion=left, homestead_land_km2=new['homestead_land_km2'],
                        middle_wet_fields_km2=fields_mw)
    return new


def rank_destinations(B, spare_km2):
    """The destinations ranked on the five measures, 1 the best and ties sharing their mean rank: spare land (most
    first), water energy with nothing recovered (least), comfort (the mean of the CM1 and GCM ranks), rain and
    lightning (the mean of their ranks, least first) and night lighting (the mean of the walking-lamp and lit-sky
    ranks, least first). Their order is by the mean of the five, ties to the dry highlands, the author's land."""
    D = B['destinations']
    keys = list(D)

    def mid(v):
        return (v[0] + v[1]) / 2.0

    def ranks(values, most_first):
        order = sorted(keys, key=lambda k: -values[k] if most_first else values[k])
        out, i = {}, 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and math.isclose(values[order[j + 1]], values[order[i]], rel_tol=1e-9,
                                                      abs_tol=1e-12):
                j += 1
            for k in order[i:j + 1]:
                out[k] = (i + j) / 2.0 + 1.0
            i = j + 1
        return out

    r = dict(spare=ranks(spare_km2, True),
             water=ranks({k: mid(D[k]['water_gw_per_billion_nothing_recovered']) for k in keys}, False),
             cm1=ranks({k: mid(D[k]['comfort_cm1_strict_hours']) for k in keys}, True),
             gcm=ranks({k: mid(D[k]['comfort_gcm_strict_hours']) for k in keys}, True),
             rain=ranks({k: mid(D[k]['rain_mm_day']) for k in keys}, False),
             flashes=ranks({k: mid(D[k]['flashes_per_km2_year']) for k in keys}, False),
             lamps=ranks({k: mid(D[k]['walking_lamp_hours']) for k in keys}, False),
             lit=ranks({k: mid(D[k]['lit_m2_per_person']) for k in keys}, False))
    measures = {}
    for k in keys:
        five = dict(spare_land=r['spare'][k], water=r['water'][k], comfort=(r['cm1'][k] + r['gcm'][k]) / 2.0,
                    rain_and_lightning=(r['rain'][k] + r['flashes'][k]) / 2.0,
                    night_lighting=(r['lamps'][k] + r['lit'][k]) / 2.0)
        five['mean'] = sum(five.values()) / 5.0
        measures[k] = five
    order = sorted(keys, key=lambda k: (measures[k]['mean'], 0 if D[k]['kind'] == 'dry' else 1, -spare_km2[k]))
    return dict(ranks=measures, order=order, spare_land_km2=spare_km2)


def destination_spare_for_ranking(B, food, moon, aloft, kept):
    """The least spare land each destination keeps over C's corners as spread, with the middle latitudes' lowland
    land counted as moved out (the spare land the move can use)."""
    least = {k: math.inf for k in DESTINATIONS}
    for combo in COMBOS:
        place = placement(B, moon - aloft, combo, kept=kept)
        place['scattered_land_km2']['middle'] = 0.0
        fb = food_block(food, CHOSEN_MIX, moon, kept)
        for e in (0, 1):
            spare = destination_spare_km2(B, fb['block'], fb['reach'], place, e, kept)
            for k in least:
                least[k] = min(least[k], max(spare[k], 0.0))
    return least


def traffic(B, aloft, kept, district='50.0'):
    """Sky towns at 10 km: their number, spacing over the sky clear of terrain outside the kept share, gathering by
    the winds, and encounters on the shared tracks by the sky-ecology study's rule: density x speed x swept width.
    Per-town rates are the busiest band's; fleet totals take the density each town meets, weighted by occupancy."""
    sky = B['sky']
    n_env = aloft * sky['envelopes_per_billion'][district]
    clear = B['moon']['sky_clear_of_terrain_km2']
    area = B['moon']['sky_clear_total_km2'] * (1.0 - kept)
    n_even = n_env / area
    gather = max(sky['gather_by_band'].values())
    busiest = max(sky['gather_by_band'], key=sky['gather_by_band'].get)
    n_busy = n_even * gather
    occ = sky['occupancy_by_band']
    n_met = n_env * sum(occ[b] ** 2 / (clear[b] * (1.0 - kept)) for b in occ)
    d_km = 2.0 * sky['envelope_radius_km']
    swept_km = 2.0 * sky['spacing_diameters'] * d_km           # centres within three diameters, either side
    km_day = JULIAN_DAY / 1e3
    drift = GUESSES['relative_drift_m_s']['value']
    duty = sky['duty']
    speeds = dict(drifting=drift, powered=[f * sky['powered_m_s'] for f in duty])
    rates = {k: [n_busy * swept_km * v * km_day for v in vs] for k, vs in speeds.items()}
    rates['uncorrelated_winds'] = n_busy * swept_km * sky['wind_10_km']['uncorrelated_relative_m_s'] * km_day
    mean_rates = {k: [n_met * swept_km * v * km_day for v in vs] for k, vs in speeds.items()}
    fleet = {k: [n_env * x / 2.0 for x in vs] for k, vs in mean_rates.items()}
    aero_d_km, margin_km = 0.1, 0.05
    aero_swept = d_km + aero_d_km + 2.0 * margin_km
    covers = [sky['aerophyte_cover_range'][0], sky['aerophyte_cover'], sky['aerophyte_cover_range'][1]]
    aero = {f'{c:g}': [c / (math.pi * (aero_d_km / 2.0) ** 2) * aero_swept * v * km_day for v in drift]
            for c in covers}
    return dict(towns=n_env, spacing_even_km=n_even ** -0.5, spacing_busiest_band_km=n_busy ** -0.5,
                busiest_band=busiest, busiest_band_gathering=gather, met_density_over_even=n_met / n_even,
                towns_under_a_storm=n_busy * math.pi * (sky['storm_width_km'] / 2) ** 2,
                conflict_distance_km=sky['spacing_diameters'] * d_km,
                encounters_per_town_day_busiest_band=rates, encounters_per_town_day_mean=mean_rates,
                encounters_fleet_day=fleet, aerophyte_passes_per_town_day=aero)


def assemble(rows_by_supply, selection):
    """The merged accounts over the selected corners."""
    if not selection:
        return None
    base = merge([rows_by_supply['fusion'][i] for i in selection])
    arr = merge([rows_by_supply['array'][i] for i in selection])
    return dict(people=base['people'], people_unplaced_billion=base['people_unplaced_billion'],
                n2o_kt_n_per_billion=base['n2o_kt_n_per_billion'], dry_highland_billion=base['dry_highland_billion'],
                dry_people_by_band=base['dry_people_by_band'], scattered_share_0_30=base['scattered_share_0_30'],
                highland_share=base['highland_share'], tower_share=base['tower_share'],
                scattered_m2_per_person=base['scattered_m2_per_person'], bands=base['bands'], land=base['land'],
                water=base['water'], phosphorus=base['phosphorus'], light=base['light'], sky=base['sky'],
                energy=dict(fusion=base['energy'], array=arr['energy']), moved=base.get('moved'))


def moon_case(B, food, moon, aloft, mix=CHOSEN_MIX):
    out = {}
    for kept in KEPT:
        corners = [(combo, e) for combo in COMBOS for e in (0, 1)]
        rows = {s: [moon_corner(B, food, moon, aloft, kept, s, c, e, mix=mix) for c, e in corners] for s in SUPPLIES}
        ok = [i for i, row in enumerate(rows['fusion']) if row['buildable']]
        overfilling = [dict(highland_share=rows['fusion'][i]['highland_share'],
                            tower_share=rows['fusion'][i]['tower_share'],
                            scattered_m2_per_person=rows['fusion'][i]['scattered_m2_per_person'],
                            end=corners[i][1], overfill_km2=rows['fusion'][i]['bands']['overfill_km2'])
                       for i in range(len(corners)) if i not in ok]
        entry = dict(corners=len(corners), buildable_corners=len(ok), overfilling_corners=overfilling,
                     full=assemble(rows, list(range(len(corners)))),
                     feasible=assemble(rows, ok) if len(ok) < len(corners) else None, remedy=None,
                     traffic={d: traffic(B, aloft, kept, d) for d in DISTRICT_T})
        if len(ok) < len(corners):
            rem = {s: [moon_corner(B, food, moon, aloft, kept, s, c, e, rule='twilit', dry_rule='poleward',
                                   wet_rule='dry_side', mix=mix) for c, e in corners] for s in SUPPLIES}
            rem_ok = [i for i, row in enumerate(rem['fusion']) if row['buildable']]
            full = assemble(rem, list(range(len(corners))))
            parts = {}
            for name, kw in (('poleward', dict(dry_rule='poleward')),
                             ('poleward_and_twilit', dict(rule='twilit', dry_rule='poleward')),
                             ('all_three', dict(rule='twilit', dry_rule='poleward', wet_rule='dry_side'))):
                parts[name] = sum(moon_corner(B, food, moon, aloft, kept, 'fusion', c, e, mix=mix, **kw)['buildable']
                                  for c, e in corners)
            entry['remedy'] = dict(rule='the dry highland towns poleward of 45 degrees first, the scattered '
                                        'communities four times as dense on the twilit and Earthlit homes, and the '
                                        'middle band\'s scattered communities on the fog-desert side of their homes '
                                        'where the fields fill the middle-wet zone',
                                   buildable_corners=len(rem_ok), buildable_corners_by_parts=parts,
                                   bands=full['bands'],
                                   dry_people_by_band=full['dry_people_by_band'],
                                   scattered_share_0_30=full['scattered_share_0_30'],
                                   land=full['land'], water=full['water'],
                                   light=dict(shielded=full['light']['shielded'], living=full['light']['living']),
                                   energy=dict(fusion=full['energy']['fusion']))
        out[kept_key(kept)] = entry
    return out


CPRIME_VARIANTS = {
    'as_placed': dict(),
    'poleward_and_twilit': dict(rule='twilit', dry_rule='poleward'),
    'homes_within_the_fields': dict(rule='twilit', dry_rule='poleward', homes_on='within_fields'),
    'to_the_dry_highlands': dict(dest='dry'),
}


def cprime_case(B, food, moon, aloft, rankings, c_case):
    """The author's C': C with the middle latitudes' lowland people moved by the ranking, homesteads staying among
    the farms there. Corners: the 8 placement combinations, the homestead density at both ends and every input at
    both ends (32). Besides the case as placed, three variants: the remedy's poleward and twilit placements, the same
    with the homesteads' homes counted within the fields they work, and the moved people sent to the dry highland
    cities instead."""
    out = {}
    for kept in KEPT:
        order = rankings[kept_key(kept)]['order']
        dry_order = [k for k in order if DESTINATIONS[k]['kind'] == 'dry']
        corners = [(combo, h, e) for combo in COMBOS for h in (0, 1) for e in (0, 1)]
        variants = {}
        for name, spec in CPRIME_VARIANTS.items():
            kw = {k: v for k, v in spec.items() if k in ('rule', 'dry_rule')}
            middle = dict(order=dry_order if spec.get('dest') == 'dry' else order,
                          homes_on=spec.get('homes_on', 'zone'))
            rows = {s: [moon_corner(B, food, moon, aloft, kept, s, c, e, mix=CHOSEN_MIX,
                                    middle=dict(middle, h_end=h), **kw) for c, h, e in corners] for s in SUPPLIES}
            ok = [i for i, row in enumerate(rows['fusion']) if row['buildable']]
            overfilling = [dict(highland_share=rows['fusion'][i]['highland_share'],
                                tower_share=rows['fusion'][i]['tower_share'],
                                scattered_m2_per_person=rows['fusion'][i]['scattered_m2_per_person'],
                                homestead_people_per_km2=GUESSES['homestead_people_per_km2']['value'][corners[i][1]],
                                end=corners[i][2], overfill_km2=rows['fusion'][i]['bands']['overfill_km2'])
                           for i in range(len(corners)) if i not in ok]
            variants[name] = dict(corners=len(corners), buildable_corners=len(ok), overfilling_corners=overfilling,
                                  full=assemble(rows, list(range(len(corners)))),
                                  feasible=assemble(rows, ok) if 0 < len(ok) < len(corners) else None)
        main = variants['as_placed']
        out[kept_key(kept)] = dict(
            corners=main['corners'], buildable_corners=main['buildable_corners'],
            overfilling_corners=main['overfilling_corners'], full=main['full'], feasible=main['feasible'],
            remedy=None, traffic=c_case[kept_key(kept)]['traffic'], destination_order=order,
            variants={k: v for k, v in variants.items() if k != 'as_placed'})
    return out


def food_menu(B, food, moon, aloft):
    """Every mix for `moon` billion residents: whether it fits by the food study's zones, how far, and its land,
    power, water, N2O and phosphorus over its feasible range, with the dry land left wild once this study's
    settlement replaces the food study's towns. The band check of moon.cases applies to the chosen mix."""
    menu = {}
    for mix in food['mixes']:
        row = {}
        for kept in KEPT:
            fb = food_block(food, mix, moon, kept)
            case = food['mixes'][mix]['cases'][str(int(moon))][kept_key(kept)]
            if fb is None:
                row[kept_key(kept)] = dict(buildable=False, overflowing=case.get('overflowing'))
                continue
            blk, reach = fb['block'], fb['reach']
            wild = []
            for e in (0, 1):
                towns_fs = food_town_km2(B, blk['people_billion']['surface'][e], e, reach)
                for combo in COMBOS:
                    place = placement(B, moon - aloft, combo)
                    used = blk['land_km2']['total'][e] - towns_fs + blk['grazing_km2'][e] + place['settled_km2']
                    wild.append(1.0 - used / (B['moon']['dry_land_km2'] * (1.0 - kept)))
            row[kept_key(kept)] = dict(
                buildable=True, both_ends=fb['both_ends'], reach=reach,
                food_land_km2=[blk['land_km2']['total'][e] - food_town_km2(B, blk['people_billion']['surface'][e],
                                                                           e, reach) for e in (0, 1)],
                grazing_km2=blk['grazing_km2'], glass_km2=blk['land_km2']['glass'],
                power_tw=blk['power_tw']['total'], from_seas_km3=blk['water_km3_year']['withdrawn_from_seas'],
                n2o_ppb_per_year=blk['n2o']['ppb_per_year'], phosphorus_lost_mt=blk['phosphorus_mt_year']['lost'],
                fish_sludge_p_mt=blk['phosphorus_mt_year']['to_fish_sludge'],
                dimming_fusion_percent=blk['heat']['fusion']['added_dimming_percent'],
                wild_share_outside_kept=[min(wild), max(wild)])
        menu[mix] = row
    return menu


# ---------------------------------------------------------------------------------------------------------------
# The tubes and the undersea cities joined.

def tubes_and_undersea(B, heritage):
    """Each undersea community with its nearest lava tube to which the land-cities inventory gives a capacity.
    Candidates within TIE_KM of the nearest distance tie, and every tied tube joins: the inventory can give two
    distinct tubes one position."""
    t = B['tubes']
    sites = {s['name']: s for s in heritage['sites']}
    capacity = {k: v for k, v in t['by_candidate'].items() if k in t['rows']}
    named = {k: v for k, v in t['rows'].items() if v['type'].startswith('named rille')}
    pairs = []
    for name, depth in zip(B['undersea']['sites'], B['undersea']['depth_m']):
        s = sites[name]
        near = sorted((gc_km(s['lat'], s['lon'], t['rows'][k]['lat'], t['rows'][k]['lon']), k) for k in capacity)
        km = near[0][0]
        tied = [k for d_km, k in near if d_km - km <= TIE_KM]
        rille = min((gc_km(s['lat'], s['lon'], v['lat'], v['lon']), k) for k, v in named.items())
        pairs.append(dict(site=name, lat=s['lat'], lon=s['lon'], depth_m=depth, km=km,
                          tubes=[capacity[k]['name'] for k in tied], tube_keys=tied,
                          tube_classes=[capacity[k]['cls'] for k in tied],
                          tube_million_small=[sum(capacity[k]['small'][i] for k in tied) for i in (0, 1)],
                          tube_million_large=[sum(capacity[k]['large'][i] for k in tied) for i in (0, 1)],
                          tied=len(tied) > 1, nearest_named_rille=t['rows'][rille[1]]['name'],
                          named_rille_km=rille[0]))
    joined = [p for p in pairs if p['km'] <= JOIN_KM]
    keys = sorted({k for p in joined for k in p['tube_keys']})
    small = [sum(capacity[k]['small'][i] for k in keys) for i in (0, 1)]
    large = [sum(capacity[k]['large'][i] for k in keys) for i in (0, 1)]
    domes = len(joined) * B['undersea']['per_dome'] / 1e6
    mapped, grav = t['totals']['mapped']['by_class'], t['totals']['with_gravity']['by_class']
    sea_coast = dict(
        mapped=dict(small=[mapped['under the sea']['small'][i] + mapped['straddles the coast']['small'][i]
                           for i in (0, 1)],
                    large=[mapped['under the sea']['large'][i] + mapped['straddles the coast']['large'][i]
                           for i in (0, 1)]),
        with_gravity=dict(small=[grav['under the sea']['small'][i] + grav['straddles the coast']['small'][i]
                                 for i in (0, 1)],
                          large=[grav['under the sea']['large'][i] + grav['straddles the coast']['large'][i]
                                 for i in (0, 1)]))
    sealing = {row['setting']: row for row in t['sealing']['settings']}
    pumping = {}
    for setting in ('the water table at the ground', 'under 500 m of sea', 'under 2,000 m of sea'):
        rings = sealing[setting]['by_ring']
        pumping[setting] = dict(seepage_m3_s=[rings['tight']['seepage_m3_s'], rings['loose']['seepage_m3_s']],
                                pumping_mw=[rings['tight']['pumping_mw'], rings['loose']['pumping_mw']])
    return dict(
        undersea_domes=len(B['undersea']['sites']), dome_residents=B['undersea']['per_dome'],
        domes_billion=B['undersea']['billion'], tie_km=TIE_KM, join_km=JOIN_KM, pairs=pairs,
        ties=[dict(site=p['site'], tubes=p['tubes']) for p in pairs if p['tied']],
        joined=dict(sites=[p['site'] for p in joined], tubes=keys, tube_million_small=small,
                    tube_million_large=large, domes_million=domes,
                    all_million=[small[0] + domes, large[1] + domes]),
        sea_and_coast_tubes_million=sea_coast,
        all_tubes_million=dict(mapped=t['totals']['mapped']['all'],
                               dry=[t['totals']['mapped']['by_class']['dry']['small'][0],
                                    t['totals']['mapped']['by_class']['dry']['large'][1]],
                               dry_with_gravity=[grav['dry']['small'][0], grav['dry']['large'][1]]),
        reference_tube_pumping=pumping)


# ---------------------------------------------------------------------------------------------------------------
# Synergies among places.

def synergies(B, cases, tubes):
    food, ph = B['food'], B['phosphorus']
    out = {}
    # (a) a submerged tube's seepage and the sea floor's sediment
    seep = {}
    ref_people = B['tubes']['reference_people_million']
    for setting, row in tubes['reference_tube_pumping'].items():
        km3 = [x * SECONDS_PER_YEAR / 1e9 for x in row['seepage_m3_s']]
        fed = [km3[i] / food['watered_fog_fields']['withdrawal_km3'][1 - i] * 1e3 for i in (0, 1)]
        km2 = [km3[i] / food['watered_fog_fields']['withdrawal_km3'][1 - i]
               * food['watered_fog_fields']['land_km2'][1 - i] for i in (0, 1)]
        p_t = [km3[i] * 1e9 * ph['near_side_sea']['mg_m3_1000_years'][i] / 1e9 for i in (0, 1)]   # t of P a year
        seep[setting] = dict(seepage_km3_year=km3, fields_watered_km2=km2, people_fed_million=fed,
                             phosphorus_t_year=p_t, pumping_mw=row['pumping_mw'])
    city_need_t = [ref_people[0] * 1e6 * ph['design_kg_per_person'][0] / 1e3,
                   ref_people[1] * 1e6 * ph['design_kg_per_person'][1] / 1e3]
    ns = ph['near_side_sea']
    out['dewatering'] = dict(reference_tube_people_million=ref_people, by_setting=seep,
                             tube_city_phosphorus_need_t_year=city_need_t,
                             near_side_sea_mg_m3_1000_years=ns['mg_m3_1000_years'])
    supplied = [ns['particulate_kg_yr'][0] / ph['design_kg_per_person'][1] / 1e6,
                ns['particulate_kg_yr'][1] / ph['design_kg_per_person'][0] / 1e6]
    sc = tubes['sea_and_coast_tubes_million']
    out['near_side_sea_sediment'] = dict(
        particulate_kt=[x / 1e6 for x in ns['particulate_kg_yr']],
        dissolved_kt=[x / 1e6 for x in ns['dissolved_kg_yr']],
        million_supplied_at_design=supplied, sea_and_coast_tubes_million=sc,
        supply_over_tube_people_mapped=[supplied[0] / sc['mapped']['large'][1], supplied[1] / sc['mapped']['small'][0]],
        supply_over_tube_people_with_gravity=[supplied[0] / sc['with_gravity']['large'][1],
                                              supplied[1] / sc['with_gravity']['small'][0]])
    # (b) the lake deltas' phosphorus for the highland and lowland farms
    out['deltas'] = {k: cases[k]['kept_30']['full']['phosphorus'] for k in cases}
    out['deltas_reach'] = ph['land_reach']
    out['sea_floor_dredging'] = dict(
        sea_mouth_particles_kt=[x / 1e6 for x in ph['sea_mouth_particles_kg_yr']],
        billion_supplied_at_design=[ph['sea_mouth_particles_kg_yr'][0] / (ph['design_kg_per_person'][1] * 1e9),
                                    ph['sea_mouth_particles_kg_yr'][1] / (ph['design_kg_per_person'][0] * 1e9)])
    # (c) fog-desert glass drinking piped sea water for the dry highland cities
    g = food['fog_glass']
    rise_m = (B['highlands']['dry_mean_height_km'] - food['fog_desert_mean_height_km']) * 1e3
    freight_gw = (1e9 * food['eaten_kg_day'] / JULIAN_DAY * MOON_SURFACE_GRAVITY * rise_m
                  / B['sky']['lift_efficiency'] / 1e9)
    per_case = {}
    for k, c in cases.items():
        dry = c['kept_30']['full']['dry_highland_billion']
        per_case[k] = dict(dry_highland_billion=dry,
                           glass_km2=[dry[0] * g['glass_km2'][0], dry[1] * g['glass_km2'][1]],
                           fog_desert_share=[dry[0] * g['glass_km2'][0] / food['fog_desert_km2'],
                                             dry[1] * g['glass_km2'][1] / food['fog_desert_km2']],
                           sea_water_km3=[dry[0] * g['from_seas_km3'][0], dry[1] * g['from_seas_km3'][1]],
                           power_tw=[dry[0] * g['power_tw'][0], dry[1] * g['power_tw'][1]],
                           pipes_tw=[dry[0] * g['pipe_tw'][0], dry[1] * g['pipe_tw'][1]],
                           n2o_kt_n=[dry[0] * g['n2o_kt_n'][0], dry[1] * g['n2o_kt_n'][1]],
                           freight_up_gw=[dry[0] * freight_gw, dry[1] * freight_gw])
    chosen_n2o = {k: c['kept_30']['full']['n2o_kt_n_per_billion'] for k, c in cases.items()}
    out['fog_desert_glass'] = dict(chosen_mix_n2o_kt_n_per_billion=chosen_n2o,
                                   per_billion=dict(glass_km2=g['glass_km2'], sea_water_km3=g['from_seas_km3'],
                                                    power_tw=g['power_tw'], pipes_tw=g['pipe_tw'],
                                                    n2o_kt_n=g['n2o_kt_n'], freight_up_gw=freight_gw,
                                                    rise_km=rise_m / 1e3),
                                   fog_desert_km2=food['fog_desert_km2'],
                                   fog_desert_mean_km_to_sea=food['fog_desert_sea_km'], by_case=per_case)
    # (d) the dark belt and the twilit high latitudes
    out['dark_belt'] = dict(
        dark_belt_km2=B['light']['dark_belt_km2'],
        summit_darkest_lux=B['light']['summit_darkest_lux'],
        highland_box_darkest_lux=B['light']['highland_box_darkest_lux'],
        polar_plain_darkest_lux=B['light']['polar_plain_darkest_lux'],
        walking_lamp_hours_by_home=B['light']['walking_hours'],
        never_below_walking_share_by_home=B['light']['never_below_walking'],
        dry_highlands_m2_per_person_30_45=B['light']['per_person']['shielded']['homes']['far_side_places'],
        metropolis_lit_km2={des: dict(
            as_a_town=[B['metropolis']['billion'] * 1e9 * x / 1e6 for x in B['light']['metropolis_town_lit'][des]],
            at_earth_area=[B['metropolis']['billion'] * 1e9 * x / 1e6
                           for x in B['light']['per_person'][des]['metropolis']]) for des in DESIGNS})
    # (e) the array's power behind Earth's food off the land
    e = B['earth']
    pm = food['power_made_kw']
    out['earth_food_off_the_land'] = dict(
        per_billion_beyond_boundaries=dict(
            power_made_tw=pm, microbial_animal_part_tw=food['microbial_animal_part_kw'],
            rectenna_km2=[pm[0] * e['rectenna_km2_per_tw'], pm[1] * e['rectenna_km2_per_tw']],
            array_collectors_km2=[pm[0] * e['collector_km2_per_tw_delivered'][0],
                                  pm[1] * e['collector_km2_per_tw_delivered'][1]],
            heat_on_earth_tw=[pm[0] * e['heat_per_w_delivered'], pm[1] * e['heat_per_w_delivered']],
            heat_on_earth_tw_thermal=[pm[0] * e['thermal_plant_heat_per_w'], pm[1] * e['thermal_plant_heat_per_w']],
            farmland_spared_km2=e['farmland_km2_per_billion'],
            share_of_night_half=[pm[0] / (B['array']['night_half_pw'] * 1e3),
                                 pm[1] / (B['array']['night_half_pw'] * 1e3)]),
        within_boundaries_billion=e['within_boundaries_billion'])
    # (f) more places whose weakness serves
    met = B['metropolis']
    dom = B['water']['domestic_km3_per_billion']
    out['summit_rain'] = dict(home_runoff_km3=met['home_runoff_km3_year'],
                              domestic_km3=[met['billion'] * dom[0], met['billion'] * dom[1]],
                              runoff_over_need=[met['home_runoff_km3_year'] / (met['billion'] * dom[1]),
                                                met['home_runoff_km3_year'] / (met['billion'] * dom[0])],
                              rain_mm_day=met['rain_mm_day'], lake_share=met['lake_share'],
                              wet_highlands_km2=B['highlands']['wet_km2'],
                              wet_highlands_rain_mm_day=B['highlands']['wet_rain_mm_day'])
    sky_l = B['water']['makeup_l_per_resident_day']
    catch = B['water']['sky_catchment_l_per_resident_per_mm']
    out['sky_rain'] = dict(catchment_l_per_resident_per_mm=catch, makeup_l_per_resident_day=sky_l,
                           rain_mm_day_that_closes_it=[sky_l[0] / catch, sky_l[1] / catch])
    t = B['transcribed']
    out['array_shield_as_water'] = dict(share_of_build_percent=t['water_shield_share_of_build'],
                                        habitat_gt_per_billion=t['array_habitat_gt_per_billion'])
    return out


# ---------------------------------------------------------------------------------------------------------------
# Earth, the array and the colonies.

def earth_block(B, billion):
    """Earth at `billion` people: food within four planetary boundaries to 10.2 billion and off the land beyond,
    its energy on the one basis (the residents' primary energy and the off-land food's electricity at fusion's input,
    2.0-3.03 W a watt) and its heat."""
    e, en = B['earth'], B['energy']
    beyond = max(billion - e['within_boundaries_billion'], 0.0)
    fed = billion - beyond
    pm = B['food']['power_made_kw']
    off = [beyond * pm[0], beyond * pm[1]]
    fusion_input = en['heat_per_w']['fusion']
    heat_array = [x * e['heat_per_w_delivered'] for x in off]
    heat_thermal = [x * e['thermal_plant_heat_per_w'] for x in off]
    k_w = e['earth_w_m2_per_tw'] * e['earth_k_per_w_m2']
    return dict(
        billion=billion, within_boundaries=billion <= e['within_boundaries_billion'], beyond_boundaries_billion=beyond,
        off_land_food_tw=off,
        off_land_rectenna_km2=[off[0] * e['rectenna_km2_per_tw'], off[1] * e['rectenna_km2_per_tw']],
        off_land_array_collectors_km2=[off[0] * e['collector_km2_per_tw_delivered'][0],
                                       off[1] * e['collector_km2_per_tw_delivered'][1]],
        off_land_heat_on_earth_tw=dict(array=heat_array, thermal=heat_thermal),
        off_land_warming_k=dict(array=[x * k_w for x in heat_array], thermal=[x * k_w for x in heat_thermal]),
        farmland_km2=[fed * e['farmland_km2_per_billion'][0], fed * e['farmland_km2_per_billion'][1]],
        farmland_share_of_today=[fed * e['farmland_share_per_billion'][0], fed * e['farmland_share_per_billion'][1]],
        energy_tw=[billion * en['band_kw'][i] + off[i] * fusion_input[i] for i in (0, 1)],
        warming_k_thermal=[billion * e['k_per_billion_thermal'][0], billion * e['k_per_billion_thermal'][1]])


def array_block(B, billion):
    a, t, en = B['array'], B['transcribed'], B['energy']
    shield = [billion * a['shield_gt_per_billion'][0], billion * a['shield_gt_per_billion'][1]]
    habitat = billion * t['array_habitat_gt_per_billion']
    energy = [billion * en['band_kw'][0], billion * en['band_kw'][1]]
    return dict(
        billion=billion, shield_gt=shield, habitat_gt=habitat, habitat_over_fleet=habitat / a['fleet_mass_gt'],
        shield_over_fleet=[shield[0] / a['fleet_mass_gt'], shield[1] / a['fleet_mass_gt']],
        stream_kg_s=[billion * t['array_growth_kg_s_per_billion'][0], billion * t['array_growth_kg_s_per_billion'][1]],
        days_of_build_stream=[shield[0] / a['matter_stream_gt_per_year'] * JULIAN_YEAR_DAYS,
                              shield[1] / a['matter_stream_gt_per_year'] * JULIAN_YEAR_DAYS],
        farm_km2=billion * a['farm_km2_per_billion'], torus_power_tw=billion * a['power_tw_per_billion'],
        energy_tw=energy, heat_on_moon_w_m2=[energy[0] * a['moon_w_m2_per_tw'], energy[1] * a['moon_w_m2_per_tw']],
        collectors_km2=[energy[0] * a['collector_km2_per_tw'][0], energy[1] * a['collector_km2_per_tw'][1]],
        commute_gw_per_million_commuters=dict(weekly=a['commute_kw']['weekly'],
                                              each_lunar_cycle=a['commute_kw']['each_lunar_cycle']))


def colonies_per_billion(B):
    """Costs per billion residents of the budding colonies from the project's figures and the two sources read."""
    en, k, c = B['energy'], B['literature']['kalpana_one'], B['literature']['city_npc']
    band = en['band_kw']
    pm = B['food']['power_made_kw']
    fusion_input = en['heat_per_w']['fusion']
    torus_farm_m2 = B['array']['farm_km2_per_billion'] * 1e6 / 1e9
    lunar_kcal = B['diet_kcal_day']
    t = B['transcribed']
    mirrors = {}
    for name in PLACES_AU:
        au = t[name][0]
        mirrors[name] = dict(au=au, w_m2=SOLAR_CONSTANT / au ** 2, w_m2_in_the_note=t[f'{name}_sunlight'],
                             collecting_area_m2=torus_farm_m2 * au ** 2,
                             mirrors_beyond_the_farm_m2=torus_farm_m2 * (au ** 2 - 1.0),
                             light_minutes_from_the_moon=[(au - 1.0) * AU / SPEED_OF_LIGHT / 60.0,
                                                          (au + 1.0) * AU / SPEED_OF_LIGHT / 60.0])
    kalpana_t = k['total_mass_t'] / k['residents']
    return dict(
        living_tw=band,
        food=dict(power_made_tw=pm, lamp_lit_agriculture_tw=k['agriculture_kw_per_resident'],
                  mars_greenhouse_m2_per_person=c['greenhouse_ha'] * 1e4 / c['people'],
                  mars_greenhouse_m2_per_person_at_lunar_diet=c['greenhouse_ha'] * 1e4 / c['people']
                  * lunar_kcal / c['diet_kcal_day'],
                  sunlit_farm_m2_per_person_at_1_au=torus_farm_m2, sunlit_farm_mirrors=mirrors),
        energy_tw=dict(primary_with_power_made_food=[band[i] + pm[i] * fusion_input[i] for i in (0, 1)],
                       mars_electricity_twe=[c['colony_kwe_per_member'], c['city_kwe_per_capita']],
                       kalpana_one=k['kw_per_resident']),
        habitat_gt=dict(stanford_torus=B['transcribed']['array_habitat_gt_per_billion'],
                        kalpana_one=kalpana_t, kalpana_one_shield_t_m2=k['shield_t_m2'],
                        kalpana_one_living_m2=k['living_m2_per_resident']),
        mars_lava_tube_district=dict(people=c['tunnel_people'], km=c['tunnel_km'], diameter_m=c['tunnel_diameter_m']),
        capacity_credited=False)


# ---------------------------------------------------------------------------------------------------------------
# The distributions.

def fmt(x):
    """A number for the generated sentences: three significant figures, thousands separated."""
    v = float(f'{x:.3g}')
    return f'{v:,.0f}' if v >= 100 else f'{v:g}'


ASSUMPTIONS = {
    'common': [
        "The Moon's ground people live mostly in the highlands: the summit metropolis's 100 million round the port, "
        'the rest of the highland people on the dry high ground poleward of 30 degrees in tower districts and '
        'highland towns, and the others in scattered communities spread over the other homes by their land.',
        "Food is the food study's power-heavy mix: vegetables and fruit in lamp-lit vertical farms, the animal part as "
        'microbial protein, the rest on rain-fed fields; the land bands are checked with its fields and the '
        'settlement together.',
        'Sky towns are the inhabited-volume reference envelopes at 10 km with 50 t a resident, powered 1-10% of the '
        'time at 5 m/s, losing 1-10% of their hydrogen a year.',
        "Residents use 2-6 kW of primary energy with industry in orbit; fusion supplies the Moon's power, its heat "
        'answered by added dimming at 0.0087-0.0139% a terawatt; electricity counts at fusion\'s input.',
        "Earth holds 9 billion, within the 10.2 billion that eat inside four planetary boundaries."],
    'A': ['The Moon keeps 2 billion on the ground and 2 billion aloft, the least the author set for the sky.',
          'The array grows to 7 billion residents and the colonies to 5 billion.'],
    'B': ['The Moon holds 3 billion on the ground and 3 billion aloft.',
          'The array holds 6 billion and the colonies 4 billion.'],
    'C': ['The Moon holds 4 billion on the ground and 4 billion aloft, the most the author set for the sky.',
          'The array holds 5 billion and the colonies 3 billion.'],
    'D': ['The Moon holds 6 billion on the ground and 2 billion aloft, keeping the sky at its least occupied.',
          'The array holds 8 billion and the colonies 5 billion.'],
    "C'": ['The Moon holds 4 billion on the ground and 4 billion aloft, as in C, the metropolis at its 100 million.',
           "The scattered communities of the 30-48.6 degree lowland homes move to the best-ranked home with room, "
           'except homesteads among the farms there at 2-20 people a km2 of field (a design guess).',
           'The array holds 5 billion and the colonies 3 billion.'],
}
UNCERTAINTIES_COMMON = [
    "Whether people grow up healthy at a sixth of Earth's gravity, and how much time at Earth's gravity a lifetime "
    'needs.',
    "The crop model's water use, which sets the rain-fed land within a factor of two, and farm yields below potential.",
    'Comfort over the land: CM1 and the GCM disagree on the air near the ground, and both put the most comfortable '
    'land on high ground poleward of 30 degrees.',
    "The sky towns' hulls, gas barriers and traffic rules: the hydrogen lost each year sets the top of the Moon's "
    'heat, and the relative drift of neighbours on a shared track is a guess.',
    "The colonies' capacity."]


VARIANT_TEXT = dict(poleward_and_twilit="with the remedy's poleward and twilit placements",
                    homes_within_the_fields=("with the remedy's poleward and twilit placements and the homesteads' "
                                             'homes within their fields'),
                    to_the_dry_highlands='with the moved people sent to the dry highland cities')


def rule_summary(k):
    """One land rule's numbers for a distribution: which corners fit, as placed, by the remedy and by C's variants, and
    the land, water, food power, heat, sky and traffic over all corners."""
    f = k['full']
    t50 = k['traffic']['50.0']
    return dict(buildable_corners=k['buildable_corners'], of=k['corners'],
                remedy_buildable_corners=k['remedy']['buildable_corners'] if k.get('remedy') else None,
                variants_buildable_corners=({n: v['buildable_corners'] for n, v in k['variants'].items()}
                                            if 'variants' in k else None),
                used_share_of_dry_land=f['land']['used_share_of_dry_land'],
                parks_share_of_dry_land=f['land']['parks_share_of_dry_land'],
                wild_share_outside_kept=f['land']['wild_share_outside_kept'],
                food_from_seas_km3=f['water']['food_from_seas_km3'], food_tw=f['energy']['fusion']['food_tw'],
                heat_tw=f['energy']['fusion']['heat_tw'],
                added_dimming_percent=f['energy']['fusion']['added_dimming_percent'],
                towns_grid_share_of_their_sky=f['sky']['towns_grid_share_of_their_sky'],
                spacing_busiest_band_km=t50['spacing_busiest_band_km'],
                encounters_fleet_day=t50['encounters_fleet_day'])


def distribution(B, key, spec, moon_cases):
    case = moon_cases[spec.get('case', f"{int(spec['moon'])}_{int(spec['aloft'])}")]
    k30, k50, kc = case['kept_30'], case['kept_50'], case['wild_30_hybrid_20']
    f30, f50 = k30['full'], k50['full']
    feas50 = k50['feasible'] or f50
    earth = earth_block(B, spec['earth'])
    array = array_block(B, spec['array'])
    col = colonies_per_billion(B)
    colonies_tw = [spec['colonies'] * x for x in col['energy_tw']['primary_with_power_made_food']]
    moon_primary = f30['energy']['fusion']['primary_tw']
    total = [earth['energy_tw'][i] + moon_primary[i] + array['energy_tw'][i] + colonies_tw[i] for i in (0, 1)]
    people_total = spec['earth'] + spec['moon'] + spec['array'] + spec['colonies']
    names = {'all_land_0_30': 'the 0-30 degree land', 'highlands_30_45': 'the 30-45 degree highlands',
             'highlands_45_90': 'the highlands poleward of 45 degrees', 'middle_wet_zone': 'the middle-wet zone',
             'fog_desert_and_polar': 'the fog desert and the polar lowland'}
    def overfill_text(k, where):
        over = [names[b] for b in LAND_BANDS if any(o['overfill_km2'][b] > 1.0 for o in k['overfilling_corners'])]
        listed = ', '.join(over[:-1]) + ' or ' + over[-1] if len(over) > 1 else ''.join(over)
        return (f"{where} its high-cost end overfills {listed} at {k['corners'] - k['buildable_corners']} of "
                f"{k['corners']} corners" if over else f'every corner fits {where}')
    half_text = overfill_text(k50, 'outside a kept half')
    compromise_text = overfill_text(kc, "under the author's compromise")
    fits = [k for k, v in kc.get('variants', {}).items() if v['buildable_corners'] == v['corners']]
    if fits and kc['buildable_corners'] < kc['corners']:
        compromise_text += ', and every corner fits ' + VARIANT_TEXT[fits[0]]
    uncertain = {
        'A': f"The array's 7 billion need {fmt(array['shield_gt'][0])}-{fmt(array['shield_gt'][1])} Gt of shield, "
             "the largest continuing material stream after the build, and the colonies' 5 billion rest on "
             'figures from designs never built.',
        'B': f"It sits between A and C: the Moon's heat and lit sky at 6 billion, and the array's "
             f"{fmt(array['shield_gt'][0])}-{fmt(array['shield_gt'][1])} Gt of shield.",
        'C': "Four billion aloft shade more of the Moon than the aerophytes' central 1% cover and put "
             f"{fmt(k30['traffic']['50.0']['towns'])} towns on the wind tracks; the Moon's heat and the dimming that "
             f'answers it are the largest of the four, and {half_text}.',
        'D': 'It asks the most of the dry high ground, the land and the array: six billion on the ground put the '
             f"most people on the dry high ground and in the scattered communities, the array's 8 billion need "
             f"{fmt(array['shield_gt'][0])}-{fmt(array['shield_gt'][1])} Gt of shield, and {half_text}.",
        "C'": "C's heat, sky and array, with the middle latitudes' lowland towns moved to the best-ranked home with "
              f"room and homesteads among the farms there; {half_text}; {compromise_text}."}
    return dict(
        name=spec['name'],
        people_billion=dict(earth=spec['earth'], moon=spec['moon'], moon_ground=spec['moon'] - spec['aloft'],
                            moon_aloft=spec['aloft'], array=spec['array'], colonies=spec['colonies'],
                            total=people_total),
        moon=dict(case=spec.get('case', f"{int(spec['moon'])}_{int(spec['aloft'])}"),
                  heat_tw=dict(fusion=f30['energy']['fusion']['heat_tw'], array=f30['energy']['array']['heat_tw']),
                  added_dimming_percent=dict(fusion=f30['energy']['fusion']['added_dimming_percent'],
                                             array=f30['energy']['array']['added_dimming_percent']),
                  heat_share_fusion=f30['energy']['fusion']['heat_share'],
                  primary_tw=moon_primary, electricity_tw=f30['energy']['fusion']['electricity_tw'],
                  resident_primary_tw=f30['energy']['fusion']['resident_primary_tw'],
                  wild_dry_land_outside_kept=dict(kept_30=f30['land']['wild_share_outside_kept'],
                                                  kept_50=f50['land']['wild_share_outside_kept'],
                                                  kept_50_feasible=feas50['land']['wild_share_outside_kept']),
                  buildable_corners=dict(kept_30=k30['buildable_corners'], kept_50=k50['buildable_corners'],
                                         of=k30['corners']),
                  remedy_buildable_corners_kept_50=k50['remedy']['buildable_corners'] if k50['remedy'] else None,
                  variants_buildable_corners_kept_50=({k: v['buildable_corners'] for k, v in k50['variants'].items()}
                                                      if 'variants' in k50 else None),
                  land_rules={kept_key(r): rule_summary(case[kept_key(r)]) for r in KEPT},
                  moved=f30['moved'],
                  settled_km2=f30['land']['settled_km2'],
                  people_by_setting=f30['people'],
                  highlands_used_share=dict(kept_30=f30['land']['highlands_used_share_outside_kept'],
                                            kept_50=f50['land']['highlands_used_share_outside_kept']),
                  dry_high_ground_settled_share=dict(kept_30=f30['land']['dry_high_ground_settled_share'],
                                                     kept_50=f50['land']['dry_high_ground_settled_share']),
                  dry_high_ground_used_share=dict(kept_30=f30['land']['dry_high_ground_used_share'],
                                                  kept_50=f50['land']['dry_high_ground_used_share']),
                  lit_sky_share_of_dry_land=dict(
                      metropolis_as_a_town=f30['light']['shielded']['share_of_dry_land_metropolis_as_a_town'],
                      metropolis_at_earth_area=f30['light']['shielded']['share_of_dry_land_metropolis_at_earth_area'],
                      living_metropolis_as_a_town=f30['light']['living']['share_of_dry_land_metropolis_as_a_town']),
                  metropolis_share_of_dark_belt=dict(
                      as_a_town=f30['light']['shielded']['metropolis_share_of_dark_belt_as_a_town'],
                      at_earth_area=f30['light']['shielded']['metropolis_share_of_dark_belt_at_earth_area']),
                  wild_sky_outside_kept_by_shade=dict(kept_30=f30['sky']['wild_sky_outside_kept_by_shade'],
                                                      kept_50=f50['sky']['wild_sky_outside_kept_by_shade']),
                  shade_share_of_moon=f30['sky']['shade_share_of_moon'],
                  towns=k30['traffic']['50.0']['towns'],
                  encounters_per_town_day_busiest_band=k30['traffic']['50.0']['encounters_per_town_day_busiest_band'],
                  encounters_fleet_day=k30['traffic']['50.0']['encounters_fleet_day'],
                  phosphorus_need_kt=f30['phosphorus']['need_kt'],
                  dredged_sediment_gt=f30['phosphorus']['dredged_sediment_gt'],
                  share_of_delta_catch=f30['phosphorus']['share_of_delta_catch'],
                  water_lift_gw=f30['water']['lift_gw'],
                  hydrogen_ppm_held=f30['energy']['fusion']['hydrogen_ppm_held'],
                  close_computing_tw=f30['energy']['fusion']['close_computing_tw']),
        earth=earth, array=array,
        colonies=dict(billion=spec['colonies'], energy_tw=colonies_tw,
                      energy_tw_lamp_lit=spec['colonies'] * col['energy_tw']['kalpana_one'],
                      habitat_gt_free_space=[spec['colonies'] * col['habitat_gt']['stanford_torus'],
                                             spec['colonies'] * col['habitat_gt']['kalpana_one']]),
        energy_tw=dict(basis='primary energy, with electricity counted at fusion\'s input of 2.0-3.0 W a watt',
                       earth=earth['energy_tw'], moon=moon_primary, array=array['energy_tw'], colonies=colonies_tw,
                       total=total,
                       share_of_night_half=[total[0] / (B['array']['night_half_pw'] * 1e3),
                                            total[1] / (B['array']['night_half_pw'] * 1e3)]),
        assumptions=ASSUMPTIONS['common'] + ASSUMPTIONS[key],
        uncertainties=UNCERTAINTIES_COMMON + [uncertain[key]])


def earth_sensitivity(B, dists):
    """Earth at 8 and 11 billion beside the distributions' 9: its food, energy and heat, and how each total moves."""
    base = earth_block(B, 9.0)
    out = {'9': base}
    for billion in EARTH_SENSITIVITY:
        blk = earth_block(B, billion)
        blk['energy_change_tw'] = [blk['energy_tw'][i] - base['energy_tw'][i] for i in (0, 1)]
        blk['totals_billion'] = {k: d['people_billion']['total'] - 9.0 + billion for k, d in dists.items()}
        out[f'{billion:g}'] = blk
    return out


# ---------------------------------------------------------------------------------------------------------------
# Checks against the bound and committed products.

def checks(B, wol, moon_cases):
    out = {}
    out['lit_sky_rule'] = dict(note='living light\'s spread mix at 4, 6 and 8 billion all on the land, recomputed '
                                    'from its per-person areas', bound=B['light']['check_spread'],
                               recomputed=lit_check_values(B, wol))
    sky = B['sky']
    n = 5.0 * sky['envelopes_per_billion']['50.0']
    out['spacing_rule'] = dict(towns=n, spacing_km=(B['moon']['sky_clear_total_km2'] * 0.7 / n) ** 0.5,
                               ways_of_living_km=wol['whole_system']['total25']['0.5']['spacing_km_whole_sky'][
                                   'kept_30'])
    chk = sky['aerophyte_passes_check']
    d = chk['diameter_m'] / 1e3
    rate = [chk['cover'] / (math.pi * (d / 2) ** 2) * (2 * d + 2 * 0.05) * v * JULIAN_DAY / 1e3
            for v in GUESSES['relative_drift_m_s']['value']]
    out['encounter_rule'] = dict(recomputed=rate, sky_ecology=chk['passes_per_day'])
    out['residents_heat_tw_per_billion'] = B['energy']['heat_tw_per_billion_fusion']
    lb = B['land_bands']
    caps = sum(lb['dry_land_km2'][k] for k in FOOD_BANDS)
    out['bands_partition_the_dry_land_km2'] = dict(bands=caps, zones=sum(lb['zones_km2'].values()),
                                                   moon_dry_land=B['moon']['dry_land_km2'])
    out['people_unplaced_billion'] = {k: v['kept_30']['full']['people_unplaced_billion'] for k, v in moon_cases.items()}
    return out


def lit_check_values(B, wol):
    """Living light's spread shares all on the land, recomputed from its per-person areas and the ways-of-living
    homes (the metropolis at Earth's area a person)."""
    pp = B['light']['per_person']['shielded']
    settings = wol['settings']
    homes = list(pp['homes'])
    land = {h: settings[h]['home_dry_land_km2'] for h in homes}
    total = sum(land.values())
    met = B['metropolis']['billion']
    out = {}
    for b in MOON_BILLION:
        vals = []
        for e in (0, 1):
            lit = met * 1e9 * pp['metropolis'][e] / 1e6
            n_surface = b - met
            for h in homes:
                lit += min(land[h], n_surface * land[h] / total * 1e9 * pp['homes'][h][e] / 1e6)
            vals.append(min(1.0, lit / B['moon']['dry_land_km2']))
        out[str(b)] = vals
    return out


# ---------------------------------------------------------------------------------------------------------------

def build(rounded=True):
    """The product; unrounded values when rounded is False (for writing the README from them once)."""
    rnd = rd if rounded else (lambda o: o)
    d = {k: load_json(KEYS[k], INPUTS[KEYS[k]]) for k in KEYS}
    b = {k: load_bound(k) for k in BOUND if BOUND[k]['schema']}
    industry_text = bound_bytes('industry_md').decode('utf-8')
    B = build_basis(d, b)
    food = d['food']
    moon_cases = {}
    for moon in MOON_BILLION:
        for aloft in ALOFT_BILLION:
            if moon - aloft < 1:
                continue
            moon_cases[f'{moon}_{aloft}'] = moon_case(B, food, moon, aloft)
    twilit = {}
    for moon, aloft in ((4, 2), (6, 3), (8, 4), (8, 2)):
        rows = [moon_corner(B, food, moon, aloft, 0.3, 'fusion', c, e, rule='twilit') for c in COMBOS for e in (0, 1)]
        twilit[f'{moon}_{aloft}'] = merge(rows)['light']
    half = {str(m): moon_cases[f'{m}_{m // 2}'] for m in MOON_BILLION}
    menu = {str(m): food_menu(B, food, m, m // 2) for m in MOON_BILLION}
    tubes = tubes_and_undersea(B, d['heritage'])
    syn = synergies(B, half, tubes)
    syn['twilit_scattering'] = {
        k: dict(spread_km2=moon_cases[k]['kept_30']['full']['light']['shielded']['scattered'],
                twilit_km2=twilit[k]['shielded']['scattered'],
                spread_share=moon_cases[k]['kept_30']['full']['light']['shielded'][
                    'share_of_dry_land_metropolis_as_a_town'],
                twilit_share=twilit[k]['shielded']['share_of_dry_land_metropolis_as_a_town'])
        for k in twilit}
    rank_spare = {kept_key(k): destination_spare_for_ranking(B, food, 8, 4, k) for k in KEPT}
    rankings = {kk: rank_destinations(B, v) for kk, v in rank_spare.items()}
    cprime = cprime_case(B, food, 8, 4, rankings, moon_cases['8_4'])
    dists = {k: distribution(B, k, spec, dict(moon_cases, cprime=cprime)) for k, spec in DISTRIBUTIONS.items()}
    chk = checks(B, d['wol'], moon_cases)
    chk['transcribed_phrases'] = transcribed_check(industry_text)
    here = {kept_key(k): dict(closed=float(k), hybrid=hybrid_share(k)) for k in KEPT}
    chk['land_rules'] = dict(note="the land rules here against the food study's", here=here,
                             food_study=food['basis']['land_rules'], equal=here == food['basis']['land_rules'])
    files = [HERE / f for f in CODE]
    product = dict(
        schema=SCHEMA,
        producer=dict(study='distributions', files={str(f.relative_to(ROOT)): digest(f) for f in files},
                      constants=constants_used(files),
                      inputs={p: digest(ROOT / p) for p in INPUTS},
                      sources=digest(HERE / 'sources.json'),
                      bound_inputs={k: dict(path_at_integration=v['path'], branch=v['branch'], commit=v['commit'],
                                            blob=v['blob'], sha256=v['sha256'], worktree=v['worktree'],
                                            read_from='in-tree',
                                            note=v['note']) for k, v in BOUND.items()}),
        evidence=EVIDENCE, reading_rule=READING_RULE, marks=MARKS,
        assumptions=dict(design_guesses=GUESSES, chosen_mix=CHOSEN_MIX, moon_billion=list(MOON_BILLION),
                         aloft_billion=list(ALOFT_BILLION), supplies=list(SUPPLIES),
                         land_rules={kept_key(k): dict(closed=float(k), hybrid=hybrid_share(k)) for k in KEPT},
                         combinations=[list(c) for c in COMBOS], twilit_homes=list(TWILIT_HOMES),
                         land_bands=list(LAND_BANDS), wet_rules=list(WET_RULES), tie_km=TIE_KM, join_km=JOIN_KM),
        literature={k: dict(value=v['value'], unit=v['unit'], source=v['source']) for k, v in LIT.items()},
        transcribed={k: dict(value=v['value'], phrase=v['phrase'], unit=v['unit']) for k, v in TRANSCRIBED.items()},
        basis=rnd(B),
        moon=dict(cases=rnd(moon_cases), twilit_scattering=rnd(twilit), cprime=rnd(cprime),
                  destination_ranking=rnd(rankings)),
        food_menu=rnd(menu),
        tubes_and_undersea=rnd(tubes),
        synergies=rnd(syn),
        earth=rnd({f'{x:g}': earth_block(B, x) for x in (8.0, 9.0, 10.0, 11.0, 12.0, 16.0)}),
        earth_sensitivity=rnd(earth_sensitivity(B, dists)),
        array=rnd({'per_billion': array_block(B, 1.0)}),
        colonies=rnd(colonies_per_billion(B)),
        distributions=rnd(dists),
        checks=rnd(chk))
    return product


def results():
    return build(rounded=True)


def main():
    product = results()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(dumps(product) + '\n')
    print(f'wrote {OUT.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
