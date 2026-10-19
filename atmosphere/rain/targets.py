"""The three target climates as the rain products read them, from the climate targets' product
(research/studies/climate_targets/results/climate_targets.json).

Until the branches meet, the product is domain/climate's corrected one, read from Git at its commit and checked
against its sha256. At the merge IN_TREE becomes True, and the rain products read the same path in the tree and
record its hash among their inputs. Either way the schema is checked.

The hot world and the lock are the product's targets 'no further dimming' and '295 K', each with its own climate. The
resort stands at 301 K, as the author set it on 18 October (domain/climate's register), and is one of the resort
band's candidates (resort_band.rows), whose row gives its 2 m air, wet bulb, land rain and share of the way between
the two runs. Its sky ladder and its air at the lowest level are the two runs' climates interpolated at that share,
the targets' own rule, which must reproduce the row's figures.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PATH = 'research/studies/climate_targets/results/climate_targets.json'
SCHEMA = 'terluna.research.climate-targets/5'
# The one switch for the merge: False reads domain/climate's product from Git at the pinned commit, True reads PATH in
# the tree.
IN_TREE = False
PIN = dict(branch='domain/climate', commit='82b6ed4757da988f6f72ef06c8a5fb0db61bfd70',
           sha256='e5c7dcf3fb733259b051978bfeb99dc932532a7290f58acff4b0f20e481d2140',
           bind='Read in-tree once the branches meet; until then from Git at the commit.')
KEYS = ('hot_world', 'resort', 'lock')
NAMES = dict(hot_world='no further dimming', lock='295 K')
RESORT_K = 301.0
# The row's figures the interpolated climate must reproduce, and where the row holds them.
ROW_CHECKS = (('land_air_2m_median_c', 'comfort'), ('wet_bulb_2m_median_c', 'comfort'), ('sea_median_c', 'comfort'),
              ('land_rain_mm_day', 'rain'))


def load():
    """The climate targets' product and its binding as a producer records it: (data, binding), the binding being
    {'inputs': {PATH: sha256}} in the tree or {'pinned_inputs': {PATH: the pin}} from Git."""
    if IN_TREE:
        raw = (ROOT / PATH).read_bytes()
        binding = dict(inputs={PATH: hashlib.sha256(raw).hexdigest()})
    else:
        raw = subprocess.run(['git', '-C', str(ROOT), 'show', f"{PIN['commit']}:{PATH}"], check=True,
                             capture_output=True).stdout
        if hashlib.sha256(raw).hexdigest() != PIN['sha256']:
            raise ValueError(f"{PATH} at {PIN['commit'][:7]} does not match its pinned sha256")
        binding = dict(pinned_inputs={PATH: dict(PIN)})
    data = json.loads(raw)
    if data.get('schema') != SCHEMA:
        raise ValueError(f'{PATH}: expected {SCHEMA}, found {data.get("schema")}')
    return data, binding


def resort_climate(data, row):
    """The resort's 2 m air from its row, and its air at the lowest level and its sky ladder from the two runs'
    climates interpolated at the row's share of the way, checked against the row."""
    a, b, f = data['runs']['design'], data['runs']['test'], row['share_of_the_way']

    def at(key):
        return a[key] + f * (b[key] - a[key])

    for key, part in ROW_CHECKS:
        if abs(at(key) - row[part][key]) > 1e-9:
            raise ValueError(f'The resort row\'s {key} does not reproduce from the runs at its share of the way')
    ladder = [dict(height_km=x['height_km'], air_c=x['air_c'] + f * (y['air_c'] - x['air_c']),
                   pressure_atm=x['pressure_atm'] + f * (y['pressure_atm'] - x['pressure_atm']))
              for x, y in zip(a['ladder'], b['ladder'])]
    return dict(ground_k=row['ground_k'], land_air_2m_median_c=row['comfort']['land_air_2m_median_c'],
                wet_bulb_2m_median_c=row['comfort']['wet_bulb_2m_median_c'], land_air_median_c=at('land_air_median_c'),
                land_rain_mm_day=row['rain']['land_rain_mm_day'], ladder=ladder)


def read(data):
    """Each target, keyed as the rain products name them: where the product holds it, its ground temperature (K), its
    share of the way from the design case to the test, its land rain and that rain with the dimming's own cut
    (mm/day, [low, high]), and its climate."""
    by_name = {t['name']: t for t in data['targets']}
    out = {}
    for key in KEYS:
        if key == 'resort':
            row = next(r for r in data['resort_band']['rows'] if r['ground_k'] == RESORT_K)
            out[key] = dict(name_in_targets=f'{RESORT_K:g} K', held_in='resort_band.rows', ground_k=row['ground_k'],
                            share_of_the_way=row['share_of_the_way'], land_rain_mm_day=row['rain']['land_rain_mm_day'],
                            land_rain_with_the_dimmings_own_cut_mm_day=row['rain']['with_the_dimmings_own_cut_mm_day'],
                            climate=resort_climate(data, row))
        else:
            t = by_name[NAMES[key]]
            out[key] = dict(name_in_targets=t['name'], held_in='targets', ground_k=t['ground_k'],
                            share_of_the_way=t['interpolated_share_of_the_way'], land_rain_mm_day=t['land_rain_mm_day'],
                            land_rain_with_the_dimmings_own_cut_mm_day=t['land_rain_with_the_dimmings_own_cut_mm_day'],
                            climate=t['climate'])
    return out
