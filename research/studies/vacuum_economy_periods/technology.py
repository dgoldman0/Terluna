"""The technology of the vacuum period, 2041-2241: each quantity of the stack as a path in a slow, a central and a fast
case, the static case that holds 2041, and the technologies that arrive later, read from period_technology.json.

A case can also move one group of the stack to another case while the rest keep theirs: 'central|power=slow' takes
the slow path for the power group's quantities and arrivals and the central path for everything else, which shows
how much each group decides on its own."""
from __future__ import annotations
from functools import lru_cache
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = ('static', 'slow', 'central', 'fast')

# The stack's quantities and arrivals in six groups, each a field of work that develops on its own.
GROUPS = {
    'launch': ('launch_USD2026_kg', 'landed_multiple', 'mass_driver'),
    'tugs_and_sources': ('tug_life_years', 'water_extraction_kWh_kg', 'water_recovered_from_rock', 'nuclear_electric'),
    'power': ('array_W_kg', 'array_local_share', 'storage_Wh_kg', 'storage_round_trip', 'fission_kg_kWe', 'fission_life_years',
              'fission_USD2024_kWe', 'beam_station_W_kg', 'fusion', 'polar_grid', 'equatorial_belt'),
    'industry': ('metal_kWh_kg', 'equipment_import_share', 'output_per_job', 'equipment_price_multiple', 'equipment_life_multiple',
                 'industry_doubling_years', 'services_output_multiple', 'built_mass_multiple', 'chips_on_moon'),
    'food_and_life_support': ('led_umol_J', 'crop_kcal_per_mol', 'microbial_diet_share', 'microbial_efficiency', 'basic_kW_resident',
                              'water_recovery'),
    'computing': ('computation_per_W_multiple',),
}
GROUP_OF = {name: group for group, names in GROUPS.items() for name in names}


def base(case):
    """The case a technology case starts from: 'central' for 'central|power=slow'."""
    return case.split('|')[0]


@lru_cache(maxsize=None)
def moves(case):
    """The groups a case moves, as {group: case}; empty for a plain case."""
    out = {}
    for move in case.split('|')[1:]:
        group, other = move.split('=')
        if group not in GROUPS or other not in CASES:
            raise ValueError('Unknown group or case in %r' % case)
        out[group] = other
    return out


def case_for(name, case):
    """The case one quantity or arrival follows: its group's case where the case moves that group, else the base case."""
    if '|' not in case:
        return case
    return moves(case).get(GROUP_OF[name], base(case))


@lru_cache(maxsize=None)
def roadmap():
    return json.loads((HERE / 'period_technology.json').read_text())


@lru_cache(maxsize=None)
def value(name, case, year):
    """One quantity in a calendar year: the 2041 value in the static case, otherwise interpolated between the stated
    years, geometrically where the path says so, and held after the last."""
    case = case_for(name, case)
    t = roadmap()
    q = t['quantities'][name]
    years = t['years']
    values = q['slow' if case == 'static' else case]
    if case == 'static' or year <= years[0]:
        return float(values[0])
    if year >= years[-1]:
        return float(values[-1])
    i = next(k for k in range(1, len(years)) if year <= years[k])
    f = (year - years[i - 1]) / (years[i] - years[i - 1])
    a, b = values[i - 1], values[i]
    if q.get('interpolate') == 'geometric' and a > 0 and b > 0:
        return float(a * (b / a) ** f)
    return float(a + (b - a) * f)


def arrival(name, case):
    """The first year a technology is available in a case, or None where it does not arrive (always in the static case)."""
    case = case_for(name, case)
    if case == 'static':
        return None
    return roadmap()['arrivals'][name]['year'][case]


def arrived(name, case, year):
    first = arrival(name, case)
    return first is not None and year >= first


def parameter(name, key, case):
    """A parameter an arrival brings: a single value, or a value by case."""
    v = roadmap()['arrivals'][name][key]
    return v[case_for(name, case)] if isinstance(v, dict) else v


class Year:
    """The technology of one calendar year in one case, read by attribute: `Year('central', 2100).array_W_kg`."""

    def __init__(self, case, year):
        self.case, self.year = case, year

    def __getattr__(self, name):
        if name.startswith('_'):
            raise AttributeError(name)
        return value(name, self.case, self.year)

    def arrived(self, name):
        return arrived(name, self.case, self.year)
