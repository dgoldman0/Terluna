"""The programme definition of the vacuum buildout's economy: its calendar, its money and its population ramps.

Every part of the coupled study reads the schedule and the population from programme.json through this module, so
that a year, a horizon or a number of residents is stated once.
"""
from __future__ import annotations
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
PATH = HERE / 'programme.json'


def load(path=PATH):
    programme = json.loads(Path(path).read_text())
    calendar, people = programme['calendar'], programme['population']
    if not calendar['first_productive_year'] < calendar['completion_year'] <= calendar['last_year']:
        raise ValueError('The calendar needs first production before completion, and completion by the last year')
    if not 0 < people['opening'] < people['founding_system'] < people['planned']:
        raise ValueError('The opening, founding and planned populations must rise in that order')
    for name, ramp in programme['ramps'].items():
        if ramp['kind'] == 'two_rates' and not calendar['first_productive_year'] < ramp['founding_system_year'] < calendar['completion_year']:
            raise ValueError('The founding system must stand between first production and completion: ' + name)
    return programme


def population(programme, ramp, year):
    """Residents housed in a year under the named ramp; the planned population from the completion year on."""
    calendar, people, spec = programme['calendar'], programme['population'], programme['ramps'][ramp]
    first, done = calendar['first_productive_year'], calendar['completion_year']
    opening, founding, planned = people['opening'], people['founding_system'], people['planned']
    kind = spec['kind']
    if kind == 'constant_rate':
        u = min(1., max(0., (year - first) / (done - first)))
        return opening * (planned / opening)**u
    if kind == 'two_rates':
        stands = spec['founding_system_year']
        if year <= stands:
            return opening * (founding / opening)**(max(0., year - first) / (stands - first))
        return founding * (planned / founding)**(min(1., (year - stands) / (done - stands)))
    if kind == 'smooth_logarithmic':
        u = max(0., min(1., (year - spec['start_year']) / (done - spec['start_year'])))
        return opening * math.exp(math.log(planned / opening) * (3 * u * u - 2 * u * u * u))
    raise ValueError('Unknown kind of ramp: ' + kind)


def schedule(programme, ramp):
    """Each year from first production to the last year, with its residents and the residents added in it."""
    calendar = programme['calendar']
    years = range(calendar['first_productive_year'], calendar['last_year'] + 1)
    rows = [dict(year=y, population=population(programme, ramp, y)) for y in years]
    for before, row in zip([None, *rows], rows):
        row['added'] = 0. if before is None else row['population'] - before['population']
    return rows


def growth_rate(programme, ramp, year):
    """The fractional growth of the population over the year that ends in the given year."""
    before = population(programme, ramp, year - 1)
    return population(programme, ramp, year) / before - 1
