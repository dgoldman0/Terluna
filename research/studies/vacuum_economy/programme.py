"""The programme definition of the vacuum buildout's economy: its calendar, its money and its population ramps.

Every part of the coupled study reads the schedule and the population from programme.json through this module, so
that a year, a horizon or a number of residents is stated once.

A ramp states the residents of each year, or, for a ramp of the kind `solved`, the rule by which they follow from
the system's own production and financing. The residents of a solved ramp are found by residents.py and given to the
definition with `solved`, which returns a copy of the programme that carries them.
"""
from __future__ import annotations
from copy import deepcopy
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
PATH = HERE / 'programme.json'
PLANNED_KINDS = ('constant_rate', 'two_rates', 'smooth_logarithmic')


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
        if ramp['kind'] == 'solved' and 'residents' in ramp:
            raise ValueError('The residents of a solved ramp come from its solution, not from the definition: ' + name)
    return programme


def planned(programme):
    """The ramps whose residents the programme states, beside those it solves."""
    return sorted(name for name, ramp in programme['ramps'].items() if ramp['kind'] in PLANNED_KINDS)


def solved(programme, ramp, residents):
    """A copy of the programme in which a solved ramp carries its residents: [[year, residents], ...] for every year
    from first production to the last year, opening on the opening population and never falling."""
    spec = programme['ramps'][ramp]
    if spec['kind'] != 'solved':
        raise ValueError('Only a solved ramp takes its residents from a solution: ' + ramp)
    calendar = programme['calendar']
    years = [int(y) for y, _ in residents]
    if years != list(range(calendar['first_productive_year'], calendar['last_year'] + 1)):
        raise ValueError('A solved ramp lists every year from first production to the last year')
    values = [float(v) for _, v in residents]
    if not math.isclose(values[0], programme['population']['opening'], rel_tol=1e-12):
        raise ValueError('A solved ramp opens on the opening population')
    if any(b < a * (1 - 1e-12) for a, b in zip(values, values[1:])):
        raise ValueError('The residents of a solved ramp never fall')
    out = deepcopy(programme)
    out['ramps'][ramp] = dict(spec, residents=[[y, v] for y, v in zip(years, values)])
    return out


def population(programme, ramp, year):
    """Residents housed in a year under the named ramp; the planned population from the completion year on, or for a
    solved ramp its last listed year's residents after that year."""
    calendar, people, spec = programme['calendar'], programme['population'], programme['ramps'][ramp]
    first, done = calendar['first_productive_year'], calendar['completion_year']
    opening, founding, planned = people['opening'], people['founding_system'], people['planned']
    kind = spec['kind']
    if kind == 'solved':
        if 'residents' not in spec:
            raise ValueError('The residents of a solved ramp come from residents.py: ' + ramp)
        points = spec['residents']
        index = min(max(0, int(year) - points[0][0]), len(points) - 1)
        if points[index][0] != year and not (year > points[-1][0] or year < points[0][0]):
            raise ValueError('A solved ramp gives residents only for whole years')
        return points[index][1]
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
