"""The inputs of the coupled study, read once: the programme, the engineering and recovery coefficients, the scenarios
of demand, the staged demand of the first pass, the financing rules and the reporting periods.

run.py and residents.py both read their inputs here, so that the files a result depends on are named in one place and
a change to how one of them reports its result leaves the other's product standing.
"""
from __future__ import annotations
from copy import deepcopy
from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path

from research.studies.enclosed_returns import markets
from . import demand as scenarios, programme as definition

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
INDUSTRY = ROOT / 'provisioning/vacuum_industry'
RETURNS = ROOT / 'research/studies/enclosed_returns'
# The source survey's estimate of the near-Earth population, with the water its pools contain, and its water by the
# velocity change of the tugs' leg from the depot.
POOL = ROOT / 'resources/source_access/results/population.json.gz'
REACH = ROOT / 'resources/source_access/results/reach.json.gz'
PERCENTILES = ('p05', 'p50', 'p95')
# The tugs freight.py chose for each scenario; run.py and residents.py read it, and freight.py does not.
FREIGHT = HERE / 'results/freight.json'
# The files the calculation reads or runs, beside those the staged order book and the financing name themselves.
READ = (HERE / 'programme.json', HERE / 'programme.py', HERE / 'inputs.py', HERE / 'scenarios.json', HERE / 'demand.py', HERE / 'mix.py',
        INDUSTRY / 'model.py', INDUSTRY / 'recovery.py', INDUSTRY / 'industry.py', INDUSTRY / 'engineering.json', INDUSTRY / 'recovery.json',
        RETURNS / 'markets.py', RETURNS / 'enclosed_inputs.json', RETURNS / 'finance.py', RETURNS / 'finance_inputs.json',
        RETURNS / 'revised.py', RETURNS / 'revised_inputs.json', RETURNS / 'bottlenecks.py', RETURNS.parent / 'complex_returns/policy.py', POOL, REACH)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# The survey's pools: the wider screen, the orbital window of the programme's intended sources (semimajor axis to
# 1.8 AU, eccentricity to 0.5, inclination to 30 degrees), and the broader near-Earth population without Don Quixote,
# whose use would need transport that reaches beyond that window.
POOLS = dict(wider='wider', broader='all_no_don_quixote')


@lru_cache(maxsize=None)
def pool():
    """The water contained in each of the source survey's two pools, at three percentiles, in kilograms.

    Each pool holds catalogued, inferred and undiscovered bodies, before recovery or any limit on access; the survey
    states them as conditional ranges, not reserves, and the wider screen lies within the broader population."""
    data = json.loads(gzip.decompress(POOL.read_bytes()))
    stated = data['scenarios']['central']['pools']
    return {name: {key: stated[key_in_survey]['total_quantiles_kg']['water_equivalent_kg'][key] for key in PERCENTILES}
            for name, key_in_survey in POOLS.items()}


def water_drawn(physical, loops, coefficients, horizons):
    """The water contained in the feed mined at the sources, to each horizon and in all, against each pool's water:
    the year the draw passes each percentile and the share of each pool's median drawn by the last year."""
    grade = coefficients['industrial']['source']['water_grade']
    if loops and loops.get('annual'):
        feed = [(r['year'], r['source_feed_kg']) for r in sorted(loops['annual'], key=lambda r: r['year'])]
    else:
        first = physical['annual'][0]['year']
        feed = [(first - 1, physical['seed']['source_feed_kg'])] + [(r['year'], r['physical']['source_feed_kg_year']) for r in physical['annual']]
    stock = pool()
    drawn, total = {}, 0.0
    past = {name: {key: None for key in levels} for name, levels in stock.items()}
    for year, kg in feed:
        total += kg * grade
        for name, levels in stock.items():
            for key, value in levels.items():
                if past[name][key] is None and total > value:
                    past[name][key] = year
        if year in horizons:
            drawn[str(year)] = total
    return dict(contained_water_drawn_kg=drawn, pool_contained_water_kg=stock, year_drawn_past=past,
                share_of_median_drawn_by_last_year={name: total / levels['p50'] for name, levels in stock.items()})


@lru_cache(maxsize=None)
def sources():
    """The programme, the engineering coefficients, the recovery measures, the demand inputs and the financing rules."""
    programme = definition.load()
    demand = markets.inputs()                                  # checks the digests of the inputs it inherits
    money = json.loads((RETURNS / 'finance_inputs.json').read_text())
    for name, pinned in money['inherited_dependencies'].items():
        if digest(RETURNS / name) != pinned:
            raise ValueError('Changed financing dependency: ' + name)
    periods = json.loads((RETURNS / 'revised_inputs.json').read_text())['periods']
    calendar, people, rules = programme['calendar'], programme['population'], demand[0]
    stated = dict(first_productive_year=calendar['first_productive_year'], completion_year=calendar['completion_year'],
                  last_year=calendar['last_year'], seed_population=people['opening'],
                  founding_population=people['founding_system'], cap_population=people['planned'])
    for key, value in stated.items():
        if rules[key] != value:
            raise ValueError('The demand inputs and the programme definition disagree on ' + key)
    if money['horizons'] != programme['horizons'] or money['rates'] != programme['discount_rates']:
        raise ValueError('The financing rules and the programme definition disagree on horizons or discount rates')
    if periods[0][0] != 0 or periods[-1][1] != calendar['last_year']:
        raise ValueError('The reporting periods do not span the programme')
    return dict(programme=programme, demand=demand, money=money, periods=periods, scenarios=scenarios.load(),
                engineering=json.loads((INDUSTRY / 'engineering.json').read_text()),
                recovery=json.loads((INDUSTRY / 'recovery.json').read_text()))


def dependencies(extra=()):
    """Each file read, with its digest: those listed here, those the staged inputs inherit and any others given."""
    s = sources()
    paths = {Path(p).resolve() for p in (*READ, *extra)}
    paths |= {(RETURNS / name).resolve() for name in s['demand'][0]['inherited_inputs']}
    paths |= {(RETURNS / name).resolve() for name in s['money']['inherited_dependencies']}
    return {str(p.relative_to(ROOT)): digest(p) for p in sorted(paths)}


@lru_cache(maxsize=None)
def reach():
    """The contained water of the source survey within each velocity change of the tugs' leg from the depot, at three
    percentiles, cumulative over the upper edges of its bins and, last, in all."""
    data = json.loads(gzip.decompress(REACH.read_bytes()))
    return dict(upper_edges_m_s=data['depot_leg_upper_edges_m_s'], within_kg=data['contained_water_within_kg'])


@lru_cache(maxsize=None)
def chosen_freight():
    """The tugs freight.py chose for each scenario, refused where any file it read has changed since."""
    product = json.loads(FREIGHT.read_text())
    stale = [path for path, value in product['dependencies'].items() if not (ROOT / path).exists() or digest(ROOT / path) != value]
    if stale:
        raise ValueError('The choice of tugs is out of date (run freight.py again): ' + ', '.join(stale))
    return product


# Freight to the sources. Under `founding` the tugs are the founding design's, chemical, on one leg of 2.5 km/s to
# every source, as in the staged studies. Under any other the tugs mine the nearest water of the source survey first,
# at the case's percentile of its water (`reach`, p50 unless stated), and fly: one named kind throughout (`chemical`,
# or an option of the engineering coefficients' transport.source_propulsion), `chosen`, the kinds freight.py chose for
# the scenario period by period, or `chosen_with_target`, its choice with the water thruster at its development target
# among the options. The first pass keeps the founding freight; a scenario's cases take `chosen` unless they state.
# Where freight.py found a period whose water no kind of tug reaches, the scenario's orders cannot be filled with the
# reach, and a case that would take that choice flies the founding freight; its record names those periods.
def stated_freight(case):
    return case.get('freight', 'founding' if case['scenario'] == 'first_pass' else 'chosen')


def unreached(case):
    """The periods whose water no kind of tug reaches, for a case that takes freight.py's choice; none otherwise."""
    stated = stated_freight(case)
    if stated not in ('chosen', 'chosen_with_target'):
        return []
    return chosen_freight()['scenarios'][case['scenario']][stated]['periods_no_kind_reaches']


def freight_of(case):
    return 'founding' if unreached(case) else stated_freight(case)


def engineering(case):
    """The engineering coefficients for one case: the model's own, with the changes the case states and its freight.

    A change names a coefficient by its path, as in `industrial.transport.earth_leg`, and gives its value. The
    coefficient must exist, so that a misspelt path fails where it is read.
    """
    data = sources()['engineering']
    stated = sources()['scenarios']['scenarios'].get(case['scenario'], {}).get('engineering', {})
    changes = {**{path: entry['value'] for path, entry in stated.items()}, **case.get('engineering', {})}
    freight = freight_of(case)
    if not changes and freight == 'founding':
        return data
    data = deepcopy(data)
    for path, value in changes.items():
        *parents, leaf = path.split('.')
        node = data
        for key in parents:
            node = node[key]
        if leaf not in node:
            raise ValueError('The engineering coefficients hold no ' + path)
        node[leaf] = value
    if freight != 'founding':
        curve = reach()
        percentile = case.get('reach', 'p50')
        data['source_reach'] = dict(upper_edges_m_s=curve['upper_edges_m_s'], within_kg=curve['within_kg'][percentile], percentile=percentile)
        if freight in ('chosen', 'chosen_with_target'):
            kind = chosen_freight()['scenarios'][case['scenario']][freight]['schedule']
        elif freight == 'chemical' or freight in data['transport']['source_propulsion']['options']:
            kind = freight
        else:
            raise ValueError('Unknown freight: ' + str(freight))
        data['industrial']['transport']['source_propulsion'] = kind
    return data
