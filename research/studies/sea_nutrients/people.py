"""People's phosphorus on the design surface: where what they lose to the water ends, lake by lake and sea by sea.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.sea_nutrients.people   # about two minutes

The phosphorus study sized what people lose, 44-180 g a person a year in the food design, and where it can be won
back, at the lake deltas. The seas' budget (run.py) found that people's phosphorus is the one input large enough to
break it. This module places that loss on the map. The atlas round's land-use screen marks the land that the chosen
distribution of 8 billion people fills with its homes, farms and pasture; the loss is spread over that land by area
and carried down the design routing as loads.py carries the land's own, in solution and on particles, to the lakes
that keep it and the seas it reaches.

It answers four things: what share of people's loss reaches each sea, and what share of it could reach one without
crossing a lake; which lakes keep it, and how few of them keep most of it, since those are where it can be won back;
what it does to the lakes' own phosphorus, lake by lake, against the lines the lake studies draw at 10 and 35 mg/m3;
and how much room the seas have left. The loss is carried per gram a person a year, so that any split between the
dissolved part and the part bound to soil can be read from it. Writes results/people.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from research.studies.sea_nutrients import loads as LO
from resources import sea_phosphorus_waters as SPW

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'people.json'
SCHEMA = 'terluna.research.sea-nutrients-people/1'
LAND_USE = ROOT / 'research' / 'studies' / 'atlas_round' / 'results' / 'land_use.json'
PHOSPHORUS = ROOT / 'resources' / 'results' / 'sea_phosphorus.json'
BUDGET = HERE / 'results' / 'sea_nutrients.json'
PEOPLE = 8.0e9                              # the chosen distribution, 4 billion on the ground and 4 aloft, fed from the ground
VARIANTS = dict(A='30% of the land wild', B='half of the land wild')
PEOPLES_LAND = 3                            # the land-use grid's kind for open land the distribution fills
DISSOLVED_G = (10.0, 30.0, 90.0)            # g a person a year reaching water in solution: cases up to the glass houses' bleed
LINES_MG_M3 = dict(nutrient_poor=10.0, rich=35.0)       # the lake studies' classes by total phosphorus
YEAR_S = LO.YEAR_S


def sha256(path) -> str:
    digest = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def by_lake(c, passing, factor):
    """What reaches each lake's exit and what the lake keeps of it, from what passes the cells where water leaves it."""
    leaving = c['leaving']
    b = c['body'][leaving]
    n = c['volume'].size
    arrives = np.bincount(b, weights=passing[leaving], minlength=n)
    kept = np.bincount(b, weights=passing[leaving] * (1.0 - factor[leaving]), minlength=n)
    return arrives, kept


def how_few(kept, area, shares=(0.5, 0.9)):
    """How many lakes, taken from the one that keeps most, keep each share of all that the lakes keep, with their area."""
    order = np.argsort(kept)[::-1]
    cum = np.cumsum(kept[order]) / max(kept.sum(), 1e-300)
    out = {}
    for s in shares:
        k = int(np.searchsorted(cum, s)) + 1
        out[f'{int(100 * s)}%'] = dict(lakes=k, area_km2=float(area[order[:k]].sum() / 1e6))
    return out


def lake_phosphorus(arrives_kg_yr, c):
    """Each spilling lake's phosphorus (mg/m3) by Vollenweider's relation: what arrives over the water leaving, times
    the share a lake of its residence passes on."""
    spills = c['spills']
    survive = np.where(spills, SPW.lake_survival_dissolved(np.where(spills, c['tau'], 0.0)), 0.0)
    return np.where(spills, 1e6 * arrives_kg_yr * survive / np.maximum(c['outflow'] * YEAR_S, 1e-30), np.nan)


def over_the_lines(tp, area):
    ok = np.isfinite(tp)
    total = area[ok].sum()
    return {name: dict(lakes=int((tp[ok] > line).sum()), area_km2=float(area[ok & (tp > line)].sum() / 1e6),
                       share_of_lake_area=float(area[ok & (tp > line)].sum() / total))
            for name, line in LINES_MG_M3.items()}


def one_variant(c, kind, natural_arrives, land, budget, log):
    receiver, fronts, factors, into, mouth, names = (c[k] for k in ('receiver', 'fronts', 'factors', 'into', 'mouth', 'names'))
    cell_area = np.broadcast_to(c['row_area'][:, None], (c['rows'], c['cols'])).ravel()
    on = (kind == PEOPLES_LAND) & c['donor']
    source = np.where(on, cell_area, 0.0)
    land_m2 = float(source.sum())
    source = source / land_m2
    n_seas = len(names)
    out = dict(peoples_land_km2=land_m2 / 1e6)
    passing = {}
    for key in ('dissolved', 'particles_central', 'runoff_that_crossed_no_lake', 'runoff'):
        p = LO.carry(receiver, fronts, source, factors[key])
        passing[key] = p
        out[key] = dict(reaches_each_sea=dict(zip(names, (np.bincount(into[mouth], weights=p[mouth] * factors[key][mouth],
                                                                      minlength=n_seas)).tolist())))
        out[key]['reaches_the_seas'] = float(sum(out[key]['reaches_each_sea'].values()))
    ends_closed = float(passing['runoff'][c['keep']].sum())
    out['ends_in_lakes_that_pass_nothing_on'] = ends_closed
    arrives, kept = by_lake(c, passing['dissolved'], factors['dissolved'])
    _, kept_particles = by_lake(c, passing['particles_central'], factors['particles_central'])
    out['dissolved']['kept_by_spilling_lakes'] = float(kept.sum())
    out['dissolved']['kept_where_water_ends'] = float(passing['dissolved'][c['keep']].sum())
    out['particles_central']['kept_by_spilling_lakes'] = float(kept_particles.sum())
    out['particles_central']['kept_where_water_ends'] = float(passing['particles_central'][c['keep']].sum())
    out['the_lakes_that_keep_it'] = dict(dissolved=how_few(kept, c['area']), particles=how_few(kept_particles, c['area']))

    # the lakes' own phosphorus, with the land's and with people's
    lakes = {}
    for end, kg in zip(('low', 'high'), land['dissolved_kg_yr']):
        natural = lake_phosphorus(natural_arrives * kg, c)
        entry = dict(the_land_alone=over_the_lines(natural, c['area']))
        for grams in DISSOLVED_G:
            with_people = lake_phosphorus(natural_arrives * kg + arrives * grams * PEOPLE / 1e3, c)
            entry[f'with {grams:g} g a person'] = over_the_lines(with_people, c['area'])
        lakes[f'the land\'s loss at its {end} end'] = entry
    out['lakes_over_the_lines'] = lakes

    # the seas: what each gram a person a year adds, and the room
    seas = {}
    for k, name in enumerate(names):
        sea = budget['seas'][name]
        share = out['dissolved']['reaches_each_sea'][name]
        per_gram = share * PEOPLE * 1e3 / (sea['area_km2'] * 1e6)                 # mg per m2 of sea a year for each g a person
        room = sea['room_beside_rivers_and_dust']['1000 years']['mg_m2_yr'][0]
        seas[name] = dict(share_of_peoples_dissolved_loss=share,
                          share_that_crossed_no_lake=out['runoff_that_crossed_no_lake']['reaches_each_sea'][name],
                          mg_m2_yr_for_each_g_a_person=per_gram,
                          room_for_a_thousand_years_mg_m2_yr=room,
                          g_a_person_that_fills_the_room=room / per_gram if per_gram > 0 else None,
                          supply_with_people_mg_m2_yr={f'{g:g} g': [sea['supply_mg_m2_yr']['rivers_and_dust'][0] + per_gram * g,
                                                                    sea['supply_mg_m2_yr']['rivers_and_dust'][1] + per_gram * g]
                                                       for g in DISSOLVED_G})
    out['seas'] = seas
    log(f'people\'s land {land_m2 / 1e12:.2f} million km2: {out["dissolved"]["reaches_the_seas"]:.4f} of the dissolved loss '
        f'reaches the seas, {out["runoff_that_crossed_no_lake"]["reaches_the_seas"]:.4f} crossing no lake')
    return out


def build(log=print):
    started = time.time()
    c = LO.prepare(log)
    use = json.loads(LAND_USE.read_text())
    grid = Path(use['product']['path'])
    if sha256(grid) != use['product']['sha256']:
        raise SystemExit(f'{grid} does not match the grid {LAND_USE} records')
    land = json.loads(PHOSPHORUS.read_text())['land_loss']
    people = json.loads(PHOSPHORUS.read_text())['people']
    budget = json.loads(BUDGET.read_text())
    natural = LO.carry(c['receiver'], c['fronts'], c['runoff'] / c['total'], c['factors']['dissolved'])
    natural_arrives, _ = by_lake(c, natural, c['factors']['dissolved'])
    variants = {}
    with np.load(grid, allow_pickle=False) as z:
        for key in VARIANTS:
            name = 'kind_' + key
            variants[key] = dict(what=VARIANTS[key], **one_variant(c, z[name].ravel(), natural_arrives, land, budget, log))
    spills = c['spills']
    return dict(
        schema=SCHEMA,
        producer=dict(module='research.studies.sea_nutrients.people', code_sha256=sha256(Path(__file__))),
        inputs={str(p.relative_to(ROOT)): sha256(p) for p in (LAND_USE, PHOSPHORUS, BUDGET, ROOT / LO.ROUTING_RECORD)},
        ignored_inputs={str(grid): use['product']['sha256']},
        evidence=('Tracers carried down the design routing from the land the land-use screen gives to people, spread '
                  'by area, with Earth\'s lake relations for what each lake keeps. The screen places kinds of land and '
                  'no town; the loss a person is the phosphorus study\'s. Not validated.'),
        reading_rule=('variants: for each wild share, the land people fill and, for a loss spread over it, the share '
                      'that reaches each sea in solution and on particles, the share that would reach a sea crossing '
                      'no lake, the share kept by lakes, how few lakes keep half and nine tenths of it, the lakes whose '
                      'phosphorus passes 10 and 35 mg/m3 with the land\'s loss alone and with each dissolved loss a '
                      'person, and for each sea what a gram a person a year adds and the grams that would fill its '
                      f'room for a thousand years. People number {PEOPLE:g}.'),
        setting=dict(people=PEOPLE, dissolved_g_a_person=list(DISSOLVED_G), lines_mg_m3=LINES_MG_M3,
                     design_loss_g_a_person=[1e3 * x for x in people['ecology_design_target']['lost_kg_per_person']],
                     eaten_kg_a_person=people['eaten_kg_per_person'], lakes=int(spills.sum()),
                     lake_area_km2=float(c['area'][spills].sum() / 1e6)),
        variants=variants, seconds=round(time.time() - started))


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    out = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(f'wrote {OUT} in {out["seconds"]} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
