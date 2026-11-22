"""The land's phosphorus through the design surface's lakes to each of its seas.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.sea_nutrients.loads   # 2-4 minutes

The phosphorus study of 11 October (resources/sea_phosphorus.md) followed the land's phosphorus through the lakes of
the atlas's drainage and found that the lakes keep 96% of it. This module repeats that screen on the design surface
as it stands since 18 November, whose lakes, seas and flows differ: it reads the design routing
(research/studies/atlas_round/design_routing.py's product on the research drive) and the redrawn floor, and carries
tracers down the routing's receivers, as resources/sea_phosphorus_waters.py does:

- the land's runoff is what each dry cell adds to the flow through it (its discharge less its donors'); phosphorus
  rides it at one concentration;
- a lake is the cells that stand at one level; its residence is its volume over the flow leaving it, and where water
  leaves it the lake passes on 1 / (1 + sqrt(tau)) of the dissolved phosphorus (Vollenweider 1976) and
  1 - exp(-beta I/C) of the particles (Siyam's form of Brune's curves), I/C slowed sixfold by lunar gravity in the
  central case; lakes that pass nothing on keep everything;
- the three upper seas spill over their rims. Where their water enters the routing again it brings the share of
  their dissolved phosphorus that a water of their residence passes on, shared among the seas it runs to by flow,
  with no further loss on the way (the design routing's record gives the flows); the main part of the near-side sea
  sends its southern part water through the neck on the same rule. Their particles stay with them.

Writes results/loads.json: for each sea the shares of the land's dissolved and particulate loss that reach it, as
shares of the whole Moon's loss, so that the phosphorus study's range of that loss gives the kilograms.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from research.studies.atlas_round import design_routing as DR
from research.studies.atlas_round import redrawn_drainage as RD
from resources import sea_phosphorus_waters as SPW
from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'loads.json'
SCHEMA = 'terluna.research.sea-nutrients-loads/1'
YEAR_S = 31557600.0
OVERFLOW_M3S = 500.0            # a dry cell that adds more than this to its flow is where an upper sea's water enters
ROUTING_RECORD = 'research/studies/atlas_round/results/design_routing.json'
PASSING_ORDER = ('moscoviense', 'orientale', 'smythii_marginis', 'near_side_main')   # each after the seas that feed it
PARTICLE_CASES = dict(low=('upper', 1.0), central=('median', MOON_SURFACE_GRAVITY / STANDARD_GRAVITY),
                      high=('lower', MOON_SURFACE_GRAVITY / STANDARD_GRAVITY))


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def kahn_fronts(receiver, donor):
    """The donors in an order in which every cell comes after all the cells that drain to it, as a list of fronts."""
    n = receiver.size
    indeg = np.bincount(receiver[donor], minlength=n)
    front = np.flatnonzero(donor & (indeg == 0))
    fronts = []
    while front.size:
        fronts.append(front)
        dest = receiver[front]
        indeg -= np.bincount(dest, minlength=n) if front.size > n // 50 else 0
        if front.size <= n // 50:
            np.subtract.at(indeg, dest, 1)
        cand = np.unique(dest)
        front = cand[(indeg[cand] == 0) & donor[cand]]
    return fronts


def carry(receiver, fronts, source, factor):
    """What passes each cell: its own source and all that its donors pass on, each cell passing on its factor's share."""
    p = source.astype(float).copy()
    for front in fronts:
        np.add.at(p, receiver[front], p[front] * factor[front])
    return p


def lake_bodies(rt, h, row_area, cols):
    """The lakes off the seas, each the cells at one level: the body of every cell (-1 off the lakes), each body's
    volume (m3) and area (m2), the flow leaving it (m3/s) and the cells where its water leaves."""
    cell, level = rt['lake_cell'].astype(np.int64), rt['lake_level_m']
    receiver, q = rt['receiver'].astype(np.int64), rt['discharge_m3s']
    values, inverse = np.unique(level, return_inverse=True)
    body = np.full(receiver.size, -1, np.int64)
    body[cell] = inverse
    a = row_area[cell // cols]
    depth = np.maximum(level.astype(float) - h.ravel()[cell], 0.0)
    volume = np.bincount(inverse, weights=a * depth, minlength=values.size)
    area = np.bincount(inverse, weights=a, minlength=values.size)
    leaving = cell[body[receiver[cell]] != inverse]
    outflow = np.bincount(body[leaving], weights=q[leaving], minlength=values.size)
    return body, volume, area, outflow, leaving


def build(log=print):
    started = time.time()
    rt = DR.load(log=log)
    meta = rt['meta']
    surf = RD.load_surface(log=lambda *a: None)
    h, row_area = surf['h'], np.asarray(surf['row_area'], float)
    rows, cols = h.shape
    names, levels = meta['seas'], meta['levels_m']
    receiver, q, sea, outlet = rt['receiver'].astype(np.int64), rt['discharge_m3s'], rt['sea'], rt['outlet']
    n = receiver.size
    donor = ~outlet
    inflow = np.bincount(receiver[donor], weights=q[donor], minlength=n)
    local = q - inflow
    dry = donor & ~rt['routing_lake']
    overflow = np.flatnonzero(dry & (local > OVERFLOW_M3S))
    runoff = np.where(dry, np.maximum(local, 0.0), 0.0)
    runoff[overflow] = 0.0
    total = float(runoff.sum())
    log(f'the land adds {total:.0f} m3/s of runoff; {overflow.size} cells take in an upper sea\'s water, {time.time() - started:.0f} s')

    body, volume, area, outflow, leaving = lake_bodies(rt, h, row_area, cols)
    spills = outflow > 0
    tau = np.where(spills, volume / np.maximum(outflow * YEAR_S, 1e-30), np.inf)
    gravity = MOON_SURFACE_GRAVITY / STANDARD_GRAVITY
    keep = rt['keeps']
    b = body[leaving]
    factors = {}
    one = np.ones(n, np.float32)
    one[keep] = 0.0
    dissolved = one.copy()
    dissolved[leaving] = np.where(spills[b], SPW.lake_survival_dissolved(np.where(spills[b], tau[b], 0.0)), 0.0)
    factors['dissolved'] = dissolved
    for case, (curve, ratio) in PARTICLE_CASES.items():
        f = one.copy()
        f[leaving] = np.where(spills[b], SPW.lake_survival_particulate(np.where(spills[b], tau[b], 1.0),
                                                                        SPW.BRUNE_BETA[curve], ratio), 0.0)
        factors[f'particles_{case}'] = f
    no_lake = one.copy()
    no_lake[leaving] = 0.0
    factors['runoff_that_crossed_no_lake'] = no_lake
    factors['runoff'] = one
    fronts = kahn_fronts(receiver, donor)
    log(f'{len(fronts)} fronts, {time.time() - started:.0f} s')
    into = np.where(donor, sea[receiver], -1)
    mouth = np.flatnonzero(into >= 0)
    arrive = {}
    for key, f in factors.items():
        p = carry(receiver, fronts, runoff, f)
        arrive[key] = np.bincount(into[mouth], weights=p[mouth] * f[mouth], minlength=len(names)) / total
        log(f'{key}: {arrive[key].sum():.4f} of the land\'s reaches the seas, {time.time() - started:.0f} s')
    ended = float(carry(receiver, fronts, runoff, factors['runoff'])[keep].sum()) / total

    # the seas that pass water on: the upper seas over their rims, and the near-side sea's main part through its neck
    area_m2 = np.broadcast_to(row_area[:, None], (rows, cols)).ravel()
    sea_volume = {name: float((area_m2[sea == k] * (levels[name] - h.ravel()[sea == k])).sum()) for k, name in enumerate(names)}
    summary = json.loads((ROOT / ROUTING_RECORD).read_text())['summary']
    weirs = summary['weirs']
    sends = {'moscoviense': weirs['released_m3s'],
             **{name: {'near_side_main': flow} for name, flow in weirs['full_seas_pass_on_m3s'].items()}}
    if summary.get('neck_spill_m3s'):
        sends['near_side_main'] = {'near_side_south': summary['neck_spill_m3s']}
    entering = float(local[overflow].sum())
    transfers = {}
    for name in PASSING_ORDER:
        to = {k: float(v) for k, v in sends.get(name, {}).items()}
        flow = sum(to.values())
        if flow:
            t = sea_volume[name] / (flow * YEAR_S)
            transfers[name] = dict(outflow_m3s=flow, residence_yr=t, passes_on=float(SPW.lake_survival_dissolved(t)), to_m3s=to)
    seas = {}
    for k, name in enumerate(names):
        seas[name] = dict(
            share_of_the_lands_dissolved_loss=float(arrive['dissolved'][k]),
            share_of_the_lands_particles={case: float(arrive[f'particles_{case}'][k]) for case in PARTICLE_CASES},
            share_of_the_lands_runoff=float(arrive['runoff'][k]),
            of_it_crossed_no_lake=float(arrive['runoff_that_crossed_no_lake'][k] / max(arrive['runoff'][k], 1e-30)),
            volume_km3=sea_volume[name] / 1e9)
    received = {name: 0.0 for name in names}
    for name in list(transfers):                      # an upper sea's own share first, then what reaches it from another
        t = transfers[name]
        t['dissolved_share_passed_on'] = (seas[name]['share_of_the_lands_dissolved_loss'] + received[name]) * t['passes_on']
        for dest, flow in t['to_m3s'].items():
            received[dest] += t['dissolved_share_passed_on'] * flow / t['outflow_m3s']
    for name in names:
        seas[name]['dissolved_share_from_other_seas'] = received[name]
    census = dict(lakes=int(volume.size), spilling=int(spills.sum()), lake_area_km2=float(area.sum() / 1e6),
                  lake_volume_km3=float(volume.sum() / 1e9), land_runoff_m3s=total,
                  runoff_ending_in_closed_lakes=ended,
                  residence_yr_by_outflow={f'p{p}': float(v) for p, v in zip((10, 50, 90), weighted_percentiles(tau[spills], outflow[spills], (10, 50, 90)))})
    return dict(
        schema=SCHEMA,
        producer=dict(module='research.studies.sea_nutrients.loads', code_sha256=sha256(Path(__file__))),
        inputs={ROUTING_RECORD: sha256(ROOT / ROUTING_RECORD)},
        evidence=('Tracers carried down the design routing\'s receivers, with Earth\'s lake relations for the phosphorus '
                  'each lake keeps (Vollenweider 1976; Siyam 2000 on Brune 1953) at a steady yearly flow and lake beds '
                  'that do not fill. Phosphorus follows the land\'s runoff at one concentration. Not validated.'),
        reading_rule=('seas: for each sea of the design surface, the shares of the whole land\'s dissolved loss and of '
                      'its particles that reach it through its own rivers, the share of the land\'s runoff that does, '
                      'how much of that crossed no lake, and the dissolved share that other seas pass on to it. '
                      'transfers: the seas that pass water on, with their residence and the share a water of that '
                      'residence passes on. The particle cases are Brune\'s upper curve at Earth\'s settling (low), '
                      'the median curve with settling slowed by lunar gravity (central) and the lower curve so slowed '
                      '(high). Multiply the shares by the land\'s loss to get kilograms.'),
        census=census, water_entering_from_upper_seas_m3s=entering, transfers=transfers, seas=seas, seconds=round(time.time() - started))


def weighted_percentiles(values, weights, percents):
    order = np.argsort(values)
    v, w = values[order], weights[order]
    return [float(np.interp(p / 100.0, np.cumsum(w) / w.sum(), v)) for p in percents]


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    out = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(f'wrote {OUT} in {out["seconds"]} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
