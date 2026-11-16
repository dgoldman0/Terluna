"""The drainage network's water at chosen years, for the globe page's years control.

    /media/projectspace/terluna-research/ocean/run.sh -m research.studies.drainage_network.stops

The network study follows the six large seas together from the atlas's state for 10,000 years with the routing's
budget closed (network.py, the case the study takes as its better guide). This module records that history at
STOP_YEARS and draws each stop's water on the 16 px/deg grid:

- each sea stands at its level that year over every cell whose joining level lies below it, so a rising sea floods
  the basins its water reaches and a falling one leaves its higher floor dry;
- a sea that is full sends its overflow down its spill path, filling each hollow on the way to the level it spills
  at; the path's other cells carry the overflow as a river with the surplus it passes on;
- every other water stays as the atlas draws it: the hollows the atlas holds at the common level and the rain-fed
  lakes, whose own course the study does not follow, and the rivers routed with every sea at the common level.

Each cell carries the bodies of water that cover it at some stop (at most two), which of them is its surface at
each stop, and whether one of them covers it at each stop; each body has its level at each stop. Year 0 is the
atlas's own water. After it a full sea stands just under its rim, and a sea within half a metre below a pass is
drawn joined over it: from year 500 the history's five-year steps leave the near-side ocean hovering at the pass,
120 m below the common level, beyond which lies a fifth of its water. The history gives each sea one level, so a
basin beyond a pass is drawn at the sea's level once the sea passes over it. The maps go to
research/runs/drainage_network/globe_stops.npz on today's grid (0 E at the left edge) and the levels, the events,
the bodies and each stop's water share to results/globe_stops.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

from research.studies.drainage_network import network as N

HERE = Path(__file__).resolve().parent
ROOT = N.ROOT
OUT = HERE / 'results' / 'globe_stops.json'
MAPS = N.RUNS / 'globe_stops.npz'
SCHEMA = 'terluna.research.drainage-network-globe-stops/1'
STOP_YEARS = (0, 250, 500, 750, 1000, 2000, 5000, 10000)   # at most eight: one bit each
FIRST_HOLLOW_ID = 7                     # ids 1-6 are the six seas in the study's order
JOIN_WITHIN_M = 0.5                     # a sea within this of a pass above it is drawn joined over it: the history's
                                        # 5-year steps and 0.1 m records leave a sea hovering at such a pass
CODE = ('research/studies/drainage_network/stops.py', *N.CODE)


def code_hash():
    h = hashlib.sha256()
    for path in CODE:
        h.update((ROOT / path).read_bytes())
    return h.hexdigest()


def runs_of(flags, values=None):
    """Start and end (exclusive) of each run of True in a 1-D boolean array; with `values`, a run also ends where
    the value changes, so that two hollows that meet on a path stay apart."""
    flags = np.asarray(flags, bool)
    out = []
    start = None
    for i, f in enumerate(flags):
        if start is not None and (not f or (values is not None and values[i] != values[i - 1])):
            out.append((start, i))
            start = None
        if f and start is None:
            start = i
    if start is not None:
        out.append((start, len(flags)))
    return out


def component_below(h2d, level, seeds):
    """The cells below `level` joined (8-neighbour) to any of the flat indices in seeds."""
    from scipy import ndimage
    below = h2d < level
    labels, _ = ndimage.label(below, structure=np.ones((3, 3), bool))
    keep = np.unique(labels.ravel()[seeds])
    keep = keep[keep > 0]
    return np.isin(labels, keep).ravel()


def overflow_bodies(geo, sea, h):
    """A full sea's overflow: the hollows its spill path fills, each as (level, cells), and the path's other cells,
    from where the water leaves the sea's rim to where it enters the sea below."""
    path, filled = geo['spill_paths'][sea]
    rim = min(geo['passes'][(sea, j)] for j in geo['seas'] if j != sea)
    leave = np.flatnonzero(geo['joined'][sea][path] >= rim - N.BELOW_RIM_M)
    if not leave.size:
        return [], np.array([], np.int64)
    river, surface = path[leave[0]:], filled[leave[0]:].astype(np.float64)
    below = np.logical_or.reduce([geo['joined'][j][river] < geo['level0'] for j in geo['seas'] if j != sea])
    if below.any():                     # the river ends where it enters the sea below
        river, surface = river[:np.argmax(below)], surface[:np.argmax(below)]
    h2d = h.reshape(geo['h'].shape)
    under = surface > h[river] + 0.5
    others = np.logical_or.reduce([geo['joined'][j] < geo['level0'] for j in geo['seas'] if j != sea])
    hollows = []
    for a, b in runs_of(under, np.round(surface, 1)):
        level = float(surface[a:b].max())
        cells = component_below(h2d, level, river[a:b])
        if (cells & others).any():      # a hollow above a sea's rim cannot reach that sea's water
            raise SystemExit(f'a hollow on {sea}\'s overflow at {level - geo["level0"]:.1f} m reaches another sea')
        hollows.append((level, cells))
    return hollows, river[~under]


def stops(geo, case, log=print):
    """Each stop's water: per cell, its bodies, which is the surface at each stop and whether it is water; per body,
    its level as drawn at each stop (relative to the common level; nan where it holds no water); and per stop, each
    sea's level in the history itself."""
    h = geo['h'].ravel().astype(np.float64)
    level0, seas, joined = geo['level0'], geo['seas'], geo['joined']
    n_stops = len(STOP_YEARS)
    six0 = np.zeros(h.size, bool)
    for s in seas:
        six0 |= joined[s] < level0
    isolated = (h < level0) & ~six0                  # held at the common level in the atlas; not followed
    records = {r['year']: r for r in case['history']['levels_m']}
    passing = {s: case['seas'][s]['with_the_seas_above']['passes_on_m3s'] for s in seas}
    bodies = [dict(id=i, name=s, kind='sea') for i, s in enumerate(seas, start=1)]
    overflow = {}
    for s in seas:                                    # the overflow of every sea that is full at some stop
        if any(records.get(y, {}).get('states', {}).get(s) == 'full' for y in STOP_YEARS[1:]):
            hollows, channel = overflow_bodies(geo, s, h)
            ids = []
            for level, cells in hollows:
                ids.append(len(bodies) + 1)
                bodies.append(dict(id=ids[-1], name=f'{s} overflow hollow {len(ids)}', kind='hollow', fed_by=s,
                                   level_m=round(level - level0, 1), cells=int(cells.sum())))
            overflow[s] = dict(hollows=list(zip(ids, hollows)), channel=channel)
            log(f'{s}: {len(hollows)} hollows on its overflow, {channel.size} channel cells')
    if len(bodies) > 255:
        raise SystemExit(f'{len(bodies)} bodies of water do not fit one byte')
    # the history rounds its levels to 0.1 m, which can lift a full sea a few centimetres over its pass: each sea
    # stands at most just under its rim, as the history's own volumes hold it
    rims = {s: min(geo['passes'][(s, j)] for j in seas if j != s) - N.BELOW_RIM_M for s in seas}
    levels = np.full((n_stops, len(bodies) + 1), np.nan, np.float32)
    ident = np.zeros((n_stops, h.size), np.uint8)
    water = np.zeros(h.size, np.uint8)
    channel_bits = np.zeros(h.size, np.uint8)
    channel_q = np.zeros(h.size, np.float32)
    shares, states_by_stop, history = [], [], []
    for k, year in enumerate(STOP_YEARS):
        if year == 0:
            rel, states = {s: 0.0 for s in seas}, {s: 'filling' for s in seas}
        else:
            rec = records[year]
            rel, states = {s: rec[s] for s in seas}, rec['states']
        history.append(dict(rel))
        surface = np.full(h.size, -np.inf)
        for i, s in enumerate(seas, start=1):
            if rel[s] is None:
                raise SystemExit(f'{s} has joined another sea by year {year}; the stops do not follow joined seas')
            lev = min(level0 + rel[s] + (JOIN_WITHIN_M if year else 0.0), rims[s])    # year 0 is the atlas's own
            covered = (joined[s] < lev) & (lev > surface)
            surface[covered], ident[k][covered] = lev, i
            levels[k, i] = lev - level0
        for s, o in overflow.items():
            if states.get(s) != 'full':
                continue
            for i, (lev, cells) in o['hollows']:
                covered = cells & (lev > surface)
                surface[covered], ident[k][covered] = lev, i
                levels[k, i] = lev - level0
            dry = o['channel'][~np.isfinite(surface[o['channel']])]
            channel_bits[dry] |= np.uint8(1 << k)
            channel_q[dry] = np.maximum(channel_q[dry], passing[s])
        wet = (ident[k] > 0) | isolated
        water[ident[k] > 0] |= np.uint8(1 << k)
        states_by_stop.append(states)
        shares.append(wet)
    return dict(ident=ident, water=water, channel_bits=channel_bits, channel_q=channel_q, levels=levels,
                bodies=bodies, states=states_by_stop, history=history, wet=shares, isolated=isolated)


def two_bodies(ident):
    """Per cell, the first and a second body that cover it at some stop, the stops at which the second is its
    surface, and how many cells more than two bodies cover."""
    n_stops, n = ident.shape
    first = np.zeros(n, np.uint8)
    for k in range(n_stops):
        take = (first == 0) & (ident[k] > 0)
        first[take] = ident[k][take]
    second = np.zeros(n, np.uint8)
    for k in range(n_stops):
        take = (ident[k] > 0) & (ident[k] != first)
        second[take] = ident[k][take]
    second_bits = np.zeros(n, np.uint8)
    extra = np.zeros(n, bool)
    for k in range(n_stops):
        on_second = (ident[k] == second) & (second > 0)
        second_bits[on_second] |= np.uint8(1 << k)
        extra |= (ident[k] > 0) & (ident[k] != first) & (ident[k] != second)
    return first, second, second_bits, int(extra.sum())


def build(log=print):
    started = time.time()
    geo = N.geometry(log)
    design = json.loads(N.OUT.read_text())['cases'][0]
    factor = design['budget']['open_water_evaporation_factor_closing_it']
    case = N.network(1.0, geo, log, evaporation=factor, record=STOP_YEARS[1:])
    drawn = stops(geo, case, log)
    first, second, second_bits, extra = two_bodies(drawn['ident'])
    rows, cols = geo['h'].shape
    area = np.repeat((geo['cell_km'] * 1e3) ** 2 * np.cos(np.radians(geo['lat'])), cols)
    moon = float(area.sum())
    N.RUNS.mkdir(parents=True, exist_ok=True)
    meta = dict(schema=SCHEMA, grid=dict(rows=rows, cols=cols, left_lon=0.0, pixels_per_degree=cols / 360.0),
                stops=list(STOP_YEARS))
    np.savez_compressed(MAPS, body_bits=drawn['water'].reshape(rows, cols), first=first.reshape(rows, cols),
                        second=second.reshape(rows, cols), second_bits=second_bits.reshape(rows, cols),
                        channel_bits=drawn['channel_bits'].reshape(rows, cols),
                        channel_q=drawn['channel_q'].reshape(rows, cols), levels=drawn['levels'],
                        metadata=json.dumps(meta))
    seas = geo['seas']
    out = dict(
        schema=SCHEMA,
        producer=dict(module='research.studies.drainage_network.stops', code_sha256=code_hash()),
        evidence=('The drainage network study\'s history of the six large seas with the routing\'s budget closed, '
                  'recorded at chosen years and drawn on LOLA\'s 16 px/deg heights above the GRAIL geoid; the '
                  'network\'s own limits hold (no water\'s own gravity or crustal sag, the design climate fixed). '
                  'A computed drawing of a computed history, not validated.'),
        reading_rule=('stops are years after the water stands as the atlas draws it. Each stop gives each sea\'s '
                      'level in the history relative to the common level (levels_m, rounded there to 0.1 m) and its '
                      'state (filling, or full and overflowing), the '
                      'events up to that year, and the share of the Moon that is water: the six seas, the hollows '
                      'their overflows fill and the hollows the atlas holds at the common level, which the study '
                      'does not follow and which stay as drawn; the rain-fed lakes above the common level are '
                      'counted apart, as the atlas counts them. bodies lists every water the maps name: ids 1-6 '
                      'the seas, the rest the hollows an overflow fills, with the level each is filled to. The maps '
                      'hold, per cell of the 16 px/deg grid with 0 E at the left edge, a bit per stop where one of '
                      'the bodies covers it (body_bits; the hollows the atlas holds at the common level are water '
                      'in every year besides), the first and second body that cover the cell at some stop (first, '
                      'second), a bit '
                      'per stop where the second is its surface (second_bits), a bit per stop for an overflow '
                      'river on dry ground (channel_bits) with the surplus it carries (channel_q, m3/s), and each '
                      'body\'s level as drawn at each stop (levels, m, nan where it holds no water): after year 0 a '
                      'sea is drawn half a metre above its level in the history, so that a sea within half a metre '
                      'below a pass is drawn over it, and at most 1 cm under its rim.'),
        case=dict(rain=1.0, open_water_evaporation_factor=factor),
        maps=str(MAPS.relative_to(ROOT)), cells_under_more_than_two_bodies=extra,
        stops=[dict(year=y,
                    levels_m=drawn['history'][k],
                    states={s: drawn['states'][k].get(s, 'filling') for s in seas},
                    events=[e for e in case['history']['events'] if e['year'] <= y],
                    water_share=round(float(area[drawn['wet'][k]].sum() / moon), 5),
                    water_km2=round(float(area[drawn['wet'][k]].sum() / 1e6)))
               for k, y in enumerate(STOP_YEARS)],
        bodies=drawn['bodies'],
        held_as_the_atlas_draws_km2=round(float(area[drawn['isolated']].sum() / 1e6)),
        seconds=round(time.time() - started))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    log(f'wrote {OUT.relative_to(ROOT)} and {MAPS} in {out["seconds"]} s')
    return out


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    build()
    return 0


if __name__ == '__main__':
    sys.exit(main())
