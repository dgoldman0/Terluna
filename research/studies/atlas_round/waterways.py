"""The design surface's waterways: each major watercourse followed whole from its mouth up its main flow, through the
lakes on its way, so that it can be named, drawn and described as one thing.

    ulimit -v 8000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.waterways   # under a minute

On this Moon a large river is mostly lakes: the design routing (design_routing.py) fills every hollow on a river's way
until it spills, so a watercourse is a chain of long lakes joined by short channels. The sustained-rivers study
(sustained_rivers.py) classes the channels reach by reach, between the lakes; this module follows each watercourse from
where it enters a sea, or a lake that passes nothing on, upstream along the donor carrying the most water at every
step, through lakes and channels alike, until its flow falls under HEAD_M3S or it reaches the place its water is put in:
the rim a full sea overflows, a cut of the works, or the near-side sea's neck.

Each waterway of MOUTH_M3S or more at its mouth is recorded with its course (the settled frame's latitudes and
longitudes, thinned), its length, how much of it runs through lakes, the lakes' area along it, its flow at the mouth
and where it is greatest, the lowest year's flow and whether any reach of it falls under the sustained-rivers
threshold in some year, the steep reaches on it (falls.py), the regions it crosses (districts.py) and what feeds its
head. A waterway's anchor is its mouth's place on today's Moon, which a redraw of the regions keeps; the names table
(place_names.json) pins its name to that anchor, in the tongue of the region holding most of its course.
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np

from climate.gcm.turned import to_turned
from research.studies.atlas_round import design_routing as DRt
from research.studies.atlas_round import district_names as DN
from research.studies.atlas_round import districts as D
from research.studies.atlas_round import sustained_rivers as SR

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'waterways.json'
NAMES = HERE / 'place_names.json'
SCHEMA = 'terluna.research.atlas-round-waterways/1'
CODE = ('research/studies/atlas_round/waterways.py',)

MOUTH_M3S = 1000.0                   # a waterway carries at least this where it ends
HEAD_M3S = 100.0                     # and is followed upstream while it carries this (the sustained rivers' threshold)
PUT_IN_M3S = 1000.0                  # a head still carrying this with no donor is water put in: an overflow or a cut
MOON_KM = D.MOON_KM
# What feeds a waterway whose head is a place the design puts water in, by the design's name for it.
FED_BY = {'salt_sea_weir': 'the overflow of {moscoviense} over its natural pass',
          'orientale': 'the overflow of {orientale}', 'smythii_marginis': 'the overflow of {smythii_marginis}',
          'humboldtianum_branch': 'a branch the works cut from {moscoviense}',
          'smythii_marginis_branch': 'a branch the works cut from {moscoviense}',
          'neck': 'the spill of {near_side_main} over its neck'}
# The routing's seas and the districts draft's anchors for them.
SEA_ANCHORS = {'moscoviense': 'sea:Mare Moscoviense', 'orientale': 'sea:Mare Orientale',
               'smythii_marginis': 'sea:Smythii–Marginis Sea', 'south_pole_aitken': 'sea:South Pole–Aitken Sea',
               'humboldtianum': 'sea:Mare Humboldtianum', 'near_side_main': 'sea:near-side sea',
               'near_side_south': 'sea:southern near-side seas'}
THIN_KM = 4.0                        # the recorded course keeps a point about every this far


def sha256(path) -> str:
    return D.sha256(path)


def rel(path) -> str:
    return D.rel(path)


def main_donors(receiver, q):
    """For every cell, the cell draining into it that carries the most water, or -1."""
    order = np.argsort(q, kind='stable')
    order = order[(receiver[order] >= 0) & (receiver[order] != order)]
    donor = np.full(receiver.size, -1, np.int64)
    donor[receiver[order]] = order                       # a later, larger donor replaces a smaller one
    return donor


def follow(mouth, donor, q, added, head=HEAD_M3S, put_in=PUT_IN_M3S):
    """The cells of a waterway from its mouth upstream along the main donors while they carry `head`, ending at a
    cell where `put_in` or more is added beyond what its donors bring: a full sea's overflow or a branch of the
    works. Returns the cells and whether the head is such a place."""
    cells = [int(mouth)]
    while True:
        if added[cells[-1]] >= put_in:
            return cells, True
        up = int(donor[cells[-1]])
        if up < 0 or q[up] < head:
            return cells, False
        cells.append(up)


def straighten(cells, lakes):
    """A course with each run of lake cells replaced by its first and last cell: the routing wanders over a lake's
    level surface, and a lake is crossed from where the water comes in to where it leaves. Returns the kept cells,
    whether the step up to each is across a lake, and the number of lakes crossed."""
    kept, across, n, i = [], [], 0, 0
    while i < len(cells):
        if not lakes[i]:
            kept.append(cells[i])
            across.append(False)
            i += 1
            continue
        j = i
        while j + 1 < len(cells) and lakes[j + 1]:
            j += 1
        kept.append(cells[i])
        across.append(False)
        if j > i:
            kept.append(cells[j])
            across.append(True)
        n += 1
        i = j + 1
    return kept, across, n


def step_km(a, b, cols, ppd=16):
    """Great-circle distance between two cells' centres."""
    la, lo = 90.0 - (a // cols + 0.5) / ppd, (a % cols + 0.5) / ppd
    lb, lob = 90.0 - (b // cols + 0.5) / ppd, (b % cols + 0.5) / ppd
    u, v = D.unit(np.array(la), np.array(lo)), D.unit(np.array(lb), np.array(lob))
    return float(MOON_KM * np.arccos(np.clip((u * v).sum(axis=-1), -1.0, 1.0)))


def thin(lat, lon, km, every=THIN_KM):
    """Points of a course about `every` km apart along it, its ends kept."""
    keep, last = [0], 0.0
    for i in range(1, len(km)):
        if km[i] - last >= every:
            keep.append(i)
            last = km[i]
    if keep[-1] != len(km) - 1:
        keep.append(len(km) - 1)
    return [[round(float(lat[i]), 3), round(float(lon[i]), 3)] for i in keep]


def build(log=print):
    started = time.time()
    routing = json.loads(DRt.OUT.read_text())
    path = Path(routing['product']['path'])
    if sha256(path) != routing['product']['sha256']:
        raise SystemExit(f'{path} no longer matches its record {DRt.OUT}')
    with np.load(path, allow_pickle=False) as z:
        meta = json.loads(str(z['metadata']))
        receiver, q = z['receiver'].astype(np.int64), z['discharge_m3s'].astype(np.float64)
        sea, wet = z['sea'], z['routing_lake'].copy()
        lake_cell = z['lake_cell']
    rows, cols = meta['grid']['rows'], meta['grid']['cols']
    in_lake = np.zeros(rows * cols, bool)
    in_lake[lake_cell] = True
    in_lake |= wet
    rotation = np.asarray(meta['frame']['rotation_today_to_turned'], float)
    donor = main_donors(receiver, q)
    cell_km2 = (MOON_KM * np.radians(1 / 16)) ** 2 * np.cos(np.radians(90.0 - (np.arange(rows) + 0.5) / 16))

    # where waterways end: a land cell draining into a sea, or the last cell of a flow that goes nowhere further
    n = rows * cols
    idx = np.arange(n)
    ok = (receiver >= 0) & (receiver != idx)
    into_sea = ok & (sea < 0) & (sea[np.where(ok, receiver, 0)] >= 0) & (q >= MOUTH_M3S)
    ends_here = ~ok & (sea < 0) & (q >= MOUTH_M3S)
    mouths = np.flatnonzero(into_sea | ends_here)
    mouths = mouths[np.argsort(-q[mouths])]
    log(f'{mouths.size} mouths of {MOUTH_M3S:.0f} m3/s or more, {time.time() - started:.0f} s')

    # the sustained rivers' cells, the falls, the cuts, the regions
    with np.load(SR.PRODUCT, allow_pickle=False) as s:
        low = np.full(n, np.nan, np.float32)
        low[s['cell']] = s['sustained_m3s']
        falters = np.zeros(n, bool)
        falters[s['cell']] = s['intermittent'] & s['major']
    falls = json.loads(D.FALLS.read_text())['reaches'] if D.FALLS.exists() else []
    fall_at = {}
    for x in falls:
        fall_at.setdefault(int(x['start']), []).append(float(x['power_gw']))
    districts = json.loads(D.OUT.read_text())
    with np.load(D.PRODUCT) as g:
        region8 = g['region']
    land = [r for r in districts['regions'] if r['kind'] == 'land']
    seas = [r for r in districts['regions'] if r['kind'] == 'sea']
    names = {r['anchor']: r['name'] for r in seas}
    region_by_number = {int(r['code'][1:]): r for r in land}

    # where the design puts water in, and what each is (the sustained-rivers study keeps the design's structure)
    put_ins = {}
    if SR.STRUCTURE.exists():
        with np.load(SR.STRUCTURE, allow_pickle=False) as st:
            for c, k, sign in zip(st['op_cell'].tolist(), st['op_kind'].tolist(), st['op_sign'].tolist()):
                if sign > 0:
                    put_ins[int(c)] = SR.OPS[int(k)]
    inflow = np.bincount(receiver[ok], weights=q[ok], minlength=n)
    added = q - inflow                                   # what a cell adds beyond its donors: its own runoff, or water put in
    taken = np.zeros(n, bool)
    waterways = []
    for mouth in mouths:
        if taken[mouth]:
            continue
        cells, put_in = follow(mouth, donor, q, added)
        fresh = [c for c in cells if not taken[c]]
        if len(fresh) < len(cells):                      # joins a waterway already followed: its own part alone
            cells, put_in = fresh, False
            if not cells or q[cells[0]] < MOUTH_M3S:
                continue
        taken[cells] = True
        whole = np.asarray(cells)
        kept, across, lakes_crossed = straighten(cells, in_lake[whole])
        arr, across = np.asarray(kept), np.asarray(across)
        seg = np.array([0.0] + [step_km(a_, b_, cols) for a_, b_ in zip(kept[:-1], kept[1:])])
        km = np.cumsum(seg)
        la, lo = 90.0 - (arr // cols + 0.5) / 16, (arr % cols + 0.5) / 16
        lat_t, lon_t = to_turned(la, lo, rotation)
        end = int(receiver[mouth]) if receiver[mouth] >= 0 and receiver[mouth] != mouth else -1
        ends_in = meta['seas'][int(sea[end])] if end >= 0 and sea[end] >= 0 else None
        head = int(arr[-1])
        feeds = 'headwaters'
        if put_in:
            kind = put_ins.get(head)
            if kind is None and put_ins:                 # the nearest place the design puts water in
                near = min(put_ins, key=lambda c: step_km(c, head, cols))
                kind = put_ins[near] if step_km(near, head, cols) < 60.0 else None
            feeds = FED_BY.get(kind, 'water the design puts in')
            for key, anchor in SEA_ANCHORS.items():
                feeds = feeds.replace('{' + key + '}', names.get(anchor, key))
        regions = {}
        for c in cells:
            num = int(region8[(c // cols) // 2, (c % cols) // 2])
            if num in region_by_number:
                regions[num] = regions.get(num, 0) + 1
        ordered = sorted(regions.items(), key=lambda t: -t[1])
        lows = low[whole]
        channel = np.isfinite(lows)
        reach_falls = [g_ for c in cells for g_ in fall_at.get(c, [])]
        anchor = SEA_ANCHORS.get(ends_in)
        w = dict(
            anchor=f'river:{la[0]:.1f},{lo[0]:.1f}',
            mouth=dict(today=[round(float(la[0]), 2), round(float(lo[0]), 2)],
                       turned=[round(float(lat_t[0]), 2), round(float(lon_t[0]), 2)]),
            head=dict(today=[round(float(la[-1]), 2), round(float(lo[-1]), 2)],
                      turned=[round(float(lat_t[-1]), 2), round(float(lon_t[-1]), 2)]),
            ends_in=names.get(anchor, ends_in) if ends_in else
            ('a waterway it joins' if end >= 0 and taken[end] else 'a lake that passes nothing on'),
            ends_in_anchor=anchor, fed_by=feeds, works=feeds.startswith('a branch the works'),
            length_km=round(float(km[-1])), lake_km=round(float(seg[across].sum())),
            channel_km=round(float(seg[~across].sum())), lakes_crossed=lakes_crossed,
            flow_mouth_m3s=round(float(q[mouth])), flow_greatest_m3s=round(float(q[whole].max())),
            flow_head_m3s=round(float(q[head])),
            lowest_year_m3s=round(float(lows[channel][0])) if channel.any() else None,
            falters_km=round(float(1.9 * falters[whole].sum())),
            steep_reaches=len(reach_falls), steep_gw=round(float(sum(reach_falls)), 1),
            regions=[dict(name=region_by_number[k]['name'], anchor=region_by_number[k]['anchor'],
                          tongue=region_by_number[k]['tongue'], share=round(v / max(sum(regions.values()), 1), 2))
                     for k, v in ordered],
            course=thin(lat_t, lon_t, km))
        waterways.append(w)
    # two waterways cannot share an anchor
    seen = {}
    for w in waterways:
        k = seen.get(w['anchor'], 0)
        seen[w['anchor']] = k + 1
        if k:
            w['anchor'] += f'#{k + 1}'
    # the names table pins a waterway's name to its anchor (place_names.json, river:<mouth>)
    table = json.loads(NAMES.read_text())['names']
    for w in waterways:
        pin = table.get(w['anchor'])
        w['name'] = pin['name'] if pin else None
        w['tongue'] = DN.TONGUES[pin['tongue']]['label'] if pin else None
        w['meaning'] = pin.get('meaning', '') if pin else ''
    log(f'{len(waterways)} waterways, {time.time() - started:.0f} s')
    return dict(
        schema=SCHEMA, producer=dict(files={f: sha256(ROOT / f) for f in CODE}),
        inputs={rel(p): sha256(p) for p in (DRt.OUT, SR.OUT, D.OUT, D.FALLS, NAMES)},
        ignored_inputs={str(path): routing['product']['sha256'], str(SR.PRODUCT): sha256(SR.PRODUCT),
                        str(D.PRODUCT): sha256(D.PRODUCT), str(SR.STRUCTURE): sha256(SR.STRUCTURE)},
        evidence=('A reading of products: the design routing\'s receivers and flows followed upstream from each mouth '
                  'along the donor carrying the most water, with the sustained-rivers study\'s lowest-year flows, the '
                  'falls and the districts draft read along each course. Flows are the design climate\'s means on its '
                  'closed budget; the climate is an estimate.'),
        reading_rule=('One entry per waterway carrying mouth_m3s or more where it ends, largest first; a tributary is '
                      'listed apart from the waterway it joins, from the join upstream. Lengths in km along the 1.9 km '
                      'cells\' centres; lake_km is the length through lakes and channel_km the rest. flow_* in m3/s; '
                      'lowest_year_m3s is the sustained discharge at the channel nearest the mouth; falters_km the '
                      'length of channel that carries under the sustained-rivers threshold in some year. fed_by is '
                      'headwaters, a full sea\'s overflow, the near-side sea\'s spill over its neck, or a branch the works cut '
                      '(works marks those); regions gives each '
                      'region\'s share of the course\'s cells. course lists latitude and longitude in the settled frame from '
                      'the mouth upstream.'),
        parameters=dict(mouth_m3s=MOUTH_M3S, head_m3s=HEAD_M3S, put_in_m3s=PUT_IN_M3S, thin_km=THIN_KM),
        waterways=waterways, seconds=round(time.time() - started))


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    record = build()
    OUT.write_text(json.dumps(record, indent=1, ensure_ascii=False) + '\n')
    print(f"wrote {OUT}: {len(record['waterways'])} waterways")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
