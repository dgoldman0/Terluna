"""The fresh seas of the design surface and the rivers that enter them: each sea's area and volume by depth on the
redrawn floor, and for each large river its mouth, its course up to the last lake it leaves, and the floor it would
run down inside the sea.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.river_oxygen.mouths   # about a minute

It reads the design routing (research/studies/atlas_round/design_routing.py's product on the research drive) and the
redrawn floor (redrawn_drainage.load_surface), and writes results/mouths.json. A mouth is a cell off the seas whose
water goes straight into a sea. The rivers kept one by one are the largest that together carry SHARE of a sea's river
water, at most PER_SEA of them; the rest are counted together as the sea's small rivers.

- The course above a mouth follows the largest donor, through the small lakes on its way, to the first lake large
  enough to bring the flow through it to its own surface's temperature, or to a higher sea. That water is the
  river's source: its level, area and mean depth are kept, with the distance and the small lakes passed, since a
  river carries its source's surface water down to the sea and loses the difference only slowly.
- The floor inside the sea is followed from the mouth down the steepest neighbour to a lowest point, and on through
  every hollow the river's water of FILL_YEARS years would fill to its sill, to the hollow that holds it longer: the
  floor's depths with distance on the first descent give the slope an underflow runs on.
- That last hollow is measured as the water below each height above its lowest point that is joined to it: its area
  and volume say how large a bottom pool a river's cold water makes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from scipy import ndimage

from research.studies.atlas_round import design_routing as DR
from research.studies.atlas_round import redrawn_drainage as RD
from shared.constants import MOON_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'mouths.json'
SCHEMA = 'terluna.research.river-oxygen-mouths/1'
SALT_SEAS = ('south_pole_aitken',)                 # the salt sea is outside this study (register, 16 November 2026)
SHARE, PER_SEA = 0.95, 20                          # the rivers kept one by one in each sea
DEPTHS_M = (0, 10, 25, 50, 100, 200, 300, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 7000)
POOL_HEIGHTS_M = (5, 10, 25, 50, 100, 200, 400, 800)
COURSE_MAX_KM = 4000.0
POOL_WINDOW_DEG = 12.0                             # half-width of the window in which a hollow is filled
RESET_M2_PER_M3S = 1.7e5                           # lake area for each m3/s through it that brings the flow to the lake's
#                                                    own surface temperature: rho c / K with K = 24 W/m2/K (waters.py)
FILL_YEARS = 10.0                                  # a hollow the river fills to its sill sooner is passed through
FILL_HEIGHTS_M = (1, 2, 5, 10, 25, 50, 100, 200, 400, 800, 1600)
NEIGHBOURS = ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1))
YEAR_S = 31557600.0


def sha256(path) -> str:
    digest = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def step_km(r0, c0, r1, c1, lat, cols):
    """Great-circle distance (km) between two cells' centres."""
    la0, la1 = np.radians(lat[r0]), np.radians(lat[r1])
    dlon = np.radians(((c1 - c0 + cols // 2) % cols - cols // 2) * 360.0 / cols)
    a = np.sin((la1 - la0) / 2) ** 2 + np.cos(la0) * np.cos(la1) * np.sin(dlon / 2) ** 2
    return 2.0 * MOON_RADIUS / 1e3 * np.arcsin(np.sqrt(a))


def seas_table(sea, depth_of, row_area, names, levels):
    out = {}
    area = np.broadcast_to(row_area[:, None], sea.shape)
    for k, name in enumerate(names):
        m = sea == k
        d, a = depth_of(k)[m].astype(float), area[m]
        below = [dict(depth_m=z, area_km2=float(a[d > z].sum() / 1e6), volume_km3=float((a[d > z] * (d[d > z] - z)).sum() / 1e9))
                 for z in DEPTHS_M if (d > z).any()]
        out[name] = dict(level_m=float(levels[name]), area_km2=float(a.sum() / 1e6), volume_km3=float((a * d).sum() / 1e9),
                         mean_depth_m=float((a * d).sum() / a.sum()), deepest_m=float(d.max()), below=below,
                         salt=name in SALT_SEAS)
    return out


def donors_of(cell, receiver, discharge, rows, cols):
    """The largest donor among a cell's neighbours and its flow."""
    r, c = divmod(cell, cols)
    best, best_q = -1, 0.0
    for dr, dc in NEIGHBOURS:
        rr, cc = r + dr, (c + dc) % cols
        if 0 <= rr < rows:
            n = rr * cols + cc
            if receiver[n] == cell and discharge[n] > best_q:
                best, best_q = n, float(discharge[n])
    return best, best_q


def course_up(cell, rt, lakes, names, levels, h, lat, rows, cols):
    """Up a river from its mouth along its largest donor, through the small lakes on its way, to the water that sets
    its temperature: the first lake large enough to bring the flow through it to its own surface's temperature
    (area at least RESET_M2_PER_M3S times the flow through it), or a higher sea. Where the neighbours no longer carry
    half the flow an overflow lands, and the course goes on from the cell that sends it. Returns the distance, the
    small lakes passed with their areas over the flow through them, and the source."""
    receiver, discharge, sea = rt['receiver'], rt['discharge_m3s'], rt['sea']
    km, here, passed, seen = 0.0, cell, [], set()
    start_h = float(h.flat[cell])
    out = dict(kind='none within the search')
    while km < COURSE_MAX_KM:
        best, best_q = donors_of(here, receiver, discharge, rows, cols)
        if best < 0 or best_q < 0.5 * discharge[here]:
            far = np.flatnonzero(receiver == here)
            far = far[far != here]
            if far.size == 0 or discharge[far].max() < 0.25 * discharge[here]:
                out = dict(kind='the course ends', ground_m=float(h.flat[here]))
                break
            best = int(far[np.argmax(discharge[far])])
            best_q = float(discharge[best])
        r, c = divmod(here, cols)
        rb, cb = divmod(best, cols)
        km += float(step_km(r, c, rb, cb, lat, cols))
        here = best
        if sea[here] >= 0:
            out = dict(kind='sea', sea=names[int(sea[here])], level_m=float(levels[names[int(sea[here])]]))
            break
        i = lakes['index'].get(here)
        if i is not None and i not in seen:
            seen.add(i)
            flow = max(float(discharge[here]), 1e-9)
            ratio = float(lakes['area_km2'][i]) * 1e6 / flow
            if ratio >= RESET_M2_PER_M3S:
                out = dict(kind='lake', level_m=float(lakes['level'][i]), area_km2=float(lakes['area_km2'][i]),
                           mean_depth_m=float(lakes['mean_depth_m'][i]), flow_through_m3s=flow)
                break
            passed.append(dict(km=km, area_km2=float(lakes['area_km2'][i]), level_m=float(lakes['level'][i]),
                               m2_per_m3s=ratio))
    r, c = divmod(here, cols)
    out.update(lat=float(lat[r]), lon_east=float((c + 0.5) * 360.0 / cols))
    return dict(km=km, fall_m=float(out.get('level_m', h.flat[here])) - start_h, small_lakes_passed=len(passed),
                small_lakes_m2_per_m3s=float(sum(x['m2_per_m3s'] for x in passed)), source=out)


def lakes_table(rt, h, row_area, cols):
    """The lakes off the seas by their levels: each lake's area and mean depth, and the lake each lake cell lies in."""
    cell, level = rt['lake_cell'], rt['lake_level_m']
    values, inverse = np.unique(level, return_inverse=True)
    a = row_area[cell // cols]
    depth = np.maximum(level.astype(float) - h.flat[cell], 0.0)
    area = np.bincount(inverse, weights=a, minlength=values.size)
    volume = np.bincount(inverse, weights=a * depth, minlength=values.size)
    return dict(index=dict(zip(cell.tolist(), inverse.tolist())), level=values, area_km2=area / 1e6,
                mean_depth_m=volume / np.maximum(area, 1.0))


def window(lowest, sea_k, sea, h, row_area, rows, cols):
    r0, c0 = divmod(lowest, cols)
    half = int(POOL_WINDOW_DEG * cols / 360.0)
    ra, rb = max(0, r0 - half), min(rows, r0 + half + 1)
    cc = (np.arange(c0 - half, c0 + half + 1)) % cols
    return (h[ra:rb][:, cc].astype(float), sea[ra:rb][:, cc] == sea_k, np.broadcast_to(row_area[ra:rb, None], (rb - ra, cc.size)),
            r0 - ra, half, ra, cc)


def last_hollow(lowest, flow_m3s, sea_k, sea, h, row_area, rows, cols):
    """From a lowest point of the floor on through each hollow the river would fill to its sill within FILL_YEARS, to
    the lowest point of the hollow that holds its water longer; the hollows passed are counted."""
    limit = flow_m3s * YEAR_S * FILL_YEARS
    passed = 0
    for _ in range(60):
        sub_h, sub_sea, area, rr, cc0, ra, cc = window(lowest, sea_k, sea, h, row_area, rows, cols)
        floor = sub_h[rr, cc0]
        moved = False
        for height in FILL_HEIGHTS_M:
            labels, _ = ndimage.label(sub_sea & (sub_h < floor + height), structure=np.ones((3, 3)))
            m = labels == labels[rr, cc0]
            deeper = m & (sub_h < floor - 0.5)
            if deeper.any():                                   # the water has run over a sill into deeper ground
                j = np.argmin(np.where(m, sub_h, np.inf))
                jr, jc = divmod(int(j), sub_h.shape[1])
                lowest, passed, moved = int((ra + jr) * cols + cc[jc]), passed + 1, True
                break
            if (area[m] * (floor + height - sub_h[m])).sum() >= limit:
                break
        if not moved:
            break
    return lowest, passed


def descent(cell, sea_k, sea, h, level, lat, rows, cols):
    """Down the steepest neighbour over the sea's floor from the sea cell a river enters to a lowest point: distance
    (km) and depth (m) along the way, thinned."""
    here, km = cell, 0.0
    dist, depth = [0.0], [float(level - h.flat[cell])]
    while True:
        r, c = divmod(here, cols)
        best, best_slope = -1, 0.0
        for dr, dc in NEIGHBOURS:
            rr, cc = r + dr, (c + dc) % cols
            if 0 <= rr < rows and sea[rr, cc] == sea_k:
                n = rr * cols + cc
                s = float(h.flat[here] - h.flat[n]) / (float(step_km(r, c, rr, cc, lat, cols)) * 1e3)
                if s > best_slope:
                    best, best_slope = n, s
        if best < 0:
            break
        rb, cb = divmod(best, cols)
        km += float(step_km(r, c, rb, cb, lat, cols))
        here = best
        dist.append(km); depth.append(float(level - h.flat[here]))
    return here, np.array(dist), np.array(depth)


def pool(lowest, sea_k, sea, h, row_area, rows, cols):
    """The hollow round a lowest point of a sea's floor: the area and volume of the water below each height above it
    that is joined to it, within a window."""
    sub_h, sub_sea, area, rr, cc0, _, _ = window(lowest, sea_k, sea, h, row_area, rows, cols)
    floor = float(h.flat[lowest])
    out = []
    for height in POOL_HEIGHTS_M:
        labels, _ = ndimage.label(sub_sea & (sub_h < floor + height), structure=np.ones((3, 3)))
        lab = labels[rr, cc0]
        if lab == 0:
            continue
        m = labels == lab
        edge = bool(m[0].any() or m[-1].any() or m[:, 0].any() or m[:, -1].any())
        out.append(dict(height_m=height, area_km2=float(area[m].sum() / 1e6),
                        volume_km3=float((area[m] * (floor + height - sub_h[m])).sum() / 1e9), reaches_window_edge=edge))
    return out


def build(log=print):
    started = time.time()
    rt = DR.load(log=log)
    meta = rt['meta']
    surf = RD.load_surface(log=lambda *a: None)
    h, lat, row_area = surf['h'], np.asarray(surf['lat'], float), np.asarray(surf['row_area'], float)
    rows, cols = h.shape
    names, levels = meta['seas'], meta['levels_m']
    sea = rt['sea'].reshape(rows, cols)
    receiver, discharge = rt['receiver'], rt['discharge_m3s']
    lakes = lakes_table(rt, h, row_area, cols)
    seas = seas_table(sea, lambda k: levels[names[k]] - h, row_area, names, levels)
    log(f'seas measured, {time.time() - started:.0f} s')
    n = receiver.size
    cells = np.arange(n)
    into = np.where(receiver != cells, rt['sea'][receiver], -1)
    is_mouth = (rt['sea'] < 0) & (into >= 0)
    rivers = {}
    for k, name in enumerate(names):
        m = np.flatnonzero(is_mouth & (into == k))
        q = discharge[m]
        order = np.argsort(q)[::-1]
        total = float(q.sum())
        seas[name].update(rivers_m3s=total, rivers_km3_yr=total * YEAR_S / 1e9, mouths=int(m.size),
                          rivers_share_of_volume_a_year=total * YEAR_S / 1e9 / seas[name]['volume_km3'])
        if name in SALT_SEAS:
            continue
        cum = np.cumsum(q[order]) / total
        keep = min(PER_SEA, int(np.searchsorted(cum, SHARE)) + 1)
        kept = []
        for i in order[:keep]:
            cell = int(m[i])
            r, c = divmod(cell, cols)
            up = course_up(cell, rt, lakes, names, levels, h, lat, rows, cols)
            first = int(receiver[cell])
            pit, dist, depth = descent(first, k, sea, h, levels[name], lat, rows, cols)
            lowest, passed = last_hollow(pit, float(q[i]), k, sea, h, row_area, rows, cols)
            lr, lc = divmod(lowest, cols)
            marks = {}
            for z in (25, 50, 100, 200, 500, 1000):
                j = int(np.argmax(depth >= z)) if (depth >= z).any() else -1
                if j > 0:
                    marks[str(z)] = dict(km=float(dist[j]), slope=float(depth[j] / (dist[j] * 1e3)))
            kept.append(dict(cell=cell, lat=float(lat[r]), lon_east=float((c + 0.5) * 360.0 / cols), flow_m3s=float(q[i]),
                             ground_m=float(h[r, c]), above_the_sea_m=float(h[r, c] - levels[name]), course=up,
                             floor=dict(first_lowest_depth_m=float(depth[-1]), first_lowest_km=float(dist[-1]),
                                        hollows_passed=passed, lowest_depth_m=float(levels[name] - h.flat[lowest]),
                                        lowest_lat=float(lat[lr]), lowest_lon_east=float((lc + 0.5) * 360.0 / cols),
                                        depth_reached=marks,
                                        pool=pool(lowest, k, sea, h, row_area, rows, cols))))
        rivers[name] = dict(kept=kept, kept_m3s=float(q[order[:keep]].sum()), small_m3s=float(q[order[keep:]].sum()),
                            small_mouths=int(m.size - keep),
                            small_median_m3s=float(np.median(q[order[keep:]][q[order[keep:]] > 0])) if (q[order[keep:]] > 0).any() else 0.0)
        log(f'{name}: {keep} rivers kept with {rivers[name]["kept_m3s"]:,.0f} of {total:,.0f} m3/s, {time.time() - started:.0f} s')
    if meta['neck'].get('joined'):
        seas['near_side_south']['through_the_neck_m3s'] = float(meta['neck']['flow_m3s'])
    return dict(seas=seas, rivers=rivers, common_level_m=float(meta['common_level_m']), routing=meta)


def write(result, out=OUT, seconds=None):
    record = dict(
        schema=SCHEMA,
        producer=dict(files={str(p.relative_to(ROOT)): sha256(p) for p in (Path(__file__).resolve(), Path(DR.__file__).resolve(),
                                                                          Path(RD.__file__).resolve())}),
        inputs={'research/studies/atlas_round/results/design_routing.json': sha256(DR.OUT)},
        ignored_inputs={str(DR.PRODUCT): sha256(DR.PRODUCT)},
        evidence=('Geometry and routing read from products: the design surface\'s seas on the redrawn floor at 16 '
                  'pixels a degree and the design routing\'s mean flows on the estimated design climate. A mouth is a '
                  'cell; a river\'s course and the floor below its mouth follow single cells of about 1.9 km, so slopes '
                  'and hollows finer than that are not seen. Nothing is measured.'),
        reading_rule=('Levels in m above the dry Moon\'s geoid, depths in m below each sea\'s own level, areas km2, '
                      'volumes km3, flows m3/s. seas[name].below lists the area and volume of the water deeper than '
                      'each depth. rivers[name].kept are the largest rivers that carry 95% of the sea\'s river water '
                      '(at most 20), each with its course up to the water that sets its temperature (course.km, '
                      'course.fall_m from that water\'s level to the mouth, the small lakes passed, and course.source: '
                      'the first lake holding at least 1.7e5 m2 for each m3/s through it, with its level, area and '
                      'mean depth, or a higher sea), '
                      'the floor below its mouth down the steepest way to a first lowest point '
                      '(floor.depth_reached[z]: the distance at which the floor first lies z m deep and the mean slope '
                      'to there), the hollows its water of ten years would fill and pass (floor.hollows_passed) and '
                      'the hollow that holds it longer (floor.lowest_depth_m and floor.pool: the joined water below '
                      'each height above its lowest point, within 12 degrees). small_m3s is the rest of the sea\'s '
                      'river water.'),
        climate=result['routing']['climate'], common_level_m=result['common_level_m'],
        levels_m=result['routing']['levels_m'], neck=result['routing']['neck'],
        seas=result['seas'], rivers=result['rivers'])
    if seconds is not None:
        record['seconds'] = seconds
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(record, indent=1) + '\n')
    return record


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--out', type=Path, default=OUT)
    args = p.parse_args(argv)
    started = time.time()
    log = lambda *a: print(*a, flush=True)                                     # noqa: E731
    write(build(log), args.out, round(time.time() - started))
    log(f'wrote {args.out}, {time.time() - started:.0f} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
