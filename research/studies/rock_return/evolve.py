"""How a region's ground changes through millions of years: the rock left where it settles, or returned by a rule.

    /media/projectspace/terluna-research/venvs/optics/bin/python -m research.studies.rock_return.evolve   # about ten minutes on four cores; wants numba

placement.py asks, for one moment, whether a rule of placement gives each cell back what it lost. This module lets
the ground answer through time. It takes squares of the design surface (placement.py writes them) and steps each
forward:

- water runs downhill and cuts the ground it crosses at a rate that rises with the land drained and the slope
  (the stream-power law, solved implicitly down the drainage; Braun and Willett 2013), scaled so that the region's
  first mean lowering is the rock module's;
- what is cut is carried to the first hollow on its way and laid down there from the bottom up, a level floor;
  what reaches the region's edge or a sea leaves;
- carriers take a share of what comes to rest at each hollow, and of what leaves at the edge or the sea, and set it
  down on the land round it by a rule: on worn ground in proportion to its wear, on worn ground up to what each
  piece lost and no more (healed), on ground worn by district, on crests, evenly, or not at all.
  Nearness is a bell on the ground, with about 86% of a place's rock set down within the range.

Two cases bracket the water. In the wet case every hollow holds a lake to its lip and spills, so that water crosses
it and cuts its outlet, and ground under the lake is not cut. In the dry case no hollow holds water: its walls are
cut to its floor and nothing leaves it. The design surface lies between.

The measures ask what the author of the design would look for: how far the ground has moved, whether the map's
forms at tens of kilometres are still its own, whether the hollows keep their depth and their roundness, and what
became of the rock. Pictures of each region under each rule go to the research drive beside the snapshots.
Writes results/evolution.json.
"""
from __future__ import annotations

import argparse
import hashlib
import heapq
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy import ndimage

try:                                                   # the kernels run compiled where numba is installed
    from numba import njit
    COMPILED = True
except ImportError:                                    # and as plain Python on small grids, as the tests use them
    COMPILED = False

    def njit(*args, **kwargs):
        if args and callable(args[0]):
            return args[0]
        return lambda f: f

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results' / 'evolution.json'
SCHEMA = 'terluna.research.rock-return-evolution/1'
DRIVE = Path('/media/projectspace/terluna-research/terrain-return/evolution')
REGIONS = DRIVE / 'regions_16ppd.npz'
PLACEMENT = HERE / 'results' / 'placement.json'
SNAPSHOTS = DRIVE / 'snapshots.npz'
AREA_EXPONENT = 0.5                  # the stream-power law's exponent on the land drained, with slope to the first power (the usual pair)
STEP_YR = 25000.0
YEARS = (1.0e6, 3.0e6, 1.0e7, 3.0e7)
RANGE_M = 20000.0                    # carriers set most rock down within this distance of where it lies (the placement module's middle range)
FAR_RANGE_M = 50000.0                # the longer range tried
DISTRICT_M = 10000.0                 # the distance over which the district rule judges wear
SURROUNDINGS_M = 5000.0              # the distance a crest is judged against
SURROUNDINGS_CELLS = 3               # the half-width, in cells, of the ground whose height range is the local relief (about 5 km)
FLOOR_ROUGHNESS_M = 0.01             # laid floors are left this rough, so that water crossing them does not follow the grid
MARGIN = 32                          # cells left out of the measures along each edge, where the open edge acts
FORM_CELLS = 8                       # forms are the ground against its mean over this half-width (about 30 km across)
HOLLOWS_FOLLOWED = 12
RULES = ('none', 'worn', 'healed', 'worn_by_district', 'crest', 'even')
CASES = dict(none=dict(rule='none'), worn=dict(rule='worn'), healed=dict(rule='healed'), worn_far=dict(rule='worn', range_m=FAR_RANGE_M),
             worn_six_tenths=dict(rule='worn', reach=0.6), worn_by_district=dict(rule='worn_by_district'),
             crest=dict(rule='crest'), even=dict(rule='even'))
WATERS = ('wet', 'dry')
EPSILON_M = 1.0e-6                   # the rise per cell that keeps water moving across a filled hollow
UNDER_M = 1.0e-3                     # ground this far under its hollow's lip counts as in the hollow


def sha256(path) -> str:
    digest = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


# ---- the kernels -------------------------------------------------------------------------------------------------
@njit(cache=True)
def flood(h, outlet, ny, nx, eps):
    """Fill every hollow to its lip from the outlets inward (priority flood; Barnes and colleagues 2014). Returns the
    filled surface, for each cell the neighbour its water goes on to, and the cells from the lowest level up."""
    n = ny * nx
    filled = h.copy()
    receiver = np.arange(n)
    order = np.empty(n, np.int64)
    done = np.zeros(n, np.bool_)
    heap = [(-1.0e300, -1)]
    for i in range(n):
        if outlet[i]:
            done[i] = True
            heapq.heappush(heap, (h[i], i))
    k = 0
    while len(heap) > 0:
        z, i = heapq.heappop(heap)
        if i < 0:
            continue
        order[k] = i
        k += 1
        y = i // nx
        x = i - y * nx
        for dy in range(-1, 2):
            yy = y + dy
            if yy < 0 or yy >= ny:
                continue
            for dx in range(-1, 2):
                xx = x + dx
                if xx < 0 or xx >= nx or (dy == 0 and dx == 0):
                    continue
                j = yy * nx + xx
                if done[j]:
                    continue
                done[j] = True
                filled[j] = max(h[j], z + eps)
                receiver[j] = i
                heapq.heappush(heap, (filled[j], j))
    return filled, receiver, order[:k]


@njit(cache=True)
def steepest(h, outlet, ny, nx, east_west, north_south):
    """Each cell's steepest way down over the bare ground; a cell with no lower neighbour, or an outlet, keeps itself."""
    n = ny * nx
    receiver = np.arange(n)
    for i in range(n):
        if outlet[i]:
            continue
        y = i // nx
        x = i - y * nx
        best = 0.0
        for dy in range(-1, 2):
            yy = y + dy
            if yy < 0 or yy >= ny:
                continue
            for dx in range(-1, 2):
                xx = x + dx
                if xx < 0 or xx >= nx or (dy == 0 and dx == 0):
                    continue
                j = yy * nx + xx
                d = math.sqrt((dy * north_south) ** 2 + (dx * east_west) ** 2)
                s = (h[i] - h[j]) / d
                if s > best:
                    best = s
                    receiver[i] = j
    return receiver


@njit(cache=True)
def lengths(receiver, nx, east_west, north_south):
    n = receiver.size
    out = np.ones(n)
    for i in range(n):
        j = receiver[i]
        if j != i:
            dy = i // nx - j // nx
            dx = (i - (i // nx) * nx) - (j - (j // nx) * nx)
            out[i] = math.sqrt((dy * north_south) ** 2 + (dx * east_west) ** 2)
    return out


@njit(cache=True)
def accumulate(receiver, order, cell_area):
    """The land that drains through each cell (m2); order lists every cell after the cell it drains to."""
    a = np.full(receiver.size, cell_area)
    for k in range(order.size - 1, -1, -1):
        i = order[k]
        j = receiver[i]
        if j != i:
            a[j] += a[i]
    return a


@njit(cache=True)
def erode(h, receiver, order, drained, length, wet, level, rate, exponent, dt):
    """The ground after one step of cutting, implicit in the height of the cell below. A cell under standing water is
    not cut, and a cell that drains into standing water is cut no lower than the water's surface."""
    new = h.copy()
    for k in range(order.size):
        i = order[k]
        j = receiver[i]
        if j == i or wet[i]:
            continue
        base = level[j] if wet[j] else new[j]
        if h[i] <= base:
            continue
        f = rate * drained[i] ** exponent * dt / length[i]
        new[i] = (h[i] + f * base) / (1.0 + f)
    return new


@njit(cache=True)
def entries(receiver, order, in_hollow, outlet):
    """For each cell, the cell where what is cut from it comes to rest: the first cell on its way down that lies in a
    hollow (itself, if it does), or the outlet where it leaves the region, or -1 if its way ends on level ground
    outside any hollow."""
    entry = np.full(receiver.size, -1)
    for k in range(order.size):
        i = order[k]
        if in_hollow[i] or outlet[i]:
            entry[i] = i
        elif receiver[i] != i:
            entry[i] = entry[receiver[i]]
    return entry


@njit(cache=True)
def fill_hollows(h, filled, cells, starts, volume, cell_area):
    """Lay each hollow's sediment down from its lowest ground up to one level, no higher than its lip. cells lists
    each hollow's cells from the lowest, hollow after hollow; returns the volume that did not fit."""
    over = 0.0
    for q in range(volume.size):
        v = volume[q]
        s0, s1 = starts[q], starts[q + 1]
        if v <= 0.0 or s1 <= s0:
            over += max(v, 0.0)
            continue
        lip = filled[cells[s0]]
        for p in range(s0, s1):
            if filled[cells[p]] < lip:
                lip = filled[cells[p]]
        count = 1
        level = h[cells[s0]]
        while v > 0.0:
            nxt = lip
            if s0 + count < s1 and h[cells[s0 + count]] < lip:
                nxt = h[cells[s0 + count]]
            need = (nxt - level) * count * cell_area
            if need >= v:
                level += v / (count * cell_area)
                v = 0.0
            else:
                v -= need
                level = nxt
                if s0 + count < s1 and level < lip:
                    count += 1
                else:
                    break
        for p in range(s0, s0 + count):
            if h[cells[p]] < level:
                h[cells[p]] = level
        over += v
    return over


# ---- one step ----------------------------------------------------------------------------------------------------
def near(field, g, distance_m):
    """A field spread over the ground round each cell: a bell of half-width distance_m / 2 in both directions on the
    ground, so that about 86% lies within distance_m. Sums are kept except for what passes the region's edge."""
    sigma = (0.5 * distance_m / g['north_south'], 0.5 * distance_m / g['east_west'])
    return ndimage.gaussian_filter(field.reshape(g['shape']), sigma=sigma, mode='constant', truncate=3.0).ravel()


def weights(rule, h, cut, g):
    """Each rule's weight per cell: where carriers that follow it would set rock down."""
    if rule == 'worn':
        return cut
    if rule == 'worn_by_district':
        return near(cut, g, DISTRICT_M)
    if rule == 'crest':
        sigma = (0.5 * SURROUNDINGS_M / g['north_south'], 0.5 * SURROUNDINGS_M / g['east_west'])
        around = ndimage.gaussian_filter(h.reshape(g['shape']), sigma=sigma, mode='nearest', truncate=3.0).ravel()
        return np.maximum(h - around, 0.0)
    if rule == 'even':
        return np.ones_like(h)
    raise ValueError(rule)


def set_down(lying, weight, g, range_m=RANGE_M):
    """Share the rock lying at each cell over the land round it, by weight and by nearness. Returns what each cell
    gets and the rock lying where no weighted land is near, cell by cell."""
    reach = near(weight, g, range_m)
    ok = (reach > 1e-300) & (lying > 0)
    per = np.where(ok, lying / np.where(ok, reach, 1.0), 0.0)
    return weight * near(per, g, range_m), np.where(ok, 0.0, lying)


def heal(lying, lost, g, range_m=RANGE_M, passes=3):
    """Carriers that fill what was worn and no more. The rock lying at each place is offered to the ground round it
    in proportion to what each cell has lost; a cell takes no more than it lost, and what it turns away is offered
    again, to the ground still short, up to passes times. Returns what each cell gets and what is left lying."""
    back = np.zeros_like(lost)
    left, want = lying.astype(np.float64).copy(), lost.astype(np.float64).copy()
    for _ in range(passes):
        reach = near(want, g, range_m)
        ok = (reach > 1e-300) & (left > 0)
        if not ok.any():
            break
        per = np.where(ok, left / np.where(ok, reach, 1.0), 0.0)
        offer = want * near(per, g, range_m)
        taken = np.where(offer > want, want / np.where(offer > 0, offer, 1.0), 1.0)      # the share of an offer a cell takes
        given = offer * taken
        back += given
        left = np.maximum(left - per * near(want * taken, g, range_m), 0.0)
        want = np.maximum(want - given, 0.0)
    return back, left


def step(h, g, water, rule, rate, dt, reach=1.0, range_m=RANGE_M):
    """One step of a region's ground. g holds the grid (shape, outlet, the cell's sides); returns the new ground and
    the step's volumes: cut, laid in hollows, set down by carriers, gone from the region."""
    ny, nx = g['shape']
    outlet, cell_area = g['outlet'], g['cell_area']
    filled, down, upward = flood(h, outlet, ny, nx, EPSILON_M)
    in_hollow = (filled - h > UNDER_M) & ~outlet
    if water == 'wet':
        receiver, order, wet = down, upward, in_hollow
    else:
        receiver = steepest(h, outlet, ny, nx, g['east_west'], g['north_south'])
        order = np.argsort(h, kind='stable')
        wet = np.zeros(h.size, np.bool_)
    drained = accumulate(receiver, order, cell_area)
    new = erode(h, receiver, order, drained, lengths(receiver, nx, g['east_west'], g['north_south']), wet, filled, rate, AREA_EXPONENT, dt)
    cut = h - new
    entry = entries(receiver, order, in_hollow, outlet)
    volume = cut * cell_area
    arrives = np.bincount(entry[entry >= 0], weights=volume[entry >= 0], minlength=h.size)
    gone = float(volume[entry < 0].sum())
    back = np.zeros_like(h)
    if rule != 'none' and reach > 0:
        # carriers work from every place the rock comes to rest, the region's edge and the sea among them: what
        # crosses the edge settles in the next hollow beyond it, within their range of this side
        land = ~wet & ~outlet
        if rule == 'healed':
            back, unplaced = heal(reach * arrives, np.where(land, volume, 0.0), g, range_m)
        else:
            back, unplaced = set_down(reach * arrives, np.where(land, weights(rule, new, cut, g), 0.0), g, range_m)
        arrives = (1.0 - reach) * arrives + unplaced
    gone += float(arrives[outlet].sum())
    labels, count = ndimage.label(in_hollow.reshape(ny, nx), structure=np.ones((3, 3), bool))
    labels = labels.ravel()
    laid = 0.0
    if count:
        inside = np.flatnonzero(labels > 0)
        per_hollow = np.bincount(labels[inside] - 1, weights=arrives[inside], minlength=count)
        cells = inside[np.lexsort((new[inside], labels[inside]))]
        starts = np.concatenate(([0], np.cumsum(np.bincount(labels[inside] - 1, minlength=count)))).astype(np.int64)
        before = new[cells].copy()
        over = fill_hollows(new, filled, cells, starts, per_hollow, cell_area)
        if g.get('rng') is not None:
            raised = cells[new[cells] > before]
            new[raised] += g['rng'].uniform(-FLOOR_ROUGHNESS_M, FLOOR_ROUGHNESS_M, raised.size)
        laid = float(per_hollow.sum() - over)
        gone += float(over)
    new = new + back / cell_area
    return new, dict(cut=float(volume.sum()), laid=laid, set_down=float(back.sum()), gone=gone)


def grid_of(region):
    ny, nx = region['ground'].shape
    outlet = np.zeros((ny, nx), bool)
    outlet[0, :] = outlet[-1, :] = outlet[:, 0] = outlet[:, -1] = True
    outlet |= region['sea']
    return dict(shape=(ny, nx), outlet=outlet.ravel(), east_west=float(region['east_west_m']), north_south=float(region['north_south_m']),
                cell_area=float(region['east_west_m'] * region['north_south_m']))


def start_of(region):
    """The ground a run starts from: the redrawn floor, with the sea's surface over the sea's cells."""
    return np.where(region['sea'], region['water'], region['ground']).astype(np.float64).ravel()


def interior(region):
    ny, nx = region['ground'].shape
    m = np.zeros((ny, nx), bool)
    m[MARGIN:ny - MARGIN, MARGIN:nx - MARGIN] = True
    return (m & ~region['sea']).ravel()


def calibrate(region, water):
    """The stream-power rate at which the region's first mean lowering is the rock module's."""
    g, h = grid_of(region), start_of(region)
    inside = interior(region)
    target = float(region['lowering'].ravel()[inside & region['eroding'].ravel()].mean())
    trial = 1.0e-7
    # the hollows' fill is not lowering: count the cut alone
    filled, down, upward = flood(h, g['outlet'], *g['shape'], EPSILON_M)
    in_hollow = (filled - h > UNDER_M) & ~g['outlet']
    if water == 'wet':
        receiver, order, wet = down, upward, in_hollow
    else:
        receiver, order, wet = steepest(h, g['outlet'], *g['shape'], g['east_west'], g['north_south']), np.argsort(h, kind='stable'), np.zeros(h.size, np.bool_)
    drained = accumulate(receiver, order, g['cell_area'])
    cut = h - erode(h, receiver, order, drained, lengths(receiver, g['shape'][1], g['east_west'], g['north_south']), wet, filled, trial, AREA_EXPONENT, 100.0)
    got = float(cut[inside & region['eroding'].ravel()].mean() / 100.0)
    return trial * target / got, target


def run(region, water, rule, years=YEARS, dt=STEP_YR, reach=1.0, rate=None, range_m=RANGE_M, seed=None):
    """A region's ground at each of the years, with the volumes moved up to each. With a seed, laid floors are left
    a centimetre rough."""
    g, h = grid_of(region), start_of(region)
    if seed is not None:
        g['rng'] = np.random.default_rng(seed)
    if rate is None:
        rate, _ = calibrate(region, water)
    marks = sorted(years)
    shots, books, total = {}, {}, dict(cut=0.0, laid=0.0, set_down=0.0, gone=0.0)
    t = 0.0
    for mark in marks:
        while t < mark - 1e-6:
            h, moved = step(h, g, water, rule, rate, min(dt, mark - t), reach, range_m)
            for k in total:
                total[k] += moved[k]
            t += min(dt, mark - t)
        shots[mark] = h.copy()
        books[mark] = dict(total)
    return shots, books, rate


# ---- the measures ------------------------------------------------------------------------------------------------
def forms(h, shape):
    """The ground against its own mean over the forms' half-width: what is left is the map's forms at that scale."""
    g = h.reshape(shape)
    return (g - ndimage.uniform_filter(g, size=2 * FORM_CELLS + 1, mode='nearest')).ravel()


def local_relief(h, shape, half=SURROUNDINGS_CELLS):
    g = h.reshape(shape)
    size = 2 * half + 1
    return (ndimage.maximum_filter(g, size=size, mode='nearest') - ndimage.minimum_filter(g, size=size, mode='nearest')).ravel()


def hollows_of(h, g):
    filled, _, _ = flood(h, g['outlet'], *g['shape'], EPSILON_M)
    depth = np.where(g['outlet'], 0.0, filled - h)
    labels, count = ndimage.label((depth > UNDER_M).reshape(g['shape']), structure=np.ones((3, 3), bool))
    return depth, labels, count


def roundness(mask, east_west, north_south):
    """How far a patch of cells is from a disc: the spread of its edge cells' distances from its centre over their
    mean (0 for a disc; a 2 by 1 ellipse gives about 0.2). None for a patch too small to say."""
    if mask.sum() < 12:
        return None
    edge = mask & ~ndimage.binary_erosion(mask, structure=np.ones((3, 3), bool))
    yy, xx = np.nonzero(mask)
    cy, cx = yy.mean() * north_south, xx.mean() * east_west
    ey, ex = np.nonzero(edge)
    r = np.hypot(ey * north_south - cy, ex * east_west - cx)
    return float(r.std() / r.mean()) if r.mean() > 0 else None


def followed_hollows(region):
    """The largest hollows of the starting ground that lie wholly inside the measured part: each one's cells, a ring
    of ground round it, its depth below the ring and the roundness of its outline at half that depth."""
    g, h = grid_of(region), start_of(region)
    depth, labels, count = hollows_of(h, g)
    inside = interior(region).reshape(g['shape'])
    sizes = ndimage.sum(np.ones(g['shape']), labels, index=np.arange(1, count + 1))
    out = []
    for q in np.argsort(sizes)[::-1][:300]:
        mask = labels == q + 1
        ring = ndimage.binary_dilation(mask, structure=np.ones((3, 3), bool), iterations=2) & ~mask
        if not (inside[mask].all() and inside[ring].all()):
            continue
        out.append(dict(mask=mask, ring=ring, cells=int(mask.sum())))
        if len(out) == HOLLOWS_FOLLOWED:
            break
    return out


def hollow_now(h, hollow, g):
    """A followed hollow on some later ground: its depth below its ring and its roundness at half depth."""
    ground = h.reshape(g['shape'])
    rim = float(np.median(ground[hollow['ring']]))
    floor = float(np.percentile(ground[hollow['mask']], 5))
    depth = rim - floor
    if depth <= 0:
        return dict(depth_m=depth, roundness=None)
    low = hollow['mask'] & (ground < floor + 0.5 * depth)
    lab, count = ndimage.label(low, structure=np.ones((3, 3), bool))
    if count > 1:
        sizes = ndimage.sum(low, lab, index=np.arange(1, count + 1))
        low = lab == int(np.argmax(sizes)) + 1
    return dict(depth_m=depth, roundness=roundness(low, g['east_west'], g['north_south']))


def measure(h, h0, region, followed, start_marks):
    g = grid_of(region)
    shape = g['shape']
    inside = interior(region)
    land = inside & region['eroding'].ravel()
    f0, f1 = forms(h0, shape)[inside], forms(h, shape)[inside]
    r0, r1 = local_relief(h0, shape)[land], local_relief(h, shape)[land]
    depth, labels, count = hollows_of(h, g)
    d0, _, _ = hollows_of(h0, g)
    marks = [hollow_now(h, hol, g) for hol in followed]
    depth_kept = [m['depth_m'] / s['depth_m'] for m, s in zip(marks, start_marks) if s['depth_m'] > 0]
    round_then = [s['roundness'] for m, s in zip(marks, start_marks) if m['roundness'] is not None and s['roundness'] is not None]
    round_now = [m['roundness'] for m, s in zip(marks, start_marks) if m['roundness'] is not None and s['roundness'] is not None]
    change = (h - h0)[inside]
    return dict(mean_change_m=float(change.mean()), rms_change_m=float(np.sqrt((change ** 2).mean())),
                share_of_ground_moved_over_100_m=float((np.abs(change) > 100.0).mean()),
                forms_kept=float(np.corrcoef(f0, f1)[0, 1]),
                forms_strength_over_start=float(f1.std() / f0.std()),
                relief_over_start=float(r1.mean() / r0.mean()),
                hollow_volume_over_start=float(depth[inside].sum() / max(d0[inside].sum(), 1e-9)),
                hollow_area_over_start=float((depth[inside] > UNDER_M).sum() / max((d0[inside] > UNDER_M).sum(), 1)),
                followed_hollows=dict(depth_over_start=dict(median=float(np.median(depth_kept)), least=float(np.min(depth_kept)), most=float(np.max(depth_kept))),
                                      roundness=dict(at_start=float(np.median(round_then)) if round_then else None,
                                                     now=float(np.median(round_now)) if round_now else None, counted=len(round_now))))


# ---- regions, runs and pictures ----------------------------------------------------------------------------------
def load_regions(names=None):
    with np.load(REGIONS, allow_pickle=False) as z:
        meta = json.loads(str(z['metadata']))
        out = {}
        for name, r in meta['regions'].items():
            if names and name not in names:
                continue
            out[name] = dict(ground=z[f'{name}__ground_m'], water=z[f'{name}__water_m'], sea=z[f'{name}__sea'],
                             eroding=z[f'{name}__eroding'], lowering=z[f'{name}__lowering_m_yr'],
                             east_west_m=r['east_west_m'], north_south_m=r['north_south_m'], centre=(r['lat_centre'], r['lon_centre']),
                             stats=r['stats'])
    return out


def one(job):
    name, water, case, region = job
    started = time.time()
    shots, books, rate = run(region, water, seed=20261122, **CASES[case])
    g, h0 = grid_of(region), start_of(region)
    followed = followed_hollows(region)
    start_marks = [hollow_now(h0, hol, g) for hol in followed]
    rows = {}
    for mark, h in shots.items():
        rows[f'{mark:g}'] = dict(measure(h, h0, region, followed, start_marks),
                                 rock={k: v / 1e9 for k, v in books[mark].items()})            # km3
    return name, water, case, dict(rate=rate, case=CASES[case], years=rows, seconds=round(time.time() - started),
                                   followed=dict(count=len(followed), cells=[hol['cells'] for hol in followed],
                                                 depth_m=[s['depth_m'] for s in start_marks],
                                                 roundness=[s['roundness'] for s in start_marks])), {f'{mark:g}': h.astype(np.float32) for mark, h in shots.items()}


def picture(name, region, water, shots, path, year='3e+07'):
    """The region's ground in relief shading: as it starts, left alone, and under each case at one year."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LightSource
    shape = region['ground'].shape
    h0 = start_of(region).reshape(shape)
    my = f'{float(year) / 1e6:g} million years'
    panels = [('as the map has it', h0), ('nothing returned, 10 million years', shots[(water, 'none')]['1e+07'].reshape(shape)),
              (f'nothing returned, {my}', shots[(water, 'none')][year].reshape(shape))]
    for case, label in (('worn', 'to worn ground, shared by wear'), ('healed', 'to worn ground, no more than it lost'),
                        ('worn_far', 'shared by wear, range 50 km'), ('worn_by_district', 'to ground worn by district'),
                        ('even', 'returned evenly'), ('crest', 'returned to crests')):
        panels.append((f'{label}, {my}', shots[(water, case)][year].reshape(shape)))
    light = LightSource(azdeg=315, altdeg=35)
    lo, hi = np.percentile(h0, [1, 99])
    wide = shape[1] * region['east_west_m'] / (shape[0] * region['north_south_m'])
    fig, axes = plt.subplots(3, 3, figsize=(5.6 * 3 * max(wide, 0.75) + 1.0, 17.5))
    extent = [0, shape[1] * region['east_west_m'] / 1e3, shape[0] * region['north_south_m'] / 1e3, 0]
    for ax, (title, h) in zip(axes.ravel(), panels):
        rgb = light.shade(np.clip(h, lo - (hi - lo), hi + (hi - lo)), cmap=plt.cm.gist_earth, vmin=lo, vmax=hi, blend_mode='overlay',
                          vert_exag=6.0, dx=region['east_west_m'], dy=region['north_south_m'])
        ax.imshow(rgb, extent=extent, aspect='equal')
        ax.set_title(title, fontsize=12)
        ax.set_xlabel('km')
    fig.suptitle(f'{name.replace("_", " ")}, {region["centre"][0]:g}° N {region["centre"][1]:g}° E: '
                 f'{"every hollow holds a lake and spills" if water == "wet" else "no hollow holds water"}', fontsize=15)
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)


def build(names=('korolev', 'southern_highlands'), workers=4, log=print):
    from multiprocessing import Pool
    started = time.time()
    regions = load_regions(names)
    small = dict(ground=np.zeros((8, 8), np.float32), water=np.full((8, 8), np.nan, np.float32), sea=np.zeros((8, 8), bool),
                 eroding=np.ones((8, 8), bool), lowering=np.zeros((8, 8), np.float32), east_west_m=1.0, north_south_m=1.0)
    for water in WATERS:                                      # compile the kernels once, before the workers start
        step(start_of(small) + np.arange(64.0), grid_of(small), water, 'even', 1.0e-7, 1.0)
    jobs = [(name, water, case, regions[name]) for name in names for water in WATERS for case in CASES]
    results, shots = {}, {}
    with Pool(workers) as pool:
        for name, water, case, row, snaps in pool.imap_unordered(one, jobs):
            results.setdefault(name, {}).setdefault(water, {})[case] = row
            shots.setdefault(name, {})[(water, case)] = snaps
            log(f'{name} {water} {case}: {row["seconds"]} s, forms kept {row["years"]["3e+07"]["forms_kept"]:.3f} at 30 million years, {time.time() - started:.0f} s')
    pictures = {}
    arrays = {}
    for name in names:
        for water in WATERS:
            path = DRIVE / f'{name}_{water}.png'
            picture(name, regions[name], water, shots[name], path)
            pictures[f'{name}, {water}'] = str(path)
            for case in CASES:
                for mark, h in shots[name][(water, case)].items():
                    arrays[f'{name}__{water}__{case}__{mark}'] = h.reshape(regions[name]['ground'].shape)
    np.savez_compressed(SNAPSHOTS, **arrays)
    placement = json.loads(PLACEMENT.read_text())
    return dict(
        schema=SCHEMA,
        producer=dict(module='research.studies.rock_return.evolve', code_sha256=sha256(Path(__file__)), compiled=COMPILED),
        inputs={'research/studies/rock_return/results/placement.json': sha256(PLACEMENT)},
        ignored_inputs={str(REGIONS): placement['product']['sha256']},
        snapshots=dict(path=str(SNAPSHOTS), sha256=sha256(SNAPSHOTS)), pictures=pictures,
        evidence=('A landscape model on squares of the design surface at 1.9 km: the stream-power law for cutting, '
                  'scaled to the rock module\'s lowering, every hollow a trap that fills from the bottom, and rules of '
                  'placement that are design guesses. No hillslope process is resolved at this grain, rain is even, '
                  'the rock has one strength, and the region\'s edge is open. It shows what kinds of change each rule '
                  'makes and how large; it predicts no place. Not validated.'),
        reading_rule=('regions: for each region, each case of water (wet: every hollow holds a lake and spills; dry: '
                      'none holds water) and each case (a rule, with the longer range or the smaller share where the case names one), at each year: mean_change_m and rms_change_m of the ground '
                      'in the measured part; forms_kept, the correlation between the starting and the later ground '
                      'once each is taken against its mean over about 30 km (1 when the forms are unchanged); '
                      'forms_strength_over_start, the spread of those forms over their starting spread; '
                      'relief_over_start, the mean height range within about 5 km over its start; hollow_volume and '
                      'hollow_area over their starts; followed_hollows, for the largest hollows of the starting map, '
                      'their depth below their rims over its start and the roundness of their outline at half depth '
                      '(0 for a disc); rock, km3 cut, laid in hollows, set down by carriers and gone from the region.'),
        setting=dict(area_exponent=AREA_EXPONENT, step_yr=STEP_YR, years=list(YEARS), range_m=RANGE_M, far_range_m=FAR_RANGE_M, district_m=DISTRICT_M,
                     surroundings_m=SURROUNDINGS_M, floor_roughness_m=FLOOR_ROUGHNESS_M, margin_cells=MARGIN, form_cells=FORM_CELLS,
                     hollows_followed=HOLLOWS_FOLLOWED, cases=CASES, waters=list(WATERS),
                     regions={name: dict(centre=regions[name]['centre'], east_west_m=regions[name]['east_west_m'],
                                         north_south_m=regions[name]['north_south_m'], **regions[name]['stats']) for name in names}),
        regions=results, seconds=round(time.time() - started))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--workers', type=int, default=4)
    args = p.parse_args(argv)
    out = build(workers=args.workers, log=lambda *a: print(*a, flush=True))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(f'wrote {OUT} in {out["seconds"]} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
