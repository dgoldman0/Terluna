"""The storm-water study's tables for the hollows' store in PlaSim's land scheme (stage (b) of the plan, 3 November 2026;
session-tools/agents/storm_water/stage_b_design.md on the research drive).

    climate/gcm/.venv/bin/python -m research.studies.storm_water.store_tables     # about a minute, after run.py

PlaSim's land scheme with the store (climate/gcm/land_store.py on domain/climate, which reads this product as a file)
holds one lumped store of the hollows' water per T21 land cell. For every dryland cell this writes
results/store_tables.json (schema terluna.research.storm-water-store-tables/1):

- per hollow age, the store's capacity and, at FILLS common fills of the cell's hollows (each holding the same share of
  its capacity), the water area and the seepage under each floor design, summed over the hollows with phase 2b's own
  geometry and seepage laws: geography's closed depressions on the 1.9 km grid with their area rising linearly between
  cells, and the crater bowls below the grid as paraboloids at half Pike's fresh depth; the clay seal everywhere, the
  seal below a 1 m stage with sandy banks above, and the seal ten times tighter (phase 2b's sensitivity case);
- per storm library, surface and hollow age, the shares of the slopes' runoff that reach the hollows, the rivers (the
  sea and the humid land) and the ground below the roots (the channel beds' and slopes' recharge), from phase 2a's
  routing; the rest soaks back into the slopes' soil and the channel beds above the roots;
- per surface, the share of a 30-minute step's rain that runs off against the rain's depth on its wetted patch and the
  bucket's relative wetness, from phase 1's Green and Ampt model at lunar gravity, each step a storm of its own at a
  steady rate on soil the bucket's wetness gives (the plant-available reading).

The capacity of a depression with more cells below its pass than phase 2b keeps (14,470 of the 94,858 whose floors lie in
the drylands) is its volume at the pass; phase 2b gives those depressions walls 1 km tall above their sixteenth cell,
which never spill, and holds their area above it at the sixteen cells' (phase2b.depression_curves), which the curve
here keeps.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

from research.studies.storm_water import extract as EX
from research.studies.storm_water import infiltration as IN
from research.studies.storm_water import phase2a as P2
from research.studies.storm_water import phase2b as B
from research.studies.storm_water import routing as RT
from research.studies.storm_water import run as R1
from shared.constants import JULIAN_YEAR_DAYS, SYNODIC_MONTH_DAYS
from shared.provenance import constants_used

OUT = R1.HERE / 'results' / 'store_tables.json'
SCHEMA = 'terluna.research.storm-water-store-tables/1'
CODE = ('research/studies/storm_water/store_tables.py', 'research/studies/storm_water/phase2b.py',
        'research/studies/storm_water/phase2a.py', 'research/studies/storm_water/routing.py',
        'research/studies/storm_water/infiltration.py', 'research/studies/storm_water/wetting.py',
        'research/studies/storm_water/run.py', 'research/studies/storm_water/extract.py')
FILLS = np.array([0.0, 1e-6, 1e-5, 1e-4, 3e-4, 1e-3, 3e-3, 0.01, 0.03, 0.1, 0.3, 1.0])     # shares of capacity
DEPTHS_MM = np.array([0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 50.0, 70.0, 100.0, 150.0,
                      200.0, 300.0])             # a step's rain on its wetted patch (mm in 30 minutes)
WETNESS = np.array([0.0, 0.005, 0.01, 0.02, 0.04, 0.08, 0.15, 0.3, 0.5, 0.75, 1.0])   # the bucket's water over wsmax
STEP_H = 0.5                                     # PlaSim's step on the design run (30 minutes)
INFILTRATION_DT_H = 1.0 / 60.0                   # phase 1's one-minute steps within a storm
FLOORS = ('clay_floored', 'sandy_banked', 'tight_seal')
TIGHT = dict(B.SEAL, k_earth_mm_h=B.SEAL['k_earth_mm_h'] / 10.0)    # phase 2b's sensitivity seal, ten times tighter
LIBRARIES, SURFACES, STATES = R1.BRACKETS, IN.MAIN, B.STATES
DEPTH_TAB = np.linspace(0.0, 30.0, 601)          # phase 2b's seepage table (simulate's depth_tab)


def geometry(v, s, A, V):
    """Stage and area of each hollow at volume v (phase 2b's geometry closure in simulate): the first segment's law
    A = g1 stage (a wedge or a paraboloid, area as the root of the volume), then the piecewise curve (on_curve)."""
    g1 = A[:, 1] / np.maximum(s[:, 1], 1e-9)
    stage = np.sqrt(2.0 * v / np.maximum(g1, 1e-12))
    area = g1 * stage
    big = v > V[:, 1]
    if big.any():
        stage[big], area[big] = B.on_curve(v[big], s[big], A[big], V[big])
    return stage, area


def seepage(dmean, area, clog, floor):
    """Seepage (m3/s) of each hollow with mean depth dmean (m) over its water area (m2) under a floor design, by phase
    2b's laws (simulate's seepage_rate): the seal under the mean depth over the whole area, or with sandy banks the
    seal below the clogged stage's area clog and the banks' sand above it; the tight seal as phase 2b's sensitivity."""
    if floor == 'tight_seal':
        return np.interp(dmean, DEPTH_TAB, B.seal_seepage(DEPTH_TAB, IN.MOON, seal=TIGHT)) * area
    q_seal = np.interp(dmean, DEPTH_TAB, B.seal_seepage(DEPTH_TAB, IN.MOON)) * area
    if floor == 'clay_floored':
        return q_seal
    q_bank = B.bank_seepage(IN.MOON)
    return np.where(area > clog, q_seal / np.maximum(area, 1e-30) * clog + q_bank * (area - clog), q_seal)


def hollows(inp, state):
    """Every hollow of the drylands for one hollow age, as one set: the grid's depressions (both ages share them) and
    one paraboloid per T21 cell and diameter bin standing for all its bowls, wall-fed and channel-fed alike. Returns
    the curves, capacities (m3), numbers, cells and kinds (0 depression, 1 bowls)."""
    keep = inp['floor_t21'] >= 0
    nb = inp['n_below'][keep]
    s, A, V, top = B.depression_curves(inp['hyps_z'][keep].astype(float), inp['hyps_a'][keep].astype(float), nb,
                                       inp['pass_m'][keep].astype(float), inp['vmax_m3'][keep])
    full = nb <= s.shape[1] - 1
    cap_d = np.where(full, top, np.maximum(inp['vmax_m3'][keep], V[np.arange(nb.size), np.minimum(nb, s.shape[1] - 1)]))
    bc = RT.bowl_classes(state)
    dry = inp['t21_dry_km2']
    ncell, nbin = dry.size, bc['number_km2'].size
    cells_b = np.repeat(np.arange(ncell), nbin)
    bins = np.tile(np.arange(nbin), ncell)
    number = bc['number_km2'][bins] * dry[cells_b]
    sb, Ab, Vb, topb = B.bowl_curves(bc['diameter_km'][bins] * 1000.0)
    curves = B.pad_curves([(s, A, V, cap_d), (sb, Ab, Vb, topb)], max(s.shape[1], sb.shape[1]))
    return dict(s=curves[0], A=curves[1], V=curves[2], cap=curves[3], number=np.r_[np.ones(cap_d.size), number],
                cell=np.r_[inp['floor_t21'][keep].astype(np.int64), cells_b],
                kind=np.r_[np.zeros(cap_d.size, np.int8), np.ones(number.size, np.int8)],
                wide=np.r_[~full, np.zeros(number.size, bool)])


def store_curve(h, ncell, dry_m2, fills=FILLS):
    """Per cell: the capacity (m over the grid's dryland area) and, at each common fill, the water area (share of that
    area) and the seepage (m/s over it) under each floor, summed over the hollows with their numbers."""
    clog = B.area_at_stage(np.full(h['cap'].size, B.BANK['clogged_stage_m']), h['s'], h['A'])
    w = h['number']
    cap = np.bincount(h['cell'], weights=w * h['cap'], minlength=ncell) / dry_m2
    area = np.zeros((ncell, fills.size))
    seep = {f: np.zeros((ncell, fills.size)) for f in FLOORS}
    for k, f in enumerate(fills):
        if f <= 0.0:
            continue
        v = f * h['cap']
        _, a = geometry(v, h['s'], h['A'], h['V'])
        dmean = np.where(a > 0.0, v / np.maximum(a, 1e-9), 0.0)
        area[:, k] = np.bincount(h['cell'], weights=w * a, minlength=ncell) / dry_m2
        for floor in FLOORS:
            seep[floor][:, k] = np.bincount(h['cell'], weights=w * seepage(dmean, a, clog, floor), minlength=ncell) / dry_m2
    return cap, area, seep


def interpolation_error(h, ncell, dry_m2, cap, area, seep):
    """The store's curve read between its fills as land_store.curve reads it (linear in the square root of the fill)
    against the hollows' own sum at the geometric middle of each pair of fills: the largest relative error of the
    area and of the clay-floored seepage over the cells with a store, area-weighted mean and maximum."""
    mids = np.sqrt(FILLS[1:-1] * FILLS[2:])
    _, a_true, s_true = store_curve(h, ncell, dry_m2, fills=np.r_[0.0, mids])
    root = np.sqrt(FILLS)
    out = {}
    for name, tab, true in (('area', area, a_true), ('seepage_clay', seep['clay_floored'], s_true['clay_floored'])):
        errs = []
        for m, f in enumerate(mids, 1):
            k = int(np.searchsorted(FILLS, f))
            w = (math.sqrt(f) - root[k - 1]) / (root[k] - root[k - 1])
            est = tab[:, k - 1] + w * (tab[:, k] - tab[:, k - 1])
            with np.errstate(invalid='ignore', divide='ignore'):
                errs.append(np.where(true[:, m] > 0.0, np.abs(est / true[:, m] - 1.0), 0.0))
        e = np.max(errs, axis=0)
        out[name] = dict(max=float(e[cap > 0].max()), area_weighted_mean=float((e * dry_m2)[cap > 0].sum() / dry_m2[cap > 0].sum()))
    return out


def routing_shares(inp, x):
    """Per library, surface and hollow age, each cell's shares of its slopes' runoff: phase 2a's fates over the runoff
    they came from, both mm a year over the cell (phase2b_inputs' fates and off). Where flows crossing a T21 border
    leave a cell's hollows, rivers and recharge above its own runoff, the three are scaled to it and the soil's share
    is nil; the product records how often and by how much."""
    years = x['pr'].shape[0] / 120.0
    yr = JULIAN_YEAR_DAYS / (12.0 * SYNODIC_MONTH_DAYS)
    out, checks = {}, {}
    for b in LIBRARIES:
        for name in SURFACES:
            rt = R1.by_tenth(inp[f'off__{b}__{name}'].astype(float), x['tenth']).sum(axis=0) / years * yr
            for s in STATES:
                f = {q: inp[f'{q}__{b}__{name}__{s}'].astype(float).sum(axis=0)
                     for q in ('bowls_inside', 'bowls_slope', 'bowls_channel', 'depressions', 'slope_recharge',
                               'bed_recharge', 'sea', 'humid_land')}
                with np.errstate(invalid='ignore', divide='ignore'):
                    share = lambda v: np.where(rt > 0.0, v / np.maximum(rt, 1e-30), 0.0)
                    hol = share(f['bowls_inside'] + f['bowls_slope'] + f['bowls_channel'] + f['depressions'])
                    riv = share(f['sea'] + f['humid_land'])
                    rch = share(f['slope_recharge'] + f['bed_recharge'])
                tot = hol + riv + rch
                over = tot > 1.0
                scale = np.where(over, 1.0 / np.maximum(tot, 1e-30), 1.0)
                hol, riv, rch = hol * scale, riv * scale, rch * scale
                soil = np.where(over, 0.0, np.maximum(1.0 - hol - riv - rch, 0.0))
                out[(b, name, s)] = np.stack([hol, riv, rch, soil], axis=1)
                a = x['area_km2']
                checks[f'{b}__{name}__{s}'] = dict(cells_scaled=int(over.sum()), largest_sum=float(tot.max()),
                                                   area_share_scaled=float(a[over].sum() / a.sum()),
                                                   hollows_share_drylands=float((hol * rt * a).sum() / (rt * a).sum()))
    return out, checks


def runoff_tables(pressure):
    """Per surface, the share of a 30-minute step's rain that runs off (DEPTHS_MM, WETNESS): the step's depth falling at
    a steady rate through the step onto soil at the effective saturation the bucket's relative wetness gives (phase 1's
    plant-available reading), in phase 1's one-minute Green and Ampt steps with its 1 mm surface store."""
    out = {}
    d = np.repeat(DEPTHS_MM[:, None], WETNESS.size, axis=1)
    for name in SURFACES:
        p = IN.parameters(name, pressure_pa=pressure)
        sat = IN.initial_saturation(p['soil'], WETNESS * IN.WSMAX_M)
        dtheta = np.repeat(IN.deficit(p, sat)[None, :], DEPTHS_MM.size, axis=0)
        off = IN.steady_rain(p, d, d / STEP_H, dtheta, 0.0, dt_h=INFILTRATION_DT_H)
        with np.errstate(invalid='ignore', divide='ignore'):
            out[name] = np.where(d > 0.0, off / np.maximum(d, 1e-30), 0.0)
    return out


def channel_beds(product_2a):
    """The channel beds' share of the slopes' runoff by class, library, surface and hollow age (phase 2a's class sums;
    the per-cell beds' take is not kept in phase 2b's inputs, so the store's tables give it within the soil's share)."""
    fates = product_2a['routing']['fates']
    return {b: {n: {s: {c: fates[b][n][s][c]['channel_bed']['share_of_runoff'] for c in R1.CLASSES + ('drylands',)}
                    for s in STATES} for n in SURFACES} for b in LIBRARIES}


def main():
    x, meta = EX.load()
    with np.load(P2.PHASE2B_INPUTS) as z:
        meta_in = json.loads(str(z['metadata']))
        inp = {k: z[k] for k in z.files if k != 'metadata'}
    ncell = x['cls'].size
    dry_km2 = inp['t21_dry_km2'].astype(float)
    dry_m2 = np.maximum(dry_km2, 1e-9) * 1e6
    stores, composition, interp = {}, {}, {}
    for s in STATES:
        h = hollows(inp, s)
        cap, area, seep = store_curve(h, ncell, dry_m2)
        stores[s] = (cap, area, seep)
        w = h['number'] * h['cap']
        total = np.bincount(h['cell'], weights=w, minlength=ncell)
        parts = dict(bowls=h['kind'] == 1, depressions_within_sixteen_cells=(h['kind'] == 0) & ~h['wide'],
                     depressions_wider=(h['kind'] == 0) & h['wide'])
        composition[s] = {k: float(w[m].sum() / w.sum()) for k, m in parts.items()}
        composition[s]['capacity_m_quantiles'] = [float(v) for v in np.quantile(cap, [0.1, 0.5, 0.9])]
        interp[s] = interpolation_error(h, ncell, dry_m2, cap, area, seep)
        print(f'store {s}: capacity median {np.median(cap):.2f} m, area at full {np.median(area[:, -1]):.3f}, '
              f'interpolation {interp[s]}', file=sys.stderr)
        del h, total
    shares, share_checks = routing_shares(inp, x)
    pressure = R1.mean_pressure_pa(x)
    tables = runoff_tables(pressure)
    product_2a = json.loads(R1.OUT.read_text())['phase2a']
    cells = []
    for c in range(ncell):
        row = dict(j=int(x['j'][c]), i=int(x['i'][c]), lat=R1.r3(x['lat'][c]), lon_e=R1.r3(x['lon'][c]),
                   unep=str(x['cls'][c]), area_km2=R1.r3(x['area_km2'][c]), dryland_km2=R1.r3(dry_km2[c]),
                   store={s: dict(capacity_m=float(stores[s][0][c]), area_share=[float(v) for v in stores[s][1][c]],
                                  seepage_m_s={f: [float(v) for v in stores[s][2][f][c]] for f in FLOORS})
                          for s in STATES},
                   shares={b: {n: {s: dict(zip(('hollows', 'rivers', 'recharge', 'soil'),
                                               (float(v) for v in shares[(b, n, s)][c])))
                                   for s in STATES} for n in SURFACES} for b in LIBRARIES})
        cells.append(row)
    product = dict(
        schema=SCHEMA, evidence=EVIDENCE, reading_rule=READING_RULE,
        basis=dict(run=meta['run'], years=meta['years'], cells=ncell, drylands_pressure_pa=pressure,
                   moon_gravity_in_earths=IN.MOON, step_h=STEP_H),
        grid=dict(nlon=int(x['grid_lon'].size), nlat=int(x['grid_lat'].size)),
        grids=dict(fill=[float(v) for v in FILLS], rain_depth_mm=[float(v) for v in DEPTHS_MM],
                   wetness=[float(v) for v in WETNESS]),
        floors=list(FLOORS), surfaces=list(SURFACES), libraries=list(LIBRARIES), hollows=list(STATES),
        floor_designs=dict(clay_floored=dict(seal=B.SEAL), sandy_banked=dict(seal=B.SEAL, bank=B.BANK),
                           tight_seal=dict(seal=TIGHT, note='phase 2b\'s sensitivity case, the loam only there')),
        surface_designs={n: IN.SURFACES[n] for n in SURFACES},
        runoff_tables={n: tables[n].tolist() for n in SURFACES},
        checks=dict(capacity_composition=composition, interpolation_between_fills=interp, routing_shares=share_checks,
                    channel_beds_share_of_runoff=channel_beds(product_2a)),
        cells=cells)
    files = {p: R1.sha256(p) for p in CODE}
    product['producer'] = dict(study='research/studies/storm_water', exporter='store_tables.py', files=files,
                               inputs={str(P2.PHASE2B_INPUTS.relative_to(R1.ROOT)): R1.sha256(P2.PHASE2B_INPUTS.relative_to(R1.ROOT)),
                                       str(EX.OUT.relative_to(R1.ROOT)): meta['npz']['sha256'],
                                       str(R1.OUT.relative_to(R1.ROOT)): R1.sha256(R1.OUT.relative_to(R1.ROOT))},
                               inputs_schema=meta_in['schema'], drainage_sources=meta_in['drainage_sources'],
                               constants=constants_used(CODE))
    head = json.dumps({k: v for k, v in product.items() if k != 'cells'}, indent=1, default=float)
    rows = ',\n'.join('  ' + json.dumps(c, separators=(',', ':')) for c in cells)
    OUT.write_text(head[:-2] + ',\n "cells": [\n' + rows + '\n ]\n}\n')
    return product


EVIDENCE = ('A model result of the offline storm-water chain on the design years (A28R_301k_design_l15_snow_s934, years '
            '10-19) and geography\'s 1.9 km depression tree, for PlaSim\'s hollows\' store (climate/gcm/land_store.py). '
            'The store\'s curves sum phase 2b\'s hollows at common fills: every hollow holding the same share of its '
            'capacity, which spreads a cell\'s water over its hollows by their capacity, where phase 2b\'s pools fill by '
            'their own inflow; the cost of that lumping is measured offline (session-tools/agents/storm_water/'
            'stage_b_report.md). The floors are designs (a clay seal scaled from the HAPEX-Sahel plateau pools, sandy '
            'banks above a 1 m stage, a seal ten times tighter), the bowls\' counts the terrain study\'s equilibrium law, '
            'their shape half Pike\'s fresh depth, and the routing shares phase 2a\'s routing under stand-in storms '
            '(Walnut Gulch\'s measured family and CM1\'s 6 km window over wet ground). The runoff tables read each '
            '30-minute step as a storm of its own at a steady rate on soil the bucket\'s wetness gives, Green and Ampt at '
            'lunar gravity on designed surfaces, so a storm\'s memory from one step to the next is lost. No PlaSim run '
            'has used these tables.')
READING_RULE = ('Each cell is a T21 point j * nlon + i of PlaSim\'s grid (north to south from 0 E). Depths, areas and '
                'seepage are per unit of the grid\'s dryland area in the cell (dryland_km2), the basis of phase 2a\'s '
                'fates and phase 2b\'s per-cell results, which PlaSim carries over the whole T21 cell as it carries the '
                'bucket. store[age].capacity_m is the hollows\' total volume (m of water); area_share and seepage_m_s '
                '(per floor) are the water area\'s share and the floors\' seepage (m/s) with every hollow at the share '
                'grids.fill of its capacity, read between fills linearly in the square root of the fill (exact for the '
                'bowls) and nil at nil volume. shares[library][surface][age] split the slopes\' runoff into the hollows, '
                'the rivers (the sea and the humid land), the ground below the roots (recharge) and the soil (the slopes\' '
                'and channel beds\' soak above the roots); they sum to 1. runoff_tables[surface][depth][wetness] is the '
                'share of a step\'s rain that runs off, bilinear in grids.rain_depth_mm (the step\'s rain on its wetted '
                'patch, mm in 30 minutes) and grids.wetness (the bucket\'s water over its 0.5 m capacity), held within '
                'both grids.')


if __name__ == '__main__':
    main()
