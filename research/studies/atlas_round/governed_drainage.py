"""The atlas round's second step, governed: the drainage on the redrawn floor with the works chosen on 16 November.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.governed_drainage     # about 15 minutes

The register's rows of 16 November divide Moscoviense's water at two weirs, feed Humboldtianum from the cut through the
ridge on Moscoviense's northern shore, and hold each upper sea at its settled level (research/decisions.md, the atlas
round's frame). This module routes the design years' rain over the redrawn floor with those works in place:

- Moscoviense, Orientale and Smythii-Marginis stand full at their settled rims, where the sill and the two shallow cuts
  hold them, and Orientale and Smythii-Marginis pass their whole surplus on over their passes;
- the weir at Moscoviense's natural pass sends the salt south sea what holds it at the common level, and the rest of
  Moscoviense's surplus leaves through the northern cut;
- past the ridge, Humboldtianum's share runs on by itself and is set so that Humboldtianum receives what holds it at the
  common level after the lakes it fills on the way take their share; the remainder is taken out of the lakes from which
  the two further cuts lead and put in past each cut, so that it reaches Smythii-Marginis, which passes it on with its
  own surplus to the near-side ocean;
- the near-side ocean is routed at the common level, and the level its balance would bring it to is reported, since
  whether it may part is decided on the climate run on the final surface.

The cuts' routes are redrawn_drainage's: from Moscoviense to each of the two seas, the route whose deepest cut is
shallowest and, among those, the one whose cuts add up to the least. Each weir is set by routing again, the share it
releases corrected by what its sea lacks or has to spare, until every held sea is within a cubic metre a second of its
need.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from geography import drainage as dr
from research.studies.atlas_round import redrawn_drainage as RD
from research.studies.drainage_network import network as NW

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'governed_drainage.json'
BULK = RD.BULK
GOVERNED_ROUTING = BULK / 'governed_routing.npz'
SCHEMA = 'terluna.research.atlas-round-governed-drainage/1'
PRODUCT_SCHEMA = 'terluna.research.atlas-round-governed-routing/1'
CODE = ('research/studies/atlas_round/governed_drainage.py', 'research/studies/atlas_round/redrawn_drainage.py',
        'research/studies/atlas_round/redrawn_surface.py', 'research/studies/drainage_network/network.py',
        'research/studies/drainage_network/routing.py', 'geography/drainage.py')
SOURCE = 'moscoviense'
SALT_SEA = 'south_pole_aitken'
HELD_BRANCH = 'humboldtianum'
REST_BRANCH = 'smythii_marginis'
FULL_SEAS = ('orientale', 'smythii_marginis')
SETTLE_LIMIT_M = 1.0                    # each weir is set until its sea's need is met within this many m3/s
ROUNDS = 10
BANDS = (0.0, 15.0, 30.0, 45.0, 60.0, 90.0)


# ---------------------------------------------------------------- a branch's take-outs and put-ins

def branch_operations(chain, search, sinks, in_sea, names, water, counted=RD.CUT_COUNTED_M):
    """What a branch of unit flow does to the routing along a CutSearch chain of (edge, water level before it) pairs:
    the water is put in where it falls out of a cut into a hollow it then runs on from by itself, and taken out of the
    lake a cut leaves from, where that lake is not the source sea. `water` marks the cells under a sea: a cut that opens
    onto a sea's water hands it the flow directly, while one that opens onto land draining to a sea puts the water in
    there to run down to it. Returns the put-ins and take-outs as (cell, +1 or -1,
    the lake's level for a take-out), the sea a cut opens straight onto (or None), the cuts deeper than `counted` as
    (edge, depth), and the cells where the water runs on by itself."""
    lev_l, wt_l, src, dst, c_to = search.lev_l, search.wt_l, search.src, search.dst, search.c_to
    ops, cuts, free, onto = [], [], [], None
    channel = False
    for e, lev in chain:
        x, y = int(src[e]), int(dst[e])
        depth = wt_l[e] - lev
        if depth >= counted:
            cuts.append((e, depth))
            if not channel:
                if in_sea[x] < 0:                    # the water leaves a lake through this cut: take it out there
                    ops.append((int(sinks[x]), -1, lev_l[x]))
                channel = True
        ly = lev if lev < lev_l[y] else lev_l[y]
        if channel and ly == lev_l[y]:               # it falls into a hollow and runs on by itself from there
            channel = False
            if water[int(c_to[e])]:
                onto = names[in_sea[y]]
            else:
                ops.append((int(c_to[e]), +1, None))
                free.append(int(c_to[e]))
    return ops, onto, cuts, free


def lake_on_the_way(r, h_flat, start, level, tolerance=0.01):
    """The first cell down the routing from `start` that lies under a full hollow's water standing at `level`: where the
    branch's water, run on from `start`, passes through the lake its next cut leaves from. None if it reaches a sea
    first."""
    for cell in RD.walk(r, start):
        if r['outlet'][cell]:
            return None
        lev = r['level'][r['label'][cell]]
        if h_flat[cell] < lev and abs(lev - level) < tolerance:
            return int(cell)
    return None


# ---------------------------------------------------------------- the governed routing

def route_with(geo, runoff, lake_net, waters, names, ops, area, lake_full):
    """The routing with flows (m3/s) put in at cells (on the ground as runoff, into a full hollow's balance where the
    cell lies under its water) and taken out of full hollows; returns the routing and each sea's supply."""
    h, lat, cell_km = geo['h'], geo['lat'], geo['cell_km']
    cols = h.shape[1]
    add_runoff = np.asarray(runoff, np.float64).copy()
    add_lake = np.asarray(lake_net, np.float64).copy()
    for cell, q in ops:
        if q == 0.0:
            continue
        i, j = divmod(int(cell), cols)
        if lake_full[cell]:
            add_lake[i, j] += q / (dr.MM_DAY * area[cell])
        elif q > 0.0:
            add_runoff[i, j] += q / (dr.MM_DAY * area[cell])
        else:
            raise ValueError(f'a take-out at cell {cell} lies under no hollow\'s water')
    r, got = NW.route_at(h, lat, cell_km, add_runoff, add_lake, waters)
    return r, dict(zip(names, got))


def govern(geo, state, runoff, lake_net, branches, source=SOURCE, salt_sea=SALT_SEA, held=HELD_BRANCH,
           rest=REST_BRANCH, full_seas=FULL_SEAS, rounds=ROUNDS, settle=SETTLE_LIMIT_M, log=print):
    """Set the weirs: the source's natural pass sends `salt_sea` its need, the branch to `held` (if any) carries what
    `held` needs after the branch's losses, the branch to `rest` the remainder, and the full seas pass their surplus on
    over their natural passes. `branches` maps a target sea to (operations per unit flow, the sea a cut opens straight
    onto). Returns the routing, the supplies and the shares."""
    started = time.time()
    names, levels, masks, loss, own, entry = (state[k] for k in ('names', 'levels', 'masks', 'loss', 'own', 'entry'))
    h = geo['h']
    rows, cols = h.shape
    area = np.repeat((geo['cell_km'] * 1e3) ** 2 * np.maximum(np.cos(np.radians(geo['lat'])), 1e-6), cols)
    waters = [(masks[k], levels[k]) for k in names]
    r0 = state['r']
    lake_full = ~r0['outlet'] & (h.ravel() < r0['level'][r0['label']])
    surplus = own[source] - loss[source]
    held_seas = [salt_sea] + ([held] if held else [])
    share = {k: max(loss[k] - own[k], 0.0) for k in held_seas}
    over = {k: max(own[k] - loss[k], 0.0) for k in full_seas}
    history = []
    for _ in range(rounds):
        if held:
            share[held] = min(share[held], max(surplus - share[salt_sea], 0.0))
        share[rest] = max(surplus - sum(share[k] for k in held_seas), 0.0)
        ops = [(entry[source], share[salt_sea])] + [(entry[k], over[k]) for k in full_seas]
        direct = {}
        for target in ([held] if held else []) + [rest]:
            unit_ops, onto = branches[target]
            ops += [(c, sign * share[target]) for c, sign in unit_ops]
            if onto is not None:
                direct[onto] = direct.get(onto, 0.0) + share[target]
        r, supply = route_with(geo, runoff, lake_net, waters, names, ops, area, lake_full)
        for k, q in direct.items():
            supply[k] += q
        lack = {k: loss[k] - supply[k] for k in held_seas}
        for k in held_seas:                          # a held sea that needs nothing keeps its weir shut
            if share[k] <= 0.0 and lack[k] < 0.0:
                lack[k] = 0.0
        new_over = {k: max(supply[k] - loss[k], 0.0) for k in full_seas}
        moved = max([abs(v) for v in lack.values()] + [abs(new_over[k] - over[k]) for k in full_seas])
        history.append(dict(shares_m3s={k: round(v, 1) for k, v in share.items()},
                            lacking_m3s={k: round(v, 1) for k, v in lack.items()}))
        if moved < settle:
            break
        for k in held_seas:
            share[k] = max(share[k] + lack[k], 0.0)
        over = new_over
        del r
    else:
        raise RuntimeError(f'the weirs did not settle: {history[-3:]}')
    log(f'weirs set in {len(history)} routings: {history[-1]}, {time.time() - started:.0f} s')
    return dict(r=r, supply=supply, share=share, over=over, history=history, direct=direct, surplus=surplus)


# ---------------------------------------------------------------- the record

def near_side_level(geo, supply, loss_field, area, rim, sea='near_side'):
    """The level the near-side ocean's balance would bring it to on a supply: where its loss to the air meets the
    supply, up to its rim (network.Table on its joining levels)."""
    level0 = geo['level0']
    table = NW.Table(geo['joined'][sea], area, loss_field, geo['h'].ravel(), rim + 1.0)
    held = table.level_for(supply, rim - NW.BELOW_RIM_M)
    a0 = table.at(level0)
    if held is None:
        a, l_, v = table.at(rim - NW.BELOW_RIM_M)
        return dict(reaches_its_rim=True, level_m=round(rim - level0, 1), area_km2=round(a / 1e6),
                    area_at_the_common_level_km2=round(a0[0] / 1e6))
    a, l_, v = table.at(held)
    return dict(reaches_its_rim=False, level_m=round(held - level0, 1), area_km2=round(a / 1e6),
                area_at_the_common_level_km2=round(a0[0] / 1e6))


def new_lakes(r_design, r_governed, area, cols, rotation):
    """The lakes the governed routing fills that the design state's left dry or closed: their area and bodies, by band
    of turned latitude."""
    from scipy import ndimage
    from climate.gcm.turned import to_turned
    fresh = r_governed['lake'] & ~r_design['lake'] & ~r_governed['outlet']
    body, n = ndimage.label(fresh.reshape(-1, cols), structure=np.ones((3, 3), bool))
    cells = np.flatnonzero(fresh)
    t_lat, _ = to_turned(90.0 - (cells // cols + 0.5) / 16.0, (cells % cols + 0.5) / 16.0, rotation)
    bands = {}
    for lo, hi in zip(BANDS[:-1], BANDS[1:]):
        m = (np.abs(t_lat) >= lo) & (np.abs(t_lat) < hi)
        if m.any():
            bands[f'{lo:g}-{hi:g}'] = round(float(area[cells[m]].sum()) / 1e6)
    return dict(km2=round(float(area[cells].sum()) / 1e6), bodies=int(n), km2_by_turned_latitude=bands)


def governed_record(geo, state, gov, branches_info, rotation, loss_field, area):
    """Each sea's balance under the works, the shares, the near-side ocean's level, the new rivers' long profiles and
    the lakes they fill."""
    r, level0, h, lat, cell_km = gov['r'], geo['level0'], geo['h'], geo['lat'], geo['cell_km']
    cols = h.shape[1]
    names, levels, loss, own = state['names'], state['levels'], state['loss'], state['own']
    surface = np.where(r['lake'], np.nan_to_num(r['lake_level'], nan=0.0), h.ravel())
    profile = lambda path: NW.long_profile(path, surface, r['discharge'], lat, cols, cell_km, level0, r['lake'])
    seas = {}
    for k in names:
        seas[k] = dict(level_m=round(levels[k] - level0, 2), area_km2=round(state['area'][k] / 1e6),
                       loss_m3s=round(loss[k]), own_rivers_m3s=round(own[k]),
                       supply_m3s=round(gov['supply'][k]), balance_m3s=round(gov['supply'][k] - loss[k]))
    rim_ns = min(geo['passes'][('near_side', j)] for j in geo['seas'] if j != 'near_side')
    seas['near_side']['would_settle'] = near_side_level(geo, gov['supply']['near_side'], loss_field, area, rim_ns)
    rivers = {}
    for k in (SOURCE,) + FULL_SEAS:
        path = RD.walk(r, state['entry'][k])
        rivers[f'{k}_over_its_natural_pass'] = dict(leaves=NW.place(state['entry'][k], cols, rotation), **profile(path))
    for target, info in branches_info.items():
        segments = []
        for cell in info['free']:                    # each stretch the branch's water runs by itself, to its sea or
            path = RD.walk(r, cell)                  # to the lake it is taken out of for the next cut
            out_at = [i for i, c in enumerate(path.tolist()) if c in set(info['take_outs'])]
            if out_at:
                path = path[:out_at[0] + 1]
                ends = 'the lake its next cut leaves from'
            else:
                end = state['sea_of'][path[-1]]
                ends = names[end] if end >= 0 else None
            segments.append(dict(starts=NW.place(cell, cols, rotation), reaches=ends,
                                 ends_at=NW.place(int(path[-1]), cols, rotation), **profile(path)))
        rivers[f'{SOURCE}_to_{target}'] = dict(cuts=info['cuts'], segments=segments, opens_onto=info['onto'])
    sent = {SALT_SEA: gov['share'][SALT_SEA], HELD_BRANCH: gov['share'][HELD_BRANCH], REST_BRANCH: gov['share'][REST_BRANCH]}
    return dict(
        moscoviense_surplus_m3s=round(gov['surplus']),
        released_m3s={k: round(v) for k, v in sent.items()},
        arrived_m3s={SALT_SEA: round(gov['supply'][SALT_SEA] - own[SALT_SEA]),
                     HELD_BRANCH: round(gov['supply'][HELD_BRANCH] - own[HELD_BRANCH]),
                     REST_BRANCH: round(gov['supply'][REST_BRANCH] - own[REST_BRANCH])},
        full_seas_pass_on_m3s={k: round(v) for k, v in gov['over'].items()},
        seas=seas, rivers=rivers, weirs_set_in_routings=len(gov['history']), weir_history=gov['history'])


def product(gov, state, geo, path=GOVERNED_ROUTING):
    """The governed routing on the routing's own budget, in the design product's form."""
    r = gov['r']
    lake = r['lake'] & ~r['outlet']
    cells = np.flatnonzero(lake).astype(np.int32)
    meta = dict(schema=PRODUCT_SCHEMA, common_level_m=geo['level0'], seas=state['names'],
                levels_m={k: state['levels'][k] for k in state['names']},
                grid=dict(pixels_per_degree=dr.PPD, rows=geo['h'].shape[0], cols=geo['h'].shape[1]),
                shares_m3s={k: round(v, 1) for k, v in gov['share'].items()},
                reading_rule=('As design_routing.npz: cells index rows from 90 N and columns from 0 E, 16 per degree, on '
                              'today\'s Moon; receiver, discharge_m3s and drained_km2 per cell, sea the index into seas, '
                              'lake_cell with lake_level_m and lake_spills. Here the weirs and cuts of 16 November divide '
                              'Moscoviense\'s water and the branches\' water runs in the rivers it joins.'))
    BULK.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, metadata=json.dumps(meta), receiver=r['receiver'].astype(np.int32),
                        discharge_m3s=r['discharge'].astype(np.float32),
                        drained_km2=(r['drained_m2'] / 1e6).astype(np.float32), sea=state['sea_of'],
                        lake_cell=cells, lake_level_m=r['lake_level'][cells].astype(np.float32),
                        lake_spills=~r['closed'][r['label'][cells]])
    return dict(path=str(path), sha256=RD.sha256(path))


def branches_for(geo, state, rotation, log=print):
    """The two branches' chains from Moscoviense (CutSearch) and their operations per unit flow, each take-out moved to
    the first cell under its lake's water."""
    started = time.time()
    levels, in_sea = RD.node_levels(geo, state)
    search = RD.CutSearch(geo, levels, in_sea, state['names'], SOURCE)
    allowed = search.shallowest(only=(HELD_BRANCH, REST_BRANCH))
    _, sinks, _ = geo['graph']
    r, h_flat = state['r'], geo['h'].ravel()
    cols = geo['h'].shape[1]
    branches, info = {}, {}
    for target in (HELD_BRANCH, REST_BRANCH):
        chain, tidy = search.route(target, allowed[target])
        ops, onto, cuts, free = branch_operations(chain, search, sinks, in_sea, state['names'], state['sea_of'] >= 0)
        placed, last_in = [], None
        for cell, sign, level in ops:
            if sign < 0:
                lake = None if last_in is None else lake_on_the_way(r, h_flat, last_in, level)
                if lake is None:
                    raise RuntimeError(f'the branch to {target} does not pass the lake at {level - geo["level0"]:.1f} m '
                                       f'its cut leaves from')
                cell = lake
            else:
                last_in = cell
            placed.append((cell, sign))
        branches[target] = (placed, onto)
        info[target] = dict(onto=onto, free=free, least_in_sum=tidy, take_outs=[c for c, sign in placed if sign < 0],
                            cuts=[dict(depth_m=round(d, 1), at=NW.place(int(search.c_to[e]), cols, rotation),
                                       pass_above_the_common_level_m=round(search.wt_l[e] - geo['level0'], 1))
                                  for e, d in cuts])
    log(f'branches: {dict((k, (len(v[0]), v[1])) for k, v in branches.items())}, {time.time() - started:.0f} s')
    return branches, info


def build(log=print):
    started = time.time()
    surf = RD.load_surface(log)
    _, _, _, _, _, runoff, lake_net, _, _ = NW.load(1.0)
    geo = RD.geometry(surf, runoff, lake_net, log)
    del runoff, lake_net
    rims = RD.rims_against_the_first_step(geo, surf['case'])
    if not all(v['agree'] for v in rims.values()):
        raise SystemExit('the rims do not reproduce the first step\'s design case')
    rotation = np.array(json.loads((ROOT / NW.FRAME).read_text())['frame']['rotation_today_to_turned'])
    cols = geo['h'].shape[1]
    area = np.repeat((geo['cell_km'] * 1e3) ** 2 * np.cos(np.radians(geo['lat'])), cols)
    drained = json.loads(RD.OUT.read_text())
    factor = drained['cases'][0]['budget']['open_water_evaporation_factor_closing_it']
    out, prod = {}, None
    branches = info = None
    for label, evaporation in (('routing_budget', 1.0), ('closed_budget', factor)):
        _, _, _, _, _, runoff, lake_net, _, _ = NW.load(1.0, evaporation)
        state = RD.design_state(geo, runoff, lake_net, log=log)
        if branches is None:                          # the cuts' routes stand on the floor, not on the budget
            branches, info = branches_for(geo, state, rotation, log)
        gov = govern(geo, state, runoff, lake_net, branches, log=log)
        loss_field = -lake_net.ravel().astype(np.float64) * dr.MM_DAY * area
        record = governed_record(geo, state, gov, info, rotation, loss_field, area)
        record['evaporation'] = evaporation
        record['lakes_the_works_fill'] = new_lakes(state['r'], gov['r'], area, cols, rotation)
        if label == 'routing_budget':
            prod = product(gov, state, geo)
        out[label] = record
        del state, gov
    inputs = [RD.OUT, ROOT / RD.SURFACE, ROOT / NW.FRAME]
    rel = lambda p: str(Path(p).relative_to(ROOT)) if str(p).startswith(str(ROOT)) else str(p)
    return dict(
        schema=SCHEMA, producer=dict(files={f: RD.sha256(ROOT / f) for f in CODE}),
        inputs={rel(p): RD.sha256(p) for p in inputs},
        ignored_inputs={rel(ROOT / NW.CLIMATOLOGY): RD.sha256(ROOT / NW.CLIMATOLOGY),
                        str(surf['grids']): surf['case']['grids']['sha256']},
        climatology=dict(run=NW.RUN, years=NW.YEARS),
        evidence=('A calculation on products: the drainage network study\'s routing of the design years\' rain over the '
                  'redrawn floor (the first step\'s settled case) with the works chosen on 16 November, each weir set by '
                  'routing again until the sea it serves receives its need. The routes of the cuts are found between 1.9 '
                  'km cells. The climate is the present design climate on the turned Moon at 8 degrees. Not validated '
                  'against any observation.'),
        reading_rule=('Levels in m above the redrawn common level; flows in m3/s; areas in km2. released is what each weir '
                      'or branch lets go and arrived what reaches its sea after the lakes on the way; the salt sea and '
                      'Humboldtianum are held at the common level, so their balance is nought within the weirs\' setting. '
                      'The near-side ocean is routed at the common level; would_settle gives the level its balance would '
                      'bring it to, up to its rim, where it would join Smythii-Marginis. rivers give the long profiles of '
                      'the overflows at their governed flows and of each branch from where its water runs on by itself; '
                      'lakes_the_works_fill are lakes the design state without the works leaves dry or closed.'),
        branches={k: dict(cuts=v['cuts'], opens_onto=v['onto'], cuts_least_in_sum=v['least_in_sum']) for k, v in info.items()},
        cases=out, governed_routing=prod, seconds=round(time.time() - started))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--out', type=Path, default=OUT)
    args = p.parse_args(argv)
    record = build(log=lambda *a: print(*a, flush=True))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(record, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')
    print(f'wrote {args.out} in {record["seconds"]} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
