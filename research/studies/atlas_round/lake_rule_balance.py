"""The design's settings made again on the corrected lake rule, and the design state they give.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python3 -m research.studies.atlas_round.lake_rule_balance   # about 6 minutes

lake_rule.py routed the design climate once with every hollow filled as one water while the seas' levels, the weirs'
shares and the factor on the open water's evaporation stayed as the present rule had set them. Here each of those
settings is made again on the corrected rule by the procedure that made it, in the design's order:

- the factor (final_drainage.budget): with every sea at the common level and the open water's evaporation as the
  climate gives it, the rivers bring the seas more than they lose; the factor on the open water's evaporation is the
  one that takes that surplus away, one step from 1 as the design's was. More of the Moon lies under lakes on the
  corrected rule, so less is left over and the factor is smaller. What the budget still leaves over at the factor is
  routed and reported for both rules;
- the weirs (governed_drainage.govern): Moscoviense, Orientale and Smythii-Marginis stand full at their rims; the
  weir at Moscoviense's pass sends the salt south sea what holds it at the common level, the northern cut carries
  what Humboldtianum needs after the lakes on its way have taken theirs, and the rest of Moscoviense's surplus goes
  down the branch to Smythii-Marginis, which passes it on with its own surplus. The cuts stay where the design put
  them. A cut can take no more than its lake can spare while standing full: where it could not, the branch carries
  on with what the lake does spare;
- the near-side sea (near_side_neck.balance, design_routing.route): each part stands where its loss to the air meets
  its rivers, the main part sending the southern seas what the neck passes while it stands over the crest and nothing
  once it stands under it; where the southern part's own rivers and the neck's flow lift it over the crest too, the
  two are one sea across the neck. The design routed once with the parts at the levels the common-level routing gave; here
  the routing and the levels are set against each other until each part's level holds on the rivers that reach it.
  A part's level can come to rest at the sill of a basin beside it, with more reaching it than it loses below the
  sill. The part is then routed with its water just over the sill and the basin under it: where the basin was a lake
  on its rivers' way, the lake's loss had been taken from the rivers already, and the part rises on from there. Where
  the part then loses more than reaches it, it stands at the sill and spills what is left over into the basin, which
  keeps the lake that spill and its own rivers hold.
  The search for a part's level takes the cells that join its water at one level together (part_tables), where
  near_side_neck.py's runs through them one by one.

The record gives each sea's level, area, loss and supply beside the design's, the lakes, the works' flows, the named
waterways' flows where they end and at their heads, and the water the seas and lakes hold. Nothing in the design is
changed: the surface, the maps and the page still stand on the present rule.

The seas' cells at other levels than the design's take a minute to find and are kept on the research drive, as are the
factor, the weirs' shares and the parts' levels once found, so a later run takes about two minutes; --restart sets
everything again on the kept cells and --refresh finds the cells again as well.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from geography import drainage as dr
from research.studies.atlas_round import design_routing as DRt
from research.studies.atlas_round import final_drainage as FD
from research.studies.atlas_round import governed_drainage as GD
from research.studies.atlas_round import lake_rule as LR
from research.studies.atlas_round import near_side_neck as N
from research.studies.atlas_round import redrawn_drainage as RD
from research.studies.atlas_round import sustained_rivers as SR
from research.studies.drainage_network import network as NW
from research.studies.drainage_network.routing import route_levels
from research.studies.south_sea_ocean import inflows as I

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'lake_rule_balance.json'
BULK = LR.BULK
WORK = BULK / 'work'
SETTING = WORK / 'balance_setting.npz'
STAGES = WORK / 'balance_stages.json'
PRODUCT = BULK / 'lake_rule_balance_16ppd.npz'
SCHEMA = 'terluna.research.atlas-round-lake-rule-balance/1'
CODE = ('research/studies/atlas_round/lake_rule_balance.py', 'research/studies/atlas_round/lake_rule.py',
        'research/studies/atlas_round/near_side_neck.py', 'research/studies/atlas_round/governed_drainage.py',
        'research/studies/atlas_round/final_drainage.py', 'research/studies/drainage_network/routing.py',
        'geography/drainage.py')
SETTING_CODE = ('research/studies/atlas_round/near_side_neck.py', 'research/studies/atlas_round/redrawn_drainage.py',
                'research/studies/drainage_network/network.py', 'research/studies/drainage_network/routing.py',
                'research/studies/south_sea_ocean/inflows.py', 'geography/crustal_loading.py', 'geography/drainage.py')
MM_DAY = dr.MM_DAY
WATERS = SR.WATERS                   # the design routing's waters: the seas, the near-side sea as its two parts
SIX = ('moscoviense', 'orientale', 'smythii_marginis', 'south_pole_aitken', 'humboldtianum', 'near_side')
MAIN, SOUTH = WATERS.index('near_side_main'), WATERS.index('near_side_south')
SOURCE, SALT_SEA, HELD, REST, FULL_SEAS = GD.SOURCE, GD.SALT_SEA, GD.HELD_BRANCH, GD.REST_BRANCH, GD.FULL_SEAS
SETTLE_M3S = GD.SETTLE_LIMIT_M       # each weir is set until its sea's need is met within this many m3/s
ROUNDS = 30                          # routings allowed for the weirs to settle
PART_ROUNDS = 20                     # routings allowed for the near-side sea's parts to settle
PART_SETTLED_M3S = 20.0              # a part has settled when its rivers and its loss at its level agree within this
PART_SETTLED_SHARE = 0.002           # or within this share of what it receives, where its shore moves by steps
PART_STILL_M = 0.005                 # or when another routing moves neither part's level by this much
PART_STILL_M3S = 0.5                 # nor the neck's flow by this much
CONTROL_LEVEL_M = 0.05               # the neck read from its record gives the design's levels back within this
LAKES_LISTED = 8                     # the largest lakes listed, and the largest that pass nothing on
COURSE_LAKES_LISTED = 3              # the largest lakes listed on the course of a waterway that loses its flow


def sha256(path) -> str:
    return LR.sha256(path)


def rel(path) -> str:
    return str(Path(path).relative_to(ROOT)) if str(path).startswith(str(ROOT)) else str(path)


# ---------------------------------------------------------------- the seas' cells at any level

def setting(d, log=print, refresh=False):
    """What the seas' cells at other levels than the design's are found from, kept on the research drive since it
    takes the floor's catchment graph twice: the sea over each cell with every sea at the common level (an index into
    SIX, -1 elsewhere; redrawn_drainage.geometry's joining levels), and near_side_neck.setting's parts of the
    near-side sea with the level at which each cell joins each part's water, the crest between them, and the floor's
    common level."""
    key = dict(structure=sha256(SR.STRUCTURE), code={f: sha256(ROOT / f) for f in SETTING_CODE})
    if SETTING.exists() and not refresh:
        with np.load(SETTING, allow_pickle=False) as z:
            meta = json.loads(str(z['metadata']))
            if meta['key'] == key:
                log(f'the seas\' cells from {SETTING}')
                return dict(meta=meta, **{k: z[k] for k in z.files if k != 'metadata'})
    started = time.time()
    st = N.setting(log)
    h, level0, lat, cell_km = st['h'], st['level0'], st['lat'], st['cell_km']
    if not np.array_equal(h, d['h']) or abs(level0 - d['level0']) > 1e-9:
        raise SystemExit('the redrawn floor is not the design structure\'s')
    rows, cols = h.shape
    flat = st['flat']
    part = np.zeros(flat.size, np.int8)
    part[st['part']['main']] = 1
    part[st['part']['south']] = 2
    joined_main, joined_south = st['joined']['main'].astype(np.float32), st['joined']['south'].astype(np.float32)
    crest, join4, surf = float(st['crest']), float(st['join4']), st['surf']
    del st
    lon = (np.arange(cols) + 0.5) / dr.PPD
    body = I.label_waters(h < level0)
    named, _ = I.name_seas(body, surf['row_area'], lat, ((lon + 180.0) % 360.0) - 180.0, surf['atlas']['bodies'])
    label, sinks, edges = RD.catchment_graph(h, lat, cell_km)
    common = np.full(flat.size, -1, np.int8)
    for i, k in enumerate(SIX):
        cells = np.flatnonzero(body.ravel() == named[k]['label'])
        seed = int(cells[np.argmin(flat[cells])])
        common[NW.joining_level(flat, label, sinks.size, edges, seed) < level0] = i
    del body, label, sinks, edges
    meta = dict(key=key, level0=float(level0), crest=crest, join4=join4, seconds=round(time.time() - started))
    WORK.mkdir(parents=True, exist_ok=True)
    arrays = dict(common=common, part=part, joined_main=joined_main, joined_south=joined_south)
    np.savez_compressed(SETTING, metadata=json.dumps(meta), **arrays)
    log(f'the seas\' cells at any level found in {meta["seconds"]} s and kept in {SETTING}')
    return dict(meta=meta, **arrays)


def design_state_waters(d, st):
    """The water over each cell (an index into WATERS) and each water's level in the state the weirs are set on: the
    upper seas full at their rims, the salt sea and Humboldtianum at the common level, and the near-side sea at the
    common level, its cells given to the part each belongs to."""
    water = np.where(d['water'] < MAIN, d['water'], -1).astype(np.int8)
    near = st['common'] == SIX.index('near_side')
    water[near & (st['part'] == 2)] = SOUTH
    water[near & (st['part'] != 2)] = MAIN
    levels = list(d['levels'][:MAIN]) + [d['level0'], d['level0']]
    return water, levels


def parted_waters(d, st, z_main, z_south):
    """The water over each cell with the near-side sea's parts at their own levels: each part's cells that join its
    water below its level (design_routing.route's)."""
    water = np.where(d['water'] < MAIN, d['water'], -1).astype(np.int8)
    water[(st['part'] == 1) & (st['joined_main'] < z_main)] = MAIN
    water[(st['part'] == 2) & (st['joined_south'] < z_south)] = SOUTH
    return water, list(d['levels'][:MAIN]) + [float(z_main), float(z_south)]


def structure_of(d, water, levels):
    h = d['h']
    outlet = (water >= 0).reshape(h.shape)
    level = np.where(water >= 0, np.asarray(levels)[np.maximum(water, 0)], 0.0).reshape(h.shape)
    return LR.hollow_structure(h, outlet, d['lat'], d['cell_km'], level)


# ---------------------------------------------------------------- the factor that closes the budget

def budget(lake, spills, delivered, sea, runoff, lake_net, precip, evap, area):
    """final_drainage.budget's figures for a routing with every sea at the common level: the rivers' surplus over the
    seas' loss, and the factor on the open water's evaporation that takes it away in one step."""
    land = ~sea & ~lake
    mm = MM_DAY
    out = dict(
        land_rain_less_evaporation_m3s=round(float(((precip - evap).ravel() * mm * area)[land].sum())),
        routing_runoff_m3s=round(float((runoff.ravel() * mm * area)[land].sum())),
        lakes_rain_less_evaporation_m3s=round(float((lake_net.ravel() * mm * area)[lake].sum())),
        seas_loss_m3s=round(float((-lake_net.ravel().astype(np.float64) * mm * area)[sea].sum())),
        delivered_to_the_seas_m3s=round(float(delivered)),
        lakes_km2=round(float(area[lake].sum()) / 1e6),
        lakes_passing_nothing_on_km2=round(float(area[lake & ~spills].sum()) / 1e6))
    out['surplus_m3s'] = out['delivered_to_the_seas_m3s'] - out['seas_loss_m3s']
    open_water = float((NW.evap_lake_of(lake_net, precip).ravel() * mm * area)[sea | lake].sum())
    out['open_water_evaporation_m3s'] = round(open_water)
    out['open_water_evaporation_factor_closing_it'] = round(1.0 + out['surplus_m3s'] / open_water, 3)
    return out


def corrected_closing(d, st, log=print):
    """The budget with every sea at the common level on the corrected rule at factor 1, the factor it gives, and what
    the budget still leaves over at that factor. Returns the factor and the two budgets."""
    started = time.time()
    h = d['h']
    sea = st['common'] >= 0
    s = LR.hollow_structure(h, sea.reshape(h.shape), d['lat'], d['cell_km'], np.where(sea, d['level0'], 0.0).reshape(h.shape))
    out = {}
    factor = 1.0
    for key in ('at_factor_1', 'at_its_factor'):
        runoff, lake_net, precip, evap = FD.load(factor)
        terms = LR.hollow_terms(s, runoff, lake_net)
        r = LR.hollow_routing(s, terms, LR.hollow_state(s, terms))
        out[key] = budget(r['lake'], r['lake_spills'], r['discharge'][sea].sum(), sea, runoff, lake_net, precip, evap, s['area'])
        out[key]['factor'] = factor
        log(f'corrected rule, every sea at the common level, factor {factor}: {out[key]}, {time.time() - started:.0f} s')
        if key == 'at_factor_1':
            factor = out[key]['open_water_evaporation_factor_closing_it']
        del r, terms
    return factor, out


def closing_factor(d, st, design_factor, log=print):
    """corrected_closing, and for the present rule what its budget leaves over at the design's factor."""
    started = time.time()
    h = d['h']
    sea = st['common'] >= 0
    factor, out = corrected_closing(d, st, log)
    runoff, lake_net, precip, evap = FD.load(design_factor)
    outlet = sea.reshape(h.shape)
    r = route_levels(h, outlet, d['lat'], d['cell_km'], runoff, lake_net, np.where(outlet, d['level0'], 0.0))
    lake = r['lake'] & ~sea
    present = budget(lake, lake & ~r['closed'][r['label']], r['discharge'][sea].sum(), sea, runoff, lake_net, precip, evap, r['area'])
    present['factor'] = design_factor
    log(f'present rule, every sea at the common level, factor {design_factor}: {present}, {time.time() - started:.0f} s')
    return factor, dict(corrected=out, present_at_the_design_factor=present)


# ---------------------------------------------------------------- the weirs

def roots_of(tree):
    """The hollow that overflows towards a sea round each hollow."""
    n_all = tree['parent'].size
    root = np.where(tree['parent'] >= 0, tree['parent'], np.arange(n_all))
    while True:
        nxt = root[root]
        if np.array_equal(nxt, root):
            return root
        root = nxt


def branch_state(s, terms, fixed, branches, root):
    """The hollows' state with flows put in at cells (`fixed`: (cell, m3/s) pairs) and the branches' water put in and
    taken out along their cuts. branches: (cells and signs along the branch, the share released into it) for each. A
    cut takes the branch's water out of the lake it leaves from, which must stand full for the cut to draw at its
    level: where the hollow round a take-out would not stand full, that take-out and all the branch carries beyond it
    are cut by what the hollow lacks. Returns the state, each branch's flow at each of its cells, and the sources."""
    label, need = s['label'], terms['need']
    flows = [[float(q)] * len(ops) for ops, q in branches]
    while True:
        sources = list(fixed)
        for (ops, _), qs in zip(branches, flows):
            sources += [(c, sign * q) for (c, sign), q in zip(ops, qs)]
        state = LR.hollow_state(s, terms, sources)
        lacks = np.bincount(root, weights=need - state['held'], minlength=root.size)
        cut = None
        for (ops, _), qs in zip(branches, flows):
            for j, (c, sign) in enumerate(ops):
                top = int(root[label[int(c)]])
                if sign < 0 and qs[j] > 0.0 and not state['full'][top] and lacks[top] > 1e-9:
                    less = min(float(lacks[top]) * (1.0 + 1e-9) + 1e-9, qs[j])   # a hair over, so the lake stands full
                    for jj in range(j, len(ops)):
                        qs[jj] = max(qs[jj] - less, 0.0)
                    cut = (c, less)
                    break
            if cut:
                break
        if cut is None:
            return state, flows, sources


def set_weirs(s, terms, water, names, works, start=None, rounds=ROUNDS, within=SETTLE_M3S, log=print):
    """governed_drainage.govern on the corrected rule: the weir at the source's pass sends the salt sea what holds it
    at its level, the branch to the held sea carries what that sea needs after the lakes on its way, the branch to the
    rest's sea the remainder of the source's surplus, and each full sea passes its surplus on at its pass. works:
    source_pass (cell), full_pass (cell by sea), held_branch and rest_branch ((cell, sign) along each). start: shares
    and overflows to begin from. The source gives no more than its surplus: a held sea whose need it cannot meet is
    left short by `unmet` (govern itself does not settle then). Returns the state, each sea's loss, own rivers and
    supply, the shares, the full seas' overflows, each branch's flow at its cells and the sources."""
    started = time.time()
    idx = {k: i for i, k in enumerate(names)}
    count = len(names)
    at_sea = water >= 0
    loss_by = -np.bincount(water[at_sea], weights=terms['net'][at_sea], minlength=count)
    own_by = LR.reaches(s, terms, LR.hollow_state(s, terms), water, count)
    loss, own = (lambda k: float(loss_by[idx[k]])), (lambda k: float(own_by[idx[k]]))
    root = roots_of(s['tree'])
    surplus = own(SOURCE) - loss(SOURCE)
    held_seas = (SALT_SEA, HELD)
    share = {k: max(loss(k) - own(k), 0.0) for k in held_seas}
    over = {k: max(own(k) - loss(k), 0.0) for k in FULL_SEAS}
    if start:
        share.update({k: float(start['share'][k]) for k in held_seas})
        over.update({k: float(start['over'][k]) for k in FULL_SEAS})
    history = []
    for _ in range(rounds):
        share[SALT_SEA] = min(share[SALT_SEA], max(surplus, 0.0))
        share[HELD] = min(share[HELD], max(surplus - share[SALT_SEA], 0.0))
        share[REST] = max(surplus - sum(share[k] for k in held_seas), 0.0)
        fixed = [(works['source_pass'], share[SALT_SEA])] + [(works['full_pass'][k], over[k]) for k in FULL_SEAS]
        branches = [(works['held_branch'], share[HELD]), (works['rest_branch'], share[REST])]
        state, flows, sources = branch_state(s, terms, fixed, branches, root)
        supply_by = LR.reaches(s, terms, state, water, count)
        lack = {k: loss(k) - float(supply_by[idx[k]]) for k in held_seas}
        room = {SALT_SEA: surplus - share[SALT_SEA], HELD: share[REST]}     # what the source could still send each
        unmet = {}
        for k in held_seas:
            if share[k] <= 0.0 and lack[k] < 0.0:    # a held sea that needs nothing keeps its weir shut
                lack[k] = 0.0
            if lack[k] > room[k] and room[k] <= within:      # the source has no more to give it
                unmet[k], lack[k] = lack[k], 0.0
        new_over = {k: max(float(supply_by[idx[k]]) - loss(k), 0.0) for k in FULL_SEAS}
        moved = max([abs(v) for v in lack.values()] + [abs(new_over[k] - over[k]) for k in FULL_SEAS])
        history.append(dict(shares_m3s={k: round(v, 1) for k, v in share.items()},
                            lacking_m3s={k: round(v, 1) for k, v in lack.items()},
                            unmet_m3s={k: round(v, 1) for k, v in unmet.items()}))
        if moved < within:
            break
        for k in held_seas:
            share[k] = max(share[k] + lack[k], 0.0)
        over = new_over
    else:
        raise RuntimeError(f'the weirs did not settle: {history[-3:]}')
    log(f'weirs set in {len(history)} settlings: {history[-1]}, {time.time() - started:.0f} s')
    return dict(state=state, loss=loss_by, own=own_by, supply=supply_by, share=share, over=over, surplus=surplus, unmet=unmet,
                flows=dict(held_branch=flows[0], rest_branch=flows[1]), sources=sources, history=history, root=root)


def govern(d, st, runoff, lake_net, start=None, log=print, routing=False):
    """The weirs set on a climate with the near-side sea at the common level, the upper seas full at their rims and
    the works where the design structure has them. Returns set_weirs' result without its state, what the gathered
    rivers bring each water, the works' cells, and the put-ins the settled weirs make: `fixed` at cells and the
    branches with their shares; with routing, also the waters and their levels and the routing's structure, terms,
    state and result."""
    names = list(WATERS)
    water, levels = design_state_waters(d, st)
    s = structure_of(d, water, levels)
    terms = LR.hollow_terms(s, runoff, lake_net)
    ops = list(zip(d['ops'], d['op_kind']))
    works = dict(source_pass=next(c for (c, _, _), k in ops if k == 'salt_sea_weir'),
                 full_pass={name: next(c for (c, _, _), k in ops if k == name) for name in FULL_SEAS},
                 held_branch=[(c, float(np.sign(q))) for (c, q, _), k in ops if k == 'humboldtianum_branch'],
                 rest_branch=[(c, float(np.sign(q))) for (c, q, _), k in ops if k == 'smythii_marginis_branch'])
    gov = set_weirs(s, terms, water, names, works, start=start, log=log)
    r = LR.hollow_routing(s, terms, gov['state'])
    supply = np.array([float(r['discharge'][water == i].sum()) for i in range(len(names))])
    if np.abs(supply - gov['supply']).max() > 1e-3:
        raise SystemExit('the rivers gathered do not bring the seas what the hollows\' overflows give them')
    least = float(r['left'].min())
    fixed = [(works['source_pass'], gov['share'][SALT_SEA])] + [(works['full_pass'][k], gov['over'][k]) for k in FULL_SEAS]
    branches = [(works['held_branch'], gov['share'][HELD]), (works['rest_branch'], gov['share'][REST])]
    out = dict(gov={k: v for k, v in gov.items() if k not in ('state', 'root')}, supply=supply, works=works,
               fixed=fixed, branches=branches, least_kept_by_a_lake_not_full=least)
    if routing:
        out.update(water=water, levels=levels, s=s, terms=terms, state=gov['state'], r=r)
    return out


# ---------------------------------------------------------------- the near-side sea's parts

def neck_from_record(st, record):
    """near_side_neck.Neck from the neck's record: the throat's section and width at each level it was cut at, and the
    sill's length. The record rounds the widths to 0.1 km; the control in build() sets what that changes."""
    level0, crest = st['meta']['level0'], st['meta']['join4']
    found = record['throat']['by_level']
    return N.Neck(crest, [crest] + [level0 + t['level_m'] for t in found], [0.0] + [float(t['section_m2']) for t in found],
                  [0.0] + [t['width_km'] * 1e3 for t in found], record['sill']['length_km'] * 1e3)


def part_tables(d, st, loss, area):
    """near_side_neck.tables on the kept parts: each part's area, loss and volume against its level. The level search
    (near_side_neck.level_for) looks for the first level at which a part's loss reaches its supply along the running
    greatest loss, which it takes cell by cell. The cells of a basin beyond a sill all join the part's water at the
    sill's level and come under it together, some losing to the air and some gaining rain, so a search cell by cell
    can stop part of the way through a basin, at a loss the part has at no level. Here the running greatest loss is
    taken level by level: the level found is a sill itself, where the basin beyond it would bring the loss over the
    supply, or a level past the basin."""
    like = dict(flat=d['h'].ravel(), level0=d['level0'], area=area,
                part=dict(main=st['part'] == 1, south=st['part'] == 2),
                joined=dict(main=st['joined_main'], south=st['joined_south']))
    tab = N.tables(like, loss)
    for t in tab.values():
        by_level(t)
    return tab


def by_level(table):
    """Set a network.Table's running greatest loss, which near_side_neck.level_for searches, level by level."""
    if table.c.size:
        last = np.flatnonzero(np.r_[table.c[1:] != table.c[:-1], True])      # the last cell of each level
        table.runmax = np.maximum.accumulate(table.loss[last])[np.searchsorted(last, np.arange(table.c.size))]
    return table


def basin_at_the_sill(st, hf, which, level):
    """The cells of the basin a part's water reaches into when it stands at `level`, a sill: the part's cells that
    join its water at exactly that level and lie below it."""
    joined = st['joined_main'] if which == 'main' else st['joined_south']
    code = 1 if which == 'main' else 2
    return np.flatnonzero((st['part'] == code) & (joined == np.float32(level)) & (hf < level))


def part_levels(tab, rivers, neck):
    """The two parts' levels on their rivers: near_side_neck.balance's steady state while the main part stands over
    the crest and sends the neck's flow south, and each part where its own loss meets its own rivers once the main
    part stands under it."""
    base = N.balance(tab, rivers, neck)
    if base is not None:
        return dict(main=base['main_level'], south=base['south_level'], neck_flow=base['flow_m3s'],
                    choked=base['choked'], main_over_the_crest=True)
    alone = {k: N.level_for(tab[k], rivers[k], neck.crest) for k in ('main', 'south')}
    if None in alone.values():
        raise RuntimeError('a part\'s own rivers would raise it over the crest while the main part cannot stand there')
    return dict(main=alone['main'], south=alone['south'], neck_flow=0.0, choked=False, main_over_the_crest=False)


def settle_parts(d, st, tab, neck, centre, runoff, lake_net, net, area, fixed, branches, rivers, start=None,
                 on_round=None, log=print):
    """The near-side sea's parts, the routing and the levels set against each other (see the module's text). tab:
    part_tables on the climate's loss; centre: the neck's throat, where the main part's spill is put in past the
    crest; fixed and branches: govern's put-ins; rivers: what the parts receive with the sea at the common level, from
    which the first levels are set; start: the levels, history and sills a run reached, to take up; on_round is given
    that state after each routing. Returns the settled levels, the waters and their levels, the routing's structure,
    terms and state, the branches' flows, what reaches each water, each part's balance, what a part spills over a
    sill, the history and the sills met."""
    started = time.time()
    names = list(WATERS)
    level0, hf, cols = d['level0'], d['h'].ravel(), d['h'].shape[1]
    z = part_levels(tab, rivers, neck)
    history, sills, tried, spills = [], [], set(), {}
    taken_up = False                                 # a settled state kept by an earlier run is routed once more
    if start:
        z, history, sills = start['z'], start['history'], start['sills']
        tried = {(k, level) for k, level in start['tried']}
        taken_up = bool(start.get('settled'))
    kept = lambda settled=False: dict(z={k: (v if isinstance(v, bool) else float(v)) for k, v in z.items()},   # noqa: E731
                                      history=history, sills=sills, tried=sorted(tried), settled=settled)
    index = dict(main=MAIN, south=SOUTH)
    for _ in range(PART_ROUNDS):
        water, levels = parted_waters(d, st, z['main'], z['south'])
        s = structure_of(d, water, levels)
        terms = LR.hollow_terms(s, runoff, lake_net)
        # While the southern part stands under the crest, the neck's flow is a river on the ground past it, put in
        # where it leaves the main part. Once both parts stand over the crest they are one sea across the neck and
        # the flow between them passes through its water.
        over_the_neck, spilt = [], 0.0
        if z['neck_flow'] > 0.0 and z['south'] <= neck.crest:
            spilt = z['neck_flow']
            over_the_neck = [(DRt.spill_cell(hf, water == MAIN, cols, centre), spilt)]
        root = roots_of(s['tree'])
        state, flows, sources = branch_state(s, terms, fixed + over_the_neck, branches, root)
        got = LR.reaches(s, terms, state, water, len(names))
        rivers = dict(main=float(got[MAIN]), south=float(got[SOUTH]) - spilt)
        loss = {k: tab[k].at(z[k])[1] for k in ('main', 'south')}
        short = dict(main=rivers['main'] - z['neck_flow'] - loss['main'], south=rivers['south'] + z['neck_flow'] - loss['south'])
        this = dict(main_level_m=round(z['main'] - level0, 2), south_level_m=round(z['south'] - level0, 2),
                    neck_flow_m3s=round(z['neck_flow'], 1), main_rivers_m3s=round(rivers['main'], 1),
                    south_rivers_m3s=round(rivers['south'], 1), main_balance_m3s=round(short['main'], 1),
                    south_balance_m3s=round(short['south'], 1))
        log(f'near-side sea: {this}, {time.time() - started:.0f} s')
        if not taken_up:
            if this in history and not any((k, z[k]) in tried for k in ('main', 'south')):   # a sill is come back to once
                raise RuntimeError(f'the near-side sea\'s parts come round to a state they left: {history[history.index(this):]}')
            history.append(this)
        taken_up = False
        if on_round:
            on_round(kept())
        # a part's shore climbs by the levels its cells join at, metres apart over a plain, so a large part is held
        # settled within a share of what it receives
        within = dict(main=max(PART_SETTLED_M3S, PART_SETTLED_SHARE * rivers['main']),
                      south=max(PART_SETTLED_M3S, PART_SETTLED_SHARE * (rivers['south'] + z['neck_flow'])))
        if all(abs(short[k]) < within[k] for k in short):
            break
        z_next = part_levels(tab, rivers, neck)
        still = max(abs(z_next[k] - z[k]) for k in ('main', 'south')) < PART_STILL_M and \
            abs(z_next['neck_flow'] - z['neck_flow']) < PART_STILL_M3S
        back = any(z_next[k] != z[k] and (k, z_next[k]) in tried for k in ('main', 'south'))
        if not still or back:                        # back: a part tried over a sill returns to it
            z = z_next
            del s, terms, state
            continue
        # The levels hold with a part left over: it stands at the sill of a basin beside it, with more reaching it
        # than it loses below the sill.
        at_sill = [k for k in ('main', 'south') if short[k] >= within[k] and basin_at_the_sill(st, hf, k, z[k]).size]
        if not at_sill or any(abs(short[k]) >= within[k] for k in ('main', 'south') if k not in at_sill):
            raise RuntimeError(f'the near-side sea\'s parts hold out of balance: {this}')
        fresh = [k for k in at_sill if (k, z[k]) not in tried]
        for k in at_sill:
            basin = basin_at_the_sill(st, hf, k, z[k])
            lowest = int(basin[np.argmin(hf[basin])])
            note = dict(part=k, sill_above_common_m=round(z[k] - level0, 2), left_over_at_it_m3s=round(short[k], 1),
                        basin_below_the_sill_km2=round(float(area[basin].sum()) / 1e6),
                        basin_loses_m3s=round(-float(net[basin].sum())), lowest_at=LR.place(lowest, cols, d['rotation']))
            if k in fresh:
                # First the part is tried just over the sill, the basin under its water: where the basin is a lake
                # on its rivers' way, the lake's loss was taken from the rivers and is now the part's own, and the
                # part rises on from there.
                tried.add((k, z[k]))
                sills.append(dict(note, then='the part\'s water tried over the sill'))
                z = dict(z, **{k: float(np.nextafter(np.float32(z[k]), np.float32(np.inf)))})
            elif not fresh:
                # Over the sill the part lost more than reached it and came back: it stands at the sill and spills
                # what is left over into the basin, which keeps the lake that spill and its own rivers hold.
                spills[names[index[k]]] = short[k]
                sills.append(dict(note, then='the part stands at the sill and spills what is left over into the basin'))
                over_the_neck = over_the_neck + [(lowest, short[k])]
        log(f'at a sill: {sills[-len(at_sill):]}')
        if fresh:
            del s, terms, state
            continue
        state, flows, sources = branch_state(s, terms, fixed + over_the_neck, branches, root)
        again = LR.reaches(s, terms, state, water, len(names))
        if np.abs(again - got).max() >= min(within.values()):
            raise RuntimeError('a basin a part spills into over its sill fills and returns the water')
        got = again
        for k in at_sill:
            short[k] -= spills[names[index[k]]]
        break
    else:
        raise RuntimeError(f'the near-side sea\'s parts did not settle: {history[-3:]}')
    if on_round:
        on_round(kept(True))
    return dict(z=z, water=water, levels=levels, s=s, terms=terms, state=state, flows=flows, sources=sources, got=got,
                rivers=rivers, short=short, spills=spills, history=history, sills=sills)


# ---------------------------------------------------------------- the record's pieces

def seas_record(names, water, levels, level0, area, net, supply, spills=None):
    """Each water's level, area, loss and supply; its balance is what it receives less what it loses and what it
    spills over a sill into a basin beside it (`spills`, by name)."""
    out = {}
    for i, k in enumerate(names):
        cells = water == i
        loss = -float(net[cells].sum())
        spilt = float((spills or {}).get(k, 0.0))
        out[k] = dict(level_above_common_m=round(levels[i] - level0, 2), area_km2=round(float(area[cells].sum()) / 1e6),
                      loses_m3s=round(loss), receives_m3s=round(float(supply[i])),
                      balance_m3s=round(float(supply[i]) - loss - spilt))
        if spilt:
            out[k]['spills_over_a_sill_m3s'] = round(spilt)
    return out


def lakes_record(lake, spills, level, hf, area_km2, total, rows, cols):
    bands = {}
    for top_lat in range(90, -90, -30):
        part = slice(int((90 - top_lat) * rows / 180) * cols, int((120 - top_lat) * rows / 180) * cols)
        bands[f'{top_lat} to {top_lat - 30}'] = round(float(area_km2[part][lake[part]].sum()))
    a = float(area_km2[lake].sum())
    return dict(area_km2=round(a), share_of_moon=round(a / total, 5),
                area_overflowing_km2=round(float(area_km2[lake & spills].sum())),
                area_not_overflowing_km2=round(float(area_km2[lake & ~spills].sum())),
                volume_km3=round(float(((level - hf) * area_km2)[lake].sum() / 1000.0)), km2_by_latitude=bands)


def stages(key):
    """What a run on the same floor, rule and climate reached, kept so that a later run need not find it again: the
    factor with its budgets, taken as it is (--restart finds it again), and the weirs' shares and the parts' levels,
    which a later run only begins from and settles again."""
    if STAGES.exists():
        kept = json.loads(STAGES.read_text())
        if kept.get('key') == key:
            return kept
    return dict(key=key)


def keep(kept):
    WORK.mkdir(parents=True, exist_ok=True)
    STAGES.write_text(json.dumps(kept, indent=1) + '\n')


# ---------------------------------------------------------------- the run

def build(log=print, refresh=False):
    started = time.time()
    d = LR.design_inputs(log)
    st = setting(d, log, refresh)
    h, level0 = d['h'], d['level0']
    rows, cols = h.shape
    hf = h.ravel()
    design_case = json.loads(FD.OUT.read_text())['cases'][DRt.BUDGET]
    design_factor = float(design_case['evaporation'])
    neck_record = json.loads(N.OUT.read_text())
    design_neck = neck_record['budgets'][DRt.BUDGET]
    names = list(WATERS)
    if not (np.array_equal(d['water'] == 3, st['common'] == 3) and np.array_equal(d['water'] == 4, st['common'] == 4)):
        raise SystemExit('the seas at the common level are not the design structure\'s')
    z_design = d['levels'][MAIN], d['levels'][SOUTH]
    water_d, _ = parted_waters(d, st, *z_design)
    if not np.array_equal(water_d, d['water']):
        raise SystemExit('the near-side sea\'s parts are not the design structure\'s')
    del water_d
    kept = stages(dict(setting=st['meta']['key'], rule=sha256(Path(LR.__file__)), climate=sha256(Path(FD.__file__))))

    # the factor
    if 'factor' in kept:
        factor, closing = kept['factor'], kept['closing']
        log(f'the factor {factor} from {STAGES}')
    else:
        factor, closing = closing_factor(d, st, design_factor, log)
        kept.update(factor=factor, closing=closing)
        keep(kept)
    log(f'the factor that closes the budget on the corrected rule: {factor} (the design\'s {design_factor}), '
        f'{time.time() - started:.0f} s')

    # the design's own state at its factor: each sea's loss, and the neck read from its record against the design's
    runoff, lake_net, _, _ = FD.load(design_factor)
    area = (((d['cell_km'] * 1000.0) ** 2 * np.maximum(np.cos(np.radians(d['lat'])), 1e-6))[:, None] * np.ones((1, cols))).ravel()
    area_km2 = area / 1e6
    total = float(area_km2.sum())
    net_design = lake_net.ravel() * MM_DAY * area
    neck = neck_from_record(st, neck_record)
    tab = part_tables(d, st, -net_design, area)
    again = part_levels(tab, dict(main=float(design_neck['supply_m3s']['main']), south=float(design_neck['supply_m3s']['south'])), neck)
    steady = design_neck['natural_neck']['steady']
    control = dict(main_level_m=round(again['main'] - level0, 2), south_level_m=round(again['south'] - level0, 2),
                   neck_flow_m3s=round(again['neck_flow']), design=dict(main_level_m=steady['main_level_m'],
                   south_level_m=steady['south_level_m'], neck_flow_m3s=steady['neck_flow_m3s']))
    log(f'the neck from its record on the design\'s rivers: {control}')
    if abs(control['main_level_m'] - steady['main_level_m']) > CONTROL_LEVEL_M or \
            abs(control['south_level_m'] - steady['south_level_m']) > 20 * CONTROL_LEVEL_M:
        raise SystemExit('the neck read from its record does not give the design\'s levels back')
    design_volume = {k: tab[k].at(z)[2] / 1e9 for k, z in zip(('main', 'south'), z_design)}
    design_seas = {k: dict(level_above_common_m=round(d['levels'][i] - level0, 2),
                           area_km2=round(float(area_km2[d['water'] == i].sum())),
                           loses_m3s=round(-float(net_design[d['water'] == i].sum())))
                   for i, k in enumerate(names)}
    del tab, net_design

    # the weirs, on the corrected rule's factor with the near-side sea at the common level
    runoff, lake_net, _, _ = FD.load(factor)
    net = lake_net.ravel() * MM_DAY * area
    g = govern(d, st, runoff, lake_net, start=kept.get('weirs'), log=log)
    gov, governed_supply, fixed, branches = g['gov'], g['supply'], g['fixed'], g['branches']
    if 'weirs' not in kept:                          # a later run begins from these and keeps how they were reached
        kept['weirs'] = dict(share=gov['share'], over=gov['over'], history=gov['history'])
        keep(kept)
    weir_history = kept['weirs']['history']
    log(f'governed: supplies {dict((k, round(float(v))) for k, v in zip(names, governed_supply))}, least kept by a lake '
        f'not full {g["least_kept_by_a_lake_not_full"]:.3g}, {time.time() - started:.0f} s')
    own_g, loss_g = gov['own'], gov['loss']
    governed = dict(
        moscoviense_surplus_m3s=round(gov['surplus']),
        released_m3s={k: round(v) for k, v in gov['share'].items()},
        arrived_m3s={k: round(float(governed_supply[names.index(k)] - own_g[names.index(k)])) for k in (SALT_SEA, HELD, REST)},
        full_seas_pass_on_m3s={k: round(v) for k, v in gov['over'].items()},
        held_seas_left_short_m3s={k: round(v) for k, v in gov['unmet'].items()},
        branch_flows_m3s={k: [round(q, 1) for q in v] for k, v in gov['flows'].items()},
        near_side_at_the_common_level=dict(
            loses_m3s=round(float(loss_g[MAIN] + loss_g[SOUTH])), own_rivers_m3s=round(float(own_g[MAIN] + own_g[SOUTH])),
            receives_m3s=round(float(governed_supply[MAIN] + governed_supply[SOUTH])),
            main_rivers_m3s=round(float(governed_supply[MAIN])), south_rivers_m3s=round(float(governed_supply[SOUTH]))),
        seas={k: dict(loses_m3s=round(float(loss_g[i])), own_rivers_m3s=round(float(own_g[i])),
                      receives_m3s=round(float(governed_supply[i]))) for i, k in enumerate(names[:MAIN])},
        weirs_set_in_settlings=len(weir_history), weir_history=weir_history)
    rivers = dict(main=float(governed_supply[MAIN]), south=float(governed_supply[SOUTH]))
    del g, gov

    # the near-side sea's parts, the routing and the levels set against each other
    tab = part_tables(d, st, -net, area)

    def on_round(state_now):
        kept['parts'] = state_now
        keep(kept)

    parts = settle_parts(d, st, tab, neck, neck_record['throat']['centre_today'], runoff, lake_net, net, area, fixed,
                         branches, rivers, start=kept.get('parts'), on_round=on_round, log=log)
    z, water, levels, s, terms, state = (parts[k] for k in ('z', 'water', 'levels', 's', 'terms', 'state'))
    flows, got, rivers, short, spills = (parts[k] for k in ('flows', 'got', 'rivers', 'short', 'spills'))
    history, sills = parts['history'], parts['sills']
    del parts

    # the rebalanced state's routing
    r = LR.hollow_routing(s, terms, state)
    q = r['discharge']
    supply = np.array([float(q[water == i].sum()) for i in range(len(names))])
    if np.abs(supply - got).max() > 1e-3:
        raise SystemExit('the rivers gathered do not bring the seas what the hollows\' overflows give them')
    b = r['budget']
    gap = b['runoff_of_dry_ground'] + b['put_in'] + b['lakes_rain_less_evaporation'] - b['reaches_the_seas'] - b['kept_by_lakes_not_full']
    if float(r['left'].min()) < -1e-6:
        raise SystemExit('a lake not full keeps less than nothing: a cut takes more than its lake holds')
    seas = seas_record(names, water, levels, level0, area, net, supply, spills)
    seas[names[MAIN]]['sends_over_the_neck_m3s'] = round(z['neck_flow'])
    seas[names[MAIN]]['balance_m3s'] = round(short['main'])
    log(f'rebalanced: {seas}, budget closes within {gap:.3g} m3/s, {time.time() - started:.0f} s')
    lake = r['lake'] & (water < 0)
    lakes = lakes_record(lake, r['lake_spills'], r['lake_level'], hf, area_km2, total, rows, cols)
    lakes.update(LR.unheld_shores(r, h, water, levels, area_km2))
    root = roots_of(r['tree'])
    lacks = np.bincount(root, weights=r['need'] - r['held'], minlength=root.size)
    take_outs = []
    for (ops_, _), qs, kind in zip(branches, flows, ('humboldtianum_branch', 'smythii_marginis_branch')):
        for (c, sign), flow in zip(ops_, qs):
            if sign < 0:
                top = int(root[r['label'][int(c)]])
                take_outs.append(dict(kind=kind, at=LR.place(int(c), cols, d['rotation']), takes_m3s=round(flow, 1),
                                      outer_pass_above_common_m=round(float(r['tree']['spill'][top] - level0), 1),
                                      stands_full=bool(r['full'][top]), overflows_besides_m3s=round(float(r['out'][top]), 1),
                                      short_of_standing_full_m3s=round(float(max(lacks[top], 0.0)), 1)))

    # beside the design
    design = DRt.load(log)
    was_lake = np.zeros(hf.size, bool)
    was_lake[design['lake_cell']] = True
    was_level = np.full(hf.size, np.nan)
    was_level[design['lake_cell']] = design['lake_level_m']
    was_spills = np.zeros(hf.size, bool)
    was_spills[design['lake_cell']] = design['lake_spills']
    design_lakes = lakes_record(was_lake, was_spills, was_level, hf, area_km2, total, rows, cols)
    was_sea = d['water'] >= 0
    now_sea = water >= 0
    changed = dict(newly_lake_km2=round(float(area_km2[lake & ~was_lake & ~was_sea].sum())),
                   no_longer_lake_km2=round(float(area_km2[was_lake & ~lake & ~now_sea].sum())),
                   sea_floor_left_dry_km2=round(float(area_km2[was_sea & ~now_sea & ~lake].sum())),
                   sea_floor_under_lakes_km2=round(float(area_km2[was_sea & ~now_sea & lake].sum())),
                   ground_newly_under_a_sea_km2=round(float(area_km2[now_sea & ~was_sea].sum())))
    volume = {k: tab[k].at(z[k])[2] / 1e9 for k in ('main', 'south')}
    water_held = dict(
        near_side_main_km3=dict(design=round(design_volume['main']), rebalanced=round(volume['main'])),
        near_side_south_km3=dict(design=round(design_volume['south']), rebalanced=round(volume['south'])),
        lakes_km3=dict(design=design_lakes['volume_km3'], rebalanced=lakes['volume_km3']))
    water_held['change_km3'] = round(sum(v['rebalanced'] - v['design'] for v in water_held.values()))
    bodies = LR.lake_bodies(r, hf, water, area_km2, was_lake, level0, cols, d['rotation'])
    lakes_closer = dict(lakes_as_bodies=bodies['count'], bodies_passing_nothing_on=bodies['passing_nothing_on'],
                        largest=bodies['largest'][:LAKES_LISTED],
                        largest_passing_nothing_on=bodies['largest_passing_nothing_on'][:LAKES_LISTED])
    record = json.loads(LR.WATERWAYS.read_text())
    ends = LR.waterway_cells(design['discharge_m3s'], cols, record)
    waterways = []
    for x, e in zip(record['waterways'], ends):
        row = dict(name=x.get('name'), anchor=x['anchor'], ends_in=x['ends_in'], fed_by=x.get('fed_by'))
        mouth, head = e['mouth'][1], e['head'][1]
        path = LR.course(design['receiver'], head, mouth)
        ashore, ashore_was = path[water[path] < 0], path[d['water'][path] < 0]
        risen = bool(ashore.size and ashore_was.size and ashore[-1] != ashore_was[-1])
        # where the river's sea has risen up its course, it meets the river above the streams that joined it on
        # the reach now under water: both flows are then read at the last cell of the course still ashore
        shore = int(ashore[-1]) if risen else mouth
        was, now = float(design['discharge_m3s'][shore]), float(q[shore])
        row['mouth'] = dict(design_m3s=round(float(design['discharge_m3s'][mouth])), rebalanced_m3s=round(now),
                            change=round(now / max(was, 1e-9) - 1.0, 4), its_sea_has_risen_to_meet_it=risen,
                            recorded_m3s=x['flow_mouth_m3s'])
        if risen:
            row['mouth']['design_where_the_sea_now_meets_it_m3s'] = round(was)
        row['head'] = dict(design_m3s=round(float(design['discharge_m3s'][head])), rebalanced_m3s=round(float(q[head])),
                           under_a_lake=bool(lake[head]), recorded_m3s=x['flow_head_m3s'])
        if now < (1.0 - LR.LOST_SHARE) * was:        # the largest lakes its design course crosses
            rows_on = np.unique(bodies['row'][path][bodies['row'][path] >= 0])
            largest = rows_on[np.argsort(-bodies['area'][rows_on])][:COURSE_LAKES_LISTED]
            row['largest_lakes_on_its_course'] = [dict(passes_water_on=bool(bodies['spills'][k]), **bodies['describe'](int(k)))
                                                  for k in largest]
        waterways.append(row)
    del bodies
    measured = json.loads(LR.OUT.read_text())
    design_supply = measured['cases']['present']['seas_receive_m3s']
    for k in names:
        design_seas[k]['receives_m3s'] = design_supply[k]
        design_seas[k]['balance_m3s'] = design_supply[k] - design_seas[k]['loses_m3s']
    spare_at_design = rivers['main'] + rivers['south'] + float(net[d['water'] == MAIN].sum() + net[d['water'] == SOUTH].sum())
    design_sent = measured['cases']['present']['all_seas_receive_m3s'] - round(sum(q_ for _, q_, _ in d['ops']))
    BULK.mkdir(parents=True, exist_ok=True)
    cells = np.flatnonzero(lake).astype(np.int32)
    np.savez_compressed(PRODUCT, metadata=json.dumps(dict(
        schema=SCHEMA + '-grid', grid=dict(pixels_per_degree=dr.PPD, rows=rows, cols=cols), common_level_m=level0,
        waters=names, levels_m=[float(v) for v in levels], open_water_evaporation_factor=factor,
        reading_rule=('cell = row * 5760 + column, rows from 90 N and columns from 0 E on today\'s Moon; water is the '
                      'index into waters of the sea over a cell, -1 off the seas; discharge_m3s the flow through each '
                      'cell; lake_cell lists the cells under lakes off the seas, with lake_level_m and lake_spills.'))),
        discharge_m3s=q, water=water, lake_cell=cells, lake_level_m=r['lake_level'][cells].astype(np.float32),
        lake_spills=r['lake_spills'][cells])
    return dict(
        schema=SCHEMA, producer=dict(files={f: sha256(ROOT / f) for f in CODE}),
        inputs={rel(p): sha256(p) for p in (DRt.OUT, FD.OUT, N.OUT, LR.OUT, LR.WATERWAYS)},
        ignored_inputs={str(SR.STRUCTURE): sha256(SR.STRUCTURE), str(DRt.PRODUCT): sha256(DRt.PRODUCT),
                        FD.CLIMATOLOGY: sha256(ROOT / FD.CLIMATOLOGY)},
        evidence=('A calculation on the design routing\'s own floor, works and climate with a lake rule that fills '
                  'each hollow as one water: the factor on the open water\'s evaporation, the weirs\' shares and the '
                  'near-side sea\'s levels are set again by the procedures that set them on the present rule, the '
                  'search for a part\'s level taking the cells that join at one level together. It is '
                  'a steady state: a lake stands wherever its supply would keep it, however long it would take to '
                  'fill. The climate is the design climate, run on a surface with the design\'s seas and lakes; no '
                  'climate has been run on the state found here. Nothing is validated against an observation.'),
        reading_rule=('Levels in m above the redrawn common level, flows in m3/s, areas in km2, volumes in km3. '
                      'factor: closing gives final_drainage.budget\'s figures with every sea at the common level, on '
                      'the corrected rule at factor 1 and at the factor that step gives, and on the present rule at '
                      'the design\'s factor; surplus_m3s at a rule\'s own factor is what its budget still leaves over. '
                      'seas: each water of the design routing (the near-side sea as its two parts) in the design and '
                      'rebalanced; loses is the open water\'s loss beyond the rain on it at the state\'s own factor, '
                      'receives what rivers, overflows and works bring it. governed is the state the weirs are set on, '
                      'the near-side sea at the common level: released is what each weir or branch lets go and arrived '
                      'what reaches its sea beyond its own rivers; branch_flows_m3s is each branch\'s flow at each of '
                      'its put-ins and take-outs in order. near_side.history is each routing of the parted sea with '
                      'the balance each part was left with; control is the neck read from its record on the design\'s '
                      'rivers beside the design\'s steady state. near_side.sills_met lists each sill of a basin beside '
                      'a part at which the part\'s level came to rest with water left over: the sill, what was left '
                      'over, the basin\'s ground below the sill, what the basin loses to the air, its lowest point, and '
                      'what was done (the part tried with the basin under its water, or left at the sill spilling into '
                      'the basin). near_side.to_spare_at_the_design_levels_m3s is what the rivers that reach both '
                      'parts bring beyond what the parts would lose at their design levels on the new factor. lakes '
                      'are off the seas; shore_cells_not_held and ground_below_and_joined_km2 are lake_support.py\'s '
                      'measures; lakes.closer lists the largest lakes as bodies (a hollow\'s one water) and the largest '
                      'that pass nothing on, each with its place on today\'s Moon and in the settled frame. changed '
                      'sets the rebalanced lakes and seas against the design\'s cell by cell. '
                      'land_and_lakes_send_the_seas_m3s is what reaches the seas less what the works and overflows '
                      'carry from sea to sea. waterways gives each named waterway\'s flow where it ends, at its design '
                      'mouth cell, and at its design head cell; where its sea has risen up its design course, the '
                      'rebalanced flow is read at the last cell of that course not under the sea, the design\'s '
                      'flow at that cell is given beside it, and change sets the one against the other. For a '
                      'waterway that loses more than a twentieth of its flow, the largest lakes its design course '
                      'crosses are listed.'),
        factor=dict(design=design_factor, rebalanced=factor, closing=closing),
        seas=dict(design=design_seas, rebalanced=seas),
        seas_share_of_moon=dict(design=round(float(area_km2[was_sea].sum()) / total, 5),
                                rebalanced=round(float(area_km2[now_sea].sum()) / total, 5)),
        land_share_of_moon=dict(design=round(float(area_km2[~was_sea & ~was_lake].sum()) / total, 5),
                                rebalanced=round(float(area_km2[~now_sea & ~lake].sum()) / total, 5)),
        governed=governed,
        near_side=dict(main_over_the_crest=bool(z['main_over_the_crest']), neck_flow_m3s=round(z['neck_flow']),
                       crest_above_common_m=round(st['meta']['join4'] - level0, 2), history=history, control=control,
                       sills_met=sills,
                       left_out_of_balance_m3s={k: round(v, 1) for k, v in short.items()},
                       to_spare_at_the_design_levels_m3s=round(spare_at_design)),
        works=dict(design={k: design_case[k] for k in ('moscoviense_surplus_m3s', 'released_m3s', 'arrived_m3s',
                                                        'full_seas_pass_on_m3s')},
                   rebalanced={k: governed[k] for k in ('moscoviense_surplus_m3s', 'released_m3s', 'arrived_m3s',
                                                        'full_seas_pass_on_m3s')},
                   branch_flows_m3s={'humboldtianum_branch': [round(v, 1) for v in flows[0]],
                                     'smythii_marginis_branch': [round(v, 1) for v in flows[1]]},
                   take_outs=take_outs),
        lakes=dict(design=design_lakes, rebalanced=lakes, closer=lakes_closer), changed=changed, water_held=water_held,
        land_and_lakes_send_the_seas_m3s=dict(design=design_sent, rebalanced=round(b['reaches_the_seas'] - b['put_in'])),
        budget_m3s={k: round(v, 1) for k, v in b.items()}, budget_gap_m3s=round(gap, 3),
        waterways=waterways, parameters=dict(settle_m3s=SETTLE_M3S, part_settled_m3s=PART_SETTLED_M3S,
                                             lakes_listed=LAKES_LISTED, course_lakes_listed=COURSE_LAKES_LISTED),
        product=dict(path=str(PRODUCT), sha256=sha256(PRODUCT)), seconds=round(time.time() - started))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--refresh', action='store_true', help='find the seas\' cells again and start the settings afresh')
    p.add_argument('--restart', action='store_true', help='start the settings afresh on the kept seas\' cells')
    args = p.parse_args(argv)
    log = lambda *a: print(*a, flush=True)                                     # noqa: E731
    if (args.refresh or args.restart) and STAGES.exists():
        STAGES.unlink()
    record = build(log, args.refresh)
    OUT.write_text(json.dumps(record, indent=1, ensure_ascii=False) + '\n')
    log(f'wrote {OUT} in {record["seconds"]} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
