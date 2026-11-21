"""The atlas round's first step: the Moon's surface redrawn under its water.

    ulimit -v 14000000; nice -n 10 python -m research.studies.atlas_round.redrawn_surface     # about half an hour

The atlas fills LOLA's floor to one level of the dry Moon's geoid. The round's frame (research/decisions.md, "The
atlas round's frame, 16 November 2026") redraws that surface on the water's own gravity and the crust's sag and
slow creep, with the seas covering 28% at the common level and Moscoviense, Orientale and Smythii-Marginis filled to
their rims on top. This module solves it:

- the water: every cell below one common level is water, the level set so that those cells cover 28% of the surface
  (the atlas's storage geometry); each upper sea stands at its rim, the lowest pass from its water to another sea's
  (the drainage network study's rule), over all the land its water joins below the rim; the lakes are a fixed load
  at the turned Moon's drainage depths wherever no sea covers them;
- the Moon's response: geography's crustal-loading study's elastic load Love numbers at every degree, the response
  within years, plus the creep of the present Moon's mantle from the geology lane's viscoelastic response: its
  change in h' and k' at degrees 2-50 between the load's arrival and a chosen time, under the central, stiff or soft
  temperature profile;
- the sea-level equation iterated as geography.crustal_loading.solve iterates it, with the upper seas held at their
  rims, and the rims found again on the deformed heights until none moves by more than 5 cm. A case whose rim keeps
  flipping between two passes, each made the lower by the load of filling to the other, is reported with the range
  of its last three rounds and flagged as unsettled.

The atlas stands at the crust's settled state, barring erosion: LOLA's landforms as they are, under the sag and creep
the water has finished causing. Under the present Moon's central temperatures the creep grows by about the same
amount every tenfold of time and does not settle within the geology lane's billion years; that lane holds the present
thermal state, which it takes to describe the Moon for about 100 million years before the cooling interior stiffens
the mantle and the sag already made is kept. The design case is therefore the central profile 100 million years
after the water arrives. The central and stiff profiles at a billion years give the most the present temperatures
could add; the soft profile, the warm mantle carried as a stated risk, settles within about 10,000 years and is taken
at 10 million. The elastic response and the central profile at 1,000 years give the first years and the first
millennium, through which the shores move toward the settled state. The 'atlas' case deforms nothing and reproduces
the drainage network study's rims on the dry Moon.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np

from geography import crustal_loading as CL
from geography import fetch_inputs
from geography import love_numbers as ln
from geography import topography as tp
from research.studies.drainage_network import network as NW
from research.studies.south_sea_ocean import inflows as I
from shared.constants import MOON_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'redrawn_surface.json'
BULK = Path('/media/projectspace/terluna-research/atlas-round/redrawn_surface')
SCHEMA = 'terluna.research.atlas-round-redrawn-surface/1'
ATLAS = 'geography/results/atlas.json'
LAKES = 'geography/turned_corrected/products/drainage_28pct_16ppd.npz'
RESPONSE = Path('/media/projectspace/terluna-research/geology/interior-response/interior_response_love.npz')
CODE = ('research/studies/atlas_round/redrawn_surface.py', 'geography/crustal_loading.py', 'geography/love_numbers.py',
        'geography/topography.py', 'geography/drainage.py', 'research/studies/drainage_network/network.py',
        'research/studies/drainage_network/routing.py', 'research/studies/south_sea_ocean/inflows.py')
UPPER = ('moscoviense', 'orientale', 'smythii_marginis')
SHARE = CL.SHARE
FINE_PPD = CL.FINE_PPD
FACTOR = CL.TRANSFORM_FACTOR
RIM_TOLERANCE_M = 0.05                  # the rims are found again until none moves by more than this
OUTER_ROUNDS = 8
CASES = {                               # the profile of the geology lane's response and the years since the water arrived
    'atlas': None,
    'elastic': ('central', 0.0),
    'central, 1,000 years': ('central', 1.0e3),
    'central, 100 million years': ('central', 1.0e8),
    'central, 1 billion years': ('central', 1.0e9),
    'stiff, 1 billion years': ('stiff', 1.0e9),
    'soft, 10 million years': ('soft', 1.0e7),
}
DESIGN = 'central, 100 million years'   # the settled state, barring erosion


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


# ---------------------------------------------------------------- the Moon's response

def creep_increment(table, profile, years):
    """The change in the load Love numbers h' and k' since the load arrived, at each of the table's degrees, after
    `years` under one profile of the geology lane's response, log-linear in time between its tabulated times."""
    t = np.asarray(table['times_years'], float)
    h = np.asarray(table[f'{profile}__h_load'], float)
    k = np.asarray(table[f'{profile}__k_load'], float)
    if years <= 0.0:
        return np.zeros(h.shape[0]), np.zeros(k.shape[0])
    if not t[1] <= years <= t[-1]:
        raise ValueError(f'{years} years lies outside the response table ({t[1]:g}-{t[-1]:g})')
    j = int(np.searchsorted(t, years))
    if t[j] == years:
        hy, ky = h[:, j], k[:, j]
    else:
        w = (math.log(years) - math.log(t[j - 1])) / (math.log(t[j]) - math.log(t[j - 1]))
        hy, ky = h[:, j - 1] + w * (h[:, j] - h[:, j - 1]), k[:, j - 1] + w * (k[:, j] - k[:, j - 1])
    return hy - h[:, 0], ky - k[:, 0]


def love_numbers(elastic, table, profile, years, nmax):
    """The elastic load Love numbers to degree nmax with the creep's change added at the table's degrees."""
    h, k = (np.array(x[:nmax + 1], float) for x in elastic)
    dh, dk = creep_increment(table, profile, years)
    deg = np.asarray(table['degrees'], int)
    keep = deg <= nmax
    h[deg[keep]] += dh[keep]
    k[deg[keep]] += dk[keep]
    return h, k


# ---------------------------------------------------------------- the water

def water_depth(eff, level, regions):
    """The water over each cell: up to the common level wherever the ground lies below it, and up to an upper sea's
    rim over the land its water joins below the rim. `regions` is a list of (mask, rim)."""
    depth = np.maximum(level - eff, 0.0)
    for mask, rim in regions:
        depth[mask] = np.maximum(rim - eff[mask], 0.0)
    return depth


TIED_M = 1.0                            # seas whose passes lie within this of the rim are joined at the rim


def rims_from(joined, seeds, seas, upper):
    """Each upper sea's rim, the lowest level at which its water joins another sea's, and every sea it joins there:
    past its pass the water runs down to several seas at the same level, so the sea its overflow reaches is the
    routing's to say (the drainage network study's spill path), not the rim's."""
    out = {}
    for k in upper:
        passes = {j: float(joined[k][seeds[j]]) for j in seas if j != k}
        level = min(passes.values())
        out[k] = (level, sorted(j for j, p in passes.items() if p <= level + TIED_M))
    return out


def geometry(eff32, level, grid, bodies, log=print):
    """The six seas at the common level on the given heights, each upper sea's rim and the land its water joins below
    the rim (the drainage network study's catchment graph and flood, run on these heights)."""
    started = time.time()
    lat, lon, cell_km, row_area = grid['lat'], grid['lon'], grid['cell_km'], grid['row_area']
    body = I.label_waters(eff32 < level)
    named, _ = I.name_seas(body, row_area, lat, ((lon + 180.0) % 360.0) - 180.0, bodies)
    seas = list(I.SEAS.values())
    flat = eff32.ravel()
    seeds = {}
    for k in seas:
        cells = np.flatnonzero(body.ravel() == named[k]['label'])
        seeds[k] = int(cells[np.argmin(flat[cells])])
    label, n_nodes, edges = NW.catchment_graph(eff32, lat, cell_km)
    joined = {k: NW.joining_level(flat, label, n_nodes, edges, seeds[k]) for k in UPPER}
    del label, edges
    rims = rims_from(joined, seeds, seas, UPPER)
    regions = {k: (joined[k] < rims[k][0] - NW.BELOW_RIM_M).reshape(eff32.shape) for k in UPPER}
    del joined
    log(f'  geometry: {n_nodes:,} catchments, rims ' + ', '.join(f'{k} {rims[k][0] - level:.2f} m'
                                                                 for k in UPPER) + f', {time.time() - started:.0f} s')
    return dict(seas=seas, seeds=seeds, rims=rims, regions=regions, body=body, named=named, catchments=n_nodes)


def lake_load(height, cells, levels, covered):
    """The lakes' weight (kg/m2) at their drainage depths, except where a sea covers them."""
    load = np.zeros(height.size)
    load[cells] = CL.WATER_DENSITY * np.maximum(levels - height.ravel()[cells], 0.0)
    load[covered.ravel()] = 0.0
    return load.reshape(height.shape)


# ---------------------------------------------------------------- the sea-level equation with the upper seas at their rims

def solve(height, row_area, response, lakes, level, regions, dev, iterations=60):
    """geography.crustal_loading.solve's iteration with the common level holding 28% of the surface and each upper sea
    at its rim: returns the level, the effective heights' change U - N (coarse grid), U, N and the history."""
    history = []
    for _ in range(iterations):
        eff = height + CL.upsample(dev, FACTOR)
        level = CL.level_for_area(eff, row_area, SHARE, level)
        load = CL.WATER_DENSITY * water_depth(eff, level, regions)
        del eff
        load += lakes
        U, N = response(CL.block_mean(load, FACTOR))
        del load
        change = float(np.max(np.abs(U - N - dev)))
        dev = U - N
        history.append(round(change, 4))
        if change < CL.TOLERANCE_M:
            break
    else:
        raise RuntimeError(f'Sea-level equation did not converge: {history[-5:]}')
    return dict(level=level, dev=dev, U=U, N=N, history=history)


# ---------------------------------------------------------------- measures

def strait_seeds(height):
    """The straits' basin seeds and the near-side maria's, as geography.crustal_loading.run places them."""
    seeds = {name: (CL.deepest_near(height, *spec['a']), CL.deepest_near(height, *spec['b']))
             for name, spec in CL.STRAITS.items() if spec['a']}
    rows, cols = CL.box_slices(dict(lat=(-86.0, -70.0), lon=(96.0, 161.0)))
    i, j = np.unravel_index(np.argmin(height[rows, cols]), height[rows, cols].shape)
    rows2, cols2 = CL.box_slices(dict(lat=(-70.0, -30.0), lon=(150.0, 230.0)))
    i2, j2 = np.unravel_index(np.argmin(height[rows2, cols2]), height[rows2, cols2].shape)
    seeds['Tenebrarum sub-basin to the South Pole-Aitken Sea'] = ((rows.start + i, cols.start + j),
                                                                   (rows2.start + i2, cols2.start + j2))
    mares = {k: CL.deepest_near(height, *v) for k, v in CL.NEARSIDE_MARIA.items()}
    return seeds, mares


def place(cell, cols):
    i, j = divmod(int(cell), cols)
    return [round(90.0 - (i + 0.5) / FINE_PPD, 3), round(((j + 0.5) / FINE_PPD + 180.0) % 360.0 - 180.0, 3)]


def measures(height, eff, level, level0, geo, grid, U, straits, mares, lakes_volume_m3, atlas):
    """What a redrawn surface holds: the levels, each sea's water, the passes and sills, the shore and the totals."""
    row_area, cell_km = grid['row_area'], grid['cell_km']
    rows, cols = height.shape
    total_area = float(row_area.sum()) * cols
    regions = [(geo['regions'][k], geo['rims'][k][0] - NW.BELOW_RIM_M) for k in UPPER]
    depth = water_depth(eff, level, regions)
    by_row = lambda x: float(row_area @ x.sum(axis=1))
    storage = eff < level
    seas = {}
    body = geo['body']
    for k in geo['seas']:
        if k in UPPER:
            mask, surface = geo['regions'][k], geo['rims'][k][0] - NW.BELOW_RIM_M
        else:
            mask, surface = body == geo['named'][k]['label'], level
        d = np.where(mask, np.maximum(surface - eff, 0.0), 0.0)
        seed = geo['seeds'][k]
        entry = dict(surface_above_the_common_level_m=round(surface - level, 2),
                     area_km2=round(by_row(d > 0) / 1e6), volume_km3=round(by_row(d) / 1e9),
                     greatest_depth_m=round(float(d.max()), 1), deepest_place=place(int(np.argmax(d)), cols),
                     floor_fall_at_the_deepest_m=round(float(-U.ravel()[int(np.argmax(d))]), 1))
        if k in UPPER:
            at_common = np.where(mask & storage, level - eff, 0.0)
            entry.update(joins_at_its_rim=geo['rims'][k][1], rise_volume_km3=round((by_row(d) - by_row(np.maximum(at_common, 0.0))) / 1e9),
                         seed=place(seed, cols))
        seas[k] = entry
        del d, mask
    sea_volume = by_row(depth)
    out = dict(level_change_m=round(level - level0, 2),
               common_level_share=round(by_row(storage) / total_area, 5),
               water_share=round(by_row(depth > 0) / total_area, 5),
               seas=seas,
               sills={k: v['sill_depth_m'] for k, v in CL.strait_sills(eff, level, straits).items()},
               nearside=CL.nearside_join(eff, level, mares),
               shore=CL.shore_measures(height, eff, level0, level, row_area, cell_km),
               seas_volume_km3=round(sea_volume / 1e9), lakes_volume_km3=round(lakes_volume_m3 / 1e9),
               water_volume_km3=round((sea_volume + lakes_volume_m3) / 1e9),
               water_beyond_the_atlas_km3=round((sea_volume + lakes_volume_m3) / 1e9 - atlas['water_km3']),
               water_beyond_the_atlas_global_layer_m=round(((sea_volume + lakes_volume_m3) / 1e9 - atlas['water_km3']) * 1e9
                                                           / (4.0 * math.pi * MOON_RADIUS ** 2), 2))
    del depth, storage
    return out


# ---------------------------------------------------------------- the run

def load_inputs(nmax, log=print):
    started = time.time()
    if fetch_inputs.restore()['status'] != 'PASS':
        raise SystemExit('BLOCKED: topography inputs missing; run python -m geography.fetch_inputs --download')
    height, _, g = tp.height_above_geoid(FINE_PPD, 200)
    rows, cols = height.shape
    atlas = json.loads((ROOT / ATLAS).read_text())
    level0 = float(atlas['sea_level_m'])
    grid = dict(lat=90.0 - (np.arange(rows) + 0.5) / FINE_PPD, lon=(np.arange(cols) + 0.5) / FINE_PPD,
                cell_km=MOON_RADIUS / 1000.0 * np.radians(1.0 / FINE_PPD), row_area=g.cell_area_m2[:, 0].copy())
    with np.load(ROOT / LAKES, allow_pickle=False) as d:
        cells, levels = d['lake_cell'], d['lake_level_m'].astype(float)
        meta = json.loads(d['metadata'].item())
    if meta['grid']['rows'] != rows or meta['grid']['cols'] != cols:
        raise ValueError('The lakes grid does not match the heights')
    filt = CL.Filter(rows // FACTOR, cols // FACTOR, nmax)
    h, k, _ = ln.elastic_spectrum(ln.Interior(ln.interior_model(compressible=True)), nmax)
    with np.load(RESPONSE, allow_pickle=False) as d:
        table = {key: d[key] for key in d.files}
    atlas_lakes = float(np.sum(np.maximum(levels - height.ravel()[cells], 0.0) * np.repeat(grid['row_area'], cols)[cells]))
    atlas_seas = float(grid['row_area'] @ np.maximum(level0 - height, 0.0).sum(axis=1))
    log(f'inputs: heights {rows} x {cols}, lakes {cells.size:,} cells, elastic Love numbers to degree {nmax}, '
        f'{time.time() - started:.0f} s')
    return dict(height=height, level0=level0, grid=grid, bodies=atlas['bodies'], cells=cells, levels=levels,
                filt=filt, elastic=(h, k), table=table,
                atlas=dict(seas_km3=round(atlas_seas / 1e9), lakes_km3=round(atlas_lakes / 1e9),
                           water_km3=round((atlas_seas + atlas_lakes) / 1e9)))


def redraw(name, spec, inp, log=print):
    """One case: the surface redrawn with the upper seas at their rims, the rims found again until they settle."""
    started = time.time()
    height, level0, grid = inp['height'], inp['level0'], inp['grid']
    rows, cols = height.shape
    coarse = (rows // FACTOR, cols // FACTOR)
    dev, level, U, N = np.zeros(coarse), level0, np.zeros(coarse), np.zeros(coarse)
    response = None
    if spec is not None:
        response = CL.Response(*love_numbers(inp['elastic'], inp['table'], spec[0], spec[1], inp['filt'].nmax), inp['filt'])
    log(f'{name}:')
    geo = geometry(height.astype(np.float32), level, grid, inp['bodies'], log)
    rounds = []
    settled = response is None
    for _ in range(OUTER_ROUNDS if response is not None else 0):
        regions = [(geo['regions'][k], geo['rims'][k][0] - NW.BELOW_RIM_M) for k in UPPER]
        eff = height + CL.upsample(dev, FACTOR)
        covered = eff < level
        del eff
        for mask, _ in regions:
            covered |= mask
        lakes = lake_load(height, inp['cells'], inp['levels'], covered)
        del covered
        sol = solve(height, grid['row_area'], response, lakes, level, regions, dev)
        del lakes
        dev, level, U, N = sol['dev'], sol['level'], sol['U'], sol['N']
        new = geometry((height + CL.upsample(dev, FACTOR)).astype(np.float32), level, grid, inp['bodies'], log)
        moved = max(abs((new['rims'][k][0] - level) - (geo['rims'][k][0] - level)) for k in UPPER)
        rounds.append(dict(iterations=len(sol['history']), last_change_m=sol['history'][-1],
                           level_change_m=round(level - level0, 3), rims_moved_m=round(moved, 3),
                           rims_m={k: round(new['rims'][k][0] - level, 2) for k in UPPER}))
        log(f'  round {len(rounds)}: {len(sol["history"])} iterations, common level {level - level0:+.2f} m, '
            f'rims moved {moved:.3f} m')
        geo = new
        if moved < RIM_TOLERANCE_M:
            settled = True
            break
    if not settled:
        # A rim that flips between two passes, each made the lower by the load of filling to the other, cycles
        # instead of settling; the rounds' last three rims give its range and the case is flagged.
        log(f'  {name}: the rims did not settle in {OUTER_ROUNDS} rounds; the last three give their range')
    eff = height + CL.upsample(dev, FACTOR)
    U_fine = CL.upsample(U, FACTOR)
    covered = eff < level
    for k in UPPER:
        covered |= geo['regions'][k]
    lakes_volume = float(grid['row_area'] @ lake_load(height, inp['cells'], inp['levels'], covered).sum(axis=1)) / CL.WATER_DENSITY
    del covered
    straits, mares = strait_seeds(height)
    out = measures(height, eff, level, level0, geo, grid, U_fine, straits, mares, lakes_volume, inp['atlas'])
    out.update(response=None if spec is None else dict(profile=spec[0], years=spec[1]), settled=settled, rounds=rounds,
               catchments=geo['catchments'], seconds=round(time.time() - started))
    if not settled:
        last = rounds[-3:]
        out['unsettled_rims_m'] = {k: [min(r['rims_m'][k] for r in last), max(r['rims_m'][k] for r in last)] for k in UPPER}
    grids = dict(U=U.astype(np.float32), N=N.astype(np.float32))
    if name == DESIGN:
        regions = [(geo['regions'][k], geo['rims'][k][0] - NW.BELOW_RIM_M) for k in UPPER]
        grids['depth_m'] = water_depth(eff, level, regions).astype(np.float32)
        for k in UPPER:
            grids[f'region__{k}'] = np.packbits(geo['regions'][k].ravel())
    del eff, U_fine
    return out, grids


def run(cases=None, nmax=CL.NMAX, log=print):
    started = time.time()
    inp = load_inputs(nmax, log)
    results = {}
    BULK.mkdir(parents=True, exist_ok=True)
    for name, spec in CASES.items():
        if cases and name not in cases:
            continue
        out, grids = redraw(name, spec, inp, log)
        slug = name.replace(',', '').replace(' ', '_')
        path = BULK / f'redrawn_surface_{slug}_nmax{nmax}.npz'
        np.savez_compressed(path, metadata=json.dumps(dict(case=name, level_m=out['level_change_m'] + inp['level0'],
                                                           shape=list(inp['height'].shape), fine_ppd=FINE_PPD,
                                                           transform_ppd=CL.TRANSFORM_PPD, nmax=nmax)), **grids)
        out['grids'] = dict(path=str(path), sha256=sha256(path))
        results[name] = out
        (BULK / f'redrawn_surface_cases_so_far_nmax{nmax}.json').write_text(
            json.dumps(results, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')
        log(f'  {name}: level {out["level_change_m"]:+.2f} m, water {out["water_volume_km3"]:,} km3 '
            f'({out["water_beyond_the_atlas_km3"]:+,} beyond the atlas), {out["seconds"]} s')
    return dict(atlas=inp['atlas'], level0=inp['level0'], cases=results, nmax=nmax, seconds=round(time.time() - started))


def write(result, out=OUT):
    inputs = [ROOT / ATLAS, ROOT / 'shared/constants.json', ROOT / 'research/studies/drainage_network/results/network.json',
              ROOT / 'geography/results/crustal_loading.json', ln.SOURCES]
    ignored = [ROOT / LAKES, RESPONSE] + [fetch_inputs.path(n) for n in ('ldem_16.img', 'ldem_4.img', 'jggrx_0420a_sha.tab')]
    rel = lambda p: str(Path(p).relative_to(ROOT)) if str(p).startswith(str(ROOT)) else str(p)
    meta = dict(
        schema=SCHEMA,
        producer=dict(files={f: sha256(ROOT / f) for f in CODE}),
        inputs={rel(p): sha256(p) for p in inputs}, ignored_inputs={rel(p): sha256(p) for p in ignored},
        evidence=('A calculation on committed products and published values: the sea-level equation with migrating '
                  'shorelines (Farrell and Clark 1976; Mitrovica and Milne 2003) on LOLA at 16 px/deg with spherical '
                  'harmonics to degree %d, the load Love numbers of geography\'s crustal-loading study (elastic, '
                  'compressible, every degree) with the geology lane\'s viscoelastic change for the present Moon '
                  '(degrees 2-50, incompressible) added, and the drainage network study\'s rims found on the deformed '
                  'heights. The water arrives at once in the response; the creep\'s flow laws are Earth olivine\'s. '
                  'Not validated against any observed loading on the Moon.' % result['nmax']),
        reading_rule=('Levels and surfaces are on the deformed Moon\'s equipotential, in metres against the atlas\'s '
                      'common level on the dry Moon\'s geoid (level_change_m) or the redrawn common level '
                      '(surface_above_the_common_level_m). The common level holds 28% of the surface as the atlas '
                      'does; the three upper seas stand 1 cm under their rims over all the land their water joins. '
                      'Read the design case, the central profile 100 million years after the water arrives, as the '
                      'settled surface the round builds on, barring erosion; the central and stiff billion-year cases '
                      'as the most the present temperatures could add; the soft case as the warm mantle\'s risk; the '
                      'elastic and 1,000-year cases as the first years and the first millennium. Sill depths are below '
                      'the common level, '
                      'four-connected 16 px/deg cells. The atlas case deforms nothing. Water beyond the atlas counts '
                      'the seas and the lakes against the atlas\'s seas at its level and the same lakes.'),
        settings=dict(share=SHARE, upper_seas=list(UPPER), rim_tolerance_m=RIM_TOLERANCE_M, below_rim_m=NW.BELOW_RIM_M,
                      fine_ppd=FINE_PPD, transform_ppd=CL.TRANSFORM_PPD, tolerance_m=CL.TOLERANCE_M,
                      cases={k: (None if v is None else dict(profile=v[0], years=v[1])) for k, v in CASES.items()},
                      design=DESIGN, lakes=LAKES),
        **result)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(meta, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')
    print(f'wrote {out} in {result["seconds"]} s')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--cases', nargs='*', help='a subset of ' + ', '.join(CASES))
    parser.add_argument('--nmax', type=int, default=CL.NMAX)
    parser.add_argument('--out', type=Path, default=OUT)
    args = parser.parse_args(argv)
    write(run(args.cases, args.nmax, log=lambda *a: print(*a, flush=True)), args.out)


if __name__ == '__main__':
    main()
