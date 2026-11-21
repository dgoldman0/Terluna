"""The atlas round's settled turn on the final water: the Moon's orientation 100 million years after the water arrives,
under the load the round's surface holds.

    ulimit -v 9000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.settled_turn     # about a minute

The register (16 November) runs the round's design climate at the Moon's settled turn, the central profile 100 million
years after the water arrives, which the geology lane's orientation series put at 6.15° (geology/interior_response.py
with geology/orientation.py). That series loads the Moon with the atlas's water, the seas held at 28% on the dry floor
and the rain-fed lakes, through the crustal-loading study's rigid geoid change. The round's final water differs: the
seas redrawn at the settled state with the upper seas full at their rims, and the governed drainage's lakes, most of
the hollows below the common level dry. This module first reproduces the series' settled turn from the series' own
load, then finds the turn under the final water with the same response:

- the load: the final water's depth over each 16 px/deg cell (each sea at its level over its cells, each lake at its
  level), 1,000 kg/m3, block-averaged to the series' 4 px/deg grid, and its degree-2 inertia as orientation.py takes it;
- the response: orientation.evolve with each profile's k'_2 and tidal k_2 spectra from the interior response product,
  the load applied when the water arrives and the rotational and tidal bulge re-forming about the moving axes;
- the frame at the settled state, the central profile at 100 million years, for the turned surface.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

from geography import crustal_loading as CL
from geology import orientation as ori
from geology import viscoelastic as ve
from research.studies.atlas_round import governed_drainage as GD
from research.studies.atlas_round import redrawn_drainage as RD

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'settled_turn.json'
SCHEMA = 'terluna.research.atlas-round-settled-turn/1'
RESPONSE = ROOT / 'geology' / 'results' / 'interior_response.json'
SERIES = ROOT / 'geology' / 'results' / 'moon_orientation.json'
CRUSTAL_NPZ = Path('/media/projectspace/terluna-research/crustal-loading/crustal_loading_nmax359.npz')
CODE = ('research/studies/atlas_round/settled_turn.py', 'geology/orientation.py', 'geology/viscoelastic.py',
        'research/studies/atlas_round/redrawn_drainage.py', 'geography/crustal_loading.py')
TIMES = np.r_[0.0, np.logspace(0.0, 9.0, 91)]          # years, the series' own
SETTLED = ('central', 1.0e8)                           # the settled state
BOUNDS = {'stiff': 1.0e9, 'soft': 1.0e7}               # the first step's bounding cases
WATER_DENSITY = 1000.0


def spectra(name, response):
    """A profile's (elastic, amplitudes) of k'_2 and of the tidal k_2."""
    sp = response['spectra'][name]
    pair = lambda kind: (sp[kind]['elastic'], np.array(sp[kind]['amplitudes'], float))
    return pair('k_load'), pair('k')


def evolve(tensor, name, response):
    """The orientation through time under a load applied when the water arrives (orientation.evolve)."""
    k_load, k_tide = spectra(name, response)
    return ori.evolve(TIMES, [(lambda t, tensor=tensor: tensor, k_load)], k_tide, ve.RATES_PER_YEAR, bulge=True)


def at(rows, years):
    """The row at the time nearest `years` in log time."""
    return min(rows, key=lambda r: abs(math.log10(1.0 + r['years']) - math.log10(1.0 + years)))


def summary(row):
    return dict(years=row['years'], pole_deg=round(row['pole_deg'], 3), pole_toward_lon_deg=round(row['pole_toward_lon_deg'], 1),
                sub_earth_lat_deg=round(row['sub_earth_lat_deg'], 2), sub_earth_lon_deg=round(row['sub_earth_lon_deg'], 2),
                face_move_deg=round(row['face_move_deg'], 2), tilt_to_ecliptic_deg=round(row['tilt_to_ecliptic_deg'], 3))


def final_water(log=print):
    """The final water's load (kg/m2) on the 4 px/deg grid, with the grid's latitudes and longitudes and the mass of
    the seas and of the lakes."""
    surf = RD.load_surface(log)
    h = surf['h'].astype(float)
    rows, cols = h.shape
    with np.load(GD.GOVERNED_ROUTING, allow_pickle=False) as d:
        meta = json.loads(str(d['metadata']))
        sea, lake_cell, lake_level = d['sea'], d['lake_cell'], d['lake_level_m'].astype(float)
    depth = np.zeros(h.size)
    flat = h.ravel()
    for i, k in enumerate(meta['seas']):
        cells = np.flatnonzero(sea == i)
        depth[cells] = np.maximum(meta['levels_m'][k] - flat[cells], 0.0)
    seas_depth = depth.copy()
    depth[lake_cell] = np.maximum(lake_level - flat[lake_cell], 0.0)
    lakes_depth = depth - seas_depth
    factor = CL.FINE_PPD // CL.TRANSFORM_PPD
    lat = 90.0 - (np.arange(rows // factor) + 0.5) / CL.TRANSFORM_PPD
    lon = (np.arange(cols // factor) + 0.5) / CL.TRANSFORM_PPD
    load = lambda d: WATER_DENSITY * CL.block_mean(d.reshape(rows, cols), factor)
    return dict(all=load(depth), seas=load(seas_depth), lakes=load(lakes_depth), lat=lat, lon=lon,
                governed_sha256=RD.sha256(GD.GOVERNED_ROUTING))


def build(log=print):
    started = time.time()
    response = json.loads(RESPONSE.read_text())
    series = json.loads(SERIES.read_text())
    # the series' own load: the crustal-loading study's rigid geoid change with the seas at 28% and the rain-fed lakes
    grids = np.load(CRUSTAL_NPZ)
    lat0, lon0 = grids['lat_deg'], grids['lon_deg']
    own = ori.inertia_of_geoid(grids['rigid__share__N'].astype(float), lat0, lon0)
    reproduced = summary(at(evolve(own, SETTLED[0], response), SETTLED[1]))
    listed = at([dict(r, years=r['years']) for r in series['time_series'][SETTLED[0]]['seas_and_lakes']], SETTLED[1])
    log(f'the series\' settled turn {listed["pole_deg"]} deg, reproduced {reproduced["pole_deg"]} deg')
    water = final_water(log)
    lat, lon = water['lat'], water['lon']
    if not (np.allclose(lat, lat0) and np.allclose(lon, lon0)):
        raise SystemExit('the final water\'s grid is not the series\'')
    loads = {k: ori.inertia_of_load(water[k], lat, lon) for k in ('all', 'seas', 'lakes')}
    out = dict(final_water_mass_kg={k: v[1] for k, v in loads.items()})
    cases = {}
    for label, (tensor, _) in (('final_water', loads['all']), ('final_seas_without_the_lakes', loads['seas'])):
        rows = {name: evolve(tensor, name, response) for name in ('central', 'stiff', 'soft')}
        cases[label] = dict(settled=summary(at(rows['central'], SETTLED[1])),
                            first_years=summary(rows['central'][0]),
                            bounds={name: summary(at(rows[name], years)) for name, years in BOUNDS.items()},
                            central_through_time=[summary(at(rows['central'], y)) for y in (1.0, 1e3, 1e5, 1e7, 1e8, 1e9)])
        if label == 'final_water':
            frame = at(rows['central'], SETTLED[1])['frame']
        log(f'{label}: settled {cases[label]["settled"]["pole_deg"]} deg, bounds '
            f'{dict((k, v["pole_deg"]) for k, v in cases[label]["bounds"].items())}')
    the_series = series['time_series']
    out.update(
        schema=SCHEMA, producer=dict(files={f: RD.sha256(ROOT / f) for f in CODE}),
        inputs={str(p.relative_to(ROOT)): RD.sha256(p) for p in (RESPONSE, SERIES, RD.OUT, GD.OUT)},
        ignored_inputs={str(CRUSTAL_NPZ): RD.sha256(CRUSTAL_NPZ), str(GD.GOVERNED_ROUTING): water['governed_sha256']},
        evidence=('A calculation on products: the geology lane\'s orientation method (geology/orientation.py: degree-2 '
                  'inertia, principal axes, the bulge re-forming through the tidal k_2) with its viscoelastic response, '
                  'applied to the round\'s final water. It reproduces the lane\'s own settled turn from the lane\'s load. '
                  'Not validated against any observation; the creep\'s flow laws are Earth olivine\'s.'),
        reading_rule=('Poles in degrees from today\'s, toward the longitude given; the sub-Earth point is the turned Moon\'s '
                      'point under Earth in today\'s coordinates. settled is the central profile 100 million years after '
                      'the water arrives; bounds are the stiff profile at a billion years and the soft at ten million, as '
                      'the first step bounds the surface. frame_today_to_turned is the settled frame for the turned '
                      'surface: its rows are the turned Moon\'s axes (toward Earth, the third, north) in today\'s frame.'),
        series_settled=dict(listed=listed, reproduced=reproduced,
                            stiff_and_soft_listed={name: summary(at(the_series[name]['seas_and_lakes'], y)) for name, y in BOUNDS.items()}),
        cases=cases, frame_today_to_turned=np.asarray(frame).tolist(), seconds=round(time.time() - started))
    return out


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
