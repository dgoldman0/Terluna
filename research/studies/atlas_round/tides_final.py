"""The monthly tide of the design surface's seas: geography's tide model (geography/tides.py) on the seven seas the
design surface holds, each at its own level, on the moisture map's quarter-degree sea pixels.

    ulimit -v 9000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.tides_final   # about five minutes

geography/tides.py computes the tide of the 28% atlas's seas, all at one level. On the design surface the near-side sea
is parted, the salt sea's shores are redrawn and the upper seas stand full at their rims, so the atlas's tide nodes miss
the upper seas and the new shores. This runs the same model on the design surface's seas (register, 17 November):

- the seas: the 4 px/deg pixels where most of a pixel lies under one sea of the design surface, labelled by that sea
  (moisture_final.map_grid on globe_surface.py's product, so the tide covers every sea pixel the moisture map's coast
  meets), each with its depth below its sea's own level from the mean ground of the pixel's cells under that sea;
- the response: tides.Response, each sea's mean level held and the water's self-attraction to degree 60, under the same
  Earth and Sun forcing, 2026-2045 every three hours (tides.forcing_coefficients), with tides.constituents' lines and
  tides.node_ranges' monthly ranges. The forcing's monomials are taken at each node's place on the turned Moon, in the
  frame of the climate the moisture map reads (final_drainage.frame_record), as near_side_neck.py takes them: Earth
  stands over the turned Moon as over today's, so a sea's tide follows its turned place. tides.py's own product is for
  today's orientation;
- stations at the named features of tides.STATION_FEATURES, each at the nearest sea pixel of the design surface; the
  dynamic check (tides.dynamic_check, its forcing turned the same way) on the near-side sea's main part, the salt sea
  and Smythii-Marginis. The seiche periods are left out.

The product (atlas-round/tides_final/tides_design_4ppd.npz on the research drive) takes tides.py's grid form (schema
terluna.geography.tides-grid/1), which the moisture map reads; the record is results/tides_final.json, with the atlas's
tide of the seas it replaces beside each sea's.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from geography import tides as T
from research.studies.atlas_round import final_drainage as FD
from research.studies.atlas_round import globe_surface as GS
from research.studies.atlas_round import near_side_neck as N

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'tides_final.json'
BULK = Path('/media/projectspace/terluna-research/atlas-round/tides_final')
PRODUCT = BULK / 'tides_design_4ppd.npz'
SCHEMA = 'terluna.research.atlas-round-tides/1'
CODE = ('research/studies/atlas_round/tides_final.py', 'research/studies/atlas_round/moisture_final.py',
        'research/studies/atlas_round/globe_surface.py', 'geography/tides.py', 'geography/lunar_ephemeris.py',
        'geography/tide_dynamics.py', 'research/studies/atlas_round/final_drainage.py',
        'research/studies/atlas_round/near_side_neck.py', 'shared/constants.json')
DYNAMIC_SEAS = (1, 3, 7)             # the near-side sea's main part, the salt sea and Smythii-Marginis
ATLAS_SEAS = {'433': (1, 2), '1417': (3,), '874': (7,), '213': (4,)}   # the atlas's seas the design surface's replace


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def design_seas(log=print):
    """The design surface's seas on the atlas's quarter-degree nodes in tides.load_seas's form, with each sea's area and
    name in tides.sea_bodies's form."""
    from research.studies.atlas_round import moisture_final as MF
    surface = MF.surface_grid()
    labels = MF.map_grid(surface).astype(np.int32)
    rows, cols = surface['rows'], surface['cols']
    k = rows // labels.shape[0]
    with np.load(GS.PRODUCT, allow_pickle=False) as z:
        h = z['height_m'].astype(float)
    body = surface['body'].reshape(rows, cols)
    level = surface['level'].reshape(rows, cols).astype(float)
    r, c = labels.shape
    blocks = lambda a: a.reshape(r, k, c, k)                                        # noqa: E731
    depth = np.zeros(labels.shape)
    for s in np.unique(labels[labels > 0]):
        under = blocks(body == s)
        n = under.sum(axis=(1, 3))
        mean_ground = np.where(under, blocks(h), 0.0).sum(axis=(1, 3)) / np.maximum(n, 1)
        sea_level = np.where(under, blocks(np.nan_to_num(level)), 0.0).sum(axis=(1, 3)) / np.maximum(n, 1)
        here = labels == s
        depth[here] = (sea_level - mean_ground)[here]
    atlas = T.load_seas()
    lat, lon = atlas['lat'], atlas['lon']
    if labels.shape != (lat.size, lon.size) or not np.allclose(lat, 90.0 - (np.arange(r) + 0.5) / 4.0) \
            or not np.allclose(lon, (np.arange(c) + 0.5) / 4.0):
        raise SystemExit('the design surface\'s sea pixels are not the atlas\'s quarter-degree nodes')
    area = (T.MOON_RADIUS / 1e3 * np.radians(0.25)) ** 2 * np.cos(np.radians(lat))[:, None] * np.ones((1, c))
    bodies = {int(s): dict(label=int(s), area_km2=float(area[labels == s].sum()), water_names=[GS.BODIES[int(s)]])
              for s in np.unique(labels[labels > 0])}
    log(f'{len(bodies)} seas over {int((labels > 0).sum()):,} quarter-degree nodes; shallowest '
        f'{float(depth[labels > 0].min()):.1f} m')
    return dict(labels=labels, depth=depth, lat=lat, lon=lon, level_m=None, source_sha256=sha256(GS.PRODUCT)), bodies


def basis_responses(response, rotation):
    """tides.basis_responses with each node's monomials taken at its place on the turned Moon (rotation today to turned,
    its rows the turned axes in today's frame): the response to each monomial forcing with and without the water's
    self-attraction, the monomials' values and the iterations."""
    la, lo = np.meshgrid(np.radians(response.grid['lat']), np.radians(response.grid['lon']), indexing='ij')
    values = T.monomial_values(T.unit_points(la, lo) @ np.asarray(rotation, float).T)
    with_sal, without, iterations = [], [], []
    for k in range(len(T.MONOMIALS)):
        h, n = response.solve(values[..., k])
        with_sal.append(h)
        without.append(response.solve(values[..., k], attraction=False)[0])
        iterations.append(n)
    return np.array(with_sal), np.array(without), values, iterations


def build(log=print):
    started = time.time()
    ephemeris_check = T.check_ephemeris()
    grid, bodies = design_seas(log)
    if (grid['depth'][grid['labels'] > 0] <= 0).any():
        raise SystemExit('a sea pixel stands no deeper than its ground')
    frame_path, frame_record = FD.frame_record()
    rotation = np.asarray(frame_record['frame']['rotation_today_to_turned'], float)
    response = T.Response(grid, bodies.keys())
    basis, plain, values, iterations = basis_responses(response, rotation)
    log(f'responses with self-attraction, {time.time() - started:.0f} s')
    time_days = np.arange(0.0, T.YEARS * 365.25, T.STEP_DAYS)
    coefficients = T.forcing_coefficients(T.JD_START + time_days)
    lines, _, amplitude, residual = T.constituents(time_days, coefficients, count=T.LINES)
    wet = response.wet
    ranges = T.node_ranges(coefficients, basis, wet, time_days)
    ranges_plain = T.node_ranges(coefficients, plain, wet, time_days)
    log(f'monthly ranges over {ranges["months"]} months, {time.time() - started:.0f} s')
    labels, area = grid['labels'][wet], response.area[wet]
    atlas = json.loads(T.RESULT.read_text())['seas']
    seas = {}
    for body in sorted(bodies, key=lambda b: -bodies[b]['area_km2']):
        m = labels == body
        w = area[m]
        seas[str(body)] = dict(
            names=bodies[body]['water_names'], area_km2=bodies[body]['area_km2'], nodes=int(m.sum()),
            monthly_range_median_m=dict(area_median=T.weighted_quantile(ranges['monthly_median'][m], w, .5),
                                        area_p95=T.weighted_quantile(ranges['monthly_median'][m], w, .95),
                                        maximum=float(ranges['monthly_median'][m].max())),
            extreme_range_m=dict(area_median=T.weighted_quantile(ranges['extreme'][m], w, .5),
                                 maximum=float(ranges['extreme'][m].max())),
            self_attraction_gain=float(np.sum(w * ranges['monthly_median'][m]) / np.sum(w * ranges_plain['monthly_median'][m])),
            atlas_seas={a: dict(names=atlas[a]['names'], area_km2=atlas[a]['area_km2'],
                                monthly_range_median_m=atlas[a]['monthly_range_median_m'])
                        for a, ours in ATLAS_SEAS.items() if body in ours and a in atlas})
    stations = []
    from geography import nomenclature
    named = nomenclature.by_name(nomenclature.features())
    sea_nodes = grid['labels'] > 0
    for names in T.STATION_FEATURES.values():
        for name in names:
            row, col, offset = T.nearest(grid, sea_nodes, named[name]['lon'], named[name]['lat'])
            stations.append(T.station_record(name, int(grid['labels'][row, col]), row, col, offset, grid, coefficients,
                                             basis, amplitude, lines, time_days))
    row, col, offset = T.nearest(grid, sea_nodes, *T.HEADLAND)
    stations.append(T.station_record('Eastern Smythii headland (wave shore studies)', int(grid['labels'][row, col]), row,
                                     col, offset, grid, coefficients, basis, amplitude, lines, time_days))
    kept, monomials = T.DYNAMIC_SEAS, T.monomial_values
    T.DYNAMIC_SEAS = tuple(b for b in DYNAMIC_SEAS if b in bodies)
    T.monomial_values = lambda points: monomials(np.asarray(points) @ rotation.T)          # the sea's turned places
    try:
        dynamics = T.dynamic_check(grid, amplitude, lines)
    finally:
        T.DYNAMIC_SEAS, T.monomial_values = kept, monomials
    log(f'dynamic check, {time.time() - started:.0f} s')
    line_records = []
    for line, amp in zip(lines, amplitude):
        on_seas = (amp @ basis[:, wet])
        line_records.append(dict(line, forcing_max_m=float(np.max(np.abs(values @ amp))),
                                 sea_response_max_m=float(np.max(np.abs(on_seas)))))
    BULK.mkdir(parents=True, exist_ok=True)
    rows, cols = np.nonzero(wet)
    np.savez_compressed(PRODUCT, lat_deg=grid['lat'], lon_deg=grid['lon'], rows=rows.astype(np.int32),
                        cols=cols.astype(np.int32), labels=labels.astype(np.int32),
                        basis=basis[:, wet].astype(np.float32), monomials=np.array(T.MONOMIALS),
                        line_period_days=np.array([l['period_days'] for l in lines]),
                        line_amplitude=(amplitude @ basis[:, wet]).astype(np.complex64),
                        monthly_range_median_m=ranges['monthly_median'].astype(np.float32),
                        extreme_range_m=ranges['extreme'].astype(np.float32),
                        metadata=np.array(json.dumps(dict(
                            schema='terluna.geography.tides-grid/1',
                            bodies={str(b): GS.BODIES[b] for b in bodies},
                            reading_rule=('Water-level response on every sea node of the atlas round\'s design surface '
                                          '(rows, cols into lat_deg, lon_deg; labels the sea, an index into '
                                          'globe_surface.BODIES). The level at time t is coefficients(t) @ basis, with '
                                          'coefficients from geography.tides.forcing_coefficients; line_amplitude gives '
                                          'each constituent\'s complex amplitude (m), Re[A exp(i w t)] with t in days '
                                          'from 2026-01-01 00:00 TT.')))))
    return dict(
        schema=SCHEMA, producer=dict(files={f: sha256(ROOT / f) for f in CODE}),
        inputs={str(p.relative_to(ROOT)): sha256(p) for p in (GS.OUT, T.CHECK, T.WATER, T.RESULT, frame_path)},
        ignored_inputs={str(GS.PRODUCT): sha256(GS.PRODUCT), str(T.ATLAS_GRID.relative_to(ROOT)): sha256(T.ATLAS_GRID)},
        product=dict(path=str(PRODUCT), sha256=sha256(PRODUCT)),
        evidence=('Equilibrium tide of the atlas round\'s design surface\'s seas from Earth and the Sun, 2026-2045, every '
                  'three hours, by geography/tides.py\'s model as it stands: degree-2 and degree-3 potential, lunar Love '
                  'numbers k2 and h2 at degree 2, each sea\'s volume held, and the water\'s self-attraction to '
                  'spherical-harmonic degree 60 with the solid Moon rigid under the load. A linear shallow-water solution '
                  'on the near-side sea\'s main part, the salt sea and Smythii-Marginis checks the equilibrium for the '
                  'eight leading lines. Each sea stands at its own level, so the seas\' tides are independent of one '
                  'another apart from their attraction.'),
        reading_rule=('As geography/results/tides.json: heights are metres of water level relative to the 2026-2045 mean; '
                      'the monthly range is the highest minus the lowest level within an anomalistic month (27.55 days), '
                      'its median taken over the months; the extreme range spans the whole record. Seas are keyed by '
                      'their index into globe_surface.BODIES (1 the near-side sea\'s main part, 2 its southern seas, 3 '
                      'the salt sea, 4 Humboldtianum, 5 Moscoviense, 6 Orientale, 7 Smythii-Marginis), each with the '
                      'atlas\'s seas it replaces; statistics weight nodes by area. Stations stand at the sea pixel '
                      'nearest each named feature, offset_km from it.'),
        settings=dict(grid_deg=0.25, harmonic_degree=T.DEGREE, love_numbers=dict(k2=T.MOON_LOVE_K2, h2=T.MOON_LOVE_H2),
                      degree3='Rigid solid Moon (Love factor 1)', water_density_kg_m3=response.water_density,
                      self_attraction_factor_degree2=float(response.sal_factor[2]),
                      self_attraction_iterations=iterations[1:], span_days=float(time_days[-1]), step_days=T.STEP_DAYS,
                      start_jd_tt=T.JD_START, pixels='moisture_final.map_grid: a sea where most of a pixel lies under one sea',
                      frame=dict(record=str(frame_path.relative_to(ROOT)), sha256=sha256(frame_path),
                                 pole_move_deg=frame_record['frame'].get('pole_move_deg'),
                                 sub_earth_point=frame_record['frame'].get('sub_earth_point'),
                                 meaning='the frame of the climate the moisture map reads; the forcing is taken at '
                                         'each node\'s turned place')),
        ephemeris_check=ephemeris_check,
        constituents=dict(lines=line_records, residual_rms_fraction=residual),
        seas=seas, stations=stations, dynamics=dynamics, neck_study_south_monthly_range_m=json.loads(
            N.OUT.read_text())['budgets']['closed_budget']['tide']['natural_neck']['south_monthly_range_m'],
        seconds=round(time.time() - started))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--out', type=Path, default=OUT)
    args = p.parse_args(argv)
    record = build(log=lambda *a: print(*a, flush=True))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(record, indent=1, ensure_ascii=False) + '\n')
    print(json.dumps({k: (v['names'][0], round(v['monthly_range_median_m']['area_median'], 3),
                          round(v['monthly_range_median_m']['maximum'], 3)) for k, v in record['seas'].items()}))
    print(f'wrote {args.out} in {record["seconds"]} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
