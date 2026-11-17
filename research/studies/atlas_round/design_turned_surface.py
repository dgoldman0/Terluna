"""The climate model's surface for the run on the design surface: the design routing's seas and lakes turned into the
Moon's settled frame and averaged onto the GCM grids, as turned_surface.py does for the final surface.

    ulimit -v 9000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.design_turned_surface   # about a minute

The run at sunlight 0.930 ran on the final surface (turned_surface.py), with the near-side sea whole at the common level
and the governed lakes of the earlier design climate's rain. Its own closed budget parts the near-side sea at its neck
and gives the lakes of its own rain (design_routing.py), so the design surface's land differs from the land that run
saw. The register (17 November) runs the design climate again on the design surface: every sea at its own level (the
near-side sea parted, its main part and its southern seas below the common level; the salt sea and Humboldtianum at
the common level; the upper seas full at their rims) and the design routing's lakes, each at its own level, on the
redrawn floor's 16 px/deg heights.

The frame is the settled turn under the design surface's own water (design_settled_turn.py, 6.72°), by the round's
rule that the run stands at the settled turn under the water its surface holds; the run at 0.930 stood at the turn
under the final water with the near-side sea whole (6.97°), and the two frames differ by about a degree.

It takes turned_surface.py's steps as they stand: the lakes surface on the atlas's 4 px/deg grid (coarse), the frame
from a rotation (frame_from_rotation), the turned files with the record in turned.py's form (write_surface) and the
count of the T21 land cells below the water level (lifted). The files are named by the routing product's content
(turned.turned_name), which the runner's experiment A28D reads; the record is results/design_turned_surface_<frame>.json,
beside the final surface's.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from research.studies.atlas_round import design_routing as DR
from research.studies.atlas_round import design_settled_turn as DS
from research.studies.atlas_round import final_drainage as FD
from research.studies.atlas_round import redrawn_drainage as RD
from research.studies.atlas_round import turned_surface as TS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CODE = ('research/studies/atlas_round/design_turned_surface.py', 'research/studies/atlas_round/turned_surface.py',
        'research/studies/atlas_round/design_routing.py', 'research/studies/atlas_round/design_settled_turn.py',
        'climate/gcm/turned.py', 'climate/gcm/boundary.py', 'research/studies/atlas_round/redrawn_drainage.py')


def record_path(frame_sha):
    return HERE / 'results' / f'design_turned_surface_{frame_sha[:16]}.json'


def lakes_surface(log=print):
    """The design routing's water on the redrawn floor in boundary.lakes_product's form on the 4 px/deg grid, with the
    seas' fraction and the record's details."""
    surf = RD.load_surface(log)
    routing = DR.load(log)
    meta = routing['meta']
    h = surf['h']
    level0 = float(surf['level0'])
    if abs(level0 - meta['common_level_m']) > RD.RIM_CHECK_M:
        raise SystemExit(f'the design routing\'s common level {meta["common_level_m"]:.3f} m is not the floor\'s '
                         f'{level0:.3f} m')
    c = TS.coarse(h, routing['sea'].reshape(h.shape), routing['lake_cell'], routing['lake_level_m'].astype(float),
                  [meta['levels_m'][k] for k in meta['seas']])
    rows, cols = c['water_fraction'].shape
    lat = 90.0 - (np.arange(rows) + 0.5) / (RD.dr.PPD / TS.FACTOR)
    lon = (np.arange(cols) + 0.5) / (RD.dr.PPD / TS.FACTOR)
    weights = np.cos(np.radians(lat))[:, None] * np.ones(cols)
    share = lambda a: float((a * weights).sum() / weights.sum())                # noqa: E731
    wet = c['water_fraction']
    routing_sha = DR.sha256(DR.PRODUCT)
    prod = dict(lat_deg=lat, lon_deg=lon, water_fraction=wet, mean_height_m=c['mean_height_m'],
                mean_water_depth_m=c['mean_water_depth_m'],
                water_surface_m=np.where(wet > 0, c['water_surface_m'], level0), meta=dict(level_m=level0),
                sha256=f'{RD.sha256(surf["grids"])[:16]}+{routing_sha[:16]}', drainage_sha256=routing_sha[:16],
                routing_run=meta['climate']['run'], sea_share=share(c['sea_fraction']),
                lake_share=share(wet - c['sea_fraction']))
    details = dict(seas=meta['seas'], levels_above_common_m={k: round(v - level0, 2) for k, v in meta['levels_m'].items()},
                   common_level_m=level0, lakes_on_seas=c['lakes_on_seas'],
                   water_at_or_below_ground=c['water_at_or_below_ground'], climate=meta['climate'],
                   routing_sha256=routing_sha, grids=str(surf['grids']), grids_sha256=RD.sha256(surf['grids']))
    return prod, c['sea_fraction'], details


def design_frame():
    """The settled turn under the design surface's water (design_settled_turn.py) as a frame in turned.read_frame's form."""
    st = json.loads(DS.OUT.read_text())
    return TS.frame_from_rotation(st['frame_today_to_turned'], st['cases']['final_water']['settled']['tilt_to_ecliptic_deg'],
                                  TS.STATE, RD.sha256(DS.OUT))


def build(log=print):
    started = time.time()
    prod, sea, details = lakes_surface(log)
    frame = design_frame()
    run_path, run_record = FD.frame_record()          # the final surface's record, which the run at 0.930 ran on
    turn = np.asarray(frame['rotation']) @ np.asarray(run_record['frame']['rotation_today_to_turned']).T
    code = {f: RD.sha256(ROOT / f)[:16] for f in CODE}
    inputs = {str(p.relative_to(ROOT)): RD.sha256(p) for p in (DS.OUT, DR.OUT, ROOT / RD.SURFACE)}
    source = str(DS.OUT.relative_to(ROOT))
    path = record_path(frame['sha256'])
    record = TS.write_surface(
        prod, sea, frame, source=source, out=TS.PRODUCTS, inputs=inputs, path=path,
        producer=dict(domain='research', study='atlas_round', files=code),
        describe=dict(
            producer=' '.join(f'{k} {v}' for k, v in code.items()),
            source_product=('the atlas round\'s redrawn floor (LOLA with the settled case\'s U - N) and its design '
                            'routing (research/studies/atlas_round/design_routing.py, design_routing_16ppd.npz)'),
            file_evidence=('The atlas round\'s design surface: the redrawn floor at the crust\'s settled state with every '
                           'sea at the level the design climate\'s closed budget gives it (the near-side sea parted at '
                           'its neck, the salt sea and Humboldtianum at the common level, the upper seas full at their '
                           f'rims) and the lakes of that budget\'s routing on climate run {prod["routing_run"]}\'s rain, '
                           'carried into the Moon\'s settled frame under that water '
                           '(research/studies/atlas_round/design_settled_turn.py); heights above the settled state\'s '
                           'geoid. A storage geometry on the routing, not shorelines shaped by erosion.'),
            evidence=('A calculation: the design routing\'s seas and lakes on the redrawn floor, averaged onto the '
                      'atlas\'s 4 pixel/degree grid from 16 pixel/degree cells, carried into the settled frame by the '
                      'atlas round\'s rotation and averaged onto the GCM grids with the dry ground\'s height as it '
                      'stands. Heights and levels stand above the settled state\'s geoid and are counted from the '
                      'redrawn common level.')))
    record.update(
        ignored_inputs={str(DR.PRODUCT): details['routing_sha256'], details['grids']: details['grids_sha256']},
        water=dict(common_level_m_above_geoid=round(details['common_level_m'], 3),
                   seas_above_common_level_m=details['levels_above_common_m'], seas=details['seas'],
                   lake_cells_on_seas=details['lakes_on_seas'],
                   wet_cells_at_or_below_their_ground=details['water_at_or_below_ground'],
                   climate=details['climate']),
        ground_below_the_water_level=TS.lifted(TS.PRODUCTS / record['surface']),
        frame_source=dict(product=source, state=TS.STATE,
                          meaning=('the central profile 100 million years after the water arrives, under the design '
                                   'surface\'s water')),
        against_the_final_surface=dict(record=str(run_path.relative_to(ROOT)),
                                       frame=run_record['frame'], shares=run_record['shares'], t21=run_record['t21'],
                                       frame_turn_deg=round(float(np.degrees(np.arccos(np.clip(
                                           (np.trace(turn) - 1.0) / 2.0, -1.0, 1.0)))), 4)),
        seconds=round(time.time() - started))
    path.write_text(json.dumps(record, indent=1) + '\n')
    log(f'{path.relative_to(ROOT)}: water {record["shares"]["turned"]["water"]:.5f} turned, '
        f'{record["shares"]["today"]["water"]:.5f} today (seas {record["shares"]["today"]["seas"]:.5f}); '
        f'T21 water {record["t21"]["turned"]["water_share"]:.4f}, mask land {record["t21"]["turned"]["mask_land_share"]:.4f}; '
        f'land cells raised to the water level {record["ground_below_the_water_level"]}; {record["seconds"]} s')
    return path, record


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    build(log=lambda *a: print(*a, flush=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
