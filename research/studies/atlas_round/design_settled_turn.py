"""The Moon's settled turn under the design surface's water: settled_turn.py's calculation with the design routing's seas
and lakes as the load.

    ulimit -v 9000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.design_settled_turn   # about a minute

The register (16 November) runs the round's design climate at the settled turn under the water its surface holds.
settled_turn.py found that turn, 6.97°, under the final water with the near-side sea whole at the common level. The
design surface (design_routing.py) parts the near-side sea at its neck, its main part 66.6 m and its southern seas
166.3 m below the common level, and carries the lakes of the run at 0.930's rain, so its water is lighter, mostly on
the near side, and the turn it settles to differs. This module runs settled_turn.build as it stands with only the load
swapped: each sea at its own level over its cells and each lake at its level, from the design routing's product, at
1,000 kg/m3 on the series' 4 px/deg grid. The record, results/design_settled_turn.json, takes settled_turn.json's form,
with that record's settled turn beside its own; design_turned_surface.py builds the climate model's surface in its
frame.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from geography import crustal_loading as CL
from research.studies.atlas_round import design_routing as DR
from research.studies.atlas_round import governed_drainage as GD
from research.studies.atlas_round import redrawn_drainage as RD
from research.studies.atlas_round import settled_turn as ST

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'design_settled_turn.json'
CODE = ('research/studies/atlas_round/design_settled_turn.py', 'research/studies/atlas_round/design_routing.py',
        *ST.CODE)


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def design_water(log=print):
    """The design surface's water as settled_turn.final_water gives the final water's: the load (kg/m2) on the 4 px/deg
    grid with the grid's latitudes and longitudes, the seas' and the lakes' apart."""
    surf = RD.load_surface(log)
    h = surf['h'].astype(float)
    rows, cols = h.shape
    routing = DR.load(log)
    meta = routing['meta']
    sea, lake_cell, lake_level = routing['sea'], routing['lake_cell'], routing['lake_level_m'].astype(float)
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
    load = lambda d: ST.WATER_DENSITY * CL.block_mean(d.reshape(rows, cols), factor)     # noqa: E731
    return dict(all=load(depth), seas=load(seas_depth), lakes=load(lakes_depth), lat=lat, lon=lon,
                governed_sha256=sha256(DR.PRODUCT))


def build(log=print):
    started = time.time()
    kept = ST.final_water
    ST.final_water = design_water
    try:
        out = ST.build(log)
    finally:
        ST.final_water = kept
    final = json.loads(ST.OUT.read_text())
    r_new, r_old = np.asarray(out['frame_today_to_turned']), np.asarray(final['frame_today_to_turned'])
    angle = lambda cosine: float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))     # noqa: E731
    rel = lambda p: str(Path(p).relative_to(ROOT))                                 # noqa: E731
    inputs = {k: v for k, v in out['inputs'].items() if k != rel(GD.OUT)}
    inputs.update({rel(DR.OUT): sha256(DR.OUT), rel(ST.OUT): sha256(ST.OUT)})
    ignored = {k: v for k, v in out['ignored_inputs'].items() if k != str(GD.GOVERNED_ROUTING)}
    ignored[str(DR.PRODUCT)] = sha256(DR.PRODUCT)
    out.update(
        producer=dict(files={f: sha256(ROOT / f) for f in CODE}), inputs=inputs, ignored_inputs=ignored,
        evidence=out['evidence'].replace('applied to the round\'s final water',
                                         'applied to the design surface\'s water (design_routing.py)'),
        water='the design surface: each sea at its own level (the near-side sea parted at its neck) and the design '
              'routing\'s lakes, from design_routing_16ppd.npz',
        against_the_final_water=dict(
            record=rel(ST.OUT), settled=final['cases']['final_water']['settled'],
            final_water_mass_kg=final['final_water_mass_kg'],
            frame_turn_deg=angle((np.trace(r_new @ r_old.T) - 1.0) / 2.0),
            pole_move_deg=angle(r_new[2] @ r_old[2]), sub_earth_move_deg=angle(r_new[0] @ r_old[0])),
        seconds=round(time.time() - started))
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--out', type=Path, default=OUT)
    args = p.parse_args(argv)
    record = build(log=lambda *a: print(*a, flush=True))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(record, indent=1) + '\n')
    print(json.dumps(dict(settled=record['cases']['final_water']['settled'],
                          against=record['against_the_final_water'])))
    print(f'wrote {args.out} in {record["seconds"]} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
