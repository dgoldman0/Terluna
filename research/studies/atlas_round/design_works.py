"""The design's fixed waters and works on the redrawn floor, kept apart from any climate and any lake rule.

    ulimit -v 9000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python3 -m research.studies.atlas_round.design_works   # under a minute

The design routing stands on things the floor and the register's decisions fix, whatever rain falls and however the
lakes are filled:

- the three upper seas' cells with each full at its rim, and the salt sea's and Humboldtianum's at the common level;
- every sea's cells at the common level, on which the budget is closed;
- the near-side sea's two parts beyond its neck, the level at which each cell joins each part's water, and the crest
  between them;
- the works: the cell past Moscoviense's natural pass where the weir's share for the salt sea runs on, the cells past
  Orientale's and Smythii-Marginis's passes where their overflows run on, and the put-ins and take-outs of the two
  branches along their cuts.

Until the corrected lake rule these were found each time by the routing modules themselves (redrawn_drainage.py's
geometry and design state, governed_drainage.py's branches, near_side_neck.py's setting) and kept by
sustained_rivers.py under the design routing record's hash, so a change of the routing lost them. This module keeps
them as a product of their own, atlas-round/design_works/design_works_16ppd.npz on the research drive, with a record,
results/design_works.json, that load() checks the product against.

The product is assembled from the pieces those modules built on 18 November for the design on the present rule:
sustained_rivers.design_structure's waters and works and lake_rule_balance.setting's cells at the common level and
parts. The cells do not depend on the climate or on the lake rule: the seas' cells follow from the floor and the
levels, the passes from the floor, and the cuts' routes from the floor and the levels the hollows fill to before
they spill. The module does not yet find them from the floor alone; that needs redrawn_drainage.geometry and
design_state, governed_drainage.branches_for and near_side_neck.setting run once more (about 15 minutes), and is
left for when the floor or the works change.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from research.studies.atlas_round import final_drainage as FD
from research.studies.atlas_round import lake_rule as LR
from research.studies.atlas_round import lake_rule_balance as LB
from research.studies.atlas_round import redrawn_drainage as RD
from research.studies.atlas_round import sustained_rivers as SR

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'design_works.json'
BULK = Path('/media/projectspace/terluna-research/atlas-round/design_works')
PRODUCT = BULK / 'design_works_16ppd.npz'
SCHEMA = 'terluna.research.atlas-round-design-works/1'
WATERS = LB.WATERS
FIXED = WATERS[:LB.MAIN]             # the waters whose cells and levels the design fixes
OPS = SR.OPS


def rel(path) -> str:
    return str(Path(path).relative_to(ROOT)) if str(path).startswith(str(ROOT)) else str(path)


def build(log=print):
    """The product from the kept structure and setting of the design on the present rule, and its record."""
    started = time.time()
    d = LR.design_inputs(log)
    st = LB.setting(d, log)
    with np.load(SR.STRUCTURE, allow_pickle=False) as z:
        structure = json.loads(str(z['metadata']))
        op_cell, op_kind, op_sign = z['op_cell'], z['op_kind'], z['op_sign']
    water = np.where(d['water'] < LB.MAIN, d['water'], -1).astype(np.int8)
    meta = dict(schema=SCHEMA + '-grid', grid=dict(pixels_per_degree=RD.dr.PPD, rows=d['h'].shape[0], cols=d['h'].shape[1]),
                common_level_m=d['level0'], cell_km=d['cell_km'], waters=list(WATERS), fixed=list(FIXED),
                levels_m=[float(v) for v in d['levels'][:LB.MAIN]], seas_at_the_common_level=list(LB.SIX),
                crest_m=st['meta']['crest'], crest_four_connected_m=st['meta']['join4'], ops=list(OPS),
                reading_rule=('cell = row * 5760 + column, rows from 90 N and columns from 0 E on today\'s Moon. water: '
                              'the index into fixed of the sea over a cell with the upper seas at their rims and the '
                              'salt sea and Humboldtianum at the common level, -1 elsewhere, each at levels_m above the '
                              'dry Moon\'s geoid. common: the index into seas_at_the_common_level of the sea over a cell '
                              'with every sea at the common level. part: 1 for the near-side sea\'s main part, 2 for '
                              'its southern part, 0 elsewhere; joined_main and joined_south the level at which a cell '
                              'joins each part\'s water. op_cell, op_kind (an index into ops) and op_sign: the works\' '
                              'put-ins (+1) and take-outs (-1), in the order the branches run.'))
    BULK.mkdir(parents=True, exist_ok=True)
    keep = op_kind != OPS.index('neck')                  # the neck's spill is found from the parts' levels each time
    np.savez_compressed(PRODUCT, metadata=json.dumps(meta), water=water, common=st['common'], part=st['part'],
                        joined_main=st['joined_main'], joined_south=st['joined_south'], lat=d['lat'],
                        op_cell=op_cell[keep], op_kind=op_kind[keep], op_sign=op_sign[keep])
    cols = d['h'].shape[1]
    record = dict(
        schema=SCHEMA, producer=dict(files={'research/studies/atlas_round/design_works.py': RD.sha256(Path(__file__))}),
        inputs={rel(ROOT / RD.SURFACE): RD.sha256(ROOT / RD.SURFACE)},
        assembled_from={str(SR.STRUCTURE): RD.sha256(SR.STRUCTURE), str(LB.SETTING): RD.sha256(LB.SETTING),
                        'structure_key': structure['key'], 'setting_key': st['meta']['key']},
        evidence=('The design\'s fixed waters and works on the redrawn floor, assembled from the structure the routing '
                  'modules built for the design on the present lake rule (sustained_rivers.design_structure, '
                  'lake_rule_balance.setting): cells that follow from the floor, the seas\' levels and the cuts chosen '
                  'on 16 November, and depend on no climate and no lake rule. Not yet found again from the floor '
                  'alone by this module.'),
        reading_rule='Levels in m above the redrawn common level; the product\'s own rule is in its metadata.',
        levels_above_common_m={k: round(float(v) - d['level0'], 2) for k, v in zip(FIXED, d['levels'])},
        crest_above_common_m=round(st['meta']['join4'] - d['level0'], 2),
        works=[dict(kind=OPS[k], sign=int(s), at=LR.place(int(c), cols, d['rotation']))
               for c, k, s in zip(op_cell[keep], op_kind[keep], op_sign[keep])],
        product=dict(path=str(PRODUCT), sha256=RD.sha256(PRODUCT)), seconds=round(time.time() - started))
    return record


def load(log=print):
    """The works as the routing takes them: `d` with the floor (h, lat, cell_km, level0), the fixed waters and their
    levels, the works' cells as (cell, sign, False) with their kinds, and the frame's rotation; `st` with the seas at
    the common level, the near-side sea's parts and joining levels, and the crest. The product must be the record's."""
    record = json.loads(OUT.read_text())
    if not PRODUCT.exists() or RD.sha256(PRODUCT) != record['product']['sha256']:
        raise SystemExit(f'{PRODUCT} is not the product {rel(OUT)} records')
    surf = RD.load_surface(log)
    with np.load(PRODUCT, allow_pickle=False) as z:
        meta = json.loads(str(z['metadata']))
        a = {k: z[k] for k in z.files if k != 'metadata'}
    h = surf['h']
    if abs(float(surf['level0']) - meta['common_level_m']) > RD.RIM_CHECK_M or list(h.shape) != [meta['grid']['rows'], meta['grid']['cols']]:
        raise SystemExit('the works were not kept on this floor')
    _, frame = FD.frame_record()
    levels = list(meta['levels_m']) + [meta['common_level_m']] * 2
    d = dict(h=h, lat=a['lat'], cell_km=float(meta['cell_km']), level0=float(meta['common_level_m']), water=a['water'],
             levels=levels, ops=[(int(c), float(s), False) for c, s in zip(a['op_cell'], a['op_sign'])],
             op_kind=[meta['ops'][k] for k in a['op_kind']],
             rotation=np.asarray(frame['frame']['rotation_today_to_turned'], float))
    st = dict(common=a['common'], part=a['part'], joined_main=a['joined_main'], joined_south=a['joined_south'],
              meta=dict(level0=float(meta['common_level_m']), crest=meta['crest_m'], join4=meta['crest_four_connected_m']))
    log(f'the design\'s waters and works from {PRODUCT}')
    return d, st


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    record = build(log=lambda *a: print(*a, flush=True))
    OUT.write_text(json.dumps(record, indent=1, ensure_ascii=False) + '\n')
    print(f'wrote {OUT} in {record["seconds"]} s: {record["levels_above_common_m"]}, {len(record["works"])} works\' cells')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
