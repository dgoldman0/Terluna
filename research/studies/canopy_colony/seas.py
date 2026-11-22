"""The seven seas of the design surface as the colony study needs them: each sea's area by depth on the redrawn
floor, its latitudes, and the river water it receives.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.canopy_colony.seas   # about a minute

It reads the design routing (research/studies/atlas_round/design_routing.py's product on the research drive) and the
redrawn floor (redrawn_drainage.load_surface) and writes results/seas.json. A sea's river water is what the cells
off the seas send straight into it. The south sea is the salt sea (register, 16 November 2026); the other six are
fresh.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from research.studies.atlas_round import design_routing as DR
from research.studies.atlas_round import redrawn_drainage as RD

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'seas.json'
SCHEMA = 'terluna.research.canopy-colony-seas/1'
SALT_SEAS = ('south_pole_aitken',)
DEPTHS_M = (0, 10, 20, 50, 100, 150, 200, 250, 300, 400, 500, 600, 750, 850, 1000, 1250, 1500, 2000, 2500, 3000, 4000,
            5000, 6000, 7000)
YEAR_S = 31557600.0


def sha256(path) -> str:
    digest = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def sea_table(depth, area, lat):
    """One sea from its cells' depths (m), areas (m2) and latitudes: its area and volume, the area of water deeper than
    each of DEPTHS_M, and the latitudes it spans with the share of its area poleward of 45 degrees."""
    below = [dict(depth_m=z, area_km2=float(area[depth > z].sum() / 1e6)) for z in DEPTHS_M if (depth > z).any()]
    return dict(area_km2=float(area.sum() / 1e6), volume_km3=float((area * depth).sum() / 1e9),
                mean_depth_m=float((area * depth).sum() / area.sum()), deepest_m=float(depth.max()), below=below,
                latitudes=[float(lat.min()), float(lat.max())],
                mean_latitude=float((area * lat).sum() / area.sum()),
                share_poleward_of_45=float(area[np.abs(lat) > 45.0].sum() / area.sum()))


def build(log=print):
    started = time.time()
    rt = DR.load(log=log)
    meta = rt['meta']
    surf = RD.load_surface(log=lambda *a: None)
    h, lat, row_area = surf['h'], np.asarray(surf['lat'], float), np.asarray(surf['row_area'], float)
    rows, cols = h.shape
    names, levels = meta['seas'], meta['levels_m']
    sea = rt['sea'].reshape(rows, cols)
    area = np.broadcast_to(row_area[:, None], sea.shape)
    lats = np.broadcast_to(lat[:, None], sea.shape)
    receiver, discharge = rt['receiver'], rt['discharge_m3s']
    into = np.where(receiver != np.arange(receiver.size), rt['sea'][receiver], -1)
    is_mouth = (rt['sea'] < 0) & (into >= 0)
    seas = {}
    for k, name in enumerate(names):
        m = sea == k
        entry = sea_table((levels[name] - h[m]).astype(float), area[m], lats[m])
        rivers = float(discharge[is_mouth & (into == k)].sum())
        entry.update(level_m=float(levels[name]), salt=name in SALT_SEAS, rivers_m3s=rivers,
                     rivers_m_yr=rivers * YEAR_S / (entry['area_km2'] * 1e6))
        seas[name] = entry
        log(f'{name}: {entry["area_km2"]:.0f} km2, mean {entry["mean_depth_m"]:.0f} m, rivers {rivers:.0f} m3/s')
    return dict(
        schema=SCHEMA,
        producer=dict(module='research.studies.canopy_colony.seas', code_sha256=sha256(Path(__file__))),
        inputs={'research/studies/atlas_round/results/design_routing.json':
                sha256(ROOT / 'research/studies/atlas_round/results/design_routing.json')},
        evidence=('Measured on the design surface\'s routing and redrawn floor at 16 pixels per degree; the rivers are '
                  'the design climate\'s estimate as the routing carries them. Not a survey of any real water.'),
        reading_rule=('seas: each of the design surface\'s seven seas at its own level, with the area of water deeper '
                      'than each depth (below), the latitudes of its cells on the design surface, and the river water '
                      'that cells off the seas send straight into it (m3/s, and as metres of water over the sea a '
                      'year). salt marks the salt sea.'),
        climate=meta.get('climate'), levels_m=levels, seas=seas, seconds=round(time.time() - started))


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    out = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(f'wrote {OUT} in {out["seconds"]} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
