"""The deep lakes' columns through time, for what their anoxic depths gather.

    python -m research.studies.lake_anoxia.profiles        # about six minutes on one core

The lake model (biosphere/ecology/lake.py) runs unchanged. Each deep crater lake of the waters' oxygen study (352 m,
three bands, three targets) starts from that study's own spun-up column, cached in research/runs/lake_oxygen/cases,
and runs its target's ten years of lunar months once more, which repeats the settled cycle; the lake taken to the
1,925 m deepest point that half the lake area exceeds runs five passes of the record, fifty years, from its start
oxygenated as if just overturned. Each run saves the column daily (temperature, the sixteen oxygen columns and their
debt of reduced substances, the water's age since surface contact) with the mixed layer's depth at every step to
research/runs/lake_anoxia/<case>.npz; run.py reads them.
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

from biosphere.ecology import lake as L
from research.studies.lake_oxygen import run as R

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = ROOT / 'research' / 'runs' / 'lake_anoxia'
TARGETS = ('lock', 'resort', 'hot')
SITES = (0, 1, 2)
DEEPEST_PASSES = 5
EVERY_DEEP = 12                       # profile every 12 steps of 7150 s, about a day


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def cached_case(st, name, **case):
    path = R.RUNS / 'cases' / f'{name}_{R.cache_key(st, name, **case)}.json'
    if not path.exists():
        raise SystemExit(f'{path} is missing; run research.studies.lake_oxygen.run first')
    return json.loads(path.read_text()), path


def capture(lake, forcing, state, every):
    steps = len(forcing)
    at = list(range(0, steps, every))
    out = lake.run(forcing, state, profiles_at=at)
    keys = sorted(out['profiles'])
    prof = out['profiles']
    return dict(step=np.array(keys), t=np.array([prof[k]['t'] for k in keys], dtype=np.float32),
                oxygen=np.array([prof[k]['oxygen'] for k in keys], dtype=np.float32),
                debt=np.array([prof[k]['debt'] for k in keys], dtype=np.float32),
                age=np.array([prof[k]['age'] for k in keys], dtype=np.float32),
                mixed_depth=out['mixed_depth'].astype(np.float32), surface_t=out['surface_t'].astype(np.float32),
                saturation=out['saturation'].astype(np.float32), mean_o2=out['mean_o2'].astype(np.float32))


def main(argv=None) -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    st = R.Setting()
    scen = R.scenarios()
    meta = dict(schema='terluna.research.lake-anoxia-profiles/1', cases={},
                producer=dict(files={str(p.resolve().relative_to(ROOT)): sha256(p)
                                     for p in (Path(__file__), Path(L.__file__), Path(R.__file__))},
                              inputs=st.inputs),
                scenarios=[s.name for s in scen])
    tic = time.time()
    for target in TARGETS:
        for site_i in SITES:
            saved, path = cached_case(st, 'deep', kind='deep', target=target, site_i=site_i, years=None)
            site = R.site_info(st, site_i)
            w, forcing, mu = R.forcing_for(st, target, site, 'lock', None, R.STEPS_PER_MEAN_DEEP)
            bsn = R.basin('deep', site)
            lake = L.Lake(bsn, R.physics(st, 'deep'), scen, st.air)
            lake.calibrate_production(forcing)
            state = R.state_from(saved['final_state'], bsn, scen)
            rec = capture(lake, forcing, state, EVERY_DEEP)
            name = f'deep_{target}_{site["band"]}'
            np.savez_compressed(OUT / f'{name}.npz', z=bsn.z, area=bsn.area, mu=mu.astype(np.float32),
                                dt=forcing.dt, pressure=forcing.pressure.astype(np.float32), **rec)
            meta['cases'][name] = dict(kind='deep', target=target, band=site['band'], from_case=path.name,
                                       steps=len(forcing), dt=forcing.dt, every=EVERY_DEEP,
                                       npz_sha256=sha256(OUT / f'{name}.npz'))
            print(f'  {name} done at {time.time() - tic:.0f} s', flush=True)
    for target in TARGETS:
        base, _ = cached_case(st, 'deep', kind='deep', target='lock' if target == 'lock' else target, site_i=0,
                              years=None)
        site = R.site_info(st, 0)
        w, forcing, mu = R.forcing_for(st, target, site, 'lock', None, R.STEPS_PER_MEAN_DEEP)
        bsn = R.basin('deep', site, R.AREA_MEDIAN_DEEPEST_M)
        lake = L.Lake(bsn, R.physics(st, 'deep'), scen, st.air)
        lake.calibrate_production(forcing)
        state = R.state_from(base['final_state'], bsn, scen)
        state.oxygen[:] = float(L.oxygen_saturation(state.t.mean(), R.DESIGN_PRESSURE_PA / 101325.0, st.salinity,
                                                    st.o2_fraction))
        state.debt[:], state.age[:], state.solute[:] = 0.0, 0.0, 0.0
        passes = []
        for k in range(DEEPEST_PASSES):
            passes.append(capture(lake, forcing, state, EVERY_DEEP * 30))
        rec = {key: np.concatenate([p[key] for p in passes]) for key in passes[0] if key != 'step'}
        rec['step'] = np.concatenate([p['step'] + k * len(forcing) for k, p in enumerate(passes)])
        name = f'deepest_{target}'
        np.savez_compressed(OUT / f'{name}.npz', z=bsn.z, area=bsn.area, dt=forcing.dt,
                            pressure=forcing.pressure.astype(np.float32), **rec)
        meta['cases'][name] = dict(kind='deepest', target=target, band=site['band'], passes=DEEPEST_PASSES,
                                   steps=len(forcing) * DEEPEST_PASSES, dt=forcing.dt, every=EVERY_DEEP * 30,
                                   npz_sha256=sha256(OUT / f'{name}.npz'))
        print(f'  {name} done at {time.time() - tic:.0f} s', flush=True)
    (OUT / 'profiles.json').write_text(json.dumps(meta, indent=1) + '\n')
    print('wrote', OUT / 'profiles.json', f'in {time.time() - tic:.0f} s')
    return 0


if __name__ == '__main__':
    sys.exit(main())
