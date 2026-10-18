"""A 20-50 m lake held meromictic: what its isolated deep water does for fifty years at 301 K.

A closed lake turned soda, or a lake fed at depth by dense springs, keeps a denser layer below a chemocline that the
nights and the cold spells cannot mix through. The lake model carries one dissolved tracer that sets density by
layer, the solutes of respiration (lake.py's solute, at lake_oxygen's 0.39e-3 kg/m3 per g/m3 of oxygen respired);
here the deep layers start with an added step of it worth step_kg_m3 of density below step_depth_m. Above the
chemocline the model lets the solutes leave at the surface with the gas, so the mixed layer stays fresh as it would
under a fresh inflow at the surface, the way rivers hold Earth's ectogenic meromictic lakes; below it the step stays
unless the mixing entrains it. The column starts from the ladder's own lake of the same depth, band and trophic state,
oxygenated as if just overturned (columns.py's start), and runs five passes of the resort's record, saved daily as
columns.py saves its passes, so the same readers take it. Each case records its step so the step can be taken off
the solute to leave respiration's own.
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np

from biosphere.ecology import lake as L
from research.studies.lake_kinds import columns as C
from research.studies.lake_kinds import setting as S
from research.studies.lake_oxygen import run as R


def code_files():
    return C.code_files() + (Path(__file__),)


def key(st: S.Setting301, case: dict) -> str:
    h = hashlib.sha256()
    for p in code_files():
        h.update(C.sha256(p).encode())
    h.update(json.dumps(st.inputs, sort_keys=True).encode())
    h.update(json.dumps(dict(ground_k=st.ground_k, **case), sort_keys=True).encode())
    return h.hexdigest()[:20]


def run(st: S.Setting301, name: str, case: dict, start):
    """case: band, depth_m, coupled, steps_per_mean, passes, start (the ladder case it begins from), step_kg_m3 and
    step_depth_m. Saves and resumes pass by pass as columns.run does, and returns the same summary."""
    C.CASES.mkdir(parents=True, exist_ok=True)
    k = key(st, case)
    final = C.CASES / f'{name}_{k}.json'
    if final.exists():
        return json.loads(final.read_text())
    tic = time.time()
    site = R.site_info(st, case['band'])
    w, forcing, mu = R.forcing_for(st, 'resort', site, 'lock', None, case['steps_per_mean'])
    bsn = C.basin(site, case['depth_m'])
    scen = R.scenarios()
    names = [s.name for s in scen]
    ph = R.physics(st, 'deep', eta=case.get('eta', 'central'), coupled=case['coupled'])
    lake = L.Lake(bsn, ph, scen, st.air)
    lake.calibrate_production(forcing)
    every = max(1, int(round(C.PROFILE_S / forcing.dt)))
    step_solute = case['step_kg_m3'] / ph.solute_density             # g/m3 of the solute tracer
    below = bsn.centre > case['step_depth_m']
    done = []
    while (C.CASES / f'{name}_{k}_pass{len(done):02d}.json').exists():
        done.append(C.CASES / f'{name}_{k}_pass{len(done):02d}.npz')
    if done:
        state = C.state_from_dict(json.loads(done[-1].with_suffix('.json').read_text())['state'])
    else:
        state = R.state_from(start, bsn, scen)
        state.oxygen[:] = C.saturated(st, float(state.t.mean()))
        state.debt[:], state.age[:], state.solute[:] = 0.0, 0.0, 0.0
        state.solute[below] = step_solute
    spun = dict(passes=0, settled=None, start_c=None, step_solute_g_m3=step_solute,
                step_layers=int(below.sum()))
    for p in range(len(done), case['passes']):
        h0, o0 = L.heat_content(bsn, state), L.oxygen_content(bsn, state)
        rec, out = C.capture(lake, forcing, state, every)
        budget = dict(R.budget_errors(bsn, state, h0, o0, out, len(forcing) * forcing.dt),
                      produced=out['produced'].tolist(), respired=out['respired'].tolist(),
                      oxygen_in=out['oxygen_in'].tolist())
        series = {s: np.asarray(out[s], dtype=float) for s in C.SERIES + C.COLUMN_SERIES}
        path = C.CASES / f'{name}_{k}_pass{p:02d}.npz'
        np.savez_compressed(path, z=bsn.z, area=bsn.area, mu=mu, dt=forcing.dt,
                            pressure=forcing.pressure.astype(np.float32), every=every, **rec, **series)
        path.with_suffix('.json').write_text(json.dumps(dict(pass_index=p, budget=budget, spin_up=spun,
                                                             interfaces={i: int(out[i]) for i in C.INTERFACES},
                                                             state=C.state_dict(bsn, state),
                                                             npz_sha256=C.sha256(path))) + '\n')
        del rec, out, series
    passes = [np.load(C.CASES / f'{name}_{k}_pass{p:02d}.npz') for p in range(case['passes'])]
    metas = [json.loads((C.CASES / f'{name}_{k}_pass{p:02d}.json').read_text()) for p in range(case['passes'])]
    out = {s: np.concatenate([z[s] for z in passes]).astype(float) for s in C.SERIES + C.COLUMN_SERIES}
    mu_all = np.concatenate([z['mu'] for z in passes]).astype(float)
    res = R.summarise('deep', out, mu_all, bsn, scen, forcing.dt)
    res['deep_mixing'] = C.mixing_record(ph, metas, out['dissipation'])
    days = len(mu_all) * forcing.dt / 86400.0
    res['budget'] = dict(heat_relative_error=max(abs(m['budget']['heat_relative_error']) for m in metas),
                         heat_error_w_m2=max(abs(m['budget']['heat_error_w_m2']) for m in metas),
                         oxygen_relative_error=max(m['budget']['oxygen_relative_error'] for m in metas),
                         oxygen_error_g_m2_day=max(m['budget']['oxygen_error_g_m2_day'] for m in metas))
    res['case'] = dict(case, name=name, key=k, band_name=site['band'], depth_m=bsn.depth, layers=bsn.n,
                       area_km2=round(bsn.area[0] / 1e6, 1), mean_depth_m=round(float(bsn.volume.sum() / bsn.area[0]),
                                                                                2),
                       step_s=forcing.dt, profile_every_steps=every, years=round(days / 365.25, 2),
                       ground_k=st.ground_k, share_of_the_way=st.share['resort'], step_solute_g_m3=step_solute,
                       seconds=round(time.time() - tic, 1), scenarios=names)
    res['files'] = {f'{name}_{k}_pass{p:02d}.npz': metas[p]['npz_sha256'] for p in range(case['passes'])}
    res['final_state'] = metas[-1]['state']
    for z in passes:
        z.close()
    final.write_text(json.dumps(res) + '\n')
    return res
