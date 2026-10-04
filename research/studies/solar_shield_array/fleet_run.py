"""Propagate a simultaneous 3-D population with independent initial states."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp

from shared import constants as K
from shared.provenance import constants_used, constants_changed
from protection.dynamics.ephemeris import Ephemeris
from protection.dynamics.cycling import CyclingEnvironment
from protection.dynamics.fleet import initial_fleet, fleet_acceleration
from .cycling_search import ROOT, HERE, SCENARIO, CENTRAL
from .cycling_feedback import RAW_SOURCES

RUN = ROOT/'research/runs/solar_shield_array/fleet'
SOURCES = RAW_SOURCES+['protection/dynamics/fleet.py',
                       'research/studies/solar_shield_array/fleet_run.py']


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def identities():
    return {p: digest(ROOT/p) for p in SOURCES}


def environment(days):
    eph = Ephemeris(epoch=SCENARIO['epoch_tdb'])
    samples = eph.sample(np.arange(-3600, days*K.JULIAN_DAY+7200, 1800.))
    eph.close()
    return CyclingEnvironment(samples, .05, CENTRAL, 'filter_only')


def read_raw(name):
    path = RUN/(name+'.npz')
    meta = json.loads(path.with_suffix('.json').read_text())
    if meta['producer']['source_hashes'] != identities():
        raise ValueError('Fleet dynamics source changed')
    if constants_changed(meta['producer']['constants']) or digest(path) != meta['sha256']:
        raise ValueError('Fleet raw constants or trace hash changed')
    return np.load(path), meta


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--name', required=True)
    p.add_argument('--radii-km', type=float, nargs='+', default=[15000.])
    p.add_argument('--planes', type=int, default=12)
    p.add_argument('--phases', type=int, default=48)
    p.add_argument('--days', type=float, default=30.)
    p.add_argument('--side-km', type=float, default=10.)
    p.add_argument('--arrangement', choices=['planar', 'isotropic'], default='isotropic')
    p.add_argument('--phase-offset', type=float, default=0.)
    p.add_argument('--output-step', type=float, default=600.)
    p.add_argument('--max-step', type=float, default=1200.)
    p.add_argument('--rtol', type=float, default=2e-10)
    a = p.parse_args(); config = vars(a)
    if Path(a.name).name != a.name:
        raise ValueError('Name must be a filename stem')
    RUN.mkdir(parents=True, exist_ok=True)
    env = environment(a.days)
    initial, labels = initial_fleet(env, np.asarray(a.radii_km)*1000,
                                    a.planes, a.phases, a.arrangement, a.phase_offset)
    count = len(initial); identity = identities(); before = time.monotonic()
    print(json.dumps(dict(stage='propagating', count=count, **config)), flush=True)
    times = np.linspace(0, a.days*K.JULIAN_DAY, int(np.ceil(a.days*K.JULIAN_DAY/a.output_step))+1)
    def derivative(t, flat):
        state = flat.reshape(count, 6)
        acc = fleet_acceleration(env, t, state, 4*K.MOON_RADIUS, a.side_km*1000)
        return np.column_stack([state[:, 3:], acc]).ravel()
    sol = solve_ivp(derivative, (times[0], times[-1]), initial.ravel(), t_eval=times,
                    method='DOP853', rtol=a.rtol,
                    atol=np.tile([.001]*3+[1e-7]*3, count), max_step=a.max_step)
    if not sol.success:
        raise RuntimeError(sol.message)
    states = sol.y.T.reshape(len(times), count, 6)
    ranges = np.linalg.norm(states[..., :3], axis=-1)
    if ranges.min() <= 4.01*K.MOON_RADIUS or ranges.max() >= 300e6:
        raise ValueError('Member crossed the radial exclusion; do not publish as completed')
    path = RUN/(a.name+'.npz')
    np.savez_compressed(path, t=times, state=states, labels=labels, initial=initial)
    if identities() != identity:
        raise ValueError('Sources changed during propagation')
    meta = dict(schema='terluna.research.solar-shield-fleet-raw/1', config=config,
        producer=dict(source_hashes=identity, constants=constants_used(identity)),
        sha256=digest(path), completed=True, count=count,
        range_km=[float(ranges.min()/1000), float(ranges.max()/1000)],
        evaluations=sol.nfev, elapsed_s=time.monotonic()-before,
        force_scope='Independent ideal spectral centre-ray SRP; intersecting shadow tubes require removal or coupled-force repropagation',
        initial_state_sha256=hashlib.sha256(initial.tobytes()).hexdigest())
    path.with_suffix('.json').write_text(json.dumps(meta, indent=2)+'\n')
    print(json.dumps(meta), flush=True)


if __name__ == '__main__':
    main()
