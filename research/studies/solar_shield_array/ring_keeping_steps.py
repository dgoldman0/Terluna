"""Radius steps of the kept rings at matched phase (diagnostic of stage 2c).

Rings slide past one another, so their representatives sit at different
phases of an orbit whose radius and osculating elements swing with the tide.
This reads the kept run's checkpoints and, at each mid-service, carries each
ring's inner neighbour along its orbit to the ring's own phase under the tide
(ring_keeping.gravity_shift), then records the radius step between them.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.ring_keeping_steps
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np

from shared.provenance import constants_used
from protection.dynamics.ring_keeping import gravity_shift
from .ring_bundle_run import environment, sun_angle
from .ring_keeping_run import DESIGN, OUT as KEPT, RUN

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE/'results/ring_keeping_steps.json'
FILES = ['research/studies/solar_shield_array/ring_keeping_steps.py', 'protection/dynamics/ring_keeping.py',
         'research/studies/solar_shield_array/ring_bundle_run.py', 'research/studies/solar_shield_array/ring_keeping_run.py',
         'protection/dynamics/cycling.py', 'protection/dynamics/ephemeris.py']


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def matched_steps(env, t, states):
    steps = []
    for m in range(len(states)-1):
        h = np.cross(states[m, :3], states[m, 3:])
        rate = np.linalg.norm(h)/np.dot(states[m, :3], states[m, :3])
        ahead = np.arctan2(np.dot(h/np.linalg.norm(h), np.cross(states[m, :3], states[m+1, :3])),
                           np.dot(states[m, :3], states[m+1, :3]))
        matched = gravity_shift(env, t, states[m], ahead/rate)
        steps.append(float(np.linalg.norm(states[m+1, :3])-np.linalg.norm(matched[:3])))
    return steps


def main():
    cpu = time.process_time()
    kept = json.loads(KEPT.read_text())
    env = environment(DESIGN['orbits']*2.0+2)
    files = sorted(RUN.glob(f"orbit_*_{kept['run_key']}.npz"), key=lambda p: int(p.name.split('_')[1]))
    rows = []
    for path in files:
        z = np.load(path, allow_pickle=False)
        t, s = z['t'], z['s']
        angles = np.array([sun_angle(env, tt, ss) for tt, ss in zip(t, s)])
        crossing = [i for i in range(len(t)-1) if angles[i] < 0 <= angles[i+1]]
        i = crossing[0] if crossing else len(t)//2
        steps = matched_steps(env, t[i], s[i])
        rows.append(dict(orbit=int(path.name.split('_')[1]), time_h=float(t[i]/3600.), steps_m=steps,
                         interior_min_m=float(min(steps[1:-1])), interior_max_m=float(max(steps[1:-1])),
                         bottom_edge_m=steps[0], top_edge_m=steps[-1]))
    out = dict(schema='terluna.research.ring-keeping-steps/1',
               producer=dict(files={f: digest(ROOT/f) for f in FILES},
                             constants=constants_used([f for f in FILES if f.endswith('.py')]),
                             inputs={'results/ring_keeping.json': digest(KEPT)}),
               evidence=('Radius steps between adjacent kept rings at one mid-service per orbit, each inner '
                         'neighbour carried to the outer ring\'s phase under the environment\'s gravity; read from '
                         'the kept run\'s checkpoints, which ring_keeping.json identifies by run key.'),
               reading_rule=('Designed step 1.0 km. Pairs 0-1 and 17-18 are the bundle\'s edges; the others are '
                             'interior.'),
               designed_step_m=DESIGN['radius_step_m'], orbits=rows, cpu_s=time.process_time()-cpu)
    OUT.write_text(json.dumps(out, indent=1)+'\n')
    for r in rows:
        print(f"orbit {r['orbit']:2d} h {r['time_h']:6.1f}: interior {r['interior_min_m']:6.0f}..{r['interior_max_m']:6.0f} m; "
              f"bottom edge {r['bottom_edge_m']:6.0f} m; top edge {r['top_edge_m']:6.0f} m")


if __name__ == '__main__':
    main()
