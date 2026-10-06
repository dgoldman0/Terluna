"""Kept ring bundle (stage 2c of relative_orbits.md).

The continuous rings of ring_bundle_run, with the radius step widened from
0.6 to 1.0 km, fly twelve orbits with continuous low-thrust keeping (ring_keeping.ContinuousKeeping): each ring's
representative is held to its starting semi-major-axis step, eccentricity and
plane relative to its neighbour toward the middle ring, compared at the same
phase; the middle ring flies free. Rings slide, so the along-track element
stays free.
Mutual shadows push the rings unequally by up to about a tenth of the sail
force as they slide, which opens kilometres of eccentricity per orbit, so the
keeping acts continuously with a time constant of thirty minutes under the
0.001 m/s^2 actuator limit. The runner records the velocity change each ring
spends per orbit and the same clearance and service-overlap checks as stage 2b.
The same keeping also flies three orbits under gravity alone, which separates
the tide's share of the keeping cost from the light's.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.ring_keeping_run
"""
from __future__ import annotations

import hashlib
import inspect
import json
import resource
import time
from pathlib import Path

import numpy as np

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics import ring_bundle as rb
from protection.dynamics.ephemeris import DEFAULT_KERNEL
from protection.dynamics.relative_orbit import rtn_frame
from protection.dynamics.ring_keeping import ContinuousKeeping, geometric_elements
from .ring_bundle_run import DESIGN as BUNDLE, FILES as BUNDLE_FILES, HOUR, SEED, GravityOnly, analyse, environment, initial

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUN = ROOT/'research/runs/solar_shield_array/ring_keeping'
OUT = HERE/'results/ring_keeping.json'
FILES = ['research/studies/solar_shield_array/ring_keeping_run.py', 'protection/dynamics/ring_keeping.py'] + BUNDLE_FILES
DYNAMICS = [f for f in FILES if not f.endswith('_run.py')]
DESIGN = dict(orbits=12, radius_step_m=1000., tau_s=1800., limit_m_s2=1e-3, sample_s=BUNDLE['sample_s'], rtol=BUNDLE['rtol'],
              atol=BUNDLE['atol'], max_step_s=BUNDLE['max_step_s'], copy_reach=BUNDLE['copy_reach'], gravity_only_orbits=3)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_key(reps, origin, label=''):
    """Checkpoint key from the settings and sources that shape each orbit's dynamics."""
    settings = {k: DESIGN[k] for k in ('radius_step_m', 'tau_s', 'limit_m_s2', 'sample_s', 'rtol', 'atol',
                                       'max_step_s', 'copy_reach')}
    settings.update({k: BUNDLE[k] for k in ('rings', 'per_ring', 'pitch_m', 'tilt_deg')})
    return hashlib.sha256(json.dumps(settings, sort_keys=True).encode()+reps.tobytes()+origin.tobytes()+
                          ''.join(digest(ROOT/f) for f in DYNAMICS).encode()+label.encode()).hexdigest()[:16]


def kept(env, reps, model, origin, orbits=None, label=''):
    """Propagate with continuous keeping; checkpoint each orbit; return samples, control use and period.

    A label names a model other than the full one; its checkpoints also record the control's radial,
    along-track and normal parts.
    """
    RUN.mkdir(parents=True, exist_ok=True)
    mid = len(reps)//2
    _, a0, n0, _ = geometric_elements(reps[mid], reps, origin)
    period = 2*np.pi/n0
    keeper = ContinuousKeeping(model, reps, mid, origin, DESIGN['tau_s'], DESIGN['limit_m_s2'])
    key = run_key(reps, origin, label)
    times, samples, spent, peak, cpu, t0, current = [], [], [], [], 0., 0., reps.copy()
    parts = []
    for k in range(DESIGN['orbits'] if orbits is None else orbits):
        path = RUN/(f'{label}_orbit_{k}_{key}.npz' if label else f'orbit_{k}_{key}.npz')
        t1 = (k+1)*period
        if path.exists():
            z = np.load(path, allow_pickle=False)
        else:
            start = time.process_time()
            at, _ = rb.propagate(keeper.acceleration, current, t0, t1, rtol=DESIGN['rtol'], atol=DESIGN['atol'],
                                 max_step=DESIGN['max_step_s'])
            grid = np.arange(t0, t1, 60.)
            grid = np.r_[grid, t1] if grid[-1] < t1 else grid
            vectors = [(keeper.control(x, y), y) for x, y in ((x, at(x)) for x in grid)]
            control = np.array([np.linalg.norm(c, axis=1) for c, _ in vectors])
            extra = {}
            if label:
                rtn = np.array([[np.abs(c[i]@rtn_frame(y[i])[j]) for j in range(3)] for c, y in vectors
                                for i in range(len(y))]).reshape(len(grid), len(reps), 3)
                extra = dict(spent_rtn=np.trapezoid(rtn, grid, axis=0))
            sample = np.arange(t0, t1, DESIGN['sample_s'])
            sample = np.r_[sample, t1] if sample[-1] < t1 else sample
            tmp = path.with_name(path.stem+'.tmp.npz')
            np.savez(tmp, t=sample, s=np.stack([at(x) for x in sample]), spent=np.trapezoid(control, grid, axis=0),
                     peak=control.max(axis=0), cpu=time.process_time()-start, **extra)
            tmp.rename(path)
            z = np.load(path, allow_pickle=False)
        times.append(z['t'] if k == 0 else z['t'][1:]); samples.append(z['s'] if k == 0 else z['s'][1:])
        spent.append(z['spent']); peak.append(z['peak']); cpu += float(z['cpu'])
        if 'spent_rtn' in z:
            parts.append(z['spent_rtn'])
        current, t0 = z['s'][-1], t1
    return (np.concatenate(times), np.concatenate(samples), np.array(spent), np.array(peak), period, cpu, key,
            np.array(parts))


def main():
    wall, cpu0 = time.monotonic(), time.process_time()
    env = environment(DESIGN['orbits']*2.0+2)
    states, ring, slot, ratio, built = initial(env, DESIGN['radius_step_m'])
    tilt = np.radians(BUNDLE['tilt_deg'])
    reps = states[slot == (BUNDLE['per_ring']-1)//2]
    lags = BUNDLE['pitch_m']/np.linalg.norm(reps[:, 3:], axis=1)
    model = rb.ContinuousRings(env, tilt, np.full(len(reps), ratio), lags, reach=DESIGN['copy_reach'])
    origin = env.at(0.)['positions']['sun']
    origin = origin/np.linalg.norm(origin)
    t, s, spent, peak, period, cpu_run, key, _ = kept(env, reps, model, origin)
    # The gravity-only model lives in the bundle runner, so its source names its checkpoints.
    label = 'gravity_'+hashlib.sha256(inspect.getsource(GravityOnly).encode()).hexdigest()[:8]
    _, _, g_spent, g_peak, _, g_cpu, g_key, g_parts = kept(env, reps, GravityOnly(env, tilt), origin,
                                                            DESIGN['gravity_only_orbits'], label)
    g_mean = g_spent.mean(axis=0)
    a0 = np.linalg.norm(reps.mean(axis=0)[:3])
    n, v = np.sqrt(K.MOON_GM/a0**3), np.sqrt(K.MOON_GM/a0)
    checks = analyse(env, t, s, np.arange(len(reps)), np.zeros(len(reps), int), tilt, n, v, copies=model.copies,
                     interior=lambda r, j: (r > 0) & (r < len(reps)-1))
    per_orbit = spent.sum(axis=0)/len(spent)
    out = dict(schema='terluna.research.ring-keeping/1',
               producer=dict(files={f: digest(ROOT/f) for f in FILES},
                             constants=constants_used([f for f in FILES if f.endswith('.py')]),
                             inputs={'results/joint_seed.json': digest(SEED), 'de440s.bsp': digest(DEFAULT_KERNEL)}),
               evidence=('The continuous rings of ring_bundle.json with a 1.0 km radius step with '
                         'continuous low-thrust keeping over twelve orbits '
                         'of each ring\'s semi-major-axis step, eccentricity vector and plane about its phase-matched '
                         'neighbour toward the middle ring, which flies free, measured with exact elements; time '
                         'constant thirty minutes, 0.001 m/s^2 limit. Control use is the time '
                         'integral of the control acceleration sampled every minute. Forces, copies and checks are '
                         'those of stage 2b; the interior clearance leaves out the two edge rings. The gravity-only '
                         'case keeps the same rings under gravity alone (finite-square quadrature, no light) and '
                         'splits its control into radial, along-track and normal parts, each integrated in '
                         'magnitude.'),
               reading_rule=('A kept bundle passes when no sample falls below 150 m and the service overlaps stay '
                             'positive; the velocity spent per ring per orbit is the per-tile keeping cost of that '
                             'ring. Ring 18 is the bundle\'s unshadowed top edge.'),
               design=DESIGN, bundle=BUNDLE, run_key=key, period_h=period/HOUR,
               keeping=dict(mean_m_s_per_orbit_by_ring=per_orbit.tolist(),
                            rings_mean_m_s_per_orbit=float(np.mean(np.delete(per_orbit, len(reps)//2))),
                            ring_orbit_max_m_s=float(spent.max()), peak_acceleration_m_s2=float(peak.max()),
                            top_ring_mean_m_s_per_orbit=float(per_orbit[-1]),
                            per_orbit_m_s=spent.round(6).tolist()),
               checks=checks, propagation_cpu_s=cpu_run,
               gravity_only=dict(orbits=DESIGN['gravity_only_orbits'], run_key=g_key, label=label,
                                 mean_m_s_per_orbit_by_ring=g_mean.tolist(),
                                 rings_mean_m_s_per_orbit=float(np.mean(np.delete(g_mean, len(reps)//2))),
                                 per_orbit_m_s=g_spent.round(6).tolist(), peak_acceleration_m_s2=float(g_peak.max()),
                                 radial_along_normal_m_s_per_orbit_by_ring=g_parts.mean(axis=0).round(6).tolist(),
                                 propagation_cpu_s=g_cpu),
               resources=dict(cpu_s=time.process_time()-cpu0, wall_s=time.monotonic()-wall,
                              max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, numerical_threads=1))
    OUT.write_text(json.dumps(out, indent=1)+'\n')
    print(json.dumps(dict(keeping={k: v for k, v in out['keeping'].items() if k in
                                   ('rings_mean_m_s_per_orbit', 'ring_orbit_max_m_s', 'peak_acceleration_m_s2',
                                    'top_ring_mean_m_s_per_orbit')},
                          clearance={k: v for k, v in checks['clearance'].items() if not k.startswith('series')},
                          interior={k: v for k, v in checks['interior_clearance'].items() if not k.startswith('series')},
                          gravity_only={k: v for k, v in out['gravity_only'].items() if k.startswith('rings_mean')},
                          service=checks['service'], resources=out['resources']), indent=1))


if __name__ == '__main__':
    main()
