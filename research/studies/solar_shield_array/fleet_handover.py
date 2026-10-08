"""Minute-cadence coverage around a fleet's weakest sampled date."""
import argparse
import json
import numpy as np

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.fleet import fleet_acceleration, coverage_snapshot
from .fleet_run import read_raw, environment, identities, digest
from .cycling_search import HERE
from .cycling_analysis import compact_series_json
from .fleet_analyze import summarize_coverage


def hermite_state(t, times, states):
    i = min(int(np.searchsorted(times, t, side='right'))-1, len(times)-2)
    i = max(0, i)
    dt = times[i+1]-times[i]; s = (t-times[i])/dt
    q0, v0 = states[i, :, :3], states[i, :, 3:]
    q1, v1 = states[i+1, :, :3], states[i+1, :, 3:]
    q = (2*s**3-3*s*s+1)*q0+(s**3-2*s*s+s)*dt*v0+(-2*s**3+3*s*s)*q1+(s**3-s*s)*dt*v1
    v = ((6*s*s-6*s)*q0+(-6*s*s+6*s)*q1)/dt+(3*s*s-4*s+1)*v0+(3*s*s-2*s)*v1
    return np.column_stack([q, v])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--product', required=True)
    p.add_argument('--name', required=True)
    p.add_argument('--step', type=float, default=60.)
    p.add_argument('--hours', type=float, default=12.)
    a = p.parse_args(); path = HERE/'results'/a.product
    case = json.loads(path.read_text())['cases'][a.name]
    coarse = case['selections']['nominal']['records']
    worst = int(np.argmin([r['ray_coverage'] for r in coarse]))
    raw, meta = read_raw(a.name); t = raw['t']; states = raw['state']
    centre = case['coverage_times_days'][worst]*K.JULIAN_DAY
    begin = np.clip(centre-a.hours*1800, t[0], max(t[0], t[-1]-a.hours*3600))
    times = np.linspace(begin, min(begin+a.hours*3600, t[-1]), int(a.hours*3600/a.step)+1)
    env = environment(meta['config']['days']); side = meta['config']['side_km']*1000
    records = []
    for time in times:
        state = hermite_state(time, t, states)
        d = fleet_acceleration(env, time, state, 4*K.MOON_RADIUS, side, diagnostics=True)
        records.append(coverage_snapshot(state, d['normal'], env.at(time), side, 4*K.MOON_RADIUS,
            sun_count=32, rotation=.291))
    sources = {**identities(), 'research/studies/solar_shield_array/fleet_handover.py': digest(__file__),
        'research/studies/solar_shield_array/fleet_analyze.py': digest(HERE/'fleet_analyze.py'),
        'research/studies/solar_shield_array/cycling_analysis.py': digest(HERE/'cycling_analysis.py')}
    product = dict(schema='terluna.research.solar-shield-fleet-handover/1',
        producer=dict(source_hashes=sources, constants=constants_used(sources), input_product_sha256=digest(path),
                      input_product=a.product, raw_sha256=meta['sha256']),
        case=a.name, parameters=vars(a), times_days=(times/K.JULIAN_DAY).tolist(),
        summary=summarize_coverage(records, times), records=records,
        evidence='Actual simultaneous phase histories, cubic-Hermite positions and velocities, finite-Sun polygon coverage every minute',
        reading_rule='A local time-resolution audit of the nominal geometric probe; collision, mutual-force and finite-extent exclusions from the parent case still apply. Finite sampling provides no continuous coverage certificate.')
    (HERE/'results/fleet_handover.json').write_text(compact_series_json(product))
    print(json.dumps(product['summary']), flush=True)


if __name__ == '__main__':
    main()
