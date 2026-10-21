"""Number, mass and projected-area distributions of snow in the aerial winter.

Run with one numerical worker: python -m atmosphere.snow.population.
The old aggregation product remains a baseline. This calculation retains its
population equations, varies the admitted geometry, and exports distributions
inside the cloud as well as the snow crossing the freezing level. It does not
turn an assumed turbulence intensity into an empirically known breakup law.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np

from atmosphere.snow import aggregation as ag, flakes as sf
from atmosphere.snow.habits import GEOMETRIES
from shared.provenance import constants_used

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'atmosphere/snow/results/population.json'
CACHE = ROOT / 'research/runs/snow_population'
SCHEMA = 'terluna.atmosphere.snow-population/1'
FILES = ('atmosphere/snow/population.py', 'atmosphere/snow/aggregation.py',
         'atmosphere/snow/habits.py', 'atmosphere/snow/flakes.py',
         'atmosphere/rain/drops.py', 'atmosphere/rain/targets.py',
         'atmosphere/radiative_convective/thermodynamics.py')
# Keep an explicit open-ended bin for the inherited model's uncontrolled rare
# tail. Its moments are retained; infinity is serialized as a null upper edge.
# Neither this bin nor any raw maximum is a physically admitted flake size.
SIZE_EDGES = np.r_[np.geomspace(1e-5, 10., 97), np.inf]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def histogram(diameter, concentration, mass, area, velocity, edges=SIZE_EDGES):
    """Bin extensive moments, preserving number, mass, area and settling flux.

    Input concentration is number/m3; area is projected m2 per particle. The
    representative diameter is concentration-weighted. Consumers must use the
    separately conserved area moment, not reconstruct it from that diameter.
    """
    d, n, m, a, v = [np.asarray(x, dtype=float) for x in (diameter, concentration, mass, area, velocity)]
    if not (d.shape == n.shape == m.shape == a.shape == v.shape):
        raise ValueError('Particle arrays must have the same shape')
    if any(np.any(~np.isfinite(x)) for x in (d, n, m, a, v)):
        raise ValueError('Particle arrays must be finite')
    if any(np.any(x < 0) for x in (d, n, m, a, v)):
        raise ValueError('Particle moments cannot be negative')
    if d.size and (d.min() < edges[0] or d.max() >= edges[-1]):
        raise ValueError(f'Distribution exceeds explicit size grid ({d.min():g}..{d.max():g} m); expand it, do not clip')
    bins = np.searchsorted(edges, d, side='right') - 1
    total = lambda w: np.bincount(bins, weights=w, minlength=len(edges)-1)
    return np.array([total(n), total(n*d), total(n*m), total(n*a), total(n*v), total(n*m*v)])


def distribution(moments, edges=SIZE_EDGES):
    """JSON representation; empty bins retain geometric centres and zero moments."""
    n, nd, nm, na, nv, nmv = np.asarray(moments)
    divide = lambda x: np.divide(x, n, out=np.zeros_like(n), where=n > 0)
    d = divide(nd)
    finite_upper = np.where(np.isfinite(edges[1:]), edges[1:], edges[:-1]*4)
    d[n == 0] = np.sqrt(edges[:-1] * finite_upper)[n == 0]
    return dict(size_edges_m=[float(x) if np.isfinite(x) else None for x in edges],
                open_ended_last_bin=bool(np.isinf(edges[-1])),
                diameter_m=d.tolist(), number_m3=n.tolist(),
                projected_area_m2=divide(na).tolist(), mass_kg_m3=nm.tolist(),
                number_flux_m2_s=nv.tolist(), mass_flux_kg_m2_s=nmv.tolist(),
                fall_speed_m_s=divide(nv).tolist(),
                totals=dict(number_m3=float(n.sum()), mass_kg_m3=float(nm.sum()),
                            projected_area_m_inv=float(na.sum()), number_flux_m2_s=float(nv.sum()),
                            mass_flux_kg_m2_s=float(nmv.sum())))


class CloudSamples:
    """Time-sample the actual super-particle population in temperature layers."""
    def __init__(self, world, geometry, top_c, layers=15, sample_s=900.0):
        self.world, self.geometry = world, geometry
        self.z = np.linspace(world.z_freeze, world.height_of(sf.T0 + top_c), layers + 1)
        self.mid = 0.5 * (self.z[:-1] + self.z[1:])
        self.states = [world.state(z) for z in self.mid]
        self.moments = np.zeros((layers, 6, len(SIZE_EDGES)-1))
        self.count = 0
        self.sample_s, self.next_time = sample_s, None

    def __call__(self, mass, number, height, time, dt):
        if self.next_time is not None and time + 1e-6 < self.next_time:
            return
        self.next_time = time + self.sample_s
        self.count += 1
        layer = np.clip(np.searchsorted(self.z, height, side='right') - 1, 0, len(self.mid)-1)
        for j, state in enumerate(self.states):
            take = layer == j
            if not take.any():
                continue
            m = mass[take]
            d = self.geometry.size(m)
            n = number[take] / (self.z[j+1] - self.z[j])
            v, _ = ag.speeds(d, m, self.world.gravity, state, self.geometry.area_ratio)
            self.moments[j] += histogram(d, n, m, self.geometry.area(d), v)

    def export(self):
        if not self.count:
            raise ValueError('No post-fill cloud samples')
        return [dict(altitude_m=float(z), temperature_k=float(s['t']), pressure_pa=float(s['p']),
                     depth_m=float(self.z[i+1]-self.z[i]), samples=self.count,
                     **distribution(self.moments[i] / self.count))
                for i, (z, s) in enumerate(zip(self.mid, self.states))]


def arriving_distribution(result, world, geometry):
    mass, number = result['masses'], result['numbers']
    d = geometry.size(mass)
    state = world.state(world.z_freeze)
    v, _ = ag.speeds(d, mass, world.gravity, state, geometry.area_ratio)
    # Convert a crossing flux into instantaneous concentration with each
    # particle's own terminal velocity, not a median or a single size.
    n = number / result['seconds'] / v
    return dict(altitude_m=world.z_freeze, temperature_k=sf.T0, pressure_pa=world.p_freeze,
                **distribution(histogram(d, n, mass, geometry.area(d), v)))


class ArrivalHistory:
    """Arrival flux during cloud startup, in fixed three-hour windows.

    This addresses how long the prescribed cloud must operate before its large
    flakes arrive. It does not simulate the cloud's decay or evaporation.
    """
    def __init__(self, world, geometry, window_s=10800.):
        self.world, self.geometry, self.window_s = world, geometry, window_s
        self.state = world.state(world.z_freeze)
        self.moments = {}
        self.maximum = {}

    def __call__(self, mass, number, time, dt):
        k = int(time // self.window_s)
        d = self.geometry.size(mass)
        v, _ = ag.speeds(d, mass, self.world.gravity, self.state, self.geometry.area_ratio)
        moment = histogram(d, number / v, mass, self.geometry.area(d), v)
        self.moments[k] = self.moments.get(k, np.zeros_like(moment)) + moment
        if (number > 0).any():
            self.maximum[k] = max(self.maximum.get(k, 0.), float(d[number > 0].max()))

    def export(self):
        if not self.moments:
            return []
        # The last, possibly incomplete window is excluded.
        return [dict(start_hours=k*self.window_s/3600., end_hours=(k+1)*self.window_s/3600.,
                     maximum_sampled_diameter_m=self.maximum.get(k),
                     **distribution(self.moments.get(k, np.zeros((6, len(SIZE_EDGES)-1))) / self.window_s))
                for k in range(max(self.moments))]


def large_tail(result, geometry):
    """Retain rare model outliers separately from the one-per-hour statistic."""
    d = geometry.size(result['masses'])
    n, m = result['numbers'], result['masses']
    live = n > 0
    total_mass = float(np.dot(m, n))
    return dict(maximum_sampled_diameter_m=float(d[live].max()) if live.any() else None,
                minimum_sampled_diameter_m=float(d[live].min()) if live.any() else None,
                thresholds=[dict(diameter_m=threshold,
                                 arrivals_per_m2_hour=float(n[d >= threshold].sum()*3600/result['seconds']),
                                 arriving_mass_share=float(np.dot(n[d >= threshold], m[d >= threshold])/total_mass)
                                 if total_mass else 0.) for threshold in (.15, .6, 1., 10.)],
                evidence='Raw model tail under extrapolated aggregate geometry and no large-aggregate fracture. '
                         'The largest sampled particle is not a physical upper limit or a validated flake.')


def run_case(spec, world, geometry, binding):
    collector = CloudSamples(world, geometry, spec['top_c'])
    arrivals = ArrivalHistory(world, geometry)
    started, next_report = time.monotonic(), 0.
    def progress(model_time, particles):
        nonlocal next_report
        if model_time >= next_report:
            print(f"  {model_time/3600.:.0f} model h, {particles} particles, "
                  f"{time.monotonic()-started:.0f} wall s", flush=True)
            next_report = model_time + 21600.
    result = ag.column(world, spec['top_c'], spec['crystals_per_litre'] * 1000,
                       spec['snowfall_set_mm_h'], geometry=geometry, observer=collector,
                       seed=spec['seed'], bins=spec.get('bins', 30), inject=spec.get('inject', 8),
                       max_step_s=spec.get('max_step_s', 30.),
                       fill=spec.get('fill', 1.5), gather=spec.get('gather', 1.5),
                       arrival_observer=arrivals, progress=progress)
    ledger = result['mass_ledger_kg_m2']
    denominator = ledger['injected'] + ledger['deposited']
    ledger['relative_residual'] = ledger['residual'] / denominator if denominator else 0.0
    return dict(**spec, geometry_parameters=geometry.to_dict(),
                summary=ag.summary(result, mass_to_size=geometry.size), large_tail=large_tail(result, geometry),
                mass_ledger_kg_m2=ledger,
                freezing_level=arriving_distribution(result, world, geometry), layers=collector.export(),
                startup_arrivals=arrivals.export(),
                inputs=binding.get('inputs', {}))


def cases():
    """Central inputs plus geometry, random-seed, Earth and resolution controls."""
    rows = []
    def add(world, rate, number=100., geometry='bf95_dmax_dendrites_area', seed=0, suffix='', **kw):
        rows.append(dict(id=f'{world}_{rate:g}_{number:g}_{geometry}_s{seed}{suffix}',
                         world=world, snowfall_set_mm_h=rate, crystals_per_litre=number,
                         geometry=geometry, seed=seed, top_c=-15., **kw))
    for rate in (0.3, 1., 3.):
        for number in (10., 100.):
            add('moon', rate, number)
    for rate in (1., 3.):
        add('moon', rate, geometry='legacy_bf95')
        add('moon', rate, geometry='bf95_dmax')
        add('earth', rate)
        for seed in (1, 2):
            add('moon', rate, seed=seed)
    # A coupled refinement of layers, particle sampling and step size tests
    # numerical sensitivity; it is not an empirical validation.
    add('moon', 3., suffix='_refined', bins=60, inject=12, max_step_s=15.)
    for geometry in ('bf95_dmax_needles_area', 'bf95_dmax_plates_area'):
        add('moon', 3., geometry=geometry)
    return rows


def analyse(selected=None):
    moon, binding = sf.moon_column()
    worlds = dict(moon=moon, earth=sf.earth_column())
    files = {p: sha(ROOT/p) for p in FILES}
    constants = constants_used(FILES)
    rows = []
    CACHE.mkdir(parents=True, exist_ok=True)
    for spec in cases():
        if selected and spec['id'] not in selected:
            continue
        identity = dict(spec=spec, files=files, inputs=binding['inputs'], constants=constants)
        key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
        cache = CACHE / f'{spec["id"]}.json'
        previous = json.loads(cache.read_text()) if cache.exists() else None
        if previous and previous['key'] == key:
            row = previous['result']
            print(f'cached {spec["id"]}', flush=True)
        else:
            print(f'computing {spec["id"]}', flush=True)
            row = run_case(spec, worlds[spec['world']], GEOMETRIES[spec['geometry']], binding)
            cache.write_text(json.dumps(dict(key=key, result=row), separators=(',', ':')) + '\n')
        rows.append(row)
        print(spec['id'], {k: round(v, 4) if isinstance(v, float) else v
                          for k, v in row['summary'].items()}, flush=True)
    return dict(schema=SCHEMA, runs=rows,
                producer=dict(domain='atmosphere', files=files, inputs=binding['inputs'],
                              constants=constants),
                evidence=('Conditional one-dimensional populations in a persistent liquid-saturated cloud. '
                          'The condensation supply is imposed; the same collision/fragmentation laws as the '
                          'earlier study remain. Geometry corrections, seed spread and a resolution control are '
                          'separate numerical evidence. No resolved turbulence or large-aggregate mechanical '
                          'fracture is introduced into these populations.'),
                reading_rule=('Each run retains number/mass/area moments in logarithmic maximum-dimension bins. '
                              'freezing_level is a crossing population converted to concentration by individual '
                              'fall speeds; layers are post-fill time-sampled concentrations in finite layers. '
                              'projected_area_m2 is per particle, weighted by number concentration; multiply by '
                              'number_m3 for area density. mass_kg_m3 is already a concentration. summary size '
                              'percentiles are weighted by arriving mass, not number; the hourly extreme is a '
                              'cumulative arrival-rate threshold, not a hard maximum. large_tail records raw rare outliers; '
                              'the final size bin is open above 10 m (null upper edge) and retains its full moments. '
                              'The export grid is not a physically admitted size range. No local snowfall '
                              'frequency or habit proportions are inferred. crystals_per_litre normalizes '
                              'the injected number flux using the fall speed of a 0.2 mm crystal; the injected '
                              '0.1–0.4 mm spectrum is not a prescribed size-resolved concentration boundary.'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', action='append', help='Case id; partial output remains in ignored run directory')
    args = parser.parse_args()
    product = analyse(args.case)
    if not product['runs']:
        raise ValueError('No matching cases')
    path = CACHE/'partial.json' if args.case else OUT
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(product, indent=1) + '\n')
    print(f'wrote {path.relative_to(ROOT)}', flush=True)


if __name__ == '__main__':
    main()
