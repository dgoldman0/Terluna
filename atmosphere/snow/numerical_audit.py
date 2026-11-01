"""Independent representation-area audit of the inherited snow population solver.

Run sequentially: python -m atmosphere.snow.numerical_audit
The imported simulator is never edited or monkeypatched. A new function sharing
its bytecode receives a private globals dictionary with a diagnostic collide
wrapper. This preserves the inherited equations and random sequence at A=1.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import time
import types
from pathlib import Path

import numpy as np

from atmosphere.snow import aggregation as ag, flakes as sf
from atmosphere.snow.habits import GEOMETRIES
from shared.provenance import constants_used

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'atmosphere/snow/results/numerical_audit.json'
SCHEMA = 'terluna.atmosphere.snow-numerical-audit/1'
FILES = ('atmosphere/snow/numerical_audit.py', 'atmosphere/snow/aggregation.py',
         'atmosphere/snow/habits.py', 'atmosphere/snow/flakes.py',
         'atmosphere/rain/drops.py', 'atmosphere/rain/targets.py',
         'atmosphere/radiative_convective/thermodynamics.py')
GEOMETRY = 'bf95_dmax_dendrites_area'
EDGES = np.r_[0., np.geomspace(1e-5, 10., 97), np.inf]
SPEC = dict(top_c=-2., number_per_m3=1e4, snowfall_mm_h=3., bins=4,
            inject=4, seed=0, fill=1., gather=1., max_step_s=30.)
SOURCE = dict(
    id='shima_2009', doi='10.1002/qj.441',
    url='https://doi.org/10.1002/qj.441',
    accessed_text='https://arxiv.org/html/physics/0701103v2',
    reading_extent='Sections 4.1.1, 4.1.4 and 5.1.3; not a full-paper validation.',
    relevant_claims=('Multiplicity is a positive integer; exhausted pairs use an integer-preserving split. '
                     'Coalescence efficiency enters the collision probability. The rigorous timestep '
                     'criterion is p much less than one; multiple collisions extend the numerical scheme.'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def private_column(collision):
    """Use the actual column bytecode with only a private collide binding."""
    namespace = dict(ag.column.__globals__, collide=collision)
    fn = types.FunctionType(ag.column.__code__, namespace, ag.column.__name__,
                            ag.column.__defaults__, ag.column.__closure__)
    fn.__kwdefaults__ = ag.column.__kwdefaults__
    return fn


class CollisionDiagnostics:
    """Replay only random draws needed to diagnose p and gt on an RNG copy.

    Counters describe sampled superparticle pairs, not real-particle encounter
    probabilities. No synthetic diagnostic draw reaches the simulator's RNG.
    """
    def __init__(self, area_m2):
        self.area = area_m2
        self.model_time = 0.
        self.calls = self.candidates = self.hits = self.multiple = self.capped = 0
        self.p_over_one = self.split_pairs = self.tiny_remainders = 0
        self.max_p = self.max_gt = self.max_gamma = 0.
        self.minimum_weight_per_m2 = None
        self.max_batch = None
        self.fragment_batches = self.fragment_batches_retained = 0
        self.fragments_made_per_m2 = self.fragments_discarded_per_m2 = 0.
        self.threshold_events = {}
        self.largest_sticking_batch = None
        self.small_remainder_examples = []

    def probe(self, m, xi, u, d, idx, volume, temperature, dt, rng, area_ratio):
        self.calls += 1
        n = idx.size
        if n < 2:
            return
        trial = copy.deepcopy(rng)
        order = trial.permutation(idx)
        j, k = order[:n//2], order[n//2:2*(n//2)]
        kernel = area_ratio * math.pi / 4. * (d[j] + d[k])**2 * np.abs(u[j]-u[k])
        p = np.maximum(xi[j], xi[k]) * kernel * dt / volume * (n*(n-1)/2)/(n//2)
        gamma = np.floor(p) + (trial.random(p.size) < p-np.floor(p))
        self.candidates += p.size
        self.p_over_one += int(np.count_nonzero(p > 1.))
        self.max_p = max(self.max_p, float(p.max()))
        self.max_gamma = max(self.max_gamma, float(gamma.max()))
        hit = gamma > 0
        self.hits += int(hit.sum())
        j, k, gamma, p = j[hit], k[hit], gamma[hit], p[hit]
        if not j.size:
            return
        swap = xi[j] < xi[k]
        a, b = np.where(swap, k, j), np.where(swap, j, k)
        limit = np.floor(xi[a]/xi[b])
        gt = np.minimum(gamma, limit)
        self.capped += int(np.count_nonzero(gamma > limit))
        self.multiple += int(np.count_nonzero(gt > 1))
        good = gt >= 1
        a, b, gt, p = a[good], b[good], gt[good], p[good]
        if not gt.size:
            return
        biggest = int(np.argmax(gt))
        if self.max_batch is None or gt[biggest] > self.max_gt:
            q = biggest
            self.max_gt = float(gt[q])
            self.max_batch = dict(model_time_h=self.model_time/3600., temperature_c=temperature-sf.T0,
                                  p=float(p[q]), gt=float(gt[q]),
                                  donor_weight_per_m2=float(xi[a[q]]/self.area),
                                  collector_weight_per_m2=float(xi[b[q]]/self.area),
                                  donor_mass_kg=float(m[a[q]]), collector_mass_kg=float(m[b[q]]),
                                  donor_diameter_m=float(d[a[q]]), collector_diameter_m=float(d[b[q]]))
        sticking = trial.random(a.size) < math.exp(.09*(temperature-sf.T0))
        rest = xi[a[sticking]]-gt[sticking]*xi[b[sticking]]
        self.split_pairs += int(np.count_nonzero(rest <= 0))
        scale = xi[a[sticking]]
        small = (rest > 0) & (rest <= 32*np.finfo(float).eps*scale)
        self.tiny_remainders += int(np.count_nonzero(small))
        if small.any() and len(self.small_remainder_examples) < 8:
            q = int(np.flatnonzero(sticking)[np.flatnonzero(small)[0]])
            self.small_remainder_examples.append(dict(model_time_h=self.model_time/3600.,
                donor_weight_per_m2=float(xi[a[q]]/self.area),
                collector_weight_per_m2=float(xi[b[q]]/self.area), gt=float(gt[q]),
                positive_remainder_per_m2=float((xi[a[q]]-gt[q]*xi[b[q]])/self.area)))
        if sticking.any():
            q = int(np.flatnonzero(sticking)[np.argmax(gt[sticking])])
            if self.largest_sticking_batch is None or gt[q] > self.largest_sticking_batch['gt']:
                self.largest_sticking_batch = dict(model_time_h=self.model_time/3600., gt=float(gt[q]),
                    p=float(p[q]), donor_weight_per_m2=float(xi[a[q]]/self.area),
                    collector_weight_per_m2=float(xi[b[q]]/self.area),
                    donor_mass_kg=float(m[a[q]]),collector_mass_kg=float(m[b[q]]))
        return dict(a=a,b=b,gt=gt,p=p,sticking=sticking,
                    ma=m[a].copy(),mb=m[b].copy(),xa=xi[a].copy(),xb=xi[b].copy())

    def __call__(self, m, xi, u, d, idx, volume, temperature, dt, rng, breakup, m_frag,
                 area_ratio=sf.AREA_RATIO, mass_to_size=sf.flake_size):
        actual_volume = volume*self.area
        batch = self.probe(m, xi, u, d, idx, actual_volume, temperature, dt, rng, area_ratio)
        made = ag.collide(m, xi, u, d, idx, actual_volume, temperature, dt, rng, breakup, m_frag,
                          area_ratio=area_ratio, mass_to_size=mass_to_size)
        # Inspect positive-weight packets after the actual, unchanged collision.
        # Threshold events diagnose numerical history; dimensions remain raw
        # extrapolations and are never admitted as physical aggregate sizes.
        live = idx[xi[idx] > 0]
        if live.size:
            post_d = mass_to_size(m[live])
            largest = int(np.argmax(post_d))
            particle = int(live[largest])
            for threshold in (.15,.6,1.,10.,1e3,1e6,1e9):
                if threshold not in self.threshold_events and post_d[largest] >= threshold:
                    event = dict(threshold_m=threshold,model_time_h=self.model_time/3600.,
                        temperature_c=temperature-sf.T0,diameter_m=float(post_d[largest]),
                        particle_mass_kg=float(m[particle]),weight_per_m2=float(xi[particle]/self.area),
                        represented_mass_kg_m2=float(m[particle]*xi[particle]/self.area),
                        layer_mass_kg_m2=float(np.dot(m[idx],xi[idx])/self.area))
                    if batch is not None:
                        involved = np.flatnonzero((batch['a']==particle)|(batch['b']==particle))
                        if involved.size:
                            q=int(involved[0])
                            event['pair_before_update']={k:float(batch[k][q]) for k in ('gt','p','ma','mb')}
                            event['pair_before_update'].update(
                                donor_weight_per_m2=float(batch['xa'][q]/self.area),
                                collector_weight_per_m2=float(batch['xb'][q]/self.area),
                                sticking=bool(batch['sticking'][q]),
                                observed_packet_was_collector=bool(batch['b'][q]==particle))
                    self.threshold_events[threshold]=event
        positive = xi[idx][xi[idx] > 0]/self.area
        if positive.size:
            value = float(positive.min())
            self.minimum_weight_per_m2 = value if self.minimum_weight_per_m2 is None else min(value, self.minimum_weight_per_m2)
        if made > 0:
            self.fragment_batches += 1
            self.fragment_batches_retained += int(made >= 1.)
            self.fragments_made_per_m2 += made/self.area
            if made < 1.:
                self.fragments_discarded_per_m2 += made/self.area
        return made

    def export(self):
        return dict(collision_calls=self.calls, candidate_pairs=self.candidates, positive_gamma_pairs=self.hits,
                    p_above_one_pairs=self.p_over_one, multiple_collision_pairs=self.multiple,
                    gamma_capped_pairs=self.capped, exhausted_sticking_pairs=self.split_pairs,
                    roundoff_scale_positive_remainders=self.tiny_remainders,
                    maximum_p=self.max_p, maximum_gamma=self.max_gamma, maximum_gt=self.max_gt,
                    first_threshold_events=list(self.threshold_events.values()),
                    largest_sticking_batch=self.largest_sticking_batch,
                    roundoff_scale_remainder_examples=self.small_remainder_examples,
                    minimum_positive_weight_per_m2=self.minimum_weight_per_m2, largest_batch=self.max_batch,
                    fragment_batches=self.fragment_batches, fragment_batches_retained=self.fragment_batches_retained,
                    fragments_made_per_m2=self.fragments_made_per_m2,
                    fragments_discarded_per_m2=self.fragments_discarded_per_m2,
                    reading_rule='Counts weight each sampled superparticle pair once. The remainder diagnostic '
                                 'flags 0 < residual <= 32 machine-epsilon times donor weight; it does not prove '
                                 'that a physically intended rare population should be removed.')


def scaled_column(world, area_m2=1., observer=None, arrival_observer=None, instrument=True,
                  progress=None, **kwargs):
    """Represent area A while preserving the physical density/supply inputs.

    The inherited column internally assumes one square metre. Scaling its input
    fluxes and its collision volume together emulates A square metres; observers
    and returned extensive quantities are normalized back to one square metre.
    """
    if not math.isfinite(area_m2) or area_m2 <= 0:
        raise ValueError('represented area must be finite and positive')
    diagnostics = CollisionDiagnostics(area_m2)
    def plain(m, xi, u, d, idx, volume, *args, **kw):
        return ag.collide(m, xi, u, d, idx, volume*area_m2, *args, **kw)
    runner = private_column(diagnostics if instrument else plain)
    args = dict(kwargs)
    args['number_per_m3'] *= area_m2
    args['snowfall_mm_h'] *= area_m2
    def cloud(m, xi, z, t, dt):
        if observer is not None:
            observer(m, xi/area_m2, z, t, dt)
    def arriving(m, xi, t, dt):
        if arrival_observer is not None:
            arrival_observer(m, xi/area_m2, t, dt)
    def reporting(t, n):
        diagnostics.model_time = t
        if progress is not None:
            progress(t,n,diagnostics)
    result = runner(world, observer=cloud, arrival_observer=arriving, progress=reporting, **args)
    result['numbers'] = result['numbers']/area_m2
    result['fragments'] /= area_m2
    result['mass_ledger_kg_m2'] = {k:v/area_m2 for k,v in result['mass_ledger_kg_m2'].items()}
    return result, diagnostics.export() if instrument else None


def moments(d, n, geometry):
    """Extensive number, mass and area moments in common open-ended bins."""
    q = np.searchsorted(EDGES, d, side='right')-1
    return np.array([np.bincount(q, weights=w, minlength=len(EDGES)-1)
                     for w in (n, n*geometry.mass(d), n*geometry.area(d))])


def packed(value):
    return dict(number_m3=value[0].tolist(), mass_kg_m3=value[1].tolist(),
                projected_area_m_inv=value[2].tolist(),
                totals={key:float(row.sum()) for key,row in zip(
                    ('number_m3','mass_kg_m3','projected_area_m_inv'), value)})


class Samples:
    def __init__(self, depth, geometry, interval_s=300.):
        self.depth, self.geometry, self.interval = depth, geometry, interval_s
        self.next_time = None
        self.count = 0
        self.moments = np.zeros((3,len(EDGES)-1))
        self.maximum_diameter = 0.

    def __call__(self, m, xi, z, t, dt):
        if self.next_time is not None and t < self.next_time:
            return
        self.next_time = t+self.interval
        self.count += 1
        d = self.geometry.size(m)
        if d.size:
            self.maximum_diameter = max(self.maximum_diameter, float(d.max()))
        self.moments += moments(d, xi/self.depth, self.geometry)

    def export(self):
        return dict(samples=self.count, interval_s=self.interval,
                    maximum_sampled_diameter_m=self.maximum_diameter,
                    **packed(self.moments/max(self.count,1)))


def arithmetic_examples():
    """Call the actual collision routine for an exhausted fractional pair."""
    class FixedRandom:
        def permutation(self, idx):
            return np.array(idx, copy=True)
        def random(self, size):
            return np.full(size, .5)
    pairs = []
    for area in (1.,2.,256.):
        m, xi = np.array([1e-9,2e-9]), np.full(2,1.5*area)
        d, u = np.array([.001,.002]), np.array([.1,.3])
        kernel = sf.AREA_RATIO*math.pi/4.*(d.sum())**2*abs(u[1]-u[0])
        before = float(np.dot(m,xi)/area)
        ag.collide(m,xi,u,d,np.arange(2),1.5*kernel*area,sf.T0,1.,FixedRandom(),False,1e-12)
        pairs.append(dict(area_m2=area, normalized_weights=(xi/area).tolist(),
                          resulting_masses_kg=m.tolist(), mass_before_kg_m2=before,
                          mass_after_kg_m2=float(np.dot(m,xi)/area), positive_packets=int((xi>0).sum())))
    fragments = [dict(area_m2=area, made_count=.5*area,
                      retained_count_per_m2=.5 if .5*area >= 1 else 0.,
                      discarded_count_per_m2=0. if .5*area >= 1 else .5) for area in (1.,2.,256.)]
    return dict(exhausted_pair=pairs, fragment_threshold=fragments,
                batch_sticking=dict(efficiency=.5, encounters=100,
                                    inherited_sticking_count_variance=100**2*.5*.5,
                                    independent_encounter_count_variance=100*.5*.5,
                                    qualification='Conditional counting example. Not an alternative solver or a '
                                                  'measurement of bias in the full population.'))


def compare(base, other):
    out = {}
    for location in ('cloud','freezing_level'):
        a,b = base[location],other[location]
        values = {}
        for key in ('number_m3','mass_kg_m3','projected_area_m_inv'):
            x,y = np.asarray(a[key]),np.asarray(b[key])
            sx,sy = float(x.sum()),float(y.sum())
            values[key] = dict(relative_total_change=(sy-sx)/sx if sx else None,
                               maximum_binned_cdf_difference=float(np.max(np.abs(np.cumsum(x)/sx-np.cumsum(y)/sy)))
                               if sx and sy else None,
                               l1_difference_over_baseline_total=float(np.abs(x-y).sum()/sx) if sx else None)
        out[location] = values
    return out


def run_case(world, area, breakup, spec=SPEC, progress=None):
    geometry = GEOMETRIES[GEOMETRY]
    depth = world.height_of(sf.T0+spec['top_c'])-world.z_freeze
    samples = Samples(depth,geometry)
    started = time.monotonic()
    result,diagnostics = scaled_column(world,area_m2=area,geometry=geometry,breakup=breakup,
                                       observer=samples,progress=progress,**spec)
    d = geometry.size(result['masses'])
    u,_ = ag.speeds(d,result['masses'],world.gravity,world.state(world.z_freeze),geometry.area_ratio)
    freezing = packed(moments(d,result['numbers']/result['seconds']/u,geometry))
    ledger = result['mass_ledger_kg_m2']
    ledger['relative_residual'] = ledger['residual']/(ledger['injected']+ledger['deposited'])
    return dict(id=f'area_{area:g}_breakup_{int(breakup)}',area_m2=area,breakup=breakup,
                wall_seconds=time.monotonic()-started,depth_m=depth,
                summary=ag.summary(result,mass_to_size=geometry.size),
                maximum_arriving_diameter_m=float(d.max()) if d.size else None,
                arriving_packet_count=int(d.size), mass_ledger_kg_m2=ledger,
                cloud=samples.export(),freezing_level=freezing,collision_diagnostics=diagnostics)


def production_context():
    """Snapshot existing production results as context; never run or modify them."""
    main = ROOT/'atmosphere/snow/results/population.json'
    paths = [main] if main.exists() else sorted((ROOT/'research/runs/snow_population').glob('*.json'))
    rows = []
    for path in paths:
        if path.name == 'partial.json':
            continue
        data = json.loads(path.read_text())
        records = data.get('runs',[data['result']] if 'result' in data else [])
        for row in records:
            tail = row.get('large_tail',{})
            maximum = tail.get('maximum_sampled_diameter_m')
            if maximum is not None:
                rows.append(dict(id=row['id'],source=str(path.relative_to(ROOT)),sha256=sha(path),
                                 maximum_arriving_diameter_m=maximum,
                                 tail_thresholds=tail.get('thresholds',[])))
    return dict(records=rows,reading_rule='Existing model outputs copied for context only; the short-column '
                'audit does not trace their extreme particles and cannot establish the cause of their tails. '
                'Ignored cache files may not be present in another checkout; their hashes identify the read bytes.')


def heavy_replay(world):
    spec=dict(SPEC,top_c=-15.,number_per_m3=1e5,bins=30,inject=8,fill=1.5,gather=1.5)
    path=ROOT/'research/runs/snow_population/moon_3_100_bf95_dmax_dendrites_area_s0.json'
    reference=json.loads(path.read_text())['result'] if path.exists() else None
    reference_identity=dict(path=str(path.relative_to(ROOT)),sha256=sha(path)) if reference else None
    started=time.monotonic()
    next_report=0.
    def progress(t,n,diag):
        nonlocal next_report
        if t>=next_report:
            print(f'  heavy {t/3600.:.0f} model h, {n} packets, '
                  f'{time.monotonic()-started:.1f} wall s, max gt={diag.max_gt:g}',flush=True)
            next_report=t+21600.
    row=run_case(world,1.,True,spec,progress=progress)
    row['id']='full_heavy100_unchanged_replay'
    row['physical_spec']=spec
    row['reference']=reference_identity
    if reference is not None:
        checks=dict(summary_exactly_equal=row['summary']==reference['summary'],
                    mass_ledger_exactly_equal=row['mass_ledger_kg_m2']==reference['mass_ledger_kg_m2'],
                    maximum_arriving_diameter_exactly_equal=row['maximum_arriving_diameter_m']
                    ==reference['large_tail']['maximum_sampled_diameter_m'])
        row['unchanged_replay_checks']=checks
        if not all(checks.values()):
            raise ValueError(f'Instrumented replay differs from cached reference: {checks}')
    return row


def analyse(areas=(1.,2.,256.), include_heavy=True):
    world,binding = sf.moon_column()
    rows=[]
    for breakup in (False,True):
        for area in areas:
            print(f'auditing area={area:g}, breakup={breakup}',flush=True)
            row=run_case(world,area,breakup)
            rows.append(row)
            print(f"  {row['wall_seconds']:.1f} s; {row['arriving_packet_count']} arrivals; "
                  f"max gt={row['collision_diagnostics']['maximum_gt']:.0f}",flush=True)
    comparisons=[]
    for breakup in (False,True):
        group=[r for r in rows if r['breakup']==breakup]
        baseline=next(r for r in group if r['area_m2']==1.)
        for other in group:
            if other is not baseline:
                comparisons.append(dict(baseline=baseline['id'],compared=other['id'],
                                        **compare(baseline,other)))
    heavy=heavy_replay(world) if include_heavy else None
    return dict(schema=SCHEMA,producer=dict(domain='atmosphere',files={p:sha(ROOT/p) for p in FILES},
                inputs=binding['inputs'],constants=constants_used(FILES)),
                evidence='Numerical audit of inherited exploratory, unvalidated column populations. '
                         'Arithmetic counterexamples establish nonhomogeneity in represented area; paired '
                         'short-column runs quantify sensitivity for one seed and one shallow prescribed cloud. '
                         'The unchanged full heavy replay, when included, traces collision arithmetic and '
                         'extreme-tail timing without revising the physical model.',
                reading_rule='All extensive particle numbers, mass ledgers and observers are normalized to '
                             'one square metre. Cloud moments are equally spaced post-fill samples averaged '
                             'over the whole prescribed depth. Freezing-level fluxes use individual settling '
                             'speeds to infer concentration. Number, mass and area are retained separately; '
                             'open size bins are bookkeeping, not physically admitted particle dimensions. '
                             'Paired-seed differences after representation-dependent branching include later '
                             'Monte Carlo divergence and are not an ensemble estimate of systematic bias.',
                control=dict(world='moon',geometry=GEOMETRY,physical_spec=SPEC,areas_m2=list(areas),
                             workers=1,instrumentation='Private column globals; independent deep-copied RNG probes.'),
                size_edges_m=[float(x) if np.isfinite(x) else None for x in EDGES],
                arithmetic=arithmetic_examples(),runs=rows,comparisons=comparisons,heavy_replay=heavy,
                production_context=production_context(),sources=[SOURCE])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--areas',type=float,nargs='+',default=[1.,2.,256.])
    parser.add_argument('--output',type=Path,default=OUT)
    parser.add_argument('--skip-heavy',action='store_true',help='Only the six short area controls')
    args=parser.parse_args()
    if 1. not in args.areas:
        parser.error('areas must include the A=1 baseline')
    product=analyse(tuple(args.areas),include_heavy=not args.skip_heavy)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(product,indent=1,allow_nan=False)+'\n')
    print(f'wrote {args.output}',flush=True)


if __name__=='__main__':
    main()
