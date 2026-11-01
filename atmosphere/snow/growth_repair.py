"""Run explicit limited-growth column controls beside the frozen donor reference.

This changes mass-space deposition transport, keeping the physical laws and
boundary forcing fixed. Published terrestrial component and population tests
are distinct from the generic Earth column here; see earth_benchmarks.md.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from atmosphere.snow import flakes as sf, repair
from atmosphere.snow.column_experiments import Kinetics, solve_column
from atmosphere.snow.habits import GEOMETRIES
from shared.provenance import constants_used

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT/'research/runs/snow_growth_repair'
OUT = ROOT/'atmosphere/snow/results/growth_repair.json'
FILES = ('atmosphere/snow/growth_repair.py', 'atmosphere/snow/column_experiments.py',
         'atmosphere/snow/limited_growth.py', 'atmosphere/snow/limited_growth_batch.py',
         'atmosphere/snow/sectional.py',
         'atmosphere/snow/aggregation.py', 'atmosphere/snow/habits.py',
         'atmosphere/snow/flakes.py', 'atmosphere/rain/drops.py',
         'atmosphere/rain/targets.py', 'atmosphere/radiative_convective/thermodynamics.py')


def identity(binding):
    def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    inputs={p:sha(ROOT/p) for p in binding.get('inputs',{})}
    if inputs != binding.get('inputs',{}):
        raise RuntimeError('Loaded column input identity changed')
    return dict(files={p:sha(ROOT/p) for p in FILES}, inputs=inputs, constants=constants_used(FILES))


def cases():
    base=dict(world='moon', top_c=-15., number_per_m3=1e5, snowfall_mm_h=3.,
              layers=30, mass_bins=96, minimum_diameter_m=2e-5,
              maximum_diameter_m=10., max_step_s=30., cfl=.3,
              fill=1.5, gather=1.5, breakup=True)
    return [dict(base,id='earth_limited_m144',world='earth',snowfall_mm_h=1.,mass_bins=144),
            dict(base,id='moon_limited_m96'),dict(base,id='moon_limited_m144',mass_bins=144)]


def run_case(spec):
    world,binding=sf.moon_column() if spec['world']=='moon' else (sf.earth_column(),{})
    before=identity(binding); geometry=GEOMETRIES['bf95_dmax_dendrites_area']
    start=time.monotonic(); next_report=0.
    def progress(t, steps):
        nonlocal next_report
        if t>=next_report:
            print(f"{spec['id']}: {t/3600:.1f} model h, {steps} steps, {time.monotonic()-start:.1f} wall s",flush=True)
            next_report=t+6*3600.
    args={k:v for k,v in spec.items() if k not in ('id','world')}
    raw=solve_column(world,geometry=geometry,kinetics=Kinetics('legacy'),
                     growth_scheme='limited',progress=progress,**args)
    if identity(binding)!=before:
        raise RuntimeError('Producer changed during run; no result admitted')
    masses=raw['masses']; d=raw['diameters_m']; n=raw['arriving_numbers_m2']; weights=masses*n
    total=weights.sum(); cdf=np.cumsum(weights)/total if total else np.zeros_like(weights)
    def quantile(p): return float(d[min(np.searchsorted(cdf,p),d.size-1)]) if total else None
    summary=dict(mass_median_m=quantile(.5),mass_d90_m=quantile(.9),mass_d99_m=quantile(.99),
                 snowfall_mm_h=float(total*3600/raw['seconds']),
                 arriving_mass_share_above={str(x):float(weights[d>=x].sum()/total) if total else 0.
                                           for x in (.005,.01,.15,.6,1.,10.)})
    return dict(id=spec['id'],specification=spec,producer=before,
                settings=raw['settings'],summary=summary,
                mass_ledger_kg_m2=raw['mass_ledger_kg_m2'],domain_adequacy=raw['domain_adequacy'],
                diagnostics=dict(raw['diagnostics'],wall_seconds=time.monotonic()-start),
                masses_kg=masses.tolist(),diameters_m=d.tolist(),arriving_numbers_m2=n.tolist(),
                gathering_seconds=raw['seconds'],
                evidence='Limited mass-space growth flux; unchanged inherited collision and physical laws. '
                         'Generic Earth/Moon forcing, not a reproduction of an observed storm or published population simulation.',
                interpretation=('Finite mass-domain test passed; independent grid/time/height convergence still required.'
                                if raw['domain_adequacy']['passed'] else
                                'Finite mass-domain test failed; surviving quantiles are diagnostic and unsuitable for snow or optics forecasts.'))


def assemble():
    rows=[]; missing=[]; references=[]
    for spec in cases():
        path=CACHE/(spec['id']+'.json')
        if not path.exists():
            missing.append(dict(id=spec['id'],status='not_performed')); continue
        row=json.loads(path.read_text()); binding=sf.moon_column()[1] if spec['world']=='moon' else {}
        if row['producer']!=identity(binding) or row['specification']!=spec:
            raise RuntimeError('Stale limited-growth control: '+spec['id'])
        rows.append(row)
    for case in ('earth_moderate_m144','moon_heavy_m96','moon_heavy_m144'):
        path=repair.CACHE/(case+'.json')
        if not path.exists(): continue
        row=json.loads(path.read_text()); binding=sf.moon_column()[1] if row['specification']['world']=='moon' else {}
        if row['producer']!=repair.identity(binding): raise RuntimeError('Stale checkpoint reference: '+case)
        references.append({key:row[key] for key in ('id','specification','producer','summary',
                                                   'mass_ledger_kg_m2','domain_adequacy','diagnostics')})
    result=dict(schema='terluna.atmosphere.snow-growth-repair/1',cases=rows,references=references,
                unperformed=missing,
                reading_rule='Compare matched Earth or lunar cases at the same mass grid and forcing. '
                             'The generic Earth column is a transfer control, not external validation. '
                             'Upper-domain losses are numerical exits, never a physical flake-size ceiling.',
                remaining=['Mass-grid, timestep and vertical convergence of any usable column distribution.',
                           'Matched published terrestrial population benchmark; see earth_benchmarks.md.',
                           'Evolving aggregate geometry, interception applicability, whole-aggregate fracture and moist trajectories.'])
    OUT.write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',choices=[x['id'] for x in cases()]); parser.add_argument('--assemble',action='store_true')
    args=parser.parse_args()
    if args.case:
        CACHE.mkdir(parents=True,exist_ok=True)
        row=run_case(next(x for x in cases() if x['id']==args.case))
        (CACHE/(args.case+'.json')).write_text(json.dumps(row,indent=2)+'\n')
        print(json.dumps({key:row[key] for key in ('id','summary','domain_adequacy','mass_ledger_kg_m2')}),flush=True)
    if args.assemble: print(json.dumps({'cases':[x['id'] for x in assemble()['cases']]}))

if __name__=='__main__': main()
