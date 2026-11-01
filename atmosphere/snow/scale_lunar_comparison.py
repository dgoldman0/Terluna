"""Run or assemble the complete lunar forcing comparison, with at most two workers.

SCALE keeps its native particle physics except in the explicitly labelled
sticking-law sensitivity. The prescribed atmosphere, source and duration match
the in-house scenario. Matching those inputs does not equate the physical laws
or establish empirical validity. Failed and incomplete attempts remain visible.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse
import json
import hashlib
import math
import numpy as np
from atmosphere.snow import scale_lunar as l, scale_aggregation as boxes
from atmosphere.snow import event_column as custom, flakes

NATIVE_CASES=[
 dict(name='lunar_prescribed_p2',packets=2,capacity_per_layer=128),
 dict(name='lunar_prescribed_p4',packets=4,capacity_per_layer=128),
 dict(name='lunar_prescribed_p8',packets=8,capacity_per_layer=192),
 dict(name='lunar_prescribed_dt30',packets=4,dt=30.,capacity_per_layer=192),
 dict(name='lunar_prescribed_seed2',packets=4,seed=2,capacity_per_layer=128),
 dict(name='lunar_prescribed_z60',packets=4,layers=60,capacity_per_layer=64),
 dict(name='lunar_reservoir_p2',packets=2,mode='reservoir',capacity_per_layer=256),
 dict(name='lunar_thermal_p2',packets=2,mode='thermal',capacity_per_layer=256),
 dict(name='lunar_sticking_p4',packets=4,capacity_per_layer=128,collection='inhouse'),
 dict(name='lunar_sticking_dt15',packets=1,dt=15.,capacity_per_layer=128,collection='inhouse'),
]
CUSTOM_CASES=[dict(area_m2=area,seed=seed) for area in (.0004,.0016) for seed in (20261021,20261022)]
OUT=l.ROOT/'atmosphere/snow/results/scale_lunar_comparison.json'


def run_native(spec):
    d=l.CACHE/'runs'/spec['name']
    if d.exists():return l.analyze(d)
    d=l.configure(**spec);l.execute(d);return l.analyze(d)


def native_group(case):
    columns=case['horizontal_columns']
    if any(not x['packets'] for x in columns):return dict(status='insufficient_arrivals')
    fields=['snowfall_mm_h','mass_d50_m','mass_d90_m','mass_d99_m']
    values=np.array([[x['snowfall_mm_h'],*x['dimension_mass_quantiles_m'][1:]] for x in columns])
    return dict(status='complete',replicates=len(columns),mean=dict(zip(fields,values.mean(axis=0).tolist())),
       standard_error=dict(zip(fields,(values.std(axis=0,ddof=1)/math.sqrt(len(values))).tolist())),
       mean_mass_cdf=np.mean([x['cdf'] for x in columns],axis=0).tolist(),
       interpretation='Eight horizontally isolated columns. Standard errors describe this artificial ensemble; they exclude shared physical-model error.')


def comparison(a,b,label,screen=True):
    if a['status']!='complete' or b['status']!='complete':return dict(label=label,status='incomplete',passed=False)
    changes={k:b['mean'][k]/a['mean'][k]-1 for k in a['mean']}
    cdf=float(np.max(np.abs(np.array(a['mean_mass_cdf'])-b['mean_mass_cdf'])))
    return dict(label=label,status='complete',relative_changes=changes,maximum_mass_cdf_difference=cdf,
      passed=(cdf<=.05 and all(abs(x)<=.05 for x in changes.values())) if screen else None,
      screen='At most 5% flux/d50/d90/d99 change and 0.05 mass-CDF distance; descriptive practical screen, not confidence or physical validity.' if screen else 'Physical sensitivity or seed variation; no pass/fail criterion.')


def assemble():
    cases={};groups={};missing=[]
    for spec in NATIVE_CASES:
        name=spec['name'];d=l.CACHE/'runs'/name
        try:
            case=l.analyze(d)
            cases[name]=case;groups[name]=native_group(case)
        except (FileNotFoundError,ValueError) as e:
            missing.append(dict(name=name,error=str(e)))
    pairs=[('lunar_prescribed_p2','lunar_prescribed_p4','particle sampling, 2 to 4 packets/minute',True),
           ('lunar_prescribed_p4','lunar_prescribed_p8','particle sampling, 4 to 8 packets/minute',True),
           ('lunar_prescribed_p8','lunar_prescribed_dt30','60 to 30 s at fixed 8 packets/minute',True),
           ('lunar_prescribed_p4','lunar_prescribed_z60','30 to 60 layers at fixed column depth and injection sampling',True),
           ('lunar_prescribed_p4','lunar_prescribed_seed2','independent seed',False),
           ('lunar_prescribed_p2','lunar_reservoir_p2','finite vapor reservoir with temperature restored',False),
           ('lunar_reservoir_p2','lunar_thermal_p2','allow native latent warming and melting; no external cooling',False),
           ('lunar_prescribed_p4','lunar_sticking_p4','replace 0.1 ice collection efficiency by inherited temperature law',False),
           ('lunar_sticking_p4','lunar_sticking_dt15','sticking sensitivity, 60 to 15 s at fixed 4 packets/minute',True)]
    contrasts=[]
    for left,right,label,screen in pairs:
        r=comparison(groups.get(left,dict(status='incomplete')),groups.get(right,dict(status='incomplete')),label,screen)
        r.update(left=left,right=right);contrasts.append(r)
    custom_cases=[]
    for spec in CUSTOM_CASES:
        p=l.CACHE/'custom_comparison'/f"ballistic_{spec['area_m2']}_{spec['seed']}.json"
        if p.exists():custom_cases.append(json.loads(p.read_text()))
        else:missing.append(dict(name=p.name,error='not completed'))
    files=['atmosphere/snow/scale_lunar.py','atmosphere/snow/scale_lunar_hook.f90',
           'atmosphere/snow/scale_lunar_comparison.py','atmosphere/snow/scale_build.py',
           'atmosphere/snow/scale_collision_observer.py','atmosphere/snow/scale_column.py',
           'atmosphere/snow/scale_aggregation.py']
    product=dict(schema='terluna.atmosphere.scale-lunar-comparison/1',
      producer={p:l.base.digest(l.ROOT/p) for p in files},expected_native_cases=NATIVE_CASES,
      expected_custom_cases=CUSTOM_CASES,cases=cases,groups=groups,comparisons=contrasts,
      custom_ballistic_without_chipping=custom_cases,missing_or_failed=missing,
      original_product=dict(path='atmosphere/snow/results/event_column.json',sha256=l.base.digest(l.ROOT/'atmosphere/snow/results/event_column.json')),
      evidence='Full-duration numerical sensitivity and refinement calculations. Matched lunar forcing, different microphysical closures. No lunar observation, contact-fracture validation or predicted storm optics.',
      reading_rule='Use the eight-column mean/standard-error summaries for comparisons with original means of realization quantiles; pooled arrival quantiles are retained separately. Native packet refinement is distinct from custom individual-particle statistical-volume refinement. The timestep comparison holds injected packets per minute fixed while refining their release timing. Reservoir experiments start with 48.34 kg/m2 vapor and impose a source ceiling, rather than an instantaneous deposition ceiling. Native thermal experiment has no imposed radiative/dynamical cooling. Neither moisture allowance is an observed storm supply.')
    boxes.save(OUT,product)
    return product


def run_custom(spec):
    w,b=flakes.moon_column();root=l.CACHE/'custom_comparison';root.mkdir(exist_ok=True)
    p=root/f"ballistic_{spec['area_m2']}_{spec['seed']}.json"
    if p.exists():return json.loads(p.read_text())
    r=custom.column(w,**spec,max_step_s=60.,breakup=False,settling_scheme='ballistic',wall_limit_s=1800.)
    r['producer']=custom.identity(b);r['status']='complete';boxes.save(p,r);return r


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['run','assemble'])
    p.add_argument('--workers',type=int,choices=(1,2),default=1);p.add_argument('--case',action='append')
    args=p.parse_args()
    if args.action=='run':
        jobs=[(run_native,s) for s in NATIVE_CASES if args.case is None or s['name'] in args.case]
        if args.case is None:jobs += [(run_custom,s) for s in CUSTOM_CASES]
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            tasks={pool.submit(fn,spec):spec for fn,spec in jobs}
            for future in as_completed(tasks):
                try:future.result();print('complete',tasks[future],flush=True)
                except Exception as e:print('failed',tasks[future],repr(e),flush=True)
    result=assemble();print('completed native',len(result['cases']),'missing',result['missing_or_failed'])
