"""Factor gravity, cloud depth and water sharing in the pinned SCALE column.

This changes forcing only and reuses the identified native binaries. Diagonal
Earth/lunar controls preserve gH and the pressure/temperature path; the two
crossed controls deliberately freeze nonhydrostatic air to isolate mechanisms.
No case here is an observed Earth snowstorm or a validated lunar prediction.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import argparse
import json
import math
import numpy as np
from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from atmosphere.snow import scale_lunar as lunar
from atmosphere.snow import scale_lunar_comparison as comparison

ROOT = lunar.ROOT
OUT = ROOT / 'atmosphere/snow/results/scale_size_diagnosis.json'
DEPTH_RATIO = MOON_SURFACE_GRAVITY / STANDARD_GRAVITY
CASES = [
    dict(name='size_moon_baseline'),
    dict(name='size_earth_g_deep', gravity=STANDARD_GRAVITY),
    dict(name='size_moon_g_shallow', depth_factor=DEPTH_RATIO),
    dict(name='size_earth_g_shallow', gravity=STANDARD_GRAVITY, depth_factor=DEPTH_RATIO),
    dict(name='size_earth_g_shallow_dt15', gravity=STANDARD_GRAVITY, depth_factor=DEPTH_RATIO, dt=15., packets=1),
    dict(name='size_moon_no_aggregation', collisions=False),
    dict(name='size_earth_no_aggregation', gravity=STANDARD_GRAVITY, depth_factor=DEPTH_RATIO, collisions=False),
    dict(name='size_moon_seeds_tenth', source_factor=.1),
    dict(name='size_moon_water_ten', water_factor=10.),
    dict(name='size_moon_water_ten_dt15', water_factor=10., dt=15., packets=1),
    dict(name='size_earth_shallow_dt15_p3', gravity=STANDARD_GRAVITY, depth_factor=DEPTH_RATIO, dt=15., packets=3),
    dict(name='size_earth_shallow_dt15_p6', gravity=STANDARD_GRAVITY, depth_factor=DEPTH_RATIO, dt=15., packets=6),
    dict(name='size_earth_shallow_dt5_p1', gravity=STANDARD_GRAVITY, depth_factor=DEPTH_RATIO, dt=5., packets=1),
]
REFERENCES = ['lunar_prescribed_p4', 'lunar_sticking_dt15']


def configure(name, gravity=MOON_SURFACE_GRAVITY, depth_factor=1., source_factor=1.,
              water_factor=1., collisions=True, dt=60., packets=4, cache=lunar.CACHE):
    values = (gravity, depth_factor, source_factor, water_factor)
    if not all(math.isfinite(x) and x > 0 for x in values):
        raise ValueError('All physical factors must be finite and positive')
    d = lunar.configure(name, packets=packets, dt=dt, collisions=collisions,
                        capacity_per_layer=192, cache=cache)
    s = json.loads((d / 'settings.json').read_text())
    old_depth = s['depth_m']
    depth = old_depth * depth_factor
    dz = depth / s['layers']
    for fname in ('init.conf', 'run.conf'):
        text = (d / fname).read_text()
        changes = dict(CONST_GRAV=gravity, DZ=dz, BUFFER_DZ=dz)
        if fname == 'run.conf':
            changes.update(sdm_zupper=depth, sdm_sdnmlvol=s['dx_m']**2*dz)
        for key, value in changes.items():
            text = lunar.boxes.set_value(text, key, value)
        (d / fname).write_text(text)
    # The hook imposes identical T,p,e at identical fractional cloud heights.
    # Only sounding heights change; its pressure/theta/water values stay exact.
    lines = (d / 'env.txt').read_text().splitlines()
    for i in range(1, len(lines)):
        fields = lines[i].split()
        fields[0] = format(float(fields[0])*depth_factor, '.17g')
        lines[i] = ' '.join(fields)
    (d / 'env.txt').write_text('\n'.join(lines)+'\n')
    forcing = (d / 'lunar_forcing.txt').read_text().splitlines()
    values = list(map(float, forcing[1].split()))
    values[0] = depth
    values[2] *= source_factor
    values[3] *= water_factor / depth_factor
    forcing[1] = ' '.join(format(x, '.17g') for x in values)
    (d / 'lunar_forcing.txt').write_text('\n'.join(forcing)+'\n')
    s.update(depth_m=depth, gravity_m_s2=gravity,
             experiment=dict(depth_factor=depth_factor, source_factor=source_factor,
                             water_factor=water_factor, collisions=collisions,
                             gravity_ratio=gravity/MOON_SURFACE_GRAVITY,
                             hydrostatic_gH_ratio=gravity/MOON_SURFACE_GRAVITY*depth_factor,
                             imposed_number_flux_m2_s=values[2],
                             imposed_deposition_ceiling_kg_m2_s=values[3]*depth,
                             maximum_deposition_per_seed_kg=values[3]*depth/values[2]),
             forcing='Factorial mechanism control derived from original lunar forcing; same T,p,e versus fractional height. Constant seed number flux and integrated deposition ceiling except explicit source/water sensitivities.',
             diagnostic_producer_sha256=lunar.base.digest(Path(__file__)))
    lunar.boxes.save(d/'settings.json', s)
    lunar.boxes.save(d/'inputs.json', {p.name:lunar.base.digest(p) for p in d.iterdir() if p.is_file() and p.name!='inputs.json'})
    return d


def population(d, case):
    s = case['settings']
    a = np.fromfile(d/'lunar_arrivals.bin', dtype='>f8').reshape(-1, len(lunar.ARRIVAL_FIELDS))
    a = a[(a[:,0] >= s['gather_start_s']) & (a[:,8] == 10)]
    area = s['horizontal_cells'] * s['dx_m']**2
    period = s['duration_s'] - s['gather_start_s']
    count = a[:,2].sum()
    mass = (a[:,2]*a[:,3]).sum()
    nmono = (a[:,2]*a[:,7]).sum()
    injected = s.get('experiment', {}).get('imposed_number_flux_m2_s', s['original_settings']['injected_number_flux_m2_s'])
    ceiling = s.get('experiment', {}).get('imposed_deposition_ceiling_kg_m2_s', s['original_settings']['snowfall_mm_h']/3600)
    return dict(number_flux_m2_s=float(count/area/period),
        mass_flux_kg_m2_s=float(mass/area/period),
        mean_mass_kg=float(mass/count), mean_monomers=float(nmono/count),
        monomer_flux_m2_s=float(nmono/area/period),
        monomer_flux_to_injected=float(nmono/area/period/injected),
        number_survival_fraction=float(count/area/period/injected),
        mass_per_arriving_monomer_kg=float(mass/nmono),
        deposition_ceiling_per_injected_seed_kg=float(ceiling/injected),
        monomers_mass_quantiles=lunar.base.weighted_quantile(a[:,7], a[:,2]*a[:,3]),
        near_sphere_mass_fraction=float((a[:,2]*a[:,3])[(a[:,6]>.8)&(a[:,6]<1.25)].sum()/mass),
        maximum_arrived_dimension_m=float(a[:,4].max()))



def weighted_cdf_distance(x, weights_x, y, weights_y):
    """Exact empirical mass-CDF supremum at every observed diameter."""
    def ordered(a,w):
        a=np.asarray(a);w=np.asarray(w)
        if len(a)!=len(w) or not len(a) or np.any(w<0) or w.sum()<=0:
            raise ValueError('Invalid weighted population')
        order=np.argsort(a)
        return a[order],np.r_[0.,np.cumsum(w[order])/w.sum()]
    x,cx=ordered(x,weights_x);y,cy=ordered(y,weights_y)
    points=np.unique(np.r_[x,y])
    delta=np.abs(cx[np.searchsorted(x,points,side='right')]-cy[np.searchsorted(y,points,side='right')])
    i=int(np.argmax(delta))
    return dict(distance=float(delta[i]),diameter_m=float(points[i]))


def exact_column_cdf(left,right):
    """Compare equally weighted column CDFs, matching native_group's estimand."""
    def population_for(name):
        path=lunar.CACHE/'runs'/name
        settings=json.loads((path/'settings.json').read_text())
        a=np.fromfile(path/'lunar_arrivals.bin',dtype='>f8').reshape(-1,len(lunar.ARRIVAL_FIELDS))
        a=a[a[:,0]>=settings['gather_start_s']]
        w=a[:,2]*a[:,3]
        for cell in range(settings['horizontal_cells']):
            mask=a[:,1]==cell
            total=w[mask].sum()
            if total<=0:raise ValueError('Empty gathering column')
            w[mask]/=total*settings['horizontal_cells']
        return a[:,4],w
    return weighted_cdf_distance(*population_for(left),*population_for(right))


def assemble():
    results = {}
    missing = []
    for name in REFERENCES+[c['name'] for c in CASES]:
        d = lunar.CACHE/'runs'/name
        try:
            # Keep original reference analyses byte-identical; new cases analyzed here.
            c = json.loads((d/'analysis.json').read_text())
            results[name] = dict(group=comparison.native_group(c), population=population(d,c),
                arrivals=c['arrivals'], final=c['final'], observer=c['observer'],
                time_blocks=c['gathering_time_blocks'], settings=c['settings'],
                execution=c['execution'], water_relative_error=c['max_relative_water_residual'],
                input_sha256=json.loads((d/'inputs.json').read_text()),
                analysis_sha256=lunar.base.digest(d/'analysis.json'))
        except (FileNotFoundError, ValueError) as e:
            missing.append(dict(name=name, error=str(e)))
    pairs = [('lunar_prescribed_p4','size_moon_baseline','inactive capacity and random-stream control',True),('size_earth_g_shallow','size_moon_baseline','hydrostatically scaled Earth to lunar cloud',False),
             ('lunar_prescribed_p4','size_earth_g_deep','increase gravity at fixed deep profile',False),
             ('size_moon_g_shallow','size_earth_g_shallow','increase gravity at fixed shallow profile',False),
             ('size_moon_g_shallow','lunar_prescribed_p4','deepen cloud at lunar gravity',False),
             ('size_earth_g_shallow','size_earth_g_deep','deepen cloud at Earth gravity',False),
             ('size_earth_g_shallow','size_earth_g_shallow_dt15','Earth shallow timestep refinement',True),
             ('size_earth_g_shallow_dt15','size_earth_shallow_dt15_p3','Earth shallow injection sampling refinement',True),
             ('size_earth_shallow_dt15_p3','size_earth_shallow_dt15_p6','Earth shallow 12 to 24 packets per minute',True),
             ('size_earth_shallow_dt15_p3','size_earth_shallow_dt5_p1','Earth shallow 15 to 5 s at fixed 12 packets per minute',True),
             ('size_earth_no_aggregation','size_moon_no_aggregation','Earth to Moon without aggregation',False),
             ('size_moon_no_aggregation','lunar_prescribed_p4','enable lunar aggregation',False),
             ('lunar_prescribed_p4','size_moon_seeds_tenth','reduce lunar seed flux tenfold',False),
             ('lunar_prescribed_p4','size_moon_water_ten','increase lunar water ceiling tenfold',False),
             ('size_moon_water_ten','size_moon_water_ten_dt15','higher-water timestep and collision-cap control',True)]
    contrasts = []
    for left,right,label,screen in pairs:
        if left in results and right in results:
            r = comparison.comparison(results[left]['group'],results[right]['group'],label,screen)
            r.update(left=left,right=right)
            r['query_grid_mass_cdf_difference']=r['maximum_mass_cdf_difference']
            r['cdf_distance_method']='sparse_query_grid'
            if screen:
                exact=exact_column_cdf(left,right)
                r['exact_mass_cdf_difference']=exact['distance']
                r['maximum_mass_cdf_difference']=exact['distance']
                r['cdf_distance_method']='exact_equal_column_empirical'
                r['exact_mass_cdf_difference_at_m']=exact['diameter_m']
                r['query_grid_screen_passed']=r['passed']
                r['passed']=r['passed'] and exact['distance']<=.05
                r['screen']='At most 5% flux/d50/d90/d99 change and 0.05 exact empirical mass-CDF distance; equal weighting of eight column-normalized distributions. Earlier sparse-query distance retained separately.'
            contrasts.append(r)
    result = dict(schema='terluna.atmosphere.scale-size-diagnosis/1',
        evidence='Executed controlled numerical attribution within the original SCALE physics. No observed Earth population reproduction or validation of lunar size limits.',
        producer={p:lunar.base.digest(ROOT/p) for p in ('atmosphere/snow/scale_size_diagnosis.py','atmosphere/snow/scale_lunar.py','atmosphere/snow/scale_lunar_hook.f90','atmosphere/snow/scale_lunar_comparison.py')},
        expected_cases=CASES, references=REFERENCES, cases=results, comparisons=contrasts, missing_or_failed=missing,
        reading_rule='Same native binary and T,p,e path in all new forcing controls; the separately retained stronger-sticking reference uses its identified modified binary. Lunar/deep and Earth/shallow preserve gH. Crossed controls are imposed nonhydrostatic profiles with dynamics disabled. Water and seed flux are fixed except separately labeled sensitivities. Eight-column sampling errors exclude physical model error.')
    lunar.boxes.save(OUT,result)
    return result


def run(spec):
    d=lunar.CACHE/'runs'/spec['name']
    if not d.exists(): d=configure(**spec)
    saved=json.loads((d/'settings.json').read_text())
    for key, default in [('depth_factor',1.),('source_factor',1.),('water_factor',1.),('collisions',True)]:
        if saved['experiment'][key] != spec.get(key,default): raise ValueError('Cached experiment differs: '+key)
    for key, option, default in [('gravity_m_s2','gravity',MOON_SURFACE_GRAVITY),('dt_s','dt',60.),('packets_per_step_per_horizontal_cell','packets',4)]:
        if saved[key] != spec.get(option,default): raise ValueError('Cached experiment differs: '+key)
    if not (d/'analysis.json').exists():
        if not (d/'execution.json').exists(): lunar.execute(d)
        lunar.analyze(d)
    return spec['name']


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['run','assemble']);p.add_argument('--case',action='append')
    p.add_argument('--workers',type=int,choices=[1,2],default=1)
    args=p.parse_args()
    if args.action=='run':
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures=[pool.submit(run,c) for c in CASES if args.case is None or c['name'] in args.case]
            for future in as_completed(futures): print('complete',future.result(),flush=True)
    r=assemble()
    print('completed',len(r['cases']),'missing',len(r['missing_or_failed']),flush=True)
