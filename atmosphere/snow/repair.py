"""Run and compare conservative sectional controls for the audited snow equations.

Cases are cached separately from the historical stochastic populations. The
mass-domain exit and independent refinements remain visible; no output is called
an empirically validated snowfall distribution. Use --case ID for one worker,
then --assemble to export all completed, identity-matching cases.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np
from atmosphere.snow import aggregation as ag, flakes as sf, sectional
from atmosphere.snow.habits import GEOMETRIES
from atmosphere.snow.population import histogram, distribution, arriving_distribution
from shared.provenance import constants_used

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / 'research/runs/snow_sectional'
OUT = ROOT / 'atmosphere/snow/results/repair.json'
SCHEMA = 'terluna.atmosphere.snow-sectional-population/1'
FILES = ('atmosphere/snow/repair.py','atmosphere/snow/sectional.py',
         'atmosphere/snow/population.py','atmosphere/snow/aggregation.py',
         'atmosphere/snow/habits.py','atmosphere/snow/flakes.py',
         'atmosphere/rain/drops.py','atmosphere/rain/targets.py',
         'atmosphere/radiative_convective/thermodynamics.py')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def identity(binding):
    inputs={p:sha(ROOT/p) for p in binding.get('inputs',{})}
    if inputs != binding.get('inputs',{}):
        raise RuntimeError('Input bytes no longer match the loaded column binding')
    return dict(files={p:sha(ROOT/p) for p in FILES}, inputs=inputs,
                constants=constants_used(FILES))


def cases():
    base = dict(world='moon', top_c=-15., crystals_per_litre=100., snowfall_mm_h=3.,
                geometry='bf95_dmax_dendrites_area', layers=30, mass_bins=96,
                minimum_diameter_m=2e-5, maximum_diameter_m=10., max_step_s=30.,
                cfl=.3, fill=1.5, gather=1.5, breakup=True)
    rows=[]
    def add(id, **changes): rows.append(dict(base, id=id, **changes))
    add('moon_heavy_m96')
    add('moon_moderate_m144',snowfall_mm_h=1.,mass_bins=144)
    add('moon_heavy_m144',mass_bins=144)
    add('moon_heavy_m216',mass_bins=216)
    add('moon_heavy_m144_dt15',mass_bins=144,max_step_s=15.)
    add('moon_heavy_m144_cfl015',mass_bins=144,cfl=.15)
    add('moon_heavy_m144_z60',mass_bins=144,layers=60)
    add('moon_heavy_m144_domain_extended',mass_bins=144,extra_upper_bins=26)
    add('earth_moderate_m144',world='earth',snowfall_mm_h=1.,mass_bins=144)
    return rows


def run_case(spec):
    world,binding=sf.moon_column() if spec['world']=='moon' else (sf.earth_column(),{})
    geometry=GEOMETRIES[spec['geometry']]
    before=identity(binding)
    args={key:spec[key] for key in ('layers','mass_bins','minimum_diameter_m','maximum_diameter_m',
                                  'max_step_s','cfl','fill','gather','breakup')}
    if spec.get('extra_upper_bins'):
        start=float(geometry.mass(args['minimum_diameter_m']))
        ratio=(float(geometry.mass(args['maximum_diameter_m']))/start)**(1/(args['mass_bins']-1))
        args['mass_grid']=start*ratio**np.arange(args['mass_bins']+spec['extra_upper_bins'])
    wall=time.monotonic(); next_report=0.
    def progress(t, steps):
        nonlocal next_report
        if t>=next_report:
            print(f"{spec['id']}: {t/3600:.1f} model h, {steps} steps, {time.monotonic()-wall:.1f} wall s",flush=True)
            next_report=t+6*3600.
    raw=sectional.column(world,spec['top_c'],spec['crystals_per_litre']*1000.,spec['snowfall_mm_h'],
                         geometry=geometry,progress=progress,**args)
    if identity(binding)!=before:
        raise RuntimeError('Producer or input identity changed during calculation; result not admitted')
    masses=raw['masses']; numbers=raw['arriving_numbers_m2']; diameters=geometry.size(masses)
    weights=masses*numbers
    cdf=np.cumsum(weights)/weights.sum() if weights.sum() else np.zeros_like(weights)
    def quantile(q):
        return float(diameters[min(np.searchsorted(cdf,q),diameters.size-1)]) if weights.sum() else None
    legacy_interface=dict(masses=masses,numbers=numbers,seconds=raw['seconds'])
    boundary=arriving_distribution(legacy_interface,world,geometry)
    beta=2*boundary['totals']['projected_area_m_inv']
    layer_products=[]
    for k,(z,n) in enumerate(zip(raw['layer_height_m'],raw['time_mean_number_m3'])):
        state=world.state(z)
        layer_products.append(dict(altitude_m=float(z),temperature_k=float(state['t']),pressure_pa=float(state['p']),
            **distribution(histogram(diameters,n,masses,geometry.area(diameters),raw['velocity_m_s'][k]))))
    arrival_rate=numbers*3600/raw['seconds']
    cumulative_rate=np.cumsum(arrival_rate[::-1])[::-1]
    eligible=diameters[cumulative_rate>=1.]
    summary=dict(mass_median_m=quantile(.5),mass_d90_m=quantile(.9),mass_d99_m=quantile(.99),
                 largest_hourly_m=float(eligible.max()) if eligible.size else None,
                 snowfall_mm_h=float(weights.sum()*3600/raw['seconds']),
                 flakes_per_m2_hour=float(arrival_rate.sum()),
                 mass_weighted_log_diameter=float(np.dot(weights,np.log(diameters))/weights.sum()) if weights.sum() else None,
                 extinction_m_inv=beta,beam_5_percent_distance_m=-math.log(.05)/beta if beta else None,
                 arriving_mass_share_above={str(d):float(weights[diameters>=d].sum()/weights.sum()) if weights.sum() else 0.
                                           for d in (.15,.6,1.,10.)})
    native={k:raw[k] for k in ('settings','mass_ledger_kg_m2','extensive_mass_ledger_kg','domain_exit_numbers_m2',
              'gather_domain_exit_mass_kg_m2','domain_adequacy','diagnostics','crossing_hours')}
    native['diagnostics']['wall_seconds']=time.monotonic()-wall
    return dict(id=spec['id'],specification=spec,**native,summary=summary,
                numerical_status='domain_passed_refinement_pending' if raw['domain_adequacy']['passed']
                                 else 'finite_domain_inadequate',
                arriving_mass_pivots_kg=masses.tolist(),arriving_diameter_pivots_m=diameters.tolist(),
                arriving_numbers_m2=numbers.tolist(),gathering_seconds=raw['seconds'],
                freezing_level=boundary,layers=layer_products,producer=before,
                evidence='Deterministic fixed-pivot solution of the inherited mean collision, chipping, deposition '
                         'and settling laws; no integer/fractional packet mixing or batched sticking. '
                         'Finite mass domain, mass-grid diffusion, time and spatial errors require separate controls. '
                         'Whole-aggregate fracture, evolving morphology and realistic moist trajectories remain absent.',
                reading_rule='Quantiles are discrete pivot quantiles and inherit mass-grid resolution. '
                             'Mass-domain exits are excluded from snowfall and explicitly accounted for. '
                             'Optical extinction uses Qext=2 and the exact discretized projected-area moment; '
                             'the 5% distance is homogeneous direct-beam attenuation, not visual range.')


def comparisons(rows):
    by_id={row['id']:row for row in rows}
    baseline=by_id.get('moon_heavy_m144')
    if baseline is None: return []
    result=[]
    for row in rows:
        if not row['id'].startswith('moon_heavy') or row is baseline: continue
        metrics={}
        for key in ('snowfall_mm_h','mass_median_m','mass_d90_m','mass_d99_m','extinction_m_inv','beam_5_percent_distance_m'):
            a,b=baseline['summary'][key],row['summary'][key]
            metrics[key]=dict(reference=a,value=b,relative_difference=(b/a-1) if a and b is not None else None)
        def mass_cdf(r,grid):
            d=np.asarray(r['arriving_diameter_pivots_m']); w=np.asarray(r['arriving_mass_pivots_kg'])*r['arriving_numbers_m2']
            cumulative=np.r_[0,np.cumsum(w)/w.sum()]
            return cumulative[np.searchsorted(d,grid,side='right')]
        grid=np.unique(np.r_[baseline['arriving_diameter_pivots_m'],row['arriving_diameter_pivots_m']])
        result.append(dict(reference=baseline['id'],case=row['id'],metrics=metrics,
            maximum_mass_cdf_difference=float(np.max(np.abs(mass_cdf(row,grid)-mass_cdf(baseline,grid))))))
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',choices=[s['id'] for s in cases()])
    parser.add_argument('--assemble',action='store_true')
    args=parser.parse_args()
    CACHE.mkdir(parents=True,exist_ok=True)
    if args.case:
        spec=next(s for s in cases() if s['id']==args.case)
        output=run_case(spec)
        (CACHE/(args.case+'.json')).write_text(json.dumps(output,indent=1,default=float)+'\n')
        print(json.dumps(dict(id=args.case,summary=output['summary'],domain=output['domain_adequacy'],ledger=output['mass_ledger_kg_m2'])),flush=True)
    if args.assemble:
        rows=[]; missing=[]
        for spec in cases():
            path=CACHE/(spec['id']+'.json')
            if not path.exists(): missing.append(spec['id']); continue
            row=json.loads(path.read_text())
            _,binding=sf.moon_column() if spec['world']=='moon' else (sf.earth_column(),{})
            if row['specification']!=spec or row['producer']!=identity(binding):
                raise RuntimeError(f'Stale case identity: {spec["id"]}')
            rows.append(row)
        OUT.write_text(json.dumps(dict(schema=SCHEMA,runs=rows,unperformed_cases=missing,
            comparisons=comparisons(rows),evidence='Numerical repair and independent controls; adequacy judged '
            'from recorded domain and refinement comparisons. No empirical lunar snow validation.'),indent=1)+'\n')
        print(f'Assembled {len(rows)} completed cases; {len(missing)} unperformed',flush=True)
    if not args.case and not args.assemble: parser.error('Choose --case or --assemble')


if __name__=='__main__': main()
