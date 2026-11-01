"""Separate seed-shape and growth-package sensitivities for 2018 humidity cases.

These are declared, unfitted seed shapes, since the observed initial shapes were
not reported per case. Fixed maximum seed span means mass changes with thickness.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import re
from atmosphere.snow import growth_humidity_validation as humidity
from atmosphere.snow import mcsnow_ventilation_controls as controls

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CACHE=ROOT/'research/runs/snow_scale/growth_humidity_history'
OUT=HERE/'results/growth_humidity_history.json'
SHAPES=(.1,.3,1.)
MODES=('native','both')
CASE_IDS=('C','D','G')
DT=.0125


def producer():
    return {'module_sha256':humidity.sha(__file__),'humidity_harness':humidity.producer()}


def configuration(row,mode,growth_mode,phi,ice_density,dt=DT):
    if phi not in SHAPES:raise ValueError('Undeclared seed shape')
    if growth_mode not in (1,2):raise ValueError('Unknown growth package')
    config=humidity.configuration(row,mode,10.,dt,ice_density)
    config.update(growth_mode=growth_mode,initial_polar_m=config['initial_equatorial_m']*phi)
    config['initial_mass_kg']=math.pi/6*config['initial_equatorial_m']**2*config['initial_polar_m']*ice_density
    return config


def branch_diagnostic(path,result,growth_mode):
    """Distinguish exact IGF2 threshold events from sparsely observed native traces."""
    traces=result['growth_trace'];first=next((i for i,r in enumerate(traces) if r['branch_criterion_met']),None)
    bracket=None
    if first is not None:
        bracket=[0. if first==0 else traces[first-1]['evaluation_elapsed_s'],traces[first]['evaluation_elapsed_s']]
    record=dict(first_saved_active_branch_source_time_s=traces[first]['evaluation_elapsed_s'] if first is not None else None,
                first_source_branch_crossing_time_bracket_s=bracket,
                source_sampling='first step and one-minute intervals; observed first-positive time is not exact crossing',
                any_saved_branch_criterion_met=first is not None)
    if growth_mode==2:
        rows,errors=controls.native.read_csv(path/'branch_diagnostic.csv',('elapsed_s','branch_threshold_ratio','branch_criterion_met'))
        if errors:raise ValueError('Malformed branching event file')
        record['native_maximum_radius_200um_events']=rows
        record['first_maximum_radius_200um_event_time_s']=next((r['elapsed_s'] for r in rows if r['branch_criterion_met']),None)
        record['event_semantics']='Post-step maximum radius crosses 200 micrometres; native pre-step branch test acts on following step.'
    return record


def assemble(do_run=False,build_dir=humidity.BUILD,cache=CACHE,output=OUT):
    build_dir=Path(build_dir);cache=Path(cache);output=Path(output)
    data=humidity.observations();byid={r['id']:r for r in data['cases']};own=producer();build=controls.verify_build(build_dir)
    if do_run:cache.mkdir(parents=True,exist_ok=True)
    match=re.search(r'\brhoi\s*=\s*([0-9.]+)_wp',(build_dir/'source/src/mo_atmo_types.f90').read_text(),re.I)
    if not match:raise ValueError('Native ice density unavailable')
    density=float(match[1]);runs=[]
    for id_ in CASE_IDS:
        row=byid[id_]
        for growth_mode in (1,2):
            for mode in MODES:
                for phi,dt in [(phi,DT) for phi in SHAPES]+[(phi,DT/4) for phi in SHAPES if (id_=='G' and phi in (.1,.3)) or (id_=='D' and phi==.1)]:
                    config=configuration(row,mode,growth_mode,phi,density,dt)
                    metadata={'experiment':'Takahashi 2018 seed-shape/growth-package sensitivity','datum_id':id_,
                        'source_pdf_sha256':humidity.SOURCE_SHA,'observations_sha256':humidity.sha(humidity.DATA),
                        'initial_aspect_ratio':phi,'initial_full_equatorial_span_m':config['initial_equatorial_m'],
                        'initial_full_polar_span_m':config['initial_polar_m'],'initial_mass_kg':config['initial_mass_kg'],
                        'initial_density_kg_m3':density,'growth_mode':growth_mode,'mode':mode,
                        'initial_state_status':'Unreported observed initial state; declared spheroid sensitivity, not fitted.',
                        'mass_convention':'Fixed 10-micrometre maximum span with full spheroid volume and bulk density; mass varies with thickness.',
                        'nominal_temperature_C':row['nominal_temperature_C'],'vapor_density_excess_g_m3':row['vapor_density_excess_g_m3']}
                    name=f'{id_}_igf{growth_mode}_{mode}_phi{phi:g}'+(f'_dt{dt:g}' if dt!=DT else '')
                    path=cache/name.replace('.','p')
                    if not path.exists():
                        if not do_run:raise FileNotFoundError(path)
                        controls.run_configuration(build_dir,path,config,metadata=metadata)
                    result=humidity.verified_record(path,config,metadata)
                    if result['build_record_sha256']!=humidity.sha(build_dir/'build.json'):raise ValueError('History sensitivity build differs')
                    record=humidity.compact(path,result,row,mode,10.,dt,0.)
                    record.update(growth_mode=growth_mode,initial_aspect_ratio=phi,initial_state=metadata,
                                  branching=branch_diagnostic(path,result,growth_mode),
                                  source_growth_factor_history=[{k:r[k] for k in ('evaluation_elapsed_s','aspect_ratio','inherent_growth_factor','effective_growth_factor','deposition_density_kg_m3','branch_criterion_met','ventilation_vapor')} for r in result['growth_trace']])
                    runs.append(record)
        if do_run:print('completed seed-shape case',id_,flush=True)
    if own!=producer() or build!=controls.verify_build(build_dir):raise RuntimeError('History sensitivity identity changed')
    refinements=[]
    for fine in (r for r in runs if r['dt_s']==DT/4):
        coarse=next(r for r in runs if r['id']==fine['id'] and r['mode']==fine['mode'] and r['growth_mode']==fine['growth_mode'] and r['initial_aspect_ratio']==fine['initial_aspect_ratio'] and r['dt_s']==DT)
        comparison=humidity.compare_refinement(coarse,fine)
        early_a=[r for r in coarse['trajectory'] if r['elapsed_s']<=300]
        early_b=[r for r in fine['trajectory'] if r['elapsed_s']<=300]
        comparison['first_five_minutes_diameter_max_relative_change']=max(abs(a['diameter_m']/b['diameter_m']-1) for a,b in zip(early_a,early_b))
        refinements.append(dict(id=fine['id'],mode=fine['mode'],growth_mode=fine['growth_mode'],initial_aspect_ratio=fine['initial_aspect_ratio'],**comparison))
    result=dict(schema='terluna.atmosphere.growth-humidity-history/1',producer=own,
        source_observations_file=str(humidity.DATA.relative_to(ROOT)),source_observations_sha256=humidity.sha(humidity.DATA),
        build_record_sha256=humidity.sha(build_dir/'build.json'),counts={'runs':len(runs),'complete':sum(r['complete'] for r in runs),'base_runs':36,'refinement_pairs':len(refinements)},refinements=refinements,
        evidence='Unfitted initial-shape and original/revised growth-package sensitivities against additional 2018 evidence; no model replacement adopted.',
        limits=['Observed initial geometry is unknown; these three aspect ratios are not observed per-case inputs or confidence bounds.',
            'Initial maximum span is fixed at10micrometres; dense spheroid mass scales with initial aspect ratio.',
            'IGF1 and IGF2 change intrinsic growth tables and branching rules together; differences cannot isolate the threshold alone.',
            'Temperatures remain nominal groups; vapor density excess and site mean pressure are prescribed.',
            'Spheroid growth does not predict secondary habit; observed mass and thickness are unavailable.',
            'Base runs use dt0.0125s; twelve selected thin-seed runs are additionally refined to0.003125s. Late failed trajectories remain failure diagnostics.'],runs=runs)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');return result


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run',action='store_true');a=parser.parse_args()
    print(json.dumps(assemble(a.run)['counts'],indent=2))
if __name__=='__main__':main()
