"""Six starting arrangements, two retained turns, and independent burn vectors."""
import json
import resource
import signal
import time

import numpy as np
from scipy.interpolate import CubicHermiteSpline

from shared.provenance import constants_used, constants_changed
from protection.dynamics.active_formation import solar_frame
from protection.dynamics.energy_coast import coast_command
from protection.dynamics.pattern_departure import expansion_maps, expansion_state
from protection.dynamics.packing_control import (packing_offsets, response_maps, matrices,
    controls_state, measured_base, separation_cuts, crossing_cuts, independent_impulses)
from protection.dynamics.pattern_return import encounter_scan, crossing_witnesses, beam_and_rate_screen
from protection.dynamics.natural_pattern import local_coverage
from protection.dynamics.optical import length
from .patterns_run import reference_shooting, affine_state
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE, SCENARIO

RUN=ROOT/'research/runs/solar_shield_array/packing'


def parents(names):
    values={n:json.loads((HERE/'results'/n).read_text()) for n in names};sources={}
    for p in values.values():
        for path,h in p['producer']['source_hashes'].items():
            if digest(ROOT/path)!=h:raise ValueError('Source changed: '+path)
        if constants_changed(p['producer']['constants']):raise ValueError('Constants changed')
        sources.update(p['producer']['source_hashes'])
    return values,sources


def load_raw(entry):
    if digest(ROOT/entry['raw_path'])!=entry['raw_sha256']:raise ValueError('Raw changed')
    return np.load(ROOT/entry['raw_path'])


def refine(env,base,maps,command,original,windows,trust=1.5):
    start,end=maps.t[[0,-1]];coarse=np.arange(start,end+1.,300.);fine=np.arange(start,end+1.,60.)
    controls=original.copy();cuts={};records=[];passed=False
    def predict(t,u):return controls_state(base,maps,t,u,original)
    for iteration in range(4):
        grid=coarse if iteration==0 else fine
        state=predict(grid,controls)
        roots=crossing_witnesses(grid,state,command,lambda t,ids:predict([t],controls)[0,ids],retain=96)
        additions=separation_cuts(grid,state,command,maps,controls)
        additions+=crossing_cuts(roots['witnesses'],command,maps,controls)
        for cut in additions:
            key=(round(cut['t'],4),cut['i'],cut['j'],cut['axis'],cut['sign'])
            cuts[key]=cut
        controls,fit=independent_impulses(controls,list(cuts.values()),windows,trust=trust)
        row=dict(iteration=iteration,optimization=fit);records.append(row)
        if not fit['success']:break
        state=predict(fine,controls)
        scan=encounter_scan(fine,state,command,clearance=350.,retain=8)
        roots=crossing_witnesses(fine,state,command,lambda t,ids:predict([t],controls)[0,ids],retain=32)
        beam=beam_and_rate_screen(env,fine,state,command)
        row.update(encounters=scan,crossings=roots,beam_and_rate=beam)
        passed=bool(not scan['sampled_pair_violations'] and not roots['crossing_near_encounters']
            and beam['minimum_sampled_beam_margin_deg']>0 and beam['maximum_sampled_rate_deg_s']<=.1)
        if passed:break
    return controls,predict(fine,controls),fine,dict(proposal_pass=passed,iterations=records,
        mean_delta_v_m_s=float(length(controls).sum(axis=1).mean()),
        maximum_tile_delta_v_m_s=float(length(controls).sum(axis=1).max()),
        burn_count=int(np.count_nonzero(length(controls)>1e-8)))


def main():
    signal.alarm(1200);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();RUN.mkdir(parents=True,exist_ok=True)
    names=['patterns.json','energy.json','energy_replay.json','energy_screen.json']
    p,sources=parents(names)
    for path in ['protection/dynamics/packing_control.py','research/studies/solar_shield_array/packing_screen.py']:
        sources[path]=digest(ROOT/path)
    env=environment(1.);reference,gradient=reference_shooting(env,21600.)
    service=load_raw(p['patterns.json']['replay']);initial=service['state'][0];terminal=service['state'][-1]
    frame=solar_frame(env.at(0.));basis=solar_frame(env.at(21600.))
    windows=np.array([[21600.,28800.],[32400.,39600.]])
    config_rows=[];cases=[];out=dict(schema='terluna.research.packing-screen/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),
            inputs={n:digest(HERE/'results'/n) for n in names}),
        epoch_tdb=SCENARIO['epoch_tdb'],start_s=21600.,end_s=43200.,windows_s=windows.tolist(),
        count=361,objective='Sparse L1 individual-vector proposal; report Euclidean impulse. Coupled service/departure and actual torque required.',
        force_scope='Measured earlier departure defects plus common gravity response. Changed packing alters shadow force and requires nonlinear replay.',
        accepted_fleet=False,accepted_departure=False,recurring_cost_measured=False,
        collection_capacity_W=None,delivered_power_W=None,configs=config_rows,cases=cases)
    def save():
        out.update(elapsed_s=time.monotonic()-before,max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        (HERE/'results/packing_screen.json').write_text(json.dumps(out,indent=2)+'\n')
    baseline=reference.sol(0.)[:6]+np.c_[packing_offsets(4000.)@frame.T,packing_offsets(4000.)@(gradient*.95).T]
    out['baseline_initial_reconstruction_error_m']=float(length(baseline[:,:3]-initial[:,:3]).max())
    if np.max(abs(baseline-initial))>1e-5:raise ValueError('Baseline initial reconstruction changed')
    entries=[(24,p['energy.json']['runs'][0]),(41,p['energy_replay.json']['runs'][0])]
    templates={}
    for number,entry in entries:
        raw=load_raw(entry);spec=p['energy_screen.json']['cases'][number]['schedule']
        command=coast_command(env,21600.,43200.,**spec)
        oldmaps=expansion_maps(env,terminal,21600.,43200.,7200.)
        maps=response_maps(env,lambda t,m=oldmaps:m.sol(t)[:6],21600.,43200.,windows,basis)
        original=np.zeros((361,2,3));original[:,0]=raw['impulses']@basis
        previous=lambda t,m=oldmaps,u=raw['impulses']:expansion_state(m,terminal,u,t)
        base,defect=measured_base(previous,maps,raw['t'],raw['state'])
        templates[number]=(command,maps,base,original,spec,defect)
    for depth in [4000.,8000.,12000.]:
        for swapped in [False,True]:
            offsets=packing_offsets(depth,swapped);new_initial=reference.sol(0.)[:6]+np.c_[offsets@frame.T,offsets@(gradient*.95).T]
            service_t=np.linspace(0.,21600.,13)
            affine=np.array([affine_state(reference,t,offsets,gradient*.95,frame) for t in service_t])
            coverage=[local_coverage(y,env.at(t),reference.sol(t)[:6]) for t,y in zip(service_t,affine)]
            config=dict(index=len(config_rows),depth_m=depth,swapped=swapped,pitch_m=8500.,rows=19,velocity_gradient_scale=.95,
                minimum_proposal_coverage=min(c['minimum_source_coverage'] for c in coverage),coverage=coverage,
                initial_relative_speed_mean_m_s=float(length(new_initial[:,3:]-reference.sol(0.)[3:6]).mean()))
            config_rows.append(config)
            delta=(new_initial-initial)@reference.sol(21600.)[6:].reshape(6,6).T
            for number,(command,maps,base,original,spec,defect) in templates.items():
                changed=lambda t,b=base,m=maps,d=delta:b(t)+np.einsum('tij,nj->tni',matrices(m,t)[:,:,:6],d)
                if config['minimum_proposal_coverage']<1-1e-9:
                    cases.append(dict(config=config['index'],schedule_case=number,proposal_pass=False,rejected='service_coverage'));save();continue
                controls,state,times,result=refine(env,changed,maps,command,original,windows)
                result.update(config=config['index'],schedule_case=number,schedule=spec,
                    nominal_attitude_m_s=p['energy_screen.json']['cases'][number]['nominal_attitude_m_s'],
                    earlier_measured_force_defect_m=defect)
                result['nominal_total_m_s']=result['mean_delta_v_m_s']+result['nominal_attitude_m_s']
                path=RUN/f"proposal_{config['index']}_{number}.npz"
                np.savez_compressed(path,t=times,state=state,controls=controls,basis=basis,initial=new_initial,windows=windows)
                result.update(raw_path=str(path.relative_to(ROOT)),raw_sha256=digest(path));cases.append(result);save()
                print(json.dumps(dict(stage='packing_proposal',config=config['index'],schedule=number,
                    passed=result['proposal_pass'],nominal_m_s=result['nominal_total_m_s'],elapsed_s=time.monotonic()-before)),flush=True)
    out['ranked_proposals']=sorted([i for i,c in enumerate(cases) if c['proposal_pass']],key=lambda i:cases[i]['nominal_total_m_s'])
    out['completed']=True
    if any(digest(ROOT/path)!=h for path,h in sources.items()):raise ValueError('Source changed')
    save()


if __name__=='__main__':main()
