"""Compact physical-prefix comparisons and final budget, without new trajectories."""
import json
import time
import numpy as np
from protection.dynamics.packing_control import response_maps,matrices
from .coupled_common import HERE,ROOT,RUN,inputs,write,read_raw,spline,state_at,digest,length
from .coupled_stable import StableExperiment


def main():
    cpu=time.process_time();wall=time.monotonic();ex=StableExperiment()
    names=['benchmark','optimization_attempt','refinement_attempt','optimization','validation_attempt','validation','cycle_probe','restart']
    products={n:json.loads((HERE/'results'/('coupled_'+n+'.json')).read_text()) for n in names}
    _,_,producer=inputs(['research/studies/solar_shield_array/coupled_cycle.py','research/studies/solar_shield_array/coupled_stable.py',
        'research/studies/solar_shield_array/coupled_audit.py','protection/dynamics/packing_control.py'])
    producer['inputs'].update({'coupled_'+n+'.json':digest(HERE/'results'/('coupled_'+n+'.json')) for n in names})
    b=products['benchmark'];base=read_raw(b['evaluations'][0]);full=read_raw(b['evaluations'][2]);half=read_raw(b['evaluations'][3])
    path=spline(base);maps=response_maps(ex.env,lambda t:state_at(path,t).mean(axis=0),21600.,base['t'][-1],
        np.array([[21600.,43200.]]),np.eye(3))
    gravitational=np.einsum('tij,nj->tni',matrices(maps,base['t'],1)[:,:,6:],ex.basis[:,:,0])
    response=(full['state']-base['state'])/.2;observed=half['state']-base['state']
    ce=length((observed-.1*response)[:,:,:3]);ge=length((observed-.1*gravitational)[:,:,:3]);pm=base['t']<=43200
    event=min(r['clearance_events_s'][0] for r in b['evaluations'] if r['clearance_events_s']);common=base['t']<event
    service=read_raw(b['service']);v=products['validation'];fine=read_raw(v['runs'][-1]);cmd=ex.command(np.array(v['selected_x']))
    t=v['terminal']['time_s'];i,j=v['terminal']['clearance_witness']['pair'];frame=cmd(t).as_matrix()
    by=state_at(path,t);fy=fine['state'][-1]
    coarse=v['coarse'];coarse_energy=b['service']['electrical_energy_J']+coarse['electrical_energy_J']
    coarse_light=b['service']['first_intercept_J']+coarse['first_intercept_J']
    out=dict(schema='terluna.research.coupled-audit/1',producer=producer,completed=True,
        accepted_return=False,accepted_cycle=False,accepted_fleet=False,
        baseline_twelve_hour_energy_J=b['service']['electrical_energy_J']+float(np.trapezoid(base['power_W'][pm].sum(axis=1),base['t'][pm])),
        physical_prefix_response=dict(end_s=43200.,coupled_maximum_position_error_m=float(ce[pm].max()),
            gravity_maximum_position_error_m=float(ge[pm].max()),
            coupled_rms_position_error_m=float(np.sqrt(np.mean(ce[pm]**2))),gravity_rms_position_error_m=float(np.sqrt(np.mean(ge[pm]**2)))),
        common_pre_event_response=dict(last_sample_s=float(base['t'][common][-1]),first_event_s=event,
            coupled_maximum_position_error_m=float(ce[common].max()),gravity_maximum_position_error_m=float(ge[common].max())),
        final_encounter=dict(time_s=t,pair=[i,j],baseline_body_displacement_m=((by[i,:3]-by[j,:3])@frame).tolist(),
            changed_body_displacement_m=((fy[i,:3]-fy[j,:3])@frame).tolist()),
        accounting_comparison=dict(relative_control_energy_difference=abs(coarse_energy-v['account']['electrical_energy_J'])/v['account']['electrical_energy_J'],
            relative_interception_difference=abs(coarse_light-v['account']['first_intercept_J'])/v['account']['first_intercept_J'],
            scope='Joint changes in time/source quadrature and the 1.86 s stopping-time difference; not separately isolated convergence orders.'),
        scientific_conclusion='Coupled local responses improve prediction; eight-mode restoration does not produce a safe repair. Full-cycle boundary debt, coverage, arrival controls and supply remain unresolved.')
    write('coupled_audit.json',out,cpu,wall)
    print(out['physical_prefix_response'],out['accounting_comparison'],out['final_encounter'],flush=True)


if __name__=='__main__':main()
