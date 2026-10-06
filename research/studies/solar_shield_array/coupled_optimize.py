"""Bounded coupled finite-difference SQP with replay-tested trust regions."""
import signal
import time
import json
import numpy as np

from protection.dynamics.coupled_control import (
    contact_set, contact_values, contact_jacobian, secant_update, trust_decision, solve_step)
from protection.dynamics.packing_control import response_maps, matrices
from .coupled_common import (HERE, ROOT, RUN, DIMENSION, SPECIFIC, limit, inputs,
    write, save, read_raw, spline, state_at, digest, length)
from .coupled_cycle import CycleExperiment, CycleDebt


def main():
    limit(); signal.alarm(1200)
    cpu=time.process_time();wall=time.monotonic()
    bench=json.loads((HERE/'results/coupled_benchmark.json').read_text())
    used=bench['resources']['cpu_s']+20.
    ceiling=1800.-500.-used
    ex=CycleExperiment()
    _,_,producer=inputs(['protection/dynamics/coupled_control.py',
        'research/studies/solar_shield_array/coupled_cycle.py',
        'research/studies/solar_shield_array/coupled_optimize.py'])
    producer['inputs']['coupled_benchmark.json']=digest(HERE/'results/coupled_benchmark.json')
    service=read_raw(bench['service']);base=read_raw(bench['evaluations'][0])
    initial=service['state'][-1];horizon=float(base['t'][-1]);x=np.zeros(DIMENSION)
    out=dict(schema='terluna.research.coupled-optimization/1',producer=producer,
        completed=False,accepted_return=False,accepted_cycle=False,accepted_fleet=False,
        cpu_ceiling_s=ceiling,validation_reserve_cpu_s=500.,derivatives=[],iterations=[],
        horizon_s=horizon, cycle_end_s=ex.end,
        basis=dict(names=['row_parity','column_parity','parity_product','row_gradient',
            'pair_325_347','pair_338_341','turn_duration','burn_end'],
            impulse_scale_m_s=.1,timing_scale_hours=.5,
            controls_repeated_at_next_departure=True,
            fixed=['initial service','tile attitudes within each date','arrival date','identity assignment','arrival packing target'],
            excludes='Arbitrary tile histories, independent normals, arrival arcs and reassignments; terminal tail transport is an explicitly approximate initializer.'))
    # Baseline tail continues through contacts solely for cycle residuals.
    row,tail,_=ex.evaluate(x,base['state'][-1],horizon,ex.end,step=240.)
    row.update(save('baseline_cycle_tail',**tail));row['restoration_only']=True;out['tail']=row
    print('cycle tail',row['cpu_s'],flush=True)
    debt=CycleDebt(ex,base,tail,service);initial_debt=debt.summary(base,x)
    out['initial_cycle_debt']=initial_debt;terminal_limit=initial_debt['normalized_rms']
    write('coupled_optimization.json',out,cpu,wall)
    # Exact repeated command equals the benchmark throughout its local horizon.
    frame_error=float(abs(ex.command(x)(base['t']).as_matrix()-base['frames']).max())
    out['benchmark_command_maximum_matrix_difference']=frame_error
    if frame_error>1e-12:raise ValueError('Changed preserved prefix command')
    step=.2;perturb=[]
    for j in range(DIMENSION):
        if time.process_time()-cpu+55 > ceiling:raise RuntimeError('Derivative budget exhausted; retain checkpoint')
        if j==0:
            raw=read_raw(bench['evaluations'][2]);row=dict(bench['evaluations'][2]);row['reused_benchmark']=True
        elif j==7:
            # At zero impulse burn timing has identically zero sensitivity.
            raw={k:np.array(v,copy=True) for k,v in base.items()};raw['x'][j]=step
            row=dict(name='zero_control_burn_timing',cpu_s=0.,analytic_zero=True)
        else:
            p=x.copy();p[j]+=step
            row,raw,_=ex.evaluate(p,initial,21600.,horizon)
            row.update(save('derivative_'+str(j),**raw))
        perturb.append(raw);row['parameter']=j;out['derivatives'].append(row)
        write('coupled_optimization.json',out,cpu,wall)
        print('derivative',j,row['cpu_s'],flush=True)
    state_jac=np.stack([(p['state']-base['state'])/step for p in perturb],axis=-1)
    frame_jac=np.stack([(p['frames']-base['frames'])/step for p in perturb],axis=-1)
    att_jac=np.stack([(p['attitude_acceleration']-base['attitude_acceleration'])/step for p in perturb],axis=-1)
    # Like-for-like defect-corrected gravity map, sharing the coupled intercept.
    windows=ex.controls(x)[1][:1]
    maps=response_maps(ex.env,lambda t:state_at(spline(base),t).mean(axis=0),21600.,horizon,windows,np.eye(3))
    gravity=np.einsum('tij,nj->tni',matrices(maps,base['t'],1)[:,:,6:],ex.basis[:,:,0])
    half=read_raw(bench['evaluations'][3]);observed=half['state']-base['state']
    coupled_error=length((observed-.1*state_jac[:,:,:,0])[:,:,:3])
    gravity_error=length((observed-.1*gravity)[:,:,:3])
    out['response_comparison']=dict(physical_held_out_impulse_m_s=.01,
        coupled_maximum_position_error_m=float(coupled_error.max()),
        gravity_only_maximum_position_error_m=float(gravity_error.max()),
        coupled_rms_position_error_m=float(np.sqrt(np.mean(coupled_error**2))),
        gravity_only_rms_position_error_m=float(np.sqrt(np.mean(gravity_error**2))),
        scope='Both start on the same coupled baseline; only the control response differs. All evaluated members retained.')
    # Carry future attitude demand explicitly as an unexecuted restoration cost.
    mask=tail['t']<=ex.arrival+21600
    future_energy=float(np.trapezoid(tail['power_W'][mask].sum(axis=1),tail['t'][mask]))
    service_energy=bench['service']['electrical_energy_J']
    radius=1.;best=base;best_row=bench['evaluations'][0]
    checkpoint=dict(x=x.tolist(),trust_radius=radius,horizon_s=horizon,
        target_arrival_s=ex.arrival,cycle_end_s=ex.end,terminal_rms_limit=terminal_limit,
        continuation='Resume coupled local Jacobian/merit search; extend accepted prefix to actual next service and repeated departure. Tail proxy is not acceptance.')
    for iteration in range(6):
        if time.process_time()-cpu+65>ceiling:
            out['stopped_reason']='Reserved independent-validation budget';break
        cuts=contact_set(best['state'],best['frames'])
        gaps=contact_values(best['state'],best['frames'],cuts)
        gap_jac=contact_jacobian(best['state'],best['frames'],state_jac,frame_jac,cuts)
        terminal0,_=debt.residual(best,x);terminal_jac=[]
        for j in range(DIMENSION):
            p={**best,'state':best['state']+state_jac[:,:,:,j]*.001}
            xx=x.copy();xx[j]+=.001
            terminal_jac.append((debt.residual(p,xx)[0]-terminal0)/.001)
        terminal_jac=np.stack(terminal_jac,axis=-1)
        times=best['t'];prefix_mask=times<=43200
        def energy(d):
            attitude=np.maximum(0,best['attitude_acceleration']+np.einsum('tnc,c->tn',att_jac,d))
            power=(attitude+length(ex.thrust(x+d,times)))*ex.mass*SPECIFIC
            return (service_energy+future_energy+float(np.trapezoid(power.sum(axis=1),times)),
                service_energy+float(np.trapezoid(power[prefix_mask].sum(axis=1),times[prefix_mask])))
        def power_margin(d):
            attitude=np.maximum(0,best['attitude_acceleration']+np.einsum('tnc,c->tn',att_jac,d)).max(axis=0)
            u,win=ex.controls(x+d)
            maximum=length(u[:,0])*2/(win[0,1]-win[0,0])
            # Triangle bound combines each member's separate attitude/translation peaks.
            return 1-(attitude+maximum)*ex.mass*SPECIFIC/ex.capacity
        terminal=lambda d:terminal0+terminal_jac@d
        d,fit=solve_step(x,radius,gaps,gap_jac,energy,power_margin,terminal,terminal_limit)
        record=dict(iteration=iteration,base_x=x.tolist(),radius=radius,contacts=len(cuts),fit=fit)
        out['iterations'].append(record)
        if not fit['success']:
            radius/=2;record['rejected_reason']='SQP subproblem infeasible or unconverged'
            write('coupled_optimization.json',out,cpu,wall);continue
        proposal=x+d
        def merit(g,e,r):return e/1e14+.2*np.mean(r**2)+100*max(0.,(400-float(g.min()))/1000)
        before=merit(gaps,energy(np.zeros(DIMENSION))[0],terminal0)
        predicted=merit(gaps+gap_jac@d,energy(d)[0],terminal(d))
        row,raw,_=ex.evaluate(proposal,initial,21600.,horizon)
        row.update(save('trial_'+str(iteration),**raw));record['replay']=row
        actual_gaps=contact_values(raw['state'],raw['frames'],cuts)
        actual_terminal,_=debt.residual(raw,proposal)
        actual_energy=service_energy+future_energy+row['electrical_energy_J']
        actual=merit(actual_gaps,actual_energy,actual_terminal)
        defect=raw['state']-best['state']-np.einsum('tnic,c->tni',state_jac,d)
        decision=trust_decision(before,predicted,actual,radius,float(length(defect[:,:,:3]).max()))
        prefix_energy=service_energy+float(np.trapezoid(raw['power_W'][prefix_mask].sum(axis=1),times[prefix_mask]))
        decision['hard_limits_pass']=bool(row['overloaded_members']==0 and prefix_energy<=51.859e12 and actual_energy<=1e14
            and np.sqrt(np.mean(actual_terminal**2))<=terminal_limit*(1+1e-6))
        newcuts=contact_set(raw['state'],raw['frames'])
        newgaps=contact_values(raw['state'],raw['frames'],newcuts)
        # Compare new branches as well, so improving retained cuts cannot hide a new worse contact.
        new_merit=merit(newgaps,actual_energy,actual_terminal)
        decision['new_contact_merit']=float(new_merit)
        if not decision['hard_limits_pass'] or new_merit>=before:
            decision['accepted']=False;decision['next_radius']=radius/2
        record.update(decision=decision,step=d.tolist(),before_merit=float(before),predicted_merit=float(predicted),
            actual_merit=float(actual),predicted_minimum_signed_gap_m=float((gaps+gap_jac@d).min()),
            actual_minimum_signed_gap_m=float(newgaps.min()),prefix_energy_J=prefix_energy,
            estimated_through_next_service_J=actual_energy,cycle_debt=debt.summary(raw,proposal))
        # Every propagated direction supplies a secant, even a rejected trial.
        state_jac=secant_update(state_jac,d,raw['state']-best['state'])
        frame_jac=secant_update(frame_jac,d,raw['frames']-best['frames'])
        att_jac=secant_update(att_jac,d,raw['attitude_acceleration']-best['attitude_acceleration'])
        if decision['accepted']:
            x=proposal;best=raw;best_row=row
        radius=decision['next_radius']
        checkpoint.update(x=x.tolist(),trust_radius=radius,last_iteration=iteration,
                          best_raw_path=best_row['raw_path'],best_raw_sha256=best_row['raw_sha256'])
        checkpoint.update(save('optimizer_state',state_jac=state_jac,frame_jac=frame_jac,attitude_jac=att_jac,
            cuts=cuts,gaps=gaps,gap_jac=gap_jac,terminal_jac=terminal_jac,x=x,t=times))
        out['checkpoint']=checkpoint;out['selected']=best_row
        write('coupled_optimization.json',out,cpu,wall)
        print('iteration',iteration,decision,'gap',float(newgaps.min()),'cpu',time.process_time()-cpu,flush=True)
        if decision['accepted'] and newgaps.min()>=399. and row['minimum_sampled_distance_m']>=399.:
            out['stopped_reason']='Local restoration achieved; extend and independently validate';break
        if radius<.0625:
            out['stopped_reason']='Trust region too small for remaining local model';break
    out['selected']=best_row;out['selected_x']=x.tolist();out['final_cycle_debt']=debt.summary(best,x)
    out['completed']=True;write('coupled_optimization.json',out,cpu,wall)


if __name__=='__main__':main()
