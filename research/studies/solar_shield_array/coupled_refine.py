"""Continue the retained coupled model with stable timing and a shorter initializer."""
import json
import time
import signal
import numpy as np
from protection.dynamics.coupled_control import contact_set,contact_values,contact_jacobian,secant_update,trust_decision,solve_step
from .coupled_common import HERE, SPECIFIC, limit, inputs, write, save, read_raw, digest, length
from .coupled_cycle import CycleDebt, CycleExperiment
from .coupled_stable import StableExperiment


def truncate(raw,end):
    n=len(raw['t']);mask=raw['t']<=end
    return {k:(v[mask] if v.ndim and len(v)==n and k not in ['basis'] else v.copy()) for k,v in raw.items()}


def main():
    limit();signal.alarm(370);cpu=time.process_time();wall=time.monotonic()
    old=json.loads((HERE/'results/coupled_optimization_attempt.json').read_text())
    bench=json.loads((HERE/'results/coupled_benchmark.json').read_text())
    previous_charge=bench['resources']['cpu_s']+old['resources']['cpu_s']+90.+20.
    ceiling=1800-500-previous_charge
    ex=StableExperiment();original=CycleExperiment();horizon=14.*3600.;x=np.zeros(8)
    service=read_raw(bench['service']);full_base=read_raw(bench['evaluations'][0]);base=truncate(full_base,horizon)
    tail=read_raw(old['tail']);tail={k:(np.r_[full_base[k][full_base['t']>=horizon],v[1:]]
        if k in ['t','state','frames','power_W','attitude_acceleration','first_intercept_W','redirected_W','translation','distance_m'] else v)
        for k,v in tail.items()}
    debt=CycleDebt(ex,base,tail,service);terminal_limit=debt.summary(base,x)['normalized_rms']
    _,_,producer=inputs(['research/studies/solar_shield_array/coupled_cycle.py','research/studies/solar_shield_array/coupled_stable.py',
        'research/studies/solar_shield_array/coupled_refine.py','protection/dynamics/coupled_control.py'])
    producer['inputs'].update({n:digest(HERE/'results'/n) for n in ['coupled_benchmark.json','coupled_optimization_attempt.json']})
    out=dict(schema='terluna.research.coupled-optimization/1',producer=producer,completed=False,
        accepted_return=False,accepted_cycle=False,accepted_fleet=False,iterations=[],horizon_s=horizon,cycle_end_s=ex.end,
        previous_numerical_charge_s=previous_charge,cpu_ceiling_s=ceiling,validation_reserve_cpu_s=500.,
        interrupted_attempt_charge=dict(completed_cpu_s=old['resources']['cpu_s'],additional_partial_allowance_s=90.,
            reason='Interrupted during third trial after detecting near-duplicate attitude knots; charge a conservative 90 s for unsaved work.'),
        basis=old['basis'],response_comparison=old['response_comparison'],tail=old['tail'],
        initial_cycle_debt=debt.summary(base,x),derivatives=old['derivatives'],
        method_changes=['Coalesce attitude knots within 1 ms before spline construction.',
            'Shorten restoration horizon from 14.5 to 14 h; retain complete cycle terminal target.',
            'Rebuild new contact rows but judge new conflicts by physical clearance event progression; a band-limited signed-gap minimum is not an invariant merit between different contact sets.'])
    derivative=[];frame_errors=[];rate_errors=[]
    for j in range(8):
        p=np.zeros(8);p[j]=.2
        if j==7:raw={k:v.copy() for k,v in base.items()}
        else:raw=truncate(read_raw(old['derivatives'][j]),horizon)
        # Reuse only after checking the complete commanded frame/rate/acceleration
        # history at integration-scale dates for baseline and finite probes.
        dates=np.r_[base['t'],(base['t'][1:]+base['t'][:-1])/2]
        f=ex.command(p);g=original.command(p)
        frame_errors.append(float(abs(f(dates).as_matrix()-g(dates).as_matrix()).max()))
        rate_errors.append([float(abs(f(dates,k)-g(dates,k)).max()) for k in [1,2]])
        derivative.append(raw)
    if max(frame_errors)>1e-12 or np.max(rate_errors)>1e-12:raise ValueError('Cannot reuse old coupled derivatives')
    out['derivative_reuse_check']=dict(maximum_frame_element_difference=max(frame_errors),
        maximum_rate_or_acceleration_difference=float(np.max(rate_errors)),
        scope='No tolerance or force change; baseline and the eight finite derivative commands match. Failed near-grid trial commands are repaired and receive new replay.')
    state_jac=np.stack([(r['state']-base['state'])/.2 for r in derivative],axis=-1)
    frame_jac=np.stack([(r['frames']-base['frames'])/.2 for r in derivative],axis=-1)
    att_jac=np.stack([(r['attitude_acceleration']-base['attitude_acceleration'])/.2 for r in derivative],axis=-1)
    service_energy=bench['service']['electrical_energy_J'];m=tail['t']<=ex.arrival+21600
    future=float(np.trapezoid(tail['power_W'][m].sum(axis=1),tail['t'][m]))
    times=base['t'];pm=times<=43200;radius=.5;best=base
    base_energy=float(np.trapezoid(base['power_W'].sum(axis=1),times))
    best_row=dict(bench['evaluations'][0],end_s=horizon,intended_end_s=horizon,electrical_energy_J=base_energy)
    best_row.update(save('refined_base',**base))
    checkpoint=dict(x=x.tolist(),trust_radius=radius,horizon_s=horizon,target_arrival_s=ex.arrival,cycle_end_s=ex.end,
        terminal_rms_limit=terminal_limit,continuation='Restore with coupled_restart.restore(); relinearize coupled responses when extending the horizon. Full cycle remains target.')
    write('coupled_optimization.json',out,cpu,wall)
    for iteration in range(4):
        if time.process_time()-cpu+65>ceiling:out['stopped_reason']='Validation reserve';break
        cuts=contact_set(best['state'],best['frames']);gaps=contact_values(best['state'],best['frames'],cuts)
        gj=contact_jacobian(best['state'],best['frames'],state_jac,frame_jac,cuts)
        t0=debt.residual(best,x)[0];tj=[]
        for j in range(8):
            p={**best,'state':best['state']+.001*state_jac[:,:,:,j]};xx=x.copy();xx[j]+=.001
            tj.append((debt.residual(p,xx)[0]-t0)/.001)
        tj=np.stack(tj,axis=-1)
        def energy(d):
            att=np.maximum(0,best['attitude_acceleration']+np.einsum('tnc,c->tn',att_jac,d))
            power=(att+length(ex.thrust(x+d,times)))*ex.mass*SPECIFIC
            return service_energy+future+float(np.trapezoid(power.sum(axis=1),times)),service_energy+float(np.trapezoid(power[pm].sum(axis=1),times[pm]))
        def power_margin(d):
            att=np.maximum(0,best['attitude_acceleration']+np.einsum('tnc,c->tn',att_jac,d)).max(axis=0)
            u,w=ex.controls(x+d);a=length(u[:,0])*2/(w[0,1]-w[0,0])
            return 1-(att+a)*ex.mass*SPECIFIC/ex.capacity
        d,fit=solve_step(x,radius,gaps,gj,energy,power_margin,lambda d:t0+tj@d,terminal_limit,margin=250.)
        rec=dict(iteration=iteration,base_x=x.tolist(),radius=radius,contacts=len(cuts),fit=fit);out['iterations'].append(rec)
        if not fit['success']:
            radius/=2;write('coupled_optimization.json',out,cpu,wall);continue
        xx=x+d;row,raw,_=ex.evaluate(xx,service['state'][-1],21600.,horizon)
        row.update(save('refined_trial_'+str(iteration),**raw));rec['replay']=row
        actual_gap=contact_values(raw['state'],raw['frames'],cuts);actual_t=debt.residual(raw,xx)[0]
        def merit(g,e,r):return e/1e14+.2*np.mean(r**2)+100*max(0,(250-float(g.min()))/1000)
        before=merit(gaps,energy(np.zeros(8))[0],t0);pred=merit(gaps+gj@d,energy(d)[0],t0+tj@d)
        total=service_energy+future+row['electrical_energy_J'];actual=merit(actual_gap,total,actual_t)
        defect=raw['state']-best['state']-np.einsum('tnic,c->tni',state_jac,d)
        decision=trust_decision(before,pred,actual,radius,float(length(defect[:,:,:3]).max()))
        prefix=service_energy+float(np.trapezoid(raw['power_W'][pm].sum(axis=1),times[pm]))
        old_event=best_row['clearance_events_s'][0] if best_row.get('clearance_events_s') else horizon
        new_event=row['clearance_events_s'][0] if row['clearance_events_s'] else horizon
        hard=bool(row['overloaded_members']==0 and prefix<=51.859e12 and total<=1e14 and np.sqrt(np.mean(actual_t**2))<=terminal_limit*(1+1e-6))
        if not hard or new_event<old_event-1.:
            decision['accepted']=False;decision['next_radius']=radius/2
        decision.update(hard_limits_pass=hard,first_floor_event_change_s=new_event-old_event)
        rec.update(decision=decision,step=d.tolist(),before_merit=float(before),predicted_merit=float(pred),actual_merit=float(actual),
            actual_minimum_retained_signed_gap_m=float(actual_gap.min()),prefix_energy_J=prefix,
            estimated_through_next_service_J=total,cycle_debt=debt.summary(raw,xx))
        state_jac=secant_update(state_jac,d,raw['state']-best['state']);frame_jac=secant_update(frame_jac,d,raw['frames']-best['frames'])
        att_jac=secant_update(att_jac,d,raw['attitude_acceleration']-best['attitude_acceleration'])
        if decision['accepted']:x=xx;best=raw;best_row=row
        radius=decision['next_radius']
        checkpoint.update(x=x.tolist(),trust_radius=radius,last_iteration=iteration,best_raw_path=best_row['raw_path'],best_raw_sha256=best_row['raw_sha256'])
        checkpoint.update(save('refined_optimizer_state',state_jac=state_jac,frame_jac=frame_jac,attitude_jac=att_jac,
            cuts=cuts,gaps=gaps,gap_jac=gj,terminal_jac=tj,x=x,t=times))
        out.update(checkpoint=checkpoint,selected=best_row,selected_x=x.tolist(),final_cycle_debt=debt.summary(best,x))
        write('coupled_optimization.json',out,cpu,wall)
        print('refine',iteration,decision,'gap',actual_gap.min(),'cpu',time.process_time()-cpu,flush=True)
        if decision['accepted'] and actual_gap.min()>=249 and row['minimum_sampled_distance_m']>=249:
            out['stopped_reason']='Local initializer restored; event-stopped extension next';break
    out.update(completed=True,selected=best_row,selected_x=x.tolist(),final_cycle_debt=debt.summary(best,x))
    write('coupled_optimization.json',out,cpu,wall)


if __name__=='__main__':main()
