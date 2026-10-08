"""One bounded rebuild retaining all newly exposed contact branches."""
import json
import signal
import time
import numpy as np
from protection.dynamics.coupled_control import contact_set,contact_values,contact_jacobian,secant_update,trust_decision,solve_step
from .coupled_common import HERE,SPECIFIC,limit,inputs,write,save,read_raw,digest,length
from .coupled_cycle import CycleDebt
from .coupled_stable import StableExperiment


def main():
    limit();cpu=time.process_time();wall=time.monotonic()
    old=json.loads((HERE/'results/coupled_refinement_attempt.json').read_text())
    bench=json.loads((HERE/'results/coupled_benchmark.json').read_text())
    charged=old['previous_numerical_charge_s']+old['resources']['cpu_s'];allowance=1800-500-charged
    signal.alarm(int(allowance))
    ex=StableExperiment();base=read_raw(old['selected']);service=read_raw(bench['service']);model=read_raw(old['checkpoint'])
    fullbase=read_raw(bench['evaluations'][0]);tail=read_raw(old['tail']);horizon=old['horizon_s']
    tail={k:(np.r_[fullbase[k][fullbase['t']>=horizon],v[1:]] if k in
        ['t','state','frames','power_W','attitude_acceleration','first_intercept_W','redirected_W','translation','distance_m'] else v) for k,v in tail.items()}
    debt=CycleDebt(ex,base,tail,service);x=np.array(old['selected_x']);radius=.5
    sets=[contact_set(base['state'],base['frames'])]
    for r in old['iterations']:
        if 'replay' in r:
            raw=read_raw(r['replay']);sets.append(contact_set(raw['state'],raw['frames']))
    cuts=np.unique(np.r_[*sets],axis=0)
    # Every retained branch is evaluated/linearized at the actual current base.
    gaps=contact_values(base['state'],base['frames'],cuts)
    sj=model['state_jac'];fj=model['frame_jac'];aj=model['attitude_jac']
    gj=contact_jacobian(base['state'],base['frames'],sj,fj,cuts)
    t0=debt.residual(base,x)[0];tj=[]
    for j in range(8):
        p={**base,'state':base['state']+.001*sj[:,:,:,j]};xx=x.copy();xx[j]+=.001
        tj.append((debt.residual(p,xx)[0]-t0)/.001)
    tj=np.stack(tj,axis=-1)
    times=base['t'];pm=times<=43200.;m=tail['t']<=ex.arrival+21600
    future=float(np.trapezoid(tail['power_W'][m].sum(axis=1),tail['t'][m]));service_energy=bench['service']['electrical_energy_J']
    def energy(d):
        att=np.maximum(0,base['attitude_acceleration']+np.einsum('tnc,c->tn',aj,d))
        p=(att+length(ex.thrust(x+d,times)))*ex.mass*SPECIFIC
        return service_energy+future+float(np.trapezoid(p.sum(axis=1),times)),service_energy+float(np.trapezoid(p[pm].sum(axis=1),times[pm]))
    def power_margin(d):
        att=np.maximum(0,base['attitude_acceleration']+np.einsum('tnc,c->tn',aj,d)).max(axis=0)
        u,w=ex.controls(x+d);a=length(u[:,0])*2/(w[0,1]-w[0,0])
        return 1-(att+a)*ex.mass*SPECIFIC/ex.capacity
    _,_,producer=inputs(['research/studies/solar_shield_array/coupled_cycle.py','research/studies/solar_shield_array/coupled_stable.py',
        'research/studies/solar_shield_array/coupled_contacts.py','protection/dynamics/coupled_control.py'])
    producer['inputs'].update({n:digest(HERE/'results'/n) for n in ['coupled_benchmark.json','coupled_refinement_attempt.json']})
    out={**old,'producer':producer,'completed':False,'iterations':[],
        'previous_numerical_charge_s':charged,'cpu_ceiling_s':allowance,
        'retained_contacts_before':len(sets[0]),'retained_contacts_after':len(cuts),
        'method_changes':old['method_changes']+['Union all exposed contact branches from four rejected coupled replays; relinearize each at the current base.']}
    d,fit=solve_step(x,radius,gaps,gj,energy,power_margin,lambda d:t0+tj@d,
        old['checkpoint']['terminal_rms_limit'],margin=250.)
    rec=dict(iteration=0,base_x=x.tolist(),radius=radius,contacts=len(cuts),fit=fit,step=d.tolist());out['iterations'].append(rec)
    out['checkpoint']={**old['checkpoint'],'last_iteration':0}
    write('coupled_optimization.json',out,cpu,wall)
    if fit['success']:
        row,raw,_=ex.evaluate(x+d,service['state'][-1],21600.,horizon)
        row.update(save('retained_contacts_trial',**raw));rec['replay']=row
        ag=contact_values(raw['state'],raw['frames'],cuts);at=debt.residual(raw,x+d)[0]
        def merit(g,e,r):return e/1e14+.2*np.mean(r**2)+100*max(0,(250-float(g.min()))/1000)
        before=merit(gaps,energy(np.zeros(8))[0],t0);pred=merit(gaps+gj@d,energy(d)[0],t0+tj@d)
        total=service_energy+future+row['electrical_energy_J'];actual=merit(ag,total,at)
        defect=raw['state']-base['state']-np.einsum('tnic,c->tni',sj,d)
        decision=trust_decision(before,pred,actual,radius,float(length(defect[:,:,:3]).max()))
        prefix=service_energy+float(np.trapezoid(raw['power_W'][pm].sum(axis=1),times[pm]))
        old_event=old['selected']['clearance_events_s'][0] if old['selected'].get('clearance_events_s') else horizon
        new_event=row['clearance_events_s'][0] if row['clearance_events_s'] else horizon
        hard=bool(row['overloaded_members']==0 and prefix<=51.859e12 and total<=1e14 and
            np.sqrt(np.mean(at**2))<=old['checkpoint']['terminal_rms_limit']*(1+1e-6))
        if not hard or new_event<old_event-1.:
            decision['accepted']=False;decision['next_radius']=radius/2
        decision.update(hard_limits_pass=hard,first_floor_event_change_s=new_event-old_event)
        rec.update(decision=decision,before_merit=float(before),predicted_merit=float(pred),actual_merit=float(actual),
            actual_minimum_retained_signed_gap_m=float(ag.min()),prefix_energy_J=prefix,
            estimated_through_next_service_J=total,cycle_debt=debt.summary(raw,x+d))
        sj=secant_update(sj,d,raw['state']-base['state']);fj=secant_update(fj,d,raw['frames']-base['frames'])
        aj=secant_update(aj,d,raw['attitude_acceleration']-base['attitude_acceleration'])
        if decision['accepted']:
            x=x+d;out.update(selected=row,selected_x=x.tolist(),final_cycle_debt=debt.summary(raw,x))
        out['checkpoint'].update(x=x.tolist(),trust_radius=decision['next_radius'],best_raw_path=out['selected']['raw_path'],best_raw_sha256=out['selected']['raw_sha256'])
        print('retained contacts',len(cuts),decision,'gap',float(ag.min()),'cpu',time.process_time()-cpu,flush=True)
    out['checkpoint'].update(save('retained_contacts_optimizer_state',state_jac=sj,frame_jac=fj,attitude_jac=aj,
        cuts=cuts,gaps=gaps,gap_jac=gj,terminal_jac=tj,x=x,t=times))
    out['completed']=True;out['stopped_reason']='One retained-contact rebuild; preserve independent validation allocation'
    write('coupled_optimization.json',out,cpu,wall)


if __name__=='__main__':main()
