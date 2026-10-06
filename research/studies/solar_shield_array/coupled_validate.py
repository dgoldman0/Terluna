"""Event-stopped extension and independent tightened replay of selected controls."""
import json
import signal
import time
import numpy as np

from protection.dynamics.square_distance import closest_squares
from protection.dynamics.natural_pattern import local_coverage
from protection.dynamics.pattern_return import beam_and_rate_screen
from protection.dynamics.service_search import cadence_account
from .coupled_common import (HERE, ROOT, SPECIFIC, limit, inputs, write, save,
    read_raw, spline, state_at, digest, length)
from .coupled_stable import StableExperiment


def audit_clearance(sol,command,start,end,point_error):
    cache={};minimum=[float('inf'),None];lower=[];unresolved=[]
    def sample(t):
        if t not in cache:
            y=sol.sol(t).reshape(361,6)
            distance,pair,projection=closest_squares(y[:,:3],command(t).as_matrix())
            cache[t]=(distance,float(length(y[:,3:]-y[:,3:].mean(axis=0)).max()))
            if distance<minimum[0]:minimum[:]=[distance,dict(time_s=t,pair=pair.tolist(),projection_m=projection.tolist())]
        return cache[t]
    def guard(a,b):
        da,va=sample(a);db,vb=sample(b);dt=b-a
        bound=min(da,db)-(max(va,vb)+10000/np.sqrt(2)*np.deg2rad(.1))*dt-.15*dt**2/2-2*point_error
        if bound<300 and dt>.5:
            mid=(a+b)/2;guard(a,mid);guard(mid,b)
        else:
            lower.append(bound)
            if bound<100:unresolved.append([a,b,bound])
    dates=np.unique(np.r_[np.arange(start,end,60.),end])
    for a,b in zip(dates[:-1],dates[1:]):guard(float(a),float(b))
    return dict(start_s=start,end_s=end,minimum_sampled_distance_m=minimum[0],witness=minimum[1],
        conditional_interval_lower_bound_m=min(lower),point_error_allowance_m=point_error,
        unresolved_intervals=unresolved,adaptive_evaluations=len(cache),minimum_step_s=.5,
        acceleration_envelope_m_s2=.15,attitude_rate_envelope_deg_s=.1,
        passed=not unresolved and min(lower)>=100)


def main():
    limit();signal.alarm(650)
    cpu=time.process_time();wall=time.monotonic()
    bench=json.loads((HERE/'results/coupled_benchmark.json').read_text())
    opt=json.loads((HERE/'results/coupled_optimization.json').read_text())
    charged=opt['previous_numerical_charge_s']+opt['resources']['cpu_s']
    if 1800-charged<320:raise RuntimeError('Insufficient remaining independent-replay budget')
    ex=StableExperiment();x=np.array(opt['selected_x']);validated_case='selected_optimizer_state'
    if not any(r.get('decision',{}).get('accepted',False) for r in opt['iterations']):
        replayed=[r for r in opt['iterations'] if 'replay' in r]
        if replayed:
            x=np.array(replayed[-1]['replay']['x']);validated_case='last_rejected_contact_rebuild'
    cmd=ex.command(x)
    _,_,producer=inputs(['research/studies/solar_shield_array/coupled_cycle.py',
                        'research/studies/solar_shield_array/coupled_stable.py',
                        'research/studies/solar_shield_array/coupled_validate.py'])
    producer['inputs'].update({name:digest(HERE/'results'/name) for name in ['coupled_benchmark.json','coupled_optimization.json']})
    out=dict(schema='terluna.research.coupled-validation/1',producer=producer,completed=False,
        accepted_return=False,accepted_next_service=False,accepted_cycle=False,accepted_fleet=False,
        selected_x=x.tolist(),best_optimizer_x=opt['selected_x'],validated_case=validated_case,
        remaining_numerical_budget_s=1800-charged,
        initial_service_preserved=True,replacement_service_required_from_s=21600.,handovers_executed=0,
        held_array_comparison_TW=238.35,target_TW=23.835,capacity_relaxation_tiles=17.29e6,
        cadence=cadence_account(ex.arrival,mass_kg=float(ex.mass.sum())),runs=[])
    coarse_service=read_raw(bench['service'])
    coarse_row,coarse,coarse_sol=ex.evaluate(x,coarse_service['state'][-1],21600.,ex.arrival,stop=True)
    coarse_row.update(save('selected_coarse_execution',**coarse));out['coarse']=coarse_row
    write('coupled_validation.json',out,cpu,wall)
    print('coarse extension',coarse_row['end_s']/3600,coarse_row['cpu_s'],flush=True)
    fine_service_row,fs,fs_sol=ex.evaluate(x,ex.arrays['initial_epoch_state'],0.,21600.,fine=True,stop=True,step=60.)
    fine_service_row.update(save('selected_fine_service',**fs));out['runs'].append(fine_service_row)
    write('coupled_validation.json',out,cpu,wall)
    fine_row,fine,fine_sol=ex.evaluate(x,fs['state'][-1],21600.,ex.arrival,fine=True,stop=True,step=60.)
    fine_row.update(save('selected_fine_execution',**fine));out['runs'].append(fine_row)
    write('coupled_validation.json',out,cpu,wall)
    print('fine extension',fine_row['end_s']/3600,fine_row['cpu_s'],flush=True)
    comparisons=[]
    for old,new,sol,raw in [(bench['service'],fine_service_row,fs_sol,coarse_service),(coarse_row,fine_row,fine_sol,coarse)]:
        end=min(old['end_s'],new['end_s']);dates=np.unique(np.r_[np.arange(new['start_s'],end,30.),end])
        states=sol.sol(dates).T.reshape(-1,361,6);other=state_at(spline(raw),dates)
        error=float(length(states[:,:,:3]-other[:,:,:3]).max())+new['reconstruction_error_m']+old['reconstruction_error_m']
        comparisons.append(dict(maximum_position_with_reconstruction_m=error,
            maximum_velocity_difference_m_s=float(length(states[:,:,3:]-other[:,:,3:]).max()),
            end_time_difference_s=abs(old['end_s']-new['end_s'])))
    out['comparisons']=comparisons
    out['service_clearance']=audit_clearance(fs_sol,cmd,0.,21600.,comparisons[0]['maximum_position_with_reconstruction_m'])
    finish=fine_row['end_s'];requested=[43200.,float(opt['horizon_s'])]
    # The final stop at exactly 100 m cannot itself pass an error allowance.
    if finish<ex.arrival:requested.append(max(21600.,finish-120.))
    out['prefix_clearance']=[]
    for end in sorted(set(t for t in requested if t<=finish and t>21600)):
        out['prefix_clearance'].append(audit_clearance(fine_sol,cmd,21600.,end,comparisons[1]['maximum_position_with_reconstruction_m']))
    out['stopping_clearance']=audit_clearance(fine_sol,cmd, max(21600.,finish-120.),finish,comparisons[1]['maximum_position_with_reconstruction_m'])
    coverage=[]
    for t in np.linspace(0.,21600.,49):
        state=fs_sol.sol(t).reshape(361,6)
        # Preserve the previously assigned receiver, from the reproduced baseline.
        receiver=state_at(spline(coarse_service),t).mean(axis=0)
        coverage.append(local_coverage(state,ex.env.at(t),receiver,limb_count=16,interior=8,rotation=.137))
    out['initial_service_coverage']=dict(dates=49,sources_per_date=25,clear_side_m=9890.,
        minimum_source_coverage=min(r['minimum_source_coverage'] for r in coverage),
        passed=all(r['minimum_source_coverage']>=1-1e-9 for r in coverage),
        spectral_model='Unchanged inherited ideal filter and source-bound stored spectral products; angular material response remains outside this dynamics model.')
    for row,raw in [(fine_service_row,fs),(fine_row,fine)]:
        row['beam_and_rate']=beam_and_rate_screen(ex.env,raw['t'],raw['state'],cmd)
    time_all=np.r_[fs['t'],fine['t'][1:]]
    powers=np.r_[fs['power_W'],fine['power_W'][1:]]
    light=np.r_[fs['first_intercept_W'],fine['first_intercept_W'][1:]]
    att=np.r_[fs['attitude_acceleration'],fine['attitude_acceleration'][1:]]
    translation=np.r_[fs['translation'],fine['translation'][1:]]
    energy=float(np.trapezoid(powers.sum(axis=1),time_all));mask=time_all<=43200
    out['account']=dict(duration_s=finish,electrical_energy_J=energy,
        twelve_hour_energy_J=float(np.trapezoid(powers[mask].sum(axis=1),time_all[mask])) if finish>=43200 else None,
        mean_electrical_W=energy/finish,
        first_intercept_J=float(np.trapezoid(light.sum(axis=1),time_all)),
        redirected_J=float(np.trapezoid(light.sum(axis=1),time_all))*ex.env.central_fraction,
        attitude_energy_J=float(np.trapezoid((att*ex.mass*SPECIFIC).sum(axis=1),time_all)),
        translation_energy_J=float(np.trapezoid((length(translation)*ex.mass*SPECIFIC).sum(axis=1),time_all)),
        ideal_exhaust_kg=energy*1.4/30000**2,propulsion_waste_heat_J=.3*energy,
        maximum_individual_power_W=float(powers.max()),minimum_installed_to_peak=float(np.min(ex.capacity/powers.max(axis=0))),
        overloaded_members=int(np.count_nonzero(powers.max(axis=0)>ex.capacity)),per_member_peak_W=powers.max(axis=0).tolist(),
        per_member_installed_W=ex.capacity.tolist(),mass_kg=float(ex.mass.sum()),hardware_changed=False,
        actual_generation_W=None,actual_bus_delivery_W=None,actual_storage_state_J=None,
        actual_delivery_losses_J=None,actual_unmet_demand_J=None,full_cycle_energy_J=None,
        supply_evidence='No generator/storage placement specified. Per-date demand and optical source histories are stored; electricity availability and storage chronology remain unknown.',
        **save('executed_chronology',t=time_all,power_W=powers,first_intercept_W=light,
            redirected_W=light*ex.env.central_fraction,attitude_acceleration=att,translation=translation))
    out['terminal']=dict(time_s=finish,state_m_m_s=fine['state'][-1].tolist(),
        attitude_xyzw=cmd(finish).as_quat().tolist(),body_angular_velocity_rad_s=cmd(finish,1).tolist(),
        state_frame='Moon-relative ICRF; actual terminal state, no reset',
        clearance_witness=out['stopping_clearance']['witness'])
    out['safe_prefix_extended']=any(r['passed'] and r['end_s']>13.799825*3600 for r in out['prefix_clearance'])
    out['next_service_executed']=False;out['subsequent_departure_executed']=False
    out['completed']=True;write('coupled_validation.json',out,cpu,wall)
    print('account TJ',energy/1e12,'prefix12',out['account']['twelve_hour_energy_J']/1e12,
        'safe extension',out['safe_prefix_extended'],'cpu',out['resources']['cpu_s'],flush=True)


if __name__=='__main__':main()
