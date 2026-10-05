"""Read measured histories and execute a bounded EM infrastructure screen."""
import json
from pathlib import Path
import resource
import signal
import time
import numpy as np
from scipy.interpolate import CubicHermiteSpline, interp1d
from scipy.spatial.transform import Rotation, RotationSpline
from scipy.optimize import brentq

from shared import constants as K
from shared.provenance import constants_used
from engineering.electromagnetic import (square, mutual, field_and_potential,
    coil_hardware, link, collector, storage, cable, plasma, rigidity, radiator)
from protection.dynamics.attitude_load import rigid_square_load
from protection.dynamics.cycle_control import return_command
from protection.dynamics.pattern_return import smooth_arc
from .packing_run import full_command
from .packing_screen import parents, load_raw
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE, SCENARIO

RUN=ROOT/'research/runs/solar_shield_array/electromagnetic'
GROUP=np.array([160,161,162,179,180,181,198,199,200])


def plain(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value))


def integrate(t,value):return np.trapezoid(value,t,axis=0)


def main():
    before=time.monotonic();signal.alarm(900)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    RUN.mkdir(parents=True,exist_ok=True)
    names=['closure_validation.json','closure_sizing.json','closure_budget.json',
           'packing_validation.json','repair_trial_local_refined_fine.json','repair_budget.json']
    p,sources=parents(names)
    for file in ['engineering/electromagnetic.py','research/studies/solar_shield_array/electromagnetic_run.py',
                 'research/studies/solar_shield_array/electromagnetic_sources.json']:
        sources[file]=digest(ROOT/file)
    out=dict(schema='terluna.research.electromagnetic-feasibility/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),
            inputs={n:digest(HERE/'results'/n) for n in names}),
        evidence='Executed vacuum line integrals, analytic engineering sensitivities and saved-trajectory accounting. No new trajectory, plasma solution, accepted design or recurrence.',
        accepted_trajectory=False,accepted_return=False,accepted_cycle=False,accepted_fleet=False,
        original_target_W=23.835e12,new_target_adopted=False,
        raw_inputs={},group_members=GROUP,stages=[],force_samples=[],pair_samples=[],group_fields=[],
        units='SI; numeric fields use their suffix. Currents are ampere-turns.',
        reading_rule='Compare engineering cases only at their stated scope; no mass, force or collection change has been propagated.')
    env=environment(3.)
    command=full_command(env,p['packing_validation.json']['schedule'])
    ret=p['repair_trial_local_refined_fine.json']
    retcmd=return_command(env,ret['run']['start_s'],ret['arrival_s'],ret['arrival_turn_s'])
    trial=next(a for a in p['repair_budget.json']['trials'] if a['source']=='repair_trial_local_refined_fine.json')
    stages=[('service',p['closure_validation.json']['runs'][0],command,
             p['closure_budget.json']['datasets']['hardware_service']['evaluations'][0],None),
            ('departure',p['closure_validation.json']['runs'][1],command,
             p['closure_budget.json']['datasets']['hardware_departure']['evaluations'][0],None),
            ('return',ret['run'],retcmd,trial['evaluations'][0],trial['grazing_refinement'][-1])]
    ratio=np.array(p['closure_sizing.json']['per_member_total_to_optical_mass_ratio'])
    mass=ratio*5e6;capacity=(ratio-1.2)*5e6*300
    prop=SCENARIO['propulsion'];w_per_n=prop['exhaust_velocity_m_s']/(2*prop['efficiency']*np.cos(np.deg2rad(prop['cant_deg'])))
    all_peak=np.zeros(361);loaded=[];snapshots=[];last=None
    for label,entry,cmd,load_entry,opt_entry in stages:
        raw=load_raw(entry);load=load_raw(load_entry);opt=load_raw(opt_entry) if opt_entry else load
        for e in [entry,load_entry,*([opt_entry] if opt_entry else [])]:out['raw_inputs'][e['raw_path']]=e['raw_sha256']
        if last is not None:assert np.max(abs(raw['state'][0]-last))<1e-6
        last=raw['state'][-1]
        t=np.unique(np.r_[np.arange(raw['t'][0],raw['t'][-1],2.),raw['t'][-1]])
        nominal=rigid_square_load(cmd,t)['torque_per_mass']
        disturbed=interp1d(load['t'],load['gravity_torque_N_m_kg']+load['radiation_torque_optical_N_m_kg']/ratio[None,:,None],axis=0)(t)
        torque=(nominal[:,None,:]-disturbed)*mass[None,:,None]
        shape=np.array([smooth_arc(ti,raw['windows']) for ti in t])
        acceleration=np.einsum('tk,nkj->tnj',shape,raw['controls'])
        force=acceleration*mass[None,:,None]
        translation=np.linalg.norm(force,axis=2)
        attitude=2*abs(torque).sum(axis=2)/10000.
        demand=(translation+attitude)*w_per_n
        # Conditional generation only: previous 50% capture, 40% conversion.
        optical=interp1d(opt['t'],opt['redirected_band_optical_W'],axis=0)(t)
        generation=.2*optical
        peak=demand.max(axis=0);all_peak=np.maximum(all_peak,peak)
        assert abs(float(integrate(t,demand.sum(axis=1)))/load_entry['cost']['propulsion_energy_J']-1)<1e-7
        s=dict(stage=label,start_s=float(t[0]),end_s=float(t[-1]),
            propulsion_energy_J=float(integrate(t,demand.sum(axis=1))),
            translation_energy_J=float(integrate(t,translation.sum(axis=1))*w_per_n),
            attitude_energy_J=float(integrate(t,attitude.sum(axis=1))*w_per_n),
            per_tile_peak_W=peak,peak_total_W=float(demand.sum(axis=1).max()),
            group_peak_total_W=float(demand[:,GROUP].sum(axis=1).max()),
            group_sum_individual_peaks_W=float(peak[GROUP].sum()),
            group_energy_J=float(integrate(t,demand[:,GROUP].sum(axis=1))),
            group_max_translation_force_N=float(translation[:,GROUP].max()),
            group_max_attitude_torque_N_m=float(np.linalg.norm(torque[:,GROUP],axis=2).max()),
            mean_conditional_generation_W=float(integrate(t,generation.sum(axis=1))/(t[-1]-t[0])))
        if label=='return':
            deficit=np.maximum(demand-generation,0)
            pooled=np.maximum(demand.sum(axis=1)-generation.sum(axis=1),0)
            s.update(net_energy_deficit_J=float(integrate(t,(demand-generation).sum(axis=1))),
                sum_local_positive_deficit_J=float(integrate(t,deficit.sum(axis=1))),
                pooled_positive_deficit_J=float(integrate(t,pooled)),
                peak_pooled_deficit_W=float(pooled.max()),
                group_deficit_J=float(integrate(t,deficit[:,GROUP].sum(axis=1))),
                per_tile_deficit_J=integrate(t,deficit),per_tile_peak_deficit_W=deficit.max(axis=0),
                pooled_translational_force_lower_bound_impulse_N_s=float(integrate(t,np.linalg.norm(force.sum(axis=1),axis=1))),
                sum_translation_impulse_N_s=float(integrate(t,translation.sum(axis=1))),
                group_external_force_lower_bound_impulse_N_s=float(integrate(t,np.linalg.norm(force[:,GROUP].sum(axis=1),axis=1))),
                group_sum_translation_impulse_N_s=float(integrate(t,translation[:,GROUP].sum(axis=1))),
                generation_scope='Conditional 0.5 optical capture times 0.4 conversion, no implemented collector or routing.')
            return_demand=demand[:,GROUP].copy();return_time=t.copy()
            loaded.append((t,demand[:,GROUP].copy(),generation[:,GROUP].copy()))
        else:loaded.append((t,demand[:,GROUP].copy(),generation[:,GROUP].copy()))
        path=RUN/(label+'_group.npz')
        np.savez_compressed(path,t=t,member_id=GROUP,propulsion_W=demand[:,GROUP],
            translational_force_N=force[:,GROUP],attitude_torque_N_m=torque[:,GROUP],
            redirected_optical_W=optical[:,GROUP],conditional_generation_W=generation[:,GROUP],
            total_propulsion_W=demand.sum(axis=1),total_conditional_generation_W=generation.sum(axis=1))
        s['raw_output']=dict(path=str(path.relative_to(ROOT)),sha256=digest(path));out['stages'].append(s)
        spline=CubicHermiteSpline(raw['t'],raw['state'][:,:,:3],raw['state'][:,:,3:])
        dates=([0.,21600.] if label=='service' else [32400.,43200.] if label=='departure' else [46800.,50400.,float(t[-1])])
        for ti in dates:
            q=spline(ti);v=spline(ti,1);frame=cmd(ti).as_matrix()
            local=(q-q[180])@frame
            pairs=local[GROUP][:,None]-local[GROUP][None,:]
            dist=np.linalg.norm(pairs,axis=2);np.fill_diagonal(dist,np.inf)
            relative_v=v[GROUP][:,None]-v[GROUP][None,:]
            inertial_r=q[GROUP][:,None]-q[GROUP][None,:]
            rate=np.linalg.norm(np.cross(inertial_r,relative_v),axis=2)/dist**2
            snapshots.append((ti,local,frame))
            acc=np.einsum('k,nkj->nj',smooth_arc(ti,raw['windows']),raw['controls'])
            out['force_samples'].append(dict(time_s=ti,group_min_centre_m=float(dist.min()),
                group_max_centre_m=float(dist[np.isfinite(dist)].max()),
                max_pair_tracking_rate_rad_s=float(rate.max()),
                max_correction_force_N=float(np.linalg.norm(acc*mass[:,None],axis=1).max()),
                group_required_force_N=(acc*mass[:,None])[GROUP],group_positions_body_m=local[GROUP],
                group_mass_kg=mass[GROUP]))
    out['reference']=dict(total_mass_kg=float(mass.sum()),group_mass_kg=float(mass[GROUP].sum()),
        installed_propulsion_capacity_W=float(capacity.sum()),per_tile_capacity_W=capacity,
        all_stage_peak_W=all_peak,overloaded_members=int(np.count_nonzero(all_peak>capacity)),
        group_overloaded_members=int(np.count_nonzero(all_peak[GROUP]>capacity[GROUP])),
        extra_actuator_mass_kg=float(np.maximum(1.25*all_peak-capacity,0).sum()/300),
        group_extra_actuator_mass_kg=float(np.maximum(1.25*all_peak[GROUP]-capacity[GROUP],0).sum()/300),
        caveat='Installed propulsion/power allocation is not measured generator output or receiver capacity.')
    retstage=out['stages'][-1]
    out['storage_cases']=[dict(wh_kg=wh,**storage(retstage['pooled_positive_deficit_J'],retstage['peak_pooled_deficit_W'],wh_kg=wh)) for wh in [150,200,400]]
    peak_group=all_peak[GROUP]
    # Local generation at all-stage per-member peak, plus one minute interruption buffer.
    local_collectors=[collector(1.25*x) for x in peak_group]
    local_storage=[storage(60*x,x) for x in peak_group]
    out['local_generation']=dict(per_tile=local_collectors,
        group_collector_mass_kg=sum(x['installed_mass_kg'] for x in local_collectors),
        group_area_m2=sum(x['area_m2'] for x in local_collectors),
        group_60s_buffer_mass_kg=sum(x['battery_mass_kg']+x['heat_rejection']['mass_kg'] for x in local_storage),
        group_bus_conductors_kg=sum(cable(x)['conductor_mass_kg'] for x in peak_group),
        group_extra_actuator_mass_kg=out['reference']['group_extra_actuator_mass_kg'],
        group_propellant_saved_kg=0.,scope='Peak-sized supplemental sunlight generation and 60 s ride-through; no credit for the unbuilt redirected-band collector. No propulsion saving is claimed.')
    out['link_cases']=[link(power,distance,kind) for kind in ['microwave','laser']
        for power in [30e6,100e6] for distance in [15000.,100000.,1000000.]]
    out['bus_cases']=[cable(100e6,voltage=v) for v in [1e4,1e5,1e6]]
    out['hub_power']=[]
    for kind in ['microwave','laser']:
        links=[link(x,100000.,kind) for x in peak_group]
        eta=links[0]['dc_efficiency'];peak_source=max(x['group_peak_total_W'] for x in out['stages'])/eta
        generation=collector(1.25*peak_source)
        out['hub_power'].append(dict(kind=kind,links=links,generation=generation,
            group_link_mass_kg=sum(x['link_installed_mass_kg'] for x in links),
            total_new_power_mass_kg=generation['installed_mass_kg']+sum(x['link_installed_mass_kg'] for x in links)+out['local_generation']['group_60s_buffer_mass_kg']+out['local_generation']['group_bus_conductors_kg'],
            mean_source_bus_W=float(sum(x['group_energy_J'] for x in out['stages'])/out['stages'][-1]['end_s']/eta),
            scope='One source hub with nine independently peak-rated links; generation sized from simultaneous demand plus 25%. 100 km is a placement sensitivity, not a propagated hub orbit.'))
    # Physical square coils at measured geometries. No dipole expansion here.
    base=square();current=1e5;coil=coil_hardware(current);out['coil_reference']=coil
    for ti,local,frame in snapshots:
        pairlist=[(180,int(i)) for i in GROUP if i!=180]+[(325,347),(151,186)]
        for i,j in pairlist:
            d=local[j]-local[i];target=square(centre=d)
            low=mutual(base,target,32);high=mutual(base,target,64)
            rel=np.linalg.norm(high['force_N_A2']-low['force_N_A2'])/max(np.linalg.norm(high['force_N_A2']),1e-30)
            n=64
            while rel>1e-4 and n<512:
                low=high;n*=2;high=mutual(base,target,n)
                rel=np.linalg.norm(high['force_N_A2']-low['force_N_A2'])/max(np.linalg.norm(high['force_N_A2']),1e-30)
            reverse=mutual(target,base,n)
            torque_balance=high['torque_N_m_A2']+reverse['torque_N_m_A2']+np.cross(d,high['force_N_A2'])
            coefficient=np.linalg.norm(high['force_N_A2'])
            out['pair_samples'].append(dict(time_s=ti,pair=[i,j],displacement_body_m=d,
                force_at_100kA_N=high['force_N_A2']*current**2,
                torque_on_target_at_100kA_N_m=high['torque_N_m_A2']*current**2,
                mutual_H=high['mutual_H'],quadrature_nodes_per_edge=n,relative_refinement=rel,
                action_reaction_relative=np.linalg.norm(high['force_N_A2']+reverse['force_N_A2'])/max(coefficient,1e-30),
                angular_momentum_relative=np.linalg.norm(torque_balance)/max(coefficient*np.linalg.norm(d),1e-30),
                ampere_turns_for_100N=float(np.sqrt(100/coefficient))))
        # All 361 loops energized identically versus only the selected nine.
        net=np.zeros((9,3));torq=np.zeros((9,3));outside=np.zeros((9,3));worst_refinement=0.
        group_set=set(GROUP.tolist())
        for k,i in enumerate(GROUP):
            for j in range(361):
                if i==j:continue
                source=square(centre=local[j]-local[i]);low=mutual(source,base,24)
                m=mutual(source,base,48);n=48
                delta=np.linalg.norm(m['force_N_A2']-low['force_N_A2'])/max(np.linalg.norm(m['force_N_A2']),1e-30)
                while delta>1e-4 and n<192:
                    low=m;n*=2;m=mutual(source,base,n)
                    delta=np.linalg.norm(m['force_N_A2']-low['force_N_A2'])/max(np.linalg.norm(m['force_N_A2']),1e-30)
                worst_refinement=max(worst_refinement,delta)
                f=m['force_N_A2']*current**2;tq=m['torque_N_m_A2']*current**2
                if j in group_set:net[k]+=f;torq[k]+=tq
                else:outside[k]+=f
        out['group_fields'].append(dict(time_s=ti,group_net_force_N=net,group_torque_N_m=torq,
            force_from_other_352_N=outside,sum_group_internal_force_N=net.sum(axis=0),
            maximum_pair_force_relative_refinement=worst_refinement,
            scope='24/48 nodes with refinement through 192 as needed; all coils parallel at +100 kA. Cross-talk diagnostic, not independently commandable pair forces.'))
    out['coil_scaling']=[]
    for distance in [10000.,15000.,30000.,100000.,300000.]:
        a=mutual(base,square(centre=[0,0,distance]),64)
        coefficient=np.linalg.norm(a['force_N_A2'])
        for force in [10.,100.,1000.,5000.]:
            amps=np.sqrt(force/coefficient);hardware=coil_hardware(amps)
            out['coil_scaling'].append(dict(coaxial_distance_m=distance,required_pair_force_N=force,
                **hardware,force_dipole_ratio=float(coefficient/(3*K.VACUUM_PERMEABILITY*1e16/(2*np.pi*distance**4)))))
    # Circuit near-field screens for a measured nearest group neighbor at hour 12.
    candidates=[x for x in out['pair_samples'] if x['time_s']==43200 and x['pair'][0]==180]
    nearest=min(candidates,key=lambda x:np.linalg.norm(x['displacement_body_m']))
    M=abs(nearest['mutual_H']);L=coil['inductance_H'];omega=2*np.pi*100
    amps=np.sqrt(100e6/(omega*M))
    out['induction']=dict(pair=nearest['pair'],mutual_H=M,self_H=L,coupling=M/L,
        frequency_Hz=100.,delivered_W=100e6,rms_ampere_turns=float(amps),
        per_coil_reactive_VA=float(omega*L*amps**2),per_coil_voltage_V=float(omega*L*amps),
        required_Q_for_10kW_cold_loss_per_coil=float(omega*L*amps**2/10000),
        warm_copper_loss_W=float(amps**2*1.68e-8*40000/(amps/1e8)),
        pair_hardware=[coil_hardware(amps)],
        scope='Ideal resonant RMS power P=omega M I1 I2; 10 kW cold loss requires 2 MW input per coil at COP .005. No measured HTS AC-loss or tuning performance. One winding hardware row must be doubled for the pair.')
    # Hub square, 100 km on a side: extended field and force, not a central dipole.
    hub=square(100000.,centre=[0,0,100000.]);hf=mutual(hub,base,64)
    out['magnetic_hub']=dict(side_m=100000.,normal_offset_m=100000.,
        coil=coil_hardware(1e6,perimeter=400000.,bundle_radius=1.),
        force_on_100kA_tile_N=hf['force_N_A2']*1e11,
        hub_reaction_N=-hf['force_N_A2']*1e11,
        scope='One-million ampere-turn hub and 100kA tile; reaction moves hub. A stationary hub needs paid external force; a coasting hub accumulates state debt.')
    # Solar-wind diagnostic contours from extended square fields.
    out['plasma_conditions']=[plasma(5,400,10,5),plasma(20,800,20,15),plasma(100,800,30,30)]
    out['field_pressure_diagnostics']=[]
    for side,amps,label in [(10000.,1e5,'tile_100kA'),(100000.,1e6,'hub_1MA')]:
        loop=square(side);moment=amps*side**2
        for pressure in [2e-9,20e-9,100e-9]:
            required=np.sqrt(2*K.VACUUM_PERMEABILITY*pressure)
            def strength(r,axis):return np.linalg.norm(field_and_potential(np.array(axis)*r,loop)[0][0])*amps
            axial=brentq(lambda r:strength(r,[0,0,1])-required,side*.001,side*1000)
            equatorial=brentq(lambda r:strength(r,[1,0,0])-required,side*.8,side*1000)
            wind_force=pressure*np.pi*equatorial**2
            out['field_pressure_diagnostics'].append(dict(label=label,pressure_Pa=pressure,required_B_T=required,
                axial_distance_m=axial,equatorial_distance_m=equatorial,
                absorbing_disk_momentum_N=wind_force,perfect_reversal_upper_N=2*wind_force,
                moment_A_m2=moment,imf_5nT_max_torque_N_m=moment*5e-9,
                scope='Vacuum |B| pressure crossings only. Circular cross-section momentum estimate assumes a barrier, not proven plasma exclusion or a protected downstream wake.'))
    out['particle_bending']=[dict(energy_eV=e,one_radian_BL_T_m=rigidity(e)) for e in [835.,1e6,1e8,1e9]]
    out['wind_force_area']=[dict(force_N=f,pressure_Pa=pressure,
        minimum_reversal_radius_m=float(np.sqrt(f/(2*np.pi*pressure)))) for f in [100.,1000.,5000.]
        for pressure in [2e-9,20e-9,100e-9]]
    # Wind/plasma interactions of a closed coil: uniform E/B gives torque, no net DC Lorentz thrust.
    out['coil_storage_comparison']=dict(all_361_stored_energy_J=361*coil['stored_energy_J'],
        fraction_of_return_deficit=361*coil['stored_energy_J']/retstage['net_energy_deficit_J'],
        current_needed_per_tile_A=float(np.sqrt(2*retstage['net_energy_deficit_J']/(361*L))))
    out['replacement']=[dict(lifetime_years=life,kg_s_per_million_kg=1e6/(life*K.JULIAN_YEAR_DAYS*K.JULIAN_DAY),
        manufacture_W_per_million_kg=1e6*1e8/(life*K.JULIAN_YEAR_DAYS*K.JULIAN_DAY)) for life in [5.,20.,100.]]
    out['resources']=dict(elapsed_s=time.monotonic()-before,max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
        raw_output_MiB=sum(f.stat().st_size for f in RUN.glob('*.npz'))/2**20,stage_limit_s=900,total_limit_s=1200,threads=1,address_space_GiB=2)
    out['completed']=True
    if any(digest(ROOT/f)!=h for f,h in sources.items()):raise ValueError('Source changed during run')
    (HERE/'results/electromagnetic.json').write_text(json.dumps(out,default=plain,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(resources=out['resources'],overloads=out['reference']['overloaded_members'],
        group_collector_mass_kg=out['local_generation']['group_collector_mass_kg'],
        return_deficit_TJ=retstage['net_energy_deficit_J']/1e12),default=plain),flush=True)


if __name__=='__main__':main()
