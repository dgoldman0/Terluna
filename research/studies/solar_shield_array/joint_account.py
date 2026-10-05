"""Full executed sequence account and optimistic local-supply rejection test."""
import time,signal,resource
import numpy as np
from scipy.interpolate import CubicSpline,interp1d
from protection.dynamics.eclipse_parallel import load
from protection.dynamics.attitude_load import rigid_square_load
from protection.dynamics.cycling import ray_geometry,visible_sun
from protection.dynamics.optical import length
from engineering.electromagnetic import collector,storage,cable,radiator
from engineering.shield_storage import required_capacity,dispatch
from shared import constants as K
from .joint_common import *
from .joint_bridge import EXTRA


def main():
    signal.alarm(220);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();cpu=time.process_time()
    extra=[*EXTRA,'research/studies/solar_shield_array/joint_account.py','engineering/electromagnetic.py','engineering/shield_storage.py']
    p,producer=setup(extra)
    name='joint_corridor_108_fine.json';parent=json.loads((HERE/'results'/name).read_text());producer['inputs'][name]=digest(HERE/'results'/name)
    inverse=json.loads((HERE/'results/joint_bridge.json').read_text())['cases'][0];inv=load_raw(inverse)
    producer['inputs']['joint_bridge.json']=digest(HERE/'results/joint_bridge.json')
    controls=CubicSpline(inv['t'],inv['acceleration'],axis=0)
    env=environment(3.);arrival=parent['arrival_s'];end=parent['runs'][-1]['end_s'];cmd=command(env,end,arrival)
    ratio=inv['mass_ratio'];mass=ratio*5e6;specific=30000/(2*.7*np.cos(np.pi/4));capacity=(ratio-1.2)*5e6*300
    out=dict(schema='terluna.research.joint-account/1',producer=producer,stages=[],accepted_return=False,
        accepted_cycle=False,accepted_fleet=False,actual_generation_W=None,actual_delivered_bus_W=None,
        actual_delivered_export_W=None,full_operating_cycle_cost_J=None,
        scope='Measured ideal-command sequence cost. Installed actuator limits and supply are separate failed gates. Conditional collector/storage estimates are not installed, and are not added to the optical first-intercept total.')
    combined_t=[];combined_p=[];combined_y=[];peak=np.zeros(361);cost_total=0.;optical_total=0.;att_total=np.zeros(361);dv_total=att_total.copy()
    for entry in parent['runs']:
        raw=load_raw(entry);path=spline(raw);evaluations=[]
        for suns,step,grid in [(8,120.,16),(16,240.,32)]:
            ts=np.unique(np.r_[np.arange(entry['start_s'],entry['end_s'],step),entry['end_s']])
            optical=[];radiation=[];gravity=[]
            for t in ts:
                y=state_at(path,t);a=load(env,t,y,cmd(t).as_matrix(),suns=suns,torque=True,grid=grid)
                optical.append(a['optical']['first_intercept_bolometric_equivalent_W']);radiation.append(a['radiation_torque_N_m']);gravity.append(a['gravity_torque_per_mass'])
            optical=np.array(optical);disturbance=np.array(gravity)+np.array(radiation)/mass[None,:,None]
            dense=np.unique(np.r_[np.arange(ts[0],ts[-1],30.),ts[-1]])
            nominal=rigid_square_load(cmd,dense)['torque_per_mass']
            torques=nominal[:,None,:]-interp1d(ts,disturbance,axis=0)(dense)
            couple=2*abs(torques).sum(axis=-1)/10000
            u=np.array([controls(t) if 21600<t<arrival else np.zeros((361,3)) for t in dense])
            powers=(couple+length(u))*mass*specific
            dv=np.trapezoid(length(u),dense,axis=0);att=np.trapezoid(couple,dense,axis=0)
            optical_E=float(np.trapezoid(optical.sum(axis=1),ts));energy=float(np.trapezoid(powers.sum(axis=1),dense))
            item=dict(suns=suns,load_step_s=step,partial_eclipse_area_nodes_per_axis=grid,
                propulsion_energy_J=energy,translation_mean_m_s=float(dv.mean()),attitude_mean_m_s=float(att.mean()),
                first_intercept_J=optical_E,redirected_band_J=optical_E*env.central_fraction,
                mean_first_intercept_W=optical_E/(ts[-1]-ts[0]),mean_redirected_W=optical_E*env.central_fraction/(ts[-1]-ts[0]),
                maximum_individual_power_W=float(powers.max()),sum_individual_peak_W=float(powers.max(axis=0).sum()),
                maximum_torque_N_m=float(length(torques*mass[None,:,None]).max()),
                maximum_command_acceleration_m_s2=float(length(u).max()))
            evaluations.append(item)
            if suns==16:
                peak=np.maximum(peak,powers.max(axis=0));cost_total+=energy;optical_total+=optical_E;att_total+=att;dv_total+=dv
                sl=slice(None) if not combined_t else slice(1,None)
                combined_t.extend(dense[sl]);combined_p.extend(powers[sl]);combined_y.extend(path(dense[sl]))
                saved=RUN/('account_'+entry['stage']+'.npz')
                np.savez_compressed(saved,t=dense,power_W=powers,attitude_equivalent_acceleration=couple,
                    translation_acceleration=u,optical_times=ts,first_intercept_W=optical,
                    redirected_band_W=optical*env.central_fraction,gravity_torque_per_mass=np.array(gravity),radiation_torque_N_m=np.array(radiation))
                item.update(raw_path=str(saved.relative_to(ROOT)),raw_sha256=digest(saved))
        out['stages'].append(dict(stage=entry['stage'],start_s=entry['start_s'],end_s=entry['end_s'],evaluations=evaluations,
            time_source_relative_energy_difference=abs(evaluations[0]['propulsion_energy_J']-energy)/energy,
            time_source_relative_optical_difference=abs(evaluations[0]['first_intercept_J']-optical_E)/optical_E))
        write('joint_account.json',out,before,cpu);print('account',entry['stage'],time.monotonic()-before,flush=True)
    times=np.array(combined_t);powers=np.array(combined_p);positions=np.array(combined_y)
    out['measured_sequence']=dict(duration_s=float(times[-1]),propulsion_energy_J=cost_total,
        initial_service_and_one_return_energy_J=cost_total-out['stages'][-1]['evaluations'][-1]['propulsion_energy_J'],
        mean_power_W=cost_total/times[-1],translation_mean_m_s=float(dv_total.mean()),attitude_mean_m_s=float(att_total.mean()),
        first_intercept_J=optical_total,redirected_band_J=optical_total*env.central_fraction,
        ideal_exhaust_kg=2*.7*cost_total/30000**2,propulsion_waste_heat_J=(1-.7)*cost_total,
        overloaded_members=int(np.count_nonzero(peak>capacity)),minimum_installed_to_peak=float(np.min(capacity/peak)),
        sum_required_individual_peak_W=float(peak.sum()),simultaneous_peak_W=float(powers.sum(axis=1).max()),
        per_member_peak_W=peak.tolist(),per_member_installed_W=capacity.tolist(),
        added_actuator_mass_for_25_percent_margin_kg=float(np.maximum(1.25*peak-capacity,0).sum()/300),
        resized_actuator_mass_propagated=False,constant_mass_approximation=True)
    # Optimistic additional local collectors: fresh unshadowed sunlight, except
    # measured finite-Sun body eclipses. Inter-array shading and actual placement
    # can only make this specified supply screen worse; no combined-light sum.
    designs=[collector(1.25*x) for x in peak];areas=np.array([x['area_m2'] for x in designs]);gen=[];eclipse=[]
    for t,q in zip(times,positions):
        rays=ray_geometry(q,env.at(t));visibility=visible_sun(rays['sun'],rays['earth'],rays['moon'])
        incident=areas*K.SOLAR_CONSTANT*(K.AU/length(rays['sun']))**2*visibility
        gen.append(incident*.30*.95*.99)
        if visibility.min()<.999999:eclipse.append(float(t))
    gen=np.array(gen);buffers=[storage(60*x,x) for x in peak]
    nominal=np.array([x['battery_mass_kg']*200*3600 for x in buffers]);rating=np.array([x['battery_mass_kg']*1000 for x in buffers])
    short=dispatch(times,powers,gen,nominal,rating);required=required_capacity(times,powers,gen)
    minimum_mass=np.maximum(required/(200*3600),peak/1000);full=dispatch(times,powers,gen,minimum_mass*200*3600,minimum_mass*1000)
    extra=RUN/'local_supply_upper_bound.npz'
    np.savez_compressed(extra,t=times,load_W=powers,generation_bus_upper_W=gen,
        short_buffer_state_J=short['state_J'],eclipse_sized_state_J=full['state_J'],
        short_buffer_nominal_J=nominal,required_nominal_J=required)
    emitted=gen/(.30*.95*.99)
    pvheat=emitted*(.9-.30);pmadloss=emitted*.30*.05;lineloss=emitted*.30*.95*.01
    out['local_supply_upper_bound']=dict(assumptions=dict(PV_efficiency=.30,PMAD_efficiency=.95,cable_efficiency=.99,
        array_W_kg=300,battery_Wh_kg=200,battery_W_kg=1000,charge_efficiency=.95,discharge_efficiency=.95,usable_depth=.8,
        no_mutual_collector_shading=True,collector_placement_validated=False,collector_mass_and_photon_forces_replayed=False),
        generation_hardware_mass_kg=float(sum(x['installed_mass_kg'] for x in designs)),collector_area_m2=float(areas.sum()),
        incident_photon_absorption_force_N=float(sum(x['absorbed_photon_force_N'] for x in designs)),
        buffer_60s_design_mass_kg=float(sum(x['battery_mass_kg'] for x in buffers)),
        buffer_60s_unserved_J=float(short['unserved_J'].sum()),first_unserved_s=short['first_unserved_s'],
        eclipse_first_last_s=[min(eclipse),max(eclipse)] if eclipse else [],
        optimistic_eclipse_sized_battery_mass_kg=float(minimum_mass.sum()),
        required_nominal_storage_J=float(required.sum()),sized_unserved_J=float(full['unserved_J'].sum()),
        additional_conductor_mass_kg=float(sum(cable(1.25*x)['conductor_mass_kg'] for x in peak)),
        maximum_battery_radiator_mass_kg=float(sum(radiator(x*(1/.95-1))['mass_kg'] for x in peak)),
        available_bus_upper_energy_J=float(np.trapezoid(gen.sum(axis=1),times)),
        PV_waste_heat_J=float(np.trapezoid(pvheat.sum(axis=1),times)),
        PMAD_loss_J=float(np.trapezoid(pmadloss.sum(axis=1),times)),
        cable_loss_J=float(np.trapezoid(lineloss.sum(axis=1),times)),
        storage_loss_J=float(full['storage_loss_J'].sum()),curtailed_bus_J=float(full['curtailed_J'].sum()),
        initial_stored_J=float(full['state_J'][0].sum()),final_stored_J=float(full['state_J'][-1].sum()),
        maximum_energy_balance_residual_J=float(abs(full['balance_residual_J']).max()),
        raw_path=str(extra.relative_to(ROOT)),raw_sha256=digest(extra))
    out['completed']=True;write('joint_account.json',out,before,cpu)
if __name__=='__main__':main()
