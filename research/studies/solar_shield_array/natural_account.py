"""Measured optical, electrical and exhaust accounts for executed prefixes only."""
import signal,resource,time
import numpy as np
from protection.dynamics.service_search import free_command
from protection.dynamics.eclipse_parallel import load
from protection.dynamics.attitude_load import rigid_square_load
from protection.dynamics.pattern_return import smooth_arc
from protection.dynamics.optical import length
from .joint_common import spline,state_at
from .natural_common import *


def main():
    signal.alarm(150);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();cpu=time.process_time()
    arrays,producer=setup(['protection/dynamics/eclipse_parallel.py','protection/dynamics/eclipse_parallel.cpp','research/studies/solar_shield_array/natural_account.py'])
    env=environment(3.);ratio=arrays['mass_ratio'];mass=5e6*ratio;capacity=(ratio-1.2)*5e6*300;specific=30000/(1.4*np.cos(np.pi/4))
    out=dict(schema='terluna.research.natural-account/1',producer=producer,completed=False,cases=[],accepted_return=False,accepted_cycle=False,accepted_fleet=False,
        actual_generation_W=None,actual_bus_delivery_W=None,actual_storage_state_J=None,full_cycle_energy_J=None,
        scope='Complete modeled accounts of the executed prefixes. No future coast, next service, generator, resized actuator or fleet demand is credited. Propellant depletion is not propagated.')
    for name in ['natural_face_refined_coarse.json','natural_free_slow_fine.json']:
        p=HERE/'results'/name;parent=json.loads(p.read_text());producer['inputs'][name]=digest(p)
        tilt=0. if parent['case']=='face_refined' else 90.;cmd=free_command(env,parent['arrival_s']+21600,tilt=tilt,arrival=parent['arrival_s'],turn_start=6.,turn_hours=6.)
        stages=[];whole=0.;first=0.;redirect=0.;peak=np.zeros(361);prefix12=0.;dv=np.zeros(361);att=dv.copy();max_torque=0.
        for entry in parent['runs']:
            raw=np.load(ROOT/entry['raw_path']);assert digest(ROOT/entry['raw_path'])==entry['raw_sha256'];path=spline(raw);controls=raw['controls'];windows=raw['windows'];evaluations=[]
            for suns,step in [(8,120.),(16,60.)]:
                ts=np.unique(np.r_[np.arange(entry['start_s'],entry['end_s'],step),entry['end_s']]);powers=[];optical=[];attitudes=[];translations=[];torques=[]
                nominal=rigid_square_load(cmd,ts)['torque_per_mass']
                for k,t in enumerate(ts):
                    y=state_at(path,t);light=load(env,t,y,cmd(t).as_matrix(),suns=suns,torque=True,grid=32)
                    torque=nominal[k]-light['gravity_torque_per_mass']-light['radiation_torque_N_m']/mass[:,None]
                    couple=2*abs(torque).sum(axis=1)/10000
                    u=np.einsum('k,nki->ni',smooth_arc(t,windows),controls) if len(windows) else np.zeros((361,3))
                    powers.append((couple+length(u))*mass*specific);optical.append(light['optical']['first_intercept_bolometric_equivalent_W'])
                    attitudes.append(couple);translations.append(length(u));torques.append(float(length(torque*mass[:,None]).max()))
                powers=np.array(powers);optical=np.array(optical);attitudes=np.array(attitudes);translations=np.array(translations)
                energy=float(np.trapezoid(powers.sum(axis=1),ts));light_energy=float(np.trapezoid(optical.sum(axis=1),ts))
                item=dict(suns=suns,step_s=step,electrical_energy_J=energy,first_intercept_J=light_energy,redirected_band_J=light_energy*env.central_fraction)
                if suns==16:
                    whole+=energy;first+=light_energy;redirect+=light_energy*env.central_fraction;peak=np.maximum(peak,powers.max(axis=0));max_torque=max(max_torque,max(torques))
                    dv+=np.trapezoid(translations,ts,axis=0);att+=np.trapezoid(attitudes,ts,axis=0)
                    if ts[0]<43200:
                        stop=min(ts[-1],43200);tt=np.r_[ts[ts<stop],stop];pp=np.interp(tt,ts,powers.sum(axis=1));prefix12+=float(np.trapezoid(pp,tt))
                    item.update(save('account_'+parent['case']+'_'+entry['stage'],t=ts,power_W=powers,first_intercept_W=optical,redirected_W=optical*env.central_fraction,attitude_equivalent_acceleration=attitudes,translation_acceleration=translations))
                evaluations.append(item)
            stages.append(dict(stage=entry['stage'],start_s=entry['start_s'],end_s=entry['end_s'],evaluations=evaluations,
                relative_energy_difference=abs(evaluations[0]['electrical_energy_J']-energy)/energy,relative_optical_difference=abs(evaluations[0]['first_intercept_J']-light_energy)/light_energy))
        duration=parent['runs'][-1]['end_s'];row=dict(case=parent['case'],duration_s=duration,stages=stages,executed_electrical_energy_J=whole,executed_mean_power_W=whole/duration,
            first_intercept_J=first,redirected_band_J=redirect,translation_mean_m_s=float(dv.mean()),attitude_mean_m_s=float(att.mean()),ideal_exhaust_kg=whole*1.4/30000**2,
            propulsion_waste_heat_J=.3*whole,maximum_torque_N_m=max_torque,maximum_individual_power_W=float(peak.max()),minimum_installed_to_peak=float(np.min(capacity/peak)),overloaded_members=int(np.count_nonzero(peak>capacity*(1+1e-9))),
            first_twelve_hours_completed=duration>=43200,first_twelve_hours_energy_J=prefix12 if duration>=43200 else None,per_member_peak_W=peak.tolist(),per_member_installed_W=capacity.tolist(),
            complete_return=False,actual_supply_validated=False)
        out['cases'].append(row);write('natural_account.json',out,before,cpu)
        print(parent['case'],whole/1e12,row['minimum_installed_to_peak'],time.monotonic()-before,flush=True)
    out['completed']=True;write('natural_account.json',out,before,cpu)
if __name__=='__main__':main()
