"""Per-tile optical and actuator account for the cheaper stopped execution."""
import time,resource,signal
import numpy as np
from scipy.interpolate import interp1d
from protection.dynamics.eclipse_parallel import load
from protection.dynamics.attitude_load import rigid_square_load
from protection.dynamics.pattern_return import smooth_arc
from .joint_common import *
from .joint_bridge import EXTRA


def main():
    signal.alarm(40);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();cpu=time.process_time();p,producer=setup([*EXTRA,'research/studies/solar_shield_array/joint_failed_account.py'])
    name='joint_early_same_coarse.json';parent=json.loads((HERE/'results'/name).read_text());producer['inputs'][name]=digest(HERE/'results'/name)
    screen=json.loads((HERE/'results/joint_screen.json').read_text());proposal=load_raw(screen['cases'][0]);producer['inputs']['joint_screen.json']=digest(HERE/'results/joint_screen.json')
    env=environment(3);ratio=proposal['mass_ratio'];mass=5e6*ratio;capacity=(ratio-1.2)*5e6*300;specific=30000/(2*.7*np.cos(np.pi/4))
    cmd=command(env,parent['arrival_s']+21600,parent['arrival_s']);stages=[];all_peaks=[]
    out=dict(schema='terluna.research.joint-failed-account/1',producer=producer,accepted_cycle=False,accepted_return=False,
        actual_generation_W=None,actual_bus_supply_W=None,full_cycle_energy_J=None,stages=stages,
        scope='Complete executed prefix only; no fleet extrapolation or independent tighter replay of this rejected candidate.')
    for entry in parent['runs']:
        raw=load_raw(entry);s=spline(raw);ts=np.unique(np.r_[np.arange(entry['start_s'],entry['end_s'],120),entry['end_s']])
        light=[];dist=[]
        for t in ts:
            a=load(env,t,state_at(s,t),cmd(t).as_matrix(),suns=16,torque=True)
            light.append(a['optical']['first_intercept_bolometric_equivalent_W'])
            dist.append(a['gravity_torque_per_mass']+a['radiation_torque_N_m']/mass[:,None])
        dense=np.unique(np.r_[np.arange(ts[0],ts[-1],10.),ts[-1]])
        torque=rigid_square_load(cmd,dense)['torque_per_mass'][:,None,:]-interp1d(ts,dist,axis=0)(dense)
        attitude=2*abs(torque).sum(axis=-1)/10000
        u=np.array([np.einsum('k,nki->ni',smooth_arc(t,proposal['windows']),proposal['controls']) if t>=21600 else np.zeros((361,3)) for t in dense])
        power=(attitude+length(u))*mass*specific;peak=power.max(axis=0);all_peaks.append(peak)
        light=np.array(light);path=RUN/('failed_account_'+entry['stage']+'.npz')
        np.savez_compressed(path,t=dense,power_W=power,optical_times=ts,first_intercept_W=light,redirected_band_W=light*env.central_fraction)
        energy=float(np.trapezoid(power.sum(axis=1),dense));optical=float(np.trapezoid(light.sum(axis=1),ts))
        stages.append(dict(stage=entry['stage'],start_s=entry['start_s'],end_s=entry['end_s'],energy_J=energy,
            first_intercept_J=optical,redirected_band_J=optical*env.central_fraction,
            mean_translation_m_s=float(np.trapezoid(length(u),dense,axis=0).mean()),
            mean_attitude_m_s=float(np.trapezoid(attitude,dense,axis=0).mean()),
            raw_path=str(path.relative_to(ROOT)),raw_sha256=digest(path)))
    peak=np.maximum.reduce(all_peaks)
    out.update(completed=True,executed_energy_J=sum(s['energy_J'] for s in stages),
        executed_first_intercept_J=sum(s['first_intercept_J'] for s in stages),
        overloaded_members=int(np.count_nonzero(peak>capacity)),minimum_individual_capacity_to_peak=float(np.min(capacity/peak)),
        per_member_peak_W=peak.tolist(),ended_at_clearance_floor=True)
    write('joint_failed_account.json',out,before,cpu)
if __name__=='__main__':main()
