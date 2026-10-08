"""Modest independent row-band sail pulses, screened before long returns."""
import signal,resource,time
import numpy as np
from scipy.spatial.transform import Rotation,RotationSpline
from protection.dynamics.service_search import free_command
from protection.dynamics.service_reduced import reduced_path,reduced_force
from protection.dynamics.oriented_clearance import closest_oriented
from protection.dynamics.oriented_tiles import photon_load
from protection.dynamics.attitude_load import rigid_square_load
from protection.dynamics.tile_torque import gravity_torque
from protection.dynamics.optical import length
from .natural_common import *


def main():
    signal.alarm(180);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();cpu=time.process_time()
    arrays,producer=setup(['protection/dynamics/service_reduced.py','protection/dynamics/oriented_clearance.py','protection/dynamics/oriented_tiles.py','research/studies/solar_shield_array/natural_groups.py'])
    env=environment(2.);ratio=arrays['mass_ratio'];mass=5e6*ratio;capacity=(ratio-1.2)*5e6*300
    labels=(np.arange(361)//19//2)%4;times=np.arange(21600.,86400.+1,120.)
    baseline=free_command(env,86400.,turn_start=6.,turn_hours=6.)
    pulse=np.sin(np.pi*(times-21600)/(86400-21600))**4
    out=dict(schema='terluna.research.natural-groups/1',producer=producer,completed=False,cases=[],accepted_return=False,accepted_cycle=False,
        scope='Four row-band normals, smooth sail pulses and zero electric translation. Centre-force proposal dynamics; finite-square contact is exact at samples and local refinement. Exact oriented shadow/torque snapshots audit pre-contact states. No grouped coupled trajectory or next passage accepted.')
    for amplitude in [5.,15.]:
        offsets=np.linspace(-amplitude,amplitude,4)
        commands=[RotationSpline(times,baseline(times)*Rotation.from_rotvec(np.c_[np.deg2rad(a)*pulse,np.zeros((len(times),2))])) for a in offsets]
        def frames(t):return np.array([c(t).as_matrix() for c in commands])[labels]
        class Command:
            def __call__(self,t):return Rotation.from_matrix(frames(t))
        sol=reduced_path(env,arrays['first_service_end_state'],21600,86400,Command(),mass_ratio=ratio)
        dates=np.arange(21600.,86400.+1,300.);ys=sol.sol(dates).T.reshape(-1,361,6);first=None;minimum=np.inf
        for k,(t,y) in enumerate(zip(dates,ys)):
            distance,pair=closest_oriented(y[:,:3],frames(t));minimum=min(minimum,distance)
            if distance<100:
                lo=dates[k-1];hi=t
                for _ in range(13):
                    mid=(lo+hi)/2;d,_=closest_oriented(sol.sol(mid).reshape(361,6)[:,:3],frames(mid))
                    if d<100:hi=mid
                    else:lo=mid
                d,pair=closest_oriented(sol.sol(hi).reshape(361,6)[:,:3],frames(hi))
                first=dict(time_s=float(hi),distance_m=d,pair=pair,bracket_width_s=float(hi-lo));break
        finish=first['time_s'] if first else 86400.
        ts=np.unique(np.r_[np.arange(21600.,finish,300.),finish]);states=sol.sol(ts).T.reshape(-1,361,6)
        nominal=np.array([rigid_square_load(c,ts)['torque_per_mass'] for c in commands]);powers=[];force_checks=[];optical=[]
        # Keep the only expensive exact-shadow evaluations before the first failure.
        for t in np.linspace(21600.,finish,3):
            y=sol.sol(t).reshape(361,6);exact=photon_load(env,t,y,frames(t),suns=8)
            force_checks.append(dict(time_s=float(t),first_intercept_W=float(exact['optical']['first_intercept_bolometric_equivalent_W'].sum()),maximum_photon_torque_N_m=float(length(exact['radiation_torque_N_m']).max())))
        for k,t in enumerate(ts):
            torques=np.empty((361,3))
            for g,c in enumerate(commands):
                chosen=labels==g;torques[chosen]=nominal[g,k]-gravity_torque(env,t,states[k,chosen,:3],c(t).as_matrix())
            powers.append(2*abs(torques).sum(axis=1)/10000*mass*30000/(1.4*np.cos(np.pi/4)))
        powers=np.array(powers)
        row=dict(amplitude_deg=amplitude,group_offsets_deg=offsets.tolist(),first_clearance_failure=first,
            estimated_prefix_energy_J=2.6480316453385626e12+float(np.trapezoid(powers.sum(axis=1),ts)),
            estimated_overloaded_members=int(np.count_nonzero(powers.max(axis=0)>capacity)),
            maximum_rate_deg_s=float(max(np.rad2deg(length(c(ts,1))).max() for c in commands)),shadow_snapshots=force_checks,
            estimated_attitude_scope='Rigid command and finite-area gravity; shadow torque recorded at three snapshots, not integrated. Energy is not a complete executed account.',
            **save('groups_'+str(int(amplitude)),t=ts,state=states,mass_ratio=ratio,frames=np.array([frames(t) for t in ts]),labels=labels,estimated_power_W=powers))
        out['cases'].append(row);write('natural_groups.json',out,before,cpu)
        print(amplitude,first,row['estimated_prefix_energy_J']/1e12,time.monotonic()-before,flush=True)
    out['completed']=True;write('natural_groups.json',out,before,cpu)
if __name__=='__main__':main()
