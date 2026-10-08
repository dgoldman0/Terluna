"""Four-group early steering proposal and exact oriented-geometry checks.

One 12-hour off-service screen, starting from the actual service terminal.
Normals vary by original depth label; electric controls are the retained early
endpoint proposal. Unshadowed dynamics generate the proposal. Exact shadow
forces are measured on selected states; this family receives no acceptance.
"""
import time,signal,resource
import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation,RotationSpline
from scipy.spatial import cKDTree
from protection.dynamics.natural_pattern import finite_gravity,retarded_sun
from protection.dynamics.oriented_tiles import photon_load,square_distance
from protection.dynamics.attitude_load import rigid_square_load
from protection.dynamics.fleet import sun_points,beam_margins
from protection.dynamics.cycling import ray_geometry
from protection.dynamics.collection import require_uneclipsed_tiles
from protection.dynamics.optical import length,unit
from .joint_common import *
from shared import constants as K


def main():
    signal.alarm(160);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();cpu=time.process_time()
    p,producer=setup(['protection/dynamics/oriented_tiles.py','research/studies/solar_shield_array/joint_grouped.py'])
    env=environment(3.);ratio=np.array(p['closure_sizing.json']['per_member_total_to_optical_mass_ratio'])
    parent=json.loads((HERE/'results/joint_screen.json').read_text());case=parent['cases'][0];raw=load_raw(case)
    producer['inputs']['joint_screen.json']=digest(HERE/'results/joint_screen.json')
    initial=load_raw(p['closure_validation.json']['runs'][0])['state'][-1]
    labels=((np.arange(361)//19)%2*2+(np.arange(361)%19)%2)
    baseline=command(env,case['arrival_s']+21600,case['arrival_s'])
    times=np.arange(21600.,64800.+1,120.);x=(times-21600)/(64800-21600)
    pulse=np.sin(np.pi*x)**4
    commands=[RotationSpline(times,baseline(times)*Rotation.from_rotvec(np.c_[np.deg2rad(angle)*pulse,np.zeros((len(times),2))])) for angle in [-10.,-10/3,10/3,10.]]
    def frames(t):return np.array([c(t).as_matrix() for c in commands])[labels]
    def free_force(t,y):
        f=frames(t);q=y[:,:3];sample=env.at(t);rays=ray_geometry(q,sample);require_uneclipsed_tiles(rays)
        coefficient=np.zeros(361)
        for source in sun_points(retarded_sun(sample),8,.317):
            ray=source-q;c=np.sum(-unit(ray)*f[:,:,2],axis=1)
            coefficient+=c*abs(c)*(K.AU/length(ray))**2/8
        sail=2*env.central_fraction*K.SOLAR_CONSTANT/(K.SPEED_OF_LIGHT*env.sigma)*(9890/10000)**2*coefficient[:,None]*f[:,:,2]/ratio[:,None]
        return finite_gravity(env,t,q,f[:,:,0],f[:,:,1])+sail,sail
    def thrust(t):return np.einsum('k,nki->ni',smooth_arc(t,raw['windows']),raw['controls'])
    def rhs(t,y):
        s=y.reshape(-1,6);return np.c_[s[:,3:],free_force(t,s)[0]+thrust(t)].ravel()
    sol=solve_ivp(rhs,[21600.,64800.],initial.ravel(),method='DOP853',max_step=240.,rtol=2e-11,
        atol=np.tile([1e-5]*3+[1e-9]*3,361),dense_output=True)
    dates=np.arange(21600.,64800.+1,600.);ys=sol.sol(dates).T.reshape(-1,361,6)
    events=[];minimum=np.inf;witness=None;force_checks=[]
    for t,y in zip(dates,ys):
        f=frames(t);pairs=cKDTree(y[:,:3]).query_pairs(10000*np.sqrt(2)+500,output_type='ndarray')
        local=np.inf;which=None
        for i,j in pairs:
            d=square_distance(y[i,:3],f[i],y[j,:3],f[j])
            if d<local:local=d;which=[int(i),int(j)]
            if d<1e-7:break
        if local<minimum:minimum=local;witness=dict(time_s=float(t),pair=which)
        if local<100:events.append(dict(time_s=float(t),distance_m=float(local),pair=which))
    # Do not propagate the rejected proposal into a claimed coupled encounter.
    # These snapshots quantify shadow-force error before the first screen hit.
    first_bad=events[0]['time_s'] if events else 64800.
    for t in np.linspace(21600,min(first_bad,43200.),3):
        y=sol.sol(t).reshape(-1,6);f=frames(t)
        exact=photon_load(env,t,y,f,suns=4);_,unshadowed=free_force(t,y)
        actual=exact['reflected_force_N']/(5e6*ratio[:,None])
        basef=baseline(t).as_matrix();same=photon_load(env,t,y,np.broadcast_to(basef,(361,3,3)),suns=4)
        from protection.dynamics.fast_parallel import load
        parallel=load(env,t,y,basef,suns=4)
        force_checks.append(dict(time_s=float(t),maximum_shadow_acceleration_defect_m_s2=float(length(actual-unshadowed).max()),
            maximum_photon_torque_N_m=float(length(exact['radiation_torque_N_m']).max()),
            first_intercept_W=float(exact['optical']['first_intercept_bolometric_equivalent_W'].sum()),
            parallel_limit_force_error_N=float(length(same['reflected_force_N']-parallel['reflected_force_N']).max())))
    nominal=np.array([rigid_square_load(c,times)['force_per_mass'] for c in commands])
    base_nom=rigid_square_load(baseline,times)['force_per_mass']
    path=RUN/'grouped_steering_proposal.npz'
    np.savez_compressed(path,t=dates,state=ys,labels=labels,frames=np.array([frames(t) for t in dates]),mass_ratio=ratio)
    out=dict(schema='terluna.research.joint-grouped-screen/1',producer=producer,completed=True,
        accepted_return=False,accepted_cycle=False,accepted_fleet=False,
        group_angles_deg=[-10.,-10/3,10/3,10.],groups=4,screen_duration_h=12,
        minimum_sampled_surface_distance_m=float(minimum),witness=witness,
        samples_below_clearance=len(events),first_sample_below_clearance=events[0] if events else None,
        nominal_mean_attitude_impulse_m_s=float(np.mean(np.trapezoid(nominal,times,axis=1)[labels])),
        common_nominal_attitude_impulse_m_s=float(np.trapezoid(base_nom,times)),
        maximum_rate_deg_s=float(max(np.rad2deg(length(c(times,1))).max() for c in commands)),
        shadow_force_checks=force_checks,raw_path=str(path.relative_to(ROOT)),raw_sha256=digest(path),
        scope='No mutual shadows during reduced propagation; exact independent-normal shadows, force/moment and finite contact used in audit. Sampled collisions reject this specific group schedule; no independent-normal coupled trajectory accepted.')
    write('joint_grouped.json',out,before,cpu);print(out['resources'],events[:1],flush=True)
if __name__=='__main__':main()
