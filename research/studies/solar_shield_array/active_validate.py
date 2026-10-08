"""Read-back validation and conditional active-control budget scales."""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.active_formation import Formation,solar_frame,coverage_window
from protection.dynamics.fleet import square_vertices
from protection.dynamics.active_retime import UV,WEIGHTS
from .fleet_run import environment,identities,digest
from .cycling_search import HERE,ROOT,SCENARIO
from .cycling_analysis import compact_series_json


def checked_raw(case):
    path=ROOT/case['raw_path']
    if digest(path)!=case['raw_sha256']:raise ValueError('Raw trace hash mismatch')
    return np.load(path)


def main():
    filenames=['active_formation.json','active_formation_replay.json','active_acquisition.json','active_retime.json','active_gap_scale.json']
    products={n:json.loads((HERE/'results'/n).read_text()) for n in filenames}
    for p in products.values():
        for n,h in p['producer']['source_hashes'].items():
            if digest(ROOT/n)!=h:raise ValueError('Product source changed: '+n)
    base=products[filenames[0]]['cases']['0.0005']; replay=products[filenames[1]]['cases']['0.0005']
    raw=checked_raw(base); repeated=checked_raw(replay)
    if not np.array_equal(raw['t'],repeated['t']):raise ValueError('Replay grids differ')
    difference=np.linalg.norm(raw['state'][...,:3]-repeated['state'][...,:3],axis=-1)
    acq=products['active_acquisition.json']; acquisition=checked_raw(acq)
    env=environment(1); duration=6*3600
    form=Formation(env,None,duration)
    basis=solar_frame(env.at(0));u,b=basis[:,:2].T
    if np.dot(np.cross(u,b),env.frames(0)[:,2])>0:b=-b
    radius=15e6;phase=-np.pi/8
    initial_c=np.r_[radius*(u*np.cos(phase)+b*np.sin(phase)),
        np.sqrt(K.MOON_GM/radius)*(-u*np.sin(phase)+b*np.cos(phase))]
    def cf(t,y):return np.r_[y[3:],env.gravity(t,y[:3])+form.isolated_sail(t,y[None,:])[0]]
    centre=solve_ivp(cf,[0,duration],initial_c,dense_output=True,method='DOP853',rtol=2e-12,
        atol=[1e-5]*3+[1e-9]*3,max_step=30.)
    if not centre.success:raise RuntimeError(centre.message)
    geometry=[]
    # Refine the claimed completed state and the end. The initial mean-ray
    # statistic retains its recorded 32-direction resolution. Full-Sun area
    # for the initially perforated formation is not claimed converged.
    for minute in [50,180]:
        i=int(np.argmin(abs(acquisition['t']-minute*60)));t=acquisition['t'][i]
        state=acquisition['state'][i]; sample=env.at(t);normal=-form.frames(t)[:,0]
        geometry.append(dict(minutes=float(t/60),**coverage_window(state,normal,sample,10000.,centre.sol(t),suns=128)))
        print(json.dumps(dict(stage='solar_refinement',**geometry[-1])),flush=True)
    residuals=[]
    for index in [0,len(raw['t'])//2,len(raw['t'])-1]:
        t=raw['t'][index];state=raw['state'][index];frame=form.frames(t);sample=env.at(t)
        vertices=square_vertices(state[:,:3],np.broadcast_to(-frame[:,0],state[:,:3].shape),
            np.broadcast_to(sample['positions']['sun'],state[:,:3].shape),10000.)
        a=(vertices[:,1]-vertices[:,0])/10000;b=(vertices[:,3]-vertices[:,0])/10000
        points=state[:,None,:3]+5000*(UV[None,:,:1]*a[:,None,:]+UV[None,:,1:]*b[:,None,:])
        averaged=np.einsum('j,ijk->ik',WEIGHTS,env.gravity(t,points))
        error=np.linalg.norm(averaged-env.gravity(t,state[:,:3]),axis=-1)
        residuals.append(dict(hours=float(t/3600),maximum_m_s2=float(error.max())))
    # Conditional scaling is bookkeeping only. The reference area supplies no
    # spatial population or assignment and is never counted as achieved coverage.
    aperture_area=np.pi*(4*K.MOON_RADIUS)**2
    base_mass=aperture_area*.05
    scaling={}
    for name,c in products['active_retime.json']['cases'].items():
        if not c['replay_agreement_within_50m']:continue
        per_kg=c['mean_electric_power_over_maneuver_W']/c['optical_mass_kg']
        scaling[name]=dict(mean_power_if_one_aperture_area_of_mass_performs_this_maneuver_W=per_kg*base_mass,
            mean_power_if_entire_3072_probe_area_inventory_performs_it_W=per_kg*1.536e14,
            peak_power_hardware_fraction_of_tile_mass=c['peak_electric_power_W']*SCENARIO['propulsion']['peak_margin_factor']/SCENARIO['propulsion']['specific_power_W_kg']/c['optical_mass_kg'])
    idx=np.flatnonzero(raw['t']<=acquisition['t'][-1])
    preset_dv=float(np.trapezoid(raw['control_m_s2'][idx].mean(axis=1),raw['t'][idx]))
    sources={**identities(),'protection/dynamics/active_formation.py':digest(ROOT/'protection/dynamics/active_formation.py'),
        'protection/dynamics/active_retime.py':digest(ROOT/'protection/dynamics/active_retime.py'),
        'research/studies/solar_shield_array/active_validate.py':digest(__file__),
        'research/studies/solar_shield_array/cycling_analysis.py':digest(HERE/'cycling_analysis.py')}
    result=dict(schema='terluna.research.active-control-validation/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),
            input_products={n:digest(HERE/'results'/n) for n in filenames}),
        evidence='Independent replay comparison, finer solar coverage and finite-extent force residuals; conditional power scaling kept separate from fleet evidence',
        formation_max_replay_position_difference_m=float(difference.max()),
        formation_replay_50m_gate=bool(difference.max()<50),
        acquisition_128_sun_coverage=geometry,formation_finite_extent_gravity_residuals=residuals,
        acquisition_mean_delta_v_m_s=acq['mean_delta_v_m_s'],preset_same_interval_mean_delta_v_m_s=preset_dv,
        conditional_inventory_scales=dict(aperture_area_m2=aperture_area,aperture_base_mass_kg=base_mass,
            cases=scaling,reading_rule='Assumes each stated mass executes one identical maneuver over the stated days. This is not an achieved fleet, an operating duty cycle, available electrical supply or delivered power. No uncovered-area fraction is used to infer moved mass.'),
        global_coverage_demonstrated=False,delivered_power_W=None)
    (HERE/'results/active_validation.json').write_text(compact_series_json(result))
    print(json.dumps({k:v for k,v in result.items() if k!='producer'}),flush=True)


if __name__=='__main__':main()
