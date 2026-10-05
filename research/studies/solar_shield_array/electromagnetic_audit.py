"""Small static allocation, beam blockage and engineering sensitivity checks.

No new trajectories. Snapshot optimization is a restricted local search, not
an optimum certificate or a fuel-saving maneuver.
"""
import json
import resource
import signal
import time
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import least_squares
from scipy.integrate import cumulative_trapezoid
from scipy.spatial.transform import Rotation, RotationSpline
from engineering.electromagnetic import (square, mutual, coil_hardware, collector,
    link, storage, radiator, field_and_potential)
from shared import constants as K
from shared.provenance import constants_used
from .electromagnetic_run import GROUP, RUN, plain
from .cycling_search import ROOT,HERE,SCENARIO
from .fleet_run import digest


def blocked_segment(source, target, centres, frame, receiver, radius):
    """Conservative circular beam tube versus intervening physical squares.

    Constant radius = max(tx,rx); both tiles and radiators need full spectral
    transparency later. Here every physical square is opaque to a link.
    """
    d=(target-source)@frame;q=(centres-source)@frame
    if abs(d[2])<1e-12:return []
    fraction=q[:,2]/d[2]
    xy=fraction[:,None]*d[None,:2]-q[:,:2]
    # Square Minkowski expanded by a circumscribed square around beam disk.
    mask=(fraction>1e-7)&(fraction<1-1e-7)&(abs(xy).max(axis=1)<5000+radius)
    mask[receiver]=False
    return np.flatnonzero(mask).tolist()


def main():
    before=time.monotonic();signal.alarm(300)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    parent_path=HERE/'results/electromagnetic.json';p=json.loads(parent_path.read_text())
    sources=dict(p['producer']['source_hashes'])
    sources['research/studies/solar_shield_array/electromagnetic_audit.py']=digest(__file__)
    for f,h in sources.items():
        if digest(ROOT/f)!=h:raise ValueError('Source changed: '+f)
    out=dict(schema='terluna.research.electromagnetic-audit/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),inputs={'electromagnetic.json':digest(parent_path)}),
        accepted_trajectory=False,accepted_cycle=False,accepted_fleet=False,
        evidence='Executed static current allocations, conservative sampled beam blockage and analytic sensitivities; no propagation or plasma exclusion.',
        current_allocations=[],beam_geometry=[],sensitivities={})
    prop=SCENARIO['propulsion'];w_per_n=prop['exhaust_velocity_m_s']/(2*prop['efficiency']*np.cos(np.deg2rad(prop['cant_deg'])))
    rng=np.random.default_rng(20261005)
    departure_load=np.load(ROOT/p['stages'][1]['raw_output']['path'])
    peak_date=float(departure_load['t'][np.argmax(departure_load['propulsion_W'].sum(axis=1))])
    for ti,label,background in [(peak_date,'departure',0.),(32400.,'departure',0.),(46800.,'return',0.),
                                (peak_date,'departure',50e-9),(46800.,'return',50e-9)]:
        n=len(GROUP)
        path=ROOT/p['stages'][1 if label=='departure' else 2]['raw_output']['path']
        if digest(path)!=p['stages'][1 if label=='departure' else 2]['raw_output']['sha256']:raise ValueError('Changed group loads')
        z=np.load(path);idx=int(np.argmin(abs(z['t']-ti)));assert abs(z['t'][idx]-ti)<.01
        # Required force is already recorded in inertial coordinates; rotate it using parent raw attitude.
        trajectory=ROOT/('research/runs/solar_shield_array/closure/closure_fine_1.npz' if label=='departure' else 'research/runs/solar_shield_array/return_repair/repair_trial_local_refined_fine.npz')
        assert digest(trajectory)==p['raw_inputs'][str(trajectory.relative_to(ROOT))]
        raw=np.load(trajectory)
        frame=RotationSpline(raw['t'],Rotation.from_matrix(raw['frames']))(ti).as_matrix()
        q=CubicHermiteSpline(raw['t'],raw['state'][:,:,:3],raw['state'][:,:,3:])(ti)
        centres=(q[GROUP]-q[180])@frame
        required_f=z['translational_force_N'][idx]@frame
        required_t=z['attitude_torque_N_m'][idx]
        force=np.zeros((n,n,3));torque=np.zeros((n,n,3));inductance=np.eye(n)*coil_hardware(1e5)['inductance_H']
        for i in range(n):
            for j in range(n):
                if i==j:continue
                v=mutual(square(centre=centres[j]-centres[i]),square(),64)
                force[i,j]=v['force_N_A2'];torque[i,j]=v['torque_N_m_A2'];inductance[i,j]=v['mutual_H']
        def response(currents):
            pair=currents[:,None]*currents[None,:]
            external=np.cross(np.c_[np.zeros((n,2)),currents*1e8],np.array([0.,background,0.]))
            return np.sum(pair[:,:,None]*force,axis=1),np.sum(pair[:,:,None]*torque,axis=1)+external
        def engine_cost(f,t):return (np.linalg.norm(f,axis=1)+2*abs(t).sum(axis=1)/10000).sum()*w_per_n
        base=engine_cost(required_f,required_t)
        for cap in [1e5,3e5]:
            def residual(u):
                f,t=response(u*cap)
                return np.r_[(f-required_f).ravel(),(2*(t-required_t)/10000).ravel()]/1000
            # Include the always-available zero-current choice; least-squares
            # force error and the actual L1 torque/propulsion cost differ.
            zero=dict(success=True,status=0,nfev=0,residual_engine_W=float(base),currents_A=np.zeros(n),label='zero_current')
            records=[zero];best=zero
            for start in [np.ones(n)*.05,*rng.uniform(-.9,.9,(7,n))]:
                fit=least_squares(residual,start,bounds=(-np.ones(n),np.ones(n)),max_nfev=200,
                                  ftol=1e-9,xtol=1e-9,gtol=1e-9)
                amps=fit.x*cap;f,t=response(amps)
                engine=engine_cost(required_f-f,required_t-t)
                row=dict(success=bool(fit.success),status=int(fit.status),nfev=int(fit.nfev),residual_engine_W=float(engine),currents_A=amps)
                records.append(row)
                if best is None or engine<best['residual_engine_W']:best=row
            f,t=response(best['currents_A']);cold=9*coil_hardware(cap)['refrigerator_input_W']
            out['current_allocations'].append(dict(time_s=ti,current_cap_A=cap,
                background_body_field_T=[0.,background,0.],
                reference_engine_W=float(base),reference_required_force_N=required_f,reference_required_torque_N_m=required_t,
                best=best,all_starts=records,cryo_input_W=cold,
                residual_engine_plus_cryo_W=best['residual_engine_W']+cold,
                group_coil_mass_kg=9*coil_hardware(cap)['installed_mass_kg'],
                coil_supply_mass_kg=collector(cold)['installed_mass_kg'],
                net_internal_force_N=f.sum(axis=0),
                external_magnetic_torque_N_m=np.cross(np.c_[np.zeros((n,2)),best['currents_A']*1e8],np.array([0.,background,0.])),
                internal_total_torque_about_origin_N_m=(t-np.cross(np.c_[np.zeros((n,2)),best['currents_A']*1e8],np.array([0.,background,0.]))+np.cross(centres,f)).sum(axis=0),
                required_external_torque_about_origin_N_m=(required_t+np.cross(centres,required_f)).sum(axis=0),
                inductance_min_eigenvalue_H=float(np.linalg.eigvalsh(inductance).min()),
                stored_coupled_energy_J=float(.5*best['currents_A']@inductance@best['currents_A']),
                scope='Eight local least-squares starts plus zero-current choice; one parallel square winding per tile, other 352 coils off. Optional uniform 50 nT body-y field is an orientation sensitivity, not a lunar field trace. It exchanges angular momentum with its external source; its small spatial-gradient force is omitted. 54 force/torque components, nine current controls. No switching, acquisition, flexible modes, changed mass or motion.'))
    # Seven measured geometries; fixed offset is a siting probe, not an orbit.
    stages=['closure/closure_fine_0.npz','closure/closure_fine_1.npz','return_repair/repair_trial_local_refined_fine.npz']
    peak=np.array(p['reference']['all_stage_peak_W'])[GROUP]
    for case in p['force_samples']:
        ti=case['time_s'];rawpath=ROOT/'research/runs/solar_shield_array'/stages[0 if ti<=21600 else 1 if ti<=43200 else 2]
        assert digest(rawpath)==p['raw_inputs'][str(rawpath.relative_to(ROOT))]
        raw=np.load(rawpath);q=CubicHermiteSpline(raw['t'],raw['state'][:,:,:3],raw['state'][:,:,3:])(ti)
        frame=raw['frames'][np.argmin(abs(raw['t']-ti))];centre=q[GROUP].mean(axis=0)
        for offset in [(0,0,1e5),(0,0,-1e5),(1e5,0,1e5),(0,1e5,1e5),(0,0,3e5)]:
            source=centre+np.array(offset)@frame.T;rows=[]
            for k,i in enumerate(GROUP):
                l=link(peak[k],float(np.linalg.norm(q[i]-source)),'microwave')
                blockers=blocked_segment(source,q[i],q,frame,int(i),max(l['receiver_diameter_m']/2,50.))
                rows.append(dict(tile=int(i),range_m=l['range_m'],blockers=blockers,
                    transmission_if_each_crossing_99percent=.99**len(blockers),
                    transmission_if_each_crossing_999percent=.999**len(blockers)))
            out['beam_geometry'].append(dict(time_s=ti,hub_offset_body_m=offset,
                clear_links=sum(not r['blockers'] for r in rows),links=rows))
    spectrum_path=ROOT/'protection/spectra/shield_transmission.json'
    spectrum=json.loads(spectrum_path.read_text())
    out['producer']['optical_product']=dict(path=str(spectrum_path.relative_to(ROOT)),sha256=digest(spectrum_path))
    # Coherent thin, lossless vacuum/silica/vacuum slab plus a first-order
    # absorption sensitivity using NIST's upper measured loss-tangent bound.
    frequency=5.8e9;thickness=1e-5;index=np.sqrt(3.87)
    interface=((index-1)/(index+1))**2
    phase=2*np.pi*index*thickness*frequency/K.SPEED_OF_LIGHT
    lossless=1/(1+4*interface/(1-interface)**2*np.sin(phase)**2)
    absorption=2*np.pi*index*.005*thickness*frequency/K.SPEED_OF_LIGHT
    out['film_crossings']=dict(laser_normal_incidence_T1064=float(np.interp(1064.,spectrum['wavelength_nm'],spectrum['transmission']['titania_stack'])),
        microwave_bare_silica_thickness_m=thickness,microwave_index=index,
        microwave_lossless_slab_T=float(lossless),microwave_first_order_absorption=float(absorption),
        scope='Crossings are potential obstructions, not proven link failures. Bare 10 micrometre silica can transmit microwaves; full coated assembly, conductors, supports, incidence, phase and heating require characterization. Laser T is stored normal-incidence spectral design only; oblique or support crossings need their own model.',
        primary_source='https://www.nist.gov/publications/measuring-permittivity-fused-silica-planar-wafer-structures-325-ghz')
    # Storage chronological usable-capacity calculation with ideal common bus.
    stage=p['stages'][-1];z=np.load(ROOT/stage['raw_output']['path'])
    t=z['t'];net=z['total_propulsion_W']-z['total_conditional_generation_W']
    chemistry=np.maximum(net,0)/.95+np.minimum(net,0)*.95
    integral=np.r_[0.,cumulative_trapezoid(chemistry,t)]
    drawdown=integral-np.minimum.accumulate(integral)
    required=float(drawdown.max())
    out['chronological_storage']=dict(usable_chemical_energy_J=required,
        nominal_capacity_at_80_percent_depth_J=required/.8,
        battery_mass_at_200Wh_kg=required/(.8*200*3600),
        peak_discharge_W=float(net.max()),
        scope='Ideal pooled bus with free routing; start fully charged, spill energy above capacity, 95% each-way efficiency. Actual links add losses and availability constraints.')
    c=p['coil_reference'];out['sensitivities']['coil_cooling']=[dict(heat_leak_W_m2=h,cop=cop,
        **coil_hardware(1e5,heat_leak=h,cop=cop)) for h,cop in [(.01,.01),(.05,.005),(.1,.002)]]
    out['sensitivities']['radiators']=[dict(temperature_K=temp,**radiator(100e6,temperature_K=temp)) for temp in [300.,350.,400.]]
    out['sensitivities']['collector_specific_power']=[dict(specific_power_W_kg=power,
        mass_kg=sum(collector(x['bus_W'],specific_power=power)['installed_mass_kg'] for x in p['local_generation']['per_tile'])) for power in [100.,300.,1000.]]
    out['sensitivities']['group_size']=[dict(tiles=n,identical_100kA_coils_mass_kg=n*c['installed_mass_kg'],
        cryo_W=n*c['refrigerator_input_W'],pair_count=n*(n-1)//2,hardware_scope='Coil accounting only; no independent control or overlap-free protection implied.') for n in [2,9,25,361]]
    out['force_budget_limits']=dict(all361_net_force_impulse_fraction=stage['pooled_translational_force_lower_bound_impulse_N_s']/stage['sum_translation_impulse_N_s'],
        group_net_force_impulse_fraction=stage['group_external_force_lower_bound_impulse_N_s']/stage['group_sum_translation_impulse_N_s'],
        scope='Triangle-inequality external force lower bounds for those existing controls, ignoring angular momentum and controllability. Complement is not achievable EM savings.')
    # Explicit regional scale comparison with existing lunar magnets.
    distance=78000e3;moment=4*3.75e20
    vacuum_axis=K.VACUUM_PERMEABILITY*moment/(2*np.pi*distance**3)
    out['existing_architecture']=dict(four_anchor_total_moment_A_m2=moment,
        optical_reference_distance_m=distance,axial_upper_field_T=vacuum_axis,
        max_torque_on_100kA_tile_N_m=1e13*vacuum_axis,
        scale_gradient_force_N=3*1e13*vacuum_axis/distance,
        scope='Far-field maximum for aligned lunar anchors at the old 78,000 km optical distance, not actual time-varying field along this orbit. Lunar anchors do not supply arbitrary local relative forces.')
    actual=[]
    for name in stages:
        raw=np.load(ROOT/'research/runs/solar_shield_array'/name)
        r=np.linalg.norm(raw['state'][:,:,:3],axis=2)
        actual.append((float(r.min()),float(r.max())))
    low=min(x[0] for x in actual);high=max(x[1] for x in actual)
    equatorial=K.VACUUM_PERMEABILITY*moment/(4*np.pi*high**3)
    axial=K.VACUUM_PERMEABILITY*moment/(2*np.pi*low**3)
    out['existing_architecture'].update(actual_moon_range_m=[low,high],
        central_dipole_orientation_field_envelope_T=[equatorial,axial],
        maximum_100kA_tile_torque_scale_N_m=1e13*axial,
        gradient_force_scale_N=3*1e13*axial/low,
        actual_scope='Measured 14.9–15.2 thousand km range, central moment orientation envelope. Four displaced anchors, lunar attitude and plasma distortions need the real field trace. This background may matter for torque even when its translational gradient force is small.')
    # BL through a path across a coil opening. No trajectory deflection or dose inference.
    from scipy.integrate import quad_vec
    paths=[]
    for side,amps in [(10000.,1e5),(100000.,1e6)]:
        value,error=quad_vec(lambda x:field_and_potential([[x,0,.1*side]],square(side))[0][0]*amps,-20*side,20*side,epsabs=1e-10)
        absolute,_=quad_vec(lambda x:np.linalg.norm(field_and_potential([[x,0,.1*side]],square(side))[0][0]*amps),-20*side,20*side,epsabs=1e-10)
        paths.append(dict(side_m=side,current_A=amps,path_height_m=.1*side,integrated_B_vector_T_m=value,
            integrated_absolute_B_T_m=float(absolute),quadrature_error_T_m=float(error),
            scope='Straight transverse path across opening, finite +/-20-side integration. Field reversals matter. Particle path changes, cusps and plasma not solved.'))
    out['field_integrals']=paths
    out['resources']=dict(elapsed_s=time.monotonic()-before,max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
        total_numerical_charge_s=10+p['resources']['elapsed_s']+time.monotonic()-before,threads=1,address_space_GiB=2)
    out['completed']=True
    (HERE/'results/electromagnetic_audit.json').write_text(json.dumps(out,default=plain,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(resources=out['resources'],allocations=[{k:r[k] for k in ['time_s','current_cap_A','reference_engine_W','residual_engine_plus_cryo_W']} for r in out['current_allocations']],storage=out['chronological_storage']),default=plain))


if __name__=='__main__':main()
