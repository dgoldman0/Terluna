"""Recurring global finite-tile schedule: geometry, return safety and cost."""
import argparse
from dataclasses import asdict
import json
import time

import numpy as np
from scipy.stats import qmc

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.active_global import Traffic, Layout, RayIndex, sun_source, pole_collisions
from protection.dynamics.optical import length
from .fleet_run import environment, digest, identities
from .cycling_search import ROOT, HERE, SCENARIO
from .cycling_analysis import compact_series_json


def sources():
    paths=['protection/dynamics/active_global.py','protection/dynamics/active_formation.py',
           'research/studies/solar_shield_array/active_global_run.py',
           'research/studies/solar_shield_array/cycling_analysis.py']
    return {**identities(),**{p:digest(ROOT/p) for p in paths}}


def disk_points(power,seed=483):
    u=qmc.Sobol(2,scramble=True,seed=seed).random_base2(power)
    r=np.sqrt(u[:,0]);a=2*np.pi*u[:,1]
    return np.column_stack([r*np.cos(a),r*np.sin(a)])


def coverage(fleet,times,receiver_power=13,sun_power=4):
    receiver=disk_points(receiver_power)*fleet.layout.protected_radius
    # Include the exact outer boundary separately; area points cannot test it.
    angle=np.arange(1024)*2*np.pi/1024
    rim=fleet.layout.protected_radius*np.column_stack([np.cos(angle),np.sin(angle)])
    xy=np.vstack([receiver,rim]); suns=disk_points(sun_power,117)
    rows=[]; witnesses=[]
    for ti in times:
        hits=[]
        for si,sun in enumerate(suns):
            src=sun_source(fleet.env.at(ti),fleet.frames(ti),sun)
            found=RayIndex(fleet,ti,src).intercept(xy)
            hits.append(found)
            for k in np.flatnonzero(~found)[:16]:
                witnesses.append(dict(time_s=float(ti),sun_xy=sun.tolist(),
                    receiver_xy_m=xy[k].tolist(),rim=bool(k>=len(receiver))))
        h=np.array(hits)
        rows.append(dict(day=float(ti/K.JULIAN_DAY),ray_coverage=float(h[:,:len(receiver)].mean()),
            all_sampled_sun_receiver_fraction=float(np.all(h[:,:len(receiver)],axis=0).mean()),
            rim_ray_coverage=float(h[:,len(receiver):].mean()),
            missed_rays=int(np.count_nonzero(~h))))
        print(json.dumps(dict(stage='coverage',**rows[-1])),flush=True)
    return dict(dates=len(times),solar_disk_directions=len(suns),area_receivers=len(receiver),
        rim_receivers=len(rim),total_ray_tests=len(times)*len(suns)*len(xy),
        mean_ray_coverage=float(np.trapezoid([r['ray_coverage'] for r in rows],times)/(times[-1]-times[0])) if len(times)>1 else rows[0]['ray_coverage'],
        minimum_ray_coverage=min(r['ray_coverage'] for r in rows),
        minimum_rim_coverage=min(r['rim_ray_coverage'] for r in rows),
        all_time_all_source_coverage_certified=False,
        scope='Direct intersections with individual identified finite squares, at actual common ephemeris dates. Receiver-plane geometry includes finite Sun and first-order moving-tile light time; sampled geometry is not a continuous coverage theorem.',
        series=rows,miss_witnesses=witnesses)


def budget(fleet,times,plane_stride=16,phase_count=128):
    n=fleet.layout.planes
    if n%plane_stride:raise ValueError('Plane stride must divide the population')
    # Midpoint representatives of equal-size plane blocks, both radial parities.
    planes=np.arange(n//plane_stride)*plane_stride+plane_stride//2
    phases=(np.arange(phase_count)+.371)*2*np.pi/phase_count
    j,par,ph=np.meshgrid(planes,[0,1],phases,indexing='ij')
    shape=j.shape;j=j.ravel();par=par.ravel();ph=ph.ravel()
    weight=fleet.count[j]/phase_count/2*plane_stride
    weight*=fleet.inventory/weight.sum()
    series=[];per_mean=np.zeros((2,*shape[:2]));per_peak=np.zeros(shape[:2])
    per_mean_b=np.zeros(shape[:2]);minimum_margin=np.full(2,np.inf);peak_normal=0.
    # Trapezoid weights on nonuniform dates, explicitly numeric.
    tw=np.zeros(len(times));tw[:-1]+=np.diff(times)/2;tw[1:]+=np.diff(times)/2;tw/=times[-1]-times[0]
    for i,ti in enumerate(times):
        lo,hi,d=fleet.power_bounds(ti,j,ph,par)
        minimum_margin=np.minimum(minimum_margin,d['margins'].min(axis=0));peak_normal=max(peak_normal,float(d['tilt'].max()))
        per_mean[0]+=tw[i]*lo.reshape(shape).mean(axis=2)
        per_mean[1]+=tw[i]*hi.reshape(shape).mean(axis=2)
        per_peak=np.maximum(per_peak,hi.reshape(shape).max(axis=2))
        series.append(dict(day=float(ti/K.JULIAN_DAY),
            mean_acceleration_lower_m_s2=float(np.dot(weight,lo)/weight.sum()),
            mean_acceleration_upper_m_s2=float(np.dot(weight,hi)/weight.sum()),
            peak_acceleration_upper_m_s2=float(hi.max()),
            fraction_requiring_over_1um_s2=float(weight[lo>1e-6].sum()/weight.sum())))
        if i%12==0:print(json.dumps(dict(stage='budget',**series[-1])),flush=True)
    means=np.array([[r['mean_acceleration_lower_m_s2'],r['mean_acceleration_upper_m_s2']] for r in series])
    mean=tw@means
    count_weight=fleet.count[planes,None]/2*plane_stride*np.ones((1,2))
    count_weight*=fleet.inventory/count_weight.sum()
    mass=fleet.inventory*fleet.layout.side**2*fleet.env.sigma
    p=SCENARIO['propulsion'];cant=np.cos(np.deg2rad(p['cant_deg']))
    watts_per_newton=p['exhaust_velocity_m_s']/(2*p['efficiency']*cant)
    mean_power=mass*mean*watts_per_newton
    flow=mass*mean/(p['exhaust_velocity_m_s']*cant)
    dry_fraction=per_peak*watts_per_newton*p['peak_margin_factor']/p['specific_power_W_kg']
    fuel_fraction=per_mean[1]*p['propellant_buffer_days']*K.JULIAN_DAY/(p['exhaust_velocity_m_s']*cant)
    coefficient=dry_fraction+fuel_fraction
    if np.max(coefficient)>=1:raise ValueError('Propulsion hardware/reserve mass closure diverges')
    mass_factor=1/(1-coefficient)
    base_tile_mass=fleet.layout.side**2*fleet.env.sigma
    loaded_mass=float(np.sum(count_weight*mass_factor)*base_tile_mass)
    loaded_power_upper=float(np.sum(count_weight*mass_factor*per_mean[1])*base_tile_mass*watts_per_newton)
    per_cycle=per_mean* (2*np.pi/np.abs(fleet.omega[planes]))[None,:,None]
    return dict(plane_stride=plane_stride,phase_count_per_parity=phase_count,dates=len(times),
        members_in_force_quadrature=len(j),mean_acceleration_bounds_m_s2=mean.tolist(),
        mean_power_bounds_W=mean_power.tolist(),propellant_rate_bounds_kg_s=flow.tolist(),
        thirty_day_energy_bounds_J=(mean_power*(times[-1]-times[0])).tolist(),
        thirty_day_propellant_bounds_kg=(flow*(times[-1]-times[0])).tolist(),
        mean_tile_delta_v_bounds_m_s=(mean*(times[-1]-times[0])).tolist(),
        cycle_delta_v_bounds_m_s=[float(per_cycle[0].min()),float(per_cycle[1].max())],
        mean_peak_distributed_acceleration_m_s2=float(np.average(per_peak,weights=count_weight)),
        maximum_sampled_acceleration_m_s2=float(per_peak.max()),
        installed_distributed_power_upper_W=float(np.sum(count_weight*per_peak)*base_tile_mass*watts_per_newton*p['peak_margin_factor']),
        added_distributed_power_hardware_upper_kg=float(np.sum(count_weight*dry_fraction)*base_tile_mass),
        seven_day_propellant_reserve_upper_kg=float(np.sum(count_weight*fuel_fraction)*base_tile_mass),
        loaded_control_only_mass_upper_kg=loaded_mass,
        loaded_control_only_mean_power_upper_W=loaded_power_upper,
        loaded_control_only_propellant_rate_upper_kg_s=loaded_power_upper*2*p['efficiency']/p['exhaust_velocity_m_s']**2,
        maximum_mass_multiplier=float(mass_factor.max()),
        minimum_beam_margins_deg=np.rad2deg(minimum_margin).tolist(),
        maximum_normal_tilt_deg=float(np.rad2deg(peak_normal)),
        scope='Nonlinear force bracket over every possible mutual illumination fraction [0,1]; spatial/phase quadrature of explicit common-date paths. Not an illumination simulation or energy optimum.',
        mass_closure='Conservative per-plane/per-parity feedback closure of the assumed 300 W/kg power hardware plus seven-day propellant. Thrust per loaded mass remains within the zero-to-full photon-force convex bracket. Collector/habitat/thermal/plume hardware is not added.',
        series=series)


def census(fleet,t):
    """Sum all actual integer slots at one date; no phase quadrature."""
    sums=np.zeros(2);peak=0.;margins=np.full(2,np.inf);tilt=0.;processed=0
    for first in range(0,fleet.layout.planes,8):
        jj=[];kk=[]
        for j in range(first,min(first+8,fleet.layout.planes)):
            jj.append(np.full(fleet.count[j],j,int));kk.append(np.arange(fleet.count[j]))
        j=np.concatenate(jj);k=np.concatenate(kk)
        phase=2*np.pi*k/fleet.count[j]+fleet.phase[j]+fleet.omega[j]*t
        lo,hi,d=fleet.power_bounds(t,j,phase,k%2)
        sums+=[lo.sum(),hi.sum()];peak=max(peak,float(hi.max()));processed+=len(j)
        margins=np.minimum(margins,d['margins'].min(axis=0));tilt=max(tilt,float(d['tilt'].max()))
        if first%256==0:print(json.dumps(dict(stage='integer_census',processed=processed,total=fleet.inventory)),flush=True)
    return dict(day=t/K.JULIAN_DAY,integer_tiles_evaluated=processed,
        mean_acceleration_bounds_m_s2=(sums/processed).tolist(),
        maximum_acceleration_m_s2=peak,minimum_beam_margins_deg=np.rad2deg(margins).tolist(),
        maximum_normal_tilt_deg=float(np.rad2deg(tilt)))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--days',type=float,default=30.)
    p.add_argument('--coverage-step-hours',type=float,default=12.)
    p.add_argument('--receiver-power',type=int,default=13)
    p.add_argument('--sun-power',type=int,default=4)
    p.add_argument('--census',action='store_true')
    p.add_argument('--output',default='active_global.json')
    a=p.parse_args();before=sources();start=time.monotonic();env=environment(a.days)
    fleet=Traffic(env,a.days);compact=Traffic(env,a.days,Layout(radial_step=1000.))
    times=np.linspace(0,a.days*K.JULIAN_DAY,int(a.days*24/a.coverage_step_hours)+1)
    optical=coverage(fleet,times,a.receiver_power,a.sun_power)
    budget_times=np.linspace(0,a.days*K.JULIAN_DAY,int(a.days*4)+1)
    costs=budget(fleet,budget_times)
    rejected=dict(layout=asdict(compact.layout),inventory=compact.inventory,
        separation=compact.separation_certificate(),
        intersections=[pole_collisions(compact,t) for t in (0,7.5*K.JULIAN_DAY,15*K.JULIAN_DAY)],
        scope='Compact comparison is rejected for actual finite-square intersections; it is not used as an operational coverage baseline')
    periods=2*np.pi/np.abs(fleet.omega)
    duty=np.arcsin(fleet.layout.protected_radius/fleet.r0)/np.pi
    result=dict(schema='terluna.research.active-global-traffic/1',
        producer=dict(source_hashes=before,constants=constants_used(before)),
        parameters=vars(a),layout=asdict(fleet.layout),
        inventory=fleet.inventory,optical_base_mass_kg=fleet.inventory*fleet.layout.side**2*env.sigma,
        radius_range_km=[float(fleet.r0.min()/1000),float((fleet.r0.max()+fleet.layout.radial_step)/1000)],
        cycle_period_days=[float(periods.min()/K.JULIAN_DAY),float(periods.max()/K.JULIAN_DAY)],
        phase_slots_per_plane=[int(fleet.count.min()),int(fleet.count.max())],
        front_passages_per_day=float(np.sum(fleet.count/periods)*K.JULIAN_DAY),
        mean_central_source_service_fraction=float(np.average(duty,weights=fleet.count)),
        nominal_handover_interval_s=[float(np.min(periods/fleet.count)),float(np.max(periods/fleet.count))],
        coverage=optical,propulsion_assumptions=SCENARIO['propulsion'],budget=costs,
        separation=fleet.separation_certificate(),compact_rejected=rejected,
        integer_census=census(fleet,15*K.JULIAN_DAY) if a.census else None,
        evidence='Prescribed simultaneous finite-square traffic with actual-date inverse dynamics; complete return paths and a global separation bound. Not a freely propagated fleet or bounded-slew attitude/control demonstration.',
        acceptance=dict(global_spatial_sampling_complete=True,return_paths_defined=True,
            global_four_lunar_radius_coverage_demonstrated=False,
            continuous_finite_sun_coverage_certified=False,all_member_ode_replay=False,
            bounded_slew_attitudes_demonstrated=False,energy_optimum=False,
            delivered_electrical_power_demonstrated_W=None),
        limitations=['Coverage samples do not exclude arbitrarily small holes between receiver/source/date samples',
            'Earth pointing has no bounded-slew realization; ideal angles, flexible membranes, thruster and plume geometry remain unimplemented',
            'Mutual shadow changes are bracketed in translational force, not simulated; multiple reflected encounters are omitted',
            'Power brackets describe this prescribed family, not a lower bound for all orbital choreography',
            'The full census is one date; month-long costs are spatial/temporal quadrature, with separate convergence checks',
            'Collector/habitat loads and delivered power remain separate'],elapsed_s=time.monotonic()-start)
    if before!=sources():raise ValueError('Source changed during run')
    (HERE/'results'/a.output).write_text(compact_series_json(result))
    print(json.dumps(dict(stage='saved',output=a.output,elapsed_s=result['elapsed_s'])),flush=True)


if __name__=='__main__':main()
