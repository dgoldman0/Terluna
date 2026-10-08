"""Measure what fails when global prescribed-circle holding is removed."""
import json
import argparse
import hashlib
import inspect
import time
import numpy as np
from scipy.integrate import solve_ivp

from shared import constants as K
from shared.provenance import constants_used,constants_changed
from protection.dynamics.active_global import Traffic,sun_source,pole_collisions
from protection.dynamics.natural_traffic import NaturalTraffic,NaturalRayIndex
from protection.dynamics.optical import length
from .active_global_run import sources,disk_points
from .fleet_run import environment,digest
from .cycling_search import ROOT,HERE,SCENARIO
from .cycling_analysis import compact_series_json


def identities():
    return {**sources(),**{p:digest(ROOT/p) for p in
        ['protection/dynamics/natural_traffic.py',
         'research/studies/solar_shield_array/natural_traffic_run.py']}}


def replay(f):
    rng=np.random.default_rng(24426)
    j=np.r_[rng.integers(0,f.layout.planes,128),[0,1,15,16,17,1135,1150,1151,1152,1153,1154,1168,2287,2302,2303]]
    k=np.array([rng.integers(f.count[x]) for x in j])
    initial=f.initial_states(j,2*np.pi*k/f.count[j]+f.phase[j],k%2)
    def rhs(t,y):
        y=y.reshape(-1,6)
        return np.column_stack([y[:,3:],f.env.gravity(t,y[:,:3])]).ravel()
    # Check the entire seven-day cache, but admit optical diagnostics only on
    # the independently validated three-day prefix. Preserve the later failure.
    t=np.unique(np.r_[np.linspace(0,3*K.JULIAN_DAY,31),
        [1379.,21617.,43291.,86001.,172981.,302003.,430017.,604800.]])
    sol=solve_ivp(rhs,[0,t[-1]],initial.ravel(),method='DOP853',t_eval=t,
        max_step=300.,rtol=2e-12,atol=np.tile([1e-5]*3+[1e-9]*3,len(j)))
    if not sol.success:raise RuntimeError(sol.message)
    actual=sol.y.T.reshape(len(t),len(j),6)
    predicted=np.array([f.path(ti,j,k)[0] for ti in t])
    difference=length(predicted[:,:,:3]-actual[:,:,:3])
    prefix=t<=3*K.JULIAN_DAY
    return dict(integer_members=len(j),max_step_s=300.,rtol=2e-12,
        maximum_position_difference_m=float(difference.max()),
        maximum_velocity_difference_m_s=float(length(predicted[:,:,3:]-actual[:,:,3:]).max()),
        times_s=t.tolist(),position_max_by_date_m=difference.max(axis=1).tolist(),
        seven_day_interpolation_within_50m=bool(difference.max()<50),
        geometry_horizon_days=3.,
        geometry_prefix_maximum_position_difference_m=float(difference[prefix].max()),
        interpolation_within_50m=bool(difference[prefix].max()<50))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--reuse-raw',action='store_true');args=ap.parse_args()
    before=identities();start=time.monotonic();env=environment(30);f=NaturalTraffic(env)
    print(json.dumps(dict(stage='natural_propagation',days=7,phase_nodes=128,plane_stride=16)),flush=True)
    raw=ROOT/'research/runs/solar_shield_array/active/natural_traffic.npz'
    if args.reuse_raw:
        meta=json.loads(raw.with_suffix('.json').read_text())
        signature=hashlib.sha256(inspect.getsource(NaturalTraffic.propagate).encode()).hexdigest()
        initial_signature=hashlib.sha256(inspect.getsource(NaturalTraffic.initial_states).encode()).hexdigest()
        if meta['source_hashes']!=sources() or meta['propagation_function_sha256']!=signature or meta.get('initial_function_sha256')!=initial_signature or constants_changed(meta['constants']) or digest(raw)!=meta['raw_sha256']:
            raise ValueError('Natural propagation cache identity mismatch')
        cache=np.load(raw);f.times=cache['t'];f.states=cache['state'];f.grid_j=cache['plane_nodes']
        f.phase_nodes=meta['phase_nodes'];f.phases=np.arange(f.phase_nodes)*2*np.pi/f.phase_nodes
        n=f.layout.planes
        f.segments=[f.grid_j[f.grid_j<n//2],f.grid_j[f.grid_j==n//2],f.grid_j[f.grid_j>n//2]]
        f.grid_shape=f.states.shape[1:-1];f.nfev=meta['nfev'];f.cached_time=None
    else:
        f.propagate()
        np.savez_compressed(raw,t=f.times,state=f.states,plane_nodes=f.grid_j)
        src=sources()
        meta=dict(schema='terluna.research.natural-traffic-raw/1',source_hashes=src,constants=constants_used(src),
            propagation_function_sha256=hashlib.sha256(inspect.getsource(NaturalTraffic.propagate).encode()).hexdigest(),
            initial_function_sha256=hashlib.sha256(inspect.getsource(NaturalTraffic.initial_states).encode()).hexdigest(),
            raw_sha256=digest(raw),days=7,phase_nodes=128,plane_stride=16,max_step_s=1200,rtol=2e-11,nfev=f.nfev)
        raw.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(json.dumps(dict(stage='propagated',nfev=f.nfev,grid_members=int(np.prod(f.grid_shape)),elapsed_s=time.monotonic()-start)),flush=True)
    verification=replay(f);print(json.dumps(dict(stage='independent_replay',**verification)),flush=True)
    if not verification['interpolation_within_50m']:
        raise ValueError('Natural interpolation exceeds placement allowance; do not publish geometry')
    receiver=disk_points(12,945)*f.layout.protected_radius
    a=np.arange(512)*2*np.pi/512
    rim=f.layout.protected_radius*np.column_stack([np.cos(a),np.sin(a)])
    xy=np.vstack([receiver,rim]);suns=disk_points(3,457)
    dates=[0.,.25,.5,1.,2.,3.];rows=[];witnesses=[];collision=[]
    for day in dates:
        t=day*K.JULIAN_DAY;hits=[]
        for s in suns:
            source=sun_source(env.at(t),f.frames(t),s)
            h=NaturalRayIndex(f,t,source).intercept(xy);hits.append(h)
            for k in np.flatnonzero(~h)[:8]:
                witnesses.append(dict(day=day,sun_xy=s.tolist(),receiver_xy_m=xy[k].tolist(),rim=bool(k>=len(receiver))))
        h=np.array(hits)
        row=dict(day=day,ray_coverage=float(h[:,:len(receiver)].mean()),
            full_sampled_sun_receiver_fraction=float(np.all(h[:,:len(receiver)],axis=0).mean()),
            rim_ray_coverage=float(h[:,len(receiver):].mean()),missed_rays=int(np.count_nonzero(~h)))
        rows.append(row);print(json.dumps(dict(stage='natural_coverage',**row)),flush=True)
        if day in (0.,1.,3.):
            c=pole_collisions(f,t,half_width=200000.);collision.append(c)
            c['scope']='Exact static square intersections in identities selected around the rejected schedule\'s nominal return clock; natural motion can carry them away from that pole. Not a search of all encounters.'
            print(json.dumps(dict(stage='natural_collisions',**c)),flush=True)
    budget_rows=[]
    planes=np.arange(0,f.layout.planes,16)+8
    j,par,ph=np.meshgrid(planes,[0,1],(np.arange(128)+.371)*2*np.pi/128,indexing='ij')
    j=j.ravel();par=par.ravel();ph=ph.ravel();weights=f.count[j]
    scale=f.inventory*f.layout.side**2*env.sigma
    p=SCENARIO['propulsion'];factor=p['exhaust_velocity_m_s']/(2*p['efficiency']*np.cos(np.deg2rad(p['cant_deg'])))
    for day in dates:
        _,hi,d=f.power_bounds(day*K.JULIAN_DAY,j,ph,par)
        power=float(np.average(hi,weights=weights)*scale*factor)
        budget_rows.append(dict(day=day,photon_cancellation_power_upper_W=power,
            minimum_beam_margins_deg=np.rad2deg(d['margins'].min(axis=0)).tolist()))
    mean_power=float(np.trapezoid([x['photon_cancellation_power_upper_W'] for x in budget_rows],np.array(dates)*K.JULIAN_DAY)/(dates[-1]*K.JULIAN_DAY))
    out=dict(schema='terluna.research.natural-traffic-release/1',
        producer=dict(source_hashes=before,constants=constants_used(before)),
        inventory=f.inventory,optical_base_mass_kg=scale,
        propagation=dict(days=7,grid_members=int(np.prod(f.grid_shape)),plane_stride=16,
            phase_nodes=128,max_step_s=1200.,rtol=2e-11,nfev=f.nfev,raw_path=str(raw.relative_to(ROOT)),raw_sha256=digest(raw)),
        independent_replay=verification,coverage=dict(horizon_days=dates[-1],area_receivers=len(receiver),rim_receivers=len(rim),
            solar_directions=len(suns),series=rows,miss_witnesses=witnesses),
        nominal_return_subset_intersections=collision,
        photon_cancellation=dict(mean_power_upper_W=mean_power,
            propellant_rate_upper_kg_s=mean_power*2*p['efficiency']/p['exhaust_velocity_m_s']**2,series=budget_rows),
        acceptance=dict(global_four_lunar_radius_coverage_demonstrated=False,
            low_energy_gap_repair_demonstrated=False,delivered_electrical_power_demonstrated_W=None),
        evidence='Common-epoch gravitational IVPs from the same tile positions with lunar circular velocities, interpolated to actual integer tiles with independent-member replay. Ideal photon-force cancellation only; no coverage or collision correction is applied.',
        limitations=['The photon-cancellation bound is not a complete fleet operating budget',
            'The seven-day interpolation failed the placement gate; only the separately checked three-day prefix is used for geometry',
            'Recorded intersections reject uncorrected natural traffic; the prescribed-family radial separation proof no longer applies',
            'Eight solar directions and finite receiver/date samples do not certify continuous coverage',
            'Neither gap repair nor collision-avoidance thrust is determined by this release diagnostic',
            'Collector/habitat loads, bounded-slew attitudes, plumes and delivered power remain separate'],
        elapsed_s=time.monotonic()-start)
    if identities()!=before:raise ValueError('Sources changed during run')
    (HERE/'results/natural_traffic.json').write_text(compact_series_json(out))
    print(json.dumps(dict(stage='natural_saved',elapsed_s=out['elapsed_s'],photon_cancellation_upper_TW=mean_power/1e12)),flush=True)


if __name__=='__main__':main()
