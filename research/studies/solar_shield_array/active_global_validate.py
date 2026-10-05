"""Focused handovers, force convergence and feedback checks for global traffic."""
import argparse
import json
import time

import numpy as np
from scipy.integrate import solve_ivp

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.active_global import Traffic, RayIndex, sun_source
from protection.dynamics.cycling import ray_geometry, visible_sun
from protection.dynamics.optical import length,unit
from .active_global_run import sources, disk_points, budget
from .fleet_run import environment,digest
from .cycling_search import ROOT,HERE,SCENARIO
from .cycling_analysis import compact_series_json


def handover(fleet,start,duration=120.,step=.5):
    """Frozen broad-phase supersets, exact integer tiles at every new date.

    The extra 1 km envelope covers the <=120 s change of ribbon light-time
    curves; this is validated against freshly rebuilt full indices below.
    No tile position/normal or ephemeris force is copied from the midpoint.
    """
    import shapely
    R=fleet.layout.protected_radius
    angle=np.arange(12)*np.pi/6
    xy=np.vstack([[0.,0.],[R/2,0.],[0.,-R/2],R*np.column_stack([np.cos(angle),np.sin(angle)])])
    suns=disk_points(3,671)
    # Explicit limb sources complement equal-area interior quadrature.
    suns=np.vstack([suns,[[1.,0.],[0,1.],[-1.,0.],[0,-1.]]])
    middle=start+duration/2;groups=[];all_ids=set()
    for s in suns:
        source=sun_source(fleet.env.at(middle),fleet.frames(middle),s)
        index=RayIndex(fleet,middle,source)
        pair=index.tree.query(shapely.points(xy),predicate='dwithin',distance=index.buffer+1000.)
        for k,point in enumerate(xy):
            lines=np.unique(pair[1,pair[0]==k]//32)
            ids=set()
            for line in lines:
                j=index.j[line];par=index.parity[line]
                radius=fleet.r0[j]+par*fleet.layout.radial_step
                local=index.source_local;w=np.r_[0.,point];ray=unit(local-w)
                wd=np.dot(w,ray);lam=-wd+np.sqrt(wd*wd+radius**2-np.dot(w,w))
                q=w+lam*ray
                q+=(fleet.env.at(middle)['moon_v']@fleet.frames(middle))*(q[0]/K.SPEED_OF_LIGHT)
                transverse=q[1]*np.cos(fleet.alpha[j])+q[2]*np.sin(fleet.alpha[j])
                phase=np.arcsin(np.clip(transverse/fleet.r0[j],-.9,.9))
                near=int(np.round(((phase-fleet.phase[j]-fleet.omega[j]*middle)*fleet.count[j]/(2*np.pi)-par)/2))*2+par
                width=int(np.ceil(duration/2*abs(fleet.omega[j])*fleet.count[j]/(4*np.pi)))+4
                ids.update((int(j),int((near+2*l)%fleet.count[j])) for l in range(-width,width+1))
            all_ids.update(ids);groups.append((s,k,ids))
    ids=np.array(sorted(all_ids));lookup={tuple(x):i for i,x in enumerate(ids)}
    groups=[(s,k,np.array([lookup[x] for x in sorted(v)],int)) for s,k,v in groups]
    half=fleet.layout.clear_side/2;hist=[];miss=[];min_margins=np.full(2,np.inf)
    normals=[];max_slew=0.;last_n=None
    times=start+np.arange(round(duration/step)+1)*step
    for ti in times:
        d=fleet.facets(ti,ids[:,0],ids[:,1]);frame=fleet.frames(ti);sample=fleet.env.at(ti)
        q=d['state'][:,:3];tau=np.maximum(q@frame[:,0],0)/K.SPEED_OF_LIGHT
        q=(q-(d['state'][:,3:]+sample['moon_v'])*tau[:,None])@frame
        n=d['normal']@frame;aa=d['a']@frame;bb=d['b']@frame
        min_margins=np.minimum(min_margins,d['margins'].min(axis=0))
        if last_n is not None:
            change=np.arctan2(length(np.cross(last_n,d['normal'])),np.abs(np.sum(last_n*d['normal'],axis=1)))
            max_slew=max(max_slew,float(change.max()/step))
        last_n=d['normal']
        assigned=[]
        for s,k,candidates in groups:
            source=sun_source(sample,frame,s)@frame;w=np.r_[0.,xy[k]];ray=unit(source-w)
            den=n[candidates]@ray
            distance=np.sum((q[candidates]-w)*n[candidates],axis=1)/np.where(np.abs(den)>1e-12,den,1e-100)
            p=w+distance[:,None]*ray-q[candidates]
            yes=(distance>0)&(np.abs(np.sum(p*aa[candidates],axis=1))<=half)&(np.abs(np.sum(p*bb[candidates],axis=1))<=half)
            if np.any(yes):assigned.append(int(candidates[np.flatnonzero(yes)[0]]))
            else:
                assigned.append(-1)
                if len(miss)<32:miss.append(dict(time_s=float(ti),receiver_xy_m=xy[k].tolist(),sun_xy=s.tolist()))
        hist.append(assigned)
    hist=np.asarray(hist);switches=(hist[1:]!=hist[:-1])&(hist[1:]>=0)&(hist[:-1]>=0)
    # Independent broad-phase reconstruction at both ends and midpoint.
    fresh=[]
    for ti in (times[0],middle,times[-1]):
        for s in suns:
            src=sun_source(fleet.env.at(ti),fleet.frames(ti),s)
            fresh.extend(RayIndex(fleet,ti,src).intercept(xy).tolist())
    return dict(start_day=start/K.JULIAN_DAY,duration_s=duration,cadence_s=step,
        receivers=len(xy),solar_directions=len(suns),unique_integer_tiles=len(ids),
        ray_tests=int(hist.size),uncovered_ray_samples=int(np.count_nonzero(hist<0)),
        observed_responsibility_changes=int(switches.sum()),
        minimum_receiver_source_handover_changes=int(switches.sum(axis=0).min()),
        rebuilt_indices_all_covered=bool(np.all(fresh)),
        minimum_beam_margins_deg=np.rad2deg(min_margins).tolist(),
        maximum_sampled_facet_normal_rate_rad_s=max_slew,miss_witnesses=miss,
        scope='Explicit local temporal handovers across the whole receiver disk, including its rim and solar limb; this does not certify every receiver between the 12-hour global samples')


def feedback_replay(fleet,start=15*K.JULIAN_DAY,days=8.,step=120.):
    """Independent ODE with unknown, rapidly varying bounded photon force.

    Nominal feedforward cancels centre gravity at the target, without using
    actual-state gravity or actual mutual illumination. The disturbance is a
    specified stress signal in [0,1], not a claimed fleet shadow solution.
    """
    planes=np.array([0,288,576,864,1152,1440,1728,2016,2303])
    j=np.repeat(planes,2)
    k=np.rint((np.tile([0.,np.pi],len(planes))-fleet.phase[j]-fleet.omega[j]*start)*fleet.count[j]/(2*np.pi)).astype(int)%fleet.count[j]
    initial,_,_=fleet.path(start,j,k);rng=np.random.default_rng(2026)
    error=unit(rng.normal(size=(len(j),3)))*20.
    initial[:,:3]+=error;initial[:,3:]+=unit(rng.normal(size=(len(j),3)))*.001
    tau=600.;kp=1/tau**2;kd=2/tau;cap=.002
    def evaluate(t,state):
        d=fleet.facets(t,j,k);sample=fleet.env.at(t);target=d['state']
        control=d['acc']-fleet.env.gravity(t,target[:,:3])+kp*(target[:,:3]-state[:,:3])+kd*(target[:,3:]-state[:,3:])
        control*=np.minimum(1,cap/np.maximum(length(control),1e-100))[:,None]
        rays=ray_geometry(state[:,:3],sample)
        cos=np.abs(np.sum(-unit(rays['sun'])*d['normal'],axis=1))
        photon=(2*fleet.env.central_fraction*K.SOLAR_CONSTANT/(K.SPEED_OF_LIGHT*fleet.env.sigma)
            *(K.AU/length(rays['sun']))**2*cos**2*visible_sun(rays['sun'],rays['earth'],rays['moon'])
            *(fleet.layout.clear_side/fleet.layout.side)**2)
        lit=.5+.5*np.sin((t-start)/137+np.arange(len(j))*.731)
        sail=photon[:,None]*lit[:,None]*d['normal']
        return fleet.env.gravity(t,state[:,:3])+sail+control,control,target
    def rhs(t,y):
        state=y.reshape(-1,6);acc,_,_=evaluate(t,state)
        return np.column_stack([state[:,3:],acc]).ravel()
    times=np.arange(start,start+days*K.JULIAN_DAY+1,300.)
    sol=solve_ivp(rhs,[times[0],times[-1]],initial.ravel(),t_eval=times,method='DOP853',
        rtol=2e-11,atol=np.tile([1e-4]*3+[1e-8]*3,len(j)),max_step=step)
    if not sol.success:raise RuntimeError(sol.message)
    states=sol.y.T.reshape(len(times),len(j),6);errors=[];thrust=[]
    for t,state in zip(times,states):
        _,u,target=evaluate(t,state);errors.append(length(state[:,:3]-target[:,:3]));thrust.append(length(u))
    errors=np.array(errors);thrust=np.array(thrust)
    path=ROOT/'research/runs/solar_shield_array/active/global_feedback.npz'
    np.savez_compressed(path,t=times,state=states,planes=j,slots=k,tracking_error=errors,control=thrust)
    return dict(days=days,start_day=start/K.JULIAN_DAY,integer_tiles=len(j),
        maximum_position_error_m=float(errors.max()),final_maximum_position_error_m=float(errors[-1].max()),
        maximum_control_m_s2=float(thrust.max()),cap_m_s2=cap,response_s=tau,
        capped_samples=int(np.count_nonzero(thrust>=cap*(1-1e-10))),
        max_step_s=step,rtol=2e-11,nfev=sol.nfev,
        completed_cycles_range=[float((times[-1]-times[0])*abs(fleet.omega[j]).min()/(2*np.pi)),float((times[-1]-times[0])*abs(fleet.omega[j]).max()/(2*np.pi))],
        raw_path=str(path.relative_to(ROOT)),raw_sha256=digest(path),
        scope='Bounded PD robustness stress replay for selected actual member identities; synthetic unknown mutual illumination varies from zero to one every 14.3 min. It is not an all-member coupled optical replay.')


def finite_extent(fleet):
    j=np.repeat(np.array([0,288,576,864,1152,1728,2303]),32)
    ph=np.tile(np.arange(32)*2*np.pi/32,len(j)//32);par=np.arange(len(j))%2
    nodes,weights=np.polynomial.legendre.leggauss(3)
    maximum=0.;worst=None
    for day in (0,7.5,15,22.5,30):
        d=fleet.facets(day*K.JULIAN_DAY,j,phase=ph,parity=par)
        centre=d['state'][:,:3];average=np.zeros_like(centre)
        for ix,x in enumerate(nodes):
            for iy,y in enumerate(nodes):
                q=centre+fleet.layout.side/2*(x*d['a']+y*d['b'])
                average+=weights[ix]*weights[iy]/4*fleet.env.gravity(day*K.JULIAN_DAY,q)
        err=length(average-fleet.env.gravity(day*K.JULIAN_DAY,centre))
        if err.max()>maximum:maximum=float(err.max());worst=day
    return dict(square_gravity_quadrature='3 by 3 Gauss-Legendre',cases=len(j)*5,
        maximum_centre_force_error_m_s2=maximum,worst_day=worst,
        scope='Sampled finite-extent correction; centre force is adequate relative to the reported translational controls, not a membrane stress calculation')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--mode',choices=['geometry','budget','feedback','assemble'],required=True);a=p.parse_args()
    env=environment(30);fleet=Traffic(env)
    ids={**sources(),str(__file__).split('/Terluna/')[-1]:digest(__file__)}
    out=dict(producer=dict(source_hashes=ids,constants=constants_used(ids)))
    if a.mode=='geometry':
        out['handover_windows']=[handover(fleet,day*K.JULIAN_DAY) for day in (0,7.5,15,22.5,29.99)]
        out['handover_refinement']=handover(fleet,15*K.JULIAN_DAY,32.,.125)
        out['finite_extent']=finite_extent(fleet)
    elif a.mode=='budget':
        out['spatial_refinement']=budget(fleet,np.linspace(0,30*K.JULIAN_DAY,121),plane_stride=8,phase_count=256)
        out['temporal_refinement']=budget(fleet,np.linspace(0,30*K.JULIAN_DAY,241),plane_stride=16,phase_count=128)
    elif a.mode=='feedback':out['feedback_replay']=feedback_replay(fleet)
    else:
        base=json.loads((HERE/'results/active_global.json').read_text())
        for mode in ('geometry','budget','feedback'):
            raw=ROOT/f'research/runs/solar_shield_array/active/global_validate_{mode}.json'
            data=json.loads(raw.read_text())
            if data['producer']!=out['producer']:raise ValueError('Validation source mismatch')
            out.update({k:v for k,v in data.items() if k!='producer'})
        out['schema']='terluna.research.active-global-validation/1'
        out['producer']['input_products']={'active_global.json':digest(HERE/'results/active_global.json')}
        base_mean=np.array(base['budget']['mean_power_bounds_W'])
        out['force_convergence']=dict(
            spatial_relative_change=(np.array(out['spatial_refinement']['mean_power_bounds_W'])/base_mean-1).tolist(),
            temporal_relative_change=(np.array(out['temporal_refinement']['mean_power_bounds_W'])/base_mean-1).tolist())
        census=np.array(base['integer_census']['mean_acceleration_bounds_m_s2'])
        match=next(x for x in base['budget']['series'] if x['day']==15.)
        quad=np.array([match['mean_acceleration_lower_m_s2'],match['mean_acceleration_upper_m_s2']])
        out['force_convergence']['full_integer_census_relative_change']=(census/quad-1).tolist()
        out['acceptance']=dict(global_four_lunar_radius_coverage_demonstrated=False,
            delivered_electrical_power_demonstrated_W=None,energy_optimum=False)
        out['evidence']='Independent geometry/index tests, local temporal refinement, force quadrature and bounded feedback stress replay; not empirical validation or a continuous all-fleet optical certificate'
        (HERE/'results/active_global_validation.json').write_text(compact_series_json(out));return
    path=ROOT/f'research/runs/solar_shield_array/active/global_validate_{a.mode}.json'
    path.write_text(compact_series_json(out))
    print(json.dumps(dict(saved=str(path))),flush=True)


if __name__=='__main__':main()
