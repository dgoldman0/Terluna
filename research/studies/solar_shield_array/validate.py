"""Cadence/derivative convergence, eclipse contacts and propagated control cases."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline, CubicHermiteSpline
from scipy.optimize import brentq, minimize_scalar

from shared import constants as K
from protection.dynamics.ephemeris import Ephemeris, DEFAULT_KERNEL, InterpolatedEphemeris
from protection.dynamics.model import geometry, target, evaluate, relative_gravity
from protection.dynamics.optical import solar_visibility, optical_projection, length, covering_radius, light_time_allowance, arriving_ray_vectors
from .run import ROOT, SCENARIO, CORE_DIVERTED, constant_coeff, closure


def eclipse_geometry(ephem, t, coefficients):
    s=ephem.sample([t]); geo=geometry(s); q=target(geo,coefficients)["q"][0]
    sun,earth=arriving_ray_vectors(s["positions"]["sun"][0]-q,s["positions"]["earth"][0]-q,
                                  s["sun_v"][0],s["earth_v"][0],s["sun_a"][0],s["earth_a"][0])
    separation=np.arctan2(length(np.cross(sun,earth)),np.dot(sun,earth))
    sr=np.arcsin(K.SUN_RADIUS/length(sun)); er=np.arcsin(K.EARTH_RADIUS/length(earth))
    return float(separation-sr-er),float(separation+sr-er),float(solar_visibility(sun,earth))


def eclipses(ephem, coefficients, step=600):
    stop=SCENARIO["primary_days"]*K.JULIAN_DAY
    t=np.linspace(0,stop,int(stop/step)+1)
    samples=ephem.sample(t);case=evaluate(samples,geometry(samples),coefficients)
    hit=np.flatnonzero(case["visibility"]<1-1e-10)
    return refine_eclipses(ephem,coefficients,t,hit),samples,case


def refine_eclipses(ephem,coefficients,t,hit):
    groups=np.split(hit,np.flatnonzero(np.diff(hit)>1)+1)
    events=[]
    for group in groups:
        if not len(group):continue
        i,j=group[0],group[-1]
        if i==0 or j==len(t)-1:raise RuntimeError("Eclipse crosses study boundary")
        begin=brentq(lambda ti:eclipse_geometry(ephem,ti,coefficients)[0],t[i-1],t[i],xtol=0.01)
        end=brentq(lambda ti:eclipse_geometry(ephem,ti,coefficients)[0],t[j],t[j+1],xtol=0.01)
        # Geometric separation locates mid-eclipse even when visibility has a flat zero plateau.
        mid=minimize_scalar(lambda ti:eclipse_geometry(ephem,ti,coefficients)[1],bounds=(begin,end),method="bounded",options={"xatol":0.1}).x
        _,inner,minimum=eclipse_geometry(ephem,mid,coefficients)
        umbra_duration=0.
        if inner<0:
            ub=brentq(lambda ti:eclipse_geometry(ephem,ti,coefficients)[1],begin,mid,xtol=0.01)
            ue=brentq(lambda ti:eclipse_geometry(ephem,ti,coefficients)[1],mid,end,xtol=0.01)
            umbra_duration=(ue-ub)/3600
        event_t=np.linspace(begin,end,601)
        ss=ephem.sample(event_t); cc=evaluate(ss,geometry(ss),coefficients)
        msv,mev=arriving_ray_vectors(ss["positions"]["sun"],ss["positions"]["earth"],ss["sun_v"],ss["earth_v"],ss["sun_a"],ss["earth_a"])
        moon_visibility=solar_visibility(msv,mev)
        reqmag=length(cc["required"])
        events.append({"start_tdb_days_from_epoch":begin/K.JULIAN_DAY,"end_tdb_days_from_epoch":end/K.JULIAN_DAY,
                       "penumbra_duration_hours":(end-begin)/3600,"umbra_duration_hours":umbra_duration,
                       "minimum_solar_fraction":minimum,"mid_tdb_days_from_epoch":mid/K.JULIAN_DAY,
                       "moon_solar_fraction_max_when_tile_dark":float(moon_visibility[cc["visibility"]<1e-6].max()) if np.any(cc["visibility"]<1e-6) else None,
                       "mean_dark_holding_mm_s2":float(reqmag[cc["visibility"]<1e-6].mean()*1000) if np.any(cc["visibility"]<1e-6) else None,
                       "solar_deficit_hours":float(np.trapezoid(1-cc["visibility"],event_t)/3600)})
    return events


def nodal_eclipses(ephem,coefficients,step=600):
    """Scan a nodal span in bounded annual chunks, then refine contacts."""
    stop=SCENARIO["nodal_span_years"]*K.JULIAN_YEAR_DAYS*K.JULIAN_DAY
    t=np.linspace(0,stop,int(stop/step)+1)
    chunk=int(K.JULIAN_YEAR_DAYS*K.JULIAN_DAY/step)
    hits=[]
    for begin in range(0,len(t),chunk):
        ss=ephem.sample(t[begin:begin+chunk]);case=evaluate(ss,geometry(ss),coefficients)
        hits.extend((begin+np.flatnonzero(case["visibility"]<1-1e-10)).tolist())
    return {"years":SCENARIO["nodal_span_years"],"scan_step_s":step,
            "events":refine_eclipses(ephem,coefficients,t,np.array(hits,dtype=int))}


def propagate(ephem, coefficients, sigma=0.062, mode="feedback", start_day=0, days=30, mesh_s=600):
    begin=start_day*K.JULIAN_DAY; end=begin+days*K.JULIAN_DAY
    grid=np.linspace(begin-mesh_s,end+mesh_s,int((end-begin)/mesh_s)+3)
    samples=ephem.sample(grid); geo=geometry(samples); nominal=evaluate(samples,geo,coefficients,sigma=sigma,fraction=CORE_DIVERTED)
    reference=CubicHermiteSpline(grid,nominal["q"],nominal["v"])
    command=CubicSpline(grid,nominal["residual"])
    required=CubicSpline(grid,nominal["required"])
    cache=InterpolatedEphemeris(samples)
    y0=np.r_[reference(begin),reference(begin,1)]
    if mode in ("feedback","feedforward"):
        y0[:3]+=np.array([1.,0.,0.]);y0[3:]+=np.array([0.,0.001,0.])
    feedback_time=3600.
    def rhs(t,y):
        q,v=y[:3],y[3:];state=cache.at(t)
        sun,earth=arriving_ray_vectors(state["positions"]["sun"]-q,state["positions"]["earth"]-q,
                                      state["sun_v"],state["earth_v"],state["sun_a"],state["earth_a"])
        visible=solar_visibility(sun,earth)
        optical=optical_projection(required(t),sun,sigma,CORE_DIVERTED,visible)[0]
        thrust=np.zeros(3) if mode=="solar_only" else command(t)
        if mode in ("feedback","eclipse_coast"):
            thrust=thrust+(reference(t)-q)/feedback_time**2+2*(reference(t,1)-v)/feedback_time
        if mode=="eclipse_coast" and visible<0.999:
            thrust=np.zeros(3)
        return np.r_[v,relative_gravity(q,state["positions"],state["moon_a"])+optical+thrust]
    def atmosphere(t,y):return length(y[:3])-SCENARIO["protected_radii"]*K.MOON_RADIUS
    atmosphere.terminal=True;atmosphere.direction=-1
    output_t=np.linspace(begin,end,int(days*24*6)+1)
    solution=solve_ivp(rhs,(begin,end),y0,method="DOP853",rtol=2e-11,
                       atol=np.array([0.001]*3+[1e-8]*3),max_step=mesh_s,
                       t_eval=output_t,events=atmosphere)
    delta=solution.y[:3].T-reference(solution.t)
    error=length(delta)
    all_shadow_margin=[]
    aperture=float(nominal["covering_radius"].max()+SCENARIO["formation_margin_m"])
    for ti,q in zip(solution.t,solution.y[:3].T):
        state=cache.at(ti);u=state["positions"]["sun"]/length(state["positions"]["sun"])
        d=np.dot(q,u); offset=length(q-d*u)
        needed=(covering_radius(d,length(state["positions"]["sun"]),SCENARIO["protected_radii"]*K.MOON_RADIUS,offset)+light_time_allowance(d,state["moon_v"],state["sun_v"])) if d>0 else np.inf
        all_shadow_margin.append(aperture-needed)
    margin=np.array(all_shadow_margin)
    failed=np.flatnonzero(margin<0)
    return {"mode":mode,"start_tdb_day":start_day,"requested_days":days,
            "propagated_days":float((solution.t[-1]-begin)/K.JULIAN_DAY),
            "success":bool(solution.success),"integrator": "DOP853; point masses with Hermite-interpolated DE440s",
            "function_evaluations":solution.nfev,"ephemeris_mesh_s":mesh_s,
            "initial_position_perturbation_m":1 if mode in ("feedback","feedforward") else 0,
            "initial_velocity_perturbation_m_s":0.001 if mode in ("feedback","feedforward") else 0,
            "max_position_error_m":float(error.max()),"end_position_error_m":float(error[-1]),
            "minimum_full_sun_coverage_margin_m":float(margin.min()),
            "first_coverage_failure_hours":float((solution.t[failed[0]]-begin)/3600) if len(failed) else None,
            "entered_protected_atmosphere":bool(len(solution.t_events[0])),
            "boundary":"Unbounded ideal electrical feedback; no plume model, sensor delay, tile elasticity or receiver/storage mass. This tests trajectory integration and controller dependence."}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=ROOT/"research/runs/solar_shield_array")
    parser.add_argument("--kernel",type=Path,default=DEFAULT_KERNEL)
    args=parser.parse_args();p=args.out/"holding.json"
    holding=json.loads(p.read_text());ephem=Ephemeris(args.kernel,SCENARIO["epoch_tdb"])
    base=constant_coeff(78000);chosen=np.array(holding["trajectory_search"]["variable_distance"]["coefficients_km"])
    result={"schema":"terluna.protection.holding-validation/1","holding_product_sha256":hashlib.sha256(p.read_bytes()).hexdigest(),
            "scope":"Numerical derivatives, cadence, ephemeris interpolation, finite-Sun contacts and trajectory integration",
            "physical_or_engineering_validation":False,"cadence":{},"propagation":{}}
    print("eclipse contacts",flush=True)
    for label,c in [("baseline",base),("variable_distance",chosen)]:
        events,samples,case=eclipses(ephem,c)
        result[label+"_eclipses"]=events
        result["cadence"][label]={"step_600_s_closure":closure(case,samples)}
    print("19-year eclipse scan",flush=True)
    result["nodal_eclipses"]=nodal_eclipses(ephem,chosen)
    t=np.linspace(0,365.25*K.JULIAN_DAY,500)
    result["derivative_convergence"]={"max_10_vs_30_s_m_s2":float(length(ephem.acceleration("moon",t,10)-ephem.acceleration("moon",t,30)).max()),
                                      "max_90_vs_30_s_m_s2":float(length(ephem.acceleration("moon",t,90)-ephem.acceleration("moon",t,30)).max())}
    sigma=holding["trajectory_search"]["variable_distance"]["fine_closure"]["total_areal_mass_kg_m2"]
    for mode,days in [("feedback",30),("feedforward",30),("solar_only",2)]:
        print("propagate",mode,flush=True)
        result["propagation"][mode]=propagate(ephem,chosen,sigma=sigma,mode=mode,days=days)
    result["propagation"]["feedback_fine"]=propagate(ephem,chosen,sigma=sigma,days=3,mesh_s=120)
    result["propagation"]["feedback_coarse_short"]=propagate(ephem,chosen,sigma=sigma,days=3)
    if result["variable_distance_eclipses"]:
        event=result["variable_distance_eclipses"][0]
        first=event["start_tdb_days_from_epoch"]-0.25
        result["propagation"]["eclipse_with_power"]=propagate(ephem,chosen,sigma=sigma,start_day=first,days=2)
        result["propagation"]["eclipse_without_power"]=propagate(ephem,chosen,sigma=sigma,mode="eclipse_coast",start_day=first,days=2)
    if result["nodal_eclipses"]["events"]:
        worst=max(result["nodal_eclipses"]["events"],key=lambda x:x["umbra_duration_hours"])
        first=worst["start_tdb_days_from_epoch"]-0.25
        result["worst_eclipse"]=worst
        result["propagation"]["worst_eclipse_with_power"]=propagate(ephem,chosen,sigma=sigma,start_day=first,days=2)
        result["propagation"]["worst_eclipse_without_power"]=propagate(ephem,chosen,sigma=sigma,mode="eclipse_coast",start_day=first,days=2)
    ephem.close()
    result["producer"]={"runner_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                        "source_hashes":holding["producer"]["source_hashes"]}
    (args.out/"validation.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("wrote validation.json",flush=True)


if __name__=="__main__":main()
