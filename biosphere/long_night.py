"""Periodic resource budgets; functional requirements, not organism validation.

Carbon units are day-reference-maintenance equivalents per fixed structural
biomass. Parameters are explicit hypothetical traits, not assigned to species.
A proof of periodic feasibility for this model is in research/findings.md.
Oxygen concentrations use g/m3 (numerically mg/L); exchange is first-order.
"""
from __future__ import annotations
from collections import deque
from dataclasses import dataclass
import math
import numpy as np

from shared.constants import SYNODIC_MONTH_DAYS as PERIOD_DAYS


def validate_trace(net,dt):
    a=np.asarray(net,dtype=float)
    if a.ndim != 1 or len(a)<2 or not np.all(np.isfinite(a)):
        raise ValueError('Finite one-dimensional periodic trace required')
    if not math.isfinite(dt) or dt<=0: raise ValueError('Positive dt required')
    return a


def periodic_storage_requirement(net,dt):
    """Necessary and sufficient finite-capacity condition for fixed net production.

    dS/dt=net(t), surplus can be spilled at S=C, initial phase can be chosen.
    Need integral(net)>=0 and C>=maximum cumulative deficit on any circular
    interval of length <= one period. The O(N) deque tests *all* such intervals.
    No rate constraints, reserve-dependent physiology, mortality or growth feedback.
    """
    a=validate_trace(net,dt); n=len(a)
    increments=a*dt; total=float(increments.sum())
    z=np.r_[0.,np.cumsum(np.tile(increments,2))]
    candidates=deque([0]); deficit=0.; start=end=0
    for j in range(1,2*n+1):
        while candidates and candidates[0]<j-n: candidates.popleft()
        d=float(z[candidates[0]]-z[j])
        if d>deficit: deficit=d;start=candidates[0];end=j
        while candidates and z[candidates[-1]]<=z[j]: candidates.pop()
        candidates.append(j)
    tol=1e-11*max(1.,float(np.abs(increments).sum()))
    return dict(cycle_net=total,cycle_balance_feasible=total>=-tol,
        minimum_capacity=deficit if total>=-tol else None,
        worst_single_cycle_deficit=deficit,
        deficit_duration_days=(end-start)*dt,
        deficit_start_phase=(start%n)/n,deficit_end_phase=(end%n)/n)


def simulate_store(net,dt,capacity,cycles=4,initial=None):
    a=validate_trace(net,dt)
    if not math.isfinite(capacity) or capacity<0 or cycles<1:
        raise ValueError('Nonnegative capacity and positive cycles required')
    store=capacity if initial is None else float(initial)
    if not 0<=store<=capacity: raise ValueError('Initial reserve exceeds capacity')
    s0=store; surplus=unmet=produced=consumed=0.;minimum=store;trace=[]
    for c in range(cycles):
        cycle_unmet=cycle_surplus=0.
        for v in a:
            amount=float(v*dt)
            if amount>=0:
                accepted=min(amount,capacity-store);store+=accepted
                surplus+=amount-accepted;cycle_surplus+=amount-accepted;produced+=amount
            else:
                demand=-amount;actual=min(store,demand);store-=actual
                consumed+=actual;unmet+=demand-actual;cycle_unmet+=demand-actual
            minimum=min(minimum,store)
        trace.append(dict(cycle=c+1,end_reserve=store,unmet=cycle_unmet,spill=cycle_surplus))
    return dict(initial_reserve=s0,final_reserve=store,minimum_reserve=minimum,
        unmet_demand=unmet,spill=surplus,positive_net_input=produced,
        delivered_net_demand=consumed,
        ledger_residual=s0+produced-consumed-surplus-store,cycles=trace)


def carbon_trace(period_days=PERIOD_DAYS,steps=1440,mean_light_assimilation=3.,
                 night_demand_fraction=.25,day_temp_k=288.,night_temp_k=288.,
                 q10=2.,demand_reference_k=288.,light_trace=None,temperature_trace=None,light_reference=None):
    """Net carbon rate at fixed biomass and maintenance.

    Light is a half sine, normalized to mean 1 during illuminated half. A is
    usable assimilate *after biosynthetic/conversion overhead*. Growth is zero;
    positive spill is unallocated surplus, never claimed as population growth.
    Night suppression and Q10 are distinct scenario assumptions.
    Optional traces allow a climate case's irradiance and temperature to be used.
    """
    if period_days<=0 or steps<8 or mean_light_assimilation<0 or not 0<=night_demand_fraction<=1 or q10<=0:
        raise ValueError('Invalid carbon scenario')
    phase=(np.arange(steps)+.5)/steps
    light=np.maximum(0.,np.sin(2*np.pi*phase)) if light_trace is None else np.asarray(light_trace,float)
    if light.shape!=(steps,) or np.any(light<0) or not np.all(np.isfinite(light)): raise ValueError('Invalid light')
    day=light>1e-12
    if light_reference is not None and (not math.isfinite(light_reference) or light_reference<=0):
        raise ValueError('Positive light reference required')
    norm=light_reference if light_reference is not None else max(float(light[day].mean()) if np.any(day) else 0.,1e-12)
    normalized=light/norm
    temp=np.where(day,day_temp_k,night_temp_k) if temperature_trace is None else np.asarray(temperature_trace,float)
    if temp.shape!=(steps,) or not np.all(np.isfinite(temp)) or np.any(temp<=0): raise ValueError('Invalid temperature')
    respiration=q10**((temp-demand_reference_k)/10)*np.where(day,1.,night_demand_fraction)
    production=mean_light_assimilation*normalized
    return production-respiration,period_days/steps,dict(production=production,respiration=respiration,
                                                       temperature_K=temp,light=normalized)


def oxygen_dark_end(initial,saturation,respiration,exchange_per_day,dark_days):
    if min(initial,saturation,respiration,exchange_per_day,dark_days)<0:
        raise ValueError('Oxygen inputs must be nonnegative')
    if exchange_per_day==0:return initial-respiration*dark_days
    k=exchange_per_day
    return saturation-respiration/k+(initial-saturation+respiration/k)*math.exp(-k*dark_days)


def oxygen_step(initial,production,respiration,k,saturation,capacity,dt):
    """Exact constant-forcing step, with recorded upper outgassing/lower shortfall.

    Capacity is an imposed dissolved-oxygen ceiling; this models immediate extra
    outgassing at that ceiling, not a measured bubble-nucleation law.
    """
    eq=saturation+(production-respiration)/k
    beta=-math.expm1(-k*dt)/k
    end=eq+(initial-eq)*math.exp(-k*dt)
    integral=eq*dt+(initial-eq)*beta
    spill=unmet=0.
    if end>capacity or end<0:
        bound=capacity if end>capacity else 0.
        ratio=(bound-eq)/(initial-eq)
        hit=max(0.,min(dt,-math.log(max(ratio,1e-300))/k))
        segment=eq*hit+(initial-eq)*(-math.expm1(-k*hit)/k)
        integral=segment+bound*(dt-hit)
        balance_rate=production-respiration+k*(saturation-bound)
        spill=max(balance_rate,0)*(dt-hit)
        unmet=max(-balance_rate,0)*(dt-hit)
        end=bound
    exchange=k*(saturation*dt-integral)
    return end,exchange,spill,unmet


def oxygen_periodic(period_days=PERIOD_DAYS,steps=1440,saturation_g_m3=8.,
                    respiration_g_m3_day=.5,mean_day_production_g_m3_day=1.5,
                    exchange_per_day=.1,threshold_g_m3=2.,capacity_g_m3=9.6):
    """Periodic mixed oxygen box with explicit inventory ceiling and shortfalls.

    Rate constants, saturation, ceiling and demand are scenario inputs. The cap
    avoids treating arbitrarily supersaturated water as usable oxygen storage.
    A prescribed aerobic demand that exceeds supply is recorded as unmet.
    """
    from scipy.optimize import brentq
    vals=[period_days,saturation_g_m3,respiration_g_m3_day,mean_day_production_g_m3_day,
          threshold_g_m3,capacity_g_m3,exchange_per_day]
    if not all(math.isfinite(v) and v>=0 for v in vals) or period_days<=0 or capacity_g_m3<=0 or steps<8 or exchange_per_day<=0:
        raise ValueError('Invalid oxygen inputs')
    dt=period_days/steps; ph=(np.arange(steps)+.5)/steps
    daylight=np.maximum(0,np.sin(2*np.pi*ph))
    production=mean_day_production_g_m3_day*daylight/daylight[daylight>0].mean()
    k=exchange_per_day
    def advance(y,collect=False):
        values=[y];exchange=spill=unmet=0.
        for p in production:
            y,e,s,u=oxygen_step(y,float(p),respiration_g_m3_day,k,saturation_g_m3,capacity_g_m3,dt)
            exchange+=e;spill+=s;unmet+=u
            if collect:values.append(y)
        return y,exchange,spill,unmet,values
    initial=brentq(lambda y:advance(y)[0]-y,0,capacity_g_m3,xtol=1e-12)
    y,exchange,spill,unmet,trace=advance(initial,True)
    minimum=min(trace);maximum=max(trace)
    ledger=y-initial-float(production.sum())*dt+respiration_g_m3_day*period_days-exchange+spill-unmet
    return dict(minimum_g_m3=float(minimum),maximum_g_m3=float(maximum),
        threshold_g_m3=threshold_g_m3,threshold_satisfied=minimum>=threshold_g_m3,
        prescribed_aerobic_demand_feasible=unmet<1e-9,
        capacity_g_m3=capacity_g_m3,outgassing_above_cap_g_m3=float(spill),
        unmet_respiration_g_m3=float(unmet),
        periodic_residual=float(y-initial),mass_balance_residual=float(ledger),
        exchange_per_day=k,respiration_g_m3_day=respiration_g_m3_day,
        trace_g_m3=np.asarray(trace))
