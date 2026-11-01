"""Positive finite-time Poisson capture with bounded cell reconstruction.

Constant-rate replenished bath only. Actual donor masses are added before
projection to three moments; positive quadrature is split at destination faces.
Poisson tails are omitted without renormalization, with analytic moment bounds.
An absorbing upper mass boundary records the first crossing and stops capture.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import time

import numpy as np
from scipy.special import gammaln
from scipy.stats import poisson

from atmosphere.snow import collision_entropy as entropy
from atmosphere.snow import collision_moments as moments


@dataclass
class PoissonRule:
    counts: np.ndarray
    weights: np.ndarray
    omitted_count_moments: np.ndarray  # analytical upper bounds on the omitted count moments
    normalization_residual: float


def poisson_tail_bound(threshold,mean,*,upper):
    """Chernoff bound on P[K>=threshold] or P[K<=threshold].

    Long-double evaluation and upward conversion retain a positive bound even
    when ordinary probability evaluation underflows. This is a mathematical
    tail bound evaluated in floating point, not interval arithmetic.
    """
    if mean == 0:return np.longdouble(float(threshold<=0) if upper else float(threshold>=0))
    if upper and threshold<=mean:return np.longdouble(1)
    if not upper and threshold<0:return np.longdouble(0)
    if not upper and threshold>=mean:return np.longdouble(1)
    k=np.longdouble(threshold);mu=np.longdouble(mean)
    exponent=-mu if k==0 else -mu+k*(1+np.log(mu/k))
    value=np.exp(exponent)
    return np.nextafter(value,np.longdouble(np.inf))


def poisson_rule(mean, tail_tolerance=1e-16, maximum_terms=100000):
    """Mode-centered recurrence; no rescaling of the retained probabilities."""
    if (not math.isfinite(mean) or mean < 0 or not 0 < tail_tolerance < .01
            or maximum_terms < 1):
        raise ValueError('Finite nonnegative mean and a probability tolerance in (0,.01) required')
    if mean == 0:
        return PoissonRule(np.array([0]),np.ones(1),np.zeros(3),0.)
    width=8*math.sqrt(mean)+16
    while True:
        lo=max(0,int(math.floor(mean-width)));hi=int(math.ceil(mean+width))
        if hi-lo+1 > maximum_terms:
            raise RuntimeError('Poisson mixture exceeds term budget')
        omitted=poisson_tail_bound(lo-1,mean,upper=False)+poisson_tail_bound(hi+1,mean,upper=True)
        if omitted <= tail_tolerance:
            break
        width*=1.5
    k=np.arange(lo,hi+1,dtype=np.int64);mode=int(mean)
    if mode < 32:
        log_peak=-mean+mode*math.log(mean)-float(gammaln(mode+1))
    else:
        x=(mean-mode)/mode
        correction=1/(12*mode)-1/(360*mode**3)+1/(1260*mode**5)-1/(1680*mode**7)+1/(1188*mode**9)
        log_peak=-.5*math.log(2*math.pi*mode)+mode*(math.log1p(x)-x)-correction
    w=np.empty(len(k),dtype=np.longdouble);center=mode-lo
    w[center]=np.exp(np.longdouble(log_peak))
    if center:
        ratio=np.arange(mode,lo,-1,dtype=np.longdouble)/np.longdouble(mean)
        w[:center]=w[center]*np.cumprod(ratio)[::-1]
    if mode < hi:
        ratio=np.longdouble(mean)/np.arange(mode+1,hi+1,dtype=np.longdouble)
        w[center+1:]=w[center]*np.cumprod(ratio)
    first=np.longdouble(mean)*(poisson_tail_bound(lo-2,mean,upper=False)+poisson_tail_bound(hi,mean,upper=True))
    second=np.longdouble(mean)**2*(poisson_tail_bound(lo-3,mean,upper=False)+poisson_tail_bound(hi-1,mean,upper=True))+first
    actual_tail=poisson.cdf(lo-1,mean)+poisson.sf(hi,mean)
    return PoissonRule(k,np.asarray(w,float),np.array([omitted,first,second],dtype=np.longdouble),
                       float(w.sum()+np.longdouble(actual_tail)-1))


def number_cdf(grid, state, queries):
    """Unnormalized number CDF with atoms compared in physical mass.

    Convert continuous queries through the same pivot coordinate used to store
    the moments. Comparing atoms after an independent boundary-coordinate
    subtraction can move an exact jump by one floating-point unit.
    """
    query = np.asarray(queries, float)
    total = np.zeros_like(query)
    for cell, measure in enumerate(grid.distributions(state)):
        if measure is None:
            continue
        if measure.atoms is not None:
            x, weight = measure.atoms
            mass = grid.pivots[cell] + grid.widths[cell]*(x-grid.offsets[cell])
            partial = sum(w*(query >= m) for m,w in zip(mass,weight))
        else:
            coordinate = (query-grid.pivots[cell])/grid.widths[cell]+grid.offsets[cell]
            partial = measure.cdf(coordinate)
        total += state[cell,0]*partial
    return total


def first_exit_count(mass,donor,upper):
    if mass >= upper:
        return 0
    k=max(1,int(math.ceil((upper-mass)/donor)))
    while mass+k*donor < upper:k+=1
    while k and mass+(k-1)*donor >= upper:k-=1
    return k


def _density(measure,z):
    A,B=measure.A,measure.B
    if B < 0 and measure.lower <= -A/(2*B) <= measure.upper:
        peak=-A/(2*B);relative=B*(z-peak)**2;peak_value=-.25*A*A/B
    else:
        peak=measure.upper if (measure.upper-measure.lower)*(A+B*(measure.upper+measure.lower)) >= 0 else measure.lower
        relative=(z-peak)*(A+B*(z+peak));peak_value=peak*(A+B*peak)
    return np.exp(relative+(peak_value-measure.logZ))


class ConstantCapture:
    """Finite-time capture operator; grid must supply entropy measures."""
    def __init__(self,grid,donor_mass,rate_s,*,quadrature_order=16,tail_tolerance=1e-16,
                 maximum_terms=100000,maximum_intervals=2000000):
        if (not math.isfinite(donor_mass) or donor_mass <= 0 or not math.isfinite(rate_s)
                or rate_s < 0 or quadrature_order < 4 or quadrature_order > 64
                or int(quadrature_order)!=quadrature_order or not 0 < tail_tolerance < .01):
            raise ValueError('Positive donor, nonnegative rate, order4..64 and a tail tolerance required')
        self.grid=grid;self.donor=donor_mass;self.rate=rate_s
        self.order=int(quadrature_order);self.tolerance=tail_tolerance
        self.maximum_terms=maximum_terms;self.maximum_intervals=maximum_intervals
        self.stats=dict(macrosteps=0,source_cells=0,positive_quadrature_nodes=0,
                        maximum_poisson_mean=0.,maximum_mixture_terms=0,
                        maximum_probability_normalization_residual=0.)

    def step(self,state,dt):
        if not math.isfinite(dt) or dt < 0:
            raise ValueError('Finite nonnegative macrostep required')
        grid=self.grid;source=np.asarray(state,float)
        measures=grid.distributions(source)
        if dt == 0 or self.rate == 0:
            return source.copy(),dict(supplied_mass=0.,upper_exit=np.zeros(3)),np.zeros(3),np.zeros(3)
        rule=poisson_rule(self.rate*dt,self.tolerance,self.maximum_terms)
        output=np.zeros_like(source,dtype=np.longdouble)
        exit_moments=np.zeros(3,dtype=np.longdouble);supplied=np.longdouble(0)
        upper=grid.bounds[-1];d=self.donor
        gx,gw=moments._legendre(self.order)
        self.stats['macrosteps']+=1
        self.stats['maximum_poisson_mean']=max(self.stats['maximum_poisson_mean'],self.rate*dt)
        self.stats['maximum_mixture_terms']=max(self.stats['maximum_mixture_terms'],len(rule.counts))
        self.stats['maximum_probability_normalization_residual']=max(self.stats['maximum_probability_normalization_residual'],abs(rule.normalization_residual))
        for cell,measure in enumerate(measures):
            if measure is None:continue
            self.stats['source_cells']+=1
            population=np.longdouble(source[cell,0])
            unit=np.zeros_like(output);unit_exit=np.zeros(3,dtype=np.longdouble)
            unit_input=np.longdouble(0)
            if measure.atoms is not None:
                for x,atom_weight in zip(*measure.atoms):
                    mass=grid.pivots[cell]+grid.widths[cell]*(x-grid.offsets[cell])
                    crossing=first_exit_count(mass,d,upper)
                    keep=rule.counts < crossing
                    products=mass+d*rule.counts[keep];weights=atom_weight*rule.weights[keep]
                    dest=np.searchsorted(grid.bounds,products,side='right')-1
                    u=(products-grid.pivots[dest])/grid.widths[dest]
                    for power in range(3):np.add.at(unit[:,power],dest,weights.astype(np.longdouble)*u**power)
                    unit_input+=np.sum(weights.astype(np.longdouble)*(d*rule.counts[keep]))
                    # Restoring the complete physical-exit probability is exact;
                    # any remaining omitted retained trajectories stay omitted.
                    if crossing <= rule.counts[-1]:
                        probability=atom_weight*poisson.sf(crossing-1,self.rate*dt)
                        exited=mass+d*crossing
                        unit_exit+=np.longdouble(probability)*np.array([1.,exited,exited**2],dtype=np.longdouble)
                        unit_input+=np.longdouble(probability)*(d*crossing)
            else:
                origin=grid.pivots[cell]+grid.widths[cell]*(measure.mean-grid.offsets[cell])
                scale=grid.widths[cell]*measure.sd
                edges=entropy._partition(measure.lower,measure.upper,measure.A,measure.B)
                # Each record is one positive integral with a known destination.
                low=[];high=[];shifts=[];probabilities=[];destinations=[]
                def add_intervals(a,b,shift,probability,destination=None):
                    if b <= a or probability == 0:return
                    cuts=edges[(edges>a)&(edges<b)]
                    if destination is None:
                        faces=(grid.bounds-origin-shift)/scale
                        cuts=np.r_[cuts,faces[(faces>a)&(faces<b)]]
                    parts=np.unique(np.r_[a,cuts,b])
                    lo=parts[:-1];hi=parts[1:]
                    dest=(np.searchsorted(grid.bounds,origin+shift+scale*.5*(lo+hi),side='right')-1
                          if destination is None else np.full(len(lo),destination,int))
                    low.extend(lo);high.extend(hi);shifts.extend([shift]*len(lo))
                    probabilities.extend([probability]*len(lo));destinations.extend(dest)
                    if len(low)>self.maximum_intervals:raise RuntimeError('Capture overlap interval budget exceeded')
                for count,probability in zip(rule.counts,rule.weights):
                    shift=d*count
                    hi=min(measure.upper,(upper-origin-shift)/scale)
                    add_intervals(measure.lower,hi,shift,probability)
                # A particle exits on capture K when its initial mass is in
                # [upper-K*d,upper-(K-1)*d). K=0 is impossible for a density.
                first=max(1,int(math.floor((upper-grid.bounds[cell+1])/d)))
                last=min(int(rule.counts[-1]),int(math.ceil((upper-grid.bounds[cell])/d))+1)
                for count in range(first,last+1):
                    lo=max(measure.lower,(upper-count*d-origin)/scale)
                    hi=min(measure.upper,(upper-(count-1)*d-origin)/scale)
                    add_intervals(lo,hi,d*count,poisson.sf(count-1,self.rate*dt),grid.pivots.size)
                low=np.asarray(low);high=np.asarray(high);shift=np.asarray(shifts)
                probability=np.asarray(probabilities);destination=np.asarray(destinations)
                # Chunking bounds peak memory without changing the quadrature.
                for start in range(0,len(low),8192):
                    stop=min(start+8192,len(low));sl=slice(start,stop)
                    half=.5*(high[sl]-low[sl]);mid=.5*(high[sl]+low[sl])
                    z=mid[:,None]+half[:,None]*gx
                    weight=half[:,None]*gw*_density(measure,z)*probability[sl,None]
                    product=origin+scale*z+shift[sl,None]
                    dest=destination[sl];inside=dest < grid.pivots.size
                    if (dest < 0).any():raise ValueError('Positive captures cannot cross lower boundary')
                    w=weight.astype(np.longdouble)
                    if inside.any():
                        u=(product[inside]-grid.pivots[dest[inside],None])/grid.widths[dest[inside],None]
                        for power in range(3):np.add.at(unit[:,power],dest[inside],np.sum(w[inside]*u**power,axis=1))
                    if (~inside).any():
                        for power in range(3):unit_exit[power]+=np.sum(w[~inside]*product[~inside]**power)
                    unit_input+=np.sum(w*shift[sl,None])
                    self.stats['positive_quadrature_nodes']+=z.size
            output+=population*unit
            exit_moments+=population*unit_exit;supplied+=population*unit_input
        # Separate analytical Poisson-tail bounds from quadrature/storage errors.
        n,m1,m2=np.asarray(grid.moments(source),dtype=np.longdouble);q0,q1,q2=rule.omitted_count_moments
        omitted=np.array([n*q0,m1*q0+n*d*q1,m2*q0+2*d*m1*q1+n*d*d*q2],dtype=np.longdouble)
        stored=np.asarray(output,float)
        errors=np.abs(output-stored.astype(np.longdouble));c=grid.pivots;h=grid.widths
        storage=np.array([errors[:,0].sum(),np.sum(c*errors[:,0]+h*errors[:,1]),
            np.sum(c*c*errors[:,0]+2*np.abs(c*h)*errors[:,1]+h*h*errors[:,2])])
        grid.distributions(stored)
        return stored,dict(supplied_mass=float(supplied),upper_exit=np.asarray(exit_moments,float)),omitted,storage


def integrate(operator,initial_state,times_s,*,max_step_s=60.,time_limit_s=600.):
    times=np.asarray(times_s,float)
    if (times.ndim!=1 or not len(times) or not np.isfinite(times).all() or (times<0).any()
            or (np.diff(times)<0).any() or not math.isfinite(max_step_s) or max_step_s<=0
            or not math.isfinite(time_limit_s) or time_limit_s<=0):
        raise ValueError('Ordered finite times and positive finite step/time budgets required')
    grid=operator.grid;state=np.asarray(initial_state,float).copy();initial=grid.moments(state)
    ledger=dict(supplied_mass=0.,collision_number_loss=0.,lower_exit=np.zeros(3),upper_exit=np.zeros(3))
    tail=np.zeros(3,dtype=np.longdouble);storage=np.zeros(3,dtype=np.longdouble)
    snapshots=[];t=0.;started=time.monotonic()
    def propagate(bound,dt):
        n,m1,m2=bound.copy();mu=operator.rate*dt;d=operator.donor
        return np.array([n,m1+d*mu*n,m2+2*d*mu*m1+d*d*(mu*mu+mu)*n],dtype=np.longdouble)
    for target in times:
        while t < target:
            if time.monotonic()-started>time_limit_s:
                raise RuntimeError(f'Finite-jump wall time limit reached at model time {t}s')
            dt=min(max_step_s,target-t)
            state,increment,missing,rounded=operator.step(state,dt)
            tail=propagate(tail,dt)+missing
            storage=propagate(storage,dt)+rounded
            ledger['supplied_mass']+=increment['supplied_mass'];ledger['upper_exit']+=increment['upper_exit']
            t=float(target) if dt==target-t else t+dt
        value=grid.moments(state);exits=ledger['upper_exit']
        snapshots.append(dict(time_s=float(target),state=state.copy(),moments=value,variance=grid.variance(state),
            ledgers={k:np.copy(v) if isinstance(v,np.ndarray) else v for k,v in ledger.items()},
            number_residual=float(value[0]+exits[0]-initial[0]),
            mass_residual=float(value[1]+exits[1]-ledger['supplied_mass']-initial[1]),
            poisson_tail_moment_upper_bounds=[float(np.nextafter(float(x),np.inf)) if x>0 else 0. for x in tail],
            storage_rounding_moment_bounds=[float(np.nextafter(float(x),np.inf)) if x>0 else 0. for x in storage]))
    return dict(snapshots=snapshots,diagnostics=dict(operator.stats),wall_seconds=time.monotonic()-started,
        reading_rule='Constant replenished bath; Poisson truncation bounds exclude reconstruction and quadrature errors. Storage bounds are propagated under the exact constant-rate transition, not a bound on nonlinear amplification by reconstruction. No renormalization, population floor, donor depletion or physical snow validation.')
