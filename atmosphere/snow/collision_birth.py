"""Positive finite-time size-dependent capture on actual-mass birth ladders.

Uniformization evolves each source quadrature mass before projecting moments.
A finite birth ladder has an explicit hitting-probability bound; its last state
is not silently admitted as a physical collector. Physical exits are absorbing.
This is a maintained bath, not closed aggregation or a moisture budget.
"""
from __future__ import annotations

import math
import time
import numpy as np
from scipy.stats import poisson
from atmosphere.snow import collision_entropy as entropy
from atmosphere.snow import collision_moments as moments
from atmosphere.snow import collision_jumps as jumps


class StepBudgetError(RuntimeError):
    pass


class BirthCapture:
    def __init__(self,grid,kernel,donor_mass,bath_number,*,quadrature_order=16,
                 maximum_poisson_mean=64.,event_capacity=256,tail_tolerance=1e-24):
        if (not np.isfinite([donor_mass,bath_number,maximum_poisson_mean]).all()
                or donor_mass<=0 or bath_number<=0 or maximum_poisson_mean<=0
                or not 4<=quadrature_order<=64 or int(quadrature_order)!=quadrature_order
                or not 16<=event_capacity<=4096 or int(event_capacity)!=event_capacity
                or not 0<tail_tolerance<.01):
            raise ValueError('Positive bath parameters and supported quadrature/event budgets required')
        self.grid=grid;self.kernel=kernel;self.donor=donor_mass;self.bath=bath_number
        self.order=int(quadrature_order);self.capacity=int(event_capacity)
        self.maximum_mean=maximum_poisson_mean;self.tolerance=tail_tolerance
        self.stats=dict(macrosteps=0,reduced_steps=0,quadrature_sources=0,kernel_evaluations=0,
                        uniformization_iterations=0,maximum_uniformization_mean=0.,
                        maximum_ladder_hitting_probability_bound=0.,
                        maximum_transition_number_residual=0.,minimum_step_s=math.inf,maximum_step_s=0.,
                        maximum_source_relative_moment_error=0.,
                        maximum_source_weight_underflow_number_bound=0.,
                        maximum_storage_moment_rounding_bounds=[0.,0.,0.])

    def _sources(self,state):
        grid=self.grid;gx,gw=moments._legendre(self.order)
        mass=[];weight=[];underflow=np.longdouble(0)
        for cell,measure in enumerate(grid.distributions(state)):
            if measure is None:continue
            number=np.longdouble(state[cell,0])
            if measure.atoms is not None:
                x,w=measure.atoms
                mass.extend(grid.pivots[cell]+grid.widths[cell]*(x-grid.offsets[cell]))
                weight.extend(number*w.astype(np.longdouble));continue
            origin=grid.pivots[cell]+grid.widths[cell]*(measure.mean-grid.offsets[cell])
            scale=grid.widths[cell]*measure.sd
            base=entropy._partition(measure.lower,measure.upper,measure.A,measure.B)
            shifted=(grid.bounds[:,None]-origin-self.donor*np.arange(self.capacity+1))/scale
            cuts=shifted[(shifted>measure.lower)&(shifted<measure.upper)]
            edges=np.unique(np.r_[base,cuts])
            half=.5*np.diff(edges);mid=.5*(edges[1:]+edges[:-1])
            z=mid[:,None]+half[:,None]*gx
            w=(half[:,None]*gw*jumps._density(measure,z)).astype(np.longdouble)
            m=origin+scale*z
            # Record a conservative quadrature-weight underflow bound. It
            # covers underflow in density evaluation and weight multiplication;
            # it is separate from the error of the quadrature formula itself.
            underflow+=number*np.longdouble(np.nextafter(0.,1.))*(measure.upper-measure.lower+w.size)
            keep=w.ravel()>0
            mass.extend(m.ravel()[keep]);weight.extend(number*w.ravel()[keep])
        mass=np.asarray(mass,float);weight=np.asarray(weight,dtype=np.longdouble)
        target=np.asarray(grid.moments(state),dtype=np.longdouble)
        represented=np.array([np.sum(weight*mass**r) for r in range(3)])
        if np.any(target):
            error=float(np.max(abs(represented-target)/np.maximum(abs(target),np.finfo(np.longdouble).tiny)))
            self.stats['maximum_source_relative_moment_error']=max(self.stats['maximum_source_relative_moment_error'],error)
        bound=float(np.nextafter(float(underflow),np.inf)) if underflow>0 else 0.
        self.stats['maximum_source_weight_underflow_number_bound']=max(self.stats['maximum_source_weight_underflow_number_bound'],bound)
        return mass,weight

    def step(self,state,dt):
        if not math.isfinite(dt) or dt<=0:raise ValueError('Positive finite step required')
        grid=self.grid;m,w=self._sources(state);cap=self.capacity;d=self.donor
        if not len(m):return np.zeros_like(state),dict(supplied_mass=0.,upper_exit=np.zeros(3)),np.zeros(3),np.zeros(3)
        output=np.zeros_like(state,dtype=np.longdouble);exit_moments=np.zeros(3,dtype=np.longdouble)
        supplied=np.longdouble(0);truncated=np.zeros(3,dtype=np.longdouble);projection=np.zeros(3,dtype=np.longdouble)
        max_mean=0.;steps=0
        # Chunk the quadrature sources to limit memory; all sources retain the
        # same requested dt, and an over-budget source rejects the whole trial.
        for start in range(0,len(m),256):
            mass=m[start:start+256];number=w[start:start+256]
            counts=np.arange(cap+1)
            ladder=mass[:,None]+d*counts
            crossing=np.array([min(cap,jumps.first_exit_count(x,d,grid.bounds[-1])) for x in mass])
            stopped=np.minimum(counts[None,:],crossing[:,None])
            physical=mass+d*crossing>=grid.bounds[-1]
            rates=self.bath*moments._kernel(self.kernel,ladder[:,:cap],d)
            rates=np.where(counts[None,:cap]>=crossing[:,None],0.,rates)
            nu=rates.max(axis=1);mu=nu*dt;max_mean=max(max_mean,float(mu.max()))
            self.stats['kernel_evaluations']+=rates.size
            if max_mean>self.maximum_mean:raise StepBudgetError('Uniformization work budget requires a smaller macrostep')
            # All uniformized paths up to R jumps lie inside the evaluated
            # ladder. Poisson tails and reaching its artificial cap are distinct.
            R=int(math.ceil(max_mean+8*math.sqrt(max_mean)+16))
            while jumps.poisson_tail_bound(R+1,max_mean,upper=True)>self.tolerance:R+=4
            if R>=cap:raise StepBudgetError('Poisson truncation reaches the finite ladder capacity')
            probability=np.divide(rates,nu[:,None],out=np.zeros_like(rates),where=nu[:,None]>0)
            v=np.zeros((len(mass),R+1));v[:,0]=1.
            coefficient=np.exp(-mu);transition=np.zeros_like(v);transition[:,0]=coefficient
            for r in range(1,R+1):
                # Update right to left algebraically using the old array.
                previous=v[:,:r].copy()
                v[:,1:r+1]=v[:,1:r+1]*(1-probability[:,1:r+1])+previous*probability[:,:r]
                v[:,0]*=1-probability[:,0]
                coefficient*=mu/r
                transition[:,:r+1]+=coefficient[:,None]*v[:,:r+1]
            actual_tail=poisson.sf(R,mu)
            self.stats['maximum_transition_number_residual']=max(self.stats['maximum_transition_number_residual'],float(np.max(abs(transition.sum(axis=1)+actual_tail-1))))
            k=stopped[:,:R+1];products=mass[:,None]+d*k
            population=number[:,None]*transition.astype(np.longdouble)
            dest=np.searchsorted(grid.bounds,products,side='right')-1
            inside=dest<grid.pivots.size
            index=dest[inside];u=(products[inside]-grid.pivots[index])/grid.widths[index]
            for power in range(3):np.add.at(output[:,power],index,population[inside]*u**power)
            outside=~inside
            for power in range(3):exit_moments[power]+=np.sum(population[outside]*products[outside]**power)
            supplied+=np.sum(population*(d*k))
            for i in range(len(mass)):
                omitted=jumps.poisson_tail_bound(R+1,mu[i],upper=True)
                hit=(np.longdouble(0) if physical[i] else jumps.poisson_tail_bound(cap,mu[i],upper=True))
                finite_max=mass[i]+d*crossing[i]
                physical_max=grid.bounds[-1]+d
                truncated+=number[i]*omitted*np.array([1,finite_max,finite_max**2],dtype=np.longdouble)
                projection+=number[i]*hit*np.array([1,physical_max,physical_max**2],dtype=np.longdouble)
                self.stats['maximum_ladder_hitting_probability_bound']=max(self.stats['maximum_ladder_hitting_probability_bound'],float(hit))
            steps+=R
        stored=np.asarray(output,float);grid.distributions(stored)
        error=np.abs(output-stored.astype(np.longdouble));c=grid.pivots;h=grid.widths
        bounds=[np.sum(error[:,0]),np.sum(c*error[:,0]+h*error[:,1]),
                np.sum(c*c*error[:,0]+2*np.abs(c*h)*error[:,1]+h*h*error[:,2])]
        self.stats['maximum_storage_moment_rounding_bounds']=[max(old,float(np.nextafter(float(x),np.inf)) if x>0 else 0.) for old,x in zip(self.stats['maximum_storage_moment_rounding_bounds'],bounds)]
        self.stats['macrosteps']+=1;self.stats['quadrature_sources']+=len(m)
        self.stats['uniformization_iterations']+=steps
        self.stats['maximum_uniformization_mean']=max(self.stats['maximum_uniformization_mean'],max_mean)
        self.stats['minimum_step_s']=min(self.stats['minimum_step_s'],dt)
        self.stats['maximum_step_s']=max(self.stats['maximum_step_s'],dt)
        return stored,dict(supplied_mass=float(supplied),upper_exit=np.asarray(exit_moments,float)),truncated,projection


def integrate(operator,initial_state,times_s,*,max_step_s=1.,time_limit_s=600.):
    times=np.asarray(times_s,float)
    if (times.ndim!=1 or not len(times) or not np.isfinite(times).all() or (times<0).any()
            or (np.diff(times)<0).any() or not math.isfinite(max_step_s) or max_step_s<=0
            or not math.isfinite(time_limit_s) or time_limit_s<=0):
        raise ValueError('Ordered finite times and positive finite budgets required')
    grid=operator.grid;state=np.asarray(initial_state,float).copy();initial=grid.moments(state)
    ledger=dict(supplied_mass=0.,collision_number_loss=0.,lower_exit=np.zeros(3),upper_exit=np.zeros(3))
    error_number=np.longdouble(0);projection_number=np.longdouble(0)
    # After any later capture, stopped masses cannot exceed this physical exit
    # bound. This propagation bounds omitted transition probability, not the
    # amplification of reconstruction or quadrature errors.
    cap=grid.bounds[-1]+operator.donor
    snapshots=[];t=0.;suggested=max_step_s;started=time.monotonic()
    for target in times:
        while t<target:
            if time.monotonic()-started>time_limit_s:raise RuntimeError(f'Birth-update wall time limit reached at model time {t}s after {operator.stats["macrosteps"]} macrosteps')
            dt=min(suggested,max_step_s,target-t)
            while True:
                if t+dt==t:raise FloatingPointError('Macrostep below time resolution')
                try:next_state,increment,missing,projection=operator.step(state,dt);break
                except StepBudgetError:dt*=.5;operator.stats['reduced_steps']+=1
            state=next_state;error_number+=missing[0];projection_number+=projection[0]
            for key in increment:ledger[key]+=increment[key]
            t=float(target) if dt==target-t else t+dt;suggested=min(max_step_s,dt*1.25)
        value=grid.moments(state);exits=ledger['upper_exit']
        snapshots.append(dict(time_s=float(target),state=state.copy(),moments=value,variance=grid.variance(state),
            ledgers={k:np.copy(v) if isinstance(v,np.ndarray) else v for k,v in ledger.items()},
            number_residual=float(value[0]+exits[0]-initial[0]),mass_residual=float(value[1]+exits[1]-ledger['supplied_mass']-initial[1]),
            uniformization_tail_moment_upper_bounds=[float(np.nextafter(float(error_number*cap**r),np.inf)) if error_number else 0. for r in range(3)],
            finite_ladder_moment_upper_bounds=[float(np.nextafter(float(projection_number*cap**r),np.inf)) if projection_number else 0. for r in range(3)]))
    diagnostics=dict(operator.stats)
    if not diagnostics['macrosteps']:diagnostics['minimum_step_s']=0.
    return dict(snapshots=snapshots,diagnostics=diagnostics,wall_seconds=time.monotonic()-started,
        reading_rule='State-dependent birth rates on actual masses; replenished bath. Uniformization and finite-ladder probability/moment bounds exclude source reconstruction, quadrature and floating-point errors. Bounds propagated using the absorbing physical-domain exit mass. Closed aggregation and moisture depletion remain unimplemented.')
