"""Bounded quadratic maximum-entropy reconstruction for collision moments.

This independent reconstruction candidate retains the frozen three-moment event
bookkeeping and time integration. The reference measure is uniform particle
number per unit mass. Positive quadratic curvature is allowed; exact atomic
limits need no finite-density fit. No physical kernel or mass floor is changed.
"""
from __future__ import annotations

from collections import OrderedDict, namedtuple
import math
import time

import numpy as np
from scipy.special import log_ndtr, dawsn

from atmosphere.snow import collision_moments as base

BathCapture = base.BathCapture
CollisionMoments = base.CollisionMoments
RealizabilityError = base.RealizabilityError
QuadratureError = base.QuadratureError


def _normal_log_probability(lower, upper):
    lower,upper = np.broadcast_arrays(np.asarray(lower,float),np.asarray(upper,float))
    result = np.full_like(lower,-np.inf)
    valid = upper > lower
    left = valid & (upper <= 0)
    right = valid & (lower >= 0)
    middle = valid & ~(left|right)
    with np.errstate(divide='ignore',invalid='ignore'):
        result[left] = log_ndtr(upper[left])+np.log(-np.expm1(log_ndtr(lower[left])-log_ndtr(upper[left])))
        result[right] = log_ndtr(-lower[right])+np.log(-np.expm1(log_ndtr(-upper[right])-log_ndtr(-lower[right])))
        result[middle] = np.log(-np.expm1(log_ndtr(lower[middle]))-np.exp(log_ndtr(-upper[middle])))
    return result


def _partition(lower,upper,A,B):
    """Resolve interior and endpoint peaks without deleting any interval."""
    width = upper-lower
    points = [lower,upper]
    curvature_scale = math.sqrt(abs(B))
    for endpoint,direction in ((lower,1),(upper,-1)):
        scale = 1/max(abs(A+2*B*endpoint),curvature_scale,1/width)
        points.extend(endpoint+direction*scale*k for k in (.25,1,2,4,8,16,32,64))
    if B < 0:
        mode = -A/(2*B)
        scale = math.sqrt(-.5/B)
        points.extend(mode+scale*k for k in (-32,-16,-8,-4,-2,-1,0,1,2,4,8,16,32))
    return np.array(sorted(set(p for p in points if lower <= p <= upper)))


def _positive_rule(lower,upper,A,B,order=24):
    """Positive integral rule; exponent differences use the peak as origin."""
    edges = _partition(lower,upper,A,B)
    half = .5*np.diff(edges)
    midpoint = .5*(edges[1:]+edges[:-1])
    gx,gw = base._legendre(order)
    nodes = (midpoint[:,None]+half[:,None]*gx).ravel()
    weights = (half[:,None]*gw).ravel()
    if B < 0 and lower <= -A/(2*B) <= upper:
        peak = -A/(2*B)
        relative = B*(nodes-peak)**2
        peak_value = -.25*A*A/B
    else:
        difference = (upper-lower)*(A+B*(upper+lower))
        peak = upper if difference >= 0 else lower
        relative = (nodes-peak)*(A+B*(nodes+peak))
        peak_value = peak*(A+B*peak)
    if not np.isfinite(relative).all() or not math.isfinite(peak_value):
        raise QuadratureError('Quadratic entropy exponent exceeds finite precision')
    weights *= np.exp(relative)
    total = float(weights.sum())
    if total <= 0 or not math.isfinite(total):
        raise QuadratureError('Quadratic entropy normalization failed')
    weights /= total
    return nodes,weights,peak_value+math.log(total)


def _statistics(lower,upper,A,B,order=24):
    nodes,weights,logZ = _positive_rule(lower,upper,A,B,order)
    mean = float(weights@nodes)
    second = float(weights@nodes**2)
    centered1 = nodes-mean
    centered2 = nodes**2-second
    hessian = np.array([[weights@centered1**2,weights@(centered1*centered2)],
                        [weights@(centered1*centered2),weights@centered2**2]])
    return np.array([mean,second-1]),hessian,logZ-B,nodes,weights,logZ


def _fit_from_initial(mean,variance,initial=None):
    """Fit the dual, then independently admit the reconstructed moments.

    The optimizer target is tighter than the quadrature admission criterion.
    A stalled optimizer may still supply an admissible reconstruction; its
    termination status and independently refined residual remain explicit.
    """
    sd = math.sqrt(variance)
    lower,upper = -mean/sd,(1-mean)/sd
    parameters = np.array([0.,-.5] if initial is None else initial,dtype=float)
    maximum_residual = math.inf
    iterations = 0
    status = 'iteration_limit'
    for iteration in range(80):
        gradient,hessian,dual,nodes,weights,logZ = _statistics(lower,upper,*parameters)
        maximum_residual = float(np.max(abs(gradient)))
        iterations = iteration+1
        if maximum_residual <= 2e-11:
            status = 'optimizer_target'
            break
        try:
            step = np.linalg.solve(hessian,-gradient)
        except np.linalg.LinAlgError:
            status = 'singular_hessian'
            break
        if not np.isfinite(step).all():
            status = 'nonfinite_direction'
            break
        directional = float(gradient@step)
        if directional >= 0:
            status = 'non_descent_direction'
            break
        accepted = False
        for halves in range(50):
            fraction = .5**halves
            trial = parameters+fraction*step
            try:
                trial_g,_,trial_dual,*_ = _statistics(lower,upper,*trial)
            except (QuadratureError,FloatingPointError,OverflowError):
                continue
            # The second condition admits residual improvement once dual
            # differences themselves approach floating-point roundoff.
            if (trial_dual <= dual+1e-4*fraction*directional
                    or (abs(trial_dual-dual) <= 64*base.EPS*max(1.,abs(dual))
                        and np.max(abs(trial_g)) < maximum_residual)):
                parameters = trial
                accepted = True
                break
        if not accepted:
            status = 'line_search_stall'
            break
    # Recompute both rules at the final parameters, including a step taken on
    # the last optimizer iteration. No stale normalization enters admission.
    coarse = _statistics(lower,upper,*parameters)
    refined = _statistics(lower,upper,*parameters,order=48)
    residual = float(np.max(abs(refined[0])))
    normalization_difference = float(abs(refined[-1]-coarse[-1]))
    if (not math.isfinite(residual) or not math.isfinite(normalization_difference)
            or residual > 2e-10 or normalization_difference > 2e-10):
        raise QuadratureError(
            f'Entropy inverse moments fail quadrature refinement after {status}: '
            f'residual {residual:.4g}, log normalization difference {normalization_difference:.4g}')
    return (float(parameters[0]),float(parameters[1]),float(refined[-1]),iterations,
            residual,status,float(np.max(abs(coarse[0]))),normalization_difference)


_FIT_CACHE = OrderedDict()
_FIT_CACHE_COUNTS = dict(hits=0,misses=0)
_FitCacheInfo = namedtuple('FitCacheInfo','hits misses maxsize currsize')


def _fit(mean,variance,initial=None):
    """Cache by exact target moments, independently of an optional seed."""
    key=(mean,variance)
    if key in _FIT_CACHE:
        _FIT_CACHE_COUNTS['hits'] += 1
        _FIT_CACHE.move_to_end(key)
        return _FIT_CACHE[key]
    _FIT_CACHE_COUNTS['misses'] += 1
    warm=initial is not None
    fallback=False
    if warm:
        try:
            if len(initial) != 2 or not np.isfinite(initial).all():
                raise QuadratureError('Invalid warm-start parameters')
            result=_fit_from_initial(mean,variance,initial)
        except (ArithmeticError,ValueError):
            fallback=True
            result=_fit_from_initial(mean,variance)
    else:
        result=_fit_from_initial(mean,variance)
    result=(*result,warm,fallback)
    _FIT_CACHE[key]=result
    if len(_FIT_CACHE) > 8192:
        _FIT_CACHE.popitem(last=False)
    return result


def _fit_cache_clear():
    _FIT_CACHE.clear()
    _FIT_CACHE_COUNTS.update(hits=0,misses=0)


def _fit_cache_info():
    return _FitCacheInfo(**_FIT_CACHE_COUNTS,maxsize=8192,currsize=len(_FIT_CACHE))


_fit.cache_clear=_fit_cache_clear
_fit.cache_info=_fit_cache_info


def _gaussian_from_positive_samples(nodes,weights,order,lower,upper):
    """Compress a positive rule by orthogonalizing sampled polynomials.

    A reorthogonalized Lanczos basis avoids a nearly singular raw-moment
    matrix for densities concentrated near both endpoints. The compressed
    multiplication operator is symmetric; its spectrum lies in the sampled
    support. The same moment and support admission gates remain explicit.
    """
    nodes,weights=np.asarray(nodes,float),np.asarray(weights,float)
    if (nodes.ndim != 1 or nodes.shape != weights.shape or nodes.size == 0
            or not np.isfinite(nodes).all() or not np.isfinite(weights).all()
            or (weights < 0).any() or weights.sum() <= 0):
        raise QuadratureError('Finite positive sampled measure required')
    total=float(weights.sum())
    first=np.sqrt(weights/total)
    basis=[first]
    for _ in range(1,order):
        residual=nodes*basis[-1]
        # Two orthogonalization passes retain distinct narrow endpoint modes
        # even when the next recurrence coefficient is very small.
        for _ in range(2):
            for vector in basis:
                residual -= (vector@residual)*vector
        norm=float(np.linalg.norm(residual))
        if norm == 0:
            break  # Exact finite sampled support has exhausted its rank.
        if not math.isfinite(norm):
            raise QuadratureError('Nonfinite sampled recurrence coefficient')
        basis.append(residual/norm)
    vectors=np.column_stack(basis)
    projected=vectors.T@(nodes[:,None]*vectors)
    projected=.5*(projected+projected.T)
    try:
        q,eigenvectors=np.linalg.eigh(projected)
    except np.linalg.LinAlgError as exc:
        raise QuadratureError('Sampled quadrature eigensystem failed') from exc
    w=total*eigenvectors[0]**2
    tolerance=256*base.EPS*max(1.,abs(lower),abs(upper))
    if (not np.isfinite(q).all() or not np.isfinite(w).all() or (w <= 0).any()
            or q.min() < lower-tolerance or q.max() > upper+tolerance):
        raise QuadratureError('Sampled quadrature leaves its support')
    expected=np.array([weights@nodes**k for k in range(2*order)])
    reproduced=np.array([w@q**k for k in range(2*order)])
    if np.max(abs(reproduced-expected)/np.maximum(1.,abs(expected))) > 2e-9:
        raise QuadratureError('Sampled quadrature does not reproduce moments')
    return q,w


class EntropyMeasure(base.BetaMeasure):
    """Maximum entropy on [0,1] with prescribed mean and variance."""
    def __init__(self,mean,variance,adjustment=0.,initial=None):
        super().__init__(mean,variance,adjustment)
        self.sd = math.sqrt(variance)
        if self.atoms is None:
            try:
                (self.A,self.B,self.logZ,self.fit_iterations,self.fit_residual,
                 self.fit_status,self.fit_optimizer_residual,
                 self.fit_normalization_difference,self.fit_warm_start,
                 self.fit_cold_fallback) = _fit(mean,variance,initial)
            except QuadratureError as exc:
                raise QuadratureError(f'{exc}; mean={mean:.17g}, variance={variance:.17g}') from exc
            self.lower,self.upper = -mean/self.sd,(1-mean)/self.sd
        else:
            self.A=self.B=self.logZ=0.
            self.fit_iterations=0
            self.fit_residual=0.
            self.fit_status='atomic'
            self.fit_optimizer_residual=0.
            self.fit_normalization_difference=0.
            self.fit_warm_start=False
            self.fit_cold_fallback=False

    def _partial(self,x,need_first=True):
        """Normalized cumulative number and first z moment."""
        query = np.asarray(x,float)
        z = (np.clip(query,0.,1.)-self.mean)/self.sd
        L,U,A,B = self.lower,self.upper,self.A,self.B
        if B < -1e-9:
            location = -A/(2*B)
            scale = math.sqrt(-.5/B)
            lo,hi,cut = (L-location)/scale,(U-location)/scale,(z-location)/scale
            log_total = _normal_log_probability(lo,hi)
            probability = np.exp(_normal_log_probability(lo,cut)-log_total)
            first = np.zeros_like(z)
            if need_first and not (lo > 12 or hi < -12):
                logleft = -.5*lo*lo-.5*math.log(2*math.pi)-log_total
                logcut = -.5*cut*cut-.5*math.log(2*math.pi)-log_total
                maximum = np.maximum(logleft,logcut)
                difference = np.exp(maximum)*(np.exp(logleft-maximum)-np.exp(logcut-maximum))
                first = location*probability+scale*difference
            elif need_first:
                for index in np.ndindex(z.shape):
                    if z[index] > L:
                        nodes,weights,log_part = _positive_rule(L,float(z[index]),A,B,32)
                        first[index] = math.exp(log_part-self.logZ)*(weights@nodes)
        elif B > 1e-9:
            root = math.sqrt(B)
            qlo,qcut = (A+2*B*L)/(2*root),(A+2*B*z)/(2*root)
            flo = np.exp(L*(A+B*L)-self.logZ)
            fcut = np.exp(z*(A+B*z)-self.logZ)
            probability = (fcut*dawsn(qcut)-flo*dawsn(qlo))/root
            first = (fcut-flo-A*probability)/(2*B) if need_first else np.zeros_like(z)
            if need_first and abs(A)/root > 24:
                for index in np.ndindex(z.shape):
                    if z[index] > L:
                        nodes,weights,log_part = _positive_rule(L,float(z[index]),A,B,32)
                        first[index] = math.exp(log_part-self.logZ)*(weights@nodes)
        else:
            # Near-linear curvature avoids subtracting two O(1/B) terms.
            # Positive Gauss integration over each query remains inexpensive
            # for this uncommon branch; no distribution is replaced by one
            # with different parameters.
            probability = np.zeros_like(z)
            first = np.zeros_like(z)
            for index in np.ndindex(z.shape):
                if z[index] > L:
                    nodes,weights,log_part = _positive_rule(L,float(z[index]),A,B,32)
                    probability[index] = math.exp(log_part-self.logZ)
                    first[index] = probability[index]*(weights@nodes)
        probability = np.where(query <= 0,0.,np.where(query >= 1,1.,probability))
        first = np.where((query <= 0)|(query >= 1),0.,first)
        cumulative_mass = self.mean*probability+self.sd*first
        interior = (query > 0)&(query < 1)
        invalid = interior & ((probability < 0)|(probability > 1)|~np.isfinite(probability))
        if need_first:
            tolerance = 128*base.EPS*np.maximum(1.,abs(cumulative_mass))
            invalid |= interior & ((cumulative_mass < -tolerance)
                                   |(cumulative_mass > query*probability+tolerance)
                                   |~np.isfinite(first))
        for index in np.ndindex(z.shape):
            if invalid[index]:
                nodes,weights,log_part = _positive_rule(L,float(z[index]),A,B,48)
                probability[index] = math.exp(log_part-self.logZ)
                first[index] = probability[index]*(weights@nodes)
        if not np.isfinite(probability).all() or not np.isfinite(first).all():
            raise QuadratureError('Nonfinite entropy cumulative moments')
        cumulative_mass = self.mean*probability+self.sd*first
        if (np.min(probability) < -1e-12 or np.max(probability) > 1+1e-12
                or (need_first and (np.min(cumulative_mass) < -1e-12
                    or np.any(cumulative_mass > np.clip(query,0.,1.)*probability+1e-12)))):
            raise QuadratureError('Entropy cumulative moments leave their support bounds')
        return probability,first

    def cdf(self,x):
        if self.atoms is not None:
            return super().cdf(x)
        return self._partial(x,need_first=False)[0]

    def partial_first(self,x):
        if self.atoms is not None:
            return super().partial_first(x)
        probability,first = self._partial(x)
        return self.mean*probability+self.sd*first

    def quadrature(self,order=2,lower=0.,upper=1.,*,stats=None):
        if not 2 <= order <= 6 or not 0 <= lower <= upper <= 1:
            raise ValueError('Ordered unit bounds and quadrature order 2..6 required')
        if lower == upper:
            return np.empty(0),np.empty(0)
        if self.atoms is not None:
            return super().quadrature(order,lower,upper,stats=stats)
        lo,hi = (lower-self.mean)/self.sd,(upper-self.mean)/self.sd
        nodes,weights,log_part = _positive_rule(lo,hi,self.A,self.B,24)
        probability = math.exp(log_part-self.logZ)
        if probability == 0:
            return np.empty(0),np.empty(0)
        center = float(weights@nodes)
        variance = float(weights@(nodes-center)**2)
        if variance <= 0:
            raise QuadratureError('Conditional entropy variance is unresolved')
        sd = math.sqrt(variance)
        centered = (nodes-center)/sd
        q,w = _gaussian_from_positive_samples(
            centered,weights,order,(lo-center)/sd,(hi-center)/sd)
        return self.mean+self.sd*(center+sd*q),probability*w


class MomentGrid(base.MomentGrid):
    """Frozen local moments, with an explicit subnormal-storage error bound.

    State is never modified. The extra cone tolerance propagates one-half of
    the smallest float spacing through N,S1,S2 normalization. It applies only
    where floating storage cannot resolve a moment defect. Reconstructed
    population and physical moment bounds are recorded, including trial stages.
    """
    def __init__(self,pivots):
        super().__init__(pivots)
        self._previous_measure = [None]*self.pivots.size
        self.fit_diagnostics = dict(
            validating_calls=0,
            reconstructed_finite_density_cells=0,
            reconstructed_atomic_cells=0,
            admitted_optimizer_stalls=0,
            warm_start_reconstructions=0,
            cold_start_fallback_reconstructions=0,
            maximum_independently_refined_residual=0.,
            maximum_normalization_log_difference=0.,
            maximum_optimizer_iterations=0,
            optimizer_termination_counts={},
            optimizer_target=2e-11,
            independent_admission_tolerance=2e-10,
            reading_rule='Residuals are standardized mean and second-moment errors under independently refined positive quadrature. Counts include cache reuse, repeated validations and rejected trial stages, not independent fits or physical events. Warm-start and fallback counts describe the cached fit used for each reconstruction. Optimizer stalls are admitted only by the unchanged moment and normalization-refinement criterion.')
        self.roundoff_diagnostics = dict(
            validating_calls_with_storage_reconstruction=0,
            storage_reconstructed_cells=0,
            maximum_affected_number=0.,
            maximum_affected_mass_upper_bound=0.,
            maximum_affected_second_moment_upper_bound=0.,
            maximum_absolute_reconstruction_first_moment_change=0.,
            maximum_absolute_reconstruction_second_moment_change=0.,
            maximum_storage_number_rounding_bound=0.,
            maximum_storage_first_moment_rounding_bound=0.,
            maximum_storage_second_moment_rounding_bound=0.,
            reading_rule='Bounds use the number and mass units of this grid; counts include repeated validations and rejected trial stages, not independent physical events. Stored populations and moments are unchanged. Positive bounds are rounded upward when converted to float.')

    @staticmethod
    def _upper_float(value):
        result = float(value)
        if np.longdouble(result) < value:
            result = float(np.nextafter(result,np.inf))
        return result

    def _storage_admissible(self,state):
        state = np.asarray(state,float)
        if state.shape != (self.pivots.size,3) or not np.isfinite(state).all():
            raise RealizabilityError('Finite state with shape (cells,3) required')
        if (state[:,0] < 0).any():
            raise RealizabilityError('Negative population')
        half_spacing = np.longdouble(np.nextafter(0.,1.))/2
        affected = []
        result = []
        try:
            for i,(number,first,second) in enumerate(state):
                if number == 0:
                    if first != 0 or second != 0:
                        raise RealizabilityError('Empty cell has nonzero moments')
                    result.append(None)
                    continue
                n = np.longdouble(number)
                u,q = np.longdouble(first)/n,np.longdouble(second)/n
                mean = u+np.longdouble(self.offsets[i])
                variance = q-u*u
                mean_roundoff = 64*base.EPS*(abs(u)+abs(self.offsets[i]))
                variance_roundoff = 64*base.EPS*(abs(q)+u*u)
                denominator = n-half_spacing
                mean_storage = half_spacing*(1+abs(u))/denominator
                second_storage = half_spacing*(1+abs(q))/denominator
                variance_storage = second_storage+2*abs(u)*mean_storage+mean_storage**2
                mean_tolerance = mean_roundoff+mean_storage
                variance_tolerance = variance_roundoff+variance_storage
                if mean < -mean_tolerance or mean > 1+mean_tolerance:
                    raise RealizabilityError(f'Cell {i} mean leaves support beyond storage uncertainty')
                bounded_mean = min(np.longdouble(1),max(np.longdouble(0),mean))
                maximum = bounded_mean*(1-bounded_mean)
                if variance < -variance_tolerance or variance > maximum+variance_tolerance+mean_tolerance:
                    raise RealizabilityError(f'Cell {i} variance leaves support beyond storage uncertainty')
                bounded_variance = min(maximum,max(np.longdouble(0),variance))
                # Preserve the frozen normal-range roundoff convention. A
                # feasible positive variance is never erased merely because
                # the larger storage uncertainty could also contain zero.
                if abs(variance) <= variance_roundoff:
                    bounded_variance = np.longdouble(0)
                elif abs(variance-maximum) <= variance_roundoff:
                    bounded_variance = maximum
                adjustment = max(abs(mean-bounded_mean),abs(variance-bounded_variance))
                uses_storage = (abs(mean-bounded_mean) > mean_roundoff
                                or abs(variance-bounded_variance) > variance_roundoff)
                if uses_storage:
                    unew = bounded_mean-np.longdouble(self.offsets[i])
                    qnew = unew*unew+bounded_variance
                    c,h = np.longdouble(self.pivots[i]),np.longdouble(self.widths[i])
                    upper = np.longdouble(self.bounds[i+1])
                    affected.append(dict(
                        maximum_affected_number=n,
                        maximum_affected_mass_upper_bound=n*upper,
                        maximum_affected_second_moment_upper_bound=n*upper*upper,
                        maximum_absolute_reconstruction_first_moment_change=abs(n*h*(unew-u)),
                        maximum_absolute_reconstruction_second_moment_change=abs(n*(2*c*h*(unew-u)+h*h*(qnew-q))),
                        maximum_storage_number_rounding_bound=half_spacing,
                        maximum_storage_first_moment_rounding_bound=(c+h)*half_spacing,
                        maximum_storage_second_moment_rounding_bound=(c*c+2*abs(c*h)+h*h)*half_spacing))
                mean_float = float(bounded_mean)
                variance_float = (mean_float*(1-mean_float) if bounded_variance == maximum
                                  else float(bounded_variance))
                result.append(base.BetaMeasure(mean_float,variance_float,float(adjustment)))
            return result
        finally:
            if affected:
                self.roundoff_diagnostics['validating_calls_with_storage_reconstruction'] += 1
                self.roundoff_diagnostics['storage_reconstructed_cells'] += len(affected)
                for key in affected[0]:
                    bound = sum((r[key] for r in affected),np.longdouble(0))
                    self.roundoff_diagnostics[key] = max(self.roundoff_diagnostics[key],self._upper_float(bound))

    def distributions(self,state):
        try:
            admissible = super().distributions(state)
        except RealizabilityError:
            admissible = self._storage_admissible(state)
        self.fit_diagnostics['validating_calls'] += 1
        reconstructed = []
        for cell,d in enumerate(admissible):
            if d is None:
                reconstructed.append(None)
                continue
            previous=self._previous_measure[cell]
            initial=None
            if previous is not None and d.atoms is None:
                # z_old=a+b*z_new transforms the same physical density into
                # the new target's centered/scaled coordinate system.
                a=(d.mean-previous.mean)/previous.sd
                b=math.sqrt(d.variance)/previous.sd
                initial=(b*(previous.A+2*previous.B*a),previous.B*b*b)
            measure = EntropyMeasure(d.mean,d.variance,d.adjustment,initial=initial)
            self._previous_measure[cell]=measure if measure.atoms is None else None
            reconstructed.append(measure)
            stats = self.fit_diagnostics
            if measure.atoms is not None:
                stats['reconstructed_atomic_cells'] += 1
                continue
            stats['reconstructed_finite_density_cells'] += 1
            stats['warm_start_reconstructions'] += int(measure.fit_warm_start)
            stats['cold_start_fallback_reconstructions'] += int(measure.fit_cold_fallback)
            stats['maximum_independently_refined_residual'] = max(
                stats['maximum_independently_refined_residual'],measure.fit_residual)
            stats['maximum_normalization_log_difference'] = max(
                stats['maximum_normalization_log_difference'],measure.fit_normalization_difference)
            stats['maximum_optimizer_iterations'] = max(
                stats['maximum_optimizer_iterations'],measure.fit_iterations)
            status = measure.fit_status
            stats['optimizer_termination_counts'][status] = stats['optimizer_termination_counts'].get(status,0)+1
            if status != 'optimizer_target':
                stats['admitted_optimizer_stalls'] += 1
        return reconstructed


def integrate(operator, initial_state, times_s, *, max_step_s=30., cfl=.3, maximum_steps=200000,
              time_limit_s=None, progress=None):
    """SSPRK2 with stable averaging and the frozen trial-stage/ledger sequence.

    Both Euler stages must be realizable; no population clipping occurs.
    Roundoff-sized moment-cone reconstruction adjustments leave conserved state
    untouched and are reported. The net number CFL is supplementary; moment
    realizability also limits steps in a high-rate bath with constant number.
    """
    times = np.asarray(times_s, float)
    if (times.ndim != 1 or not times.size or not np.isfinite(times).all()
            or (times < 0).any() or (np.diff(times) < 0).any()
            or not math.isfinite(max_step_s) or max_step_s <= 0 or not 0 < cfl <= 1):
        raise ValueError('Ordered finite times, positive step and CFL in (0,1] required')
    grid, state = operator.grid, np.asarray(initial_state, float).copy()
    initial = grid.moments(state)
    ledger = dict(supplied_mass=0., collision_number_loss=0.,
                  lower_exit=np.zeros(3), upper_exit=np.zeros(3))
    stats = dict(steps=0, rejected_steps=0, minimum_step_s=math.inf,
                 maximum_step_s=0., maximum_roundoff_reconstruction_adjustment=0.)
    snapshots, t, suggested = [], 0., max_step_s
    started = time.monotonic()
    if time_limit_s is not None and (not math.isfinite(time_limit_s) or time_limit_s <= 0):
        raise ValueError('Positive finite wall time limit required')
    def validate(value):
        distributions = grid.distributions(value)
        adjustment = max((d.adjustment for d in distributions if d is not None), default=0.)
        stats['maximum_roundoff_reconstruction_adjustment'] = max(
            stats['maximum_roundoff_reconstruction_adjustment'], adjustment)
    validate(state)
    for target in times:
        while t < target:
            if stats['steps'] >= maximum_steps:
                raise RuntimeError(f'Maximum step budget exhausted at model time {t:.17g}s')
            if time_limit_s is not None and time.monotonic()-started >= time_limit_s:
                raise RuntimeError(f'Wall time limit reached at model time {t:.17g}s after {stats["steps"]} steps')
            rhs0, rate0 = operator.evaluate(state)
            losing = rhs0[:,0] < 0
            # An unrepresentably long decay time is an infinite step bound.
            with np.errstate(over='ignore'):
                limit = float(np.min(-state[losing,0]/rhs0[losing,0])) if losing.any() else math.inf
            dt = min(suggested, max_step_s, target-t, cfl*limit)
            accepted = False
            for _ in range(80):
                if t+dt == t or dt <= 0:
                    raise FloatingPointError('Realizable step fell below floating-point time resolution')
                stage = state+dt*rhs0
                try:
                    validate(stage)
                    rhs1, rate1 = operator.evaluate(stage)
                    second = stage+dt*rhs1
                    validate(second)
                    # Two separately rounded subnormal halves can change a
                    # stationary state. This equivalent convex average keeps
                    # identical inputs identical at every storage magnitude.
                    candidate = state+.5*(second-state)
                    validate(candidate)
                    accepted = True
                    break
                except RealizabilityError:
                    stats['rejected_steps'] += 1
                    dt *= .5
            if not accepted:
                raise FloatingPointError('Unable to construct realizable SSP stages')
            state = candidate
            for key in ledger:
                ledger[key] += .5*dt*(rate0[key]+rate1[key])
            t = float(target) if dt == target-t else t+dt
            suggested = min(max_step_s, dt*1.6)
            stats['steps'] += 1
            stats['minimum_step_s'] = min(stats['minimum_step_s'], dt)
            stats['maximum_step_s'] = max(stats['maximum_step_s'], dt)
            if progress is not None:
                progress(t, dict(stats))
        moments = grid.moments(state)
        exits = ledger['lower_exit']+ledger['upper_exit']
        snapshots.append(dict(time_s=float(target), state=state.copy(), moments=moments,
                              variance=grid.variance(state),
                              ledgers={k: np.copy(v) if isinstance(v,np.ndarray) else v for k,v in ledger.items()},
                              number_residual=float(moments[0]+exits[0]+ledger['collision_number_loss']-initial[0]),
                              mass_residual=float(moments[1]+exits[1]-ledger['supplied_mass']-initial[1])))
    stats.update(operator.stats)
    stats['ssp_average'] = 'state + 0.5*(second_euler_stage-state)'
    if not stats['steps']:
        stats['minimum_step_s'] = 0.
    return dict(snapshots=snapshots, diagnostics=stats)
