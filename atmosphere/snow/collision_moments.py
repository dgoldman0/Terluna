"""Bounded local moment closure for aggregation, without mass-pivot projection.

Every cell stores N, integral(u n), integral(u**2 n), u=(m-pivot)/width.
A bounded Beta distribution supplies unresolved shape. Collisions use actual
quadrature masses and route exact product moments into destination cells. This
is an aggregation-only closure experiment: no deposition, chipping or column
transport is included, and matching three moments does not establish tail shape.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math
import time

import numpy as np
from scipy.integrate import quad_vec
from scipy.linalg import eigh, eigh_tridiagonal
from scipy.special import betainc, betaincc, betaln

EPS = np.finfo(float).eps


@lru_cache(maxsize=8)
def _legendre(order):
    return np.polynomial.legendre.leggauss(order)



class RealizabilityError(ValueError):
    """A state is outside the bounded nonnegative moment cone."""


class QuadratureError(FloatingPointError):
    """A positive supported quadrature could not be established."""


def _gaussian_from_moments(moment, order, lower, upper):
    """Gaussian quadrature in an already centered/scaled coordinate."""
    moment = np.asarray(moment, dtype=float)
    h0 = np.array([[moment[i+j] for j in range(order)] for i in range(order)])
    h1 = np.array([[moment[i+j+1] for j in range(order)] for i in range(order)])
    try:
        nodes, vectors = eigh(h1, h0, check_finite=False)
    except np.linalg.LinAlgError as exc:
        raise QuadratureError('Conditional moment matrix is not positive') from exc
    weights = (moment[:order] @ vectors)**2
    tolerance = 256*EPS*max(1., abs(lower), abs(upper))
    if (not np.isfinite(nodes).all() or not np.isfinite(weights).all()
            or (weights <= 0).any() or nodes.min() < lower-tolerance
            or nodes.max() > upper+tolerance):
        raise QuadratureError('Conditional quadrature leaves its support')
    reproduced = np.array([weights @ nodes**k for k in range(2*order)])
    if np.max(np.abs(reproduced-moment)/np.maximum(1., np.abs(moment))) > 2e-9:
        raise QuadratureError('Conditional quadrature does not reproduce moments')
    return nodes, weights


@dataclass
class BetaMeasure:
    """Unit-number bounded distribution on x in [0,1]."""
    mean: float
    variance: float
    adjustment: float = 0.

    @property
    def kappa(self):
        return self.mean*(1-self.mean)/self.variance-1

    @property
    def atoms(self):
        if self.variance == 0:
            return np.array([self.mean]), np.ones(1)
        if self.variance == self.mean*(1-self.mean):
            return np.array([0., 1.]), np.array([1-self.mean, self.mean])
        return None

    def cdf(self, x):
        x = np.asarray(x, dtype=float)
        atoms = self.atoms
        if atoms is not None:
            return sum(w*(x >= p) for p, w in zip(*atoms))
        k = self.kappa
        return betainc(self.mean*k, (1-self.mean)*k, np.clip(x, 0., 1.))

    def partial_first(self, x):
        """Integral x' dF(x') over x' <= x."""
        x = np.asarray(x, dtype=float)
        atoms = self.atoms
        if atoms is not None:
            return sum(w*p*(x >= p) for p, w in zip(*atoms))
        k = self.kappa
        return self.mean*betainc(self.mean*k+1, (1-self.mean)*k,
                                np.clip(x, 0., 1.))

    def _probability(self, lower, upper):
        k = self.kappa
        a, b = self.mean*k, (1-self.mean)*k
        fl, fu = float(betainc(a,b,lower)), float(betainc(a,b,upper))
        sl, su = float(betaincc(a,b,lower)), float(betaincc(a,b,upper))
        return fu-fl if fu <= sl else sl-su

    def _full(self, order):
        """Stable Jacobi recurrence, centered on the mean and scaled by SD."""
        mu, k, sd = self.mean, self.kappa, math.sqrt(self.variance)
        a, b = mu*k, (1-mu)*k
        diagonal = np.zeros(order)
        off = np.ones(order-1)
        for n in range(1, order):
            diagonal[n] = (-2*n*(2*mu-1)/(k+2*n)
                           *((k+n-1)/(k+2*n-2))/sd)
            if n > 1:
                off[n-1] = math.sqrt((n/(2*n+k-2))*((n+a-1)/(2*n+k-2))
                                     *((n+b-1)/(2*n+k-1))*((n+k-2)/(2*n+k-3)))/sd
        z, vectors = eigh_tridiagonal(diagonal, off, check_finite=False)
        return mu+sd*z, vectors[0]**2

    def _centered_integrals(self, lower, upper, degree):
        """Truncated moments about the full mean; endpoint recurrence."""
        mu, k = np.longdouble(self.mean), np.longdouble(self.kappa)
        a, b = float(mu*k), float((1-mu)*k)
        probability = np.longdouble(self._probability(lower, upper))
        def g(x):
            if x == 0 or x == 1:
                return np.longdouble(0.)
            return np.longdouble(math.exp(a*math.log(x)+b*math.log1p(-x)-betaln(a, b)))
        gl, gu = g(lower), g(upper)
        yl, yu = np.longdouble(lower)-mu, np.longdouble(upper)-mu
        result = [probability, (gl-gu)/k]
        for r in range(2, degree+1):
            result.append(((r-1)*mu*(1-mu)*result[r-2]
                           +(r-1)*(1-2*mu)*result[r-1]
                           -(yu**(r-1)*gu-yl**(r-1)*gl))/(k+r-1))
        return np.array(result)

    def _positive_integrals(self, lower, upper, degree, probability, stats=None):
        """Adaptive centered positive-density integration for difficult segments.

        The integration coordinate spans the segment; points around the Beta
        peak prevent a narrow peak in a broad interval from being missed. The
        normalization uses the independent incomplete-Beta probability. No
        probability cutoff or population floor enters this fallback.
        """
        center = min(upper, max(lower, self.mean))
        scale = min(.5*(upper-lower), math.sqrt(self.variance))
        a, b = self.mean*self.kappa, (1-self.mean)*self.kappa
        lognorm = float(betaln(a, b))+math.log(probability)-math.log(scale)
        zlower, zupper = (lower-center)/scale, (upper-center)/scale
        endpoint = 'left' if lower == 0 and a < 1 else ('right' if upper == 1 and b < 1 else None)
        if endpoint is not None:
            width = upper-lower
            exponent = a if endpoint == 'left' else b
            endpoint_norm = exponent*math.log(width)-math.log(exponent)-betaln(a,b)-math.log(probability)
        def coordinate_and_density(t):
            if np.isscalar(t):
                if endpoint == 'left':
                    x = width*t**(1/a)
                    return (x-center)/scale, math.exp(endpoint_norm+(b-1)*math.log1p(-x))
                if endpoint == 'right':
                    y = width*t**(1/b)
                    return ((1-center)-y)/scale, math.exp(endpoint_norm+(a-1)*math.log1p(-y))
                x = center+scale*t
                if not 0 < x < 1:
                    return t,0.
                return t,math.exp((a-1)*math.log(x)+(b-1)*math.log1p(-x)-lognorm)
            if endpoint == 'left':
                x = width*t**(1/a)
                return (x-center)/scale, np.exp(endpoint_norm+(b-1)*np.log1p(-x))
            if endpoint == 'right':
                y = width*t**(1/b)
                return ((1-center)-y)/scale, np.exp(endpoint_norm+(a-1)*np.log1p(-y))
            x = center+scale*t
            valid = (x > 0) & (x < 1)
            with np.errstate(divide='ignore', invalid='ignore'):
                value = np.exp((a-1)*np.log(x)+(b-1)*np.log1p(-x)-lognorm)
            return t, np.where(valid,value,0.)
        sd = math.sqrt(self.variance)
        points = [(self.mean+v*sd-center)/scale
                  for v in (-32, -16, -8, -4, -2, -1, 0, 1, 2, 4, 8, 16, 32)]
        points += list((-32., -16., -8., -4., -2., -1., 0., 1., 2., 4., 8., 16., 32.))
        points = sorted(set(p for p in points if zlower < p < zupper))
        if endpoint is not None:
            physical = [center+scale*p for p in points]
            points = ([(x/width)**a for x in physical] if endpoint == 'left'
                      else [((1-x)/width)**b for x in physical])
            points = sorted(set(p for p in points if 0 < p < 1))
            zlower, zupper = 0., 1.
        # Vectorized positive quadrature usually resolves these smooth,
        # transformed pieces cheaply. Admit it only after independent order
        # refinement of normalization and directly centered moments; otherwise
        # retain the adaptive integration below.
        edges = np.array([zlower,*points,zupper])
        half = .5*np.diff(edges)
        midpoint = .5*(edges[1:]+edges[:-1])
        previous = None
        probability_tolerance = max(2e-8,4*float(np.spacing(probability)/probability))
        for order in (16,32,64,128):
            gx,gw = _legendre(order)
            parameter = (midpoint[:,None]+half[:,None]*gx).ravel()
            z,density_value = coordinate_and_density(parameter)
            weight = (half[:,None]*gw).ravel()*density_value
            total = float(weight.sum())
            if total <= 0 or not math.isfinite(total):
                break
            weight /= total
            local_mean = float(weight@z)
            local_variance = float(weight@(z-local_mean)**2)
            if local_variance <= 0:
                break
            local_sd = math.sqrt(local_variance)
            standardized = np.array([weight@((z-local_mean)/local_sd)**r
                                     for r in range(degree+1)])
            signature = np.r_[total,local_mean,local_variance,standardized[3:]]
            if previous is not None:
                denominator = np.r_[1.,max(1.,abs(local_mean)),local_variance,
                                     np.maximum(1.,abs(standardized[3:]))]
                if (np.max(abs(signature-previous)/denominator) < 2e-10
                        and abs(total-1) <= probability_tolerance):
                    standardized[:3] = [1.,0.,1.]
                    if stats is not None:
                        stats['refined_gauss_integrations'] += 1
                    return center+scale*local_mean,scale*local_sd,standardized
            previous = signature
        if stats is not None:
            stats['adaptive_segment_integrations'] += 1
        def raw_integrand(t):
            z, density = coordinate_and_density(t)
            return density*z**np.arange(3)
        value, error = quad_vec(raw_integrand, zlower, zupper,
                                epsabs=2e-11, epsrel=2e-11, points=points, limit=400)
        # Subnormal CDF probabilities have an explicitly quantized relative
        # precision. This tolerance follows float spacing, never a mass floor.
        probability_tolerance = max(2e-8, 4*float(np.spacing(probability)/probability))
        if (not np.isfinite(value).all() or error > 2e-8*max(1.,np.max(abs(value)))
                or abs(value[0]-1) > probability_tolerance):
            raise QuadratureError('Adaptive truncated-Beta integration failed')
        value /= value[0]
        mean_local = value[1]
        variance = value[2]-mean_local*mean_local
        if variance <= 0:
            raise QuadratureError('Truncated distribution lost positive variance')
        sd_local = math.sqrt(variance)
        # Integrate centered powers directly. Binomial subtraction of broad
        # interval moments can destroy orders3..7 for a narrow conditional peak.
        powers = np.arange(degree+1)
        def centered_integrand(t):
            z, density = coordinate_and_density(t)
            return density*((z-mean_local)/sd_local)**powers
        centered, error = quad_vec(centered_integrand, zlower, zupper,
            epsabs=2e-10, epsrel=2e-11, points=points, limit=400)
        if not np.isfinite(centered).all() or error > 2e-8*max(1.,np.max(abs(centered))):
            raise QuadratureError('Centered truncated-Beta integration failed')
        centered /= centered[0]
        if abs(centered[1]) > 2e-8 or abs(centered[2]-1) > 2e-8:
            raise QuadratureError('Conditional moments disagree between integration passes')
        centered[:3] = [1.,0.,1.]
        return center+scale*mean_local, scale*sd_local, centered

    def quadrature(self, order=2, lower=0., upper=1., *, stats=None):
        if not (0 <= lower <= upper <= 1) or not 2 <= order <= 6:
            raise ValueError('Ordered unit bounds and quadrature order 2..6 required')
        if lower == upper:
            return np.empty(0), np.empty(0)
        atoms = self.atoms
        if atoms is not None:
            nodes, weights = atoms
            # Adjacent segments own their lower endpoint, but only the final
            # interval owns x=1. This avoids duplicating endpoint atoms.
            keep = (nodes >= lower) & ((nodes < upper) | ((upper == 1) & (nodes == 1)))
            return nodes[keep], weights[keep]
        if lower == 0 and upper == 1:
            return self._full(order)
        probability = self._probability(lower, upper)
        if probability < 0 or not math.isfinite(probability):
            raise QuadratureError('Invalid incomplete-Beta segment probability')
        if probability == 0:
            return np.empty(0), np.empty(0)
        # Higher truncated moments and narrow intervals amplify recurrence
        # cancellation. Integrate them directly in a centered local coordinate.
        fallback = order > 2 or upper-lower < 2*math.sqrt(self.variance)
        if not fallback:
            raw = self._centered_integrals(lower, upper, 2*order-1)/probability
            shift = raw[1]
            variance = raw[2]-shift*shift
            if variance <= 0:
                fallback = True
            else:
                sd = np.sqrt(variance)
                mean = np.longdouble(self.mean)+shift
                moments = np.array([sum(math.comb(r,j)*(-shift)**(r-j)*raw[j]
                                         for j in range(r+1))/sd**r
                                    for r in range(2*order)], dtype=float)
                moments[:3] = [1., 0., 1.]
                try:
                    nodes, weights = _gaussian_from_moments(
                        moments, order, float((lower-mean)/sd), float((upper-mean)/sd))
                    return np.asarray(mean+sd*nodes, float), probability*weights
                except QuadratureError:
                    fallback = True
        if fallback:
            if stats is not None:
                stats['numerical_segment_integrations'] += 1
            mean, sd, moments = self._positive_integrals(lower, upper, 2*order-1, probability, stats)
            nodes, weights = _gaussian_from_moments(moments, order,
                                                   (lower-mean)/sd, (upper-mean)/sd)
            return mean+sd*nodes, probability*weights


class MomentGrid:
    """Fixed cells with three local moments; masses retain continuous values."""
    def __init__(self, pivots):
        self.pivots = np.asarray(pivots, dtype=float).copy()
        if (self.pivots.ndim != 1 or self.pivots.size < 2
                or not np.isfinite(self.pivots).all() or (self.pivots <= 0).any()
                or (np.diff(self.pivots) <= 0).any()):
            raise ValueError('At least two strictly increasing positive pivots required')
        lower = self.pivots[0]**2/self.pivots[1]
        upper = self.pivots[-1]**2/self.pivots[-2]
        self.bounds = .5*(np.r_[lower, self.pivots]+np.r_[self.pivots, upper])
        self.widths = np.diff(self.bounds)
        self.offsets = (self.pivots-self.bounds[:-1])/self.widths

    def from_particles(self, masses, numbers):
        masses, numbers = np.broadcast_arrays(np.asarray(masses, float), np.asarray(numbers, float))
        if (not np.isfinite(masses).all() or not np.isfinite(numbers).all()
                or (numbers < 0).any() or (masses < self.bounds[0]).any()
                or (masses >= self.bounds[-1]).any()):
            raise ValueError('Finite nonnegative populations inside the mass domain required')
        result = np.zeros((self.pivots.size, 3))
        index = np.searchsorted(self.bounds, masses, side='right')-1
        u = (masses-self.pivots[index])/self.widths[index]
        for k in range(3):
            np.add.at(result[:, k], index, numbers*u**k)
        self.distributions(result)
        return result

    def distributions(self, state):
        state = np.asarray(state, float)
        if state.shape != (self.pivots.size, 3) or not np.isfinite(state).all():
            raise RealizabilityError('Finite state with shape (cells,3) required')
        if (state[:,0] < 0).any():
            raise RealizabilityError('Negative population')
        output = []
        for i, (number, first, second) in enumerate(state):
            if number == 0:
                if first != 0 or second != 0:
                    raise RealizabilityError('Empty cell has nonzero moments')
                output.append(None)
                continue
            local_mean = np.longdouble(first)/number
            local_second = np.longdouble(second)/number
            variance = local_second-local_mean*local_mean
            mean = local_mean+np.longdouble(self.offsets[i])
            mean_tolerance = 64*EPS*(abs(local_mean)+abs(self.offsets[i]))
            variance_tolerance = 64*EPS*(abs(local_second)+local_mean*local_mean)
            if mean < -mean_tolerance or mean > 1+mean_tolerance:
                raise RealizabilityError(f'Cell {i} mean leaves support')
            bounded_mean = min(np.longdouble(1), max(np.longdouble(0), mean))
            maximum = bounded_mean*(1-bounded_mean)
            if variance < -variance_tolerance or variance > maximum+variance_tolerance+mean_tolerance:
                raise RealizabilityError(f'Cell {i} variance leaves support: {float(variance):.17g}')
            bounded_variance = min(maximum, max(np.longdouble(0), variance))
            if abs(variance) <= variance_tolerance:
                bounded_variance = np.longdouble(0)
            elif abs(variance-maximum) <= variance_tolerance:
                bounded_variance = maximum
            adjustment = max(abs(float(mean-bounded_mean)), abs(float(variance-bounded_variance)))
            mean_float = float(bounded_mean)
            variance_float = float(bounded_variance)
            if bounded_variance == maximum:
                variance_float = mean_float*(1-mean_float)
            measure = BetaMeasure(mean_float, variance_float, adjustment)
            if measure.atoms is None and not math.isfinite(measure.kappa):
                raise RealizabilityError('Beta concentration exceeds floating-point range')
            output.append(measure)
        return output

    def moments(self, state):
        n, u, u2 = np.asarray(state).T
        c, h = self.pivots, self.widths
        return np.array([n.sum(), np.sum(c*n+h*u),
                         np.sum(c*c*n+2*c*h*u+h*h*u2)])

    def variance(self, state):
        n = np.asarray(state)[:,0]
        total = n.sum()
        if total == 0:
            return 0.
        mean = self.moments(state)[1]/total
        measures = self.distributions(state)
        return float(sum(n[i]*((self.pivots[i]+self.widths[i]*(d.mean-self.offsets[i])-mean)**2
                               +self.widths[i]**2*d.variance)
                         for i,d in enumerate(measures) if d is not None)/total)

    def cdf(self, state, queries, weighting='number'):
        if weighting not in ('number', 'mass'):
            raise ValueError('CDF weighting must be number or mass')
        queries = np.asarray(queries, float)
        output = np.zeros_like(queries)
        for i,d in enumerate(self.distributions(state)):
            if d is None:
                continue
            x = (queries-self.bounds[i])/self.widths[i]
            partial = d.cdf(x)
            if weighting == 'mass':
                partial = self.bounds[i]*partial+self.widths[i]*d.partial_first(x)
            output += state[i,0]*partial
        total = self.moments(state)[0 if weighting == 'number' else 1]
        return output/total if total else output


def _kernel(kernel, left, right):
    values = np.asarray(kernel(left, right), float)
    values = np.broadcast_to(values, np.broadcast_shapes(np.shape(left), np.shape(right)))
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError('Kernel values must be finite and nonnegative')
    return values


class _BaseOperator:
    def __init__(self, grid, kernel, quadrature_order):
        if not 2 <= quadrature_order <= 6:
            raise ValueError('Quadrature order 2..6 required')
        self.grid, self.kernel, self.quadrature_order = grid, kernel, quadrature_order
        self.stats = dict(numerical_segment_integrations=0, adaptive_segment_integrations=0,
                          refined_gauss_integrations=0, kernel_evaluations=0)
        self._symmetry_checked = False

    def _segments(self, measure, cell, donor):
        grid = self.grid
        cuts = (grid.bounds-donor-grid.bounds[cell])/grid.widths[cell]
        cuts = np.r_[0., cuts[(cuts > 0) & (cuts < 1)], 1.]
        nodes, weights = [], []
        for lower, upper in zip(cuts[:-1], cuts[1:]):
            x, w = measure.quadrature(self.quadrature_order, float(lower), float(upper), stats=self.stats)
            nodes.extend(x)
            weights.extend(w)
        x, weights = np.asarray(nodes), np.asarray(weights)
        return grid.pivots[cell]+grid.widths[cell]*(x-grid.offsets[cell]), weights

    def _assemble(self, state, outer, inner, donor, collector, weight, bath):
        grid = self.grid
        rhs = np.zeros_like(state)
        ledger = dict(supplied_mass=0., collision_number_loss=0.,
                      lower_exit=np.zeros(3), upper_exit=np.zeros(3))
        if not len(weight):
            return rhs, ledger
        outer, inner = np.asarray(outer,int), np.asarray(inner,int)
        donor, collector, weight = map(np.asarray, (donor, collector, weight))
        rates = weight*_kernel(self.kernel, donor, collector)
        if not self._symmetry_checked:
            reverse = _kernel(self.kernel, collector, donor)
            if not np.allclose(reverse, _kernel(self.kernel, donor, collector), rtol=1e-12, atol=0):
                raise ValueError('Collision kernel must be symmetric')
            self._symmetry_checked = True
        self.stats['kernel_evaluations'] += rates.size
        product = donor+collector
        destination = np.searchsorted(grid.bounds, product, side='right')-1
        inside = (destination >= 0) & (destination < grid.pivots.size)
        same = inside & (destination == inner)
        u = (collector-grid.pivots[inner])/grid.widths[inner]
        increment = donor/grid.widths[inner]
        # Stable replacement of a collector that remains in its cell.
        np.add.at(rhs[:,1], inner[same], rates[same]*increment[same])
        np.add.at(rhs[:,2], inner[same], rates[same]*(2*u[same]*increment[same]+increment[same]**2))
        separate = ~same
        for k in range(3):
            np.add.at(rhs[:,k], inner[separate], -rates[separate]*u[separate]**k)
        crossing = inside & separate
        dest = destination[crossing]
        up = (product[crossing]-grid.pivots[dest])/grid.widths[dest]
        for k in range(3):
            np.add.at(rhs[:,k], dest, rates[crossing]*up**k)
        if bath:
            ledger['supplied_mass'] = float(rates@donor)
        else:
            ud = (donor-grid.pivots[outer])/grid.widths[outer]
            for k in range(3):
                np.add.at(rhs[:,k], outer, -rates*ud**k)
            ledger['collision_number_loss'] = float(rates.sum())
        for label, selector in (('lower_exit', destination < 0),
                                ('upper_exit', destination >= grid.pivots.size)):
            ledger[label] = np.array([np.sum(rates[selector]*product[selector]**k) for k in range(3)])
        return rhs, ledger


class BathCapture(_BaseOperator):
    """Collector evolution in a prescribed replenished monodisperse bath."""
    def __init__(self, grid, kernel, bath_mass, bath_number, quadrature_order=2):
        super().__init__(grid, kernel, quadrature_order)
        if not math.isfinite(bath_mass) or not math.isfinite(bath_number) or min(bath_mass,bath_number) <= 0:
            raise ValueError('Positive finite bath mass and number required')
        self.bath_mass, self.bath_number = bath_mass, bath_number

    def evaluate(self, state):
        outer, inner, donor, collector, weight = [], [], [], [], []
        for j, measure in enumerate(self.grid.distributions(state)):
            if measure is None:
                continue
            mass, w = self._segments(measure, j, self.bath_mass)
            count = len(mass)
            outer.extend([0]*count)
            inner.extend([j]*count)
            donor.extend([self.bath_mass]*count)
            collector.extend(mass)
            weight.extend(state[j,0]*self.bath_number*w)
        return self._assemble(state, outer, inner, donor, collector, weight, True)


class CollisionMoments(_BaseOperator):
    """Closed-box Smoluchowski aggregation with a local bounded-Beta closure."""
    def __init__(self, grid, kernel, quadrature_order=2):
        super().__init__(grid, kernel, quadrature_order)

    def evaluate(self, state):
        distributions = self.grid.distributions(state)
        active = [i for i,d in enumerate(distributions) if d is not None]
        outer, inner, donor, collector, weight = [], [], [], [], []
        for ii, i in enumerate(active):
            x, w = distributions[i].quadrature(self.quadrature_order)
            masses = self.grid.pivots[i]+self.grid.widths[i]*(x-self.grid.offsets[i])
            for j in active[ii:]:
                factor = .5 if i == j else 1.
                for mi, wi in zip(masses, w):
                    mj, wj = self._segments(distributions[j], j, mi)
                    count = len(mj)
                    outer.extend([i]*count)
                    inner.extend([j]*count)
                    donor.extend([mi]*count)
                    collector.extend(mj)
                    weight.extend(factor*state[i,0]*state[j,0]*wi*wj)
        return self._assemble(state, outer, inner, donor, collector, weight, False)


def integrate(operator, initial_state, times_s, *, max_step_s=30., cfl=.3, maximum_steps=200000,
              time_limit_s=None, progress=None):
    """SSPRK2 with trial-stage realizability rejection and explicit ledgers.

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
                    candidate = .5*state+.5*second
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
    if not stats['steps']:
        stats['minimum_step_s'] = 0.
    return dict(snapshots=snapshots, diagnostics=stats)
