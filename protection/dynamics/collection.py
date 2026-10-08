"""An optical-energy inventory for parallel finite tiles at actual dates.

First-intercept power counts a source ray once, on its first clear aperture.
It is a bolometric geometric equivalent, not the power absorbed by a filter.
The separately reported redirected band uses the dynamics' ideal spectral
fraction. Electrical conversion, collection of the outgoing beam and delivery
require an explicit device and routing model; this module grants no such credit.
"""
import numpy as np

from shared import constants as K
from .cycling import ray_geometry
from .fleet import sun_points
from .natural_pattern import retarded_sun
from .parallel_shadow import source_illumination
from .optical import length, unit


def parallel_interception(q, frame, sources, band_fraction, clear_side=9890., weights=None):
    """Per-tile W from a specified, unocculted finite-source quadrature.

    Shared square planes permit exact sampled-source shadow unions. As in
    the parent force model, irradiance and incidence are evaluated at tile
    centres. Both faces count; their separate ledgers expose any one-sided
    collector assumption. No transmitted-band energy is added downstream.
    """
    q, frame, sources = np.asarray(q), np.asarray(frame), np.asarray(sources)
    if q.ndim != 2 or q.shape[1] != 3 or not len(q) or not np.isfinite(q).all():
        raise ValueError('Finite tile centres must have shape (members, 3)')
    if frame.shape != (3, 3) or not np.allclose(frame.T@frame, np.eye(3), atol=1e-10):
        raise ValueError('A shared orthonormal square frame is required')
    if sources.ndim != 2 or sources.shape[1] != 3 or not len(sources):
        raise ValueError('Finite source points must have shape (sources, 3)')
    if not np.isfinite(sources).all() or not 0 <= band_fraction <= 1 or clear_side <= 0:
        raise ValueError('Invalid source, spectral fraction or aperture')
    # source_illumination uses a 10 km physical-square neighbour bound.
    if clear_side > 10000.:
        raise ValueError('This neighbour reduction supports apertures up to 10 km')
    weights = np.full(len(sources), 1/len(sources)) if weights is None else np.asarray(weights)
    if weights.shape != (len(sources),) or np.any(weights < 0) or not np.isclose(weights.sum(), 1.):
        raise ValueError('Nonnegative source weights must sum to one')
    n = frame[:, 2]
    unshadowed = np.zeros(len(q)); first = np.zeros(len(q))
    front = np.zeros(len(q)); back = np.zeros(len(q)); force = np.zeros((len(q), 3))
    for source, weight in zip(sources, weights):
        incident = source-q
        distance = length(incident)
        if np.any(distance <= clear_side):
            raise ValueError('Source lies within the finite-tile neighbourhood')
        cosine = unit(incident)@n
        power = weight*K.SOLAR_CONSTANT*(K.AU/distance)**2*clear_side**2*abs(cosine)
        visible = source_illumination(q, frame, source, clear_side)
        unique = power*visible
        unshadowed += power; first += unique
        front += np.where(cosine >= 0, unique, 0.)
        back += np.where(cosine < 0, unique, 0.)
        # Incoming minus specular outgoing momentum, on either illuminated face.
        force -= (2*band_fraction/K.SPEED_OF_LIGHT*unique*cosine)[:, None]*n
    return dict(unshadowed_aperture_W=unshadowed,
        first_intercept_bolometric_equivalent_W=first,
        front_first_intercept_W=front, back_first_intercept_W=back,
        redirected_band_optical_W=first*band_fraction,
        mutual_shadow_equivalent_loss_W=unshadowed-first,
        reflected_force_N=force)


def require_uneclipsed_tiles(rays, side=10000.):
    """Reject body eclipses, including a conservative whole-tile envelope.

    Inflating the apparent Sun and each body by the tile's circumsphere
    encloses their directions from any point on the tile. Overlapping cones
    require the future joint source/area body-plus-array visibility integral.
    Multiplying independent body and array visibility fractions is unsafe.
    """
    radius = side/np.sqrt(2)
    sun = rays['sun']; ds = length(sun)
    if np.any(ds <= K.SUN_RADIUS+radius):
        raise ValueError('Tile intersects the solar exclusion')
    solar_angle = np.arcsin((K.SUN_RADIUS+radius)/ds)
    for name, body_radius in [('earth', K.EARTH_RADIUS), ('moon', K.MOON_RADIUS)]:
        body = rays[name]; distance = length(body)
        if np.any(distance <= body_radius+radius):
            raise ValueError('Tile intersects a solid-body exclusion')
        separation = np.arctan2(length(np.cross(sun, body)), np.sum(sun*body, axis=-1))
        angle = np.arcsin((body_radius+radius)/distance)
        possible = (distance-body_radius < ds+K.SUN_RADIUS) & (separation <= solar_angle+angle)
        if np.any(possible):
            raise ValueError('Body eclipse requires joint source/area collection and force accounting')


def uneclipsed_collection(env, t, state, frame, suns=16, rotation=.317, side=10000.):
    """Collection potential on the existing uneclipsed parallel-pattern model."""
    sample = env.at(t)
    require_uneclipsed_tiles(ray_geometry(state[:, :3], sample), side)
    sources = sun_points(retarded_sun(sample), suns, rotation)
    return parallel_interception(state[:, :3], frame, sources, env.central_fraction, side-110.)


def integrate_inventory(times, power):
    """Integrate a (time, member) W ledger, retaining per-member joules."""
    times, power = np.asarray(times), np.asarray(power)
    if (times.ndim != 1 or len(times) < 2 or np.any(np.diff(times) <= 0) or
            power.ndim != 2 or power.shape[0] != len(times) or
            not np.isfinite(times).all() or not np.isfinite(power).all() or np.any(power < 0)):
        raise ValueError('Power requires finite, nonnegative samples at strictly increasing dates')
    per_member = np.trapezoid(power, times, axis=0)
    total = power.sum(axis=1)
    return dict(energy_J=float(per_member.sum()),
        mean_W=float(per_member.sum()/(times[-1]-times[0])),
        sampled_minimum_W=float(total.min()), sampled_maximum_W=float(total.max()),
        per_member_energy_J=per_member.tolist())
