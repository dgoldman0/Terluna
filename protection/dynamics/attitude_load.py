"""Ideal rigid-square attitude load and an explicit zero-net-force actuator.

Three independent couples act at opposite edge midpoints, with lever arm
equal to the tile side. This is a propulsion accounting option, not a chosen
structural/actuator design. Wheels, sail torque and flexible modes need their
own hardware and dynamics models.
"""
import numpy as np

from shared import constants as K
from .optical import length, unit
from .cycling import ray_geometry


def rigid_square_load(command, times, side=10000.):
    """Euler torque per unit mass in the square's rotating body frame."""
    omega, alpha = command(times, 1), command(times, 2)
    inertia = side**2/12*np.array([1., 1., 2.])
    torque = alpha*inertia+np.cross(omega, omega*inertia)
    # Each torque component uses an opposite pair F, -F separated by side.
    force_per_mass = 2*np.abs(torque).sum(axis=-1)/side
    return dict(omega=omega, alpha=alpha, torque_per_mass=torque,
        force_per_mass=force_per_mass,
        kinetic_energy_per_mass=.5*np.sum(omega**2*inertia, axis=-1),
        angular_momentum_per_mass=omega*inertia)


def disturbance_torque_bound(env, t, state, frame, side=10000.):
    """Pointwise norm envelope for gravity-gradient and reflected-band torque.

    A square lies in its circumsphere. Point-mass gravity has gradient norm
    at most 2 GM/d^3 there; a uniform square has mean squared radius side²/6.
    Reflected momentum lies along the shared normal. Its largest permitted
    local cosine bounds either illuminated face, allowing arbitrary shadows.
    """
    q = state[:, :3]; sample = env.at(t); radius = side/np.sqrt(2)
    gradient = 2*K.MOON_GM/(length(q)-radius)**3
    for name, gm in zip(env.names, env.gm):
        gradient += 2*gm/(length(sample['positions'][name]-q)-radius)**3
    gravity = side**2/6*gradient
    rays = ray_geometry(q, sample); distance = length(rays['sun'])
    cone = np.arcsin((K.SUN_RADIUS+radius)/(distance-radius))
    cosine = np.minimum(1., abs(unit(rays['sun'])@frame[:, 2])+np.sin(cone))
    force = (2*env.central_fraction*K.SOLAR_CONSTANT/(K.SPEED_OF_LIGHT*env.sigma)*
        (K.AU/(distance-radius))**2*cosine**2*((side-110.)/side)**2)
    radiation = radius*force
    # L1 of an unknown torque vector is at most sqrt(3) times its norm.
    couple_force = 2*np.sqrt(3)*(gravity+radiation)/side
    return gravity, radiation, couple_force
