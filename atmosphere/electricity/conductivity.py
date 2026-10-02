"""The air's small ions and conductivity: ion pairs made at rate q are lost by recombination with each other and by
attachment to aerosol particles and cloud droplets, so that in steady state

    q = alpha n^2 + S n,   S = sum_j beta_j Z_j,

with n the density of either sign, alpha the ion-ion recombination coefficient and S the attachment rate to particles
of concentration Z_j and attachment coefficient beta_j. The conductivity is sigma = e n (mu+ + mu-), each mobility
scaling inversely with the air's number density, and charge in the air relaxes with time constant eps0 / sigma.
"""
from __future__ import annotations
import numpy as np

ELEMENTARY_CHARGE = 1.602176634e-19          # C
EPSILON_0 = 8.8541878128e-12                 # F/m


def ion_density(q, alpha, sink=0.0):
    """Steady small-ion density (m-3) of each sign for production q (m-3 s-1), recombination alpha (m3 s-1) and
    attachment rate sink (s-1)."""
    q, alpha, sink = (np.asarray(a, float) for a in (q, alpha, sink))
    return 2.0 * q / (sink + np.sqrt(sink ** 2 + 4.0 * alpha * q))      # the positive root, stable when sink >> alpha n


def mobility(mu_standard, number_density, standard_number_density):
    """A small ion's mobility (m2 V-1 s-1) at the air's number density, from its reduced mobility at standard
    conditions."""
    return mu_standard * standard_number_density / np.asarray(number_density, float)


def conductivity(n, mu_positive, mu_negative):
    """Conductivity (S/m) from equal densities n (m-3) of positive and negative small ions."""
    return ELEMENTARY_CHARGE * np.asarray(n, float) * (np.asarray(mu_positive, float) + np.asarray(mu_negative, float))


def relaxation_time(sigma):
    """The time (s) over which charge in air of conductivity sigma (S/m) leaks away: eps0 / sigma."""
    return EPSILON_0 / np.asarray(sigma, float)
