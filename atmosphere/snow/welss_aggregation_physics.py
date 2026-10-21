"""Selected McSnow aggregation closures for the declared Earth experiment.

Pure ice, habit0, agggeo1, phi=1, sticking option1, no random speed perturbation.
This adapter supplies component equations; it is not a coupled population solver
and does not alter the frozen lunar column or the Earth growth reference.

Derived from McSnow revision 1e428d7e901cbc5b91f2040e575f736ec5963a74:
Copyright (c) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss.
Authors: Slavko Brdar, Axel Seifert, Christoph Siewert, Markus Karrer,
Jan-Niklas Welss. SPDX-License-Identifier: BSD-3-Clause.
See references/mcsnow_bsd_3_clause.txt and welss_aggregation_physics.md.
"""
from __future__ import annotations

import numpy as np

from atmosphere.snow import boehm
from atmosphere.snow import welss_column as w


def geometry(mass, monomer):
    """Source mass geometry and monomer/aggregate capacitance [m].

    monomer is a boolean flag for the source condition mm<1.01. Actual monomer
    count has no other effect in this selected diagnostic-geometry experiment.
    The source's rounding-sensitive void-volume sign remains unchanged.
    """
    flag = np.asarray(monomer)
    if flag.dtype.kind != 'b':
        raise ValueError('monomer must be a boolean flag')
    mass, flag = np.broadcast_arrays(np.asarray(mass, float), flag)
    result = w.geometry(mass)
    porous_ratio = np.where(flag, 1 / w.PI, .25)
    result['C'] = result['D'] * np.where(result['void'] > 0, porous_ratio, .5)
    return result


def microphysics(mass, monomer, state):
    """Frozen Earth growth equations with the source aggregate C branch."""
    shape = geometry(mass, monomer)
    mass = np.broadcast_to(np.asarray(mass, float), shape['D'].shape)
    result = w.microphysics(mass, state)
    result['growth'] = result['growth'] * (shape['C'] / result['C'])
    result['C'] = shape['C']
    return result


def sticking(temperature_k):
    """Exact intervals of McSnow stick_effi_tempOnly option1."""
    t = np.asarray(temperature_k, float)
    if not np.isfinite(t).all() or (t <= 0).any():
        raise ValueError('temperature must be finite and positive')
    tc = t - w.SOURCE_CONSTANTS['T_3']
    return np.select([tc >= 0., tc >= -4., tc >= -6., tc >= -9.,
                      tc >= -12.5, tc >= -17., tc >= -20.],
                     [1., .1, .6, .1, .4, 1., .4], default=.25)


def kernel(mass_left, mass_right, state):
    """Broadcast hydrodynamic aggregation kernel [m3/s] in shared air.

    For phi=1, ckt0 uses outer radii D/2, not sqrt(projected area/pi).
    Rough-ice internal speed width SIG=.25 enters interception, while the outer
    encounter rate uses deterministic terminal-speed differences. SIG is not a
    measured ambient turbulent velocity. Raw source Ec is never clipped to1.
    """
    left, right = np.broadcast_arrays(np.asarray(mass_left, float),
                                      np.asarray(mass_right, float))
    a, b = w.geometry(left), w.geometry(right)
    p = boehm.particle_properties(a['D'], left, w.SOURCE_CONSTANTS['grav'], state,
                                 a['q'], sigma=.25)
    q = boehm.particle_properties(b['D'], right, w.SOURCE_CONSTANTS['grav'], state,
                                 b['q'], sigma=.25)
    efficiency = boehm.collision_efficiency(p, q, state)
    area = w.PI / 4 * (a['D'] + b['D'])**2
    relative_speed = np.abs(p['v'] - q['v'])
    stick = sticking(state['T'])
    return dict(K=efficiency * stick * area * relative_speed, Ec=efficiency,
                sticking=stick, encounter_area_m2=area,
                relative_speed_m_s=relative_speed)
