"""Held screen with a zoned aperture (stage 4 of relative_orbits.md).

The published held-screen benchmark gives every tile 50 g/m², the climate stack.
Only the central window, whose rays reach the solid Moon, has to deliver the
climate spectrum; the annulus out to four lunar radii has to stop the UV, and
its rays continue past the Moon toward Earth near eclipse-season new moons, so
it passes visible light. This runner repeats the published area quadrature on
the selected moving trajectory (holding.json, variable_distance) with a base
areal mass per zone and a chosen redirected fraction for annulus tiles, with
the same per-tile power, propellant and mass closure.

The first case repeats the published 50 g/m² result with the published optics
as a check. Annulus tiles that pass visible light redirect about as much light
as the climate stack (the stack's 13.8%); the annulus film's own UV
transmission comes from the protection domain.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.zoned_aperture
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.ephemeris import DEFAULT_KERNEL, Ephemeris
from protection.dynamics.model import geometry, relative_gravity, target
from protection.dynamics.optical import (arriving_ray_vectors, axial_optical, covering_radius, length,
                                         light_time_allowance, optical_projection, solar_visibility)
from .aperture import disk_nodes
from .run import CORE_DIVERTED, SCENARIO

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
HOLDING = HERE/'results/holding.json'
OUT = HERE/'results/zoned_aperture.json'
FILES = ['research/studies/solar_shield_array/zoned_aperture.py', 'research/studies/solar_shield_array/aperture.py',
         'research/studies/solar_shield_array/run.py', 'protection/dynamics/model.py', 'protection/dynamics/optical.py',
         'protection/dynamics/ephemeris.py']


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def zoned_aperture(samples, geo, coefficients, core_sigma, annulus_sigma, annulus_fraction, rings=4, angles=12,
                   protected_radii=4, fields=False):
    """evaluate_aperture of aperture.py with a base areal mass per zone and a chosen annulus redirected fraction.

    With fields, also returns each node's position and velocity relative to the Moon, the propulsive
    acceleration it needs (the residual after the optical bound), its closed areal mass and its weight.
    """
    trajectory = target(geo, coefficients)
    distance = trajectory['coords'][:, 0]
    offset = length(trajectory['coords'][:, 1:])
    sun_distance = length(samples['positions']['sun'])
    ray_margin = light_time_allowance(distance, samples['moon_v'], samples['sun_v'])
    radius = float((covering_radius(distance, sun_distance, protected_radii*K.MOON_RADIUS, offset)+ray_margin).max()
                   + SCENARIO['formation_margin_m'])
    core_radii = covering_radius(distance, sun_distance, K.MOON_RADIUS)+ray_margin
    core_radius = min(float((core_radii+offset).max()), radius)
    nodes, weights = disk_nodes(radius, core_radius, rings, angles)
    in_core = length(nodes) <= core_radius
    displacement = np.einsum('pi,nij->npj', nodes, geo['frame'][:, 1:])
    q = trajectory['q'][:, None, :]+displacement
    acc = trajectory['a'][:, None, :]+np.einsum('pi,nij->npj', nodes, geo['frame_dd'][:, 1:])
    positions = {k: v[:, None, :] for k, v in samples['positions'].items()}
    required = acc-relative_gravity(q, positions, samples['moon_a'][:, None, :])
    sun, earth = arriving_ray_vectors(positions['sun']-q, positions['earth']-q, samples['sun_v'][:, None, :],
                                      samples['earth_v'][:, None, :], samples['sun_a'][:, None, :],
                                      samples['earth_a'][:, None, :])
    visibility = solar_visibility(sun, earth)
    absolute_xy = trajectory['coords'][:, None, 1:]+nodes[None, :, :]
    # Rays that can reach the solid Moon keep the climate spectrum; the others redirect the annulus fraction.
    diverted = np.where(length(absolute_xy) <= core_radii[:, None], CORE_DIVERTED, annulus_fraction)
    prop = SCENARIO['propulsion']
    ve, eta, kappa = prop['exhaust_velocity_m_s'], prop['efficiency'], prop['specific_power_W_kg']
    cant = np.cos(np.deg2rad(prop['cant_deg']))
    buffer = prop['propellant_buffer_days']*K.JULIAN_DAY
    base = np.where(in_core, core_sigma, annulus_sigma)
    sigma = base.copy()
    t = samples['t']
    for _ in range(100):
        residual = optical_projection(required, sun, sigma, diverted, visibility)[1]
        mag = length(residual)
        mean = np.trapezoid(mag, t, axis=0)/(t[-1]-t[0])
        peak = mag.max(axis=0)*prop['peak_margin_factor']
        mdot_area = sigma*mean/(ve*cant)
        power_area = mdot_area*ve*ve/(2*eta)
        peak_area = sigma*peak*ve/(2*eta*cant)
        new_sigma = base+peak_area/kappa+mdot_area*buffer
        if np.max(abs(new_sigma-sigma)) < 1e-11:
            break
        sigma = new_sigma
    else:
        raise RuntimeError('Per-tile mass closure did not converge')
    area = np.pi*radius*radius
    core_weight = float(weights[in_core].sum())
    result = dict(core_sigma_kg_m2=core_sigma, annulus_sigma_kg_m2=annulus_sigma, annulus_redirected_fraction=annulus_fraction,
                  aperture_radius_km=radius/1e3, core_radius_km=core_radius/1e3, core_area_fraction=core_weight,
                  mean_power_TW=float(area*(power_area@weights)/1e12),
                  installed_design_peak_power_TW=float(area*(peak_area@weights)/1e12),
                  propellant_kg_s=float(area*(mdot_area@weights)),
                  base_optical_mass_kg=float(area*(base@weights)), total_mass_kg=float(area*(sigma@weights)),
                  power_hardware_mass_kg=float(area*(peak_area@weights)/kappa),
                  mean_power_core_TW=float(area*(power_area[in_core]@weights[in_core])/1e12),
                  mean_power_annulus_TW=float(area*(power_area[~in_core]@weights[~in_core])/1e12))
    if not fields:
        return result
    velocity = trajectory['v'][:, None, :]+np.einsum('pi,nij->npj', nodes, geo['frame_d'][:, 1:])
    return result, dict(t=t, nodes=nodes, weights=weights, in_core=in_core, q=q, v=velocity, residual=residual,
                        sigma=sigma, aperture_radius_m=radius)


def main():
    wall, cpu0 = time.monotonic(), time.process_time()
    holding = json.loads(HOLDING.read_text())
    coefficients = np.array(holding['trajectory_search']['variable_distance']['coefficients_km'])
    days = SCENARIO['primary_days']
    eph = Ephemeris(DEFAULT_KERNEL, SCENARIO['epoch_tdb'])
    coarse = eph.sample(np.linspace(0, days*K.JULIAN_DAY, int(days*8)+1))
    eph.close()
    geo = geometry(coarse)
    published = holding['aperture_quadrature']['variable_distance']
    cases = [('published_check', .05, .05, 1.0)]
    cases += [('visible_passing_annulus_50g', .05, .05, CORE_DIVERTED)]
    cases += [(f'annulus_{int(a*1000)}g', .05, a, CORE_DIVERTED) for a in (.02, .01, .005)]
    cases += [(f'core_26g_annulus_{int(a*1000)}g', .026, a, CORE_DIVERTED) for a in (.01, .005)]
    results = []
    for label, core, annulus, fraction in cases:
        results.append(dict(label=label, **zoned_aperture(coarse, geo, coefficients, core, annulus, fraction)))
        print(f"{label:28s} {results[-1]['mean_power_TW']:8.2f} TW  {results[-1]['propellant_kg_s']:9.0f} kg/s  "
              f"mass {results[-1]['total_mass_kg']:.3e} kg  core {results[-1]['mean_power_core_TW']:.2f} TW")
    out = dict(schema='terluna.research.zoned-aperture/1',
               producer=dict(files={f: digest(ROOT/f) for f in FILES},
                             constants=constants_used([f for f in FILES if f.endswith('.py')]),
                             inputs={'results/holding.json': digest(HOLDING), 'de440s.bsp': digest(DEFAULT_KERNEL),
                                     'scenario.json': digest(HERE/'scenario.json')}),
               evidence=('Area quadrature of independently held tiles on the published selected moving trajectory, '
                         'DE440s inverse dynamics sampled every three hours over the primary year, the ideal '
                         'optical-control bound with per-zone redirected fractions, and the published propulsion '
                         'closure (30 km/s, 70%, 45 degree cant, 300 W/kg, seven-day buffer).'),
               reading_rule=('published_check repeats aperture_quadrature.variable_distance of holding.json. Other cases '
                             'change only the base areal mass per zone and the annulus redirected fraction.'),
               published_mean_power_TW=published['mean_power_TW'], core_redirected_fraction=CORE_DIVERTED,
               cases=results, resources=dict(cpu_s=time.process_time()-cpu0, wall_s=time.monotonic()-wall))
    OUT.write_text(json.dumps(out, indent=1)+'\n')


if __name__ == '__main__':
    main()
