"""Area quadrature of independently held tiles, including per-tile mass closure."""
from __future__ import annotations
import numpy as np
from shared import constants as K
from protection.dynamics.ephemeris import BODY_GM
from protection.dynamics.model import target, relative_gravity
from protection.dynamics.optical import length, covering_radius, solar_visibility, optical_projection, axial_optical, light_time_allowance, arriving_ray_vectors
from .run import SCENARIO, CORE_DIVERTED


def disk_nodes(radius, core_radius, rings=4, angles=12):
    """Equal-area radial midpoints, with the central climate window resolved."""
    coords, weights = [], []
    for lo, hi in ((0, core_radius), (core_radius, radius)):
        if hi <= lo:
            continue
        for j in range(rings):
            r = np.sqrt(lo*lo+(hi*hi-lo*lo)*(j+0.5)/rings)
            for angle in (np.arange(angles)+0.5)*2*np.pi/angles:
                coords.append([r*np.cos(angle), r*np.sin(angle)])
                weights.append((hi*hi-lo*lo)/(radius*radius*rings*angles))
    return np.array(coords), np.array(weights)


def evaluate_aperture(samples, geo, coefficients, mode="ideal", base_sigma=0.05,
                      protected_radii=4, rings=4, angles=12):
    trajectory = target(geo, coefficients)
    distance = trajectory["coords"][:, 0]
    offset = length(trajectory["coords"][:, 1:])
    sun_distance = length(samples["positions"]["sun"])
    ray_margin=light_time_allowance(distance,samples["moon_v"],samples["sun_v"])
    radius = float((covering_radius(distance, sun_distance, protected_radii*K.MOON_RADIUS, offset)+ray_margin).max()+SCENARIO["formation_margin_m"])
    core_radii = covering_radius(distance, sun_distance, K.MOON_RADIUS)+ray_margin
    core_radius = min(float((core_radii+offset).max()), radius)
    nodes, weights = disk_nodes(radius, core_radius, rings, angles)
    displacement = np.einsum("pi,nij->npj", nodes, geo["frame"][:, 1:])
    q = trajectory["q"][:, None, :]+displacement
    acc = trajectory["a"][:, None, :]+np.einsum("pi,nij->npj", nodes, geo["frame_dd"][:, 1:])
    positions = {k:v[:, None, :] for k,v in samples["positions"].items()}
    required = acc-relative_gravity(q, positions, samples["moon_a"][:, None, :])
    sun,earth=arriving_ray_vectors(positions["sun"]-q,positions["earth"]-q,samples["sun_v"][:,None,:],samples["earth_v"][:,None,:],samples["sun_a"][:,None,:],samples["earth_a"][:,None,:])
    visibility = solar_visibility(sun, earth)
    absolute_xy = trajectory["coords"][:, None, 1:]+nodes[None, :, :]
    # Any tile whose rays can reach the solid Moon must preserve the climate spectrum.
    diverted = np.where(length(absolute_xy) <= core_radii[:, None], CORE_DIVERTED, 1.0)
    prop = SCENARIO["propulsion"]
    ve, eta, kappa = prop["exhaust_velocity_m_s"], prop["efficiency"], prop["specific_power_W_kg"]
    cant = np.cos(np.deg2rad(prop["cant_deg"]))
    buffer = prop["propellant_buffer_days"]*K.JULIAN_DAY
    sigma = np.full(len(nodes), base_sigma)
    t = samples["t"]
    for iteration in range(100):
        if mode == "none":
            residual = required
        elif mode == "ideal":
            residual = optical_projection(required, sun, sigma, diverted, visibility)[1]
        elif mode in ("absorb", "reflect"):
            coefficient = (1 if mode == "absorb" else 2)*diverted
            residual = required-axial_optical(sun, sigma, coefficient, visibility)
        else:
            raise ValueError(mode)
        mag = length(residual)
        mean = np.trapezoid(mag, t, axis=0)/(t[-1]-t[0])
        peak = mag.max(axis=0)*prop["peak_margin_factor"]
        mdot_area = sigma*mean/(ve*cant)
        power_area = mdot_area*ve*ve/(2*eta)
        peak_area = sigma*peak*ve/(2*eta*cant)
        new_sigma = base_sigma+peak_area/kappa+mdot_area*buffer
        if np.max(abs(new_sigma-sigma)) < 1e-11:
            break
        sigma = new_sigma
    else:
        raise RuntimeError("Per-tile mass closure did not converge")
    area = np.pi*radius*radius
    mean_required = np.trapezoid(length(required)@weights, t)/(t[-1]-t[0])
    mean_residual = mean@weights
    return {"optical_mode": mode, "quadrature_nodes": len(nodes), "rings_per_zone": rings,
            "azimuths": angles, "aperture_radius_km": radius/1000,
            "aperture_area_m2": area, "protected_radii": protected_radii,
            "mean_required_mm_s2": float(mean_required*1000),
            "area_mean_residual_mm_s2": float(mean_residual*1000),
            "mean_power_TW": float(area*(power_area@weights)/1e12),
            "installed_design_peak_power_TW": float(area*(peak_area@weights)/1e12),
            "propellant_kg_s": float(area*(mdot_area@weights)),
            "gyr_propellant_kg": float(area*(mdot_area@weights)*K.JULIAN_DAY*K.JULIAN_YEAR_DAYS*1e9),
            "total_mass_kg": float(area*(sigma@weights)),
            "base_optical_mass_kg": area*base_sigma,
            "power_hardware_mass_kg": float(area*(peak_area@weights)/kappa),
            "buffer_mass_kg": float(area*(mdot_area@weights)*buffer),
            "areal_mass_min_max_kg_m2": [float(sigma.min()), float(sigma.max())],
            "solar_reduction_percent": float((1-mean_residual/mean_required)*100),
            "storage_included": False,
            "ray_flight_allowance_max_km":float(ray_margin.max()/1000),
            "coverage_rule": "Worst-epoch full-Sun aperture plus source/target ray-flight motion allowance and 50 m margin; independently held area nodes. Central rays retain the climate spectrum."}
