"""Conditional partial-gravity process screens, not demonstrated product gains.

The casting model uses isolated-sphere Stokes motion at constant viscosity.
Published material inputs and hypothetical process constraints remain separate.
The spin model gives a rotating liquid's hydrostatic paraboloid and load ratios;
it does not predict cured optical quality, equipment costs or export revenue.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from numbers import Real
from pathlib import Path

from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
INPUT_PATH = HERE / "partial_gravity_inputs.json"
INPUT_SCHEMA = "terluna.research.partial-gravity-inputs/1"


def _number(name, value, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite number")
    value = float(value)
    if not math.isfinite(value) or value < 0 or (positive and value == 0):
        raise ValueError(f"{name} must be finite and {'positive' if positive else 'nonnegative'}")
    return value


def _digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_inputs(path=INPUT_PATH):
    data = json.loads(Path(path).read_text())
    if data.get("schema") != INPUT_SCHEMA:
        raise ValueError("Unexpected partial-gravity input schema")
    return data


def stokes_speed(density_contrast, radius, viscosity, gravity):
    """Magnitude for a dilute rigid sphere or immobilized spherical bubble."""
    drho = _number("density contrast", density_contrast)
    a = _number("radius", radius, positive=True)
    mu = _number("viscosity", viscosity, positive=True)
    g = _number("gravity", gravity)
    result = 2 * drho * g * a * a / (9 * mu)
    if not math.isfinite(result):
        raise ValueError("Stokes speed exceeds finite arithmetic")
    return result


def _parameters(material, process):
    rho = _number("fluid density", material["fluid_density_kg_m3"], positive=True)
    rho_p = _number("particle density", material["particle_density_kg_m3"], positive=True)
    if rho_p <= rho:
        raise ValueError("This settling screen requires particles denser than the liquid")
    mu = _number("dynamic viscosity", material["dynamic_viscosity_Pa_s"], positive=True)
    particle_radius = _number("particle diameter", material["particle_diameter_m"], positive=True) / 2
    bubble_radius = _number("bubble radius", process["bubble_radius_m"], positive=True)
    depth = _number("escape depth", process["escape_depth_m"], positive=True)
    drift = _number("particle drift tolerance", process["particle_drift_limit_m"], positive=True)
    duration = _number("duration", process["liquid_interval_s"], positive=True)
    sigma = _number("surface tension", process["surface_tension_N_m"], positive=True)
    return rho, rho_p, mu, particle_radius, bubble_radius, depth, drift, duration, sigma


def casting_screen(material, process, gravity, *, duration_s=None):
    """Free-surface arrival and particle drift; arrival does not prove venting.

    Gas density is neglected. Both phases use an immobile-interface Stokes law.
    The diameter-based Bond number is rho*g*(2*a_bubble)^2/sigma.
    Ideal zero gravity is a reference limit, not a residual-acceleration model.
    """
    rho, rho_p, mu, ap, ab, depth, drift, duration, sigma = _parameters(material, process)
    g = _number("gravity", gravity)
    if duration_s is not None:
        duration = _number("duration", duration_s, positive=True)
    vp = stokes_speed(rho_p - rho, ap, mu, g)
    vb = stokes_speed(rho, ab, mu, g)
    particle_motion, bubble_motion = vp * duration, vb * duration
    bubble_pass = bubble_motion >= depth
    particle_pass = particle_motion <= drift
    return dict(
        gravity_m_s2=g, liquid_interval_s=duration,
        particle_speed_m_s=vp, bubble_speed_m_s=vb,
        particle_displacement_m=particle_motion, bubble_displacement_m=bubble_motion,
        bubble_surface_arrival_time_s=depth / vb if vb else None,
        particle_drift_limit_time_s=drift / vp if vp else None,
        latest_bubble_birth_time_s=duration - depth / vb if vb else None,
        particle_reynolds=rho * vp * (2 * ap) / mu,
        bubble_reynolds=rho * vb * (2 * ab) / mu,
        bubble_bond_diameter=rho * g * (2 * ab) ** 2 / sigma,
        bubble_capillary_number=mu * vb / sigma,
        bubble_reaches_surface_before_immobilization=bubble_pass,
        particle_drift_within_limit=particle_pass,
        both_constraints_pass=bubble_pass and particle_pass,
    )


def gravity_window(material, process):
    """Gravity interval for one specified, constant-viscosity liquid interval."""
    _, _, _, _, _, depth, drift, duration, _ = _parameters(material, process)
    unit = casting_screen(material, process, 1.0)
    lower = depth / (unit["bubble_speed_m_s"] * duration)
    upper = drift / (unit["particle_speed_m_s"] * duration)
    return dict(
        minimum_for_bubble_arrival_m_s2=lower,
        maximum_for_particle_uniformity_m_s2=upper,
        interval_nonempty=lower <= upper,
        bubble_to_particle_speed_ratio=unit["bubble_speed_m_s"] / unit["particle_speed_m_s"],
        minimum_required_displacement_ratio=depth / drift,
    )


def spin_paraboloid(gravity, focal_length, aperture_radius):
    """z(r)=omega^2*r^2/(2*g)=r^2/(4*f) for uniform axial gravity.

    The speed-squared load proxy compares identical geometry and material.
    It omits hydrostatic, thermal, capillary, wall and solidification stresses.
    No finite-focus gravitational paraboloid is defined at zero axial gravity.
    """
    g = _number("gravity", gravity)
    f = _number("focal length", focal_length, positive=True)
    radius = _number("aperture radius", aperture_radius, positive=True)
    if g == 0:
        return dict(gravity_m_s2=g, hydrostatic_paraboloid_defined=False,
                    angular_speed_rad_s=None, rotational_speed_rpm=None,
                    rim_speed_m_s=None, rim_centripetal_acceleration_m_s2=None,
                    rotation_load_proxy_m2_s2=None, edge_sag_m=None)
    omega = math.sqrt(g / (2 * f))
    return dict(
        gravity_m_s2=g, hydrostatic_paraboloid_defined=True,
        angular_speed_rad_s=omega, rotational_speed_rpm=omega * 60 / (2 * math.pi),
        rim_speed_m_s=omega * radius,
        rim_centripetal_acceleration_m_s2=omega * omega * radius,
        rotation_load_proxy_m2_s2=(omega * radius) ** 2,
        edge_sag_m=radius * radius / (4 * f),
    )


def calculate(settings):
    material = settings["paper_tabulated_material_inputs"]["values"]
    process = settings["illustrative_casting_process"]["values"]
    gravities = dict(earth=STANDARD_GRAVITY, moon=MOON_SURFACE_GRAVITY,
                     ideal_microgravity=0.0)
    casting = {name: casting_screen(material, process, g) for name, g in gravities.items()}
    duration = process["liquid_interval_s"]
    matched_duration = duration * MOON_SURFACE_GRAVITY / STANDARD_GRAVITY
    matched_mu = material["dynamic_viscosity_Pa_s"] * STANDARD_GRAVITY / MOON_SURFACE_GRAVITY
    altered_material = dict(material, dynamic_viscosity_Pa_s=matched_mu)
    long_duration = settings["long_interval_sensitivity"]["duration_s"]
    long_process = dict(process, liquid_interval_s=long_duration)
    spin_inputs = settings["illustrative_spin_cast_geometry"]
    spin = {name: spin_paraboloid(g, spin_inputs["focal_length_m"],
                                 spin_inputs["aperture_radius_m"])
            for name, g in gravities.items()}
    return dict(
        gravity_window=gravity_window(material, process), casting=casting,
        earth_matching_controls=dict(
            same_viscosity_immobilize_after_s=matched_duration,
            time_control=casting_screen(material, process, STANDARD_GRAVITY,
                                        duration_s=matched_duration),
            same_duration_dynamic_viscosity_Pa_s=matched_mu,
            viscosity_control=casting_screen(altered_material, process, STANDARD_GRAVITY),
            condition="Both are mathematical controls. Earlier immobilization or a viscosity change must preserve filling, chemistry and final properties; no feasible recipe is established.",
        ),
        long_interval_sensitivity=dict(
            duration_s=long_duration,
            gravity_window=gravity_window(material, long_process),
            cases={name: casting_screen(material, long_process, g) for name, g in gravities.items()},
            evidence=settings["long_interval_sensitivity"]["evidence"],
        ),
        spin_paraboloid=dict(
            geometry=spin_inputs, cases=spin,
            moon_over_earth_angular_speed=math.sqrt(MOON_SURFACE_GRAVITY / STANDARD_GRAVITY),
            moon_over_earth_inertial_load_proxy=MOON_SURFACE_GRAVITY / STANDARD_GRAVITY,
            limitation="Same geometry/material comparison only; no cured optical quality, complete stress budget or manufacturing-cost benefit established.",
        ),
    )


def results(input_path=INPUT_PATH):
    settings = load_inputs(input_path)
    return dict(
        schema="terluna.research.partial-gravity/1",
        producer=dict(
            study="world_economy",
            files={"research/studies/world_economy/partial_gravity.py": _digest(__file__),
                   "research/studies/world_economy/partial_gravity_inputs.json": _digest(input_path)},
            constants=constants_used([__file__]),
        ),
        evidence=settings["evidence"], sources=settings["sources"],
        material_input_scope=settings["paper_tabulated_material_inputs"]["evidence"],
        hypothetical_process_scope=settings["illustrative_casting_process"]["evidence"],
        reading_rule="SI units. Ideal constant-viscosity Stokes and hydrostatic screens; null times mean no buoyant arrival at zero gravity. A fixed-cycle pass is not a demonstrated advantage over an optimized terrestrial or orbital process. No economic premium is assigned.",
        assumptions=settings["assumptions"], **calculate(settings),
    )


def main():
    root = next((p for p in HERE.parents if (p / "shared/constants.json").is_file()), None)
    default = (root / "research/runs/world_economy/partial_gravity.json"
               if root else HERE / "partial_gravity_result.json")
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=default)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results(), indent=2, allow_nan=False) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
