"""Atmosphere-loading and sealed-construction screens; no forecast or cost model.

Pressure clocks assume retained gas and a specified globally averaged column.
Thin-column pressure is an inventory-equivalent bound, not a local plume model.
Shell outputs omit practical construction allowances and are not safe designs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from numbers import Real
from pathlib import Path

from shared.constants import (
    JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_GM, MOON_RADIUS,
    MOON_SURFACE_GRAVITY, STANDARD_GRAVITY,
)
from shared.provenance import ROOT, constants_used

HERE = Path(__file__).resolve().parent
INPUT_PATH = HERE / "industrial_physics_inputs.json"
INPUT_SCHEMA = "terluna.research.industrial-physics-inputs/1"
SCHEMA = "terluna.research.industrial-physics/1"


def _number(name, value, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be numeric")
    value = float(value)
    if not math.isfinite(value) or value < 0 or (positive and value == 0):
        raise ValueError(f"{name} must be finite and {'positive' if positive else 'nonnegative'}")
    return value


def _digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _finite_tree(value):
    if isinstance(value, dict):
        for item in value.values():
            _finite_tree(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            _finite_tree(item)
    elif isinstance(value, float) and not math.isfinite(value):
        raise ValueError("Nonfinite arithmetic result")


def load_inputs(path=INPUT_PATH, root=ROOT):
    data = json.loads(Path(path).read_text())
    if data.get("schema") != INPUT_SCHEMA:
        raise ValueError("Unexpected industrial physics input schema")
    selected = {}
    for key, spec in data["products"].items():
        source = Path(root) / spec["path"]
        if _digest(source) != spec["sha256"]:
            raise ValueError(f"Input changed; review and rebind explicitly: {source}")
        product = json.loads(source.read_text())
        if spec.get("schema") and product.get("schema") != spec["schema"]:
            raise ValueError(f"Unexpected source schema: {source}")
        values = {}
        for field in spec["fields"]:
            value = product
            for component in field.split("."):
                value = value[component]
            values[field] = value
        selected[key] = dict(**spec, values=values)
    _finite_tree(data)
    return data, selected


def loading_screen(radius_m, gravity_m_s2, final_mass_kg,
                   reference_surface_pressure_Pa, end_year, start_year,
                   thresholds_Pa, year_s):
    """Two linear retained-mass clocks, both ending at end_year.

    The calibrated clock preserves the reference column's pressure/mass ratio.
    The thin clock takes all mass at surface gravity/area. For a spherical,
    hydrostatic atmosphere it is an upper pressure bound at fixed total mass.
    Neither predicts local pressure during a release or the exospheric phase.
    """
    radius = _number("radius", radius_m, positive=True)
    gravity = _number("gravity", gravity_m_s2, positive=True)
    final_mass = _number("final mass", final_mass_kg, positive=True)
    reference_pressure = _number("reference pressure", reference_surface_pressure_Pa, positive=True)
    end = _number("end year", end_year, positive=True)
    start = _number("start year", start_year)
    seconds = _number("seconds per year", year_s, positive=True)
    if start >= end:
        raise ValueError("Loading must start before its end")
    thresholds = [_number("threshold", p, positive=True) for p in thresholds_Pa]
    if len(set(thresholds)) != len(thresholds) or not thresholds:
        raise ValueError("Pressure thresholds must be nonempty and unique")
    area = 4 * math.pi * radius**2
    thin_mass_per_pa = area / gravity
    thin_final_pressure = final_mass / thin_mass_per_pa
    if reference_pressure > thin_final_pressure * (1 + 1e-12):
        raise ValueError("Reference pressure exceeds the spherical hydrostatic mass bound")
    calibrated_mass_per_pa = final_mass / reference_pressure
    duration = end - start
    rows = []
    for pressure in thresholds:
        row = dict(threshold_Pa=pressure)
        for name, mass_per_pa in (
            ("thin_column", thin_mass_per_pa),
            ("reference_column_scaled", calibrated_mass_per_pa),
        ):
            mass = pressure * mass_per_pa
            fraction = mass / final_mass
            reached = fraction <= 1
            elapsed = duration * fraction if reached else None
            row[name] = dict(
                retained_mass_kg=mass, fraction_of_final_mass=fraction,
                reached_by_end=reached,
                elapsed_loading_years=elapsed,
                elapsed_loading_seconds=elapsed * seconds if reached else None,
                crossing_year=start + elapsed if reached else None,
            )
        rows.append(row)
    result = dict(
        start_year=start, end_year=end, loading_years=duration,
        final_mass_kg=final_mass, retained_loading_kg_s=final_mass/(duration*seconds),
        surface_area_m2=area, gravity_m_s2=gravity,
        thin_mass_per_Pa_kg=thin_mass_per_pa,
        reference_mass_per_Pa_kg=calibrated_mass_per_pa,
        thin_equivalent_final_pressure_Pa=thin_final_pressure,
        reference_final_pressure_Pa=reference_pressure,
        spherical_mass_factor=calibrated_mass_per_pa/thin_mass_per_pa,
        thresholds=rows,
    )
    _finite_tree(result)
    return result


def tensile_shell_screen(radius_m, internal_pressure_Pa, external_pressure_Pa,
                         density_kg_m3, allowable_stress_Pa):
    """Ideal isotropic membrane: full sphere and long-cylinder barrel only.

    Zero differential gives zero pressure-only material, not a massless building.
    An external-overpressure case requires a separate buckling calculation.
    """
    radius = _number("radius", radius_m, positive=True)
    inside = _number("internal pressure", internal_pressure_Pa)
    outside = _number("external pressure", external_pressure_Pa)
    density = _number("density", density_kg_m3, positive=True)
    stress = _number("allowable stress", allowable_stress_Pa, positive=True)
    if outside > inside:
        raise ValueError("External overpressure is outside the tensile-shell screen")
    delta = inside - outside
    sphere_t = delta * radius / (2 * stress)
    cylinder_t = delta * radius / stress
    volume = 4 * math.pi * radius**3 / 3
    result = dict(
        radius_m=radius, internal_pressure_Pa=inside, external_pressure_Pa=outside,
        differential_pressure_Pa=delta,
        sphere=dict(thickness_m=sphere_t,
                    shell_mass_kg=4*math.pi*radius**2*sphere_t*density,
                    enclosed_volume_m3=volume,
                    shell_mass_per_volume_kg_m3=1.5*density*delta/stress),
        long_cylinder_barrel=dict(thickness_m=cylinder_t,
                    shell_mass_per_barrel_volume_kg_m3=2*density*delta/stress),
        thin_wall_ratio_sphere=sphere_t/radius,
        thin_wall_ratio_cylinder=cylinder_t/radius,
        thin_wall_approximation_applicable=cylinder_t/radius <= 0.05,
    )
    _finite_tree(result)
    return result


def external_sphere_screen(radius_m, differential_pressure_Pa, density_kg_m3,
                           elastic_modulus_Pa, poisson_ratio, knockdown_factor):
    """Classical perfect-sphere elastic buckling and an assumed reduction.

    The assumed reduction is a sensitivity, not a safety factor or a code check.
    Seams, ports, imperfections, supports, plasticity and other loads are omitted.
    """
    radius = _number("radius", radius_m, positive=True)
    pressure = _number("external differential", differential_pressure_Pa)
    density = _number("density", density_kg_m3, positive=True)
    modulus = _number("elastic modulus", elastic_modulus_Pa, positive=True)
    nu = _number("Poisson ratio", poisson_ratio)
    factor = _number("knockdown factor", knockdown_factor, positive=True)
    if nu >= 0.5 or factor > 1:
        raise ValueError("Require 0 <= Poisson ratio < 0.5 and 0 < knockdown <= 1")
    perfect_ratio = math.sqrt(pressure*math.sqrt(3*(1-nu*nu))/(2*modulus))
    adjusted_ratio = perfect_ratio / math.sqrt(factor)
    result = dict(radius_m=radius, differential_pressure_Pa=pressure,
        assumed_knockdown_factor=factor,
        perfect_shell_thickness_m=radius*perfect_ratio,
        assumed_reduced_capacity_thickness_m=radius*adjusted_ratio,
        perfect_shell_mass_per_volume_kg_m3=3*density*perfect_ratio,
        assumed_reduced_capacity_mass_per_volume_kg_m3=3*density*adjusted_ratio,
        ideal_membrane_compressive_stress_Pa=pressure/(2*adjusted_ratio) if pressure else 0,
        limitation="Elastic equilibrium/buckling geometry only; neither thickness is a construction prescription.")
    _finite_tree(result)
    return result


def handling_screen(mass_kg, gravity_m_s2, reference_gravity_m_s2,
                    lift_height_m, acceleration_m_s2, speed_m_s, friction_coefficient):
    mass = _number("mass", mass_kg, positive=True)
    gravity = _number("gravity", gravity_m_s2, positive=True)
    reference = _number("reference gravity", reference_gravity_m_s2, positive=True)
    height = _number("lift height", lift_height_m)
    acceleration = _number("acceleration", acceleration_m_s2)
    speed = _number("speed", speed_m_s)
    friction = _number("friction coefficient", friction_coefficient)
    result = dict(mass_kg=mass, lunar_weight_N=mass*gravity,
        reference_weight_N=mass*reference, weight_ratio=gravity/reference,
        lunar_lift_energy_J=mass*gravity*height,
        reference_lift_energy_J=mass*reference*height,
        inertial_force_both_N=mass*acceleration,
        kinetic_energy_both_J=0.5*mass*speed**2,
        lunar_weight_only_friction_limit_N=friction*mass*gravity,
        reference_weight_only_friction_limit_N=friction*mass*reference)
    _finite_tree(result)
    return result


def ring_energy_screen(gm_m3_s2, surface_radius_m, circular_orbit_radius_m):
    gm = _number("GM", gm_m3_s2, positive=True)
    surface = _number("surface radius", surface_radius_m, positive=True)
    orbit = _number("orbit radius", circular_orbit_radius_m, positive=True)
    if orbit <= surface:
        raise ValueError("Orbit must be above the surface")
    potential = gm*(1/surface - 1/orbit)
    kinetic = gm/(2*orbit)
    result = dict(circular_orbit_radius_m=orbit,
        potential_energy_increase_J_kg=potential,
        final_circular_kinetic_energy_J_kg=kinetic,
        mechanical_energy_increase_J_kg=potential+kinetic,
        surface_escape_energy_J_kg=gm/surface,
        limitation="Starts stationary in a nonrotating lunar frame; excludes drag, propulsion losses, launch/capture hardware, traffic geometry and Earth/Sun perturbations. This energy is not an electricity bill or export revenue.")
    _finite_tree(result)
    return result


def beam_span_screen(gravity_m_s2, reference_gravity_m_s2):
    """Same beam section/material and uniform self-weight alone.

    Absolute deflection is proportional to g*L^4; bending stress to g*L^2.
    A relative-deflection criterion, live load or other constraint changes this.
    """
    gravity = _number("gravity", gravity_m_s2, positive=True)
    reference = _number("reference gravity", reference_gravity_m_s2, positive=True)
    ratio = reference/gravity
    result = dict(equal_absolute_sag_span_ratio=ratio**0.25,
                  equal_bending_stress_span_ratio=ratio**0.5,
                  load_scope="Uniform self-weight only; fixed section, material and boundary type. No machinery, pressure, thermal, dynamic, buckling or whole-building cost inference.")
    _finite_tree(result)
    return result


def capillary_screen(liquid_density_kg_m3, surface_tension_N_m, gravity_m_s2,
                     lengths_m, head_height_m):
    density = _number("liquid density", liquid_density_kg_m3, positive=True)
    tension = _number("surface tension", surface_tension_N_m, positive=True)
    gravity = _number("gravity", gravity_m_s2, positive=True)
    head = _number("head height", head_height_m)
    lengths = [_number("length", value, positive=True) for value in lengths_m]
    if not lengths or len(set(lengths)) != len(lengths):
        raise ValueError("Capillary lengths must be nonempty and unique")
    result = dict(capillary_length_m=math.sqrt(tension/(density*gravity)),
                  hydrostatic_head_Pa=density*gravity*head,
                  bond_numbers=[dict(length_m=value, bond_number=density*gravity*value**2/tension)
                                for value in lengths],
                  limitation="Chosen liquid and geometry; Bond number compares gravity and capillarity, not wetting, flow, finished quality or a validated process window.")
    _finite_tree(result)
    return result


def ideal_pump_work(volume_m3, initial_pressure_Pa, final_pressure_Pa):
    """Reversible isothermal evacuation of rigid volume to initial ambient P.

    Real pumps, leaks, outgassing and conductance require more work/time.
    This is not an estimate of installed pumping equipment or consumption.
    """
    volume = _number("volume", volume_m3, positive=True)
    initial = _number("initial pressure", initial_pressure_Pa, positive=True)
    final = _number("final pressure", final_pressure_Pa)
    if final > initial:
        raise ValueError("Evacuation final pressure must not exceed initial pressure")
    fractional_drop = (initial-final)/initial
    if final == 0:
        work = volume*initial
    elif fractional_drop < 0.001:
        # Avoid subtracting nearly equal terms for a very small pressure drop.
        d = fractional_drop
        work = volume*initial*d*d*(0.5+d*(1/6+d*(1/12+d/20)))
    else:
        work = volume*(initial-final*(1+math.log(initial/final)))
    result = dict(volume_m3=volume, initial_pressure_Pa=initial, final_pressure_Pa=final,
                  ideal_reversible_work_J=work)
    _finite_tree(result)
    return result


def batch_pumpdown_screen(process_minutes, pumpdown_minutes):
    process = _number("useful batch process time", process_minutes, positive=True)
    pumping = _number("pumpdown time", pumpdown_minutes)
    result = dict(process_minutes=process, pumpdown_minutes=pumping,
        original_cycle_minutes=process+pumping, without_pumpdown_minutes=process,
        ideal_capacity_ratio_if_all_pumpdown_removed=(process+pumping)/process,
        original_cycle_fraction_removed=pumping/(process+pumping),
        limitation="One repeatedly vented chamber, sequential pumpdown and process. Omits loading, parallel chambers, load locks, continuous operation and local heads; no whole-factory productivity gain.")
    _finite_tree(result)
    return result


def deposition_head_screen(volume_rate_m3_h, density_kg_m3, duty_fraction,
                           annual_target_kg, year_s):
    rate = _number("volume rate", volume_rate_m3_h, positive=True)
    density = _number("density", density_kg_m3, positive=True)
    duty = _number("duty fraction", duty_fraction, positive=True)
    target = _number("annual target", annual_target_kg, positive=True)
    seconds = _number("seconds per year", year_s, positive=True)
    if duty > 1:
        raise ValueError("Duty fraction must be <= 1")
    hourly = rate*density
    annual = hourly*duty*seconds/3600
    heads = target/annual
    result = dict(volume_rate_m3_h=rate, mass_rate_kg_h=hourly,
        assumed_duty_fraction=duty, annual_target_kg=target,
        deposited_mass_per_head_kg_year=annual,
        continuous_head_equivalents=heads, minimum_integer_heads=math.ceil(heads),
        limitation="Historical demonstrated deposition rate and chosen utilization/target. Before rejection, finishing, feedstock production and support equipment; neither a future maximum rate nor saleable-product capacity.")
    _finite_tree(result)
    return result


def accelerator_screen(speed_m_s, acceleration_m_s2):
    speed = _number("launch speed", speed_m_s)
    acceleration = _number("acceleration", acceleration_m_s2, positive=True)
    result = dict(speed_m_s=speed, acceleration_m_s2=acceleration,
        constant_acceleration_length_m=speed**2/(2*acceleration),
        acceleration_duration_s=speed/acceleration,
        limitation="Starts at rest and neglects gravity along the track, drag, trajectory shaping, packaging and capture. No assumption that fragile components survive this acceleration.")
    _finite_tree(result)
    return result


def results(input_path=INPUT_PATH, root=ROOT):
    data, products = load_inputs(input_path, root)
    mass = products["air_inventory"]["values"]["inventories.air_kg"]
    if mass != products["historical_air"]["values"]["inventory.mass_kg"]:
        raise ValueError("Pressure calibration and air inventory do not share a final mass")
    pressure = products["historical_air"]["values"]["constants.pressure_Pa"]
    ring_radius = products["ring_layout"]["values"]["design.middle_radius_km"] * 1000
    schedule = data["loading_scenarios"]
    material = data["illustrative_material"]
    enclosures = data["enclosures"]
    liquid = data["illustrative_liquid"]
    pump = data["pumpdown_example"]
    deposition = data["historical_deposition_scale"]
    launch = data["illustrative_accelerator"]
    shells = []
    for inside in enclosures["internal_pressures_Pa"]:
        for radius in enclosures["radii_m"]:
            shells.append(dict(
                against_vacuum=tensile_shell_screen(radius, inside, 0,
                    material["density_kg_m3"], material["allowable_stress_Pa"]),
                matched_ambient=tensile_shell_screen(radius, inside, inside,
                    material["density_kg_m3"], material["allowable_stress_Pa"]),
                vacuum_chamber_after_air=external_sphere_screen(radius, inside,
                    material["density_kg_m3"], material["elastic_modulus_Pa"],
                    material["poisson_ratio"], enclosures["illustrative_buckling_knockdown"]),
            ))
    files = [Path(__file__), Path(input_path)]
    def relative(path):
        try:
            return str(path.resolve().relative_to(Path(root).resolve()))
        except ValueError:
            return path.name
    output = dict(schema=SCHEMA,
        producer=dict(study="world_economy", files={relative(p):_digest(p) for p in files},
                      constants=constants_used([Path(__file__)])),
        evidence=data["evidence"], sources=data["sources"], products=products,
        assumptions=data["assumptions"], reading_rule=data["reading_rule"],
        loading=[loading_screen(MOON_RADIUS, MOON_SURFACE_GRAVITY, mass, pressure,
            schedule["end_year"], start, schedule["thresholds_Pa"], JULIAN_DAY*JULIAN_YEAR_DAYS)
            for start in schedule["start_years"]],
        material_inputs=material, enclosure_cases=shells,
        process_inputs={key:data[key] for key in ("illustrative_liquid", "pumpdown_example",
            "historical_deposition_scale", "illustrative_accelerator", "beam_span_assumptions")},
        uplift_ballast=[dict(differential_pressure_Pa=p,
            lunar_deadweight_kg_m2=p/MOON_SURFACE_GRAVITY,
            earth_deadweight_kg_m2=p/STANDARD_GRAVITY,
            lunar_depth_m=p/(MOON_SURFACE_GRAVITY*enclosures["illustrative_ballast_density_kg_m3"]),
            scope="Roof uplift balanced solely by deadweight per projected area; not the net load on a fully closed hull.")
            for p in enclosures["internal_pressures_Pa"]],
        handling=handling_screen(gravity_m_s2=MOON_SURFACE_GRAVITY,
            reference_gravity_m_s2=STANDARD_GRAVITY, **data["handling"]),
        ring_transfer=ring_energy_screen(MOON_GM, MOON_RADIUS, ring_radius),
        beam_spans=beam_span_screen(MOON_SURFACE_GRAVITY, STANDARD_GRAVITY),
        capillary={name:capillary_screen(liquid["liquid_density_kg_m3"],
            liquid["surface_tension_N_m"], gravity, liquid["lengths_m"], liquid["head_height_m"])
            for name, gravity in (("moon", MOON_SURFACE_GRAVITY), ("earth", STANDARD_GRAVITY))},
        ideal_pumping=ideal_pump_work(pump["chamber_volume_m3"], pump["initial_pressure_Pa"],
                                    pump["ideal_work_final_pressure_Pa"]),
        batch_pumpdown=[dict(source_working_pressure_Pa=case["working_pressure_Pa"],
            **batch_pumpdown_screen(process, case["pumpdown_minutes"]))
            for process in pump["illustrative_process_minutes"] for case in pump["reported_pumpdown_cases"]],
        deposition_heads=[deposition_head_screen(rate, deposition["illustrative_density_kg_m3"],
            deposition["assumed_duty_fraction"], deposition["illustrative_annual_target_kg"],
            JULIAN_DAY*JULIAN_YEAR_DAYS) for rate in deposition["reported_volume_rates_m3_h"]],
        accelerators=[accelerator_screen(launch["speed_m_s"], acceleration)
                      for acceleration in launch["accelerations_m_s2"]],
    )
    _finite_tree(output)
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", type=Path, default=INPUT_PATH)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path, default=ROOT / "research/runs/world_economy/industrial_physics.json")
    args = parser.parse_args()
    product = results(args.inputs, args.root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(product, indent=2, allow_nan=False)+"\n")
    print(args.output)


if __name__ == "__main__":
    main()
