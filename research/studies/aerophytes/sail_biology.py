"""Requirements arithmetic for a living, controllable hanging sail.

This reads the saved sailing product; it does not solve flight, validate that
product's assumed lift/drag coefficients, or select a viable anatomy. Scenario
rates, dimensions and material fractions are choices, not measured reef traits.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from shared.constants import MOON_GM, MOON_RADIUS, SYNODIC_MONTH_DAYS, JULIAN_YEAR_DAYS

ROOT = Path(__file__).resolve().parents[3]
STUDY = Path('research/studies/aerophytes')
PARENT = STUDY / 'results/aerophytes.json'
EXPECTED_PARENT_SHA256 = '066cb18739272578220b264f1382652d8e57e86f215e625f69f102f0306b6fd3'


def _positive(**values):
    if any(not math.isfinite(v) or v <= 0 for v in values.values()):
        raise ValueError(f'finite positive values required: {tuple(values)}')


def _nonnegative(**values):
    if any(not math.isfinite(v) or v < 0 for v in values.values()):
        raise ValueError(f'finite nonnegative values required: {tuple(values)}')


def propagation(path_m, speed_m_s):
    """Travel only: no detection, synapse, processing or actuator time included."""
    _nonnegative(path=path_m)
    _positive(speed=speed_m_s)
    return dict(path_m=path_m, speed_m_s=speed_m_s,
                one_way_s=path_m/speed_m_s, round_trip_s=2*path_m/speed_m_s)


def potential_difference(depth_m, top_height_m, gm_m3_s2, radius_m):
    """Energy per kg raised from a fixed lower height to a fixed upper height."""
    _nonnegative(depth=depth_m, top_height=top_height_m)
    _positive(gm=gm_m3_s2, radius=radius_m)
    if depth_m > top_height_m:
        raise ValueError('the lower point must remain above sea level')
    upper = radius_m + top_height_m
    lower = upper - depth_m
    return gm_m3_s2*depth_m/(lower*upper)


def liquid_head(depth_m, top_height_m, density_kg_m3, gm_m3_s2, radius_m):
    """Static incompressible column; friction, osmotic and external-air heads omitted.

    Splitting this into pumping stages does not remove its total lifting work.
    """
    _positive(density=density_kg_m3)
    specific = potential_difference(depth_m, top_height_m, gm_m3_s2, radius_m)
    return dict(depth_m=depth_m, density_kg_m3=density_kg_m3,
                static_head_pa=density_kg_m3*specific, upward_work_j_per_kg=specific)


def retrieval(depth_m, paid_line_m, payload_kg, reel_speed_m_s,
              top_height_m, gm_m3_s2, radius_m):
    """Fixed-upper-height payload-only lifting floor and constant payout travel time.

    This excludes raising the line itself, drag, friction, actuator mass, storage,
    losses and changing reef altitude. It does not model a biological winch.
    """
    _nonnegative(payload=payload_kg)
    _positive(line=paid_line_m, reel_speed=reel_speed_m_s)
    if paid_line_m < depth_m:
        raise ValueError('paid line cannot be shorter than its vertical drop')
    seconds = paid_line_m/reel_speed_m_s
    energy = payload_kg*potential_difference(depth_m, top_height_m, gm_m3_s2, radius_m)
    return dict(depth_m=depth_m, paid_line_m=paid_line_m, payload_kg=payload_kg,
                reel_speed_m_s=reel_speed_m_s, travel_s=seconds,
                payload_gravitational_energy_j=energy,
                payload_mean_lifting_power_w=energy/seconds)


def panel_actuation(load_n_m2, area_m2, force_arm_m, tendon_arm_m, angle_deg, seconds):
    """Illustrative constant resisting torque about one hinge; not sail aerodynamics.

    load_n_m2 is an assumed resultant per area, not automatically dynamic pressure.
    Actual torque depends on the aerodynamic centre, geometry and deformation.
    Mechanical work excludes inertia, friction, holding and metabolic efficiency.
    """
    _nonnegative(load=load_n_m2, force_arm=force_arm_m, angle=angle_deg)
    _positive(area=area_m2, tendon_arm=tendon_arm_m, seconds=seconds)
    force = load_n_m2*area_m2
    torque = force*force_arm_m
    work = torque*math.radians(angle_deg)
    return dict(load_n_m2=load_n_m2, area_m2=area_m2, force_arm_m=force_arm_m,
                tendon_arm_m=tendon_arm_m, angle_deg=angle_deg, seconds=seconds,
                resultant_n=force, hinge_torque_nm=torque,
                tendon_force_n=torque/tendon_arm_m, mechanical_work_j=work,
                mechanical_power_w=work/seconds)


def equal_tether_split(count):
    """Same total section and length in independent, equally loaded round lines.

    Geometric projected-width ratio only: no wakes, new force balance, sheaths,
    joints, safety redundancy, taper or stronger bundle material included.
    """
    if isinstance(count, bool) or not isinstance(count, int) or count < 1:
        raise ValueError('count must be a positive integer')
    return dict(count=count, diameter_ratio_each=1/math.sqrt(count),
                total_section_ratio=1., total_projected_width_ratio=math.sqrt(count))


def wet_load(area_m2, water_kg_m2):
    """Additional retained liquid, separate from dry tissue and existing water."""
    _positive(area=area_m2)
    _nonnegative(water=water_kg_m2)
    return dict(water_kg_m2=water_kg_m2, added_water_kg=water_kg_m2*area_m2)


def renewal(area_m2, dry_kg_m2, replacements_per_year, carbon_fraction,
            phosphorus_fraction, phosphorus_recovery):
    """Material requirements only; no photosynthetic or metabolic budget is closed.

    New-tissue carbon is incorporated carbon, not total synthesis cost or new
    external carbon demand. P loss assumes the unrecovered material leaves the
    reef; material retained by its food web would instead remain internal.
    """
    _positive(area=area_m2)
    _nonnegative(dry_mass=dry_kg_m2, replacements=replacements_per_year)
    for v in (carbon_fraction, phosphorus_fraction, phosphorus_recovery):
        if not math.isfinite(v) or not 0 <= v <= 1:
            raise ValueError('material fractions must be in [0, 1]')
    if carbon_fraction + phosphorus_fraction > 1:
        raise ValueError('carbon and phosphorus cannot exceed the dry mass')
    stock = area_m2*dry_kg_m2
    annual = stock*replacements_per_year
    return dict(area_m2=area_m2, dry_kg_m2=dry_kg_m2,
                replacements_per_year=replacements_per_year,
                carbon_fraction=carbon_fraction, phosphorus_fraction=phosphorus_fraction,
                phosphorus_recovery=phosphorus_recovery, dry_inventory_kg=stock,
                new_dry_tissue_kg_year=annual,
                incorporated_carbon_kg_year=annual*carbon_fraction,
                lost_phosphorus_kg_year=annual*phosphorus_fraction*(1-phosphorus_recovery))


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def evaluate(root=ROOT, scenario_path=None, source_path=None):
    root = Path(root)
    scenario_path = Path(scenario_path or root/STUDY/'sail_biology_scenarios.json')
    source_path = Path(source_path or root/STUDY/'sail_biology_sources.json')
    actual_hash = _sha(root/PARENT)
    if actual_hash != EXPECTED_PARENT_SHA256:
        raise ValueError('parent product changed; review its evidence before rebinding')
    parent = json.loads((root/PARENT).read_text())
    if parent['schema'] != 'terluna.research.aerophytes/1':
        raise ValueError('unexpected parent product schema')
    s = json.loads(scenario_path.read_text())
    if s.get('schema') != 'terluna.research.aerophyte-sail-biology-scenarios/1':
        raise ValueError('unexpected scenario schema')
    gravity = dict(gm_m3_s2=MOON_GM, radius_m=MOON_RADIUS)
    area = math.pi*s['body_radius_m']**2
    wing_area = s['wing_share']*area
    rows = [r for r in parent['sailing']['tethered_sails']
            if r['body_radius_m'] == s['body_radius_m'] and r['wing_share'] == s['wing_share']
            and r['lower_km'] == s['sail_height_km'] and r['body_km'] == s['body_height_km']]
    if not rows:
        raise ValueError('no matching parent sailing cases')
    inherited = []
    for r in rows:
        inherited.append(dict(body_drag_coefficient=r['body_drag_coefficient'],
            allowable_mpa=r['allowable_mpa'], wing_force_n=r['wing_force_n'],
            mean_resultant_n_m2=r['wing_force_n']/wing_area,
            bottom_ballast_kg=r['ballast_kg_per_projected_m2']*area,
            tether_kg=r['tether_kg_per_projected_m2']*area,
            tether_diameter_m=r['tether_diameter_m'], paid_line_m=r['tether_length_m'],
            depth_m=r['wing_depth_m']))
    depth = (s['body_height_km']-s['sail_height_km'])*1000
    top_height = s['body_height_km']*1000
    loads = [min(r['mean_resultant_n_m2'] for r in inherited),
             max(r['mean_resultant_n_m2'] for r in inherited), s['stress_case_resultant_n_m2']]
    def name(path):
        try:
            return str(Path(path).relative_to(root))
        except ValueError:
            return str(path)
    inputs = [root/PARENT, scenario_path, source_path]
    return dict(schema='terluna.research.aerophyte-sail-biology/1',
        producer=dict(path=str(STUDY/'sail_biology.py'), sha256=_sha(__file__),
                      inputs={name(p): _sha(p) for p in inputs}),
        evidence='Requirements arithmetic and saved conditional sailing loads. No biological material, actuator, navigation, growth or survival validation.',
        reading_rule='The parent coefficients and force balance are not re-solved. Signal rates and anatomy are scenarios. Payload retrieval work is a lower bound at fixed reef altitude; panel work is a specified constant-torque example. New-tissue carbon is not a complete carbon cost. Recovery only changes export if lost material leaves the reef.',
        constants_used=dict(MOON_GM=MOON_GM, MOON_RADIUS=MOON_RADIUS,
                            SYNODIC_MONTH_DAYS=SYNODIC_MONTH_DAYS, JULIAN_YEAR_DAYS=JULIAN_YEAR_DAYS),
        scenario_limits=s['limits'],
        projected_area_m2=area, wing_area_m2=wing_area, inherited_cases=inherited,
        signals=[propagation(length, speed) for length in s['signal_path_m'] for speed in s['signal_speed_m_s']],
        liquid_columns=[liquid_head(d, top_height, s['liquid_density_kg_m3'], **gravity)
                        for d in [s['local_hydraulic_stage_m'], depth]],
        retrieval=[retrieval(r['depth_m'], r['paid_line_m'], r['bottom_ballast_kg'], speed,
                            top_height, **gravity)
                   for r in inherited for speed in s['reel_speed_m_s']],
        panel_actuation=[panel_actuation(load, **s['panel']) for load in loads],
        tether_split=[equal_tether_split(n) for n in s['tether_count']],
        wet_loads=[wet_load(wing_area, water) for water in s['water_loading_kg_m2']],
        renewal=[renewal(wing_area, dry_mass, rate, **s['tissue'])
                 for dry_mass in s['dry_sail_kg_m2']
                 for rate in [s['slow_replacement_per_year'], 1.,
                              JULIAN_YEAR_DAYS/SYNODIC_MONTH_DAYS]])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'research/runs/aerophytes/sail_biology.json')
    args = parser.parse_args()
    output = evaluate()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False)+'\n')
    print(args.output)


if __name__ == '__main__':
    main()
