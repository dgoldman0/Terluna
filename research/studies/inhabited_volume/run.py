"""Inhabited-volume requirements; closed-form screens, never a capacity forecast.

Run with python -m research.studies.inhabited_volume.run. Cross-lane inputs are
schema-checked products. Historical population.py and its output stay intact.
"""
from __future__ import annotations
import hashlib
import json
import math
from pathlib import Path
from shared.constants import MOON_GM, MOON_RADIUS, JULIAN_YEAR_DAYS, JULIAN_DAY
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.inhabited-volume/1'
OUT = HERE / 'results/inhabited_volume.json'
INPUTS = {
    'shared/scenarios/population.json': 'terluna.scenario.population/1',
    'research/studies/sky_ships/results/sky_ships.json': 'terluna.research.sky-ships/1',
    'provisioning/results/population.json': 'terluna.provisioning.population/1',
}
# Hypothetical district payloads, separate from envelope and flight hardware.
# Wet stocks include stored water/food; buildings include floors, partitions and furnishings support.
DISTRICTS = {
    'light20': dict(building_t=10.5, people_furnishings_t=2.5, wet_stocks_t=3., utilities_local_stocks_t=4.),
    'central50': dict(building_t=28., people_furnishings_t=4., wet_stocks_t=8., utilities_local_stocks_t=10.),
    'heavy100': dict(building_t=59.5, people_furnishings_t=5.5, wet_stocks_t=15., utilities_local_stocks_t=20.),
}
LIFT_SHARES = dict(district_payload=.5, envelope_flight_structure=.25, control_flight_services=.1, reserve=.15)
GAS_PURITY = .98  # sky_ships/run.py PURITY, already pinned by the consumed air product
ROAMING_DUTIES = (.01,.05,.10)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def load_product(path, schema):
    data = json.loads(Path(path).read_text())
    if data.get('schema') != schema:
        raise ValueError(f'Unexpected schema for {path}')
    return data


def validate_cases(shared):
    for name, case in shared['cases'].items():
        parts = [case[k] for k in ('earth', 'lunar_surface', 'lunar_aerial', 'orbital', 'elsewhere')]
        if min(parts) < 0 or not math.isclose(sum(parts), case['total']):
            raise ValueError(f'Population not conserved in {name}')
    if any(not 0 <= f <= 1 for f in shared['aerial_share_sensitivity']):
        raise ValueError('Aerial share out of range')


def envelope(air, people=0., radius_m=500., tonnes_person=50., payload_share=.5, gas='hydrogen'):
    if radius_m <= 0 or tonnes_person <= 0 or not 0 < payload_share <= 1 or people < 0:
        raise ValueError('Nonphysical envelope inputs')
    volume = 4 * math.pi / 3 * radius_m**3
    area = math.pi * radius_m**2
    lift = air[f'lift_{gas}_kg_m3']
    gross = volume * lift
    payload = gross * payload_share
    residents = payload / (tonnes_person * 1e3)
    count = people / residents
    gas_mass = volume * (air['density_kg_m3'] - lift)
    contaminant_mass = volume*air['density_kg_m3']*(1-GAS_PURITY)
    component_mass = gas_mass-contaminant_mass
    reserve = LIFT_SHARES['reserve']
    hardware = LIFT_SHARES['envelope_flight_structure']+LIFT_SHARES['control_flight_services']
    nominal = 1-reserve
    return dict(height_km=air['height_km'], gas=gas, radius_m=radius_m, volume_m3=volume,
                projected_area_km2=area/1e6, gross_supported_t=gross/1e3,
                district_payload_t=payload/1e3, other_and_reserve_t=(gross-payload)/1e3,
                gas_mass_t=gas_mass/1e3, neutral_trim_gas_t=gas_mass*nominal/1e3 if payload_share==.5 else None,
                neutral_trim_ballonet_air_t=air['density_kg_m3']*volume*reserve/1e3 if payload_share==.5 else None,
                residents=residents, fleet_equivalent_count=count,
                fleet_count_rounded_up=math.ceil(count), people=people, district_t_per_person=tonnes_person,
                fleet_district_gt=people*tonnes_person/1e9,
                fleet_other_and_reserve_gt=count*(gross-payload)/1e12,
                fleet_envelope_control_hardware_gt=count*gross*hardware/1e12 if payload_share==.5 else None,
                fleet_unallocated_reserve_gt=count*gross*reserve/1e12 if payload_share==.5 else None,
                fleet_gas_gt=count*gas_mass/1e12,
                fleet_lifting_species_gt=count*component_mass/1e12,
                fleet_gas_entrained_air_gt=count*contaminant_mass/1e12,
                fleet_nominal_lifting_species_gt=count*component_mass*nominal/1e12 if payload_share==.5 else None,
                fleet_neutral_trim_gas_gt=count*gas_mass*nominal/1e12 if payload_share==.5 else None,
                sum_projected_area_km2=count*area/1e6,
                sum_projected_share_moon=count*area/(4*math.pi*MOON_RADIUS**2),
                spacing_3_diameters_square_grid_km2=count*(6*radius_m)**2/1e6,
                spacing_3_diameters_square_grid_share_moon=count*(6*radius_m)**2/(4*math.pi*MOON_RADIUS**2))


def drag(air, radius_m, relative_m_s, cd=.47, efficiency=.7):
    if relative_m_s < 0 or cd < 0 or not 0 < efficiency <= 1:
        raise ValueError('Nonphysical drag inputs')
    force = .5*air['density_kg_m3']*cd*math.pi*radius_m**2*relative_m_s**2
    return dict(relative_air_speed_m_s=relative_m_s, cd=cd, propulsive_efficiency=efficiency,
                drag_n=force, thrust_power_w=force*relative_m_s,
                input_power_w=force*relative_m_s/efficiency)


def water(people, litres_day, recovery, height_km, pump_efficiency=.7):
    if people < 0 or litres_day < 0 or not 0 <= recovery <= 1 or height_km < 0 or not 0 < pump_efficiency <= 1:
        raise ValueError('Nonphysical water inputs')
    # 1 kg/L is an explicitly assumed fresh-water conversion, not a seawater property.
    kg_s = people*litres_day*(1-recovery)/JULIAN_DAY
    potential = MOON_GM*(1/MOON_RADIUS - 1/(MOON_RADIUS+height_km*1e3))
    return dict(people=people, domestic_litres_person_day=litres_day, recovery=recovery, height_km=height_km,
                gross_km3_year=people*litres_day*JULIAN_YEAR_DAYS/1e12,
                makeup_km3_year=people*litres_day*(1-recovery)*JULIAN_YEAR_DAYS/1e12,
                makeup_kg_s=kg_s, ideal_lift_j_kg=potential,
                ideal_lift_kwh_m3=potential/3600,
                pump_efficiency=pump_efficiency, makeup_lift_tw=kg_s*potential/pump_efficiency/1e12)


def results():
    data = {p: load_product(ROOT/p, s) for p, s in INPUTS.items()}
    shared, ships, pop = [data[p] for p in INPUTS]
    validate_cases(shared)
    air = {a['height_km']: a for a in ships['air']}
    floor = pop['affluence']['living']['metropolis_floor_m2']
    area = 4*math.pi*MOON_RADIUS**2/1e6
    heat_case = pop['heat']['moon']['european']['industry_in_orbit_thermal']
    districts = {name: dict(**d, total_t_per_person=sum(d.values()),
                           building_kg_per_m2=d['building_t']*1e3/floor) for name,d in DISTRICTS.items()}
    # Reference population comes from the shared case, never an independent allocation.
    reference_people = shared['cases']['total25']['lunar_aerial']*1e9
    lift = [envelope(air[h], reference_people, tonnes_person=d['total_t_per_person'], gas=g)
            for h in (0.,10.,20.,40.) for d in districts.values() for g in ('hydrogen','helium')]
    central = envelope(air[10.], reference_people)
    stationary = []
    for speed in (0.,5.,10.,20.,30.):
        row = drag(air[10.],500.,speed)
        row.update(fleet_input_tw=row['input_power_w']*central['fleet_equivalent_count']/1e12,
                   input_kw_per_person=row['input_power_w']/central['residents']/1e3)
        stationary.append(row)
    roaming = []
    for speed in (5.,10.):
        active = drag(air[10.],500.,speed)
        for duty in ROAMING_DUTIES:
            roaming.append(dict(relative_air_speed_m_s=speed,powered_duty_fraction=duty,
                mean_input_kw_per_person=active['input_power_w']*duty/central['residents']/1e3,
                fleet_mean_input_tw=active['input_power_w']*duty*central['fleet_equivalent_count']/1e12,
                one_event_energy_mwh={str(hours):active['input_power_w']*hours/1e6 for hours in (1,6,24)}))
    cases = {}
    for name, c in shared['cases'].items():
        lunar = (c['lunar_surface']+c['lunar_aerial'])*1e9
        ground, aloft = c['lunar_surface']*1e9, c['lunar_aerial']*1e9
        held = envelope(air[10.],aloft)
        heat = lunar*heat_case['kw_per_person']/1e9
        cases[name] = dict(population_billion=c, lunar_population_billion=lunar/1e9,
            lunar_total_floor_km2=lunar*floor/1e6, aerial_floor_km2=aloft*floor/1e6,
            ground_residential_footprint_km2_at_40000_per_km2=ground/40000.,
            all_lunar_ground_counterfactual_km2_at_40000_per_km2=lunar/40000.,
            food_cropland_km2=[lunar*x/1e6 for x in pop['food']['moon_m2_per_person_planning']],
            cropland_fraction_dry=[lunar*x/1e6/pop['food']['moon_dry_land_km2'] for x in pop['food']['moon_m2_per_person_planning']],
            wild_aerial_food_credited_people=0.,
            runoff_m3_person_year=pop['water']['moon_runoff_km3']*1e9/lunar,
            local_human_heat_tw=heat, human_heat_w_m2=heat*1e6/area,
            extra_dimming_percent=[lunar/1e9*x for x in heat_case['dimming_percent_per_billion']],
            aerial_reference=held,
            aerial_domestic_water=[water(aloft,l,r,h) for l in (50.,150.) for r in (.9,.99) for h in (10.,20.)],
            ground_analogue_economy_stock_gt={key:lunar/1e9*value for key,value in pop['materials']['stocks_t_per_person'].items()},
            aerial_share_sensitivity=[dict(aerial_share=f,aerial_floor_km2=lunar*f*floor/1e6,
                    reference=envelope(air[10.],lunar*f)) for f in shared['aerial_share_sensitivity']])
    return dict(schema=SCHEMA, producer=dict(study='inhabited_volume',
                files={str(Path(__file__).relative_to(ROOT)):digest(__file__)},
                constants=constants_used([__file__]), inputs={p:digest(ROOT/p) for p in INPUTS}),
        evidence='Closed-form geometry, Archimedean lift, drag and service-account screens on committed model products. No habitat structure, population capacity, weather simulation, ecosystem or sustainable allocation is established.',
        reading_rule='Shares are fractions; _percent is percentage points. Gt=1e12 kg. Sum projected areas ignores overlap and is not radiative-transfer shade or ecological occupancy. Grid area is an imposed single-layer spacing diagnostic. Lift-gas stock is the 98%-species / 2%-air mixture, separate from supported mass; component species are recorded separately. Historical economy stock proxies overlap district stocks and must not be added to them.',
        assumptions=dict(reference_case='total25',reference_height_km=10.,reference_radius_m=500.,
            floor_m2_per_person=floor,lift_allocation_shares=LIFT_SHARES,districts=districts,
            gas_density_rule='Ambient density minus sky_ships tabulated gross supported-mass density gives the 98%-volume lifting species plus 2% entrained-air mixture. Component masses subtract that air, using rounded tables; not pure-hydrogen inventory.',
            lifting_gas_purity_by_volume=GAS_PURITY,
            habitation_direction=shared['habitation_direction'],
            roaming_rule='Co-moving drift is the default for much floating habitation, per author 2026-10-09. Powered5/10m/s repositioning at 1/5/10% duty is a sensitivity, not a weather-derived safe operating schedule. Event energies at 1/6/24 hours are reserve requirements, not an installed battery specification.',
            reserve_rule='15% of gross lift remains unallocated contingency, not manufactured hardware. At nominal load, an ideal air ballonet occupies 15% of the envelope and lift gas 85% to give neutral buoyancy; full-gas stock is capacity inventory, with spare gas stored externally. Actual trim/pressure/compartment design unresolved.',
            pump_water_density_kg_l=1.,pump_efficiency=.7,
            drag_rule='Unstreamlined sphere Cd=0.47 and 70% propulsive efficiency are sensitivities, not validated city aerodynamics; input power excludes treatment and transport. Co-moving drift has zero relative speed in this screen; tethering transfers static load to ground.',
            human_heat_rule='Inherited 2.778kW/person thermal supply with industry in orbit; excludes city-specific stationkeeping, pumping and their generation losses. Dimming is an inherited broad climate-response screen, not a new climate result.'),
        air=[air[h] for h in (0.,10.,20.,40.)],lift_cases=lift,
        radius_sensitivity=[envelope(air[10.],reference_people,radius_m=r) for r in (250.,500.,1000.)],
        held_city_drag=stationary, roaming_control_sensitivity=roaming,
        water_cases=[water(reference_people,l,r,h) for l in (50.,150.) for r in (.9,.99) for h in (10.,20.)],
        cases=cases,
        food_sensitivity=dict(wild_aerial_food_credited_people=0,
            note='Ecology and resources separately evaluate the 1% projected-cover, 1000gC/m2/year, 5% harvest case. Its gross calorie equivalents are not dependable production or a complete diet; no capacity here is increased by it.'),
        architecture_assessment=dict(
            ground='Residential density translates to land directly; farms, conservation, traffic and regional water remain separate.',
            supported_tethered='Elevated floors need actual structures. A tether restrains buoyancy and drag; a tether from below does not lift a heavy city. Tower-supported capacity requires tower and foundation sizing.',
            drifting_buoyant='Author default for much floating habitation: roam with the air, budgeting navigation, collision and storm avoidance, rendezvous, rescue and occasional repositioning; safe routes and duty cycles unresolved.',
            held_buoyant='Archimedean mass allowances are conditional; holding power rises with relative wind cubed, while gust and tether loads rise with its square.'),
        unresolved=['Complete buoyancy and load path, gas purity/leakage/fire separation, rain/icing loads and compartment-loss resilience.',
                    'Altitude oxygen availability, chronic partial-gravity health, evacuation, accessibility and rescue.',
                    'Regional wind shear and turbulence, water/food inventories through lunar night, actual routes and service access.',
                    'Agricultural water and phosphorus closure, shared surface-aerial sunlight and climate/biological aerosol effects.',
                    'Ecological exclusion corridors, projected opacity, artificial light and acceptable occupancy; no zones or allowable sky share selected.'])


def main():
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(results(),indent=2)+'\n')
    print(OUT.relative_to(ROOT))

if __name__ == '__main__':
    main()
