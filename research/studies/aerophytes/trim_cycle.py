"""Conditional fixed-altitude mass/trim budgets, not a flight controller.

Supported mass excludes lifting gas and ambient-air ballonet contents. The gas
lift coefficient already subtracts the lifting gas's own weight. Every other
actual mass change, including reserve and metabolic water, belongs in the signed
ledger. Capacity alone is not a mass flow or a control actuator.
"""
from __future__ import annotations

import math

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS
from research.studies.aerophytes.biology import fermentation_feed


def _nonnegative(**values):
    if any(not math.isfinite(v) or v < 0 for v in values.values()):
        raise ValueError(f'Finite nonnegative inputs required: {tuple(values)}')


def _positive(**values):
    if any(not math.isfinite(v) or v <= 0 for v in values.values()):
        raise ValueError(f'Finite positive inputs required: {tuple(values)}')


def net_lift_per_hydrogen_kg(air_density_kg_m3, lifting_mix_density_kg_m3,
                             hydrogen_mole_fraction=.98):
    """Constant-p,T lifting mixture, with its remaining fraction ambient air.

    Removal/addition changes mixture volume at this composition; separation and
    purity control are not solved. This is not constant-volume gas replacement.
    """
    _positive(air_density=air_density_kg_m3, mixture_density=lifting_mix_density_kg_m3)
    if not 0 < hydrogen_mole_fraction <= 1:
        raise ValueError('Hydrogen mole fraction must lie in (0,1]')
    h2_density = lifting_mix_density_kg_m3-(1-hydrogen_mole_fraction)*air_density_kg_m3
    if not 0 < h2_density <= lifting_mix_density_kg_m3 < air_density_kg_m3:
        raise ValueError('Mixture densities and composition are incompatible')
    return (air_density_kg_m3-lifting_mix_density_kg_m3)/h2_density


def trim_ledger(net_lift_kg_per_hydrogen_kg, *,
                humidity_water_uptake_kg=0., rain_water_kg=0., other_water_intake_kg=0.,
                evaporation_kg=0., drainage_kg=0., water_dump_kg=0.,
                metabolic_water_change_kg=0., dry_mass_change_kg=0.,
                other_supported_mass_change_kg=0., hydrogen_in_kg=0.,
                hydrogen_leak_kg=0., hydrogen_vent_kg=0., initial_residual_kg=0.):
    """Signed step: upward residual = net-gas-lift change - supported-mass change.

    Inflows/outflows are nonnegative magnitudes; dry, metabolic-water and other
    mass changes are signed. Initial residual is net support minus actual load.
    Water lost by drainage/dumping must not also be entered as evaporation.
    This is one common pressure/temperature state, not a vertical trajectory.
    """
    _positive(lift_coefficient=net_lift_kg_per_hydrogen_kg)
    _nonnegative(humidity=humidity_water_uptake_kg, rain=rain_water_kg,
                 intake=other_water_intake_kg, evaporation=evaporation_kg,
                 drainage=drainage_kg, dump=water_dump_kg, hydrogen_in=hydrogen_in_kg,
                 leak=hydrogen_leak_kg, vent=hydrogen_vent_kg)
    signed = (metabolic_water_change_kg, dry_mass_change_kg,
              other_supported_mass_change_kg, initial_residual_kg)
    if any(not math.isfinite(v) for v in signed):
        raise ValueError('Signed changes must be finite')
    water = (humidity_water_uptake_kg+rain_water_kg+other_water_intake_kg
             +metabolic_water_change_kg-evaporation_kg-drainage_kg-water_dump_kg)
    supported = water+dry_mass_change_kg+other_supported_mass_change_kg
    h2 = hydrogen_in_kg-hydrogen_leak_kg-hydrogen_vent_kg
    gas_lift = net_lift_kg_per_hydrogen_kg*h2
    residual = initial_residual_kg+gas_lift-supported
    return dict(water_mass_change_kg=water, dry_mass_change_kg=dry_mass_change_kg,
                supported_mass_change_kg=supported, hydrogen_mass_change_kg=h2,
                net_gas_lift_change_kg=gas_lift, initial_residual_kg=initial_residual_kg,
                final_upward_residual_kg=residual,
                ledger_residual_kg=residual-initial_residual_kg-gas_lift+supported,
                reading_rule='Positive residual tends upward. Gas weight is already inside net gas lift; do not subtract it again.')


def vent_for_neutrality(final_upward_residual_kg, net_lift_kg_per_hydrogen_kg,
                        available_hydrogen_kg=None):
    """Additional vent OR refill after all other simultaneous flows are netted.

    This is one optional fixed-altitude strategy. Actual gas inventory limits
    venting; expansion capacity and production rate limit refill separately.
    No simultaneous vent/refill is imposed by this algebraic control step.
    """
    _positive(lift_coefficient=net_lift_kg_per_hydrogen_kg)
    if not math.isfinite(final_upward_residual_kg):
        raise ValueError('Finite residual required')
    vent = max(final_upward_residual_kg, 0)/net_lift_kg_per_hydrogen_kg
    refill = max(-final_upward_residual_kg, 0)/net_lift_kg_per_hydrogen_kg
    feasible = None
    if available_hydrogen_kg is not None:
        _nonnegative(inventory=available_hydrogen_kg)
        feasible = vent <= available_hydrogen_kg
    return dict(additional_hydrogen_vent_kg=vent, additional_hydrogen_refill_kg=refill,
                vent_net_lift_loss_kg=vent*net_lift_kg_per_hydrogen_kg,
                vent_inventory_feasible=feasible,
                residual_after_requested_control_kg=final_upward_residual_kg+
                    (refill-vent)*net_lift_kg_per_hydrogen_kg,
                rate_capacity_and_flight_control_solved=False)


def water_inventory_band(initial_free_water_kg, minimum_free_water_kg,
                         maximum_free_water_kg, *, inflow_kg_day=0.,
                         outflow_kg_day=0., interval_days=1.):
    """Actual controllable free water, not total wet tissue or spare lift capacity.

    Constant-flow inventory only; rates may be netted from fog/rain, transpiration
    and drainage. Endpoint and time-to-bound diagnostics do not silently clamp a
    violated water budget. Minimum includes water unavailable for routine trim.
    """
    _nonnegative(initial=initial_free_water_kg, minimum=minimum_free_water_kg,
                 maximum=maximum_free_water_kg, inflow=inflow_kg_day,
                 outflow=outflow_kg_day, duration=interval_days)
    if not minimum_free_water_kg <= initial_free_water_kg <= maximum_free_water_kg:
        raise ValueError('Initial water must lie inside the prescribed band')
    rate = inflow_kg_day-outflow_kg_day
    final = initial_free_water_kg+rate*interval_days
    return dict(initial_free_water_kg=initial_free_water_kg,
                actual_usable_water_for_dump_kg=initial_free_water_kg-minimum_free_water_kg,
                unfilled_water_capacity_kg=maximum_free_water_kg-initial_free_water_kg,
                net_water_change_kg=rate*interval_days, unconstrained_final_free_water_kg=final,
                required_extra_water_kg=max(minimum_free_water_kg-final,0),
                required_extra_drainage_kg=max(final-maximum_free_water_kg,0),
                days_to_minimum=(initial_free_water_kg-minimum_free_water_kg)/(-rate) if rate < 0 else None,
                days_to_maximum=(maximum_free_water_kg-initial_free_water_kg)/rate if rate > 0 else None,
                inventory_band_feasible=minimum_free_water_kg <= final <= maximum_free_water_kg)


def hydrogen_refill_budget(hydrogen_kg, interval_days, projected_area_m2=1.,
                           hydrogen_energy_j_kg=120e6,
                           fermentation_yield_mol_h2_per_mol_glucose=2.1):
    """Cost of an actual refill, with alternative photon and fermentation ledgers.

    120 MJ/kg is the declared rounded LHV basis, not process input energy. Divide
    chemical power by consistent solar conversion efficiency, or use the measured
    fermentative feed yield. Do not charge both pathways for the same gas. Annual
    equivalents assume this mass transfer repeats every interval, not a forecast.
    """
    _nonnegative(hydrogen=hydrogen_kg)
    _positive(duration=interval_days, area=projected_area_m2, energy=hydrogen_energy_j_kg)
    feed = fermentation_feed(hydrogen_kg, fermentation_yield_mol_h2_per_mol_glucose)
    annual = JULIAN_YEAR_DAYS/interval_days/projected_area_m2
    energy = hydrogen_kg*hydrogen_energy_j_kg
    return dict(hydrogen_refill_kg=hydrogen_kg, hydrogen_chemical_energy_j=energy,
                mean_hydrogen_chemical_w_m2=energy/(interval_days*JULIAN_DAY*projected_area_m2),
                repeated_cycle_hydrogen_kg_m2_year=hydrogen_kg*annual,
                fermentation_glucose_feed_kg=feed['glucose_feed_kg'],
                fermentation_feed_carbon_kg=feed['feed_carbon_kg'],
                repeated_cycle_fermentation_carbon_kg_m2_year=feed['feed_carbon_kg']*annual,
                repeated_cycle_fermentation_glucose_kg_m2_year=feed['glucose_feed_kg']*annual,
                equipment_water_heat_and_product_separation_included=False)


def retained_hydrogen_compression(payload_mass_loss_kg, initial_gas_volume_m3,
                                  air_density_kg_m3, initial_gas_pressure_pa,
                                  ambient_pressure_pa, original_design_pressure_pa,
                                  available_extra_pressure_pa=0.):
    """Retain every gas molecule while reducing gas displacement isothermally.

    Either the outer bladder shrinks or an ambient-pressure air ballonet fills
    the lost gas volume. The latter needs a pressure-bearing gas/air partition.
    Fixed gas mass means deltaV=payload_loss/rho_air, NOT loss/(rho_air-rho_gas).
    The inner gas pressure rises; a freely moving equal-pressure partition cannot
    implement this state. Work figures omit inefficiency, hardware and heat flow.
    """
    _nonnegative(loss=payload_mass_loss_kg, original_design=original_design_pressure_pa,
                 available_extra=available_extra_pressure_pa)
    _positive(volume=initial_gas_volume_m3, density=air_density_kg_m3,
              initial_pressure=initial_gas_pressure_pa, ambient_pressure=ambient_pressure_pa)
    if initial_gas_pressure_pa < ambient_pressure_pa:
        raise ValueError('Initial gas pressure must not be below ambient')
    delta_v = payload_mass_loss_kg/air_density_kg_m3
    final_v = initial_gas_volume_m3-delta_v
    common = dict(required_volume_reduction_m3=delta_v, final_gas_volume_m3=final_v,
                  retained_gas_mass_change_kg=0., geometry_feasible=final_v > 0)
    if final_v <= 0:
        return dict(**common, final_gas_pressure_pa=None, added_pressure_pa=None,
                    required_design_pressure_pa=None, compression_work_on_gas_j=None,
                    minimum_work_above_ambient_j=None, extra_pressure_allowance_sufficient=False)
    ratio = initial_gas_volume_m3/final_v
    final_p = initial_gas_pressure_pa*ratio
    extra_p = final_p-initial_gas_pressure_pa
    gross_work = initial_gas_pressure_pa*initial_gas_volume_m3*math.log(ratio)
    net_work = gross_work-ambient_pressure_pa*delta_v
    return dict(**common, final_gas_pressure_pa=final_p, added_pressure_pa=extra_p,
                required_design_pressure_pa=original_design_pressure_pa+extra_p,
                compression_work_on_gas_j=gross_work,
                minimum_work_above_ambient_j=max(net_work,0.),
                extra_pressure_allowance_sufficient=extra_p <= available_extra_pressure_pa,
                actual_added_ballonet_air_kg=payload_mass_loss_kg,
                ballonet_air_mass_applies_only_to_fixed_outer_hull=True)


def ballonet_air_trim(payload_mass_loss_kg, outer_volume_m3, initial_gas_volume_m3,
                      air_density_kg_m3, ambient_pressure_pa, initial_pressure_pa,
                      original_design_pressure_pa, available_extra_pressure_pa=0.):
    """Alternative: fixed outer hull, equal pressure in gas and air compartments.

    Add an actual air mass equal to lost payload, retaining the H2 mixture. The
    whole interior pressure rises: delta p=delta m*p_ambient/(rho_air*V_outer).
    Ballonet air is compressed; treating it as unchanged ambient-density air
    would violate its pressure equilibrium. A reversible isothermal compressor
    admits air from ambient against the increasing hull pressure.
    """
    _nonnegative(loss=payload_mass_loss_kg, original_design=original_design_pressure_pa,
                 available_extra=available_extra_pressure_pa)
    _positive(volume=outer_volume_m3, gas_volume=initial_gas_volume_m3,
              density=air_density_kg_m3, ambient_pressure=ambient_pressure_pa,
              initial_pressure=initial_pressure_pa)
    if initial_gas_volume_m3 > outer_volume_m3 or initial_pressure_pa < ambient_pressure_pa:
        raise ValueError('Gas volume and initial pressure must fit the fixed hull')
    extra_p = payload_mass_loss_kg*ambient_pressure_pa/(air_density_kg_m3*outer_volume_m3)
    final_p = initial_pressure_pa+extra_p
    gas_final = initial_gas_volume_m3*initial_pressure_pa/final_p
    primitive = lambda p: p*math.log(p/ambient_pressure_pa)-p
    work = outer_volume_m3*(primitive(final_p)-primitive(initial_pressure_pa))
    return dict(actual_added_ballonet_air_kg=payload_mass_loss_kg,
                final_gas_volume_m3=gas_final, final_ballonet_volume_m3=outer_volume_m3-gas_final,
                final_internal_pressure_pa=final_p, added_pressure_pa=extra_p,
                required_design_pressure_pa=original_design_pressure_pa+extra_p,
                minimum_isothermal_pumping_work_j=max(work,0.),
                extra_pressure_allowance_sufficient=extra_p <= available_extra_pressure_pa,
                displacement_change_m3=0., retained_hydrogen_mass_change_kg=0.,
                total_mass_change_after_air_trim_kg=0.)


def evaluate(air_density_kg_m3, lifting_mix_density_kg_m3, pressure_pa,
             hydrogen_mole_fraction=.98):
    """Small conditional one-projected-m² comparisons at a caller-supplied state."""
    lift = net_lift_per_hydrogen_kg(air_density_kg_m3,lifting_mix_density_kg_m3,hydrogen_mole_fraction)
    evaporation_cases = []
    for water in (.01,.1,1.,5.):
        ledger = trim_ledger(lift,evaporation_kg=water)
        control = vent_for_neutrality(ledger['final_upward_residual_kg'],lift)
        refill = hydrogen_refill_budget(control['additional_hydrogen_vent_kg'],1.)
        evaporation_cases.append(dict(water_cycle_kg_m2=water,ledger=ledger,control=control,
                                      next_refill_cost=refill))
    # Initial neutral payload of20kg/m², not the maximum hull capacity.
    gas_volume = 20./(air_density_kg_m3-lifting_mix_density_kg_m3)
    return dict(net_supported_kg_per_hydrogen_kg=lift,
                daily_evaporation_then_refill=evaporation_cases,
                free_water_band=water_inventory_band(5.,1.,9.,outflow_kg_day=1.,interval_days=7.),
                compression_for_one_kg_loss_from_twenty_kg_payload=
                    retained_hydrogen_compression(1.,gas_volume,air_density_kg_m3,
                        pressure_pa,pressure_pa,400.),
                reading_rule='Conditional fixed-altitude strategies. Actual net mass changes, water availability, gas geometry, rates and flight are required.')
