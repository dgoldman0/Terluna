"""Conditional organism carbon, storage and replacement ledgers.

All area-normalised inputs refer to one square metre of projected organism area.
The functions do not predict photosynthesis, hydrogen physiology, mortality or
bud shape. Supply those traits explicitly. Baseline NPP is already net of host
respiration; use gross_carbon_budget to audit its required gross production.
"""
from __future__ import annotations

import math

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, SYNODIC_MONTH_DAYS


def _nonnegative(**values):
    if any(not math.isfinite(v) or v < 0 for v in values.values()):
        raise ValueError(f'finite nonnegative inputs required: {tuple(values)}')


def _fraction(name, value, *, positive=False):
    if not math.isfinite(value) or not (0 < value <= 1 if positive else 0 <= value <= 1):
        raise ValueError(f'{name} must be a fraction in {"(0,1]" if positive else "[0,1]"}')


def gross_carbon_budget(gpp_c_kg_m2_year, living_dry_kg_m2,
                       maintenance_glucose_kg_kg_day,
                       construction_glucose_kg_kg_dry=1.39,
                       carbon_fraction_dry=0.45, day_fraction=0.5,
                       night_maintenance_fraction=0.25,
                       glucose_carbon_fraction=0.4):
    """GPP -> maintenance + construction respiration + new biomass carbon.

    Living structural tissue excludes inert envelope and stored carbohydrate.
    Rates are glucose mass per living dry mass per day, not carbon fractions.
    Construction cost includes both retained tissue carbon and growth respiration.
    All maintenance-free assimilate is allocated to tissue in this upper screen.
    Community heterotrophs and hydrogen production are not included here.
    """
    _nonnegative(gpp=gpp_c_kg_m2_year, tissue=living_dry_kg_m2,
                 maintenance=maintenance_glucose_kg_kg_day,
                 construction=construction_glucose_kg_kg_dry)
    for name, value in [('carbon', carbon_fraction_dry), ('glucose carbon', glucose_carbon_fraction)]:
        _fraction(name, value, positive=True)
    _fraction('day', day_fraction)
    _fraction('night maintenance', night_maintenance_fraction)
    construction_c = construction_glucose_kg_kg_dry * glucose_carbon_fraction
    if construction_c < carbon_fraction_dry:
        raise ValueError('construction substrate supplies less carbon than the new tissue contains')
    annual_glucose_maintenance = (living_dry_kg_m2 * maintenance_glucose_kg_kg_day *
        JULIAN_YEAR_DAYS * (day_fraction + (1 - day_fraction) * night_maintenance_fraction))
    maintenance_c = annual_glucose_maintenance * glucose_carbon_fraction
    available_c = gpp_c_kg_m2_year - maintenance_c
    new_dry = max(available_c, 0) / construction_c
    growth_respiration = new_dry * (construction_c - carbon_fraction_dry)
    npp = gpp_c_kg_m2_year - maintenance_c - growth_respiration
    return dict(gpp_c_kg_m2_year=gpp_c_kg_m2_year,
        maintenance_c_kg_m2_year=maintenance_c,
        growth_respiration_c_kg_m2_year=growth_respiration,
        new_dry_kg_m2_year=new_dry, npp_c_kg_m2_year=npp,
        structural_maintenance_feasible=available_c >= 0,
        maintenance_deficit_c_kg_m2_year=max(-available_c, 0),
        carbon_residual_kg_m2_year=gpp_c_kg_m2_year-maintenance_c-growth_respiration-npp,
        accounting='NPP includes host maintenance and construction respiration; do not subtract them again.')


def dark_storage(maintenance_c_kg_m2_day, dark_days,
                 community_c_kg_m2_day=0, reserve_decay_per_day=0,
                 reserve_carbon_fraction=4/9, reserve_recovery_fraction=1,
                 final_usable_reserve_c_kg_m2=0):
    """Minimum storage with a fixed tissue-maintenance floor.

    Usable reserve obeys dS/dt=-F-k*S; structural biomass is never cannibalised.
    F includes constant host and dependent-community demand. k acts only on
    stored reserve. Reserve does not become new maintenance-bearing tissue.
    Recovery is an explicitly imposed usability fraction; 4/9 is starch C mass.
    No gas-production substrate is included unless added to the fixed demand.
    """
    _nonnegative(host=maintenance_c_kg_m2_day, dark=dark_days,
                 community=community_c_kg_m2_day, decay=reserve_decay_per_day,
                 final=final_usable_reserve_c_kg_m2)
    _fraction('reserve carbon', reserve_carbon_fraction, positive=True)
    _fraction('recovery', reserve_recovery_fraction, positive=True)
    demand = maintenance_c_kg_m2_day + community_c_kg_m2_day
    k = reserve_decay_per_day
    if k == 0:
        usable_initial = final_usable_reserve_c_kg_m2 + demand * dark_days
    else:
        usable_initial = (final_usable_reserve_c_kg_m2 * math.exp(k * dark_days)
                          + demand * math.expm1(k * dark_days) / k)
    stored_c = usable_initial / reserve_recovery_fraction
    return dict(minimum_stored_c_kg_m2=stored_c,
                minimum_reserve_dry_kg_m2=stored_c / reserve_carbon_fraction,
                initial_usable_c_kg_m2=usable_initial,
                maintenance_floor_c_kg_m2_day=demand,
                floor_consumption_c_kg_m2=demand * dark_days,
                reserve_decay_loss_c_kg_m2=usable_initial-final_usable_reserve_c_kg_m2-demand*dark_days,
                final_usable_c_kg_m2=final_usable_reserve_c_kg_m2,
                reading_rule='Reserve is additional dry payload, separate from the living structural mass.')


def net_surplus_budget(npp_c_kg_m2_year, hydrogen_cost_w_m2=0,
                       turnover_dry_kg_m2_year=0, community_food_c_kg_m2_year=0,
                       harvest_c_kg_m2_year=0, reproduction_c_kg_m2_year=0,
                       carbon_fraction_dry=0.45, biomass_energy_j_kg_dry=18e6):
    """Allocate a baseline NPP budget without charging baseline respiration twice.

    hydrogen_cost_w_m2 is an already converted assimilate/biomass opportunity
    cost, not raw H2 chemical output. Conversion, separation and pumping must
    be charged by the caller; direct photolysis uses the separate area function.
    The energetic equivalent is a requirements screen, not a biochemical yield.
    Community food means all allocated food carbon, not only animal biomass.
    """
    if not math.isfinite(npp_c_kg_m2_year):
        raise ValueError('finite NPP required')
    _nonnegative(hydrogen=hydrogen_cost_w_m2, turnover=turnover_dry_kg_m2_year,
                 community=community_food_c_kg_m2_year, harvest=harvest_c_kg_m2_year,
                 reproduction=reproduction_c_kg_m2_year)
    _fraction('carbon', carbon_fraction_dry, positive=True)
    if not math.isfinite(biomass_energy_j_kg_dry) or biomass_energy_j_kg_dry <= 0:
        raise ValueError('positive biomass energy required')
    h2_c = hydrogen_cost_w_m2 * JULIAN_DAY * JULIAN_YEAR_DAYS * carbon_fraction_dry / biomass_energy_j_kg_dry
    turnover_c = turnover_dry_kg_m2_year * carbon_fraction_dry
    allocated = h2_c + turnover_c + community_food_c_kg_m2_year + harvest_c_kg_m2_year + reproduction_c_kg_m2_year
    residual = npp_c_kg_m2_year - allocated
    return dict(hydrogen_opportunity_c_kg_m2_year=h2_c,
        turnover_c_kg_m2_year=turnover_c, total_allocated_c_kg_m2_year=allocated,
        remaining_c_kg_m2_year=residual, allocation_feasible=residual >= 0,
        reading_rule='All host baseline respiration is already inside NPP; night storage is an inventory, not a second annual sink.')


def photosynthetic_hydrogen_allocation(hydrogen_chemical_w_m2, incident_w_m2,
                                      light_to_hydrogen_efficiency,
                                      baseline_gpp_c_kg_m2_year):
    """Exclusive projected-area/time allocation for H2 and carbon fixation.

    Incident power and efficiency must have the same spectral and duty basis.
    For example a measured efficiency under pulsed PAR cannot be multiplied
    directly by all-solar annual-mean irradiance. The remaining GPP still pays
    the maintenance of the complete organism, including the H2-producing cells.
    """
    _nonnegative(hydrogen=hydrogen_chemical_w_m2, incident=incident_w_m2,
                 gpp=baseline_gpp_c_kg_m2_year)
    _fraction('light-to-H2', light_to_hydrogen_efficiency, positive=True)
    full_output = incident_w_m2 * light_to_hydrogen_efficiency
    required_fraction = hydrogen_chemical_w_m2 / full_output if full_output else (0 if hydrogen_chemical_w_m2 == 0 else None)
    assigned_fraction = min(required_fraction, 1) if required_fraction is not None else 1
    return dict(required_hydrogen_area_time_fraction=required_fraction,
        area_time_feasible=required_fraction is not None and required_fraction <= 1,
        remaining_gpp_c_kg_m2_year=baseline_gpp_c_kg_m2_year*(1-assigned_fraction),
        unmet_hydrogen_w_m2=max(hydrogen_chemical_w_m2-full_output, 0),
        reading_rule='Exclusive H2 allocation gives no simultaneous carbon credit for the same photons.')


def fermentation_feed(hydrogen_kg, mol_h2_per_mol_glucose=2.1,
                      hydrogen_molar_mass_kg_mol=0.002016,
                      glucose_molar_mass_kg_mol=0.180156,
                      glucose_carbon_fraction=0.4):
    """Gross sugar committed to a fermentative H2 pathway, not carbon destroyed.

    2.1 mol/mol is one measured reactor yield; 4 is the acetate-route ceiling.
    Residual acetate/butyrate and CO2 have explicit downstream fates. No recovery
    credit, biological pump, gas purification or nutrient closure is presumed.
    Chemical molar masses are declared stoichiometric inputs.
    """
    _nonnegative(hydrogen=hydrogen_kg)
    for value in (mol_h2_per_mol_glucose, hydrogen_molar_mass_kg_mol, glucose_molar_mass_kg_mol):
        if not math.isfinite(value) or value <= 0:
            raise ValueError('positive stoichiometric inputs required')
    _fraction('glucose carbon', glucose_carbon_fraction, positive=True)
    glucose = hydrogen_kg / hydrogen_molar_mass_kg_mol / mol_h2_per_mol_glucose * glucose_molar_mass_kg_mol
    return dict(glucose_feed_kg=glucose, feed_carbon_kg=glucose*glucose_carbon_fraction,
        kg_glucose_per_kg_hydrogen=glucose_molar_mass_kg_mol/(hydrogen_molar_mass_kg_mol*mol_h2_per_mol_glucose),
        reading_rule='Gross diverted sugar carbon, not all respired or exported; no double credit for the fermentation products.')


def bud_attachment_budget(parent_margin_kg, bud_mass_kg, bud_net_gas_lift_kg):
    """Combined static mass gate during budding; parent margin after all reserves.

    Net gas lift already subtracts gas mass. Parent must support every negative
    margin on the developmental trajectory. Independent trim/structure, nights,
    storms and reproductive maturation remain separate conditions.
    """
    _nonnegative(parent=parent_margin_kg, bud=bud_mass_kg, lift=bud_net_gas_lift_kg)
    bud_margin = bud_net_gas_lift_kg - bud_mass_kg
    return dict(bud_lift_margin_kg=bud_margin,
        parent_bud_combined_margin_kg=parent_margin_kg+bud_margin,
        attached_mass_feasible=parent_margin_kg+bud_margin >= 0,
        independent_mass_feasible=bud_margin >= 0)


def reproduction_budget(parent_area_m2, available_c_kg_m2_year, bud_dry_kg,
                        bud_hydrogen_carbon_cost_kg=0, bud_reserve_c_kg=0,
                        reproductive_lifetime_years=1, recruitment_probability=1,
                        carbon_fraction_dry=0.45,
                        construction_glucose_kg_kg_dry=None,
                        glucose_carbon_fraction=0.4):
    """Continuous-allocation lower bound on bud time and expected replacements.

    With construction cost None, supply net-primary-production carbon already
    net of construction respiration. Otherwise supply assimilate after only
    maintenance: construction respiration is charged here. Inputs must exclude
    the adult's upkeep, turnover, H2 leakage replacement and other allocations.
    Recruitment includes reaching reproductive adulthood, not merely detaching.
    No mortality distribution, increasing bud photosynthesis or age structure is
    inferred; expected replacements >=1 is necessary, not ecological stability.
    """
    _nonnegative(area=parent_area_m2, dry=bud_dry_kg, hydrogen=bud_hydrogen_carbon_cost_kg,
                 reserve=bud_reserve_c_kg, lifetime=reproductive_lifetime_years)
    if not math.isfinite(available_c_kg_m2_year):
        raise ValueError('finite available carbon required')
    _fraction('recruitment', recruitment_probability)
    _fraction('carbon', carbon_fraction_dry, positive=True)
    _fraction('glucose carbon', glucose_carbon_fraction, positive=True)
    if construction_glucose_kg_kg_dry is None:
        construction_c = bud_dry_kg * carbon_fraction_dry
        stage = 'net primary production; host construction respiration already paid'
    else:
        _nonnegative(construction=construction_glucose_kg_kg_dry)
        if construction_glucose_kg_kg_dry*glucose_carbon_fraction < carbon_fraction_dry:
            raise ValueError('construction cost below tissue carbon')
        construction_c = bud_dry_kg * construction_glucose_kg_kg_dry * glucose_carbon_fraction
        stage = 'post-maintenance assimilate; construction respiration paid here'
    cost_c = construction_c + bud_hydrogen_carbon_cost_kg + bud_reserve_c_kg
    if cost_c <= 0:
        raise ValueError('positive bud cost required')
    annual_c = max(available_c_kg_m2_year, 0) * parent_area_m2
    buds_per_year = annual_c / cost_c
    expected = buds_per_year * reproductive_lifetime_years * recruitment_probability
    return dict(carbon_budget_stage=stage, bud_carbon_cost_kg=cost_c,
        minimum_years_per_bud=cost_c/annual_c if annual_c else None,
        maximum_buds_per_year=buds_per_year,
        expected_adult_recruits_per_lifetime=expected,
        lifetime_replacement_possible=expected >= 1,
        minimum_reproductive_years=1/(buds_per_year*recruitment_probability) if buds_per_year*recruitment_probability else None)


def evaluate():
    """Small reproducible sensitivities; root runner owns output and provenance."""
    dark_days = SYNODIC_MONTH_DAYS / 2
    return dict(
        evidence='Analytic necessary conditions. No aerophyte physiology, H2 organ, offspring or sustainable community is validated.',
        maintenance=[dict(living_dry_kg_m2=mass, maintenance_glucose_kg_kg_day=rate,
                          **gross_carbon_budget(gpp, mass, rate))
                     for mass in (0.1, 1., 10.) for rate in (0.001, 0.003, 0.01, 0.031)
                     for gpp in (1., 3., 5.)],
        dark_storage=[dict(living_dry_kg_m2=mass, maintenance_glucose_kg_kg_day=rate,
                           night_maintenance_fraction=night,
                           **dark_storage(mass*rate*0.4*night, dark_days))
                      for mass in (0.1, 1., 10.) for rate in (0.001, 0.003, 0.01, 0.031)
                      for night in (0.1, 0.25, 1.)],
        net_npp_hydrogen=[dict(hydrogen_cost_w_m2=power, **net_surplus_budget(1, power))
                         for power in (0., 0.1, 1., 2.)],
        fermentation=[dict(mol_h2_per_mol_glucose=y, **fermentation_feed(1, y))
                      for y in (1.63, 2.1, 2.32, 4.)],
        photo_allocation=[dict(assumed_incident_w_m2=100., assumed_efficiency=eta,
                               hydrogen_chemical_w_m2=p,
                               **photosynthetic_hydrogen_allocation(p, 100, eta, 3))
                          for eta in (0.005, 0.01, 0.02) for p in (0.1, 1., 3.)])
