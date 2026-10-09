"""Coupled requirements for hypothetical photosynthetic floaters.

Run python -m research.studies.floater_viability.run. No design is certified by
passing these algebraic gates. Independent trait sweeps are requirements, not
evidence that one material or organism jointly supplies those traits.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import hashlib
import json
import math
from pathlib import Path

from shared.constants import GAS_CONSTANT, JULIAN_DAY, JULIAN_YEAR_DAYS, SYNODIC_MONTH_DAYS
from shared.provenance import constants_used
from research.studies.floater_viability import (biology, envelope, environment, gas_biology, growth,
    mechanics, navigation, photoperiod, resources, sailing, storms, structure, trim_cycle, water)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.floater-viability/2'
YEAR_S = JULIAN_DAY * JULIAN_YEAR_DAYS
H2_LHV_J_KG = 120e6  # rounded chemical-energy basis of the explicitly assumed solar efficiency
CONSTRUCTION_C_PER_KG = 1.39 * .4  # glucose construction analogue; material-specific chemistry unresolved
OUT = HERE / 'results/floater_viability.json'
# Design guesses for the coupled giants, each given once with its range in the notes.
COMPARTMENTS_ROUND = 8          # least number of gas cells in a round body (4-16)
SPARE_TRIM_PA = 300.            # superpressure a round body's hull carries to bank rain (100-1000 Pa)
SPARE_TRIM_PA_COLONY = 100.     # the same for a colony module, which drinks as it loses (100-300 Pa)
SKIN_WATER_KG_M2 = .5           # least water held per m² of outer skin, the fire barrier (0.1-2 kg/m²)
RADIANT_FRACTION = .1           # share of a hydrogen fire's heat radiated (0.03-0.16 in the literature)
FIRE_BURN_S = 30.               # time a torn gas cell burns out (30-120 s); the dose does not depend on it
FIRE_DISTANCE_RADII = 2.        # neighbour skin's distance from a burning cell's centre, in cell radii (1-2)
DUMPABLE_WATER_KG_M2 = 4.       # free water above the 1 kg/m² floor of the 5 kg/m² store
SIZE_DIAMETERS = (20., 30., 40., 50., 60., 80., 100., 150., 200., 300., 400., 500., 700., 1000., 1500., 2000.,
                  2500., 3000., 4000., 5000.)
MODULE_DIAMETERS = (20., 30., 40., 50., 60., 80., 100., 150., 200., 300.)
MODULE_RADIUS_M = 30.           # the colony's founding module, 60 m across, the smallest that passes at 10 km


@dataclass(frozen=True)
class Traits:
    """Requirements tuple: values are not a demonstrated joint phenotype."""
    name: str = 'fixed_dark_half_cycle_stress'
    architecture: str = 'sphere'
    allowable_mpa: float = 5.
    material_density_kg_m3: float = 1500.
    structural_minimum_um: float = 25.
    barrier_um: float = 100.
    permeability_barrer: float = .173
    partition_structural_kg_m2: float = .05
    structural_dry_fraction: float = .9
    maximum_lift_utilization: float = .7
    living_dry_kg_m2: float = 1.
    community_dry_host_fraction: float = .1
    living_water_fraction: float = .9
    free_water_kg_m2: float = 5.
    reserve_water_fraction: float = .5
    gpp_c_kg_m2_year: float = 2.
    maintenance_glucose_kg_kg_day: float = .003
    night_maintenance_fraction: float = .25
    consumer_food_c_kg_m2_year: float = .1
    inert_turnover_per_year: float = .1
    living_turnover_per_year: float = .5
    incident_cycle_mean_w_m2: float = 200.
    assumed_full_solar_to_h2_lhv_efficiency: float = .01
    fermentation_h2_mol_per_glucose_mol: float = 2.1
    bottom_overpressure_pa: float = 20.
    relative_gust_m_s: float = 10.
    relative_humidity: float = .6
    leaf_above_air_k: float = 5.
    co2_ppm: float = 400.
    leaf_fraction_host_respiration: float = .5
    annual_dark_fraction: float = .5
    longest_dark_days: float = SYNODIC_MONTH_DAYS/2
    dark_hydrogen_strategy: str = 'fermentation'
    lobe_radius_m: float = 2.
    tendon_allowable_mpa: float = 100.
    tendon_turnover_per_year: float | None = None
    partition_count: float = 1.         # great-circle walls inside the hull, for their mass
    gas_air_partitions: float = 1.      # of those, the walls with hydrogen on one side and air on the other
    fitting_share: float = 0.
    seam_mass_share: float = 0.
    cell_area_factor: float = 1.
    network_kg_m2: float = 0.
    skin_water_kg_m2: float = 0.        # water held per m² of outer skin (fire barrier), carried as mass


def night_ballast_requirement(hydrogen_loss_kg, air_density, lifting_mix_density,
                             hydrogen_mole_fraction=.98):
    """Constant-mixture first-order lift loss; purity maintenance is unsolved."""
    if hydrogen_loss_kg < 0 or not 0 < hydrogen_mole_fraction <= 1:
        raise ValueError('Nonnegative loss and physical composition required')
    h2_density = lifting_mix_density - (1-hydrogen_mole_fraction)*air_density
    if not 0 < h2_density < lifting_mix_density < air_density:
        raise ValueError('Inconsistent densities')
    return hydrogen_loss_kg*(air_density-lifting_mix_density)/h2_density


def contracting_pure_hydrogen_bladder(radius_m, flux_mol_m2_s, pressure_pa,
                                     temperature_k, time_s):
    """Independent constant-p,T spherical contraction; pure gas, no air ingress.

    dr/dt=-JRT/p. Fixed initial wetted area overestimates loss for this shape.
    It is not a dynamics model for the partitioned,98%-H2 reference organism.
    """
    if min(radius_m, pressure_pa, temperature_k) <= 0 or min(flux_mol_m2_s,time_s) < 0:
        raise ValueError('Invalid gas/geometry inputs')
    radius = max(0.,radius_m-flux_mol_m2_s*GAS_CONSTANT*temperature_k/pressure_pa*time_s)
    return dict(final_radius_m=radius, inventory_fraction_remaining=(radius/radius_m)**3)


def parent_funded_reproduction(available_assimilate_c_kg_m2_year,
                               full_surface_gpp_c_kg_m2_year,
                               bud_construction_and_reserve_c_kg_m2,
                               bud_h2_energy_j_m2, full_surface_h2_w_m2,
                               existing_h2_area_fraction):
    """Parent's fixed collecting area funds one same-sized daughter.

    Annual available assimilate already pays current gas, maintenance, turnover
    and consumers. Diverting a fraction y of parent area to daughter gas removes
    GPP0*y carbon each year. T=(Cbud+GPP0*Egas/Pgas/year)/annual_available.
    Daughter photosynthesis, developmental upkeep, attachment mechanics, gas
    processing equipment and mortality are omitted. This is a conditional
    allocation time, not a universal lower bound or a developmental simulation.
    """
    if min(full_surface_gpp_c_kg_m2_year,bud_construction_and_reserve_c_kg_m2,
           bud_h2_energy_j_m2,existing_h2_area_fraction) < 0 or full_surface_h2_w_m2 <= 0:
        raise ValueError('Invalid construction inputs')
    fill_years = bud_h2_energy_j_m2/full_surface_h2_w_m2/YEAR_S
    if available_assimilate_c_kg_m2_year <= 0 or existing_h2_area_fraction >= 1:
        return dict(parent_funded_years=None, full_surface_gas_fill_years=fill_years,
                    daughter_gas_area_fraction=None, total_gas_area_fraction=None,
                    allocations_close=False)
    years = (bud_construction_and_reserve_c_kg_m2 + full_surface_gpp_c_kg_m2_year*fill_years)/available_assimilate_c_kg_m2_year
    gas_share = fill_years/years if years else 0.
    residual = years*(available_assimilate_c_kg_m2_year-full_surface_gpp_c_kg_m2_year*gas_share)-bud_construction_and_reserve_c_kg_m2
    return dict(parent_funded_years=years,full_surface_gas_fill_years=fill_years,
                daughter_gas_area_fraction=gas_share,
                total_gas_area_fraction=existing_h2_area_fraction+gas_share,
                carbon_residual_kg_m2=residual,
                allocations_close=abs(residual) < 1e-8 and existing_h2_area_fraction+gas_share <= 1+1e-12)


def coupled_case(radius_m, air, traits=Traits()):
    """Full-capacity pressure sizing, neutral actual gas inventory, joint ledgers.

    Two separate barrier layers: exterior and one great-circle partition. Full
    exterior plus partition area is a conservative gas-contact sensitivity. It
    is not claimed to be the actual area of a solved ballonet geometry.
    Daytime replacement uses exclusive solar allocation; nighttime replacement
    commits fermentable substrate with no coproduct recovery credit. Pure-H2
    separation and suitable organs are assumed interfaces and remain unbuilt.
    """
    t=traits; disc=math.pi*radius_m**2; volume=4*disc*radius_m/3
    if (t.cell_area_factor < 1 or t.partition_count < 0 or not 0 <= t.gas_air_partitions <= t.partition_count
            or min(t.fitting_share,t.seam_mass_share,t.network_kg_m2,t.skin_water_kg_m2) < 0):
        raise ValueError('Invalid colony, wall or joint traits')
    area=disc*t.cell_area_factor
    if not 0 <= t.annual_dark_fraction <= 1 or t.longest_dark_days < 0:
        raise ValueError('Invalid light regime')
    if (t.annual_dark_fraction == 0) != (t.longest_dark_days == 0):
        raise ValueError('Zero-dark duty and duration must agree')
    if t.dark_hydrogen_strategy not in ('fermentation','water_ballast'):
        raise ValueError('Unknown dark hydrogen strategy')
    rho=air['density_kg_m3']; delta=air['lift_hydrogen_kg_m3']; gas_density=rho-delta
    p=air['pressure_atm']*101325.; temp=air['temperature_c']+273.15
    gravity=mechanics.gravity_at_height(air['height_km']*1000)
    loads=mechanics.pressure_loads(radius_m,rho,gas_density,gravity,t.bottom_overpressure_pa,t.relative_gust_m_s)
    barrier_mass_area=t.barrier_um*1e-6*t.material_density_kg_m3
    if t.architecture == 'sphere':
        skin=mechanics.spherical_skin(radius_m,loads['design_differential_pressure_pa'],
            t.allowable_mpa*1e6,t.material_density_kg_m3,t.structural_minimum_um*1e-6,barrier_mass_area)
        exterior=skin['area_m2']
    elif t.architecture == 'lobed':
        skin=mechanics.lobed_skin_tendons(radius_m,loads['design_differential_pressure_pa'],
            min(t.lobe_radius_m,radius_m),t.allowable_mpa*1e6,t.material_density_kg_m3,
            t.tendon_allowable_mpa*1e6,t.material_density_kg_m3,t.structural_minimum_um*1e-6,barrier_mass_area)
        exterior=skin['film_area_m2']
    else:
        raise ValueError('Unknown architecture')
    partition=mechanics.partition_mass(radius_m,barrier_mass_area+t.partition_structural_kg_m2,t.partition_count)
    tendon=skin.get('tendon_mass_kg',0.)
    fittings=t.fitting_share*tendon
    seams=t.seam_mass_share*(skin['total_structure_mass_kg']-tendon)
    network=t.network_kg_m2*area
    structure=skin['total_structure_mass_kg']+partition+fittings+seams+network
    skin_water=t.skin_water_kg_m2*exterior
    # Walls between gas cells hold hydrogen on both sides and lose none; only the outer skin and the
    # gas-air partitions pass hydrogen out of the body.
    wetted=exterior+t.gas_air_partitions*disc
    flux=envelope.permeation_flux(t.permeability_barrer*envelope.BARRER_SI,t.barrier_um*1e-6,.98*p,0.)
    leak=envelope.measured_flux_budget(hydrogen_flux_mol_m2_s=flux,gas_wetted_area_m2=wetted,projected_area_m2=area)
    annual_h2_kg_m2=leak['hydrogen_loss_kg_day']*JULIAN_YEAR_DAYS/area
    night_days=t.longest_dark_days
    dark_h2_kg_m2=leak['hydrogen_loss_kg_day']*night_days/area
    # A shared H2-dark/carbon-deficit duty is imposed, not inferred from lux.
    fermented_share=t.annual_dark_fraction if t.dark_hydrogen_strategy=='fermentation' else 0.
    night_feed=biology.fermentation_feed(annual_h2_kg_m2*fermented_share,t.fermentation_h2_mol_per_glucose_mol)
    dark_feed=biology.fermentation_feed(dark_h2_kg_m2,t.fermentation_h2_mol_per_glucose_mol)
    daily_feed=biology.fermentation_feed(leak['hydrogen_loss_kg_day']/area if t.dark_hydrogen_strategy=='fermentation' else 0.,t.fermentation_h2_mol_per_glucose_mol)
    daytime_h2_chemical_w_m2=annual_h2_kg_m2*(1-fermented_share)*H2_LHV_J_KG/YEAR_S
    photo=biology.photosynthetic_hydrogen_allocation(daytime_h2_chemical_w_m2,
        t.incident_cycle_mean_w_m2,t.assumed_full_solar_to_h2_lhv_efficiency,t.gpp_c_kg_m2_year)
    gross=biology.gross_carbon_budget(photo['remaining_gpp_c_kg_m2_year'],t.living_dry_kg_m2,
        t.maintenance_glucose_kg_kg_day,day_fraction=1-t.annual_dark_fraction,
        night_maintenance_fraction=t.night_maintenance_fraction)
    host_night_c_day=t.living_dry_kg_m2*t.maintenance_glucose_kg_kg_day*t.night_maintenance_fraction*.4
    community_c_day=t.consumer_food_c_kg_m2_year/JULIAN_YEAR_DAYS
    storage=biology.dark_storage(host_night_c_day+daily_feed['feed_carbon_kg'],
        night_days,community_c_kg_m2_day=community_c_day)
    stored_wet=storage['minimum_reserve_dry_kg_m2']/(1-t.reserve_water_fraction)
    community_dry=t.living_dry_kg_m2*t.community_dry_host_fraction
    tissue_wet=(t.living_dry_kg_m2+community_dry)/(1-t.living_water_fraction)
    actual_mass=structure+skin_water+area*(tissue_wet+t.free_water_kg_m2+stored_wet)
    trim=mechanics.neutral_trim(radius_m,rho,gas_density,actual_mass)
    utilization=trim['required_lifting_mix_volume_fraction']
    structural_dry=structure/area*t.structural_dry_fraction
    tendon_dry=(tendon+fittings+network)/area*t.structural_dry_fraction  # ties renew as tendon fibre
    tendon_turnover=t.inert_turnover_per_year if t.tendon_turnover_per_year is None else t.tendon_turnover_per_year
    turnover_dry=((structural_dry-tendon_dry)*t.inert_turnover_per_year+tendon_dry*tendon_turnover
                  +t.living_dry_kg_m2*t.living_turnover_per_year)
    turnover_assimilate=turnover_dry*CONSTRUCTION_C_PER_KG
    available=(photo['remaining_gpp_c_kg_m2_year']-gross['maintenance_c_kg_m2_year']
               -night_feed['feed_carbon_kg']-t.consumer_food_c_kg_m2_year-turnover_assimilate)
    out=dict(radius_m=radius_m,diameter_m=2*radius_m,height_km=air['height_km'],traits=asdict(t),
        pressure=loads,skin=skin,partition_mass_kg=partition,gas_contact_area_m2=wetted,
        mass=dict(structure_kg=structure,structural_dry_kg_m2=structural_dry,tendon_dry_kg_m2=tendon_dry,
            fittings_kg=fittings,seams_kg=seams,network_kg=network,partitions_kg=partition,skin_water_kg=skin_water,
            living_and_community_wet_kg_m2=tissue_wet,free_water_kg_m2=t.free_water_kg_m2,
            stored_reserve_wet_kg_m2=stored_wet,actual_supported_mass_kg=actual_mass,
            supported_kg_m2=actual_mass/area,**trim),
        hydrogen=dict(**leak,annual_loss_kg_m2=annual_h2_kg_m2,dark_loss_kg_m2=dark_h2_kg_m2,
            full_solar_efficiency_lhv=t.assumed_full_solar_to_h2_lhv_efficiency,
            daytime_replacement_chemical_w_m2=daytime_h2_chemical_w_m2,**photo),
        carbon=dict(gpp_before_h2_c_kg_m2_year=t.gpp_c_kg_m2_year,
            gpp_after_h2_c_kg_m2_year=photo['remaining_gpp_c_kg_m2_year'],
            host_maintenance_c_kg_m2_year=gross['maintenance_c_kg_m2_year'],
            night_h2_gross_substrate_c_kg_m2_year=night_feed['feed_carbon_kg'],
            night_h2_glucose_kg_m2_year=night_feed['glucose_feed_kg'],
            consumer_food_c_kg_m2_year=t.consumer_food_c_kg_m2_year,
            turnover_dry_kg_m2_year=turnover_dry,turnover_construction_c_kg_m2_year=turnover_assimilate,
            remaining_assimilate_c_kg_m2_year=available,
            stage='Post-maintenance assimilate after current turnover, consumers and H2; future construction respiration unpaid.'),
        storage=storage,
        gates=dict(mass_and_reserve=utilization<=t.maximum_lift_utilization,
                   solar_h2_allocation=photo['area_time_feasible'],positive_construction_surplus=available>0,
                   passive_unfrozen_mean_air=air['temperature_c']>0,
                   assembled_material_demonstrated=False,water_route_closed=False,
                   hydrogen_organs_and_purity_closed=False,flight_and_lifecycle_validated=False))
    net_leaf=max(0.,photo['remaining_gpp_c_kg_m2_year']-
        t.leaf_fraction_host_respiration*gross['maintenance_c_kg_m2_year'])
    leaf_temp=max(5.,air['temperature_c'])+t.leaf_above_air_k
    water=environment.water_for_assimilation(net_leaf*1000,air['temperature_c'],leaf_temp,
        t.relative_humidity,p,co2_ppm=t.co2_ppm)
    night_cuticle=environment.transpiration(air['temperature_c'],max(5.,air['temperature_c']),.9,p,.0016)
    out['water']=dict(**water,night_cuticle=night_cuticle,
        solar_splitting_feed_kg_m2_year=leak['water_feed_kg_day']*JULIAN_YEAR_DAYS/area*(1-fermented_share),
        night_substrate_chemical_water_unclosed=True,collection_and_rainfall_unestablished=True,
        note='Net-leaf diffusion demand; fixed Ci/Ca. All light-derived carbon shares one water flow. Leaf heating and RH are scenarios.')
    ballast=night_ballast_requirement(dark_h2_kg_m2,rho,gas_density)
    out['solar_only_night_alternative']=dict(required_water_ballast_release_kg_m2=ballast,
        free_water_share=ballast/t.free_water_kg_m2 if t.free_water_kg_m2 else None,
        meaning='If no dark H2 synthesis and all other supported mass is constant: fixed-initial-area loss demand at maintained composition. Evaporation, intake and metabolic mass changes must be netted before any actual ballast release. Beyond the inventory this ceases to be a physical loss prediction. Water must be replaced. No dynamics or purity control.')
    if t.dark_hydrogen_strategy=='water_ballast':
        out['gates']['dark_ballast_stock']=ballast<=t.free_water_kg_m2
        out['water']['ballast_replenishment_kg_m2_year']=night_ballast_requirement(
            annual_h2_kg_m2*t.annual_dark_fraction,rho,gas_density)
        out['water']['ballast_cycle_note']='Gas lost during darkness is restored during light; discharged liquid ballast must be collected again. Fixed initial-area flux is conservative only for specified contraction geometry, and composition/control remain unsolved.'
    out['reproduction']=None
    if trim['neutral_trim_possible']:
        h2=envelope.hydrogen_budget(gas_volume_m3=volume*utilization,gas_wetted_area_m2=wetted,
            projected_area_m2=area,pressure_pa=p,temperature_k=temp,hydrogen_mole_fraction=.98,
            permeability_mol_m_m2_s_pa=t.permeability_barrer*envelope.BARRER_SI,thickness_m=t.barrier_um*1e-6)
        actual_h2=trim['actual_hydrogen_mass_kg']
        out['hydrogen'].update(inventory_kg=actual_h2,
            ideal_pvt_inventory_kg=h2['hydrogen_inventory_kg'],
            ideal_pvt_vs_inherited_density_fraction=h2['hydrogen_inventory_kg']/actual_h2-1,
            initial_turnover_days=actual_h2/leak['hydrogen_loss_kg_day'] if leak['hydrogen_loss_kg_day'] else None,
            initial_inflation_minimum_work_j=actual_h2/envelope.HYDROGEN_MOLAR_MASS_KG*envelope.WATER_SPLITTING_GIBBS_J_MOL,
            inventory_convention='Rounded inherited density mixture used consistently for supported mass, H2 inventory, ballast and reproduction; ideal pVT discrepancy explicitly recorded.')
        frac=dark_h2_kg_m2*area/actual_h2
        out['solar_only_night_alternative'].update(initial_inventory_fraction_lost=frac,
            initial_flux_approximation_small_loss=frac<=.05,
            initial_flux_demand_exceeds_inventory=frac>=1,
            inventory_capped_ballast_loss_kg_m2=min(ballast,actual_mass/area))
        bud_c=(structural_dry+t.living_dry_kg_m2+community_dry)*CONSTRUCTION_C_PER_KG+storage['minimum_stored_c_kg_m2']
        repro=parent_funded_reproduction(available,t.gpp_c_kg_m2_year,bud_c,
            actual_h2/area*H2_LHV_J_KG,
            t.incident_cycle_mean_w_m2*t.assumed_full_solar_to_h2_lhv_efficiency,
            photo['required_hydrogen_area_time_fraction'])
        repro['minimum_parent_reproductive_years_at_half_recruitment']=(repro['parent_funded_years']*2 if repro['parent_funded_years'] else None)
        repro['bud_construction_and_reserve_c_kg_m2']=bud_c
        out['reproduction']=repro
        b=trim_cycle.net_lift_per_hydrogen_kg(rho,gas_density)
        daily_leak=leak['hydrogen_loss_kg_day']/area
        simultaneous=trim_cycle.trim_ledger(b,other_water_intake_kg=1.,evaporation_kg=1.,
            hydrogen_leak_kg=daily_leak,hydrogen_in_kg=daily_leak)
        out['trim_comparisons']=dict(
            balanced_daily_water_and_gas_flow=simultaneous,
            one_kg_water_loss_and_unreplaced_daily_gas_loss=trim_cycle.trim_ledger(
                b,evaporation_kg=1.,hydrogen_leak_kg=daily_leak),
            retained_gas_one_kg_water_loss=trim_cycle.retained_hydrogen_compression(
                1.,volume*utilization/area,rho,p,p,loads['design_differential_pressure_pa']),
            fixed_hull_air_ballonet_one_kg_water_loss=trim_cycle.ballonet_air_trim(
                1.,volume/area,volume*utilization/area,rho,p,p,loads['design_differential_pressure_pa']),
            boundary='One projected square metre. Conditional mass-flow and alternative pressure diagnostics, not controller mass or energy already included in the design. Net intake timing, water buffering, metabolic product fates and flight remain unresolved.')
    return out


def gas_alternatives(air):
    """Same-pressure/temperature ideal gases; permeability/metabolism not compared."""
    p=air['pressure_atm']*101325.; temp=air['temperature_c']+273.15; rho=air['density_kg_m3']
    species=[dict(gas=n,molar_mass_kg_mol=m,density_kg_m3=p*m/GAS_CONSTANT/temp,
                  net_lift_kg_m3=rho-p*m/GAS_CONSTANT/temp)
             for n,m in [('hydrogen',.002016),('helium',.0040026),('methane',.016043)]]
    return dict(height_km=air['height_km'],pure_gases=species,
        warm_air=[dict(above_ambient_k=d,net_lift_kg_m3=rho*d/(temp+d)) for d in (5,10,20)],
        boundary='No gas production/purification, heat balance or common permeability assumed. Helium requires an external stock; methane has a carbon cycle.')


def design_families(height_km, structural):
    """Assembled round giants and colony modules at one height, from the structure component.

    Every hull carries its storm gust and the spare superpressure its water strategy needs:
    round bodies bank storm rain (SPARE_TRIM_PA), colony modules drink as they lose
    (SPARE_TRIM_PA_COLONY). Every outer skin holds at least SKIN_WATER_KG_M2 of water.
    """
    film = structure.assembled_allowable(20e6, .7, 10., 1.5)
    tendon = structure.assembled_allowable(100e6, .8, 100., 1.5)
    strong = structure.assembled_allowable(300e6, .8, 100., 1.5)
    gust = next(r for r in structural['design_gusts'] if r['height_km'] == height_km)
    base = replace(Traits(), allowable_mpa=film['allowable_pa']/1e6, relative_gust_m_s=gust['design_gust_m_s'],
                   skin_water_kg_m2=SKIN_WATER_KG_M2)
    round_ = replace(base, name='assembled_round', architecture='lobed',
        bottom_overpressure_pa=Traits().bottom_overpressure_pa+SPARE_TRIM_PA,
        tendon_allowable_mpa=tendon['allowable_pa']/1e6, tendon_turnover_per_year=1/100,
        partition_count=1+COMPARTMENTS_ROUND/2, gas_air_partitions=1., fitting_share=structure.end_fitting_share(.05),
        seam_mass_share=structure.seam_mass_share(.05, 2.))
    collagen = structure.COLLAGEN_SUSTAINED_PA*.8/1.5
    return dict(film=film, tendon=tendon, strong_tendon=strong, collagen_tendon_pa=collagen, gust=gust, base=base,
                round=round_, round_strong=replace(round_, name='assembled_round_strong_fibre',
                                                   tendon_allowable_mpa=strong['allowable_pa']/1e6),
                round_collagen=replace(round_, name='assembled_round_collagen_tendon',
                                       tendon_allowable_mpa=collagen/1e6))


def colony_module(radius_m, families, spare_pa=SPARE_TRIM_PA_COLONY):
    """A sphere module in a one-layer hexagonal raft, with its share of the raft's ties.

    The ties are tendon fibre, sized at the century allowable and renewed once a century.
    The film carries the round bodies' seam and ripstop allowance (2.5%), which also
    arrests its tears.
    """
    base, tendon, gust = families['base'], families['tendon'], families['gust']
    raft = structure.hex_raft(radius_m, 1.)
    ties = structure.colony_network(gust['design_pressure_pa'], 2*radius_m, tendon['allowable_pa'], base.material_density_kg_m3)
    return replace(base, name='colony_module', architecture='sphere',
                   bottom_overpressure_pa=Traits().bottom_overpressure_pa+spare_pa,
                   seam_mass_share=structure.seam_mass_share(.05, 2.),
                   cell_area_factor=raft['cell_over_disc'], network_kg_m2=ties, tendon_turnover_per_year=1/100)


def fire_skin_water(case, cells, distance_radii=FIRE_DISTANCE_RADII):
    """Water per m² of skin that the radiant heat of one burning gas cell dries at a neighbour cell.

    The torn cell's hydrogen burns in air outside the body; its heat radiates from a
    point at the cell's centre to a neighbour cell's skin distance_radii cell radii away
    (2 by default; 1 for skin that adjoins the burning cell, four times the dose). Walls
    between cells hold hydrogen on both sides and cannot burn from within.
    """
    hydrogen = case['hydrogen'].get('inventory_kg')
    if not hydrogen:
        return None
    cell_radius = case['radius_m']/cells**(1/3)
    dose = storms.neighbour_dose(hydrogen/cells, FIRE_BURN_S, RADIANT_FRACTION, distance_radii*cell_radius)['dose_j_m2']
    return dose/environment.LATENT_J_KG


def settle_case(diameter, air, traits, cells=None, fire_distance_radii=FIRE_DISTANCE_RADII):
    """A body whose gas cells and skin water answer damage and fire.

    cells=None is a round body: the fewest gas cells, at least COMPARTMENTS_ROUND and
    even, such that losing one is answered by dropping DUMPABLE_WATER_KG_M2 of free
    water. A fixed count (1 for a colony module) is kept. In both, the outer skin holds
    the larger of SKIN_WATER_KG_M2 and the water a burning neighbour cell would dry,
    carried as mass. The cells are slack gas cells inside the pressure-bearing hull: a
    lost cell's volume fills with ballonet air at the hull's pressure, so no wall carries
    the crown pressure.
    """
    least = traits.skin_water_kg_m2

    def walls(n):
        return (1+n/2) if cells is None else traits.partition_count

    def settled(n):
        water_ = least
        for _ in range(100):
            t = replace(traits, skin_water_kg_m2=water_, partition_count=walls(n))
            c = coupled_case(diameter/2, air, t)
            if not c['mass']['neutral_trim_possible']:
                return c, t, False
            need = max(least, fire_skin_water(c, n, fire_distance_radii))
            if abs(need-water_) <= 1e-9*max(1., need):
                return c, t, True
            water_ = need
        raise RuntimeError('skin water did not settle')

    def needed(case):
        n = max(COMPARTMENTS_ROUND, math.ceil(case['mass']['supported_kg_m2']/DUMPABLE_WATER_KG_M2))
        return n + n % 2

    if cells is not None:
        c, t, _ = settled(cells)
        return c, t, cells
    # Walls and skin water only add mass, so the fewest walls and least water bound the count from below.
    floor = coupled_case(diameter/2, air, replace(traits, partition_count=walls(COMPARTMENTS_ROUND)))
    n = needed(floor)
    for _ in range(5000):
        c, t, ok = settled(n)
        if not ok or needed(c) <= n:
            return c, t, max(n, needed(c))
        n += 2
    raise RuntimeError('gas cells did not settle')


def size_gates(case, traits, air, cells, spare_pa, ripstop_m=4., opening_ratio=.05,
               fire_distance_radii=FIRE_DISTANCE_RADII):
    """Hard gates for one coupled case, with damage, tear, fire and trim requirements beside them.

    Gates: mass with its 30% reserve, positive construction carbon, the solar H2
    share, and staying aloft after losing one cell by dropping free water.
    Beside them: how fast an arrested tear empties one cell, the skin water that
    stops a burning cell igniting its neighbour (carried as mass), and the water
    swing the hull's spare superpressure absorbs.
    """
    radius = case['radius_m']
    area = math.pi*radius**2*traits.cell_area_factor
    supported = case['mass']['supported_kg_m2']
    loss = structure.compartment_loss(supported, DUMPABLE_WATER_KG_M2, cells)
    gas_fraction = case['mass']['required_lifting_mix_volume_fraction']
    hydrogen = case['hydrogen'].get('inventory_kg')
    fire_water = fire_skin_water(case, cells, fire_distance_radii)
    crown = case['pressure']['crown_static_pressure_pa']
    gas_density = air['density_kg_m3']-air['lift_hydrogen_kg_m3']
    tear = structure.tear_outflow(ripstop_m, opening_ratio, crown, gas_density)
    cell_gas = 4/3*math.pi*radius**3*min(gas_fraction, 1.)/cells
    column = 4/3*radius/traits.cell_area_factor
    swing = water.allowed_swing(spare_pa, column, air['density_kg_m3'], air['pressure_atm']*101325.)
    gates = dict(mass=case['gates']['mass_and_reserve'], carbon=case['gates']['positive_construction_surplus'],
                 solar_h2=case['gates']['solar_h2_allocation'], one_compartment_loss=loss['stays_aloft'])
    return dict(diameter_m=2*radius, gas_fraction=gas_fraction, supported_kg_m2=supported,
                structural_dry_kg_m2=case['mass']['structural_dry_kg_m2'],
                tendon_dry_kg_m2=case['mass']['tendon_dry_kg_m2'],
                surplus_c_kg_m2_year=case['carbon']['remaining_assimilate_c_kg_m2_year'],
                hydrogen_kg_m2=hydrogen/area if hydrogen else None,
                hydrogen_replacement_c_kg_m2_year=(case['carbon']['gpp_before_h2_c_kg_m2_year']
                    -case['carbon']['gpp_after_h2_c_kg_m2_year']+case['carbon']['night_h2_gross_substrate_c_kg_m2_year']),
                turnover_construction_c_kg_m2_year=case['carbon']['turnover_construction_c_kg_m2_year'],
                compartments=cells, compartment_hydrogen_kg=hydrogen/cells if hydrogen else None,
                arrested_tear_m=ripstop_m, tear_outflow_m3_s=tear['volume_m3_s'],
                compartment_empties_hours=cell_gas/tear['volume_m3_s']/3600 if tear['volume_m3_s'] else None,
                fire_skin_water_kg_m2=fire_water, skin_water_carried_kg_m2=traits.skin_water_kg_m2,
                skin_water_kg_per_projected_m2=case['mass']['skin_water_kg']/area,
                spare_superpressure_pa=spare_pa, design_pressure_pa=case['pressure']['design_differential_pressure_pa'],
                gas_column_m=column, water_swing_kg_m2_at_spare=swing,
                gates=gates, all_gates=all(gates.values()))


def carbon_stop(rows, surplus_at):
    """Diameter where the surplus first reaches zero, by bisection on the model between the grid rows
    that bracket it; None if the grid never reaches it."""
    ordered = sorted(rows, key=lambda r: r['diameter_m'])
    for a, b in zip(ordered, ordered[1:]):
        if a['surplus_c_kg_m2_year'] > 0 >= b['surplus_c_kg_m2_year']:
            lo, hi = a['diameter_m'], b['diameter_m']
            while hi/lo-1 > 2e-4:
                mid = math.sqrt(lo*hi)
                lo, hi = (mid, hi) if surplus_at(mid) > 0 else (lo, mid)
            return math.sqrt(lo*hi)
    return None


def surviving_range(rows, stop=None):
    """The first contiguous run of diameters that passes every gate, and where carbon runs out."""
    ordered = sorted(rows, key=lambda r: r['diameter_m'])
    passing = [r['diameter_m'] for r in ordered if r['all_gates']]
    if not passing:
        return dict(smallest_m=None, largest_m=None, passing_m=[], carbon_stop_m=stop, largest_is_grid_edge=None)
    diameters = [r['diameter_m'] for r in ordered]
    run_ = []
    for d in diameters[diameters.index(passing[0]):]:
        if d in passing:
            run_.append(d)
        else:
            break
    return dict(smallest_m=run_[0], largest_m=run_[-1], passing_m=passing, carbon_stop_m=stop,
                largest_is_grid_edge=run_[-1] == diameters[-1])


def size_limits(air, structural):
    """Round giants and colony modules at 10 and 20 km through every gate.

    At 10 km two sensitivities carry the skin water for a neighbour skin that adjoins a
    burning cell (one cell radius, four times the dose) at the reference light.
    """
    out = {}
    for h in (10, 20):
        fam = design_families(h, structural)
        rows, stops = {}, {}
        cases = []
        for key in ('round_collagen', 'round', 'round_strong'):
            for gpp in (2., 5.):
                cases.append((key if gpp == 2. else f'{key}_gpp5', replace(fam[key], gpp_c_kg_m2_year=gpp),
                              None, SPARE_TRIM_PA, FIRE_DISTANCE_RADII))
        if h == 10:
            cases.append(('round_fire_at_one_radius', fam['round'], None, SPARE_TRIM_PA, 1.))
        for label, spare in (('colony_module', SPARE_TRIM_PA_COLONY), ('colony_module_spare_300_pa', SPARE_TRIM_PA)):
            for gpp in (2., 5.):
                if label != 'colony_module' and gpp != 2.:
                    continue
                cases.append((label if gpp == 2. else f'{label}_gpp5', (spare, gpp), 1, spare, FIRE_DISTANCE_RADII))
        if h == 10:
            cases.append(('colony_module_fire_at_one_radius', (SPARE_TRIM_PA_COLONY, 2.), 1, SPARE_TRIM_PA_COLONY, 1.))
        for name, traits, cells, spare, radii in cases:
            module = cells == 1

            def build(d, traits=traits, module=module, radii=radii):
                t0 = (replace(colony_module(d/2, fam, traits[0]), gpp_c_kg_m2_year=traits[1]) if module else traits)
                return settle_case(d, air[h], t0, cells=1 if module else None, fire_distance_radii=radii)

            rows[name] = []
            for d in (MODULE_DIAMETERS if module else SIZE_DIAMETERS):
                c, t, n = build(d)
                row = size_gates(c, t, air[h], n, spare, fire_distance_radii=radii)
                if module:
                    row['gates']['one_compartment_loss'] = True  # one module of many; the colony drops water for it
                    row['all_gates'] = all(row['gates'].values())
                    row['network_kg_m2'] = t.network_kg_m2
                rows[name].append(row)
            stops[name] = carbon_stop(rows[name], lambda d, build=build: build(d)[0]['carbon']['remaining_assimilate_c_kg_m2_year'])
        out[f'{h}_km'] = dict(design_gust_m_s=fam['gust']['design_gust_m_s'],
            film_allowable_mpa=fam['film']['allowable_pa']/1e6,
            tendon_allowable_mpa=fam['tendon']['allowable_pa']/1e6,
            strong_tendon_allowable_mpa=fam['strong_tendon']['allowable_pa']/1e6,
            collagen_tendon_allowable_mpa=fam['collagen_tendon_pa']/1e6,
            spare_superpressure_pa=dict(round=SPARE_TRIM_PA, colony_module=SPARE_TRIM_PA_COLONY),
            rows=rows, ranges={k: surviving_range(v, stops[k]) for k, v in rows.items()})
    return out


def _build(case, traits, radius_m, cost):
    area = math.pi*radius_m**2*traits.cell_area_factor
    return growth.construction_per_area(case['hydrogen']['inventory_kg']/area,
        case['mass']['structural_dry_kg_m2'], traits.living_dry_kg_m2*(1+traits.community_dry_host_fraction),
        case['storage']['minimum_stored_c_kg_m2'], cost)


def growth_tables(air, structural, nursery):
    """Ages of round giants from a 40-m juvenile and of colonies from one founding module."""
    fam = design_families(10, structural)
    routes = dict(photolysis=None, fermentation_2_1=2.1, fermentation_4_0=4.)
    targets = (100., 500., 1000.)  # radii for 200 m, 1 km and 2 km bodies
    radii = [20., 25., 30., 40., 50., 60., 75., 100., 125., 150., 200., 250., 300., 400., 500., 600., 750., 1000.,
             1250., 1500.]
    out = dict(round=[], colony=[], start_diameter_m=40., colony_start='one founding module',
               target_diameters_m=[2*r for r in targets])
    for gpp in (2., 5.):
        for route, yield_ in routes.items():
            cost = (growth.photo_carbon_per_hydrogen(gpp, Traits().incident_cycle_mean_w_m2,
                                                     Traits().assumed_full_solar_to_h2_lhv_efficiency)
                    if yield_ is None else growth.fermentation_carbon_per_hydrogen(yield_))
            for fibre in ('round', 'round_strong'):
                rows = []
                for r in radii:
                    c, t, n = settle_case(2*r, air[10], replace(fam[fibre], gpp_c_kg_m2_year=gpp))
                    if not c['mass']['neutral_trim_possible']:
                        continue
                    build = _build(c, t, r, cost)
                    rows.append(dict(radius_m=r, construction_c_kg_m2=build['total_c_kg_m2'],
                                     gas_share=build['gas_share'], compartments=n,
                                     surplus_c_kg_m2_year=c['carbon']['remaining_assimilate_c_kg_m2_year']))
                for share in (1., .5):
                    ages = growth.round_body_ages(rows, 20., targets, share)
                    out['round'].append(dict(fibre=fibre, gpp_c_kg_m2_year=gpp, hydrogen_route=route,
                                             carbon_per_kg_h2=cost, **ages, rows=rows if share == 1. else None))
            m0 = replace(colony_module(MODULE_RADIUS_M, fam), gpp_c_kg_m2_year=gpp)
            c, m, _ = settle_case(2*MODULE_RADIUS_M, air[10], m0, cells=1)
            build = _build(c, m, MODULE_RADIUS_M, cost)
            start = 2*math.sqrt(3)*MODULE_RADIUS_M**2
            for share in (1., .5):
                flat = growth.flat_body_ages(build['total_c_kg_m2'], c['carbon']['remaining_assimilate_c_kg_m2_year'],
                                             start, [math.pi*r*r for r in targets], share)
                out['colony'].append(dict(gpp_c_kg_m2_year=gpp, hydrogen_route=route, carbon_per_kg_h2=cost,
                                          growth_share=share, module_diameter_m=2*MODULE_RADIUS_M, start_area_m2=start,
                                          construction=build, surplus_c_kg_m2_year=c['carbon']['remaining_assimilate_c_kg_m2_year'],
                                          **flat))
    out['nursery'] = nursery
    return out


def canopy_nursery(air, structural, growth_eval):
    """Juveniles filled in a megaforest crown by host sugar and their own light, net of their gas loss."""
    fam = design_families(10, structural)
    rows = []
    designs = [('round', d, fam['round'], None) for d in (40., 50., 60.)]
    designs.append(('colony_module', 2*MODULE_RADIUS_M, colony_module(MODULE_RADIUS_M, fam), 1))
    for design, d, traits, cells in designs:
        c, t, n = settle_case(d, air[10], traits, cells=cells)
        if not c['mass']['neutral_trim_possible']:
            continue
        hydrogen = c['hydrogen']['inventory_kg']
        loss = c['hydrogen']['hydrogen_loss_kg_day']*JULIAN_YEAR_DAYS
        own = .5*math.pi*(d/2)**2*Traits().incident_cycle_mean_w_m2*Traits().assumed_full_solar_to_h2_lhv_efficiency*YEAR_S/H2_LHV_J_KG
        for host in growth_eval['hosts']:
            if host['tree_height_m'] < 300:
                continue
            for y in (2.1, 4.):
                fill = growth.nursery_fill_years(hydrogen, y, host['diverted_c_kg_year'], own, loss)
                rows.append(dict(design=design, juvenile_diameter_m=d, hydrogen_kg=hydrogen,
                                 passes_mass_gate=c['gates']['mass_and_reserve'],
                                 host_height_m=host['tree_height_m'],
                                 host_npp_case=host['npp_case'], diverted_share=host['diverted_share'],
                                 diverted_c_kg_year=host['diverted_c_kg_year'], mol_h2_per_mol_glucose=y,
                                 own_hydrogen_kg_year=own, gas_loss_kg_year=loss, years=fill['years'],
                                 years_on_own_light_only=hydrogen/(own-loss) if own > loss else None))
    big, _, _ = settle_case(200., air[10], fam['round'])
    canopy = growth_eval['canopy']['0']
    gas = big['hydrogen']['inventory_kg']
    full = mechanics.neutral_trim(100., air[10]['density_kg_m3'], air[10]['density_kg_m3']-air[10]['lift_hydrogen_kg_m3'], 0.)['full_hydrogen_mass_kg']
    hectare_years = []
    for label, kg in (('trimmed_200m', gas), ('full_200m', full)):
        for y in (4., 2.1, 1.63):
            carbon = kg*growth.fermentation_carbon_per_hydrogen(y)
            for npp_name, npp in growth_eval['host_npp_cases'].items():
                hectare_years.append(dict(gas=label, hydrogen_kg=kg, mol_h2_per_mol_glucose=y, sugar_c_kg=carbon,
                                          npp_case=npp_name, npp_c_kg_m2_year=npp,
                                          hectare_years=carbon/(npp*1e4)))
    return dict(juveniles=rows, two_hundred_metre_gas=hectare_years,
                canopy_gross_c_kg_m2_year=canopy['gpp_c_kg_m2_year'], canopy_npp_c_kg_m2_year=canopy['npp_c_kg_m2_year'])


def lifespans(structural_storms, module_radius_m=MODULE_RADIUS_M):
    """Lightning strikes over centuries: survival of round bodies and module losses in colonies.

    A strike that ignites a round body either burns one gas cell (the body drops
    water and lives) or spreads through a shared, dry skin (it dies): the two
    bounds. In a colony each strike burns the modules it reaches, which regrow.
    """
    rows = []
    for r in structural_storms['lightning']['strikes']:
        rate = r['strikes_per_year']
        rows.append(dict(place=r['place'], flashes_per_km2_year=r['flashes_per_km2_year'],
                         body_radius_m=r['body_radius_m'], striking_distance_m=r['striking_distance_m'],
                         strikes_per_year=rate,
                         strikes_in_300_years=r['strikes_in_300_years'],
                         survive_300_years_if_fire_spreads=growth.survival(300., rate),
                         survive_1000_years_if_fire_spreads=growth.survival(1000., rate),
                         strikes_per_century=100*rate))
    return rows


def run():
    ships=json.loads((ROOT/'research/studies/sky_ships/results/sky_ships.json').read_text())
    if ships['schema']!='terluna.research.sky-ships/1': raise ValueError('Unexpected air input')
    air={a['height_km']:a for a in ships['air']}
    reference=Traits()
    families=[reference,
        replace(reference,name='low_barrier_requirement',barrier_um=25,permeability_barrer=.039),
        replace(reference,name='higher_permeability_requirement',barrier_um=100,permeability_barrer=10.),
        replace(reference,name='thick_higher_permeability_requirement',barrier_um=1000,permeability_barrer=10.)]
    diameters=[2.,5.,10.,20.,30.,40.,60.,100.,150.,200.,300.,500.,1000.]
    scans=[]
    for f in families:
        for architecture in ('sphere','lobed'):
            for height in (10,20,40):
                for living in (.1,1.,10.):
                    t=replace(f,architecture=architecture,living_dry_kg_m2=living)
                    rows=[]
                    for d in diameters:
                        c=coupled_case(d/2,air[height],t)
                        rows.append(dict(diameter_m=d,gas_fraction=c['mass']['required_lifting_mix_volume_fraction'],
                            required_h2_solar_share=c['hydrogen']['required_hydrogen_area_time_fraction'],
                            carbon_surplus_c_kg_m2_year=c['carbon']['remaining_assimilate_c_kg_m2_year'],
                            mass_gate=c['gates']['mass_and_reserve'],
                            coupled_mass_solar_carbon_gates=(c['gates']['mass_and_reserve'] and c['gates']['solar_h2_allocation'] and c['gates']['positive_construction_surplus']),
                            parent_funded_years=c['reproduction']['parent_funded_years'] if c['reproduction'] else None))
                    scans.append(dict(traits=asdict(t),height_km=height,rows=rows,
                        boundary='Passing three algebraic conditions leaves material, water, weather, chemistry, development and lifecycle unresolved.'))
    cases=[coupled_case(d/2,air[10],f) for f in families for d in (10.,40.,100.,300.)]
    light_regimes=[]
    for f in families:
        for name,fraction,duration in [('half_cycle_stress',.5,SYNODIC_MONTH_DAYS/2),
            ('short_deficit_sensitivity',.1,3.),('wind_matched_sun_following_sensitivity',0.,0.)]:
            light_regimes.append(coupled_case(20.,air[10],replace(f,
                name=f.name+'__'+name,annual_dark_fraction=fraction,longest_dark_days=duration)))
    dark_ballast_cases=[coupled_case(d/2,air[10],replace(f,dark_hydrogen_strategy='water_ballast'))
                       for f in families for d in (40.,100.)]
    sensitivities={}
    for field,values in [('allowable_mpa',(.1,1.,5.,20.)),
        ('living_dry_kg_m2',(.1,1.,10.)),('maintenance_glucose_kg_kg_day',(.001,.003,.01,.031)),
        ('assumed_full_solar_to_h2_lhv_efficiency',(.001,.005,.01,.03)),
        ('gpp_c_kg_m2_year',(1.,2.,5.)),('free_water_kg_m2',(0.,5.,25.)),
        ('inert_turnover_per_year',(.01,.1,1.)),('relative_gust_m_s',(0.,5.,10.,20.,30.))]:
        sensitivities[field]=[coupled_case(20.,air[10],replace(reference,**{field:v})) for v in values]
    structural_scans=[]
    for strength in (1.,5.,10.):
        t=replace(reference,allowable_mpa=strength)
        rows=[]
        for d in diameters:
            c=coupled_case(d/2,air[10],t)
            rows.append(dict(diameter_m=d,gas_fraction=c['mass']['required_lifting_mix_volume_fraction'],
                mass_gate=c['gates']['mass_and_reserve']))
        structural_scans.append(dict(traits=asdict(t),height_km=10,rows=rows))
    motion_coupling=[]
    for d in (40.,100.):
        c=coupled_case(d/2,air[10],replace(reference,annual_dark_fraction=0,longest_dark_days=0))
        for cd in (.47,.05):
            for mismatch in (0.,.25,.5,1.,2.):
                motion=navigation.drag_power(mismatch,air[10]['density_kg_m3'],cd,.25)
                carbon_cost=motion['input_w_m2']*YEAR_S/40e6
                available=c['carbon']['remaining_assimilate_c_kg_m2_year']-carbon_cost
                repro=parent_funded_reproduction(available,reference.gpp_c_kg_m2_year,
                    c['reproduction']['bud_construction_and_reserve_c_kg_m2'],
                    c['hydrogen']['inventory_kg']/(math.pi*(d/2)**2)*H2_LHV_J_KG,
                    reference.incident_cycle_mean_w_m2*reference.assumed_full_solar_to_h2_lhv_efficiency,
                    c['hydrogen']['required_hydrogen_area_time_fraction'])
                motion_coupling.append(dict(diameter_m=d,drag_coefficient=cd,mismatch_m_s=mismatch,
                    available_before_motion_c_kg_m2_year=c['carbon']['remaining_assimilate_c_kg_m2_year'],
                    motion=motion,motion_carbon_equivalent_kg_m2_year=carbon_cost,
                    available_after_motion_c_kg_m2_year=available,reproduction=repro,
                    boundary='Continuous useful light and imposed vector mismatch. Actual route, propulsor mass/efficiency and latitude-dependent GPP remain unestablished.'))
    envelope_comparators=envelope.reference_thickness_scan()
    for row in envelope_comparators:
        row['gas_area_per_projected_area']=4.
        row['hydrogen_chemical_w_m2_projected']=row['hydrogen_flux_mol_m2_s']*4*envelope.HYDROGEN_MOLAR_MASS_KG*H2_LHV_J_KG
    structural=structure.evaluate(round_spare_pa=SPARE_TRIM_PA,module_spare_pa=SPARE_TRIM_PA_COLONY)
    growth_eval=growth.evaluate()
    sailing_eval=sailing.evaluate()
    periods={f'equator_riding_{z:g}_km_mean_wind':r['relative_period_days'] for z in (10.,20.)
             for r in sailing_eval['day_length'] if r['latitude_deg']==0. and r['height_km']==z}
    water_eval=water.evaluate(relative_periods_days=periods)
    storms_eval=storms.evaluate()
    nursery=canopy_nursery(air,structural,growth_eval)
    limits=size_limits(air,structural)
    ages=growth_tables(air,structural,nursery)
    components=('mechanics','envelope','biology','environment','navigation','photoperiod','trim_cycle',
                'structure','growth','water','storms','sailing','gas_biology','resources')
    producers=[HERE/'run.py']+[HERE/f'{name}.py' for name in components]
    sources=[HERE/f'{name}_sources.json' for name in components if (HERE/f'{name}_sources.json').exists()]
    modules=(environment,photoperiod,navigation,structure,growth,water,storms,sailing,resources)
    inputs=list(dict.fromkeys([str(p.relative_to(ROOT)) for p in sources]+[f for m in modules for f in m.INPUT_FILES]))
    def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
    out=dict(schema=SCHEMA,
        producer=dict(files={str(p.relative_to(ROOT)):sha(p) for p in producers},
            inputs={p:sha(ROOT/p) for p in inputs},constants=constants_used(producers)),
        evidence='Coupled analytic requirements using inherited air and measured material/physiology comparators. No integrated organism, material trait combination, water route, flight controller or life cycle demonstrated.',
        reading_rule='Diameter and radius both explicit. Organism-area units use projected piR²; envelope area is separate. Mass includes hydrated tissue, inert structure, free water and night reserve. Carbon is post-maintenance assimilate when labelled; no reuse of NPP for baseline respiration.',
        assumptions=dict(reference=asdict(reference),hydrogen_lhv_j_kg=H2_LHV_J_KG,
            construction_assimilate_c_per_kg_dry=CONSTRUCTION_C_PER_KG,
            hydrogen_night_mode='Day solar splitting and dark glucose fermentation use each explicitly imposed light regime. Shared carbon/H2-dark duty is a sensitivity. Separation, organs and pumping unsolved.',
            light_regimes='The14.77-day dark interval is a stress test. Short-deficit and Sun-following cases are evaluated separately at fixed supplied GPP and cycle-mean irradiance to isolate timing costs, not predict their latitude-dependent productivity.',
            material_combinations='Permeability cases bracket measured conditions; strength, density, thickness, turnover, humidity protection and biological synthesis are independently imposed. No family is an observed membrane.'),
        environment=environment.evaluate(),biology=biology.evaluate(),gas_alternatives=gas_alternatives(air[10]),
        photoperiod=photoperiod.evaluate(),navigation=navigation.evaluate(),motion_coupling=motion_coupling,
        trim_cycle=trim_cycle.evaluate(air[10]['density_kg_m3'],
            air[10]['density_kg_m3']-air[10]['lift_hydrogen_kg_m3'],air[10]['pressure_atm']*101325.),
        envelope_comparators=envelope_comparators,cases=cases,size_scans=scans,structural_size_scans=structural_scans,
        light_regimes=light_regimes,dark_ballast_cases=dark_ballast_cases,sensitivities=sensitivities,
        structure=structural,size_limits=limits,growth=growth_eval,canopy_nursery=nursery,ages=ages,
        water=water_eval,storms=storms_eval,lifespans=lifespans(storms_eval),sailing=sailing_eval,
        gas_biology=gas_biology.evaluate(),resources=resources.evaluate(),
        design_guesses=dict(compartments_round=COMPARTMENTS_ROUND,spare_trim_pa=SPARE_TRIM_PA,
            spare_trim_pa_colony=SPARE_TRIM_PA_COLONY,skin_water_kg_m2=SKIN_WATER_KG_M2,
            radiant_fraction=RADIANT_FRACTION,fire_burn_s=FIRE_BURN_S,fire_distance_radii=FIRE_DISTANCE_RADII,
            dumpable_water_kg_m2=DUMPABLE_WATER_KG_M2,module_radius_m=MODULE_RADIUS_M,
            ranges=dict(compartments_round=[4,16],spare_trim_pa=[100.,1000.],spare_trim_pa_colony=[100.,300.],
                        skin_water_kg_m2=[.1,2.],radiant_fraction=[.03,.16],fire_burn_s=[30.,120.],
                        fire_distance_radii=[1.,2.],
                        joint_efficiency=[.5,.9],film_renewal_years=[1.,30.],
                        tendon_renewal_years=[30.,300.],tendon_coupon_mpa=[100.,1000.]),
            gas_cells='Slack gas cells inside the pressure-bearing hull: a lost cell fills with ballonet air at the hull pressure, so walls carry no crown pressure.'),
        no_claims=['self-sustaining species','guaranteed harvest','resolved gas purity','population carrying capacity','storm climatology over centuries'])
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    return out


if __name__=='__main__':
    run()
