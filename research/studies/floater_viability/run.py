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
from research.studies.floater_viability import biology, envelope, environment, mechanics, navigation, photoperiod, trim_cycle

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.floater-viability/1'
YEAR_S = JULIAN_DAY * JULIAN_YEAR_DAYS
H2_LHV_J_KG = 120e6  # rounded chemical-energy basis of the explicitly assumed solar efficiency
CONSTRUCTION_C_PER_KG = 1.39 * .4  # glucose construction analogue; material-specific chemistry unresolved
OUT = HERE / 'results/floater_viability.json'


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
    t=traits; area=math.pi*radius_m**2; volume=4*area*radius_m/3
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
            min(2.,radius_m),t.allowable_mpa*1e6,t.material_density_kg_m3,
            100e6,t.material_density_kg_m3,t.structural_minimum_um*1e-6,barrier_mass_area)
        exterior=skin['film_area_m2']
    else:
        raise ValueError('Unknown architecture')
    partition=mechanics.partition_mass(radius_m,barrier_mass_area+t.partition_structural_kg_m2)
    structure=skin['total_structure_mass_kg']+partition
    wetted=exterior+area
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
    actual_mass=structure+area*(tissue_wet+t.free_water_kg_m2+stored_wet)
    trim=mechanics.neutral_trim(radius_m,rho,gas_density,actual_mass)
    utilization=trim['required_lifting_mix_volume_fraction']
    structural_dry=structure/area*t.structural_dry_fraction
    turnover_dry=structural_dry*t.inert_turnover_per_year+t.living_dry_kg_m2*t.living_turnover_per_year
    turnover_assimilate=turnover_dry*CONSTRUCTION_C_PER_KG
    available=(photo['remaining_gpp_c_kg_m2_year']-gross['maintenance_c_kg_m2_year']
               -night_feed['feed_carbon_kg']-t.consumer_food_c_kg_m2_year-turnover_assimilate)
    out=dict(radius_m=radius_m,diameter_m=2*radius_m,height_km=air['height_km'],traits=asdict(t),
        pressure=loads,skin=skin,partition_mass_kg=partition,gas_contact_area_m2=wetted,
        mass=dict(structure_kg=structure,structural_dry_kg_m2=structural_dry,
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
    producers=[HERE/f for f in ('run.py','mechanics.py','envelope.py','biology.py','environment.py','navigation.py','photoperiod.py','trim_cycle.py')]
    sources=[HERE/f'{name}_sources.json' for name in ('mechanics','envelope','biology','environment','navigation','photoperiod')]
    inputs=list(dict.fromkeys([str(p.relative_to(ROOT)) for p in sources]+list(environment.INPUT_FILES)+list(photoperiod.INPUT_FILES)+list(navigation.INPUT_FILES)))
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
        no_claims=['safe size range','self-sustaining species','guaranteed harvest','resolved gas purity','storm survival','population carrying capacity'])
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    return out


if __name__=='__main__':
    run()
