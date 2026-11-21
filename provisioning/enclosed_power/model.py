"""Vacuum-phase surface electricity alternatives on a common service basis.

The dispatch is a cyclic energy screen, not a terrain illumination, reactor,
grid-reliability or whole-economy simulation. All costs are conditional.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_pinned(entry):
    path = ROOT / entry['path']
    if digest(path) != entry['sha256']:
        raise ValueError(f"Input hash mismatch: {entry['path']}")
    if path.suffix == '.gz':
        with gzip.open(path, 'rt') as f:
            return json.load(f)
    return json.loads(path.read_text())


def profile(period_h, samples, site='ordinary', polar_lit=.85, polar_dark_h=100):
    phase = (np.arange(samples) + .5) / samples
    if site == 'ordinary':
        sun = (phase < .5).astype(float)
    elif site == 'polar_sensitivity':
        # Two synthetic shadows, not an observed site or an ephemeris.
        second_dark = (1 - polar_lit) * period_h - polar_dark_h
        if not 0 <= second_dark <= polar_dark_h:
            raise ValueError('Invalid synthetic polar outage pattern')
        dark = abs(phase - .45) < polar_dark_h / period_h / 2
        dark |= abs(phase - .90) < second_dark / period_h / 2
        sun = (~dark).astype(float)
    else:
        raise ValueError(site)
    return sun, period_h / samples


def shifted_load(sun, shift):
    """Move a fraction of dark-time demand to daylight, preserving total energy."""
    if not 0 <= shift <= 1:
        raise ValueError('Shift must be in [0,1]')
    lit = float(sun.mean())
    return np.where(sun > 0, 1 + shift * (1 - lit) / lit, 1 - shift)


def dispatch(sun, dt, load, nuclear_delivered_kW, solar_derate, grid_eff, round_trip):
    """One kW mean load; storage state is deliverable electrical energy.

    At the minimum PV size, integrated charging after round-trip losses exactly
    supplies deficits. The SOC range supplies the cyclic store, with no energy
    borrowed from the separately retained critical reserve.
    """
    if not 0 < round_trip <= 1 or not 0 < solar_derate <= 1 or not 0 < grid_eff <= 1:
        raise ValueError('Invalid efficiency')
    def increments(pv):
        surplus = pv * sun * solar_derate * grid_eff + nuclear_delivered_kW - load
        return np.where(surplus >= 0, surplus * round_trip, surplus), surplus
    lo, hi = 0., 1.
    if increments(0)[0].sum() > 1e-9:
        raise ValueError('Nuclear surplus requires a separate load-following case')
    while increments(hi)[0].sum() < 0:
        hi *= 2
    for _ in range(65):
        mid = (lo + hi) / 2
        if increments(mid)[0].sum() < 0:
            lo = mid
        else:
            hi = mid
    step, surplus = increments(hi)
    cumulative = np.r_[0., np.cumsum(step) * dt]
    return {'PV_nameplate_kW': hi, 'storage_output_kWh': float(np.ptp(cumulative)),
            'charge_kW': float(max(0., surplus.max())),
            'discharge_kW': float(max(0., -surplus.min())),
            'cycle_residual_kWh': float(cumulative[-1]),
            'cycle_load_kWh': float(load.sum() * dt),
            'solar_to_bus_kWh': float((hi * sun * solar_derate * grid_eff).sum() * dt),
            'nuclear_to_bus_kWh': nuclear_delivered_kW * len(sun) * dt,
            'storage_loss_kWh': float(np.maximum(surplus, 0).sum() * dt * (1 - round_trip)),
            'actual_lit_fraction': float(sun.mean()),
            'load_peak_kW': float(load.max())}


def emitting_area(heat_kW, temperature, constants, common):
    return heat_kW * 1000 / (common['radiator_emissivity'] * constants['physics']['stefan_boltzmann_W_m2_K4'] * temperature**4 * common['radiator_effective_fraction'])


def architecture(name, tech, cfg, eng, constants, *, mean_MW=1000,
                 firm_ratio=1.5, site='ordinary', shift=0, beam_useful=.2):
    if name not in {'beamed_legacy', 'solar_battery', 'solar_RFC', 'solar_thermal', 'hybrid_50', 'fission'}:
        raise ValueError(name)
    if mean_MW <= 0 or firm_ratio < 1 or not 0 < beam_useful <= 1:
        raise ValueError('Invalid capacity or beam utilization')
    c = cfg['common']; scale = mean_MW * 1000
    period = constants['moon']['synodic_month_days'] * constants['time']['julian_day_s'] / 3600
    sun, dt = profile(period, c['samples'], site, c['polar_lit_fraction'], c['polar_max_dark_hours'])
    load = shifted_load(sun, shift)
    # Changed production schedules need at least their new daytime peak.
    required_peak = max(firm_ratio, float(load.max()))
    components = []
    def add(label, mass, price, life, destination='surface'):
        components.append({'name': label, 'mass_kg': mass, 'supply_USD2026': price,
                           'life_years': life, 'destination': destination})
    def generator(label, mass, price, destination='surface'):
        share = c['generator_active_fraction']
        add(label + '_active', mass * share, price * share, c['generator_active_life_years'], destination)
        add(label + '_structure', mass * (1-share), price * (1-share), c['generator_structure_life_years'], destination)
    nuclear = 0.; details = {}
    if name == 'beamed_legacy':
        p = eng['industrial']['power']
        dc = required_peak / (p['RF_efficiency'] * p['beam_capture'] * p['receiver_efficiency'] * p['end_of_life_factor'] * beam_useful)
        mass = dc * 1000 / p['space_system_W_kg'] * tech['beam_mass_multiplier']
        aperture_floor = p['stations'] * math.pi * (p['aperture_m']/2)**2 * c['beam_aperture_floor_kg_m2'] / scale
        # Optimistic lower bound, not an extra duplicate subsystem charge.
        mass = max(mass, aperture_floor)
        price_kg = p['space_USD2026_W'] * p['space_system_W_kg'] * tech['beam_price_multiplier'] / tech['beam_mass_multiplier']
        generator('orbital_power', mass, mass * price_kg, 'orbit')
        receiver_mass = required_peak * 1000 * eng['demand']['power_hardware_kg_W']
        generator('ground_receiver', receiver_mass, required_peak * 1000 * eng['industrial']['prices']['ground_receiver_USD2026_W'])
        details = {'solar_DC_kW': dc, 'aperture_floor_kg_per_mean_kW': aperture_floor,
                   'aperture_floor_binding': mass <= aperture_floor * (1+1e-10),
                   'beam_useful_capacity_fraction': beam_useful,
                   'receiver_conversion_heat_kW_at_mean': 1/p['receiver_efficiency']-1,
                   'no_shield_credit': True}
    elif name == 'solar_thermal':
        eta = tech['thermal_engine_efficiency']
        retention = tech['thermal_store_retention']
        optical = tech['thermal_collector_efficiency'] * c['thermal_collector_derate']
        # Dispatch variable here is incident solar kW on collecting area.
        # Equivalent electrical SOC charges after the assumed cycle heat loss.
        d = dispatch(sun,dt,load,0,optical*eta,c['grid_efficiency'],retention)
        area = d['PV_nameplate_kW'] * 1000 / constants['sun']['solar_constant_W_m2']
        engine = required_peak / c['grid_efficiency']
        gross_heat = d['storage_output_kWh'] / (eta*c['grid_efficiency']*retention)
        media = gross_heat / tech['thermal_store_output_kWhth_kg']
        store_hardware = media * tech['thermal_store_hardware_kg_kg_media']
        generator('thermal_collectors',area*tech['thermal_collector_kg_m2'],area*tech['thermal_collector_USD2026_m2'])
        generator('heat_engine',engine*tech['thermal_engine_kg_kW'],engine*tech['thermal_engine_USD2026_kW'])
        add('thermal_media',media,media*tech['thermal_media_USD2026_kg'],None,'local_bulk')
        add('thermal_store_hardware',store_hardware,store_hardware*tech['thermal_store_hardware_USD2026_kg'],c['storage_tank_life_years'])
        rejected = engine*(1/eta-1)
        rad = emitting_area(rejected,tech['thermal_radiator_K'],constants,c)
        add('heat_engine_radiator',rad*c['radiator_kg_m2'],rad*c['radiator_USD2026_m2'],c['generator_active_life_years'])
        details = d | {'collector_area_m2':area,'gross_stored_heat_kWhth':gross_heat,
                       'thermal_media_kg':media,'thermal_store_hardware_kg':store_hardware,
                       'heat_engine_rated_kW':engine,'heat_engine_radiator_m2':rad,
                       'thermal_store_retention':retention,
                       'PV_nameplate_kW':None,
                       'incident_solar_nameplate_kW':d['PV_nameplate_kW'],
                       'meaning':'Cycle-equivalent retention screen, not a transient conduction or insulation design; storage medium is local bulk, not imported equipment.'}
    else:
        nuclear = {'fission': 1., 'hybrid_50': .5}.get(name, 0.)
        rated_nuclear = nuclear * required_peak / c['grid_efficiency']
        if nuclear:
            generator('fission_system', rated_nuclear * tech['nuclear_kg_kW'], rated_nuclear * tech['nuclear_USD2026_kW'])
            heat = rated_nuclear * (1/tech['nuclear_thermal_efficiency'] - 1)
            area = emitting_area(heat, tech['nuclear_radiator_K'], constants, c)
            details.update(nuclear_rated_kW=rated_nuclear, nuclear_conversion_heat_kW=heat,
                           nuclear_radiator_m2=area,
                           nuclear_radiator_mass_within_system_kg=area*c['radiator_kg_m2'])
        if name != 'fission':
            battery = name == 'solar_battery'
            efficiency = tech['battery_round_trip'] if battery else tech['FC_HHV_efficiency'] * tech['electrolyser_HHV_efficiency']
            d = dispatch(sun, dt, load, nuclear, tech['PV_output_derate'], c['grid_efficiency'], efficiency)
            generator('surface_PV', d['PV_nameplate_kW'] * tech['PV_kg_kW'], d['PV_nameplate_kW'] * tech['PV_USD2026_kW'])
            details.update(d)
            details['PV_area_m2'] = d['PV_nameplate_kW'] * 1000 / (constants['sun']['solar_constant_W_m2'] * c['PV_efficiency'])
            peak_store = max(d['discharge_kW'], (1-nuclear) * required_peak)
            if battery:
                add('battery_pack', d['storage_output_kWh']/tech['battery_usable_kWh_kg'], d['storage_output_kWh']*tech['battery_USD2026_usable_kWh'], c['battery_life_years'])
                rating = max(peak_store, d['charge_kW'])
                add('battery_converter', rating*tech['battery_converter_kg_kW'], rating*tech['battery_converter_USD2026_kW'], c['generator_active_life_years'])
                heat = max(peak_store*(1/math.sqrt(efficiency)-1), d['charge_kW']*(1-math.sqrt(efficiency)))
                area = emitting_area(heat, c['battery_temperature_K'], constants, c)
            else:
                hhv = eng['industrial']['source']['hydrogen_HHV_kWh_kg']
                hydrogen = d['storage_output_kWh'] / (hhv * tech['FC_HHV_efficiency'])
                oxygen = hydrogen * eng['industrial']['source']['oxygen_to_hydrogen_mass_ratio']
                water = hydrogen + oxygen
                tanks = hydrogen*tech['H2_tank_kg_kg'] + oxygen*tech['O2_tank_kg_kg'] + water*tech['water_tank_kg_kg']
                add('RFC_inventory', water, water*c['water_inventory_USD2026_kg'], None)
                add('RFC_tanks', tanks, tanks*tech['tank_USD2026_kg'], c['storage_tank_life_years'])
                add('RFC_fuel_cells', peak_store*tech['FC_kg_kW'], peak_store*tech['FC_USD2026_kW'], c['RFC_active_life_years'])
                add('RFC_electrolysers', d['charge_kW']*tech['electrolyser_kg_kW'], d['charge_kW']*tech['electrolyser_USD2026_kW'], c['RFC_active_life_years'])
                add('RFC_balance_of_plant', peak_store*tech['RFC_BOP_kg_kW'], peak_store*tech['RFC_BOP_USD2026_kW'], c['RFC_active_life_years'])
                heat = max(peak_store*(1/tech['FC_HHV_efficiency']-1), d['charge_kW']*(1-tech['electrolyser_HHV_efficiency']))
                area = emitting_area(heat, c['RFC_temperature_K'], constants, c)
                details.update(RFC_round_trip=efficiency, hydrogen_kg=hydrogen, oxygen_kg=oxygen,
                               water_equivalent_inventory_kg=water, tank_mass_kg=tanks,
                               initial_charge_input_kWh=d['storage_output_kWh']/efficiency)
            add('storage_thermal', area*c['radiator_kg_m2'], area*c['radiator_USD2026_m2'], c['generator_active_life_years'])
            details['storage_radiator_m2'] = area
            if not battery:
                store_mass = sum(x['mass_kg'] for x in components if x['name'].startswith('RFC_') or x['name']=='storage_thermal')
                details['RFC_complete_output_kWh_kg'] = d['storage_output_kWh']/store_mass if store_mass else 0
    return {'architecture': name, 'mean_MW':mean_MW, 'firm_to_mean':required_peak,
            'site':site,'night_load_shift':shift,'components_per_mean_kW':components,
            'mass_kg_per_mean_kW':sum(x['mass_kg'] for x in components),
            'supply_USD2026_per_mean_kW':sum(x['supply_USD2026'] for x in components),
            'nuclear_mean_share':nuclear,'details':details,
            'boundary':'Incremental generating and night-cycle storage systems. Common critical reserve, use-side heat and ordinary local grid retained separately. Fixed qualification, detailed installation, outage probabilities and whole-industry feedback unpriced.'}


def cost(case, cfg, constants, freight, imported):
    if not 0 <= imported <= 1:
        raise ValueError('Imported mass fraction must be in [0,1]')
    c=cfg['common']; comp=case['components_per_mean_kW']
    def delivery(item):
        dest=item['destination']
        if dest=='local_bulk':
            return item['mass_kg']*c['thermal_media_delivered_local_USD2026_kg']
        local=c['local_surface_delivery_USD2026_kg'] if dest=='surface' else c['local_orbit_delivery_USD2026_kg']
        return item['mass_kg']*(imported*freight[dest]+(1-imported)*local)
    capital=sum(x['supply_USD2026']+delivery(x) for x in comp)
    yearly=c['annual_upkeep_fraction']*case['supply_USD2026_per_mean_kW']
    for x in comp:
        loss = c['annual_thermal_media_replacement_fraction'] if x['name']=='thermal_media' else c['annual_RFC_inventory_loss_fraction'] if x['name']=='RFC_inventory' else 0
        yearly += loss * (x['supply_USD2026'] + delivery(x))
    annual_MWh=constants['time']['julian_year_days']*constants['time']['julian_day_s']/3600/1000
    yearly+=case['nuclear_mean_share']/c['grid_efficiency']*annual_MWh*c['nuclear_fuel_operations_USD2026_MWh']
    values={}
    for horizon in c['horizons_years']:
        pv=0.
        for x in comp:
            for y in range(0,horizon,x['life_years'] or horizon):
                pv+=(x['supply_USD2026']+delivery(x))/(1+c['discount_rate'])**y
        annuity=sum((1+c['discount_rate'])**(-y) for y in range(1,horizon+1))
        pv+=yearly*annuity
        values[str(horizon)]={'cost_PV_USD2026_per_mean_kW':pv,
                             'levelized_USD2026_MWh':pv/(annual_MWh*annuity)}
    return {'installed_supply_and_delivery_USD2026_per_mean_kW':capital,
            'annual_upkeep_and_fuel_USD2026_per_mean_kW':yearly,
            'imported_mass_fraction':imported, 'freight_USD2026_kg':freight,
            'imported_hardware_and_fluid_kg_per_mean_kW':sum(x['mass_kg'] for x in comp if x['destination']!='local_bulk')*imported,
            'local_hardware_and_fluid_kg_per_mean_kW':sum(x['mass_kg'] for x in comp if x['destination']!='local_bulk')*(1-imported),
            'always_local_bulk_kg_per_mean_kW':sum(x['mass_kg'] for x in comp if x['destination']=='local_bulk'),
            'lifecycle':values}


def run(cfg, data):
    eng=data['engineering']; constants=data['constants']; common=cfg['common']
    names=['beamed_legacy','solar_battery','solar_RFC','solar_thermal','hybrid_50','fission']
    comparisons=[]
    for tech_name,tech in cfg['technologies'].items():
        for name in names:
            case=architecture(name,tech,cfg,eng,constants,mean_MW=common['reference_mean_MW'],firm_ratio=common['reference_firm_to_mean'])
            case['technology_case']=tech_name
            case['costs']=[{'freight_case':f,**cost(case,cfg,constants,freight,imp)}
                           for f,freight in common['freight_scenarios'].items() for imp in common['mass_import_fractions']]
            comparisons.append(case)
    tech=cfg['technologies']['reference_assumptions']; freight=common['freight_scenarios']['reference']
    sensitivities=[]
    for site in ['ordinary','polar_sensitivity']:
        for shift in common['shift_fractions']:
            for name in ['solar_RFC','solar_battery','solar_thermal']:
                case=architecture(name,tech,cfg,eng,constants,site=site,shift=shift)
                case['costs']=[cost(case,cfg,constants,freight,imp) for imp in [1,.1]]
                sensitivities.append(case)
    beam_sens=[]
    for use in common['beam_useful_capacity_fractions']:
        case=architecture('beamed_legacy',tech,cfg,eng,constants,beam_useful=use)
        case['costs']=[cost(case,cfg,constants,freight,imp) for imp in [1,.1]]
        beam_sens.append(case)
    beam_curves=[]
    for use in np.linspace(.2,.6,41):
        case=architecture('beamed_legacy',tech,cfg,eng,constants,beam_useful=float(use))
        beam_curves.append({'useful_fraction':float(use),'costs':[
            {'imported_mass_fraction':imp,
             'installed_supply_and_delivery_USD2026_per_mean_kW':cost(case,cfg,constants,freight,imp)['installed_supply_and_delivery_USD2026_per_mean_kW']}
            for imp in [1,.1]]})
    aperture_sens=[]
    for density in common['beam_aperture_floor_sensitivities_kg_m2']:
        trial=deepcopy(cfg);trial['common']['beam_aperture_floor_kg_m2']=density
        case=architecture('beamed_legacy',tech,trial,eng,constants,mean_MW=1)
        case['aperture_floor_kg_m2']=density;case['cost']=cost(case,trial,constants,freight,1)
        aperture_sens.append(case)
    # Direct initial supply+delivery crossover, at the common large-system load.
    beamA=next(x for x in beam_sens if x['details']['beam_useful_capacity_fraction']==.2)
    beamB=next(x for x in beam_sens if x['details']['beam_useful_capacity_fraction']==.4)
    nuclear=next(x for x in comparisons if x['architecture']=='fission' and x['technology_case']=='reference_assumptions')
    thresholds=[]
    for imported in [1,.1]:
        ca=cost(beamA,cfg,constants,freight,imported)['installed_supply_and_delivery_USD2026_per_mean_kW']
        cb=cost(beamB,cfg,constants,freight,imported)['installed_supply_and_delivery_USD2026_per_mean_kW']
        cn=cost(nuclear,cfg,constants,freight,imported)['installed_supply_and_delivery_USD2026_per_mean_kW']
        slope=(ca-cb)/(1/.2-1/.4); fixed=ca-slope/.2
        rated=nuclear['details']['nuclear_rated_kW']
        hardware=tech['nuclear_USD2026_kW']*rated
        thresholds.append({'imported_mass_fraction':imported,
            'beam_useful_capacity_fraction_to_match_reference_fission':slope/(cn-fixed),
            'nuclear_supply_USD2026_kW_to_match_inherited_beam':(ca-(cn-hardware))/rated,
            'reading_rule':'Initial supply and delivery only; fixed qualification, site construction and lifecycle differences are not included in this crossover.'})
    pilots=[]
    for mean in common['pilot_mean_MW']:
        for name in ['beamed_legacy','solar_RFC']:
            case=architecture(name,tech,cfg,eng,constants,mean_MW=mean)
            case['cost']=cost(case,cfg,constants,freight,1)
            pilots.append(case)
    # Separate intermittent service: one MW during daylight, zero night supply.
    # The cost routine's kW-mean denominator is not used for these one-off stocks.
    daytime_pv_kW=1000/(tech['PV_output_derate']*common['grid_efficiency'])
    daytime_mass=daytime_pv_kW*tech['PV_kg_kW']
    daytime_supply=daytime_pv_kW*tech['PV_USD2026_kW']
    daytime={'lit_output_MW':1.,'cycle_mean_MW':.5,'PV_nameplate_MW':daytime_pv_kW/1000,
             'equipment_mass_kg':daytime_mass,
             'supply_and_delivery_USD2026':daytime_supply+daytime_mass*freight['surface'],
             'nighttime_electric_supply_MW':0.,
             'reading_rule':'Different service from the continuous comparisons: daylight-only output at an ordinary site. Survival heaters, short buffers and night standby must be added for actual equipment.'}
    physical=data['industry']['cases']
    reference=next(x for x in physical if x['name']==cfg['reference_case'])
    checkpoints=[]; mappings=[]
    for case in physical:
        if case['variant']!='combined':continue
        for row in case['physical_case']['annual']:
            if row['year']!=200 and not(case['name']==cfg['reference_case'] and row['year'] in [8,32,100,220]):continue
            mean=(eng['demand']['occupied_mean_W']*row['population']/eng['normalization']['population']+row['physical']['lunar_process_W'])/1e6
            firm=row['physical']['ground_firm_W']/1e6
            rec={'case':case['name'],'year':row['year'],'mean_MW':mean,'firm_MW':firm,
                 'population':row['population'],'maximum_shift_if_all_current_lunar_processing_moves':row['physical']['lunar_process_W']/1e6/mean}
            mappings.append(rec)
            if case['name']!=cfg['reference_case']:continue
            outcomes=[]
            for name in names:
                a=architecture(name,tech,cfg,eng,constants,mean_MW=mean,firm_ratio=firm/mean)
                cc=cost(a,cfg,constants,freight,1)
                outcomes.append({'architecture':name,
                                 'capacity_supply_and_delivery_USD2041':cc['installed_supply_and_delivery_USD2026_per_mean_kW']*mean*1000*eng['normalization']['USD2026_to_USD2041'],
                                 'system_mass_kg':a['mass_kg_per_mean_kW']*mean*1000,
                                 'imported_hardware_and_fluid_kg':cc['imported_hardware_and_fluid_kg_per_mean_kW']*mean*1000,
                                 'local_bulk_kg':cc['always_local_bulk_kg_per_mean_kW']*mean*1000})
            checkpoints.append(rec|{'comparative_capacity_costs':outcomes})
    first=reference['physical_case']['annual'][0]
    stock=first['asset_stocks_kg']['lunar_beam_constellation']
    quoted=stock*eng['templates']['lunar_beam_constellation']['price']*eng['normalization']['USD2026_to_USD2041']
    area_moon=4*math.pi*constants['moon']['radius_m']**2
    terminal=next(x for x in checkpoints if x['year']==200)
    return {'schema':'terluna.provisioning.enclosed-power/1','evidence':cfg['evidence'],
            'boundary':cfg['boundary'],'accounting':cfg['accounting'],
            'daytime_productive_pilot':daytime,'comparisons':comparisons,'load_and_site_sensitivities':sensitivities,
            'beam_utilization_sensitivities':beam_sens,'beam_cost_curves':beam_curves,'small_scale_beam_floor':pilots,
            'beam_aperture_floor_sensitivities':aperture_sens,'initial_cost_crossovers':thresholds,
            'global_context':{'year200_reference_mean_MW':terminal['mean_MW'],
                'year200_fission_total_heat_W_m2_global_equivalent':terminal['mean_MW']*1e6/common['grid_efficiency']/tech['nuclear_thermal_efficiency']/area_moon,
                'equatorial_circumference_km':2*math.pi*constants['moon']['radius_m']/1000,
                '100km_equatorial_longitude_solar_time_difference_hours':100e3/(2*math.pi*constants['moon']['radius_m'])*constants['moon']['synodic_month_days']*24,
                'meaning':'Area-averaged context, not a climate or local temperature calculation. No atmospheric or shield credit is used.'},
            'physical_load_mappings':mappings,'reference_checkpoints':checkpoints,
            'baseline_opening_stock_audit':{'surviving_beam_hardware_kg':stock,
               'at_inherited_equipment_quote_USD2041':quoted,
               'ratio_to_total_seed_outlay':quoted/reference['physical_case']['seed']['cost']['total_USD2041'],
               'reading_rule':'A surviving stock valued at the inherited purchase quote, not a reconstruction of gross seed spending or an avoidable cash total. Delivery, receivers and attrition are additional.'},
            'period_hours':constants['moon']['synodic_month_days']*24}


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,default=HERE/'results');args=ap.parse_args()
    inp=HERE/'enclosed_power_inputs.json'; cfg=json.loads(inp.read_text()); data={k:read_pinned(v) for k,v in cfg['inputs'].items()}
    result=run(cfg,data);result['dependencies']=cfg['inputs']|{'power_inputs':{'path':str(inp.relative_to(ROOT)),'sha256':digest(inp)}}
    result['producer']={'path':str(Path(__file__).resolve().relative_to(ROOT)),'sha256':digest(__file__)}
    args.out.mkdir(parents=True,exist_ok=True);(args.out/'power_comparison.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for case in result['comparisons']:
        if case['technology_case']=='reference_assumptions':
            co=next(c for c in case['costs'] if c['freight_case']=='reference' and c['imported_mass_fraction']==1)
            print(case['architecture'],round(case['mass_kg_per_mean_kW'],1),round(co['installed_supply_and_delivery_USD2026_per_mean_kW']/1000,2),round(co['lifecycle']['100']['levelized_USD2026_MWh']))


if __name__=='__main__':main()
