"""Comparable source-assumption accounts and explicit electricity sensitivities."""
import json
import resource
import signal
import time

import numpy as np
from scipy.interpolate import CubicHermiteSpline

from shared import constants as K
from shared.provenance import constants_used
from engineering.shield_power import routed_power, optical_requirement, exhaust_rate
from protection.dynamics.natural_pattern import retarded_sun
from protection.dynamics.optical import length
from .packing_screen import load_raw
from .repair_screen import source_parents, write_product
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE, SCENARIO


def main():
    before=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    names=['holding.json','closure_budget.json','closure_validation.json','closure_sizing.json','pilot_budget.json']
    p,sources=source_parents(names,['engineering/shield_power.py',
        'research/studies/solar_shield_array/repair_trade.py'])
    holding=p['holding.json'];fraction=holding['core_diverted_fraction']
    benchmark=holding['aperture_quadrature']['variable_distance_finer']
    area=benchmark['aperture_area_m2'];moon=np.pi*K.MOON_RADIUS**2
    extra=SCENARIO['climate_flux_before_dimming_W_m2']*(1-SCENARIO['climate_sunlight_scale'])
    rows=[]
    for label,a in [('lunar_disk',moon),('held_78000km_aperture',holding['aperture_quadrature']['baseline']['aperture_area_m2']),
                    ('held_selected_aperture',area),('same_361_clear_tiles',361*9890.**2)]:
        rows.append(dict(label=label,area_m2=a,face_on_first_intercept_W=a*K.SOLAR_CONSTANT,
            extra_five_percent_dimming_W=a*extra,
            combined_filter_and_dimming_band_W=a*K.SOLAR_CONSTANT*fraction))
    # A fixed climate window enclosing the annual source/target cone is a
    # conservative allocation screen. Holding's area model permits additional
    # annular authority at each date; its ideal-force output is not harvested power.
    core=holding['trajectory_search']['variable_distance']['fine_summary']['core_area_fraction']
    broad=area*K.SOLAR_CONSTANT*((1-core)+core*fraction)
    out=dict(schema='terluna.research.shield-collection-trade/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),
            inputs={n:digest(HERE/'results'/n) for n in names}),
        accepted_cycle=False,accepted_fleet=False,new_target_adopted=False,
        original_target_W=23.835e12,full_cycle_energy_J=None,full_cycle_average_W=None,
        actual_fleet_collection_W=None,actual_fleet_delivered_W=None,
        one_AU_face_on_source_comparison=rows,
        normalization='Same 1361 W/m2 irradiance, face-on clear collecting area, no eclipse or mutual shading. The original 586 TW applies only to the lunar-disk extra dimming. These are source-assumption normalizations, not annual harvested-energy simulations.',
        broad_static_allocation=dict(climate_window_area_fraction=core,optically_available_W=broad,
            scope='Climate-window rejected band plus all photons outside a fixed enclosing climate window. The held dynamics allowed full annular optical authority. Collection of that light and compatibility with the ideal thrust are unverified.'),
        benchmark=dict(propulsion_W=benchmark['mean_power_TW']*1e12,total_mass_kg=benchmark['total_mass_kg'],
            propellant_kg_s=benchmark['propellant_kg_s'],storage_included=False),
        hardware_note='Moving trajectories carry their propagated per-member 20% fixed allocation and power mass; held benchmark carries its own power mass and seven-day propellant buffer. Neither is a collector bill of materials.')
    env=environment(3.);faceon=0.;intercept=0.;band=0.;stages=[]
    for k,label in enumerate(['hardware_service','hardware_departure']):
        entry=p['closure_validation.json']['runs'][k];raw=load_raw(entry)
        s=CubicHermiteSpline(raw['t'],raw['state'][:,:,:3],raw['state'][:,:,3:])
        ev=p['closure_budget.json']['datasets'][label]['evaluations'][-1]
        t=load_raw(ev)['t']
        power=np.array([np.sum(K.SOLAR_CONSTANT*(K.AU/length(retarded_sun(env.at(ti))-s(ti)))**2*9890.**2) for ti in t])
        f=float(np.trapezoid(power,t));a=ev['optical']['first_intercept_bolometric_equivalent_W'];b=ev['optical']['redirected_band_optical_W']
        faceon+=f;intercept+=a['energy_J'];band+=b['energy_J']
        stages.append(dict(label=label,duration_s=float(t[-1]-t[0]),first_intercept_energy_J=a['energy_J'],
            first_intercept_mean_W=a['mean_W'],redirected_band_energy_J=b['energy_J'],redirected_band_mean_W=b['mean_W'],
            clear_film_faceon_energy_J=f))
    local_fraction=intercept/faceon
    out['measured_twelve_hour_passage']=dict(stages=stages,first_intercept_energy_J=intercept,
        redirected_band_energy_J=band,faceon_unshadowed_energy_J=faceon,
        first_interception_fraction_of_same_clear_film=local_fraction,
        first_intercept_average_W=intercept/43200.,redirected_band_average_W=band/43200.,
        scope='All 361 members, six service hours plus six departure hours at actual dates. No full-cycle or inter-pattern factor has been measured.')
    n=p['pilot_budget.json']['rows'][0]['tiles'];full_band=n*9890.**2*K.SOLAR_CONSTANT*fraction
    full_light=n*9890.**2*K.SOLAR_CONSTANT
    ratio=np.array(p['closure_sizing.json']['per_member_total_to_optical_mass_ratio'])
    total_mass=n*.05*10000.**2*ratio.mean();optical_mass=n*.05*10000.**2
    static_same=area*K.SOLAR_CONSTANT*fraction
    out['capacity_relaxation_sensitivity']=dict(tiles=n,physical_area_m2=n*10000.**2,clear_area_m2=n*9890.**2,
        full_faceon_first_intercept_W=full_light,full_faceon_redirected_band_W=full_band,
        nominal_optical_mass_kg=optical_mass,nominal_installed_mass_kg=total_mass,
        nominal_extra_hardware_mass_kg=total_mass-optical_mass,
        fraction_needed_to_match_static_same_band=static_same/full_band,
        fraction_needed_to_match_static_broad_allocation=broad/full_band,
        passage_only_fraction_if_repeated_every_two_days=local_fraction/4,
        inter_pattern_retention_needed_for_two_day_passage_to_match_same_band=static_same/(full_band*local_fraction/4),
        scope='Define effective fraction f as cycle-mean projected first-intercept clear area divided by N times clear area, at 1 AU. It includes duty, orientation, body occultation and inter-pattern overlap. N is a capacity relaxation; neither placement nor f is established. Twelve-hour passage times N is only a thought experiment.')
    efficiencies=[dict(label='low',capture=.25,conversion=.20,delivery=.70),
                  dict(label='middle',capture=.50,conversion=.40,delivery=.85),
                  dict(label='high',capture=.80,conversion=.60,delivery=.95)]
    out['conversion_cases']=efficiencies;out['target_sensitivities']=[]
    prop=SCENARIO['propulsion']
    for load in [23.835e12,50e12,100e12,238.35e12]:
        fuel=exhaust_rate(load,prop['exhaust_velocity_m_s'],prop['efficiency'])
        row=dict(assumed_propulsion_W=load,benchmark_share=load/(benchmark['mean_power_TW']*1e12),
            exhaust_kg_s=fuel,exhaust_kg_year=fuel*K.JULIAN_DAY*K.JULIAN_YEAR_DAYS,conversion=[])
        for e in efficiencies:
            kwargs={k:e[k] for k in ['capture','conversion','delivery']}
            required=optical_requirement(100e12,**kwargs,propulsion_W=load)
            additional=(load-23.835e12)/(e['capture']*e['conversion'])
            row['conversion'].append(dict(label=e['label'],required_redirected_band_W=required,
                effective_fraction_for_100TW=required/full_band,
                extra_band_to_offset_higher_load_W=additional,
                extra_fraction_to_offset_higher_load=additional/full_band,
                effective_fraction_to_match_static_uniform_net=(static_same+(load-benchmark['mean_power_TW']*1e12)/(e['capture']*e['conversion']))/full_band,
                effective_fraction_to_match_static_broad_net=(broad+(load-benchmark['mean_power_TW']*1e12)/(e['capture']*e['conversion']))/full_band))
        out['target_sensitivities'].append(row)
    out['delivered_sensitivities']=[]
    out['static_electrical_sensitivities']=[dict(optical_allocation=label,conversion_case=e['label'],
        **routed_power(power,e['capture'],e['conversion'],e['delivery'],benchmark['mean_power_TW']*1e12))
        for label,power in [('uniform_band',static_same),('broad_annular_allowance',broad)] for e in efficiencies]
    for f in [.001,.01,.05,.10,.15]:
        for load in [23.835e12,50e12,100e12]:
            for e in efficiencies:
                for other in [0.,10e12]:
                    r=routed_power(full_band*f,e['capture'],e['conversion'],e['delivery'],load,other)
                    out['delivered_sensitivities'].append(dict(effective_fraction=f,conversion_case=e['label'],**r))
    held_dry=benchmark['total_mass_kg']-benchmark['buffer_mass_kg']
    year=K.JULIAN_DAY*K.JULIAN_YEAR_DAYS
    out['replacement_sensitivities']=[dict(lifetime_years=y,installed_replacement_kg_s=total_mass/(y*year),
        optical_replacement_kg_s=optical_mass/(y*year),held_dry_replacement_kg_s=held_dry/(y*year)) for y in [10,30,100]]
    out['material_trade']=dict(held_dry_hardware_mass_kg=held_dry,moving_installed_mass_kg=total_mass,
        comparison='Identical whole-hardware service lives and complete disposal, no recycling; exclude held propellant buffer from replacement to avoid charging exhaust twice. Moving inventory and hardware remain hypothetical.',
        rows=[dict(assumed_propulsion_W=r['assumed_propulsion_W'],
            same_lifetime_break_even_years=(total_mass-held_dry)/(benchmark['propellant_kg_s']-r['exhaust_kg_s'])/year
                if r['exhaust_kg_s']<benchmark['propellant_kg_s'] else None,
            annual_material_account=[dict(lifetime_years=y,
                moving_exhaust_plus_replacement_kg_s=r['exhaust_kg_s']+total_mass/(y*year),
                held_exhaust_plus_replacement_kg_s=benchmark['propellant_kg_s']+held_dry/(y*year)) for y in [10,30,100]])
            for r in out['target_sensitivities']])
    out['sensitivity_limits']='Efficiencies, f, collector layout and other loads are assumptions. Thermal-control electric power is included only in the stated other-load scenarios; collector heat and transfer loss are separate. Added collector/route mass and optical momentum have not been propagated, so these calculations accept no generation architecture or revised operating target.'
    out['completed']=True;write_product('repair_trade.json',out,sources,before)
    print(json.dumps(dict(local_fraction=local_fraction,static_uniform_band_PW=static_same/1e15,
        static_broad_PW=broad/1e15,capacity_band_PW=full_band/1e15,elapsed_s=time.monotonic()-before)),flush=True)


if __name__=='__main__':main()
