"""Aerophyte water and heat requirements on committed lunar environmental products.

Pure component functions plus evaluate(root). No climate, weather trajectory,
altitude light field, or physiological acclimation is simulated here.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

from shared.constants import (JULIAN_DAY, JULIAN_YEAR_DAYS,
                              STEFAN_BOLTZMANN, SYNODIC_MONTH_DAYS)

ROOT = Path(__file__).resolve().parents[3]
INPUT_FILES = (
    'research/studies/sky_ships/results/sky_ships.json',
    'climate/results/crm/ring_ring_equator.json',
    'biosphere/ecology/results/people.json',
    'illumination/surface_light/results/surface_light.json',
    'biosphere/canopy/results/plant.json',
    'biosphere/canopy/results/canopy.json',
    'research/studies/aerophytes/environment_sources.json',
)
WATER_MOLAR_KG = 0.01801528
CARBON_MOLAR_KG = 0.012011
WATER_CP_J_KG_K = 4180.0  # rounded liquid-water engineering approximation
LATENT_J_KG = 2.45e6       # constant near room temperature; not an ice calculation


def saturation_pa(temperature_c):
    """Buck 1996 saturation over LIQUID water, -40..80 C; no pressure enhancement.

    Negative temperatures describe hypothetical supercooled water, not proof that
    tissues remain unfrozen. Equation documented by NCAR EOL (Vömel 2025).
    """
    if not -40 <= temperature_c <= 80:
        raise ValueError('liquid-water approximation restricted to -40..80 C')
    return 611.21 * math.exp((18.678-temperature_c/234.5)*temperature_c/(257.14+temperature_c))


def series_conductance(leaf_mol_m2_s, boundary_mol_m2_s):
    """Conductances on the SAME one-sided leaf reference area, in series."""
    if leaf_mol_m2_s < 0 or boundary_mol_m2_s <= 0:
        raise ValueError('nonnegative leaf and positive boundary conductance required')
    return leaf_mol_m2_s*boundary_mol_m2_s/(leaf_mol_m2_s+boundary_mol_m2_s)


def transpiration(air_c, leaf_c, relative_humidity, pressure_pa,
                  leaf_conductance_mol_m2_s, boundary_conductance_mol_m2_s=0.5):
    """Vapour loss per ONE-SIDED leaf area; wet interior at saturation.

    g_total * (e_leaf-e_air)/P, dilute-vapour Fick approximation. Below dewpoint
    loss is zero; condensation is flagged, never silently credited as usable water.
    """
    if not 0 <= relative_humidity <= 1 or pressure_pa <= 0:
        raise ValueError('physical humidity and pressure required')
    es = saturation_pa(leaf_c)
    ea = relative_humidity*saturation_pa(air_c)
    if max(es, ea) >= pressure_pa:
        raise ValueError('vapour pressure must be below total pressure')
    g = series_conductance(leaf_conductance_mol_m2_s, boundary_conductance_mol_m2_s)
    mol_s = g*max(0.0, es-ea)/pressure_pa
    kg_s = WATER_MOLAR_KG*mol_s
    return dict(vpd_pa=es-ea, total_conductance_mol_m2_s=g,
                water_kg_m2_day=kg_s*JULIAN_DAY, latent_w_m2=kg_s*LATENT_J_KG,
                condensation_possible=es < ea)


def carbon_water_requirement(npp_g_c_m2_year, air_c, leaf_c, relative_humidity,
                             pressure_pa, co2_ppm=400., ci_ca=0.7,
                             assimilation_to_npp=2.):
    """C3 stomatal-diffusion requirement, per projected organism area.

    E/A=1.6*VPD/(p_CO2*(1-Ci/Ca)). A is integrated net LEAF CO2 assimilation,
    not gross carboxylation. Its ratio to whole-organism NPP includes non-leaf
    maintenance and growth costs. Boundary layer, cuticular loss and CAM omitted.
    """
    if npp_g_c_m2_year < 0 or not 0 <= relative_humidity <= 1 or pressure_pa <= 0:
        raise ValueError('nonnegative production and physical air required')
    if co2_ppm <= 0 or not 0 <= ci_ca < 1 or assimilation_to_npp < 1:
        raise ValueError('positive CO2, Ci/Ca<1 and assimilation/NPP>=1 required')
    vpd = max(0., saturation_pa(leaf_c)-relative_humidity*saturation_pa(air_c))
    mol_water_mol_c = 1.6*vpd/(co2_ppm*1e-6*pressure_pa*(1-ci_ca))
    water = npp_g_c_m2_year/1000/CARBON_MOLAR_KG*assimilation_to_npp*mol_water_mol_c*WATER_MOLAR_KG
    return dict(water_kg_m2_year=water, water_kg_m2_lunar_cycle=water*SYNODIC_MONTH_DAYS/JULIAN_YEAR_DAYS,
                mol_water_per_mol_leaf_net_c=mol_water_mol_c,
                assimilation_to_npp=assimilation_to_npp, vpd_pa=vpd)


def water_for_assimilation(net_leaf_g_c_m2_year, air_c, leaf_c, relative_humidity,
                          pressure_pa, co2_ppm=400., ci_ca=.7):
    """Direct coupling for annual NET leaf CO2 uptake per organism footprint.

    If the caller supplies canopy GPP, subtract leaf respiration in its carbon
    account first. Whole-organism NPP must not be substituted without restoring
    stem/root and growth costs; microbial/CAM pathways require their own models.
    """
    out=carbon_water_requirement(net_leaf_g_c_m2_year,air_c,leaf_c,relative_humidity,
                                pressure_pa,co2_ppm,ci_ca,assimilation_to_npp=1.)
    out.pop('assimilation_to_npp')
    return dict(net_leaf_g_c_m2_year=net_leaf_g_c_m2_year,air_c=air_c,leaf_c=leaf_c,
                relative_humidity=relative_humidity,pressure_pa=pressure_pa,
                co2_ppm=co2_ppm,ci_ca=ci_ca,**out)


def water_cycle(day_loss_kg_m2_day, night_loss_kg_m2_day, leaf_area_ratio=1.,
                rainless_days=SYNODIC_MONTH_DAYS, water_store_kg_m2=9., usable_share=0.2):
    """Equal day/night and specified no-refill duration; no temporal weather claim.

    Water store is per projected organism area; only usable_share is expendable.
    Total tissue water cannot all be consumed while living tissues remain hydrated.
    """
    if min(day_loss_kg_m2_day, night_loss_kg_m2_day, leaf_area_ratio,
           rainless_days, water_store_kg_m2) < 0 or not 0 <= usable_share <= 1:
        raise ValueError('nonnegative stocks/flows and a usable fraction required')
    daily = (day_loss_kg_m2_day+night_loss_kg_m2_day)*leaf_area_ratio/2
    return dict(mean_loss_kg_m2_day=daily, demand_kg_m2=daily*rainless_days,
                usable_store_kg_m2=water_store_kg_m2*usable_share,
                autonomy_days=None if daily == 0 else water_store_kg_m2*usable_share/daily,
                required_total_store_kg_m2=None if usable_share == 0 else daily*rainless_days/usable_share)


def fog_collection(liquid_g_m3, relative_speed_m_s, capture_fraction,
                   collector_area_ratio=1., fog_duty=1.):
    """Collected existing DROPLETS per projected area; no condensation heat charged.

    Ambient liquid content, relative flow and area must be jointly attainable.
    Aerodynamic work, leaf wetting, salinity and biological uptake remain external.
    """
    if min(liquid_g_m3, relative_speed_m_s, collector_area_ratio) < 0:
        raise ValueError('nonnegative fog inputs required')
    if not 0 <= capture_fraction <= 1 or not 0 <= fog_duty <= 1:
        raise ValueError('fractions outside [0,1]')
    flux = liquid_g_m3/1000*relative_speed_m_s*capture_fraction*collector_area_ratio
    return dict(in_fog_kg_m2_hour=flux*3600, mean_kg_m2_day=flux*JULIAN_DAY*fog_duty)


def dew_collection(rejection_w_m2, leaf_c, air_c, relative_humidity):
    """Energy ceiling only: condensed vapour releases latent heat.

    Rejection is NET cooling left after sensible/radiative gains. No mass-transfer
    supply or area benefit is assumed. Condensation requires a sub-dewpoint surface.
    """
    if rejection_w_m2 < 0 or not 0 <= relative_humidity <= 1:
        raise ValueError('positive cooling and physical humidity required')
    possible = saturation_pa(leaf_c) < relative_humidity*saturation_pa(air_c)
    return dict(below_dewpoint=possible,
                energy_ceiling_kg_m2_day=rejection_w_m2*JULIAN_DAY/LATENT_J_KG if possible else 0.)


def rain_load(rain_mm_h, duration_h, catch_share=1., drainage_kg_m2_h=0.):
    """Freshwater 1 mm =1 kg/m²; horizontal footprint, initially empty catchment."""
    if min(rain_mm_h, duration_h, drainage_kg_m2_h) < 0 or not 0 <= catch_share <= 1:
        raise ValueError('nonnegative rain/drainage required')
    incoming = rain_mm_h*catch_share
    return dict(incoming_kg_m2_h=incoming,
                accumulated_kg_m2=max(0., incoming-drainage_kg_m2_h)*duration_h,
                minimum_steady_drainage_kg_m2_h=incoming)


def leaf_heat_requirement(air_c,radiation_c,target_c,relative_humidity,pressure_pa,
                          leaf_conductance_mol_m2_s,sensible_w_m2_k=10.,emissivity=.96):
    """Net absorbed/heater power needed at target leaf temperature, same area as leaf_equilibrium."""
    if sensible_w_m2_k <= 0 or not 0 < emissivity <= 1:
        raise ValueError('positive sensible transfer and physical emissivity required')
    water=transpiration(air_c,target_c,relative_humidity,pressure_pa,leaf_conductance_mol_m2_s)
    radiation=2*emissivity*STEFAN_BOLTZMANN*((target_c+273.15)**4-(radiation_c+273.15)**4)
    sensible=sensible_w_m2_k*(target_c-air_c)
    return dict(target_c=target_c,net_heat_required_w_m2=radiation+sensible+water['latent_w_m2'],
                radiative_loss_w_m2=radiation,sensible_loss_w_m2=sensible,**water)


def leaf_equilibrium(air_c, radiation_c, absorbed_w_m2, relative_humidity, pressure_pa,
                     leaf_conductance_mol_m2_s, sensible_w_m2_k=10., emissivity=.96):
    """Steady two-face leaf heat screen, all terms per one-sided leaf area.

    Each face sees the same effective IR radiation temperature; this is a stated
    boundary, NOT the temperature of deep space. h combines both faces; it is an
    imposed sensitivity, not derived from the aerophyte's bulk diameter or wind.
    Dew, heat conduction to body, biochemical heat and freezing are not solved.
    """
    if absorbed_w_m2 < 0 or sensible_w_m2_k <= 0 or not 0 < emissivity <= 1:
        raise ValueError('physical thermal inputs required')

    def terms(t):
        water = transpiration(air_c,t,relative_humidity,pressure_pa,leaf_conductance_mol_m2_s)
        rad = 2*emissivity*STEFAN_BOLTZMANN*((t+273.15)**4-(radiation_c+273.15)**4)
        sensible = sensible_w_m2_k*(t-air_c)
        residual = absorbed_w_m2-rad-sensible-water['latent_w_m2']
        return residual, rad, sensible, water

    lo, hi = -40., 80.
    if terms(lo)[0] < 0 or terms(hi)[0] > 0:
        raise ValueError('thermal root outside liquid-water calculation range')
    for _ in range(70):
        mid = (lo+hi)/2
        if terms(mid)[0] > 0: lo = mid
        else: hi = mid
    temperature = (lo+hi)/2
    residual, rad, sensible, water = terms(temperature)
    return dict(leaf_c=temperature, radiative_loss_w_m2=rad,
                sensible_loss_w_m2=sensible, energy_residual_w_m2=residual,
                liquid_state_unestablished=temperature < 0,
                dew_heat_omitted=water['condensation_possible'], **water)


def thermal_store(water_kg_m2, allowable_cooling_k, net_loss_w_m2,
                  duration_days=SYNODIC_MONTH_DAYS/2):
    """Sensible heat in liquid water only; no freezing/insulation/metabolic credit."""
    if min(water_kg_m2,allowable_cooling_k,net_loss_w_m2,duration_days) < 0:
        raise ValueError('nonnegative thermal inputs required')
    energy = water_kg_m2*WATER_CP_J_KG_K*allowable_cooling_k
    demand = net_loss_w_m2*duration_days*JULIAN_DAY
    return dict(stored_heat_j_m2=energy, demanded_heat_j_m2=demand,
                covers_interval=energy >= demand,
                autonomy_days=None if net_loss_w_m2 == 0 else energy/net_loss_w_m2/JULIAN_DAY,
                required_water_kg_m2=None if allowable_cooling_k == 0 else demand/(WATER_CP_J_KG_K*allowable_cooling_k))


def evaluate(root=ROOT):
    """Return JSON-compatible requirements and inherited environment facts."""
    root=Path(root)
    data=[json.loads((root/p).read_text()) for p in INPUT_FILES[:-1]]
    ships, ring, people, light, plant, canopy=data
    expected=('terluna.research.sky-ships/1','terluna.climate.crm-ring/1',
              'terluna.biosphere.ecology-people/1','terluna.illumination.surface-light/1',
              'terluna.biosphere.plant-carbon-cycle/1',
              'terluna.biosphere.canopy-photosynthesis/1')
    for d,s in zip(data,expected):
        if d.get('schema') != s:
            raise ValueError(f'Unexpected environment schema {d.get("schema")}, expected {s}')
    air=[a for a in ships['air'] if a['height_km'] in (0,10,20,25,30,40)]
    ten=next(a for a in air if a['height_km']==10)
    temperature=ten['temperature_c']; pressure=ten['pressure_atm']*101325
    evaporation=[]; carbon=[]
    for a in air:
        if a['height_km'] not in (10,20,40): continue
        for rh in (.6,.9,.98):
            for warming in (0.,5.):
                # Negative mean air needs an explicit warmed, liquid tissue scenario.
                leaf=max(a['temperature_c'],5.)+warming
                for g in (.0016,.0096,.05,.1,.3):
                    evaporation.append(dict(height_km=a['height_km'],relative_humidity=rh,
                        leaf_c=leaf,leaf_conductance_mol_m2_s=g,
                        **transpiration(a['temperature_c'],leaf,rh,a['pressure_atm']*101325,g)))
                for ratio in (1.,2.,3.):
                    for co2 in (400.,1000.,2000.):
                        carbon.append(dict(height_km=a['height_km'],relative_humidity=rh,leaf_c=leaf,
                            npp_g_c_m2_year=1000.,air_c=a['temperature_c'],pressure_pa=a['pressure_atm']*101325,
                            co2_ppm=co2,co2_partial_pa=co2*1e-6*a['pressure_atm']*101325,
                            **carbon_water_requirement(1000,a['temperature_c'],leaf,rh,a['pressure_atm']*101325,
                                                       co2_ppm=co2,assimilation_to_npp=ratio)))
    day=transpiration(temperature,temperature+5,.6,pressure,.1)
    night=transpiration(temperature,temperature,.9,pressure,.0016)
    return dict(
        evidence='Environment requirements on inherited global-mean air and one equatorial CM1 ring. No aerophyte trajectory, measured lunar humidity/LWC, aerial radiative field, viable water cycle or thermal physiology is established.',
        reading_rule='Water kg per projected organism m² unless leaf area is stated. Conductance and leaf heat use ONE-SIDED leaf area; leaf-area ratio converts to organism footprint. carbon_water_requirements scales whole-organism NPP to net leaf uptake; gross_to_net_leaf_water subtracts an assumed leaf-respiration share from GPP.',
        assumptions=dict(latent_j_kg=LATENT_J_KG,water_heat_capacity_j_kg_k=WATER_CP_J_KG_K,
            boundary_conductance_mol_m2_s=.5,relative_humidity_grid=[.6,.9,.98],
            liquid_water_g_m3=[.05,.1,.5],effective_ir_offset_k=[0,-10],
            note='RH,LWC,leaf heating,conductance,area,duty,capture,h and IR temperature are independent requirements sensitivities, not a jointly validated organism/weather scenario.'),
        air=air,mean_freezing_height_km=people['sky']['freezing_height_km'],
        co2=dict(inherited_canopy_partial_pa=canopy['settings']['co2_pa'],
                 baseline_mixing_ppm=400.,sensitivity_mixing_ppm=[400.,1000.,2000.],
                 note='400ppm is the inherited column and canopy convention; higher surface pressure raises lunar CO2 partial pressure. The1000/2000ppm alternatives are physiological sensitivities, not adopted climate compositions.'),
        cloud_fractions={k:v['cloud_fraction_ring'] for k,v in people['sky']['bands'].items()},
        ring_water=ring['water_budget'],ring_storms=ring['storms'],ring_storm_tracks=ring['storm_tracks'],
        light=dict(ground_overhead=light['cases']['moon_1.2atm/design']['summaries']['90.0']['0.1'],
                   ground_twilight=plant['twilight']['moon'],
                   note='Ground spectra and plant twilight are context only. No altitude light gain or darkness reduction is used to validate aerophyte production.'),
        transpiration=evaporation,carbon_water_requirements=carbon,
        gross_to_net_leaf_water=[dict(gross_g_c_m2_year=gpp,net_leaf_share=f,
            **water_for_assimilation(gpp*f,temperature,temperature+5,.6,pressure,co2_ppm=co2))
            for gpp in (2000.,5000.) for f in (.7,.9) for co2 in (400.,1000.,2000.)],
        storage=dict(reference_day=day,reference_night=night,
            cases=[dict(leaf_area_ratio=l,water_store_kg_m2=w,
                        **water_cycle(day['water_kg_m2_day'],night['water_kg_m2_day'],l,
                                      water_store_kg_m2=w)) for l in (1.,3.) for w in (9.,90.)]),
        fog=[dict(liquid_g_m3=c,relative_speed_m_s=v,capture=.3,fog_duty=d,
                  **fog_collection(c,v,.3,fog_duty=d)) for c in (.05,.1,.5) for v in (0.,1.,3.) for d in (.01,.1,.5)],
        dew=[dict(net_rejection_w_m2=q,leaf_c=5.,air_c=temperature,relative_humidity=.6,
                  **dew_collection(q,5.,temperature,.6)) for q in (10.,50.)],
        rain=[dict(rate_mm_h=r,duration_h=h,**rain_load(r,h))
              for r in (ring['storms']['rain_max_mm_h']['median'],ring['storms']['rain_max_mm_h']['p90'],
                        ring['storms']['rain_max_mm_h']['max']) for h in (.1,1.)],
        thermal=[dict(air_c=temperature,relative_humidity=.6,absorbed_w_m2=q,
                      sensible_w_m2_k=h,leaf_conductance_mol_m2_s=g,
                      **leaf_equilibrium(temperature,temperature,q,.6,pressure,g,h))
                 for q in (100.,300.,600.) for h in (2.,10.,30.) for g in (.0016,.1)],
        night_thermal=[dict(air_c=a['temperature_c'],radiation_c=a['temperature_c']-10,
                            **leaf_equilibrium(a['temperature_c'],a['temperature_c']-10,0,.6,
                                               a['pressure_atm']*101325,.0016,10.))
                       for a in air if a['height_km'] in (10,20,40)],
        high_air_liquid_heat=[dict(height_km=a['height_km'],sensible_w_m2_k=h,
             **leaf_heat_requirement(a['temperature_c'],a['temperature_c'],5.,.6,
                                     a['pressure_atm']*101325,.0016,h))
             for a in air if a['height_km']==40 for h in (2.,10.,30.)],
        heat_storage=[dict(water_kg_m2=w,net_loss_w_m2=q,allowable_cooling_k=10.,
                           **thermal_store(w,10.,q)) for w in (9.,90.) for q in (.1,1.,10.)],
        unresolved=['Lagrangian RH,cloud-liquid/ice,rain and shear histories along controllable biological paths.',
                    'Leaf/substrate optics and IR sky field at altitude; heat and gas exchange under actual leaf geometry.',
                    'Closed-stomata permeability at cold/UV-filtered conditions, water-uptake rate, desiccation/ice tolerance.',
                    'Joint carbon-water-lift budget, rain shedding and drainage, water compartment failures, cloud depletion.',
                    'GCM programme remains paused; no new climate run made.'])
