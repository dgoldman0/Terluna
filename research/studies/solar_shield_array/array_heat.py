"""Heat in the solar shield and habitat array: where it arises, how hot it is, what recovering work from it could give
and what rejecting it costs.

A closed-form screen for the planned array-industry branch (array_industry.md). It reads the integrated ledger (the
ring fleet's inventory, its film upkeep and the electricity cases), the attitude schemes (the momentum store) and the
protection domain's films (the sunlight the tiles absorb). It computes the films' thermal emissivity from the
protection domain's optical tables, and takes the rest from the assumptions in DESIGN, each sourced in
array_industry_sources.json. It runs no dynamics and no device model.

- tiles: the sunlight the ring fleet's films absorb, the temperature they run at, and the temperature difference
  across a film;
- collectors: the heat of converting light for the ledger's electricity cases;
- film plant: remaking the fleet's film every 10-20 years at glass-melting temperature, and what thermophotovoltaic
  cells facing the hot glass could return;
- computing: collector and radiator area and radiator mass per terawatt at a processor's rejection temperature;
- the shield's own loads: reflectivity trim by device type, and the momentum store by rotor size and rim speed;
- placement: how much of a terawatt used in orbit or on the Moon ends up in the Moon's heat budget.

    python -m research.studies.solar_shield_array.array_heat
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from shared import constants as K
from shared.provenance import constants_used
from engineering.heat_recovery import radiator_area
from protection import model as optics
from protection.spectra import short_wave as sw

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE/'results/array_heat.json'
FILES = ['research/studies/solar_shield_array/array_heat.py', 'engineering/heat_recovery.py']
LEDGER = HERE/'results/integrated_ledger.json'
ATTITUDE = HERE/'results/attitude_schemes.json'
FILMS = ROOT/'protection/spectra/annulus_film.json'
HELD = HERE/'results/held_eccentricity.json'
SOURCES = HERE/'array_industry_sources.json'
SCHEMA = 'terluna.research.array-heat/1'
YEAR = K.JULIAN_YEAR_DAYS*K.JULIAN_DAY

DESIGN = dict(
    orbit_radius_km=20000.,
    # The tiles: the comparison's coated annulus film and the window stack, titania and silica in micrometres.
    films=dict(annulus=[.1, 4.], window=[1., 10.]),
    infrared_um=[2., 125.],                 # the silica (Franta) and titania (Siefke) tables' reach
    lit_factor=[.5, 2/math.pi],             # orbit-mean projection of a radial-facing face on the Sun; shadows lower it
    silica_conductivity_W_m_K=1.,           # fused silica near 150-250 K (about 1.4 at room temperature)
    thermoradiative_measured_W_m2=9.41e-3,  # radchenkov_2025: best measured, a 292.5 K diode facing a 150 K emitter
    # Collectors converting light for the electricity cases: a multijunction array's solar absorptance and infrared
    # emittance (spectrolab_xtj), radiating from both faces.
    collector_absorptance=.88,
    collector_emissivity=[.85, .85],        # front, back
    # The film plant: the heat to raise silica into the melt from room temperature (2.0-2.3 MJ/kg at 2,000-2,200 K,
    # janaf_silica) up to evaporating it with a plant's losses (13.4-15 MJ/kg as the floor, janaf_silica).
    # Thermophotovoltaic cells facing the plant's hot surfaces capture a share of their emission (an assumption) and
    # convert 30-40% of it (lapotin_2022, tervo_2022, with cells cooled to 25 C; warmer cells convert less).
    film_specific_energy_MJ_kg=[2., 20.],
    film_hot_K=2000.,
    tpv_captured_share=[.3, .6],
    tpv_efficiency=[.3, .4],
    cell_K=330.,
    # Computing: processors rejecting at about 330 K from two-faced panels of 3 kg/m^2 (NASA's goal, el_genk_2025) to
    # 14 kg/m^2 (the ISS's panels, iss_atcs); collectors at 20-30%; present accelerators (nvidia_h100).
    radiator_K=330.,
    radiator_emissivity=.9,
    radiator_kg_m2=[3., 14.],
    collector_conversion=[.2, .3],
    accelerator_flop_per_J=1.4e12,
    # The shield's reflectivity trim: up to a quarter of the redirected band's reflectivity (photon_control,
    # attitude_schemes), from devices over a share of each tile. Three classes: liquid-crystal films that need power to
    # hold their mirror state, as IKAROS flew (1-5 W/m^2, ma_2017, pdlc_films); electrochromic films with memory,
    # switched (10-80 J/m^2 for a quarter's change, derived from coloration efficiency; window products take
    # 0.2-2 kJ/m^2 for a full tint, electrochromic_windows); and display-type bistable media (about 1-5 J/m^2 per
    # update, lin_2024). Switched 2-20 times an orbit.
    trim_fraction=.25,
    trim_device_share=[.25, .89],
    trim_held_W_m2=[1., 5.],
    trim_electrochromic_J_m2=[10., 80.],
    trim_bistable_J_m2=[1., 5.],
    trim_switches_per_orbit=[2., 20.],
    photovoltaic_on_tile=.25,
    # The momentum store: rotors spinning about axes in the tile's plane, their effective radius held by the
    # clearance to the next ring, cycled once an orbit.
    store_rotor_radius_m=[100., 400.],
    store_rim_speed_m_s=[300., 1000.],
    store_round_trip=[.85, .95],
)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def span(values):
    return [float(min(values)), float(max(values))]


def fleet(ledger):
    """The nestable ring fleet at 20,000 km with the stack's radius layout (both clearances)."""
    row = {r['radius_km']: r for r in ledger['ring_fleet']['by_radius']}[20000.]
    layout = row['stack_layout']['by_clearance']
    factors = [layout[c]['coupled']['tiles_over_one_radius'] for c in sorted(layout)]
    tiles = [row['tiles']*row['plane_motion']['ring_factor']*f for f in factors]
    mass = [row['matched_annulus']['silica_4_um']['oversized_optical_mass_Gt']*f for f in factors]
    side = ledger['requirements']['tile_side_m']
    years = ledger['requirements']['replacement_years']
    return dict(tiles=span(tiles), area_km2=span([t*side**2/1e6 for t in tiles]), optical_mass_Gt=span(mass),
                tile_mass_kg=float(np.mean(mass)*1e12/np.mean(tiles)),
                film_upkeep_Gt_per_year=[min(mass)/max(years), max(mass)/min(years)],
                window_ring_share=row['window_ring_share'])


def planck(w_um, temperature_K):
    w = np.asarray(w_um)*1e-6
    return w**-5/np.expm1(K.PLANCK*K.SPEED_OF_LIGHT/(w*K.BOLTZMANN*temperature_K))


def emissivity(titania_um, silica_um, temperature_K):
    """Normal thermal emissivity of a coated film from its sunward and shaded faces (1 - R - T, Planck-weighted over
    the tables' 2-125 um at the film's temperature), layered as the protection domain's films are."""
    w = np.geomspace(*DESIGN['infrared_um'], 4000)
    silica, titania = sw.indices(w*1e3)
    ar = sw.stored_ar_layers()
    porous, mixed = optics.mix_n(np.ones_like(silica), silica, .5), optics.mix_n(silica, titania, .5)
    layers = [porous, silica, mixed, titania, mixed, silica, porous]
    thickness = [ar[0], ar[1], ar[2], titania_um, ar[3], silica_um, ar[4]]
    weight = planck(w, temperature_K)
    faces = []
    for order in (1, -1):
        r, t = optics.stack_rt(w, layers[::order], thickness[::order])
        faces.append(float(np.trapezoid((1-r-t)*weight, w)/np.trapezoid(weight, w)))
    return faces, float(sum(ar)+titania_um+silica_um)


def tile(absorptance, titania_um, silica_um):
    """A tile facing the Sun at the sunward crossing: absorbed sunlight, the temperature at which both faces radiate it
    away, and the temperature difference that carries the shaded face's share through the film."""
    q = absorptance*K.SOLAR_CONSTANT
    t = 180.
    for _ in range(30):
        (front, back), thick = emissivity(titania_um, silica_um, t)
        t = (q/((front+back)*K.STEFAN_BOLTZMANN))**.25
    through = q*back/(front+back)
    across = through*thick*1e-6/DESIGN['silica_conductivity_W_m_K']
    return dict(solar_absorptance=absorptance, absorbed_W_m2=q, emissivity_sunward=front, emissivity_shaded=back,
                thickness_um=thick, temperature_K=t, across_film_K=across, carnot_across_film=across/t,
                carnot_work_W_m2=through*across/t)


def tiles(films, inventory):
    window = films['window_stack']['bolometric']['A']
    annulus = next(f for f in films['films'] if f['coated'] and f['titania_um'] == DESIGN['films']['annulus'][0]
                   and f['silica_um'] == DESIGN['films']['annulus'][1])['bolometric']['A']
    rows = dict(annulus=tile(annulus, *DESIGN['films']['annulus']), window=tile(window, *DESIGN['films']['window']))
    share = inventory['window_ring_share']
    mean = share*window+(1-share)*annulus
    area = np.array(inventory['area_km2'])*1e6
    lit = np.array(DESIGN['lit_factor'])
    intercepted = [float(area.min()*K.SOLAR_CONSTANT*lit.min()), float(area.max()*K.SOLAR_CONSTANT*lit.max())]
    carnot = max(r['carnot_work_W_m2'] for r in rows.values())
    return dict(films=rows, mean_solar_absorptance=mean,
                intercepted_PW=[v/1e15 for v in intercepted], absorbed_PW=[v*mean/1e15 for v in intercepted],
                carnot_across_all_films_TW=float(area.max()*carnot/1e12),
                thermoradiative_measured_everywhere_TW=float(area.max()*DESIGN['thermoradiative_measured_W_m2']/1e12))


def collectors(ledger):
    """The heat of the electricity cases' collectors, facing the Sun, and what the best measured thermoradiative
    diodes would return from their backs."""
    a = DESIGN['collector_absorptance']
    emit = sum(DESIGN['collector_emissivity'])*K.STEFAN_BOLTZMANN
    rows = []
    for case in ledger['electricity']['cases']:
        d, eta = case['demand_TW'], case['conversion']
        area = case['ring_fleet']['collector_km2_for_demand']*1e6
        rows.append(dict(demand_TW=d, conversion=eta, collector_km2=area/1e6, heat_TW=d*(a/eta-1),
                         heat_per_electric=a/eta-1,
                         temperature_K=((a-eta)*K.SOLAR_CONSTANT/emit)**.25,
                         thermoradiative_measured_GW=area*DESIGN['thermoradiative_measured_W_m2']/1e9))
    return rows


def film_plant(inventory):
    """Remaking the fleet's film: the plant's power, what thermophotovoltaic cells facing the hot glass could return,
    and the panel area that rejects the cells' heat."""
    upkeep = np.array(inventory['film_upkeep_Gt_per_year'])*1e12/YEAR            # kg/s
    energy = np.array(DESIGN['film_specific_energy_MJ_kg'])*1e6
    power = [float(upkeep.min()*energy.min()), float(upkeep.max()*energy.max())]
    share, eff = DESIGN['tpv_captured_share'], DESIGN['tpv_efficiency']
    recovered = [power[0]*share[0]*eff[0], power[1]*share[1]*eff[1]]
    cell_heat = [power[0]*share[0]*(1-eff[1]), power[1]*share[1]*(1-eff[0])]
    panel = [radiator_area(h, DESIGN['cell_K'], DESIGN['radiator_emissivity'])['emitting_area_m2']/2 for h in cell_heat]
    return dict(upkeep_kg_s=span(upkeep), power_TW=[p/1e12 for p in power],
                recovered_TW=[r/1e12 for r in recovered], recovered_share=[share[0]*eff[0], share[1]*eff[1]],
                cell_heat_TW=[h/1e12 for h in cell_heat], cell_panel_km2=[p/1e6 for p in panel],
                melt_radiation_W_m2=K.STEFAN_BOLTZMANN*DESIGN['film_hot_K']**4)


def computing():
    """One terawatt of processors: the collectors that power it and the two-faced panels that reject its heat."""
    rad = radiator_area(1e12, DESIGN['radiator_K'], DESIGN['radiator_emissivity'])
    panel = rad['emitting_area_m2']/2
    return dict(radiator_panel_km2_per_TW=panel/1e6, radiator_net_flux_W_m2=rad['net_flux_W_m2'],
                radiator_Mt_per_TW=[panel*m/1e9 for m in DESIGN['radiator_kg_m2']],
                collector_km2_per_TW=[1e12/(e*K.SOLAR_CONSTANT)/1e6 for e in DESIGN['collector_conversion'][::-1]],
                accelerator_flop_per_s_per_TW=DESIGN['accelerator_flop_per_J']*1e12)


def orbit_period_s():
    a = DESIGN['orbit_radius_km']*1e3
    return 2*math.pi*math.sqrt(a**3/K.MOON_GM)


def trim(inventory, films, held):
    """Reflectivity trim over a share of every tile: the fleet's power for each class of device, and the share of a
    tile that photovoltaics would need to supply it. Its duty: the trim changes the redirected
    band's reflection, which gives 2R of a film's normal pressure 2R + A, so it reaches trim 2R/(2R + A) of the push;
    holding the stack's eccentricities takes the share of the eccentricity drive held_eccentricity.json records."""
    area = np.array(inventory['area_km2'])*1e6
    share = np.array(DESIGN['trim_device_share'])
    devices = [float(area.min()*share.min()), float(area.max()*share.max())]
    period = orbit_period_s()
    n = DESIGN['trim_switches_per_orbit']
    supply = DESIGN['photovoltaic_on_tile']*K.SOLAR_CONSTANT*float(np.mean(DESIGN['lit_factor']))
    classes = dict(held=DESIGN['trim_held_W_m2'],
                   electrochromic=[DESIGN['trim_electrochromic_J_m2'][0]*n[0]/period,
                                   DESIGN['trim_electrochromic_J_m2'][1]*n[1]/period],
                   bistable=[DESIGN['trim_bistable_J_m2'][0]*n[0]/period, DESIGN['trim_bistable_J_m2'][1]*n[1]/period])
    out = {}
    for name, density in classes.items():
        out[name] = dict(W_m2=list(density), fleet_TW=[density[0]*devices[0]/1e12, density[1]*devices[1]/1e12],
                         photovoltaic_share_of_tile=[density[0]*share.min()/supply, density[1]*share.max()/supply])
    out['device_km2'] = [d/1e6 for d in devices]
    window = films['window_stack']['bolometric']
    annulus = next(f for f in films['films'] if f['coated'] and f['titania_um'] == DESIGN['films']['annulus'][0]
                   and f['silica_um'] == DESIGN['films']['annulus'][1])['bolometric']
    out['reach_share_of_push'] = {name: DESIGN['trim_fraction']*2*f['R']/(2*f['R']+f['A'])
                                  for name, f in (('annulus', annulus), ('window', window))}
    out['eccentricity_holding_share'] = {c: dict(max=v['share']['max'], mean=v['share']['mean'],
                                                 rings_over_a_tenth=v['share']['rings_over']['0.1'])
                                         for c, v in held['by_clearance'].items()}
    out['each_mW_m2_TW'] = [d*1e-3/1e12 for d in devices]
    return out


def store(attitude, inventory):
    """The momentum store as rotors spinning about axes in the tile's plane: kinetic energy and rotor mass for the
    stored momentum, and the mean power its round-trip losses draw if it fills and empties once an orbit."""
    h = attitude['schemes']['balanced_400_m']['store_N_m_s']['trim_0.25']['median']
    axes = attitude['schemes']['balanced_400_m']['required_orbit_mean_N_m']
    period = orbit_period_s()
    rows = []
    for r in DESIGN['store_rotor_radius_m']:
        for v in DESIGN['store_rim_speed_m_s']:
            energy, mass = .5*h*v/r, h/(v*r)
            loss = [(1-e)*energy/period for e in DESIGN['store_round_trip'][::-1]]
            rows.append(dict(rotor_radius_m=r, rim_speed_m_s=v, energy_GJ=energy/1e9, rotor_t=mass/1e3,
                             share_of_tile_mass=mass/inventory['tile_mass_kg'],
                             fleet_rotor_Gt=[mass*t/1e12 for t in inventory['tiles']],
                             tile_W=loss, fleet_GW=[loss[0]*inventory['tiles'][0]/1e9, loss[1]*inventory['tiles'][1]/1e9]))
    return dict(momentum_N_m_s=h, required_orbit_mean_N_m=axes, orbit_period_days=period/K.JULIAN_DAY, rotors=rows)


def placement(ledger):
    """The share of a terawatt released in orbit that the Moon intercepts, and the heat a terawatt adds to the Moon's
    global-mean budget when it is used there."""
    ratio = K.MOON_RADIUS/(DESIGN['orbit_radius_km']*1e3)
    intercepted = (1-math.sqrt(1-ratio**2))/2
    per_tw = 1e12/(4*math.pi*K.MOON_RADIUS**2)
    sunlight = ledger['requirements']['climate_top_of_air_W_m2']/4
    return dict(moon_share_of_isotropic_heat=intercepted, moon_W_m2_per_TW_used_on_moon=per_tw,
                moon_W_m2_per_TW_released_in_orbit=per_tw*intercepted, design_sunlight_global_mean_W_m2=sunlight,
                share_of_sunlight_per_TW_used_on_moon=per_tw/sunlight)


def main():
    ledger, attitude, films, held = (json.loads(p.read_text()) for p in (LEDGER, ATTITUDE, FILMS, HELD))
    inventory = fleet(ledger)
    out = dict(schema=SCHEMA,
               producer=dict(files={f: digest(ROOT/f) for f in FILES}, constants=constants_used(FILES)),
               inputs={str(p.relative_to(ROOT)): digest(p) for p in (LEDGER, ATTITUDE, FILMS, HELD, SOURCES)},
               evidence=('Closed-form screen on the named products and the assumptions in design, each sourced in '
                         'array_industry_sources.json; the films\' thermal emissivity comes from the protection '
                         'domain\'s optical tables (normal incidence, 2-125 um). No dynamics, device or plant model.'),
               reading_rule=('Powers in W, TW or PW as named; areas km2 of panel footprint where named panel, else of '
                             'tile or collector; temperatures K. Ranges pair the low ends and the high ends of the '
                             'assumptions. Tiles are taken facing the Sun at the sunward crossing for their temperature '
                             'and over the orbit for the fleet totals. The trim and store rows bound device classes and '
                             'rotor sizes for the planned branch to choose between.'),
               design=DESIGN, fleet=inventory, tiles=tiles(films, inventory), collectors=collectors(ledger),
               film_plant=film_plant(inventory), computing=computing(), trim=trim(inventory, films, held),
               store=store(attitude, inventory), placement=placement(ledger))
    out['checks'] = dict(
        collector_area_matches_ledger=all(math.isclose(r['collector_km2']*1e6*r['conversion']*K.SOLAR_CONSTANT,
                                                       r['demand_TW']*1e12, rel_tol=1e-9) for r in out['collectors']),
        rotor_energy_is_half_m_v2=all(math.isclose(r['energy_GJ']*1e9, .5*r['rotor_t']*1e3*r['rim_speed_m_s']**2,
                                                   rel_tol=1e-12) for r in out['store']['rotors']),
        tiles_colder_than_collectors=max(f['temperature_K'] for f in out['tiles']['films'].values())
        < min(r['temperature_K'] for r in out['collectors']))
    OUT.write_text(json.dumps(out, indent=1)+'\n')
    print(json.dumps(dict(fleet=inventory, tiles={k: {kk: round(vv, 6) if isinstance(vv, float) else vv
                                                       for kk, vv in v.items()} for k, v in out['tiles']['films'].items()},
                          absorbed_PW=out['tiles']['absorbed_PW'], checks=out['checks']), indent=1))


if __name__ == '__main__':
    main()
