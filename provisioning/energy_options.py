"""The Moon's energy supply on a common basis: beamed power from the ring fleet against its alternatives.

    python -m provisioning.energy_options
    # -> provisioning/results/energy_options.json, in under a minute

Every option is set against the loads the domain synthesis's array case sizes (C' with the food blend, 8 billion on
the Moon) and judged per terawatt delivered on the same measures:
- the heat each delivered watt leaves on the Moon beyond the sunlight it would have absorbed anyway, as W/m2 and as
  points of the shield's dimming (a point is 1% of the hot world's absorbed sunlight, as the heat-levers study and
  the climate targets count it), with what a point costs in rain and growth;
- the coupling between that dimming and the tiles' power: the light a point of dimming intercepts and the electricity
  collectors make from it;
- land and sea, storage through the 354-hour night, mass and materials, reliability and maturity.

The air's absorption is the power screen's ITU-R P.676-13 model (provisioning/power.py), here run through the moon
column at 295 and 301 K, through a 1.0 atm column for comparison, and from receivers raised to 40 km. Sunlight at the
ground, cloud and the low winds come from the design run's own output (provisioning/energy_extract.py ->
results/energy_climate.json). Everything else is arithmetic on committed products and on published, recalled or
estimated values, each labelled in the product with its kind and source.

Closed-form arithmetic and one-dimensional radiative transfer; nothing here designs or simulates a plant, a beam, a
store or a grid.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

from atmosphere.radiative_convective import thermodynamics as th
from provisioning import power as P
from shared.constants import (AVOGADRO, ELEMENTARY_CHARGE, GAS_CONSTANT, JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_GM,
                              MOON_RADIUS, MOON_SURFACE_GRAVITY, STANDARD_ATMOSPHERE, SYNODIC_MONTH_DAYS)
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SCHEMA = 'terluna.provisioning.energy-options/1'
OUT = HERE / 'results' / 'energy_options.json'
SOURCES = HERE / 'energy_options_sources.json'
CODE = ('provisioning/energy_options.py', 'provisioning/power.py', 'atmosphere/radiative_convective/thermodynamics.py')

STUDIES = 'research/studies'
CLIMATE = 'provisioning/results/energy_climate.json'
POWER = 'provisioning/results/power.json'
SCREEN = 'provisioning/results/first_screen.json'
LEDGERS = f'{STUDIES}/domain_synthesis/results/ledgers.json'
LEVERS = f'{STUDIES}/heat_levers/results/heat_levers.json'
WATER = 'shared/scenarios/water.json'
RINGS = P.RINGS
INPUTS = (CLIMATE, POWER, SCREEN, LEDGERS, LEVERS, WATER) + tuple(RINGS.values())

TW = 1e12
HOURS_PER_YEAR = JULIAN_YEAR_DAYS * 24.0
MONTH_H = SYNODIC_MONTH_DAYS * 24.0
MOON_AREA = 4.0 * math.pi * MOON_RADIUS ** 2

# ---------------------------------------------------------------------------------------------------------------
# Values this screen takes from outside the products, each with its kind: published (source and how far read in
# energy_options_sources.json), recalled (from general knowledge, not checked this session), estimate (a stated
# approximation), design guess (a value chosen for the comparison).

V = dict(
    collector_conversion=dict(value=(0.20, 0.30), kind='project',
                              basis="the climate targets' and heat-levers study's 20-30% for collectors on dimmed light"),
    rectenna_peak_w_m2=dict(value=230.0, kind='project', basis="the 1978 reference's peak, as the synthesis fills the "
                                                                "ground's patches"),
    receiver_buffer=dict(value=1.25, kind='project', basis="the power screen's ground closed round each patch to the "
                                                           "public 10 W/m2, a quarter more"),
    rectenna_areal_kg_m2=dict(value=(5.0, 20.0), kind='design guess',
                              basis='a ground rectenna of mesh, dipoles and footings; the 1978 reference used steel '
                                    'panels on continuous concrete footings, not quantified here'),
    sps_kg_per_kw=dict(value=10.0, kind='published',
                       basis="the 1978 reference satellite, about 5e4 t for 5 GW delivered (OTA 1981 summary, from a "
                             "search snippet; not read in full); thin-film tiles would be lighter"),
    line_raised=dict(value=(0.97, 0.98), kind='estimate',
                     basis='a raised receiver: conductor down its tether and a converter, 2-3% in all'),
    raised_heights_km=dict(value=(20.0, 30.0, 40.0), kind='design guess', basis='above the storms (tops 22 km, the '
                                                                                '0 C level 27-28 km) and the reefs'),
    pv_module=dict(value=(0.22, 0.30), kind='published/estimate',
                   basis="22%: commercial wafer modules from 2024 production (Fraunhofer ISE Photovoltaics Report "
                         "2025); 30%: an estimate for commercial tandem modules"),
    pv_derate=dict(value=0.85, kind='estimate', basis='temperature, inverter, wiring and soiling together'),
    pv_tilt_gain=dict(value=(1.0, 1.15), kind='estimate',
                      basis="tilting or tracking against a flat panel, small in the deep air's diffuse light"),
    ground_cover_ratio=dict(value=(0.4, 0.7), kind='estimate', basis='panel area over the land a plant occupies'),
    pv_albedo=dict(value=(0.05, 0.10), kind='estimate', basis='glass-faced modules with their frames and gaps'),
    cloud_reflect_up=dict(value=(0.3, 0.6), kind='estimate',
                          basis="the share of light reflected up from the ground that a cloud above sends back down"),
    surface_kernel_clear=dict(value=0.34, kind='project',
                              basis="the heat-levers study's clear-sky kernel at the equator: the Moon reflects 0.34 "
                                    "more of the arriving sunlight for each unit of surface albedo"),
    pv_t_per_mwp=dict(value=(100.0, 130.0), kind='published/estimate',
                      basis='modules about 55 t/MWp (Mirletz et al. 2022: 97 Mt for 1.75 TW), mounting steel 60-67 '
                            't/MWp and aluminium 6-8 t/MWp (search summaries); Earth mounting, lighter at a sixth g'),
    hydrogen_round_trip=dict(value=0.40, kind='project', basis='Headley & Schoenung 2020, in the power screen'),
    fuel_cell=dict(value=0.5, kind='project', basis="the power screen's electricity per unit of hydrogen's heating "
                                                    "value"),
    hydrogen_kwh_kg=dict(value=33.3, kind='published', basis="hydrogen's lower heating value"),
    hydrogen_store_density=dict(value=15.0, kind='recalled', basis='kg/m3 of hydrogen near 20 MPa and 300 K'),
    pumped_round_trip=dict(value=(0.75, 0.80), kind='recalled', basis="Earth's pumped hydro"),
    pumped_head_m=dict(value=3740.0, kind='project', basis="Korolev's head, the first screen"),
    battery_round_trip=dict(value=0.88, kind='recalled', basis='lithium-ion packs'),
    battery_kwh_kg=dict(value=(0.15, 0.3), kind='project', basis="the power screen's pack specific energy"),
    iron_air_round_trip=dict(value=(0.40, 0.45), kind='published', basis='Form Energy target for its 100-hour '
                                                                         'iron-air system (secondary sources)'),
    grid_long_loss=dict(value=(0.09, 0.19), kind='project',
                        basis='the first screen: HVDC 2-3.5% per 1,000 km and 1% a converter over 4,000-5,500 km '
                              '(US EIA 2018)'),
    grid_local_loss=dict(value=0.02, kind='estimate', basis='local lines and one converter'),
    grid_line_gw=dict(value=12.0, kind='published', basis='Changji-Guquan, 1,100 kV, 12 GW over more than 3,000 km '
                                                         '(Hitachi Energy; in sources.json)'),
    grid_distance_km=dict(value=(4000.0, 5500.0), kind='project', basis="the first screen's distances to the day side"),
    grid_al_t_per_km=dict(value=30.0, kind='estimate',
                          basis='aluminium for both poles of a 12 GW line at 1,100 kV losing 2.5% per 1,000 km'),
    roof_usable_share=dict(value=(0.3, 0.5), kind='estimate', basis='share of roofs and pavements that carries panels'),
    grid_peak_transfer=dict(value=(0.5, 0.7), kind='estimate',
                            basis="the share of the ground's load the lines carry at once, the night side's loads "
                                  "with C''s far-side lean"),
    grid_row_km=dict(value=(0.05, 0.1), kind='estimate', basis='right of way a line'),
    aerial_transmission=dict(value=(0.80, 0.90), kind='estimate',
                             basis='sunlight reaching panels 30 km up, above the cloud and most of the vapour'),
    aerial_albedo_below=dict(value=0.40, kind='project', basis="the design run's planetary albedo (high-cloud "
                                                               "product), what panels above the cloud cover"),
    aerial_escape=dict(value=0.9, kind='estimate', basis='share of light reflected at 30 km that leaves the air'),
    laser_wall_plug=dict(value=(0.5, 0.7), kind='published',
                         basis='fibre lasers above 50% (IPG) and diode lasers to 71% at 780 nm (2023)'),
    laser_capture=dict(value=0.9, kind='estimate', basis="beam spread and the dense air's turbulence"),
    laser_cell=dict(value=(0.40, 0.65), kind='published/estimate',
                    basis='68.9% for a GaAs cell under 858 nm laser light (Fraunhofer ISE 2021), about 20% '
                          'laser-to-electric at the receiver in the 2025 DARPA POWER record; arrays at high flux '
                          'between them'),
    rayleigh_earth_850nm=dict(value=0.017, kind='recalled', basis="Earth's Rayleigh optical depth at 850 nm"),
    rayleigh_earth_550nm=dict(value=0.097, kind='recalled', basis="Earth's Rayleigh optical depth at 550 nm"),
    optical_other_depth=dict(value=(0.05, 0.30), kind='estimate',
                             basis="water vapour's window absorption and aerosol in the deep, moist column"),
    optical_escape=dict(value=(0.0, 0.5), kind='estimate', basis='share of the light lost from the beam that leaves '
                                                                  'the Moon instead of warming it'),
    rx_5_8=dict(value=0.85, kind='estimate', basis='the 2.45 GHz chain kept at 5.8 GHz; 82% is published for a 5.8 '
                                                   'GHz rectenna, 90.6% for Brown 1977 at 2.45 GHz'),
    thermal_efficiency_fission=dict(value=(0.33, 0.45), kind='recalled', basis='light-water to advanced reactors'),
    heat_per_w_fusion=dict(value=(2.0, 3.0), kind='project', basis="the domain synthesis's fusion plant"),
    mev_per_fission=dict(value=200.0, kind='recalled', basis='energy a fission of U-233'),
    thorium_ppm=dict(value=(4.0, 12.0), kind='published', basis='Procellarum KREEP Terrane, Lunar Prospector '
                                                               '(Lawrence et al.; search summaries)'),
    mev_dt=dict(value=17.6, kind='recalled', basis='energy a D-T fusion'),
    heat_flow_mw_m2=dict(value=(14.0, 21.0), kind='published',
                         basis='Apollo 15 and 17 (Langseth et al. 1976); Siegler & Smrekar 2014 put the internal '
                               'flow at 9-13 mW/m2'),
    geothermal_efficiency=dict(value=(0.10, 0.15), kind='recalled', basis='binary plants on low-grade heat'),
    otec_w_m2=dict(value=(0.023, 0.047), kind='published/estimate',
                   basis="Earth's 7-14 TW (Rajagopalan & Nihous 2013) over about 3e14 m2 of warm ocean, carried to "
                         "the Moon's seas by area alone; their deep water's temperature is open"),
    hydro_developable=dict(value=(0.3, 0.6), kind='estimate', basis="share of the rivers' gross power a fleet of "
                                                                    "plants could take"),
    turbine=dict(value=0.85, kind='recalled', basis='hydro turbine and generator'),
    earth_wind_realized_w_m2=dict(value=0.5, kind='published', basis='US wind farms in 2016 (Miller & Keith 2018)'),
    earth_land_wind_power_density_w_m2=dict(value=(323.0, 410.0), kind='published',
                                            basis="the Global Wind Atlas's mean power density over its whole area "
                                                  "and over the windiest half (Badger et al.; search summary)"),
    moon_gust_factor=dict(value=(1.0, 3.0), kind='estimate',
                          basis="instantaneous over 3-day-mean wind cubes; the design run's means lack the gusts"),
    food_kcal=dict(value=2500.0, kind='project', basis="the food studies' diet"),
    food_wet_kg=dict(value=(0.55, 2.0), kind='estimate', basis='a person-day shipped dry to fresh'),
    phosphorus_kg_person_year=dict(value=0.44, kind='project', basis="README finding 11: the metropolis passes 44 kt "
                                                                     "of P a year for 100 million"),
    computing_kept=dict(value=0.25, kind='project', basis="the power screen's millisecond control kept on the Moon"),
    aloft_share=dict(value=0.5, kind='project', basis="C': 4 of the Moon's 8 billion aloft"),
)
FREQS_GHZ = (1.0, 2.45, 5.8, 10.0)
HEIGHTS_KM = (0.0, 4.0, 10.0, 20.0, 30.0, 40.0)
COLUMNS = (('1.2 atm, 295 K', 121590.0, 295.0), ('1.2 atm, 301 K', 121590.0, 301.0),
           ('1.0 atm, 301 K', STANDARD_ATMOSPHERE, 301.0))
CRITICAL_DAYS = 7.0

EVIDENCE = ("Arithmetic on committed products, one-dimensional gaseous absorption (ITU-R P.676-13 Annex 1, the power "
            "screen's model, reproduced against ITU-R Study Group 3's Earth values) through prescribed clear-sky "
            "columns, and means of the design run's own sunlight, cloud and winds. Published values are labelled with "
            "their sources and how far each was read; recalled values, estimates and design guesses are labelled as "
            "such. No plant, store, beam, grid or line is designed or simulated, and nothing here is validated against "
            "a built system.")
READING_RULE = ("Pairs are [low-cost end, high-cost end] of C' with the food blend, as the domain synthesis gives them; "
                "where an option's own range is wide its low end goes with the low-cost loads and its high end with "
                "the high-cost loads, so a pair spans both. 'Per TW' means per terawatt delivered to the loads. Heat "
                "is what the option leaves on the Moon beyond the sunlight it would have absorbed anyway, at an "
                "efficacy of 1 as the climate targets take it; a point of dimming is 1% of the hot world's absorbed "
                "sunlight, 225.4 W/m2. Land in km2; storage in TWh delivered; masses in Mt or Gt.")


def sha256(path) -> str:
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def load(path):
    return json.loads((ROOT / path).read_text())


sig = P.sig


def val(key):
    return V[key]['value']


def lohi(x):
    return (float(x[0]), float(x[1])) if isinstance(x, (tuple, list)) else (float(x), float(x))


# ---------------------------------------------------------------------------------------------------------------
# The loads and the basis.

def loads(ledgers, power):
    """Delivered TW by component: the ledgers' array case is intake at the Earth-service heat per watt."""
    array = ledgers['heat']['array']
    factor = power['absorption']['link']['heat_per_w_delivered_without_air']
    parts = {k: [v[0] / factor, v[1] / factor] for k, v in array['heat_by_component_tw'].items()}
    parts['close_computing'] = [v / factor for v in array['close_computing_heat_tw']]
    aloft = val('aloft_share')
    sky_keys = ('sky_control', 'sky_hydrogen_remade', 'sky_trips_and_food_lift')
    sky = [sum(parts[k][e] for k in sky_keys) + aloft * parts['residents'][e] for e in (0, 1)]
    ground = [parts['food'][e] + parts['water'][e] + (1 - aloft) * parts['residents'][e]
              + parts['close_computing'][e] for e in (0, 1)]
    total = [sky[e] + ground[e] for e in (0, 1)]
    return dict(delivered_tw=parts, sky_tw=sky, ground_tw=ground, total_tw=total, intake_factor=factor,
                rule="the synthesis's array heat by component over the Earth-service 1.376 W a delivered watt; the "
                     "towns take their control, hydrogen, trips and the aloft half of the residents")


def basis(levers, power):
    absorbed = levers['basis']['absorbed_sunlight_w_m2']
    point_tw = 0.01 * absorbed * MOON_AREA / TW
    light = levers['shield']['light_tw_per_point']
    conv = lohi(val('collector_conversion'))
    return dict(
        moon_area_m2=MOON_AREA, moon_w_m2_per_tw=TW / MOON_AREA,
        point=dict(absorbed_sunlight_w_m2=absorbed, w_m2=0.01 * absorbed, heat_tw=point_tw,
                   light_at_the_shield_tw=light, electricity_tw=[light * conv[0], light * conv[1]],
                   growth_percent=levers['shield']['growth_through_light_percent_per_point'],
                   land_rain_percent=levers['shield']['land_rain_percent_per_point']),
        absorbing_dimmer_glow_return=levers['shield']['absorbing_dimmer_glow_return'],
        resort_dimming_points=levers['basis']['targets_dimming_in_all_percent']['resort, 301 K'],
        resort_dimming_light_tw=[light * p for p in levers['basis']['targets_dimming_in_all_percent']['resort, 301 K']],
        resort_dimming_electricity_tw=[light * levers['basis']['targets_dimming_in_all_percent']['resort, 301 K'][0]
                                       * conv[0],
                                       light * levers['basis']['targets_dimming_in_all_percent']['resort, 301 K'][1]
                                       * conv[1]],
        resort_dimming_electricity_at_the_lower_conversion_tw=[
            light * p * conv[0] for p in levers['basis']['targets_dimming_in_all_percent']['resort, 301 K']],
        beams_now_points=levers['beams']['all_beams_points'],
        tile_electricity_tw=power['synthesis_check']['dimmer_against_tile_electricity']['needed_tw'],
        planned_dimmer_electricity_tw=power['synthesis_check']['dimmer_against_tile_electricity']['dimmer_tw'])


def per_tw(h, b):
    """What a delivered TW costs at h W of heat a delivered watt: points, W/m2, rain and growth."""
    pts = h / b['point']['heat_tw']
    ret = b['absorbing_dimmer_glow_return']
    return dict(heat_w_per_w=h, w_m2=h * b['moon_w_m2_per_tw'], points=pts,
                points_if_the_dimmer_absorbs=[pts / (1 - ret['rolled']), pts / (1 - ret['as_built'])],
                land_rain_percent=[pts * b['point']['land_rain_percent'][0], pts * b['point']['land_rain_percent'][1]],
                growth_percent=[pts * b['point']['growth_percent'][0], pts * b['point']['growth_percent'][1]])


# ---------------------------------------------------------------------------------------------------------------
# The air: absorption by frequency, receiver height and column.

def column(dry_pa, ts_k):
    air = th.earthlike_air(dry_pa, P.COLUMN['co2_ppm'])
    spec = th.ColumnSpec(th.MOON, air, ts_k, P.COLUMN['stratosphere_k'], th.manabe_wetherald(),
                         troposphere_levels=P.COLUMN['levels'][0], upper_levels=P.COLUMN['levels'][1])
    col = th.build_column(spec)
    z = col['z_m'] / 1e3
    lnp, lnx = np.log(col['p_pa']), np.log(np.maximum(col['x_h2o'], 1e-30))

    def profile(h_km):
        p = np.exp(np.interp(h_km, z, lnp)) / 100.0
        e = np.exp(np.interp(h_km, z, lnx)) * p
        return np.interp(h_km, z, col['t_k']), p - e, e
    return dict(profile=profile, air=air, top_km=float(z[-1]), surface_pressure_pa=float(col['surface_pressure_pa']),
                air_column_t_m2=float(col['layer_mass_kg_m2'].sum() / 1e3))


def air_table(elevations):
    out = {}
    for name, dry, ts in COLUMNS:
        col = column(dry, ts)
        o2 = col['air'].fractions['O2'] / P.EARTH_O2
        n2 = (col['air'].fractions['N2'] / P.EARTH_N2) ** 2
        rows = {}
        for f in FREQS_GHZ:
            for h in HEIGHTS_KM:
                edges = P.layer_edges(h, col['top_km'])
                for e in elevations:
                    dry_db, wet_db = P.slant_path(f, e, col['profile'], MOON_RADIUS / 1e3, edges, o2, n2)
                    rows[f'{f:g}_ghz_{h:g}_km_{e:g}_deg'] = dict(db=dry_db + wet_db, vapour_share=wet_db / (dry_db + wet_db),
                                                                 share_lost=float(P.loss(dry_db + wet_db)))
        out[name] = dict(surface_pressure_pa=col['surface_pressure_pa'], air_column_t_m2=col['air_column_t_m2'],
                         o2_mole_fraction=col['air'].fractions['O2'], rows=rows)
    return out


# ---------------------------------------------------------------------------------------------------------------
# Beams.

def beam_option(name, share_lost, chain, b, rx=None, line=1.0, receiver_land=True, notes=''):
    """One beam variant per delivered watt: its heat on the Moon by part, the tiles' electricity, the dimming its heat
    takes, the light and electricity that dimming gives the fleet back, and the closed ground of its receivers."""
    cap, rxe, bus, tx = chain['capture'], rx or chain['rx'], chain['bus'], chain['tx']
    h = 1.0 / (cap * (1 - share_lost) * rxe * bus * line)
    parts = dict(delivered=1.0, bus=1.0 / bus - 1.0, rectenna=(1.0 / rxe - 1.0) / bus,
                 line=(1.0 / line - 1.0) / (rxe * bus) if line < 1 else 0.0,
                 spill=(1.0 / cap - 1.0) / (rxe * bus * line),
                 air=share_lost / (1 - share_lost) / (cap * rxe * bus * line))
    tile = h / tx
    cost = per_tw(h, b)
    light = cost['points'] * b['point']['light_at_the_shield_tw']
    elec = [light * c for c in lohi(val('collector_conversion'))]
    intake = 1.0 / (cap * rxe * bus * line)                  # arriving at the receiver plane, spill included
    land = intake / val('rectenna_peak_w_m2') * TW / 1e6 * val('receiver_buffer') if receiver_land else 0.0
    return dict(name=name, share_lost=share_lost, heat_parts_w_per_w=parts, tile_electricity_w_per_w=tile,
                link_dc_efficiency=1.0 / tile, **cost,
                dimming_light_tw_per_tw=light, dimming_electricity_tw_per_tw=elec,
                self_supply_share=[e / tile for e in elec], ground_closed_km2_per_tw=land, notes=notes)


def pairwise(name, shares, chain, b, **kw):
    """An option over its own range: [low heat, high heat] from the smaller and larger share lost."""
    lines = lohi(kw.pop('line', 1.0))
    return [beam_option(name, shares[0], chain, b, line=lines[1], **kw),
            beam_option(name, shares[1], chain, b, line=lines[0], **kw)]


def beams(power, table, b):
    """The chosen system as the power screen gives it, then every variant for the ground's beams, whose loads the
    alternatives compete for; the towns keep 2.45 GHz to their crowns in every variant."""
    chain = power['absorption']['link']['chain']
    lw = power['absorption']['load_weighted']
    link = power['absorption']['link']
    e_g = power['reserves']['storms']['fades']['p99_rate']['ground']['operating']['elevation_deg']
    e_t = power['reserves']['storms']['fades']['p99_rate']['towns']['operating']['elevation_deg']
    e_op = round(e_g, 2)
    base, now = '1.2 atm, 295 K', '1.2 atm, 301 K'
    db = (lw['2.45_ghz']['ground_best_tile_db'], lw['2.45_ghz']['ground_db'])      # highest tiles, any usable tile
    db58 = (lw['5.8_ghz']['ground_best_tile_db'], lw['5.8_ghz']['ground_db'])
    rows = table[base]['rows']

    def carried(dbs, f, h, col=base, ref_f=2.45):
        ratio = (table[col]['rows'][f'{f:g}_ghz_{h:g}_km_{e_op:g}_deg']['db']
                 / rows[f'{ref_f:g}_ghz_0_km_{e_op:g}_deg']['db'])
        return [float(P.loss(d * ratio)) for d in dbs]
    opts = dict(
        chosen_system=[beam_option('2.45 GHz as chosen, towns and ground together', s, chain, b,
                                   notes='the power screen: beams from the highest tiles and from any above 30 degrees')
                       for s in (link['with_the_best_tiles']['share_of_beams_lost'][0],
                                 link['share_of_beams_lost'][0])],
        towns_2_45_ghz=pairwise('2.45 GHz to the towns at 10 km', [float(P.loss(lw['2.45_ghz']['towns_best_tile_db'])),
                                                                   float(P.loss(lw['2.45_ghz']['towns_db']))], chain,
                                b, receiver_land=False, notes="on the towns' own crowns"),
        ground_2_45_ghz=pairwise('2.45 GHz to the ground, as chosen', [float(P.loss(d)) for d in db], chain, b),
        ground_no_air=pairwise('the same link with no air', (0.0, 0.0), chain, b,
                               notes='the floor any beam meets: spill, rectennas and bus'),
        ground_2_45_ghz_at_301_k=pairwise('2.45 GHz to the ground through the design climate\'s 301 K column',
                                          carried(db, 2.45, 0.0, now), chain, b),
        ground_2_45_ghz_through_1_atm=pairwise('2.45 GHz to the ground through 1.0 atm at 301 K (for comparison; the '
                                               'air is a design condition)', carried(db, 2.45, 0.0, '1.0 atm, 301 K'),
                                               chain, b),
        ground_1_ghz=pairwise('1 GHz to the ground', carried(db, 1.0, 0.0), chain, b,
                              notes='apertures 2.45 times as wide for the same receiver'),
        ground_5_8_ghz=pairwise('5.8 GHz to the ground', [float(P.loss(d)) for d in db58], chain, b,
                                rx=val('rx_5_8'), notes='storms stop it: the 90th-percentile storm takes 68-97%'))
    for h in val('raised_heights_km'):
        opts[f'raised_{h:g}_km_2_45_ghz'] = pairwise(
            f'2.45 GHz to receivers raised to {h:g} km, lines down', carried(db, 2.45, h), chain, b,
            line=val('line_raised'), receiver_land=False,
            notes='the ground beneath stays open; tethers stand in the flight bands')
    opts['raised_30_km_5_8_ghz'] = pairwise('5.8 GHz to receivers raised to 30 km, lines down',
                                            carried(db58, 5.8, 30.0, ref_f=5.8), chain, b, rx=val('rx_5_8'),
                                            line=val('line_raised'), receiver_land=False,
                                            notes='above the storms, so 5.8 GHz holds through them')
    return opts, dict(ground_operating_elevation_deg=e_g, towns_operating_elevation_deg=e_t,
                      rule='variants carry the power screen\'s load-weighted ground loss by the ratio of losses at the '
                           'ground beams\' mean operating elevation (an estimate); pairs are [highest tiles, any tile '
                           'above 30 degrees] and, for raised receivers, the better and worse line')


# ---------------------------------------------------------------------------------------------------------------
# Optical links.

def optical(power, b, climate):
    col = power['absorption']['column']['air_column_t_m2']
    earth_col = STANDARD_ATMOSPHERE / 9.80665 / 1e3
    tau_r = val('rayleigh_earth_850nm') * col / earth_col
    other = lohi(val('optical_other_depth'))
    e = math.radians(61.23)
    t = [math.exp(-(tau_r + o) / math.sin(e)) for o in other]                 # clear-sky transmission, both ends
    wp, cell, esc = lohi(val('laser_wall_plug')), lohi(val('laser_cell')), lohi(val('optical_escape'))
    cap, bus = val('laser_capture'), power['absorption']['link']['chain']['bus']
    # Low heat: clear air, the best cell, half the lost light leaving; high heat the reverse.
    h_lo = (1 - esc[1] * (1 - t[0])) / (cap * t[0] * cell[1] * bus)
    h_hi = (1 - esc[0] * (1 - t[1])) / (cap * t[1] * cell[0] * bus)
    tile = [1 / (wp[1] * cap * t[0] * cell[1] * bus), 1 / (wp[0] * cap * t[1] * cell[0] * bus)]
    cloud = climate['whole']['land']['cloud']
    tau_550 = val('rayleigh_earth_550nm') * col / earth_col
    return dict(rayleigh_depth_850nm=tau_r, rayleigh_depth_550nm=tau_550,
                direct_sunbeam_rayleigh_only_overhead=math.exp(-tau_550),
                clear_transmission_at_61_deg=t,
                heat_w_per_w=[h_lo, h_hi], tile_electricity_w_per_w=tile,
                points=[h_lo / b['point']['heat_tw'], h_hi / b['point']['heat_tw']],
                share_of_time_cloud_blocks_a_ground_receiver=cloud,
                notes="cloud stops the beam: the design run covers the land's sky 36% of the time on average")


# ---------------------------------------------------------------------------------------------------------------
# Local solar: sunlight at the ground, the night, land, albedo and stores.

def solve_store(shape, eta):
    """For a constant unit load fed by c*shape (mean 1) and a store of round trip eta: the c that balances, the share
    the store carries and its energy per unit load over the cycle (in units of the cycle's mean load-hours)."""
    s = np.asarray(shape, float) / np.mean(shape)

    def residual(c):
        return np.mean(np.minimum(c * s, 1.0)) + eta * np.mean(np.maximum(c * s - 1.0, 0.0)) - 1.0
    lo, hi = 1.0, 1.0 / eta * 1.5 + 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if residual(mid) < 0:
            lo = mid
        else:
            hi = mid
    c = 0.5 * (lo + hi)
    deficit = np.maximum(1.0 - c * s, 0.0)
    # The store's size: the largest run of deficit, the cycle wrapped once.
    run, best = 0.0, 0.0
    for x in np.concatenate([deficit, deficit]):
        run = run + x if x > 0 else 0.0
        best = max(best, run)
    best = min(best, deficit.sum())
    return dict(generation_over_load=c, store_share=float(np.mean(deficit)),
                store_cycle_fraction=float(best / s.size))


def load_weights(power, bands):
    """Each band's share of the ground's load, the power screen's groups spread over their latitudes by area."""
    groups = power['grid']['groups']
    w = {k: [0.0, 0.0] for k in bands}
    for g in groups.values():
        lo, hi = g['lat_deg']
        for k, band in bands.items():
            a, c = max(lo, band['lat_deg'][0]), min(hi, band['lat_deg'][1])
            if c <= a:
                continue
            share = ((math.sin(math.radians(c)) - math.sin(math.radians(a)))
                     / (math.sin(math.radians(hi)) - math.sin(math.radians(lo))))
            for e in (0, 1):
                w[k][e] += share * g['intake_tw'][e]
    tot = [sum(v[e] for v in w.values()) for e in (0, 1)]
    return {k: [v[0] / tot[0], v[1] / tot[1]] for k, v in w.items()}


def local_solar(climate, power, b, beams_h, screen):
    bands = climate['bands']
    weights = load_weights(power, bands)
    whole = climate['whole']
    module, derate, tilt = lohi(val('pv_module')), val('pv_derate'), lohi(val('pv_tilt_gain'))
    gcr = lohi(val('ground_cover_ratio'))
    land_albedo = climate['run']['land_albedo']
    pv_alb = lohi(val('pv_albedo'))
    d_alb = (land_albedo - pv_alb[1], land_albedo - pv_alb[0])
    band0 = bands['0']
    t_clear = band0['land_ghi_clear_w_m2'] / band0['land_toa_w_m2']
    t_up_clear = val('surface_kernel_clear') / t_clear
    reflect = lohi(val('cloud_reflect_up'))
    cloud = whole['land']['cloud']
    t_up = (t_up_clear * (1 - cloud * reflect[1]), t_up_clear * (1 - cloud * reflect[0]))
    stores = dict(hydrogen=lohi(val('hydrogen_round_trip')), pumped_water=lohi(val('pumped_round_trip')),
                  batteries=lohi(val('battery_round_trip')), iron_air=lohi(val('iron_air_round_trip')))
    by_band = {}
    for k, band in bands.items():
        shape = band['day_profile_ghi_w_m2']
        rows = {}
        for name, (r_lo, r_hi) in stores.items():
            a, c = solve_store(shape, r_hi), solve_store(shape, r_lo)
            rows[name] = dict(generation_over_load=[a['generation_over_load'], c['generation_over_load']],
                              store_share=[a['store_share'], c['store_share']],
                              store_twh_per_tw=[a['store_cycle_fraction'] * MONTH_H,
                                                c['store_cycle_fraction'] * MONTH_H])
        s = np.asarray(shape, float) / np.mean(shape)
        remote = float(np.mean(np.maximum(1.0 - s, 0.0)))       # the share a grid brings from elsewhere, no store
        by_band[k] = dict(lat_deg=band['lat_deg'], ghi_w_m2=band['land_ghi_w_m2'],
                          ghi_clear_w_m2=band['land_ghi_clear_w_m2'], toa_w_m2=band['land_toa_w_m2'],
                          ground_over_top=band['land_ghi_w_m2'] / band['land_toa_w_m2'],
                          load_weight=weights[k], stores=rows, grid_remote_share=remote)
    # The ground's loads keep C''s composition at both ends, so one weighting serves both.
    ghi = sum(by_band[k]['ghi_w_m2'] * weights[k][0] for k in bands)
    remote = sum(by_band[k]['grid_remote_share'] * weights[k][0] for k in bands)
    long_loss, local_loss = lohi(val('grid_long_loss')), val('grid_local_loss')
    t_per = lohi(val('pv_t_per_mwp'))
    routes = {}
    # Per delivered TW: panels, land, heat, light; [better, worse] with the better or worse module, tilt, spacing and
    # store or line loss.
    for name in ('global_grid', 'hydrogen', 'pumped_water', 'batteries', 'iron_air'):
        if name == 'global_grid':
            f = [1 - remote * long_loss[0] - (1 - remote) * local_loss,
                 1 - remote * long_loss[1] - (1 - remote) * local_loss]
            store = [0.0, 0.0]
        else:
            f = [1.0 / sum(by_band[k]['stores'][name]['generation_over_load'][i] * weights[k][0] for k in bands)
                 for i in (0, 1)]
            store = [sum(by_band[k]['stores'][name]['store_twh_per_tw'][i] * weights[k][0] for k in bands)
                     for i in (0, 1)]
        gain = [module[1] * derate * tilt[1] * f[0], module[0] * derate * tilt[0] * f[1]]   # delivered per W of light
        panels = [TW / (ghi * g) / 1e6 for g in gain]
        peak = [panels[0] * 1e6 * module[1] * 1000 / TW, panels[1] * 1e6 * module[0] * 1000 / TW]
        routes[name] = dict(delivered_w_m2_of_panel=[ghi * g for g in gain], delivered_over_generated=f,
                            panel_km2_per_tw=panels, land_km2_per_tw=[panels[0] / gcr[1], panels[1] / gcr[0]],
                            store_twh_per_tw=store,
                            heat_w_per_w=[d_alb[0] * t_up[0] / gain[0], d_alb[1] * t_up[1] / gain[1]],
                            surface_light_w_per_w=[1 / g for g in gain],
                            peak_capacity_tw_per_tw=peak,
                            pv_mass_gt_per_tw=[peak[0] * 1e6 * t_per[0] / 1e9, peak[1] * 1e6 * t_per[1] / 1e9])
    # The grid round the Moon: lines from the day side carrying the share the night and low sun leave.
    transfer = lohi(val('grid_peak_transfer'))
    dist = lohi(val('grid_distance_km'))
    lines_km = [transfer[0] * TW / (val('grid_line_gw') * 1e9) * dist[0],
                transfer[1] * TW / (val('grid_line_gw') * 1e9) * dist[1]]
    row = lohi(val('grid_row_km'))
    grid = dict(remote_share=remote, peak_transfer_tw_per_tw=list(transfer), line_km_per_tw=lines_km,
                aluminium_mt_per_tw=[x * val('grid_al_t_per_km') / 1e6 for x in lines_km],
                corridor_km2_per_tw=[lines_km[0] * row[0], lines_km[1] * row[1]])
    # The store's own scale per delivered TW (hydrogen, water, batteries).
    h2 = routes['hydrogen']['store_twh_per_tw']
    fc, lhv = val('fuel_cell'), val('hydrogen_kwh_kg')
    head = val('pumped_head_m')
    pw = routes['pumped_water']['store_twh_per_tw']
    bt = routes['batteries']['store_twh_per_tw']
    kwh_kg = lohi(val('battery_kwh_kg'))
    stores_scale = dict(
        hydrogen_mt_per_tw=[x * 1e9 / (fc * lhv) / 1e9 for x in h2],
        hydrogen_volume_km3_per_tw_at_20_mpa=[x * 1e9 / (fc * lhv) / val('hydrogen_store_density') / 1e9 for x in h2],
        pumped_water_km3_per_tw=[x * 1e12 * 3600 / (1000.0 * MOON_SURFACE_GRAVITY * head * val('turbine')) / 1e9
                                 for x in pw],
        battery_gt_per_tw=[bt[0] * 1e9 / kwh_kg[1] / 1e12, bt[1] * 1e9 / kwh_kg[0] / 1e12])
    # Beams' dimming against panels: the sunlight each takes from the ground per delivered watt.
    ground_over_top = whole['ghi_w_m2'] / whole['toa_w_m2']
    beam_light = [h * b['point']['light_at_the_shield_tw'] / b['point']['heat_tw'] * ground_over_top for h in beams_h]
    # Panels above the cloud: they cover a bright planet and darken it.
    t30, alb_below, esc = lohi(val('aerial_transmission')), val('aerial_albedo_below'), val('aerial_escape')
    f_h2 = routes['hydrogen']['delivered_over_generated']
    f_grid = routes['global_grid']['delivered_over_generated']
    aerial = dict(hydrogen=[(alb_below - pv_alb[1]) * esc / (module[1] * derate * f_h2[0]),
                            (alb_below - pv_alb[0]) * esc / (module[0] * derate * f_h2[1])],
                  global_grid=[(alb_below - pv_alb[1]) * esc / (module[1] * derate * f_grid[0]),
                               (alb_below - pv_alb[0]) * esc / (module[0] * derate * f_grid[1])])
    sea = whole['water']['albedo']
    floating = [(sea - pv_alb[1]) * t_up[0], (sea - pv_alb[0]) * t_up[1]]
    rooftops = levers_roofs(routes)
    return dict(
        whole=dict(toa_w_m2=whole['toa_w_m2'], ghi_w_m2=whole['ghi_w_m2'], ghi_clear_w_m2=whole['ghi_clear_w_m2'],
                   ground_over_top=ground_over_top, ground_over_top_clear=whole['ghi_clear_w_m2'] / whole['toa_w_m2'],
                   land_cloud=cloud),
        load_weighted_ghi_w_m2=ghi, by_band=by_band,
        albedo=dict(land=land_albedo, pv=list(pv_alb), change=list(d_alb), clear_sky_ground_over_top_equator=t_clear,
                    escape_clear=t_up_clear, escape_all_sky=list(t_up)),
        routes=routes, grid=grid, stores=stores_scale,
        sunlight_taken_from_the_ground_w_per_w=dict(beams_dimming=beam_light,
                                                    panels_global_grid=routes['global_grid']['surface_light_w_per_w'],
                                                    panels_hydrogen=routes['hydrogen']['surface_light_w_per_w']),
        aerial_above_the_cloud=dict(heat_w_per_w=aerial, transmission=list(t30),
                                    rule='panels at 30 km over a planet of albedo 0.40, through the night on '
                                         'hydrogen or a grid'),
        floating_on_the_seas=dict(albedo_change_times_escape=floating,
                                  rule="the design run's open water against panels of 0.05-0.10: near zero, the "
                                       "sign open"),
        rooftops=rooftops,
        first_screen_metropolis_panels_km2=screen['night']['level_panels_km2_for_the_mean'])


def levers_roofs(routes=None):
    roofs = [260000.0, 740000.0]
    out = dict(roof_and_pavement_km2=roofs, kind='project',
               basis="the heat-levers study's roofs and pavements of the settled land; a usable share of it carries "
                     "panels stored through the night in hydrogen")
    if routes is not None:
        use = lohi(val('roof_usable_share'))
        w = routes['hydrogen']['delivered_w_m2_of_panel']
        out['delivered_tw'] = [roofs[0] * 1e6 * use[0] * w[1] / TW, roofs[1] * 1e6 * use[1] * w[0] / TW]
    return out


# ---------------------------------------------------------------------------------------------------------------
# The Moon's other flows and the thermal plants.

def other_flows(climate, screen, power, dry_km2, sea_share, rings):
    rivers = screen['waters']['rivers']['gross_upper_bound_tw']
    dev = lohi(val('hydro_developable'))
    tides = [sum(v['one_way_gw'] for v in screen['waters']['seas'].values()) / 1e3,
             sum(v['both_ways_gw'] for v in screen['waters']['seas'].values()) / 1e3]
    flow = lohi(val('heat_flow_mw_m2'))
    geo_eff = lohi(val('geothermal_efficiency'))
    heat_tw = [f * 1e-3 * MOON_AREA / TW for f in flow]
    otec = lohi(val('otec_w_m2'))
    seas_km2 = sea_share * MOON_AREA / 1e6
    land = climate['whole']['land']
    p_s = power['absorption']['column']['surface_pressure_pa']
    air = th.earthlike_air(P.COLUMN['dry_pressure_pa'], P.COLUMN['co2_ppm'])
    t_air = climate['whole']['land']['air_k_weighted_by_sunlight']
    rho = p_s / (GAS_CONSTANT / air.molar_mass * t_air)
    density = 0.5 * rho * land['wind_cube_m3_s3']
    earth = lohi(val('earth_land_wind_power_density_w_m2'))
    gust = lohi(val('moon_gust_factor'))
    ratio = [density * gust[0] / earth[1], density * gust[1] / earth[0]]
    realized = [val('earth_wind_realized_w_m2') * r for r in ratio]
    dry_land_m2 = dry_km2 * 1e6
    cm1 = {}
    for name, ring in rings.items():
        rows = {}
        for where, stats, height in (('land_10_m', ring['wind_10m_m_s']['land'], 0.01),
                                     ('water_10_m', ring['wind_10m_m_s']['water'], 0.01),
                                     ('1_km', ring['wind_aloft_m_s']['1_km'], 1.0)):
            k, c = weibull(stats['median'], stats['p90'])
            fit = k >= 0.8                       # a near-calm ring's percentiles fit no useful Weibull
            cube = c ** 3 * math.gamma(1 + 3 / k) if fit else None
            # An upper bound from the percentiles alone: each slice of time at its slice's top speed.
            bound = (0.5 * stats['median'] ** 3 + 0.4 * stats['p90'] ** 3 + 0.09 * stats['p99'] ** 3
                     + 0.01 * stats['max'] ** 3)
            rows[where] = dict(median_m_s=stats['median'], p90_m_s=stats['p90'], p99_m_s=stats['p99'],
                               max_m_s=stats['max'], weibull_shape=k if fit else None,
                               weibull_scale_m_s=c if fit else None,
                               power_density_w_m2=0.5 * rho * cube if fit else None,
                               power_density_bound_w_m2=0.5 * rho * bound)
        cm1[name] = dict(latitude_deg=ring['latitude_deg'], **rows)
    return dict(
        hydro=dict(gross_tw=rivers, developable_tw=[rivers * dev[0] * val('turbine'), rivers * dev[1] * val('turbine')],
                   heat_w_per_w=0.0),
        tides=dict(tw=tides, heat_w_per_w=0.0),
        geothermal=dict(heat_flow_tw=heat_tw, electricity_if_all_taken_tw=[heat_tw[0] * geo_eff[0],
                                                                           heat_tw[1] * geo_eff[1]],
                        heat_w_per_w_mining_stored_heat=[1 / geo_eff[1], 1 / geo_eff[0]]),
        otec=dict(seas_km2=seas_km2, tw=[seas_km2 * 1e6 * otec[0] / TW, seas_km2 * 1e6 * otec[1] / TW],
                  heat_w_per_w=0.0),
        wind=dict(air_density_kg_m3=rho, mean_speed_m_s=land['wind_m_s'], share_below_cut_in=land['share_below_cut_in'],
                  power_density_w_m2=density, against_earth=ratio, realized_w_m2=realized,
                  km2_per_tw=[TW / r / 1e6 for r in realized[::-1]],
                  stress_dissipation_w_m2=land['stress_dissipation_w_m2'],
                  dissipation_over_dry_land_tw=land['stress_dissipation_w_m2'] * dry_land_m2 / TW,
                  cm1_rings=cm1,
                  cm1_rule="the cloud-resolving rings' instantaneous 3-hourly winds (the earlier design case): a "
                           "Weibull through each median and 90th percentile where its shape passes 0.8, and an upper "
                           "bound with each slice of time at its slice's top speed, at the design air's density",
                  heat_w_per_w=0.0))


def weibull(median, p90):
    """Shape and scale of the Weibull distribution through a median and a 90th percentile."""
    if median <= 0 or p90 <= median:
        return 0.0, 0.0
    k = math.log(math.log(10.0) / math.log(2.0)) / math.log(p90 / median)
    return k, median / math.log(2.0) ** (1.0 / k)


def thermal_plants(b):
    eff = lohi(val('thermal_efficiency_fission'))
    j_per_kg = val('mev_per_fission') * 1e6 * ELEMENTARY_CHARGE * AVOGADRO / 0.233
    year = HOURS_PER_YEAR * 3600.0
    fissioned = [TW * year / e / j_per_kg / 1e3 for e in eff[::-1]]           # t a year per TWe
    th_ppm = lohi(val('thorium_ppm'))
    ore = [fissioned[0] / (th_ppm[1] * 1e-6) / 1e9, fissioned[1] / (th_ppm[0] * 1e-6) / 1e9]   # Gt of rock a year
    fus = lohi(val('heat_per_w_fusion'))
    j_dt = val('mev_dt') * 1e6 * ELEMENTARY_CHARGE
    d_t = [TW * year * h / j_dt * 2.014e-3 / AVOGADRO / 1e3 for h in fus]                # t of deuterium a year
    li6 = [TW * year * h / j_dt * 6.015e-3 / AVOGADRO / 1e3 for h in fus]
    return dict(
        fission=dict(heat_w_per_w=[1 / eff[1], 1 / eff[0]], points=[1 / eff[1] / b['point']['heat_tw'],
                                                                     1 / eff[0] / b['point']['heat_tw']],
                     fissioned_t_per_twe_year=fissioned, thorium_rock_gt_per_twe_year=ore),
        fusion=dict(heat_w_per_w=list(fus), points=[f / b['point']['heat_tw'] for f in fus],
                    deuterium_t_per_twe_year=d_t, lithium6_t_per_twe_year=li6))


# ---------------------------------------------------------------------------------------------------------------
# Loads moved to orbit.

def orbit(levers, load, b, beams_h, ring_km):
    food = levers['beams']
    carbon = food['food_carbon']
    people = carbon['people']
    moved = carbon['moved_diet_share']
    losses = carbon['loss_factor']
    # As the heat-levers study counts the carbon: what is eaten at the low end, with all the food study's losses
    # falling after arrival at the high end.
    arrive = [moved[0], moved[1] * losses[1]]
    w_per_person = val('food_kcal') * 4184.0 / JULIAN_DAY
    food_energy = [people * w_per_person * arrive[e] / TW for e in (0, 1)]
    r_ring = ring_km * 1e3
    fall = MOON_GM * (1 / MOON_RADIUS - 1 / r_ring) + MOON_GM / (2 * r_ring)
    wet = lohi(val('food_wet_kg'))
    mass = [people * wet[e] * arrive[e] / JULIAN_DAY for e in (0, 1)]
    fall_tw = [m * fall / TW for m in mass]
    arriving = [food_energy[e] + fall_tw[e] for e in (0, 1)]
    removed = food['food_tw']
    computing = [load['delivered_tw']['close_computing'][e] * (1 - val('computing_kept')) * beams_h[e]
                 for e in (0, 1)]
    phosphorus = [people * arrive[e] * val('phosphorus_kg_person_year') / 1e9 for e in (0, 1)]
    per_w = [(food_energy[e] + mass[e] * fall / TW) / load['delivered_tw']['food'][e] for e in (0, 1)]
    return dict(
        food=dict(beams_heat_removed_tw=removed, food_energy_arriving_tw=food_energy, fall_tw=fall_tw,
                  fall_mj_per_kg=fall / 1e6, mass_kg_s=mass, net_heat_saved_tw=[removed[e] - arriving[e] for e in (0, 1)],
                  points_saved=[(removed[e] - arriving[e]) / b['point']['heat_tw'] for e in (0, 1)],
                  heat_on_the_moon_w_per_w_of_food_power_moved=per_w,
                  carbon_arriving_gt_yr=carbon['arriving_gt_c_yr'],
                  carbon_times_the_managed_import=carbon['times_the_managed_import'],
                  phosphorus_arriving_mt_yr=phosphorus),
        computing=dict(moved_share=1 - val('computing_kept'), beams_heat_removed_tw=computing,
                       points_saved=[c / b['point']['heat_tw'] for c in computing]),
        heat_in_orbit_w_m2_per_tw=levers['basis']['ring_w_m2_per_tw'])


# ---------------------------------------------------------------------------------------------------------------
# Mixes against the loads.

def mixes(load, b, beam_opts, solar, plants, orb, power, ledgers):
    sky, ground, total = load['sky_tw'], load['ground_tw'], load['total_tw']
    link = power['absorption']['link']
    h_t, h_g = link['towns_heat_per_w'], link['ground_heat_per_w']
    pt = b['point']['heat_tw']
    resort = b['resort_dimming_points']
    land_free = ledgers['land']['outside_kept_and_parks_km2']
    used = ledgers['land']['with_heritage']['1_kg']['used_km2']
    rec = beam_opts['ground_2_45_ghz'][0]['ground_closed_km2_per_tw']
    critical = power['reserves']['ring_fault']['energy_twh']
    out = {}

    def row(name, heat, land, store, notes):
        pts = [heat[e] / pt for e in (0, 1)]
        return dict(name=name, heat_tw=heat, w_m2=[x * b['moon_w_m2_per_tw'] for x in heat], points=pts,
                    share_of_the_resorts_dimming=[pts[0] / resort[1], pts[1] / resort[0]],
                    land_rain_percent=[pts[0] * b['point']['land_rain_percent'][0],
                                       pts[1] * b['point']['land_rain_percent'][1]],
                    growth_percent=[pts[0] * b['point']['growth_percent'][0], pts[1] * b['point']['growth_percent'][1]],
                    ground_km2=land, share_of_the_land_outside_kept_and_parks=[land[e] / land_free for e in (0, 1)],
                    against_the_land_now_used=[land[e] / used[e] for e in (0, 1)], store_twh=store, notes=notes)
    out['A'] = row('Beams as chosen: 2.45 GHz to the towns and the ground',
                   [sky[e] * h_t + ground[e] * h_g for e in (0, 1)],
                   [ground[e] * rec for e in (0, 1)], critical,
                   'a fault across the ring leaves the critical loads on a week of stores')
    f = orb['food']
    c = orb['computing']
    ground_b = [ground[e] - load['delivered_tw']['food'][e]
                - load['delivered_tw']['close_computing'][e] * c['moved_share'] for e in (0, 1)]
    out['B'] = row('Beams, with the beamed food and three quarters of close computing made at the ring radius',
                   [sky[e] * h_t + ground_b[e] * h_g + f['food_energy_arriving_tw'][e] + f['fall_tw'][e]
                    for e in (0, 1)],
                   [ground_b[e] * rec for e in (0, 1)], critical,
                   "food's carbon and phosphorus arrive with it and need return paths")
    grid = solar['routes']['global_grid']
    lines = solar['grid']
    h2 = solar['routes']['hydrogen']
    out['C'] = row('Beams to the towns; the ground on local solar with a grid round the Moon',
                   [sky[0] * h_t + ground[0] * grid['heat_w_per_w'][0],
                    sky[1] * h_t + ground[1] * grid['heat_w_per_w'][1]],
                   [ground[0] * (grid['land_km2_per_tw'][0] + lines['corridor_km2_per_tw'][0]),
                    ground[1] * (grid['land_km2_per_tw'][1] + lines['corridor_km2_per_tw'][1])],
                   [0.0, 0.0], 'lines carry the night side from the day side; the panels darken the ground')
    out['C']['line_km'] = [ground[e] * lines['line_km_per_tw'][e] for e in (0, 1)]
    out['D'] = row('Beams to the towns; the ground on local solar with hydrogen through the night',
                   [sky[0] * h_t + ground[0] * h2['heat_w_per_w'][0],
                    sky[1] * h_t + ground[1] * h2['heat_w_per_w'][1]],
                   [ground[0] * h2['land_km2_per_tw'][0], ground[1] * h2['land_km2_per_tw'][1]],
                   [ground[0] * h2['store_twh_per_tw'][0], ground[1] * h2['store_twh_per_tw'][1]],
                   'a store of about 300 TWh for each TW, refilled every lunar day')
    out['D']['hydrogen_mt'] = [ground[e] * solar['stores']['hydrogen_mt_per_tw'][e] for e in (0, 1)]
    fus = plants['fusion']['heat_w_per_w']
    fis = plants['fission']['heat_w_per_w']
    hot = [min(fus[0], fis[0]), max(fus[1], fis[1])]
    out['E'] = row('Beams to the towns; the ground on fission or fusion plants',
                   [sky[0] * h_t + ground[0] * hot[0], sky[1] * h_t + ground[1] * hot[1]],
                   [0.0, 0.0], [0.0, 0.0], "plants' land is small; their heat is two to three times what they deliver")
    out['F'] = row('Fusion for all, the planning assumption of 27 September',
                   [total[0] * fus[0], total[1] * fus[1]], [0.0, 0.0], [0.0, 0.0],
                   'a conditional technology; the towns would carry plants or take lines')
    return out


def options_table(b, beam_opts, opt, solar, flows, plants, orb, power):
    """Every option per delivered TW on the common measures, each figure with its evidence kind."""
    pt = b['point']['heat_tw']
    crit = power['reserves']['ring_fault']
    week_per_tw = [crit['energy_twh'][e] / crit['delivered_total_tw'][e] for e in (0, 1)]
    sps = val('sps_kg_per_kw')
    rect = lohi(val('rectenna_areal_kg_m2'))

    def heat(h):
        h = lohi(h)
        return dict(heat_w_per_w=list(h), points_per_tw=[h[0] / pt, h[1] / pt])

    def beam(key, kind):
        o = beam_opts[key]
        land = o[0]['ground_closed_km2_per_tw']
        return dict(**heat([o[0]['heat_w_per_w'], o[1]['heat_w_per_w']]), kind=kind,
                    ground_km2_per_tw=land,
                    rectenna_gt_per_tw=([land * 1e6 / val('receiver_buffer') * r / 1e12 for r in rect]
                                        if land else 0.0),
                    tile_electricity_w_per_w=[o[0]['tile_electricity_w_per_w'], o[1]['tile_electricity_w_per_w']],
                    notes=o[0]['notes'])
    r = solar['routes']
    rows = dict(
        beams_2_45_ghz=dict(**beam('ground_2_45_ghz', 'project and computed'),
                            store_twh_per_tw=week_per_tw, orbit_mt_per_tw_at_1978_reference=sps * 1e9 / 1e9,
                            failure='a fault across the ring is the one common mode; storms cost 1-23% of a beam',
                            maturity='rectennas measured at 82-91%; tens of kW over a few km on the ground; '
                                     'milliwatts in space (2023); no system built'),
        beams_5_8_ghz=dict(**beam('ground_5_8_ghz', 'project and computed'),
                           failure='storms stop it', maturity='as 2.45 GHz'),
        beams_1_ghz=dict(**beam('ground_1_ghz', 'computed (P.676 at its lower limit)'),
                         failure='as 2.45 GHz', maturity='as 2.45 GHz; apertures 2.45 times as wide'),
        beams_raised_30_km=dict(**beam('raised_30_km_2_45_ghz', 'computed on estimates'),
                                failure='a tether or its platform; lines in the flight bands',
                                maturity='tethered aerostats carry kilowatts to megawatts on Earth (recalled)'),
        beams_raised_30_km_5_8_ghz=dict(**beam('raised_30_km_5_8_ghz', 'computed on estimates'),
                                        failure='as raised 2.45 GHz; above the storms', maturity='as above'),
        optical=dict(**heat(opt['heat_w_per_w']), kind='estimate', ground_km2_per_tw='small',
                     tile_electricity_w_per_w=opt['tile_electricity_w_per_w'],
                     failure=f"cloud blocks a ground receiver {opt['share_of_time_cloud_blocks_a_ground_receiver']:.0%} "
                             f"of the time", maturity='800 W over 8.6 km (2025); 68.9% cells under laser light'),
        food_made_in_orbit=dict(heat_w_per_w=orb['food']['heat_on_the_moon_w_per_w_of_food_power_moved'],
                                points_per_tw=[x / pt for x in orb['food']['heat_on_the_moon_w_per_w_of_food_power_moved']],
                                kind='computed on the heat-levers product', ground_km2_per_tw=0.0,
                                failure='food security rides on the array and on arrivals',
                                maturity='controlled-environment farming and microbial protein exist on Earth; '
                                         'nothing at orbital scale'),
        local_solar_grid=dict(**heat(r['global_grid']['heat_w_per_w']), kind='computed on the design run and estimates',
                              panel_km2_per_tw=r['global_grid']['panel_km2_per_tw'],
                              ground_km2_per_tw=[r['global_grid']['land_km2_per_tw'][e]
                                                 + solar['grid']['corridor_km2_per_tw'][e] for e in (0, 1)],
                              store_twh_per_tw=0.0, line_km_per_tw=solar['grid']['line_km_per_tw'],
                              pv_gt_per_tw=r['global_grid']['pv_mass_gt_per_tw'],
                              failure='distributed; a cut in the lines darkens the night side',
                              maturity='PV and 1,100 kV HVDC are commercial; no grid girdles a world'),
        local_solar_hydrogen=dict(**heat(r['hydrogen']['heat_w_per_w']), kind='computed on the design run and estimates',
                                  panel_km2_per_tw=r['hydrogen']['panel_km2_per_tw'],
                                  ground_km2_per_tw=r['hydrogen']['land_km2_per_tw'],
                                  store_twh_per_tw=r['hydrogen']['store_twh_per_tw'],
                                  hydrogen_mt_per_tw=solar['stores']['hydrogen_mt_per_tw'],
                                  pv_gt_per_tw=r['hydrogen']['pv_mass_gt_per_tw'],
                                  failure='distributed; the store must carry a 354-hour night',
                                  maturity='electrolysers and fuel cells are commercial; no store of this size'),
        local_solar_pumped_water=dict(**heat(r['pumped_water']['heat_w_per_w']), kind='computed on estimates',
                                      ground_km2_per_tw=r['pumped_water']['land_km2_per_tw'],
                                      store_twh_per_tw=r['pumped_water']['store_twh_per_tw'],
                                      water_km3_per_tw=solar['stores']['pumped_water_km3_per_tw'],
                                      failure='the reservoirs swing hundreds of metres every lunar day',
                                      maturity='pumped hydro is commercial; no reservoir near this size'),
        local_solar_batteries=dict(**heat(r['batteries']['heat_w_per_w']), kind='computed on estimates',
                                   ground_km2_per_tw=r['batteries']['land_km2_per_tw'],
                                   store_twh_per_tw=r['batteries']['store_twh_per_tw'],
                                   battery_gt_per_tw=solar['stores']['battery_gt_per_tw'],
                                   failure='distributed', maturity='commercial at hours, not at 300 TWh a TW'),
        panels_above_the_cloud=dict(**heat(solar['aerial_above_the_cloud']['heat_w_per_w']['global_grid']),
                                    kind='estimate', failure='as aerial platforms', maturity='none at scale'),
        floating_panels=dict(heat_w_per_w=solar['floating_on_the_seas']['albedo_change_times_escape'],
                             kind='estimate', failure='as local solar',
                             maturity='floating PV is commercial on reservoirs'),
        wind=dict(heat_w_per_w=0.0, kind='computed on the design run', ground_km2_per_tw=flows['wind']['km2_per_tw'],
                  failure="the design run's low winds", maturity='commercial on Earth'),
        hydro=dict(heat_w_per_w=0.0, kind='project', total_tw=flows['hydro']['developable_tw'], maturity='commercial'),
        otec=dict(heat_w_per_w=0.0, kind='estimate from Earth by area', total_tw=flows['otec']['tw'],
                  maturity='pilot plants on Earth (recalled)'),
        tides=dict(heat_w_per_w=0.0, kind='project', total_tw=flows['tides']['tw'], maturity='commercial'),
        geothermal=dict(heat_w_per_w=[0.0, flows['geothermal']['heat_w_per_w_mining_stored_heat'][1]],
                        kind='published heat flow', total_tw=flows['geothermal']['electricity_if_all_taken_tw'],
                        maturity='commercial on Earth\'s hot spots; the Moon has none'),
        fission=dict(**heat(plants['fission']['heat_w_per_w']), kind='recalled efficiency',
                     ground_km2_per_tw='small', fissioned_t_per_twe_year=plants['fission']['fissioned_t_per_twe_year'],
                     failure='plant by plant', maturity='commercial; thorium breeding demonstrated (recalled)'),
        fusion=dict(**heat(plants['fusion']['heat_w_per_w']), kind='project (conditional technology)',
                    ground_km2_per_tw='small', failure='plant by plant',
                    maturity='target gain 4.13 at NIF (2025); ITER research 2034, D-T 2039; no plant has made '
                             'net electricity'))
    return rows


def trades(b, beam_opts, solar, orb, dry_km2):
    """What each way of easing the beams' heat buys and costs: points of dimming saved per TW it serves, and for local
    solar the extra ground per point saved, against the ground's beams ([better, worse] each)."""
    pt = b['point']['heat_tw']
    beam = beam_opts['ground_2_45_ghz']
    h_b = [beam[0]['heat_w_per_w'], beam[1]['heat_w_per_w']]
    land_b = beam[0]['ground_closed_km2_per_tw']
    out = {}
    for name in ('global_grid', 'hydrogen'):
        r = solar['routes'][name]
        corridor = solar['grid']['corridor_km2_per_tw'] if name == 'global_grid' else [0.0, 0.0]
        saved = [(h_b[0] - r['heat_w_per_w'][0]) / pt, (h_b[1] - r['heat_w_per_w'][1]) / pt]
        extra = [r['land_km2_per_tw'][0] + corridor[0] - land_b, r['land_km2_per_tw'][1] + corridor[1] - land_b]
        out[f'local_solar_{name}'] = dict(points_saved_per_tw=saved, extra_ground_km2_per_tw=extra,
                                          ground_km2_per_point_saved=[extra[0] / saved[0], extra[1] / saved[1]],
                                          share_of_dry_land_per_point_saved=[extra[0] / saved[0] / dry_km2,
                                                                             extra[1] / saved[1] / dry_km2])
    for h in val('raised_heights_km'):
        r = beam_opts[f'raised_{h:g}_km_2_45_ghz']
        out[f'raised_{h:g}_km'] = dict(points_saved_per_tw=[(h_b[0] - r[0]['heat_w_per_w']) / pt,
                                                            (h_b[1] - r[1]['heat_w_per_w']) / pt],
                                       ground_opened_km2_per_tw=land_b)
    food = orb['food']['heat_on_the_moon_w_per_w_of_food_power_moved']
    out['food_made_in_orbit'] = dict(points_saved_per_tw_of_food_power=[(h_b[0] - food[0]) / pt,
                                                                        (h_b[1] - food[1]) / pt],
                                     ground_opened_km2_per_tw=land_b)
    out['rule'] = ("against the ground's 2.45 GHz beams; a point of dimming costs the heat-levers study's 0.55-1.38% of "
                   "the land's rain and 0.07-0.30% of the growth through light")
    return out


# ---------------------------------------------------------------------------------------------------------------

def results():
    data = {p: load(p) for p in INPUTS}
    climate, power, screen, ledgers, levers, water = (data[p] for p in INPUTS[:6])
    rings = {k: data[v] for k, v in RINGS.items()}
    if climate['schema'] != 'terluna.provisioning.energy-climate/1':
        raise ValueError('energy_climate.json is not the extract this screen reads')
    load_ = loads(ledgers, power)
    b = basis(levers, power)
    e_g = power['reserves']['storms']['fades']['p99_rate']['ground']['operating']['elevation_deg']
    e_t = power['reserves']['storms']['fades']['p99_rate']['towns']['operating']['elevation_deg']
    elevations = (90.0, round(e_g, 2), round(e_t, 2), 30.0)
    table = air_table(elevations)
    beam_opts, geometry = beams(power, table, b)
    beams_h = [beam_opts['chosen_system'][e]['heat_w_per_w'] for e in (0, 1)]
    solar = local_solar(climate, power, b, beams_h, screen)
    flows = other_flows(climate, screen, power, ledgers['land']['dry_land_km2'], water['surface_water_share'], rings)
    plants = thermal_plants(b)
    orb = orbit(levers, load_, b, beams_h, ledgers['array_power']['ring_middle_radius_km'])
    opt = optical(power, b, climate)
    mix = mixes(load_, b, beam_opts, solar, plants, orb, power, ledgers)
    table_ = options_table(b, beam_opts, opt, solar, flows, plants, orb, power)
    trade = trades(b, beam_opts, solar, orb, ledgers['land']['dry_land_km2'])
    sources = json.loads(SOURCES.read_text())
    body = sig(dict(
        assumptions={k: dict(value=list(v['value']) if isinstance(v['value'], tuple) else v['value'], kind=v['kind'],
                             basis=v['basis']) for k, v in V.items()},
        loads=load_, basis=b, air=dict(columns=table, geometry=geometry,
                                       rule='dB and share lost from each receiver height to the column top; the '
                                            'vapour share is the water lines and continuum over the whole path'),
        beams=beam_opts, optical=opt, local_solar=solar, other_flows=flows, thermal_plants=plants, orbit=orb,
        options=table_, trades=trade, mixes=mix, sources_count=len(sources['sources'])), 4)
    return dict(schema=SCHEMA,
                producer=dict(domain='provisioning', files={p: sha256(p) for p in CODE},
                              inputs={p: sha256(p) for p in INPUTS},
                              sources={str(SOURCES.relative_to(ROOT)): sha256(SOURCES.relative_to(ROOT))},
                              constants=constants_used(list(CODE))),
                evidence=EVIDENCE, reading_rule=READING_RULE, **body)


def dumps(product) -> str:
    return json.dumps(product, indent=1) + '\n'


def main(argv=None) -> int:
    product = results()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(dumps(product))
    for k, m in product['mixes'].items():
        print(k, m['name'], m['heat_tw'], m['points'], m['ground_km2'])
    return 0


if __name__ == '__main__':
    sys.exit(main())
