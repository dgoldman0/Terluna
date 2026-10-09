"""A first screen of what the Open Moon's committed products already say about provisioning.

    python -m provisioning.screen
    # -> provisioning/results/first_screen.json

Six results, each read from committed products with the literature values named in sources.json:
- **Heat.** What the Moon may take from human use against the heat the ring fleet's own tiles radiate onto it. The
  tiles absorb 65-83 PW (solar_shield_array/results/array_heat.json) and face the Moon, so each face turned to the
  Moon sends it the share (R/r)^2 of its emission. Layers of tiles that hide one another from the Moon bound it from
  below. Warming comes from the GCM's response, 1.9 K per 1% of sunlight (climate/gcm), and 1.0-1.6 W/m2 per K.
- **Food land.** The canopy model's day fruit (biosphere/canopy/results/fruit.json) as food energy per square metre
  and the land a person needs, against the land round the summit port (geography atlas less its rain-fed lakes).
- **The night.** Carrying the summit metropolis through a 354-hour night on Korolev pumped storage (summit_tower).
- **Tides and rivers.** Tidal power per area against Earth's, the seas' theoretical tidal power, and an upper bound
  on the rivers' gross power.
- **Transit.** Accelerations people can stand in, scaled to lunar gravity, and the distances they set.
- **Fusion fuel.** The deuterium in the seas and lakes at Earth's ocean D/H.
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

from shared.constants import (AVOGADRO, MOON_RADIUS, MOON_SURFACE_GRAVITY, SPEED_OF_LIGHT, STANDARD_GRAVITY, STEFAN_BOLTZMANN,
                              SYNODIC_MONTH_DAYS)
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SCHEMA = 'terluna.provisioning.first-screen/1'
HEAT = ROOT / 'research' / 'studies' / 'solar_shield_array' / 'results' / 'array_heat.json'
LAYOUT = ROOT / 'research' / 'studies' / 'solar_shield_array' / 'results' / 'ring_layout.json'
FRUIT = ROOT / 'biosphere' / 'canopy' / 'results' / 'fruit.json'
TOWER = ROOT / 'research' / 'studies' / 'summit_tower' / 'results' / 'summit_tower.json'
ATLAS_GRID = ROOT / 'geography' / 'products' / 'atlas_28pct_4ppd.npz'
ATLAS = ROOT / 'geography' / 'results' / 'atlas.json'
DRAINAGE = ROOT / 'geography' / 'results' / 'drainage.json'
DRAINAGE_GRID = ROOT / 'geography' / 'products' / 'drainage_28pct_16ppd.npz'
TIDES = ROOT / 'geography' / 'results' / 'tides.json'
OUT = HERE / 'results' / 'first_screen.json'
INPUTS = (HEAT, LAYOUT, FRUIT, TOWER, ATLAS_GRID, ATLAS, DRAINAGE, DRAINAGE_GRID, TIDES)

# Heat.
GCM_K_PER_PERCENT_SUNLIGHT = 1.9          # climate/gcm/README.md, the dimmer branches
GCM_W_M2_PER_K = (1.0, 1.6)               # research/studies/atmospheric_co2/README.md, climate/gcm compare --settle
DESIGN_TOP_OF_AIR_W_M2 = 1173.0           # the 5% dimmer's sunlight at the top of the air (climate/gcm)
INFRARED_EFFICACY = (0.5, 1.0)            # absorbed aloft, the tiles' infrared may warm the ground less than sunlight
HEAT_MIRROR_EMISSIVITY = (0.03, 0.1)      # thermal emissivity of a Moon-facing heat mirror: silver stacks to oxide films
STRONG_EMITTER = 0.85                     # an outer face made a strong thermal emitter
INDUSTRY_TW = (100.0, 1000.0)             # orbital computing and industry, for scale
DATA_CENTRES_TW = 0.047                   # the world's data centres in 2024, 415 TWh (IEA 2025)
# Food.
FRUIT_KCAL_PER_G = 0.30                   # the day fruit has today's watermelon's make-up (USDA: 30 kcal per 100 g)
PERSON_KCAL_PER_DAY = 2500.0
CLEAR_SKY_SHARE = 0.89                    # the GCM's ground gets 89% of clear-sky sunlight (biosphere/canopy)
LOSS_FACTOR = (1.3, 1.5)                  # harvest, storage and household losses (FAO)
DIET_FACTOR = (1.5, 2.5)                  # a varied diet with protein and fats against the energy-only fruit
SUMMIT = (5.4125, 201.3665)               # deg N, deg E (decisions register)
METROPOLIS_PEOPLE = 100e6                 # decisions register, 2026-09-27
METROPOLIS_RADIUS_KM = 32.0
KW_PER_PERSON = 2.0                       # habitation/summit_metropolis brief
RADII_KM = (100.0, 300.0, 500.0)
# The night.
EARTH_PUMPED_STORAGE_TWH = 9.0            # world pumped-hydro energy capacity, of order (IHA)
# Tides and rivers.
SEAWATER = 1025.0
TIDAL_CYCLE_DAYS = 27.4                   # the forcing lines M' (27.55 d) and F (27.21 d); geography/README.md
EARTH_SEMIDIURNAL_H = 12.42
FRESH_WATER = 1000.0
# Transit.
STANDING_ACCELERATION_EARTH = (1.0, 1.5)  # m/s2 that standing riders take in service on Earth
LINE_SPEED = 30.0                         # m/s
# Fusion fuel.
VSMOW_D_PER_H = 155.76e-6
MEV_PER_D_CATALYSED = 43.2 / 6.0          # catalysed D-D burns six deuterons for 43.2 MeV
YEARS = 1e9

EVIDENCE = ('A screen on committed products: the shield study\'s array heat and ring layout, the canopy model\'s '
            'fruit, the summit tower\'s storage and panels, the atlas, drainage and tide products. Closed-form '
            'arithmetic with the stated literature values; no climate, grid, farm or reactor model is run.')
READING_RULE = ('heat in PW and W/m2 averaged over the Moon, warming in K; food in kcal per m2 per day and m2 per '
                'person; land in km2; night storage in TWh and km3; power in GW; transit in m/s2 and km; fusion '
                'fuel in kg, J and TW sustained over 1e9 years.')


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def two_layers(e_out, e_moon_outer, e_moon_inner):
    """Share of two stacked layers' heat that leaves toward the Moon, with each layer's light absorbed equally.

    The layers are grey and opaque in the thermal infrared, parallel and close: the outer layer radiates to space
    through e_out and to the inner layer through e_moon_outer, and the inner layer to the outer one through e_out
    and to the Moon through e_moon_inner. Returns the share and each layer's emitted flux per unit of absorbed flux.
    """
    k = 1.0 / (1.0 / e_moon_outer + 1.0 / e_out - 1.0)
    outer, inner = np.linalg.solve([[e_out + k, -k], [-k, e_moon_inner + k]], [1.0, 1.0])
    return e_moon_inner * inner / 2.0, (outer, inner)


def heat_mirror(films, window):
    """The Moon's share of the fleet's heat with a heat mirror on Moon-facing faces, for one layer and for two.

    About two layers overlap in the Moon's view, and a ring hidden behind another passes part of that ring's heat on
    inward. The mirror goes on every tile, or on the layer nearest the Moon alone; the outer faces stay as they are
    or become strong emitters. Temperatures are the hottest layer's, with the films' own absorbed sunlight.
    """
    weight = {n: (window if n == 'window' else 1 - window) for n in films}
    def mean(fn):
        return sum(weight[n] * fn(f) for n, f in films.items())
    def hottest(fn):
        return max(round(float((max(fn(f)[1]) * f['absorbed_W_m2'] / STEFAN_BOLTZMANN) ** 0.25)) for f in films.values())
    as_is = mean(lambda f: two_layers(f['emissivity_sunward'], f['emissivity_shaded'], f['emissivity_shaded'])[0])
    out = dict(two_layers_as_is=round(as_is, 3))
    for e in HEAT_MIRROR_EMISSIVITY:
        for outer_name, e_out in (('outer_faces_as_now', None), ('outer_faces_strong', STRONG_EMITTER)):
            eo = lambda f: e_out or f['emissivity_sunward']
            one = mean(lambda f: e / (e + eo(f)))
            every = lambda f: two_layers(eo(f), e, e)
            nearest = lambda f: two_layers(eo(f), f['emissivity_shaded'], e)
            out[f'mirror_{e:g}_{outer_name}'] = dict(
                one_layer_share=round(one, 3),
                every_tile=dict(share=round(mean(lambda f: every(f)[0]), 3), cut=round(as_is / mean(lambda f: every(f)[0]), 1),
                                hottest_layer_K=hottest(every)),
                nearest_layer_only=dict(share=round(mean(lambda f: nearest(f)[0]), 3), cut=round(as_is / mean(lambda f: nearest(f)[0]), 1),
                                        hottest_layer_K=hottest(nearest)))
    return out


def industry(array_heat, fleet, window, absorbed, on_moon, sunlight_global, area):
    """What computing and industry can do against the glow, and what they add to the Moon's heat.

    The glow is a fixed share of the heat the films absorb, so keeping a watt of it off the Moon means taking that
    share's inverse out of the films and rejecting it where the Moon does not see it. Computing and industry draw
    their power from light the fleet collects; their heat reaches the Moon by where it is released. The window rings
    carry the climate stack round their whole orbits and dim the Moon only while they cross the window, so light
    they absorb there heats the films far beyond the light it keeps from the Moon.
    """
    tiles, comp, placement = array_heat['tiles'], array_heat['computing'], array_heat['placement']
    share = [w * area / (a * 1e15) for w, a in zip(on_moon, absorbed)]
    absorbed_tw = [a * 1e3 for a in absorbed]
    per_tw_km2 = [c + comp['radiator_panel_km2_per_TW'] for c in comp['collector_km2_per_TW']]
    fleet_km2 = float(np.mean(fleet['area_km2']))
    window_per_moon = [window * i * 1e15 / (sunlight_global * area) for i in tiles['intercepted_PW']]
    return dict(
        glow_share_of_absorbed=[round(x, 5) for x in share],
        fleet_w_per_w_kept_off_moon=[round(1 / share[1]), round(1 / share[0])],
        glow_w_m2_per_percent_absorbed=[round(on_moon[0] / (100 * tiles['mean_solar_absorptance']), 2),
                                        round(on_moon[1] / (100 * tiles['mean_solar_absorptance']), 2)],
        computing_tw_to_halve_glow=[round(absorbed_tw[0] / 2, -2), round(absorbed_tw[1] / 2, -2)],
        data_centres_tw=DATA_CENTRES_TW,
        glow_cut_by_industry_tw={f'{tw:g}': [round(tw / absorbed_tw[1], 5), round(tw / absorbed_tw[0], 5)] for tw in INDUSTRY_TW},
        harvest_from_films_tw=dict(carnot_across_films=round(tiles['carnot_across_all_films_TW'], 2),
                                   thermoradiative_measured_everywhere=round(tiles['thermoradiative_measured_everywhere_TW'], 1)),
        fleet_area_share_for_industry_tw={f'{tw:g}': [round(tw * k / fleet_km2, 5) for k in per_tw_km2] for tw in INDUSTRY_TW},
        moon_w_m2_from_industry_tw_released_in_orbit={f'{tw:g}': round(tw * placement['moon_W_m2_per_TW_released_in_orbit'], 4) for tw in INDUSTRY_TW},
        moon_w_m2_from_industry_tw_used_on_moon={f'{tw:g}': round(tw * placement['moon_W_m2_per_TW_used_on_moon'], 2) for tw in INDUSTRY_TW},
        window_rings_intercept_per_moon_sunlight=[round(x, 1) for x in window_per_moon],
        absorbing_dimmer_gives_back=[round(x * g, 3) for x, g in zip(window_per_moon, share)])


def heat(array_heat, layout):
    films, fleet, design = array_heat['tiles']['films'], array_heat['fleet'], array_heat['design']
    window = fleet['window_ring_share']
    rings = layout['flown']['matched']['rings']
    radii = [min(r['radius_km'] for r in rings), max(r['radius_km'] for r in rings)]
    tilt = np.radians(max(abs(r['tilt_deg']) for r in rings))
    absorbed = array_heat['tiles']['absorbed_PW']
    # Each tile's absorbed light leaves through both faces in proportion to their emissivities. The tiles stay
    # radial-facing, so the face turned to the Moon is the same one all orbit: the shaded face at the sunward crossing.
    one_layer = sum((window if n == 'window' else 1 - window) * f['emissivity_shaded'] / (f['emissivity_sunward'] + f['emissivity_shaded'])
                    for n, f in films.items())
    # A ring that hides another from the Moon absorbs that ring's glow and sends it on through both its faces, so the
    # stack as a whole divides its heat between its innermost and outermost faces: evenly for a deep stack, by one
    # layer's emissivities for a single one. The share toward the Moon lies between the two.
    shares = (0.5, one_layer)
    view = [(MOON_RADIUS / 1e3 / r) ** 2 for r in radii[::-1]]           # (R/r)^2, outer ring first
    area = 4 * np.pi * MOON_RADIUS ** 2
    on_moon = [a * 1e15 * s * v / area for a, s, v in zip(absorbed, shares, view)]
    # The rings share one node line, so seen from the Moon they fill the lune of tilts within the largest tilt:
    # a share 2*tilt/pi of the sphere at their radius.
    band_km2 = 2 * tilt / np.pi * 4 * np.pi * np.mean(radii) ** 2
    layers = np.mean(fleet['area_km2']) / band_km2
    sunlight_global = DESIGN_TOP_OF_AIR_W_M2 / 4
    k_per_w = GCM_K_PER_PERCENT_SUNLIGHT / (0.01 * sunlight_global)    # K per W/m2 of top-of-air sunlight
    lo, hi = on_moon
    tw = 1e12 / area
    return dict(
        tiles_absorbed_pw=absorbed, moon_facing_share=dict(deep_stack=shares[0], one_layer=round(one_layer, 3)),
        ring_radius_km=[round(r) for r in radii], view_factor=[round(v, 5) for v in view],
        layers_in_the_moons_view=round(float(layers), 2), fleet_infrared_w_m2=[round(lo, 2), round(hi, 2)],
        heat_mirror=heat_mirror(films, window),
        industry=industry(array_heat, fleet, window, absorbed, on_moon, sunlight_global, area),
        share_of_design_sunlight=[round(lo / sunlight_global, 4), round(hi / sunlight_global, 4)],
        against_dimmer_w_m2=round(sunlight_global / 0.95 * 0.05, 1),
        warming_k_if_like_sunlight=[round(lo * k_per_w, 1), round(hi * k_per_w, 1)],
        warming_k_range=[round(lo * INFRARED_EFFICACY[0] / GCM_W_M2_PER_K[1], 1), round(hi * INFRARED_EFFICACY[1] / GCM_W_M2_PER_K[0], 1)],
        human_use=dict(w_m2_per_tw=round(tw, 4), k_per_tw=[round(tw / GCM_W_M2_PER_K[1], 4), round(tw / GCM_W_M2_PER_K[0], 4)],
                       tw_for_one_kelvin=[round(GCM_W_M2_PER_K[0] / tw), round(GCM_W_M2_PER_K[1] / tw)],
                       metropolis_gw=METROPOLIS_PEOPLE * KW_PER_PERSON / 1e6,
                       metropolis_w_m2_moon=round(METROPOLIS_PEOPLE * KW_PER_PERSON * 1e3 / area, 5),
                       metropolis_w_m2_local=round(METROPOLIS_PEOPLE * KW_PER_PERSON * 1e3 / (np.pi * (METROPOLIS_RADIUS_KM * 1e3) ** 2), 1),
                       fleet_infrared_as_tw_used=[round(lo / tw), round(hi / tw)]))


def land_round_summit(grid, drainage_grid):
    lat, lon = grid['lat_deg'], grid['lon_deg']
    assert lat[0] > lat[-1] and abs(lon[0]) < 1, 'atlas rows from the north, columns from 0 E'
    land = (grid['water_label'] == 0).astype(float)
    lakes16 = np.zeros(2880 * 5760, bool)
    lakes16[drainage_grid['lake_cell']] = True
    lakes = lakes16.reshape(720, 4, 1440, 4).mean(axis=(1, 3))
    dry = np.clip(land - lakes, 0, 1)
    la, lo = np.radians(lat)[:, None], np.radians(lon)[None, :]
    cell = (np.radians(0.25) ** 2) * np.cos(la) * (MOON_RADIUS / 1e3) ** 2 * np.ones_like(lo)
    s = np.radians(SUMMIT)
    dist = MOON_RADIUS / 1e3 * np.arccos(np.clip(np.sin(la) * np.sin(s[0]) + np.cos(la) * np.cos(s[0]) * np.cos(lo - s[1]), -1, 1))
    return {f'{r:g}': dict(dry_land_km2=round(float(np.sum(cell * dry * (dist <= r))), -2),
                           lakes_km2=round(float(np.sum(cell * lakes * (dist <= r))), -2)) for r in RADII_KM}


def food(fruit, round_summit):
    per_cycle_days = SYNODIC_MONTH_DAYS
    sites = {}
    for side, bands in fruit['results'].items():
        for band in bands:
            for strategy in ('earth_like', 'idle_quarter'):
                shares = fruit['results'][side][band][strategy]['day_fruit']['best']['fed']['shares']
                for key in ('0.5', '0.6'):
                    kg = shares[key]['harvest_kg']
                    kcal = kg * 1000 * FRUIT_KCAL_PER_G / per_cycle_days
                    sites[f'{side}/{band}/{strategy}/{key}'] = dict(harvest_kg_per_cycle=round(kg, 2), kcal_m2_day=round(kcal, 1),
                                                                    m2_per_person=round(PERSON_KCAL_PER_DAY / kcal, 0))
    best = min(v['m2_per_person'] for k, v in sites.items() if k.startswith('far/0'))
    worst = max(v['m2_per_person'] for k, v in sites.items() if k.startswith('far/0'))
    planning = [best / CLEAR_SKY_SHARE * LOSS_FACTOR[0] * DIET_FACTOR[0], worst / CLEAR_SKY_SHARE * LOSS_FACTOR[1] * DIET_FACTOR[1]]
    need = [p * METROPOLIS_PEOPLE / 1e6 for p in (best, worst)]
    need_planning = [p * METROPOLIS_PEOPLE / 1e6 for p in planning]
    metropolis_km2 = np.pi * METROPOLIS_RADIUS_KM ** 2
    return dict(sites=sites, far_side_equator_m2_per_person=[best, worst],
                planning_m2_per_person=[round(x) for x in planning],
                metropolis_cropland_km2=dict(potential=[round(x, -2) for x in need], planning=[round(x, -2) for x in need_planning]),
                metropolis_area_km2=round(metropolis_km2), land_round_summit_km2=round_summit,
                planning_share_of_dry_land=round_shares(need_planning, round_summit))


def round_shares(need, round_summit):
    return {r: [round(n / v['dry_land_km2'], 2) for n in need] for r, v in round_summit.items()}


def night(tower):
    storage = tower['night_storage']
    hours = storage['night_hours']
    energy_twh = METROPOLIS_PEOPLE * KW_PER_PERSON * 1e3 * hours / 1e12
    water_km3 = energy_twh * 1e9 / storage['delivered_kwh_per_m3'] / 1e9
    panels_km2 = METROPOLIS_PEOPLE * KW_PER_PERSON * 1e3 / tower['sunlight']['horizontal_panel_w_m2'] / 1e6
    return dict(night_hours=round(hours, 1), metropolis_night_twh=round(energy_twh, 1),
                times_earths_pumped_storage=round(energy_twh / EARTH_PUMPED_STORAGE_TWH, 1),
                korolev_water_km3=round(water_km3, 1), head_m=storage['head_m'],
                level_panels_km2_for_the_mean=round(panels_km2, -2))


def waters(tides, atlas, drainage, grid):
    cycle = TIDAL_CYCLE_DAYS * 86400
    per_area_ratio = MOON_SURFACE_GRAVITY / STANDARD_GRAVITY * EARTH_SEMIDIURNAL_H * 3600 / cycle
    seas = {}
    for key, sea in tides['seas'].items():
        a, h = float(sea['area_km2']) * 1e6, sea['monthly_range_median_m']['area_median']
        p = 0.5 * SEAWATER * MOON_SURFACE_GRAVITY * a * h ** 2 / cycle
        if a > 1e11:
            seas[', '.join(sea['names'][:2]) or key] = dict(area_km2=int(float(sea['area_km2'])), range_m=round(h, 2),
                                                           one_way_gw=round(p / 1e9, 1), both_ways_gw=round(2 * p / 1e9, 1))
    heights = grid['height_m']
    land = grid['water_label'] == 0
    lat = np.radians(grid['lat_deg'])[:, None] * np.ones_like(heights)
    mean_height = float(np.sum((heights - atlas['sea_level_m']) * land * np.cos(lat)) / np.sum(land * np.cos(lat)))
    q = drainage['rivers']['discharge_to_sea_m3s']
    return dict(tidal_power_per_area_moon_over_earth=round(per_area_ratio, 5), seas=seas,
                rivers=dict(discharge_m3s=q, mean_land_height_above_sea_m=round(mean_height),
                            gross_upper_bound_tw=round(FRESH_WATER * MOON_SURFACE_GRAVITY * q * mean_height / 1e12, 2)))


def transit():
    scale = MOON_SURFACE_GRAVITY / STANDARD_GRAVITY
    a = [x * scale for x in STANDING_ACCELERATION_EARTH]
    return dict(standing_acceleration_m_s2=dict(earth=STANDING_ACCELERATION_EARTH, moon=[round(x, 3) for x in a]),
                stop_from_m_s=LINE_SPEED, stopping_km=dict(earth=[round(LINE_SPEED ** 2 / (2 * x) / 1e3, 2) for x in STANDING_ACCELERATION_EARTH[::-1]],
                                                           moon=[round(LINE_SPEED ** 2 / (2 * x) / 1e3, 2) for x in a[::-1]]),
                tilt_of_apparent_vertical_deg_at_1_m_s2=dict(earth=round(np.degrees(np.arctan(1 / STANDARD_GRAVITY)), 1),
                                                              moon=round(np.degrees(np.arctan(1 / MOON_SURFACE_GRAVITY)), 1)))


def fusion(atlas, drainage):
    water = atlas['water_mass_kg'] + drainage['lakes']['volume_km3'] * 1e12
    deuterons = 2 * water / 0.018015 * AVOGADRO * VSMOW_D_PER_H
    joules = deuterons * MEV_PER_D_CATALYSED * 1.602176634e-13
    return dict(water_kg=float(f'{water:.3g}'), deuterium_kg=float(f'{deuterons / AVOGADRO * 2.014e-3:.3g}'),
                energy_j=float(f'{joules:.3g}'), tw_for_1e9_years=round(joules / (YEARS * 3.15576e7) / 1e12, 1),
                note='The seas and lakes are imported; the operation can deliver more deuterium with any water.')


def latency():
    return dict(round_trip_s_to_ring_radius=round(2 * 20000e3 / SPEED_OF_LIGHT, 3))


def main(argv=None) -> int:
    array_heat, layout, fruit, tower, atlas, drainage, tides = (json.loads(p.read_text()) for p in (HEAT, LAYOUT, FRUIT, TOWER, ATLAS, DRAINAGE, TIDES))
    grid, dgrid = np.load(ATLAS_GRID), np.load(DRAINAGE_GRID)
    files = ['provisioning/screen.py']
    product = dict(
        schema=SCHEMA,
        producer=dict(domain='provisioning', files={f: digest(ROOT / f) for f in files},
                      inputs={str(p.relative_to(ROOT)): digest(p) for p in INPUTS}, constants=constants_used(files)),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        heat=heat(array_heat, layout),
        food=food(fruit, land_round_summit(grid, dgrid)),
        night=night(tower),
        waters=waters(tides, atlas, drainage, grid),
        transit=transit(),
        fusion_fuel=fusion(atlas, drainage),
        latency=latency())
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    print(json.dumps({k: v for k, v in product.items() if k not in ('producer', 'evidence', 'reading_rule')}, indent=1)[:7000])
    return 0


if __name__ == '__main__':
    sys.exit(main())
