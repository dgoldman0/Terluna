"""Where and when light limits growth on the Open Moon, and what the planned dimming costs the biosphere through it.

    python -m research.studies.light_limitation.run
    # -> results/light_limitation.json

On 13 October 2026 the author asked where and when light is the factor that limits growth, if anywhere: it bears on
the gas pair, and it is worth knowing in general. Light limits growth where nothing else does first, so this screen
takes each part of the living Moon in turn and reads from the project's products what limits it there:

- the land, band by band: the share of the canopy's light-limited growth the rain supports, and the cap the
  sunlight sets on what watered land can transpire (the food systems study, from ecology's crop model and the design
  climate); the GCM's own soil-water stress beside them; and the wet margins where water gathers (ecology's scenario);
- the waters: the nutrient-poor open seas and the river-fed coasts of the living-water designs, and the lakes, which
  keep nearly all the phosphorus the land sends them;
- the sky: the aerophytes, whose growth follows their light in the sky ecology;
- the cycle: the canopy's response to the light through the lunar day and the night's store (canopy light response).

Each place's light elasticity, the share of its growth lost for each share of light taken, comes from the canopy
light response where the canopy model applies and from the literature named beside it elsewhere. Weighted by each
place's production it gives the share of the biosphere's growth that each per cent of dimming costs through light,
and the cost of the dimming the plan now carries. The dimming also weakens the rain, which this screen sets beside
the light from Earth's models; the Moon's own response needs its GCM.
"""
from __future__ import annotations
import hashlib
import json
import math
import sys
from pathlib import Path
import numpy as np

from shared.constants import JULIAN_YEAR_DAYS, SYNODIC_MONTH_DAYS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'light_limitation.json'
SCHEMA = 'terluna.research.light-limitation/1'
INPUTS = {
    'biosphere/canopy/results/light_response.json': 'terluna.biosphere.canopy-light-response/1',
    'research/studies/food_systems/results/food_systems.json': 'terluna.research.food-systems/1',
    'research/studies/food_systems/results/food_blend.json': 'terluna.research.food-blend/1',
    'biosphere/ecology/results/first_screen.json': None,
    'biosphere/living_water/waters.json': 'terluna.biosphere.living-water/1',
    'research/studies/sky_ecology/results/sky_ecology.json': 'terluna.research.sky-ecology/1',
    'research/studies/domain_synthesis/results/ledgers.json': 'terluna.research.domain-synthesis-ledgers/1',
}
# The GCM design case's climatology, a generated product outside Git, bound by its hash as the plant model binds it.
CLIMATOLOGY = 'climate/gcm/products/climatology_A28_dim5_moon.npz'
CLIMATOLOGY_SCHEMA = 'terluna.climate.gcm-climatology/1'
CYCLES_PER_YEAR = JULIAN_YEAR_DAYS / SYNODIC_MONTH_DAYS
STRATEGIES = ('earth_like', 'idle_quarter')       # the two evergreen ways through the night, the potential's two ends
CANOPY_LATITUDE = {'0-15': '0', '15-30': '30', '30-45': '30', '45-60': '60', '60-75': '60', '75-90': '60'}
# ExoPlaSim's land surface (landmod.f90, as the design run sets it): evaporation runs below its potential where the
# soil holds less than drhsfull = 0.4 of the field capacity wsmax = 0.5 m.
GCM_STRESS_M = 0.4 * 0.5
# The waters, design guesses on Earth's waters: net primary production of Earth's oligotrophic open ocean and of its
# shelves and upwelling coasts (Field et al. 1998), the living-water designs' coasts as a share of the seas, and the
# lakes', which keep nearly all the land's phosphorus, from mesotrophic to eutrophic lakes.
OPEN_SEA_G_C_M2_YR = (50.0, 150.0)
COAST_G_C_M2_YR = (300.0, 1000.0)
COAST_SHARE_OF_SEAS = (0.05, 0.15)
LAKE_G_C_M2_YR = (100.0, 400.0)
# Light elasticities outside the canopy model. A nutrient-poor sea: the light term of the vertically generalised
# production model, E0/(E0 + 4.1) with E0 in mol photons m-2 per day (Behrenfeld & Falkowski 1997), whose elasticity
# is 4.1/(E0 + 4.1), at the Moon's daylight surface light (by_sun_height). A self-shaded, nutrient-replete water column
# (the lakes and the coasts): Talling's (1957) areal production, proportional to ln(2 E0/Ek), with Ek 100-300 umol
# m-2 s-1 at fixed biomass, up to 1 where the biomass grows with the light.
VGPM_E = 4.1
EK_UMOL = (100.0, 300.0)
# Earth's balanced dimming (GeoMIP G1, twelve models): an abrupt quadrupling of CO2 offset by solar reductions of
# 3.5-5.0%, which held the global temperature within 0.3 K and cut global precipitation by 4.5 +/- 1.3%, and land
# evaporation by about 10% with CO2 closing the stomata (Tilmes et al. 2013, Table 1 and section 4).
G1_SOLAR_REDUCTION = (3.8, 4.0, 4.7, 4.1, 4.3, 4.5, 3.9, 4.1, 3.5, 5.0, 4.7, 4.0)
G1_PRECIPITATION = (4.5 - 1.3, 4.5 + 1.3)
# The fleet's glow, answered by rolling both halves of each orbit outside the service arc: its added dimming on
# research/array-industry (8be6d1e, solar_shield_array/glow.md), which this branch does not carry.
GLOW_DIMMING_PERCENT = 2.64
EVIDENCE = ("A screen of committed products and stated guesses. The land is the food systems study's band basis "
            "(ecology's crop model, the design climate's rain and the sunlight's cap on transpiration) and ecology's "
            "wet-margin scenario, with the GCM's soil water beside it; the canopy's light response is the plant "
            "model's; the waters' production and the elasticities outside the canopy model are Earth's, named in "
            "sources.json; the rain's response to dimming is Earth's GeoMIP G1. Production is the canopy model's "
            "growth, clear-sky potential times the rain's or the sunlight's share; nothing is simulated here.")
READING_RULE = ("Pairs are [low, high]. Production is Gt C a year of growth (new tissue after respiration, the canopy "
                "model's growth and the waters' and sky's net primary production); areas km2; elasticity the share "
                "of growth lost per share of light taken, 1 where growth follows the light and 0 where it ignores it. "
                "The land's potential takes the plant model's two evergreen strategies at the equator, Earth-like and "
                "idling at a quarter, as its two ends. limit names what holds growth where it is: water, light (with "
                "the sunlight's energy for transpiration where that binds first) or nutrients. The weighted "
                "elasticity's low end pairs the light-limited places' low production and elasticity with the rest's "
                "high production, and its high end the reverse.")


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path, schema=None):
    data = json.loads((ROOT / path).read_text())
    if schema and data.get('schema') != schema:
        raise ValueError(f'{path}: unexpected schema {data.get("schema")}')
    return data


def rnd(x, sig=3):
    if isinstance(x, dict):
        return {k: rnd(v, sig) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [rnd(v, sig) for v in x]
    if isinstance(x, bool) or not isinstance(x, float) or x == 0.0 or not math.isfinite(x):
        return x
    return round(x, sig - 1 - int(math.floor(math.log10(abs(x)))))


def pair(f):
    return [f(0), f(1)]


def gt_c(km2, g_m2_yr):
    return km2 * 1e6 * g_m2_yr / 1e15


def gcm_soil_water():
    """The design climate's land soil water by band of latitude: the share of land below the GCM's own stress."""
    raw = (ROOT / CLIMATOLOGY).read_bytes()
    z = np.load(ROOT / CLIMATOLOGY)
    meta = json.loads(str(z['metadata']))
    if meta['schema'] != CLIMATOLOGY_SCHEMA or meta['run'] != 'A28_dim5_moon':
        raise ValueError('unexpected climatology')
    lat, land, water, rain = z['lat'].astype(float), z['lsm'] > 0.5, z['mrso_m'], z['pr_mm_day']
    w = np.cos(np.radians(lat))[:, None] * np.ones(water.shape)
    out = {}
    for lo, hi in ((0, 15), (15, 30), (30, 45), (45, 60), (60, 75), (75, 90)):
        rows = ((np.abs(lat) >= lo) & (np.abs(lat) < hi))[:, None] & land
        a = float((w * rows).sum())
        out[f'{lo}-{hi}'] = dict(
            land_share_below_stress=float((w * rows * (water < GCM_STRESS_M)).sum() / a),
            mean_soil_water_m=float((w * rows * water).sum() / a), mean_rain_mm_day=float((w * rows * rain).sum() / a))
    return out, hashlib.sha256(raw).hexdigest()


def land(food, response, screen, moon_km2):
    """The land band by band: the rain-fed growth, the light-limited potential with water ample, and the wet margins
    where water gathers, with what limits each and how its growth follows the light."""
    rain, bands = food['basis']['rain'], food['basis']['zones']['bands']
    cloud = food['basis']['diet']['cloud_share']
    resp = response['bands']
    g_eq = [resp['0'][s]['base']['growth'] * CYCLES_PER_YEAR for s in STRATEGIES]     # g C m-2 a year
    growth_e = {b: [resp[b][s]['elasticity']['0.05']['growth'] for s in STRATEGIES] for b in resp}
    margins_km2 = [s * moon_km2 for s in screen['trace_gases']['methane']['wetland_share_scenario']]
    wet_belt = [b for b in bands if b['abs_lat_deg'][1] <= 30.0]
    belt_km2 = sum(b['dry_land_km2'] for b in wet_belt)
    out, rows = {}, []
    for b in bands:
        key = f"{b['abs_lat_deg'][0]:g}-{b['abs_lat_deg'][1]:g}"
        area, gpp, cap = b['dry_land_km2'], rain['gpp_factor'][key], rain['watered_yield_cap'][key]
        support = rain['support_scaled'][key]
        light_share = gpp * cloud                          # the canopy's growth under the design clouds
        energy_binds = [cap[e] < light_share for e in (0, 1)]
        e_light = growth_e[CANOPY_LATITUDE[key]]
        watered_e = [1.0 if energy_binds[e] else e_light[e] for e in (0, 1)]
        share = [min(light_share, cap[e]) for e in (0, 1)]
        margin = [m * area / belt_km2 if b in wet_belt else 0.0 for m in margins_km2]
        rain_fed = pair(lambda e: gt_c(area - margin[e], g_eq[e] * gpp * support[e]))
        out[key] = dict(dry_land_km2=area, wet_share=rain['wet_share'][key], rain_supports=support,
                        light_potential_share_of_equator=light_share, sunlight_cap_share_of_equator=cap,
                        watered_limit=['energy for transpiration' if x else 'light' for x in energy_binds],
                        watered_elasticity=watered_e, wet_margin_km2=margin, rain_fed_gt_c=rain_fed,
                        wet_margin_gt_c=pair(lambda e: gt_c(margin[e], g_eq[e] * share[e])),
                        potential_gt_c=pair(lambda e: gt_c(area, g_eq[e] * share[e])))
        rows.append(out[key])
    margin_e = [sum(r['wet_margin_gt_c'][e] * r['watered_elasticity'][e] for r in rows) /
                sum(r['wet_margin_gt_c'][e] for r in rows) for e in (0, 1)]
    totals = dict(rain_fed_gt_c=pair(lambda e: sum(r['rain_fed_gt_c'][e] for r in rows)),
                  wet_margin_gt_c=pair(lambda e: sum(r['wet_margin_gt_c'][e] for r in rows)),
                  potential_gt_c=pair(lambda e: sum(r['potential_gt_c'][e] for r in rows)),
                  wet_margin_km2=margins_km2, wet_margin_elasticity=margin_e, cloud_share=cloud,
                  equator_growth_g_c_m2_yr=g_eq)
    totals['rain_fed_share_of_potential'] = pair(lambda e: (totals['rain_fed_gt_c'][e] + totals['wet_margin_gt_c'][e])
                                                 / totals['potential_gt_c'][e])
    return dict(bands=out, totals=totals)


def daylight(response, latitude_deg):
    """Photosynthetic light at the surface through the sunlit half of the cycle at a latitude, umol m-2 s-1, from
    the light by Sun height (the solved sky's low Sun included)."""
    s = response['by_sun_height']
    h, par = np.array(s['elevation_deg']), np.array(s['par_umol_m2_s'])
    hour = np.linspace(-90.0, 90.0, 3601)
    elevation = np.degrees(np.arcsin(np.cos(np.radians(latitude_deg)) * np.cos(np.radians(hour))))
    return np.interp(elevation, h, par)


def water_column_elasticity(light):
    """The light elasticity of a nutrient-poor sea, from the production model's light term weighted by production
    through the day, and of a self-shaded, nutrient-replete column at fixed biomass, from Talling's ln(2 E0/Ek) at
    each Ek."""
    e0 = light * 86400.0 / 1e6                              # as a day's photons at that rate, mol m-2
    f = e0 / (e0 + VGPM_E)
    sea = float(np.sum(f * VGPM_E / (e0 + VGPM_E)) / np.sum(f))
    column = []
    for ek in EK_UMOL:
        lit = light > ek / 2.0
        prod = np.log(2.0 * light[lit] / ek)
        column.append(float(np.sum(prod * (1.0 / prod)) / np.sum(prod)))
    return sea, column


def waters(screen, living, response, moon_km2, sea_share):
    lakes_km2 = screen['trace_gases']['methane']['lake_share'] * moon_km2
    seas_km2 = sea_share * moon_km2
    by_lat = {lat: water_column_elasticity(daylight(response, lat)) for lat in (0.0, 60.0)}
    sea_e = sorted(v[0] for v in by_lat.values())
    column_fixed = min(min(v[1]) for v in by_lat.values())
    coast = [s * seas_km2 for s in COAST_SHARE_OF_SEAS]
    return dict(
        lakes=dict(km2=lakes_km2, gt_c=[gt_c(lakes_km2, g) for g in LAKE_G_C_M2_YR], limit='light',
                   elasticity=[column_fixed, 1.0],
                   basis='the lakes keep 96% of the phosphorus the land sends them (food systems, sea phosphorus)'),
        open_sea=dict(km2=pair(lambda e: seas_km2 - coast[1 - e]),
                      gt_c=pair(lambda e: gt_c(seas_km2 - coast[1 - e], OPEN_SEA_G_C_M2_YR[e])),
                      limit='nutrients', elasticity=sea_e,
                      chlorophyll_mg_m3=living['waters']['open_sea']['chlorophyll_mg_m3']),
        coasts=dict(km2=coast, gt_c=pair(lambda e: gt_c(coast[e], COAST_G_C_M2_YR[e])), limit='light',
                    elasticity=[column_fixed, 1.0],
                    chlorophyll_mg_m3=living['waters']['productive_coast']['chlorophyll_mg_m3']),
        elasticity_by_latitude={f'{k:g}': dict(open_sea=v[0], column_fixed_biomass=v[1]) for k, v in by_lat.items()},
        daylight_mean_umol_m2_s={f'{lat:g}': float(daylight(response, lat).mean()) for lat in (0.0, 60.0)})


def sky(sky_ecology, response):
    totals = sky_ecology['food_web']['totals']
    lo = min(min(t['reference_light_biomass_mt_c_year']) for t in totals if t['cover'] == min(r['cover'] for r in totals))
    hi = max(max(t['reference_light_biomass_mt_c_year']) for t in totals if t['cover'] == max(r['cover'] for r in totals))
    thin = response['bands']['0']['by_lai']['1']['elasticity']['0.05']
    whole = response['bands']['0']['earth_like']['elasticity']['0.05']['growth']
    return dict(gt_c=[lo / 1e3, hi / 1e3], covers=[min(t['cover'] for t in totals), max(t['cover'] for t in totals)],
                limit='light', elasticity=[thin, whole],
                basis="aerophytes at the sky ecology's reference light, 0.1-2% cover, across its classes")


def cycle(response):
    """When light limits: the deep canopy's response by the Sun's height and over the cycle's parts, and the night's
    store."""
    s = response['by_sun_height']
    at = {f'{h:g}': e for h, e in zip(s['elevation_deg'], s['elasticity']['0.05'])}
    parts = {b: response['bands'][b]['by_part_at_5_percent'] for b in response['bands']}
    store = {b: [response['bands'][b][st]['elasticity']['0.05']['store'] for st in STRATEGIES]
             for b in response['bands']}
    gpp = {b: [response['bands'][b][st]['elasticity']['0.05']['gpp'] for st in STRATEGIES][0]
           for b in response['bands']}
    growth = {b: [response['bands'][b][st]['elasticity']['0.05']['growth'] for st in STRATEGIES]
              for b in response['bands']}
    by_lai = {b: {k: v['elasticity']['0.05'] for k, v in response['bands'][b]['by_lai'].items()}
              for b in response['bands']}
    return dict(by_sun_height=at, by_part=parts, gpp=gpp, growth=growth, store=store, by_lai=by_lai)


def ledger(land_, waters_, sky_):
    """Production by place, what limits it and how it follows the light; the weighted elasticity."""
    places = {
        'rain_fed_land': dict(gt_c=land_['totals']['rain_fed_gt_c'], limit='water', elasticity=[0.0, 0.0]),
        'wet_margins': dict(gt_c=land_['totals']['wet_margin_gt_c'], limit='light',
                            elasticity=land_['totals']['wet_margin_elasticity']),
        'lakes': {k: waters_['lakes'][k] for k in ('gt_c', 'limit', 'elasticity')},
        'coasts': {k: waters_['coasts'][k] for k in ('gt_c', 'limit', 'elasticity')},
        'open_sea': {k: waters_['open_sea'][k] for k in ('gt_c', 'limit', 'elasticity')},
        'sky': {k: sky_[k] for k in ('gt_c', 'limit', 'elasticity')},
    }
    lit = [k for k, v in places.items() if v['limit'] == 'light']
    rest = [k for k in places if k not in lit]
    total = pair(lambda e: sum(v['gt_c'][e] for v in places.values()))

    def weighted(e_lit, e_rest):
        num = sum(places[k]['gt_c'][e_lit] * places[k]['elasticity'][e_lit] for k in lit) + \
            sum(places[k]['gt_c'][e_rest] * places[k]['elasticity'][e_rest] for k in rest)
        den = sum(places[k]['gt_c'][e_lit] for k in lit) + sum(places[k]['gt_c'][e_rest] for k in rest)
        return num / den
    shares = {k: pair(lambda e: v['gt_c'][e] / total[e]) for k, v in places.items()}
    light_share = [sum(places[k]['gt_c'][0] for k in lit) / (sum(places[k]['gt_c'][0] for k in lit)
                                                             + sum(places[k]['gt_c'][1] for k in rest)),
                   sum(places[k]['gt_c'][1] for k in lit) / (sum(places[k]['gt_c'][1] for k in lit)
                                                             + sum(places[k]['gt_c'][0] for k in rest))]
    return dict(places=places, total_gt_c=total, share_of_production=shares, light_limited_share=light_share,
                weighted_elasticity=[weighted(0, 1), weighted(1, 0)])


def dimming(ledgers):
    gases = ledgers['gases']['held']['added_dimming_percent']
    heat = ledgers['heat']['array']['with_close_computing']['added_dimming_percent']
    total = [gases + heat[0] + GLOW_DIMMING_PERCENT, gases + heat[1] + GLOW_DIMMING_PERCENT]
    return dict(gases_percent=gases, people_heat_percent=heat, glow_percent=GLOW_DIMMING_PERCENT, total_percent=total,
                basis=('the gases at the held pair and people\'s heat with power from the array and close computing '
                       '(domain synthesis); the glow with both halves of each orbit rolled outside the service arc '
                       '(research/array-industry)'))


def rain_pathway(ledger_, dimming_):
    """Earth's balanced dimming beside the light: the rain that a dimming answering a warming takes on Earth, and what
    the same share would cost the rain-fed land. It applies to the dimming that answers warming in the air: the
    gases', and at most the glow's, part of which reaches the ground through the infrared window. People's heat is
    released at the ground, where the dimming that answers it takes about as much energy from the surface as the heat
    adds, so it is left out."""
    mean_solar = sum(G1_SOLAR_REDUCTION) / len(G1_SOLAR_REDUCTION)
    per_percent = [G1_PRECIPITATION[0] / mean_solar, G1_PRECIPITATION[1] / mean_solar]
    water = sorted(ledger_['share_of_production']['rain_fed_land'])
    in_air = [dimming_['gases_percent'], dimming_['gases_percent'] + dimming_['glow_percent']]
    per = [per_percent[0] * water[0], per_percent[1] * water[1]]
    return dict(g1_solar_reduction_percent=[min(G1_SOLAR_REDUCTION), max(G1_SOLAR_REDUCTION)],
                g1_mean_solar_reduction_percent=mean_solar, g1_precipitation_cut_percent=list(G1_PRECIPITATION),
                precipitation_cut_per_percent_dimming=per_percent, rain_fed_share_of_production=water,
                dimming_answering_warming_in_the_air_percent=in_air,
                rain_cut_if_the_moon_matched_percent=[per_percent[0] * in_air[0], per_percent[1] * in_air[1]],
                growth_cut_per_percent_dimming_if_the_moon_matched=per,
                growth_cut_at_the_planned_dimming_percent=[per[0] * in_air[0], per[1] * in_air[1]])


def results():
    data = {p: load(p, s) for p, s in INPUTS.items()}
    response = data['biosphere/canopy/results/light_response.json']
    food = data['research/studies/food_systems/results/food_systems.json']
    screen = data['biosphere/ecology/results/first_screen.json']
    moon_km2 = food['basis']['zones']['moon_km2']
    soil, climatology_sha = gcm_soil_water()
    land_ = land(food, response, screen, moon_km2)
    waters_ = waters(screen, data['biosphere/living_water/waters.json'], response, moon_km2,
                     food['basis']['zones']['water_share'])
    sky_ = sky(data['research/studies/sky_ecology/results/sky_ecology.json'], response)
    ledger_ = ledger(land_, waters_, sky_)
    dimming_ = dimming(data['research/studies/domain_synthesis/results/ledgers.json'])
    blend = data['research/studies/food_systems/results/food_blend.json']['blend']['8']['wild_30_hybrid_20']
    food_ = dict(weights=blend['weights'], fields_km2=blend['fields_km2'], land_km2=blend['land_km2'],
                 irrigation_tw=blend['power_tw']['irrigation'], lamps_tw=blend['power_tw']['lamps'],
                 rain_used_km3_year=blend['water_km3_year']['rain_used'])
    cost = [ledger_['weighted_elasticity'][e] * dimming_['total_percent'][e] for e in (0, 1)]
    return dict(
        schema=SCHEMA,
        producer=dict(files={'research/studies/light_limitation/run.py': sha256(__file__)},
                      inputs={p: sha256(ROOT / p) for p in INPUTS},
                      ignored_inputs={CLIMATOLOGY: climatology_sha}),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        guesses=dict(open_sea_g_c_m2_yr=OPEN_SEA_G_C_M2_YR, coast_g_c_m2_yr=COAST_G_C_M2_YR,
                     coast_share_of_seas=COAST_SHARE_OF_SEAS, lake_g_c_m2_yr=LAKE_G_C_M2_YR, ek_umol_m2_s=EK_UMOL,
                     gcm_stress_soil_water_m=GCM_STRESS_M),
        land=rnd(land_), gcm_soil_water=rnd(soil), food=rnd(food_), waters=rnd(waters_), sky=rnd(sky_),
        cycle=rnd(cycle(response)), ledger=rnd(ledger_), dimming=rnd(dimming_),
        cost_through_light=rnd(dict(per_percent_dimming=ledger_['weighted_elasticity'],
                                    at_the_planned_dimming_percent=cost)),
        rain_pathway=rnd(rain_pathway(ledger_, dimming_)))


def dumps(product) -> str:
    return json.dumps(product, indent=2, allow_nan=False) + '\n'


def main(argv=None) -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(dumps(results()))
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
