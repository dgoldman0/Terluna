"""Whole-plant carbon over the lunar cycle: what is left for growth, and the store the night needs, for ways
of getting through the 354-hour night.

    python -m biosphere.canopy.plant        # writes results/plant.json (about a minute, one thread)

The base stand (model.py: leaf area index 5) is followed through one solar cycle at each latitude band,
with its leaves at the land air temperature of the chosen climate through the lunar day, on the near and
far sides (climate/gcm/climatology.py, run A28_dim5, years 15-24, composited by hour angle).

Light. Above 30 degrees of Sun height the canopy gets the surface-light product's spectra. Lower down the
two-stream product falls short of the spherical sky atlas (illumination/sky), so its diffuse light is
scaled until the total matches the atlas; and after sunset, when the Moon's tall sky stays sunlit far above
the ground, the canopy gets the atlas's twilight illuminance with the spectrum of the diffuse light at
1 degree, down to 35 degrees below the horizon. The atlas uses unfiltered sunlight, so its twilight is
reduced by the shield's share at 1 degree.

Carbon. Maintenance respiration is the leaves' own (the leaf model's) plus that of stems and roots, which
follows Q10 = 2 and is set so that the same stand under Earth's sky, with the same temperatures over a
24-hour day, turns half of its gross photosynthesis into new tissue (carbon-use efficiency 0.5; Waring et
al. 1998). New tissue costs a quarter of its carbon again in growth respiration. The smallest store that
carries the plant through the cycle is the largest cumulative deficit of photosynthesis less maintenance
(long_night.periodic_storage_requirement).

Strategies: respiring in the dark as by day (earth_like); slowing the maintenance of leaves, stems and roots
to a half, a quarter or a tenth whenever photosynthesis no longer covers it (idle_*); and shedding the
canopy when photosynthesis at dusk no longer covers its upkeep and growing it again from sunrise
(regrow_*), which saves the leaves' night respiration but pays for a new canopy every cycle (90 g of dry
leaf per m2 of leaf, the PROSPECT default, 47% carbon, plus growth respiration) with the leaf area rising
linearly over the first days of daylight.
"""
from __future__ import annotations
from dataclasses import replace
import hashlib
import json
import math
from pathlib import Path
import numpy as np

from shared.constants import SYNODIC_MONTH_DAYS
from biosphere.long_night import periodic_storage_requirement
from biosphere.canopy import model as canopy

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / 'results'
SCHEMA = 'terluna.biosphere.plant-carbon-cycle/1'
CLIMATOLOGY = ROOT / 'climate' / 'gcm' / 'products' / 'climatology_A28_dim5.npz'
CLIMATOLOGY_SCHEMA = 'terluna.climate.gcm-climatology/1'
CLIMATE_YEARS = (15, 24)
ATLASES = dict(moon=ROOT / 'illumination' / 'sky' / 'data' / 'moon_no_ozone_atlas.npz',
               earth=ROOT / 'illumination' / 'sky' / 'data' / 'earth_atlas.npz')
ATLAS_CASES = dict(moon='moon_1.2atm/unfiltered', earth='earth_control/unfiltered')
BANDS = {'0': (0.0, 15.0), '30': (20.0, 40.0), '60': (50.0, 70.0)}      # |latitude| of the land rows
SIDES = ('near', 'far')
Q10 = 2.0
CUE_EARTH = 0.5
CUE_RANGE = (0.4, 0.6)
GROWTH_RESPIRATION = 0.25
LEAF_DRY_G_M2 = 90.0
CARBON_FRACTION = 0.47
RAMP_DAYS = 5.0
RAMP_RANGE = (2.0, 10.0)
LAIS = (1.0, 2.0, 3.0, 4.0, 5.0)
TEMPS_C = np.arange(10.0, 36.0, 1.0)
TWILIGHT_DEG = (-35.0, -30.0, -25.0, -20.0, -18.0, -15.0, -12.0, -10.0, -8.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0)
SCALED_BELOW_DEG = 30.0
DARK_UMOL = 1.0
STEPS = 5760
MOON_DAY_S = SYNODIC_MONTH_DAYS * 86400.0
GRAMS = canopy.GRAMS_C_PER_UMOL
LUX_Y = np.array([0.2126, 0.7152, 0.0722])
STRATEGIES = dict(earth_like=(1.0, False), idle_half=(0.5, False), idle_quarter=(0.25, False),
                  idle_tenth=(0.1, False), regrow_leaves=(1.0, True), regrow_leaves_idle_quarter=(0.25, True))


def temperatures():
    """Land air temperature (C) through the lunar day for each side and latitude band, and the climatology."""
    raw = CLIMATOLOGY.read_bytes()
    data = np.load(CLIMATOLOGY)
    meta = json.loads(str(data['metadata']))
    if meta['schema'] != CLIMATOLOGY_SCHEMA or 'tas_land_sun_near' not in data.files:
        raise ValueError('Regenerate the climatology: python -m climate.gcm.climatology A28_dim5:15-24')
    if meta['run'] != 'A28_dim5' or tuple(meta['years']) != CLIMATE_YEARS:
        raise ValueError(f"Expected run A28_dim5, years {CLIMATE_YEARS}; found {meta['run']} {meta['years']}")
    lat, hour = data['lat'].astype(float), data['hour_angle_deg'].astype(float)
    traces = {}
    for side in SIDES:
        field = data[f'tas_land_sun_{side}'].astype(float) - 273.15
        for band, (lo, hi) in BANDS.items():
            rows = (np.abs(lat) >= lo) & (np.abs(lat) < hi)
            f = field[rows]
            w = np.cos(np.radians(lat[rows]))[:, None] * np.isfinite(f)
            traces[side, band] = np.nansum(np.nan_to_num(f) * w, axis=0) / np.maximum(w.sum(axis=0), 1e-300)
    return hour, traces, hashlib.sha256(raw).hexdigest()[:16]


def along_cycle(hour_centres, values, steps=STEPS):
    """A periodic hour-angle trace at the steps of one cycle (from midnight; noon at the middle)."""
    hour = (np.arange(steps) + 0.5) / steps * 360.0 - 180.0
    xp = np.concatenate([hour_centres - 360.0, hour_centres, hour_centres + 360.0])
    return np.interp(hour, xp, np.tile(values, 3))


def light():
    """Photon spectra at the ground over a black ground, per world, at Sun heights from -35 to 90 degrees.

    Returns the wavelengths, per world dict(elevations, direct, diffuse, reflectance_from_below, par), and the
    hashes of the sources.
    """
    lam, worlds, light_sha = canopy.surface_light()
    product = json.loads(canopy.LIGHT.read_bytes())
    atlas_rows = {w: {r['sun_deg']: r['total_ratio'] for r in product['checks']['atlas_comparison'][ATLAS_CASES[w]]['rows']}
                  for w in worlds}
    par = (lam > canopy.PAR_NM[0]) & (lam < canopy.PAR_NM[1])
    shas = dict(surface_light=light_sha)
    out = {}
    for w, src in worlds.items():
        case = product['cases'][canopy.WORLDS[w]['case']]['summaries']
        bare = product['cases'][ATLAS_CASES[w]]['summaries']
        heights = src['heights']
        direct, diffuse = src['direct'].copy(), src['diffuse'].copy()
        for i, h in enumerate(heights):
            if h < SCALED_BELOW_DEG:
                s = case[str(h)]['0.1']
                total, beam = s['illuminance_lux'], s['direct_illuminance_lux']
                diffuse[i] *= (atlas_rows[w][h] * total - beam) / (total - beam)
        atlas = np.load(ATLASES[w])
        shas[f'atlas_{w}'] = hashlib.sha256(ATLASES[w].read_bytes()).hexdigest()[:16]
        lux = (atlas['direct_horizontal'] + atlas['diffuse']) @ LUX_Y
        low = int(np.argmin(np.abs(heights - 1.0)))
        s1 = case['1.0']['0.1']
        shield = s1['illuminance_lux'] / bare['1.0']['0.1']['illuminance_lux']
        per_lux = src['diffuse'][low] / (s1['illuminance_lux'] - s1['direct_illuminance_lux'])
        twilight = np.array([per_lux * shield * float(np.interp(e, atlas['suns'], lux)) for e in TWILIGHT_DEG])
        order = np.argsort(heights)
        elevations = np.concatenate([TWILIGHT_DEG, heights[order]])
        out[w] = dict(elevations=elevations,
                      direct=np.vstack([np.zeros_like(twilight), direct[order]]),
                      diffuse=np.vstack([twilight, diffuse[order]]),
                      reflectance_from_below=src['reflectance_from_below'])
        out[w]['par'] = (out[w]['direct'] + out[w]['diffuse'])[:, par].sum(axis=1)
    return lam, out, shas


class Tables:
    """Canopy gross photosynthesis by leaf area, leaf temperature and Sun height (umol m-2 s-1 of ground)."""

    def __init__(self):
        lam, self.light, self.shas = light()
        heights = self.light['moon']['elevations']
        positive = heights[heights > 0]
        optics = canopy.Optics(lam, positive)
        self.gpp = {}
        for world, src in self.light.items():
            co2 = canopy.WORLDS[world]['co2_pa']
            for lai in (LAIS if world == 'moon' else (5.0,)):
                stand = replace(canopy.Stand(), lai=lai)
                resp = optics(stand)
                direct_resp = [resp['direct'][int(np.argmin(np.abs(positive - max(e, positive.min()))))] for e in heights]
                self.gpp[world, lai] = np.array([[canopy.rates(stand, (direct_resp[i], resp['diffuse']), lam,
                                                               src['direct'][i], src['diffuse'][i],
                                                               src['reflectance_from_below'], co2, t_c=t)['gpp_umol_m2_s']
                                                  for i in range(heights.size)] for t in TEMPS_C])
        self.elevations = heights

    def photosynthesis(self, world, lai, elevation, temp):
        """Gross photosynthesis at each step, interpolated in Sun height, temperature and leaf area."""
        k = np.clip(np.searchsorted(TEMPS_C, temp) - 1, 0, TEMPS_C.size - 2)
        w = np.clip((temp - TEMPS_C[k]) / (TEMPS_C[k + 1] - TEMPS_C[k]), 0.0, 1.0)
        steps = np.arange(elevation.size)

        def at(table):
            by_t = np.array([np.interp(elevation, self.elevations, row, left=0.0) for row in table])
            return (1 - w) * by_t[k, steps] + w * by_t[k + 1, steps]
        lai = np.broadcast_to(np.asarray(lai, dtype=float), elevation.shape)
        if np.all(lai == 5.0):
            return at(self.gpp[world, 5.0])
        grid = np.array([0.0, *LAIS])
        values = np.vstack([np.zeros(elevation.size)] + [at(self.gpp[world, a]) for a in LAIS])
        j = np.clip(np.searchsorted(grid, lai) - 1, 0, grid.size - 2)
        f = (lai - grid[j]) / (grid[j + 1] - grid[j])
        return (1 - f) * values[j, steps] + f * values[j + 1, steps]

    def par(self, world, elevation):
        return np.interp(elevation, self.elevations, self.light[world]['par'], left=0.0)


def leaf_respiration(lai, temp):
    """Leaf respiration of a stand of the given leaf area (umol m-2 s-1 of ground), in proportion to leaf area."""
    base = np.array([canopy.respiration(canopy.Stand(), t) for t in TEMPS_C])
    return np.interp(temp, TEMPS_C, base) * np.asarray(lai, dtype=float) / 5.0


def cycle(tables, world, latitude_deg, temp, solar_day_s, stem_root_25, night=1.0, regrow=False,
          ramp_days=RAMP_DAYS, idle_when='deficit', twilight=True):
    """Carbon over one solar cycle (g C per m2 of ground) for one way of getting through the night."""
    hour = (np.arange(STEPS) + 0.5) / STEPS * 360.0 - 180.0
    elevation = np.degrees(np.arcsin(math.cos(math.radians(latitude_deg)) * np.cos(np.radians(hour))))
    up = elevation > 0
    dt = solar_day_s / STEPS
    stem_root = stem_root_25 * Q10 ** ((temp - 25.0) / 10.0)
    lit = up if not twilight else np.ones(STEPS, dtype=bool)
    full_canopy = tables.photosynthesis(world, 5.0, elevation, temp) * lit
    upkeep = leaf_respiration(5.0, temp) + stem_root
    if regrow:
        # Leaves fall when the afternoon's photosynthesis no longer covers their upkeep and grow again from sunrise.
        afternoon = hour > 0
        falls = hour[afternoon & (full_canopy < upkeep)].min() if np.any(afternoon & (full_canopy < upkeep)) else 90.0
        days_since_sunrise = (hour + 90.0) / 360.0 * solar_day_s / 86400.0
        leaves_on = (hour >= -90.0) & (hour < falls)
        lai = np.where(leaves_on, 5.0 * np.clip(days_since_sunrise / ramp_days, 0.0, 1.0), 0.0)
    else:
        lai = np.full(STEPS, 5.0)
    gpp = tables.photosynthesis(world, lai, elevation, temp) * lit
    full = leaf_respiration(lai, temp) + stem_root
    idle = (gpp < full) if idle_when == 'deficit' else (tables.par(world, elevation) < DARK_UMOL)
    slow = np.where(idle, night, 1.0)
    r_leaf = leaf_respiration(lai, temp) * slow
    r_stem_root = stem_root * slow
    build = np.zeros(STEPS)
    if regrow:
        cost = 5.0 * LEAF_DRY_G_M2 * CARBON_FRACTION * (1 + GROWTH_RESPIRATION) / GRAMS     # umol C per m2
        growing = leaves_on & (lai < 5.0)
        build[growing] = cost / (growing.sum() * dt)
    net = gpp - r_leaf - r_stem_root - build
    store = periodic_storage_requirement(net, dt)
    total = lambda x: float(np.sum(x) * dt * GRAMS)
    growth = max(float(store['cycle_net'] * GRAMS), 0.0) / (1 + GROWTH_RESPIRATION)
    par = tables.par(world, elevation)
    return dict(gpp=total(gpp), twilight_gpp=total(gpp * ~up),
                leaf_respiration=total(r_leaf), stem_root_respiration=total(r_stem_root),
                respiration_in_the_dark=total((r_leaf + r_stem_root) * ~up), new_canopy=total(build),
                growth=growth, growth_per_24h=growth * 86400.0 / solar_day_s,
                carbon_use_efficiency=growth / total(gpp), balance_closes=bool(store['cycle_balance_feasible']),
                store=float(store['worst_single_cycle_deficit'] * GRAMS),
                deficit_hours=float(store['deficit_duration_days'] / 3600.0),
                idle_hours=float(np.sum(idle) * dt / 3600.0) if night < 1.0 else 0.0,
                dark_hours=float(np.sum(par < DARK_UMOL) * dt / 3600.0))


def calibrate(tables, temp_trace, cue=CUE_EARTH):
    """Stem and root maintenance at 25 C (umol m-2 s-1) giving the Earth stand the carbon-use efficiency cue."""
    probe = cycle(tables, 'earth', 0.0, temp_trace, 86400.0, 0.0)
    unit = cycle(tables, 'earth', 0.0, temp_trace, 86400.0, 1.0)
    needed = probe['gpp'] * (1 - (1 + GROWTH_RESPIRATION) * cue) - probe['leaf_respiration']
    if needed <= 0:
        raise ValueError(f'Leaf respiration alone exceeds the budget for a carbon-use efficiency of {cue}')
    return needed / unit['stem_root_respiration']


def twilight_table(tables):
    """Photosynthetic light after sunset and the hours at the equator until it falls below given levels."""
    rate = 360.0 / MOON_DAY_S * 3600.0                      # degrees of hour angle per hour
    out = {}
    for world, src in tables.light.items():
        e = src['elevations']
        dusk = e <= 0
        per_hour = rate if world == 'moon' else 15.0
        below = {}
        down, light_ = e[dusk][::-1], src['par'][dusk][::-1]          # from the horizon downward
        for level in (100.0, 50.0, 10.0, 1.0):
            if light_[0] < level:
                below[str(level)] = 0.0
                continue
            k = int(np.argmax(light_ < level)) if np.any(light_ < level) else light_.size - 1
            f = (math.log(level) - math.log(light_[k - 1])) / (math.log(light_[k]) - math.log(light_[k - 1]))
            below[str(level)] = float(-(down[k - 1] + f * (down[k] - down[k - 1])) / per_hour)
        out[world] = dict(elevation_deg=e[dusk].tolist(), par_umol_m2_s=src['par'][dusk].tolist(),
                          hours_after_sunset_above=below)
    return out


EVIDENCE = ('A carbon budget for one stand over one solar cycle: the canopy model\'s photosynthesis, with the light '
            'below 30 degrees and after sunset taken from the spherical sky atlas (unfiltered sunlight on an '
            'exponential column, scaled by the shield), at the chosen climate\'s land air temperatures through the '
            'lunar day (3-day GCM means composited by hour angle); leaf respiration from the leaf model; stem and root '
            'maintenance set from an Earth carbon-use efficiency of 0.5 (0.4-0.6 tested) with Q10 = 2. The '
            'strategies are scenarios: whether leaves can idle without dark-induced senescence, how fast metabolism '
            'can slow and recover, and how fast a canopy can regrow are not modelled. Clear sky; no water limits, '
            'nutrients, reproduction or losses to herbivores.')
READING_RULE = ('Carbon is g C per m2 of ground over one solar cycle (the Moon: the synodic month; Earth: a 24-hour '
                'day with the same temperature trace). growth is new tissue after maintenance and growth respiration; '
                'store is the smallest reserve that carries the plant through the cycle and deficit_hours how long it '
                'runs down; twilight_gpp is photosynthesis with the Sun below the horizon; dark_hours has less than '
                '1 umol m-2 s-1 of photosynthetic light. results[side][band][strategy]; bands are |latitude| ranges of '
                'land. earth is the calibration stand. twilight gives the light after sunset and, at the equator, '
                'the hours after sunset it stays above each level. sensitivity holds the equatorial near side for '
                'other Earth carbon-use efficiencies, canopy regrowth times, idling only in the dark, and every '
                'strategy with the light cut off at sunset (without_twilight).')


def run():
    hour_centres, traces, climate_sha = temperatures()
    tables = Tables()
    equator = along_cycle(hour_centres, traces['near', '0'])
    stem_root_25 = calibrate(tables, equator)
    earth = cycle(tables, 'earth', 0.0, equator, 86400.0, stem_root_25)
    results, climate = {}, {}
    for side in SIDES:
        for band in BANDS:
            temp = along_cycle(hour_centres, traces[side, band])
            day = np.abs(hour_centres) < 90
            climate.setdefault(side, {})[band] = dict(
                day_mean_c=float(traces[side, band][day].mean()), night_mean_c=float(traces[side, band][~day].mean()),
                coldest_c=float(traces[side, band].min()), warmest_c=float(traces[side, band].max()),
                by_hour_angle_c=traces[side, band].tolist())
            results.setdefault(side, {})[band] = {
                name: cycle(tables, 'moon', float(band), temp, MOON_DAY_S, stem_root_25, night=n, regrow=r)
                for name, (n, r) in STRATEGIES.items()}
    sensitivity = dict(carbon_use_efficiency={}, regrowth_days={})
    for cue in CUE_RANGE:
        r25 = calibrate(tables, equator, cue)
        sensitivity['carbon_use_efficiency'][str(cue)] = dict(
            stem_root_25=r25, **{name: cycle(tables, 'moon', 0.0, equator, MOON_DAY_S, r25, night=n, regrow=r)
                                 for name, (n, r) in STRATEGIES.items()})
    for ramp in RAMP_RANGE:
        sensitivity['regrowth_days'][str(ramp)] = cycle(tables, 'moon', 0.0, equator, MOON_DAY_S, stem_root_25,
                                                        regrow=True, ramp_days=ramp)
    sensitivity['idle_quarter_only_in_the_dark'] = cycle(tables, 'moon', 0.0, equator, MOON_DAY_S, stem_root_25,
                                                         night=0.25, idle_when='dark')
    sensitivity['without_twilight'] = {name: cycle(tables, 'moon', 0.0, equator, MOON_DAY_S, stem_root_25, night=n,
                                                   regrow=r, twilight=False)
                                       for name, (n, r) in STRATEGIES.items()}
    digest = lambda p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()[:16]
    return dict(
        schema=SCHEMA,
        producer=dict(domain='biosphere', files={p: digest(p) for p in (
            'biosphere/canopy/plant.py', 'biosphere/canopy/model.py', 'biosphere/canopy/radiation.py',
            'biosphere/canopy/leaf.py', 'biosphere/canopy/leaf_optics.py', 'biosphere/long_night.py')},
            inputs=dict(**{f'{k}_sha256': v for k, v in tables.shas.items()}, climatology_sha256=climate_sha,
                        climatology_schema=CLIMATOLOGY_SCHEMA, climatology_run='A28_dim5',
                        climatology_years=CLIMATE_YEARS)),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        units=dict(carbon='g C m-2 of ground per solar cycle', temperature='C', light='umol m-2 s-1 of 400-700 nm'),
        settings=dict(q10=Q10, earth_carbon_use_efficiency=CUE_EARTH, growth_respiration=GROWTH_RESPIRATION,
                      leaf_dry_g_m2=LEAF_DRY_G_M2, carbon_fraction=CARBON_FRACTION, regrowth_days=RAMP_DAYS,
                      stem_root_respiration_25c=stem_root_25, bands_abs_latitude=BANDS,
                      atlas_scaling_below_deg=SCALED_BELOW_DEG, twilight_deg=TWILIGHT_DEG, dark_umol=DARK_UMOL,
                      strategies={k: dict(dark_share=v[0], regrow_leaves=v[1]) for k, v in STRATEGIES.items()},
                      stand=canopy.Stand().__dict__),
        twilight=twilight_table(tables), climate=climate, earth=earth, results=results, sensitivity=sensitivity)


def main():
    product = canopy._round(run())
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / 'plant.json').write_text(json.dumps(product, indent=1) + '\n')
    print('wrote', RESULTS / 'plant.json')


if __name__ == '__main__':
    main()
