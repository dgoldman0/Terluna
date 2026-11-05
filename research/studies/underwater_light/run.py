"""Light under the Open Moon's seas and lakes through the lunar day and night.

    nice -n 10 env OPENBLAS_NUM_THREADS=1 python3 -m research.studies.underwater_light.run     # 1-2 minutes, 0.7 GB

The first run builds the floors' depth tables from geography's heights (seconds, 1.5 GB) and caches them on the drive.

The solved spherical sky (research/runs/optical_comfort/spherical, the sea-appearance study's source) gives the direct
beam and diffuse light on the water in each spectral channel for every solar elevation from -90 to 90 degrees,
including the long twilight with the Sun below the horizon. The near side's earthlight is the same per-source response
at the Earth's elevation, weighted by the calendar's Earth spectrum at full phase (illumination/calendar), as the
calendar's own export weights it. Through a flat surface (Fresnel) and down the column (the domain's water optics with
fresh water's scattering and ultraviolet absorption, Lee et al.'s attenuation per stream; optics.py) the light becomes
spectral downwelling irradiance with depth, photosynthetic photons, illuminance and a fish rod's photon catch.

The Sun's arc is the present axis's at equinox (declination 0; the 1.54 degree seasonal swing is a sensitivity at
75 degrees), so latitude stays a parameter: a place moved by the turned Moon reads its new latitude, and the lit
floors are summed twice, at today's latitudes and at the turned Moon's (the elastic frame set for the second set of
simulations), over today's depths and lakes. Each water's
chlorophyll runs from its dawn value at sunrise to its dusk value at sunset and back through the night, mixed in the
logarithm of the light. The ultraviolet comes from the surface-light product (1 nm, Sun 1-90 degrees, two-stream).
The floor's depths are geography's: the 16 pixel/degree heights against the atlas's sea level, as geography/drainage.py
builds them, and the drainage product's rain-fed lakes. The Earth control runs the same waters under Earth's sky.

writes results/underwater_light.json (schema terluna.research.underwater-light/1) and the depth tables and light
tables to /media/projectspace/terluna-research/underwater-light/.
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np

from biosphere.ecology.screen import govardovskii
from illumination.water_column import model as domain
from research.studies.underwater_light import optics
from shared.constants import AVOGADRO, MOON_EQUATOR_TO_ECLIPTIC_DEG, PLANCK, SPEED_OF_LIGHT, SYNODIC_MONTH_DAYS
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SPHERICAL = ROOT / 'research/runs/optical_comfort/spherical'
SKY = {'moon': SPHERICAL / 'moon_1.2atm_standard.npz', 'earth': SPHERICAL / 'earth_control_standard.npz'}
CALENDAR = ROOT / 'illumination/calendar/results/transfer.json'
SURFACE = ROOT / 'illumination/surface_light/results/surface_light.json'
WATERS = ROOT / 'biosphere/living_water/waters.json'
ATLAS = ROOT / 'geography/results/atlas.json'
DRAINAGE = ROOT / 'geography/products/drainage_28pct_16ppd.npz'
TOPOGRAPHY_INPUTS = [ROOT / 'geography/inputs' / n for n in ('ldem_16.img', 'ldem_4.img', 'jggrx_0420a_sha.tab')]
SEA_APPEARANCE = ROOT / 'research/studies/sea_appearance/results/water_colours.json'
BULK = Path('/media/projectspace/terluna-research/underwater-light')
OUT = HERE / 'results' / 'underwater_light.json'
SCHEMA = 'terluna.research.underwater-light/1'

DEPTH = np.unique(np.r_[np.arange(0.0, 20.0, 0.02), np.arange(20.0, 200.0, 0.2), np.arange(200.0, 2000.01, 2.0)])
ELEVATION = np.arange(-90.0, 90.01, 0.5)
STEPS = 1418                                  # half-hour steps through the 708.7-hour lunar day
PERIOD_H = SYNODIC_MONTH_DAYS * 24.0
LATITUDES = np.arange(0.0, 75.01, 5.0)        # reported places
AREA_LATITUDES = np.arange(0.0, 90.01, 2.0)   # for the floor, read per 16 px/deg row
BANDS = ((0, 15), (15, 30), (30, 45), (45, 60), (60, 75), (75, 90))
# The turned Moon: its north axis in today's body frame (x toward 0 N 0 E, z north), the third row of the elastic
# frame's rotation, the orientation set for the second set of simulations (geology/results/moon_orientation.json on
# domain/geology, 45c3d817; research/studies/moon_orientation on main). Its north pole lies at 81.95 N, 167.14 W.
TURNED_NORTH = np.array([-0.13658877239395323, -0.03117532763335884, 0.990137165347744])
TURNED_POLE_DEG = (81.94628788358153, -167.14292426123106)
LATITUDE_BINS = np.arange(-90.0, 90.01, 0.25)  # the floor by latitude in a frame, cells binned at their own latitudes
PPD = 16
DEPTH_EDGES = np.arange(0.0, 1000.01, 0.5)    # the floors' depth histograms (LOLA heights come in 0.5 m steps)
EVENING_DEG = (0.0, -6.0, -12.0, -18.0, -24.0, -30.0, -40.0)
PROFILE_M = (0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 50.0, 75.0, 100.0, 150.0, 200.0)
EARTH_ELEVATIONS = (90.0, 25.0)               # near side: under the sub-Earth point, and the development site's Earth
ROD_PEAK_NM = 480.0                           # a mesopelagic fish's rod pigment (assumption, 470-500 nm)

# Thresholds. Benthic minima are lunar-day means of photosynthetic photons at the floor, in mol per m2 per Earth day,
# Earth's reference points (sources.json); charophytes' are percentages of sub-surface light at their depth limits,
# converted with the Earth control's clear-sky equinoctial day-mean surface light at 40 degrees.
BENTHIC = {'seagrass': 5.1, 'coral_reef': 1.41, 'microphytobenthos': 0.4}
CHAROPHYTE_PERCENT = (3.0, 9.0)
CHAROPHYTE_LATITUDE = 40.0
EUPHOTIC_DOSE = 0.4                           # mol photons m-2 d-1, the microphytobenthos' minimum (Gattuso et al. 2006)
DARK_PAR = 1.0                                # umol m-2 s-1: "dark" for the hours count
COMPENSATION = {'microphytobenthos': 2.8, 'seagrass': 28.1}   # instantaneous community compensation, umol m-2 s-1
FEEDING_LUX = (0.01, 0.1, 1.0)                # visual feeding thresholds of fishes (Ono et al. 2010's review)
DEEP_LIMIT_M = 1000.0                         # daylight's limit for deep-sea fishes' eyes (Warrant and Locket 2004)
# A design water for the deep, nutrient-poor crater lakes (a design guess of this study, labelled as such).
CRATER_LAKE = dict(chlorophyll_mg_m3=dict(best=0.3, range=[0.1, 1.0]),
                   dissolved_organic_440_per_m=dict(best=0.05, range=[0.02, 0.3]),
                   fines_g_m3=dict(best=0.2, range=[0.0, 1.0]), soil='developed',
                   basis=('Design guess: a rain-fed crater lake on a small catchment, as poor in nutrients as the open '
                          'sea with a little more plankton, and more dissolved organic matter and fines from its '
                          'developed-soil shores. Between the open sea and the productive coast of waters.json.'))
# Earth's clearest open ocean, for the deep-sea threshold only: the open sea's lower bounds, at sea salinity.
CLEAREST_OCEAN = dict(chlorophyll_mg_m3=0.03, cdom_440=0.01, fines_g_m3=0.0, salinity_g_kg=35.0)
UMOL_DAY = 86400.0 / 1e6                      # umol m-2 s-1 to mol m-2 per Earth day


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def photons_per_watt(wavelength_nm):
    """umol photons s-1 per W at each wavelength."""
    return np.asarray(wavelength_nm, float) * 1e-9 / (PLANCK * SPEED_OF_LIGHT) / AVOGADRO * 1e6


# ==== Light on the water ===========================================================================================
def sky(world):
    """Direct and diffuse light on the water (W m-2 per channel) on the study's elevation grid, with the per-source
    response kept for other sources."""
    with np.load(SKY[world], allow_pickle=False) as z:
        a, sa = np.degrees(z['a']), np.degrees(z['sa'])
        diffuse, direct = z['moments'][0, :, :, 5].copy(), z['beam'][0, :, :, 1].copy()
        energy, xyz, wavelength, band = z['energy'].copy(), z['xyz'].copy(), z['wavelength'].copy(), z['band'].copy()
    on = lambda grid, table, e: np.array([np.interp(e, grid, table[:, c]) for c in range(table.shape[1])]).T
    photons = photons_per_watt(wavelength)
    weights = dict(par=photons * ((band[:, 0] >= 400) & (band[:, 1] <= 700)), lux=xyz[:, 1] / energy,
                   rod=photons * govardovskii(wavelength, ROD_PEAK_NM))
    return dict(world=world, wavelength=wavelength, band=band, energy=energy, weights=weights,
                direct=on(sa, direct, ELEVATION) * energy, diffuse=on(a, diffuse, ELEVATION) * energy,
                response=lambda e: (on(sa, direct, np.atleast_1d(e))[0], on(a, diffuse, np.atleast_1d(e))[0]))


def earthlight(light, elevation_deg, phase_deg=0.0):
    """The Earth's direct and diffuse light on the water (W m-2 per channel) at an elevation and phase: the calendar's
    Earth band irradiance shared over each band's channels as the shielded Sun's energy is (illumination/calendar)."""
    t = json.loads(CALENDAR.read_text())
    spectrum = t['earth']['spectrum']
    edges = np.asarray(spectrum['band_edges_nm'])
    phases = np.asarray(t['phase_deg'])
    earth = np.array([np.interp(phase_deg, phases, col) for col in np.asarray(spectrum['earth_irradiance_w_m2_by_phase']).T])
    index = np.array([int(np.flatnonzero((edges == b).all(1))[0]) for b in light['band']])
    shielded = np.bincount(index, light['energy'], len(edges))
    energy = light['energy'] * (earth / shielded)[index]
    direct, diffuse = light['response'](elevation_deg)
    calendar_lux = (np.interp(elevation_deg, t['direct_elevation_deg'], t['earth']['direct_lux'][0])
                    + np.interp(elevation_deg, t['diffuse_elevation_deg'], t['earth']['diffuse_lux'][0]))
    return direct * energy, diffuse * energy, float(calendar_lux)


# ==== Light down the column ========================================================================================
def table(light, column):
    """Photosynthetic photons (umol m-2 s-1), illuminance (lux) and rod catch (umol m-2 s-1, rod-weighted) at every
    elevation and depth: an array (elevations, 3, depths)."""
    w = np.stack([light['weights'][k] for k in ('par', 'lux', 'rod')], 1)
    out = np.tensordot(light['diffuse'][:, :, None] * w[None], column.diffuse, axes=([1], [0]))
    for k, e in enumerate(ELEVATION):
        if e > 0 and np.any(light['direct'][k] > 0):
            out[k] += np.tensordot(light['direct'][k][:, None] * w, column.direct(90.0 - e), axes=([0], [0]))
    return out


def at(tab, elevation_deg):
    """A table read at an elevation, linear between the grid's half degrees."""
    x = np.clip((np.asarray(elevation_deg, float) - ELEVATION[0]) / 0.5, 0, len(ELEVATION) - 1.000001)
    i = np.floor(x).astype(int)
    f = (x - i).reshape(x.shape + (1,) * (tab.ndim - 1))
    return (1 - f) * tab[i] + f * tab[i + 1]


def first_below(depth, values, level):
    """The first depth at which a profile falls to a level, log-linear between grid depths; nan if never."""
    v = np.asarray(values, float)
    if v[0] <= level:
        return 0.0
    k = np.flatnonzero(v <= level)
    if not len(k):
        return float('nan')
    j = k[0]
    y0, y1 = np.log(max(v[j - 1], 1e-300)), np.log(max(v[j], 1e-300))
    return float(depth[j - 1] + (depth[j] - depth[j - 1]) * (y0 - np.log(level)) / (y0 - y1))


def day(lat, declination=0.0):
    """Hour angle (degrees from noon), elevation and the chlorophyll phase (0 at dawn, 1 at dusk) through one day."""
    hour = (np.arange(STEPS) + 0.5) / STEPS * 360.0 - 180.0
    e = optics.sun_elevation(lat, hour, declination)
    c = np.nan_to_num(-np.tan(np.radians(lat)) * np.tan(np.radians(declination)))
    h0 = float(np.degrees(np.arccos(np.clip(c, -1, 1))))
    after = np.where(hour > h0, hour - h0, hour + 360.0 - h0)          # degrees since sunset
    f = np.where(np.abs(hour) <= h0, (hour + h0) / max(2 * h0, 1e-9), 1 - after / max(360.0 - 2 * h0, 1e-9))
    return hour, e, np.clip(f, 0, 1)


def through_day(dawn, dusk, lat, declination=0.0, quantity=0):
    """A quantity at every step and depth (steps, depths), the chlorophyll's phase mixed in the logarithm."""
    _, e, f = day(lat, declination)
    lo, hi = at(dawn[:, quantity], e), at(dusk[:, quantity], e)
    return np.exp((1 - f)[:, None] * np.log(np.maximum(lo, 1e-300)) + f[:, None] * np.log(np.maximum(hi, 1e-300)))


def twilight_share(surface, lat):
    """The share of a day's light on the water that arrives with the Sun below the horizon."""
    e = day(lat)[1]
    light = np.interp(e, ELEVATION, surface)
    return float(light[e <= 0].sum() / light.sum())


def hours_below(series, level):
    return float(np.sum(series < level) * PERIOD_H / STEPS)


# ==== The floors ===================================================================================================
def depth_tables():
    """Per 16 px/deg row, cumulative counts of sea and lake cells shallower than each edge, with the cell areas, and the
    same floors' areas binned at each cell's latitude today and on the turned Moon; cached on the drive under the
    inputs' hashes."""
    identity = {str(p.relative_to(ROOT)): sha256(p) for p in (ATLAS, DRAINAGE, ROOT / 'geography/topography.py',
                                                              *TOPOGRAPHY_INPUTS)}
    layout = 'cumulative counts per row with totals, and areas by latitude today and on the turned Moon, 2026-11-04'
    cache = BULK / 'depths_16ppd.npz'
    if cache.exists():
        with np.load(cache, allow_pickle=False) as z:
            if ('layout' in z.files and str(z['layout']) == layout and json.loads(str(z['identity'])) == identity
                    and np.array_equal(z['edges'], DEPTH_EDGES)):
                return {k: z[k].copy() for k in z.files if k not in ('identity', 'layout')} | dict(identity=identity)
    from geography import topography as tp
    atlas = json.loads(ATLAS.read_text())
    level = float(atlas['sea_level_m'])
    h, _, grid = tp.height_above_geoid(PPD, atlas['grid']['geoid_degree'])
    rows, cols = h.shape
    area = grid.cell_area_m2[:, 0] / 1e6
    sea = np.zeros((rows, len(DEPTH_EDGES)), np.int64)
    sea_all = (h < level).sum(1)
    for r in range(rows):
        d = level - h[r][h[r] < level]
        sea[r] = np.searchsorted(np.sort(d), DEPTH_EDGES, side='left')
    d = np.load(DRAINAGE, allow_pickle=False)
    cells = d['lake_cell'].astype(np.int64)
    lake_depth = d['lake_level_m'].astype(float) - h.ravel()[cells]
    lake_row = cells // cols
    lake = np.zeros_like(sea)
    order = np.argsort(lake_row, kind='stable')
    starts = np.searchsorted(lake_row[order], np.arange(rows + 1))
    for r in range(rows):
        lake[r] = np.searchsorted(np.sort(lake_depth[order[starts[r]:starts[r + 1]]]), DEPTH_EDGES, side='left')
    out = dict(lat=grid.lat_deg, area_km2=area, sea=sea, lake=lake, edges=DEPTH_EDGES, sea_all=sea_all,
               lake_all=np.bincount(lake_row, minlength=rows),
               sea_cells=np.array([int((h < level).sum())]), lake_negative=np.array([int((lake_depth < 0).sum())]),
               sea_level_m=np.array([level]), bin_lat=0.5 * (LATITUDE_BINS[1:] + LATITUDE_BINS[:-1]))
    # The same floors binned at each cell's own latitude: today's, which checks the binning against the rows, and the
    # turned Moon's. Lake cells below their floor count from the first edge on.
    sea_row, sea_col = np.nonzero(h < level)
    places = dict(sea=(sea_row, sea_col, level - h[sea_row, sea_col]),
                  lake=(lake_row, cells % cols, np.maximum(lake_depth, 0.0)))
    for frame, north in (('today', np.array([0.0, 0.0, 1.0])), ('turned', TURNED_NORTH)):
        for kind, (r, c, depth) in places.items():
            lat = turned_latitude(grid.lat_deg[r], grid.lon_deg[c], north)
            out[f'{kind}_{frame}'], out[f'{kind}_{frame}_all'] = binned_areas(lat, depth, area[r])
    BULK.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, identity=json.dumps(identity), layout=layout, **out)
    return out | dict(identity=identity)


def turned_latitude(lat_deg, lon_deg, north=TURNED_NORTH):
    """Latitude (degrees) in the frame whose north axis is `north`, a unit vector in today's body frame."""
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    s = np.cos(lat) * (north[0] * np.cos(lon) + north[1] * np.sin(lon)) + north[2] * np.sin(lat)
    return np.degrees(np.arcsin(np.clip(s, -1.0, 1.0)))


def binned_areas(lat_deg, depth_m, area_km2):
    """km2 of floor shallower than each depth edge in each latitude bin, and each bin's whole floor."""
    hist, _, _ = np.histogram2d(lat_deg, depth_m, bins=[LATITUDE_BINS, DEPTH_EDGES], weights=area_km2)
    shallower = np.concatenate([np.zeros((hist.shape[0], 1)), np.cumsum(hist, axis=1)], axis=1)
    return shallower, np.histogram(lat_deg, bins=LATITUDE_BINS, weights=area_km2)[0]


def by_band(lat, area):
    out = {f'{lo}-{hi}': float(area[(lat >= lo) & (lat < hi + (hi == 90))].sum()) for lo, hi in BANDS}
    out['all'] = float(area.sum())
    return out


def lit_area(tables, kind, limit_by_lat, frame=None):
    """km2 of floor shallower than the limit depth at each row's latitude, by band; limit_by_lat(|lat|) in m. With a
    frame ('today' or 'turned'), the floor binned at each cell's latitude in that frame."""
    if frame is not None:
        lat = np.abs(tables['bin_lat'])
        limit = np.clip(limit_by_lat(lat), 0, DEPTH_EDGES[-1])
        return by_band(lat, np.array([np.interp(z, DEPTH_EDGES, a) for z, a in zip(limit, tables[f'{kind}_{frame}'])]))
    lat = np.abs(tables['lat'])
    limit = np.clip(limit_by_lat(lat), 0, DEPTH_EDGES[-1])
    counts = np.array([np.interp(z, DEPTH_EDGES, c) for z, c in zip(limit, tables[kind])])
    return by_band(lat, counts * tables['area_km2'])


def floor_by_band(tables, kind, depth_m, frame=None):
    """km2 of floor shallower than depth_m (all of it if None) by band, per row or binned in a frame."""
    if frame is not None:
        lat = np.abs(tables['bin_lat'])
        return by_band(lat, tables[f'{kind}_{frame}_all'] if depth_m is None else
                       np.array([np.interp(depth_m, DEPTH_EDGES, a) for a in tables[f'{kind}_{frame}']]))
    lat = np.abs(tables['lat'])
    counts = (tables[f'{kind}_all'] if depth_m is None else
              np.array([np.interp(depth_m, DEPTH_EDGES, c) for c in tables[kind]])) * tables['area_km2']
    return by_band(lat, counts)


# ==== Ultraviolet ==================================================================================================
def ultraviolet(product, case, water, depth, top_nm=400.0):
    """UV-A, UV-B and DNA-weighted irradiance (W m-2) against depth at the product's solar elevations; top_nm cuts
    the band short, to isolate the action spectrum's measured part."""
    c = product['cases'][case]
    w = np.asarray(product['centres_nm'])
    keep = (w >= 280) & (w < top_nm)
    r = np.asarray(product['columns'][c['column']]['reflectance_from_below'])
    column = optics.Column(w[keep], depth, **water)
    weights = np.stack([(w[keep] >= 315).astype(float), (w[keep] < 315).astype(float), optics.dna_damage(w[keep])], 1)
    suns = sorted(c['spectra'], key=float)
    rows = []
    for s in suns:
        spec = c['spectra'][s]
        direct = np.asarray(spec['direct_w_m2_nm'])
        total = (direct + np.asarray(spec['diffuse_black_w_m2_nm'])) / (1 - 0.1 * r)
        below = column.below(direct[keep], (total - direct)[keep], 90.0 - float(s))
        rows.append(weights.T @ below)
    return np.array([float(s) for s in suns]), np.array(rows)          # (suns), (suns, 3, depths)


def uv_day_mean(suns, rows, lat):
    """The lunar-day mean of the ultraviolet rows: linear in elevation, falling to nothing from 1 degree to 0."""
    _, e, _ = day(lat)
    nodes = np.r_[0.0, suns]
    values = np.concatenate([np.zeros((1,) + rows.shape[1:]), rows])
    x = np.clip(e, 0, 90)
    j = np.clip(np.searchsorted(nodes, x) - 1, 0, len(nodes) - 2)
    f = ((x - nodes[j]) / (nodes[j + 1] - nodes[j]))[:, None, None]
    series = (1 - f) * values[j] + f * values[j + 1]
    return np.where((e > 0)[:, None, None], series, 0.0).mean(0)


# ==== The study ====================================================================================================
def water_cases(record):
    ratio = record['light_cycle']['chlorophyll_dusk_over_dawn']['best']
    waters = dict(record['waters']) | {'crater_lake': CRATER_LAKE}
    out = {}
    for name, wtr in waters.items():
        chl = wtr['chlorophyll_mg_m3']['best']
        base = dict(cdom_440=wtr['dissolved_organic_440_per_m']['best'], fines_g_m3=wtr['fines_g_m3']['best'],
                    soil=wtr['soil'])
        out[name] = {phase: base | dict(chlorophyll_mg_m3=chl * s) for phase, s in
                     (('dawn', 2 / (1 + ratio)), ('mean', 1.0), ('dusk', 2 * ratio / (1 + ratio)))}
    return out, ratio


def rounded(x, digits=4):
    if isinstance(x, dict):
        return {k: rounded(v, digits) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [rounded(v, digits) for v in x]
    if isinstance(x, (float, np.floating)):
        return None if not np.isfinite(x) else float(f'{float(x):.{digits}g}')
    return x


def build():
    start = time.time()
    record = json.loads(WATERS.read_text())
    cases, ratio = water_cases(record)
    lights = {world: sky(world) for world in SKY}
    calendar = json.loads(CALENDAR.read_text())
    tables = depth_tables()
    product = json.loads(SURFACE.read_text())
    surface_par = {}
    checks = {}

    # The twilight's diffuse light against the calendar's converged point-source solve (lux).
    moon = lights['moon']
    cal = np.interp(ELEVATION, calendar['diffuse_elevation_deg'], calendar['solar']['diffuse_lux'])
    mine = moon['diffuse'] @ moon['weights']['lux']
    band = (ELEVATION >= -40) & (ELEVATION <= 90)
    checks['twilight_diffuse_vs_calendar_max_relative'] = float(np.max(np.abs(mine[band] / cal[band] - 1)))

    # The deep-sea fishes' daylight limit: the rod catch at 1000 m in Earth's clearest ocean, overhead Sun.
    earth = lights['earth']
    clear = optics.Column(earth['wavelength'], DEPTH, **CLEAREST_OCEAN)
    k = int(np.flatnonzero(ELEVATION == 90.0)[0])
    rod = clear.below(earth['direct'][k], earth['diffuse'][k], 0.0).T @ earth['weights']['rod']
    deep_threshold = float(np.interp(DEEP_LIMIT_M, DEPTH, rod))

    # Charophytes: percentages of the Earth control's clear-sky day-mean surface light at 40 degrees.
    _, e40, _ = day(CHAROPHYTE_LATITUDE)
    earth_surface = (earth['direct'] + earth['diffuse']) @ earth['weights']['par']
    earth_mean_40 = float(np.interp(e40, ELEVATION, earth_surface).mean())
    thresholds = dict(BENTHIC) | {f'charophyte_{p:g}pct': p / 100 * earth_mean_40 * UMOL_DAY for p in CHAROPHYTE_PERCENT}

    zones, floors, dark, spectra, eyes, seasonal, profiles, twilight = {}, {}, {}, {}, {}, {}, {}, {}
    floors_turned, binning, bulk = {}, [], {}
    for world, light in lights.items():
        surface = (light['direct'] + light['diffuse']) @ light['weights']['par']
        surface_par[world] = {f'{lat:g}': dict(noon=float(np.interp(90 - lat, ELEVATION, surface)),
                                               day_mean=float(np.interp(day(lat)[1], ELEVATION, surface).mean()))
                              for lat in LATITUDES}
        for p in surface_par[world].values():
            p['day_mean_mol_m2_d'] = p['day_mean'] * UMOL_DAY
        zones[world], spectra[world], profiles[world] = {}, {}, {}
        twilight[world] = {f'{lat:g}': twilight_share(surface, lat) for lat in (0.0, 45.0, 75.0)}
        for name, phases in cases.items():
            tabs = {ph: table(light, optics.Column(light['wavelength'], DEPTH, **kw)) for ph, kw in phases.items()}
            bulk[f'{world}_{name}_mean'] = tabs['mean'].astype(np.float32)
            rows = {}
            for lat in LATITUDES:
                noon = at(tabs['mean'], 90 - lat)
                par = through_day(tabs['dawn'], tabs['dusk'], lat)
                mean = par.mean(0)
                z_dose = first_below(DEPTH, mean * UMOL_DAY, EUPHOTIC_DOSE)
                pick = [int(np.argmin(np.abs(DEPTH - z))) for z in PROFILE_M]
                profiles[world].setdefault(name, {})[f'{lat:g}'] = dict(noon=noon[0][pick].tolist(), day_mean=mean[pick].tolist())
                rows[f'{lat:g}'] = dict(
                    euphotic_noon_1pct_m=first_below(DEPTH, noon[0] / noon[0][0], 0.01),
                    euphotic_noon_1pct_dawn_m=first_below(DEPTH, at(tabs['dawn'], 90 - lat)[0] / at(tabs['dawn'], 90 - lat)[0][0], 0.01),
                    euphotic_noon_1pct_dusk_m=first_below(DEPTH, at(tabs['dusk'], 90 - lat)[0] / at(tabs['dusk'], 90 - lat)[0][0], 0.01),
                    euphotic_dose_m=z_dose,
                    surface_day_mean_par_below=float(mean[0]),
                    twilight_feeding_end_m={f'{x:g}_lux': first_below(DEPTH, noon[1], x) for x in FEEDING_LUX},
                    disphotic_end_m=first_below(DEPTH, noon[2], deep_threshold),
                    benthic_limit_m={t: first_below(DEPTH, mean * UMOL_DAY, v) for t, v in thresholds.items()})
                if world == 'moon' and name in ('open_sea', 'productive_coast', 'crater_lake'):
                    i = np.searchsorted(DEPTH, z_dose) if np.isfinite(z_dose) else 0
                    dark.setdefault(name, {})[f'{lat:g}'] = dict(
                        hours_sun_below_horizon=float(np.sum(day(lat)[1] <= 0) * PERIOD_H / STEPS),
                        hours_surface_below_1umol=hours_below(par[:, 0], DARK_PAR),
                        hours_below_microphytobenthos_compensation_at_dose_depth=hours_below(
                            par[:, i], COMPENSATION['microphytobenthos']))
            zones[world][name] = rows
            if world == 'moon':
                means = [through_day(tabs['dawn'], tabs['dusk'], lat).mean(0) * UMOL_DAY for lat in AREA_LATITUDES]
                lim = {t: np.array([first_below(DEPTH, m, v) for m in means]) for t, v in thresholds.items()}
                kinds = ('lake',) if name == 'crater_lake' else ('sea', 'lake') if name == 'productive_coast' else ('sea',)
                for kind in kinds:
                    limits = {t: (lambda lat, z=np.nan_to_num(v, nan=DEPTH_EDGES[-1]): np.interp(lat, AREA_LATITUDES, z))
                              for t, v in lim.items()}
                    floors.setdefault(kind, {})[name] = {t: lit_area(tables, kind, f) for t, f in limits.items()}
                    floors[kind][name]['limit_depth_m_by_latitude'] = {t: dict(zip([f'{x:g}' for x in AREA_LATITUDES[::5]],
                                                                                   v[::5].tolist())) for t, v in lim.items()}
                    floors_turned.setdefault(kind, {})[name] = {t: lit_area(tables, kind, f, 'turned')
                                                                for t, f in limits.items()}
                    # Binned at today's latitudes, the floor must repeat the rows' areas.
                    for t, f in limits.items():
                        rows_, binned = floors[kind][name][t], lit_area(tables, kind, f, 'today')
                        binning.extend(abs(binned[b] / rows_[b] - 1) for b in rows_ if rows_[b] > 1000.0)
                if name in ('open_sea', 'productive_coast', 'crater_lake'):
                    seasonal[name] = {f'{d:+.2f}': first_below(DEPTH, through_day(tabs['dawn'], tabs['dusk'], 75.0, d).mean(0)
                                                                * UMOL_DAY, EUPHOTIC_DOSE)
                                      for d in (-MOON_EQUATOR_TO_ECLIPTIC_DEG, 0.0, MOON_EQUATOR_TO_ECLIPTIC_DEG)}
            # Spectral downwelling irradiance at noon, by 10 nm band, at chosen depths (W m-2 per band).
            column = optics.Column(light['wavelength'], DEPTH, **phases['mean'])
            edges = np.unique(light['band'], axis=0)
            index = np.array([int(np.flatnonzero((edges == b).all(1))[0]) for b in light['band']])
            pick = [int(np.argmin(np.abs(DEPTH - z))) for z in (0, 1, 2, 5, 10, 20, 50, 100, 200)]
            spectra[world][name] = {}
            for lat in (0.0, 45.0):
                k = int(np.flatnonzero(ELEVATION == 90 - lat)[0])
                below = column.below(light['direct'][k], light['diffuse'][k], lat)[:, pick]
                spectra[world][name][f'{lat:g}'] = {f'{DEPTH[j]:g}': np.bincount(index, below[:, n], len(edges)).tolist()
                                                     for n, j in enumerate(pick)}
            spectra[world]['band_edges_nm'] = edges.tolist()
            # Eyes: the evening's light and the near side's full Earth.
            ev = {}
            for e in EVENING_DEG:
                prof = at(tabs['mean'], e)
                hours = {f'{lat:g}': float(np.degrees(np.arccos(np.clip(np.sin(np.radians(e)) / np.cos(np.radians(lat)), -1, 1)))
                                           - 90.0) * PERIOD_H / 360.0 for lat in (0.0, 45.0)}
                ev[f'{e:g}'] = dict(surface_lux=float(((light['direct'] + light['diffuse']) @ light['weights']['lux'])[
                                        int(np.flatnonzero(ELEVATION == e)[0])]),
                                    hours_after_sunset=hours,
                                    feeding_m={f'{x:g}_lux': first_below(DEPTH, prof[1], x) for x in FEEDING_LUX},
                                    deep_limit_m=first_below(DEPTH, prof[2], deep_threshold))
            eyes.setdefault(world, {})[name] = dict(evening=ev)
            if world == 'moon':
                full = {}
                for ee in EARTH_ELEVATIONS:
                    d_e, f_e, cal_lux = earthlight(light, ee)
                    below = column.below(d_e, f_e, 90.0 - ee)
                    lux, rodc = below.T @ light['weights']['lux'], below.T @ light['weights']['rod']
                    full[f'{ee:g}'] = dict(surface_lux=float((d_e + f_e) @ light['weights']['lux']),
                                           calendar_lux=cal_lux,
                                           feeding_m={f'{x:g}_lux': first_below(DEPTH, lux, x) for x in FEEDING_LUX},
                                           deep_limit_m=first_below(DEPTH, rodc, deep_threshold))
                eyes[world][name]['full_earth'] = full
    # The check against the sea-appearance study: the 1% depth of photosynthetic photons, Sun 45 degrees up.
    old = json.loads(SEA_APPEARANCE.read_text())
    ref = {c['water']: c['moon']['depth_one_percent_photons_m'] for c in old['cases']
           if c['phase'] == 'mean' and c['soil'] == 'developed'}
    checks['sea_appearance_45deg'] = {}
    for name in ref:
        mean = cases[name]['mean']
        light = lights['moon']
        k = int(np.flatnonzero(ELEVATION == 45.0)[0])
        variants = {'this_study': optics.Column(light['wavelength'], DEPTH, **mean),
                    'seawater_scattering': optics.Column(light['wavelength'], DEPTH, salinity_g_kg=38.4, **mean)}
        row = dict(sea_appearance=ref[name])
        for label, col in variants.items():
            below = col.below(light['direct'][k], light['diffuse'][k], 45.0).T @ light['weights']['par']
            row[label] = first_below(DEPTH, below / below[0], 0.01)
        # Their method: all light through K_d at 45 degrees, no surface, their 400-700 nm channel rule.
        iop = domain.water(light['wavelength'], mean['chlorophyll_mg_m3'], mean['cdom_440'], mean['fines_g_m3'],
                           mean['soil'])
        kd = optics.kd(iop, 45.0)
        inside = (light['wavelength'] >= 400) & (light['wavelength'] <= 700)
        ph = ((light['direct'][k] + light['diffuse'][k]) * photons_per_watt(light['wavelength']))[inside]
        prof = (ph[None, :] * np.exp(-np.outer(DEPTH, kd[inside]))).sum(1)
        row['their_method_reproduced'] = first_below(DEPTH, prof / prof[0], 0.01)
        checks['sea_appearance_45deg'][name] = row

    # Ultraviolet at depth, the Moon's shielded sky and Earth's, mean chlorophyll.
    uv = {}
    uv_depth = np.unique(np.r_[np.arange(0.0, 20.0, 0.02), np.arange(20.0, 200.0, 0.2)])
    for world, case in (('moon', 'moon_1.2atm/design'), ('earth', 'earth_control/unfiltered')):
        uv[world] = {}
        for name in ('open_sea', 'productive_coast', 'river_mouth', 'crater_lake'):
            variants = {'central': {}, 'no_pure_water_uv_absorption': dict(uv_absorption=False),
                        'cdom_uv_slope_0.025': dict(cdom_uv_slope=0.025)} if name == 'open_sea' else {'central': {}}
            uv[world][name] = {}
            for label, extra in variants.items():
                suns, rows = ultraviolet(product, case, cases[name]['mean'] | extra, uv_depth)
                suns340, rows340 = ultraviolet(product, case, cases[name]['mean'] | extra, uv_depth[:1], top_nm=340.0)
                dna340 = {f'{lat:g}': np.interp(90 - lat, suns340, rows340[:, 2, 0]) for lat in (0.0, 45.0)}
                out = {}
                for lat in (0.0, 45.0):
                    mean = uv_day_mean(suns, rows, lat)
                    noon = np.array([np.interp(90 - lat, suns, rows[:, q, j]) for q in range(3)
                                     for j in range(len(uv_depth))]).reshape(3, -1)
                    out[f'{lat:g}'] = dict(
                        noon_surface_below_w_m2=dict(uva=float(noon[0, 0]), uvb=float(noon[1, 0]), dna=float(noon[2, 0])),
                        day_mean_surface_below_w_m2=dict(uva=float(mean[0, 0]), uvb=float(mean[1, 0]), dna=float(mean[2, 0])),
                        dna_day_mean_10pct_m=first_below(uv_depth, mean[2] / mean[2, 0], 0.1),
                        dna_day_mean_1pct_m=first_below(uv_depth, mean[2] / mean[2, 0], 0.01),
                        uva_day_mean_10pct_m=first_below(uv_depth, mean[0] / mean[0, 0], 0.1),
                        uva_day_mean_1pct_m=first_below(uv_depth, mean[0] / mean[0, 0], 0.01),
                        uvb_day_mean_1pct_m=(first_below(uv_depth, mean[1] / mean[1, 0], 0.01) if mean[1, 0] > 0
                                             else float('nan')),
                        dna_below_340nm_noon_surface_below_w_m2=float(dna340[f'{lat:g}']),
                        profile_day_mean_w_m2={f'{z:g}': [float(mean[q, int(np.argmin(np.abs(uv_depth - z)))])
                                                          for q in range(3)] for z in (0, 1, 2, 5, 10, 20, 50)})
                uv[world][name][label] = out
    # The DNA weighting's spectral make-up just below the surface, Sun overhead.
    w = np.asarray(product['centres_nm'])
    share = {}
    for world, case in (('moon', 'moon_1.2atm/design'), ('earth', 'earth_control/unfiltered')):
        spec = product['cases'][case]['spectra']['90.0']
        r = np.asarray(product['columns'][product['cases'][case]['column']]['reflectance_from_below'])
        total = (np.asarray(spec['direct_w_m2_nm']) + np.asarray(spec['diffuse_black_w_m2_nm'])) / (1 - 0.1 * r)
        dna = total * optics.dna_damage(w) * (w >= 280) * (w < 400)
        share[world] = {f'{lo}-{hi}': float(dna[(w >= lo) & (w < hi)].sum() / dna.sum())
                        for lo, hi in ((280, 315), (315, 340), (340, 400))}
    # The twilight's share of the day's 360-400 nm light on the water (spherical sky).
    uva_channels = (moon['band'][:, 1] <= 400)
    e0 = day(0.0)[1]
    uva = np.interp(e0, ELEVATION, (moon['direct'] + moon['diffuse'])[:, uva_channels].sum(1))
    checks['twilight_share_of_360_400nm_day_dose_equator'] = float(uva[e0 < 1].sum() / uva.sum())

    BULK.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(BULK / 'light_tables.npz', depth_m=DEPTH, elevation_deg=ELEVATION,
                        quantities=np.array(['par_umol_m2_s', 'lux', 'rod_umol_m2_s']), **bulk)
    shallow = {kind: floor_by_band(tables, kind, 30.0) for kind in ('sea', 'lake')}
    totals = {kind: floor_by_band(tables, kind, None) for kind in ('sea', 'lake')}
    shallow_turned = {kind: floor_by_band(tables, kind, 30.0, 'turned') for kind in ('sea', 'lake')}
    totals_turned = {kind: floor_by_band(tables, kind, None, 'turned') for kind in ('sea', 'lake')}
    checks['floor_binned_today_vs_rows_max_relative'] = float(max(binning))
    files = [Path(__file__), HERE / 'optics.py', ROOT / 'illumination/water_column/model.py',
             ROOT / 'illumination/water_column/fresh_water.py',
             ROOT / 'biosphere/ecology/screen.py', ROOT / 'geography/topography.py']
    inputs = {str(p.relative_to(ROOT)): sha256(p) for p in (*SKY.values(), CALENDAR, SURFACE, WATERS, SEA_APPEARANCE)}
    inputs |= tables['identity']
    inputs |= {f"illumination/water_column/inputs/{f['name']}": f['sha256'] for f in
               json.loads((ROOT / 'illumination/water_column/inputs.json').read_text())['files']}
    return dict(
        schema=SCHEMA,
        producer=dict(lane='research', runner=str(Path(__file__).relative_to(ROOT)),
                      files={str(p.relative_to(ROOT)): sha256(p) for p in files}, constants=constants_used(files)),
        inputs=inputs,
        evidence=('A calculation on committed products: the solved spherical sky (clear, ground albedo 0.1), the '
                  "calendar's Earth spectrum, the surface-light product's ultraviolet, the biosphere's design guesses "
                  "for the waters and this study's design guess for the crater lakes, through Earth's ocean optics "
                  "with fresh water's scattering and ultraviolet absorption; geography's 16 px/deg depths. Earth's "
                  'light requirements are reference points, and organisms are engineered. No organism, survival '
                  'through the night, wave focusing, inelastic scattering or real reef shape is modelled.'),
        reading_rule=('Depths in metres below a flat surface; light just below the surface is depth 0. PAR is 400-700 '
                      'nm photons in umol m-2 s-1; doses in mol m-2 per Earth day (umol m-2 s-1 x 0.0864). zones'
                      '[world][water][latitude]: euphotic_noon_1pct_m (1% of the PAR just below the surface at local '
                      'noon, mean chlorophyll; dawn and dusk beside it), euphotic_dose_m (lunar-day mean PAR falls to '
                      f'{EUPHOTIC_DOSE} mol m-2 d-1), twilight_feeding_end_m (noon illuminance falls to the fish feeding '
                      'thresholds), disphotic_end_m (noon rod catch falls to the deep-sea limit; aphotic below), '
                      'benthic_limit_m (lunar-day mean PAR at the floor falls to each threshold). floor[sea|lake][water]'
                      '[threshold][band] is km2 of floor lit above the threshold, bands of |latitude| in degrees, today\'s '
                      'latitudes; floor_turned is the same floor with each cell at its latitude on the turned Moon '
                      '(turned_frame), today\'s depths and lakes. For other places on the turned Moon read each at its '
                      'new latitude. Declination 0 throughout except seasonal_75.'),
        units=dict(depth='m', par='umol m-2 s-1', dose='mol m-2 d-1 (Earth day)', illuminance='lux',
                   rod='umol m-2 s-1 weighted by an A1 template at 480 nm', ultraviolet='W m-2',
                   dna='W m-2 weighted by the DNA-damage action spectrum, 1 at 300 nm', area='km2'),
        earth_terms_fixed=dict(
            pure_water_scattering=f'scaled from 38.4 to {optics.SEA_SALINITY} g/kg: x{optics.scattering_factor(optics.SEA_SALINITY):.4f}',
            pure_water_absorption='IOCCG (2018) Table 1.1 composite below 420 nm in place of Smith and Baker (1981) '
                                  'below 380 nm and Pope and Fry 380-415 nm',
            refractive_index=optics.WATER_INDEX,
            dissolved_organic_uv='domain slope 0.0176 per nm carried into the ultraviolet; 0.025 below 350 nm as a '
                                 'sensitivity'),
        settings=dict(depth_grid_m='0-20 by 0.02, 20-200 by 0.2, 200-2000 by 2', elevation_step_deg=0.5,
                      steps_per_lunar_day=STEPS, lunar_day_h=PERIOD_H, sky_streams=8, rod_peak_nm=ROD_PEAK_NM,
                      chlorophyll_dusk_over_dawn=ratio, declination_deg=0.0,
                      seasonal_declination_deg=MOON_EQUATOR_TO_ECLIPTIC_DEG),
        thresholds=dict(benthic_mol_m2_d=thresholds, euphotic_dose_mol_m2_d=EUPHOTIC_DOSE,
                        earth_clear_sky_day_mean_par_40deg_umol_m2_s=earth_mean_40,
                        feeding_lux=FEEDING_LUX, deep_limit=dict(depth_in_clearest_ocean_m=DEEP_LIMIT_M,
                                                                 water=CLEAREST_OCEAN, rod_umol_m2_s=deep_threshold),
                        dark_par_umol_m2_s=DARK_PAR, compensation_umol_m2_s=COMPENSATION),
        waters={name: dict(phases={ph: {k: v for k, v in kw.items()} for ph, kw in ph_.items()},
                           status='design guess of this study' if name == 'crater_lake' else 'biosphere/living_water/waters.json')
                for name, ph_ in cases.items()} | {'crater_lake_record': CRATER_LAKE},
        surface_par=surface_par, zones=zones, floor=floors, floor_shallower_than_30m_km2=shallow, floor_total_km2=totals,
        floor_turned=floors_turned, floor_shallower_than_30m_turned_km2=shallow_turned,
        floor_total_turned_km2=totals_turned,
        turned_frame=dict(north_axis_in_todays_frame=TURNED_NORTH.tolist(), north_pole_lat_lon_deg=list(TURNED_POLE_DEG),
                          latitude_bin_deg=float(LATITUDE_BINS[1] - LATITUDE_BINS[0]),
                          source=('the elastic frame with the seas and lakes in place, set for the second set of '
                                  'simulations: geology/results/moon_orientation.json on domain/geology (45c3d817), '
                                  'research/studies/moon_orientation on main')),
        depth_table_notes=dict(sea_cells=int(tables['sea_cells'][0]), lake_cells_below_their_floor=int(tables['lake_negative'][0]),
                               sea_level_m=float(tables['sea_level_m'][0])),
        dark_hours=dark, seasonal_75=seasonal, ultraviolet=uv, dna_weighting_share_overhead_surface=share,
        par_profiles_umol_m2_s=dict(depth_m=list(PROFILE_M), by_world=profiles),
        twilight_share_of_day_par=twilight,
        eyes=eyes, spectra_noon_w_m2_per_band=spectra, checks=checks,
        elapsed_s=time.time() - start)


def main():
    # Four significant figures for the results; the provenance keeps its values exact.
    product = {k: v if k in ('producer', 'inputs') else rounded(v) for k, v in build().items()}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1, allow_nan=False) + '\n')
    print(json.dumps(product['checks'], indent=1))
    for world in product['zones']:
        for name, rows in product['zones'][world].items():
            r = rows['0']
            print(world, name, 'noon 1%', r['euphotic_noon_1pct_m'], 'dose', r['euphotic_dose_m'],
                  'disphotic', r['disphotic_end_m'])


if __name__ == '__main__':
    main()
