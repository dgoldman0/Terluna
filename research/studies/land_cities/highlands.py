"""Highland and tall cities: the high ground on the 28% atlas, its air, rain and lightning, how tall buildings get
at a sixth of Earth's gravity, the land and water the highland people take, and how they get about.

Grid sums on the atlas, closed-form physics and the engineering lattice model (engineering/towers/lattice.py) run on
building-sized towers. No climate, flow or structural analysis runs here.
"""
from __future__ import annotations

import math

import numpy as np

from engineering.towers import lattice as lt
from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY, SYNODIC_MONTH_DAYS

from research.studies.land_cities import grid

HIGHLAND_KM = 4.0                       # the ways-of-living study's highlands: ground 4 km and more above the sea
HEIGHT_BANDS_KM = ((4.0, 5.0), (5.0, 6.0), (6.0, 7.0), (7.0, 8.0), (8.0, 9.0), (9.0, 10.0), (10.0, 13.0))
LOW_BANDS_KM = ((0.0, 1.0), (1.0, 2.0), (2.0, 3.0), (3.0, 4.0))
LAT_BANDS_DEG = ((0.0, 15.0), (15.0, 30.0), (30.0, 45.0), (45.0, 60.0), (60.0, 90.01))
COMFORT_POLEWARD_DEG = 30.0             # both climate models put the most comfortable land poleward of 30 degrees
REGISTER_LAPSE_K_PER_KM = 1.0           # register: the highlands' air about 1 C cooler per km (CM1, corrected)
SUMMIT_RADII_KM = (32.0, 50.0, 100.0, 200.0, 300.0)
LUNAR_DAY_H = SYNODIC_MONTH_DAYS * 24.0
WATER_DENSITY = 1000.0                  # kg/m3, fresh water (the seas and lakes are fresh)
SECONDS_PER_YEAR = JULIAN_YEAR_DAYS * JULIAN_DAY

# Building towers for the lattice model: a square floor plate on every storey, the frame inside a closed facade.
STOREY_M = 4.0                          # floor to floor
PLATE_M = 50.0                          # side of the floor plate, 2,500 m2 a storey
STOREYS = (30, 60, 100, 200, 400, 600)
CLOSED_FACADE_SOLIDITY = 1e-3           # the frame stands inside the facade; the facade takes the wind
FLOOR_MASS_KG_M2 = 500.0                # the lattice model's floors: structure, fit-out, facade share and occupants
BASE_SHARES = np.geomspace(1.0 / 40.0, 0.5, 40)
EXPONENTS = (1.0, 1.5, 2.0, 3.0)
SWAY_LIMIT = 1.0 / 500.0
FLUCTUATING_SWAY = 0.3                  # the summit tower study's gust share of the steady sway, for acceleration
TOWER_GROUND_KM = 6.0                   # towers stand on highland ground 6 km above the sea (air from the design run)
EARTH_DESIGN_GUST = 50.0                # the summit tower study's Earth-like design gust
EARTH_SERVICE_SHARE = 22.5 / 32.4       # Earth's service gust in the same proportion to its design gust as the Moon's
# Lateral-system steel of designed tube, diagrid and hexagrid towers on Earth, kg per m2 of floor, by storeys
# (Kheyroddin et al. 2018: 9-15 psf at 48 storeys and 20-44 psf at 96).
DESIGNED_LATERAL_KG_M2 = {48: (44.0, 72.0), 96: (98.0, 215.0)}
DESIGNED_LATERAL_SOURCE = 'kheyroddin_2018'
COVERAGES = (0.1, 0.2, 0.3)             # share of a district's ground under towers (the tower district: a fifth)
FLOOR_M2_PER_PERSON = 70.0              # the metropolis brief's floor for home, work and services
# Densities for the people on the ground (people per km2 of settled land).
DENSITIES = {'towns': 10_000.0, 'metropolis': 40_000.0, 'tower_district': 85_700.0}
GROUND_BILLIONS = (2.0, 4.0, 6.0)
HIGHLAND_SHARES = (0.5, 0.75, 1.0)      # share of the ground's people who live in the highlands
KEPT_SHARES = (0.3, 0.5)                # wild land kept in every home (ways of living, after ecology)
RECOVERY = (0.0, 0.5, 0.9, 0.99)        # share of domestic water recovered and reused in the city
TRIP_KM = (50.0, 100.0, 300.0, 1000.0)  # trips from a highland city down to the coast and across the Moon
CLIMB_M = 6000.0                        # the climb from the coast to a highland city, for the sky boat


def r(x, n=3):
    if x is None or isinstance(x, (bool, int, str)):
        return x
    x = float(x)
    if x == 0.0 or not math.isfinite(x):
        return x
    return float(f'{x:.{n}g}')


def band_mask(at: grid.Atlas, lo_km, hi_km):
    return (at.hsl_km >= lo_km) & (at.hsl_km < hi_km)


def lat_mask(at: grid.Atlas, lo_deg, hi_deg):
    a = np.abs(at.LA)
    return (a >= lo_deg) & (a < hi_deg)


def area_mean(field, weight, mask):
    w = weight * mask
    s = float(w.sum())
    return float((field * w).sum() / s) if s > 0 else None


# ---------------------------------------------------------------------------------------------------------------
# The high ground.

def land(at: grid.Atlas, rain, evap, flashes, home):
    """Dry land by height band and latitude band, with its rain, runoff and lightning. home is the ways-of-living
    study's highland home (its highland towns' and terrain-held districts' shares)."""
    runoff = np.maximum(rain - evap, 0.0)

    def summary(mask):
        dry = float((at.dry * mask).sum())
        if dry <= 0.0:
            return dict(dry_land_km2=0.0)
        return dict(dry_land_km2=r(dry, 4), far_side_share=r(float((at.dry * mask * (at.side == 2)).sum()) / dry),
                    mean_height_above_sea_km=r(area_mean(at.hsl_km, at.dry, mask)),
                    rain_mm_day=r(area_mean(rain, at.dry, mask)),
                    runoff_km3_year=r(float((runoff * at.land * mask).sum()) * 1e-6 * JULIAN_YEAR_DAYS),
                    flashes_per_km2_year=r(area_mean(flashes, at.dry, mask)))
    high = at.hsl_km >= HIGHLAND_KM
    by_height = []
    for lo, hi in HEIGHT_BANDS_KM:
        m = band_mask(at, lo, hi)
        row = dict(height_above_sea_km=[lo, hi], **summary(m))
        row['by_latitude'] = [dict(abs_lat_deg=[a, min(b, 90.0)],
                                   dry_land_km2=r(float((at.dry * m * lat_mask(at, a, b)).sum()), 4))
                              for a, b in LAT_BANDS_DEG]
        by_height.append(row)
    by_latitude = [dict(abs_lat_deg=[a, min(b, 90.0)], **summary(high & lat_mask(at, a, b))) for a, b in LAT_BANDS_DEG]
    total = summary(high)
    pole = high & lat_mask(at, COMFORT_POLEWARD_DEG, 90.01)
    lowlands = [dict(height_above_sea_km=[lo, hi], dry_land_km2=r(float((at.dry * band_mask(at, lo, hi)).sum()), 4))
                for lo, hi in LOW_BANDS_KM]
    return dict(all_highland=total, poleward_of_30=summary(pole), equatorward_of_30=summary(high & ~pole),
                by_height=by_height, by_latitude=by_latitude, lowlands=lowlands,
                dry_land_km2=r(float(at.dry.sum()), 4),
                ways_of_living_home=dict(dry_land_km2=r(home['km2'], 4), far_side_share=home['far_side_share'],
                                         note=('the highland towns\' and terrain-held districts\' shared home: the '
                                               'ground 4 km up left once the metropolis circle, the shores of seas '
                                               'and great lakes and the wet margins have taken theirs')))


def summit_rise(at: grid.Atlas, summit, rain, flashes):
    """Dry land round the summit by distance, with its mean height, rain and lightning."""
    d = grid.great_circle_km(at.LA, at.LO, summit['lat_deg'], summit['lon_deg'])
    out = []
    for radius in SUMMIT_RADII_KM:
        m = d <= radius
        dry = float((at.dry * m).sum())
        out.append(dict(radius_km=radius, dry_land_km2=r(dry, 4),
                        lake_share=r(float((at.land * at.f_lake * m).sum()) / float((at.area * m).sum())),
                        mean_height_above_sea_km=r(area_mean(at.hsl_km, at.dry, m)),
                        rain_mm_day=r(area_mean(rain, at.dry, m)),
                        flashes_per_km2_year=r(area_mean(flashes, at.dry, m)),
                        people_at_metropolis_density_billion=r(dry * DENSITIES['metropolis'] / 1e9)))
    return out


def wettest_and_stormiest(at: grid.Atlas, summit, rain, flashes):
    """The dry-land cells with the most rain and the most flashes on the Moon, where they stand and how far from the
    summit."""
    land = at.dry > 0.0
    out = {}
    for label, field in (('wettest', rain), ('stormiest', flashes)):
        i, j = np.unravel_index(int(np.argmax(np.where(land, field, -np.inf))), field.shape)
        out[label] = dict(lat_deg=r(float(at.lat[i]), 4), lon_deg=r(float(at.lon[j]), 4),
                          height_above_sea_km=r(float(at.hsl_km[i, j])), rain_mm_day=r(float(rain[i, j])),
                          flashes_per_km2_year=r(float(flashes[i, j])),
                          km_from_summit=r(float(grid.great_circle_km(at.lat[i], at.lon[j], summit['lat_deg'],
                                                                      summit['lon_deg']))))
    return out


def comfort(flat, box, land_out):
    """Comfortable hours a lunar day on land by latitude and height: CM1's flat rings corrected for height against
    the GCM's 3-day means, and the 3-D box over the highlands at 44.7 S by height; and the GCM's land air carried to
    each highland latitude band's mean height along its own fitted rate."""
    rows = []
    for c in flat['corrected_for_height']:
        rows.append(dict(abs_lat_deg=[c['from_deg'], c['to_deg']], height_m=[c['from_m'], c['to_m']],
                         cm1_mean_height_km=r(c['mean_height_km']),
                         cm1_strict_hours=r(c['comfortable_hours']['strict']),
                         cm1_loose_hours=r(c['comfortable_hours']['loose']),
                         gcm_strict_hours=r(c['gcm']['comfortable_hours']['strict']),
                         gcm_loose_hours=r(c['gcm']['comfortable_hours']['loose']),
                         gcm_cell_height_km=r(c['gcm']['cell_height_km']),
                         cm1_air_c=r(c['air_c']), gcm_air_c=r(c['gcm']['air_c']),
                         cm1_dewpoint_c=r(c['dewpoint_c']), gcm_dewpoint_c=r(c['gcm']['dewpoint_c'])))
    rates = [dict(abs_lat_deg=[g['from_deg'], g['to_deg']], air_at_sea_level_c=r(g['air_at_sea_level_c']),
                  air_per_km_c=r(g['air_per_km_c']), dewpoint_per_km_c=r(g['dewpoint_per_km_c']))
             for g in flat['gcm_height_rates']]
    warmth = []
    for band in land_out['by_latitude']:
        z = band.get('mean_height_above_sea_km')
        if z is None:
            continue
        mid = 0.5 * (band['abs_lat_deg'][0] + band['abs_lat_deg'][1])
        g = next(x for x in flat['gcm_height_rates'] if x['from_deg'] <= mid < x['to_deg'])
        warmth.append(dict(abs_lat_deg=band['abs_lat_deg'], mean_height_km=z,
                           gcm_land_air_c=r(g['air_at_sea_level_c'] + g['air_per_km_c'] * z),
                           register_lapse_air_c=r(g['air_at_sea_level_c'] - REGISTER_LAPSE_K_PER_KM * z)))
    by_height = [dict(group=b['group'], mean_height_km=r(b['mean_height_m'] / 1e3),
                      strict_hours=r(b['comfortable_hours']['strict']), loose_hours=r(b['comfortable_hours']['loose']),
                      day_air_median_c=r(b['air_day_c']['median']), night_air_median_c=r(b['air_night_c']['median']),
                      too_hot_hours=r(b['uncomfortable_hours']['too_hot']),
                      too_cold_hours=r(b['uncomfortable_hours']['too_cold']))
                 for b in box['by_height']]
    return dict(lunar_day_hours=r(LUNAR_DAY_H, 4),
                definition=('comfortable: air at 2 m of 18-26 C with a dewpoint of at most 15 C (strict), or 17-27 C '
                            'and 16 C (loose), as the climate products define it'),
                rings_corrected_for_height=rows, gcm_land_air_by_latitude=rates, highland_warmth=warmth,
                highland_box=dict(site=box['site'], by_height=by_height, gcm_cells=box['gcm']),
                land_shares=flat['land_shares'])


# ---------------------------------------------------------------------------------------------------------------
# How tall buildings get.

def moon_air(winds, ground_km):
    zk = np.asarray(winds['height_km'], float) * 1e3
    rho = np.log(np.asarray(winds['density_mean_kg_m3'], float))
    return lambda z: np.exp(np.interp(ground_km * 1e3 + np.asarray(z, float), zk, rho))


def earth_air(z):
    """The International Standard Atmosphere's density (kg/m3) in the troposphere, as the summit tower study."""
    z = np.asarray(z, dtype=float)
    R, lapse = 287.05287, 0.0065
    t = 288.15 - lapse * np.minimum(z, 11e3)
    p = 101325.0 * (t / 288.15) ** (STANDARD_GRAVITY / (R * lapse))
    return p / (R * t)


def size_tower(storeys, world, design, service):
    """The lightest lattice frame for a tower of the given storeys, every storey a full plate, meeting the summit
    tower study's checks (sway within 1/500 of the height in the service gust, buckling factor at least 3, base no
    wider than half the height)."""
    H = storeys * STOREY_M
    floors = lt.Floors(spacing_m=STOREY_M, storeys=1, storey_height_m=STOREY_M, area_per_storey_m2=PLATE_M ** 2,
                       plan_share=1.0, mass_per_area_kg_m2=FLOOR_MASS_KG_M2)
    frame = lt.Frame(solidity=CLOSED_FACADE_SOLIDITY, top_width_m=PLATE_M)
    shares = BASE_SHARES[BASE_SHARES * H >= PLATE_M]
    if shares.size == 0 or shares[0] * H > PLATE_M:
        shares = np.concatenate([[PLATE_M / H], shares])
    S, P = np.meshgrid(shares, EXPONENTS, indexing='ij')
    res = lt.size(H, H * S.ravel(), P.ravel(), design, world, frame, floors, steps=300)
    ok = ((res['sway_per_wind2_m'] * service ** 2 <= SWAY_LIMIT * H) & (res['buckling_factor'] >= 3.0)
          & (res['base_width_m'] <= 0.5 * H))
    if not ok.any():
        return None
    i = int(np.flatnonzero(ok)[np.argmin(res['frame_mass_kg'][ok])])
    pick = {k: float(v[i]) for k, v in res.items()}
    sway = pick['sway_per_wind2_m'] * service ** 2
    from_wind = pick['base_moment_nm'] / (math.sqrt(2.0) * pick['base_width_m'])
    windward = world.gravity_m_s2 * pick['total_mass_kg'] / 4.0 - from_wind
    return dict(storeys=storeys, height_m=H, base_width_m=r(pick['base_width_m']), flare_exponent=pick['exponent'],
                floor_area_m2=r(pick['floor_area_m2']), frame_steel_t=r(pick['frame_mass_kg'] / 1e3),
                frame_steel_kg_per_m2=r(pick['frame_mass_kg'] / pick['floor_area_m2']),
                period_s=r(1.0 / pick['frequency_hz']), sway_service_m=r(sway),
                top_acceleration_m_s2=r((2.0 * math.pi * pick['frequency_hz']) ** 2 * FLUCTUATING_SWAY * sway),
                buckling_factor=r(pick['buckling_factor']),
                wind_share_of_leeward_leg=r(from_wind / pick['base_leg_force_n']),
                windward_leg_force_mn=r(windward / 1e6),
                windward_leg_in_tension=bool(windward < 0.0))


def canyon(storeys, coverage):
    """Street canyons between towers on a square grid: towers of the plate's width at the pitch that gives the
    coverage. The sky view factor at the middle of an infinitely long canyon's floor is cos(atan(2H/W)), the view
    factor of the open strip above it. The towers stand apart, so the streets' sides are open between them and the
    true view is wider: the canyon's factor is a lower bound."""
    pitch = PLATE_M / math.sqrt(coverage)
    gap = pitch - PLATE_M
    H = storeys * STOREY_M
    return dict(street_m=r(gap), height_to_width=r(H / gap),
                sky_view_factor_at_least=r(math.cos(math.atan(2.0 * H / gap))),
                open_share_of_street_side=r(gap / pitch))


def towers(winds, tower_study, per_billion_heat_tw):
    """Steel per floor for towers of given storeys on the Moon's highlands and on Earth, the wind against the weight,
    the people per km2 they hold and the heat and daylight that cap that."""
    design = tower_study['design_winds']['model_based_design_gust_m_s']
    service = tower_study['design_winds']['service_gust_m_s']
    moon = lt.World('moon', MOON_SURFACE_GRAVITY, moon_air(winds, TOWER_GROUND_KM))
    earth = lt.World('earth', STANDARD_GRAVITY, earth_air)
    earth_service = EARTH_DESIGN_GUST * EARTH_SERVICE_SHARE
    rows = []
    for n in STOREYS:
        m = size_tower(n, moon, design, service)
        e = size_tower(n, earth, EARTH_DESIGN_GUST, earth_service)
        density = {f'{c:g}': r(c * n * 1e6 / FLOOR_M2_PER_PERSON) for c in COVERAGES}
        heat = [r(c * n * 1e6 / FLOOR_M2_PER_PERSON * kw * 1e3 / 1e6) for c in (0.2,) for kw in per_billion_heat_tw]
        rows.append(dict(storeys=n, height_m=n * STOREY_M, moon=m, earth=e,
                         earth_to_moon_steel=(r(e['frame_steel_kg_per_m2'] / m['frame_steel_kg_per_m2'])
                                              if m and e else None),
                         people_per_km2_by_coverage=density,
                         heat_w_m2_at_coverage_0_2=heat,
                         canyon_at_coverage_0_2=canyon(n, 0.2)))
    # The Moon's height whose frame steel per floor equals Earth's at each height, by interpolation in the table.
    ms = np.array([row['moon']['frame_steel_kg_per_m2'] for row in rows if row['moon']])
    mh = np.array([row['height_m'] for row in rows if row['moon']])
    same = []
    for row in rows:
        if row['earth'] and row['earth']['frame_steel_kg_per_m2'] <= ms.max():
            h = float(np.interp(row['earth']['frame_steel_kg_per_m2'], ms, mh))
            same.append(dict(earth_height_m=row['height_m'], moon_height_m=r(h), ratio=r(h / row['height_m'])))
    rho_m = float(moon.air_density(0.0))
    q_m = 0.5 * rho_m * design ** 2
    q_e = 0.5 * float(earth_air(0.0)) * EARTH_DESIGN_GUST ** 2
    designed = []
    for n, (lo, hi) in DESIGNED_LATERAL_KG_M2.items():
        e = size_tower(n, earth, EARTH_DESIGN_GUST, earth_service)
        if e:
            designed.append(dict(storeys=n, model_frame_kg_per_m2=e['frame_steel_kg_per_m2'],
                                 designed_lateral_kg_per_m2=[lo, hi],
                                 model_share=[r(e['frame_steel_kg_per_m2'] / hi), r(e['frame_steel_kg_per_m2'] / lo)],
                                 source=DESIGNED_LATERAL_SOURCE))
    return dict(
        earth_frame_against_designed_towers=designed,
        basis=dict(storey_m=STOREY_M, plate_m=PLATE_M, floor_mass_kg_m2=FLOOR_MASS_KG_M2,
                   moon_ground_above_sea_km=TOWER_GROUND_KM, moon_design_gust_m_s=design, moon_service_gust_m_s=service,
                   earth_design_gust_m_s=EARTH_DESIGN_GUST, earth_service_gust_m_s=r(earth_service),
                   material=lt.STEEL.name, sway_limit=SWAY_LIMIT, floor_m2_per_person=FLOOR_M2_PER_PERSON),
        wind=dict(moon_air_density_kg_m3=r(rho_m), moon_design_pressure_pa=r(q_m),
                  earth_air_density_kg_m3=r(float(earth_air(0.0))), earth_design_pressure_pa=r(q_e),
                  pressure_ratio=r(q_m / q_e), gravity_ratio=r(MOON_SURFACE_GRAVITY / STANDARD_GRAVITY),
                  wind_to_weight_against_earth=r((q_m / q_e) / (MOON_SURFACE_GRAVITY / STANDARD_GRAVITY))),
        self_weight_length_km=dict(moon=tower_study['materials']['steel']['self_weight_length_km_moon'],
                                   earth=tower_study['materials']['steel']['self_weight_length_km_earth']),
        by_storeys=rows, same_steel_per_floor=same)


# ---------------------------------------------------------------------------------------------------------------
# People, water and travel.

def ground_people(land_out, summit_out):
    """Settled land for the ground's people in the highlands at the planning densities, against the high ground."""
    homes = {'all_highland': land_out['all_highland']['dry_land_km2'],
             'ways_of_living_home': land_out['ways_of_living_home']['dry_land_km2'],
             'poleward_of_30': land_out['poleward_of_30']['dry_land_km2']}
    rows = []
    for g in GROUND_BILLIONS:
        for s in HIGHLAND_SHARES:
            for name, dens in DENSITIES.items():
                km2 = g * s * 1e9 / dens
                rows.append(dict(ground_billion=g, highland_share=s, density=name, people_per_km2=dens,
                                 settled_km2=r(km2),
                                 share_of=dict({k: r(km2 / v) for k, v in homes.items()},
                                               **{f'all_highland_outside_kept_{int(k * 100)}':
                                                  r(km2 / (homes['all_highland'] * (1 - k))) for k in KEPT_SHARES})))
    return dict(homes_km2=homes, densities=DENSITIES, rows=rows,
                summit_100km=dict(dry_land_km2=summit_out[2]['dry_land_km2'],
                                  people_at_40000_billion=summit_out[2]['people_at_metropolis_density_billion']))


def lift_energy(km3_year, lift_m, efficiency):
    """Mean power (GW) to lift a yearly volume of water by lift_m at a pump efficiency."""
    return km3_year * 1e9 * WATER_DENSITY * MOON_SURFACE_GRAVITY * lift_m / efficiency / SECONDS_PER_YEAR / 1e9


def weighted_pct(values, weights, qs):
    order = np.argsort(values, kind='stable')
    v, w = values[order], weights[order]
    cw = np.cumsum(w) / w.sum()
    return [float(v[min(np.searchsorted(cw, q), v.size - 1)]) for q in qs]


def water(at: grid.Atlas, land_out, domestic_km3_per_billion, efficiency):
    """Lifting domestic water to the highlands from the nearest standing water (the seas, atlas lakes or rain-fed
    lakes) and from the sea itself, with and without recovery, and how many people each latitude band's own runoff
    would supply."""
    km, nearest = at.km_to_water()
    level = at.water_level_m().ravel()[nearest]
    lift_km = np.maximum(at.height_m - level, 0.0) / 1e3
    high = (at.hsl_km >= HIGHLAND_KM) & (at.dry > 0)
    per_km = lift_energy(1.0, 1000.0, efficiency)          # GW per km3 a year per km of lift
    lo, hi = domestic_km3_per_billion
    nearest_rows = []
    for a, b in LAT_BANDS_DEG:
        m = high & lat_mask(at, a, b)
        if not m.any():
            continue
        lifts = weighted_pct(lift_km[m], at.dry[m], (0.5, 0.9))
        dist = weighted_pct(km[m], at.dry[m], (0.5, 0.9))
        sea_lift = weighted_pct(at.hsl_km[m], at.dry[m], (0.5,))[0]
        nearest_rows.append(dict(abs_lat_deg=[a, min(b, 90.0)], lift_km_p50_p90=[r(x) for x in lifts],
                                 km_to_water_p50_p90=[r(x) for x in dist], lift_from_sea_km_p50=r(sea_lift),
                                 gw_per_billion_from_nearest_p50=[r(per_km * v * lifts[0]) for v in (lo, hi)],
                                 gw_per_billion_from_sea_p50=[r(per_km * v * sea_lift) for v in (lo, hi)]))
    m = high
    all_lift = weighted_pct(lift_km[m], at.dry[m], (0.5, 0.9))
    sea_lift = weighted_pct(at.hsl_km[m], at.dry[m], (0.5,))[0]
    recovery = {f'{k:g}': dict(from_nearest_gw=[r(per_km * v * all_lift[0] * (1 - k)) for v in (lo, hi)],
                               from_sea_gw=[r(per_km * v * sea_lift * (1 - k)) for v in (lo, hi)])
                for k in RECOVERY}
    local = []
    for row in land_out['by_latitude']:
        runoff = row.get('runoff_km3_year') or 0.0
        local.append(dict(abs_lat_deg=row['abs_lat_deg'], runoff_km3_year=runoff,
                          billion_supplied_no_recovery=[r(runoff / v) for v in (hi, lo)],
                          billion_supplied_90_recovered=[r(runoff / (v * 0.1)) for v in (hi, lo)]))
    return dict(domestic_km3_per_billion=[r(lo), r(hi)], pump_efficiency=efficiency,
                gw_per_km3_year_per_km_of_lift=r(per_km), kwh_per_m3_per_km_of_lift=r(
                    WATER_DENSITY * MOON_SURFACE_GRAVITY * 1e3 / efficiency / 3.6e6),
                highland_lift_km_p50_p90=[r(x) for x in all_lift], highland_lift_from_sea_km_p50=r(sea_lift),
                by_latitude=nearest_rows, per_billion_by_recovery=recovery, local_runoff=local)


def access(fleet):
    """Trips between highland cities and the coasts by sky boat, sky ferry and regional ship. The sky boat's battery
    holds its usable energy over the reserve; the fleet sizes it for the boat's range in cruise with its hover, so a
    climb comes out of that range: the range with a climb is the range less the cruise the climb's energy buys."""
    boat = next(x for x in fleet['electric']['moon'] if x['name'].startswith('sky boat'))
    tech = fleet['electric']['technology']
    eff = tech['efficiency']
    ferry = next(x for x in fleet['sky_ferries'] if x['length_m'] == 250.0)
    regional = next(x for x in fleet['regional']['ships'] if x['length_m'] == 350.0)
    climb_kwh_per_km = boat['mass_kg'] * MOON_SURFACE_GRAVITY * 1e3 / eff / 3.6e6
    cruise_kwh_per_km = boat['cruise_kw'] / (boat['cruise_m_s'] * 3.6)
    mission_kwh = boat['battery_kg'] * tech['battery_wh_kg'] * tech['usable'] / tech['reserve'] / 1e3
    allowance_kwh = mission_kwh - boat['range_km'] * cruise_kwh_per_km
    climb_kwh = CLIMB_M / 1e3 * climb_kwh_per_km
    range_with_climb = boat['range_km'] - climb_kwh / cruise_kwh_per_km
    trips = []
    for km in TRIP_KM:
        trips.append(dict(km=km,
                          sky_boat=dict(kwh_with_climb=r(km * cruise_kwh_per_km + climb_kwh + allowance_kwh),
                                        hours=r(km / (boat['cruise_m_s'] * 3.6)),
                                        within_range_with_climb=km <= range_with_climb),
                          sky_ferry=dict(wh_per_passenger=r(km * ferry['wh_per_passenger_km']),
                                         hours=r(km / (ferry['cruise_m_s'] * 3.6)),
                                         within_route=km <= ferry['route_km']),
                          regional_ship=dict(wh_per_passenger=r(km * regional['wh_per_passenger_km']),
                                             hours=r(km / (regional['cruise_m_s'] * 3.6)),
                                             within_route=km <= regional['route_km'])))
    return dict(sky_boat=dict(mass_kg=boat['mass_kg'], seats=boat['seats'], cruise_m_s=boat['cruise_m_s'],
                              range_km=boat['range_km'], cruise_kwh_per_km=r(cruise_kwh_per_km),
                              climb_kwh_per_km_of_height=r(climb_kwh_per_km), mission_kwh=r(mission_kwh),
                              hover_allowance_kwh=r(allowance_kwh), climb_kwh=r(climb_kwh),
                              range_with_climb_km=r(range_with_climb)),
                sky_ferry=dict(length_m=ferry['length_m'], passengers=ferry['passengers'],
                               line_per_hour_per_direction=ferry['line_per_hour_per_direction'],
                               wh_per_passenger_km=ferry['wh_per_passenger_km'], cruise_m_s=ferry['cruise_m_s'],
                               route_km=ferry['route_km']),
                regional_ship=dict(length_m=regional['length_m'], passengers=regional['passengers'],
                                   wh_per_passenger_km=regional['wh_per_passenger_km'],
                                   cruise_m_s=regional['cruise_m_s'], route_km=regional['route_km']),
                climb_m_for_trips=CLIMB_M, trips=trips)


def scattered(settings, ground_billion=GROUND_BILLIONS, highland_shares=(0.5, 0.75)):
    """The communities beyond the highlands: each other surface home's land, and how the ground's people beyond
    the highlands spread over them in proportion to land (the ways-of-living study's spread rule)."""
    names = ('coastal_and_lake_towns', 'wetland_edge_towns', 'canopy_towns', 'farming_country', 'nearside_towns',
             'limb_towns', 'far_side_places', 'high_latitude_towns')
    homes = {k: float(settings[k]['home_dry_land_km2']) for k in names}
    total = sum(homes.values())
    rows = []
    for g in ground_billion:
        for s in highland_shares:
            rest = g * (1.0 - s)
            per_km2 = rest * 1e9 / total
            rows.append(dict(ground_billion=g, highland_share=s, scattered_billion=r(rest),
                             people_per_km2_of_home=r(per_km2),
                             town_share_of_home_at_10000_per_km2=r(per_km2 / DENSITIES['towns']),
                             by_home_billion={k: r(rest * v / total) for k, v in homes.items()}))
    return dict(homes_km2={k: r(v, 4) for k, v in homes.items()}, total_km2=r(total, 4), rows=rows,
                undersea_people=settings['undersea_communities']['design']['residents_all_candidates'])
