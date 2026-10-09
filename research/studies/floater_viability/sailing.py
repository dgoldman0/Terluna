"""Sailing and shape: tethered wings in shear, sea drogues, streamlining and the Sun's pace.

Pure functions plus evaluate(). Winds are eastward-positive (toward later local
time). The global-mean profile is the design GCM's layer-mean eastward wind;
the CM1 rings give the wind by local hour at 0, 45 and 80 degrees.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from shared.constants import JULIAN_DAY, MOON_GM, MOON_RADIUS, SYNODIC_MONTH_DAYS
from research.studies.floater_viability.photoperiod import solar_height, sun_follow_speed

ROOT = Path(__file__).resolve().parents[3]
RINGS = ('ring_equator', 'ring_45n_lsw', 'ring_45s_lsw', 'ring_80n_lsw', 'ring_80s_lsw')
INPUT_FILES = (
    'climate/gcm/products/climatology_A28_dim5_moon.npz',
    'climate/results/gcm/global_winds_A28_dim5_moon.json',
    'research/studies/sky_ships/results/sky_ships.json',
    *(f'climate/results/crm/ring_{name}.json' for name in RINGS),
    'research/studies/floater_viability/sailing_sources.json',
)
WATER_DENSITY_KG_M3 = 1025.  # seawater, shared/scenarios/waves.json convention
METRES_PER_DEGREE = math.pi*MOON_RADIUS/180


def _positive(**values):
    if any(not math.isfinite(v) or v <= 0 for v in values.values()):
        raise ValueError(f'finite positive inputs required: {tuple(values)}')


def _nonnegative(**values):
    if any(not math.isfinite(v) or v < 0 for v in values.values()):
        raise ValueError(f'finite nonnegative inputs required: {tuple(values)}')


# --- the wind profile and the length of a floater's day -------------------------

def sigma_heights(sigma, surface_pressure_pa, heights_km, pressures_pa, extrapolate_below=False):
    """Height of each sigma level, by log-pressure interpolation in a mean profile.

    With extrapolate_below, a level below the profile's lowest point takes the
    log-pressure slope of the two lowest points.
    """
    lp = [math.log(p) for p in pressures_pa]
    out = []
    for s in sigma:
        target = math.log(s*surface_pressure_pa)
        if target > lp[0]:
            if extrapolate_below:
                slope = (heights_km[1]-heights_km[0])/(lp[1]-lp[0])
                out.append(heights_km[0]+(target-lp[0])*slope)
            else:
                out.append(heights_km[0])
            continue
        if target < lp[-1]:
            out.append(None)
            continue
        for i in range(1, len(lp)):
            if target >= lp[i]:
                f = (lp[i-1]-target)/(lp[i-1]-lp[i])
                out.append(heights_km[i-1]+f*(heights_km[i]-heights_km[i-1]))
                break
    return out


def relative_solar_period(latitude_deg, height_m, eastward_wind_m_s):
    """Days for a body riding a steady eastward wind to see the Sun return.

    A body at fixed ground position sees one synodic month; moving east at u
    against the Sun's westward pace v makes the local time advance (v+u)/v
    times faster. Westward wind faster than the Sun runs local time backward.
    """
    v = sun_follow_speed(latitude_deg, height_m)
    if v == 0:
        return None
    rate = v + eastward_wind_m_s
    if rate == 0:
        return None
    return SYNODIC_MONTH_DAYS*v/abs(rate)


def hour_attractors(hour_angles_deg, eastward_m_s, latitude_deg, height_m):
    """Local hours where the wind holds a drifting body at the Sun's pace.

    dH/dt is proportional to v_sun + u(H). A crossing of zero with u falling
    through it is stable: bodies converge there. Returns each crossing with its
    solar elevation, stability and the e-folding time of approach in days.
    """
    v = sun_follow_speed(latitude_deg, height_m)
    radius = (MOON_RADIUS+height_m)*math.cos(math.radians(latitude_deg))
    hs = list(hour_angles_deg)
    us = list(eastward_m_s)
    n = len(hs)
    out = []
    for i in range(n):
        j = (i+1) % n
        h0, h1 = hs[i], hs[j] + (360 if j == 0 else 0)
        a, b = v+us[i], v+us[j]
        if a == 0 or a*b < 0:
            f = 0. if a == 0 else a/(a-b)
            h = h0 + f*(h1-h0)
            h = (h+180) % 360 - 180
            slope = (us[j]-us[i])/math.radians(h1-h0)  # m/s per radian of hour angle
            stable = slope < 0
            efold = abs(radius/slope)/JULIAN_DAY if slope else None
            out.append(dict(hour_angle_deg=h, solar_elevation_deg=solar_height(latitude_deg, h),
                            stable=stable, efold_days=efold, sun_pace_m_s=v))
    return out


def circuit_period(hour_angles_deg, eastward_m_s, latitude_deg, height_m, steps=7200):
    """Days for a body riding the ring's hour-varying wind to come round to the same hour.

    The circuit integral of (R cos lat) dH/(v_sun + u(H)); None where the wind
    holds the body (v_sun + u reaches zero somewhere on the circuit).
    """
    v = sun_follow_speed(latitude_deg, height_m)
    radius = (MOON_RADIUS+height_m)*math.cos(math.radians(latitude_deg))
    hs = list(hour_angles_deg)
    ext_h = [h-360 for h in hs] + hs + [h+360 for h in hs]
    ext_u = list(eastward_m_s)*3
    total = 0.
    for i in range(steps):
        h = -180 + (i+.5)*360/steps
        rate = v + float(np.interp(h, ext_h, ext_u))
        if rate <= 0:
            return None
        total += radius*math.radians(360/steps)/rate
    return total/JULIAN_DAY


def hour_dwell(hour_angles_deg, eastward_m_s, latitude_deg, height_m, half_width_deg=10., westward_scale=1.):
    """Days a drifting body spends within half_width of the hour where it drifts slowest.

    Integrates dt = (R cos(lat)) dH/(v_sun + u(H)) across the window with u
    interpolated in hour angle; None where the body is held (a stable crossing).
    westward_scale multiplies the westward (negative) winds, a test of how far the
    result rests on the strength of the westward branch.
    """
    v = sun_follow_speed(latitude_deg, height_m)
    radius = (MOON_RADIUS+height_m)*math.cos(math.radians(latitude_deg))
    hs = list(hour_angles_deg)
    us = [u*westward_scale if u < 0 else u for u in eastward_m_s]
    rates = [v+u for u in us]
    k = min(range(len(hs)), key=lambda i: abs(rates[i]))
    centre = hs[k]
    ext_h = [h-360 for h in hs] + hs + [h+360 for h in hs]
    ext_u = us*3
    steps = 400
    total = 0.
    for i in range(steps):
        h = centre - half_width_deg + (i+.5)*2*half_width_deg/steps
        rate = v + float(np.interp(h, ext_h, ext_u))
        if rate <= 0:
            return dict(centre_hour_deg=centre, slowest_drift_m_s=rates[k], dwell_days=None,
                        solar_elevation_deg=solar_height(latitude_deg, centre))
        total += radius*math.radians(2*half_width_deg/steps)/rate
    return dict(centre_hour_deg=centre, slowest_drift_m_s=rates[k], half_width_deg=half_width_deg,
                dwell_days=total/JULIAN_DAY, solar_elevation_deg=solar_height(latitude_deg, centre))


# --- sailing between two moving media --------------------------------------------

def _drag(rho, cda, r):
    m = math.hypot(*r)
    return (.5*rho*cda*m*r[0], .5*rho*cda*m*r[1])


def tethered_sail(body_drag_area_m2, wing_area_m2, delta_wind_m_s, body_density_kg_m3,
                  wing_density_kg_m3, lift_coefficient=1., wing_drag_coefficient=.1,
                  tether_length_m=0., tether_diameter_m=0., tether_drag_coefficient=1.2,
                  tack=1):
    """Steady velocity of a body and a wing or drogue in a second medium.

    Frame: the body's own air. The other medium moves at delta_wind along +x.
    Body drag on -u, wing drag and lift on (delta - u), and a vertical tether's
    drag through the linear shear between them balance. Returns u, the
    system's velocity relative to the body's air: along the shear and across it.
    """
    _positive(body=body_drag_area_m2, rho_b=body_density_kg_m3, rho_w=wing_density_kg_m3)
    _nonnegative(wing=wing_area_m2, dw=delta_wind_m_s, cl=lift_coefficient, cd=wing_drag_coefficient,
                 tether=tether_length_m, diameter=tether_diameter_m)
    if tack not in (1, -1):
        raise ValueError('tack is +1 or -1')
    dw = (delta_wind_m_s, 0.)
    tether_cda = tether_drag_coefficient*tether_diameter_m*tether_length_m
    nodes = (0., .25, .5, .75, 1.)
    weights = (1/12, 1/3, 1/6, 1/3, 1/12)  # composite Simpson on four panels

    def force(u):
        fb = _drag(body_density_kg_m3, body_drag_area_m2, (-u[0], -u[1]))
        rw = (dw[0]-u[0], dw[1]-u[1])
        mw = math.hypot(*rw)
        q = .5*wing_density_kg_m3*wing_area_m2*mw
        fw = (q*(wing_drag_coefficient*rw[0] - tack*lift_coefficient*rw[1]),
              q*(wing_drag_coefficient*rw[1] + tack*lift_coefficient*rw[0]))
        ft = [0., 0.]
        if tether_cda:
            rho_t = .5*(body_density_kg_m3+wing_density_kg_m3)
            for z, w in zip(nodes, weights):
                r = (z*dw[0]-u[0], z*dw[1]-u[1])
                d = _drag(rho_t, tether_cda, r)
                ft[0] += w*d[0]
                ft[1] += w*d[1]
        return (fb[0]+fw[0]+ft[0], fb[1]+fw[1]+ft[1]), fw

    if delta_wind_m_s == 0 or wing_area_m2 == 0:
        return dict(along_m_s=0., across_m_s=0., speed_m_s=0., wing_force_n=0.)
    ratio = math.sqrt(wing_density_kg_m3*wing_area_m2*max(lift_coefficient, wing_drag_coefficient)
                      / (body_density_kg_m3*body_drag_area_m2))
    u = [delta_wind_m_s*ratio/(1+ratio)*.5, tack*delta_wind_m_s*ratio/(1+ratio)*.5]
    for _ in range(200):
        f, _ = force(u)
        if math.hypot(*f) < 1e-9*max(1., body_drag_area_m2):
            break
        h = 1e-6*max(1., delta_wind_m_s)
        f1, _ = force((u[0]+h, u[1]))
        f2, _ = force((u[0], u[1]+h))
        j = ((f1[0]-f[0])/h, (f2[0]-f[0])/h, (f1[1]-f[1])/h, (f2[1]-f[1])/h)
        det = j[0]*j[3]-j[1]*j[2]
        if det == 0:
            raise ValueError('singular sailing balance')
        du = ((j[3]*f[0]-j[1]*f[1])/det, (-j[2]*f[0]+j[0]*f[1])/det)
        step = 1.
        while step > 1e-4:
            trial = (u[0]-step*du[0], u[1]-step*du[1])
            if math.hypot(*force(trial)[0]) < math.hypot(*f):
                u = list(trial)
                break
            step /= 2
        else:
            break
    f, fw = force(u)
    return dict(along_m_s=u[0], across_m_s=u[1], speed_m_s=math.hypot(*u),
                wing_force_n=math.hypot(*fw), residual_n=math.hypot(*f))


def tether_diameter(force_n, allowable_pa):
    """Round tether that carries `force_n` at its allowable stress."""
    _nonnegative(force=force_n)
    _positive(allowable=allowable_pa)
    return math.sqrt(4*force_n/(math.pi*allowable_pa))


def hanging_depth(horizontal_force_n, wing_weight_n, line_weight_n_m, length_m):
    """Vertical depth of a tether's lower end below its top.

    Constant horizontal tension H; vertical tension W + w s at arc length s above
    the wing. D = (sqrt(H² + (W + wL)²) - sqrt(H² + W²))/w.
    """
    h, wb, w, length = horizontal_force_n, wing_weight_n, line_weight_n_m, length_m
    if w == 0:
        return length*wb/math.hypot(h, wb) if wb else 0.
    return (math.hypot(h, wb+w*length)-math.hypot(h, wb))/w


def hanging_tether(horizontal_force_n, length_m, allowable_pa, density_kg_m3, gravity_m_s2, depth_share=.9):
    """Lightest tether and wing ballast that hang the wing depth_share x L below the body.

    The tether carries its own weight, the wing's weighted ballast and the horizontal
    pull; its top tension sqrt(H² + (W + wL)²) may not exceed the allowable stress
    times its section. Scans the section in steps of 0.01 decade, then refines around
    the lightest design in steps of 0.0002 decade; None if no section serves.
    """
    _positive(force=horizontal_force_n, length=length_m, allowable=allowable_pa, density=density_kg_m3,
              gravity=gravity_m_s2)
    if not 0 < depth_share < 1:
        raise ValueError('depth share in (0, 1)')
    h, length, sig = horizontal_force_n, length_m, allowable_pa
    target = depth_share*length

    def design(a):
        w = density_kg_m3*gravity_m_s2*a
        room = (sig*a)**2 - h*h
        if room <= 0:
            return None
        cap = math.sqrt(room) - w*length
        if cap < 0 or hanging_depth(h, cap, w, length) < target:
            return None
        if hanging_depth(h, 0., w, length) >= target:
            wb = 0.
        else:
            lo, hi = 0., cap
            for _ in range(100):
                mid = (lo+hi)/2
                lo, hi = (lo, mid) if hanging_depth(h, mid, w, length) >= target else (mid, hi)
            wb = hi
        mass = (wb + w*length)/gravity_m_s2
        return dict(area_m2=a, diameter_m=math.sqrt(4*a/math.pi), tether_kg=w*length/gravity_m_s2,
                    ballast_kg=wb/gravity_m_s2, hanging_kg=mass, top_tension_n=math.hypot(h, wb+w*length),
                    depth_m=hanging_depth(h, wb, w, length), ballast_over_pull=wb/h,
                    free_hanging_limit_m=sig/(density_kg_m3*gravity_m_s2))

    best, best_k = None, None
    for k in range(1, 401):
        d = design(h/sig*10**(k/100))
        if d is not None and (best is None or d['hanging_kg'] < best['hanging_kg']):
            best, best_k = d, k
    if best is None:
        return None
    for j in range(-50, 51):
        d = design(h/sig*10**(best_k/100 + j/5000))
        if d is not None and d['hanging_kg'] < best['hanging_kg']:
            best = d
    return best


def meridional_reach(across_m_s, days=SYNODIC_MONTH_DAYS):
    """Degrees of latitude a steady cross-wind speed carries a body in `days`."""
    return abs(across_m_s)*days*JULIAN_DAY/METRES_PER_DEGREE


def sail_case(body_radius_m, body_cd, wing_share, delta_wind_m_s, tether_length_m,
              body_density, wing_density, allowable_pa=28e6, tether_density=1500., cl=1., wing_cd=.1,
              gravity_m_s2=1.62, depth_share=.9):
    """A sphere-equivalent body, its wing sized as a share of its projected area, and a hanging tether.

    The tether is sized for its own weight, the wing's ballast and the body's pull
    so that the wing hangs depth_share of the tether's length below the body; the
    balance is solved again with the sized tether's drag until the section settles.
    """
    area = math.pi*body_radius_m**2
    d = 0.
    hang = None
    for _ in range(8):
        result = tethered_sail(body_cd*area, wing_share*area, delta_wind_m_s, body_density, wing_density,
                               cl, wing_cd, tether_length_m, d)
        pull = .5*body_density*body_cd*area*result['speed_m_s']**2
        if not tether_length_m or pull == 0:
            break
        hang = hanging_tether(pull, tether_length_m, allowable_pa, tether_density, gravity_m_s2, depth_share)
        if hang is None:
            break
        if abs(hang['diameter_m']-d) <= 1e-6*hang['diameter_m']:
            break
        d = hang['diameter_m']
    per_area = (lambda key: hang[key]/area if hang else None)
    return dict(body_radius_m=body_radius_m, body_drag_coefficient=body_cd, wing_share=wing_share,
                delta_wind_m_s=delta_wind_m_s, tether_length_m=tether_length_m, allowable_mpa=allowable_pa/1e6,
                tether_diameter_m=hang['diameter_m'] if hang else None,
                tether_kg_per_projected_m2=per_area('tether_kg'), ballast_kg_per_projected_m2=per_area('ballast_kg'),
                hanging_kg_per_projected_m2=per_area('hanging_kg'),
                wing_depth_m=hang['depth_m'] if hang else None,
                ballast_over_pull=hang['ballast_over_pull'] if hang else None,
                free_hanging_limit_m=allowable_pa/(tether_density*gravity_m_s2),
                degrees_latitude_per_lunar_day=meridional_reach(result['across_m_s']), **result)


def evaluate(root=ROOT):
    """Global-mean shear, day lengths by height, hour attractors and sailing reach."""
    root = Path(root)
    clim = np.load(root/INPUT_FILES[0], allow_pickle=False)
    winds = json.loads((root/INPUT_FILES[1]).read_text())
    ships = json.loads((root/INPUT_FILES[2]).read_text())
    meta = json.loads(str(clim['metadata']))
    if meta.get('schema') != 'terluna.climate.gcm-climatology/1' or winds.get('schema') != 'terluna.climate.gcm-global-winds/1':
        raise ValueError('unexpected wind inputs')
    sigma = [float(s) for s in clim['sigma']]
    ua = [float(u) for u in clim['ua_layer_mean']]
    surface = meta['configuration']['pressure_pa']
    # The wind product's 0 km pressure is not hydrostatic with the levels above it (they agree
    # with rho g dz to 1% from 2.5 km up), so heights come from the 2.5 km level upward and the
    # lowest sigma level is placed by extrapolating that profile's log pressure.
    keep = [i for i, z in enumerate(winds['height_km']) if z >= 2.5]
    heights = sigma_heights(sigma, surface, [winds['height_km'][i] for i in keep],
                            [winds['pressure_pa'][i] for i in keep], extrapolate_below=True)
    profile = [dict(sigma=s, height_km=h, eastward_m_s=u) for s, h, u in zip(sigma, heights, ua)]
    resolved = sorted((p for p in profile if p['height_km'] is not None), key=lambda p: p['height_km'])
    zs = [p['height_km'] for p in resolved]
    us = [p['eastward_m_s'] for p in resolved]

    def wind_at(z):
        return float(np.interp(z, zs, us))

    shear = [dict(lower_km=a, upper_km=b, delta_m_s=wind_at(b)-wind_at(a),
                  per_km=(wind_at(b)-wind_at(a))/(b-a)) for a, b in ((zs[0], 10.), (5., 10.), (10., 20.), (10., 30.), (20., 40.))]
    westernmost = min(us)
    periods = [dict(latitude_deg=lat, height_km=z, eastward_m_s=wind_at(z),
                    sun_pace_m_s=sun_follow_speed(lat, z*1000),
                    relative_period_days=relative_solar_period(lat, z*1000, wind_at(z)),
                    longest_geometric_night_days=relative_solar_period(lat, z*1000, wind_at(z))/2)
               for lat in (0., 30., 45., 60.) for z in (5., 10., 20., 30., 40.)]
    rings = {}
    for name, path in zip(RINGS, INPUT_FILES[3:3+len(RINGS)]):
        ring = json.loads((root/path).read_text())
        if ring.get('schema') != 'terluna.climate.crm-ring/1':
            raise ValueError('unexpected ring')
        c = ring['circulation']
        lat = ring['latitude_deg']
        rows = []
        for z, row in zip(c['heights_km'], c['eastward_m_s']):
            if z > 40:
                continue
            mean = sum(row)/len(row)
            rows.append(dict(height_km=z, hour_mean_eastward_m_s=mean, westmost_m_s=min(row),
                             westmost_hour_deg=c['hour_angle_deg'][row.index(min(row))],
                             sun_pace_m_s=sun_follow_speed(lat, z*1000),
                             relative_period_days=circuit_period(c['hour_angle_deg'], row, lat, z*1000),
                             relative_period_days_at_hour_mean_wind=relative_solar_period(lat, z*1000, mean),
                             attractors=hour_attractors(c['hour_angle_deg'], row, lat, z*1000),
                             dwell=hour_dwell(c['hour_angle_deg'], row, lat, z*1000),
                             dwell_if_westward_weaker={f'{x:g}': hour_dwell(c['hour_angle_deg'], row, lat, z*1000,
                                                                          westward_scale=x)['dwell_days']
                                                       for x in (.95, .9, .8)}))
        rings[name] = dict(latitude_deg=lat, heights=rows)
    rho = {a['height_km']: a['density_kg_m3'] for a in ships['air']}
    g10 = MOON_GM/(MOON_RADIUS+10000.)**2
    sails = []
    for allowable in (28e6, 92e6):
        for radius in (20., 100., 1000.):
            for cd in (.47, .05):
                for share in (.01, .1, .3):
                    for low, high in ((zs[0], 10.), (5., 10.)):
                        dw = wind_at(high)-wind_at(low)
                        sails.append(dict(lower_km=low, body_km=high, **sail_case(
                            radius, cd, share, abs(dw), (high-low)*1000/.9, rho[10],
                            float(np.interp(low, [0, 10], [rho[0], rho[10]])), allowable_pa=allowable,
                            gravity_m_s2=g10)))
    sea = []
    for wind in (2., 4., 6.):
        for radius in (5., 20.):
            for share in (.001, .01, .05):
                r = tethered_sail(.47*math.pi*radius**2, share*math.pi*radius**2, wind, rho[0],
                                  WATER_DENSITY_KG_M3, 1., .1)
                sea.append(dict(wind_m_s=wind, body_radius_m=radius, drogue_share=share,
                                drift_along_wind_m_s=wind-r['along_m_s'], across_m_s=r['across_m_s'],
                                course_off_downwind_deg=math.degrees(math.atan2(abs(r['across_m_s']), wind-r['along_m_s']))))
    return dict(
        evidence='Design GCM layer-mean eastward wind (10 years of 3-day means), the CM1 rings composited by local hour at 0, +-45 and +-80 degrees, and steady force balances for a body with a tethered wing or a sea drogue. The rings are two-dimensional: air cannot converge on their storms from the sides, which likely strengthens the low westward branch, so the sunset dwell is given with weaker branches beside it. Wing and drogue sizes, lift and drag coefficients and tether strength are stated design choices; tethers carry their own weight and the wing\'s ballast.',
        reading_rule='Eastward positive, toward later local time. A body moving east sees a shorter solar cycle; one moving west at the Sun\'s pace holds its hour. Sail speeds are relative to the body\'s own air; across is perpendicular to the shear. A latitude and solar-hour wind product from another branch will replace the single global-mean profile and the five rings.',
        global_mean_profile=profile, shear=shear, westernmost_mean_wind_m_s=westernmost,
        day_length=periods, rings=rings, tethered_sails=sails, sea_drogue=sea,
        streamlining=dict(sphere_cd=.47, streamlined_cd=.05, drag_ratio=.47/.05,
                          sail_speed_gain_small_wing=math.sqrt(.47/.05)),
        unresolved=['Winds by latitude and local hour at every height (another branch builds the product).',
                    'North-south winds: the rings carry none, so the mean meridional drift is open.',
                    'Wing and tether dynamics, tether fluttering and the drag of a real fibre bundle.',
                    'Biological control of a kilometres-long tether and its wing.'])
