"""The wind as a vector over the whole Moon, in the ground's frame and the Sun's, and the solar days it gives a body
that drifts with it, written as a data product for other lanes.

    climate/gcm/.venv/bin/python -m climate.gcm.zonal_winds A28_dim5_moon:20-29
    # -> results/gcm/zonal_winds_A28_dim5_moon.json, with routes and finer Sun-frame fields in the matching .npz

The subsolar point moves west over the Moon at v_sun = 2 pi (R + z) cos(lat) / T_synodic, 4.30 m/s at the equator
10 km up. A body drifting east at u passes through the local solar hours at u + v_sun and lives a solar day
T' = T_synodic v_sun / |u + v_sun|: u = 0 gives the 29.53-day lunar day, an eastward wind shortens it, a westward
wind slower than v_sun lengthens it, u = -v_sun holds the hour, and a faster westward wind turns the Sun back so that
it rises in the west.

For the chosen model years it places every layer of every 3-day output mean at its height above the model's sea
level, as global_winds.py does, and interpolates the eastward and northward wind linearly to fixed heights, 0-90 km
in 2.5 km steps. Heights below a cell's ground are left out of the statistics; between the ground and the lowest
layer's midpoint the lowest layer's wind stands. It gathers, weighted by cell area where cells are pooled:
the wind's mean and percentiles by latitude and height; the same in the Sun's frame, by local solar hour angle
(degrees east of the subsolar point, from climatology.subsolar_longitudes, in 36 bins), latitude and height; how
often the wind holds the hour (u within 0.5 or 1 m/s of -v_sun), blows westward, or turns the Sun back; the solar
day of a body drifting with the local wind; the vertical shear over 2.5, 5, 10 and 15 km; and how much of the wind's
variation about its zonal mean the Sun's frame, the geography and the hour at each place explain.

It then carries passive parcels horizontally at fixed heights in the Sun's frame, through the steady mean flow of
that frame and through the 3-day means in sequence, and records the time each takes to pass once through all solar
hours, its routes, how long it stays within one 30-degree span of hours, its time in sunlight, and where parcels
gather; and it finds the points where the steady flow stops and what kind of point each is.
"""
from __future__ import annotations
import hashlib
import json
import re
import sys
from pathlib import Path
import numpy as np

from climate.gcm.climatology import HOUR_BINS, subsolar_longitudes
from climate.gcm.global_winds import SEA_TOLERANCE_M, to_heights, weighted_percentiles
from climate.gcm.site_winds import RESULTS, layer_heights
from shared.constants import JULIAN_DAY, MOON_GM, MOON_RADIUS, SYNODIC_MONTH_DAYS

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = 'terluna.climate.gcm-zonal-winds/1'
HEIGHTS_KM = np.arange(0.0, 90.01, 2.5)
SUN_RATE = 2.0 * np.pi / (SYNODIC_MONTH_DAYS * JULIAN_DAY)     # rad/s: the subsolar point's westward turn
HOLD_M_S = (0.5, 1.0)
SHEAR_KM = (2.5, 5.0, 10.0, 15.0)
PERCENTILES = (5, 10, 25, 50, 75, 90, 95)
SUN_PERCENTILES = (10, 50, 90)
DAY_PERCENTILES = (10, 25, 50, 75, 90)
DAY_CLASS_EDGES_D = (0.0, 7.0, 15.0, SYNODIC_MONTH_DAYS, 60.0, 180.0, np.inf)
U_BINS = np.arange(-40.0, 60.0001, 0.05)
V_BINS = np.arange(-30.0, 30.0001, 0.05)
SUN_U_BINS = np.arange(-40.0, 60.0001, 0.25)
SUN_V_BINS = np.arange(-30.0, 30.0001, 0.25)
SHEAR_BINS = np.arange(0.0, 40.0001, 0.05)
DAY_BINS_D = np.geomspace(0.1, 1.0e5, 601)
LATITUDE_BANDS_DEG = ((0.0, 20.0), (20.0, 40.0), (40.0, 60.0), (60.0, 75.0), (75.0, 90.0))
HEIGHT_BANDS_KM = ((0.0, 5.0), (5.0, 15.0), (15.0, 30.0), (30.0, 45.0), (45.0, 60.0), (60.0, 90.0))
TRAJECTORY_KM = (2.5, 5.0, 7.5, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0)
RELEASE_LAT_DEG = np.arange(-80.0, 80.1, 10.0)
RELEASE_HOUR_DEG = np.arange(-165.0, 166.0, 30.0)
STEP_S = 3600.0
RECORD_S = 12 * 3600.0
STEADY_YEARS = 5
ROUTE_DAYS = 60
DWELL_SPAN_DEG = 30.0
DWELL_DAYS = (1, 2, 3, 5, 7, 10, 15, 20, 30, 45, 60, 90, 120, 180, 270, 360, 540, 720, 1080, 1440, 1800, 2520, 3240)
OCCUPANCY_HOURS, OCCUPANCY_LATS = 12, 18
GATHER_DAYS = (0, 30, 60, 90, 180, 357, 730, 1095, 1787, 2502, 3575)
POSITION_EVERY_DAYS = 90
EVIDENCE = ('ExoPlaSim 3-day output means at T21 (about 170 km cells) on 10 sigma layers over the whole Moon, '
            'interpolated in height to fixed heights above the model\'s sea level. Averaging over 3 days removes '
            'storms, gusts and turbulence: the winds, days and shears are those of the large-scale flow, and each '
            'mean spans 36 degrees of the Sun\'s path. The model\'s ground is its cells\' mean elevation. The layer '
            'midpoints stand near 0.9, 4, 10, 18, 29, 42, 56, 74 and 98 km above the ground, so a shear over less than '
            'their spacing is the linear share of the shear between two layers, and above 74 km the wind is '
            'interpolated toward the layer near 98 km, the last below the model\'s top layer. Parcels are carried '
            'horizontally at fixed heights by the 3-day mean wind alone, as a body that holds its height and does not '
            'steer would be.')
READING_RULE = (
    "Heights are above the model's sea level (height_km), latitudes the model's rows (lat_deg, north first) and hour "
    'angles degrees east of the subsolar point (hour_angle_deg, bin centres, 0 noon, -90 sunrise, +90 sunset). Arrays'
    ' run [height][lat], [height][lat][hour] or, for shear, [dz][height][lat] with dz from shear_dz_km and the shear '
    'from height_km to height_km + dz. Statistics at a height use only the samples whose ground lies below it '
    '(coverage); null marks no samples. u is eastward, v northward, in m/s, 3-day means. sun_speed_m_s is v_sun, the '
    "subsolar point's westward speed at that latitude and height. The Sun-relative eastward speed of a drifting body "
    'is u + v_sun; hold_0_5 and hold_1 are the shares of samples with |u + v_sun| under 0.5 and 1 m/s, westward the '
    "share with u < 0, reversed the share with u + v_sun < 0 (the Sun rises in the west). solar_day gives the day T' "
    "= T_synodic v_sun / |u + v_sun| (Earth days) of a body drifting with each sample's wind: its percentiles over "
    "both directions, and the shares of samples in each class of day_class_edges_days for the Sun's ordinary motion "
    '(normal) and turned back (reversed), [class][height][lat]. Percentiles come from histograms (0.05 m/s bins; 0.25'
    " m/s in the Sun's frame). variance gives, by height, the shares of the variance of u and v about their mean by "
    'latitude and height that the mean by hour (Sun frame), by place (geography) and by hour at each place explain, '
    'adjusted for the number of cells, area-weighted. bands pool the arrays over bands of height and absolute '
    'latitude, weighted by area. trajectories: parcels released at release_lat_deg x release_hour_deg at each of '
    "heights_km are carried horizontally in the Sun's frame by the Sun-relative wind; where the ground stands above a"
    ' parcel it rides the lowest layer, about 0.9 km above the ground. steady uses the mean Sun-frame flow for '
    "steady_years; sequence uses the 3-day means in order, each placed under the Sun's fitted track and interpolated "
    "linearly in time in the Sun's frame, released at the start of each model year and followed one year (first_year)"
    ' and, for the first release, through all years (long_run). Per-parcel arrays run [height][release lat][release '
    'hour]; first_year summarises its releases [height][release lat]. first_pass_days is the time to move a full 360 '
    "degrees of hour from the start (null if not within the run), first_pass_sign +1 for the Sun's ordinary way and "
    "-1 turned back; mean_day_days the run's length over its net circuits, negative with the Sun turned back; "
    'dwell_days the longest of dwell_lengths_days over which the parcel stays within a 30-degree span of hours; '
    "sun_up_share its share of time with the Sun above the horizon; light its mean cosine of the Sun's zenith angle "
    '(0 at night; 1/pi for a fixed place on the equator) and light_ratio that over the mean fixed places at the '
    'latitudes it passed would get; sunsets and western_sunrises count its sunsets and the times the Sun rose again '
    'in the west, steps passing close to a pole left out (near_pole_steps). occupancy is the share of parcel time, '
    'weighted by the area each release stands for, in 30-degree sectors of hour and 10-degree bands of latitude, '
    '[height][lat band][hour sector], over the second half of the steady and long runs and over the whole first year '
    'of the releases, beside area_share, the share a uniform spread would give. gathering gives, [day][height], the '
    'median distance from each parcel to its nearest neighbour released at the same height (and, in first_year, the '
    'same time) and the share with a neighbour within one model cell (cell_km). fixed_points lists where the steady '
    'flow stops at each trajectory height: a sink gathers parcels, a source spreads them, a saddle does both; '
    'e_folding_days is the slower rate at which parcels near it close in or spread. The .npz beside this file holds '
    '[height][lat][hour] arrays: the Sun-frame percentiles and spread (sun_u_p10_m_s ... sun_v_sd_m_s), shares '
    '(sun_hold_0_5, sun_westward), sample counts (sun_samples) and, [dz][height][lat][hour], the mean shear '
    '(shear_sun_*); steady_u_m_s and steady_v_m_s [trajectory height][lat][hour], the flow the steady routes follow; '
    "the daily routes of the first route_days (steady_routes, sequence_routes) and the long run's positions every "
    'position_every_days (long_run_positions), each [record][parcel][hour, lat] in hundredths of a degree, parcels in'
    ' the order height, release lat, release hour.')


def sun_speed(lat_deg, height_m):
    """Westward speed (m/s) of the subsolar point under a body at that latitude and height."""
    return SUN_RATE * (MOON_RADIUS + np.asarray(height_m, dtype=float)) * np.cos(np.radians(lat_deg))


def day_rate(u, v_sun):
    """Solar days per lunar day of a body drifting east at u: 1 at rest, 0 holding the hour, negative with the Sun
    turned back."""
    return (np.asarray(u, dtype=float) + v_sun) / v_sun


def solar_day(u, v_sun):
    """Length (Earth days) of the solar day of a body drifting east at u, under a subsolar point moving west at
    v_sun; infinite where the drift holds the hour."""
    with np.errstate(divide='ignore'):
        return SYNODIC_MONTH_DAYS * np.asarray(v_sun, dtype=float) / np.abs(np.asarray(u, dtype=float) + v_sun)


def hour_angles(sub, lon):
    """Hour angle (degrees east of the subsolar point, -180..180) of every column at every output, and its bin, as
    climatology.sun_composite counts them."""
    hour = (np.asarray(lon)[None, :] - np.asarray(sub)[:, None] + 180.0) % 360.0 - 180.0
    return hour, np.clip(((hour + 180.0) / 360.0 * HOUR_BINS).astype(int), 0, HOUR_BINS - 1)


def hour_centres():
    return (np.arange(HOUR_BINS) + 0.5) / HOUR_BINS * 360.0 - 180.0


def sun_track(sub, step_deg):
    """The Sun's longitude at each output on a line of the given westward step, fitted to subsolar longitudes that
    the grid rounds to its columns."""
    n = np.arange(np.size(sub))
    start = np.degrees(np.angle(np.exp(1j * np.radians(np.asarray(sub) + step_deg * n)).mean()))
    return start - step_deg * n


def ground_heights(ps_pa, psl_pa, t0_k, q0, gas_constant):
    """Each sample's ground height (m) above the model's sea level, from the sea-level and surface pressures at the
    lowest layer's virtual temperature, as global_winds.py finds it."""
    from atmosphere.radiative_convective.thermodynamics import R_VAPOUR
    tv0 = t0_k * (1.0 + q0 * (R_VAPOUR / gas_constant - 1.0))
    phi = gas_constant * tv0 * np.log(np.maximum(psl_pa, ps_pa) / ps_pa)
    return MOON_GM / (MOON_GM / MOON_RADIUS - phi) - MOON_RADIUS


def read_year(path, gas_constant, keep_km=TRAJECTORY_KM):
    """One model year: u and v at the fixed heights (height, time, lat, lon), NaN below the ground; the same at the
    heights keep_km with the lowest layer's wind standing below the ground; the ground's height; the subsolar
    longitudes."""
    import netCDF4
    with netCDF4.Dataset(path) as d:
        lat = np.asarray(d['lat'][:], dtype=float)
        lon = np.asarray(d['lon'][:], dtype=float)
        sigma = np.asarray(d['lev'][:], dtype=float)
        up = np.argsort(sigma)[::-1]                         # lowest layer first
        field = lambda name: np.moveaxis(np.asarray(d[name][:], dtype=float)[:, up], 1, -1).reshape(-1, sigma.size)
        u, v, t, q = (field(n) for n in ('ua', 'va', 'ta', 'hus'))
        ps = np.asarray(d['ps'][:], dtype=float).reshape(-1) * 100.0      # hPa, as labelled
        psl = np.asarray(d['psl'][:], dtype=float).reshape(-1)            # Pa, although labelled hPa
        czen = np.asarray(d['czen'][:], dtype=float)
    shape = (HEIGHTS_KM.size, czen.shape[0], lat.size, lon.size)
    above, _ = layer_heights(sigma[up], ps, t, q, gas_constant, MOON_GM, MOON_RADIUS)
    ground = ground_heights(ps, psl, t[:, 0], q[:, 0], gas_constant)
    z = above + ground[:, None]
    hm = HEIGHTS_KM * 1e3
    uz, vz = to_heights(z, u, hm).reshape(shape), to_heights(z, v, hm).reshape(shape)
    keep = [int(np.argmin(np.abs(HEIGHTS_KM - k))) for k in keep_km]
    fill_u, fill_v = uz[keep].astype(np.float32), vz[keep].astype(np.float32)
    under = (hm[:, None] < ground[None, :] - SEA_TOLERANCE_M).reshape(shape)
    uz[under] = np.nan
    vz[under] = np.nan
    return dict(lat=lat, lon=lon, u=uz, v=vz, fill_u=fill_u, fill_v=fill_v,
                ground=ground.reshape(shape[1:]).astype(np.float32), sub=subsolar_longitudes(czen, lat, lon))


def bin_index(x, edges):
    return np.clip(((x - edges[0]) / (edges[1] - edges[0])).astype(np.int64), 0, edges.size - 2)


class Totals:
    """Running sums over the samples of the chosen years."""

    def __init__(self, nlat, nlon, ntraj):
        nh, nb, nd, nk = HEIGHTS_KM.size, HOUR_BINS, len(SHEAR_KM), len(HOLD_M_S)
        z = lambda *s, dtype=float: np.zeros(s, dtype=dtype)
        self.nlat, self.nlon, self.outputs = nlat, nlon, 0
        self.n, self.su, self.sv, self.suu, self.svv = (z(nh, nlat) for _ in range(5))
        self.hist_u, self.hist_v = z(nh, nlat, U_BINS.size - 1), z(nh, nlat, V_BINS.size - 1)
        self.hold, self.west, self.back = z(nk, nh, nlat), z(nh, nlat), z(nh, nlat)
        self.day = z(2, nh, nlat, DAY_BINS_D.size - 1)
        self.sn, self.ssu, self.ssv, self.ssuu, self.ssvv = (z(nh, nlat, nb) for _ in range(5))
        self.shold, self.swest, self.sback = z(nk, nh, nlat, nb), z(nh, nlat, nb), z(nh, nlat, nb)
        self.shist_u = z(nh, nlat, nb, SUN_U_BINS.size - 1, dtype=np.int32)
        self.shist_v = z(nh, nlat, nb, SUN_V_BINS.size - 1, dtype=np.int32)
        self.gn, self.gsu, self.gsv = (z(nh, nlat, nlon) for _ in range(3))
        self.pn, self.psu, self.psv = (z(nh, nlat, nlon, nb) for _ in range(3))
        self.dn, self.dsu, self.dsv, self.dss = (z(nd, nh, nlat) for _ in range(4))
        self.dhist = z(nd, nh, nlat, SHEAR_BINS.size - 1)
        self.sdn, self.sdsu, self.sdsv, self.sdss = (z(nd, nh, nlat, nb) for _ in range(4))
        self.tn, self.tsu, self.tsv = (z(ntraj, nlat, nb) for _ in range(3))

    def add(self, year):
        u, v, lat = year['u'], year['v'], year['lat']
        nh, times, nlat, nlon = u.shape
        nb = HOUR_BINS
        self.outputs += times
        _, k = hour_angles(year['sub'], year['lon'])                         # (time, lon)
        shape = u.shape
        cell = np.broadcast_to((np.arange(nh)[:, None] * nlat + np.arange(nlat)[None, :])[:, None, :, None], shape)
        hour = np.broadcast_to(k[None, :, None, :], shape)
        col = np.broadcast_to(np.arange(nlon)[None, None, None, :], shape)
        vsun = np.broadcast_to(sun_speed(lat[None, :], HEIGHTS_KM[:, None] * 1e3)[:, None, :, None], shape)
        ok = np.isfinite(u)
        c, kk, j, uu, vv, vs = cell[ok], hour[ok], col[ok], u[ok], v[ok], vsun[ok]
        us = uu + vs
        nc = nh * nlat
        add = lambda total, index, weights=None: total.__iadd__(
            np.bincount(index, weights, minlength=total.size).reshape(total.shape))
        for total, w in ((self.n, None), (self.su, uu), (self.sv, vv), (self.suu, uu * uu), (self.svv, vv * vv),
                         (self.west, (uu < 0.0).astype(float)), (self.back, (us < 0.0).astype(float))):
            add(total, c, w)
        for m, limit in enumerate(HOLD_M_S):
            add(self.hold[m], c, (np.abs(us) < limit).astype(float))
        add(self.hist_u, c * (U_BINS.size - 1) + bin_index(uu, U_BINS))
        add(self.hist_v, c * (V_BINS.size - 1) + bin_index(vv, V_BINS))
        with np.errstate(divide='ignore'):
            day = SYNODIC_MONTH_DAYS * vs / np.abs(us)
        b = np.clip(np.searchsorted(DAY_BINS_D, day) - 1, 0, DAY_BINS_D.size - 2)
        add(self.day, ((us < 0.0) * nc + c) * (DAY_BINS_D.size - 1) + b)
        s = c * nb + kk
        for total, w in ((self.sn, None), (self.ssu, uu), (self.ssv, vv), (self.ssuu, uu * uu),
                         (self.ssvv, vv * vv), (self.swest, (uu < 0.0).astype(float)),
                         (self.sback, (us < 0.0).astype(float))):
            add(total, s, w)
        for m, limit in enumerate(HOLD_M_S):
            add(self.shold[m], s, (np.abs(us) < limit).astype(float))
        self.shist_u += np.bincount(s * (SUN_U_BINS.size - 1) + bin_index(uu, SUN_U_BINS),
                                    minlength=self.shist_u.size).reshape(self.shist_u.shape).astype(np.int32)
        self.shist_v += np.bincount(s * (SUN_V_BINS.size - 1) + bin_index(vv, SUN_V_BINS),
                                    minlength=self.shist_v.size).reshape(self.shist_v.shape).astype(np.int32)
        g = c * nlon + j
        for total, w in ((self.gn, None), (self.gsu, uu), (self.gsv, vv)):
            add(total, g, w)
        p = g * nb + kk
        for total, w in ((self.pn, None), (self.psu, uu), (self.psv, vv)):
            add(total, p, w)
        del c, kk, j, uu, vv, vs, us, day, b, s, g, p
        for d, dz in enumerate(SHEAR_KM):
            m = int(round(dz / (HEIGHTS_KM[1] - HEIGHTS_KM[0])))
            du, dv = u[m:] - u[:-m], v[m:] - v[:-m]
            okd = np.isfinite(du) & np.isfinite(dv)
            cd = (cell[:nh - m] + 0)[okd]
            kd = hour[:nh - m][okd]
            du, dv = du[okd], dv[okd]
            ds = np.hypot(du, dv)
            for total, w in ((self.dn[d], None), (self.dsu[d], du), (self.dsv[d], dv), (self.dss[d], ds)):
                total[:nh - m] += np.bincount(cd, w, minlength=(nh - m) * nlat).reshape(nh - m, nlat)
            self.dhist[d, :nh - m] += np.bincount(cd * (SHEAR_BINS.size - 1) + bin_index(ds, SHEAR_BINS),
                                                  minlength=(nh - m) * nlat * (SHEAR_BINS.size - 1)
                                                  ).reshape(nh - m, nlat, -1)
            sd = cd * nb + kd
            for total, w in ((self.sdn[d], None), (self.sdsu[d], du), (self.sdsv[d], dv), (self.sdss[d], ds)):
                total[:nh - m] += np.bincount(sd, w, minlength=(nh - m) * nlat * nb).reshape(nh - m, nlat, nb)
        nt = year['fill_u'].shape[0]
        tc = np.broadcast_to((np.arange(nt)[:, None] * nlat + np.arange(nlat)[None, :])[:, None, :, None] * nb
                             + k[None, :, None, :], (nt, times, nlat, nlon)).ravel()
        for total, w in ((self.tn, None), (self.tsu, year['fill_u'].ravel().astype(float)),
                         (self.tsv, year['fill_v'].ravel().astype(float))):
            add(total, tc, w)


def share_explained(n, s, total_n, total_s, total_q):
    """Share of the sum of squares about each row's mean (rows: the leading axes of total_*) that the means of the
    cells along the remaining axes of n and s explain, adjusted for the number of cells."""
    axes = tuple(range(total_n.ndim, n.ndim))
    with np.errstate(invalid='ignore', divide='ignore'):
        sst = total_q - total_s ** 2 / total_n
        ssb = np.where(n > 0, s ** 2 / np.maximum(n, 1), 0.0).sum(axis=axes) - total_s ** 2 / total_n
        cells = (n > 0).sum(axis=axes)
        ssb_adj = ssb - (cells - 1) * (sst - ssb) / np.maximum(total_n - cells, 1)
    return sst, ssb_adj, ssb


def r(a, n=2):
    """Rounded nested lists with null for non-finite values and 0 for zeros."""
    a = np.asarray(a, dtype=float)
    out = np.round(a, n).astype(object)
    out[np.round(a, n) == 0] = 0
    out[~np.isfinite(a)] = None
    return out.tolist()


def summarise(T: Totals, lat) -> tuple[dict, dict]:
    """The product's statistics from the running sums, and the finer arrays for the .npz."""
    nh, nb = HEIGHTS_KM.size, HOUR_BINS
    vsun = sun_speed(lat[None, :], HEIGHTS_KM[:, None] * 1e3)
    with np.errstate(invalid='ignore', divide='ignore'):
        n = np.where(T.n > 0, T.n, np.nan)
        sn = np.where(T.sn > 0, T.sn, np.nan)
        up = weighted_percentiles(T.hist_u, U_BINS, PERCENTILES)
        vp = weighted_percentiles(T.hist_v, V_BINS, PERCENTILES)
        ground = dict(coverage=r(T.n / (T.outputs * T.nlon), 3), u_mean_m_s=r(T.su / n), v_mean_m_s=r(T.sv / n),
                      **{f'u_p{p}_m_s': r(up[i]) for i, p in enumerate(PERCENTILES)},
                      **{f'v_p{p}_m_s': r(vp[i]) for i, p in enumerate(PERCENTILES)})
        sun_u, sun_v = T.ssu / sn, T.ssv / sn
        sun = dict(u_mean_m_s=r(sun_u), v_mean_m_s=r(sun_v))
        holding = dict(hold_0_5=r(T.hold[0] / n, 3), hold_1=r(T.hold[1] / n, 3), westward=r(T.west / n, 3),
                       reversed=r(T.back / n, 3),
                       sun_frame=dict(hold_1=r(T.shold[1] / sn, 3), reversed=r(T.sback / sn, 3)))
        both = T.day.sum(axis=0)
        dp = day_percentiles(both, DAY_PERCENTILES)
        classes = day_classes(T.day)
        solar = dict(day_class_edges_days=[round(float(e), 2) if np.isfinite(e) else None for e in DAY_CLASS_EDGES_D],
                     **{f'p{p}_days': r(dp[i], 1) for i, p in enumerate(DAY_PERCENTILES)},
                     normal=r(classes[0] / n[None], 3), reversed=r(classes[1] / n[None], 3))
        nd = len(SHEAR_KM)
        dn = np.where(T.dn > 0, T.dn, np.nan)
        spct = weighted_percentiles(T.dhist, SHEAR_BINS, (10, 50, 90))
        shear = dict(du_mean_m_s=r(T.dsu / dn), dv_mean_m_s=r(T.dsv / dn), speed_mean_m_s=r(T.dss / dn),
                     speed_p10_m_s=r(spct[0]), speed_p50_m_s=r(spct[1]), speed_p90_m_s=r(spct[2]))
        variance = variance_shares(T, lat)
        sdn = np.where(T.sdn > 0, T.sdn, np.nan)
        su_p = weighted_percentiles(T.shist_u, SUN_U_BINS, SUN_PERCENTILES)
        sv_p = weighted_percentiles(T.shist_v, SUN_V_BINS, SUN_PERCENTILES)
        arrays = dict(
            sun_u_mean_m_s=sun_u, sun_v_mean_m_s=sun_v,
            sun_u_sd_m_s=np.sqrt(np.maximum(T.ssuu / sn - sun_u ** 2, 0.0)),
            sun_v_sd_m_s=np.sqrt(np.maximum(T.ssvv / sn - sun_v ** 2, 0.0)),
            **{f'sun_u_p{p}_m_s': su_p[i] for i, p in enumerate(SUN_PERCENTILES)},
            **{f'sun_v_p{p}_m_s': sv_p[i] for i, p in enumerate(SUN_PERCENTILES)},
            sun_hold_0_5=T.shold[0] / sn, sun_westward=T.swest / sn,
            sun_samples=T.sn, shear_sun_du_mean_m_s=T.sdsu / sdn, shear_sun_dv_mean_m_s=T.sdsv / sdn,
            shear_sun_speed_mean_m_s=T.sdss / sdn)
    bands = band_summary(T, lat, vsun)
    stats = dict(sun_speed_m_s=r(vsun, 3), ground_frame=ground, sun_frame=sun, holding=holding, solar_day=solar,
                 shear_dz_km=list(SHEAR_KM), shear=shear, variance=variance, bands=bands)
    return stats, arrays


def day_percentiles(hist, percentiles):
    """Percentiles (days) of the solar day from histograms on DAY_BINS_D (last axis), log-interpolated."""
    log_edges = np.log(DAY_BINS_D)
    p = weighted_percentiles(hist, log_edges, percentiles)
    return np.exp(p)


def day_classes(day_hist):
    """Sample counts (direction, class, ...) in the classes of DAY_CLASS_EDGES_D from the fine day histograms; a
    class edge falls between fine bins at its nearest bin edge."""
    edges = np.searchsorted(DAY_BINS_D, np.minimum(DAY_CLASS_EDGES_D, DAY_BINS_D[-1]))
    edges[0], edges[-1] = 0, DAY_BINS_D.size - 1
    return np.stack([day_hist[..., a:b].sum(axis=-1) for a, b in zip(edges[:-1], edges[1:])], axis=1)


def variance_shares(T: Totals, lat) -> dict:
    """By height: shares of the variance of u and v about their mean by latitude and height that the Sun-frame
    mean, the geographic mean and the mean by hour at each place explain, area-weighted, with the variation's rms
    and the rms left about the Sun-frame mean."""
    w = np.cos(np.radians(lat))[None, :]
    out = {}
    for name, ts, tq, sun_s, geo_n, geo_s, place_s in (('u', T.su, T.suu, T.ssu, T.gn, T.gsu, T.psu),
                                                       ('v', T.sv, T.svv, T.ssv, T.gn, T.gsv, T.psv)):
        sst, sun_b, sun_raw = share_explained(T.sn, sun_s, T.n, ts, tq)
        _, geo_b, _ = share_explained(geo_n, geo_s, T.n, ts, tq)
        _, place_b, _ = share_explained(T.pn, place_s, T.n, ts, tq)
        with np.errstate(invalid='ignore', divide='ignore'):
            per = lambda a: np.where(T.n > 0, a / np.maximum(T.n, 1), 0.0)
            total = (w * per(sst)).sum(axis=1)
            share = lambda b: (w * per(b)).sum(axis=1) / total
            out[name] = dict(share_sun_frame=r(share(sun_b), 3), share_geography=r(share(geo_b), 3),
                             share_hour_at_place=r(share(place_b), 3),
                             rms_m_s=r(np.sqrt(total / (w * (T.n > 0)).sum(axis=1)), 2),
                             rms_after_sun_frame_m_s=r(np.sqrt((w * per(sst - sun_raw)).sum(axis=1)
                                                               / (w * (T.n > 0)).sum(axis=1)), 2))
    return out


def band_summary(T: Totals, lat, vsun) -> list:
    """The arrays pooled over bands of height and absolute latitude, each sample weighted by its cell's area."""
    out = []
    centres = hour_centres()
    for h0, h1 in HEIGHT_BANDS_KM:
        hs = (HEIGHTS_KM >= h0) & (HEIGHTS_KM <= h1)
        for a0, a1 in LATITUDE_BANDS_DEG:
            rows = (np.abs(lat) >= a0) & (np.abs(lat) < a1)
            w = np.cos(np.radians(lat))[None, :] * hs[:, None] * rows[None, :]       # (height, lat)
            total = (w * T.n).sum()
            if total <= 0:
                continue
            mean = lambda a: float((w * a).sum() / total)
            pool = lambda hist: (w[..., None] * hist).sum(axis=(0, 1))
            upct = weighted_percentiles(pool(T.hist_u)[None], U_BINS, (10, 50, 90))[:, 0]
            dpct = day_percentiles(pool(T.day.sum(axis=0))[None], DAY_PERCENTILES)[:, 0]
            classes = day_classes(np.stack([pool(T.day[0]), pool(T.day[1])]))
            sw = w[..., None]                                                        # (height, lat, hour)
            hold_by_hour = (sw * T.shold[1]).sum(axis=(0, 1)) / np.maximum((sw * T.sn).sum(axis=(0, 1)), 1e-300)
            back_by_hour = (sw * T.sback).sum(axis=(0, 1)) / np.maximum((sw * T.sn).sum(axis=(0, 1)), 1e-300)
            best = int(np.argmax(hold_by_hour))
            band = dict(height_km=[h0, h1], abs_lat_deg=[a0, a1],
                        coverage=round(float((w * T.n).sum() / (w.sum() * T.outputs * T.nlon)), 3),
                        sun_speed_m_s=[round(float(vsun[hs][:, rows].min()), 2), round(float(vsun[hs][:, rows].max()), 2)],
                        u_mean_m_s=round(mean(T.su), 2), v_mean_m_s=round(mean(T.sv), 2),
                        u_p10_m_s=round(float(upct[0]), 2), u_p50_m_s=round(float(upct[1]), 2),
                        u_p90_m_s=round(float(upct[2]), 2),
                        hold_0_5=round(mean(T.hold[0]), 3), hold_1=round(mean(T.hold[1]), 3),
                        westward=round(mean(T.west), 3), reversed=round(mean(T.back), 3),
                        normal_classes=[round(float(x / total), 3) for x in classes[0]],
                        reversed_classes=[round(float(x / total), 3) for x in classes[1]],
                        **{f'day_p{p}_days': round(float(dpct[i]), 1) for i, p in enumerate(DAY_PERCENTILES)},
                        hold_1_peak_hour_deg=round(float(centres[best]), 1),
                        hold_1_peak_share=round(float(hold_by_hour[best]), 3),
                        hold_1_by_hour=r(hold_by_hour, 3), reversed_by_hour=r(back_by_hour, 3))
            for d, dz in enumerate(SHEAR_KM):
                ok = hs & (HEIGHTS_KM + dz <= HEIGHTS_KM[-1] + 1e-9)
                wd = np.cos(np.radians(lat))[None, :] * ok[:, None] * rows[None, :]
                if (wd * T.dn[d]).sum() <= 0:
                    continue
                sp = weighted_percentiles((wd[..., None] * T.dhist[d]).sum(axis=(0, 1))[None], SHEAR_BINS,
                                          (50, 90))[:, 0]
                key = str(dz).replace('.', '_')
                band[f'shear_{key}km_p50_m_s'] = round(float(sp[0]), 2)
                band[f'shear_{key}km_p90_m_s'] = round(float(sp[1]), 2)
                band[f'shear_{key}km_du_mean_m_s'] = round(float((wd * T.dsu[d]).sum() / (wd * T.dn[d]).sum()), 2)
            out.append(band)
    return out


def east_north(lat_deg, lon_deg):
    """Unit vectors east and north, each (3, ...), at points on the sphere (degrees)."""
    la, lo = np.broadcast_arrays(np.radians(lat_deg), np.radians(lon_deg))
    east = np.stack([-np.sin(lo), np.cos(lo), np.zeros_like(lo)])
    north = np.stack([-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)])
    return east, north


def cartesian(u, v, lat_deg, lon_deg):
    """Horizontal vectors (u east, v north, arrays (..., lat, lon)) as Cartesian components (..., 3, rows, lon) on
    rows of ascending latitude, with a row at each pole holding the mean vector of the nearest row; returned with
    the rows' latitudes."""
    order = np.argsort(lat_deg)
    lat = np.asarray(lat_deg, dtype=float)[order]
    u, v = np.asarray(u, dtype=float)[..., order, :], np.asarray(v, dtype=float)[..., order, :]
    east, north = east_north(lat[:, None], np.asarray(lon_deg, dtype=float)[None, :])
    w = u[..., None, :, :] * east + v[..., None, :, :] * north
    cap = lambda row: np.repeat(row.mean(axis=-1, keepdims=True), w.shape[-1], axis=-1)
    return np.concatenate([cap(w[..., :1, :]), w, cap(w[..., -1:, :])], axis=-2), np.concatenate([[-90.0], lat, [90.0]])


def interpolate(xyz, rows_deg, lon0_deg, step_deg, lat_deg, lon_deg, lead=()):
    """East and north components at points (lat_deg, lon_deg) of a Cartesian vector field xyz (..., 3, rows,
    columns) on rows_deg (ascending) and periodic columns lon0 + k step: bilinear in the components, projected on
    the local east and north. lead holds an index array for each leading axis, picking each point's field."""
    ncol = xyz.shape[-1]
    x = ((np.asarray(lon_deg, dtype=float) - lon0_deg) / step_deg) % ncol
    j0 = np.floor(x).astype(int)
    fx = x - j0
    j0 %= ncol
    j1 = (j0 + 1) % ncol
    i0 = np.clip(np.searchsorted(rows_deg, lat_deg) - 1, 0, rows_deg.size - 2)
    fy = np.clip((lat_deg - rows_deg[i0]) / (rows_deg[i0 + 1] - rows_deg[i0]), 0.0, 1.0)
    comp = np.arange(3)[:, None]
    pick = lambda i, j: xyz[tuple(a[None, :] for a in lead) + (comp, i[None, :], j[None, :])]
    w = ((1 - fy) * (1 - fx) * pick(i0, j0) + (1 - fy) * fx * pick(i0, j1)
         + fy * (1 - fx) * pick(i0 + 1, j0) + fy * fx * pick(i0 + 1, j1))
    east, north = east_north(lat_deg, lon_deg)
    return (w * east).sum(axis=0), (w * north).sum(axis=0)


def position(lat_deg, hour_deg):
    la, ho = np.radians(lat_deg), np.radians(hour_deg)
    return np.stack([np.cos(la) * np.cos(ho), np.cos(la) * np.sin(ho), np.sin(la)], axis=-1)


def coordinates(p):
    return np.degrees(np.arcsin(np.clip(p[:, 2], -1.0, 1.0))), np.degrees(np.arctan2(p[:, 1], p[:, 0]))


def integrate(velocity, lat0, hour0, radius_m, seconds, t0=0.0, step_s=STEP_S, record_s=RECORD_S, over=None):
    """Carry parcels in the Sun's frame for the given seconds with the classical fourth-order Runge-Kutta scheme on
    the sphere. velocity(t, lat, hour) gives each parcel's Sun-relative east and north speed (m/s) and radius_m the
    radius of its sphere; t0 is each parcel's start time (s). Returns the unwrapped hour angle and the latitude at
    every record, the time of the first full pass through 360 degrees of hour (days, NaN if none) and its sign
    (+1 the Sun's ordinary way), the shares of time with the Sun up and, with over(t, lat, hour), over higher ground,
    the mean cosine of the Sun's zenith angle (0 at night) and that of fixed places at the latitudes passed
    (cos(lat) / pi), the sunsets and the times the Sun rose again in the west (the hour falling back through +90
    degrees), and which parcels moved more than 90 degrees of hour in one step, passing close to a pole; such steps
    count no sunset or western sunrise."""
    lat0, hour0 = np.asarray(lat0, dtype=float), np.asarray(hour0, dtype=float)
    p = position(lat0, hour0)
    t = np.zeros_like(lat0) + t0
    n = lat0.size

    def rate(time, q):
        la, ho = coordinates(q)
        us, vs = velocity(time, la, ho)
        east, north = east_north(la, ho)
        return ((us * east + vs * north) / radius_m).T

    steps, every, dt = int(round(seconds / step_s)), int(round(record_s / step_s)), step_s
    lat, hour = coordinates(p)
    unwrapped, start = hour0.copy(), hour0.copy()
    first, sign = np.full(n, np.nan), np.zeros(n)
    sun_up, light, fixed_light, above = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
    sunsets, western = np.zeros(n), np.zeros(n)
    jumped = np.zeros(n, dtype=bool)
    hours, lats = [unwrapped.astype(np.float32)], [lat.astype(np.float32)]
    for s in range(steps):
        k1 = rate(t, p)
        k2 = rate(t + dt / 2, p + dt / 2 * k1)
        k3 = rate(t + dt / 2, p + dt / 2 * k2)
        k4 = rate(t + dt, p + dt * k3)
        p = p + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        p /= np.linalg.norm(p, axis=1, keepdims=True)
        t = t + dt
        lat, new = coordinates(p)
        dh = (new - hour + 180.0) % 360.0 - 180.0
        jump = np.abs(dh) > 90.0
        jumped |= jump
        before = np.abs(unwrapped - start)
        turns = np.floor((unwrapped + dh - 90.0) / 360.0) - np.floor((unwrapped - 90.0) / 360.0)
        sunsets += np.where(jump, 0.0, np.maximum(turns, 0.0))
        western += np.where(jump, 0.0, np.maximum(-turns, 0.0))
        unwrapped, hour = unwrapped + dh, new
        moved = np.abs(unwrapped - start)
        crossed = np.isnan(first) & (moved >= 360.0)
        if crossed.any():
            frac = (360.0 - before[crossed]) / np.maximum(moved[crossed] - before[crossed], 1e-12)
            first[crossed] = (s + frac) * dt / JULIAN_DAY
            sign[crossed] = np.sign(unwrapped[crossed] - start[crossed])
        cz = np.cos(np.radians(lat)) * np.cos(np.radians(hour))
        sun_up += cz > 0.0
        light += np.maximum(cz, 0.0)
        fixed_light += np.cos(np.radians(lat)) / np.pi
        if over is not None:
            above += over(t, lat, hour)
        if (s + 1) % every == 0:
            hours.append(unwrapped.astype(np.float32))
            lats.append(lat.astype(np.float32))
    return dict(hours=np.array(hours), lats=np.array(lats), first_pass_days=first, first_pass_sign=sign,
                sun_up_share=sun_up / steps, light=light / steps, fixed_light=fixed_light / steps,
                sunsets=sunsets, western_sunrises=western,
                over_ground_share=None if over is None else above / steps, jumped=jumped,
                end_hour=unwrapped, end_lat=lat)


def longest_dwell(hours, record_days, span_deg=DWELL_SPAN_DEG, lengths=DWELL_DAYS):
    """For each parcel (columns of hours: unwrapped hour angles at records record_days apart), the longest of the
    lengths (days) over which it stays within a span of span_deg of hours; 0 if none."""
    from scipy.ndimage import maximum_filter1d, minimum_filter1d
    out = np.zeros(hours.shape[1])
    for length in lengths:
        size = int(round(length / record_days)) + 1
        if size > hours.shape[0]:
            break
        half = size // 2
        spread = (maximum_filter1d(hours, size, axis=0) - minimum_filter1d(hours, size, axis=0))
        held = (spread[half:hours.shape[0] - (size - 1 - half)] <= span_deg).any(axis=0)
        out[held] = length
    return out


def fixed_points(us, vs, lat_deg, hour_deg, radius_m):
    """Where a Sun-relative flow (us east, vs north, arrays (lat, hour) on evenly spaced, periodic hours) stops,
    found cell by cell on its bilinear interpolant, with the kind of each point from the flow's gradient there: a
    sink gathers parcels, a source spreads them, a saddle gathers along one direction and spreads along the other.
    e_folding_days is 1 / the smaller magnitude of the real parts of the gradient's eigenvalues."""
    order = np.argsort(lat_deg)
    lat = np.asarray(lat_deg, dtype=float)[order]
    us, vs = np.asarray(us, dtype=float)[order], np.asarray(vs, dtype=float)[order]
    nh = us.shape[1]
    step = 360.0 / nh
    bil = lambda c, x, y: c[0] * (1 - x) * (1 - y) + c[1] * x * (1 - y) + c[2] * (1 - x) * y + c[3] * x * y
    dbx = lambda c, x, y: (c[1] - c[0]) * (1 - y) + (c[3] - c[2]) * y
    dby = lambda c, x, y: (c[2] - c[0]) * (1 - x) + (c[3] - c[1]) * x
    found, seen = [], []
    for i in range(lat.size - 1):
        for j in range(nh):
            jj = (j + 1) % nh
            a = np.array([us[i, j], us[i, jj], us[i + 1, j], us[i + 1, jj]])
            b = np.array([vs[i, j], vs[i, jj], vs[i + 1, j], vs[i + 1, jj]])
            if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
                continue
            if a.min() > 0 or a.max() < 0 or b.min() > 0 or b.max() < 0:
                continue
            for x, y in ((0.5, 0.5), (0.1, 0.1), (0.9, 0.1), (0.1, 0.9), (0.9, 0.9)):
                for _ in range(50):
                    jac = np.array([[dbx(a, x, y), dby(a, x, y)], [dbx(b, x, y), dby(b, x, y)]])
                    if abs(np.linalg.det(jac)) < 1e-14:
                        break
                    dx, dy = np.linalg.solve(jac, [bil(a, x, y), bil(b, x, y)])
                    x, y = x - dx, y - dy
                    if abs(dx) + abs(dy) < 1e-12:
                        break
                if not (-1e-9 <= x <= 1 + 1e-9 and -1e-9 <= y <= 1 + 1e-9):
                    continue
                if abs(bil(a, x, y)) + abs(bil(b, x, y)) > 1e-8:
                    continue
                la = lat[i] + y * (lat[i + 1] - lat[i])
                ho = (hour_deg[0] + (j + x) * step + 180.0) % 360.0 - 180.0
                if any(abs(a0 - la) < 1e-3 and abs((h0 - ho + 180.0) % 360.0 - 180.0) < 1e-3 for a0, h0 in seen):
                    continue
                seen.append((la, ho))
                mx = radius_m * np.cos(np.radians(la)) * np.radians(step)
                my = radius_m * np.radians(lat[i + 1] - lat[i])
                grad = np.array([[dbx(a, x, y) / mx, dby(a, x, y) / my], [dbx(b, x, y) / mx, dby(b, x, y) / my]])
                trace, det = np.trace(grad), np.linalg.det(grad)
                kind = 'saddle' if det < 0 else ('sink' if trace < 0 else 'source')
                rates = np.abs(np.linalg.eigvals(grad).real)
                found.append(dict(hour_deg=round(float(ho), 1), lat_deg=round(float(la), 1), kind=kind,
                                  spiral=bool(trace ** 2 < 4 * det),
                                  divergence_per_day=round(float(trace * JULIAN_DAY), 3),
                                  e_folding_days=round(float(1.0 / (rates.min() * JULIAN_DAY)), 1)
                                  if rates.min() > 0 else None))
    return found


def releases(n_heights):
    lat0, hour0 = np.meshgrid(RELEASE_LAT_DEG, RELEASE_HOUR_DEG, indexing='ij')
    return (np.repeat(np.arange(n_heights), lat0.size), np.tile(lat0.ravel(), n_heights),
            np.tile(hour0.ravel(), n_heights))


def occupancy(res, layer, lat0, n_heights, second_half=True):
    """Share of parcel time in sectors of hour and bands of latitude, by height, each parcel weighted by the area
    its release stands for; over the second half of the run."""
    h, la = res['hours'], res['lats']
    if second_half:
        h, la = h[h.shape[0] // 2:], la[la.shape[0] // 2:]
    sector = np.clip(((((h + 180.0) % 360.0)) / 360.0 * OCCUPANCY_HOURS).astype(int), 0, OCCUPANCY_HOURS - 1)
    band = np.clip(((la + 90.0) / 180.0 * OCCUPANCY_LATS).astype(int), 0, OCCUPANCY_LATS - 1)
    w = np.broadcast_to(np.cos(np.radians(lat0))[None, :], h.shape)
    out = np.zeros((n_heights, OCCUPANCY_LATS, OCCUPANCY_HOURS))
    lay = np.broadcast_to(layer[None, :], h.shape)
    np.add.at(out, (lay.ravel(), band.ravel(), sector.ravel()), w.ravel())
    return out / out.sum(axis=(1, 2), keepdims=True)


def area_shares():
    edges = np.radians(np.linspace(-90.0, 90.0, OCCUPANCY_LATS + 1))
    return 0.5 * (np.sin(edges[1:]) - np.sin(edges[:-1])) / OCCUPANCY_HOURS


def gathering(res, groups, radius_m, cell_km, record_days, days=GATHER_DAYS):
    """How close parcels drifting together come: at each of the days within the run, the median distance (km) from
    each parcel to its nearest neighbour in its group (a height, or a height and a release) and the share of parcels
    with a neighbour within one model cell, by height. groups[i] is (height index, group label) of parcel i."""
    heights = np.unique(groups[:, 0])
    out = dict(days=[], nearest_km_p50=[], within_cell_share=[])
    for day in days:
        k = int(round(day / record_days))
        if k >= res['hours'].shape[0]:
            break
        q = position(res['lats'][k].astype(float), res['hours'][k].astype(float))
        near, share = [], []
        for h in heights:
            dist, cell = [], []
            for g in np.unique(groups[groups[:, 0] == h, 1]):
                m = (groups[:, 0] == h) & (groups[:, 1] == g)
                angle = np.arccos(np.clip(q[m] @ q[m].T, -1.0, 1.0))
                np.fill_diagonal(angle, np.inf)
                d = angle.min(axis=1) * radius_m[m] / 1e3
                dist.append(d)
            d = np.concatenate(dist)
            near.append(float(np.median(d)))
            share.append(float(np.mean(d < cell_km)))
        out['days'].append(day)
        out['nearest_km_p50'].append(near)
        out['within_cell_share'].append(share)
    return dict(days=out['days'], nearest_km_p50=r(out['nearest_km_p50'], 1),
                within_cell_share=r(out['within_cell_share'], 3))


def per_parcel(res, run_days, lat0, hour0, record_days, n_heights):
    """Per-parcel results reshaped [height][release lat][release hour]."""
    shape = (n_heights, RELEASE_LAT_DEG.size, RELEASE_HOUR_DEG.size)
    net = res['end_hour'] - hour0
    circuits = np.abs(net) / 360.0
    with np.errstate(divide='ignore', invalid='ignore'):
        mean_day = np.where(circuits >= 1.0, np.sign(net) * run_days / circuits, np.nan)
    out = dict(first_pass_days=r(res['first_pass_days'].reshape(shape), 1),
               first_pass_sign=r(res['first_pass_sign'].reshape(shape), 0),
               mean_day_days=r(mean_day.reshape(shape), 1),
               dwell_days=r(longest_dwell(res['hours'], record_days).reshape(shape), 0),
               sun_up_share=r(res['sun_up_share'].reshape(shape), 3), light=r(res['light'].reshape(shape), 3),
               light_ratio=r((res['light'] / res['fixed_light']).reshape(shape), 3),
               sunsets=r(res['sunsets'].reshape(shape), 0), western_sunrises=r(res['western_sunrises'].reshape(shape), 0),
               end_lat_deg=r(res['end_lat'].reshape(shape), 1),
               end_hour_deg=r(((res['end_hour'] + 180.0) % 360.0 - 180.0).reshape(shape), 1),
               near_pole_steps=int(res['jumped'].sum()))
    if res['over_ground_share'] is not None:
        out['over_ground_share'] = r(res['over_ground_share'].reshape(shape), 3)
    return out


def release_summary(res, n_releases, n_heights, record_days):
    """First-year results of every release, summarised by height and release latitude over the releases and the
    release hours."""
    shape = (n_releases, n_heights, RELEASE_LAT_DEG.size, RELEASE_HOUR_DEG.size)
    first = np.where(np.isnan(res['first_pass_days']), np.inf, res['first_pass_days']).reshape(shape)
    sign = res['first_pass_sign'].reshape(shape)
    dwell = longest_dwell(res['hours'], record_days).reshape(shape)
    pool = lambda a: np.moveaxis(a, 0, -2).reshape(n_heights, RELEASE_LAT_DEG.size, -1)
    fp = np.percentile(pool(first), (10, 50, 90), axis=-1)
    return dict(first_pass_p10_days=r(fp[0], 1), first_pass_p50_days=r(fp[1], 1), first_pass_p90_days=r(fp[2], 1),
                no_pass_share=r((~np.isfinite(pool(first))).mean(axis=-1), 3),
                reversed_share=r((pool(sign) < 0).mean(axis=-1), 3),
                dwell_p50_days=r(np.median(pool(dwell), axis=-1), 0), dwell_max_days=r(pool(dwell).max(axis=-1), 0),
                sun_up_share_p50=r(np.median(pool(res['sun_up_share'].reshape(shape)), axis=-1), 3),
                light_ratio_p50=r(np.median(pool((res['light'] / res['fixed_light']).reshape(shape)), axis=-1), 3),
                western_sunrise_share=r((pool(res['western_sunrises'].reshape(shape)) > 0).mean(axis=-1), 3),
                western_sunrises_per_sunset=r(pool(res['western_sunrises'].reshape(shape)).sum(axis=-1)
                                              / np.maximum(pool(res['sunsets'].reshape(shape)).sum(axis=-1), 1), 3),
                over_ground_share_mean=r(pool(res['over_ground_share'].reshape(shape)).mean(axis=-1), 3),
                near_pole_steps=int(res['jumped'].sum()))


def route_array(res, every):
    """Daily routes of the first ROUTE_DAYS as int16 hundredths of a degree: (record, parcel, [hour, lat])."""
    n = ROUTE_DAYS * every + 1
    h = (res['hours'][:n:every] + 180.0) % 360.0 - 180.0
    return np.round(np.stack([h, res['lats'][:n:every]], axis=-1) * 100.0).astype(np.int16)


def trajectories(steady_u, steady_v, fill_u, fill_v, ground, track, lat, lon, dt_out, outputs_per_year):
    """Parcels carried through the steady Sun-frame flow and through the 3-day means in sequence."""
    nt = len(TRAJECTORY_KM)
    heights_m = np.asarray(TRAJECTORY_KM) * 1e3
    layer, lat0, hour0 = releases(nt)
    radius = MOON_RADIUS + heights_m[layer]
    year_s = outputs_per_year * dt_out
    record_days = RECORD_S / JULIAN_DAY
    every = int(round(JULIAN_DAY / RECORD_S))
    hours = hour_centres()
    xyz, rows = cartesian(steady_u, steady_v, lat, hours)
    cell_km = MOON_RADIUS * np.radians(lon[1] - lon[0]) / 1e3
    alone = np.stack([layer, np.zeros_like(layer)], axis=1)

    def steady(t, la, ho):
        u, v = interpolate(xyz, rows, hours[0], 360.0 / HOUR_BINS, la, ho, lead=(layer,))
        return u + SUN_RATE * radius * np.cos(np.radians(la)), v

    print('steady routes ...', flush=True)
    res = integrate(steady, lat0, hour0, radius, STEADY_YEARS * year_s)
    steady_out = dict(years=STEADY_YEARS, **per_parcel(res, STEADY_YEARS * year_s / JULIAN_DAY, lat0, hour0,
                                                       record_days, nt),
                      occupancy=r(occupancy(res, layer, lat0, nt), 4),
                      gathering=gathering(res, alone, radius, cell_km, record_days))
    routes = dict(steady_routes=route_array(res, every))
    del res
    fixed = {f'{k:g}': fixed_points(steady_u[i] + sun_speed(lat, heights_m[i])[:, None], steady_v[i], lat, hours,
                                    MOON_RADIUS + heights_m[i]) for i, k in enumerate(TRAJECTORY_KM)}

    n_out = fill_u.shape[0]
    seq = np.empty((n_out, nt, 3, lat.size + 2, lon.size), dtype=np.float32)
    for n in range(n_out):
        seq[n], seq_rows = cartesian(fill_u[n], fill_v[n], lat, lon)
    step = lon[1] - lon[0]
    ground_mean = ground.mean(axis=0)

    def frame(t):
        x = t / dt_out - 0.5
        n0 = np.clip(np.floor(x).astype(int), 0, n_out - 2)
        return n0, np.clip(x - n0, 0.0, 1.0)

    def make(lay, rad):
        def velocity(t, la, ho):
            n0, w = frame(t)
            u0, v0 = interpolate(seq, seq_rows, lon[0], step, la, ho + track[n0], lead=(n0, lay))
            u1, v1 = interpolate(seq, seq_rows, lon[0], step, la, ho + track[n0 + 1], lead=(n0 + 1, lay))
            return (1 - w) * u0 + w * u1 + SUN_RATE * rad * np.cos(np.radians(la)), (1 - w) * v0 + w * v1

        def over(t, la, ho):
            n0, w = frame(t)
            lo = ho + track[n0] + w * (track[n0 + 1] - track[n0])
            i = np.abs(lat[None, :] - la[:, None]).argmin(axis=1)
            j = np.round(((lo - lon[0]) % 360.0) / step).astype(int) % lon.size
            return ground_mean[i, j] > heights_m[lay]
        return velocity, over

    n_releases = n_out // outputs_per_year
    rel = np.repeat(np.arange(n_releases), layer.size)
    lay_all = np.tile(layer, n_releases)
    velocity, over = make(lay_all, np.tile(radius, n_releases))
    print(f'routes through the 3-day means: {n_releases} releases ...', flush=True)
    res = integrate(velocity, np.tile(lat0, n_releases), np.tile(hour0, n_releases), np.tile(radius, n_releases),
                    year_s, t0=rel * year_s, over=over)
    first_year = release_summary(res, n_releases, nt, record_days)
    first_year['occupancy'] = r(occupancy(res, lay_all, np.tile(lat0, n_releases), nt, second_half=False), 4)
    first_year['gathering'] = gathering(res, np.stack([lay_all, rel], axis=1), np.tile(radius, n_releases), cell_km,
                                        record_days)
    del res
    velocity, over = make(layer, radius)
    print('the first release through all years ...', flush=True)
    res = integrate(velocity, lat0, hour0, radius, n_releases * year_s, over=over)
    long_run = dict(years=n_releases, **per_parcel(res, n_releases * year_s / JULIAN_DAY, lat0, hour0, record_days, nt),
                    occupancy=r(occupancy(res, layer, lat0, nt), 4),
                    gathering=gathering(res, alone, radius, cell_km, record_days))
    routes['sequence_routes'] = route_array(res, every)
    k = int(round(POSITION_EVERY_DAYS / record_days))
    h = (res['hours'][::k] + 180.0) % 360.0 - 180.0
    routes['long_run_positions'] = np.round(np.stack([h, res['lats'][::k]], axis=-1) * 100.0).astype(np.int16)
    out = dict(heights_km=list(TRAJECTORY_KM), release_lat_deg=RELEASE_LAT_DEG.tolist(),
               release_hour_deg=RELEASE_HOUR_DEG.tolist(), step_s=STEP_S, record_s=RECORD_S,
               dwell_span_deg=DWELL_SPAN_DEG, dwell_lengths_days=list(DWELL_DAYS), route_days=ROUTE_DAYS,
               cell_km=round(float(cell_km), 1), position_every_days=POSITION_EVERY_DAYS,
               occupancy_hour_edges_deg=np.linspace(-180, 180, OCCUPANCY_HOURS + 1).tolist(),
               occupancy_lat_edges_deg=np.linspace(-90, 90, OCCUPANCY_LATS + 1).tolist(),
               area_share=r(area_shares(), 4), fixed_points=fixed, steady=steady_out,
               sequence=dict(releases=n_releases, first_year=first_year, long_run=long_run))
    return out, routes


def zonal_winds(folder: str, first: int, last: int, gas_constant: float, dt_out: float):
    from climate.gcm.exoplasim_run import RUNS
    totals, fills_u, fills_v, grounds, subs, paths = None, [], [], [], [], []
    for year in range(first, last + 1):
        path = RUNS / folder / 'model' / f'MOST.{year:05d}.nc'
        y = read_year(path, gas_constant)
        if totals is None:
            lat, lon = y['lat'], y['lon']
            totals = Totals(lat.size, lon.size, len(TRAJECTORY_KM))
        totals.add(y)
        fills_u.append(y['fill_u'].transpose(1, 0, 2, 3))
        fills_v.append(y['fill_v'].transpose(1, 0, 2, 3))
        grounds.append(y['ground'])
        subs.append(y['sub'])
        paths.append(path)
        outputs_per_year = y['sub'].size
        del y
        print(f'year {year} gathered', flush=True)
    stats, arrays = summarise(totals, lat)
    sub = np.concatenate(subs)
    step_deg = 360.0 * dt_out / (SYNODIC_MONTH_DAYS * JULIAN_DAY)
    track = sun_track(sub, step_deg)
    residual = (sub - track + 180.0) % 360.0 - 180.0
    with np.errstate(invalid='ignore', divide='ignore'):
        steady_u, steady_v = totals.tsu / totals.tn, totals.tsv / totals.tn
    traj, routes = trajectories(steady_u, steady_v, np.concatenate(fills_u), np.concatenate(fills_v),
                                np.concatenate(grounds), track, lat, lon, dt_out, outputs_per_year)
    sun = dict(synodic_month_days=SYNODIC_MONTH_DAYS, output_interval_s=dt_out,
               track_step_deg=round(step_deg, 5), track_start_lon_deg=round(float(track[0]), 3),
               track_residual_max_deg=round(float(np.abs(residual).max()), 3))
    axes = dict(height_km=r(HEIGHTS_KM, 1), lat_deg=r(lat, 3), hour_angle_deg=r(hour_centres(), 1))
    arrays.update(routes, steady_u_m_s=steady_u, steady_v_m_s=steady_v)
    return dict(axes=axes, sun=sun, **stats, trajectories=traj, samples=int(totals.n.sum())), arrays, paths


INNER = re.compile(r'\[([^\[\]{}"]*)\]')


def compact(product) -> str:
    """JSON with every innermost list of numbers on one line."""
    text = json.dumps(product, indent=1)
    return INNER.sub(lambda m: '[' + ','.join(p.strip() for p in m.group(1).split(',')) + ']'
                     if m.group(1).strip() else '[]', text) + '\n'


def main(argv=None) -> int:
    args = list(argv or sys.argv[1:])
    folder_out = RESULTS
    if '--out' in args:                                      # trial runs write elsewhere
        k = args.index('--out')
        folder_out = Path(args[k + 1])
        del args[k:k + 2]
    if len(args) != 1 or ':' not in args[0]:
        raise SystemExit('Give a run and year range, for example A28_dim5_moon:20-29 [--out folder]')
    folder, span = args[0].split(':')
    first, last = (int(v) for v in span.split('-'))
    from atmosphere.radiative_convective.thermodynamics import earthlike_air
    from climate.gcm.exoplasim_run import RUNS
    from shared.provenance import constants_used
    progress = json.loads((RUNS / folder / 'progress.json').read_text())
    model = progress['configuration']['model']
    dt_out = model['steps_per_lunar_day'] / model['writes_per_lunar_day'] * model['timestep_min'] * 60.0
    air = earthlike_air(progress['configuration']['pressure_pa'], 400.0)
    data, arrays, paths = zonal_winds(folder, first, last, air.gas_constant, dt_out)
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    files = ['climate/gcm/zonal_winds.py', 'climate/gcm/global_winds.py', 'climate/gcm/site_winds.py',
             'climate/gcm/climatology.py']
    inputs = {f'climate/gcm/runs/{folder}/model/{p.name}': digest(p) for p in paths}
    product = dict(schema=SCHEMA, run=folder, years=[first, last], configuration=progress['configuration'],
                   producer=dict(domain='climate', files={f: digest(ROOT / f) for f in files}, inputs=inputs,
                                 constants=constants_used(files)),
                   air_gas_constant_j_kg_k=round(air.gas_constant, 4), evidence=EVIDENCE, reading_rule=READING_RULE,
                   **data)
    folder_out.mkdir(parents=True, exist_ok=True)
    out = folder_out / f'zonal_winds_{folder}.json'
    out.write_text(compact(product))
    meta = dict(schema=SCHEMA, run=folder, years=[first, last], json=out.name, reading_rule=READING_RULE,
                evidence=EVIDENCE)
    store = {k: (v.astype(np.int32) if k.endswith('samples') else v if v.dtype.kind == 'i' else
                 v.astype(np.float32) if k.startswith('steady_') else v.astype(np.float16))
             for k, v in arrays.items() if not k.startswith(('sun_u_mean', 'sun_v_mean'))}
    np.savez_compressed(out.with_suffix('.npz'), metadata=json.dumps(meta), **store)
    w = np.cos(np.radians(np.array(data['axes']['lat_deg'])))
    for k, z in enumerate(data['axes']['height_km']):
        u = np.array(data['ground_frame']['u_mean_m_s'][k], dtype=float)
        hold = np.array(data['holding']['hold_1'][k], dtype=float)
        back = np.array(data['holding']['reversed'][k], dtype=float)
        ok = np.isfinite(u)
        print(f"{z:5.1f} km  u {np.nansum(u * w) / w[ok].sum():6.2f} m/s  hold(1 m/s) "
              f"{np.nansum(hold * w) / w[ok].sum():.3f}  reversed {np.nansum(back * w) / w[ok].sum():.3f}  "
              f"median day at the equator {data['solar_day']['p50_days'][k][15]} d")
    print(out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
