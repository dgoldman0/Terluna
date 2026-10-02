"""What the CM1 storms hold for charging: graupel, cloud ice and snow, the supercooled cloud water among them, and
the charging zone where they meet.

    climate/gcm/.venv/bin/python -m climate.crm.mixed_phase_analysis box_0e --from-day 29.5

writes ../results/crm/mixed_phase_<case>.json.

The Morrison scheme carries the mass and number of each ice type and gives each an exponential size distribution of
fixed density (cloud ice 500, snow 100 and graupel 400 kg/m3), its slope held within the scheme's bounds. Its fall
speeds V = a D^b (rho_su/rho)^c take the host gravity through a = a_Earth (g/g_Earth)^((b+1)/3), as cm1_run.py
patches them; the build with Earth's fall speeds keeps Earth's. From the snapshots (every 3 model hours) over the
span the analysis gives:
- the charging zone, the air from 0 to -40 C where graupel, cloud ice or snow, and supercooled cloud water each
  reach a threshold mixing ratio: its depth per unit of ground by temperature band, how often it is present, its
  course through the lunar day, and the width and depth of its patches;
- in the zone, by temperature band, the mass, number, mass-weighted diameter and mass-weighted fall speed of
  graupel, cloud ice and snow, the cloud water, the air's density and the updrafts;
- the graupel's bulk residence time aloft: its mass above the melting level over the rate it falls through it;
- how deep the storms reach: the cloud tops and the heights of the 0 and -40 C levels.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from math import gamma
from pathlib import Path
import sys
import numpy as np
from climate.crm import ring_analysis as ra
from climate.crm import ring_comfort as rcf
from climate.crm import ring_crossings as rc

RESULTS = ra.RESULTS
SCHEMA = 'terluna.climate.crm-mixed-phase/1'
G_EARTH = 9.81                                    # the gravity Morrison's fall-speed laws were fitted for
RHOSU = 85000.0 / (287.04 * 273.15)               # Morrison's reference air density for fall speeds (kg/m3)
RD, CP = 287.04, 1005.7
# Morrison's ice types with ihail = 0 (morrison.F): mass and number fields, V = a D^b (m/s, D in m) with the air
# density correction (rho_su/rho)^c, particle density (kg/m3), the bounds on the slope lambda (1/m) and the cap on the
# mass-weighted fall speed (m/s, before the gravity and air density factors).
SPECIES = {
    'ice': dict(q='qi', n='nci', a=700.0, b=1.0, c=0.35, rho=500.0, lam=(1.0 / 350.0e-6, 1.0 / 1.0e-6), cap=1.2),
    'snow': dict(q='qs', n='ncs', a=11.72, b=0.41, c=0.54, rho=100.0, lam=(1.0 / 2000.0e-6, 1.0 / 10.0e-6), cap=1.2),
    'graupel': dict(q='qg', n='ncg', a=19.3, b=0.37, c=0.54, rho=400.0, lam=(1.0 / 2000.0e-6, 1.0 / 20.0e-6), cap=20.0)}
THRESHOLDS = (1.0e-6, 1.0e-5, 1.0e-4)             # kg/kg that graupel, ice or snow, and cloud water must each reach
MAIN = 1.0e-5
WARMEST_C, COLDEST_C = 0.0, -40.0
BANDS_C = ((0.0, -10.0), (-10.0, -20.0), (-20.0, -30.0), (-30.0, -40.0))
CLOUD_TOP_KG_KG = 1.0e-5                          # condensate that counts toward the cloud top
HOUR_BINS = 12
EVIDENCE = ('CM1 r22.0 at lunar gravity with Morrison two-moment microphysics (graupel), its fall speeds scaled for '
            'the host gravity; the particles are the scheme\'s bulk populations with exponential size distributions '
            'and fixed densities, sampled every 3 model hours. Patch widths are limited by the grid spacing, and '
            'quantities that change within hours (a storm cell\'s life, particle residence) are sampled coarsely.')
READING_RULE = ('The charging zone is the air from 0 to -40 C holding at least the threshold mixing ratio of graupel, '
                'of cloud ice and snow together, and of cloud water. Depths per unit of ground are the zone\'s volume '
                'over the domain\'s ground area, averaged over the span\'s snapshots; a ring\'s columns stand for a '
                'unit width. Band values are volume-weighted means over the zone cells in each band; diameters and fall '
                'speeds are mass-weighted means of each cell\'s distribution, weighted again by the cell\'s mass. '
                'Residence is the mean graupel mass above the melting level over the mean graupel flux through the '
                'lowest level colder than 0 C.')


def fall_gravity(record: dict) -> float:
    """The gravity the build's fall speeds follow: Earth's in the build with Earth's fall speeds, the host's otherwise."""
    from climate.crm.cm1_run import BUILD_PATCHES, EARTH_FALL_PATCH
    label = record['build']['label']
    return G_EARTH if EARTH_FALL_PATCH in BUILD_PATCHES.get(label, []) else float(record['build']['gravity_m_s2'])


def gravity_factor(species: dict, g_fall: float) -> float:
    """How the patch scales a species' fall-speed coefficient: (g/g_Earth)^((b+1)/3)."""
    return (g_fall / G_EARTH) ** ((species['b'] + 1.0) / 3.0)


def slope(q, n, species: dict):
    """The exponential distribution's slope lambda (1/m), held within the scheme's bounds, and the number (per kg)
    that goes with it, as Morrison adjusts it; nan where there is no mass."""
    q, n = np.asarray(q, float), np.asarray(n, float)
    ok = (q > 0.0) & (n > 0.0)
    lam = np.full(q.shape, np.nan)
    lam[ok] = (np.pi * species['rho'] * n[ok] / q[ok]) ** (1.0 / 3.0)
    lam = np.clip(lam, *species['lam'])
    return lam, np.where(ok, lam ** 3 * q / (np.pi * species['rho']), 0.0)


def fall_speed(lam, rho_air, species: dict, g_fall: float, moment: str = 'mass'):
    """Mass- or number-weighted fall speed (m/s) of an exponential distribution with slope lam, at air density
    rho_air, for fall speeds at gravity g_fall, with the scheme's cap on the mass-weighted speed."""
    k = gravity_factor(species, g_fall) * (RHOSU / np.asarray(rho_air, float)) ** species['c']
    b = species['b']
    if moment == 'mass':
        return np.minimum(k * species['a'] * gamma(4.0 + b) / 6.0 / lam ** b, species['cap'] * k)
    return k * species['a'] * gamma(1.0 + b) / lam ** b


def speed_at(d, rho_air, species: dict, g_fall: float):
    """Fall speed (m/s) of a particle of diameter d (m)."""
    return gravity_factor(species, g_fall) * (RHOSU / np.asarray(rho_air, float)) ** species['c'] * species['a'] * \
        np.asarray(d, float) ** species['b']


def zone_mask(t_c, qc, qg, qi, qs, threshold: float):
    """The charging zone: 0 to -40 C with graupel, cloud ice and snow together, and cloud water at the threshold."""
    return (t_c < WARMEST_C) & (t_c >= COLDEST_C) & (qc >= threshold) & (qg >= threshold) & (qi + qs >= threshold)


def periodic_labels(mask):
    """Connected patches of a periodic 2-D (or 1-D) mask, joined across the edges: an integer label per cell, 0
    outside the patches."""
    from scipy import ndimage
    mask = np.atleast_2d(mask)
    labels, count = ndimage.label(mask)
    parent = list(range(count + 1))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for first, last in ((labels[:, 0], labels[:, -1]), (labels[0, :], labels[-1, :])):
        for a, b in zip(first, last):
            if a and b:
                ra_, rb_ = find(a), find(b)
                if ra_ != rb_:
                    parent[max(ra_, rb_)] = min(ra_, rb_)
    roots = np.array([find(a) for a in range(count + 1)])
    return roots[labels]


def hour_angles(record: dict, geo: dict, t_s: float, day_s: float, ncol: int):
    """Local hour angle (degrees, 0 = noon) of every column: the box's site for a box, each column's for a ring."""
    cfg = record['configuration']
    if cfg.get('kind') == 'box':
        h = cfg['start_hour_angle_deg'] + cfg['site']['lon_deg'] + 360.0 * t_s / day_s
        return np.full(ncol, (h + 180.0) % 360.0 - 180.0)
    return ra.hour_angle(geo, t_s, day_s)


def crossing_height(zh, t_c, value: float):
    """The height (m) where a temperature profile first falls through value going up, by linear interpolation;
    nan if it never does."""
    below = np.flatnonzero(t_c < value)
    if not below.size or below[0] == 0:
        return float('nan')
    k = below[0]
    return float(zh[k - 1] + (zh[k] - zh[k - 1]) * (t_c[k - 1] - value) / (t_c[k - 1] - t_c[k]))


def stats(values) -> dict:
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if not v.size:
        return dict(n=0)
    return dict(n=int(v.size), mean=float(v.mean()), p10=float(np.percentile(v, 10)), median=float(np.median(v)),
                p90=float(np.percentile(v, 90)), max=float(v.max()))


def analyse(name: str, from_day: float, to_day: float = np.inf) -> dict:
    from climate.crm.cm1_run import solar_day_s
    case = ra.RUNS / name
    record = json.loads((case / 'case.json').read_text())
    geo = ra.case_geometry(case)
    zh, zw = geo['zh'], geo['zw']
    dz = np.diff(zw)
    grid = record['grid']
    nx, ny = grid['nx'], grid.get('ny', 1)
    ncol = nx * ny
    output_s, day_s = record['configuration']['output_s'], solar_day_s()
    g_fall = fall_gravity(record)
    nxy, layout = ra.grads_layout(case / 'cm1out_s.ctl')
    offsets, offset = {}, 0
    for field, levels in layout:
        offsets[field] = (offset, levels)
        offset += levels * nxy
    numbers = rcf.snapshots(case, output_s, from_day * 86400.0, to_day * 86400.0)

    depth = {t: [] for t in THRESHOLDS}                                    # zone depth per unit of ground, per snapshot
    band_depth = {t: np.zeros(len(BANDS_C)) for t in THRESHOLDS}
    hour_sum, hour_count = np.zeros(HOUR_BINS), np.zeros(HOUR_BINS)
    keys = ('vol', 'rho', 'p', 'z', 'qc', 'w', 'mg', 'ng', 'mi', 'ni', 'ms', 'ns', 'dg', 'di', 'ds', 'vg', 'vi', 'vs',
            'dv_i', 'vol_i', 'dv_s', 'vol_s')
    acc = {k: np.zeros(len(BANDS_C)) for k in keys}
    w_zone, zone_columns, widths, depths_m = [], 0, [], []
    tops, z0c, z40c = [], [], []
    p_sum, t_sum, rho_sum = np.zeros(len(zh)), np.zeros(len(zh)), np.zeros(len(zh))
    mass_aloft, flux_melt, present = [], [], 0
    for n in numbers:
        raw = np.memmap(case / f'cm1out_t{n:06d}_s.dat', dtype='<f4', mode='r')

        def get(field):
            o, levels = offsets[field]
            return np.array(raw[o:o + levels * nxy], float).reshape(levels, nxy)
        th, prs, qv = get('th'), get('prs'), get('qv')
        t = th * (prs / 1.0e5) ** (RD / CP)
        t_c = t - 273.15
        rho = prs / (RD * t * (1.0 + 0.61 * qv))
        qc, qg, qi, qs = get('qc'), get('qg'), get('qi'), get('qs')
        w = get('winterp')
        cond = qc + get('qr') + qg + qi + qs
        cloudy = cond >= CLOUD_TOP_KG_KG
        if cloudy.any():
            tops.append(float(zh[np.max(np.where(cloudy.any(axis=1)))]))
        mean_t = t_c.mean(axis=1)
        p_sum += prs.mean(axis=1)
        t_sum += mean_t
        rho_sum += rho.mean(axis=1)
        z0c.append(crossing_height(zh, mean_t, 0.0))
        z40c.append(crossing_height(zh, mean_t, -40.0))
        t_s = (n - 1) * output_s
        hours = hour_angles(record, geo, t_s, day_s, ncol)
        for thr in THRESHOLDS:
            zone = zone_mask(t_c, qc, qg, qi, qs, thr)
            col_depth = (zone * dz[:, None]).sum(axis=0)
            depth[thr].append(float(col_depth.mean()))
            for b, (warm, cold) in enumerate(BANDS_C):
                band_depth[thr][b] += float(((zone & (t_c < warm) & (t_c >= cold)) * dz[:, None]).sum() / ncol)
            if thr != MAIN:
                continue
            k = np.clip(np.digitize(hours, np.linspace(-180.0, 180.0, HOUR_BINS + 1)) - 1, 0, HOUR_BINS - 1)
            np.add.at(hour_sum, k, col_depth)
            np.add.at(hour_count, k, 1.0)
            if not zone.any():
                continue
            present += 1
            lam, num, mix = {}, {}, {}
            for key, sp in SPECIES.items():
                mix[key] = np.where(zone, get(sp['q']), 0.0)
                lam[key], num[key] = slope(mix[key], np.where(zone, get(sp['n']), 0.0), sp)
            vol = dz[:, None] * np.ones((1, nxy))
            for b, (warm, cold) in enumerate(BANDS_C):
                m = zone & (t_c < warm) & (t_c >= cold)
                if not m.any():
                    continue
                v = vol[m]
                r = rho[m]
                acc['vol'][b] += v.sum()
                acc['rho'][b] += (r * v).sum()
                acc['p'][b] += (prs[m] * v).sum()
                acc['z'][b] += ((zh[:, None] * np.ones((1, nxy)))[m] * v).sum()
                acc['qc'][b] += (r * qc[m] * v).sum()
                acc['w'][b] += (w[m] * v).sum()
                speeds = {}
                for key, mk, nk, dk, vk in (('graupel', 'mg', 'ng', 'dg', 'vg'), ('ice', 'mi', 'ni', 'di', 'vi'),
                                            ('snow', 'ms', 'ns', 'ds', 'vs')):
                    sp = SPECIES[key]
                    mass = r * mix[key][m]                                  # kg/m3
                    ok = mass > 0.0
                    lm = np.where(ok, lam[key][m], 1.0)
                    speeds[key] = np.where(ok, fall_speed(lm, r, sp, g_fall), np.nan)
                    acc[mk][b] += (mass * v).sum()
                    acc[nk][b] += (r * num[key][m] * v).sum()
                    acc[dk][b] += (np.where(ok, 4.0 / lm, 0.0) * mass * v).sum()
                    acc[vk][b] += (np.where(ok, speeds[key], 0.0) * mass * v).sum()
                for other, dk, vk in (('ice', 'dv_i', 'vol_i'), ('snow', 'dv_s', 'vol_s')):
                    both = np.isfinite(speeds['graupel']) & np.isfinite(speeds[other])
                    acc[dk][b] += (np.abs(speeds['graupel'] - speeds[other])[both] * v[both]).sum()
                    acc[vk][b] += v[both].sum()
            w_zone.append(w[zone])
            columns = zone.any(axis=0)
            zone_columns += int(columns.sum())
            labels = periodic_labels(columns.reshape(ny, nx))
            top = np.where(zone, zw[1:, None], -np.inf).max(axis=0)
            bottom = np.where(zone, zw[:-1, None], np.inf).min(axis=0)
            for patch in np.unique(labels[labels > 0]):
                cells = (labels == patch).ravel()
                area = cells.sum() * grid['dx_m'] * (grid['dx_m'] if ny > 1 else 1.0)
                widths.append(2.0 * np.sqrt(area / np.pi) if ny > 1 else cells.sum() * grid['dx_m'])
                depths_m.append(float(top[cells].max() - bottom[cells].min()))
        # graupel above the melting level and its flux through the lowest level colder than 0 C
        lam_g, _ = slope(qg, get('ncg'), SPECIES['graupel'])
        vg_all = np.where(qg > 0.0, fall_speed(lam_g, rho, SPECIES['graupel'], g_fall), 0.0)
        cold = t_c < 0.0
        k0 = np.argmax(cold, axis=0)
        has = cold.any(axis=0)
        cols = np.arange(nxy)
        mass_aloft.append(float(np.where(cold, rho * qg * dz[:, None], 0.0).sum() / ncol))
        flux_melt.append(float(np.where(has, rho[k0, cols] * qg[k0, cols] * vg_all[k0, cols], 0.0).sum() / ncol))

    count = max(1, len(numbers))
    vol = acc['vol']
    with np.errstate(invalid='ignore', divide='ignore'):
        bands = []
        for b, (warm, cold) in enumerate(BANDS_C):
            v = vol[b]
            if v <= 0.0:
                bands.append(dict(warm_c=warm, cold_c=cold, depth_m=0.0))
                continue
            row = dict(warm_c=warm, cold_c=cold, depth_m=float(band_depth[MAIN][b] / count),
                       height_km=float(acc['z'][b] / v / 1000.0), pressure_hpa=float(acc['p'][b] / v / 100.0),
                       air_density_kg_m3=float(acc['rho'][b] / v), cloud_water_g_m3=float(acc['qc'][b] / v * 1000.0),
                       updraft_m_s=float(acc['w'][b] / v))
            for key, mk, nk, dk, vk in (('graupel', 'mg', 'ng', 'dg', 'vg'), ('ice', 'mi', 'ni', 'di', 'vi'),
                                        ('snow', 'ms', 'ns', 'ds', 'vs')):
                row[key] = dict(mass_g_m3=float(acc[mk][b] / v * 1000.0), number_per_litre=float(acc[nk][b] / v / 1000.0),
                                diameter_mm=float(acc[dk][b] / acc[mk][b] * 1000.0) if acc[mk][b] > 0 else None,
                                fall_speed_m_s=float(acc[vk][b] / acc[mk][b]) if acc[mk][b] > 0 else None)
            row['graupel_speed_over_ice_m_s'] = float(acc['dv_i'][b] / acc['vol_i'][b]) if acc['vol_i'][b] > 0 else None
            row['graupel_speed_over_snow_m_s'] = float(acc['dv_s'][b] / acc['vol_s'][b]) if acc['vol_s'][b] > 0 else None
            bands.append(row)
    w_all = np.concatenate(w_zone) if w_zone else np.zeros(0)
    melt = float(np.mean(flux_melt)) if flux_melt else 0.0
    return dict(
        schema=SCHEMA, case=name, purpose=record['configuration'].get('purpose'),
        build=dict(label=record['build']['label'], gravity_m_s2=record['build']['gravity_m_s2'], fall_speeds_at_m_s2=g_fall),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        span_days=[(numbers[0] - 1) * output_s / 86400.0, (numbers[-1] - 1) * output_s / 86400.0] if numbers else None,
        snapshots=len(numbers), columns=ncol, dx_m=grid['dx_m'], kind=record['configuration'].get('kind', 'ring'),
        fall_speed_factor={k: gravity_factor(sp, g_fall) for k, sp in SPECIES.items()},
        mean_profile=dict(z_km=(zh / 1000.0).tolist(), pressure_hpa=(p_sum / count / 100.0).tolist(),
                          temperature_c=(t_sum / count).tolist(), air_density_kg_m3=(rho_sum / count).tolist()),
        storm_depth=dict(cloud_top_km=stats(np.array(tops) / 1000.0), zero_c_km=float(np.nanmean(z0c) / 1000.0),
                         minus_40_c_km=float(np.nanmean(z40c) / 1000.0)),
        zone=dict(threshold_kg_kg=MAIN, share_of_snapshots=present / count,
                  share_of_columns=zone_columns / (count * ncol),
                  depth_m={f'{t:g}': float(np.mean(depth[t])) for t in THRESHOLDS},
                  depth_by_snapshot_m=stats(depth[MAIN]),
                  by_hour_angle=[dict(hour_angle_deg=float(-165.0 + 30.0 * h),
                                      depth_m=float(hour_sum[h] / hour_count[h]) if hour_count[h] else None)
                                 for h in range(HOUR_BINS)],
                  patches=dict(count=len(widths), width_km=stats(np.array(widths) / 1000.0),
                               depth_km=stats(np.array(depths_m) / 1000.0)),
                  updraft_m_s=stats(w_all), updraft_share_above_5_m_s=float((w_all > 5.0).mean()) if w_all.size else None,
                  bands=bands),
        graupel_aloft=dict(mass_kg_m2=float(np.mean(mass_aloft)) if mass_aloft else 0.0,
                           flux_through_melting_level_g_m2_h=melt * 3.6e6,
                           residence_h=float(np.mean(mass_aloft) / melt / 3600.0) if melt > 0 else None))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('case')
    parser.add_argument('--from-day', type=float, default=29.5)
    parser.add_argument('--to-day', type=float, default=np.inf)
    parser.add_argument('--output', type=Path, help='where to write the summary (default: the results folder)')
    args = parser.parse_args(argv)
    result = analyse(args.case, args.from_day, args.to_day)
    result['producer'] = dict(domain='climate', files={f'crm/{Path(__file__).name}':
                                                       hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]})
    path = args.output or RESULTS / f'mixed_phase_{args.case}.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rc.rounded(result), indent=1) + '\n')
    z = result['zone']
    print(f'{path}: {result["snapshots"]} snapshots; charging zone {z["depth_m"][f"{MAIN:g}"]:.0f} m deep per unit of '
          f'ground, present in {z["share_of_snapshots"]:.0%} of snapshots; graupel residence '
          f'{result["graupel_aloft"]["residence_h"]} h')
    return 0


if __name__ == '__main__':
    sys.exit(main())
