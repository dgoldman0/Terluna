"""The GCM beside a CM1 ring by height and local time: the air near the ground against the layers above it, the
ground, and the heat and vapour the ground gives the air, over the ring's land and its sea; python -m
climate.crm.gcm_profiles ring_45n_corrected --from-day 29.5.

The GCM's 3-day means over the rows of the ring's reference band, binned by local time into ten 36-degree bins as
ring_analysis bins its GCM composites, and the ring's snapshots after --from-day in the same bins. Each model's air is
averaged by mass over the GCM's lowest layers, their spans in heights above the ground taken from the GCM's own columns.
The GCM's land stands at its cells' own heights and the ring's at sea level, so over land the two compare through
differences in height above the ground: the air near the ground against the layers above it, the ground against the
air, and the fluxes. The 3-day means smooth the GCM's day a little more than the ring's binned snapshots."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'crm'
SCHEMA = 'terluna.climate.crm-gcm-profiles/1'
BINS = 10                          # 36 degrees of hour angle, as ring_analysis.gcm_equator_composites
LAYERS = 4                         # the GCM's lowest four layers, about 0-2, 2-7, 7-14 and 14-24 km up
RD, CP, LV = 287.04, 1005.7, 2.501e6


def bin_index(hour_deg):
    """Local-time bin (0 to BINS - 1) of hour angles in degrees, -180 to 180, noon 0."""
    return np.clip(((np.asarray(hour_deg, float) + 180.0) / (360.0 / BINS)).astype(int), 0, BINS - 1)


def interface_heights(ta, hus, sigma, gravity):
    """Heights (m) above the ground of PlaSim's layer interfaces for its lowest LAYERS layers, ground first: each
    layer's thickness rd Tv / g ln(sigma below / sigma above) from its virtual temperature. ta and hus are
    (..., lev), the top level first, as PlaSim writes them."""
    from climate.crm.cm1_run import plasim_half_levels
    edges = plasim_half_levels(sigma)                                   # top (0) first, the ground (1) last
    tv = np.asarray(ta, float) * (1.0 + 0.608 * np.asarray(hus, float))
    z = [np.zeros(tv.shape[:-1])]
    for m in range(LAYERS):
        level = len(sigma) - 1 - m                                       # the layers from the ground up
        z.append(z[-1] + RD * tv[..., level] / gravity * np.log(edges[level + 1] / edges[level]))
    return np.stack(z, axis=-1)


def layer_means(temp, prs, zh, spans):
    """Mass-weighted mean temperature of each column over each span of heights above the ground: temp and prs are
    (nz, nx) on heights zh (nz,); spans holds LAYERS + 1 interface heights. Returns (LAYERS, nx)."""
    mass = np.abs(np.gradient(np.asarray(prs, float), axis=0))
    out = []
    for m in range(LAYERS):
        inside = (zh >= spans[m]) & (zh < spans[m + 1])
        w = mass[inside]
        out.append((np.asarray(temp, float)[inside] * w).sum(axis=0) / w.sum(axis=0))
    return np.array(out)


def gcm_composite(run, years, latitude, band):
    """The GCM's 3-day means by local time over the land and the sea-level sea of the rows within `band` degrees of
    the latitude: ground, air near the ground and the lowest layers' temperatures (C), the layers' interface
    heights (m above the ground) and the ground's heat and vapour fluxes, upward positive (W/m2)."""
    import netCDF4
    from climate.crm.cm1_run import planet
    from climate.gcm.climatology import subsolar_longitudes
    gravity = planet()['gravity_m_s2']
    folder = HERE.parent / 'gcm' / 'runs' / run / 'model'
    names = ('ta', 'hus', 'ts', 'tas', 'hfss', 'hfls', 'czen')
    sums = {s: dict(ground_c=np.zeros(BINS), air_c=np.zeros(BINS), sensible_w_m2=np.zeros(BINS),
                    latent_w_m2=np.zeros(BINS), layers_c=np.zeros((BINS, LAYERS)),
                    heights_m=np.zeros((BINS, LAYERS + 1)))
            for s in ('land', 'water')}
    count = {s: np.zeros(BINS) for s in sums}
    for year in range(years[0], years[1] + 1):
        with netCDF4.Dataset(folder / f'MOST.{year:05d}.nc') as d:
            f = {k: np.asarray(d[k][:], float) for k in names}
            lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
            sigma = np.asarray(d['lev'][:], float)
            lsm = np.asarray(d['lsm'][:], float)[0] > 0.5
        topography = np.loadtxt(folder.parent / 'inputs' / 'moon_topography.sra', skiprows=1)
        height = topography.ravel().reshape(lsm.shape) / gravity
        rows = (np.abs(lat - latitude) < band)[:, None]
        sub = subsolar_longitudes(f['czen'], lat, lon)
        k = bin_index((lon[None, :] - sub[:, None] + 180.0) % 360.0 - 180.0)          # (time, lon)
        for surface, cells in (('land', rows & lsm), ('water', rows & ~lsm & (height < 1.0))):
            jj, ii = np.nonzero(cells)
            kk = k[:, ii]                                                               # (time, cells)
            ta, hus = (f[n][:, :, jj, ii].transpose(0, 2, 1) for n in ('ta', 'hus'))      # (time, cells, lev)
            z = interface_heights(ta, hus, sigma, gravity)
            values = dict(ground_c=f['ts'][:, jj, ii] - 273.15, air_c=f['tas'][:, jj, ii] - 273.15,
                          sensible_w_m2=-f['hfss'][:, jj, ii], latent_w_m2=-f['hfls'][:, jj, ii])
            for b in range(BINS):
                m = kk == b
                count[surface][b] += m.sum()
                for key, v in values.items():
                    sums[surface][key][b] += v[m].sum()
                sums[surface]['layers_c'][b] += (ta[m][:, ::-1][:, :LAYERS] - 273.15).sum(axis=0)
                sums[surface]['heights_m'][b] += z[m].sum(axis=0)
        land_height = float(height[rows & lsm].mean())
    out = {s: {key: (v / (count[s][:, None] if v.ndim == 2 else count[s])).tolist() for key, v in sums[s].items()}
           for s in sums}
    return dict(run=run, years=list(years), band_deg=band, rows_deg=[float(v) for v in lat[rows[:, 0]]],
                land_height_m=land_height, **out)


def ring_composite(name, from_day, spans):
    """The ring's snapshots after from_day by local time over its land and its sea in the same bins: ground and air
    at 2 m (C), each column's mass-weighted temperature over each GCM layer's span above the ground (spans: surface ->
    LAYERS + 1 heights in m), the ground's heat and vapour fluxes (W/m2, upward positive) and the mixed layer (km)."""
    from climate.crm import ring_analysis as ra
    from climate.crm.cm1_run import solar_day_s
    case = ra.RUNS / name
    geo = ra.case_geometry(case)
    tap = geo['record']['configuration']['output_s']
    outputs = sorted(int(p.name[8:14]) for p in case.glob('cm1out_t*_s.dat'))
    use = [n for n in outputs if (n - 1) * tap >= from_day * 86400.0]
    if not use:
        raise RuntimeError(f'no output after day {from_day}')
    zh, land, day_s = geo['zh'], geo['land'], solar_day_s()
    sums = {s: dict(ground_c=np.zeros(BINS), air_c=np.zeros(BINS), sensible_w_m2=np.zeros(BINS),
                    latent_w_m2=np.zeros(BINS), mixed_layer_km=np.zeros(BINS), layers_c=np.zeros((BINS, LAYERS)))
            for s in ('land', 'water')}
    count = {s: np.zeros(BINS) for s in sums}
    for n in use:
        d = ra.read_snapshot(case, n)
        k = bin_index(ra.hour_angle(geo, (n - 1) * tap, day_s))
        temp = d['th'] * (d['prs'] / 1.0e5) ** (RD / CP)
        rho = ra.air_density(d)
        values = dict(ground_c=d['tsk'] - 273.15, air_c=d['t2'] - 273.15, sensible_w_m2=d['thflux'] * rho * CP,
                      latent_w_m2=d['qvflux'] * rho * LV, mixed_layer_km=d['hpbl'] / 1000.0)
        for surface, sel in (('land', land), ('water', ~land)):
            layers = layer_means(temp[:, sel], d['prs'][:, sel], zh, spans[surface]) - 273.15    # (LAYERS, columns)
            kk = k[sel]
            for b in range(BINS):
                m = kk == b
                count[surface][b] += m.sum()
                for key, v in values.items():
                    sums[surface][key][b] += v[sel][m].sum()
                sums[surface]['layers_c'][b] += layers[:, m].sum(axis=1)
    out = {s: {key: (v / (count[s][:, None] if v.ndim == 2 else count[s])).tolist() for key, v in sums[s].items()}
           for s in sums}
    return dict(case=name, from_day=from_day, snapshots=len(use), **out)


def contrasts(side) -> dict:
    """Differences in height above the ground, by local time: the air near the ground and the lowest layer against
    the second layer (about 2-7 km up), and the ground against the air near it."""
    layers = np.asarray(side['layers_c'])
    air, ground = np.asarray(side['air_c']), np.asarray(side['ground_c'])
    return dict(air_minus_layer2_k=(air - layers[:, 1]).tolist(),
                layer1_minus_layer2_k=(layers[:, 0] - layers[:, 1]).tolist(),
                ground_minus_air_k=(ground - air).tolist())


def compare(name, from_day) -> dict:
    from climate.crm import ring_analysis as ra
    record = json.loads((ra.RUNS / name / 'case.json').read_text())
    if record.get('path'):
        raise ValueError('a circle of latitude only: a tilted ring crosses many GCM rows')
    gcm = gcm_composite(record['gcm']['run'], record['gcm']['years'], record['latitude_deg'],
                        record['reference']['band_deg'])
    spans = {s: np.asarray(gcm[s]['heights_m']).mean(axis=0) for s in ('land', 'water')}
    ring = ring_composite(name, from_day, spans)
    for side in (gcm, ring):
        for s in ('land', 'water'):
            side[s].update(contrasts(side[s]))
    return dict(
        schema=SCHEMA, case=name, latitude_deg=record['latitude_deg'],
        hour_angle_deg=list((np.arange(BINS) + 0.5) * 360.0 / BINS - 180.0),
        layer_spans_km={s: [round(float(v) / 1000.0, 3) for v in spans[s]] for s in spans},
        evidence=('The GCM\'s 3-day means and the ring\'s snapshots of the same forcing, side by side by local time '
                  'and height above the ground; neither model is independent of the other, since the ring is held to '
                  'the GCM\'s air from 8-16 km up and takes its vertical wind. The GCM\'s land stands at its cells\' '
                  'own heights (land_height_m) and the ring\'s at sea level, so the land\'s temperatures compare '
                  'through the contrasts in height above the ground, the sea\'s directly.'),
        reading_rule=('Ten 36-degree bins of local hour angle (noon 0, afternoon positive), centres in hour_angle_deg. '
                      'Temperatures in C; layers_c are the mass-weighted means over the GCM\'s lowest four layers, '
                      'spans in layer_spans_km above the ground. Fluxes in W/m2, from the ground to the air positive. '
                      'mixed_layer_km is CM1\'s boundary-layer depth.'),
        gcm=gcm, ring=ring)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('case')
    parser.add_argument('--from-day', type=float, default=29.5)
    args = parser.parse_args(argv)
    out = compare(args.case, args.from_day)
    out['producer'] = dict(domain='climate', files={f'crm/{p.name}': hashlib.sha256(p.read_bytes()).hexdigest()[:16]
                                                    for p in (Path(__file__), HERE / 'ring_analysis.py')})
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / f'profiles_{args.case}.json').write_text(json.dumps(out, indent=1) + '\n')
    for s in ('land', 'water'):
        print(f'\n{s} (GCM land {out["gcm"]["land_height_m"]:.0f} m up)' if s == 'land' else f'\n{s}')
        print('  hour angle           ' + ''.join(f'{h:7.0f}' for h in out['hour_angle_deg']))
        for label, side in (('GCM', out['gcm'][s]), ('CM1', out['ring'][s])):
            for key in ('air_c', 'ground_c', 'air_minus_layer2_k', 'layer1_minus_layer2_k', 'sensible_w_m2',
                        'latent_w_m2'):
                print(f'  {label} {key:18s}' + ''.join(f'{v:7.1f}' for v in side[key]))
    return 0


if __name__ == '__main__':
    sys.exit(main())
