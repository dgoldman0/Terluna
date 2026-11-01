#!/usr/bin/env python3
"""How much of its falling rain the GCM evaporates, against CM1's lunar rain at each box's place, and the constant
that matches each.

    climate/gcm/.venv/bin/python -m climate.gcm.rain_evaporation A28R_301k_cold:12-19 A28R_301k_cold_moist:0-5 \\
        A28R_301k_cold_g02:0-4 A28R_301k_cold_g04:0-4                     # a minute: results/gcm/rain_evaporation.json

PlaSim rains out what its large-scale and convective schemes condense at once, and then (rainmod's mkrain) each
unsaturated layer the precipitation falls through evaporates a share gamma of its saturation deficit each time step,
corrected for the latent cooling and at most the precipitation there; gamma, 0.01 at T21, was set for Earth and knows
nothing of the drops' size or fall. A run's evaporated share is measured two ways from its output: the rain formed,
the column integral of the large-scale and convective heating (codes 27 and 28, written before mkrain) over the latent
heat of vaporisation, against the rain reaching the ground; and the residual of the summed physics tendency (code 21)
less every listed process, mkrain's cooling, over the rain formed. The heating codes are snapshots at each output time
and the rain a mean over it, so only means over many outputs compare; shares conditioned on raining columns do not.
The runs at one gamma are pooled, their sums of rain formed and at the ground added.

Each CM1 box is set beside the GCM at its place: the small box at 0 degrees east (its rain test of 20 October) and its
twin on the design climate beside the GCM's equatorial land within 10 degrees of latitude and 20 of longitude of it,
the box over the Near-side Sea beside the GCM's sea at sea level within 11 and 20 degrees of 5 S, 45 W. A box's share
is the share of the water leaving its cloud that evaporates on the way down over its second lunar day
(box_analysis.py), and the constant that gives it at the box's place is interpolated in log gamma between the pooled
runs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE.parent / 'results' / 'gcm'
CRM = HERE.parent / 'results' / 'crm'
OUT = RESULTS / 'rain_evaporation.json'
SCHEMA = 'terluna.climate.gcm-rain-evaporation/2'
CP_J_KG_K = 1005.0
LV_J_KG = 2.501e6
PLASIM_GAMMA = 0.01
# The GCM at each box's place: centre, half-widths (degrees) and the surface taken.
REGIONS = {'cm1_box': dict(lat_deg=0.0, lon_deg=0.0, half_lat_deg=10.0, half_lon_deg=20.0, surface='land'),
           'near_side_sea': dict(lat_deg=-5.0, lon_deg=-45.0, half_lat_deg=11.0, half_lon_deg=20.0, surface='sea')}
# The CM1 boxes and the place each is set against; the small box's control keeps CM1's own rain limits.
BOXES = (('box_0e_small_lunar_rain_day29', 'cm1_box', 'the small box at 0 E, the lunar rain limits, forced by the GCM '
                                                       'before its corrections'),
         ('box_0e_design', 'cm1_box', 'the land at 0 E on the design climate, forced by the GCM alone, the lunar rain '
                                      'limits and rain law'),
         ('box_near_side_sea', 'near_side_sea', 'the Near-side Sea at 5 S, 45 W on the design climate, forced by the '
                                                'GCM alone, the lunar rain limits and rain law'))
CONTROL = 'box_0e_small_rain_control_day29'
LISTED = ('dt_nudge', 'dt_sfc', 'dt_vdif', 'dt_lw', 'dt_sw', 'dt_lsp', 'dt_conv', 'dt_dca')


def box_share(case: str, product: str | None = None):
    """A CM1 box's share of its falling water evaporated, from its box_analysis product, and the product's path; None
    while the box has no product."""
    path = CRM / f'box_{product or case}.json'
    if not path.exists():
        return None, None
    d = json.loads(path.read_text())['microphysics']
    return d[case]['share_of_falling_water_evaporated'], path


def region_mask(lat, lon, land, height, region):
    rows = np.abs(lat - region['lat_deg']) <= region['half_lat_deg']
    cols = np.abs((lon - region['lon_deg'] + 180.0) % 360.0 - 180.0) <= region['half_lon_deg']
    cells = rows[:, None] & cols[None, :]
    return cells & (land if region['surface'] == 'land' else ~land & (height < 1.0))


def measure(run, years):
    """A run's sums over the years, globally, over land and sea and at each region, of the rain formed and the rain at
    the ground (kg/m2/s times area weight), and globally of mkrain's residual cooling, with its shares."""
    import netCDF4
    cfg = json.loads((HERE / 'runs' / run / 'progress.json').read_text())['configuration']
    gravity = cfg['planet']['gravity']
    gamma = cfg.get('rain_evaporation', {}).get('gamma', PLASIM_GAMMA)
    folder = HERE / 'runs' / run
    sums, residual = {}, [0.0, 0.0]
    for year in years:
        with netCDF4.Dataset(folder / 'model' / f'MOST.{year:05d}.nc') as d:
            sigma = np.asarray(d['lev'][:], float)
            ps = np.asarray(d['ps'][:], float) * 100.0
            half = [0.0]
            for s in sigma:
                half.append(2.0 * s - half[-1])
            dp = np.diff(np.array(half))[None, :, None, None] * ps[:, None]
            heat = {k: np.asarray(d[k][:], float) for k in ('dt_diab', 'dt_lsp', 'dt_conv') + LISTED}
            formed = ((heat['dt_lsp'] + heat['dt_conv']) * CP_J_KG_K * dp / gravity).sum(axis=1) / LV_J_KG
            listed = sum(heat[k] for k in LISTED)
            cooling = ((heat['dt_diab'] - listed) * CP_J_KG_K * dp / gravity).sum(axis=1) / LV_J_KG
            ground = np.asarray(d['pr'][:], float) * 1000.0                  # m/s of water to kg/m2/s
            lat = np.asarray(d['lat'][:], float)
            lon = np.asarray(d['lon'][:], float)
            land = np.asarray(d['lsm'][:], float)[0] > 0.5
        height = np.loadtxt(folder / 'inputs' / 'moon_topography.sra', skiprows=1).ravel().reshape(land.shape) / gravity
        w = np.cos(np.radians(lat))[None, :, None] * np.ones(ps.shape)
        masks = dict(global_=np.ones(land.shape, bool), land=land, sea=~land,
                     **{k: region_mask(lat, lon, land, height, r) for k, r in REGIONS.items()})
        for key, mask in masks.items():
            s = sums.setdefault(key.rstrip('_'), [0.0, 0.0, 0.0])
            m = np.broadcast_to(mask, ps.shape)
            s[0] += float((formed * w * m).sum())
            s[1] += float((ground * w * m).sum())
            s[2] += float((w * m).sum())
        residual[0] += float((-cooling * w).sum())
        residual[1] += float(w.sum())
    return dict(run=run, years=list(years), gamma=gamma, sums=sums, residual=residual, shares=shares(sums, residual))


def shares(sums, residual=None):
    rows = {k: dict(formed_mm_day=86400.0 * f / a, ground_mm_day=86400.0 * g / a, evaporated_share=1.0 - g / f)
            for k, (f, g, a) in sums.items()}
    if residual:
        rows['global']['evaporated_share_by_residual'] = ((residual[0] / residual[1])
                                                          / (sums['global'][0] / sums['global'][2]))
    return rows


def pool(measured):
    """The runs at each gamma together: their sums added, with their years, by gamma in increasing order."""
    out = {}
    for m in measured:
        p = out.setdefault(m['gamma'], dict(gamma=m['gamma'], runs=[], sums={}, residual=[0.0, 0.0]))
        p['runs'].append(dict(run=m['run'], years=m['years']))
        for k, v in m['sums'].items():
            p['sums'][k] = [a + b for a, b in zip(p['sums'].get(k, [0.0, 0.0, 0.0]), v)]
        p['residual'] = [a + b for a, b in zip(p['residual'], m['residual'])]
    return [dict(gamma=g, runs=p['runs'], years=sum(len(r['years']) for r in p['runs']),
                 shares=shares(p['sums'], p['residual'])) for g, p in sorted(out.items())]


def fit(pooled, region, target):
    """gamma at which the region's evaporated share reaches the target, interpolated in log gamma between the pooled
    runs; None outside them or without a target."""
    if target is None:
        return None
    pts = sorted((p['gamma'], p['shares'][region]['evaporated_share']) for p in pooled)
    g, s = np.array([p[0] for p in pts]), np.array([p[1] for p in pts])
    if not (np.all(np.diff(s) > 0) and s.min() <= target <= s.max()):
        return None
    return float(np.exp(np.interp(target, s, np.log(g))))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('cases', nargs='+', metavar='RUN:FIRST-LAST')
    args = parser.parse_args(argv)
    measured = []
    for case in args.cases:
        run, span = case.rsplit(':', 1)
        first, _, last = span.partition('-')
        measured.append(measure(run, range(int(first), int(last or first) + 1)))
    pooled = pool(measured)
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    cm1, inputs = {}, {}
    for case, region, label in BOXES:
        share, path = box_share(case)
        cm1[case] = dict(region=region, about=label, share_of_falling_water_evaporated=share,
                         gamma=fit(pooled, region, share))
        if path:
            inputs[str(path.relative_to(ROOT))] = digest(path)
    control, path = box_share(CONTROL, CONTROL)
    if path:
        inputs[str(path.relative_to(ROOT))] = digest(path)
    out = dict(schema=SCHEMA, regions=REGIONS, cm1=cm1, cm1_limits_control=control,
               measured=[{k: m[k] for k in ('run', 'years', 'gamma', 'shares')} for m in measured], pooled=pooled,
               evidence=('Means over GCM output of the rain its schemes form and the rain at the ground, and of '
                         'mkrain\'s residual cooling, pooled over the runs at each gamma; each CM1 share from one box '
                         'over one lunar day with the lunar rain limits, one realization each.'),
               reading_rule=('pooled holds, for each gamma, the evaporated shares (1 less the rain at the ground over '
                             'the rain formed) globally, over land and sea and at each region; cm1 holds each box\'s '
                             'share, the region it is set against and the gamma interpolated in log between the pooled '
                             'runs at which that region evaporates the box\'s share, null outside them or before the '
                             'box has run.'),
               producer=dict(domain='climate', files={'gcm/rain_evaporation.py': digest(__file__)}, inputs=inputs))
    RESULTS.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    for p in pooled:
        s = p['shares']
        print(f"gamma {p['gamma']} ({p['years']} years): evaporated {s['global']['evaporated_share']:.3f} (residual "
              f"{s['global']['evaporated_share_by_residual']:.3f}), land {s['land']['evaporated_share']:.3f}, sea "
              f"{s['sea']['evaporated_share']:.3f}, " + ', '.join(f"{k} {s[k]['evaporated_share']:.3f}" for k in REGIONS))
    for case, row in cm1.items():
        print(f"{case}: share {row['share_of_falling_water_evaporated']}, gamma {row['gamma']}")
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
