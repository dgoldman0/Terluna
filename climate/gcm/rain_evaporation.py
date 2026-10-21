#!/usr/bin/env python3
"""How much of its falling rain the GCM evaporates, against CM1's lunar rain, and the constant that matches it.

    climate/gcm/.venv/bin/python -m climate.gcm.rain_evaporation A28R_301k_cold:18-19 A28R_301k_cold_g02:1 \\
        A28R_301k_cold_g04:1                                        # seconds: results/gcm/rain_evaporation.json

PlaSim rains out what its large-scale and convective schemes condense at once, and then (rainmod's mkrain) each
unsaturated layer the precipitation falls through evaporates a share gamma of its saturation deficit each time step,
corrected for the latent cooling and at most the precipitation there; gamma, 0.01 at T21, was set for Earth and knows
nothing of the drops' size or fall. A run's evaporated share is measured two ways from its output: the rain formed,
the column integral of the large-scale and convective heating (codes 27 and 28, written before mkrain) over the latent
heat of vaporisation, against the rain reaching the ground; and the residual of the summed physics tendency (code 21)
less every listed process, mkrain's cooling, over the rain formed. The heating codes are snapshots at each output time
and the rain a mean over it, so only means over many outputs compare; shares conditioned on raining columns do not.
CM1's small box at 0 degrees east (climate/crm, its rain test of 20 October) evaporates 48% of the water leaving its
cloud with the rain limits carried to lunar gravity and 86% with CM1's own; the GCM is measured on equatorial land
within 10 degrees of latitude and 20 of longitude of the box. The constant that gives the box's 48% there is
interpolated in log gamma between the runs given.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE.parent / 'results' / 'gcm'
OUT = RESULTS / 'rain_evaporation.json'
SCHEMA = 'terluna.climate.gcm-rain-evaporation/1'
CP_J_KG_K = 1005.0
LV_J_KG = 2.501e6
CM1 = {'lunar_limits': ROOT / 'climate' / 'results' / 'crm' / 'box_box_0e_small_lunar_rain_day29.json'}
BOX = dict(latitude_deg=10.0, longitude_deg=20.0)       # half-widths about 0 degrees east on the equator
LISTED = ('dt_nudge', 'dt_sfc', 'dt_vdif', 'dt_lw', 'dt_sw', 'dt_lsp', 'dt_conv', 'dt_dca')


def cm1_shares():
    """CM1's share of the falling water evaporated with the lunar limits and with CM1's own, as its products hold them."""
    d = json.loads(CM1['lunar_limits'].read_text())['microphysics']
    return dict(lunar_limits=d['box_0e_small_lunar_rain_day29']['share_of_falling_water_evaporated'],
                cm1_limits=d['box_0e_small_rain_control_day29']['share_of_falling_water_evaporated'])


def measure(run, years):
    """A run's evaporated share of its falling precipitation over the years: globally, over land and sea and at the
    CM1 box's place, by the rain formed against the rain at the ground, and globally by mkrain's residual cooling."""
    import netCDF4
    cfg = json.loads((HERE / 'runs' / run / 'progress.json').read_text())['configuration']
    gravity = cfg['planet']['gravity']
    gamma = cfg.get('rain_evaporation', {}).get('gamma', 0.01)
    sums = {}
    residual = [0.0, 0.0]
    for year in years:
        with netCDF4.Dataset(HERE / 'runs' / run / 'model' / f'MOST.{year:05d}.nc') as d:
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
            land = np.asarray(d['lsm'][:], float) > 0.5
            land = np.broadcast_to(land, ps.shape) if land.ndim == 2 else land
        w = np.cos(np.radians(lat))[None, :, None] * np.ones(ps.shape)
        lon_e = np.where(lon > 180.0, lon - 360.0, lon)
        box = ((np.abs(lat)[None, :, None] <= BOX['latitude_deg']) &
               (np.abs(lon_e)[None, None, :] <= BOX['longitude_deg'])) & land
        for key, mask in (('global', np.ones(ps.shape, bool)), ('land', land), ('sea', ~land), ('cm1_box', box)):
            s = sums.setdefault(key, [0.0, 0.0, 0.0])
            s[0] += float((formed * w * mask).sum())
            s[1] += float((ground * w * mask).sum())
            s[2] += float((w * mask).sum())
        residual[0] += float((-cooling * w).sum())
        residual[1] += float(w.sum())
    rows = {k: dict(formed_mm_day=86400.0 * f / a, ground_mm_day=86400.0 * g / a, evaporated_share=1.0 - g / f)
            for k, (f, g, a) in sums.items()}
    rows['global']['evaporated_share_by_residual'] = (residual[0] / residual[1]) / (sums['global'][0] / sums['global'][2])
    return dict(run=run, years=list(years), gamma=gamma, shares=rows)


def fit(measured, target):
    """gamma at which the CM1 box's evaporated share reaches the target, interpolated in log gamma."""
    pts = sorted((m['gamma'], m['shares']['cm1_box']['evaporated_share']) for m in measured)
    g, s = np.array([p[0] for p in pts]), np.array([p[1] for p in pts])
    if not (s.min() <= target <= s.max()):
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
    cm1 = cm1_shares()
    gamma = fit(measured, cm1['lunar_limits'])
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    out = dict(schema=SCHEMA, cm1=cm1, measured=measured, gamma_for_cm1_lunar_limits=gamma,
               evidence=('Means over GCM output of the rain its schemes form and the rain at the ground, and of '
                         'mkrain\'s residual cooling; CM1\'s share from one small box over one lunar day with the '
                         'lunar rain limits. The fit matches the box\'s place, which sits at the global mean in the '
                         'design years.'),
               reading_rule=('measured holds each run\'s gamma and its evaporated shares (1 less the rain at the '
                             'ground over the rain formed) globally, over land and sea and at the CM1 box\'s place; '
                             'gamma_for_cm1_lunar_limits is the gamma interpolated in log between the runs at which '
                             'the box\'s place evaporates CM1\'s share with the lunar limits, or null outside them.'),
               producer=dict(domain='climate', files={'gcm/rain_evaporation.py': digest(__file__)},
                             inputs={str(p.relative_to(ROOT)): dict(path=str(p.relative_to(ROOT)), sha256=digest(p))
                                     for p in CM1.values()}))
    RESULTS.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    for m in measured:
        s = m['shares']
        print(f"{m['run']} {m['years'][0]}-{m['years'][-1]}, gamma {m['gamma']}: evaporated {s['global']['evaporated_share']:.3f} "
              f"(residual {s['global']['evaporated_share_by_residual']:.3f}), land {s['land']['evaporated_share']:.3f}, "
              f"sea {s['sea']['evaporated_share']:.3f}, CM1 box {s['cm1_box']['evaporated_share']:.3f}; rain formed "
              f"{s['global']['formed_mm_day']:.2f}, at the ground {s['global']['ground_mm_day']:.2f} mm/day")
    print(f"CM1: {cm1}; gamma for the lunar limits' share at the box: {gamma}")
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
