#!/usr/bin/env python3
"""What the GCM's 3-day means miss of its eddies: the sideways heat and vapour transport at a ring's row from a run
that writes means every 6.5 hours, against the same outputs averaged in blocks of eleven, the 3-day means the
runner writes as standard.

    climate/gcm/.venv/bin/python -m climate.crm.gcm_eddies A28_corrected_970_eddies 0-1 --latitude 45 --block 11

Writes results/crm/gcm_eddies_<run>.json. Each level's heating is also given in W/m², through the level's own mass:
the surface pressure times its depth in sigma over the Moon's gravity, times the air's heat capacity.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from climate.crm import ring_heat_budget as hb

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'crm'
SCHEMA = 'terluna.climate.crm-gcm-eddies/1'
CP = 1004.0                     # J/kg/K, dry air
DAY = 86400.0


def sigma_depths(sigma):
    """Depths in sigma of the GCM's layers (the top first) whose full levels sit midway between half levels running
    from 0 at the top to 1 at the ground."""
    half = [1.0]
    for s in sigma[::-1]:
        half.append(2.0 * s - half[-1])
    half = np.array(half[::-1])
    return half[1:] - half[:-1], half


def compare(run, years, latitude, block, surface_pa, gravity):
    """The run's sideways transport from its outputs as written and from block means of them, by level."""
    fine = hb.gcm_sideways(run, years, latitude, None)
    coarse = hb.gcm_sideways(run, years, latitude, None, block=block)
    depth, half = sigma_depths(np.asarray(fine['sigma']))
    per_k_day = CP * surface_pa * depth / gravity / DAY                   # W/m2 for 1 K/day in each layer
    levels = []
    for k, s in enumerate(fine['sigma']):
        written, blocks = fine['eddies_k_day'][k], coarse['eddies_k_day'][k]
        levels.append(dict(
            sigma=s, sigma_span=[float(half[k]), float(half[k + 1])], w_m2_per_k_day=float(per_k_day[k]),
            eddies_k_day=dict(written=written, blocks=blocks, missed=written - blocks),
            eddies_w_m2=dict(written=written * per_k_day[k], blocks=blocks * per_k_day[k],
                             missed=(written - blocks) * per_k_day[k]),
            mean_meridional_k_day=fine['mean_meridional_k_day'][k],
            eddies_vapour_g_kg_day=dict(written=fine['vapour']['eddies_g_kg_day'][k],
                                        blocks=coarse['vapour']['eddies_g_kg_day'][k])))
    lowest = levels[::-1]
    return dict(schema=SCHEMA, run=run, years=list(years), latitude_deg=latitude, rows_deg=fine['rows_deg'],
                block=block, surface_pa=surface_pa, levels=levels,
                lowest_two_w_m2={k: sum(row['eddies_w_m2'][k] for row in lowest[:2]) for k in ('written', 'blocks', 'missed')},
                evidence=('The GCM\'s zonal-eddy heat and vapour flux convergence at the row and the rows either side, '
                          'from each output\'s deviations from its zonal mean, averaged over the run; the block means '
                          'hold the same data as coarser means would, so the difference is what those means miss. '
                          'Covariances faster than the outputs themselves stay out of both.'),
                reading_rule=('K/day of temperature at each GCM level (the top first), and W/m2 through each level\'s '
                              'mass; positive warms the air. The mean meridional advection is the same either way.'))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('run')
    parser.add_argument('years', help='first-last model years, for example 0-1')
    parser.add_argument('--latitude', type=float, default=45.0)
    parser.add_argument('--block', type=int, default=11)
    args = parser.parse_args(argv)
    from climate.crm.cm1_run import planet
    years = tuple(int(x) for x in args.years.split('-'))
    record = json.loads((HERE.parent / 'gcm' / 'runs' / args.run / 'progress.json').read_text())
    import ast
    cfg = record['configuration']
    cfg = ast.literal_eval(cfg) if isinstance(cfg, str) else cfg
    out = compare(args.run, years, args.latitude, args.block, float(cfg['pressure_pa']), planet()['gravity_m_s2'])
    out['producer'] = dict(domain='climate', files={f'crm/{p.name}': hashlib.sha256(p.read_bytes()).hexdigest()[:16]
                                                    for p in (Path(__file__), HERE / 'ring_heat_budget.py')})
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / f'gcm_eddies_{args.run}.json').write_text(json.dumps(out, indent=1) + '\n')
    print(f'rows {out["rows_deg"]}, block {args.block}')
    for row in out['levels'][::-1][:4]:
        e = row['eddies_k_day']
        print(f'sigma {row["sigma_span"][1]:.3f}-{row["sigma_span"][0]:.3f}: eddies {e["written"]:+.3f} as written, '
              f'{e["blocks"]:+.3f} in block means, {e["missed"]:+.3f} missed K/day ({row["w_m2_per_k_day"]:.1f} W/m2 '
              f'per K/day); mean meridional {row["mean_meridional_k_day"]:+.3f}')
    print('lowest two levels, W/m2:', {k: round(v, 2) for k, v in out['lowest_two_w_m2'].items()})
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
