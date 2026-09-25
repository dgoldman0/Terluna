"""Balance surface temperatures from equilibrium columns solved at several surface temperatures.

    python -m atmosphere.middle_atmosphere.balance      # writes results/balance.json

For each case solved at 288 K and at other surface temperatures (the *_ts<T>
cases of run.py), the top-of-atmosphere surplus of the clear column is
interpolated linearly between solved states to where surplus + C = 0, for an
assumed cloud effect C. Each state has its own ozone, stratosphere and
tropopause, so the balance includes how the upper atmosphere responds; C is an
assumption, not a result.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
CLOUD_EFFECTS = (0.0, -10.0, -20.0, -30.0, -40.0)
BASES = ('earth_control', 'moon_1.2atm_titania_stack', 'moon_1.0atm_titania_stack', 'moon_1.2atm_edge_200nm',
         'moon_1.0atm_edge_200nm', 'moon_1.2atm_edge_220nm', 'moon_1.2atm_edge_230nm', 'moon_1.2atm_edge_240nm',
         'moon_1.2atm_edge_310nm', 'moon_1.2atm_none', 'moon_1.2atm_edge_200nm_no_n2o', 'moon_1.0atm_edge_200nm_no_n2o',
         'moon_1.2atm_edge_200nm_slow_mixing', 'moon_1.0atm_edge_200nm_slow_mixing',
         'moon_1.2atm_edge_200nm_wet_stratosphere')


def series(results, base):
    out = []
    for path in sorted((results / 'cases').glob(f'{base}*.json')):
        d = json.loads(path.read_text())
        if d['case'] == base or d['case'].startswith(f'{base}_ts'):
            out.append((d['surface_temperature_k'], d['net_toa'], d['ozone_du'], d['tropopause_k']))
    return sorted(out)


def crossing(ts, net):
    """Surface temperatures where net falls through zero with warming (stable balances)."""
    found = []
    for i in range(len(ts) - 1):
        if net[i] >= 0 > net[i + 1]:
            found.append(float(ts[i] + (ts[i + 1] - ts[i]) * net[i] / (net[i] - net[i + 1])))
    return found


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--results', type=Path, default=HERE / 'results')
    args = parser.parse_args(argv)
    record = dict(schema='terluna.atmosphere.balance-temperatures/1',
                  evidence=' '.join(part.replace('\n', ' ') for part in __doc__.split('\n\n')[1:]),
                  reading_rule=('balance_k[C] lists the surface temperatures, within the solved range, at which the '
                                'clear column plus cloud effect C (W/m^2) balances; an empty list means the balance '
                                'lies outside the solved range.'),
                  cases={})
    for base in BASES:
        pts = series(args.results, base)
        if len(pts) < 2:
            continue
        ts, net, du, trop = (np.array(v) for v in zip(*pts))
        record['cases'][base] = dict(
            surface_temperature_k=ts.tolist(), net_toa_w_m2=net.round(3).tolist(), ozone_du=du.round(1).tolist(),
            tropopause_k=trop.round(1).tolist(),
            net_slope_w_m2_k=float(np.polyfit(ts, net, 1)[0]),
            balance_k={f'{c:g}': crossing(ts, net + c) for c in CLOUD_EFFECTS})
        print(base, {c: [round(x, 1) for x in v] for c, v in record['cases'][base]['balance_k'].items()},
              'slope', round(record['cases'][base]['net_slope_w_m2_k'], 2), flush=True)
    (args.results / 'balance.json').write_text(json.dumps(record, indent=1) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
