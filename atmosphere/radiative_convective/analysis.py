"""Read the inverse-climate table: balance temperatures, runaway limits and sensitivities.

    python -m atmosphere.radiative_convective.analysis [--results DIR]

Balance is where ASR + C = OLR for an assumed, stated cloud radiative effect C
(W/m^2, negative when clouds cool). C is an input here, not a prediction: Earth's
observed global mean is about -20 W/m^2 (CERES), and a slowly rotating Moon could
differ in either direction.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
CLOUD_EFFECTS = (0.0, -10.0, -20.0, -30.0, -40.0, -60.0)


def load(results: Path):
    with open(results / 'inverse_climate.csv', newline='') as handle:
        rows = list(csv.DictReader(handle))
    for r in rows:
        for k, v in r.items():
            if k not in ('key', 'sweep', 'scenario', 'planet', 'humidity'):
                r[k] = float(v) if v not in ('', None) else math.nan
    return rows


def curve(rows, sweep, scenario, field):
    pts = sorted((r['ts_k'], r[field]) for r in rows if r['sweep'] == sweep and r['scenario'] == scenario
                 and not math.isnan(r[field]))
    return np.array([p[0] for p in pts]), np.array([p[1] for p in pts])


def balance(ts, net):
    """Surface temperatures where net (ASR + C - OLR) crosses zero from above (stable balances)."""
    out = []
    for i in range(len(ts) - 1):
        a, b = net[i], net[i + 1]
        if a >= 0 > b:
            out.append(float(ts[i] + (ts[i + 1] - ts[i]) * a / (a - b)))
    return out


def summarise(rows):
    names = sorted({r['scenario'] for r in rows})
    summary = {}
    for name in names:
        s = {}
        ts, olr = curve(rows, 'earthlike_humidity', name, 'olr')
        _, asr = curve(rows, 'earthlike_humidity', name, 'asr')
        if len(ts) >= 2 and len(asr) == len(ts):
            s['earthlike_humidity'] = dict(
                ts_k=ts.tolist(), olr=olr.round(2).tolist(), asr=asr.round(2).tolist(),
                balance_k={f'{c:g}': balance(ts, asr + c - olr) for c in CLOUD_EFFECTS},
                olr_slope_w_m2_k=float(np.polyfit(ts, olr, 1)[0]),
                net_slope_w_m2_k=float(np.polyfit(ts, olr - asr, 1)[0]))
            _, dlr = curve(rows, 'earthlike_humidity', name, 'surface_down_lw')
            _, sup = curve(rows, 'earthlike_humidity', name, 'surface_up_lw')
            s['earthlike_humidity']['surface_net_lw_loss'] = (sup - dlr).round(2).tolist()
        ts, olr = curve(rows, 'saturated', name, 'olr')
        _, asr = curve(rows, 'saturated', name, 'asr')
        if len(ts):
            i = int(np.argmax(olr))
            s['saturated'] = dict(ts_k=ts.tolist(), olr=olr.round(2).tolist(),
                                  asr=asr.round(2).tolist() if len(asr) == len(ts) else None,
                                  max_olr=float(olr[i]), max_olr_at_k=float(ts[i]),
                                  olr_at_hottest=float(olr[-1]))
            _, strat = curve(rows, 'saturated', name, 'stratospheric_h2o')
            s['saturated']['stratospheric_h2o'] = strat.tolist()
        co2 = sorted((r['co2_ppm'], r['olr'], r['asr']) for r in rows
                     if r['scenario'] == name and r['humidity'] == 'manabe_wetherald'
                     and r['stratosphere_k'] == 200.0 and r['ts_k'] == 290.0)
        if len(co2) >= 2:
            c, o, a = (np.array(v) for v in zip(*co2))
            s['co2_at_290k'] = dict(co2_ppm=c.tolist(), olr=o.round(3).tolist(), asr=a.round(3).tolist(),
                                    olr_change_per_doubling=float(np.polyfit(np.log2(c), o, 1)[0]))
        summary[name] = s
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--results', type=Path, default=HERE / 'results')
    args = parser.parse_args(argv)
    summary = summarise(load(args.results))
    (args.results / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
