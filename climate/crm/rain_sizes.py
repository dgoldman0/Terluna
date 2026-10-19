"""Where a CM1 run's rain sits against Morrison's limits on its size: the rain's mass by mean diameter, the shares beyond
the breakup cutoff and the equilibrium and held at the cap, and how thin the rain is that the cap holds.

    python -m climate.crm.rain_sizes box_0e_small --from-day 29.5        # about a minute for the small box
    # -> climate/results/crm/rain_sizes_<case>.json

Morrison's rain has an exponential spectrum whose mean diameter is 1/lambda = (q_r / (pi rho_w N_r))^(1/3), from the
rain's mass and number mixing ratios that every 3-D output holds (CM1's rho_w, 997 kg/m3). Each grid cell with rain
above 1e-8 kg/kg, where the scheme's self-collection and breakup act, counts with its rain's mass. The case's limits are
those its namelist and build set (cm1_run.rain_limits): CM1's own, or those of the rain-limits patch.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

from climate.crm import cm1_run, ring_analysis as ra

RUNS, RESULTS = ra.RUNS, ra.RESULTS
SCHEMA = 'terluna.climate.crm-rain-sizes/1'
RHO_W = 997.0                         # Morrison's RHOW
SMALLEST_RAIN = 1.0e-8                # kg/kg: below it the scheme neither self-collects nor breaks rain
EDGES_MM = (0.0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0)
CAP_TOLERANCE = 0.005                 # a cell within 0.5% of the cap is held by it
EVIDENCE = ('CM1\'s own state: rain mass and number in every 3-D output of the span, read as Morrison\'s exponential '
            'spectrum. It says where the scheme puts the rain\'s sizes under its limits; it measures no drop.')
READING_RULE = ('Diameters in mm. mass_share gives the share of the rain\'s mass (grid cells with rain above 1e-8 kg/kg, '
                'weighted by rain mass per square metre) whose mean diameter 1/lambda lies in each bin of edges_mm. '
                'limits are the case\'s own (cm1_run.rain_limits). at_cap is the share within 0.5% of the cap on the '
                'mean diameter; held_cells describes those cells and all rain cells by rain content (g/m3) and height '
                '(km above the model\'s ground).')


def namelist_settings(case: Path) -> dict:
    """The rain-limits entries of a case's namelist, CM1's own defaults where it has none."""
    path = case / 'namelist.input' if (case / 'namelist.input').exists() else case / 'namelist.template'
    text = path.read_text()
    settings = dict(brkup_db=0.0, brkup_gscale=0, brkup_rhor=1.0, brkup_rhoexp=0.0)
    for key in settings:
        try:
            settings[key] = type(settings[key])(float(cm1_run.namelist_value(text, 'param3', key)))
        except (AttributeError, ValueError):
            pass
    return settings


def census(name: str, from_day: float, stride: int = 1) -> dict:
    case = RUNS / name
    record = json.loads((case / 'case.json').read_text())
    gravity = record['build']['gravity_m_s2']
    limits = cm1_run.rain_limits(gravity, namelist_settings(case))
    scale = record.get('time_scale', 1.0)
    output_s = record['configuration']['output_s'] * scale
    numbers = sorted(int(p.name[8:14]) for p in case.glob('cm1out_t*_s.dat'))
    first = int(round(from_day * 86400.0 * scale / output_s)) + 1
    used = [n for n in numbers if n >= first][::stride]
    if not used:
        raise RuntimeError(f'{name} has no output after day {from_day}')
    zw = np.loadtxt(case / 'input_grid_z')
    zh, dz = 0.5 * (zw[1:] + zw[:-1]), np.diff(zw)
    edges = np.asarray(EDGES_MM) * 1e-3
    mass = np.zeros(len(EDGES_MM))                        # the last bin runs from the last edge up
    cap = limits['mean_diameter_cap_m']
    weights, sizes, contents, heights = [], [], [], []
    held_contents, held_heights = [], []
    path_total = 0.0
    for n in used:
        s = ra.read_snapshot(case, n)
        qr, nr = s['qr'], s['ncr']
        temperature = s['th'] * (s['prs'] / 1e5)**(287.04 / 1005.7)
        rho = s['prs'] / (287.04 * temperature * (1.0 + 0.608 * s['qv']))
        rain = (qr > SMALLEST_RAIN) & (nr > 0)
        mean = (qr[rain] / (np.pi * RHO_W * nr[rain]))**(1.0 / 3.0)
        content = (qr * rho)[rain]                                        # kg/m3
        w = (qr * rho * dz[:, None])[rain]                                 # kg/m2
        path_total += float((qr * rho * dz[:, None]).sum() / qr.shape[1])
        np.add.at(mass, np.clip(np.searchsorted(edges, mean, side='right') - 1, 0, len(EDGES_MM) - 1), w)
        z = np.broadcast_to(zh[:, None], qr.shape)[rain]
        held = mean >= cap * (1.0 - CAP_TOLERANCE)
        weights.append(w)
        sizes.append(mean)
        contents.append(content)
        heights.append(z)
        held_contents.append(content[held])
        held_heights.append(z[held])
    w, d = np.concatenate(weights), np.concatenate(sizes)
    total = w.sum()
    order = np.argsort(d)
    cumulative = np.cumsum(w[order]) / total
    median = float(d[order][np.searchsorted(cumulative, 0.5)])
    held_c, held_z = np.concatenate(held_contents), np.concatenate(held_heights)
    all_c, all_z = np.concatenate(contents), np.concatenate(heights)

    def spread(values, unit):
        if values.size == 0:
            return None
        return dict(median=float(np.median(values) * unit), p10=float(np.percentile(values, 10) * unit),
                    p90=float(np.percentile(values, 90) * unit))

    return dict(
        schema=SCHEMA, case=name, purpose=record.get('purpose') or record['configuration'].get('purpose'),
        build=record['build']['label'], executable_sha256=record['build']['executable_sha256'],
        evidence=EVIDENCE, reading_rule=READING_RULE, span_days=[used[0] * output_s / 86400.0 / scale,
                                                                 used[-1] * output_s / 86400.0 / scale],
        snapshots=len(used), stride=stride, limits={k: v * (1e3 if k.endswith('_m') else 1.0) for k, v in limits.items()},
        edges_mm=list(EDGES_MM), mass_share=[float(m / total) for m in mass],
        mass_weighted_mean_diameter_mm=float(np.sum(w * d) / total * 1e3), mass_median_mean_diameter_mm=median * 1e3,
        beyond_cutoff=float(w[d > limits['breakup_cutoff_m']].sum() / total),
        beyond_equilibrium=float(w[d > limits['equilibrium_m']].sum() / total),
        at_cap=float(w[d >= cap * (1.0 - CAP_TOLERANCE)].sum() / total),
        rain_water_path_kg_m2=path_total / len(used),
        held_cells=dict(count=int(held_c.size), share_of_rain_cells=float(held_c.size / all_c.size),
                        rain_g_m3=spread(held_c, 1e3), height_km=spread(held_z, 1e-3)),
        rain_cells=dict(count=int(all_c.size), rain_g_m3=spread(all_c, 1e3), height_km=spread(all_z, 1e-3)))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('case')
    parser.add_argument('--from-day', type=float, default=29.5)
    parser.add_argument('--stride', type=int, default=1, help='take every stride-th output')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args(argv)
    result = census(args.case, args.from_day, args.stride)
    result['producer'] = dict(domain='climate', files={f'crm/{Path(f).name}': hashlib.sha256(Path(f).read_bytes()).hexdigest()[:16]
                                                       for f in (__file__, ra.__file__, cm1_run.__file__)})
    path = args.output or RESULTS / f'rain_sizes_{args.case}.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=1) + '\n')
    print(f"{path}: {result['snapshots']} snapshots; {100 * result['beyond_cutoff']:.0f}% of the rain's mass beyond the "
          f"cutoff, {100 * result['at_cap']:.2f}% at the cap; mass-weighted 1/lambda "
          f"{result['mass_weighted_mean_diameter_mm']:.2f} mm")
    return 0


if __name__ == '__main__':
    sys.exit(main())
