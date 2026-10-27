"""The heat and vapour a CM1 box carries through given heights, split by the way it carries them (27 October 2026; the
record is the section on the single column with CM1's own convection in climate/gcm/README.md).

    climate/gcm/.venv/bin/python -m climate.crm.box_lid     # ../results/crm/box_lid_near_side_sea_night.json, seconds

CM1's domain diagnostics write, at every w level, the domain-mean vertical flux of potential temperature and of vapour
in three parts (domaindiag.F, scalar_flux):
- resolved (thfr, qvfr): w times the scalar's departure from its level's mean, the scalar carried to the w level by
  sixth-order interpolation;
- the advection's implicit diffusion (thfd, qvfd): w times the difference between the fifth-order interpolation the
  advection uses and the sixth-order one, the part of the resolved advection that diffuses;
- subgrid (thfs, qvfs): -khv dS/dz, the boundary-layer scheme's own flux. The lunar boxes run CM1's simple boundary
  layer (ipbl 2) with no LES subgrid model, so khv is that scheme's diffusivity, l^2 sqrt(max(S^2 - N^2/Pr, small)).

Each output is one time step's mean over the box's columns. lid_split takes the w levels nearest given heights and
averages each part over the outputs of a span, in W/m2 (rho cp times the Exner function times the potential
temperature's flux; rho L times the vapour's), with the scheme's diffusivity there; grouped, when asked, by quarter of
the lunar day. The module's own product splits the held-night boxes' fluxes at the top of their boundary layer, 6.45
and 6.99 km, over the second of their two days.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from climate.crm import box_budget as bb
from climate.crm import cm1_output as co
from climate.crm import ring_analysis as ra

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'crm'
OUT = RESULTS / 'box_lid_near_side_sea_night.json'
SCHEMA = 'terluna.climate.crm-box-lid/1'
DAY, HOUR = 86400.0, 3600.0
LV = bb.XLV                                   # J/kg, CM1's latent heat of vaporization
LID_M = (6450.0, 6990.0)                      # the w levels at the top of the sea box's night boundary layer
NIGHTS = {'6 km': 'box_near_side_sea_night_6km', '3 km': 'box_near_side_sea_night_3km',
          '2 km': 'box_near_side_sea_night_2km'}
NIGHT_SPAN_H = (24.0, 48.0)                   # the held nights' second day, hours from each run's start
QUARTERS = (('sunrise-noon', -90.0, 0.0), ('noon-sunset', 0.0, 90.0), ('sunset-midnight', 90.0, 180.0),
            ('midnight-sunrise', -180.0, -90.0))
PARTS = ('resolved', 'diffusion', 'subgrid')


def quarter_name(hour_deg: float) -> str:
    h = (float(hour_deg) + 180.0) % 360.0 - 180.0
    return next(name for name, a, b in QUARTERS if a <= h < b)


def output_fluxes(s: dict, w: dict) -> dict:
    """One domain diagnostic's fluxes at its w levels (ground first), upward, W/m2, by part, with khv (m2/s): s and w
    are its scalar-level and w-level files as box_budget.read_diag reads them."""
    prs = np.array(w['prsinterp'], float)
    prs[0] = s['psfc']                                  # the lowest w level is the ground
    rho = np.asarray(w['rhointerp'], float)
    to_heat = rho * co.CP * (np.maximum(prs, 1.0) / co.P00) ** (co.RD / co.CP)
    to_vapour = rho * LV
    out = dict(mtime=float(s['mtime']), khv=np.asarray(w['khv'], float))
    for part, suffix in zip(PARTS, 'rds'):
        out[f'heat_{part}'] = to_heat * np.asarray(w[f'thf{suffix}'], float)
        out[f'vapour_{part}'] = to_vapour * np.asarray(w[f'qvf{suffix}'], float)
    return out


def diag_files(case: Path) -> list:
    """Every domain diagnostic of a case as (scalar-level file, w-level file, their layouts), over CM1's output
    configurations: a run on from a restart that writes a different number of output fields starts a new
    configuration, whose files carry a suffix (.02) and are numbered from 1 (cm1.F, nconfig)."""
    import re
    layouts, out = {}, []
    for path in sorted(case.glob('cm1out_diag_*_s*.dat')):
        m = re.match(r'cm1out_diag_(\d+)_s(\.\d+)?\.dat$', path.name)
        if not m:
            continue
        suffix = m.group(2) or ''
        if suffix not in layouts:
            layouts[suffix] = tuple(ra.grads_layout(case / f'cm1out_diag_{k}{suffix}.ctl')[1] for k in 'sw')
        out.append((path, case / f'cm1out_diag_{m.group(1)}_w{suffix}.dat', *layouts[suffix]))
    return out


def read_case(name: str, first_s: float = -np.inf, last_s: float = np.inf) -> dict:
    """A case's w levels (m) and the fluxes (output_fluxes) of every domain diagnostic after first_s up to last_s (model
    seconds), in time order, each with its local hour angle (degrees from noon) where the case runs a lunar day."""
    from climate.crm.cm1_run import solar_day_s
    from climate.crm.ring_crossings import local_hour_angle
    case = ra.RUNS / name
    record = json.loads((case / 'case.json').read_text())
    cfg = record['configuration']
    zw = np.asarray(ra.case_geometry(case)['zw'], float)
    held = record.get('held') is not None                # a held hour: the clock that times the Sun is stopped
    rows = []
    for s_path, w_path, lay_s, lay_w in diag_files(case):
        s = bb.read_diag(case, s_path, lay_s)
        if not first_s + 1.0 < s['mtime'] <= last_s + 1.0:
            continue
        row = output_fluxes(s, bb.read_diag(case, w_path, lay_w))
        row['hour'] = (record['held']['hour_angle_deg'] if held else
                       local_hour_angle(cfg['start_hour_angle_deg'], cfg['site']['lon_deg'], s['mtime'], solar_day_s()))
        rows.append(row)
    rows.sort(key=lambda r: r['mtime'])
    return dict(case=name, zw_m=zw, rows=rows, l_inf_m=bb.namelist_value((case / 'namelist.input').read_text(), 'l_inf'))


def split(rows, zw, heights_m) -> list:
    """The parts at the w levels nearest heights_m, averaged over rows (output_fluxes' dicts): for heat and vapour the
    resolved, diffusion and subgrid parts and their total (W/m2, upward), the resolved advection's sum (resolved plus
    diffusion), the subgrid part's share of the mean total and the share of the subgrid part with the diffusion, and the
    outputs' spread (standard deviation) of the total and of the subgrid part; and the scheme's mean diffusivity khv."""
    zw = np.asarray(zw, float)
    out = []
    for z in heights_m:
        k = int(np.argmin(np.abs(zw - z)))
        entry = dict(height_m=float(zw[k]), outputs=len(rows), khv_m2_s=float(np.mean([r['khv'][k] for r in rows])))
        for q in ('heat', 'vapour'):
            parts = {p: np.array([r[f'{q}_{p}'][k] for r in rows]) for p in PARTS}
            total = sum(parts.values())
            mean = {p: float(v.mean()) for p, v in parts.items()}
            share = (lambda x: x / float(total.mean())) if total.mean() != 0.0 else (lambda x: None)
            entry[q] = dict(mean, total=float(total.mean()), resolved_advection=mean['resolved'] + mean['diffusion'],
                            subgrid_share=share(mean['subgrid']),
                            subgrid_and_diffusion_share=share(mean['subgrid'] + mean['diffusion']),
                            total_spread=float(total.std()), subgrid_spread=float(parts['subgrid'].std()))
        out.append(entry)
    return out


def lid_split(name: str, heights_m=LID_M, first_s: float = -np.inf, last_s: float = np.inf,
              by_quarter: bool = False) -> dict:
    """A case's fluxes split by part at the w levels nearest heights_m over the outputs after first_s up to last_s
    (model seconds): over the whole span, or by quarter of the lunar day."""
    data = read_case(name, first_s, last_s)
    out = dict(case=name, l_inf_m=data['l_inf_m'], span_s=[None if not np.isfinite(first_s) else first_s,
                                                              None if not np.isfinite(last_s) else last_s])
    if by_quarter:
        out['quarters'] = {}
        for qname, _, _ in QUARTERS:
            rows = [r for r in data['rows'] if quarter_name(r['hour']) == qname]
            out['quarters'][qname] = split(rows, data['zw_m'], heights_m) if rows else None
    else:
        out['levels'] = split(data['rows'], data['zw_m'], heights_m)
    return out


def analyse() -> dict:
    grids = {label: lid_split(name, LID_M, NIGHT_SPAN_H[0] * HOUR, NIGHT_SPAN_H[1] * HOUR)
             for label, name in NIGHTS.items()}
    return dict(schema=SCHEMA, heights_asked_m=list(LID_M), span_hours=list(NIGHT_SPAN_H), grids=grids,
                evidence=('CM1\'s own domain diagnostics of the three held-night boxes (one realization each, 6, 3 and '
                          '2 km grids, CM1\'s simple boundary-layer scheme on at every spacing), each output one time '
                          'step\'s mean over the box\'s columns, half-hourly; no new runs.'),
                reading_rule=('grids: by spacing, the w levels nearest the heights asked, the parts averaged over the '
                              'outputs after hour 24 up to hour 48 of each run. heat and vapour in W/m2, upward '
                              'positive: resolved (thfr, qvfr), diffusion (the advection\'s implicit diffusion, thfd, '
                              'qvfd), subgrid (the boundary-layer scheme, thfs, qvfs), total their sum, '
                              'resolved_advection resolved plus diffusion; subgrid_share the subgrid part over the mean '
                              'total, subgrid_and_diffusion_share the subgrid part and the diffusion together over it; '
                              'total_spread and subgrid_spread the standard deviations over the outputs. '
                              'khv_m2_s: the scheme\'s mean diffusivity at the level. heat = rho cp (p/p0)^(R/cp) times '
                              'the potential temperature\'s flux; vapour = rho L times the vapour\'s, CM1\'s constants.'),
                producer=dict(domain='climate', files={f'crm/{p.name}': hashlib.sha256(p.read_bytes()).hexdigest()[:16]
                                                       for p in (Path(__file__), HERE / 'box_budget.py')}))


def compact(x, digits: int = 4):
    if isinstance(x, dict):
        return {str(k): compact(v, digits) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [compact(v, digits) for v in x]
    if isinstance(x, np.ndarray):
        return compact(x.tolist(), digits)
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (int, np.integer)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        return None if not np.isfinite(float(x)) else float(f'{float(x):.{digits}g}')
    return x


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--output', type=Path, default=OUT)
    args = parser.parse_args(argv)
    out = analyse()
    args.output.write_text(json.dumps(compact(out), indent=1) + '\n')
    for label, g in out['grids'].items():
        for lev in g['levels']:
            h, v = lev['heat'], lev['vapour']
            print(f"{label} {lev['height_m'] / 1000:.2f} km ({lev['outputs']} outputs): heat total {h['total']:.1f} W/m2 "
                  f"(resolved {h['resolved']:.1f}, diffusion {h['diffusion']:.1f}, subgrid {h['subgrid']:.1f}; "
                  f"shares {h['subgrid_share']:.2f}, {h['subgrid_and_diffusion_share']:.2f}); vapour total "
                  f"{v['total']:.1f} (resolved {v['resolved']:.1f}, diffusion {v['diffusion']:.1f}, subgrid "
                  f"{v['subgrid']:.1f}; shares {v['subgrid_share']:.2f}, {v['subgrid_and_diffusion_share']:.2f}); "
                  f"khv {lev['khv_m2_s']:.1f} m2/s")
    print(f'wrote {args.output}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
