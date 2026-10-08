"""A grounded tower in the electrified storms: the potential the storms' charge sets up where a tower's top would stand.

    python -m atmosphere.electricity.tower
    # -> atmosphere/electricity/results/tower_exposure.json

A grounded conductor holds its top at the ground's potential, so it meets the air at its top height with the whole
potential difference the storms' charge sets up there. When that difference passes a leader's inception potential,
the tip launches an upward leader and the tower starts its own lightning. The module solves the same difference
equations as the electrified CM1's field solver (climate/crm/fortran/terluna_lightning.F: second-order differences,
periodic sides, a stretched vertical grid, phi = 0 at the flat ground and at the 150 km top) for each saved charge
field of the canonical run box_0e_elec_corrected and its two ten-minute windows. Every column stands for a place the
tower could stand, so the statistics run over all columns and outputs. The potential is the storm's own, before the
tower changes the field; the tower's induced charge and the leaders themselves are left out.

Thresholds: the leader-crossing potentials adopted for the lunar storms (climate/crm README, "The leader's
crossing": a tip leading the air by 0.225 MV for a positive leader, 0.4 MV for a negative one) and Rizk's critical
potential for stable upward leaders from tall free-standing structures at sea-level density, 1.556 MV / (1 + 3.89/h).
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

from shared.constants import SPEED_OF_LIGHT, VACUUM_PERMEABILITY

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCHEMA = 'terluna.atmosphere.tower-exposure/1'
RUNS = ROOT / 'climate' / 'crm' / 'runs'
CASES = ('box_0e_elec_corrected', 'box_0e_elec_corrected_storms_first', 'box_0e_elec_corrected_storms_second')
OUT = HERE / 'results' / 'tower_exposure.json'
EPS0 = 1.0 / (VACUUM_PERMEABILITY * SPEED_OF_LIGHT ** 2)
TOP_HEIGHTS_KM = (24.0, 35.6)          # the port's crown above flat ground, and above sea level with the summit beneath
THRESHOLDS_MV = dict(positive_leader=0.225, negative_leader=0.4, rizk_sea_level=1.556)
SECONDS_PER_YEAR = 365.25 * 86400.0
EVIDENCE = ('The electrified CM1 box\'s saved charge (two lunar days every 3 hours, and two 1.5-day windows every 10 '
            'minutes) on its flat, periodic 6-km grid, with the potential solved again from it. A tower\'s top meets '
            'the storm\'s potential at its height; the tower\'s own charge, the summit\'s terrain beneath the storms '
            'and the leaders\' growth are outside this screen.')
READING_RULE = ('Potentials in MV relative to the ground. For each top height, share_of_time_above[threshold] is the '
                'share of column-outputs where |phi| at that height passes the threshold; episodes_per_year counts, in '
                'the ten-minute windows, how often a column\'s |phi| rises through it, per column per year scaled by the '
                'storm days the windows cover; percentiles describe |phi| over all columns and outputs.')


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def potential(q, zh, zf, dx, dy):
    """phi (V) from charge density q (C/m3, shape nk x nj x ni), the electrified CM1's difference equations."""
    nk, nj, ni = q.shape
    dz = np.diff(zf)
    below = np.concatenate([[2.0 * (zh[0] - zf[0])], np.diff(zh)])
    above = np.concatenate([np.diff(zh), [2.0 * (zf[-1] - zh[-1])]])
    lo, up = 1.0 / (dz * below), 1.0 / (dz * above)
    kx = 2.0 * np.pi * np.fft.fftfreq(ni)
    ky = 2.0 * np.pi * np.fft.fftfreq(nj)
    lam = ((2.0 * np.cos(kx)[None, :] - 2.0) / dx ** 2 + ((2.0 * np.cos(ky)[:, None] - 2.0) / dy ** 2 if nj > 1 else 0.0))
    rhs = np.fft.fft2(-q / EPS0, axes=(1, 2))
    diag = (-(lo + up))[:, None, None] + lam[None, :, :]
    diag[0] -= lo[0]                                                   # ghost below the ground holds -phi
    diag[-1] -= up[-1]                                                 # ghost above the top holds -phi
    sub, sup = lo[1:], up[:-1]                                         # coefficients of phi[k-1] and phi[k+1]
    c, d = np.zeros_like(rhs), np.zeros_like(rhs)                      # Thomas algorithm, every mode at once
    c[0], d[0] = sup[0] / diag[0], rhs[0] / diag[0]
    for k in range(1, nk):
        m = diag[k] - sub[k - 1] * c[k - 1]
        if k < nk - 1:
            c[k] = sup[k] / m
        d[k] = (rhs[k] - sub[k - 1] * d[k - 1]) / m
    phi = np.zeros_like(rhs)
    phi[-1] = d[-1]
    for k in range(nk - 2, -1, -1):
        phi[k] = d[k] - c[k] * phi[k + 1]
    return np.real(np.fft.ifft2(phi, axes=(1, 2)))


def rizk_mv(height_m, density_ratio):
    return THRESHOLDS_MV['rizk_sea_level'] * density_ratio / (1.0 + 3.89 / height_m)


def case_fields(case: Path):
    """Each saved output's potential at the top heights for every column, with the output times (s)."""
    from climate.crm import elec_analysis as ea, ring_analysis as ra
    record = json.loads((case / 'case.json').read_text())
    grid = record['grid']
    zf = np.loadtxt(case / 'input_grid_z')
    zh = ea.heights_km(case / 'cm1out_s.ctl') * 1e3
    ni, nj, dx = grid['nx'], grid['ny'], grid['dx_m']
    outputs = sorted(int(p.name[8:14]) for p in case.glob('cm1out_t*_s.dat'))
    out, rho_at = [], {}
    for n in outputs:
        s = ea.snapshot_charge(case, n, dx, dx, zf)
        q = s['net'].reshape(zh.size, nj, ni)
        phi = potential(q, zh, zf, dx, dx)
        row = {}
        for h in TOP_HEIGHTS_KM:
            k = np.searchsorted(zh, h * 1e3)
            w = (h * 1e3 - zh[k - 1]) / (zh[k] - zh[k - 1])
            row[h] = ((1 - w) * phi[k - 1] + w * phi[k]).ravel() / 1e6          # MV at every column
            rho_at[h] = float(np.mean((1 - w) * s['rho'][k - 1] + w * s['rho'][k]))
        out.append(row)
    return outputs, out, rho_at, float(record['configuration']['output_s'])


def summarise(series, density_ratio, covered_days):
    stats = {}
    for h in TOP_HEIGHTS_KM:
        a = np.abs(np.array([r[h] for r in series]))                   # outputs x columns
        thresholds = dict(THRESHOLDS_MV, rizk_at_height=rizk_mv(h * 1e3, density_ratio[h]))
        thresholds.pop('rizk_sea_level')
        above = {k: round(float((a > v).mean()), 5) for k, v in thresholds.items()}
        rises = {k: int(((a[1:] > v) & (a[:-1] <= v)).sum()) for k, v in thresholds.items()}
        per_column_year = {k: round(v / a.shape[1] / (covered_days / 365.25), 3) for k, v in rises.items()}
        stats[f'{h:g}'] = dict(thresholds_mv={k: round(v, 3) for k, v in thresholds.items()},
                               share_of_time_above=above, episodes_per_column_per_year=per_column_year,
                               percentiles_mv={f'p{p}': round(float(np.percentile(a, p)), 4) for p in (50, 90, 99, 99.9)},
                               max_mv=round(float(a.max()), 3), air_density_kg_m3=round(density_ratio[h] * 1.225, 3))
    return stats


def main(argv=None) -> int:
    result = {}
    for name in CASES:
        case = RUNS / name
        if not case.is_dir():
            result[name] = dict(missing=True)
            continue
        outputs, series, rho, step_s = case_fields(case)
        days = (len(outputs) - 1) * step_s / 86400.0
        result[name] = dict(outputs=len(outputs), output_s=step_s, days=round(days, 2))
        density_ratio = {h: rho[h] / 1.225 for h in rho}
        result[name]['by_top_height_km'] = summarise(series, density_ratio, days)
        print(name, len(outputs), 'outputs', flush=True)
    product = dict(schema=SCHEMA, producer=dict(domain='atmosphere', files={'electricity/tower.py': digest(__file__)},
                                                runs=list(CASES)),
                   evidence=EVIDENCE, reading_rule=READING_RULE, cases=result)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    print(OUT)
    return 0


if __name__ == '__main__':
    sys.exit(main())
