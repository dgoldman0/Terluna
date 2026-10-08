"""The Open Moon's air with its lightning: the design column's photochemistry with the nitrogen oxides its storms make.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.joint_synthesis.air_chemistry [--processes N]
    # -> results/air_chemistry.json, with profiles in results/air_chemistry/

The middle-atmosphere model (atmosphere/middle_atmosphere) solves the design column (1.2 atm, the titania stack,
288 K) with no source of nitrogen oxides, so its odd nitrogen stays near zero. Its chemistry takes a column source of
lightning NO spread through the troposphere by air mass (chemistry.solve's lightning_no_per_cm2_s). This study
supplies that source from the electrified box's nitrogen yields (research/studies/atmospheric_electricity) carried
over the whole Moon by the flash climatology (atmosphere/electricity/results/occurrence_A28_dim5_moon.json), and
reads what it does to ozone near the ground and aloft, to the surface ultraviolet, to the odd nitrogen, and to NO and
O at 0.3 Pa, where the loss response's infrared cooling (atmosphere/loss_response/infrared.py) reads them. The source
enters through the solver's existing argument, so the domain's files and their products stay as they are.
"""
from __future__ import annotations
import argparse
import csv
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

from atmosphere.middle_atmosphere import chemistry as ch, equilibrium as eq, fetch_inputs, run as mrun
from atmosphere.radiative_convective import fetch_inputs as rc_inputs
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.joint-synthesis-air-chemistry/1'
BASE = 'moon_1.2atm_titania_stack'
MIDDLE = ROOT / 'atmosphere' / 'middle_atmosphere' / 'results'
OCCURRENCE = ROOT / 'atmosphere' / 'electricity' / 'results' / 'occurrence_A28_dim5_moon.json'
NITROGEN = ROOT / 'research' / 'studies' / 'atmospheric_electricity' / 'results' / 'nitrogen_oxides.json'
OUT = HERE / 'results' / 'air_chemistry.json'
PROFILES = HERE / 'results' / 'air_chemistry'
KG_N_KM2_YR_TO_MOLECULES_CM2_S = 1e3 / 14.0067 * 6.02214076e23 / 1e10 / (365.25 * 86400.0)
SPECIES = ('O3', 'NO', 'NO2', 'HNO3', 'OH', 'HO2', 'H2', 'O')
EVIDENCE = ('One-dimensional, global- and diurnal-mean photochemistry of the design column (atmosphere/middle_'
            'atmosphere) with a column source of lightning NO spread through the troposphere by air mass. The source '
            'is one electrified box\'s nitrogen yield carried over the Moon by its convective rain; the chemistry '
            'carries O, H and N families with H2 at 0.53 ppm and no methane, carbon monoxide or organic compounds, '
            'which on Earth add to ozone made from nitrogen oxides.')
READING_RULE = ('Mixing ratios as mole fractions (ppb where labelled). source_molecules_cm2_s is the column NO '
                'source; cases compare with the base column at the same surface temperature. near_ground is the '
                'lowest layer; at_35km is the storms\' charging level; at_0_3pa the base the loss response reads. '
                'uv_index values are the model\'s clear-sky global-mean and subsolar indices.')


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def sources() -> dict:
    """Moon-mean column NO sources (molecules/cm2/s) from the box's yields and the flash climatology."""
    occ = json.loads(OCCURRENCE.read_text())
    nit = json.loads(NITROGEN.read_text())['estimate']
    ratios = [c['ratio_to_box_moon_mean'] for by_sea in occ['cases'].values() for c in by_sea.values()]
    central = occ['cases']['linear']['1']['ratio_to_box_moon_mean']
    out = {}
    for route in ('per_flash', 'per_joule', 'per_metre'):
        low, high = nit[route]['kg_n_per_km2_yr_range']
        out[route] = dict(box_kg_n_km2_yr=nit[route]['kg_n_per_km2_yr'], box_range=[low, high],
                          moon_kg_n_km2_yr=nit[route]['kg_n_per_km2_yr'] * central,
                          moon_range_kg_n_km2_yr=[low * min(ratios), high * max(ratios)])
    s = lambda kg: kg * KG_N_KM2_YR_TO_MOLECULES_CM2_S
    return dict(flash_ratio_central=central, flash_ratio_range=[min(ratios), max(ratios)], routes=out,
                cases={'per_joule': s(out['per_joule']['moon_kg_n_km2_yr']),
                       'per_joule_high': s(out['per_joule']['moon_range_kg_n_km2_yr'][1])})


def with_lightning(rate):
    """chemistry.solve with a lightning source, for equilibrium's calls."""
    base = ch.solve
    def solve(*args, **kwargs):
        kwargs.setdefault('lightning_no_per_cm2_s', rate)
        return base(*args, **kwargs)
    return base, solve


def read_profile(path) -> dict:
    rows = list(csv.DictReader(open(path, newline='')))
    return {k: np.array([float(r[k]) if r[k] not in ('', None) else np.nan for r in rows]) for k in rows[0]}


def levels(profile) -> dict:
    """Mixing ratios near the ground, at the storms' charging level and at the loss response's base."""
    p, z = profile['p_pa'], profile['z_km']
    order = np.argsort(p)
    at_p = lambda name, pa: float(np.exp(np.interp(np.log(pa), np.log(p[order]), np.log(np.maximum(profile[name][order], 1e-30)))))
    order_z = np.argsort(z)
    at_z = lambda name, km: float(np.interp(km, z[order_z], profile[name][order_z]))
    ground = int(np.argmax(p))
    return dict(near_ground={s: float(profile[s][ground]) for s in SPECIES},
                at_35km={s: at_z(s, 35.0) for s in SPECIES},
                at_0_3pa={s: at_p(s, 0.3) for s in ('O', 'NO')})


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--processes', type=int, default=0, help='worker processes for line-by-line sunlight')
    parser.add_argument('--only', nargs='*', help='source cases to run (per_joule, per_joule_high)')
    args = parser.parse_args(argv)
    if args.processes:
        eq.SW_CONFIG = replace(eq.optics.SHORTWAVE, processes=args.processes)
    for restore in (fetch_inputs.restore, rc_inputs.restore):
        if restore()['status'] != 'PASS':
            print('BLOCKED: middle-atmosphere inputs missing', file=sys.stderr)
            return 1
    src = sources()
    base_case = mrun.cases()[BASE]
    base_profile = read_profile(MIDDLE / 'profiles' / f'{BASE}.csv')
    base_summary = json.loads((MIDDLE / 'cases' / f'{BASE}.json').read_text())
    PROFILES.mkdir(parents=True, exist_ok=True)
    previous = json.loads(OUT.read_text()) if OUT.is_file() else {}
    results = dict(previous.get('cases', {}))
    solver = eq.Solver()
    for label, rate in src['cases'].items():
        name = f'{BASE}_lightning_{label}'
        if args.only and label not in args.only:
            continue
        if name in results and (PROFILES / f'{name}.csv').is_file():
            continue                                   # restartable: a finished case is kept
        original, patched = with_lightning(rate)
        eq.ch.solve = patched
        try:
            start = time.time()
            res = solver.solve(replace(base_case, name=name))
            summary = mrun.summarise(res, time.time() - start)
            mrun._write_csv(PROFILES / f'{name}.csv', mrun.profile_rows(res, solver))
        finally:
            eq.ch.solve = original
        profile = read_profile(PROFILES / f'{name}.csv')
        results[name] = dict(source_molecules_cm2_s=rate, summary={k: summary[k] for k in mrun.FIELDS if k in summary},
                             **levels(profile))
        OUT.write_text(json.dumps(dict(schema=SCHEMA, cases=results), indent=1) + '\n')     # checkpoint
        print(f"{name}: source {rate:.3g} /cm2/s  O3 {summary['ozone_du']:.3g} DU  near-ground O3 "
              f"{results[name]['near_ground']['O3'] * 1e9:.3g} ppb  NO+NO2 {sum(results[name]['near_ground'][s] for s in ('NO', 'NO2')) * 1e9:.3g} ppb  "
              f"UV index {summary['uv_index_subsolar']:.3g}  ({summary['seconds']:.0f} s)", flush=True)
    files = {str(Path(f).relative_to(ROOT)): digest(f) for f in
             (__file__, eq.__file__, ch.__file__, mrun.__file__)}
    product = dict(schema=SCHEMA,
                   producer=dict(study='joint_synthesis', files=files, constants=constants_used(files),
                                 inputs={str(p.relative_to(ROOT)): digest(p) for p in
                                         (OCCURRENCE, NITROGEN, MIDDLE / 'profiles' / f'{BASE}.csv', MIDDLE / 'cases' / f'{BASE}.json')}),
                   evidence=EVIDENCE, reading_rule=READING_RULE, base_case=BASE, sources=src,
                   base=dict(summary={k: base_summary[k] for k in mrun.FIELDS if k in base_summary}, **levels(base_profile)),
                   cases=results,
                   not_run={label: 'not run in this product' for label in src['cases'] if f'{BASE}_lightning_{label}' not in results})
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    print(OUT)
    return 0


if __name__ == '__main__':
    sys.exit(main())
