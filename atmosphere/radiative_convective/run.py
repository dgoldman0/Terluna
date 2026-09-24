"""Scenario sweeps of the inverse climate calculation, written as a data product.

    python -m atmosphere.radiative_convective.run                # every sweep
    python -m atmosphere.radiative_convective.run --only earthlike_humidity saturated
    python -m atmosphere.radiative_convective.run --out DIR --list

Rows are appended as they finish, so an interrupted run resumes where it stopped.
The results file carries the schema, the producing code and input hashes, and
the evidence statement consumers must keep with the numbers.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import time

from atmosphere.radiative_convective import climate as cl, thermodynamics as th, fetch_inputs

HERE = Path(__file__).resolve().parent
SCHEMA = 'terluna.atmosphere.inverse-climate/1'
EVIDENCE = ('One-dimensional clear-sky calculations for prescribed columns: a saturated moist adiabat to a fixed '
            'stratospheric temperature, a prescribed relative-humidity profile, no clouds, no ozone, CH4 or N2O. '
            'Outgoing longwave radiation (OLR) and absorbed solar radiation (ASR) are global means for the stated '
            'surface temperature; they are not a solved climate, a cloud prediction or a habitability result. '
            'Where ASR equals OLR the column could be in balance; clouds shift that balance by an amount this '
            'calculation does not determine.')

EARTH_1 = dict(planet=th.EARTH, dry_pressure_pa=101325.0)
MOON_1 = dict(planet=th.MOON, dry_pressure_pa=101325.0)
MOON_12 = dict(planet=th.MOON, dry_pressure_pa=121590.0)
BASES = (('earth_1.0atm', EARTH_1), ('moon_1.0atm', MOON_1), ('moon_1.2atm', MOON_12))


def sweeps():
    """Named sweeps: (scenario, surface temperature, include solar)."""
    out = {}
    rows = []
    for name, base in BASES:
        for ts in (260.0, 270.0, 280.0, 290.0, 300.0, 310.0, 320.0):
            rows.append((cl.Scenario(name, co2_ppm=400.0, **base), ts, True))
    out['earthlike_humidity'] = rows
    rows = []
    for name, base in BASES:
        for ts in (250.0, 270.0, 290.0, 300.0, 310.0, 320.0, 330.0, 340.0, 350.0, 360.0, 370.0):
            rows.append((cl.Scenario(name, co2_ppm=400.0, humidity='saturated', **base), ts, True))
    out['saturated'] = rows
    rows = []
    for name, base in BASES:
        for co2 in (150.0, 280.0, 1000.0):
            rows.append((cl.Scenario(name, co2_ppm=co2, **base), 290.0, True))
    out['carbon_dioxide'] = rows
    rows = []
    for name, base in BASES[1:]:
        for strat in (150.0, 250.0):
            rows.append((cl.Scenario(name, co2_ppm=400.0, stratosphere_k=strat, **base), 290.0, False))
    out['stratosphere'] = rows
    # Sunlight filtered by each spectral shield: the columns of the first two sweeps,
    # solar part only (outgoing longwave does not depend on the sunlight).
    out['shields'] = [(sc, ts, True) for sweep in ('earthlike_humidity', 'saturated') for sc, ts, _ in out[sweep]]
    return out


FIELDS = ['key', 'sweep', 'scenario', 'planet', 'dry_pressure_pa', 'co2_ppm', 'humidity', 'stratosphere_k',
          'surface_albedo', 'shield', 'ts_k', 'surface_pressure_pa', 'olr', 'asr', 'albedo', 'net_toa', 'incident_sw',
          'surface_down_lw', 'surface_up_lw', 'surface_sw', 'atmosphere_sw', 'sw_longward_of_grid',
          'water_column_kg_m2', 'air_column_kg_m2', 'tropopause_pa', 'tropopause_height_km', 'stratospheric_h2o',
          'shield_product', 'seconds']


def key(sc: cl.Scenario, ts: float, shield: str = 'none') -> str:
    return (f'{sc.name}|{sc.humidity}|co2={sc.co2_ppm:g}|strat={sc.stratosphere_k:g}|albedo={sc.surface_albedo:g}'
            f'|ts={ts:g}|shield={shield}')


def _read(table):
    """Rows of the results table, migrating rows written before shields existed."""
    if not table.is_file():
        return []
    with open(table, newline='') as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        header = reader.fieldnames or []
    if 'shield' not in header:
        for r in rows:
            r['shield'] = 'none'
            r['key'] = r['key'] + '|shield=none'
        _write(table, rows)
    return rows


def _write(table, rows):
    with open(table, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, extrasaction='ignore')
        writer.writeheader()
        for r in rows:
            writer.writerow({f: r.get(f, '') for f in FIELDS})


RUN_MODULES = ('climate.py', 'thermodynamics.py', 'optics.py', 'longwave.py', 'shortwave.py', 'spectroscopy.py',
               'fetch_inputs.py', 'run.py')


def producer():
    digest = {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()[:16] for name in RUN_MODULES}
    inputs = hashlib.sha256((HERE / 'inputs.json').read_bytes()).hexdigest()[:16]
    return dict(domain='atmosphere', model='atmosphere/radiative_convective', files=digest, inputs_manifest=inputs)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out', type=Path, default=HERE / 'results')
    parser.add_argument('--only', nargs='*')
    parser.add_argument('--list', action='store_true')
    args = parser.parse_args(argv)
    plan = sweeps()
    if args.list:
        for name, rows in plan.items():
            print(name, len(rows))
        return 0
    if fetch_inputs.restore()['status'] != 'PASS':
        print('BLOCKED: spectroscopic inputs missing; run fetch_inputs --download', file=sys.stderr)
        return 1
    args.out.mkdir(parents=True, exist_ok=True)
    table = args.out / 'inverse_climate.csv'
    rows = _read(table)
    done = {r['key'] for r in rows}
    for sweep, plan_rows in plan.items():
        if args.only and sweep not in args.only:
            continue
        for sc, ts, solar in plan_rows:
            shields = cl.SHIELDS[1:] if sweep == 'shields' else ('none',)
            wanted = [s_ for s_ in shields if key(sc, ts, s_) not in done]
            if not wanted:
                continue
            start = time.time()
            reuse = next((r for r in rows if r['key'] == key(sc, ts, 'none') and r.get('olr')), None)
            results = cl.evaluate(sc, ts, solar=solar, longwave_too=reuse is None, shields=wanted)
            for result in results:
                if reuse is not None:
                    for f in ('olr', 'surface_down_lw', 'surface_up_lw'):
                        result[f] = float(reuse[f])
                    if 'asr' in result:
                        result['net_toa'] = result['asr'] - result['olr']
                result.update(key=key(sc, ts, result['shield']), sweep=sweep,
                              seconds=round((time.time() - start) / len(results), 1))
                rows.append({f: result.get(f, '') for f in FIELDS})
                done.add(result['key'])
                print(f"{sweep:18s} {sc.name:13s} {sc.humidity:17s} co2={sc.co2_ppm:6g} Ts={ts:5.1f} "
                      f"shield={result['shield']:13s} OLR={float(result['olr']):7.2f} "
                      f"ASR={float(result.get('asr', float('nan'))):7.2f}", flush=True)
            _write(table, rows)
    meta = dict(schema=SCHEMA, producer=producer(), evidence=EVIDENCE,
                units=dict(fluxes='W m-2, global mean', pressure='Pa', temperature='K', columns='kg m-2'),
                reading_rule=('Compare ASR with OLR at the same scenario, shield and surface temperature. Positive '
                              'net_toa means the clear column would warm. shield names a transmission in '
                              'protection/spectra/shield_transmission.json (none = unfiltered sunlight). Interpolate '
                              'between surface temperatures linearly; do not extrapolate beyond the tabulated range. '
                              'Keep the evidence statement with any use.'),
                table=table.name)
    (args.out / 'inverse_climate.json').write_text(json.dumps(meta, indent=2) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
