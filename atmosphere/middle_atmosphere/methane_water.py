"""The design air's methane in the upper air behind the titania stack: the water its oxidation can leave above the
tropopause, at its two bounds, and the exobase with it.

    OPENBLAS_NUM_THREADS=1 python -m atmosphere.middle_atmosphere.methane_water      # about 40 minutes on 4 cores

The design air holds 1.9 ppm of methane and 0.1 ppm of N2O (register, 13 October), which the middle atmosphere's
cases predate: they carry Earth's 330 ppb of N2O and no methane, whose chemistry the column leaves out. Methane passes
the cold trap that holds the water entering the upper air to 0.02-0.07 ppm behind the stack, and each molecule broken
above it leaves two of water. Two bounds bracket where that hydrogen goes: methane unbroken below the 0.3 Pa base,
the water above the tropopause the cold trap's alone; or methane broken entirely as soon as it passes the tropopause,
the cold trap's water plus twice the methane, 3.8 ppm, at every level above it (an upper bound on the water, since
oxidation builds with height). The three design cases behind the exobase table (cooling in local thermodynamic
equilibrium; CO2's 15 um band out of it with all the absorbed near-infrared turned to heat; with only its collisional
share) are solved at both bounds with the design air's N2O and the curved air's light from the 301 K twilight ledger
(curved_light.py's heating ratio), and the thermal column runs from each new 0.3 Pa base with the exobase table's
heats (curved_light.exobase). Methane's own infrared bands and chemistry stay out; the water's cooling and its
photochemistry are the column's. Each solve is kept on the drive by name, so a stopped run resumes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import replace
from pathlib import Path

from atmosphere.middle_atmosphere import curved_light as cv, equilibrium as eq, escape, run as runner, chemistry as ch

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / 'results' / 'methane_water.json'
SCHEMA = 'terluna.atmosphere.methane-water/1'
LEDGER = ROOT / 'climate' / 'results' / 'gcm' / 'twilight_301k.json'
METHANE_PPM = 1.9                      # the design air's methane (register, 13 October: the middle gas pair)
N2O_PPB = 100.0                        # and its N2O
BOUNDS = ('unbroken', 'broken')


def design_case(case, bound, trap_water):
    """A design case with the design air's N2O and, broken, the methane's water above the tropopause."""
    water = None if bound == 'unbroken' else trap_water + 2.0 * METHANE_PPM * 1e-6
    return replace(case, name=f'{case.name}_design_air_{bound}', n2o_ppb=N2O_PPB, stratospheric_h2o=water)


def analyse():
    ledger = json.loads(LEDGER.read_text())
    lp, ratio = cv.curved_ratio(ledger)
    factor = cv.curved_factor(ledger)
    cases = runner.cases()
    previous = json.loads(cv.OUT.read_text())['exobase'] if cv.OUT.exists() else {}
    out = dict(schema=SCHEMA, methane_ppm=METHANE_PPM, n2o_ppb=N2O_PPB, cases={})
    for name in cv.EXOBASE_CASES:
        stored = cv.stored(name)
        trap = float(stored['summary']['stratospheric_h2o'])
        row = dict(cold_trap_water=trap, curved_light_14_october=previous.get(name, {}).get('curved'))
        for bound in BOUNDS:
            case = design_case(cases[name], bound, trap)
            solution = cv.solved(f'{case.name}_{cv.tag(lp, ratio)}', case, factor)
            row[bound] = dict(water_above_tropopause=float(solution['summary']['stratospheric_h2o']),
                              tropopause_k=float(solution['summary']['tropopause_k']),
                              t_1pa_k=float(solution['summary']['t_1pa']),
                              solve_seconds=solution['seconds'], **cv.exobase(solution, solution['summary']))
            cv.log(f'{case.name}: water above the tropopause {row[bound]["water_above_tropopause"]:.2e}, '
                   f'base {row[bound]["base_temperature_k"]:.1f} K')
        out['cases'][name] = row
    out['ledger'] = dict(path=str(LEDGER.relative_to(ROOT)), sha256=hashlib.sha256(LEDGER.read_bytes()).hexdigest()[:16])
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--processes', type=int, default=cv.PROCESSES, help='processes for the near-infrared optics')
    cv.PROCESSES = parser.parse_args(argv).processes
    out = analyse()
    out['evidence'] = ('One-dimensional, global- and diurnal-mean, clear-sky radiative-convective and photochemical '
                       'equilibrium (equilibrium.py) for the three design cases behind the titania stack at a 288 K '
                       'ground, with the curved air\'s near-infrared heating, the design air\'s N2O and the water above '
                       'the tropopause at the two bounds of methane\'s oxidation; the thermal column above the 0.3 Pa '
                       'base as escape.py runs it. Methane\'s own bands and chemistry are left out.')
    out['reading_rule'] = ('cases holds each design case: the cold trap\'s water above the tropopause as stored, the '
                           '14 October check\'s curved-light exobase rows (curved_light.json) and, for each bound, the '
                           'water held above the tropopause, the tropopause and 1 Pa temperatures, the 0.3 Pa base and '
                           'the exobase rows at the exobase table\'s heats (curved_light.exobase): exobase temperature '
                           'and radius, molecular and atomic-oxygen loss.')
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    out['producer'] = dict(domain='atmosphere', files={f'middle_atmosphere/{Path(m.__file__).name}': digest(m.__file__)
                                                         for m in (sys.modules[__name__], cv, eq, escape, runner, ch)})
    OUT.parent.mkdir(parents=True, exist_ok=True)
    tmp = OUT.with_suffix('.tmp')
    tmp.write_text(json.dumps(out, indent=1, default=float) + '\n')
    tmp.replace(OUT)
    cv.log(f'wrote {OUT.relative_to(ROOT)}')
    for name, row in out['cases'].items():
        cells = ' | '.join(f"{b}: base {row[b]['base_temperature_k']:.1f} K" for b in BOUNDS)
        print(f'{name}: {cells}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
