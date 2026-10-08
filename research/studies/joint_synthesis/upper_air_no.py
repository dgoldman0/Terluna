"""Lightning's NO in the upper air: how much NO at the thermal column's base would change the shield's loss screen.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.joint_synthesis.upper_air_no
    # -> results/upper_air_no.json

The loss response's thermal column (atmosphere/loss_response) radiates NO's 5.3 um band with the mole fraction the
middle-atmosphere case holds at 0.3 Pa, which for the design column under the titania stack is about 1e-16: that
column makes almost no NO. Under the shield the air has no oxidant to remove lightning's NO (results/air_chemistry.json),
so it can mix upward. This study holds NO at the base at 1, 10 and 100 ppb in place of the case's own value and reruns
the domain's state sweep (model.euv_sweep) for the design shield, both treatments of near-infrared heat that bound it,
quiet Sun and solar maximum, comparing the exobase temperature and the molecular loss at the swarm's standard gap
levels. Nothing in the domain is changed: the base value enters through infrared.base_mixing in this process only.
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path

from atmosphere.loss_response import infrared as ir, model as lm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.joint-synthesis-upper-air-no/1'
OUT = HERE / 'results' / 'upper_air_no.json'
NO_LEVELS = (None, 1e-9, 1e-8, 1e-7)
TREATMENTS = ('collisional', 'all_heats')
ACTIVITIES = ('quiet', 'solar_maximum')
GAPS = (1e-4, 2e-4, 3e-4)
EVIDENCE = ('The domain\'s thermal column and limb tables at traced states, with NO at the 0.3 Pa base set by hand; '
            'the solar-wind and exosphere steps that follow in the full chain are not rerun.')
READING_RULE = ('exobase_temperature_k and molecular_loss_kg_s at each gap level (UV transmission f) for each NO level; '
                '"case" is the middle-atmosphere case\'s own NO. Loss in kg/s of the thermal column\'s molecular escape.')


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def main(argv=None) -> int:
    original = ir.base_mixing
    out = {}
    for level in NO_LEVELS:
        label = 'case' if level is None else f'{level:g}'
        ir.base_mixing = original if level is None else (
            lambda case, species, base_pa, _v=level: _v if species == 'NO' else original(case, species, base_pa))
        lm._COLUMN.clear(); lm._BASE_RADIUS.clear()
        out[label] = {}
        for treatment in TREATMENTS:
            for activity in ACTIVITIES:
                rows = lm.euv_sweep('titania_stack', treatment, activity, glow=True, leaks=GAPS)
                out[label][f'{treatment}/{activity}'] = [
                    dict(gap=r['leak_fraction'], status=r['status'], exobase_temperature_k=round(r.get('exobase_temperature_k', float('nan')), 1),
                         molecular_loss_kg_s=r.get('molecular_loss_kg_s')) for r in rows]
                print(label, treatment, activity, [(r['leak_fraction'], r['status'], round(r.get('exobase_temperature_k', 0), 1),
                                                    r.get('molecular_loss_kg_s')) for r in rows], flush=True)
    ir.base_mixing = original
    files = {str(Path(p).relative_to(ROOT)): digest(p) for p in (__file__, lm.__file__, ir.__file__)}
    OUT.write_text(json.dumps(dict(schema=SCHEMA, producer=dict(study='joint_synthesis', files=files),
                                   evidence=EVIDENCE, reading_rule=READING_RULE, no_levels=[l for l in NO_LEVELS if l],
                                   sweeps=out), indent=1, default=str) + '\n')
    print(OUT)
    return 0


if __name__ == '__main__':
    sys.exit(main())
