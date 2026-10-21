"""Snow populations with shorter prescribed saturated growth paths.

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m atmosphere.snow.growth_depth

Reuses population.run_case without changing its frozen equations. The four
Moon cases keep 100/L input-flux normalization, seed 0, representative dendrite
area and 1 or 3 mm/h imposed condensation. Cloud tops at -5 and -10 C shorten
the growth path while preserving total supplied water per area per time. This
is a controlled depth sensitivity, not a cloud-frequency or habit prediction.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from atmosphere.snow import flakes as sf, population as pop
from atmosphere.snow.habits import GEOMETRIES
from shared.provenance import constants_used

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'atmosphere/snow/results/growth_depth.json'
CACHE = ROOT / 'research/runs/snow_growth_depth'
SCHEMA = 'terluna.atmosphere.snow-growth-depth/1'
FILES = ('atmosphere/snow/growth_depth.py',) + pop.FILES
GEOMETRY = 'bf95_dmax_dendrites_area'


def cases():
    return [dict(id=f'moon_{rate:g}_100_{GEOMETRY}_s0_top{abs(top):g}',
                 world='moon', snowfall_set_mm_h=rate, crystals_per_litre=100.,
                 geometry=GEOMETRY, seed=0, top_c=top)
            for top in (-5., -10.) for rate in (1., 3.)]


def cache_identity(spec, files, inputs, constants):
    identity = dict(spec=spec, files=files, inputs=inputs, constants=constants)
    return hashlib.sha256(json.dumps(identity, sort_keys=True, allow_nan=False).encode()).hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(value, indent=1, allow_nan=False) + '\n')
    temporary.replace(path)


def analyse(selected=None):
    moon, binding = sf.moon_column()
    files = {path: pop.sha(ROOT/path) for path in FILES}
    constants = constants_used(FILES)
    rows = []
    CACHE.mkdir(parents=True, exist_ok=True)
    for spec in cases():
        if selected and spec['id'] not in selected:
            continue
        key = cache_identity(spec, files, binding['inputs'], constants)
        cache = CACHE / f'{spec["id"]}.json'
        previous = json.loads(cache.read_text()) if cache.exists() else None
        if previous and previous['key'] == key:
            row = previous['result']
            print(f'cached {spec["id"]}', flush=True)
        else:
            print(f'computing {spec["id"]}', flush=True)
            row = pop.run_case(spec, moon, GEOMETRIES[spec['geometry']], binding)
            row['cloud_depth_m'] = float(moon.height_of(sf.T0 + spec['top_c']) - moon.z_freeze)
            write_json(cache, dict(key=key, result=row))
        rows.append(row)
        print(spec['id'], {name: round(value, 4) if isinstance(value, float) else value
                          for name, value in row['summary'].items()}, flush=True)
    # Reject a mixed-source product if another task changed a producer during a run.
    for path, expected in files.items():
        if pop.sha(ROOT/path) != expected:
            raise RuntimeError(f'Producer changed during growth-depth run: {path}')
    if constants_used(FILES) != constants:
        raise RuntimeError('Shared constants changed during growth-depth run')
    for path, expected in binding['inputs'].items():
        if pop.sha(ROOT/path) != expected:
            raise RuntimeError(f'Column input changed during growth-depth run: {path}')
    return dict(schema=SCHEMA, runs=rows,
                producer=dict(domain='atmosphere', files=files, inputs=binding['inputs'], constants=constants),
                evidence='Four conditional one-dimensional, persistent, liquid-saturated cloud calculations. '
                         'Only cloud-top temperature and depth change relative to the central population cases. '
                         'Total imposed water supply is fixed at each compared rate; condensation per unit depth '
                         'therefore increases in a shorter cloud. Single seed, inherited collision/fragmentation laws, '
                         'no resolved turbulence or large-aggregate mechanical fracture.',
                reading_rule='Runs have the snow-population/1 run structure, with cloud_depth_m added. '
                             'Compare top -5 and -10 C to the same-rate top -15 C population case; do not infer '
                             'cloud occurrence from these imposed depths. Summaries weight arriving mass; '
                             'number and area bins retain their original moments. The 100/L parameter normalizes '
                             'injected number flux with the speed of a 0.2 mm crystal, not a fixed size-resolved '
                             'concentration boundary. Dendrite-area geometry is held for the controlled comparison; '
                             'warmer cloud tops do not establish dendritic crystal production. These runs retain '
                             'a saturated column rather than the intermittent dry gaps in archived weather.',
                comparison=dict(top_temperatures_c=[-5., -10.], reference_top_temperature_c=-15.,
                                reference_product='atmosphere/snow/results/population.json',
                                fixed=['total imposed condensation per area per time', '100/L flux normalization',
                                       'seed 0', GEOMETRY, 'default population discretization']))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', action='append', help='Case ID; selected output stays in the ignored cache directory')
    args = parser.parse_args()
    product = analyse(args.case)
    if not product['runs']:
        raise ValueError('No matching growth-depth cases')
    path = CACHE/'partial.json' if args.case else OUT
    write_json(path, product)
    print(f'wrote {path.relative_to(ROOT)}', flush=True)


if __name__ == '__main__':
    main()
