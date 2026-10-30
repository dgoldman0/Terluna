"""McSnow's lunar column with the air below a cloud base held at a set humidity over ice, for testing PlaSim's snow
sublimation against lunar flakes (29 October 2026; the record is mcsnow_subcloud.md beside this file).

McSnow's idealised profile (mo_atmo.f90, init_atmo_ideal) holds its 'LWC region' between h1 and h2 at water saturation
and below h1 dries the air at a fixed lapse of relative humidity, 0.05% per metre (mo_atmo_types' rhlr), from which it
also takes the surface's vapour pressure for the column's dry pressure. The snow work's reference columns set h1 to 0,
the cloud reaching the freezing level, so their snow never sublimates above it. This module adds one declared change to
the corrected build (mcsnow_diffusivity.py), off by default: with atmo_nml's ssat_below at -1 or above, the air below h1
holds that supersaturation over ice, and the column's pressure stays as with its cloud to the ground (the surface vapour
pressure at water saturation, as with h1 = 0), so the cloud above h1 is the reference columns' own. Off, the build
gives the corrected build's runs byte for byte (`parity`).

The box cases repeat a reference column (its recorded configuration and seed) with h1 at the cloud base of a design box
and ssat_below at the mean humidity over ice below it (session-tools/snow_sublimation subcloud_air.py on the research
drive, by its rule). McSnow gives each cell the supply over the column's depth, per metre, so the column's own supply
left as it is gives each cloud cell its share. The first test took the moderate column (col_q1_f1_base: 1 mm/h, the
reference crystal flux, 96 hours): sub_sea_b and sub_land_b; the first pair, sub_sea and sub_land, raised the supply by
the column's depth over the cloud's, which gave each cloud cell 1.32 and 2.03 times its share, and stands as the record
of that first reading. The second test writes snapshots every 30 minutes in place of every 6 hours, which leaves the
runs themselves unchanged (the writer works on copies): sub_sea_c and sub_land_c repeat sub_sea_b and sub_land_b,
sub_sea_l and sub_land_l take the light column (col_q0.3_f1_base: 0.3 mm/h, 132 hours), and ref_q0.3_f1, ref_q1_f1 and
ref_q3_f1 rerun the three references at the reference crystal flux (0.3, 1 and 3 mm/h) with ssat_below off.

    OPENBLAS_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_subcloud build
    OPENBLAS_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_subcloud parity
    OPENBLAS_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_subcloud run --case sub_sea_c

Runs go to research/runs/snow_agent_3/subcloud.
"""
from __future__ import annotations

import argparse
import contextlib
import json
from pathlib import Path

from atmosphere.snow import mcsnow_diffusivity as dv
from atmosphere.snow import mcsnow_growth_aggregation as growagg

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CACHE = ROOT / 'research/runs/snow_agent_3/subcloud'
MARKER = 'Terluna: the air below the cloud base held over ice'
OFF = -999.0
PARITY_CASE = 'col_sat_f0.01_base'
MODERATE = 'col_q1_f1_base'
LIGHT = 'col_q0.3_f1_base'
HEAVY = 'col_q3_f1_base'
SNAP_S = 1800.0                 # s: the second test's snapshot interval
# Each box's cloud base above its freezing level and the supersaturation over ice below it, by subcloud_air.py's rule:
# the lowest height above the box's freezing level where its mean humidity over ice reaches 0.98, and the mean humidity
# over ice below it, weighted by the air's mass (box_near_side_sea and box_0e_design over days 29.5-59.0).
BOXES = {'sea': dict(h1_m=3957.0, ssat_below=0.9469 - 1.0),
         'land': dict(h1_m=8254.0, ssat_below=0.9604 - 1.0)}
# Each case: the design box whose cloud base and sub-cloud air it takes (none for a reference rerun), the reference
# column it repeats, the first pair's raised supply, and the snapshot interval when it departs from the column's.
CASES = {'sub_sea': dict(box='sea', base=MODERATE, raised=True),
         'sub_land': dict(box='land', base=MODERATE, raised=True),
         'sub_sea_b': dict(box='sea', base=MODERATE), 'sub_land_b': dict(box='land', base=MODERATE),
         'sub_sea_c': dict(box='sea', base=MODERATE, dt_snap=SNAP_S),
         'sub_land_c': dict(box='land', base=MODERATE, dt_snap=SNAP_S),
         'sub_sea_l': dict(box='sea', base=LIGHT, dt_snap=SNAP_S),
         'sub_land_l': dict(box='land', base=LIGHT, dt_snap=SNAP_S),
         'ref_q0.3_f1': dict(box=None, base=LIGHT, dt_snap=SNAP_S),
         'ref_q1_f1': dict(box=None, base=MODERATE, dt_snap=SNAP_S),
         'ref_q3_f1': dict(box=None, base=HEAVY, dt_snap=SNAP_S)}


def _replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError('Expected one occurrence of: ' + old[:70])
    return text.replace(old, new)


NML_DECLARATION = '    ssat,         & ! [1] supersaturation in the upper part of troposphere\n'
NML_LIST = '  NAMELIST /atmo_nml/ atmo_type, ssat, h0, h1, h2, nh1, nh2, area_factor, atmofile, IGF_id\n'
NML_DEFAULT = '  ssat = 0.05_wp\n'
NML_PRINT = '    PRINT *, "atmo_nml::h1: ", h1\n'
ATMO_USE = ('  USE mo_atmo_nml,     ONLY: ssat, h0, h1, h2, atmo_type, area_factor, &\n'
            '                             atmofile\n')
ATMO_BELOW = ('    atmo%pv   = 1e-2_wp * rh * psatw\n'
              '    atmo%ssat = atmo%pv / psati - 1._wp\n'
              '    atmo%rh  = atmo%pv / psatw * 100._wp\n'
              '  END IF\n'
              '\n'
              '  ! surface vapor pressure\n'
              '  pv_srf = 1e-2_wp * (100._wp-h1*rhlr) * psatw\n'
              '#endif\n')


def subcloud_sources(world):
    """The corrected build's sources with the declared change: ssat_below in atmo_nml, read by init_atmo_ideal."""
    out = dict(dv.fixed_sources(world))
    if world != 'moon':
        return out
    nml = (growagg.SOURCE / 'src/mo_atmo_nml.f90').read_text()
    nml = _replace_once(nml, NML_DECLARATION, NML_DECLARATION
                        + '    ssat_below,   & ! [1] Terluna: supersaturation over ice held below h1; off below -1\n')
    nml = _replace_once(nml, NML_LIST, NML_LIST[:-1] + ', ssat_below\n')
    nml = _replace_once(nml, NML_DEFAULT, NML_DEFAULT + f'  ssat_below = {OFF:.1f}_wp\n')
    nml = _replace_once(nml, NML_PRINT, NML_PRINT
                        + '    IF ( ssat_below >= -1._wp ) PRINT *, "atmo_nml::ssat_below: ", ssat_below\n')
    out['mo_atmo_nml.f90'] = nml
    atmo = _replace_once(out['mo_atmo.f90'], ATMO_USE, ATMO_USE[:-1] + ', ssat_below\n')
    out['mo_atmo.f90'] = _replace_once(
        atmo, ATMO_BELOW,
        '    atmo%pv   = 1e-2_wp * rh * psatw\n'
        '    atmo%ssat = atmo%pv / psati - 1._wp\n'
        '    atmo%rh  = atmo%pv / psatw * 100._wp\n'
        f'    IF ( ssat_below >= -1._wp .AND. z <= h1 ) THEN   ! {MARKER}\n'
        '      atmo%pv   = (ssat_below + 1._wp) * psati\n'
        '      atmo%ssat = ssat_below\n'
        '      atmo%rh  = atmo%pv / psatw * 100._wp\n'
        '    END IF\n'
        '  END IF\n'
        '\n'
        '  ! surface vapor pressure\n'
        '  pv_srf = 1e-2_wp * (100._wp-h1*rhlr) * psatw\n'
        '  IF ( ssat_below >= -1._wp ) pv_srf = psatw   ! the pressure as with the cloud to the ground (h1 = 0)\n'
        '#endif\n')
    return out


@contextlib.contextmanager
def subcloud():
    """growagg's build and run machinery on the changed sources, in this module's cache; the Earth build stays the old
    cache's, as under mcsnow_diffusivity.corrected."""
    saved = growagg.patched_sources, growagg.CACHE, growagg.build_path
    old_build_path = growagg.build_path

    def build_path(world):
        if world != 'earth':
            return old_build_path(world)
        cache, growagg.CACHE = growagg.CACHE, dv.OLD
        try:
            return old_build_path(world)
        finally:
            growagg.CACHE = cache

    growagg.patched_sources, growagg.CACHE, growagg.build_path = subcloud_sources, CACHE, build_path
    try:
        yield
    finally:
        growagg.patched_sources, growagg.CACHE, growagg.build_path = saved


def build():
    with subcloud():
        record = growagg.build('moon')
        path = growagg.build_path('moon')
    return dict(path=str(path.relative_to(ROOT)), executable_sha256=record['executable_sha256'])


def parity():
    """The corrected build's col_sat_f0.01_base again on the changed build with ssat_below off: every kept output
    byte for byte."""
    old = json.loads((dv.CACHE / 'runs' / PARITY_CASE / 'run.json').read_text())
    with subcloud():
        new = growagg.run('parity_' + PARITY_CASE, 'moon', old['configuration'], seed=old['seed'])
    same = {name: old['kept_outputs'].get(name) == digest for name, digest in new['kept_outputs'].items()}
    return dict(case=PARITY_CASE, returncode=new['returncode'], identical=all(same.values()) and
                set(new['kept_outputs']) == set(old['kept_outputs']), files=same, wall_s=new['wall_s'])


def case_config(case):
    """The reference column's recorded configuration with the case's cloud base and sub-cloud air, its own supply or,
    in the first pair, the supply raised by the column's depth over the cloud's; and the case's snapshot interval."""
    spec = CASES[case]
    record = json.loads((dv.CACHE / 'runs' / spec['base'] / 'run.json').read_text())
    config = json.loads(json.dumps(record['configuration']))
    if spec['box'] is not None:
        box = BOXES[spec['box']]
        depth = config['grid_nml']['dom_top']
        config['atmo_nml'].update(h1=float(box['h1_m']), ssat_below=float(box['ssat_below']))
        if spec.get('raised'):
            config['sp_nml']['water_supply'] = float(config['sp_nml']['water_supply'] * depth / (depth - box['h1_m']))
    if 'dt_snap' in spec:
        config['sp_nml']['dt_snap'] = float(spec['dt_snap'])
    return config, record['seed']


def run_case(case):
    config, seed = case_config(case)
    with subcloud():
        out = growagg.run(case, 'moon', config, seed=seed)
    return dict(case=case, returncode=out['returncode'], wall_s=out['wall_s'],
                executable_sha256=out['executable_sha256'])


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('action', choices=('build', 'parity', 'run'))
    parser.add_argument('--case', choices=tuple(CASES))
    args = parser.parse_args()
    if args.action == 'build':
        print(json.dumps(build(), indent=1))
    elif args.action == 'parity':
        print(json.dumps(parity(), indent=1))
    else:
        print(json.dumps(run_case(args.case), indent=1))


if __name__ == '__main__':
    main()
