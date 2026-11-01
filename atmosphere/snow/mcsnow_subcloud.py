"""McSnow's lunar column with the air below a cloud base held at a set humidity over ice, for testing PlaSim's snow
sublimation against lunar flakes (29 October 2026; the record is mcsnow_subcloud.md beside this file).

McSnow's idealised profile (mo_atmo.f90, init_atmo_ideal) holds its 'LWC region' between h1 and h2 at water saturation
and below h1 dries the air at a fixed lapse of relative humidity, 0.05% per metre (mo_atmo_types' rhlr), from which it
also takes the surface's vapour pressure for the column's dry pressure. The snow work's reference columns set h1 to 0,
the cloud reaching the freezing level, so their snow never sublimates above it. This module adds one declared change to
the corrected build (mcsnow_diffusivity.py), off by default: with atmo_nml's ssat_below at -1 or above, the air below h1
holds that supersaturation over ice, and the column's pressure stays as with its cloud to the ground (the surface vapour
pressure at water saturation, as with h1 = 0), so the cloud above h1 is the reference columns' own. A second declared
change, also off by default, caps the aggregates: with sp_nml's flake_cap above zero, an aggregate whose span passes it
breaks into equal pieces, each a share of its mass and monomers, until each lies within the cap, the superparticle's
multiplicity carrying the pieces so that the ice it stands for is unchanged. McSnow has no breakup law, and without one
snow of 1 mm/h from -30 C grows aggregates metres across; the cap stands in with the snow work's survivable span (the
load analogy of turbulence_survival.md: the largest Earth flakes observed, 4-6 cm, carry the load of lunar flakes of
8-15 cm). Off, the build gives the corrected build's runs byte for byte (`parity`).

The box cases repeat a reference column (its recorded configuration and seed) with h1 at the cloud base of a design box
and ssat_below at the mean humidity over ice below it (session-tools/snow_sublimation subcloud_air.py on the research
drive, by its rule). McSnow gives each cell the supply over the column's depth, per metre, so the column's own supply
left as it is gives each cloud cell its share. The first test took the moderate column (col_q1_f1_base: 1 mm/h, the
reference crystal flux, 96 hours): sub_sea_b and sub_land_b; the first pair, sub_sea and sub_land, raised the supply by
the column's depth over the cloud's, which gave each cloud cell 1.32 and 2.03 times its share, and stands as the record
of that first reading. The second test writes snapshots every 30 minutes in place of every 6 hours, which leaves the
runs themselves unchanged (the writer works on copies): sub_sea_c and sub_land_c repeat sub_sea_b and sub_land_b,
sub_sea_l and sub_land_l take the light column (col_q0.3_f1_base: 0.3 mm/h, 132 hours), and ref_q0.3_f1, ref_q1_f1 and
ref_q3_f1 rerun the three references at the reference crystal flux (0.3, 1 and 3 mm/h) with ssat_below off. The
third test adds references built as the snow work builds its columns (growagg.column_config), light snow of 0.1 mm/h
from -15 C and four rates from -30 C, the coldest top McSnow's habit table allows, with box runs on them, two at the
dry end of the boxes' sub-cloud air and one from a top halfway between.

    OPENBLAS_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_subcloud build
    OPENBLAS_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_subcloud parity
    OPENBLAS_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_subcloud run --case sub_sea_c
    OPENBLAS_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_subcloud run --case rdeep_q1_f1 --smoke --store DIR

--smoke runs the case for 6 model hours under smoke_<case>; --store DIR places a new run's directory in DIR, linked
from the cache, so that large snapshot files can sit on another disk.

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
FIRST_BUILD = 'b5c6f8c1'          # the sub-cloud build before the flake cap, whose parity run keeps its name
OFF = -999.0
PARITY_CASE = 'col_sat_f0.01_base'
MODERATE = 'col_q1_f1_base'
LIGHT = 'col_q0.3_f1_base'
HEAVY = 'col_q3_f1_base'
DEEP = 'col_deep_q0.3_f1_base'
SNAP_S = 1800.0                 # s: the second and third tests' snapshot interval
SMOKE_S = 6.0 * 3600.0          # s: a smoke run's length
CAP_M = 0.14                    # m: the flake cap, the lunar span loaded as Earth's largest observed flakes (6 cm)
CAP_LOW_M = 0.10                # m: its sensitivity, the span loaded as 4 cm Earth flakes
SEED = 20261025                 # the reference columns' seed
# Reference columns the snow work did not run, built as it builds its columns: water supply (mm/h), crystal flux over
# the reference, duration, and the top's temperature (-15 C unless set). Durations follow the snow work's, about 3^0.27
# longer for each threefold lighter rate and 1.27 times longer from -30 C than from -15 C.
COLUMNS = {'q0.1_f1': dict(water_mm_h=0.1, flux_factor=1.0, duration_h=180.0),
           'deep_q0.1_f1': dict(water_mm_h=0.1, flux_factor=1.0, duration_h=228.0, top_C=-30.0),
           'deep_q1_f1': dict(water_mm_h=1.0, flux_factor=1.0, duration_h=120.0, top_C=-30.0),
           'deep_q3_f1': dict(water_mm_h=3.0, flux_factor=1.0, duration_h=96.0, top_C=-30.0),
           'mid_q1_f1': dict(water_mm_h=1.0, flux_factor=1.0, duration_h=108.0, top_C=-22.5)}
# Each box's cloud base above its freezing level and the supersaturation over ice below it, by subcloud_air.py's rule:
# the lowest height above the box's freezing level where its mean humidity over ice reaches 0.98, and the mean humidity
# over ice below it, weighted by the air's mass (box_near_side_sea and box_0e_design over days 29.5-59.0).
# The dry variants hold the air at the flux-weighted 10th percentile of the humidity over ice that the boxes' falling
# ice meets below the cloud base (session-tools box_snow_rates.py, the same span).
BOXES = {'sea': dict(h1_m=3957.0, ssat_below=0.9469 - 1.0),
         'land': dict(h1_m=8254.0, ssat_below=0.9604 - 1.0),
         'sea_dry': dict(h1_m=3957.0, ssat_below=0.9033 - 1.0),
         'land_dry': dict(h1_m=8254.0, ssat_below=0.8909 - 1.0)}
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
         'ref_q3_f1': dict(box=None, base=HEAVY, dt_snap=SNAP_S),
         'ref_q0.1_f1': dict(box=None, column='q0.1_f1', dt_snap=SNAP_S),
         'rdeep_q0.1_f1': dict(box=None, column='deep_q0.1_f1', dt_snap=SNAP_S),
         'rdeep_q0.3_f1': dict(box=None, base=DEEP, dt_snap=SNAP_S),
         'rdeep_q1_f1': dict(box=None, column='deep_q1_f1', dt_snap=SNAP_S),
         'rdeep_q3_f1': dict(box=None, column='deep_q3_f1', dt_snap=SNAP_S),
         'sdeep_sea_q1': dict(box='sea', column='deep_q1_f1', dt_snap=SNAP_S),
         'sdeep_land_q1': dict(box='land', column='deep_q1_f1', dt_snap=SNAP_S),
         'sdeep_sea_q0.1': dict(box='sea', column='deep_q0.1_f1', dt_snap=SNAP_S),
         'sdeep_land_q0.1': dict(box='land', column='deep_q0.1_f1', dt_snap=SNAP_S),
         'sub_sea_q0.1': dict(box='sea', column='q0.1_f1', dt_snap=SNAP_S),
         'sub_land_q0.1': dict(box='land', column='q0.1_f1', dt_snap=SNAP_S),
         'sub_sea_dry': dict(box='sea_dry', base=MODERATE, dt_snap=SNAP_S),
         'sdeep_land_dry': dict(box='land_dry', column='deep_q1_f1', dt_snap=SNAP_S),
         'smid_sea_q1': dict(box='sea', column='mid_q1_f1', dt_snap=SNAP_S),
         'sub_sea_q3': dict(box='sea', base=HEAVY, dt_snap=SNAP_S),
         'sdeep_land_q3': dict(box='land', column='deep_q3_f1', dt_snap=SNAP_S)}
# The fifth test repeats the references and box runs with the aggregates capped at the survivable span (cap_), the heavy
# box runs only so, and three references at the cap's lower end (capl_); a smoke case caps at 5 mm so that its pieces
# form within the smoke's 6 hours.
CAPPED = ('ref_q0.1_f1', 'ref_q0.3_f1', 'ref_q1_f1', 'ref_q3_f1', 'rdeep_q0.1_f1', 'rdeep_q0.3_f1', 'rdeep_q1_f1',
          'rdeep_q3_f1', 'sub_sea_c', 'sub_land_c', 'sub_sea_l', 'sub_land_l', 'sdeep_sea_q1', 'sdeep_land_q1',
          'sdeep_sea_q0.1', 'sdeep_land_q0.1', 'sub_sea_q0.1', 'sub_land_q0.1', 'sub_sea_dry', 'sdeep_land_dry',
          'smid_sea_q1', 'sub_sea_q3', 'sdeep_land_q3')
CASES.update({'cap_' + name: dict(CASES[name], cap=CAP_M) for name in CAPPED})
CASES.update({'capl_' + name: dict(CASES[name], cap=CAP_LOW_M)
              for name in ('ref_q1_f1', 'rdeep_q0.3_f1', 'rdeep_q1_f1')})
CASES['cap5mm_ref_q1_f1'] = dict(CASES['ref_q1_f1'], cap=0.005)
for name in ('sub_sea_q3', 'sdeep_land_q3'):         # heavy snow runs only with the cap
    CASES.pop(name)


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


CAP_DECLARATION = '    depo_residual,    &  ! [kg] a sublimating particle below this mass is removed; 0 off (Terluna)\n'
CAP_LISTS = '            box_sub, box_settle, skip_bulk, depo_residual, depo_censor\n'
CAP_DEFAULT = '  depo_censor    = 0.0_wp\n'
CAP_USE = '                                & box_sub, box_settle, skip_bulk\n'
CAP_STEP = '            ! motion of SPs ! TODO ensure that only adding is done after sedi as it breaks the cell order\n'
CAP_ROUTINE = '''#ifndef WDROPS
    !
    ! Terluna: an aggregate whose span passes flake_cap breaks into equal pieces, each a share of its mass and
    ! monomers, until each lies within the cap; the multiplicity carries the pieces, so the ice is unchanged.
    !
    SUBROUTINE cap_flake( sp )
      TYPE(t_sp), INTENT(inout) :: sp
      INTEGER(i8) :: f
      INTEGER     :: n
      CALL m2d_particle( sp )          ! the span after this step's growth, which deposition leaves unset
      DO n = 1, 20
        IF ( sp%d <= flake_cap .OR. sp%mm < 2.0_wp ) EXIT
        f = MIN( MAX(2_i8, CEILING((sp%d/flake_cap)**2.1_wp, i8)), INT(sp%mm, i8) )
        sp%xi  = sp%xi * f
        sp%m_i = sp%m_i / REAL(f, wp)
        sp%m_r = sp%m_r / REAL(f, wp)
        sp%v_r = sp%v_r / REAL(f, wp)
        sp%m_f = sp%m_f / REAL(f, wp)
        sp%m_w = sp%m_w / REAL(f, wp)
        sp%V_i = sp%V_i / REAL(f, wp)
        sp%mm  = sp%mm / REAL(f, wp)
        CALL m2d_particle( sp )
      END DO
    END SUBROUTINE cap_flake
#endif

'''


def cap_sources(text_nml, text_loop):
    """The flake cap's declared lines in mo_sp_nml.f90 and mo_loop_array.f90."""
    nml = _replace_once(text_nml, CAP_DECLARATION, CAP_DECLARATION
                        + '    flake_cap,        &  ! [m] an aggregate passing this span breaks into pieces within it; '
                          '0 off (Terluna)\n')
    if nml.count(CAP_LISTS) != 2:
        raise ValueError('Expected the namelist and public lists once each')
    nml = nml.replace(CAP_LISTS, CAP_LISTS[:-1] + ', flake_cap\n')
    nml = _replace_once(nml, CAP_DEFAULT, CAP_DEFAULT + '  flake_cap      = 0.0_wp\n')
    loop = _replace_once(text_loop, CAP_USE, CAP_USE[:-1] + ', flake_cap\n')
    loop = _replace_once(loop, CAP_STEP, '#ifndef WDROPS\n'
                         '            IF ( flake_cap > 0.0_wp ) CALL cap_flake( sp )   ! Terluna: the survivable span\n'
                         '#endif\n' + CAP_STEP)
    loop = _replace_once(loop, '    SUBROUTINE tic( t1 )\n', CAP_ROUTINE + '    SUBROUTINE tic( t1 )\n')
    return nml, loop


def subcloud_sources(world):
    """The corrected build's sources with the declared changes: ssat_below in atmo_nml, read by init_atmo_ideal, and
    flake_cap in sp_nml, applied in the step loop."""
    out = dict(dv.fixed_sources(world))
    if world != 'moon':
        return out
    out['mo_sp_nml.f90'], out['mo_loop_array.f90'] = cap_sources(out['mo_sp_nml.f90'], out['mo_loop_array.f90'])
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
    """The corrected build's col_sat_f0.01_base again on the changed build with ssat_below and flake_cap off: every
    kept output byte for byte. The run is named by the build, so each build's parity runs afresh."""
    old = json.loads((dv.CACHE / 'runs' / PARITY_CASE / 'run.json').read_text())
    with subcloud():
        exe = growagg.verify_build('moon')['executable_sha256']
        name = 'parity_' + PARITY_CASE if exe.startswith(FIRST_BUILD) else f'parity_{exe[:12]}_{PARITY_CASE}'
        new = growagg.run(name, 'moon', old['configuration'], seed=old['seed'])
    same = {name: old['kept_outputs'].get(name) == digest for name, digest in new['kept_outputs'].items()}
    return dict(case=PARITY_CASE, returncode=new['returncode'], identical=all(same.values()) and
                set(new['kept_outputs']) == set(old['kept_outputs']), files=same, wall_s=new['wall_s'])


def case_config(case):
    """The reference column's recorded configuration with the case's cloud base and sub-cloud air, its own supply or,
    in the first pair, the supply raised by the column's depth over the cloud's; and the case's snapshot interval."""
    spec = CASES[case]
    if 'column' in spec:
        config = growagg.column_config(band=growagg.band_factors(None), snap_h=6.0, **COLUMNS[spec['column']])
        seed = SEED
    else:
        record = json.loads((dv.CACHE / 'runs' / spec['base'] / 'run.json').read_text())
        config, seed = json.loads(json.dumps(record['configuration'])), record['seed']
    if spec['box'] is not None:
        box = BOXES[spec['box']]
        depth = config['grid_nml']['dom_top']
        config['atmo_nml'].update(h1=float(box['h1_m']), ssat_below=float(box['ssat_below']))
        if spec.get('raised'):
            config['sp_nml']['water_supply'] = float(config['sp_nml']['water_supply'] * depth / (depth - box['h1_m']))
    if 'dt_snap' in spec:
        config['sp_nml']['dt_snap'] = float(spec['dt_snap'])
    if 'cap' in spec:
        config['sp_nml']['flake_cap'] = float(spec['cap'])
    return config, seed


def run_case(case, smoke=False, store=None):
    """One case, or with smoke its first 6 model hours under smoke_<case>; with store, a new run's directory is made
    in that directory and linked from the cache."""
    config, seed = case_config(case)
    name = case
    if smoke:
        name = 'smoke_' + case
        config['general_nml'].update(time_end=SMOKE_S, dt_1dprof=SMOKE_S, dt_2dprof=SMOKE_S)
    path = CACHE / 'runs' / name
    if store is not None and not path.exists() and not path.is_symlink():
        target = Path(store) / name
        target.mkdir(parents=True, exist_ok=True)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.symlink_to(target.resolve())
    with subcloud():
        out = growagg.run(name, 'moon', config, seed=seed)
    return dict(case=name, returncode=out['returncode'], wall_s=out['wall_s'],
                executable_sha256=out['executable_sha256'])


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('action', choices=('build', 'parity', 'run'))
    parser.add_argument('--case', choices=tuple(CASES))
    parser.add_argument('--smoke', action='store_true')
    parser.add_argument('--store')
    args = parser.parse_args()
    if args.action == 'build':
        print(json.dumps(build(), indent=1))
    elif args.action == 'parity':
        print(json.dumps(parity(), indent=1))
    else:
        print(json.dumps(run_case(args.case, smoke=args.smoke, store=args.store), indent=1))


if __name__ == '__main__':
    main()
