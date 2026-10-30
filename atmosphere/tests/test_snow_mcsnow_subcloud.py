"""McSnow's column with the air below a cloud base held over ice (atmosphere/snow/mcsnow_subcloud.py): the change stands
on the corrected build and touches only its declared lines, it is off by default, the Earth build is untouched, and the
cases repeat their reference columns with each box's cloud base, sub-cloud air and the cloud's supply, departing only
where declared."""
import difflib
import json

import pytest

from atmosphere.snow import mcsnow_diffusivity as dv
from atmosphere.snow import mcsnow_growth_aggregation as grow
from atmosphere.snow import mcsnow_subcloud as sc


def _sources_or_skip():
    if not (grow.SOURCE / 'src/mo_atmo_nml.f90').exists():
        pytest.skip('pinned McSnow source copy is not installed')


def test_the_change_touches_only_its_declared_lines_and_leaves_earth_alone():
    _sources_or_skip()
    new, old = sc.subcloud_sources('moon'), dv.fixed_sources('moon')
    changed = []
    for name in sorted(set(new) | set(old)):
        a = old.get(name, (grow.SOURCE / 'src' / name).read_text())
        changed += [line for line in difflib.unified_diff(a.splitlines(), new[name].splitlines(), lineterm='', n=0)
                    if line.startswith(('+', '-')) and not line.startswith(('+++', '---'))]
    expected = [
        '-                             atmofile',
        '+                             atmofile, ssat_below',
        f'+    IF ( ssat_below >= -1._wp .AND. z <= h1 ) THEN   ! {sc.MARKER}',
        '+      atmo%pv   = (ssat_below + 1._wp) * psati',
        '+      atmo%ssat = ssat_below',
        '+      atmo%rh  = atmo%pv / psatw * 100._wp',
        '+    END IF',
        '+  IF ( ssat_below >= -1._wp ) pv_srf = psatw   ! the pressure as with the cloud to the ground (h1 = 0)',
        '+    ssat_below,   & ! [1] Terluna: supersaturation over ice held below h1; off below -1',
        '-  NAMELIST /atmo_nml/ atmo_type, ssat, h0, h1, h2, nh1, nh2, area_factor, atmofile, IGF_id',
        '+  NAMELIST /atmo_nml/ atmo_type, ssat, h0, h1, h2, nh1, nh2, area_factor, atmofile, IGF_id, ssat_below',
        '+  ssat_below = -999.0_wp',
        '+    IF ( ssat_below >= -1._wp ) PRINT *, "atmo_nml::ssat_below: ", ssat_below',
    ]
    assert sorted(changed) == sorted(expected)
    assert sc.subcloud_sources('earth') == dv.fixed_sources('earth')
    assert f'  ssat_below = {sc.OFF:.1f}_wp' in new['mo_atmo_nml.f90'] and sc.OFF < -1.0      # off by default
    assert all(len(line) <= 132 for name in ('mo_atmo_nml.f90', 'mo_atmo.f90') for line in new[name].splitlines())


def test_the_cases_repeat_their_reference_columns_with_each_boxs_air():
    records = {}
    for name in {spec['base'] for spec in sc.CASES.values() if 'base' in spec} | {sc.DEEP}:
        path = dv.CACHE / 'runs' / name / 'run.json'
        if not path.exists():
            pytest.skip('the corrected reference runs are not at hand')
        records[name] = json.loads(path.read_text())
    assert {sc.MODERATE, sc.LIGHT, sc.HEAVY, sc.DEEP} <= set(records)
    # the snow work's own column builder rebuilds its recorded references, so the new columns are of their family
    for name, options in ((sc.LIGHT, dict(water_mm_h=0.3, flux_factor=1.0, duration_h=132.0)),
                          (sc.MODERATE, dict(water_mm_h=1.0, flux_factor=1.0, duration_h=96.0)),
                          (sc.HEAVY, dict(water_mm_h=3.0, flux_factor=1.0, duration_h=72.0)),
                          (sc.DEEP, dict(water_mm_h=0.3, flux_factor=1.0, duration_h=168.0, top_C=-30.0))):
        assert grow.column_config(band=grow.band_factors(None), snap_h=6.0, **options) == records[name]['configuration']
        assert records[name]['seed'] == sc.SEED
    for case, spec in sc.CASES.items():
        if 'column' in spec:
            reference = dict(seed=sc.SEED, configuration=grow.column_config(band=grow.band_factors(None), snap_h=6.0,
                                                                            **sc.COLUMNS[spec['column']]))
        else:
            reference = records[spec['base']]
        before = reference['configuration']
        config, seed = sc.case_config(case)
        assert seed == reference['seed']
        depth = config['grid_nml']['dom_top']
        supply = before['sp_nml']['water_supply']
        if spec['box'] is None:                         # a reference rerun: its own air
            assert config['atmo_nml'] == before['atmo_nml'] and 'ssat_below' not in config['atmo_nml']
        else:
            box = sc.BOXES[spec['box']]
            assert config['atmo_nml']['h1'] == box['h1_m'] and config['atmo_nml']['ssat_below'] == box['ssat_below']
            assert -1.0 <= box['ssat_below'] < 0.0 < box['h1_m'] < depth
            assert {k: v for k, v in config['atmo_nml'].items() if k not in ('h1', 'ssat_below')} == \
                {k: v for k, v in before['atmo_nml'].items() if k != 'h1'}
        if spec.get('raised'):                          # the first pair's, raised by the depth over the cloud's
            assert config['sp_nml']['water_supply'] == pytest.approx(supply * depth / (depth - box['h1_m']), rel=1e-12)
        else:                                           # each cloud cell has its column's own supply per metre
            assert config['sp_nml']['water_supply'] == supply
        assert config['sp_nml']['dt_snap'] == spec.get('dt_snap', before['sp_nml']['dt_snap'])
        rest = lambda c: {k: v for k, v in c.items() if k not in ('atmo_nml', 'sp_nml')}
        assert rest(config) == rest(before)
        assert {k: v for k, v in config['sp_nml'].items() if k not in ('water_supply', 'dt_snap')} == \
            {k: v for k, v in before['sp_nml'].items() if k not in ('water_supply', 'dt_snap')}
    # the second and third tests' runs write snapshots every 30 minutes, the reference columns every 6 hours
    assert all(spec.get('dt_snap') == sc.SNAP_S for case, spec in sc.CASES.items()
               if case not in ('sub_sea', 'sub_land', 'sub_sea_b', 'sub_land_b'))
    # the dry variants keep their box's cloud base and are drier below it
    for box in ('sea', 'land'):
        assert sc.BOXES[box + '_dry']['h1_m'] == sc.BOXES[box]['h1_m']
        assert -1.0 < sc.BOXES[box + '_dry']['ssat_below'] < sc.BOXES[box]['ssat_below'] < 0.0
    assert all(records[n]['configuration']['sp_nml']['dt_snap'] % sc.SNAP_S == 0 for n in records)
