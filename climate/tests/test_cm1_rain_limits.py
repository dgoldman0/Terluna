"""The rain-limits patch: it applies to CM1's own source on every Morrison build's series, keeps CM1's limits at its
default, scales them from its namelist as the runner's mirror does, and the runner writes the derived lunar values."""
import hashlib
import json
import math
import re

import pytest

from climate.crm import cm1_run as c
from shared.constants import (CM1_GRAVITY, CM1_RAIN_BREAKUP_CUTOFF, CM1_RAIN_BREAKUP_RATE, CM1_RAIN_FALL_CAP,
                              CM1_RAIN_MEAN_DIAMETER_CAP, MOON_SURFACE_GRAVITY)

SOURCE = c.CM1_HOME / c.SOURCE['version'] / 'src'
needs_source = pytest.mark.skipif(not (SOURCE / 'morrison.F').exists(), reason=f'CM1 is not unpacked in {SOURCE}')
PRODUCT_PATH = c.HERE.parents[1] / c.RAIN_LIMITS_PRODUCT
PRODUCT = json.loads(PRODUCT_PATH.read_text())
FILES = ('input.F', 'cm1.F', 'param.F', 'morrison.F')


def patched(series):
    """CM1's own files after a series of patches, keyed by name."""
    texts = {}
    for name, marker, edits in series:
        if name not in texts:
            texts[name] = (SOURCE / name).read_text(encoding='latin-1')
        texts[name] = c.apply_patch(texts[name], name, marker, edits)
    return texts


def rain_block(text):
    start = text.index('         ! Terluna rain limits (the header of this file)')
    return text[start:text.index('\n\n', start)]


def test_the_patch_is_one_edit_set_per_file_each_adding_its_marker():
    assert [p[0] for p in c.RAIN_LIMITS_PATCHES] == list(FILES)
    for name, marker, edits in c.RAIN_LIMITS_PATCHES:
        assert marker == c.RAIN_LIMITS_MARKER and any(marker in new for _, new, _ in edits), name
    assert c.BUILD_PATCHES['moon_omp_rain_limits'] == c.RAIN_LIMITS_PATCHES
    assert c.BUILDS['moon_omp_rain_limits'] == c.BUILDS['moon_omp']
    # No build before it takes it, so their executables and runs stay as they were; the rain-law build stands on it.
    descendants = ('moon_omp_rain_law', 'moon_omp_rain_law_pbl', 'moon_omp_rain_law_pbl_wind',
                   'moon_omp_rain_law_pbl_wind_gust')
    assert all(c.RAIN_LIMITS_PATCHES[0] not in patches for label, patches in c.BUILD_PATCHES.items()
               if label not in ('moon_omp_rain_limits', *descendants))
    for label in descendants:
        assert c.BUILD_PATCHES[label][:len(c.RAIN_LIMITS_PATCHES)] == c.RAIN_LIMITS_PATCHES
    assert not any(p in c.PATCHES for p in c.RAIN_LIMITS_PATCHES)


def test_the_header_says_what_it_changes_its_default_and_its_source():
    header = c.RAIN_LIMITS_HEADER
    assert all(line.startswith('!') and len(line) <= 120 for line in header.splitlines())
    for phrase in ('What it changes.', 'Default. brkup_gscale = 0', 'Source.', 'github.com/dgoldman0/Terluna',
                   'RAIN_LIMITS_PATCHES', 'atmosphere/rain/breakup.py', '300 um', '2300/m', '2.8 mm', 'VCAPR',
                   'sqrt(9.81/g)', '9.1*(g/9.81)**0.25', 'Straub et al. 2010'):
        assert phrase in ' '.join(header.split()), phrase
    for macro in ('TERLUNA_G', 'TERLUNA_CO2', 'TERLUNA_O2', '//', '/*', '#'):     # cpp would read these
        assert macro not in header, macro


@needs_source
def test_the_patch_applies_once_to_cm1s_own_source():
    texts = patched([*c.PATCHES, *c.RAIN_LIMITS_PATCHES])
    for name, marker, edits in c.RAIN_LIMITS_PATCHES:
        assert c.apply_patch(texts[name], name, marker, edits) == texts[name], name     # idempotent
    morrison = texts['morrison.F']
    assert morrison.startswith(c.RAIN_LIMITS_HEADER)
    assert morrison.count('dum=2.-exp(KB_r*(1./lamr(k)-dum1))') == 2 and 'exp(2300.' not in morrison
    assert 'use input, only : brkup_gscale,brkup_rhor,brkup_rhoexp,dowr,outfile' in morrison
    assert 'l_inf,brkup_db,brkup_gscale,brkup_rhor,brkup_rhoexp,ndcnst,nt_c' in texts['input.F']
    assert texts['param.F'].count('MPI_BCAST(brkup_') == 4


@needs_source
def test_the_patch_applies_on_every_morrison_builds_series():
    for label, (gravity, mode) in c.BUILDS.items():
        if label in c.BUILD_FILES:                         # the electrified builds run NSSL, not Morrison
            continue
        series = [*c.PATCHES, *[p for p in c.BUILD_PATCHES.get(label, []) if p[0] != 'Makefile'],
                  *([] if label == 'moon_omp_rain_limits' else c.RAIN_LIMITS_PATCHES)]
        texts = patched(series)
        assert c.RAIN_LIMITS_MARKER in texts['morrison.F'], label


@needs_source
def test_the_default_keeps_cm1s_own_limits():
    core = patched(c.PATCHES)['morrison.F']
    morrison = patched([*c.PATCHES, *c.RAIN_LIMITS_PATCHES])['morrison.F']
    # CM1's own lines stand as they were.
    for line in ('         LAMMAXR = 1./20.E-6\n', '         LAMMINR = 1./2800.E-6\n',
                 '         cm1db = 300.e-6  ! default cutoff diameter for breakup\n',
                 '         D_b = cm1db  ! cutoff diameter for breakup\n',
                 '         VCAPR = 9.1*(G/9.81)**((BR+1.)/3.)\n'):
        assert core.count(line) == 1 and morrison.count(line) == 1, line
    # Outside the IF on brkup_gscale = 1 the block sets only the rate, to CM1's 2300; inside, the four limits.
    block = rain_block(morrison)
    outside = re.sub(r'IF \(brkup_gscale\.EQ\.1\) THEN.*?END IF', '', block, flags=re.S)
    assignments = [line.strip() for line in outside.splitlines() if re.match(r'^\s+[A-Za-z_]\w* =', line)]
    assert assignments == ['KB_r = 2300.']
    inside = block[block.index('IF (brkup_gscale.EQ.1) THEN'):block.index('END IF')]
    for expression in ('terluna_scap = SQRT(9.81/G)', 'terluna_scol = terluna_scap*brkup_rhor**brkup_rhoexp',
                       'D_b = D_b*terluna_scol', 'KB_r = KB_r/terluna_scol', 'LAMMINR = LAMMINR/terluna_scap',
                       'VCAPR = 9.1*(G/9.81)**0.25'):
        assert expression in inside, expression
    # Undoing the patch's edits gives back the core's source: it adds nothing else.
    undone = morrison[len(c.RAIN_LIMITS_HEADER):]
    for old, new, _ in reversed(next(p for p in c.RAIN_LIMITS_PATCHES if p[0] == 'morrison.F')[2][1:]):
        undone = undone.replace(new, old)
    assert undone == core
    cm1 = patched([*c.PATCHES, *c.RAIN_LIMITS_PATCHES])['cm1.F']
    assert '      brkup_gscale = 0' in cm1 and '      brkup_rhor = 1.0\n      brkup_rhoexp = 0.0\n' in cm1


def test_the_mirror_gives_cm1s_own_limits_at_the_default_for_any_gravity():
    for gravity in (CM1_GRAVITY, MOON_SURFACE_GRAVITY, 3.7):
        for settings in ({}, dict(brkup_db=0.0, brkup_gscale=0, brkup_rhor=1.3, brkup_rhoexp=0.7)):
            limits = c.rain_limits(gravity, settings)
            assert limits['breakup_cutoff_m'] == CM1_RAIN_BREAKUP_CUTOFF
            assert limits['breakup_rate_per_m'] == CM1_RAIN_BREAKUP_RATE
            assert limits['mean_diameter_cap_m'] == CM1_RAIN_MEAN_DIAMETER_CAP
            assert limits['fall_speed_cap_m_s'] == pytest.approx(CM1_RAIN_FALL_CAP * (gravity / CM1_GRAVITY)**0.6)
            assert limits['equilibrium_m'] == pytest.approx(0.6014e-3, abs=1e-7)
    # Switched on at Earth's gravity and in Earth's air it gives CM1's own limits as well.
    earth = c.rain_limits(CM1_GRAVITY, dict(brkup_gscale=1, brkup_rhor=1.0, brkup_rhoexp=0.69))
    assert earth['breakup_cutoff_m'] == pytest.approx(CM1_RAIN_BREAKUP_CUTOFF)
    assert earth['mean_diameter_cap_m'] == pytest.approx(CM1_RAIN_MEAN_DIAMETER_CAP)
    assert earth['fall_speed_cap_m_s'] == pytest.approx(CM1_RAIN_FALL_CAP)
    with pytest.raises(ValueError):
        c.rain_limits(MOON_SURFACE_GRAVITY, dict(brkup_gscale=1, brkup_rhor=0.0, brkup_rhoexp=0.69))


def test_the_lunar_settings_carry_the_derived_values():
    assert PRODUCT['schema'] == c.RAIN_LIMITS_SCHEMA
    lunar = PRODUCT['lunar_limits']
    settings, product = c.rain_limit_settings('lunar', 1.2 * c.SEA_LEVEL_PA)
    assert product == {c.RAIN_LIMITS_PRODUCT: hashlib.sha256(PRODUCT_PATH.read_bytes()).hexdigest()}
    assert settings == dict(brkup_db=CM1_RAIN_BREAKUP_CUTOFF, brkup_gscale=1, brkup_rhor=1.2,
                            brkup_rhoexp=lunar['namelist']['brkup_rhoexp'])
    assert settings['brkup_rhoexp'] == round(lunar['alpha'], 3)
    limits = c.rain_limits(MOON_SURFACE_GRAVITY, settings)
    moon = lunar['moon']
    assert limits['breakup_cutoff_m'] * 1e3 == pytest.approx(moon['breakup_cutoff_mm'], rel=1e-12)
    assert limits['breakup_rate_per_m'] == pytest.approx(moon['breakup_rate_per_m'], rel=1e-12)
    assert limits['equilibrium_m'] * 1e3 == pytest.approx(moon['equilibrium_mm'], rel=1e-12)
    assert limits['mean_diameter_cap_m'] * 1e3 == pytest.approx(moon['mean_diameter_cap_mm'], rel=1e-12)
    assert limits['fall_speed_cap_m_s'] == pytest.approx(moon['fall_speed_cap_m_s'], rel=1e-12)
    # The patch's scale at the design pressure stands within 0.5% of the drop physics' own at the equilibrium.
    assert lunar['scale_collisional'] == pytest.approx(lunar['derived_scale_collisional'], rel=5e-3)
    assert lunar['scale_capillary'] == pytest.approx(math.sqrt(CM1_GRAVITY / MOON_SURFACE_GRAVITY))
    # CM1's own limits, written out for the control.
    assert c.rain_limit_settings('earth', 1.2 * c.SEA_LEVEL_PA) == (
        dict(brkup_db=CM1_RAIN_BREAKUP_CUTOFF, brkup_gscale=0, brkup_rhor=1.0, brkup_rhoexp=0.0), None)
    with pytest.raises(ValueError):
        c.rain_limit_settings('mars', c.SEA_LEVEL_PA)


def test_the_namelist_takes_the_entries_after_brkup_db():
    text = ' &param3\n kdiff2 = 0.0,\n brkup_db = 300.0e-6,\n ndcnst = 100.0,\n /\n &param4\n stretch_x = 0,\n /\n'
    settings, _ = c.rain_limit_settings('lunar', 121239.27)
    out = c.with_rain_limits(text, settings)
    lines = out.splitlines()
    i = lines.index(' brkup_db = 0.0003,')
    assert lines[i + 1:i + 4] == [' brkup_gscale = 1,', ' brkup_rhor = 1.1965,',
                                  f" brkup_rhoexp = {settings['brkup_rhoexp']},"]
    assert lines[i + 4] == ' ndcnst = 100.0,' and out.endswith(' &param4\n stretch_x = 0,\n /\n')
    assert c.with_rain_limits(out, settings) == out
    assert c.namelist_value(c.with_rain_limits(out, c.rain_limit_settings('earth', 1.0)[0]), 'param3',
                            'brkup_gscale') == '0'


def test_a_case_writes_its_rain_limits_from_its_sounding(tmp_path):
    (tmp_path / 'input_sounding').write_text('   1212.3927     279.4226      9.58115\n')
    (tmp_path / 'namelist.template').write_text(' &param3\n brkup_db = 300.0e-6,\n /\n')
    assert c.apply_rain_limits(tmp_path, dict(), MOON_SURFACE_GRAVITY) is None
    record = c.apply_rain_limits(tmp_path, dict(rain_limits='lunar'), MOON_SURFACE_GRAVITY)
    text = (tmp_path / 'namelist.template').read_text()
    assert c.namelist_value(text, 'param3', 'brkup_rhor') == '1.1965'
    assert record['kind'] == 'lunar' and record['namelist']['brkup_rhor'] == 1.1965
    assert record['limits'] == c.rain_limits(MOON_SURFACE_GRAVITY, record['namelist'])


def test_the_rain_test_arms_differ_in_their_rain_limits_alone():
    control, lunar = c.CASES['box_0e_small_rain_control'], c.CASES['box_0e_small_lunar_rain']
    assert (control['rain_limits'], lunar['rain_limits']) == ('earth', 'lunar')
    same = lambda case: {k: v for k, v in case.items() if k not in ('rain_limits', 'purpose')}
    assert same(control) == same(lunar)
    assert control['build'] == 'moon_omp_rain_limits' and control['inputs_from'] == 'box_0e_small'
    # Apart from the build, the rain limits and the budgets they write, they are the small box itself.
    assert {k: v for k, v in same(control).items() if k not in ('build', 'inputs_from', 'namelist')} == \
           {k: v for k, v in c.CASES['box_0e_small'].items() if k not in ('build', 'purpose')}
    assert set(control['namelist']) == {'param9'}
