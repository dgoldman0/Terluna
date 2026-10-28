"""The rain-law patch: it applies to CM1's source after the rain limits, keeps Morrison's own law and self-collection
in every expression it touches at its default, reads its table as the fall-law product's mirror does, and the runner
writes the product's table and settings into a case."""
import hashlib
import json
import math

import numpy as np
import pytest

from atmosphere.rain import fall_law
from climate.crm import cm1_run as c

SOURCE = c.CM1_HOME / c.SOURCE['version'] / 'src'
needs_source = pytest.mark.skipif(not (SOURCE / 'morrison.F').exists(), reason=f'CM1 is not unpacked in {SOURCE}')
PRODUCT_PATH = c.HERE.parents[1] / c.RAIN_LAW_PRODUCT
PRODUCT = json.loads(PRODUCT_PATH.read_text())
FILES = ('input.F', 'cm1.F', 'param.F', 'morrison.F')


def patched(series):
    texts = {}
    for name, marker, edits in series:
        if name not in texts:
            texts[name] = (SOURCE / name).read_text(encoding='latin-1')
        texts[name] = c.apply_patch(texts[name], name, marker, edits)
    return texts


def test_the_patch_is_one_edit_set_per_file_each_adding_its_marker():
    assert [p[0] for p in c.RAIN_LAW_PATCHES] == list(FILES)
    for name, marker, edits in c.RAIN_LAW_PATCHES:
        assert marker == c.RAIN_LAW_MARKER and any(marker in new for _, new, _ in edits), name
    assert c.RAIN_LAW_MARKER not in c.RAIN_LIMITS_MARKER and c.RAIN_LIMITS_MARKER not in c.RAIN_LAW_MARKER
    assert c.BUILD_PATCHES['moon_omp_rain_law'] == c.RAIN_LIMITS_PATCHES + c.RAIN_LAW_PATCHES + [c.TWILIGHT_PATCH]
    assert c.BUILDS['moon_omp_rain_law'] == c.BUILDS['moon_omp']
    assert all(c.RAIN_LAW_PATCHES[0] not in patches for label, patches in c.BUILD_PATCHES.items()
               if label not in ('moon_omp_rain_law', 'moon_omp_rain_law_pbl', 'moon_omp_rain_law_pbl_wind',
                                'moon_omp_rain_law_pbl_wind_gust'))


def test_the_header_says_what_it_changes_its_default_and_its_source():
    header = c.RAIN_LAW_HEADER
    assert all(line.startswith('!') and len(line) <= 120 for line in header.splitlines())
    for phrase in ('What it changes.', 'Default. rain_fall_law = 0 and rain_selfcoll = 0', 'Source.',
                   'github.com/dgoldman0/Terluna', 'RAIN_LAW_PATCHES', 'atmosphere/rain/fall_law.py', 'UMR', 'UNR',
                   'EPSR', 'terluna_rain_law.txt', '5.78 m3/kg/s', 'VCAPR', 'Loftus and Wordsworth 2021'):
        assert phrase in ' '.join(header.split()), phrase
    for macro in ('TERLUNA_G', 'TERLUNA_CO2', 'TERLUNA_O2', '//', '/*', '#'):
        assert macro not in header, macro


@needs_source
def test_the_patch_applies_once_after_the_rain_limits():
    texts = patched([*c.PATCHES, *c.RAIN_LIMITS_PATCHES, *c.RAIN_LAW_PATCHES])
    for name, marker, edits in c.RAIN_LAW_PATCHES:
        assert c.apply_patch(texts[name], name, marker, edits) == texts[name], name     # idempotent
    morrison = texts['morrison.F']
    assert morrison.startswith(c.RAIN_LAW_HEADER + c.RAIN_LIMITS_HEADER)
    assert morrison.count('nragg(k) = -KSC_r*dum*qr3d(k)*nr3d(k)*rho(k)') == 2 and '-5.78*dum' not in morrison
    assert morrison.count('TERLUNA_RAIN(LAMR(K),1)') == 4 and morrison.count('TERLUNA_RAIN(LAMR(K),2)') == 4
    assert morrison.count('TERLUNA_RAIN(LAMR(K),3)') == 2 and morrison.count('TERLUNA_RAIN(DLAMR,') == 2
    assert 'use input, only : rain_fall_law,rain_selfcoll' in morrison
    assert 'brkup_rhor,brkup_rhoexp,rain_fall_law,rain_selfcoll,ndcnst,nt_c' in texts['input.F']
    assert texts['param.F'].count('MPI_BCAST(rain_') == 2
    assert '      rain_fall_law = 0' in texts['cm1.F'] and '      rain_selfcoll = 0.0\n' in texts['cm1.F']
    # Without the rain limits it has nothing to stand on.
    with pytest.raises(RuntimeError):
        patched([*c.PATCHES, *c.RAIN_LAW_PATCHES])


@needs_source
def test_the_default_keeps_morrisons_own_law():
    limits = patched([*c.PATCHES, *c.RAIN_LIMITS_PATCHES])['morrison.F']
    morrison = patched([*c.PATCHES, *c.RAIN_LIMITS_PATCHES, *c.RAIN_LAW_PATCHES])['morrison.F']
    # Every expression the patch touches stands as it was, in the branch its default takes.
    for expression, count in (('            UMR = ARN(K)*CONS4/(LAMR(K)**BR)\n', 4),
                              ('            UNR = ARN(K)*CONS6/LAMR(K)**BR\n', 4),
                              (c.RAIN_LAW_SEDIMENT, 1), (c.RAIN_LAW_EPSR, 2)):
        assert limits.count(expression) == count and morrison.count(expression) == count, expression
    assert morrison.count('IF (TERLUNA_RLAW) THEN') == 11
    assert '     LOGICAL, PRIVATE ::   TERLUNA_RLAW = .FALSE.' in morrison
    # The self-collection rate is CM1's unless the namelist gives one; the table is read only when asked for.
    block = c.RAIN_LAW_BLOCK
    assert block.index('KSC_r = 5.78') < block.index('IF (rain_selfcoll.GT.0.) KSC_r = rain_selfcoll')
    assert 'TERLUNA_RLAW = .TRUE.' in block[block.index('IF (rain_fall_law.EQ.1'):block.index('ELSE IF')]
    # Undoing the patch's edits gives back the rain-limits source: it adds nothing else.
    undone = morrison[len(c.RAIN_LAW_HEADER):]
    for old, new, _ in reversed(next(p for p in c.RAIN_LAW_PATCHES if p[0] == 'morrison.F')[2][1:]):
        undone = undone.replace(new, old)
    assert undone == limits


def fortran_rain(table, lam, m):
    """TERLUNA_RAIN as the patch writes it: Fortran's INT truncates toward zero, its arrays count from 1."""
    rows = [table['ln_mass_weighted'], table['ln_number_weighted'], table['ln_ventilation']][m - 1]
    x = (math.log(lam) - table['ln_lambda_first']) / table['ln_lambda_step']
    i = min(max(int(x), 0), table['points'] - 2)
    w = min(max(x - i, 0.0), 1.0)
    return math.exp((1.0 - w) * rows[i] + w * rows[i + 1])


def test_the_fortran_reading_matches_the_products_mirror():
    table = PRODUCT['table']
    for lam in (1.0, 99.0, 100.0, 145.3, 597.0, 1000.0, 2.0e4, 5.0e4, 1.0e5, 3.0e5):
        mirror = fall_law.interpolate(table, np.array([lam]))[:, 0]
        for m in (1, 2, 3):
            assert fortran_rain(table, lam, m) == pytest.approx(mirror[m - 1], rel=1e-12), (lam, m)


def test_the_lunar_settings_carry_the_product():
    settings, text, product = c.rain_law_settings('lunar')
    assert product == {c.RAIN_LAW_PRODUCT: hashlib.sha256(PRODUCT_PATH.read_bytes()).hexdigest()}
    assert settings == PRODUCT['namelist']
    lines = text.splitlines()
    t = PRODUCT['table']
    head = lines[0].split()
    assert int(head[0]) == t['points'] == len(lines) - 1
    assert float(head[1]) == pytest.approx(t['ln_lambda_first'], rel=1e-9)
    assert float(head[2]) == pytest.approx(t['ln_lambda_step'], rel=1e-9)
    rows = np.array([[float(x) for x in line.split()] for line in lines[1:]])
    for k, key in enumerate(('ln_mass_weighted', 'ln_number_weighted', 'ln_ventilation')):
        assert np.allclose(rows[:, k], t[key], rtol=1e-9, atol=1e-12), key
    assert c.rain_law_settings('earth') == (dict(rain_fall_law=0, rain_selfcoll=0.0), None, None)
    with pytest.raises(ValueError):
        c.rain_law_settings('mars')


def test_the_namelist_takes_the_entries_after_brkup_rhoexp():
    text = ' &param3\n brkup_db = 300.0e-6,\n ndcnst = 100.0,\n /\n &param4\n stretch_x = 0,\n /\n'
    text = c.with_rain_limits(text, c.rain_limit_settings('lunar', 121239.27)[0])
    out = c.with_rain_law(text, c.rain_law_settings('lunar')[0])
    lines = out.splitlines()
    i = next(n for n, line in enumerate(lines) if line.startswith(' brkup_rhoexp ='))
    assert lines[i + 1:i + 3] == [' rain_fall_law = 1,', f" rain_selfcoll = {PRODUCT['namelist']['rain_selfcoll']},"]
    assert lines[i + 3] == ' ndcnst = 100.0,' and out.endswith(' &param4\n stretch_x = 0,\n /\n')
    assert c.with_rain_law(out, c.rain_law_settings('lunar')[0]) == out
    assert c.namelist_value(c.with_rain_law(out, c.rain_law_settings('earth')[0]), 'param3', 'rain_fall_law') == '0'


def test_a_case_writes_its_rain_law_and_its_table(tmp_path):
    (tmp_path / 'namelist.template').write_text(' &param3\n brkup_db = 300.0e-6,\n brkup_gscale = 1,\n'
                                                ' brkup_rhor = 1.1965,\n brkup_rhoexp = 0.687,\n /\n')
    assert c.apply_rain_law(tmp_path, dict(rain_limits='lunar')) is None
    with pytest.raises(ValueError):
        c.apply_rain_law(tmp_path, dict(rain_law='lunar', rain_limits='earth'))
    record = c.apply_rain_law(tmp_path, dict(rain_law='lunar', rain_limits='lunar'))
    assert record['kind'] == 'lunar' and record['table'] == c.RAIN_LAW_TABLE
    assert (tmp_path / c.RAIN_LAW_TABLE).read_text() == c.rain_law_settings('lunar')[1]
    text = (tmp_path / 'namelist.template').read_text()
    assert c.namelist_value(text, 'param3', 'rain_fall_law') == '1'
    assert float(c.namelist_value(text, 'param3', 'rain_selfcoll')) == PRODUCT['namelist']['rain_selfcoll']
