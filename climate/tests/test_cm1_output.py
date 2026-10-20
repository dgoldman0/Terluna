"""CM1's output conventions: its radiative tendencies stored at constant volume, its temperature, its equation set."""
import numpy as np
import pytest

from climate.crm import cm1_output as co
from climate.crm import rrtmg_edge


def snapshot():
    """Two levels of a lunar column's water, kg/kg, as a snapshot holds them: moist air low down, drier cloud above."""
    return dict(qv=np.array([0.016, 0.002]), qc=np.array([0.0, 1e-4]), qr=np.array([2e-5, 0.0]),
                qi=np.array([0.0, 3e-5]), qs=np.array([0.0, 1e-5]), qg=np.array([0.0, 0.0]))


def test_rrtmg_rate_undoes_what_cm1_stores_under_the_bryan_fritsch_set():
    d = snapshot()
    hr = np.array([-0.30, -0.25]) / 86400.0
    # radiation_driver.F: rcvm = cp/(cv + cvv*max(0,qv) + cpl*max(0,qc+qr) + cpi*max(0,qi+qs+qg)), lwten = hr*rcvm
    cv, cvv = 1005.7 - 287.04, 1870.0 - 461.5
    rcvm = 1005.7 / (cv + cvv * d['qv'] + 4190.0 * (d['qc'] + d['qr']) + 2106.0 * (d['qi'] + d['qs'] + d['qg']))
    stored = hr * rcvm
    assert np.allclose(co.rrtmg_rate(stored, d, 2), hr, rtol=1e-12)
    assert np.allclose(co.rrtmg_rate(stored, d, 1), stored)
    factor = co.cp_over_cvm(d)
    assert 1.35 < factor[0] < 1.37 < factor[1] < 1.40


def test_a_mixing_ratio_the_snapshot_lacks_counts_as_none():
    d = dict(qv=np.array([0.01, -1e-6]))
    assert np.allclose(co.cp_over_cvm(d), 1005.7 / (1005.7 - 287.04 + (1870.0 - 461.5) * np.array([0.01, 0.0])))


def test_the_applied_rate_is_what_the_potential_temperature_takes():
    stored = np.array([-0.4, 0.2])
    pi = np.array([1.05, 0.97])
    thraten = stored * (1005.7 - 287.04) / 1005.7 / pi                  # radiation_driver.F under eqtset 2
    assert np.allclose(co.applied_rate(stored, 2), pi * thraten)
    assert np.allclose(co.applied_rate(stored, 1), stored)              # eqtset 1: (swten + lwten)/pi
    hr = np.array([-0.3, -0.1])
    dry = dict(qv=np.zeros(2))
    assert np.allclose(co.applied_rate(hr * co.cp_over_cvm(dry), 2), hr)   # dry air takes RRTMG's own heating


def test_the_equation_set_reads_from_the_namelist(tmp_path):
    (tmp_path / 'namelist.input').write_text(' &param2\n ptype     =  5,\n eqtset    =  2,\n /\n')
    assert co.eqtset(tmp_path) == 2
    (tmp_path / 'namelist.input').write_text(' &param2\n ptype     =  5,\n /\n')
    with pytest.raises(KeyError):
        co.eqtset(tmp_path)


def test_the_temperature_takes_cm1s_own_exner_function():
    d = dict(th=np.array([300.0]), prs=np.array([121400.0]))
    assert co.temperature(d)[0] == pytest.approx(300.0 * 1.214 ** (287.04 / 1005.7), rel=1e-12)


def test_the_conventions_are_the_lunar_builds():
    """CM1's constants, its conversion after the RRTMG calls and the ground's temperature it hands RRTMG."""
    driver, constants = rrtmg_edge.BUILD / 'radiation_driver.F', rrtmg_edge.BUILD / 'constants.F'
    if not (driver.exists() and constants.exists()):
        pytest.skip('the lunar CM1 build is not at hand')
    text = driver.read_text()
    for line in ('rcvm = cp/( cv+cvv*max(0.0,qa(i,j,k,nqv))', '+cpl*max(0.0,(dum10(i,j,k)+dum11(i,j,k)))',
                 '+cpi*max(0.0,(dum9(i,j,k)+dum7(i,j,k)+dum8(i,j,k))) )', 'lwten(i,j,k) = lwten(i,j,k)*rcvm',
                 'swten(i,j,k) = swten(i,j,k)*rcvm',
                 'thraten(i,j,k) = (swten(i,j,k)+lwten(i,j,k))*cvdcp/(pi0(i,j,k)+ppi(i,j,k))',
                 'dum2(i,j,1) = tsk(i,j)'):
        assert line in text, line
    source = constants.read_text()
    defaults = source.split('!  cm1 defaults:')[1].split('ENDIF')[0]
    for name, value in (('rd', co.RD), ('cp', co.CP), ('cpv', co.CPV)):
        assert f'{name.ljust(6)} = {value}' in defaults
    for name, value in (('rv', co.RV), ('cpl', co.CPL), ('cpi', co.CPI)):
        assert f'{name.ljust(6)} = {value}' in source
    assert 'cv     = cp-rd' in source and 'cvv    = cpv-rv' in source
