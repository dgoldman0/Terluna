import shutil
import subprocess

import numpy as np
import pytest

from climate.gcm import plasim_ice_optics as pio


def test_the_radius_reads_the_table_as_reicalc_does_and_follows_its_trend_below():
    assert np.allclose(pio.ice_radius(np.arange(180.0, 275.0)), pio.RETAB)
    assert np.isclose(pio.ice_radius(189.3), 0.7 * pio.RETAB[9] + 0.3 * pio.RETAB[10])
    # Below 180 K each kelvin takes the table's first step again, so the radius meets the table at 180 K.
    step = pio.RETAB[1] / pio.RETAB[0]
    assert np.isclose(pio.ice_radius(179.999999), pio.RETAB[0], rtol=1e-6)
    assert np.isclose(pio.ice_radius(170.0) / pio.ice_radius(171.0), 1.0 / step)


def test_the_near_infrared_weights_are_the_sunlights_shares():
    from atmosphere.radiative_convective import fetch_inputs
    try:
        fetch_inputs.path('TSIS1_HSRS_stride100.csv')
    except Exception:
        pytest.skip('TSIS-1 spectrum not fetched')
    weights = pio.nir_weights()
    assert np.isclose(sum(weights), 1.0)
    assert tuple(round(w, 4) for w in weights) == pio.NIR_WEIGHTS


def test_ice_conserves_visible_light_and_reflects_more_of_it_than_plasims_cloud():
    for t, wp in ((189.0, 1.5), (217.0, 4.7)):
        for mu in (0.1, 0.5, 1.0):
            ice = pio.ice_sw(wp, t, mu)
            r1, _, r1s, t1s = pio.pifm(wp * pio.band_optics(float(pio.ice_radius(t)))[0][0],
                                       1.0 - pio.EC_C[0], pio.EC_E[0] + pio.EC_F[0] * float(pio.ice_radius(t)), mu)
            assert np.isclose(r1s + t1s, 1.0, atol=1e-4)                     # visible: scattering only
            assert ice['r2'] + ice['t2'] < 1.0 and ice['r2s'] + ice['t2s'] < 1.0
    # The coldest high cloud: about PlaSim's infrared, several times its sunlight reflection.
    liquid = pio.disk_mean(lambda m: pio.liquid_sw(1.5, m)['r1'])
    ice = pio.disk_mean(lambda m: pio.ice_sw(1.5, 189.0, m)['r1'])
    assert ice > 5.0 * liquid
    assert abs(pio.emissivity(1.5, pio.ice_kappa(189.0)) - pio.emissivity(1.5, pio.LIQUID_KAPPA_M2_G)) < 0.03
    # At 217 K the crystals are larger: less than half PlaSim's infrared for the same water.
    assert pio.emissivity(4.7, pio.ice_kappa(217.0)) < 0.5 * pio.emissivity(4.7, pio.LIQUID_KAPPA_M2_G)


STAND_IN = '''
      module radmod
      implicit none
      integer, parameter :: NHOR = 5, NLEV = 1
      real :: dcc(NHOR,NLEV), dt(NHOR,NLEV), dql(NHOR,NLEV), dp(NHOR), dsigma(NLEV), ga, ticeopt
      end module radmod
'''
DRIVER = '''
      program drive
      use radmod
      implicit none
      real :: zmu0(NHOR), zcs(NHOR), z(NHOR,7)
      logical :: losun(NHOR)
      integer :: j
      ga = {ga}
      dsigma(1) = {dsigma}
      ticeopt = {ticeopt}
      dt(:,1) = (/{t}/)
      dql(:,1) = (/{q}/)
      dp(:) = (/{p}/)
      dcc(:,1) = (/{cc}/)
      zmu0(:) = (/{mu}/)
      zcs(:) = (/{cs}/)
      losun(:) = (/{sun}/)
      z(:,:) = -1.
      call ice_cloud_sw(1,zmu0,zcs,losun,z(:,1),z(:,2),z(:,3),z(:,4),z(:,5),z(:,6))
      call ice_cloud_lw(1,z(:,7))
      do j = 1, NHOR
       print '(7es26.17)', z(j,:)
      enddo
      end program drive
'''


def test_the_fortran_patch_matches_its_python_copy(tmp_path):
    if shutil.which('gfortran') is None:
        pytest.skip('no gfortran')
    from climate.gcm import exoplasim_run as er
    ga, dsigma, ps = 1.6242188765608128, 0.085, 117000.0
    t = [189.3, 217.0, 165.0, 240.0, 200.0]                # the fourth is too warm for ice, the fifth unlit
    water = [1.5, 4.7, 0.9, 3.0, 2.0]
    q = [w / (1000.0 * ps / ga * dsigma) for w in water]
    mu, cs, cc = [0.5, 0.2, 0.9, 0.5, 0.5], [1.0, 0.6, 1.0, 1.0, 1.0], [0.3, 0.1, 0.05, 0.2, 0.4]
    sun = ['.true.'] * 4 + ['.false.']
    fmt = lambda xs: ','.join(repr(float(x)) for x in xs)
    code = STAND_IN + er.ICE_ROUTINES + DRIVER.format(ga=repr(ga), dsigma=repr(dsigma), ticeopt=repr(pio.ICE_BELOW_K),
                                                      t=fmt(t), q=fmt(q), p=fmt([ps] * 5), cc=fmt(cc), mu=fmt(mu),
                                                      cs=fmt(cs), sun=','.join(sun))
    (tmp_path / 'ice.f90').write_text(code)
    subprocess.run(['gfortran', '-fdefault-real-8', '-O0', '-o', 'ice', 'ice.f90'], cwd=tmp_path, check=True)
    out = np.loadtxt(subprocess.run(['./ice'], cwd=tmp_path, check=True, capture_output=True, text=True).stdout.split('\n'))
    for j in range(3):                                      # the cold, lit columns
        ice = pio.ice_sw(water[j], t[j], mu[j])
        want = [cs[j] * ice['r1'] + (1 - cs[j]) * ice['r1s'], ice['r1s'], cs[j] * ice['r2'] + (1 - cs[j]) * ice['r2s'],
                cs[j] * ice['t2'] + (1 - cs[j]) * ice['t2s'], ice['r2s'], ice['t2s']]
        assert np.allclose(out[j, :6], want, rtol=1e-9, atol=1e-12)
        assert np.isclose(out[j, 6], 1.0 - cc[j] * pio.emissivity(water[j], pio.ice_kappa(t[j])), rtol=1e-12)
    assert np.all(out[3] == -1.0)                           # too warm: PlaSim's own optics stay
    assert np.all(out[4, :6] == -1.0)                       # unlit: the sunlight is left alone
    assert out[4, 6] != -1.0                                # the infrared is not


def test_the_shield_weights_are_the_design_shields_sunlight():
    from atmosphere.radiative_convective import fetch_inputs
    try:
        fetch_inputs.path('TSIS1_HSRS_stride100.csv')
    except Exception:
        pytest.skip('TSIS-1 spectrum not fetched')
    weights = pio.shielded_nir_weights()
    assert np.isclose(sum(weights), 1.0)
    assert tuple(round(w, 4) for w in weights) == pio.SHIELD_NIR_WEIGHTS
