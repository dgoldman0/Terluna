"""Cloud water for PlaSim's large-scale condensation (climate/gcm/condensate.py): the fall speed's constants come from the
lunar fall law; the Python copy keeps each layer's water and energy; the patch applies once to the shared rainmod.f90
and is undone exactly; mklsp and mkcwat, compiled with ExoPlaSim's flags against a stub of PlaSim's module, give the
Python copy's tendencies, cloud water and precipitation to 1e-10, and with the switch off mklsp prints every digit
the unpatched routine prints; the patched module sources compile with PlaSim's modules; and the runner records the
options, their namelist and their diagnostics."""
import inspect
import shutil
import subprocess

import numpy as np
import pytest

from climate.gcm import condensate as cnd
from climate.gcm import exoplasim_run as er

GCM_FLAGS = ('-O3', '-cpp', '-fcheck=all', '-ffpe-trap=invalid,zero,overflow', '-ffpe-summary=none',
             '-finit-real=zero', '-fdefault-real-8')
FORTRAN = pytest.mark.skipif(shutil.which('gfortran') is None, reason='gfortran is not installed')
C = dict(ga=1.6242188765608128, gascon=288.3767, acpd=288.3767 / 0.286, ADV=0.8546, ALV=2.5008e6, ALS=2.8345e6,
         rdbrv=288.3767 / 461.51, ra1=610.78, ra2=17.2693882, ra4=35.86, TMELT=273.16, deltsec=1800.0,
         deltsec2=3600.0)
num = lambda xs: ' '.join(repr(float(x)) for x in np.atleast_1d(xs))


def source(name):
    path = er.TWILIGHT_PACKAGE / 'exoplasim' / 'plasim' / 'src' / name
    return path if path.exists() else None


def apply(text, edits, marker):
    if marker in text:
        return text
    for old, new in edits:
        assert text.count(old) == 1, old
        text = text.replace(old, new)
    return text


def undo(text, edits, marker):
    if marker not in text:
        return text
    for old, new in reversed(edits):
        assert text.count(new) == 1, new[:80]
        text = text.replace(new, old)
    assert marker not in text
    return text


def routine(text, name):
    return text[text.index(f'      subroutine {name}\n'):text.index(f'      end subroutine {name}\n')] + \
        f'      end subroutine {name}\n'


# ------------------------------------------------------------------------------------------------------ columns

def columns(n=8, nlev=15, seed=5):
    """Columns (top first) on sigma levels: warm moist air near saturation low down, some layers supersaturated by
    their tendencies, cloud water in some layers, some of them now unsaturated, cold layers above, convective rain
    falling through some."""
    rng = np.random.default_rng(seed)
    half = np.r_[np.geomspace(0.02, 0.9, nlev - 3), 0.95, 0.98, 1.0]
    sigma = 0.5 * (np.r_[0.0, half[:-1]] + half)
    dsigma = np.diff(np.r_[0.0, half])
    ps = 121500.0 + rng.uniform(-500.0, 500.0, n)
    z = -C['gascon'] * 280.0 / C['ga'] * np.log(sigma)                   # m
    t = 302.0 - 1.2e-3 * z
    t = np.maximum(t, 200.0)[None, :] + rng.uniform(-1.0, 1.0, (n, nlev))
    qs, _ = cnd.saturation(t, ps[:, None], sigma[None, :], C)
    rh = rng.uniform(0.85, 1.03, (n, nlev))
    q = rh * qs
    s = dict(dq=q, dt=np.c_[t, np.full(n, 303.0)], dp=ps, sigma=sigma, dsigma=dsigma,
             dqdt=rng.uniform(-2e-8, 4e-8, (n, nlev)), dtdt=rng.uniform(-2e-5, 1e-5, (n, nlev)),
             dqcw=np.where(rng.uniform(size=(n, nlev)) < 0.5, rng.uniform(0.0, 1.5e-3, (n, nlev)), 0.0),
             dprll=np.zeros((n, nlev)), dprsll=np.zeros((n, nlev)), dcwcond=np.zeros((n, nlev)),
             dprcl=np.where(rng.uniform(size=(n, nlev)) < 0.2, rng.uniform(0.0, 3e-9, (n, nlev)), 0.0),
             dprscl=np.zeros((n, nlev)))
    s['dt'] = s['dt'][:, :nlev]
    return s


def copy(s):
    return {k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in s.items()}


def test_the_fall_speed_comes_from_the_lunar_fall_law():
    fit = cnd.fall_speed_fit()
    assert abs(fit["v0_m_s"] - cnd.PARAMETERS["vtr0_cw"]) < 0.005
    assert abs(fit["exponent"] - cnd.PARAMETERS["vtrexp_cw"]) < 0.005
    assert fit["largest_relative_error"] < 0.02
    assert cnd.PARAMETERS['fr0_cw'] == pytest.approx(3.0 / 86400.0)


def test_the_python_copy_keeps_each_layers_water_and_energy():
    s = columns()
    before = copy(s)
    cnd.mklsp(s, C)
    out = cnd.mkcwat(s, C)
    mass = s['dp'][:, None] * s['dsigma'][None, :] / C['ga']
    vapour = (s['dqdt'] - before['dqdt']) * C['deltsec']           # the leapfrog's mean change in a step
    cloud = s['dqcw'] - before['dqcw']
    fallen = s['dprll'] * C['deltsec'] / s['dsigma'][None, :]      # kg/kg handed to the precipitation in a step
    assert np.allclose((vapour + cloud + fallen) * mass, 0.0, atol=1e-12)
    heat = C['acpd'] * (1.0 + C['ADV'] * before['dq']) * (s['dtdt'] - before['dtdt'])
    assert np.allclose(heat + C['ALV'] * (s['dqdt'] - before['dqdt']), 0.0, atol=1e-12)
    assert (s['dqcw'] >= 0.0).all()
    assert (out['evaporation'] > 0.0).any() and (out['to_precipitation'] > 0.0).any()
    assert (s['dcwcond'] > 0.0).any()


def test_below_freezing_the_condensate_falls_at_once():
    s = columns()
    s['dt'][:] = 250.0
    qs, _ = cnd.saturation(s['dt'], s['dp'][:, None], s['sigma'][None, :], C)
    s['dq'] = qs * np.random.default_rng(3).uniform(0.9, 1.05, qs.shape)   # near its own saturation
    before = copy(s)
    cnd.mklsp(s, C)
    cnd.mkcwat(s, C)
    assert (s['dqcw'] == 0.0).all()
    assert (s['dcwcond'] == 0.0).all()
    assert (s['dprll'] > 0.0).any()


def test_the_patch_applies_once_and_undoes():
    path = source('rainmod.f90')
    if path is None:
        pytest.skip('no ExoPlaSim source at hand')
    text = path.read_text()
    plain = undo(text, cnd.CLOUD_WATER_EDITS, cnd.MARKER)
    patched = apply(plain, cnd.CLOUD_WATER_EDITS, cnd.MARKER)
    assert apply(patched, cnd.CLOUD_WATER_EDITS, cnd.MARKER) == patched
    assert undo(patched, cnd.CLOUD_WATER_EDITS, cnd.MARKER) == plain
    assert 'call mkcwat' in patched and 'subroutine mkcwat' in patched


# ------------------------------------------------------------------------------------------------------ Fortran

STUB = """
      module rainmod
      integer, parameter :: NHOR={nhor}, NLEV={nlev}, NLEP={nlev}+1, NLON=1, NLAT=1, NROOT=0
      real :: dt(NHOR,NLEP), dq(NHOR,NLEP), dtdt(NHOR,NLEP), dqdt(NHOR,NLEP), dp(NHOR)
      real :: sigma(NLEV), dsigma(NLEV)
      real :: deltsec, deltsec2, ga, gascon, acpd, adv, alv, als, rdbrv, ra1, ra2, ra4, tmelt
      real :: dprll(NHOR,NLEV), dprsll(NHOR,NLEV), dprcl(NHOR,NLEV), dprscl(NHOR,NLEV)
      real :: dqcw(NHOR,NLEV), dcwcond(NHOR,NLEV)
      real :: ncw_cm3, vtr0_cw, fr0_cw, vtrexp_cw
      real :: dgp3d(NHOR,NLEV,25)
      real :: dentropy(NHOR,40), dentro(NHOR), dentrot(NHOR,NLEV), dentroq(NHOR,NLEV), dentrop(NHOR)
      real :: dentro3d(NHOR,NLEV,40), denergy(NHOR,40), dener3d(NHOR,NLEV,40)
      integer :: ndiaggp=1, ndiaggp3d=25, nentropy=0, nentro3d=0, nenergy=0, nener3d=0
      integer :: nstorain=0, nprint=0, nprhor=1, mypid=0, nud=6, mtspd=48, ncwat=1
      end module rainmod
      subroutine get_random_root(z, n)
      integer :: n
      real :: z(n)
      z(:) = 0.5
      end subroutine get_random_root
      subroutine mpgagp(a, b, n)
      integer :: n
      real :: a(*), b(*)
      end subroutine mpgagp
"""
HARNESS = """
      program harness
      use rainmod
      implicit none
      integer :: jhor, jlev
      real :: zt0(NHOR,NLEP), zq0(NHOR,NLEP)
      read(*,*) deltsec, deltsec2, ga, gascon, acpd, adv, alv, als, rdbrv, ra1, ra2, ra4, tmelt
      read(*,*) ncw_cm3, vtr0_cw, fr0_cw, vtrexp_cw, ncwat
      read(*,*) sigma, dsigma
      do jhor=1,NHOR
       read(*,*) dp(jhor)
       read(*,*) (dt(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dq(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dtdt(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dqdt(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dqcw(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dprcl(jhor,jlev), jlev=1,NLEV)
      enddo
      dprscl(:,:)=0.
      zt0(:,:)=dtdt(:,:)
      zq0(:,:)=dqdt(:,:)
      call mklsp
      {second}
      do jhor=1,NHOR
       write(*,'(*(es25.16e3))') (dtdt(jhor,jlev)-zt0(jhor,jlev), jlev=1,NLEV)
       write(*,'(*(es25.16e3))') (dqdt(jhor,jlev)-zq0(jhor,jlev), jlev=1,NLEV)
       write(*,'(*(es25.16e3))') (dqcw(jhor,jlev), jlev=1,NLEV)
       write(*,'(*(es25.16e3))') (dprll(jhor,jlev), jlev=1,NLEV)
       write(*,'(*(es25.16e3))') (dgp3d(jhor,jlev,15), jlev=1,NLEV)
       write(*,'(*(es25.16e3))') (dgp3d(jhor,jlev,24), jlev=1,NLEV)
      enddo
      end program harness
"""
OUTPUTS = ('dtdt', 'dqdt', 'dqcw', 'dprll', 'dgp15', 'dgp24')


def run_fortran(folder, code, s, ncwat, with_mkcwat):
    n, nlev = s['dq'].shape
    folder.mkdir()
    files = {'stub.f90': STUB.format(nhor=n, nlev=nlev), 'rain.f90': code,
             'harness.f90': HARNESS.replace('{second}', 'call mkcwat' if with_mkcwat else '')}
    for name, body in files.items():
        (folder / name).write_text(body)
    subprocess.run(['gfortran', *GCM_FLAGS, '-o', 'harness', *files], cwd=folder, check=True, capture_output=True)
    p = cnd.PARAMETERS
    lines = [num([C['deltsec'], C['deltsec2'], C['ga'], C['gascon'], C['acpd'], C['ADV'], C['ALV'], C['ALS'],
                  C['rdbrv'], C['ra1'], C['ra2'], C['ra4'], C['TMELT']]),
             num([p['ncw_cm3'], p['vtr0_cw'], p['fr0_cw'], p['vtrexp_cw']]) + f' {ncwat}',
             num(np.r_[s['sigma'], s['dsigma']])]
    for j in range(n):
        lines += [num(s['dp'][j])] + [num(s[k][j]) for k in ('dt', 'dq', 'dtdt', 'dqdt', 'dqcw', 'dprcl')]
    text = subprocess.run([str(folder / 'harness')], input='\n'.join(lines) + '\n', capture_output=True, text=True,
                          check=True).stdout
    rows = [np.array(r.split(), float) for r in text.strip('\n').split('\n')]
    return {k: np.array([rows[j * len(OUTPUTS) + i] for j in range(n)]) for i, k in enumerate(OUTPUTS)}


@FORTRAN
def test_the_fortran_gives_the_python_copys_cloud_water(tmp_path):
    path = source('rainmod.f90')
    if path is None:
        pytest.skip('no ExoPlaSim source at hand')
    text = apply(path.read_text(), cnd.CLOUD_WATER_EDITS, cnd.MARKER)
    code = routine(text, 'mklsp') + routine(text, 'mkcwat')
    for seed in (5, 6, 7):
        s = columns(seed=seed)
        got = run_fortran(tmp_path / f'cw{seed}', code, s, 1, True)
        py = copy(s)
        cnd.mklsp(py, C)
        out = cnd.mkcwat(py, C)
        scale = lambda a: max(np.abs(a).max(), 1e-30)
        for key, want in (('dtdt', py['dtdt'] - s['dtdt']), ('dqdt', py['dqdt'] - s['dqdt']), ('dqcw', py['dqcw']),
                          ('dprll', py['dprll'])):
            assert np.abs(got[key] - want).max() <= 1e-10 * scale(want), key
        mass = s['dp'][:, None] * s['dsigma'][None, :] / C['ga'] / 1000.0
        assert np.abs(got['dgp24'] - out['evaporation'] * mass).max() <= 1e-10 * scale(out['evaporation'] * mass)


@FORTRAN
def test_with_the_switch_off_mklsp_prints_every_digit_the_old_one_prints(tmp_path):
    path = source('rainmod.f90')
    if path is None:
        pytest.skip('no ExoPlaSim source at hand')
    text = path.read_text()
    old = routine(undo(text, cnd.CLOUD_WATER_EDITS, cnd.MARKER), 'mklsp')
    new = routine(apply(text, cnd.CLOUD_WATER_EDITS, cnd.MARKER), 'mklsp')
    assert old != new
    s = columns(seed=9)
    a = run_fortran(tmp_path / 'old', old, s, 0, False)
    b = run_fortran(tmp_path / 'new', new, s, 0, False)
    for key in OUTPUTS:
        assert np.array_equal(a[key], b[key]), key


@FORTRAN
def test_the_patched_sources_compile_with_plasims_modules(tmp_path):
    path = source('rainmod.f90')
    if path is None:
        pytest.skip('no ExoPlaSim source at hand')
    src = path.parent
    (tmp_path / 'resmod.f90').write_text('      module resmod\n      parameter(NLAT_ATM = 32)\n'
                                         '      parameter(NLEV_ATM = 15)\n      parameter(NPRO_ATM = 1)\n'
                                         '      end module resmod\n')
    mod = (src / 'plasimmod.f90').read_text()
    for name, marker, edits in er.PLASIM_SWITCH_PATCHES:
        if name == 'plasimmod.f90':
            mod = apply(mod, edits, marker)
    (tmp_path / 'plasimmod.f90').write_text(mod)
    rain = apply((src / 'rainmod.f90').read_text(), cnd.CLOUD_WATER_EDITS, cnd.MARKER)
    (tmp_path / 'rainmod.f90').write_text(rain)
    restart = (src / 'restartmod.f90').read_text()
    for name, marker, edits in er.PLASIM_SWITCH_PATCHES:
        if name == 'restartmod.f90':
            restart = apply(restart, edits, marker)
    (tmp_path / 'restartmod.f90').write_text(restart)
    for name in ('resmod.f90', 'plasimmod.f90', 'rainmod.f90', 'restartmod.f90'):
        subprocess.run(['gfortran', *GCM_FLAGS, '-c', name], cwd=tmp_path, check=True, capture_output=True)


def test_every_switch_patch_applies_once_to_the_shared_source():
    for name, marker, edits in er.PLASIM_SWITCH_PATCHES + er.HUS_MEAN_PATCHES:
        path = source(name)
        if path is None:
            pytest.skip('no ExoPlaSim source at hand')
        text = path.read_text()
        plain = undo(text, edits, marker)
        once = apply(plain, edits, marker)
        assert apply(once, edits, marker) == once
        assert undo(once, edits, marker) == plain


def test_the_runner_records_the_options_their_namelist_and_diagnostics():
    for f in (er.run, er.configuration):
        names = inspect.signature(f).parameters
        assert 'large_scale_condensate' in names and 'humidity_mean' in names
    rec = cnd.setting()
    assert rec['scheme'] == 'cloud' and rec['parameters']['ncw_cm3'] == 100.0
    assert set(cnd.namelist_values()) == {'NCW_CM3', 'VTR0_CW', 'FR0_CW', 'VTREXP_CW'}
    assert sorted(cnd.DIAGNOSTICS) == [42, 43, 44, 45] and 21 + cnd.EXTRA_SLOTS == 25
    assert er.HUS_MEAN_CODE == 233
