"""PlaSim's rain evaporated from its drops (climate/gcm/rain_drops.py): the scheme's constants come from the lunar fall
law, the breakup product and Morrison's closures; the Python copy keeps the column's water and evaporates only where the
air is unsaturated; in saturated air the drops grow toward the lunar breakup size and stop there; the patch applies once
to the shared rainmod.f90 on lines no other patch touches and is undone exactly; mkrain with the patch, compiled with
ExoPlaSim's flags against a stub of PlaSim's module, gives the Python copy's tendencies and precipitation to 1e-10 under
both birth rules, and with nevapprec 1 prints every digit the unpatched mkrain prints; the patched rainmod.f90 compiles
with PlaSim's modules; and the runner records the rule and writes its namelist."""
import inspect
import json
import math
import shutil
import subprocess

import numpy as np
import pytest

from climate.gcm import exoplasim_run as er
from climate.gcm import levels as lv
from climate.gcm import rain_drops as rd
from climate.tests.test_condensate import FORTRAN, GCM_FLAGS, apply, num, routine, source, undo

C = dict(ga=1.6242188765608128, gascon=288.3767, acpd=288.3767 / 0.286, ADV=0.8546, ALV=2.5008e6, ALS=2.8345e6,
         rdbrv=288.3767 / 461.51, ra1=610.78, ra2=17.2693882, ra4=35.86, TMELT=273.16, deltsec=1800.0,
         deltsec2=3600.0, gamma=0.01)


def saturation(t, p):
    """PlaSim's saturation mixing ratio as mkrain forms it (before its correction factor)."""
    return C['rdbrv'] * C['ra1'] * np.exp(C['ra2'] * (t - C['TMELT']) / (t - C['ra4'])) / p


def columns(n=8, seed=5, ground_freezing=False):
    """Lunar columns (top first) on the fifteen levels: a moist adiabat's slow cooling from about 300 K, freezing near
    24 km, the air 55-100% humid and some layers a little supersaturated, snow forming above the freezing level and rain
    below it in most columns, none in some."""
    rng = np.random.default_rng(seed)
    half = np.asarray(lv.half_levels('l15'), float)
    nlev = half.size
    edges = np.r_[0.0, half]
    sigma = 0.5 * (edges[1:] + edges[:-1])
    dsigma = np.diff(edges)
    ps = 121500.0 + rng.uniform(-3000.0, 1000.0, n)
    z = -C['gascon'] * 285.0 / C['ga'] * np.log(sigma)
    t = (300.5 - 1.15e-3 * z)[None, :] + rng.uniform(-1.0, 1.0, (n, nlev))
    t = np.maximum(t, 190.0)
    p = ps[:, None] * sigma[None, :]
    rh = np.where(rng.uniform(size=(n, nlev)) < 0.1, rng.uniform(1.0, 1.02, (n, nlev)),
                  rng.uniform(0.55, 0.97, (n, nlev)))
    q = rh * saturation(t, p)
    cold = t <= C['TMELT']
    raining = rng.uniform(size=n) < 0.75
    form = np.where(rng.uniform(size=(n, nlev)) < 0.6, rng.uniform(0.0, 4e-10, (n, nlev)), 0.0) * raining[:, None]
    conv = form * rng.uniform(0.3, 1.0, (n, nlev))
    large = form - conv
    ground = np.where(ground_freezing & (rng.uniform(size=n) < 0.5), 270.0, t[:, -1] + 1.0)
    return dict(dt=np.c_[t, ground], dq=q, dp=ps, sigma=sigma, dsigma=dsigma,
                dtdt=rng.uniform(-2e-5, 1e-5, (n, nlev)), dqdt=rng.uniform(-2e-8, 2e-8, (n, nlev)),
                dprcl=conv, dprscl=np.where(cold, conv, 0.0), dprll=large, dprsll=np.where(cold, large, 0.0))


def copy(s):
    return {k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in s.items()}


# ------------------------------------------------------------------------------------------------------ Python

def test_the_constants_come_from_the_fall_law_the_breakup_and_morrisons_closures():
    law = json.loads(rd.FALL_LAW.read_text())
    breakup = json.loads(rd.BREAKUP.read_text())
    assert rd.PARAMETERS['rdrop_ksc'] == law['self_collection']['lunar_rate_m3_kg_s'] == 2.919
    assert rd.PARAMETERS['rdrop_rhoexp'] == breakup['lunar_limits']['namelist']['brkup_rhoexp'] == 0.687
    assert np.array_equal(rd.PRODUCTS['table'][2], law['table']['ln_ventilation'])
    text = rd.declarations()
    for value in (rd.RHOSU, 0.78, 0.308, 997.0, 8.794e-5, 1.81, 1.496e-6, 120.0, 0.0003, 2300.0, 0.0028, 2e-05, 9.1):
        assert repr(value) in text, value
    lim = rd.limits(C['ga'])
    assert 1.0 / math.exp(lim['lnlam_min']) == pytest.approx(6.88e-3, rel=1e-3)        # the lunar cap, 6.88 mm
    assert lim['vcap'] == pytest.approx(5.805, rel=1e-3)                               # 9.1 (g/9.81)^(1/4)
    inv = rd.inverse(C['ga'])
    for d in (0.1e-3, 0.465e-3, 1.68e-3, 4e-3):                                       # the inverse reads back
        r = math.log(rd.number_per_flux(d, inv))
        assert math.exp(-rd.slope_of(r, inv)) == pytest.approx(d, rel=1e-6)


def test_the_python_copy_keeps_the_columns_water_and_evaporates_only_unsaturated_air():
    for seed in (5, 6):
        s = columns(seed=seed)
        before = copy(s)
        out = rd.mkrain(s, C)
        mass = s['dp'] / C['ga']
        formed = (s['dprcl'] + s['dprll']).sum(axis=1)                                 # sigma-weighted, per second
        ground = (out['dprc'] + out['dprl']) * 1000.0 / mass
        gone = ((out['rain_evaporation'] + out['snow_sublimation']) * s['dsigma'][None, :]).sum(axis=1)
        assert np.allclose(formed, ground + gone, rtol=1e-12, atol=1e-24)
        zq = before['dq'] + before['dqdt'] * C['deltsec2']
        zt = before['dt'][:, :-1] + before['dtdt'] * C['deltsec2']
        qs = saturation(zt, s['dp'][:, None] * s['sigma'][None, :])
        qs = np.minimum(qs, C['rdbrv'])
        qs = qs / (1.0 - (1.0 / C['rdbrv'] - 1.0) * qs)
        assert (out['rain_evaporation'][zq >= qs] == 0.0).all()
        assert (out['rain_evaporation'] >= 0.0).all() and (out['rain_evaporation'] > 0.0).any()
        rained = out['dprc'] + out['dprl'] - out['dprs'] > 0.0
        assert (out['number'][rained] > 0.0).all()


def test_in_saturated_air_the_drops_grow_toward_the_breakup_size_and_stop_there():
    s = columns(n=1, seed=8)
    s['dq'] = saturation(s['dt'][:, :-1], s['dp'][:, None] * s['sigma'][None, :]) * 1.01   # no evaporation
    s['dt'][:] = 290.0                                                              # rain all the way down
    s['dprcl'][:] = 0.0
    s['dprcl'][0, 6] = 3e-9                                                         # rain born high, 15-23 km up
    s['dprscl'][:] = s['dprll'][:] = s['dprsll'][:] = 0.0
    inv = rd.inverse(C['ga'])
    out = rd.mkrain(copy(s), C)
    flux = out['dprc'][0] * 1000.0
    size = math.exp(-rd.slope_of(math.log(out['number'][0] / flux), inv))
    scol = math.sqrt(9.81 / C['ga']) * (s['dp'][0] / 101325.0) ** rd.PARAMETERS['rdrop_rhoexp']
    equilibrium = (rd.CM1_RAIN_BREAKUP_CUTOFF + math.log(2.0) / rd.CM1_RAIN_BREAKUP_RATE) * scol
    assert rd.PARAMETERS['rdrop_melt'] < size <= equilibrium * 1.001
    assert size > 0.5 * equilibrium


# ------------------------------------------------------------------------------------------------------ Fortran

def test_the_patch_applies_once_and_undoes_on_lines_no_other_patch_touches():
    path = source('rainmod.f90')
    if path is None:
        pytest.skip('no ExoPlaSim source at hand')
    text = path.read_text()
    plain = undo(text, rd.RAIN_DROP_EDITS, rd.MARKER)
    patched = apply(plain, rd.RAIN_DROP_EDITS, rd.MARKER)
    assert apply(patched, rd.RAIN_DROP_EDITS, rd.MARKER) == patched
    assert undo(patched, rd.RAIN_DROP_EDITS, rd.MARKER) == plain
    import climate.gcm.condensate as cnd
    import climate.gcm.massflux_fortran as mff
    others = [e for _, _, edits in er.PATCHES + er.SBM_DEPTH_PATCHES + er.SBM_CAP_PATCHES + mff.MASSFLUX_PATCHES
              + cnd.CLOUD_WATER_PATCHES for e in edits]
    for old, _ in rd.RAIN_DROP_EDITS:
        assert all(old not in o and old not in n for o, n in others), old


STUB = """
      module rainmod
      integer, parameter :: NHOR={nhor}, NLEV={nlev}, NLEP={nlev}+1, NLEM={nlev}-1, NLON=1, NLAT=1, NROOT=0
      real :: dt(NHOR,NLEP), dq(NHOR,NLEP), dtdt(NHOR,NLEP), dqdt(NHOR,NLEP), dp(NHOR)
      real :: sigma(NLEV), dsigma(NLEV)
      real :: deltsec, deltsec2, ga, gascon, acpd, adv, alv, als, rdbrv, ra1, ra2, ra4, tmelt, gamma
      real :: dprll(NHOR,NLEV), dprsll(NHOR,NLEV), dprcl(NHOR,NLEV), dprscl(NHOR,NLEV)
      real :: dprc(NHOR), dprl(NHOR), dprs(NHOR)
      real :: dgp3d(NHOR,NLEV,25)
      real :: dentropy(NHOR,40), dentro(NHOR), dentrot(NHOR,NLEV), dentroq(NHOR,NLEV), dentrop(NHOR)
      real :: dentro3d(NHOR,NLEV,40), denergy(NHOR,40), dener3d(NHOR,NLEV,40)
      integer :: ndiaggp=1, nentropy=0, nentro3d=0, nenergy=0, nener3d=0, nprc=2
      integer :: nevapprec= 1
{declarations}
      end module rainmod
"""
HARNESS = """
      program harness
      use rainmod
      implicit none
      integer :: jhor, jlev
      real :: zt0(NHOR,NLEP), zq0(NHOR,NLEP)
      read(*,*) deltsec, deltsec2, ga, gascon, acpd, adv, alv, als, rdbrv, ra1, ra2, ra4, tmelt, gamma
      read(*,*) nevapprec{switch}
      read(*,*) sigma, dsigma
      do jhor=1,NHOR
       read(*,*) dp(jhor)
       read(*,*) (dt(jhor,jlev), jlev=1,NLEP)
       read(*,*) (dq(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dtdt(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dqdt(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dprcl(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dprscl(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dprll(jhor,jlev), jlev=1,NLEV)
       read(*,*) (dprsll(jhor,jlev), jlev=1,NLEV)
      enddo
      zt0(:,:)=dtdt(:,:)
      zq0(:,:)=dqdt(:,:)
      call mkrain
      do jhor=1,NHOR
       write(*,'(*(es25.16e3))') (dtdt(jhor,jlev)-zt0(jhor,jlev), jlev=1,NLEV)
       write(*,'(*(es25.16e3))') (dqdt(jhor,jlev)-zq0(jhor,jlev), jlev=1,NLEV)
       write(*,'(*(es25.16e3))') dprc(jhor), dprl(jhor), dprs(jhor)
      enddo
      end program harness
"""
SETTINGS = ', ' + ', '.join(rd.NAMELIST_NAMES)


def mkrain_code(text, patched):
    end = 'end subroutine rdropevap' if patched else 'end subroutine mkrain'
    start = text.index('      subroutine mkrain\n')
    return text[start:text.index(f'      {end}\n') + len(f'      {end}\n')]


def run_fortran(folder, code, s, nevapprec, patched, prm=rd.PARAMETERS):
    n, nlev = s['dq'].shape
    folder.mkdir()
    files = {'stub.f90': STUB.format(nhor=n, nlev=nlev, declarations=rd.declarations(prm) if patched else ''),
             'rain.f90': code, 'harness.f90': HARNESS.replace('{switch}', SETTINGS if patched else '')}
    for name, body in files.items():
        (folder / name).write_text(body)
    subprocess.run(['gfortran', *GCM_FLAGS, '-o', 'harness', *files], cwd=folder, check=True, capture_output=True)
    lines = [num([C[k] for k in ('deltsec', 'deltsec2', 'ga', 'gascon', 'acpd', 'ADV', 'ALV', 'ALS', 'rdbrv', 'ra1',
                                 'ra2', 'ra4', 'TMELT', 'gamma')]),
             f'{nevapprec} ' + (' '.join(repr(prm[k]) for k in rd.NAMELIST_NAMES) if patched else ''),
             num(np.r_[s['sigma'], s['dsigma']])]
    for j in range(n):
        lines += [num(s['dp'][j])] + [num(s[k][j]) for k in ('dt', 'dq', 'dtdt', 'dqdt', 'dprcl', 'dprscl', 'dprll',
                                                               'dprsll')]
    text = subprocess.run([str(folder / 'harness')], input='\n'.join(lines) + '\n', capture_output=True, text=True,
                          check=True).stdout
    rows = [np.array(r.split(), float) for r in text.strip('\n').split('\n')]
    return dict(dtdt=np.array([rows[3 * j] for j in range(n)]), dqdt=np.array([rows[3 * j + 1] for j in range(n)]),
                ground=np.array([rows[3 * j + 2] for j in range(n)]))


@FORTRAN
def test_the_fortran_gives_the_python_copys_rain_under_both_birth_rules(tmp_path):
    path = source('rainmod.f90')
    if path is None:
        pytest.skip('no ExoPlaSim source at hand')
    code = mkrain_code(apply(path.read_text(), rd.RAIN_DROP_EDITS, rd.MARKER), True)
    scale = lambda a: max(np.abs(a).max(), 1e-30)
    for rule in rd.BIRTH_RULES:
        prm = dict(rd.PARAMETERS, **rd.BIRTH_RULES[rule])
        for seed, freezing in ((5, False), (6, True)):
            s = columns(seed=seed, ground_freezing=freezing)
            got = run_fortran(tmp_path / f'{rule}{seed}', code, s, 2, True, prm)
            py = copy(s)
            out = rd.mkrain(py, C, prm)
            for key, want in (('dtdt', py['dtdt'] - s['dtdt']), ('dqdt', py['dqdt'] - s['dqdt']),
                              ('ground', np.c_[out['dprc'], out['dprl'], out['dprs']])):
                assert np.abs(got[key] - want).max() <= 1e-10 * scale(want), (rule, seed, key)


@FORTRAN
def test_with_nevapprec_1_the_patched_mkrain_prints_every_digit_the_old_one_prints(tmp_path):
    path = source('rainmod.f90')
    if path is None:
        pytest.skip('no ExoPlaSim source at hand')
    text = path.read_text()
    old = mkrain_code(undo(text, rd.RAIN_DROP_EDITS, rd.MARKER), False)
    new = mkrain_code(apply(text, rd.RAIN_DROP_EDITS, rd.MARKER), True)
    assert old != new
    for seed in (9, 10):
        s = columns(seed=seed, ground_freezing=seed == 10)
        a = run_fortran(tmp_path / f'old{seed}', old, s, 1, False)
        b = run_fortran(tmp_path / f'new{seed}', new, s, 1, True)
        for key in a:
            assert np.array_equal(a[key], b[key]), key


@FORTRAN
def test_the_patched_rainmod_compiles_with_plasims_modules(tmp_path):
    path = source('rainmod.f90')
    if path is None:
        pytest.skip('no ExoPlaSim source at hand')
    src = path.parent
    (tmp_path / 'resmod.f90').write_text('      module resmod\n      parameter(NLAT_ATM = 32)\n'
                                         '      parameter(NLEV_ATM = 15)\n      parameter(NPRO_ATM = 1)\n'
                                         '      end module resmod\n')
    for name in ('plasimmod.f90', 'restartmod.f90'):                      # the shared copies, as patched
        shutil.copy(src / name, tmp_path / name)
    (tmp_path / 'rainmod.f90').write_text(apply(path.read_text(), rd.RAIN_DROP_EDITS, rd.MARKER))
    for name in ('resmod.f90', 'plasimmod.f90', 'rainmod.f90', 'restartmod.f90'):
        subprocess.run(['gfortran', *GCM_FLAGS, '-c', name], cwd=tmp_path, check=True, capture_output=True)


def test_the_runner_records_the_rule_and_writes_its_namelist():
    for f in (er.run, er.configuration):
        assert inspect.signature(f).parameters['rain_drops'].default is None
    for rule, births in rd.BIRTH_RULES.items():
        rec = rd.setting(rule)
        assert rec['nevapprec'] == 2 and set(rec['namelist']) == set(rd.NAMELIST_NAMES)
        assert all(rec['namelist'][k] == v for k, v in births.items())
    with pytest.raises(ValueError):
        rd.setting('gamma')
    text = (er.HERE / 'exoplasim_run.py').read_text()
    assert "model._edit_namelist('rainmod_namelist', 'NEVAPPREC', '2')" in text
    assert 'rdr.RAIN_DROP_PATCHES if rain_drops' in text
    with pytest.raises(ValueError):
        er.configuration('A28R', 4, rain_drops='one_size', convection=dict(scheme='massflux', start=None))


@FORTRAN
def test_in_dry_air_heavy_rain_and_extreme_sizes_the_fortran_stays_finite_and_follows_the_copy(tmp_path):
    """Rain falling heavily into air 5-30% humid, born at the floor's 20 um, at the sea box's 0.465 mm or at 6.5 mm:
    the drops evaporate within a step's half, break up hard or saturate their layer, and the compiled routine, with
    ExoPlaSim's traps on invalid arithmetic, division by zero and overflow, still gives the copy's numbers."""
    path = source('rainmod.f90')
    if path is None:
        pytest.skip('no ExoPlaSim source at hand')
    code = mkrain_code(apply(path.read_text(), rd.RAIN_DROP_EDITS, rd.MARKER), True)
    scale = lambda a: max(np.abs(a).max(), 1e-30)
    for size in (2e-5, rd.CM1_SEA_BIRTH_SIZE_M, 6.5e-3):
        for join in (0, 1):
            prm = dict(rd.PARAMETERS, rdrop_melt=size, rdrop_warm=size, nrdrop_join=join)
            s = columns(seed=12, ground_freezing=True)
            rng = np.random.default_rng(13)
            s['dq'] = s['dq'] * np.where(rng.uniform(size=s['dq'].shape) < 0.5, rng.uniform(0.05, 0.3, s['dq'].shape),
                                         1.0)
            for key in ('dprcl', 'dprscl', 'dprll', 'dprsll'):
                s[key] = s[key] * 40.0
            got = run_fortran(tmp_path / f'dry{size}_{join}', code, s, 2, True, prm)
            py = copy(s)
            out = rd.mkrain(py, C, prm)
            for key, want in (('dtdt', py['dtdt'] - s['dtdt']), ('dqdt', py['dqdt'] - s['dqdt']),
                              ('ground', np.c_[out['dprc'], out['dprl'], out['dprs']])):
                assert np.isfinite(got[key]).all() and np.isfinite(want).all()
                assert np.abs(got[key] - want).max() <= 1e-10 * scale(want), (size, join, key)
